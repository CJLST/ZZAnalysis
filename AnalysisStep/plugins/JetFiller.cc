#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/EventSetup.h>

#include <DataFormats/PatCandidates/interface/Jet.h>

#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>

#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <JetMETCorrections/Objects/interface/JetCorrectionsRecord.h>

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;


class JetFiller : public edm::EDProducer {
 public:
  /// Constructor
  explicit JetFiller(const edm::ParameterSet&);
    
  /// Destructor
  ~JetFiller(){};  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  edm::EDGetTokenT<edm::View<pat::Jet> > jetToken;
  const StringCutObjectSelector<pat::Jet, true> cut;
  const std::string bTaggerName;
  const std::string jecType;
  const CutSet<pat::Jet> flags;
  edm::EDGetTokenT<edm::ValueMap<float> > qgToken;
  edm::EDGetTokenT<edm::ValueMap<float> > axis2Token;
  edm::EDGetTokenT<edm::ValueMap<int> > multToken;
  edm::EDGetTokenT<edm::ValueMap<float> > ptDToken;

};


JetFiller::JetFiller(const edm::ParameterSet& iConfig) :
  jetToken(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("src"))),
  cut(iConfig.getParameter<std::string>("cut")),
  bTaggerName(iConfig.getParameter<std::string>("bTaggerName")),
  jecType(iConfig.getParameter<std::string>("jecType")),
  flags(iConfig.getParameter<edm::ParameterSet>("flags"))
{

  qgToken = consumes<edm::ValueMap<float> >(edm::InputTag("QGTagger", "qgLikelihood"));
  axis2Token = consumes<edm::ValueMap<float> >(edm::InputTag("QGTagger", "axis2"));
  multToken = consumes<edm::ValueMap<int> >(edm::InputTag("QGTagger", "mult"));
  ptDToken = consumes<edm::ValueMap<float> >(edm::InputTag("QGTagger", "ptD"));

  produces<pat::JetCollection>();
}


void
JetFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //--- Get jets
  Handle<edm::View<pat::Jet> > jetHandle;
  iEvent.getByToken(jetToken, jetHandle);

  //--- q/g tagger
  Handle<edm::ValueMap<float> > qgHandle; 
  iEvent.getByToken(qgToken, qgHandle);
  Handle<edm::ValueMap<float> > axis2Handle; 
  iEvent.getByToken(axis2Token, axis2Handle);
  Handle<edm::ValueMap<int> > multHandle; 
  iEvent.getByToken(multToken, multHandle);
  Handle<edm::ValueMap<float> > ptDHandle; 
  iEvent.getByToken(ptDToken, ptDHandle);

  //-- JEC uncertanties 
  ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get(jecType,JetCorParColl);
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty jecUnc(JetCorPar);

  //--- Output collection
  auto_ptr<pat::JetCollection> result( new pat::JetCollection() );

  for(auto jet = jetHandle->begin(); jet != jetHandle->end(); ++jet){
    pat::Jet j(*jet);

    //--- Apply selection cut.
    if (!cut(j)) continue;

    //--- will re-embed b-tagger (so that it is chosen in one place only)
    float bTagger = j.bDiscriminator(bTaggerName);

    //--- Retrieve the q/g likelihood
    edm::RefToBase<pat::Jet> jetRef(edm::Ref<edm::View<pat::Jet> >(jetHandle, jet - jetHandle->begin()));
    float qgLikelihood = (*qgHandle)[jetRef];
    float axis2 = (*axis2Handle)[jetRef];
    int mult = (*multHandle)[jetRef];
    float ptD = (*ptDHandle)[jetRef];

    //--- Get JEC uncertainties 
    jecUnc.setJetEta(j.eta());
    jecUnc.setJetPt(j.pt());
    float jec_unc = jecUnc.getUncertainty(true);

    //--- Embed user variables
    j.addUserFloat("bTagger",bTagger);
    j.addUserFloat("qgLikelihood",qgLikelihood);
    j.addUserFloat("axis2",axis2);
    j.addUserFloat("mult",mult);
    j.addUserFloat("ptD",ptD);
    j.addUserFloat("jec_unc", jec_unc);

    //--- Embed flags (ie flags specified in the "flags" pset)
    for(CutSet<pat::Jet>::const_iterator flag = flags.begin(); flag != flags.end(); ++flag) {
      j.addUserFloat(flag->first,int((*(flag->second))(j)));
    }
    
    result->push_back(j);
  }

  iEvent.put(result);
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(JetFiller);

