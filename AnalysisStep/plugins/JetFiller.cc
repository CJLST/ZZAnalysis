#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>

#include <DataFormats/PatCandidates/interface/Jet.h>

#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>

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
  const CutSet<pat::Jet> flags;
  edm::EDGetTokenT<edm::ValueMap<float> > qgToken;

};


JetFiller::JetFiller(const edm::ParameterSet& iConfig) :
  jetToken(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("src"))),
  cut(iConfig.getParameter<std::string>("cut")),
  bTaggerName(iConfig.getParameter<std::string>("bTaggerName")),
  flags(iConfig.getParameter<edm::ParameterSet>("flags"))
{
  qgToken = consumes<edm::ValueMap<float> >(edm::InputTag("QGTagger", "qgLikelihood"));
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

    //--- Embed user variables
    j.addUserFloat("bTagger",bTagger);
    j.addUserFloat("qgLikelihood",qgLikelihood);

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

