#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/EventSetup.h>

#include <DataFormats/PatCandidates/interface/Jet.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <JetMETCorrections/Objects/interface/JetCorrectionsRecord.h>
#include <JetMETCorrections/Modules/interface/JetResolution.h>

#include <DataFormats/Math/interface/deltaR.h>

#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>

#include "TRandom3.h"

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
  int sampleType;
  int setup;
  const StringCutObjectSelector<pat::Jet, true> cut;
  bool isMC_;
  const std::string bTaggerName;
  const std::string jecType;
  bool applyJER_;
  const std::string jerType;
  const CutSet<pat::Jet> flags;
  edm::EDGetTokenT<double> rhoToken;
  edm::EDGetTokenT<edm::ValueMap<float> > qgToken;
  edm::EDGetTokenT<edm::ValueMap<float> > axis2Token;
  edm::EDGetTokenT<edm::ValueMap<int> > multToken;
  edm::EDGetTokenT<edm::ValueMap<float> > ptDToken;

  JME::JetResolution resolution_pt; 
  JME::JetResolution resolution_phi;
  JME::JetResolutionScaleFactor resolution_sf;

};


JetFiller::JetFiller(const edm::ParameterSet& iConfig) :
  jetToken(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("src"))),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  cut(iConfig.getParameter<std::string>("cut")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  bTaggerName(iConfig.getParameter<std::string>("bTaggerName")),
  jecType(iConfig.getParameter<std::string>("jecType")),
  applyJER_(iConfig.getParameter<bool>("applyJER")),
  jerType(iConfig.getParameter<std::string>("jerType")),
  flags(iConfig.getParameter<edm::ParameterSet>("flags"))
{

  rhoToken = consumes<double>(LeptonIsoHelper::getMuRhoTag(sampleType, setup));

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

  //--- rho
  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken, rhoHandle);
  double rho = *rhoHandle;

  //--- q/g tagger
  Handle<edm::ValueMap<float> > qgHandle; 
  iEvent.getByToken(qgToken, qgHandle);
  Handle<edm::ValueMap<float> > axis2Handle; 
  iEvent.getByToken(axis2Token, axis2Handle);
  Handle<edm::ValueMap<int> > multHandle; 
  iEvent.getByToken(multToken, multHandle);
  Handle<edm::ValueMap<float> > ptDHandle; 
  iEvent.getByToken(ptDToken, ptDHandle);

  //--- JEC uncertanties 
  ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get(jecType,JetCorParColl);
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty jecUnc(JetCorPar);

  //--- JER
  resolution_pt = JME::JetResolution::get(iSetup, jerType+"_pt");
  resolution_phi = JME::JetResolution::get(iSetup, jecType+"_phi");
  resolution_sf = JME::JetResolutionScaleFactor::get(iSetup, jecType);

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

    //--- JER
    j.addUserFloat("pt_nojer", j.pt());
    float pt_jer   = -1.;
    float pt_jerup = -1.;
    float pt_jerdn = -1.;
    if(isMC_ && applyJER_){
      JME::JetParameters res_parameters = {{JME::Binning::JetPt, j.pt()}, {JME::Binning::JetEta, j.eta()}, {JME::Binning::Rho, rho}};
      float res_pt  = resolution_pt .getResolution(res_parameters);
      //float res_phi = resolution_phi.getResolution(res_parameters); //not used

      JME::JetParameters sf_parameters = {{JME::Binning::JetEta, j.eta()}, {JME::Binning::Rho, rho}};
      float sf    = resolution_sf.getScaleFactor(sf_parameters);
      float sf_up = resolution_sf.getScaleFactor(sf_parameters, Variation::UP);
      float sf_dn = resolution_sf.getScaleFactor(sf_parameters, Variation::DOWN);

      //-- 'hybrid' method = scaling for well matched jets, smearing for the rest
      //- check matching to gen jets
      const reco::GenJet* genJet = j.genJet();
      bool matchedJet = genJet 
	&& ( reco::deltaR(j.eta(),j.phi(),genJet->eta(),genJet->phi()) < 0.2 ) 
	&& ( fabs(j.pt()-genJet->pt()) < 3*res_pt*j.pt() );
      if(matchedJet){
	//- apply scaling
	float gen_pt = genJet->pt();
	pt_jer   = max( 0., gen_pt + sf   *(j.pt()-gen_pt) );
	pt_jerup = max( 0., gen_pt + sf_up*(j.pt()-gen_pt) );
	pt_jerdn = max( 0., gen_pt + sf_dn*(j.pt()-gen_pt) );
      }else{
	//- apply smearing
	TRandom3 rand;
	rand.SetSeed(abs(static_cast<int>(sin(j.phi())*100000)));
	float smear = rand.Gaus(0,1.);
	float sigma   = sqrt(sf   *sf   -1.) * res_pt*j.pt();
	float sigmaup = sqrt(sf_up*sf_up-1.) * res_pt*j.pt();
	float sigmadn = sqrt(sf_dn*sf_dn-1.) * res_pt*j.pt();
	pt_jer   = max( 0., smear*sigma   + j.pt() );
	pt_jerup = max( 0., smear*sigmaup + j.pt() );
	pt_jerdn = max( 0., smear*sigmadn + j.pt() );
      } 

      j.setP4(reco::Particle::PolarLorentzVector(pt_jer, j.eta(), j.phi(), (pt_jer/j.pt())*j.mass()));
    }

    //--- Embed user variables
    j.addUserFloat("bTagger",bTagger);
    j.addUserFloat("qgLikelihood",qgLikelihood);
    j.addUserFloat("axis2",axis2);
    j.addUserFloat("mult",mult);
    j.addUserFloat("ptD",ptD);
    j.addUserFloat("jec_unc", jec_unc);
    j.addUserFloat("pt_jerup", pt_jerup);
    j.addUserFloat("pt_jerdn", pt_jerdn);

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

