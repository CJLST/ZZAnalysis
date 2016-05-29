/** \class bareZFiller
 *
 *  No description available.
 *
 *  \author N. Amapane (Torino)
 *  \author S. Bolognesi (JHU)
 *  \author C. Botta (CERN)
 *  \author S. Casasso (Torino)
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <CommonTools/Utils/interface/StringCutObjectSelector.h>
#include <CommonTools/Utils/interface/StringObjectFunction.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/Candidate/interface/ShallowCloneCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Photon.h>

#include <ZZAnalysis/AnalysisStep/interface/PhotonFwd.h>

#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>
#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>

#include <string>
#include <Math/VectorUtil.h>

using namespace std;

namespace {
//  bool recomputeIsoForFSR=false; // add "combRelIsoPFFSRCorr" for Z leptons (a' la Legacy)
}


class bareZFiller : public edm::EDProducer {
 public:
  /// Constructor
  explicit bareZFiller(const edm::ParameterSet&);
    
  /// Destructor
  ~bareZFiller(){};  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};
  
  //edm::EDGetTokenT<edm::View<reco::CompositeCandidate> > jetPairToken;
  edm::EDGetTokenT<edm::View<pat::Electron> > eleToken;
  edm::EDGetTokenT<edm::View<pat::Photon> > phToken;

//  const StringCutObjectSelector<pat::CompositeCandidate, true> preBestZSelection;
//  int sampleType;
//  int setup;

  const StringCutObjectSelector<pat::CompositeCandidate, true> cut;
  const CutSet<pat::CompositeCandidate> flags;
//  int FSRMode;
//  bool embedDaughterFloats;
};


bareZFiller::bareZFiller(const edm::ParameterSet& iConfig) :
 // jetPairToken(consumes<edm::View<reco::CompositeCandidate> >(iConfig.getParameter<edm::InputTag>("src"))),
  eleToken(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("src_electrons"))),
  phToken(consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("src_photons"))),

  //preBestZSelection(iConfig.getParameter<std::string>("bestZAmong")),
  //sampleType(iConfig.getParameter<int>("sampleType")),
  //setup(iConfig.getParameter<int>("setup")),
  cut(iConfig.getParameter<std::string>("cut")),
  flags(iConfig.getParameter<edm::ParameterSet>("flags"))
  //embedDaughterFloats(iConfig.getUntrackedParameter<bool>("embedDaughterFloats",true))
{
  produces<pat::CompositeCandidateCollection>();
}

bool greater_pt(pat::CompositeCandidate a, pat::CompositeCandidate b) {return a.pt() > b.pt(); }

void
bareZFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {  
  using namespace edm;
  using namespace std;
  using namespace reco;
  
  std::auto_ptr<pat::CompositeCandidateCollection> result(new pat::CompositeCandidateCollection);

  //LogPrint("") <<" IN JJCAND!"; 

  Handle<View<pat::Photon> > softPhotons;
  iEvent.getByToken(phToken, softPhotons);

  Handle<View<pat::Electron> > softElectrons;
  iEvent.getByToken(eleToken, softElectrons);
 
 
  std::auto_ptr<pat::CompositeCandidateCollection> e_ph_candidates(new pat::CompositeCandidateCollection);
 
//  vector<pat::CompositeCandidate> ZCandidates;
  
  for(size_t i_ph = 0; i_ph < softPhotons->size(); ++i_ph) {  //ph = softPhotons->begin(); ph != softPhotons->end(); ++ph) {
      for(size_t i_ele = 0; i_ele < softElectrons->size(); ++i_ele){ //auto ele = softElectrons->begin(); ele != softElectrons->end(); ++ele) {
          edm::Ptr<pat::Electron> electron = softElectrons->ptrAt(i_ele);
          edm::Ptr<pat::Photon> photon = softPhotons->ptrAt(i_ph);
          //photon->addUserFloat("charge", -1 * electron->charge());
          reco::CompositeCandidate z_cand = CompositeCandidate(electron->charge(), electron->p4() + photon->p4());
          z_cand.addDaughter(*electron);
          z_cand.addDaughter(*photon);

          LogPrint("") << "Built a Z candidate with mass " << z_cand.mass();
          pat::CompositeCandidate myCand(z_cand); 

          //if (!cut(myCand)) continue;    
          e_ph_candidates->push_back(myCand); 
      }
  }
  std::sort (e_ph_candidates->begin(), e_ph_candidates->end(), greater_pt);
     
    //LogPrint("") << "Photons:";
    //for(auto jet : jets) LogPrint("") << "jet pt " << jet.pt();

  //myCand.addUserFloat("is_best", 1);
 //for(auto ZCandidate : *e_ph_candidates) {

/*
 for (size_t i = 0; i < e_ph_candidates->size(); ++i) {
    pat::CompositeCandidate& ZCandidate = (*e_ph_candidates).at(i); //[i];
    //for(CutSet<pat::CompositeCandidate>::const_iterator flag = flags.begin(); flag != flags.end(); ++flag) {
    //    ZCandidate.addUserFloat(flag->first,int((*(flag->second))(ZCandidate)));
    //}

  }
*/
  LogPrint("") <<"expected size tle_Z"<<e_ph_candidates->size();
  iEvent.put(e_ph_candidates);
}



#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(bareZFiller);
