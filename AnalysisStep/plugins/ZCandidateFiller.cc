/** \class ZCandidateFiller
 *
 *  Create a collection of Zs including FSR and additional variables stored as userFloats.
 *
 *  \author N. Amapane (Torino)
 *  \author S. Bolognesi (JHU)
 *  \author C. Botta (CERN)
 *  \author S. Casasso (Torino)
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/one/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <CommonTools/Utils/interface/StringCutObjectSelector.h>
#include <CommonTools/Utils/interface/StringObjectFunction.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/Candidate/interface/ShallowCloneCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <ZZAnalysis/AnalysisStep/interface/PhotonFwd.h>

#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>
#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>

#include <string>
#include <Math/VectorUtil.h>

using namespace std;

class ZCandidateFiller : public edm::one::EDProducer<> {
 public:
  /// Constructor
  explicit ZCandidateFiller(const edm::ParameterSet&);
    
  /// Destructor
  ~ZCandidateFiller(){};  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};
  
  edm::EDGetTokenT<edm::View<reco::CompositeCandidate> > candidateToken;
  const StringCutObjectSelector<pat::CompositeCandidate, true> preBestZSelection;
  int sampleType;
  int setup;
  const CutSet<pat::CompositeCandidate> cuts;
  int FSRMode;
  bool embedDaughterFloats;
};


ZCandidateFiller::ZCandidateFiller(const edm::ParameterSet& iConfig) :
  candidateToken(consumes<edm::View<reco::CompositeCandidate> >(iConfig.getParameter<edm::InputTag>("src"))),
  preBestZSelection(iConfig.getParameter<std::string>("bestZAmong")),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  cuts(iConfig.getParameter<edm::ParameterSet>("flags")),
  embedDaughterFloats(iConfig.getUntrackedParameter<bool>("embedDaughterFloats",true))
{
  
  string mode = iConfig.getParameter<string>("FSRMode");
  if      (mode == "skip")   FSRMode = 0;
  // "Legacy" Run I mode (1) is no longer supported.
  else if (mode == "RunII")  FSRMode = 2;
  else {
    cout << "ZCandidateFiller: FSRMode " << FSRMode << " not supported" << endl;
    abort();
  }
  
  produces<pat::CompositeCandidateCollection>();
}


void
ZCandidateFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  
  using namespace edm;
  using namespace std;
  using namespace reco;
  
  auto result = std::make_unique<pat::CompositeCandidateCollection>();

  //-- Get LL candidates
  Handle<View<reco::CompositeCandidate> > LLCands;
  iEvent.getByToken(candidateToken, LLCands);

  //--- Fill user info

  const float ZmassValue = 91.1876;

  float closestZeeMassDiff = 99999.;
  float closestZmmMassDiff = 99999.;
  float closestZMassDiff   = 99999.;
  float closestLLMassDiff = 99999.;

  int bestZeeidx = -1;
  int bestZmmidx = -1;
  int bestZidx = -1;
  int bestLLidx = -1;

  //--- Loop over LL Candidates
  for(unsigned int i = 0; i < LLCands->size(); ++i) {
    const CompositeCandidate& c = (*LLCands)[i];
    pat::CompositeCandidate myCand(c); 

    if (embedDaughterFloats){  
      userdatahelpers::embedDaughterData(myCand);
    }

    int id0 = myCand.daughter(0)->pdgId();
    int id1 = myCand.daughter(1)->pdgId();
    bool OS = (id0*id1)<0;
    bool SF = abs(id0)==abs(id1);

    // ------------------------------
    // FSR recovery
    // ------------------------------
    if (FSRMode==2) { // Run II
      float mll = myCand.mass(); // pre-FSR mass
      for (int dauIdx=0; dauIdx<2; ++dauIdx) { 
        const Candidate* d = myCand.daughter(dauIdx);
        const PhotonPtrVector* gammas = userdatahelpers::getUserPhotons(d);
        if (gammas==0) continue;
        assert(gammas->size()<=1); // Must have already been preselected.
        if (gammas->size()==1){
          const pat::PFParticle* fsr = gammas->begin()->get();
          pat::PFParticle myFsr(*fsr);
          myCand.setP4(myCand.p4()+fsr->p4());
          myFsr.setPdgId(22);
          myFsr.addUserFloat("leptIdx",dauIdx);
          //        myFsr.addUserFloat("gRelIso",0.);
          myCand.addDaughter(myFsr);
        }
      }
      
      if (myCand.numberOfDaughters()>2) {
        myCand.addUserFloat("mll",mll);
      }

    } else if (FSRMode==0) { // no FSR
      myCand.addUserFloat("dauWithFSR",-1);
      myCand.addUserFloat("d0.combRelIsoPFFSRCorr",myCand.userFloat("d0.combRelIsoPF"));
      myCand.addUserFloat("d1.combRelIsoPFFSRCorr",myCand.userFloat("d1.combRelIsoPF"));
    }

    
    //--- Find "best Z" (closest to mZ) among those passing the "bestZAmong" selection (2011 PRL logic), now deprecated!!!
    if (preBestZSelection(myCand)) {
      float diffZmass = fabs(ZmassValue - myCand.mass());
      if (diffZmass < closestLLMassDiff) { // Best among any ll in the collection
        bestLLidx = i;
        closestLLMassDiff = diffZmass;
      }
      if (OS&&SF) {
        if (diffZmass < closestZMassDiff) { // Best among all OSSF pairs in the collection
          bestZidx = i;
          closestZMassDiff = diffZmass;
        }
        if (abs(id0) == 13) { 
          if (diffZmass < closestZmmMassDiff) { // Best among all mu+mu- pairs in the collection
            bestZmmidx = i;
            closestZmmMassDiff = diffZmass;
          }
        } else if (abs(id0) == 11) {
          if (diffZmass < closestZeeMassDiff) { // Best among all e+e- pairs in the collection
            bestZeeidx = i;
            closestZeeMassDiff = diffZmass;
          }
        }
      }
    }

    result->push_back(myCand);
  }

  //--- Embed isBestZ flag (must be done in a separate loop)
  for (int i = 0; i< (int)result->size(); ++i) {
    pat::CompositeCandidate& myCand = (*result)[i];    
    myCand.addUserFloat("isBestZ",  (i==bestZidx));
    myCand.addUserFloat("isBestZmm",(i==bestZmmidx));
    myCand.addUserFloat("isBestZee",(i==bestZeeidx));
    myCand.addUserFloat("isBestInColl", (i==bestLLidx));

    //--- Embed flags (ie cuts specified in the "flags" pset)
    //    We do this here so that isBestZ is available within the cuts
    for(CutSet<pat::CompositeCandidate>::const_iterator cut = cuts.begin(); cut != cuts.end(); ++cut) {
      myCand.addUserFloat(cut->first,int((*(cut->second))(myCand)));
    }
  }
  
  iEvent.put(std::move(result));
}



#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ZCandidateFiller);
