/** \class ZCandidateFiller
 *
 *  No description available.
 *
 *  $Date: 2012/10/10 22:28:52 $
 *  $Revision: 1.14 $
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
#include <ZZAnalysis/AnalysisStep/interface/PhotonFwd.h>

#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>
#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>

#include <string>
#include <Math/VectorUtil.h>

class ZCandidateFiller : public edm::EDProducer {
 public:
  /// Constructor
  explicit ZCandidateFiller(const edm::ParameterSet&);
    
  /// Destructor
  ~ZCandidateFiller(){};  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};
  
  edm::InputTag theCandidateTag;
  const StringCutObjectSelector<pat::CompositeCandidate, true> preBestZSelection;
  int sampleType;
  int setup;
  const CutSet<pat::CompositeCandidate> cuts;
  bool embedDaughterFloats;
};


ZCandidateFiller::ZCandidateFiller(const edm::ParameterSet& iConfig) :
  theCandidateTag(iConfig.getParameter<edm::InputTag>("src")),
  preBestZSelection(iConfig.getParameter<std::string>("bestZAmong")),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  cuts(iConfig.getParameter<edm::ParameterSet>("flags")),
  embedDaughterFloats(iConfig.getUntrackedParameter<bool>("embedDaughterFloats",true))
{
  produces<pat::CompositeCandidateCollection>();
}


void
ZCandidateFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  
  using namespace edm;
  using namespace std;
  using namespace reco;
  
  std::auto_ptr<pat::CompositeCandidateCollection> result(new pat::CompositeCandidateCollection);

  //-- Get LL candidates
  Handle<View<reco::CompositeCandidate> > LLCands;
  iEvent.getByLabel(theCandidateTag, LLCands);

#define USE_FSR
#ifdef USE_FSR
  // Get rho, to recompute isolation for leptons with FSR
  double rhoForMu, rhoForEle;
  {
    edm::Handle<double> rhoHandle;
    iEvent.getByLabel(LeptonIsoHelper::getMuRhoTag(sampleType, setup), rhoHandle);
    rhoForMu = *rhoHandle;
    iEvent.getByLabel(LeptonIsoHelper::getEleRhoTag(sampleType, setup), rhoHandle);
    rhoForEle = *rhoHandle;
  }
#endif

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

#define USE_FSR
#ifdef USE_FSR
    // ------------------------------
    // FSR recovery
    // ------------------------------

    //loop on the 2 daughters; apply mass cuts on (llg) and store 
    // the highest-pT and the lowest-DR assocated gamma.
    double    maxPT = -1.;
    double    minDR = 9999.;
    const pat::PFParticle* maxPTg=0;
    const pat::PFParticle* minDRg=0;
    int maxPTgLep=-1; // Index of daughter to which the above photons
    int minDRgLep=-1; // are associated

    for (int dauIdx=0; dauIdx<2; ++dauIdx) { 
      const Candidate* d = myCand.daughter(dauIdx);
      const PhotonPtrVector* gammas = userdatahelpers::getUserPhotons(d);
      if (gammas==0) continue;
      for (PhotonPtrVector::const_iterator g = gammas->begin();
	   g!= gammas->end(); ++g) {
	const pat::PFParticle* gamma = g->get();
	reco::Candidate::LorentzVector p4G = gamma->p4();
	reco::Candidate::LorentzVector p4LL = myCand.p4();
	double mLLG = (p4LL + p4G).M();
	bool movesToZPeak = (fabs(mLLG-ZmassValue) < fabs(myCand.mass()-ZmassValue));
// // Debug
// 	myCand.addUserFloat("mass",myCand.mass());
// 	myCand.addUserFloat("mLLG",mLLG);
// 	myCand.addUserFloat("movesToZPeak",movesToZPeak);
	if (movesToZPeak && mLLG<100. && mLLG>4) { // Mass cuts (4 is implicit)

	  double pt = gamma->pt();
	  if (pt>maxPT) {
	    maxPT  = pt;
	    maxPTg = gamma;
	    maxPTgLep = dauIdx;
	  }
	  
	  double dR = ROOT::Math::VectorUtil::DeltaR(gamma->momentum(),d->momentum());
	  if (dR<minDR) {
	    minDR  = dR;
	    minDRg = gamma;
	    minDRgLep = dauIdx;
	  }
	}
      } // end loop on photons
    } // end loop on daughters (leptons)
    
    // Define the selected FSR photon.
    const pat::PFParticle* fsr=0;    
    int lepWithFsr=-1; 
    if (maxPTg!=0) { // at least 1 photon selected
      if (maxPT>4) { // First case: take highest-pT
	fsr=maxPTg;
	lepWithFsr=maxPTgLep;
      } else {
	fsr=minDRg;
	lepWithFsr=minDRgLep;
      }
    }

    myCand.addUserFloat("dauWithFSR",lepWithFsr); // Index of the cand daughter with associated FSR photon

    if (fsr!=0) {
      // Add daughter and set p4.
      myCand.addUserFloat("mll",myCand.mass()); // for debug purposes
      myCand.setP4(myCand.p4()+fsr->p4());
//      myCand.addDaughter(reco::ShallowCloneCandidate(fsr->masterClone()),"FSR"); //FIXME: fsr does not have a masterClone
      pat::PFParticle myFsr(*fsr);
      myFsr.setPdgId(22); // Fix: photons that are isFromMu have abs(pdgId)=13!!!
      myCand.addDaughter(myFsr,"FSR");
    }

    // Recompute iso for leptons with FSR    
    for (int dauIdx=0; dauIdx<2; ++dauIdx) { 
      const Candidate* d = myCand.daughter(dauIdx);
      float fsrCorr = 0; // The correction to PFPhotonIso
      if (fsr!=0) {
	//	if (!fsr->isFromMuon()) { // Type 1 photons should be subtracted from muon iso cones
	double dR = ROOT::Math::VectorUtil::DeltaR(fsr->momentum(),d->momentum());
	if (dR<0.4 && ((d->isMuon() && dR > 0.01) ||
		       (d->isElectron() && (fabs((static_cast<const pat::Electron*>(d->masterClone().get()))->superCluster()->eta()) < 1.479 || dR > 0.08)))) {
	  fsrCorr = fsr->pt();
	}
	//}
      }

      float rho = ((d->isMuon())?rhoForMu:rhoForEle);
      float combRelIsoPFCorr =  LeptonIsoHelper::combRelIsoPF(sampleType, setup, rho, d, fsrCorr);
      
      string base;
      stringstream str;
      str << "d" << dauIdx << ".";
      str >> base;	  
      myCand.addUserFloat(base+"combRelIsoPFFSRCorr",combRelIsoPFCorr);
    } // end loop on daughters
#else 
    // if FSR is not activated
    myCand.addUserFloat("dauWithFSR",-1);
    myCand.addUserFloat("d0.combRelIsoPFFSRCorr",myCand.userFloat("d0.combRelIsoPF"));
    myCand.addUserFloat("d1.combRelIsoPFFSRCorr",myCand.userFloat("d1.combRelIsoPF"));

#endif

    //--- Find "best Z" (closest to mZ) among those passing the "bestZAmong" selection (2011 PRL logic)
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

    //--- Embed shortcut variables
    // Pairwise isolation
    float isoSum2 = myCand.userFloat("d0.combRelIso") + myCand.userFloat("d1.combRelIso");
    myCand.addUserFloat("isoSum2",isoSum2);
//     myCand.addUserFloat("OS",OS);
//     myCand.addUserFloat("SF",SF);

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
  
  iEvent.put(result);
}



#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ZCandidateFiller);
