/** \class FSRAnalyzer
 *
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/LuminosityBlock.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <DataFormats/Common/interface/View.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>
#include <ZZAnalysis/AnalysisStep/interface/MCHistoryTools.h>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>

#include <Math/VectorUtil.h>

#include <iostream>
#include <typeinfo>

//bool debug = true;

class FSRAnalyzer : public edm::EDAnalyzer {
 public:
  /// Constructor
  FSRAnalyzer(const edm::ParameterSet& pset);

  virtual void beginJob();
  virtual void endJob();  


  /// Destructor
  ~FSRAnalyzer();

  // Operations
  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);  

};

FSRAnalyzer::FSRAnalyzer(const edm::ParameterSet& pset){}

FSRAnalyzer::~FSRAnalyzer(){}


using namespace std;
using namespace edm;


void FSRAnalyzer::beginJob(){
}

void FSRAnalyzer::endJob(){
}


void
FSRAnalyzer::analyze(const edm::Event & event, const edm::EventSetup& eventSetup) {
  
  const float ZmassValue = 91.1876;

  Handle<View<pat::CompositeCandidate> > ZColl;
  event.getByLabel("ZCand", ZColl);

  edm::Handle<edm::View<pat::PackedCandidate> > pfCands; 
  event.getByLabel("packedPFCandidates",pfCands);

  MCHistoryTools mch(event);
  // These are all gen FSR photons coming from leptons from the H
  vector<const reco::Candidate *> genFSR = mch.genFSR();
  //  vector<const reco::Candidate *> genLep = mch.genZLeps();
  int genFinalState = mch.genFinalState();

  vector<int> nRecoFSRMatchedToGen(genFSR.size(),0);

  // NOTES
  // Since we want to check the effect of the Z mass criteria on FSR, we start from reconstructed Zs.  // To avoid possible double counting due to combinatorics we have to run on 2e2mu events skipping 
  // events with additional leptons (ie more than 2 Z candidates)
  // We then compare collected FSRs to all gen-level FSR, regardless of the matching of reco 
  // to gen leptons. In this way an FSR is considered correctly reconstructed even if it is 
  // not attached to the lepton it was originated from (which is fine since we care about 
  // the 4l mass at the end); also, in this way we can provide numer per-photon when there 
  // are 2 or more FSR from the same lepton.
  // OTOH, FSR efficiency is a bit more tricky since a gen FSR should be counted in the 
  // denominator only  when the corresponding lepton is reconstructed; to avoid ambiguities 
  // or a complicated selection I restrict to fully reconstructed events 
  // (ie exactly 2 Zs in 2e2mu events) for the time being.

  if (genFinalState==EEMM) {
    if (ZColl->size()!=2) {
      if (ZColl->size()>2) cout << "WARNING: Additional leptons!" << endl;      
      return;
    }
    

    // We now have exactly 2 Zs 
    int zzflavour=1;
    reco::Candidate::LorentzVector p4ZZ;
    for( View<pat::CompositeCandidate>::const_iterator Z = ZColl->begin(); 
	 Z != ZColl->end(); ++Z) {
      p4ZZ += (Z->daughter(0)->p4()+Z->daughter(1)->p4());
      zzflavour*=Z->daughter(0)->pdgId()*Z->daughter(1)->pdgId();
    }
    
    //also skip events where the reconstructed final state is not 2e2mu
    if (zzflavour!=11*11*13*13) {
      cout << "WARNING: 2Zs, but wrong flavour!" << endl;
      return;
    }
    
    for( View<pat::CompositeCandidate>::const_iterator Z = ZColl->begin(); 
	 Z != ZColl->end(); ++Z) {
      
      // Do not use Z->p4() as this may already include FSR, in our standard workflow
      reco::Candidate::LorentzVector p4LL= Z->daughter(0)->p4()+Z->daughter(1)->p4();

      for (int idau=0; idau<2; ++idau) {
	const reco::Candidate* l=Z->daughter(idau);
	
	int lID = l->pdgId();

	// Photons tentatively associated to this lepton (may be >1)
	const PhotonPtrVector* gammas = userdatahelpers::getUserPhotons(l);

	//	cout << "LEP " << l->pdgId() << " " << l->pt() << " " << (gammas?gammas->size():0) << endl;
	if (gammas==0) continue;


	for (PhotonPtrVector::const_iterator g = gammas->begin();
	     g!= gammas->end(); ++g) {
	  const pat::PFParticle* gamma = g->get();
	  reco::Candidate::LorentzVector p4G = gamma->p4();
	  double dR = ROOT::Math::VectorUtil::DeltaR(gamma->momentum(),l->momentum());
	  double pT = p4G.pt();
	  double mLL = p4LL.M();
	  double mLLG = (p4LL + p4G).M();
	  bool movesToZPeak = (fabs(mLLG-ZmassValue) < fabs(mLL-ZmassValue));
	  bool isFake  = true;
	  double dRGenVsReco = -1.;
	  double pTGen = -1.;
	  double etaGen = 0;
	  double phiGen = 0.;
	  int igen = MCHistoryTools::fsrMatch(g->get(),genFSR);

	  double neu, chg;
	  LeptonIsoHelper::fsrIso(gamma, pfCands, neu, chg);
	  double gRelIso = (neu + chg)/gamma->pt();	  
	  
	  if (igen>=0) {
	    dRGenVsReco = ROOT::Math::VectorUtil::DeltaR(genFSR[igen]->momentum(),gamma->momentum());
	    if (dRGenVsReco<0.3) { //Matching cut -- FIXME
	      isFake=false;
	      nRecoFSRMatchedToGen[igen]++;
	      pTGen = genFSR[igen]->pt();
	      etaGen = genFSR[igen]->eta();
	      phiGen = genFSR[igen]->phi();
	    }
	  }

	  cout << "FSR:" 
	       << event.id().run() << ":"
	       << event.id().luminosityBlock() << ":"
	       << event.id().event() << " "
	       << lID << " " // this is the ID of the reco l the photon is associated to 
	       << dR << " "  // reco FSR vs reco lep
	       << pT << " "  // reco FSR pT
	       << gamma->eta() << " " 
	       << gamma->phi() << " "
	       << gRelIso << " " // photon iso
	       << mLL  << " "
	       << mLLG << " " 
	       << p4ZZ.M() << " "  // mass of the reconstructed 4l 
	       << (p4ZZ+p4G).M() << " " // mass of the reconstructed 4l + g
	       << movesToZPeak << " "
	       << isFake << " " 
	       << dRGenVsReco << " " // dR of gen FSR to reco FSR, if not fake 
	       << pTGen << " "       // pT of gen FSR, if not fake
	       << etaGen << " "
	       << phiGen << " "
	       << endl; 

	  if (isFake && dR<0.07 && pT>20 && gRelIso<1) cout << "DEBUG: "
							    << event.id().run() << ":"
							    << event.id().luminosityBlock() << ":"
							    << event.id().event() << endl;

	}
      }
    }
    
    // Now check if any of the gen FSR has not been associated to any lepton.
    // This is meaningful only for fully reconstructed events (eg we do not care about 
    // FSR for for non reconstructed or out of acceptance leptons)

    for (unsigned j=0; j<genFSR.size(); ++j) {
      if (genFSR[j]->pt()>2. ){
	bool lost=false;
	if (nRecoFSRMatchedToGen[j]==0) {
	  lost = true; // FIXME would be nice to add DR to the gen leptton
	}
	cout << "LOST " << genFSR[j]->pt() << " "
	     << genFSR[j]->eta() << " "
	     << genFSR[j]->phi() << " "
	     << " " << lost << endl;
      }
    }
  }
}





#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(FSRAnalyzer);
