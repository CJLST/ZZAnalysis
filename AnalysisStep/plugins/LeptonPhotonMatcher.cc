/** \class LeptonPhotonMatcher
 *
 *  No description available.
 *
 *  $Date: 2012/10/17 11:32:15 $
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

#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/PFParticle.h>
#include <ZZAnalysis/AnalysisStep/interface/PhotonFwd.h>
#include <Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h>
#include <EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h>
#include <DataFormats/GeometryVector/interface/VectorUtil.h> 
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidate.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>

#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>

#include <ZZAnalysis/AnalysisStep/interface/MCHistoryTools.h>
#include <Math/VectorUtil.h>
#include <TMath.h>

#include <vector>
#include <string>



using namespace edm;
using namespace std;
using namespace reco;

class LeptonPhotonMatcher : public edm::EDProducer {
 public:
  /// Constructor
  explicit LeptonPhotonMatcher(const edm::ParameterSet&);
    
  /// Destructor
  ~LeptonPhotonMatcher(){};  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};
  PhotonPtr selectFSR(const PhotonPtrVector& photons, const reco::LeafCandidate::Vector& lepMomentum);

  edm::EDGetTokenT<pat::MuonCollection> muonToken;
  edm::EDGetTokenT<pat::ElectronCollection> electronToken;
  //edm::EDGetTokenT<pat::ElectronCollection> looseElectronToken;
  edm::EDGetTokenT<pat::PhotonCollection> tleToken;

  edm::EDGetTokenT<edm::View<pat::PFParticle> > photonToken;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate> > pfCandToken;
  int selectionMode;
  int sampleType;
  int setup;
  bool debug;
  edm::EDGetTokenT<double> rhoForMuToken;
  edm::EDGetTokenT<double> rhoForEleToken;

  float muon_iso_cut;
  float electron_iso_cut;
};


LeptonPhotonMatcher::LeptonPhotonMatcher(const edm::ParameterSet& iConfig) :
  muonToken(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
  electronToken(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronSrc"))),
  //looseElectronToken(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("looseElectronSrc"))),
  tleToken(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("tleSrc"))),
  photonToken(consumes<edm::View<pat::PFParticle> >(iConfig.getParameter<edm::InputTag>("photonSrc"))),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  debug(iConfig.getUntrackedParameter<bool>("debug",false))
{
  pfCandToken = consumes<edm::View<pat::PackedCandidate> >(edm::InputTag("packedPFCandidates"));
  rhoForMuToken = consumes<double>(LeptonIsoHelper::getMuRhoTag(sampleType, setup));
  rhoForEleToken = consumes<double>(LeptonIsoHelper::getEleRhoTag(sampleType, setup));

  string mode = iConfig.getParameter<string>("photonSel");
  
  if      (mode == "skip")        selectionMode = 0; // no FSR
  else if (mode == "passThrough") selectionMode = 1; // for debug
  else if (mode == "Legacy")      selectionMode = 2;
  else if (mode == "RunII")       selectionMode = 3;
  else {
    cout << "LeptonPhotonMatcher: mode " << mode << " not supported" << endl;
    abort();
  }
  
  muon_iso_cut = iConfig.getParameter<double>("muon_iso_cut");
  electron_iso_cut = iConfig.getParameter<double>("electron_iso_cut");


  produces<pat::MuonCollection>("muons");
  produces<pat::ElectronCollection>("electrons");
  //produces<pat::ElectronCollection>("looseElectrons");
  produces<pat::PhotonCollection>("electronstle");

}


void
LeptonPhotonMatcher::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //--- Get leptons and rho
  //  edm::Handle<pat::MuonRefVector> muonHandle;
  edm::Handle<pat::MuonCollection> muonHandle;
  iEvent.getByToken(muonToken, muonHandle);

  //  edm::Handle<pat::ElectronRefVector> electronHandle;
  edm::Handle<pat::ElectronCollection> electronHandle;
  iEvent.getByToken(electronToken, electronHandle);

//  edm::Handle<pat::ElectronCollection> looseElectronHandle;
//  iEvent.getByToken(looseElectronToken, looseElectronHandle);


  //  edm::Handle<pat::ElectronRefVector> electronHandle;
  edm::Handle<pat::PhotonCollection> tleHandle;
  iEvent.getByToken(tleToken, tleHandle);

  //--- Get the photons
  edm::Handle<edm::View<pat::PFParticle> > photonHandle;
  iEvent.getByToken(photonToken, photonHandle);

  //--- Get the PF cands
  edm::Handle<edm::View<pat::PackedCandidate> > pfCands; 
  iEvent.getByToken(pfCandToken, pfCands);

  // Output collections
  auto_ptr<pat::MuonCollection> resultMu( new pat::MuonCollection() );
  auto_ptr<pat::ElectronCollection> resultEle( new pat::ElectronCollection() );
//  auto_ptr<pat::ElectronCollection> resultLooseEle( new pat::ElectronCollection() );
  auto_ptr<pat::PhotonCollection> resultTle( new pat::PhotonCollection() );

  // Associate a vector of Ptr<Photon> to lepton pointers
  typedef map<const reco::Candidate*, PhotonPtrVector> PhotonLepMap;
  PhotonLepMap theMap;


  if (selectionMode!=0 && muonHandle->size()+electronHandle->size()>0) {
    //----------------------
    // Loop on photons
    //----------------------
    for (unsigned int i=0;i<photonHandle->size();++i) {
      
      // Get the photon as edm::Ptr
      PhotonPtr g = photonHandle->ptrAt(i);

      // Photon preselection (is currently already applied on pat::Photon collection)
      if (!(g->pt()>2. && fabs(g->eta())<2.4)) continue;

      //---------------------
      // // Supercluster veto
      //---------------------
      bool SCVeto=false;
      for (unsigned int j = 0; j< electronHandle->size(); ++j){
	const pat::Electron* e = &((*electronHandle)[j]);
	if (e->userFloat("isSIP")){
	  double dR = reco::deltaR(*(e->superCluster()), *g);
	  if ((fabs(g->eta() - e->superCluster()->eta())<0.05 && fabs(reco::deltaPhi(g->phi(), e->superCluster()->phi()))<2.) || dR<0.15) {	    
	    SCVeto=true;
  	    if (debug) cout << "SC veto: "<< g->pt() << " " << e->pt() << " " << dR << " "
			    << fabs(g->eta() - e->superCluster()->eta()) << " " << reco::deltaPhi(g->phi(), e->superCluster()->phi()) <<endl;
	    break;
	  } 
	}
      }
      if (debug) cout << "GAMMA: " << g->pt() << " " << g->eta() << " " << g->phi() << " SCVeto: " << SCVeto << endl;
      if (SCVeto) continue;


      //------------------------------------------------------
      // Get the closest lepton among those satisfying loose ID + SIP
      //------------------------------------------------------
      double dRMin(10e9);
      const reco::Candidate* closestLep = 0;

      // Loop over pat::Muon
      for (unsigned int j = 0; j< muonHandle->size(); ++j){
	//      const pat::Muon* m = ((*muonHandle)[j]).get();
	const pat::Muon* m = &((*muonHandle)[j]);
	if (! m->userFloat("isSIP")) continue; 
	double dR = ROOT::Math::VectorUtil::DeltaR(m->momentum(),g->momentum());
	if (dR>0.5) continue;
	if (dR<dRMin) {
	  dRMin = dR;
	  closestLep = m;
	}
      }//end loop over muon collection

      //---------------------
      // Loop over pat::Electron
      //---------------------
      for (unsigned int j = 0; j< electronHandle->size(); ++j){
	//      const pat::Electron* e = ((*electronHandle)[j]).get();
	const pat::Electron* e = &((*electronHandle)[j]);
	if ( ! e->userFloat("isSIP")) continue;
	double dR = ROOT::Math::VectorUtil::DeltaR(e->momentum(),g->momentum());
	if (dR>0.5) continue;
	if (dR<dRMin) {
	  dRMin = dR;
	  closestLep = e;
	}
      }//end loop over electron collection

      // Add photon to the vector that will be attached as userData for the corresponding lepton 
      if(closestLep!=0) {
	// Now that we know the closest lepton, apply Photon Selection
	bool accept = false;
	double gRelIso = 999., neu(999.), chg(999.), chgByWorstPV(999.);
	double pT = g->pt();

	if (selectionMode==1) { // passThrough: no photon selection, for FSR studies
	  accept = (dRMin<0.5 && pT>2.); 

	} else if (selectionMode==3) { // RunII
	  if (dRMin<0.5 && g->pt()>2. && dRMin/pT/pT<0.012) {
	    LeptonIsoHelper::fsrIso(&(*g), pfCands, neu, chg, chgByWorstPV);
	    gRelIso = (neu + chg)/pT;
	    if (gRelIso<1.8) accept = true;
	  }
	} else if (selectionMode==2) { // Legacy
	  if( dRMin<0.07 ){
	    if (g->pt()>2.) accept = true;
	  } else if (g->pt()>4 && dRMin<0.5 ){ // DR<0.5 is implicit, but does not hurt
	    // double relIso = g->relIso(0.5); // This is buggy, needs to recompute it.
	    LeptonIsoHelper::fsrIso(&(*g), pfCands, neu, chg, chgByWorstPV);
	    gRelIso = (neu + chg)/g->pt();
	    // For collections where this is precomputed
	    // double gRelIso2 = (g->userFloat("fsrPhotonPFIsoChHadPUNoPU03pt02") + g->userFloat("fsrPhotonPFIsoNHadPhoton03")) / g->pt();
	    if (gRelIso<1.) accept = true;
	  }
	}

	if(debug) cout << "   " << "   closest lep: " << closestLep->pdgId() << " " << closestLep->pt() <<  " gRelIso: " << gRelIso << " (ch: " << chg << " n+p: " <<  neu << " ) " << " dRMin: " << dRMin << " accept: " << accept << endl;
	if (accept) theMap[closestLep].push_back(g);
      }
    } // end of loop over photon collection
  }
  
  // Loop over muons again to write the result as userData
  PhotonPtrVector allSelFSR;
  for (unsigned int j = 0; j< muonHandle->size(); ++j){
    //    const pat::Muon* m = ((*muonHandle)[j]).get(); // Pointer to original mu
    const pat::Muon* m = &((*muonHandle)[j]);
    //---Clone the pat::Muon
    pat::Muon newM(*m);
    if (selectionMode!=0) {
      PhotonLepMap::const_iterator fsr = theMap.find(m);
      if (fsr!=theMap.end()) {
	if (selectionMode==3) { // Run II: select one per lepton; highest-pT if >4GeV, lowest-DR otherwise
	  PhotonPtr g = selectFSR(fsr->second,m->momentum());
	  PhotonPtrVector gv = {g};	  
	  newM.addUserData("FSRCandidates",gv);
	  allSelFSR.push_back(g);
	} else { //Legacy, etc.: keep all
	  newM.addUserData("FSRCandidates",fsr->second);
	}
      }
    }
    resultMu->push_back(newM);
  }

  //Loop over electrons again to write the result as userData
  for (unsigned int j = 0; j< electronHandle->size(); ++j){
    const pat::Electron* e = &((*electronHandle)[j]);
    //---Clone the pat::Electron
    pat::Electron newE(*e);
    if (selectionMode!=0) {
      PhotonLepMap::const_iterator fsr = theMap.find(e);
      if (fsr!=theMap.end()) {
	if (selectionMode==3) { // Run II: select one per lepton; highest-pT if >4GeV, lowest-DR otherwise
	  PhotonPtr g = selectFSR(fsr->second,e->momentum());
	  PhotonPtrVector gv = {g};	  
	  newE.addUserData("FSRCandidates",gv);
	  allSelFSR.push_back(g);
	} else { //Legacy, etc.: keep all
	  newE.addUserData("FSRCandidates",fsr->second);
	}
      }
    }
    resultEle->push_back(newE);
  }

  //Loop over electrons again to write the result as userData
  for (unsigned int j = 0; j< tleHandle->size(); ++j){
    const pat::Photon* e = &((*tleHandle)[j]);
    //---Clone the pat::Electron
    pat::Photon newE(*e);
    /*if (selectionMode!=0) {
      PhotonLepMap::const_iterator fsr = theMap.find(e);
      if (fsr!=theMap.end()) {
	if (selectionMode==3) { // Run II: select one per lepton; highest-pT if >4GeV, lowest-DR otherwise
	  PhotonPtr g = selectFSR(fsr->second,e->momentum());
	  PhotonPtrVector gv = {g};	  
	  newE.addUserData("FSRCandidates",gv);
	  allSelFSR.push_back(g);
	} else { //Legacy, etc.: keep all
	  newE.addUserData("FSRCandidates",fsr->second);
	}
      }
    }*/
    resultTle->push_back(newE);
  }

  //Loop over electrons again to write the result as userData
/*
  for (unsigned int j = 0; j< looseElectronHandle->size(); ++j){
    const pat::Electron* e = &((*looseElectronHandle)[j]);
    //---Clone the pat::Electron
    pat::Electron newE(*e);
    resultLooseEle->push_back(newE);
  }
*/

  //Recompute isolation of all leptons subtracting FSR from the cone (only for Run II strategy)
  if (selectionMode==3){
    double rhoForMu, rhoForEle;
    {
      edm::Handle<double> rhoHandle;
      iEvent.getByToken(rhoForMuToken, rhoHandle);
      rhoForMu = *rhoHandle;
      iEvent.getByToken(rhoForEleToken, rhoHandle);
      rhoForEle = *rhoHandle;
    }

    for (pat::MuonCollection::iterator m= resultMu->begin(); m!=resultMu->end(); ++m){
      float fsrCorr = 0; // The correction to PFPhotonIso
      for (PhotonPtrVector::const_iterator g = allSelFSR.begin();g!= allSelFSR.end(); ++g) {
	    const pat::PFParticle* gamma = g->get();
	    double dR = ROOT::Math::VectorUtil::DeltaR(gamma->momentum(),m->momentum());
	    // Check if the photon is in the lepton's iso cone and not vetoed
	    if (dR<0.3 && dR > 0.01) {
	      fsrCorr += gamma->pt();
	    }
      } 
      float combRelIsoPFCorr =  LeptonIsoHelper::combRelIsoPF(sampleType, setup, rhoForMu, *m, fsrCorr);
      m->addUserFloat("combRelIsoPFFSRCorr", combRelIsoPFCorr);
      m->addUserFloat("passCombRelIsoPFFSRCorr",combRelIsoPFCorr < muon_iso_cut);
    }

    for (pat::ElectronCollection::iterator e= resultEle->begin(); e!=resultEle->end(); ++e){
      float fsrCorr = 0; // The correction to PFPhotonIso
      for (PhotonPtrVector::const_iterator g = allSelFSR.begin();g!= allSelFSR.end(); ++g) {
	    const pat::PFParticle* gamma = g->get();
	    double dR = ROOT::Math::VectorUtil::DeltaR(gamma->momentum(),e->momentum());
	    // Check if the photon is in the lepton's iso cone and not vetoed
	    if (dR<0.3 && (fabs(e->superCluster()->eta()) < 1.479 || dR > 0.08)) {
	        fsrCorr += gamma->pt();
 	    }
      }
      float combRelIsoPFCorr =  LeptonIsoHelper::combRelIsoPF(sampleType, setup, rhoForEle, *e, fsrCorr);
      e->addUserFloat("combRelIsoPFFSRCorr", combRelIsoPFCorr);
      e->addUserFloat("passCombRelIsoPFFSRCorr",combRelIsoPFCorr < electron_iso_cut);
    }
    for (pat::PhotonCollection::iterator e= resultTle->begin(); e!=resultTle->end(); ++e){
      float fsrCorr = 0; // The correction to PFPhotonIso
      /*
      for (PhotonPtrVector::const_iterator g = allSelFSR.begin();g!= allSelFSR.end(); ++g) {
	    const pat::PFParticle* gamma = g->get();
	    double dR = ROOT::Math::VectorUtil::DeltaR(gamma->momentum(),e->momentum());
	    // Check if the photon is in the lepton's iso cone and not vetoed
	    if (dR<0.3 && (fabs(e->superCluster()->eta()) < 1.479 || dR > 0.08)) {
	        fsrCorr += gamma->pt();
 	    }
      }*/
      float combRelIsoPFCorr = LeptonIsoHelper::combRelIsoPF(sampleType, setup, rhoForEle, *e, fsrCorr);
      e->addUserFloat("combRelIsoPFFSRCorr", combRelIsoPFCorr);
      // Isolation os included in TLE ID
      e->addUserFloat("passCombRelIsoPFFSRCorr",combRelIsoPFCorr < electron_iso_cut); //LeptonIsoHelper::isoCut(&*e)); // FIXME should move this to the .py, once we drop support for the old FSR strategy
    }
//    edm::LogError("") << "About to touch loose of size : " << resultLooseEle->size();

/*
    for (pat::ElectronCollection::iterator e= resultLooseEle->begin(); e != resultLooseEle->end(); ++e){
      float fsrCorr = 0; // The correction to PFPhotonIso
      
      for (PhotonPtrVector::const_iterator g = allSelFSR.begin();g!= allSelFSR.end(); ++g) {
	    const pat::PFParticle* gamma = g->get();
	    double dR = ROOT::Math::VectorUtil::DeltaR(gamma->momentum(),e->momentum());
	    // Check if the photon is in the lepton's iso cone and not vetoed
	    if (dR<0.3 && (fabs(e->superCluster()->eta()) < 1.479 || dR > 0.08)) {
	        fsrCorr += gamma->pt();
 	    }
      }
  //    edm::LogVerbatim("") << "Loose!";
      float combRelIsoPFCorr = LeptonIsoHelper::combRelIsoPF(sampleType, setup, rhoForEle, *e, fsrCorr);
      e->addUserFloat("combRelIsoPFFSRCorr", combRelIsoPFCorr);
      // Isolation os included in TLE ID
      e->addUserFloat("passCombRelIsoPFFSRCorr",combRelIsoPFCorr < electron_iso_cut); //LeptonIsoHelper::isoCut(&*e)); // FIXME should move this to the .py, once we drop support for the old FSR strategy
    }
*/


  }
  
  //Put the result in the event
  iEvent.put(resultMu,"muons");
  iEvent.put(resultEle,"electrons");
  iEvent.put(resultTle,"electronstle");
  //iEvent.put(resultLooseEle,"looseElectrons");

}

PhotonPtr LeptonPhotonMatcher::selectFSR(const PhotonPtrVector& photons, const reco::LeafCandidate::Vector& lepMomentum){ 
  // select one photon per lepton; highest-pT if >4GeV, lowest-DR otherwise
//   PhotonPtr g = *(std::max_element(photons.begin(),photons.end(), [](const PhotonPtr& g1, const PhotonPtr& g2){return g1->pt()<g2->pt();}));
//   if (g->pt()<=4) {
//     g = *(std::min_element(photons.begin(),photons.end(), [lepMomentum](const PhotonPtr& g1, const PhotonPtr& g2){return ROOT::Math::VectorUtil::DeltaR(g1->momentum(),lepMomentum)<ROOT::Math::VectorUtil::DeltaR(g2->momentum(),lepMomentum);}));
//   }

  //Select lowest-DR/ET2
  PhotonPtr g = *(std::min_element(photons.begin(),photons.end(), [lepMomentum](const PhotonPtr& g1, const PhotonPtr& g2){return  (ROOT::Math::VectorUtil::DeltaR(g1->momentum(),lepMomentum)/g1->pt()/g1->pt())<(ROOT::Math::VectorUtil::DeltaR(g2->momentum(),lepMomentum)/g2->pt()*g2->pt());}));

  // Select highest-ET
//   PhotonPtr g = *(std::max_element(photons.begin(),photons.end(), [](const PhotonPtr& g1, const PhotonPtr& g2){return g1->pt()<g2->pt();}));
  return g;
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(LeptonPhotonMatcher);

