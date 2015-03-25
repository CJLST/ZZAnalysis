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
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>

#include <ZZAnalysis/AnalysisStep/interface/MCHistoryTools.h>
#include <Math/VectorUtil.h>
#include "TMath.h"

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

  const edm::InputTag theMuonTag_;
  const edm::InputTag theElectronTag_;
  const edm::InputTag thePhotonTag_;
  bool matchFSR;
  bool debug;
};


LeptonPhotonMatcher::LeptonPhotonMatcher(const edm::ParameterSet& iConfig) :
  theMuonTag_(iConfig.getParameter<InputTag>("muonSrc")),
  theElectronTag_(iConfig.getParameter<InputTag>("electronSrc")),
  thePhotonTag_(iConfig.getParameter<InputTag>("photonSrc")),
  matchFSR(iConfig.getParameter<bool>("matchFSR")),
  debug(iConfig.getUntrackedParameter<bool>("debug",false))
{
  produces<pat::MuonCollection>("muons");
  produces<pat::ElectronCollection>("electrons");
}


void
LeptonPhotonMatcher::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //--- Get leptons and rho
  //  edm::Handle<pat::MuonRefVector> muonHandle;
  edm::Handle<pat::MuonCollection> muonHandle;
  iEvent.getByLabel(theMuonTag_, muonHandle);

  //  edm::Handle<pat::ElectronRefVector> electronHandle;
  edm::Handle<pat::ElectronCollection> electronHandle;
  iEvent.getByLabel(theElectronTag_, electronHandle);

  //--- Get the photons
  edm::Handle<edm::View<pat::PFParticle> > photonHandle;
  iEvent.getByLabel(thePhotonTag_, photonHandle);


  // Output collections
  auto_ptr<pat::MuonCollection> resultMu( new pat::MuonCollection() );
  auto_ptr<pat::ElectronCollection> resultEle( new pat::ElectronCollection() );

  // Associate a vector of Ptr<Photon> to lepton pointers
  typedef map<const reco::Candidate*, PhotonPtrVector> PhotonLepMap;
  PhotonLepMap theMap;


  if (matchFSR) {
    //----------------------
    // Loop on photons
    //----------------------
    for (unsigned int i=0;i<photonHandle->size();++i) {
      if (muonHandle->size()+electronHandle->size()==0) continue;
      
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
	if (e->userFloat("isGood")){
	  double dR = ROOT::Math::VectorUtil::DeltaR(e->momentum(),g->momentum());
	  if ((fabs(g->eta() - e->eta())<0.05 && fabs(reco::deltaPhi(g->phi(), e->phi()))<2.) || dR<0.15) {	    
	    SCVeto=true;
  	    if (debug) cout << "SC veto: "<< g->pt() << " " << e->pt() << " " << dR << " "
			    << fabs(g->eta() - e->superCluster()->eta()) << " " << reco::deltaPhi(g->phi(), e->phi()) <<endl;
	    break;
	  } 
	}
      }
      if (debug) cout << "GAMMA: " << g->pt() << " " << g->eta() << " " << g->phi() << " SCVeto: " << SCVeto << endl;
      if (SCVeto) continue;


      //------------------------------------------------------
      // Get the closest lepton among those satisfying ID, SIP
      //------------------------------------------------------
      double dRMin(10e9);
      const reco::Candidate* closestLep = 0;

      // Loop over pat::Muon
      for (unsigned int j = 0; j< muonHandle->size(); ++j){
	//      const pat::Muon* m = ((*muonHandle)[j]).get();
	const pat::Muon* m = &((*muonHandle)[j]);
	if (! m->userFloat("isGood")) continue; 
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
	if ( ! e->userFloat("isGood")) continue;
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
	double gRelIso = 999.;
	if( dRMin<0.07 ){
	  if (g->pt()>2.) accept = true;
	} else if ( dRMin<0.5 ){ // That's implicit, but does not hurt
	  // double relIso = g->relIso(0.5); // This is buggy, needs to recompute it.
	  gRelIso = (g->userFloat("fsrPhotonPFIsoChHadPUNoPU03pt02") + g->userFloat("fsrPhotonPFIsoNHadPhoton03")) / g->pt();
	  if (g->pt()>4 && gRelIso<1.) accept = true;
	}
	if(debug) cout << "   " << "   closest lep: " << closestLep->pdgId() << " " << closestLep->pt() << " gRelIso: " << gRelIso << " dRMin: " << dRMin << " accept: " << accept << endl;
	if (accept) theMap[closestLep].push_back(g);
      }
    }// loop over photon collection
  }
  
  // Loop over muons again to write the result as userData
  for (unsigned int j = 0; j< muonHandle->size(); ++j){
    //    const pat::Muon* m = ((*muonHandle)[j]).get(); // Pointer to original mu
    const pat::Muon* m = &((*muonHandle)[j]);
    //---Clone the pat::Muon
    pat::Muon newM(*m);
    if (matchFSR) {
      PhotonLepMap::const_iterator fsr = theMap.find(m);
      if (fsr!=theMap.end()) {
	newM.addUserData("FSRCandidates",fsr->second);
      }
    }
    resultMu->push_back(newM);
  }

  //Loop over electrons again to write the result as userData
  for (unsigned int j = 0; j< electronHandle->size(); ++j){
    //const pat::Electron* e = ((*electronHandle)[j]).get();
    const pat::Electron* e = &((*electronHandle)[j]);
    //---Clone the pat::Electron
    pat::Electron newE(*e);
    if (matchFSR) {
      PhotonLepMap::const_iterator fsr = theMap.find(e);
      if (fsr!=theMap.end()) {
	newE.addUserData("FSRCandidates",fsr->second);
      }
    }
    resultEle->push_back(newE);      
  }
  
  //Put the result in the event
  iEvent.put(resultMu,"muons");
  iEvent.put(resultEle,"electrons");
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(LeptonPhotonMatcher);

