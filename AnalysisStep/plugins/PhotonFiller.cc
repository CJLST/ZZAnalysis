/** \class PhotonFiller
 *
 *  Preselect photons from all packedPFCandidates, 
 *  SuperCluster-based veto with electrons and store
 *  the surviving photons in the event for further use
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

#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/PFParticle.h>
#include <ZZAnalysis/AnalysisStep/interface/PhotonFwd.h>
//#include <EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h>
//#include <DataFormats/GeometryVector/interface/VectorUtil.h> 
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidate.h>

#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>

//#include <ZZAnalysis/AnalysisStep/interface/MCHistoryTools.h>
#include <Math/VectorUtil.h>
#include <TMath.h>

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;

class PhotonFiller : public edm::EDProducer {
 public:
  /// Constructor
  explicit PhotonFiller(const edm::ParameterSet&);
    
  /// Destructor
  ~PhotonFiller(){};  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  edm::EDGetTokenT<pat::ElectronCollection> electronToken;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate> > pfCandToken;
  int selectionMode;
  int sampleType;
  int setup;
  bool debug;
};


PhotonFiller::PhotonFiller(const edm::ParameterSet& iConfig) :
  electronToken(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronSrc"))),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  debug(iConfig.getUntrackedParameter<bool>("debug",false))
{
  pfCandToken = consumes<edm::View<pat::PackedCandidate> >(edm::InputTag("packedPFCandidates"));

  string mode = iConfig.getParameter<string>("photonSel");
  
  if      (mode == "skip")        selectionMode = 0; // no FSR
  else if (mode == "passThrough") selectionMode = 1; // for debug
  else if (mode == "Legacy")      selectionMode = 2;
  else if (mode == "RunII")       selectionMode = 3;
  else {
    cout << "PhotonFiller: mode " << mode << " not supported" << endl;
    abort();
  }
  
  produces<reco::PFCandidateCollection>(); 
}


void
PhotonFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //  edm::Handle<pat::ElectronRefVector> electronHandle;
  edm::Handle<pat::ElectronCollection> electronHandle;
  iEvent.getByToken(electronToken, electronHandle);

  //--- Get the PF cands
  edm::Handle<edm::View<pat::PackedCandidate> > pfCands; 
  iEvent.getByToken(pfCandToken, pfCands);

  // Output collections
  std::auto_ptr<reco::PFCandidateCollection> result( new reco::PFCandidateCollection );

  if (selectionMode!=0) {
    //----------------------
    // Loop on photons
    //----------------------
    for (unsigned int i=0;i<pfCands->size();++i) {
      
      // Get the candidate as edm::Ptr
      //edm::Ptr<pat::PackedCandidate> g = pfCands->ptrAt(i);
      edm::Ptr<pat::PackedCandidate> g = pfCands->ptrAt(i);
      
      // We only want photons
      if (g->pdgId()!=22) continue;

      // Photon preselection (is currently already applied on pat::PackedCandidate collection)
      if (!(g->pt()>2. && fabs(g->eta())<2.4)) continue;

      //---------------------
      // // Supercluster veto
      //---------------------
      bool SCVeto=false;
      
      if (electronHandle->size()>0) {
        for (unsigned int j = 0; j< electronHandle->size(); ++j){
          const pat::Electron* e = &((*electronHandle)[j]);
          if (e->userFloat("isSIP")){
            if (setup==2016) {
              if ((e->associatedPackedPFCandidates()).size()) {
    	        edm::RefVector < pat::PackedCandidateCollection > pfcands = e->associatedPackedPFCandidates();
    	        for ( auto itr: pfcands ) {
                  if (g.get()==&(*itr)) {
      	            SCVeto=true; 
      	            if (debug) cout << "SC veto: " << itr->eta() << " " << itr->phi() << "   " 
      	                	    << fabs(g->eta() - itr->eta()) << " " << reco::deltaPhi(g->phi(), itr->phi()) <<endl;
      	            break;
      	          }  
                }
              }
            } else {
      	      double dR = reco::deltaR(*(e->superCluster()), *g);
              if ((fabs(g->eta() - e->superCluster()->eta())<0.05 && fabs(reco::deltaPhi(g->phi(), e->superCluster()->phi()))<2.) || dR<0.15) {
    	        SCVeto=true;
    	        if (debug) cout << "SC veto: "<< g->pt() << " " << e->pt() << " " << dR << " "
    	  		        << fabs(g->eta() - e->superCluster()->eta()) << " " << reco::deltaPhi(g->phi(), e->superCluster()->phi()) <<endl;
    	        break;
    	      } 
    	    }
    	  }  
        }
      }
      if (debug) cout << "PhotonFiller: gamma:" << g->pt() << " " << g->eta() << " " << g->phi() << " SCVeto: " << SCVeto << endl;

      if (SCVeto) continue;
      
      result->push_back(reco::PFCandidate(0, g->p4(), reco::PFCandidate::gamma));
      result->back().setStatus(0);

    } // end of loop over photon collection
  }
  
  //Put the result in the event
  iEvent.put(result);
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(PhotonFiller);

