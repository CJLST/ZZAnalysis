#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>

#include <DataFormats/PatCandidates/interface/Muon.h>

// Kalman Muon Corrections 
#include "KaMuCa/Calibration/interface/KalmanMuonCalibrator.h"

#include "TLorentzVector.h"

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;


class KalmanPATMuonCorrector : public edm::EDProducer {
 public:
  /// Constructor
  explicit KalmanPATMuonCorrector(const edm::ParameterSet&);
    
  /// Destructor
  ~KalmanPATMuonCorrector(){
    delete calibrator;
  };  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  edm::EDGetTokenT<vector<pat::Muon> > muonToken_;
  string identifier_;
  bool isMC_;
  bool isSync_;

  KalmanMuonCalibrator* calibrator;

};


KalmanPATMuonCorrector::KalmanPATMuonCorrector(const edm::ParameterSet& iConfig) :
  muonToken_(consumes<vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("src"))),
  identifier_(iConfig.getParameter<string>("identifier")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  isSync_(iConfig.getParameter<bool>("isSynchronization")),
  calibrator(0)
{
  if (identifier_!="None") calibrator = new KalmanMuonCalibrator(identifier_);

  produces<pat::MuonCollection>();
}


void
KalmanPATMuonCorrector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Input collection
  edm::Handle<vector<pat::Muon> > muonHandle;
  iEvent.getByToken(muonToken_, muonHandle);
  const vector<pat::Muon> *inputMuons = muonHandle.product();

  // Output collection
  auto outputMuons = std::make_unique<vector<pat::Muon> >();

  TLorentzVector p4;

  for (unsigned i=0; i<inputMuons->size(); ++i) {

    pat::Muon mu = inputMuons->at(i);

    double oldpt=mu.pt();
    double newpt=oldpt;
    double oldpterr=mu.muonBestTrack()->ptError();
    double newpterr=oldpterr;

    if (calibrator != 0 && mu.muonBestTrackType() == 1 && oldpt<=200.) { //skip correction for passthrough mode, or if muon pt does not come from InnerTrack, or for muons above 200 GeV
      if(isMC_){
	/// ====== ON MC (correction plus smearing) =====
	double corrPt = calibrator->getCorrectedPt(oldpt, mu.eta(), mu.phi(), mu.charge());
	if(!isSync_) {
	  newpt = calibrator->smear(corrPt, mu.eta());
	} else {
	  newpt = calibrator->smearForSync(corrPt, mu.eta());
	}
	newpterr = newpt * calibrator->getCorrectedError(newpt, mu.eta(), oldpterr/newpt);
      } else {
	/// ====== ON DATA (correction only) =====
	if(mu.pt()>2.0 && fabs(mu.eta())<2.4){
	  newpt = calibrator->getCorrectedPt(oldpt, mu.eta(), mu.phi(), mu.charge());
	  newpterr = newpt * calibrator->getCorrectedError(newpt, mu.eta(), oldpterr/newpt);
	} else {
	  // keep old values
	}
      }
    }

    p4.SetPtEtaPhiM(newpt, mu.eta(), mu.phi(), mu.mass());
    mu.setP4(reco::Particle::PolarLorentzVector(p4.Pt(), p4.Eta(), p4.Phi(), mu.mass()));
    mu.addUserFloat("correctedPtError",newpterr);
    
    outputMuons->push_back(mu);
  }

  iEvent.put(std::move(outputMuons));

}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(KalmanPATMuonCorrector);

