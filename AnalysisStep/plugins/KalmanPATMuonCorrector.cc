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
  isSync_(iConfig.getParameter<bool>("isSynchronization"))
{
  calibrator = new KalmanMuonCalibrator(identifier_);

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
  vector<pat::Muon> *outputMuons = new vector<pat::Muon>();

  TLorentzVector p4;

  for (unsigned i=0; i<inputMuons->size(); ++i) {

    pat::Muon mu = inputMuons->at(i);

    double newPt, newPtError;

    if(isMC_){
      /// ====== ON MC (correction plus smearing) =====

      double corrPt = calibrator->getCorrectedPt(mu.pt(), mu.eta(), mu.phi(), mu.charge());
      double corrPtError = corrPt * calibrator->getCorrectedError(corrPt, mu.eta(), mu.bestTrack()->ptError()/corrPt );

      if(!isSync_) 
	newPt = calibrator->smear(corrPt, mu.eta());
      else
	newPt = calibrator->smearForSync(corrPt, mu.eta());
      newPtError = newPt * calibrator->getCorrectedErrorAfterSmearing(newPt, mu.eta(), corrPtError / newPt );

    }else{
      /// ====== ON DATA (correction only) =====

      newPt = calibrator->getCorrectedPt(mu.pt(), mu.eta(), mu.phi(), mu.charge());
      newPtError = newPt * calibrator->getCorrectedError(newPt, mu.eta(), mu.bestTrack()->ptError()/newPt );

    }

    p4.SetPtEtaPhiM(newPt, mu.eta(), mu.phi(), mu.mass());
    mu.setP4(reco::Particle::PolarLorentzVector(p4.Pt(), p4.Eta(), p4.Phi(), mu.mass()));
    mu.addUserFloat("correctedPtError",newPtError);
    
    outputMuons->push_back(mu);
  }

  std::auto_ptr<std::vector<pat::Muon> > ptr(outputMuons);
  iEvent.put(ptr);

}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(KalmanPATMuonCorrector);

