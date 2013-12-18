#ifndef DiLeptonKinFitter_h
#define DiLeptonKinFitter_h

#include "PhysicsTools/KinFitter/interface/TKinFitter.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>

#include <TLorentzVector.h>

class DiLeptonKinFitter : public TKinFitter {

 public:
  
  DiLeptonKinFitter( const TString&  name, const TString& title, double mass = 91.1876 );

  void set_mass( double mass ) { mass_ = mass; };

  std::pair<TLorentzVector,TLorentzVector> fit(const reco::Candidate* Z1);
  double getChi2();
  int getStatus();
  double getChi2Prob();

  double ErrMuP(const reco::Muon* Z1l); 
  double ErrMuTheta(const reco::Muon* Z1l); 
  double ErrMuPhi(const reco::Muon* Z1l); 
//   double ErrEleP(const reco::Candidate* Z1l); 
//   double ErrEleTheta(const reco::Candidate* Z1l); 
//   double ErrElePhi(const reco::Candidate* Z1l); 


 private:

  double mass_;
  TString name_;
  TString title_;

};

#endif
