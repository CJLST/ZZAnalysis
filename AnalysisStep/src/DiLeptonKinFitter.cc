#include "ZZAnalysis/AnalysisStep/interface/DiLeptonKinFitter.h"
#include "PhysicsTools/KinFitter/interface/TFitConstraintM.h"
#include "PhysicsTools/KinFitter/interface/TFitConstraintMGaus.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleMCPInvSpher.h" //for muons
#include "PhysicsTools/KinFitter/interface/TFitParticleMCSpher.h" //for electrons 

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TString.h>
#include <iostream>

using namespace std;
typedef reco::Candidate::LorentzVector LorentzVector;


DiLeptonKinFitter::DiLeptonKinFitter( const TString&  name, const TString& title, double mass ) : TKinFitter( name, title ) {

  name_ = name;
  title_ = title;
  mass_ = mass;  
}


std::pair<TLorentzVector,TLorentzVector> DiLeptonKinFitter::fit(const reco::Candidate* Z1) {

  this->reset();

  //get the Z1 leptons TLorentzVector
  const reco::Candidate* Z1Lp= Z1->daughter(0);	
  const reco::Candidate* Z1Lm= Z1->daughter(1);
  int id11 = Z1Lp->pdgId();
  int id12 = Z1Lm->pdgId();
  math::XYZTLorentzVector p11 = Z1Lp->p4();
  math::XYZTLorentzVector p12 = Z1Lm->p4();

  
  //If there are fsr photons attached to the Z1 leptons I sum them to the p4 to use in the refit (first approssimation)
  int d0FSR = (static_cast<const pat::CompositeCandidate*>(Z1->masterClone().get()))->userFloat("dauWithFSR"); // FIXME: obsolete, need to move to bew FSR coding scheme
  if (d0FSR>=0) {
    if (Z1->numberOfDaughters()!=3) cout << "ERROR: ZZCandidateFiller: problem in FSR" << endl;
    if (d0FSR==0) {
      p11 = p11 + Z1->daughter(2)->p4();
    } else if (d0FSR==1){
      p12 = p12 + Z1->daughter(2)->p4();  
    }
  }
    
  TLorentzVector pL11(p11.x(),p11.y(),p11.z(),p11.t());
  TLorentzVector pL12(p12.x(),p12.y(),p12.z(),p12.t());
  TVector3 p3L11(p11.x(),p11.y(),p11.z());
  TVector3 p3L12(p12.x(),p12.y(),p12.z());
  double pL11m = pL11.M();
  double pL12m = pL12.M();

  std::pair<TLorentzVector,TLorentzVector> lepfitpair;
  //if muons
  if (fabs(id11)==13 && fabs(id12)==13){

    const reco::Muon* Z1mup = static_cast<const reco::Muon*>(Z1Lp->masterClone().get());
    const reco::Muon* Z1mum = static_cast<const reco::Muon*>(Z1Lm->masterClone().get());

    TMatrixD m_mu1(3,3);
    TMatrixD m_mu2(3,3);
    
    m_mu1(0,0) = this->ErrMuP(Z1mup);      // 1/p
    m_mu1(1,1) = this->ErrMuTheta(Z1mup);  // theta
    m_mu1(2,2) = this->ErrMuPhi(Z1mup);    // phi
    m_mu2(0,0) = this->ErrMuP(Z1mum);      // 1/p
    m_mu2(1,1) = this->ErrMuTheta(Z1mum);  // theta
    m_mu2(2,2) = this->ErrMuPhi(Z1mum);    // phi
    
    TFitParticleMCPInvSpher *fitMu1 = new TFitParticleMCPInvSpher( "Mu1", "Mu1", &p3L11, pL11m, &m_mu1 );
    TFitParticleMCPInvSpher *fitMu2 = new TFitParticleMCPInvSpher( "Mu2", "Mu2", &p3L12, pL12m, &m_mu2 );
    
    TFitConstraintMGaus *mCons_muons = new TFitConstraintMGaus( "ZMassConstraint_muons", "ZMass-Constraint_muons", 0, 0 , mass_, 2.4952/2.355 );
    mCons_muons->addParticles1( fitMu1, fitMu2 );
    
    
    this->addMeasParticle( fitMu1 );
    this->addMeasParticle( fitMu2 );
    this->addConstraint( mCons_muons );
    

    //Set convergence criteria
    this->setMaxNbIter( 30 );
    this->setMaxDeltaS( 1e-2 );
    this->setMaxF( 1e-1 );
    this->setVerbosity(0);
    
    //Perform the fit
    TKinFitter::fit();
    
    TLorentzVector mu1_kinfit(*fitMu1->getCurr4Vec());
    TLorentzVector mu2_kinfit(*fitMu2->getCurr4Vec());
    
    
    lepfitpair.first  = mu1_kinfit;
    lepfitpair.second = mu2_kinfit;
    
    delete fitMu1; delete fitMu2; delete mCons_muons;
    
  }else{ //if electrons
    //not implemented for the moment
    
    lepfitpair.first  = pL11;
    lepfitpair.second = pL12;
  }
  
  return lepfitpair;
  
}


double DiLeptonKinFitter::getChi2( ) {
  
  //Get the Chi2
  double chi2 =  TKinFitter::getS();
  
  return chi2;
}

double DiLeptonKinFitter::getChi2Prob( ) {
  
  //Get the Chi2 and ndf
  double chi2Prob = TMath::Prob(TKinFitter::getS(),TKinFitter::getNDF());
  return chi2Prob;
}


int DiLeptonKinFitter::getStatus( ) {
  
  //Get the Chi2
  double status =  TKinFitter::getStatus();
  
  return status;
}


double DiLeptonKinFitter::ErrMuP(const reco::Muon* Z1l) {
  
  double err2;
  err2 = Z1l->innerTrack()->qoverpError();
  err2 *= err2;
  //here the error can be rescaled with Gio/Mingshui's map
  return err2;

}


double DiLeptonKinFitter::ErrMuTheta(const reco::Muon* Z1l) {
  
  double err2;
  err2 = Z1l->innerTrack()->thetaError();
  err2 *= err2;
  return err2;

}


double DiLeptonKinFitter::ErrMuPhi(const reco::Muon* Z1l) {
  
  double err2;
  err2 = Z1l->innerTrack()->phiError();
  err2 *= err2;
  return err2;
  
}






