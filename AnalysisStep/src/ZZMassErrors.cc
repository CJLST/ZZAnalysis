#include <ZZAnalysis/AnalysisStep/interface/ZZMassErrors.h>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
// #include "TrackingTools/Records/interface/TransientTrackRecord.h"
// #include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianCurvilinearToCartesian.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/CartesianTrajectoryError.h"

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>

#include <cmath>


using namespace reco;
using namespace std;


namespace {
  // This code from RecoEgamma/EgammaElectronAlgos/src/ElectronEnergyCorrector.cc,  revision=1.21
  float energyError( float E, float * par ) { 
    return sqrt( pow(par[0]/sqrt(E),2) + pow(par[1]/E,2) + pow(par[2],2) ) ; 
  }

  
  float simpleParameterizationUncertainty(const pat::Electron* electron) {
    float error = 999. ;

#if CMSSW_VERSION<500
    double ecalEnergy = electron->ecalEnergy() ;
#else
    double ecalEnergy = electron->correctedEcalEnergy() ;
#endif

    if (electron->isEB()) {
      float parEB[3] = { 5.24e-02,  2.01e-01, 1.00e-02} ;
      error =  ecalEnergy*energyError(ecalEnergy,parEB) ;
    } else if (electron->isEE()) {
      float parEE[3] = { 1.46e-01, 9.21e-01, 1.94e-03} ;
      error =  ecalEnergy*energyError(ecalEnergy,parEE) ;
    }
    return error;
  }
}

ZZMassErrors::ZZMassErrors(const edm::EventSetup& eventSetup) {
  eventSetup.get<IdealMagneticFieldRecord>().get(magfield);  
}


double ZZMassErrors::pError(const pat::Electron* ele) {
  double dp=0;
  if (! ele->ecalDriven()) {
    dp = simpleParameterizationUncertainty(ele);
  } else {
    dp = ele->p4Error(GsfElectron::P4_COMBINATION);
  }
  return dp;
}


float ZZMassErrors::massError(const reco::Candidate* Z1Lp, 
			      const reco::Candidate* Z1Lm, 
			      const reco::Candidate* Z2Lp, 
			      const reco::Candidate* Z2Lm) {
  
  //bool print=true;
  //bool debug=true;
  bool debug=false;
  
  // 3vector for the momentum of each lepton
  // use physics vector TVector3 here to get easy access to the angle between two vectors
  TVector3 Lepton1(Z1Lp->px(), Z1Lp->py(), Z1Lp->pz());
  TVector3 Lepton2(Z1Lm->px(), Z1Lm->py(), Z1Lm->pz());
  TVector3 Lepton3(Z2Lp->px(), Z2Lp->py(), Z2Lp->pz());
  TVector3 Lepton4(Z2Lm->px(), Z2Lm->py(), Z2Lm->pz());
  
  double p1 = Lepton1.Mag();
  double p2 = Lepton2.Mag();
  double p3 = Lepton3.Mag();
  double p4 = Lepton4.Mag();
  
  
  // 3x3 error matrices for the 3momentum error of each lepton
  TMatrixDSym Lepton1Error(3);
  TMatrixDSym Lepton2Error(3);
  TMatrixDSym Lepton3Error(3);
  TMatrixDSym Lepton4Error(3);
  
  const reco::Candidate* cands[4] = {Z1Lp,Z1Lm,Z2Lp,Z2Lm};
  TMatrixDSym* errors[4] = {&Lepton1Error,&Lepton2Error,&Lepton3Error,&Lepton4Error};
  
  for (int i = 0; i<4; ++i) {
    double dp = 0.;
    const pat::Muon* mu = dynamic_cast<const pat::Muon*>(cands[i]->masterClone().get());
    const pat::Electron* ele = dynamic_cast<const pat::Electron*>(cands[i]->masterClone().get());
    if (mu!=0) {
      const Track* t = mu->track().get();
      
      //     TransientTrack tt = trackBuilder->build(t);
      //     const CartesianTrajectoryError& cartErr = tt.initialFreeState().cartesianError();
      
      GlobalTrajectoryParameters gp(GlobalPoint(t->vx(), t->vy(),  t->vz()),
				    GlobalVector(t->px(),t->py(),t->pz()),
				    t->charge(),
				    magfield.product());
      JacobianCurvilinearToCartesian curv2cart(gp);
      CartesianTrajectoryError cartErr= ROOT::Math::Similarity(curv2cart.jacobian(), t->covariance());
      const AlgebraicSymMatrix66 & m = cartErr.matrix();
      
      for (int j=0; j<3; ++j){
	for (int k=0; k<3; ++k){
	  (*(errors[i]))(j,k) = m[j+3][k+3];
	}
      }
    } else {
      dp = pError(ele);
      
      double p = ele->p();
      (*(errors[i]))(0,0) = pow(dp*ele->px()/p,2);
      (*(errors[i]))(0,1) = pow(dp*ele->px()/p,2)*ele->py()/ele->px() ;
      (*(errors[i]))(0,2) = pow(dp*ele->px()/p,2)*ele->pz()/ele->px() ;
      (*(errors[i]))(1,0) = (*(errors[i]))(0,1);
      (*(errors[i]))(1,1) = pow(dp*ele->py()/p,2) ;
      (*(errors[i]))(1,2) = pow(dp*ele->py()/p,2)*ele->pz()/ele->py() ;
      (*(errors[i]))(2,0) = (*(errors[i]))(0,2);
      (*(errors[i]))(2,1) = (*(errors[i]))(1,2);
      (*(errors[i]))(2,2) = pow(dp*ele->pz()/p,2) ;
    }
    
//     if (debug) {
//       printParticle();
//       printParticle(cands[i],i);
//       if (ele!=0) cout << " ele p= " << ele->p() << " isEB: " << ele->isEB() << " classification: " << ele->classification() 
// 		       << " ecalDriven: " << ele->ecalDriven() << " dp: " << dp
// 		       << " ecalEnergyError: " << ele->ecalEnergyError()  << " p4Error: " << ele->p4Error(GsfElectron::P4_COMBINATION) << endl;
//       (*(errors[i])).Print();      
//     }
  }
  
  // 12x12 error matrix of the 4l system
  TMatrixDSym FourLeptonError(12); // this initializes to 0
  FourLeptonError.SetSub(0,0,Lepton1Error);	
  FourLeptonError.SetSub(3,3,Lepton2Error);	
  FourLeptonError.SetSub(6,6,Lepton3Error);	
  FourLeptonError.SetSub(9,9,Lepton4Error);	
  
  // Jacobian matrix for d(m^2)
  TMatrixD Jacobian(1,12);
  
  // m^2 = sum_i,j [ p_i p_j ( 1 - cos (theta(i,j)) ) ] 
  double invMassSquared;
  
  invMassSquared  = 2.*p1*p2*(1-cos(Lepton1.Angle(Lepton2))) + 2.*p1*p3*(1-cos(Lepton1.Angle(Lepton3))) + 2.*p1*p4*(1-cos(Lepton1.Angle(Lepton4)));
  invMassSquared = invMassSquared + 2.*p2*p3*(1-cos(Lepton2.Angle(Lepton3))) + 2.*p2*p4*(1-cos(Lepton2.Angle(Lepton4))); 
  invMassSquared = invMassSquared + 2.*p3*p4*(1-cos(Lepton3.Angle(Lepton4))); 
  
  // derivatives of p1,p2,p3,p4 wrt px, py and pz 
  double dp1dp1x = Lepton1.x()/p1, dp2dp2x = Lepton2.x()/p2, dp3dp3x = Lepton3.x()/p3, dp4dp4x = Lepton4.x()/p4;
  double dp1dp1y = Lepton1.y()/p1, dp2dp2y = Lepton2.y()/p2, dp3dp3y = Lepton3.y()/p3, dp4dp4y = Lepton4.y()/p4;
  double dp1dp1z = Lepton1.z()/p1, dp2dp2z = Lepton2.z()/p2, dp3dp3z = Lepton3.z()/p3, dp4dp4z = Lepton4.z()/p4;
  
  // derivatives of cos(pi,pj) wrt px, py and pz
  double dcos12dp1x = (Lepton2.x()*p1*p2 - Lepton1.x()*p2*p2*cos(Lepton1.Angle(Lepton2)))/(p1*p1*p2*p2) ; 							
  double dcos13dp1x = (Lepton3.x()*p1*p3 - Lepton1.x()*p3*p3*cos(Lepton1.Angle(Lepton3)))/(p1*p1*p3*p3) ;							
  double dcos14dp1x = (Lepton4.x()*p1*p4 - Lepton1.x()*p4*p4*cos(Lepton1.Angle(Lepton4)))/(p1*p1*p4*p4) ;							
  double dcos23dp2x = (Lepton3.x()*p2*p3 - Lepton2.x()*p3*p3*cos(Lepton2.Angle(Lepton3)))/(p2*p2*p3*p3) ;							
  double dcos24dp2x = (Lepton4.x()*p2*p4 - Lepton2.x()*p4*p4*cos(Lepton2.Angle(Lepton4)))/(p2*p2*p4*p4) ;							
  double dcos34dp3x = (Lepton4.x()*p3*p4 - Lepton3.x()*p4*p4*cos(Lepton3.Angle(Lepton4)))/(p3*p3*p4*p4) ;							
  double dcos12dp1y = (Lepton2.y()*p1*p2 - Lepton1.y()*p2*p2*cos(Lepton1.Angle(Lepton2)))/(p1*p1*p2*p2) ;							
  double dcos13dp1y = (Lepton3.y()*p1*p3 - Lepton1.y()*p3*p3*cos(Lepton1.Angle(Lepton3)))/(p1*p1*p3*p3) ;							
  double dcos14dp1y = (Lepton4.y()*p1*p4 - Lepton1.y()*p4*p4*cos(Lepton1.Angle(Lepton4)))/(p1*p1*p4*p4) ;							
  double dcos23dp2y = (Lepton3.y()*p2*p3 - Lepton2.y()*p3*p3*cos(Lepton2.Angle(Lepton3)))/(p2*p2*p3*p3) ;							
  double dcos24dp2y = (Lepton4.y()*p2*p4 - Lepton2.y()*p4*p4*cos(Lepton2.Angle(Lepton4)))/(p2*p2*p4*p4) ;							
  double dcos34dp3y = (Lepton4.y()*p3*p4 - Lepton3.y()*p4*p4*cos(Lepton3.Angle(Lepton4)))/(p3*p3*p4*p4) ;							
  double dcos12dp1z = (Lepton2.z()*p1*p2 - Lepton1.z()*p2*p2*cos(Lepton1.Angle(Lepton2)))/(p1*p1*p2*p2) ;							
  double dcos13dp1z = (Lepton3.z()*p1*p3 - Lepton1.z()*p3*p3*cos(Lepton1.Angle(Lepton3)))/(p1*p1*p3*p3) ;							
  double dcos14dp1z = (Lepton4.z()*p1*p4 - Lepton1.z()*p4*p4*cos(Lepton1.Angle(Lepton4)))/(p1*p1*p4*p4) ;							
  double dcos23dp2z = (Lepton3.z()*p2*p3 - Lepton2.z()*p3*p3*cos(Lepton2.Angle(Lepton3)))/(p2*p2*p3*p3) ;							
  double dcos24dp2z = (Lepton4.z()*p2*p4 - Lepton2.z()*p4*p4*cos(Lepton2.Angle(Lepton4)))/(p2*p2*p4*p4) ;							
  double dcos34dp3z = (Lepton4.z()*p3*p4 - Lepton3.z()*p4*p4*cos(Lepton3.Angle(Lepton4)))/(p3*p3*p4*p4) ;
  double dcos12dp2x = (Lepton1.x()*p1*p2 - Lepton2.x()*p1*p1*cos(Lepton1.Angle(Lepton2)))/(p1*p1*p2*p2) ; 							
  double dcos13dp3x = (Lepton1.x()*p1*p3 - Lepton3.x()*p1*p1*cos(Lepton1.Angle(Lepton3)))/(p1*p1*p3*p3) ;							
  double dcos14dp4x = (Lepton1.x()*p1*p4 - Lepton4.x()*p1*p1*cos(Lepton1.Angle(Lepton4)))/(p1*p1*p4*p4) ;							
  double dcos23dp3x = (Lepton2.x()*p2*p3 - Lepton3.x()*p2*p2*cos(Lepton2.Angle(Lepton3)))/(p2*p2*p3*p3) ;							
  double dcos24dp4x = (Lepton2.x()*p2*p4 - Lepton4.x()*p2*p2*cos(Lepton2.Angle(Lepton4)))/(p2*p2*p4*p4) ;							
  double dcos34dp4x = (Lepton3.x()*p3*p4 - Lepton4.x()*p3*p3*cos(Lepton3.Angle(Lepton4)))/(p3*p3*p4*p4) ;							
  double dcos12dp2y = (Lepton1.y()*p1*p2 - Lepton2.y()*p1*p1*cos(Lepton1.Angle(Lepton2)))/(p1*p1*p2*p2) ;							
  double dcos13dp3y = (Lepton1.y()*p1*p3 - Lepton3.y()*p1*p1*cos(Lepton1.Angle(Lepton3)))/(p1*p1*p3*p3) ;							
  double dcos14dp4y = (Lepton1.y()*p1*p4 - Lepton4.y()*p1*p1*cos(Lepton1.Angle(Lepton4)))/(p1*p1*p4*p4) ;							
  double dcos23dp3y = (Lepton2.y()*p2*p3 - Lepton3.y()*p2*p2*cos(Lepton2.Angle(Lepton3)))/(p2*p2*p3*p3) ;							
  double dcos24dp4y = (Lepton2.y()*p2*p4 - Lepton4.y()*p2*p2*cos(Lepton2.Angle(Lepton4)))/(p2*p2*p4*p4) ;							
  double dcos34dp4y = (Lepton3.y()*p3*p4 - Lepton4.y()*p3*p3*cos(Lepton3.Angle(Lepton4)))/(p3*p3*p4*p4) ;							
  double dcos12dp2z = (Lepton1.z()*p1*p2 - Lepton2.z()*p1*p1*cos(Lepton1.Angle(Lepton2)))/(p1*p1*p2*p2) ;							
  double dcos13dp3z = (Lepton1.z()*p1*p3 - Lepton3.z()*p1*p1*cos(Lepton1.Angle(Lepton3)))/(p1*p1*p3*p3) ;							
  double dcos14dp4z = (Lepton1.z()*p1*p4 - Lepton4.z()*p1*p1*cos(Lepton1.Angle(Lepton4)))/(p1*p1*p4*p4) ;							
  double dcos23dp3z = (Lepton2.z()*p2*p3 - Lepton3.z()*p2*p2*cos(Lepton2.Angle(Lepton3)))/(p2*p2*p3*p3) ;							
  double dcos24dp4z = (Lepton2.z()*p2*p4 - Lepton4.z()*p2*p2*cos(Lepton2.Angle(Lepton4)))/(p2*p2*p4*p4) ;							
  double dcos34dp4z = (Lepton3.z()*p3*p4 - Lepton4.z()*p3*p3*cos(Lepton3.Angle(Lepton4)))/(p3*p3*p4*p4) ;	
  
  // derivatives of g = m^2 wrt px, py and pz
  double dgdp1x = 2.*dp1dp1x*p2*(1-cos(Lepton1.Angle(Lepton2))) + 2.*dp1dp1x*p3*(1-cos(Lepton1.Angle(Lepton3))) + 2.*dp1dp1x*p4*(1-cos(Lepton1.Angle(Lepton4)));
  dgdp1x = dgdp1x - 2.*p1*p2*dcos12dp1x - 2.*p1*p3*dcos13dp1x - 2.*p1*p4*dcos14dp1x ;						
  double dgdp1y = 2.*dp1dp1y*p2*(1-cos(Lepton1.Angle(Lepton2))) + 2.*dp1dp1y*p3*(1-cos(Lepton1.Angle(Lepton3))) + 2.*dp1dp1y*p4*(1-cos(Lepton1.Angle(Lepton4)));
  dgdp1y = dgdp1y - 2.*p1*p2*dcos12dp1y - 2.*p1*p3*dcos13dp1y - 2.*p1*p4*dcos14dp1y ;							
  double dgdp1z = 2.*dp1dp1z*p2*(1-cos(Lepton1.Angle(Lepton2))) + 2.*dp1dp1z*p3*(1-cos(Lepton1.Angle(Lepton3))) + 2.*dp1dp1z*p4*(1-cos(Lepton1.Angle(Lepton4)));
  dgdp1z = dgdp1z - 2.*p1*p2*dcos12dp1z - 2.*p1*p3*dcos13dp1z - 2.*p1*p4*dcos14dp1z ;							
  
  double dgdp2x = 2.*dp2dp2x*p1*(1-cos(Lepton1.Angle(Lepton2))) + 2.*dp2dp2x*p3*(1-cos(Lepton2.Angle(Lepton3))) + 2.*dp2dp2x*p4*(1-cos(Lepton2.Angle(Lepton4)));
  dgdp2x = dgdp2x -2.*p1*p2*dcos12dp2x - 2.*p2*p3*dcos23dp2x - 2.*p2*p4*dcos24dp2x ;		
  double dgdp2y = 2.*dp2dp2y*p1*(1-cos(Lepton1.Angle(Lepton2))) + 2.*dp2dp2y*p3*(1-cos(Lepton2.Angle(Lepton3))) + 2.*dp2dp2y*p4*(1-cos(Lepton2.Angle(Lepton4)));
  dgdp2y = dgdp2y - 2.*p1*p2*dcos12dp2y - 2.*p2*p3*dcos23dp2y - 2.*p2*p4*dcos24dp2y ;		
  double dgdp2z = 2.*dp2dp2z*p1*(1-cos(Lepton1.Angle(Lepton2))) + 2.*dp2dp2z*p3*(1-cos(Lepton2.Angle(Lepton3))) + 2.*dp2dp2z*p4*(1-cos(Lepton2.Angle(Lepton4)));
  dgdp2z = dgdp2z - 2.*p1*p2*dcos12dp2z - 2.*p2*p3*dcos23dp2z - 2.*p2*p4*dcos24dp2z ;		
  
  double dgdp3x = 2.*dp3dp3x*p1*(1-cos(Lepton1.Angle(Lepton3))) + 2.*dp3dp3x*p2*(1-cos(Lepton2.Angle(Lepton3))) + 2.*dp3dp3x*p4*(1-cos(Lepton3.Angle(Lepton4)));
  dgdp3x = dgdp3x - 2.*p1*p3*dcos13dp3x - 2.*p2*p3*dcos23dp3x - 2.*p3*p4*dcos34dp3x ;		
  double dgdp3y = 2.*dp3dp3y*p1*(1-cos(Lepton1.Angle(Lepton3))) + 2.*dp3dp3y*p2*(1-cos(Lepton2.Angle(Lepton3))) + 2.*dp3dp3y*p4*(1-cos(Lepton3.Angle(Lepton4)));
  dgdp3y = dgdp3y - 2.*p1*p3*dcos13dp3y - 2.*p2*p3*dcos23dp3y - 2.*p3*p4*dcos34dp3y ;		
  double dgdp3z = 2.*dp3dp3z*p1*(1-cos(Lepton1.Angle(Lepton3))) + 2.*dp3dp3z*p2*(1-cos(Lepton2.Angle(Lepton3))) + 2.*dp3dp3z*p4*(1-cos(Lepton3.Angle(Lepton4)));
  dgdp3z = dgdp3z - 2.*p1*p3*dcos13dp3z - 2.*p2*p3*dcos23dp3z - 2.*p3*p4*dcos34dp3z ;		
  
  double dgdp4x = 2.*dp4dp4x*p1*(1-cos(Lepton1.Angle(Lepton4))) + 2.*dp4dp4x*p2*(1-cos(Lepton2.Angle(Lepton4))) + 2.*dp4dp4x*p3*(1-cos(Lepton3.Angle(Lepton4)));
  dgdp4x = dgdp4x - 2.*p1*p4*dcos14dp4x - 2.*p2*p4*dcos24dp4x - 2.*p3*p4*dcos34dp4x ;		
  double dgdp4y = 2.*dp4dp4y*p1*(1-cos(Lepton1.Angle(Lepton4))) + 2.*dp4dp4y*p2*(1-cos(Lepton2.Angle(Lepton4))) + 2.*dp4dp4y*p3*(1-cos(Lepton3.Angle(Lepton4)));
  dgdp4y = dgdp4y - 2.*p1*p4*dcos14dp4y - 2.*p2*p4*dcos24dp4y - 2.*p3*p4*dcos34dp4y ;		
  double dgdp4z = 2.*dp4dp4z*p1*(1-cos(Lepton1.Angle(Lepton4))) + 2.*dp4dp4z*p2*(1-cos(Lepton2.Angle(Lepton4))) + 2.*dp4dp4z*p3*(1-cos(Lepton3.Angle(Lepton4)));
  dgdp4z = dgdp4z - 2.*p1*p4*dcos14dp4z - 2.*p2*p4*dcos24dp4z - 2.*p3*p4*dcos34dp4z ;		
  
  Jacobian(0,0)  = dgdp1x ;	
  Jacobian(0,1)  = dgdp1y ;
  Jacobian(0,2)  = dgdp1z ;
  Jacobian(0,3)  = dgdp2x ;
  Jacobian(0,4)  = dgdp2y ;
  Jacobian(0,5)  = dgdp2z ;
  Jacobian(0,6)  = dgdp3x ;
  Jacobian(0,7)  = dgdp3y ;
  Jacobian(0,8)  = dgdp3z ;
  Jacobian(0,9)  = dgdp4x ;
  Jacobian(0,10) = dgdp4y ;
  Jacobian(0,11) = dgdp4z ;	
  
  if (debug) {
    for (int i=0; i<12; i++) cout << "Jacobian(0,"<<i<<") = " << Jacobian(0,i) << endl;
  }	
  
  // now compute dg2=Jacobian*Error*Jacobian^T
  TMatrixD C(FourLeptonError,TMatrixD::kMultTranspose,Jacobian); 
  TMatrixD dg2=Jacobian*C;
  
  // g=m2, dg^2=4m^2dm^2																															  
  double dm2 = 0. ;
  if (invMassSquared!=0.) dm2=dg2(0,0)/(4.*invMassSquared);	
  
  //cout << "Invariant mass = " << sqrt(invMassSquared)	<< " variance = " << dm2 << " error = " << sqrt(dm2) << endl;																														  
  return  sqrt(dm2);
}

