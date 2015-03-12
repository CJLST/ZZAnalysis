/** \file
 *
 *  $Date: 2013/05/13 17:10:20 $
 *  $Revision: 1.6 $
 *  \author N. Amapane
 */

#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>
#include <Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h>
//#include <EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h>
#include <ZZAnalysis/AnalysisStep/interface/CustomElectronEffectiveArea.h>

#include <iostream>

using namespace std;
using namespace edm;
using namespace pat;
using namespace reco;


int correctionType = 2; //1 = rho; 2 = dbeta;

InputTag LeptonIsoHelper::getMuRhoTag(int sampleType, int setup) {
  InputTag rhoTag;
  if (setup==2011) {
    rhoTag = InputTag("fixedGridRhoFastjetAll","");
  } else if (setup==2012) { 
    rhoTag = InputTag("fixedGridRhoFastjetAll","");
  } else if (setup==2015) { 
    rhoTag = InputTag("fixedGridRhoFastjetAll","");
  } else {
    cout << "LeptonIsoHelper: Incorrect setup: " << setup << endl;
    abort();
  }
  return rhoTag;
}

InputTag LeptonIsoHelper::getEleRhoTag(int sampleType, int setup) {
  InputTag rhoTag;
  if (setup==2011) {
    rhoTag = InputTag("fixedGridRhoFastjetAll","");
  } else if (setup==2012) {
    rhoTag = InputTag("fixedGridRhoFastjetAll","");
  } else if (setup==2015) {
    rhoTag = InputTag("fixedGridRhoFastjetAll","");
  } else {
    cout << "LeptonIsoHelper: Incorrect setup: " << setup << endl;
    abort();
  }
  return rhoTag;
}


float LeptonIsoHelper::combRelIsoPF(int sampleType, int setup, double rho, const pat::Muon& l, float fsr) {
  float PFChargedHadIso   = l.chargedHadronIso();
  float PFNeutralHadIso   = l.neutralHadronIso();
  float PFPhotonIso       = l.photonIso();
  float PFPUChargedHadIso = l.puChargedHadronIso();
    
  MuonEffectiveArea::MuonEffectiveAreaTarget EAsetup;
  if (setup==2011) {
    EAsetup = MuonEffectiveArea::kMuEAData2011;
  } else if (setup==2012) { 
    EAsetup = MuonEffectiveArea::kMuEAData2012;
  } else if (setup==2015) { 
    EAsetup = MuonEffectiveArea::kMuEAData2012;
    //FIXME: no Run II muon effective areas in the Muon package yet
  } else {
    cout << "LeptonIsoHelper: Incorrect setup: " << setup << endl;
    abort();
  }

  if (correctionType==1) {
    float EA = MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuGammaAndNeutralHadronIso04, 
						       l.eta(), EAsetup);
    return  (PFChargedHadIso + max(0., PFNeutralHadIso + PFPhotonIso - fsr - rho * EA))/l.pt();

  } else if (correctionType==2) {
    return  (PFChargedHadIso + max(0., PFNeutralHadIso + PFPhotonIso - fsr - 0.5*PFPUChargedHadIso))/l.pt();
  }
  return 0;
}


float LeptonIsoHelper::combRelIsoPF(int sampleType, int setup, double rho, const pat::Electron& l, float fsr) {
  float PFChargedHadIso   = l.chargedHadronIso();
  float PFNeutralHadIso   = l.neutralHadronIso();
  float PFPhotonIso       = l.photonIso();

  ElectronEffectiveArea::ElectronEffectiveAreaTarget EAsetup;
  if (setup==2011) {
    EAsetup = ElectronEffectiveArea::kEleEAData2011;
  } else if (setup==2012) { 
    EAsetup = ElectronEffectiveArea::kEleEAData2012;
  } else if (setup==2015) { 
    EAsetup = ElectronEffectiveArea::kEleEAPhys14MC; //FIXME: replace with EAs from data when available
  } else {
    cout << "LeptonIsoHelper: Incorrect setup: " << setup << endl;
    abort();
  }

  float EA = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04,
							     l.superCluster()->eta(), EAsetup);
  return  (PFChargedHadIso + max(0., PFNeutralHadIso + PFPhotonIso - fsr - rho * EA))/l.pt();
}


float LeptonIsoHelper::combRelIsoPF(int sampleType, int setup, double rho, const Candidate* lep, float fsr) {
  // should check if lep->hasMasterClone()?  
  if (lep->isMuon()) {
    const pat::Muon* mu = dynamic_cast<const pat::Muon*>(lep->masterClone().get());
    return combRelIsoPF(sampleType, setup, rho, *mu, fsr);
  } else if (lep->isElectron()) {
    const pat::Electron* ele = dynamic_cast<const pat::Electron*>(lep->masterClone().get());
    return combRelIsoPF(sampleType, setup, rho, *ele, fsr);    
  } else {
    cout << "ERROR: LeptonIsoHelper: unknown type" << endl;
    abort();
  }
  return 0;
}

