#ifndef ZZAnalysis_Tools_Histograms_h
#define ZZAnalysis_Tools_Histograms_h

/** \class Histograms
 *
 * A set of histograms for LLLL candidates.
 *
 *  $Date: 2012/07/04 10:45:59 $
 *  $Revision: 1.9 $
 *  \author N. Amapane - Torino
 *  \author C. Botta - Torino
 */

#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include <iostream>
#include <algorithm>
#include <vector>

#ifdef HCMSSW
#include <ZZAnalysis/AnalysisStep/interface/TFileServiceWrapper.h>
#endif

using namespace std;


class HCand{
public:
  HCand(const TString& aName_)  {
    TString N = aName_;
    TH1F::SetDefaultSumw2(kTRUE);
    
    fs = new TFileServiceWrapper();
    
    //Plots for the candidate
    hZZMass = fs->makeTH1F(N+"_hZZMass",N+"_hZZMass", 2000, 0., 1000.);
    hZ1Mass = fs->makeTH1F(N+"_hZ1Mass",N+"_hZ1Mass", 400, 0., 200.);
    hZ2Mass = fs->makeTH1F(N+"_hZ2Mass",N+"_hZ2Mass", 400, 0., 200.);
//     hZaMass = fs->makeTH1F(N+"_hZaMass",N+"_hZaMass", 400, 0., 200.);
//     hZbMass = fs->makeTH1F(N+"_hZbMass",N+"_hZbMass", 400, 0., 200.);

    //Plots for the 4 individual leptons. Each set is sorted by the plotted variable 
    hPt1 = fs->makeTH1F(N+"_hPt1",N+"_hPt1", 150, 0., 150.);
    hPt2 = fs->makeTH1F(N+"_hPt2",N+"_hPt2", 150, 0., 150.);
    hPt3 = fs->makeTH1F(N+"_hPt3",N+"_hPt3", 150, 0., 150.);
    hPt4 = fs->makeTH1F(N+"_hPt4",N+"_hPt4", 150, 0., 150.);

    hCombRelIso1 = fs->makeTH1F(N+"_hCombRelIso1",N+"_hCombRelIso1", 100, 0., 5.);
    hCombRelIso2 = fs->makeTH1F(N+"_hCombRelIso2",N+"_hCombRelIso2", 100, 0., 5.);
    hCombRelIso3 = fs->makeTH1F(N+"_hCombRelIso3",N+"_hCombRelIso3", 100, 0., 5.);
    hCombRelIso4 = fs->makeTH1F(N+"_hCombRelIso4",N+"_hCombRelIso4", 100, 0., 5.);  
    
    hSIP1 = fs->makeTH1F(N+"_hSIP1",N+"_hSIP1", 500, 0., 100.);
    hSIP2 = fs->makeTH1F(N+"_hSIP2",N+"_hSIP2", 500, 0., 100.);
    hSIP3 = fs->makeTH1F(N+"_hSIP3",N+"_hSIP3", 500, 0., 100.);
    hSIP4 = fs->makeTH1F(N+"_hSIP4",N+"_hSIP4", 500, 0., 100.);
	  
    hCosThetaStar = fs->makeTH1F(N+"_hCosThetaStar",N+"_hCosThetaStar", 100, 0., 1.);
    hCosTheta1 = fs->makeTH1F(N+"_hCosTheta1",N+"_hCosTheta1", 100, 0., 1.);
    hCosTheta2 = fs->makeTH1F(N+"_hCosTheta2",N+"_hCosTheta2", 100, 0., 1.);
    hPhi = fs->makeTH1F(N+"_hPhi",N+"_hPhi", 100, -3.14159265, 3.14159265);
    hPhi1 = fs->makeTH1F(N+"_hPhi1",N+"_hPhi1", 100, -3.14159265, 3.14159265);
    hLD = fs->makeTH1F(N+"_hLD",N+"_hLD", 60, 0., 1.);
    hLD_himass = fs->makeTH1F(N+"_hLD_himass",N+"_hLD_himass", 60, 0., 1.);
    hLD_lowmass = fs->makeTH1F(N+"_hLD_lowmass",N+"_hLD_lowmass", 60, 0., 1.);
    hPseudoLD = fs->makeTH1F(N+"_hPseudoLD",N+"_hPseudoLD", 60, 0., 1.);
  }

  HCand(TString &aName, TFile* file){
    name = aName;
    TString N = aName;

    hZZMass = (TH1F *) file->Get(N+"_hZZMass");
    hZ1Mass = (TH1F *) file->Get(N+"_hZ1Mass");
    hZ2Mass = (TH1F *) file->Get(N+"_hZ2Mass");
//     hZaMass = (TH1F *) file->Get(N+"_hZaMass");
//     hZbMass = (TH1F *) file->Get(N+"_hZbMass");

    hPt1 = (TH1F *) file->Get(N+"_hPt1");
    hPt2 = (TH1F *) file->Get(N+"_hPt2");
    hPt3 = (TH1F *) file->Get(N+"_hPt3");
    hPt4 = (TH1F *) file->Get(N+"_hPt4");

    hCombRelIso1 = (TH1F *) file->Get(N+"_hCombRelIso1");
    hCombRelIso2 = (TH1F *) file->Get(N+"_hCombRelIso2");
    hCombRelIso3 = (TH1F *) file->Get(N+"_hCombRelIso3");
    hCombRelIso4 = (TH1F *) file->Get(N+"_hCombRelIso4");
    
    hSIP1 = (TH1F *) file->Get(N+"_hSIP1");
    hSIP2 = (TH1F *) file->Get(N+"_hSIP2");
    hSIP3 = (TH1F *) file->Get(N+"_hSIP3");
    hSIP4 = (TH1F *) file->Get(N+"_hSIP4");
	  
    hCosThetaStar = (TH1F *) file->Get(N+"_hCosThetaStar");
    hCosTheta1    = (TH1F *) file->Get(N+"_hCosTheta1");
    hCosTheta2    = (TH1F *) file->Get(N+"_hCosTheta2");
    hPhi          = (TH1F *) file->Get(N+"_hPhi");
    hPhi1         = (TH1F *) file->Get(N+"_hPhi1");
    hLD           = (TH1F *) file->Get(N+"_hLD");
    hLD_himass    = (TH1F *) file->Get(N+"_hLD_himass");
    hLD_lowmass   = (TH1F *) file->Get(N+"_hLD_lowmass");
    hPseudoLD     = (TH1F *) file->Get(N+"_hPseudoLD");

  }


  ~HCand(){
    delete fs;
  }

  // Assume that vectors are sorted!
  void Fill(float ZZMass, float Z1Mass, float Z2Mass, float ZaMass, float ZbMass,
	    std::vector<float>& pt, std::vector<float>& combRelIso, 
	    std::vector<float>& SIP, float weight, // FIXME
	    float cosThetaStar = 0, float cosTheta1 = 0, float cosTheta2 = 0, float phi = 0, float phi1 = 0, float LD=0, float pseudoLD=0) {

    
    hZZMass->Fill(ZZMass, weight);
    hZ1Mass->Fill(Z1Mass, weight);
    hZ2Mass->Fill(Z2Mass, weight);
//     hZaMass->Fill(ZaMass, weight);
//     hZbMass->Fill(ZbMass, weight);

    hPt1->Fill(pt[0], weight);
    hPt2->Fill(pt[1], weight);
    hPt3->Fill(pt[2], weight);
    hPt4->Fill(pt[3], weight);

    hCombRelIso1->Fill(combRelIso[0], weight);
    hCombRelIso2->Fill(combRelIso[1], weight);
    hCombRelIso3->Fill(combRelIso[2], weight);
    hCombRelIso4->Fill(combRelIso[3], weight);

    hSIP1->Fill(SIP[0], weight);
    hSIP2->Fill(SIP[1], weight);
    hSIP3->Fill(SIP[2], weight);
    hSIP4->Fill(SIP[3], weight);
	  
    hCosThetaStar->Fill(cosThetaStar, weight);
    hCosTheta1->Fill(cosTheta1, weight);
    hCosTheta2->Fill(cosTheta2, weight);
    hPhi->Fill(normalizePhi(phi), weight);
    hPhi1->Fill(normalizePhi(phi1), weight);
    hLD->Fill(LD,weight);
    if (ZZMass<=180) {
      hLD_lowmass->Fill(LD,weight);
    } else {
      hLD_himass->Fill(LD,weight);
    }
    
    if (LD>0.5) {
      hPseudoLD->Fill(pseudoLD,weight);
    }
  }


  void NormUnit() {
    float factor = 1./hZZMass->Integral(); //FIXME??
    Scale(factor);
  }


  void Scale(float factor) {
    cout << "Histograms::Scale needs to be fixed" << endl;
    abort();
//     if (hZZMass==0) cout << "error: no histogram in " << name << endl;

//     hZZMass->Scale(factor);
//     hZ1Mass->Scale(factor);
//     hZ2Mass->Scale(factor);
//     hZaMass->Scale(factor);
//     hZbMass->Scale(factor);
 
//     hPt1->Scale(factor);
//     hPt2->Scale(factor);
//     hPt3->Scale(factor);
//     hPt4->Scale(factor);

//     hCombRelIso1->Scale(factor);
//     hCombRelIso2->Scale(factor);
//     hCombRelIso3->Scale(factor);
//     hCombRelIso4->Scale(factor);

//     hSIP1->Scale(factor);
//     hSIP2->Scale(factor);
//     hSIP3->Scale(factor);
//     hSIP4->Scale(factor);
  }
  

  
 public:
  double normalizePhi(double theValue) { 
    if( theValue > TMath::TwoPi() || theValue < -TMath::TwoPi()) {
      theValue = fmod( theValue, TMath::TwoPi());
    }
    if (theValue <= -TMath::Pi()) theValue += TMath::TwoPi();
    if (theValue > TMath::Pi()) theValue -= TMath::TwoPi();
    return theValue;
  }


  TString name;
  TFileServiceWrapper* fs;

  TH1F * hZZMass;
  TH1F * hZ1Mass;
  TH1F * hZ2Mass;
//   TH1F * hZaMass;
//   TH1F * hZbMass;

  TH1F * hPt1;
  TH1F * hPt2;
  TH1F * hPt3;
  TH1F * hPt4;

  TH1F * hCombRelIso1;
  TH1F * hCombRelIso2;
  TH1F * hCombRelIso3;
  TH1F * hCombRelIso4;

  TH1F * hSIP1;
  TH1F * hSIP2;
  TH1F * hSIP3;
  TH1F * hSIP4;
	
  TH1F * hCosThetaStar;
  TH1F * hCosTheta1;
  TH1F * hCosTheta2;
  TH1F * hPhi;
  TH1F * hPhi1;
  TH1F * hLD;
  TH1F * hLD_himass;
  TH1F * hLD_lowmass;
  TH1F * hPseudoLD;
	

};


#endif
