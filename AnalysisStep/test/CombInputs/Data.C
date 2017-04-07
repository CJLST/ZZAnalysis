#include "TDirectory.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "THStack.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TCut.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"
#include "ZZAnalysis/AnalysisStep/interface/Category.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include "ZZAnalysis/AnalysisStep/src/cConstants.cc"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooWorkspace.h"
#include "RooHist.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "TPaveText.h"
#include "RooAddPdf.h"
#include "RooBreitWigner.h"
#include "RooFitResult.h"
#include "RooFFTConvPdf.h"
#include "RooAddition.h"
#include "RooMinuit.h"
#include "Math/MinimizerOptions.h"
#include <iomanip>
#include "RooAbsCollection.h"
#include "RooWorkspace.h"

using namespace RooFit;


#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>


int Wait() {
     cout << " Continue [<RET>|q]?  "; 
     char x;
     x = getchar();
     if ((x == 'q') || (x == 'Q')) return 1;
     return 0;
}

Double_t effSigma(TH1 *hist );
Double_t effSigma(RooAbsPdf *pdf, RooRealVar *obs, Int_t nbins);
float getFitEdge(float mass, float width, bool low);
float getFitEdgeHighMass(float mass, float width, bool low);

string schannel, scategory, ssample, sselAna;

void fitSignalShapeSimul(int massBin[40], int maxMassBin, int selAna, int ch, int cat, int sample, double rangeLow[40], 
                 double rangeHigh[40],double bwSigma[40], double fitValues[7], double fitErrors[7], double covQual[1]);


void all(int selAna =-10,  int channels=-1, int categ =-10, int sample = 0 ){

  if (selAna == 0) sselAna = "Morinod";
  if (selAna == 1) sselAna = "ICHEP";
  if (selAna == 2) sselAna = "Mor17";

  if (channels == 0) schannel = "4mu";
  if (channels == 1) schannel = "4e";
  if (channels == 2) schannel = "2e2mu";

  if (categ < 0 ) scategory = "ALL";
  if (categ == 0 && selAna == 0 ) scategory = "Untagged";
  if (categ == 1 && selAna == 0 ) scategory = "VBFtagged";

  if (categ == 0 && selAna == 1 ) scategory = "Untagged";
  if (categ == 1 && selAna == 1 ) scategory = "VBF1JetTagged";
  if (categ == 2 && selAna == 1 ) scategory = "VBF2JetTagged";
  if (categ == 3 && selAna == 1 ) scategory = "VHLeptTagged";
  if (categ == 4 && selAna == 1 ) scategory = "VHHadrTagged";
  if (categ == 5 && selAna == 1 ) scategory = "ttHTagged";

  if (categ == 0 && selAna == 2 ) scategory = "UntaggedMor17";
  if (categ == 1 && selAna == 2 ) scategory = "VBF1JetTaggedMor17";
  if (categ == 2 && selAna == 2 ) scategory = "VBF2JetTaggedMor17";
  if (categ == 3 && selAna == 2 ) scategory = "VHLeptTaggedMor17";
  if (categ == 4 && selAna == 2 ) scategory = "VHHadrTaggedMor17";
  if (categ == 5 && selAna == 2 ) scategory = "ttHTaggedMor17";
  if (categ == 6 && selAna == 2 ) scategory = "VHMETTaggedMor17";

  if (sample ==1) ssample = "ggH";
  if (sample ==2) ssample = "VBFH";

  double bwSigma[40];
  int mass[40]; int id[40]; double xLow[40]; double xHigh[40];
  int maxMassBin;
  maxMassBin = 1; 

  float masses[1] = {125};
  for(int i=0;i<1;++i) {
    mass[i] = masses[i]; 
    id[i]=masses[i]; 
    xLow[i] = 105.;  
    xHigh[i] = 140.;  
  }
  // -----------------------

  double massV[40],massE[40];
  for(int i=0; i<maxMassBin;++i){
    massV[i]=mass[i];
    massE[i]=0;
  }

  double a1Val[40],a1Err[40];
  double a2Val[40],a2Err[40];
  double n1Val[40],n1Err[40];
  double n2Val[40],n2Err[40];
  double meanCBVal[40],meanCBErr[40];
  double sigmaCBVal[40],sigmaCBErr[40];
  double covQualVal[40];
  double scale3[40];
  double fitValues[12];
  double fitErrors[12];
  double covQual[1];
    
  fitSignalShapeSimul(mass,maxMassBin,selAna,channels,categ,sample,/* 10.,doSfLepton,*/xLow,xHigh,bwSigma,fitValues,fitErrors,covQual); 
     
 
  }

    void fitSignalShapeSimul(int massBin[40],int maxMassBin, int selAna, int channels, int categ, int sample,double rangeLow[40], 
                             double rangeHigh[40], double bwSigma[40], double fitValues[7], double fitErrors[7], double covQual[1])
 {
 // ------ root settings ---------
  gROOT->Reset();  
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kFALSE);
  gStyle->SetPadGridY(kFALSE);
  gStyle->SetOptStat("iourme");
  gStyle->SetOptFit(11);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  // ------------------------------ 
  ROOT::Math::MinimizerOptions::SetDefaultTolerance( 1.E-7);
 
  float m4l;
  bool useQGTagging = false;
  bool useVHMETTagged = true;
  Short_t z1flav, z2flav; 
//  float weight;
  Double_t mass4l, kd;

  Short_t ExtraZ;
  Short_t nExtraLeptons;
  Short_t nCleanedJets;  

  float ZZPt, p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, PHJ_VAJHU, p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, PWH_hadronic_VAJHU, PZH_hadronic_VAJHU, PFMET;
 
  float p_GG_SIG_ghg2_1_ghz1_1_JHUGen, p_QQB_BKG_MCFM;
  Short_t nJets;
  Short_t nBTaggedJets;
  std::vector<float> * JETQGLikeliHood = 0;
  std::vector<float> * jetpt = 0;
  std::vector<float> * jeteta = 0;
  std::vector<float> * jetphi = 0;
  std::vector<float> * jetmass = 0;
  float jetQGLL[100];
  float jetPHI[100];
  float jet30pt[10];
  float jet30eta[10];
  float jet30phi[10];
  float jet30mass[10];
  float Fisher;

  double xMin,xMax,xInit;
  xMin = rangeLow[0];
  xMax = rangeHigh[maxMassBin-1];
  cout << "Fit range: [" << xMin << " , " << xMax << "]." << endl;

  RooRealVar x("mass","m_{4l}",125.,xMin,xMax,"GeV");
  x.setBins(70);
  RooRealVar w("myW","myW",1.0,-1000.,1000.);
  RooCategory massrc("massrc","massrc");
  char tempmass[4];
  for (int i=0; i<maxMassBin; i++) {
    sprintf(tempmass,"mh%d",massBin[i]);
    massrc.defineType(tempmass,massBin[i]);
  }
  RooArgSet ntupleVarSet(x,w,massrc);
  RooDataSet dataset("mass4l","mass4l",ntupleVarSet,WeightVar("myW"));

  stringstream FileName[40];
  for (int i=0; i<maxMassBin; i++) {
    if(sample==1) FileName[i] << "root://lxcms03//data3/Higgs/170222/AllData/ZZ4lAnalysis.root";
    else {
      cout << "Wrong sample ." << endl;
      return;
    }
    cout << "Using " << FileName[i].str() << endl;
    TFile* ggFile = TFile::Open(FileName[i].str().c_str()); 
    TTree* ggTree = (TTree*) ggFile->Get("ZZTree/candTree");
    int  nentries = ggTree->GetEntries();
    
    //--- ggTree part
    ggTree->SetBranchAddress("ZZMass",&m4l);
    ggTree->SetBranchAddress("Z1Flav",&z1flav);
    ggTree->SetBranchAddress("Z2Flav",&z2flav);
 //   ggTree->SetBranchAddress("overallEventWeight",&weight);
    ggTree->SetBranchAddress("nExtraLep",&nExtraLeptons);
    ggTree->SetBranchAddress("nCleanedJets",&nJets);
    ggTree->SetBranchAddress("nCleanedJetsPt30BTagged",&nBTaggedJets);
    ggTree->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal",&p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal);
    ggTree->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal",&p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal);

    ggTree->SetBranchAddress("DiJetFisher",&Fisher);
    ggTree->SetBranchAddress("nExtraZ",&ExtraZ);
    ggTree->SetBranchAddress("nCleanedJetsPt30",&nCleanedJets);
    ggTree->SetBranchAddress("JetQGLikelihood",&JETQGLikeliHood);
    ggTree->SetBranchAddress("p_JQCD_SIG_ghg2_1_JHUGen_JECNominal",&PHJ_VAJHU);
    ggTree->SetBranchAddress("p_JVBF_SIG_ghv1_1_JHUGen_JECNominal", &p_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
    ggTree->SetBranchAddress("pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal", &pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
    ggTree->SetBranchAddress("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal", &PWH_hadronic_VAJHU);
    ggTree->SetBranchAddress("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal",&PZH_hadronic_VAJHU);
    ggTree->SetBranchAddress("p_GG_SIG_ghg2_1_ghz1_1_JHUGen",&p_GG_SIG_ghg2_1_ghz1_1_JHUGen);
    ggTree->SetBranchAddress("p_QQB_BKG_MCFM",&p_QQB_BKG_MCFM);

    ggTree->SetBranchAddress("PFMET",&PFMET);
    ggTree->SetBranchAddress("JetPt",&jetpt);
    ggTree->SetBranchAddress("JetEta",&jeteta);
    ggTree->SetBranchAddress("JetPhi",&jetphi);
    ggTree->SetBranchAddress("JetMass",&jetmass);
    ggTree->SetBranchAddress("ZZPt",&ZZPt);

   
string sschannel, sscategory;

  if (channels == 0) sschannel = "4mu";
  if (channels == 1) sschannel = "4e";
  if (channels == 2) sschannel = "2e2mu";

  if (categ == 0 && selAna == 1 ) sscategory = "Untagged";
  if (categ == 1 && selAna == 1 ) sscategory = "VBF1JetTagged";
  if (categ == 2 && selAna == 1 ) sscategory = "VBF2JetTagged";
  if (categ == 3 && selAna == 1 ) sscategory = "VHLeptTagged";
  if (categ == 4 && selAna == 1 ) sscategory = "VHHadrTagged";
  if (categ == 5 && selAna == 1 ) sscategory = "ttHTagged";

  if (categ == 0 && selAna == 2 ) sscategory = "UntaggedMor17";
  if (categ == 1 && selAna == 2 ) sscategory = "VBF1JetTaggedMor17";
  if (categ == 2 && selAna == 2 ) sscategory = "VBF2JetTaggedMor17";
  if (categ == 3 && selAna == 2 ) sscategory = "VHLeptTaggedMor17";
  if (categ == 4 && selAna == 2 ) sscategory = "VHHadrTaggedMor17";
  if (categ == 5 && selAna == 2 ) sscategory = "ttHTaggedMor17";
  if (categ == 6 && selAna == 2 ) sscategory = "VHMETTaggedMor17";


    stringstream output;
    output << "data_obs_"<< sscategory<<"_"<< sschannel <<".root";
    cout<< output.str()<<endl;
    TFile *newFile = new TFile(output.str().c_str(),"RECREATE");
    newFile->cd();
    TTree* newTree = new TTree("data_obs","data_obs");
    newTree->Branch("mass4l",&mass4l,"mass4l/D");
    newTree->Branch("kd",&kd,"kd/D");
   

    //--- rooFit part
    xInit = (double) massBin[i];
    //---------  

    for(int k=0; k<nentries; k++){
      ggTree->GetEvent(k);
         
      int nj = 0;
      for (unsigned int nj = 0; nj < JETQGLikeliHood->size(); nj++) {
	jetQGLL[nj] = (*JETQGLikeliHood)[nj];
 	}
      for(unsigned int kjet =0 ; kjet < jetphi->size(); kjet++){
	jetPHI[kjet] = (*jetphi)[kjet];
	} 
      int njet30 = 0;
      for (unsigned int ijet = 0; ijet < jetpt->size(); ijet++) { 
	if ( (*jetpt)[ijet] > 30. ) {
	  jet30pt[njet30] = (*jetpt)[ijet];      
	  jet30eta[njet30] = (*jeteta)[ijet];
	  jet30phi[njet30] = (*jetphi)[ijet];
	  jet30mass[njet30] = (*jetmass)[ijet];
	  njet30++;
	}
      }  
      int Cat = -10 ;
      if (selAna == 0) Cat = categoryMor16(nJets, p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal );	    
      if (selAna == 1) Cat = categoryIchep16(nExtraLeptons, ExtraZ, nCleanedJets, nBTaggedJets, jetQGLL, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, PHJ_VAJHU, p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, PWH_hadronic_VAJHU, PZH_hadronic_VAJHU, jetPHI, m4l, useQGTagging);
      if (selAna == 2) Cat = categoryMor17(nExtraLeptons, ExtraZ, nCleanedJets, nBTaggedJets, jetQGLL, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, PHJ_VAJHU, p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, PWH_hadronic_VAJHU, PZH_hadronic_VAJHU, jetPHI, m4l, PFMET, useVHMETTagged, useQGTagging);
 

      if (categ >= 0 && categ != Cat ) continue;
      
      if(channels==0 && z1flav*z2flav != 28561) continue;
      if(channels==1 && z1flav*z2flav != 14641) continue;
      if(channels==2 && z1flav*z2flav != 20449) continue;
     
      mass4l=m4l;
      kd = p_GG_SIG_ghg2_1_ghz1_1_JHUGen / ( p_GG_SIG_ghg2_1_ghz1_1_JHUGen + p_QQB_BKG_MCFM *getDbkgkinConstant(z1flav*z2flav,m4l) ) ;
      newTree->Fill(); 
      
//      ntupleVarSet.setCatIndex("massrc",massBin[i]);
//      ntupleVarSet.setRealValue("mass",m4l);
//      ntupleVarSet.setRealValue("myW",weight);
//      if(x.getVal()>xMin && x.getVal()<xMax)
//	dataset.add(ntupleVarSet, weight);
      //--------
     
    }
  newTree->Write("data_obs");
  newFile->Close();

  }
}

