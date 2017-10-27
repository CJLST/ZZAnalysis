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
  maxMassBin = 5; 

  float masses[5] = {120,124,125,126,130};
  for(int i=0;i<5;++i) {
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
 
      cout << "meanCB_p0 value "  << fitValues[0] << " , " << "mean_p1 value:"  << fitValues[6] << endl;
      cout << "sigmaCB_p0 value " << fitValues[1] << " , " << "sigma_p1 value:" << fitValues[7] << endl;
      cout << "a1_p0 value "      << fitValues[2] << " , " << "a1_p1 value:"    << fitValues[8]  << endl;
      cout << "n1_p0 value "      << fitValues[3] << " , " << "n1_p1 value:"    << fitValues[9]  << endl;
      cout << "a2_p0 value "      << fitValues[4] << " , " << "a2_p1 value:"    << fitValues[10]  << endl;
      cout << "n2_p0 value "      << fitValues[5] << " , " << "n2_p1 value:"    << fitValues[11]  << endl;


      string filename = "signal_shape_parametrization_13TeV_" + ssample + "_" + schannel + "_" + sselAna + "_" + scategory + "." + "yaml" ;
      ofstream outFile;
      outFile.open(filename);
      if(channels == 2)outFile<<"shape : " <<"\"RooDCBall::"<<ssample<<"_mass(mean,sigma,alpha,n,alpha2,n2)\""<< endl;
      outFile << schannel <<"    :" << endl;
      outFile <<"    mean   : " <<"'"<<fitValues[0]<<"+"<<"("<<fitValues[6] <<"+1)*(@0-125)"<<"'"<<endl;
      outFile <<"    sigma  : " <<"'"<<fitValues[1]<<"+"<<"("<<fitValues[7] <<")*(@0-125)"<<"'"<<endl;
      outFile <<"    alpha  : " <<"'"<<fitValues[2]<<"+"<<"("<<fitValues[8] <<")*(@0-125)"<<"'"<<endl;
      outFile <<"    n      : " <<"'"<<fitValues[3]<<"+"<<"("<<fitValues[9] <<")*(@0-125)"<<"'"<<endl;
      outFile <<"    alpha2 : " <<"'"<<fitValues[4]<<"+"<<"("<<fitValues[10]<<")*(@0-125)"<<"'"<<endl;
      outFile <<"    n2     : " <<"'"<<fitValues[5]<<"+"<<"("<<fitValues[11]<<")*(@0-125)"<<"'"<<endl;
      outFile << endl;
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
  
  Short_t z1flav, z2flav; 
  float weight;
  bool useQGTagging = false; 
  bool useVHMETTagged = true; 
  Short_t ExtraZ;
  Short_t nExtraLeptons;
  Short_t nCleanedJets;   

  float ZZPt, p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, PHJ_VAJHU, p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, PWH_hadronic_VAJHU, PZH_hadronic_VAJHU,p_HadWH_mavjj_JECNominal, p_HadZH_mavjj_JECNominal,PFMET;
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
    if(sample==1)      FileName[i] << "root://lxcms03//data3/Higgs/170222/ggH" << massBin[i] << "/ZZ4lAnalysis.root";
    else if(sample==2) FileName[i] << "root://lxcms03//data3/Higgs/170222/VBFH" << massBin[i] << "/ZZ4lAnalysis.root";
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
    ggTree->SetBranchAddress("overallEventWeight",&weight);
    ggTree->SetBranchAddress("nExtraLep",&nExtraLeptons);
    ggTree->SetBranchAddress("nCleanedJets",&nJets);
    ggTree->SetBranchAddress("nCleanedJetsPt30BTagged",&nBTaggedJets);
    ggTree->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSF",&nBTaggedJets);

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
    ggTree->SetBranchAddress("p_HadWH_mavjj_JECNominal", &p_HadWH_mavjj_JECNominal);
    ggTree->SetBranchAddress("p_HadZH_mavjj_JECNominal",&p_HadZH_mavjj_JECNominal);

    ggTree->SetBranchAddress("PFMET",&PFMET);
    ggTree->SetBranchAddress("JetPt",&jetpt);
    ggTree->SetBranchAddress("JetEta",&jeteta);
    ggTree->SetBranchAddress("JetPhi",&jetphi);
    ggTree->SetBranchAddress("JetMass",&jetmass);
    ggTree->SetBranchAddress("ZZPt",&ZZPt);

    //--- rooFit part
    xInit = (double) massBin[i];
    //---------  

    for(int k=0; k<nentries; k++){
      ggTree->GetEvent(k);
         
      int nj = 0;
      for (unsigned int nj = 0; nj < JETQGLikeliHood->size(); nj++) {
	jetQGLL[nj] = (*JETQGLikeliHood)[nj];
 	}
      int kj = 0;
      for (unsigned int kj = 0; kj < jetphi->size(); kj++) {
        jetPHI[kj] = (*jetphi)[nj];
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
	   if (selAna == 1) Cat = categoryIchep16(nExtraLeptons, ExtraZ, nCleanedJets, nBTaggedJets, jetQGLL, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, PHJ_VAJHU, p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, PWH_hadronic_VAJHU, PZH_hadronic_VAJHU, p_HadWH_mavjj_JECNominal, p_HadZH_mavjj_JECNominal,jetPHI, m4l, useQGTagging);
      if (selAna == 2) Cat = categoryMor17(nExtraLeptons, ExtraZ, nCleanedJets, nBTaggedJets, jetQGLL, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, PHJ_VAJHU, p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, PWH_hadronic_VAJHU, PZH_hadronic_VAJHU, p_HadWH_mavjj_JECNominal, p_HadZH_mavjj_JECNominal,jetPHI, m4l, PFMET, useVHMETTagged, useQGTagging);

      if (categ >= 0 && categ != Cat ) continue;
      
      if(channels==0 && z1flav*z2flav != 28561) continue;
      if(channels==1 && z1flav*z2flav != 14641) continue;
      if(channels==2 && z1flav*z2flav != 20449) continue;
      if (weight <= 0 ) cout << "Warning! Negative weight events" << endl;
      
      ntupleVarSet.setCatIndex("massrc",massBin[i]);
      ntupleVarSet.setRealValue("mass",m4l);
      ntupleVarSet.setRealValue("myW",weight);
      if(x.getVal()>xMin && x.getVal()<xMax)
	dataset.add(ntupleVarSet, weight);
      //--------
     
    }
  }

  cout << "dataset n entries: " << dataset.sumEntries() << endl;
  RooDataSet *dataset125 = (RooDataSet*)dataset.reduce("massrc == massrc::mh125");
  cout << "dataset 125 n entries: " << dataset125->sumEntries() << endl;

  TCanvas *c1 = new TCanvas("c1","c1",725,725);
  TCanvas *c2 = new TCanvas("c2","c2",725,725);


  //--- double CrystalBall

  RooRealVar mH("mH","MH",125.,90.,900.);

  RooRealVar mean_p0("mean_p0","mean_p0",125.,120., 130.) ;
  RooRealVar mean_p1("mean_p1","mean_p1",0, -1.5, 1.5);

  RooRealVar sigma_p0("sigma_p0","sigma_p0",1.6, 0, 30);
  RooRealVar sigma_p1("sigma_p1","sigma_p1",0,-0.5, 0.5);

  RooRealVar a1_p0("a1_p0","a1_p0", 1.46, 0.5, 5.);
  RooRealVar a1_p1("a1_p1","a1_p1", 0, -0.5, 0.5);

  RooRealVar n1_p0("n1_p0","n_p0", 1.92, 0, 10);
  RooRealVar n1_p1("n1_p1","n_p1", 0, -0.5, 0.5);

  RooRealVar a2_p0("a2_p0","a2_p0", 1.46, 1, 10);
  RooRealVar a2_p1("a2_p1","a2_p1", 0, -0.5, 0.5);

  RooRealVar n2_p0("n2_p0","n2_p0", 20., 1., 50.);
  RooRealVar n2_p1("n2_p1","n2_p1", 0, -1, 1);
  
  RooArgSet* params_Inter = new RooArgSet(RooArgList( mean_p0, sigma_p0, a1_p0, n1_p0, a2_p0, n2_p0));
  RooArgSet* params_slope = new RooArgSet(RooArgList( mean_p1, sigma_p1, a1_p1, n1_p1, a2_p1, n2_p1));
  
  if (sample != 1 || categ > 0 ) {
    if(channels==0 ) params_Inter->readFromFile("Ch0_Cat0_paraI.txt");
    if(channels==1 ) params_Inter->readFromFile("Ch1_Cat0_paraI.txt");
    if(channels==2 ) params_Inter->readFromFile("Ch2_Cat0_paraI.txt");
  }

  // Prefit at 125
  RooDoubleCB DCBall("DCBall","Double Crystal ball",x,mean_p0,sigma_p0,a1_p0,n1_p0,a2_p0,n2_p0);

  RooFitResult *fitres = (RooFitResult*)DCBall.fitTo(*dataset125,SumW2Error(1),Range(xMin,xMax),Strategy(2),NumCPU(8),Save(true));
  
  cout << "Prefit done" << endl;

  stringstream frameTitle;
  if(channels==0){frameTitle << "4#mu, m_{H} = "; }
  if(channels==1){frameTitle << "4e, m_{H} = ";}
  if(channels==2){frameTitle << "2e2#mu, m_{H} = ";}
  frameTitle << "125 GeV";

  RooPlot* xframe = x.frame() ;
  xframe->SetTitle("");
  xframe->SetName("m4lplot");
  dataset125->plotOn(xframe,DataError(RooAbsData::SumW2), MarkerStyle(kOpenCircle), MarkerSize(1.1) );
  int col ;
  if(channels==0) col=kOrange+7;
  if(channels==1) col=kAzure+2; 
  if(channels==2) col=kGreen+3; 
  DCBall.plotOn(xframe,LineColor(col));

  RooHist* hpull = xframe->pullHist();
  RooPlot* frame3 = x.frame(Title("Pull Distribution")) ;
  frame3->addPlotable(hpull,"P");

  double * Rentries = hpull->GetY();
  TH1D *hh = new TH1D("Pull Distribution","Pull Distribution",10,-5,5);
  for (int i =0; i < 70; i++ ){
  hh->Fill(Rentries[i]);}

  TLegend *legend = new TLegend(0.20,0.45,0.45,0.60,NULL,"brNDC");
  legend->SetBorderSize(     0);
  legend->SetFillColor (     0);
  legend->SetTextAlign (    12);
  legend->SetTextFont  (    42);
  legend->SetTextSize  (0.03);

  TH1F *dummyPoints = new TH1F("dummyP","dummyP",1,0,1);
  TH1F *dummyLine = new TH1F("dummyL","dummyL",1,0,1);
  TH1F *dummyLine1 = new TH1F("dummyL","dummyL",1,0,1);

  dummyPoints->SetMarkerStyle(kOpenCircle);
  dummyPoints->SetMarkerSize(1.1);
  dummyLine->SetLineColor(col);

  legend->AddEntry(dummyPoints, "Simulation", "pe");
  legend->AddEntry(dummyLine, "Parametric Model", "l");

  TPaveText *text = new TPaveText(0.15,0.90,0.77,0.98,"brNDC");
  text->AddText("CMS Simulation");
  text->SetBorderSize(0);
  text->SetFillStyle(0);
  text->SetTextAlign(12);
  text->SetTextFont(42);
  text->SetTextSize(0.03);

  TPaveText *titlet = new TPaveText(0.15,0.80,0.60,0.85,"brNDC");
  titlet->AddText(frameTitle.str().c_str());
  titlet->SetBorderSize(0);
  titlet->SetFillStyle(0);
  titlet->SetTextAlign(12);
  titlet->SetTextFont(132);
  titlet->SetTextSize(0.045);

  TPaveText *sigmat = new TPaveText(0.15,0.65,0.77,0.78,"brNDC");
  stringstream sigmaval0, sigmaval1, sigmaval2;
  sigmaval0 << fixed;
  sigmaval0 << setprecision(1);
  sigmaval0 << "m_{dCB} = " << mean_p0.getVal() << " GeV";
  sigmaval1 << fixed;
  sigmaval1 << setprecision(1);
  sigmaval1 << "#sigma_{dCB} = " << sigma_p0.getVal() << " GeV";
  
  sigmat->AddText(sigmaval0.str().c_str());
  sigmat->AddText(sigmaval1.str().c_str());
  sigmat->SetBorderSize(0);
  sigmat->SetFillStyle(0);
  sigmat->SetTextAlign(12);
  sigmat->SetTextFont(132);
  sigmat->SetTextSize(0.04);
  
  xframe->GetYaxis()->SetTitleOffset(1.5);

  c1->cd(); 
  stringstream nameFile, nameFileC, nameFilePng;
  nameFile    << "PrefitM125" << "_" << sselAna << "_" << ssample << "_" << schannel << "_"<< scategory << ".pdf";
  nameFileC   << "PrefitM125" << "_" << sselAna << "_" << ssample << "_" << schannel << "_"<< scategory << ".C";
  nameFilePng << "PrefitM125" << "_" << sselAna << "_" << ssample << "_" << schannel << "_"<< scategory << ".png";
  
  xframe->Draw();
  gPad->Update(); legend->Draw(); text->Draw(); titlet->Draw();

  c1->Print (nameFile.str().c_str());
  c1->SaveAs(nameFileC.str().c_str());
  c1->SaveAs(nameFilePng.str().c_str());

  c2->Divide(1,2);

  c2->cd(1) ;
  hh->Draw() ;

  c2->cd(2) ;
  frame3->Draw() ;
  frame3->SetMinimum(-3);
  frame3->SetMaximum(3);

  TLine *line1 = new TLine(105,0,150,0);
  line1->SetLineColor(kRed);
  line1->Draw();

  stringstream nameFilePull, nameFilePullC, nameFilePullPng;
  nameFilePull    << "Pre_PullM"    << "_" << sselAna << "_" << ssample << "_" << schannel << "_"<< scategory << ".pdf";
  nameFilePullC    << "Pre_PullM"   << "_" << sselAna << "_" << ssample << "_" << schannel << "_"<< scategory << ".C";
  nameFilePullPng   << "Pre_PullM"  << "_" << sselAna << "_" << ssample << "_" << schannel << "_"<< scategory << ".png";

  c2->Print(nameFilePull.str().c_str());
  c2->SaveAs(nameFilePullC.str().c_str());
  c2->SaveAs(nameFilePullPng.str().c_str());

  mean_p0.setConstant(kTRUE);
  sigma_p0.setConstant(kTRUE);
  a1_p0.setConstant(kTRUE);
  n1_p0.setConstant(kTRUE);
  a2_p0.setConstant(kTRUE);
  n2_p0.setConstant(kTRUE);
  RooArgSet * params = DCBall.getParameters(x);

  if (sample == 1 && categ == 0 ){ 
    if(channels==0 ) params->writeToFile("Ch0_Cat0_paraI.txt") ;
    if(channels==1 ) params->writeToFile("Ch1_Cat0_paraI.txt") ;
    if(channels==2 ) params->writeToFile("Ch2_Cat0_paraI.txt") ;
  }

  //a1_p1.setConstant(kTRUE);
  n1_p1.setConstant(kTRUE);
  //a2_p1.setConstant(kTRUE);
  n2_p1.setConstant(kTRUE);
  
  // Simultaneous fit
  RooSimultaneous rs("rs","rs",massrc);  

  RooFormulaVar* mean[40];
  RooFormulaVar* sigma[40];
  RooFormulaVar* a1[40];
  RooFormulaVar* a2[40];
  RooFormulaVar* n1[40];
  RooFormulaVar* n2[40];
  RooDoubleCB* DCBall2[40];

  for (int i=0; i<maxMassBin; i++) {
    char formulamass[200];
    int massdiff = massBin[i] - 125;
    sprintf(formulamass,"@0+(@1*%d)+%d",massdiff,massdiff);
    mean[i] = new RooFormulaVar("mean_CB","mean_CB",formulamass,RooArgList(mean_p0,mean_p1));
    sprintf(formulamass,"@0+(@1*%d)",massdiff);
    sigma[i] = new RooFormulaVar("sigma_CB","sigma_CB",formulamass,RooArgList(sigma_p0,sigma_p1));
    a1[i] = new RooFormulaVar("a1_CB","a1_CB",formulamass,RooArgList(a1_p0,a1_p1));
    n1[i] = new RooFormulaVar("n1_CB","n1_CB",formulamass,RooArgList(n1_p0,n1_p1));
    a2[i] = new RooFormulaVar("a2_CB","a2_CB",formulamass,RooArgList(a2_p0,a2_p1));
    n2[i] = new RooFormulaVar("n2_CB","n2_CB",formulamass,RooArgList(n2_p0,n2_p1));
    DCBall2[i] = new RooDoubleCB("DCBall2","Double Crystal ball 2",x,*mean[i],*sigma[i],*a1[i],*n1[i],*a2[i],*n2[i]);

    char tempmass[4];
    sprintf(tempmass,"mh%d",massBin[i]);
    rs.addPdf(*DCBall2[i], tempmass);

  }

 RooArgSet * params2 = rs.getParameters(RooArgList(x,massrc));
 if (sample != 1 ||  categ > 0 ) {
    if(channels==0 )params2->readFromFile("Ch0_Cat0_paraf.txt") ;
    if(channels==1 )params2->readFromFile("Ch1_Cat0_paraf.txt") ;
    if(channels==2 )params2->readFromFile("Ch2_Cat0_paraf.txt") ;	
  }

  RooFitResult *fitressim = (RooFitResult*)rs.fitTo(dataset,SumW2Error(1),Range(xMin,xMax),Strategy(2),NumCPU(8),Save(true));
 
  mean_p1.setConstant(kTRUE);
  sigma_p1.setConstant(kTRUE);
  a1_p1.setConstant(kTRUE);
  n1_p1.setConstant(kTRUE);
  a2_p1.setConstant(kTRUE);
  n2_p1.setConstant(kTRUE);

  if (sample == 1 && categ == 0 ) {
  
    if(channels==0 )params2->writeToFile("Ch0_Cat0_paraf.txt") ;
    if(channels==1 )params2->writeToFile("Ch1_Cat0_paraf.txt") ;
    if(channels==2 )params2->writeToFile("Ch2_Cat0_paraf.txt") ;	
  }
 
  cout << "Full fit done" << endl;

  double ChiSq[5]; 
  for (int i=0; i<maxMassBin; i++) {

    stringstream frameTitle;
    if(channels==0){frameTitle << "4#mu, m_{H} = "; }
    if(channels==1){frameTitle << "4e, m_{H} = ";}
    if(channels==2){frameTitle << "2e2#mu, m_{H} = ";}
    frameTitle << massBin[i] << " GeV";

    RooPlot* xframe2 = x.frame() ;
    xframe2->SetTitle("");
    xframe2->SetName("m4lplot");
    char tempmass[84];
    sprintf(tempmass,"massrc == massrc::mh%d",massBin[i]);
    char tempmass2[4];
    sprintf(tempmass2,"mh%d",massBin[i]);
    
    dataset.plotOn(xframe2,DataError(RooAbsData::SumW2), MarkerStyle(kOpenCircle), MarkerSize(1.1), Cut(tempmass) );
    rs.plotOn(xframe2,LineColor(col),Slice(massrc,tempmass2),ProjWData(massrc,dataset));

    RooHist* hpull = xframe2->pullHist();
    RooPlot* frame4 = x.frame(Title("Pull Distribution")) ;
    frame4->addPlotable(hpull,"P");
    cout << "Reduced chi2 = " << xframe2->chiSquare() << endl;
    ChiSq[i] = xframe2->chiSquare();

  double * Rentries = hpull->GetY();
  TH1D *hh1 = new TH1D("Pull Distribution","Pull Distribution",10,-5,5);
  for (int i =0; i < 70; i++ ){
  hh1->Fill(Rentries[i]);}
  
    // cosmetics
    TPaveText *titlet2 = new TPaveText(0.15,0.80,0.60,0.85,"brNDC");
    titlet2->AddText(frameTitle.str().c_str());
    titlet2->SetBorderSize(0);
    titlet2->SetFillStyle(0);
    titlet2->SetTextAlign(12);
    titlet2->SetTextFont(132);
    titlet2->SetTextSize(0.045);   
  
    xframe2->GetYaxis()->SetTitleOffset(1.5);
    
    c1->cd();
    stringstream nameFile, nameFileC, nameFilePng;
    nameFile    << "fitM" << massBin[i] << "_" << sselAna << "_" << ssample << "_" << schannel << "_"<< scategory << ".pdf";
    nameFileC   << "fitM" << massBin[i] << "_" << sselAna << "_" << ssample << "_" << schannel << "_"<< scategory << ".C";
    nameFilePng << "fitM" << massBin[i] << "_" << sselAna << "_" << ssample << "_" << schannel << "_"<< scategory << ".png";
  
    xframe2->Draw(); 
    gPad->Update(); legend->Draw(); text->Draw(); titlet2->Draw();
    c1->Print (nameFile.str().c_str());
    c1->SaveAs(nameFileC.str().c_str());
    c1->SaveAs(nameFilePng.str().c_str());

    c2->cd(1) ;
    hh1->Draw();

    c2->cd(2) ;
    frame4->Draw() ;
    frame4->SetMinimum(-3);
    frame4->SetMaximum(3);
    line1->Draw();

    stringstream nameFilePull, nameFilePullC, nameFilePullPng;
    nameFilePull     << "PullM" << massBin[i] << "_" << sselAna << "_" << ssample << "_" << schannel << "_"<< scategory << ".pdf";
    nameFilePullC    << "PullM" << massBin[i] << "_" << sselAna << "_" << ssample << "_" << schannel << "_"<< scategory << ".C";
    nameFilePullPng  << "PullM" << massBin[i] << "_" << sselAna << "_" << ssample << "_" << schannel << "_"<< scategory << ".png";
    c2->Print (nameFilePull.str().c_str());
    c2->SaveAs(nameFilePullC.str().c_str());
    c2->SaveAs(nameFilePullPng.str().c_str());
    }

    for (int j=0; j<5 ; j++ )
    cout<<"Chi2 for mass " << massBin[j] << " GeV: " << ChiSq[j] << endl; 

   Int_t ni = 5;
   Double_t xi[ni], yi[ni];
   for (Int_t i=0;i<ni;i++) {
      xi[i] = massBin[i];
      yi[i] = ChiSq[i];
   }
    stringstream nameGraphChi;
    nameGraphChi    << "Chi2"  << "_" << ssample << "_" << schannel << "_"<< scategory;
    TGraph *gr  = new TGraph(ni,xi,yi);
    gr->SetTitle(nameGraphChi.str().c_str());
    gr->GetXaxis()->SetTitle("Mass GeV");
    gr->GetYaxis()->SetTitle("Chi2/ndof");
    gr->SetMinimum(0);
    gr->SetMaximum(5);

    TCanvas *c3 = new TCanvas("c3","Chi2 ",200,10,600,400);
    gr->Draw();
    stringstream nameFileChi, nameFilePngChi;
    nameFileChi    << "Chi2" << "_" << sselAna << "_" << ssample << "_" << schannel << "_"<< scategory << ".pdf"; 
    nameFilePngChi    << "Chi2" << "_" << sselAna << "_" << ssample << "_" << schannel << "_"<< scategory << ".png";               
    c3->SaveAs(nameFileChi.str().c_str());
    c3->SaveAs(nameFilePngChi.str().c_str());


    if(fitValues!=0){
   
    fitValues[0] = mean_p0.getVal();
    fitValues[1] = sigma_p0.getVal();
    fitValues[2] = a1_p0.getVal();
    fitValues[3] = n1_p0.getVal();
    fitValues[4] = a2_p0.getVal();
    fitValues[5] = n2_p0.getVal();
        
    fitValues[6]  = mean_p1.getVal();
    fitValues[7]  = sigma_p1.getVal();
    fitValues[8]  = a1_p1.getVal();
    fitValues[9]  = n1_p1.getVal();
    fitValues[10] = a2_p1.getVal();
    fitValues[11] = n2_p1.getVal();

    }  
  
   }

