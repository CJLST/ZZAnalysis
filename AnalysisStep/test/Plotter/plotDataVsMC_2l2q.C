/* 
 * usage: 
 * -specify parameters at the end of this file
 * -run with:
 *   root -l plotDataVsMC.C++
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Math/DistFunc.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TROOT.h"
#include "TSpline.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"

#include "tdrstyle.C"
// #include "CMS_lumi.C"
#include "plotUtils.C"
#include "ZZAnalysis/AnalysisStep/src/kFactors.C"

#include <ZZAnalysis/AnalysisStep/src/Category.cc>
#include <ZZAnalysis/AnalysisStep/src/bitops.cc>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>

using namespace std;




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


const int nVariables = 18;
string varName[nVariables] = {
  "ZZMass",
  "ZZPt",
  "ZZMassRefit",
  "Z1Mass",
  "Z2Mass",
  "Z2Flav",
  "Jet1Pt",
  "Jet2Pt",
  "JetQG",
  "JetBtagger",
  "Lep1Pt",
  "Lep2Pt",
  "AbsCosTheta1",
  "CosTheta2",
  "CosThetaStar",
  "PhiStar",
  "Phi", 
  "NExtraJets"
};
string varXLabel[nVariables] = {
  "m_{2#font[12]{l}2q} (GeV)",
  "p_{T2#font[12]{l}2q} (GeV)",
  "m_{2#font[12]{l}2q} (GeV)",
  "m_{jj} (GeV)",
  "m_{#font[12]{l}#font[12]{l}} (GeV)",
  "flavor",
  "p_{Tj1} (GeV)",
  "p_{Tj2} (GeV)",
  "Quark gluon likelihood",
  "B tagger",
  "p_{T#font[12]{l}1} (GeV)",
  "p_{T#font[12]{l}2} (GeV)",
  "|cos(#theta_{1})|",
  "cos(#theta_{2})",
  "cos(#theta^{*})",
  "#Phi^{*}",
  "#Phi",
  "N_{extra-jets}"
};
string varYLabel[nVariables] = {
  "Events / 25 GeV",
  "Events / 10 GeV",
  "Events / 25 GeV",
  "Events / 2.5 GeV",
  "Events / 2.5 GeV",
  "Events",
  "Events / 10 GeV",
  "Events / 10 GeV",
  "Events / 20 GeV",
  "Events / 0.028",
  "Events / 0.056",
  "Events",
  "Events",
  "Events",
  "Events",
  "Events",
  "Events",
  "Events"
};
Int_t  varNbin[nVariables] = { 50, 50, 50,  44,  44,  400,  50,  50,  50,  50,  50,  50,  50, 50, 50, 50, 50, 4};
Float_t varMin[nVariables] = {  250,  0,  250,  40,  40,  -200,  0, 0, -0.2, -0.2, 0,  0, -0.2, -1.2, -1.2, -3.15, -3.15, -0.5 };
Float_t varMax[nVariables] = { 1500, 500, 1500, 150, 150, 0, 500, 500, 1.2, 1.2, 500, 500, 1.2, 1.2, 1.2 , 3.15, 3.15, 3.5 };
Bool_t varLogx[nVariables] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
Bool_t varLogy[nVariables] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


enum Process {Data=0, BulkG=1, Spin0=2, DYjets=3};
const int nProcesses = 4;
string sProcess[nProcesses] = {"Data", "BulkG", "Spin0", "DY"};
string processLabel[nProcesses] = {"Data", "G^{*}(800)#rightarrowZZ (x50)", "H_{NWA}#rightarrowZZ (x50)", "Z + 2/3/4 jets"};
Float_t scaleF[nProcesses] = {1.,50.,50.,1.};

const int nFS = 3;
string sFS[nFS] = {"ee","all","mm"};



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void plotDataVsMC_2l2q(string dirout = "test13TeV")
{
  
  float lumin = 10.0;   // ICHEP?
  setTDRStyle();
  // gStyle->SetOptStat(1111111);
  const int nDatasets = 5;
  
  TFile* inputFile[nDatasets];
  TChain* inputTree[nDatasets];
  TH1F* hCounters[nDatasets]; 
  Long64_t NGenEvt[nDatasets];
  // Float_t partialEventWeight[nDatasets];

  Float_t overallEventWeight;
  Float_t ZZMass;
  Float_t ZZMassRefit;
  Float_t ZZPt;
  Float_t Z1Mass;
  Float_t Z2Mass;
  Short_t Z2Flav;
  vector<Float_t> *LepPt = 0;
  Short_t nJets;
  Short_t nJetsBTagged;
  vector<Float_t> *JetPt = 0;
  vector<Short_t> *JetIsInZZCand = 0;
  vector<Float_t> *JetBTagger = 0;
  vector<Float_t> *JetQGLikelihood = 0;
  Float_t helcosthetaZ1;  
  Float_t helcosthetaZ2;  
  Float_t costhetastar;	  
  Float_t helphi;	  
  Float_t phistarZ1;  
  Float_t xsec;

  TH1F* h1[nVariables][nProcesses][nFS];
  for(int rs=0; rs<nFS; rs++){   //ee, mumu, or all
    for(int pr=0; pr<nProcesses; pr++){
      for(int v=0; v<nVariables; v++){
	h1[v][pr][rs] = new TH1F(Form("h1_%s_%s_%s",varName[v].c_str(),sFS[rs].c_str(),sProcess[pr].c_str()),
				 Form("h1_%s_%s_%s",varName[v].c_str(),sFS[rs].c_str(),sProcess[pr].c_str()),		
				 varNbin[v],varMin[v],varMax[v]);
      }
    }
  }

  //---------- Will loop over all datasets

  for (int d=0; d<nDatasets; d++) {
    inputTree[d] = new TChain("ZZTree/candTree");
  }
  inputTree[0]->Add("../prod/test2l2q/BulkGrav800_Chunk0/ZZ2l2qAnalysis.root");
  inputFile[0] = new TFile("../prod/test2l2q/BulkGrav800_Chunk0/ZZ2l2qAnalysis.root");
  hCounters[0] = (TH1F*)inputFile[0]->Get("ZZTree/Counters");
  NGenEvt[0] = hCounters[0]->GetBinContent(1);

  inputTree[1]->Add("../prod/test2l2q/Higgs750_Chunk0/ZZ2l2qAnalysis.root");
  inputFile[1] = new TFile("../prod/test2l2q/Higgs750_Chunk0/ZZ2l2qAnalysis.root");
  hCounters[1] = (TH1F*)inputFile[1]->Get("ZZTree/Counters");
  NGenEvt[1] = hCounters[1]->GetBinContent(1);

  ifstream list("../prod/test2l2q/goodDY.txt");
  char fileName[300]; 
  char filestring[400];
  NGenEvt[2] = 0;
  NGenEvt[3] = 0;
  NGenEvt[4] = 0;
 
  while (list >> fileName) {
    sprintf(filestring,"../prod/test2l2q/%s",fileName); 
    TFile* ftemp = new TFile(filestring);
    cout << filestring << endl;
    if (string(fileName).find("DY2") != std::string::npos) {
      inputTree[2]->Add(filestring);
      hCounters[2] = (TH1F*)ftemp->Get("ZZTree/Counters");
      NGenEvt[2] += hCounters[2]->GetBinContent(1);
    }
    else if (string(fileName).find("DY3") != std::string::npos) {
      inputTree[3]->Add(filestring);
      hCounters[3] = (TH1F*)ftemp->Get("ZZTree/Counters");
      NGenEvt[3] += hCounters[3]->GetBinContent(1);
    }
    else if (string(fileName).find("DY4") != std::string::npos) {
      inputTree[4]->Add(filestring);
      hCounters[4] = (TH1F*)ftemp->Get("ZZTree/Counters");
      NGenEvt[4] += hCounters[4]->GetBinContent(1);
    }
  }
  
  for (int d=0; d<nDatasets; d++) {
	
    inputTree[d]->SetBranchAddress("overallEventWeight", &overallEventWeight);
    inputTree[d]->SetBranchAddress("ZZMass", &ZZMass);
    inputTree[d]->SetBranchAddress("ZZMassRefit", &ZZMassRefit);
    inputTree[d]->SetBranchAddress("ZZPt", &ZZPt);
    inputTree[d]->SetBranchAddress("Z1Mass", &Z1Mass);
    inputTree[d]->SetBranchAddress("Z2Mass", &Z2Mass);
    inputTree[d]->SetBranchAddress("Z2Flav", &Z2Flav);
    inputTree[d]->SetBranchAddress("LepPt", &LepPt);
    inputTree[d]->SetBranchAddress("JetPt", &JetPt);
    inputTree[d]->SetBranchAddress("JetIsInZZCand", &JetIsInZZCand);
    inputTree[d]->SetBranchAddress("JetBTagger", &JetBTagger);
    inputTree[d]->SetBranchAddress("JetQGLikelihood", &JetQGLikelihood);
    inputTree[d]->SetBranchAddress("helcosthetaZ1",&helcosthetaZ1);  
    inputTree[d]->SetBranchAddress("helcosthetaZ2",&helcosthetaZ2);  
    inputTree[d]->SetBranchAddress("costhetastar",&costhetastar);	  	
    inputTree[d]->SetBranchAddress("helphi",	  &helphi);	  
    inputTree[d]->SetBranchAddress("phistarZ1",   &phistarZ1); 
    inputTree[d]->SetBranchAddress("xsec", &xsec);

    //---------- Process tree

    Long64_t entries = inputTree[d]->GetEntries();
    cout<<"Processing dataset "<<d<<" ("<<entries<<" entries of "<< NGenEvt[d] << ") ..."<<endl;

    for (Long64_t z=0; z<entries; ++z){

      // if(DEBUG && z>1000) break;

      inputTree[d]->GetEntry(z);

      int process;
      if (d==0) process=1;
      else if (d==1) process=2;
      else process=3;
      
      Double_t eventWeight = ( lumin * 1000 * scaleF[process] / NGenEvt[d] ) * xsec ;// * overallEventWeight ;

      //----- fill histograms

      int fsstart,fsend;
      if (abs(Z2Flav)==121) {
	fsstart=0;
	fsend=2;
      } else {
	fsstart=1;
	fsend=3;
      }

      // --cuts (change here)
      float pt1stJet = 0.;
      for (unsigned int nJet=0; nJet<JetPt->size(); nJet++) {
	if (JetIsInZZCand->at(nJet)) {
	  if (pt1stJet < JetPt->at(nJet)) pt1stJet = JetPt->at(nJet);
	}
      }
      if (pt1stJet < 150.) continue;
      
      float pt1stLep = 0.;
      for (unsigned int nLep=0; nLep<LepPt->size(); nLep++) {
	if (pt1stLep < LepPt->at(nLep)) pt1stLep = LepPt->at(nLep);
      }
      if (pt1stLep < 150.) continue;

      for(int rs=fsstart; rs<fsend; rs++){

	h1[0][process][rs]->Fill(ZZMass,eventWeight);
	// cout << process << " " << rs << " " << eventWeight << endl;
	h1[1][process][rs]->Fill(ZZPt,eventWeight);
        h1[2][process][rs]->Fill(ZZMassRefit,eventWeight);
	h1[3][process][rs]->Fill(Z1Mass,eventWeight);
        h1[4][process][rs]->Fill(Z2Mass,eventWeight);
        h1[5][process][rs]->Fill(Z2Flav,eventWeight);
        
        int nExtraJet = 0;
        float pt1stJet = 0.;
        for (unsigned int nJet=0; nJet<JetPt->size(); nJet++) {
	  if (JetIsInZZCand->at(nJet)) {
            h1[9][process][rs]->Fill(JetBTagger->at(nJet),eventWeight);
            h1[8][process][rs]->Fill(JetQGLikelihood->at(nJet),eventWeight);
	    if (pt1stJet < 0.0005) pt1stJet = JetPt->at(nJet);
	    else {
	      if (JetPt->at(nJet) > pt1stJet) {
		h1[6][process][rs]->Fill(JetPt->at(nJet),eventWeight);
                h1[7][process][rs]->Fill(pt1stJet,eventWeight);
	      } else {
		h1[7][process][rs]->Fill(JetPt->at(nJet),eventWeight);
                h1[6][process][rs]->Fill(pt1stJet,eventWeight);
	      }
	    }
	  } else nExtraJet++; 
	}
	
        float pt1stLep = 0.;
        for (unsigned int nLep=0; nLep<LepPt->size(); nLep++) {
	  if (pt1stLep < 0.0005) pt1stLep = LepPt->at(nLep);
	  else {
	    if (LepPt->at(nLep) > pt1stLep) {
	      h1[10][process][rs]->Fill(LepPt->at(nLep),eventWeight);
	      h1[11][process][rs]->Fill(pt1stLep,eventWeight);
	    } else {
	      h1[11][process][rs]->Fill(LepPt->at(nLep),eventWeight);
	      h1[10][process][rs]->Fill(pt1stLep,eventWeight);
	    }
	  }
	}
  
	h1[12][process][rs]->Fill(abs(helcosthetaZ1),eventWeight);
	h1[13][process][rs]->Fill(helcosthetaZ2,eventWeight);
        h1[14][process][rs]->Fill(costhetastar,eventWeight);
        h1[15][process][rs]->Fill(helphi,eventWeight);
        h1[16][process][rs]->Fill(phistarZ1,eventWeight);
	h1[17][process][rs]->Fill(nExtraJet,eventWeight);
      }
    }
		 
  }
  
  TCanvas c1;
  c1.cd();

  for(int rs=0; rs<nFS; rs++){   //ee, mumu, or all
    for(int v=0; v<nVariables; v++){
      h1[v][3][rs]->GetXaxis()->SetTitle(varXLabel[v].c_str());
      h1[v][3][rs]->GetYaxis()->SetTitle(varYLabel[v].c_str());
      h1[v][3][rs]->SetFillStyle(1);
      h1[v][3][rs]->SetMinimum(0.);
      h1[v][3][rs]->SetLineColor(kGreen+2);
      h1[v][3][rs]->SetFillColor(kGreen+2);
      h1[v][2][rs]->SetLineColor(kBlue-1);
      h1[v][1][rs]->SetLineColor(kRed-1);

      TLegend *legend = new TLegend(0.70,0.75,0.95,0.90,NULL,"brNDC");
      legend->SetBorderSize(     0);
      legend->SetFillColor (     0);
      legend->SetTextAlign (    12);
      legend->SetTextFont  (    42);
      legend->SetTextSize  (0.03);
      
      for(int ipr=1; ipr<nProcesses; ipr++){ 
	legend->AddEntry(h1[v][ipr][rs], processLabel[ipr].c_str() , "l");
      }
      c1.cd();

      // cout << rs << " " << h1[v][3][rs]->GetEntries() << endl;
      h1[v][3][rs]->Draw("hist");
      h1[v][2][rs]->Draw("histsame"); 
      h1[v][1][rs]->Draw("histsame");
      legend->Draw("same"); 
      c1.SaveAs(Form("~/www/graviton/%s/%s_%s.png",dirout.c_str(),varName[v].c_str(),sFS[rs].c_str()));
    }
  }

}
