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

bool useHTBinned = false;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


const int nVariables = 21;
string varName[nVariables] = {
  "ZZMass",
  "ZZPt",
  "ZZMassRefit",
  "Z1Mass",
  "Z2Mass",
  "Z1Pt",
  "Z2Pt",
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
  "NExtraJets",
  "MET"
};
string varXLabel[nVariables] = {
  "m_{2#font[12]{l}2q} (GeV)",
  "p_{T2#font[12]{l}2q} (GeV)",
  "m_{2#font[12]{l}2q} (GeV)",
  "m_{jj} (GeV)",
  "m_{#font[12]{l}#font[12]{l}} (GeV)",
  "p_{T,jj} (GeV)",
  "p_{T,#font[12]{l}#font[12]{l}} (GeV)",
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
  "N_{extra-jets}",
  "MET (GeV)"
};
string varYLabel[nVariables] = {
  "Events / 25 GeV",
  "Events / 10 GeV",
  "Events / 25 GeV",
  "Events / 2.5 GeV",
  "Events / 2.5 GeV",
  "Events / 10 GeV",
  "Events / 10 GeV",
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
  "Events",
  "Events / 6 GeV"
};
Int_t  varNbin[nVariables] = { 50, 50, 50,  44,  44, 50,50, 400,  50,  50,  50,  50,  50,  50,  50, 50, 50, 25, 25, 4, 50};
Float_t varMin[nVariables] = {  250,  0,  250,  40,  40,  90, 90, -200,  0, 0, -0.2, -0.2, 0,  0, -0.2, -1.2, -1.2, -3.15, -3.15, -0.5, 0. };
Float_t varMax[nVariables] = { 1500, 500, 1500, 150, 150, 800, 800, 0, 500, 500, 1.2, 1.2, 500, 500, 1.2, 1.2, 1.2 , 3.15, 3.15, 3.5, 300. };
Bool_t varLogx[nVariables] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
Bool_t varLogy[nVariables] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1};


enum Process {Data=0, BulkG=1, Spin0=2, DYjets=3, TTBar=4, Diboson=5};
const int nProcesses = 6;
string sProcess[nProcesses] = {"Data", "BulkG", "Spin0", "DY", "TT", "VV"};
string processLabel[nProcesses] = {"Data", "G^{*}(800)#rightarrowZZ (x50)", "H_{NWA}#rightarrowZZ (x50)", "Z + jets (HT > 100 GeV)", "ttbar", "WZ, ZZ"};
Float_t scaleF[nProcesses] = {1.,50.,50.,1.,1.,1.};

const int nFS = 3;
string sFS[nFS] = {"ee","all","mm"};

const int nType = 4;
string typeS[nType] = {"resolvedSB","mergedSB","mergedSR","resolvedSR"};

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
  const int nDatasets = 9;
  
  TFile* inputFile[nDatasets];
  TChain* inputTree[nDatasets];
  TH1F* hCounters[nDatasets]; 
  Long64_t NGenEvt[nDatasets];
  string Dsname[nDatasets] = {"BulkGrav800","Higgs750","DYHT100","DYHT200","DYHT400","DYHT600","TTBar","WZDib","ZZDib"};
  // Float_t partialEventWeight[nDatasets];

  if (!useHTBinned) {
    Dsname[2] = "DYJetsToLL";
    processLabel[3] = "Z + jets";
  }

  Float_t overallEventWeight;
  vector<Float_t> *ZZMass = 0;
  vector<Float_t> *ZZMassRefit = 0;
  vector<Float_t> *ZZPt = 0;
  vector<Float_t> *Z1Mass = 0;
  vector<Float_t> *Z2Mass = 0;
  vector<Float_t> *Z1Pt = 0;
  vector<Float_t> *Z2Pt = 0;
  vector<Short_t> *Z2Flav = 0;
  vector<Short_t> *ZZCandType = 0;
  vector<Float_t> *LepPt = 0;
  Short_t nJets;
  Short_t nJetsBTagged;
  vector<Float_t> *JetPt = 0;
  vector<bool> *JetIsInZZCand = 0;
  vector<Float_t> *JetBTagger = 0;
  vector<Float_t> *JetQGLikelihood = 0;
  vector<Float_t> *helcosthetaZ1 = 0;  
  vector<Float_t> *helcosthetaZ2 = 0;  
  vector<Float_t> *costhetastar = 0;	  
  vector<Float_t> *helphi = 0;	  
  vector<Float_t> *phistarZ1 = 0;  
  Float_t xsec;
  Float_t Met;

  TH1F* h1[nVariables][nProcesses][nFS][nType];
  for(int rs=0; rs<nFS; rs++){   //ee, mumu, or all
    for(int pr=0; pr<nProcesses; pr++){
      for(int nt=0; nt<nType; nt++){
	for(int v=0; v<nVariables; v++){
	  h1[v][pr][rs][nt] = new TH1F(Form("h1_%s_%s_%s_%s",varName[v].c_str(),sFS[rs].c_str(),typeS[nt].c_str(),sProcess[pr].c_str()),
				       Form("h1_%s_%s_%s_%s",varName[v].c_str(),sFS[rs].c_str(),typeS[nt].c_str(),sProcess[pr].c_str()),		
				       varNbin[v],varMin[v],varMax[v]);
	}
      }
    }
  }

  //---------- Will loop over all datasets
  for (int d=0; d<nDatasets; d++) {
    NGenEvt[d] = 0;
    inputTree[d] = new TChain("ZZTree/candTree");
  }

  ifstream list("./goodDatasets.txt");
  char fileName[200];
  char filestring[200];
  while (list >> fileName) {
    sprintf(filestring,"root://eoscms//eos/cms/%s",fileName); 
    TFile* ftemp = TFile::Open(filestring);
    cout << filestring << endl;
    for (int d=0; d<nDatasets; d++) {
      if (string(fileName).find(Dsname[d].c_str()) != std::string::npos) {
	inputTree[d]->Add(filestring);
	hCounters[d] = (TH1F*)ftemp->Get("ZZTree/Counters");
	NGenEvt[d] += hCounters[d]->GetBinContent(1);
      }
    }
  }

  for (int d=0; d<nDatasets; d++) {

    if (!useHTBinned && d>2 && d<6) continue;   // in this case there is just one DY
	
    inputTree[d]->SetBranchAddress("overallEventWeight", &overallEventWeight);
    inputTree[d]->SetBranchAddress("ZZMass", &ZZMass);
    inputTree[d]->SetBranchAddress("ZZMassRefit", &ZZMassRefit);
    inputTree[d]->SetBranchAddress("ZZPt", &ZZPt);
    inputTree[d]->SetBranchAddress("Z1Mass", &Z1Mass);
    inputTree[d]->SetBranchAddress("Z2Mass", &Z2Mass);
    inputTree[d]->SetBranchAddress("Z1Pt", &Z1Pt);
    inputTree[d]->SetBranchAddress("Z2Pt", &Z2Pt);
    inputTree[d]->SetBranchAddress("Z2Flav", &Z2Flav);
    inputTree[d]->SetBranchAddress("ZZCandType", &ZZCandType);
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
    inputTree[d]->SetBranchAddress("PFMET", &Met);

    //---------- Process tree

    Long64_t entries = inputTree[d]->GetEntries();
    float eff = float(entries)/float(NGenEvt[d]);
    cout<<"Processing dataset "<<d<<" ("<<entries<<" entries of "<< NGenEvt[d] << " = " << eff*100. << "%) ..."<<endl;

    for (Long64_t z=0; z<entries; ++z){

      // cout<<"Processing entry "<<z<<endl;

      inputTree[d]->GetEntry(z);

      int process;
      if (d==0) process=1;
      else if (d==1) process=2;
      else if (d>1 && d<6) process=3;
      else if (d==6) process=4;
      else process=5;

      Double_t eventWeight = ( lumin * 1000 * scaleF[process] / NGenEvt[d] ) * xsec ;// * overallEventWeight ;
      if (z == 0) cout << "cross-section = " << xsec << " pb; eventweight = " << eventWeight << endl;       

      // find leading lepton and leading jet
      float pt1stJet = 0.0001;
      float pt2ndJet = 0.0001;
      int nInJets = 0;
      int nExtraJets = 0;
      for (unsigned int nJet=0; nJet<JetPt->size(); nJet++) {
	if (JetIsInZZCand->at(nJet)) {
	  if (pt1stJet < JetPt->at(nJet)) {
            pt2ndJet = pt1stJet;
	    pt1stJet = JetPt->at(nJet);
	  } else if (pt2ndJet < JetPt->at(nJet)) {
	    pt2ndJet = JetPt->at(nJet);
	  }
	  nInJets++;
	} else
	nExtraJets++;  
      } 
      /* if (nInJets!=2) cout << "nInJets " << nInJets << endl;
	 cout << "nExtraJets " << nExtraJets << endl;  */

      float pt1stLep = 0.0001;
      float pt2ndLep = 0.0001;
      for (unsigned int nLep=0; nLep<LepPt->size(); nLep++) {
	if (pt1stLep < LepPt->at(nLep)) {
	  pt2ndLep = pt1stLep;
	  pt1stLep = LepPt->at(nLep);
	} else if (pt2ndLep < LepPt->at(nLep)) {
	  pt2ndLep = LepPt->at(nLep);
	}
      }  

    
      // --cuts (change here)
      // if (pt1stJet < 150.) continue;
      if (pt1stLep < 150.) continue;
      // if (ZZPt < 20.) continue;
      // if (Z1Mass < 70.) continue; 

      //----- fill histograms

      int fsstart,fsend;
      if (abs(Z2Flav->at(0))==121) {
	fsstart=0;
	fsend=2;
      } else {
	fsstart=1;
	fsend=3;
      }

      for(int rs=fsstart; rs<fsend; rs++){
      
	for(unsigned int theCand=0; theCand<ZZMass->size(); theCand++){
	
	  int typ = ZZCandType->at(theCand)+2;
          if (typ>2) typ--;
         
          // if (z == 0) cout << ZZMass->size() << " " << ZZPt->size() << " " << ZZMassRefit->size() << " " << Z1Mass->size() << " " << Z2Mass->size() << " " << Z2Flav->size() << " " << helcosthetaZ1->size() << " " << helcosthetaZ2->size() << " " << costhetastar->size() << " "  << helphi->size() << " " << phistarZ1->size() << " " << endl;

	  h1[0][process][rs][typ]->Fill(ZZMass->at(theCand),eventWeight);
	  // cout << process << " " << rs << " " << eventWeight << endl;
	  h1[1][process][rs][typ]->Fill(ZZPt->at(theCand),eventWeight);
   	  if (typ==0 || typ==3) h1[2][process][rs][typ]->Fill(ZZMassRefit->at(theCand),eventWeight);   // only for resolved
	  h1[3][process][rs][typ]->Fill(Z1Mass->at(theCand),eventWeight);
	  h1[4][process][rs][typ]->Fill(Z2Mass->at(theCand),eventWeight);
          h1[5][process][rs][typ]->Fill(Z1Pt->at(theCand),eventWeight);
	  h1[6][process][rs][typ]->Fill(Z2Pt->at(theCand),eventWeight);
	  h1[7][process][rs][typ]->Fill(Z2Flav->at(theCand),eventWeight);
	  
	  h1[8][process][rs][typ]->Fill(pt1stJet,eventWeight);
	  h1[9][process][rs][typ]->Fill(pt2ndJet,eventWeight);

          if (typ==0 || typ==3) {
	    for (unsigned int nJet=0; nJet<JetPt->size(); nJet++) {
	      if (JetIsInZZCand->at(nJet) && JetBTagger->at(nJet) > -0.999) {
		h1[11][process][rs][typ]->Fill(JetBTagger->at(nJet),eventWeight);
		h1[10][process][rs][typ]->Fill(JetQGLikelihood->at(nJet),eventWeight); 
	      } 
	    }
	  }
	  
	  h1[12][process][rs][typ]->Fill(pt1stLep,eventWeight);
	  h1[13][process][rs][typ]->Fill(pt2ndLep,eventWeight);
	  h1[14][process][rs][typ]->Fill(abs(helcosthetaZ1->at(theCand)),eventWeight);
	  h1[15][process][rs][typ]->Fill(helcosthetaZ2->at(theCand),eventWeight);
	  h1[16][process][rs][typ]->Fill(costhetastar->at(theCand),eventWeight);
	  h1[17][process][rs][typ]->Fill(helphi->at(theCand),eventWeight);
	  h1[18][process][rs][typ]->Fill(phistarZ1->at(theCand),eventWeight);
	  h1[19][process][rs][typ]->Fill(nExtraJets,eventWeight);
	  h1[20][process][rs][typ]->Fill(Met,eventWeight);
	  
	}
      }		 
    }
    
  }

  TCanvas c1;
  c1.cd();
  
  for(int rs=0; rs<nFS; rs++){   //ee, mumu, or all
    for(int v=0; v<nVariables; v++){
      for(int nt=0; nt<nType; nt++){
      
        h1[v][4][rs][nt]->Add(h1[v][5][rs][nt]);
        h1[v][3][rs][nt]->Add(h1[v][4][rs][nt]);

	h1[v][3][rs][nt]->GetXaxis()->SetTitle(varXLabel[v].c_str());
	h1[v][3][rs][nt]->GetYaxis()->SetTitle(varYLabel[v].c_str());
	h1[v][3][rs][nt]->SetFillStyle(1);
	h1[v][3][rs][nt]->SetMinimum(0.1);
	h1[v][3][rs][nt]->SetLineColor(kGreen+2);
	h1[v][3][rs][nt]->SetFillColor(kGreen+2);
	h1[v][4][rs][nt]->SetFillStyle(1);
	h1[v][4][rs][nt]->SetMinimum(0.1);
	h1[v][4][rs][nt]->SetLineColor(kYellow+2);
	h1[v][4][rs][nt]->SetFillColor(kYellow+2);
        h1[v][5][rs][nt]->SetFillStyle(1);
	h1[v][5][rs][nt]->SetMinimum(0.1);
	h1[v][5][rs][nt]->SetLineColor(kMagenta-3);
	h1[v][5][rs][nt]->SetFillColor(kMagenta-3);
	h1[v][2][rs][nt]->SetLineColor(kBlue-1);
        h1[v][2][rs][nt]->SetLineWidth(3);
	h1[v][1][rs][nt]->SetLineColor(kRed-1);
        h1[v][1][rs][nt]->SetLineWidth(3);
	
	TLegend *legend = new TLegend(0.70,0.75,0.95,0.90,NULL,"brNDC");
	legend->SetBorderSize(     0);
	legend->SetFillColor (     0);
	legend->SetTextAlign (    12);
	legend->SetTextFont  (    42);
	legend->SetTextSize  (0.03);
	
	for(int ipr=1; ipr<nProcesses; ipr++){ 
	  legend->AddEntry(h1[v][ipr][rs][nt], processLabel[ipr].c_str() , "l");
	}
	c1.cd();
 
        
	// cout << rs << " " << h1[v][3][rs][nt]->GetEntries() << endl;
	h1[v][3][rs][nt]->Draw("hist");
	h1[v][4][rs][nt]->Draw("histsame"); 
	h1[v][5][rs][nt]->Draw("histsame");
	h1[v][2][rs][nt]->Draw("histsame"); 
	h1[v][1][rs][nt]->Draw("histsame");
        gPad->SetLogy(varLogy[v]);
	gPad->SetLogx(varLogx[v]);

	legend->Draw("same"); 
	c1.SaveAs(Form("~/www/graviton/%s/%s_%s_%s.png",dirout.c_str(),varName[v].c_str(),sFS[rs].c_str(),typeS[nt].c_str()));
      }
    }
  }

}
