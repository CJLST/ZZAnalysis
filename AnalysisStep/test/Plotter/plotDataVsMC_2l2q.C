
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
#include "TLine.h"
#include "TSpline.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include "tdrstyle.C"
// #include "CMS_lumi.C"
#include "plotUtils.C"
#include "ZZAnalysis/AnalysisStep/src/kFactors.C"

#include <ZZAnalysis/AnalysisStep/src/Category.cc>
#include <ZZAnalysis/AnalysisStep/src/bitops.cc>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>

using namespace std;

int useHTBinned = 2;         // 0 - use simple DY inclusive
                             // 1 - use ht binned
                             // 2 - use jet binned + b-enricchement

bool enforceNarrowWidth = true;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


const int nVariables = 25;
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
  "MET",
  "Z1Tau21",
  "BDT",
  "ZjetMELA",
  "vbfMELA"
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
  "MET (GeV)",
  "#tau_{21} (J)",
  "BDT output",
  "ZjetMELA",
  "vbfMELA",
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
  "Events / 0.028",
  "Events / 0.056",
  "Events / 20 GeV",
  "Events / 20 GeV",
  "Events",
  "Events",
  "Events",
  "Events",
  "Events",
  "Events",
  "Events / 6 GeV",
  "Events / 0.04",
  "Events / 0.025",
  "Events / 0.025", 
  "Events / 0.025",
};
Int_t  varNbin[nVariables] = { 50, 50, 50,  44,  44, 50,50, 400,  50,  50,  50,  50,  50,  50,  50, 50, 50, 25, 25, 4, 50, 26, 40, 40, 40};
Float_t varMin[nVariables] = {  250,  0,  250,  40,  40,  90, 90, -200,  0, 0, -0.2, -0.2, 0,  0, -0.2, -1.2, -1.2, -3.15, -3.15, -0.5, 0., -0.05,-0.3,0.,0.};
Float_t varMax[nVariables] = { 1500, 500, 1500, 150, 150, 800, 800, 0, 500, 500, 1.2, 1.2, 500, 500, 1.2, 1.2, 1.2 , 3.15, 3.15, 3.5, 300., 1.05, 0.3, 1.,1.};
Bool_t varLogx[nVariables] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
Bool_t varLogy[nVariables] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0};

const int nMasses = 13;
string signalMasses[nMasses] = {"200","250","300","350","400","450","500","550","600","750","800","1000","2000"};

enum Process {Data=0, BulkG=1, Spin0=2, DYjets=3, TTBar=4, Diboson=5};
const int nProcesses = 6;
string sProcess[nProcesses] = {"Data", "Spin0750", "Spin01000", "DY", "TT", "VV"};
string processLabel[nProcesses] = {"Data", "H_{NWA}(750)#rightarrowZZ (x100)", "H_{NWA}(1000)#rightarrowZZ (x100)", "Z + jets", "ttbar", "WZ, ZZ"};

//Float_t scaleF[nProcesses] = {1.,20.,20.,1.,1.,1.};

// WITH NNLO k-FACTOR FOR Z+JET
Float_t scaleF[nProcesses] = {1.,100.,100.,1.231,1.,1.};

const int nFS = 3;
string sFS[nFS] = {"ee","all","mm"};

const int nType = 12;
string typeS[nType] = {"resolvedSB","mergedSB","mergedSR","resolvedSR","resolvedSBbtag","mergedSBbtag","mergedSRbtag","resolvedSRbtag","resolvedSBvbf","mergedSBvbf","mergedSRvbf","resolvedSRvbf"};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void densityHist(TH1F* hist) 
{
  for (int i=1; i<=hist->GetNbinsX(); i++) {
    float binsize = hist->GetXaxis()->GetBinUpEdge(i) - hist->GetXaxis()->GetBinLowEdge(i);
    hist->SetBinContent(i,hist->GetBinContent(i)*50./binsize);
    hist->SetBinError(i,hist->GetBinError(i)*50./binsize);
  }
  hist->GetXaxis()->SetTitle(varXLabel[0].c_str());  
  hist->GetYaxis()->SetTitle("Events / 50 GeV");
}

void plotDataVsMC_2l2q(string dirout = "test13TeV", string theNtupleFile = "./goodDatasetsWithData.txt", bool norm = false, bool CR = false, bool draw = true)
{
  
  float lumin = 2.6;   // Moriond
  setTDRStyle();
  // gStyle->SetOptStat(1111111);
  const int nDatasets = 13;          // Moriond: 11
  const int nDatasetsMC = 11;         // Moriond: 9
  
  TFile* inputFile[nDatasets];
  TChain* inputTree[nDatasets];
  TH1F* hCounters[nDatasets]; 
  Long64_t NGenEvt[nDatasets];
  Float_t NEvtNarrow[nDatasets];
  Float_t sumWeights[nDatasets];
  Int_t mass[nDatasets];

  string Dsname[nDatasets] = {"ggHiggs750","ggHiggs1000","DY1JetsToLL","DY2JetsToLL","DY3JetsToLL","DY4JetsToLL","DYBJetsToLL","DYBFiltJetsToLL","TTBar","WZDib","ZZDib","DoubleEG2016B","DoubleMu2016B"};

  // string Dsname[nDatasets] = {"BulkGrav800","Higgs750","DY1JetsToLL","DY2JetsToLL","DY3JetsToLL","DY4JetsToLL","DYBJetsToLL","DYBFiltJetsToLL","TTBar","WZDib","ZZDib","DoubleEG2016B","DoubleMu2016B"};
  
  if (useHTBinned == 0) {
    Dsname[2] = "DYJetsToLL";
  } else if (useHTBinned == 1) {
    Dsname[2] = "DYHT100";
    Dsname[3] = "DYHT200";
    Dsname[4] = "DYHT400";
    Dsname[5] = "DYHT600";
    processLabel[3] = "Z + jets (HT > 100 GeV)"; 
  }

  /// I/O to TMVA 
  float tmvaZZPt, tmvaZ2Mass, tmvaZ1Pt, tmvaZ1tau21, tmvaZ2Pt, tmvaLepPt1, tmvaLepPt2, tmvaJetPt1, tmvaJetPt2;
  float tmvaJetQGLikelihood1, tmvaJetQGLikelihood2, tmvaabshelcosthetaZ1, tmvahelcosthetaZ2, tmvacosthetastar, tmvahelphi, tmvaphistarZ1; 

  TString outfileName( "TMVAAndRoofitInputs.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TTree* outputTree[4];   // ordered as 0 = sig merged
                          //            1 = bkg merged
                          //            2 = sig resolved       
                          //            3 = bkg resolved   
  outputTree[0] = new TTree("TreeSM","signal tree merged");
  outputTree[1] = new TTree("TreeBM","background tree merged");
  outputTree[2] = new TTree("TreeSR","signal tree resolved");
  outputTree[3] = new TTree("TreeBR","background tree resolved");

  for (int t=0; t<4; t++) {
    outputTree[t]->Branch("tmvaZZPt", &tmvaZZPt);
    outputTree[t]->Branch("tmvaZ2Mass", &tmvaZ2Mass);
    outputTree[t]->Branch("tmvaZ1tau21", &tmvaZ1tau21);
    outputTree[t]->Branch("tmvaZ1Pt", &tmvaZ1Pt);
    outputTree[t]->Branch("tmvaZ2Pt", &tmvaZ2Pt);
    outputTree[t]->Branch("tmvaLepPt1", &tmvaLepPt1);
    outputTree[t]->Branch("tmvaJetPt1", &tmvaJetPt1);
    outputTree[t]->Branch("tmvaLepPt2", &tmvaLepPt2);
    outputTree[t]->Branch("tmvaJetPt2", &tmvaJetPt2);
    outputTree[t]->Branch("tmvaJetQGLikelihood1", &tmvaJetQGLikelihood1);
    outputTree[t]->Branch("tmvaJetQGLikelihood2", &tmvaJetQGLikelihood2);
    outputTree[t]->Branch("tmvaabshelcosthetaZ1",&tmvaabshelcosthetaZ1);  
    outputTree[t]->Branch("tmvahelcosthetaZ2",&tmvahelcosthetaZ2);  
    outputTree[t]->Branch("tmvacosthetastar",&tmvacosthetastar);	  	
    outputTree[t]->Branch("tmvahelphi",	  &tmvahelphi);	  
    outputTree[t]->Branch("tmvaphistarZ1",   &tmvaphistarZ1);
  }

  TMVA::Reader *readerM = new TMVA::Reader( "!Color:!Silent" );

  readerM->AddVariable("tmvaZZPt", &tmvaZZPt);
  readerM->AddVariable("tmvaZ2Mass", &tmvaZ2Mass);
  readerM->AddVariable("tmvaZ1tau21", &tmvaZ1tau21);
  // readerM->AddVariable("tmvaZ1Pt", &tmvaZ1Pt);
  // readerM->AddVariable("tmvaZ2Pt", &tmvaZ2Pt);
  /* readerM->AddVariable("tmvaLepPt1", &tmvaLepPt1);
  readerM->AddVariable("tmvaJetPt1", &tmvaJetPt1);
  readerM->AddVariable("tmvaLepPt2", &tmvaLepPt2);
  readerM->AddVariable("tmvaJetPt2", &tmvaJetPt2); */
  readerM->AddVariable("tmvaabshelcosthetaZ1",&tmvaabshelcosthetaZ1);  
  readerM->AddVariable("tmvahelcosthetaZ2",&tmvahelcosthetaZ2);  
  readerM->AddVariable("tmvacosthetastar",&tmvacosthetastar);	  	
  readerM->AddVariable("tmvahelphi",	  &tmvahelphi);	  
  readerM->AddVariable("tmvaphistarZ1",   &tmvaphistarZ1);
  readerM->BookMVA( "BDT method" , "../TMVA/weights/merged/TMVAClassification_BDT.weights.xml" ); 
  
  TMVA::Reader *readerR = new TMVA::Reader( "!Color:!Silent" );

  readerR->AddVariable("tmvaZZPt", &tmvaZZPt);
  readerR->AddVariable("tmvaZ2Mass", &tmvaZ2Mass);
  // readerR->AddVariable("tmvaZ1Pt", &tmvaZ1Pt);
  // readerR->AddVariable("tmvaZ2Pt", &tmvaZ2Pt);
  /* readerR->AddVariable("tmvaLepPt1", &tmvaLepPt1);
  readerR->AddVariable("tmvaJetPt1", &tmvaJetPt1);
  readerR->AddVariable("tmvaLepPt2", &tmvaLepPt2);
  readerR->AddVariable("tmvaJetPt2", &tmvaJetPt2);  */
  readerR->AddVariable("tmvaJetQGLikelihood1", &tmvaJetQGLikelihood1);
  readerR->AddVariable("tmvaJetQGLikelihood2", &tmvaJetQGLikelihood2);
  readerR->AddVariable("tmvaabshelcosthetaZ1",&tmvaabshelcosthetaZ1);  
  readerR->AddVariable("tmvahelcosthetaZ2",&tmvahelcosthetaZ2);  
  readerR->AddVariable("tmvacosthetastar",&tmvacosthetastar);	  	
  readerR->AddVariable("tmvahelphi",	  &tmvahelphi);	  
  readerR->AddVariable("tmvaphistarZ1",   &tmvaphistarZ1);
  readerR->BookMVA( "BDT method" , "../TMVA/weights/resolved/TMVAClassification_BDT.weights.xml" );
  /// END I/O to TMVA

  Int_t RunNumber;
  Long64_t EventNumber;
  Int_t LumiNumber;
  Float_t genEventWeight;
  Float_t genHMass;
  vector<Float_t> *ZZMass = 0;
  vector<Float_t> *ZZMassRefit = 0;
  vector<Float_t> *ZZPt = 0;
  vector<Float_t> *Z1Mass = 0;
  vector<Float_t> *Z2Mass = 0;
  vector<Float_t> *Z1Pt = 0;
  vector<Float_t> *Z1tau21 = 0;
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
  vector<Float_t> *p0plus_VAJHU = 0; 
  vector<Float_t> *pqqZJJ_VAMCFM = 0;  
  vector<Float_t> *pvbf_VAJHU_highestPTJets = 0; 
  vector<Float_t> *phjj_VAJHU_highestPTJets = 0; 
  Float_t xsec;
  Float_t Met;

  TH1F* h1[nVariables][nProcesses][nFS][nType];
  for(int rs=0; rs<nFS; rs++){   //ee, mumu, or all
    for(int pr=0; pr<nProcesses; pr++){
      for(int nt=0; nt<nType; nt++){
	for(int v=0; v<nVariables; v++){
	  if (nt<4) 
	    h1[v][pr][rs][nt] = new TH1F(Form("h1_%s_%s_%s_%s",varName[v].c_str(),sFS[rs].c_str(),typeS[nt].c_str(),sProcess[pr].c_str()),
					 Form("h1_%s_%s_%s_%s",varName[v].c_str(),sFS[rs].c_str(),typeS[nt].c_str(),sProcess[pr].c_str()),		
					 varNbin[v],varMin[v],varMax[v]);
	  else
	    h1[v][pr][rs][nt] = new TH1F(Form("h1_%s_%s_%s_%s",varName[v].c_str(),sFS[rs].c_str(),typeS[nt].c_str(),sProcess[pr].c_str()),
					 Form("h1_%s_%s_%s_%s",varName[v].c_str(),sFS[rs].c_str(),typeS[nt].c_str(),sProcess[pr].c_str()),		
					 varNbin[v]/2,varMin[v],varMax[v]);
          h1[v][pr][rs][nt]->Sumw2();
	}
      }
    }
  }
  
  TH1F* hmass[nProcesses][nType];
  Float_t binsincl[] = { 450, 475, 500, 525, 550, 575, 600, 640, 680, 720, 760, 800, 850, 900, 1000, 1100, 1200, 1400, 1600, 2000};
  Int_t  binnumincl = sizeof(binsincl)/sizeof(Float_t) - 1; 
  Float_t binsbtag[] = { 450, 500, 550, 600, 680, 760, 850, 1000, 1200, 1600, 2000};
  Int_t  binnumbtag = sizeof(binsbtag)/sizeof(Float_t) - 1; 
  for(int pr=0; pr<nProcesses; pr++){
    for(int nt=0; nt<nType; nt++){
      if (nt<4) 
	 hmass[pr][nt] = new TH1F(Form("hmass_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
				  Form("hmass_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),		
			          binnumincl,binsincl);
      else
	 hmass[pr][nt] = new TH1F(Form("hmass_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
				  Form("hmass_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),		
				  binnumbtag,binsbtag);
      hmass[pr][nt]->Sumw2();
    }
  }

  //---------- Will loop over all datasets
  for (int d=0; d<nDatasets; d++) {
    if (d<nDatasets) { 
      NGenEvt[d] = 0;      NEvtNarrow[d] = 0.;
      sumWeights[d] = 0.;  mass[d] = 0;
    }
    inputTree[d] = new TChain("ZZTree/candTree");
  }

  ifstream list(theNtupleFile.c_str());
  char fileName[200];
  char filestring[200];
  while (list >> fileName) {
    if (string(fileName).find("store") != std::string::npos) sprintf(filestring,"root://eoscms//eos/cms/%s",fileName); 
    else sprintf(filestring,"%s",fileName);
    TFile* ftemp = TFile::Open(filestring);
    cout << filestring << endl;
    for (int d=0; d<nDatasets; d++) {
      if (string(fileName).find(Dsname[d].c_str()) != std::string::npos) {
	inputTree[d]->Add(filestring);
	if (d<nDatasetsMC) {
	  hCounters[d] = (TH1F*)ftemp->Get("ZZTree/Counters");
	  NGenEvt[d] += hCounters[d]->GetBinContent(1);
          sumWeights[d] += hCounters[d]->GetBinContent(41);
	}
	if (d<3) {
	  for (int m=0; m<nMasses; m++) {
	    if (string(fileName).find(signalMasses[m].c_str()) != std::string::npos) mass[d] = atoi(signalMasses[m].c_str());	     
	  }
	} 
      }
    }
    if (string(fileName).find("160601") != std::string::npos) {
      scaleF[1] = 30.;      scaleF[2] = 30.;   // forgot 2tau2q
    }
  }

  for (int d=0; d<nDatasets; d++) {

    if (useHTBinned == 0 && d>2 && d<8) continue;   // in this case there is just one DY
    if (useHTBinned == 1 && d>5 && d<8) continue;   // in this case there are just four DY

    inputTree[d]->SetBranchAddress("RunNumber", &RunNumber);
    inputTree[d]->SetBranchAddress("EventNumber", &EventNumber);
    inputTree[d]->SetBranchAddress("LumiNumber", &LumiNumber);
    inputTree[d]->SetBranchAddress("genHEPMCweight", &genEventWeight);
    inputTree[d]->SetBranchAddress("ZZMass", &ZZMass);
    inputTree[d]->SetBranchAddress("ZZMassRefit", &ZZMassRefit);
    inputTree[d]->SetBranchAddress("ZZPt", &ZZPt);
    inputTree[d]->SetBranchAddress("Z1Mass", &Z1Mass);
    inputTree[d]->SetBranchAddress("Z1tau21", &Z1tau21);
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
    inputTree[d]->SetBranchAddress("GenHMass", &genHMass);
    inputTree[d]->SetBranchAddress("PFMET", &Met);
    inputTree[d]->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU );
    inputTree[d]->SetBranchAddress("pqqZJJ_VAMCFM", &pqqZJJ_VAMCFM );  
    inputTree[d]->SetBranchAddress("pvbf_VAJHU_highestPTJets", &pvbf_VAJHU_highestPTJets );
    inputTree[d]->SetBranchAddress("phjj_VAJHU_highestPTJets", &phjj_VAJHU_highestPTJets );

    //---------- Process tree

    Long64_t entries = inputTree[d]->GetEntries();

    int process;
    if (d==0) process=1;
    else if (d==1) process=2;
    else if (d>1 && d<8) process=3;
    else if (d==8) process=4;
    else if (d>8 && d<11) process=5;
    else process=0;

    // for synchronization
    ofstream myfile;
    int nPassMerged = 0 ;
    int nPassResol = 0 ;
    if (process == 2) myfile.open("synchronization.txt");
    
    if (process>0) {
    float eff = float(entries)/float(NGenEvt[d]);
    cout<<"Processing dataset "<<d<<" ("<<entries<<" entries of "<< NGenEvt[d] << " = " << eff*100. << "%) ..."<<endl;
    } else {
      cout<<"Processing dataset "<<d<<" ("<<entries<<" entries)"<<endl;
    }

    for (Long64_t z=0; z<entries; ++z){

      // cout<<"Processing entry "<<z<<endl;

      inputTree[d]->GetEntry(z);
      bool writeThis = (z<100 && process==2); 

      Double_t eventWeight = 1. ;
      if (process>0) eventWeight = ( lumin * 1000 * scaleF[process]  / NGenEvt[d] ) * xsec ;
      if (process > 0 && z == 0) cout << "cross-section = " << xsec << " pb; eventweight = " << eventWeight << endl;       
      if (d == 7)  eventWeight *= (NGenEvt[d]*genEventWeight/sumWeights[d]);
      // use event weight for multiple hadronization sample

      // keep only events around nominal mass!!!   
      if (enforceNarrowWidth && mass[d] > 0 && fabs(genHMass-(float)mass[d]) > 0.01*mass[d]) continue;
      NEvtNarrow[d] = NEvtNarrow[d] + 1.;
      
      // find leading jets (notice this vector also includes subjets (identified by a btagger value of -1) which must be treated apart
      float pt1stJet = 0.0001;
      float pt2ndJet = 0.0001;
      float btag1stJet = 0.;
      float btag2ndJet = 0.;
      float qglik1stJet = 0.;
      float qglik2ndJet = 0.;

      float pt1stSubjet = 0.0001;
      float pt2ndSubjet = 0.0001;
      float btag1stSubjet = 0.;
      float btag2ndSubjet = 0.;
 
      int nInJets = 0;
      int nExtraJets = 0;

      for (unsigned int nJet=0; nJet<JetPt->size(); nJet++) {
	if (JetQGLikelihood->at(nJet) > -800.) {            // real jets
 	  if (JetIsInZZCand->at(nJet) ) {         
	    if (pt1stJet < JetPt->at(nJet)) {
	      pt2ndJet = pt1stJet;
	      pt1stJet = JetPt->at(nJet);
	      btag2ndJet = btag1stJet;
	      btag1stJet = JetBTagger->at(nJet);
	      qglik2ndJet = qglik1stJet;
	      qglik1stJet = JetQGLikelihood->at(nJet);
	    } else if (pt2ndJet < JetPt->at(nJet)) {
	      pt2ndJet = JetPt->at(nJet);
              btag2ndJet = JetBTagger->at(nJet);
	      qglik2ndJet = JetQGLikelihood->at(nJet);
	    }
	    nInJets++;
	  } else nExtraJets++;  
	}
      } 

      for (unsigned int nJet=0; nJet<JetPt->size(); nJet++) {
	if (JetQGLikelihood->at(nJet) < -800.) {            // subjets
	  if (JetIsInZZCand->at(nJet) ) {         
	    if (pt1stSubjet < JetPt->at(nJet)) {
	      pt2ndSubjet = pt1stSubjet;
	      pt1stSubjet = JetPt->at(nJet);
	      btag2ndSubjet = btag1stSubjet;
	      btag1stSubjet = JetBTagger->at(nJet);
	    } else if (pt2ndSubjet < JetPt->at(nJet)) {
	      pt2ndSubjet = JetPt->at(nJet);
              btag2ndSubjet = JetBTagger->at(nJet);
	    }
	  }
	}
      } 
        
      // find leading leptons
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
   
      //----- fill histograms

      int fsstart,fsend;
      if (abs(Z2Flav->at(0))==121) {
	fsstart=0;
	fsend=2;
      } else {
	fsstart=1;
	fsend=3;
      }

      int preferType = 0;    // choose merged or resolved
      // dump for synchronization 
      if (ZZMass->size() == 1 && abs(ZZCandType->at(0)) == 1) {
	if (writeThis) myfile << RunNumber << ":" << EventNumber << ":" << LumiNumber << ":" << Z2Mass->at(0)  << ":" << (abs(Z2Flav->at(0))==121 ? "Ele:" : "Muo:")  << pt1stLep << ":" << pt2ndLep << ":" << ZZMass->at(0) << ":" << Z1Mass->at(0) << ":" << Z1tau21->at(0) << ":" << Z1Pt->at(0) << ":-1:-1:-1:-1:-1:" << Met << endl; 
	nPassMerged++;
	preferType = 1;
      }
      else if (ZZMass->size() == 1 && abs(ZZCandType->at(0)) == 2) {
	if (writeThis) myfile << RunNumber << ":" << EventNumber << ":" << LumiNumber << ":" << Z2Mass->at(0)  << ":" << (abs(Z2Flav->at(0))==121 ? "Ele:" : "Muo:")  << pt1stLep << ":" << pt2ndLep << ":-1:-1:-1:-1:" << ZZMass->at(0) << ":" << Z1Mass->at(0) << ":" << pt1stJet << ":" << pt2ndJet << ":" << btag1stJet << ":" << Met << endl;
	nPassResol++;
        preferType = 2;
      }
      else if (ZZMass->size() == 2 && abs(ZZCandType->at(0)) == 1) {
	if (writeThis) myfile << RunNumber << ":" << EventNumber << ":" << LumiNumber << ":" << Z2Mass->at(0)  << ":" << (abs(Z2Flav->at(0))==121 ? "Ele:" : "Muo:")  << pt1stLep << ":" << pt2ndLep << ":" << ZZMass->at(0) << ":" << Z1Mass->at(0) << ":" << Z1tau21->at(0) << ":" << Z1Pt->at(0) << ":" << ZZMass->at(1) << ":" << Z1Mass->at(1) << ":" << pt1stJet << ":" << pt2ndJet << ":" << btag1stJet << ":" << Met << endl;
	nPassMerged++; nPassResol++;
        if (ZZCandType->at(0) > 0) preferType = 2;
        else if (ZZCandType->at(1) > 0) preferType = 1;
        else preferType = 2;
      }
      else if (ZZMass->size() == 2 && abs(ZZCandType->at(0)) == 2) {
	if (writeThis) myfile << RunNumber << ":" << EventNumber << ":" << LumiNumber << ":" << Z2Mass->at(0)  << ":" << (abs(Z2Flav->at(0))==121 ? "Ele:" : "Muo:")  << pt1stLep << ":" << pt2ndLep << ":" << ZZMass->at(1) << ":" << Z1Mass->at(1) << ":" << Z1tau21->at(1) << ":" << Z1Pt->at(1) << ":" << ZZMass->at(0) << ":" << Z1Mass->at(0) << ":" << pt1stJet << ":" << pt2ndJet << ":" << btag1stJet << ":" << Met << endl;
	nPassMerged++; nPassResol++;
        if (ZZCandType->at(1) > 0) preferType = 2;
        else if (ZZCandType->at(0) > 0) preferType = 1;
        else preferType = 2;
      }
    
      // end dump for synchronization

      for(int rs=fsstart; rs<fsend; rs++){
      
	for(unsigned int theCand=0; theCand<Z1Mass->size(); theCand++){
	
	  if (!CR) {   // regular events 
            if (abs(ZZCandType->at(theCand)) != preferType) continue;
	    int typ = ZZCandType->at(theCand)+2;
	    if (typ>2) typ--;
	    
	    int whichTmvaTree = -1;
	    if (typ==2 && process == 2) whichTmvaTree = 0;
	    if (typ==2 && process == 3) whichTmvaTree = 1;
	    if (typ==3 && process == 2) whichTmvaTree = 2;
	    if (typ==3 && process == 3) whichTmvaTree = 3;
 
            float mela = 1./(1.+0.035*(pqqZJJ_VAMCFM->at(theCand)/p0plus_VAJHU->at(theCand)));
            float vbfmela = 1./(1.+0.06*(phjj_VAJHU_highestPTJets->at(theCand)/pvbf_VAJHU_highestPTJets->at(theCand)));

	    if ((typ==0 || typ==3) && nExtraJets > 1 && vbfmela > 0.5) typ=typ+8;
	    if ((typ==1 || typ==2) && nExtraJets > 1 && vbfmela > 0.5) typ=typ+8;             
  	    
	    if ((typ==0 || typ==3) && btag1stJet > 0.605 && btag2ndJet > 0.605) typ=typ+4;
	    if ((typ==1 || typ==2) && btag1stSubjet > 0.605 && btag2ndSubjet > 0.605) typ=typ+4;
	    
	    tmvaZZPt = (float)ZZPt->at(theCand); 
	    tmvaZ2Mass = (float)Z2Mass->at(theCand);
	    tmvaZ1Pt = (float)Z1Pt->at(theCand);
	    tmvaZ2Pt = (float)Z2Pt->at(theCand);	  
	    tmvaLepPt1 = (float)pt1stLep;
	    tmvaLepPt2 = (float)pt2ndLep;
	    tmvaJetQGLikelihood1 = (float)qglik1stJet;
	    tmvaJetQGLikelihood2 = (float)qglik2ndJet;
	    tmvaabshelcosthetaZ1 = (float)abs(helcosthetaZ1->at(theCand)); 
	    tmvahelcosthetaZ2 = (float)helcosthetaZ2->at(theCand);
	    tmvacosthetastar = (float)costhetastar->at(theCand);
	    tmvahelphi = (float)helphi->at(theCand);
	    tmvaphistarZ1 = (float)phistarZ1->at(theCand);
	    tmvaZ1tau21 = (float)Z1tau21->at(theCand); 
	    if (typ==0 || typ==3 || typ==4 || typ==7 || typ==8 || typ==11) {   // only resolved
	      tmvaJetPt1 = (float)pt1stJet;
	      tmvaJetPt2 = (float)pt2ndJet;
	    } else {
	      tmvaJetPt1 = (float)pt1stSubjet;
	      tmvaJetPt2 = (float)pt2ndSubjet;
	    }              
	    float bdt;
	    if (typ==0 || typ==3 || typ==4 || typ==7 || typ==8 || typ==11) bdt = readerR->EvaluateMVA( "BDT method" );
	    else bdt = readerM->EvaluateMVA( "BDT method" );
	    
	    // test MET and BDT cut 
	    // if (typ < 4 && bdt < 0.08) continue;
	    // if (typ > 3 && Met > 50.) continue;
	    
	    // test basic cuts 
	    // if (tmvaJetPt1 < 150.) continue;
	    // if (tmvaJetPt2 < 40.) continue;
	    // if (tmvaLepPt1 < 150.) continue;
	    // if (tmvaLepPt2 < 30.) continue;	    
	    if (tmvaZ2Mass < 60.) continue;
            if ((typ==1 || typ==2 || typ==5 || typ==6 || typ==9 || typ==10) && tmvaZ1tau21 > 0.6) continue;
            // if (mela < 0.8) continue;
	    
	    h1[0][process][rs][typ]->Fill(ZZMass->at(theCand),eventWeight);
	    h1[1][process][rs][typ]->Fill(ZZPt->at(theCand),eventWeight);
	    
	    h1[2][process][rs][typ]->Fill(ZZMassRefit->at(theCand),eventWeight);
            if (rs == 1) {
	      if (typ==0 || typ==3 || typ==4 || typ==7 || typ==8 || typ==11) hmass[process][typ]->Fill(ZZMassRefit->at(theCand),eventWeight);
	      else hmass[process][typ]->Fill(ZZMass->at(theCand),eventWeight);
            }
 
	    h1[3][process][rs][typ]->Fill(Z1Mass->at(theCand),eventWeight);
	    
	    h1[4][process][rs][typ]->Fill(Z2Mass->at(theCand),eventWeight);
	    h1[5][process][rs][typ]->Fill(Z1Pt->at(theCand),eventWeight);
	    h1[6][process][rs][typ]->Fill(Z2Pt->at(theCand),eventWeight);
	    h1[7][process][rs][typ]->Fill(Z2Flav->at(theCand),eventWeight);
	    
	    if (typ==0 || typ==3 || typ==4 || typ==7 || typ==8 || typ==11) {   // only resolved
	      h1[8][process][rs][typ]->Fill(pt1stJet,eventWeight);
	      h1[9][process][rs][typ]->Fill(pt2ndJet,eventWeight);
	    } else {
	      h1[8][process][rs][typ]->Fill(pt1stSubjet,eventWeight);
	      h1[9][process][rs][typ]->Fill(pt2ndSubjet,eventWeight);
	    }
	    
	    if (typ==0 || typ==3 || typ==4 || typ==7  || typ==8 || typ==11) {   // only resolved
	      for (unsigned int nJet=0; nJet<JetPt->size(); nJet++) {
		if (JetIsInZZCand->at(nJet) && JetQGLikelihood->at(nJet) > -800. /*remove subjets of fat jet also included in this collection! */) {
		  h1[11][process][rs][typ]->Fill(JetBTagger->at(nJet),eventWeight);
		  h1[10][process][rs][typ]->Fill(JetQGLikelihood->at(nJet),eventWeight); 
		} 
	      }
	    } else {  
	      for (unsigned int nJet=0; nJet<JetPt->size(); nJet++) {
		if (JetIsInZZCand->at(nJet) && JetQGLikelihood->at(nJet) < -800. /*subjets of fat jet also included in this collection! */) {
		  h1[11][process][rs][typ]->Fill(JetBTagger->at(nJet),eventWeight);
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

	    if (typ==1 || typ==2 || typ==5 || typ==6  || typ==9 || typ==10) h1[21][process][rs][typ]->Fill(Z1tau21->at(theCand),eventWeight);   // only merged
	    h1[22][process][rs][typ]->Fill(bdt,eventWeight);
            h1[23][process][rs][typ]->Fill(mela,eventWeight);
            if (nExtraJets > 1) h1[24][process][rs][typ]->Fill(vbfmela,eventWeight);
	    
	    if (rs == 1 && whichTmvaTree > -1) outputTree[whichTmvaTree]->Fill();
	  } else { // control region for QG (only fille some variables)

	    h1[3][process][rs][3]->Fill(Z1Mass->at(theCand),eventWeight);
	    h1[4][process][rs][3]->Fill(Z2Mass->at(theCand),eventWeight);
	    h1[5][process][rs][3]->Fill(Z1Pt->at(theCand),eventWeight);
	    h1[6][process][rs][3]->Fill(Z2Pt->at(theCand),eventWeight);
	    h1[7][process][rs][3]->Fill(Z2Flav->at(theCand),eventWeight);
	    
	    h1[8][process][rs][3]->Fill(pt1stJet,eventWeight);
	    	    
	    for (unsigned int nJet=0; nJet<JetPt->size(); nJet++) {
	      if (JetIsInZZCand->at(nJet) && JetQGLikelihood->at(nJet) > -800. /*remove subjets of fat jet also included in this collection! */) {
		h1[11][process][rs][3]->Fill(JetBTagger->at(nJet),eventWeight);
		h1[10][process][rs][3]->Fill(JetQGLikelihood->at(nJet),eventWeight); 
	      }
	    }

	    h1[12][process][rs][3]->Fill(pt1stLep,eventWeight);
	    h1[13][process][rs][3]->Fill(pt2ndLep,eventWeight);
	  }
	}
      }		 
    }
    if (mass[d] > 0 && enforceNarrowWidth) NEvtNarrow[d] = NEvtNarrow[d] / entries;
    if (process==2 && !draw) { 
      myfile.close(); 
      float eff = float(nPassMerged)/float(NGenEvt[d]);
      cout<<"Pass merged analysis = "<<nPassMerged<<"/"<< NGenEvt[d] << " (" << eff*100. << "%)"<<endl;
      eff = float(nPassResol)/float(NGenEvt[d]);
      cout<<"Pass resolved analysis = "<<nPassResol<<"/"<< NGenEvt[d] << " (" << eff*100. << "%) ..."<<endl;
      break;
    }  
  }

  if (!draw) return; 
  TCanvas c1;
  c1.cd();

  TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.35,0.95,0.97);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.02,0.95,0.35);
  pad2->Draw();
 
  
  for(int rs=0; rs<nFS; rs++){   //ee, mumu, or all
    for(int v=0; v<nVariables; v++){
      for(int nt=0; nt<nType; nt++){

	if (enforceNarrowWidth) {
          cout << "Only selecting " << NEvtNarrow[0]*100. << "% events close to nominal mass (mimic narrow width)" << endl;
	  cout << "Only selecting " << NEvtNarrow[1]*100. << "% events close to nominal mass (mimic narrow width)" << endl;
	  h1[v][1][rs][nt]->Scale(1./NEvtNarrow[0]);
	  h1[v][2][rs][nt]->Scale(1./NEvtNarrow[1]);
	}

	float normal = 1.;
	if (norm) {
	  if (h1[v][0][rs][nt]->Integral() > 0) normal = (float)h1[v][3][rs][nt]->Integral()/(float)h1[v][0][rs][nt]->Integral();
	  h1[v][0][rs][nt]->Scale(normal);
	  if (h1[v][2][rs][nt]->Integral() > 0) normal = (float)h1[v][3][rs][nt]->Integral()/(float)h1[v][2][rs][nt]->Integral();
	  h1[v][2][rs][nt]->Scale(normal);
	  if (h1[v][1][rs][nt]->Integral() > 0) normal = (float)h1[v][3][rs][nt]->Integral()/(float)h1[v][1][rs][nt]->Integral();
	  h1[v][1][rs][nt]->Scale(normal);
	}	


	pad1->cd();
      
        h1[v][4][rs][nt]->Add(h1[v][5][rs][nt]);
        h1[v][3][rs][nt]->Add(h1[v][4][rs][nt]);

        if (v==0 && rs==0) {    // do just once!
	  hmass[4][nt]->Add(hmass[5][nt]);
	}

	h1[v][3][rs][nt]->GetXaxis()->SetTitle(varXLabel[v].c_str());
	h1[v][3][rs][nt]->GetYaxis()->SetTitle(varYLabel[v].c_str());
        h1[v][0][rs][nt]->GetXaxis()->SetTitle(varXLabel[v].c_str());
	h1[v][0][rs][nt]->GetYaxis()->SetTitle(varYLabel[v].c_str());
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
        h1[v][0][rs][nt]->SetMarkerStyle(20);
        h1[v][0][rs][nt]->SetMinimum(0.1);
	
	TLegend *legend = new TLegend(0.70,0.75,0.95,0.90,NULL,"brNDC");
	legend->SetBorderSize(     0);
	legend->SetFillColor (     0);
	legend->SetTextAlign (    12);
	legend->SetTextFont  (    42);
	legend->SetTextSize  (0.03);
	
        legend->AddEntry(h1[v][0][rs][nt], processLabel[0].c_str() , "p");
	for(int ipr=1; ipr<nProcesses; ipr++){ 
          string legType = "f";
          if (ipr<3) legType = "l"; 
	  legend->AddEntry(h1[v][ipr][rs][nt], processLabel[ipr].c_str() , legType.c_str());
	}
	// c1.cd();
 
        
	// cout << rs << " " << h1[v][3][rs][nt]->GetEntries() << endl;
        if (nt==0 || nt==1 || nt==4 || nt==5 || nt==8 || nt==9 || CR) {

	  float a = h1[v][3][rs][nt]->GetMaximum();
	  float b = h1[v][0][rs][nt]->GetMaximum();
	  
	  if (a>b) {
	    h1[v][3][rs][nt]->Draw("hist");
	    h1[v][0][rs][nt]->Draw("esame");
	  } else {
	    h1[v][0][rs][nt]->Draw("e");
	    h1[v][3][rs][nt]->Draw("histsame");
	  }
	  
	  h1[v][4][rs][nt]->Draw("histsame"); 
	  h1[v][5][rs][nt]->Draw("histsame");
	  h1[v][2][rs][nt]->Draw("histsame"); 
	  h1[v][1][rs][nt]->Draw("histsame");

          h1[v][0][rs][nt]->Draw("esame");
	  
        } else {

	  h1[v][3][rs][nt]->Draw("hist");
	  h1[v][4][rs][nt]->Draw("histsame"); 
	  h1[v][5][rs][nt]->Draw("histsame");
	  h1[v][2][rs][nt]->Draw("histsame"); 
	  h1[v][1][rs][nt]->Draw("histsame");
        }

        gPad->SetLogy(varLogy[v]);
	gPad->SetLogx(varLogx[v]);

	legend->Draw("same"); 

	pad2->cd();

        gPad->SetLogy(0);
        TH1F* ratio = (TH1F*)h1[v][0][rs][nt]->Clone();
	ratio->Add(h1[v][0][rs][nt],h1[v][3][rs][nt],1,-1);
	ratio->Divide(ratio,h1[v][0][rs][nt]);
        ratio->SetMinimum(-0.5);
        ratio->SetMaximum(0.5); 
        ratio->GetYaxis()->SetTitle("(Data-MC)/Data");
        if (!(nt==0 || nt==1 || nt==4 || nt==5 || nt==8 || nt==9 || CR)) {
	  ratio->SetLineColor(kWhite);
	  ratio->SetMarkerColor(kWhite);
	}
	ratio->Draw("e");
        TLine line(h1[v][0][rs][nt]->GetXaxis()->GetBinLowEdge(1),0.,
		   h1[v][0][rs][nt]->GetXaxis()->GetBinUpEdge(h1[v][0][rs][nt]->GetNbinsX()-1),0.);
	line.SetLineColor(kRed);
        line.SetLineStyle(kDashed);
        line.Draw("same");  

	c1.SaveAs(Form("~/www/graviton/%s/%s_%s_%s.png",dirout.c_str(),varName[v].c_str(),sFS[rs].c_str(),typeS[nt].c_str()));
      }
    }
  }
  
  outputFile->cd();
  for (int t=0; t<4; t++) {
    outputTree[t]->Write();
  }
  for(int nt=0; nt<nType; nt++){
    // do histogram -> density histogram conversion
    densityHist(hmass[0][nt]);
    hmass[0][nt]->Write();
    densityHist(hmass[3][nt]);
    hmass[3][nt]->Write();
    densityHist(hmass[4][nt]);
    hmass[4][nt]->Write();
  }
  outputFile->Close(); 
  
}
