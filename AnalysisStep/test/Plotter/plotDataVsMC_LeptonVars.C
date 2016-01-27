/* 
 * usage: 
 * -specify parameters at the end of this file
 * -run with:
 *   root -l plotDataVsMC_LeptonVars.C++
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
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "plotUtils.C"

using namespace std;

#define DEBUG 0
#define DRAWLINES 0

#define DRAWDATAMCRATIO 1





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


const int nVariables = 14;
string varName[nVariables] = {
  "Nvtx",
  "Mll",
  "Pt1",
  "Pt2",
  "Eta1",
  "Eta2",
  "Phi1",
  "Phi2",
  "SIP1",
  "SIP2",
  "BDT1",
  "BDT2",
  "Iso1",
  "Iso2",
};
string varXLabel[nVariables] = {
  "number of primary vertices",
  "m_{#font[12]{l+l-}} [GeV]",
  "p_{T} (leading lepton) [GeV]",//"p_{T}^{lead.} (GeV)",
  "p_{T} (trailing lepton) [GeV]",//"p_{T}^{trail.} (GeV)",
  "#eta (leading lepton)",
  "#eta (trailing lepton)",
  "#phi (leading lepton)",
  "#phi (trailing lepton)",
  "SIP (leading lepton)",
  "SIP (trailing lepton)",
  "BDT (leading lepton)",
  "BDT (trailing lepton)",
  "isolation (leading lepton)",
  "isolation (trailing lepton)",
};
string varYLabel[nVariables] = {
  "Events",
  "Events / 1 GeV",
  "Events / 2 GeV",
  "Events / 2 GeV",
  "Events / 0.1",
  "Events / 0.1",
  "Events / 0.1",
  "Events / 0.1",
  "Events / 0.05",
  "Events / 0.05",
  "Events / 0.02",
  "Events / 0.02",
  "Events / 0.01",
  "Events / 0.01",
};
Bool_t plotThisVar[2][nVariables] = {
  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,},
  {1,1,1,1,1,1,1,1,1,1,1,1,1,1,},
};
Int_t  varNbin[nVariables] = { 50,  60, 100, 100, 60, 60, 80, 80, 90  , 90  , 100, 100, 55   , 55   , };
Float_t varMin[nVariables] = {  0,  60,   0,   0, -3, -3, -4, -4,  0  ,  0  ,  -1,  -1,  0   ,  0   , };
Float_t varMax[nVariables] = { 50, 120, 200, 200,  3,  3,  4,  4,  4.5,  4.5,   1,   1,  0.55,  0.55, };
Bool_t varLogx[nVariables] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,};
Bool_t varLogy[nVariables] = {1,1,1,1,0,0,0,0,1,1,1,1,1,1,};
Int_t restrictCountVar[nVariables] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,};
Float_t varMinFactor[nVariables] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,};
Int_t varCMSPos[nVariables] = {11,11,11,11,11,11,11,11,11,11,33,33,11,11,};
Int_t varLegPos[nVariables] = {33,33,33,33,33,33,33,33,33,33,11,11,33,33,};

const int nVarPairs = 1;
string varPairName[nVarPairs] = { "Pt1VsEta1", };
string varPairXLabel[nVarPairs] = { "p_{T} (leading lepton) [GeV]", };
string varPairYLabel[nVarPairs] = { "#eta (leading lepton)", };
Int_t  varPairXNbin[nVarPairs] = { 100, };
Float_t varPairXMin[nVarPairs] = {   0, };
Float_t varPairXMax[nVarPairs] = { 200, };
Int_t  varPairYNbin[nVarPairs] = { 60, };
Float_t varPairYMin[nVarPairs] = { -3, };
Float_t varPairYMax[nVarPairs] = {  3, };
Bool_t varPairLogx[nVarPairs] = {0,};
Bool_t varPairLogy[nVarPairs] = {0,};
Bool_t plotThisVarPair[2][nVarPairs] = {
  {0,}, 
  {1,}, 
};



enum Process {Data=0, DY=1, ttbar=2, DiBoson=3};
const int nProcesses = 4;
string sProcess[nProcesses] = {"Data", "DY", "ttbar", "DiBoson"};
string processLabel[nProcesses] = {"Data", "Z + jets", "t#bar{t}", "diboson"};
Int_t processFillColor[nProcesses] = {TColor::GetColor("#000000"), TColor::GetColor("#669966"), TColor::GetColor("#9a6666"), DRAWLINES?TColor::GetColor("#99ccff"):TColor::GetColor("#8bc5ff")};
Int_t processLineColor[nProcesses] = {TColor::GetColor("#000000"), TColor::GetColor("#003300"), TColor::GetColor("#5f3f3f"), TColor::GetColor("#000099")};
Bool_t useProcess[nProcesses] = {1,1,1,1,};

enum FinalState {fs2mu=0, fs2e=1};
const int nFinalStates = 2;
string sFinalState[nFinalStates+1] = {"2mu", "2e", "2l"};
string fsLabel[nFinalStates+1] = {"2#mu", "2e", "2#font[12]{l}"};





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void doHistograms(string inputFilePath_MC, string inputFilePath_Data, double lumi)
{

  const int nDatasets = 25;
  string datasets[nDatasets] = {
    "DoubleMu2015B",
    "DoubleEG2015B",
    "MuonEG2015B",
    "SingleEle2015B",
    "DoubleMu2015C",
    "DoubleEG2015C",
    "MuonEG2015C",
    "SingleEle2015C",
    "DoubleMu2015C_50ns",
    "DoubleEG2015C_50ns",
    "MuonEG2015C_50ns",
    "SingleEle2015C_50ns",
    "DoubleMu2015D",
    "DoubleEG2015D",
    "MuonEG2015D",
    "SingleEle2015D",
    "DoubleMu2015Dv4",
    "DoubleEG2015Dv4",
    "MuonEG2015Dv4",
    "SingleEle2015Dv4",
    "DYJetsToLL_M50",
    "DYJetsToLL_M10to50",
    "TTTo2L2Nu",
    "WWTo2L2Nu",
    "WZJets",
  };

  TFile* inputFile[nDatasets];
  TTree* inputTree[nDatasets];
  TH1F* hCounters[nDatasets];
  Long64_t NGenEvt[nDatasets];
  Double_t gen_sumWeights[nDatasets];
  Float_t partialSampleWeight[nDatasets];

  Int_t nRun;
  Long64_t nEvent;
  Int_t nLumi;
  Short_t Nvtx;
  Short_t NObsInt;
  Float_t NTrueInt;
  Float_t xsec;
  Float_t overallEventWeight;
  Short_t Zsel;
  Float_t ZMass;
  Float_t ZPt;
  Float_t ZEta;
  Float_t ZPhi;
  Short_t ZFlav;
  vector<Float_t> *LepPt = 0;
  vector<Float_t> *LepEta = 0;
  vector<Float_t> *LepPhi = 0;
  vector<Short_t> *LepLepId = 0;
  vector<Float_t> *LepSIP = 0;
  vector<Float_t> *LepTime = 0;
  vector<Bool_t> *LepisID = 0;
  vector<Float_t> *LepBDT = 0;
  vector<Char_t> *LepMissingHit = 0;
  vector<Float_t> *LepChargedHadIso = 0;
  vector<Float_t> *LepNeutralHadIso = 0;
  vector<Float_t> *LepPhotonIso = 0;
  vector<Float_t> *LepCombRelIsoPF = 0;
  vector<Float_t> *LepCombRelIsoPFPreFSR = 0;
  vector<Float_t> *fsrPt = 0;
  vector<Float_t> *fsrEta = 0; 
  vector<Float_t> *fsrPhi = 0;
  vector<Float_t> *fsrDR = 0;
  vector<Short_t> *fsrLept = 0;
  vector<Short_t> *fsrLeptID = 0;

  TH1F* h1[nVariables][nFinalStates+1][nProcesses];
  TH2F* h2[nVarPairs ][nFinalStates+1][nProcesses];
  for(int fs=0; fs<nFinalStates+1; fs++){
    for(int pr=0; pr<nProcesses; pr++){
      for(int v=0; v<nVariables; v++){
	h1[v][fs][pr] = new TH1F(Form("h1_%s_%s_%s",varName[v].c_str(),sFinalState[fs].c_str(),sProcess[pr].c_str()),
				 Form(";%s;%s",varXLabel[v].c_str(),varYLabel[v].c_str()),
				 varNbin[v],varMin[v],varMax[v]);
      }
      for(int v2=0; v2<nVarPairs; v2++){
	h2[v2][fs][pr] = new TH2F(Form("h2_%s_%s_%s",varPairName[v2].c_str(),sFinalState[fs].c_str(),sProcess[pr].c_str()),
				  Form(";%s;%s",varPairXLabel[v2].c_str(),varPairYLabel[v2].c_str()),
				  varPairXNbin[v2],varPairXMin[v2],varPairXMax[v2],
				  varPairYNbin[v2],varPairYMin[v2],varPairYMax[v2]);
      }
    }
  }
  
  int currentProcess;
  int currentFinalState;


  //---------- Will loop over all datasets

  for(int d=0; d<nDatasets; d++){

    //----- assign dataset to correct process
    currentProcess = -1;
    if(datasets[d].find("2015")!=string::npos)
      currentProcess = Data;
    if(datasets[d].find("DY")!=string::npos) 
      currentProcess = DY;
    if(datasets[d]=="WWTo2L2Nu"||datasets[d]=="WZJets")
      currentProcess = DiBoson;
    if(datasets[d]=="TTTo2L2Nu")
      currentProcess = ttbar;

    string inputFileName = string(Form("%s%s/ZZ4lAnalysis.root",(currentProcess==Data?inputFilePath_Data:inputFilePath_MC).c_str(),datasets[d].c_str()));
    inputFile[d] = TFile::Open(inputFileName.c_str());

    hCounters[d] = (TH1F*)inputFile[d]->Get("ZTree/Counters");
    gen_sumWeights[d] = (Long64_t)hCounters[d]->GetBinContent(1);
    partialSampleWeight[d] = lumi * 1000 / gen_sumWeights[d] ;

    inputTree[d] = (TTree*)inputFile[d]->Get("ZTree/candTree");

    inputTree[d]->SetBranchAddress("RunNumber",&nRun);
    inputTree[d]->SetBranchAddress("EventNumber",&nEvent);
    inputTree[d]->SetBranchAddress("LumiNumber",&nLumi);
    inputTree[d]->SetBranchAddress("Nvtx",&Nvtx);
    inputTree[d]->SetBranchAddress("NObsInt",&NObsInt);
    inputTree[d]->SetBranchAddress("NTrueInt",&NTrueInt);
    inputTree[d]->SetBranchAddress("xsec",&xsec);
    inputTree[d]->SetBranchAddress("overallEventWeight",&overallEventWeight);
    inputTree[d]->SetBranchAddress("Zsel",&Zsel);
    inputTree[d]->SetBranchAddress("ZMass",&ZMass);
    inputTree[d]->SetBranchAddress("ZPt",&ZPt);
    inputTree[d]->SetBranchAddress("ZEta",&ZEta);
    inputTree[d]->SetBranchAddress("ZPhi",&ZPhi);
    inputTree[d]->SetBranchAddress("ZFlav",&ZFlav);
    inputTree[d]->SetBranchAddress("LepPt",&LepPt);
    inputTree[d]->SetBranchAddress("LepEta",&LepEta);
    inputTree[d]->SetBranchAddress("LepPhi",&LepPhi);
    inputTree[d]->SetBranchAddress("LepLepId",&LepLepId);
    inputTree[d]->SetBranchAddress("LepSIP",&LepSIP);
    inputTree[d]->SetBranchAddress("LepTime",&LepTime);
    inputTree[d]->SetBranchAddress("LepisID",&LepisID);
    inputTree[d]->SetBranchAddress("LepBDT",&LepBDT);
    inputTree[d]->SetBranchAddress("LepMissingHit",&LepMissingHit);
    inputTree[d]->SetBranchAddress("LepChargedHadIso",&LepChargedHadIso);
    inputTree[d]->SetBranchAddress("LepNeutralHadIso",&LepNeutralHadIso);
    inputTree[d]->SetBranchAddress("LepPhotonIso",&LepPhotonIso);
    inputTree[d]->SetBranchAddress("LepCombRelIsoPF",&LepCombRelIsoPF);
    inputTree[d]->SetBranchAddress("LepCombRelIsoPFPreFSR",&LepCombRelIsoPFPreFSR);
    inputTree[d]->SetBranchAddress("fsrPt",&fsrPt);
    inputTree[d]->SetBranchAddress("fsrEta",&fsrEta);
    inputTree[d]->SetBranchAddress("fsrPhi",&fsrPhi);
    inputTree[d]->SetBranchAddress("fsrLept",&fsrLept);
    inputTree[d]->SetBranchAddress("fsrDR",&fsrDR);
    inputTree[d]->SetBranchAddress("fsrLeptId",&fsrLeptID);


    //---------- Process tree

    Long64_t entries = inputTree[d]->GetEntries();
    cout<<"Processing dataset "<<datasets[d]<<" ("<<entries<<" entries) ..."<<endl;

    for (Long64_t z=0; z<entries; ++z){

      if(DEBUG && z>1000) break;

      inputTree[d]->GetEntry(z);
      
      if( !(Zsel==30 || Zsel==50) ) continue;

      if(LepPt->size()!=2){
       	cout<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", LepPt.size() = "<<LepPt->size()<<endl; // should not happen anymore
       	continue;
      }

      Short_t l1 = (LepPt->at(0)>LepPt->at(1))?0:1; // leading lepton
      Short_t l2 = (LepPt->at(0)>LepPt->at(1))?1:0; // trailing lepton

      if(LepPt->at(l1)<20.) continue;
      if(LepPt->at(l2)<10.) continue;

      Double_t eventWeight = partialSampleWeight[d] * xsec * overallEventWeight ;

      //----- find final state

      currentFinalState = -1;
      if(ZFlav==-121)
	currentFinalState = fs2e;
      else if(ZFlav==-169)
	currentFinalState = fs2mu;
      else
	cerr<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", ZFlav="<<ZFlav<<endl;


      //----- fill histograms

      Float_t varVal[nVariables] = {
	(Float_t)Nvtx,
	ZMass,
	LepPt->at(l1),
	LepPt->at(l2),
	LepEta->at(l1),
	LepEta->at(l2),
	LepPhi->at(l1),
	LepPhi->at(l2),
	LepSIP->at(l1),
	LepSIP->at(l2),
	LepBDT->at(l1),
	LepBDT->at(l2),
	LepCombRelIsoPF->at(l1),
	LepCombRelIsoPF->at(l2),
      };
      Float_t varPairVal[nVarPairs][2] = {
	{ LepPt->at(l1), LepEta->at(l1) },
      };

      for(int v=0; v<nVariables; v++)
	h1[v][currentFinalState][currentProcess]->Fill(varVal[v],(currentProcess==Data)?1.:eventWeight);
      for(int v2=0; v2<nVarPairs; v2++)
	h2[v2][currentFinalState][currentProcess]->Fill(varPairVal[v2][0],varPairVal[v2][1],(currentProcess==Data)?1.:eventWeight);

    } // end for entries

  } // end for datasets


  //---------- Fill 'inclusive' counters
  for(int pr=0; pr<nProcesses; pr++){
    for(int fs=0; fs<nFinalStates; fs++){
      for(int v=0; v<nVariables; v++)
	h1[v][nFinalStates][pr]->Add(h1[v][fs][pr]);
      for(int v2=0; v2<nVarPairs; v2++)
	h2[v2][nFinalStates][pr]->Add(h2[v2][fs][pr]);
    }
  }
  

  //---------- Write histograms to a ROOT file
  TFile* fOutHistos = new TFile(Form("histos_plotZDataVsMC_%.3ffb-1.root",lumi),"recreate");
  fOutHistos->cd();
  for(int fs=0; fs<nFinalStates+1; fs++){
    for(int pr=0; pr<nProcesses; pr++){
      for(int v=0; v<nVariables; v++){
	h1[v][fs][pr]->Write(h1[v][fs][pr]->GetName());
	delete h1[v][fs][pr];
      }
      for(int v2=0; v2<nVarPairs; v2++){
	h2[v2][fs][pr]->Write(h2[v2][fs][pr]->GetName());
	delete h2[v2][fs][pr];
      }
    }
  }
  
  fOutHistos->Close();
  delete fOutHistos;

}





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void DrawDataMC(TCanvas* c, TH1F** h, int v, string lumiText, Bool_t logX = false, Bool_t logY = false) {

  //----- prepare canvas
  TStyle* style = (TStyle*)gROOT->GetStyle("tdrStyle")->Clone();
  style->cd();
  c->cd();
  if(logX) c->SetLogx();
  if(logY) c->SetLogy();

  //----- prepare MC histograms
  TH1F* hStacks[nProcesses];
  bool first = true;
  int previous = -1;
  for(int pr=nProcesses-1; pr>=1; pr--){
    if(useProcess[pr]){
      hStacks[pr] = (TH1F*)h[pr]->Clone();
      if(!first) hStacks[pr]->Add(hStacks[previous]);
      hStacks[pr]->SetFillColor(processFillColor[pr]);
      hStacks[pr]->SetLineColor(DRAWLINES?processLineColor[pr]:processFillColor[pr]);
      if(!DRAWLINES) hStacks[pr]->SetLineColorAlpha(processFillColor[pr],0.);
      hStacks[pr]->GetXaxis()->SetTitleOffset(1.2);
      hStacks[pr]->GetYaxis()->SetTitleOffset(1.4);
      hStacks[pr]->GetXaxis()->SetLabelFont(42);
      hStacks[pr]->GetYaxis()->SetLabelFont(42);
      hStacks[pr]->GetXaxis()->SetTitleFont(42);
      hStacks[pr]->GetYaxis()->SetTitleFont(42);
      if(logX) hStacks[pr]->GetXaxis()->SetMoreLogLabels();
      if(logX) hStacks[pr]->GetXaxis()->SetNoExponent();
      first = false;
      previous = pr;
    }
  }
  int idxSumMC = previous;

  //----- prepare data graph
  h[0]->SetMarkerStyle(20);
  h[0]->SetMarkerColor(kBlack);
  h[0]->SetMarkerSize(1.);
  TGraphAsymmErrors* gData = getDataGraph(h[0]);
  gData->SetMarkerSize(0.5);

  //----- adjust Y axis
  const int npoints = gData->GetN();
  Float_t gDataErrorBarUp[npoints];
  for(int i=0; i<npoints; i++) gDataErrorBarUp[i] = gData->GetY()[i] + gData->GetEYhigh()[i] ;
  Float_t cmax = TMath::Max( (Float_t)hStacks[idxSumMC]->GetMaximum(), (Float_t)TMath::MaxElement(npoints,gDataErrorBarUp) );
  cmax *= logY ? 30. : 1.5 ;
  hStacks[idxSumMC]->SetMaximum(cmax);
  Float_t cminlog = hStacks[idxSumMC]->GetMaximum() / 20000000.;// varMinFactor[v];
  if(logY) hStacks[idxSumMC]->SetMinimum(0.3);//(cminlog);

  //----- prepare legend
  int legPos = varLegPos[v];
  float legLeft = (legPos==11) ? 0.2 : 0.75 ;
  float legWidth = 0.2;
  float legUp = 0.94;
  if(!DRAWDATAMCRATIO) legUp -= 0.05;
  if(!DRAWDATAMCRATIO && legPos==varCMSPos[v]) legUp -= 0.12;
  float legHeight = 0.15;
  if(DRAWDATAMCRATIO) legHeight /= 0.72;
  TLegend* lgd = new TLegend(legLeft,legUp-legHeight,legLeft+legWidth,legUp);
  lgd->SetFillStyle(0);
  lgd->SetBorderSize(0);
  lgd->AddEntry(h[0],processLabel[0].c_str(),"p");
  for(int pr=1; pr<nProcesses; pr++)
    if(useProcess[pr])
      lgd->AddEntry(hStacks[pr],processLabel[pr].c_str(),"f");

  //----- draw everything
  if(!DRAWDATAMCRATIO){

    //--- no Data/MC graph -> draw on main pad
    for(int pr=1; pr<nProcesses; pr++){
      if(useProcess[pr]){
	hStacks[pr]->Draw( pr==idxSumMC ? "HIST" : "HIST SAME" );
	if(!logY) hStacks[pr]->GetYaxis()->SetLabelSize(0.025);
      }
    }
    gData->Draw("P");
    lgd->Draw();
    gPad->RedrawAxis();

  }else{

    //--- with Data/MC graph -> need 2 pads
    TPad* pad1 = new TPad("pad1", "pad1", 0., 0.28, 1., 0.95);
    TPad* pad2 = new TPad("pad2", "pad2", 0., 0.13, 1., 0.28);
    if(logX){pad1->SetLogx(); pad2->SetLogx();}
    if(logY) pad1->SetLogy();
    pad1->SetMargin(0.16,0.02,0.03,0.);
    pad2->SetMargin(0.16,0.02,0.,0.);

    //--- dummy histogram to get the X axis right
    TH1F* hBlank = (TH1F*)hStacks[idxSumMC]->Clone();
    hBlank->Reset();
    hBlank->Draw();
    hBlank->GetYaxis()->SetLabelSize(0.);

    pad1->Draw();
    pad2->Draw();

    //--- main pad
    pad1->cd();
    for(int pr=1; pr<nProcesses; pr++){
      if(useProcess[pr]){
	hStacks[pr]->GetYaxis()->SetTitleOffset(1.);
	hStacks[pr]->GetYaxis()->SetTitleSize(0.06);
	hStacks[pr]->GetYaxis()->SetLabelSize(0.06);
	if(!logY) hStacks[pr]->GetYaxis()->SetLabelSize(0.04);
	hStacks[pr]->Draw( pr==idxSumMC ? "HIST" : "HIST SAME" );
      }
    }
    gData->Draw("P");
    hStacks[idxSumMC]->GetXaxis()->SetLabelSize(0.);
    lgd->Draw();
    pad1->RedrawAxis();
    
    //--- Data/MC pad
    pad2->cd();
    pad2->SetGridy();
    TH1F* hOne = (TH1F*)h[0]->Clone();
    hOne->Reset();
    for(int i=1; i<=hOne->GetNbinsX(); i++) hOne->SetBinContent(i,1.);
    hOne->SetLineColor(kGray);
    hOne->Draw();
    hOne->GetYaxis()->SetRangeUser(0.45,1.55);
    hOne->GetYaxis()->SetNdivisions(206);
    hOne->GetYaxis()->SetTitle("Data/MC");
    hOne->GetYaxis()->SetTitleOffset(0.3);
    hOne->GetYaxis()->SetTitleSize(0.2);
    hOne->GetYaxis()->SetLabelSize(0.17);
    hOne->GetYaxis()->SetTitleFont(42);
    hOne->GetYaxis()->SetLabelFont(42);
    hOne->GetYaxis()->CenterTitle();
    hOne->GetXaxis()->SetTickSize(0.1);
    hOne->GetYaxis()->SetTickSize(0.02);
    TGraphAsymmErrors* gRatio = getDataOverMCGraph(gData,hStacks[idxSumMC]);
    gRatio->Draw("P");
    pad2->RedrawAxis();

  }

  //----- print official CMS labels and lumi
  writeExtraText = true;
  extraText  = "Preliminary";
  lumi_sqrtS = lumiText + " (13 TeV)";
  //CMS_lumi( c, 0, varCMSPos[v] );
  CMS_lumi( c, 0, DRAWDATAMCRATIO ? 0 : varCMSPos[v] );

}


void DrawMC2D(TCanvas* c, TH2F** h2, int v2, string lumiText, Bool_t logX = false, Bool_t logY = false) {

  //----- prepare canvas
  TStyle* style = (TStyle*)gROOT->GetStyle("tdrStyle")->Clone();
  style->SetPadTopMargin(0.05);
  style->SetPadBottomMargin(0.19);
  style->SetPadLeftMargin(0.12);
  style->SetPadRightMargin(0.12);
  style->cd();
  setColZGradient_2();
  c->cd();
  c->UseCurrentStyle();
  if(logX) c->SetLogx();
  if(logY) c->SetLogy();

  //----- prepare MC histogram
  TH2F* h2Stacked = 0;
  bool first = true;
  for(int pr=nProcesses-1; pr>=1; pr--){
    if(useProcess[pr]){
      if(first) h2Stacked = (TH2F*)h2[pr]->Clone();
      else h2Stacked->Add(h2[pr]);
      first = false;
    }
  }
  h2Stacked->GetXaxis()->SetTitleOffset(1.2);
  h2Stacked->GetYaxis()->SetTitleOffset(1.4);
  h2Stacked->GetXaxis()->SetLabelFont(42);
  h2Stacked->GetYaxis()->SetLabelFont(42);
  h2Stacked->GetXaxis()->SetTitleFont(42);
  h2Stacked->GetYaxis()->SetTitleFont(42);
  if(logX) h2Stacked->GetXaxis()->SetMoreLogLabels();
  if(logX) h2Stacked->GetXaxis()->SetNoExponent();
  //h2Stacked->SetMinimum(-1e-20); // avoid white bins

  //----- draw everything
  h2Stacked->Draw("COLZ");
  gPad->RedrawAxis();

  //----- adjust color axis
  c->Update();
  TPaletteAxis* pal = (TPaletteAxis*)h2Stacked->GetListOfFunctions()->FindObject("palette");
  pal->SetX1NDC(0.9);
  pal->SetX2NDC(0.93);
  h2Stacked->GetZaxis()->SetLabelSize(0.035);
  h2Stacked->GetZaxis()->SetLabelFont(42);

  //----- print official CMS labels and lumi
  writeExtraText = true;
  extraText  = "Preliminary";
  lumi_sqrtS = lumiText + " (13 TeV)";
  CMS_lumi( c, 0, 0 );

}


void doPlots(string outputDirectory, int variableList, int varPairList, double lumi, string lumiText)
{

  setTDRStyle();
  gStyle->SetGridColor(kGray);
  //TGaxis::SetMaxDigits(3);

  //---------- retrieve histograms from the ROOT file
  TH1F* h1[nVariables][nFinalStates+1][nProcesses];
  TH2F* h2[nVarPairs ][nFinalStates+1][nProcesses];
  TFile* fInHistos = TFile::Open(Form("histos_plotZDataVsMC_%.3ffb-1.root",lumi));
  for(int fs=0; fs<nFinalStates+1; fs++){
    for(int pr=0; pr<nProcesses; pr++){
      for(int v=0; v<nVariables; v++){
	h1[v][fs][pr] = (TH1F*)fInHistos->Get( Form("h1_%s_%s_%s",varName[v].c_str(),sFinalState[fs].c_str(),sProcess[pr].c_str()) );
	if(fs==fs2e ){
	  h1[v][fs][pr]->GetXaxis()->SetTitle(ReplaceString(h1[v][fs][pr]->GetXaxis()->GetTitle(),"lepton","electron").c_str());
	  h1[v][fs][pr]->GetXaxis()->SetTitle(ReplaceString(h1[v][fs][pr]->GetXaxis()->GetTitle(),"l+l-","e^{+}e^{-}").c_str());
	}
	if(fs==fs2mu){
	  h1[v][fs][pr]->GetXaxis()->SetTitle(ReplaceString(h1[v][fs][pr]->GetXaxis()->GetTitle(),"lepton","muon").c_str());
	  h1[v][fs][pr]->GetXaxis()->SetTitle(ReplaceString(h1[v][fs][pr]->GetXaxis()->GetTitle(),"l+l-","#mu^{+}#mu^{-}").c_str());
	}
	if(restrictCountVar[v]) restrictXAxis(h1[v][fs][pr],restrictCountVar[v]);
      }
      for(int v2=0; v2<nVarPairs; v2++){
	h2[v2][fs][pr] = (TH2F*)fInHistos->Get( Form("h2_%s_%s_%s",varPairName[v2].c_str(),sFinalState[fs].c_str(),sProcess[pr].c_str()) );
      }
    }
  }

  //---------- do the plots (1 canvas per variable)

  string canvasName;
  TCanvas* c1[nVariables][nFinalStates+1];
  TCanvas* c2[nVarPairs ][nFinalStates+1];
  for(int fs=0; fs<nFinalStates; fs++){
    for(int v=0; v<nVariables; v++){
      if(!plotThisVar[variableList][v]) continue;
      canvasName = string(Form("c1_%s_%s",varName[v].c_str(),sFinalState[fs].c_str()));
      c1[v][fs] = new TCanvas(canvasName.c_str(),canvasName.c_str(),500,500);
      DrawDataMC(c1[v][fs],h1[v][fs],v,lumiText,varLogx[v],varLogy[v]);
      SaveCanvas(outputDirectory,c1[v][fs]);
    }
    for(int v2=0; v2<nVarPairs; v2++){
      if(!plotThisVarPair[varPairList][v2]) continue;
      canvasName = string(Form("c2_%s_%s",varPairName[v2].c_str(),sFinalState[fs].c_str()));
      c2[v2][fs] = new TCanvas(canvasName.c_str(),canvasName.c_str(),500,500);
      DrawMC2D(c2[v2][fs],h2[v2][fs],v2,lumiText,varPairLogx[v2],varPairLogy[v2]);
      SaveCanvas(outputDirectory,c2[v2][fs]);
    }
  }

}





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void plotDataVsMC_LeptonVars(bool redoHistograms = true) {

  // --------------- inputs ---------------

  // Define input/output location
  string inputPathMC   = "";
  string inputPathData = "";
  string outputPath = "PlotsDataVsMC_LeptonVars/";

  // Define the luminosity

  // all 50ns (2015B + 2015C_50ns) : 71.52/pb
  // 25ns 2015D (no v4) : 578.3/pb
  // 25ns 2015D + 2015Dv4 (Oct. 17th JSON) : 832.31/pb
  // all 25ns (2015C + 2015D + 2015Dv4) (Nov. 13th Silver JSON) : 2.46/fb

  float lumi = 2.6;
  string lumiText = "2.6 fb^{-1}";

  // Choose a list of 1D plots
  int variableList = 1;

  // Choose a list of 2D plots
  int varPairList = 1;


  // --------------- processing ---------------

  gSystem->Exec(("mkdir -p "+outputPath).c_str());

  // Prepare the histograms and store them in a ROOT file (to be done only once, it can take a few minutes)
  if(redoHistograms)
    doHistograms(inputPathMC, inputPathData, lumi);

  // Do the plots (+- instantaneous)
  doPlots(outputPath, variableList, varPairList, lumi, lumiText);

}

