#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TRandom3.h"
#include "Math/DistFunc.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TF1.h"
#include "TLatex.h"
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
#include "TEfficiency.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include "tdrstyle.C"
// #include "CMS_lumi.C"
#include "plotUtils.C"
#include "ZZAnalysis/AnalysisStep/src/kFactors.C"

#include <ZZAnalysis/AnalysisStep/src/Category.cc>
#include <ZZAnalysis/AnalysisStep/src/bitops.cc>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>
#include <ZZAnalysis/AnalysisStep/test/Plotter/fit_functions.C>

using namespace std;

int useHTBinned = 2;         // 0 - use simple DY inclusive
                             // 1 - use ht binned
                             // 2 - use jet binned + b-enricchement

bool enforceNarrowWidth = true;
bool unblind = true;

int onlyOneLep = 1;          // 0 - ee
                             // 1 - all leptons
                             // 2 - mumu

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

const int nVariables = 27;
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
  "JetQGProduct",
  "ZjetMELAspin0",
  "vbfMELA",
  "ZZMasshighMELA",
  "ZjetMELAspin2",
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
  "Quark gluon likelihood product",
  "ZjetMELAspin0",
  "vbfMELA",
  "m_{2#font[12]{l}2q} (GeV) for ZjetMELAspin0 > 0.5",
  "ZjetMELAspin2",
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
  "Events / 25 GeV",
  "Events / 0.025",
};
Int_t  varNbin[nVariables] = { 70, 50, 70,  56,  44, 50,50, 400,  50,  50,  50,  50,  50,  50,  50, 50, 50, 25, 25, 4, 50, 78, 50, 40, 82, 50, 40};
Float_t varMin[nVariables] = {  250,  0,  250,  40,  40,  90, 90, -200,  0, 0, -0.2, -0.2, 0,  0, -0.2, -1.2, -1.2, 0., 0., -0.5, 0., -0.05,-0.2,0.,-1.05, 250, 0.};
Float_t varMax[nVariables] = { 2000, 500, 2000, 180, 150, 800, 800, 0, 500, 500, 1.2, 1.2, 500, 500, 1.2, 1.2, 1.2 , 3.15, 3.15, 3.5, 300., 1.05, 1.2, 1.,1., 2000.,1.};

Bool_t varLogx[nVariables] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
Bool_t varLogy[nVariables] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0};

enum Process {Data=0, DYjets=1, TTBar=2}; 
const int nProcesses = 3;
string sProcess[nProcesses] = {"Data", "DY", "TT"}; 
string processLabel[nProcesses] = {"Data e#mu", "Z + jets", "t#bar{t},WW"};

// WITH NNLO k-FACTOR FOR Z+JET
// 1.238*0.95
Float_t scaleF[nProcesses] = {1.,1.1761,1.};

const int nFS = 3;
string sFS[nFS] = {"ee","all","mm"};

const int nType = 12;
string typeS[nType] = {"resolvedSB","mergedSB","mergedSR","resolvedSR","resolvedSBbtag","mergedSBbtag","mergedSRbtag","resolvedSRbtag","resolvedSBvbf","mergedSBvbf","mergedSRvbf","resolvedSRvbf"};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

bool passAdditionalCuts(float tmvaZ2Mass, float tmvaZ1Pt, float tmvaZ2Pt, float tmvaZ1tau21,  bool merged) {
  if (tmvaZ2Mass < 60. || tmvaZ2Mass>120.) return false;
  if (tmvaZ1Pt < 100. || tmvaZ2Pt < 100.) return false;   // TEST!
  if (merged && tmvaZ1tau21 > 0.6) return false;
  //if (zzmass < 400.) return false;
 
  return true;
}

TGraphAsymmErrors* doBkgEstGraph(int n, TH1F* central, TH1F* up, TH1F* down, bool isRatio) {
  float theX[n];
  float theY[n];
  float theEX[n];
  float theEYup[n];
  float theEYdown[n];
  for(int ibin=1; ibin<=central->GetNbinsX(); ibin++){
    theX[ibin-1] = central->GetXaxis()->GetBinCenter(ibin);
    theY[ibin-1] = central->GetBinContent(ibin);
    if (isRatio) theY[ibin-1] = 0.;
    theEX[ibin-1] = central->GetXaxis()->GetBinCenter(ibin) - central->GetXaxis()->GetBinLowEdge(ibin);
    theEYup[ibin-1] = up->GetBinContent(ibin)-central->GetBinContent(ibin);
    if (isRatio) theEYup[ibin-1] /= central->GetBinContent(ibin);
    theEYdown[ibin-1] = central->GetBinContent(ibin)-down->GetBinContent(ibin);
    if (isRatio) theEYdown[ibin-1] /= central->GetBinContent(ibin);
  }
  TGraphAsymmErrors* tg = new TGraphAsymmErrors(n,theX,theY,theEX,theEX,theEYup,theEYdown); 
  tg->SetLineColor(kBlue+2);
  tg->SetMarkerColor(kBlue+2);
  tg->SetFillColor(kBlue+2);
  tg->SetFillStyle(3004);
  tg->SetLineWidth(3);
  tg->SetMarkerStyle(25);
  return tg;
}
  
float deltaPhi(float phi1, float phi2) 
{
  float delt = phi1-phi2;
  while (delt >= 3.1416) delt -= 6.2832;
  while (delt < -3.1416) delt += 6.2832;
  return delt;
}

void densityHist(TH1F* hist) 
{
  for (int i=1; i<=hist->GetNbinsX(); i++) {
    float binsize = hist->GetXaxis()->GetBinUpEdge(i) - hist->GetXaxis()->GetBinLowEdge(i);
    hist->SetBinContent(i,hist->GetBinContent(i)*50./binsize);
    hist->SetBinError(i,hist->GetBinError(i)*50./binsize);
  }
  hist->GetXaxis()->SetTitle("m_{ZZ} (GeV)");  
  hist->GetYaxis()->SetTitle("Events / 50 GeV");
}

float getDVBF2jetsConstant(float ZZMass){
  float par[9]={
    1.876,
    -55.488,
    403.32,
    0.3906,
    80.8,
    27.7,
    -0.06,
    54.97,
    309.96
  };
  float kappa =
    pow(1.-atan((ZZMass-par[1])/par[2])*2./TMath::Pi(), par[0])
    + par[3]*exp(-pow((ZZMass-par[4])/par[5], 2))
    + par[6]*exp(-pow((ZZMass-par[7])/par[8], 2));
  float constant = kappa/(1.-kappa);
  return constant;
}

float getDZjjspin0Constant(float ZZMass){
  float constant = 0.035*(3.05+(0.0005*ZZMass)-(2.73/(1.+exp((ZZMass-2500.)/258.26)))) ;
  return constant;
}

float getDZjjspin2Constant(float ZZMass){
  return 0.14;
}

bool debug_ = false;

void plotDataVsMC_2l2qTemplate_emu(string dirout = "ReminiAOD2016_mZZ500_EMu", string theNtupleFile = "./goodDatasets_Moriond2017_Reminiaod_Emu.txt", bool sync = false, bool CR = false, int draw = 1, bool isRefitSB = false, bool useDyNLO = false)
{ 

  // draw = 0 : do not draw, but save inputs for fits
  //      = 1 : draw only PAS-style plots
  //      = 2 : draw all
  // CR = false : SB/SR events
  //      true  : Z+1jet CR to check q/g tagger performance

  useDyNLO = false;
  draw = 2;

  if(useDyNLO) dirout+="_nloDY"; 
  cout<<"dirout "<<dirout<<endl;

  const float lumin = 35.8;   // ICHEP total
  setTDRStyle();

  const int nDatasets = 17;
  const int nDatasetsMC = 8;  

  // cut on mZZ before plotting
  float minMZZresolved = 200.0;
  float minMZZmerged = 200.0;

  if(draw>0){
       minMZZresolved = 500.0;
       minMZZmerged = 500.0;
  }

  // range of mZZ that stores in the output rootfile
  const float minMZZ = 500;
  const float maxMZZ = 3500;

  string mZZtype = "";
  if(isRefitSB)  mZZtype="Refit";
  else mZZtype="Reco";

  TFile* inputFile[nDatasets];
  TChain* inputTree[nDatasets];
  TH1F* hCounters[nDatasets]; 
  Long64_t NGenEvt[nDatasets];
  Float_t NEvtNarrow[nDatasets];
  Float_t sumWeights[nDatasets];
  Int_t mass[nDatasets];

  string Dsname[nDatasets] = {"DY1Jet","DY2Jet","DY3Jet","DY4Jet","DYBFiltJets","DYBJets","TTBar","WW","MuEG2016Bv1","MuEG2016Bv2","MuEG2016C","MuEG2016D","MuEG2016E","MuEG2016F","MuEG2016G","MuEG2016Hv2","MuEG2016Hv3"};

  // ZPT reweighting for LO MC
  TFile* fKzpt = new TFile("ZPTCorrection/Kfac_pTZ.root");
  TH1D* NLO = (TH1D*)fKzpt->Get("NLO");// This NLO is 2015 measured results
  TH1D* LO = (TH1D*)fKzpt->Get("LO");// This LO is the LO*1.238 where 1.238 in inclusive NLO/LO K factor
  TH1D* Kzpt_ratio = (TH1D*)NLO->Clone("Kzpt_ratio");
  Kzpt_ratio->Divide(LO);
  TGraphErrors* gKzpt_ratio = new TGraphErrors(Kzpt_ratio);

  // ZPT reweighting for NLO DY
  TFile* fdyzpt = new TFile("ZPTCorrection/UnfoldingOutputZPt.root");
  TH1D* hdyzptdt = (TH1D*)fdyzpt->Get("hUnfold");
  TH1D* hdyzptmc = (TH1D*)fdyzpt->Get("hTruth");
  TH1D* hdyzpt_dtmc_ratio = (TH1D*)hdyzptdt->Clone("hdyzpt_dtmc_ratio");
  hdyzpt_dtmc_ratio->Divide(hdyzptmc);
  TGraphErrors* gdyzpt_dtmc_ratio = new TGraphErrors(hdyzpt_dtmc_ratio);

  // EWK NOT AVAILABLE
  TFile* fewk = new TFile("EWKCorrection/kfactor.root");
  TH1D* hEWK = (TH1D*)fewk->Get("kfactor_ewk");
  TGraphErrors* gKewk = new TGraphErrors(hEWK); 

  if (useHTBinned == 0) {
    Dsname[2] = "DYJetsToLL";
  } else if (useHTBinned == 1) {
    Dsname[2] = "DYHT100";
    Dsname[3] = "DYHT200";
    Dsname[4] = "DYHT400";
    Dsname[5] = "DYHT600";
    processLabel[3] = "Z + jets (HT > 100 GeV)"; 
  }

  // trigger weights
  TFile *f_data_se  = new TFile("trigeff/eff_data_2d_se.root","READ");
  TFile *f_data_de1 = new TFile("trigeff/eff_data_2d_de1.root","READ");
  TFile *f_data_de2 = new TFile("trigeff/eff_data_2d_de2.root","READ");
  TFile *f_data_sm  = new TFile("trigeff/eff_data_2d_sm.root","READ");
  TFile *f_data_dm1 = new TFile("trigeff/eff_data_2d_dm1.root","READ");
  TFile *f_data_dm2 = new TFile("trigeff/eff_data_2d_dm2.root","READ");
  
  TCanvas *c_data_se = (TCanvas*)f_data_se->Get("c1");
  TEfficiency *eff_data_se = (TEfficiency*)c_data_se->GetPrimitive("den_se_2d_clone");
  TCanvas *c_data_de1 = (TCanvas*)f_data_de1->Get("c1");
  TEfficiency *eff_data_de1 = (TEfficiency*)c_data_de1->GetPrimitive("den_de1_2d_clone");
  TCanvas *c_data_de2 = (TCanvas*)f_data_de2->Get("c1");
  TEfficiency *eff_data_de2 = (TEfficiency*)c_data_de2->GetPrimitive("den_de2_2d_clone");
  TCanvas *c_data_sm = (TCanvas*)f_data_sm->Get("c1");
  TEfficiency *eff_data_sm = (TEfficiency*)c_data_sm->GetPrimitive("den_sm_2d_clone");
  TCanvas *c_data_dm1 = (TCanvas*)f_data_dm1->Get("c1");
  TEfficiency *eff_data_dm1 = (TEfficiency*)c_data_dm1->GetPrimitive("den_dm1_2d_clone");
  TCanvas *c_data_dm2 = (TCanvas*)f_data_dm2->Get("c1");
  TEfficiency *eff_data_dm2 = (TEfficiency*)c_data_dm2->GetPrimitive("den_dm2_2d_clone");

  TFile *f_mc_se  = new TFile("trigeff/eff_mc_2d_se.root","READ");
  TFile *f_mc_de1 = new TFile("trigeff/eff_mc_2d_de1.root","READ");
  TFile *f_mc_de2 = new TFile("trigeff/eff_mc_2d_de2.root","READ");
  TFile *f_mc_sm  = new TFile("trigeff/eff_mc_2d_sm.root","READ");
  TFile *f_mc_dm1 = new TFile("trigeff/eff_mc_2d_dm1.root","READ");
  TFile *f_mc_dm2 = new TFile("trigeff/eff_mc_2d_dm2.root","READ");
  
  TCanvas *c_mc_se = (TCanvas*)f_mc_se->Get("c2");
  TEfficiency *eff_mc_se = (TEfficiency*)c_mc_se->GetPrimitive("den_se_2d_clone");
  TCanvas *c_mc_de1 = (TCanvas*)f_mc_de1->Get("c2");
  TEfficiency *eff_mc_de1 = (TEfficiency*)c_mc_de1->GetPrimitive("den_de1_2d_clone");
  TCanvas *c_mc_de2 = (TCanvas*)f_mc_de2->Get("c2");
  TEfficiency *eff_mc_de2 = (TEfficiency*)c_mc_de2->GetPrimitive("den_de2_2d_clone");
  TCanvas *c_mc_sm = (TCanvas*)f_mc_sm->Get("c2");
  TEfficiency *eff_mc_sm = (TEfficiency*)c_mc_sm->GetPrimitive("den_sm_2d_clone");
  TCanvas *c_mc_dm1 = (TCanvas*)f_mc_dm1->Get("c2");
  TEfficiency *eff_mc_dm1 = (TEfficiency*)c_mc_dm1->GetPrimitive("den_dm1_2d_clone");
  TCanvas *c_mc_dm2 = (TCanvas*)f_mc_dm2->Get("c2");
  TEfficiency *eff_mc_dm2 = (TEfficiency*)c_mc_dm2->GetPrimitive("den_dm2_2d_clone");
  // end trigger weights

  float tau21bin[25] = {
    -0.00769231,
    0.0346154,
    0.0769231,
    0.119231,
    0.161538,
    0.203846,
    0.246154,
    0.288462,
    0.330769,
    0.373077,
    0.415385,
    0.457692,
    0.5,
    0.542308,
    0.584615,
    0.626923,
    0.669231,
    0.711538,
    0.753846,
    0.796154,
    0.838461,
    0.880769,
    0.923077,
    0.965385,
    1.00769} ;

  float tau21corr[24] = {0,
    0      ,
    0,
    0.290173  ,
    -0.377686 ,
    0.0977722  ,
    0.412889   ,
    0.322422   ,
    -0.169658  ,
    -0.272581  ,
    -0.384607  ,
    0.109909   ,
    -0.123669  ,
    -0.158826  ,
    0.0764657  ,
    0.155703   ,
    -0.0750966 ,
    -0.0484963 ,
    0.239322   ,
    0.107814   ,
    -1.40918   ,
    0      ,
    0     ,
    0      };

  /// I/O to TMVA 
  float tmvaZZPt, tmvaZ2Mass, tmvaZ1Pt, tmvaZ1tau21, tmvaZ2Pt, tmvaLepPt1, tmvaLepPt2, tmvaJetPt1, tmvaJetPt2;
  float tmvaJetQGLikelihood1, tmvaJetQGLikelihood2, tmvaabshelcosthetaZ1, tmvahelcosthetaZ2, tmvacosthetastar, tmvahelphi, tmvaphistarZ1; 

  TString postfix = "njetLoDY";
  if(useDyNLO) postfix = "nloDY";

  TString outfileName( "TMVAAndRoofitInputs_MuEG_"+postfix+".root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  ofstream mengch;
  mengch.open("mengch.txt");

  char filestring[400];
  float xMin[nType],xMax[nType];

  Int_t RunNumber;
  Long64_t EventNumber;
  Int_t LumiNumber;
  Float_t genEventWeight;
  Float_t overallEventWeight;
  Float_t Nvtx;
  Float_t genHMass;
  vector<Short_t> *ZZsel = 0;
  vector<Float_t> *ZZMass = 0;
  vector<Float_t> *ZZMass_up = 0;
  vector<Float_t> *ZZMass_dn = 0;
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
  vector<Float_t> *LepEta = 0;
  vector<Float_t> *LepPhi = 0;
  Short_t nJets;
  Short_t nJetsBTagged;
  vector<Float_t> *JetPt = 0;
  vector<Float_t> *JetPhi = 0;
  vector<Float_t> *JetEta = 0;
  vector<bool> *JetIsInZZCand = 0;
  vector<Float_t> *JetBTagger = 0;
  vector<Float_t> *JetIsBtaggedWithSF = 0;
  vector<Float_t> *JetQGLikelihood = 0;
  vector<Float_t> *helcosthetaZ1 = 0;  
  vector<Float_t> *helcosthetaZ2 = 0;  
  vector<Float_t> *costhetastar = 0;	  
  vector<Float_t> *helphi = 0;	  
  vector<Float_t> *phistarZ1 = 0;
  vector<Float_t> *p0plus_VAJHU = 0; 
  vector<Float_t> *p2bplus_VAJHU = 0; 
  vector<Float_t> *pqqZJJ_VAMCFM = 0;  
  vector<Float_t> *pvbf_VAJHU_highestPTJets = 0; 
  vector<Float_t> *phjj_VAJHU_highestPTJets = 0; 
  Float_t xsec;
  Float_t Met;
  Float_t GenLep1Pt;
  Float_t GenLep1Phi;
  Float_t GenLep2Pt;
  Float_t GenLep2Phi;
  Int_t GenLep1Id,GenLep2Id;

  TH1F* h1[nVariables][nProcesses][nFS][nType];
  TH2F* h2_mjjVSdR[nProcesses][nFS][nType];

  for(int rs=0; rs<nFS; rs++){   //ee, mumu, or all

    for(int pr=0; pr<nProcesses; pr++){

      for(int nt=0; nt<nType; nt++){

	for(int v=0; v<nVariables; v++){

          if (nt == 0 || nt == 3)
            h2_mjjVSdR[pr][rs][nt] = new TH2F(Form("h2_mjjVSdR_%s_%s_%s",
                                                    sFS[rs].c_str(),typeS[nt].c_str(),sProcess[pr].c_str()),
                                                   "mjjVSdR", 80, 0, 4, 56, 40, 180);
          else
            h2_mjjVSdR[pr][rs][nt] = new TH2F(Form("h2_mjjVSdR_%s_%s_%s",
                                                   sFS[rs].c_str(),typeS[nt].c_str(),sProcess[pr].c_str()),
                                                   "mjjVSdR", 40, 0, 4, 28, 40, 180);

	  if (nt == 0 || ( nt == 3 && (v!=0 && v!=2)) ) 
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
  
  TH1F* hmass[nProcesses][nType+7];
  TH1F* hmassRefit[nProcesses][nType+4];
  TH1F* hmass_up[nProcesses][nType+4];
  TH1F* hmass_down[nProcesses][nType+4];

  Float_t binsincl[] = { 500, 513, 525, 537, 550, 563, 575, 587, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 825, 850, 875, 900, 950, 1000, 1050, 1100, 1150, 1200, 1300, 1400, 1500, 1600, 2000, 2400, 3000};
  Int_t  binnumincl = sizeof(binsincl)/sizeof(Float_t) - 1; 
  Float_t binsbtag[] = { 500, 550, 600, 650, 700, 750, 800, 900, 1000, 1100, 1200, 1400, 1600, 2000, 2400, 3000};
  Int_t  binnumbtag = sizeof(binsbtag)/sizeof(Float_t) - 1; 
  TString hmasstags[nType+7] = {"","","untagged, merged jet","untagged, resolved jets",
				"","","b-tagged, merged jet","b-tagged, resolved jets",
				"","","VBF-tagged, merged jet","VBF-tagged, resolved jets",
				"","","","","all categories, resolved jets","all categories, resolved jets","all categories, resolved jets"};

  for(int pr=0; pr<nProcesses; pr++){

    for(int nt=0; nt<nType+7; nt++){

      if (nt == 0 || nt == 3) {        // mZZ high stats

	if (draw > 0) { hmass[pr][nt] = new TH1F(Form("hmass_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
					         Form("hmass_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),		
			 		         binnumincl,binsincl);
 
                        hmassRefit[pr][nt] = new TH1F(Form("hmassRefit_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
                                                      Form("hmassRefit_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
                                                      binnumincl,binsincl);

        }
	else{ hmass[pr][nt] = new TH1F(Form("hmass_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
				       Form("hmass_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),		
				       int((maxMZZ-minMZZ)/5), minMZZ, maxMZZ);

              hmassRefit[pr][nt] = new TH1F(Form("hmassRefit_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
                                            Form("hmassRefit_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
                                            int((maxMZZ-minMZZ)/5), minMZZ, maxMZZ);

        } 
	if (draw > 0) hmass_up[pr][nt] = new TH1F(Form("hmass_up_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
						  Form("hmass_up_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),			                   	               binnumincl,binsincl);
	else hmass_up[pr][nt] = new TH1F(Form("hmass_up_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
					 Form("hmass_up_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),		
					 int((maxMZZ-minMZZ)/5), minMZZ, maxMZZ);
	if (draw > 0) hmass_down[pr][nt] = new TH1F(Form("hmass_down_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
						    Form("hmass_down_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),		   					         binnumincl,binsincl);
	else  hmass_down[pr][nt] = new TH1F(Form("hmass_down_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
					    Form("hmass_down_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),	
					    int((maxMZZ-minMZZ)/5), minMZZ, maxMZZ);
      } else if (nt<nType) { // mZZ low stats
	 if (draw > 0)  { hmass[pr][nt] = new TH1F(Form("hmass_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
						   Form("hmass_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),		
						   binnumbtag,binsbtag);
                          hmassRefit[pr][nt] = new TH1F(Form("hmassRefit_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
                                                        Form("hmassRefit_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
                                                        binnumbtag,binsbtag);
         }
	 else{ 
              hmass[pr][nt] = new TH1F(Form("hmass_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
	   			       Form("hmass_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),		
				       int((maxMZZ-minMZZ)/5), minMZZ, maxMZZ);

              hmassRefit[pr][nt] = new TH1F(Form("hmassRefit_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
                                            Form("hmassRefit_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
                                            int((maxMZZ-minMZZ)/5), minMZZ, maxMZZ);

         }
	 if (draw > 0) hmass_up[pr][nt] = new TH1F(Form("hmass_up_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
						   Form("hmass_up_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),								 binnumbtag,binsbtag);
	 else hmass_up[pr][nt] = new TH1F(Form("hmass_up_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
					  Form("hmass_up_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),		
					  int((maxMZZ-minMZZ)/5), minMZZ, maxMZZ);
	 if (draw > 0 ) hmass_down[pr][nt] = new TH1F(Form("hmass_down_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
						      Form("hmass_down_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),								   binnumbtag,binsbtag);
	 else  hmass_down[pr][nt] = new TH1F(Form("hmass_down_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),
					     Form("hmass_down_%s_%s",typeS[nt].c_str(),sProcess[pr].c_str()),	
					     int((maxMZZ-minMZZ)/5), minMZZ, maxMZZ);
      } else if (nt<nType+4) { // mZZ non b-tagged
	hmass[pr][nt] = new TH1F(Form("hmass_%snobtag_%s",typeS[nt-nType].c_str(),sProcess[pr].c_str()),
				 Form("hmass_%snobtag_%s",typeS[nt-nType].c_str(),sProcess[pr].c_str()),		
				 binnumincl,binsincl);
        hmassRefit[pr][nt] = new TH1F(Form("hmassRefit_%snobtag_%s",typeS[nt-nType].c_str(),sProcess[pr].c_str()),
                                 Form("hmassRefit_%snobtag_%s",typeS[nt-nType].c_str(),sProcess[pr].c_str()),
                                 binnumincl,binsincl);

        hmass_up[pr][nt] = new TH1F(Form("hmass_up_%snobtag_%s",typeS[nt-nType].c_str(),sProcess[pr].c_str()),
				    Form("hmass_up_%snobtag_%s",typeS[nt-nType].c_str(),sProcess[pr].c_str()),		
				    binnumincl,binsincl);
	hmass_down[pr][nt] = new TH1F(Form("hmass_down_%snobtag_%s",typeS[nt-nType].c_str(),sProcess[pr].c_str()),
				      Form("hmass_down_%snobtag_%s",typeS[nt-nType].c_str(),sProcess[pr].c_str()),		
				      binnumincl,binsincl);
      } else if (nt == nType+4)          // MELA spin-0
	hmass[pr][nt] = new TH1F(Form("hmelaspin0_resolvedSR_%s",sProcess[pr].c_str()),
				 Form("hmelaspin0_resolvedSR_%s",sProcess[pr].c_str()),		
				 20,0.,1.);
      else if (nt == nType+5)         // MELA spin-2
	hmass[pr][nt] = new TH1F(Form("hmelaspin2_resolvedSR_%s",sProcess[pr].c_str()),
				 Form("hmelaspin2_resolvedSR_%s",sProcess[pr].c_str()),		
				 20,0.,1.);
      else if (nt == nType+6)         // VBF-tagged
	hmass[pr][nt] = new TH1F(Form("hmelavbf_allSR_%s",sProcess[pr].c_str()),
				 Form("hmelavbf_allSR_%s",sProcess[pr].c_str()),		
				 42,-1.05,1.05);
 
      hmass[pr][nt]->Sumw2();
    }
  }
  
  TF1* ffit[nType];
  TF1* ffitup[nType];
  TF1* ffitdown[nType];
  
  TH1F* hbkg[nType];
  TH1F* hbkg_up[nType];
  TH1F* hbkg_down[nType];

  for(int nt=0; nt<nType; nt++){
    hbkg[nt] = (TH1F*)hmass[0][nt]->Clone();
    hbkg_up[nt] = (TH1F*)hmass[0][nt]->Clone();
    hbkg_down[nt] = (TH1F*)hmass[0][nt]->Clone();
  }

  if (unblind) {
    for(int nt=0; nt<nType; nt++){
      ifstream parInput;
      float parValue[20];  string parName;  
      xMin[nt] = minMZZ;
      xMax[nt] = maxMZZ;
      if (nt == 3 || nt == 6 || nt == 10 || nt == 7 || nt == 11) {   // MODIFY HERE
        if (nt == 3) {
	  parInput.open("/afs/cern.ch/user/t/tocheng/public/ZZ2l2q_Moriond/ReminiAOD/resolved_untagged.txt");
	  //cout << "Reading resolved untagged fit file" << endl;

	}
	else if (nt == 7) {
	  parInput.open("/afs/cern.ch/user/t/tocheng/public/ZZ2l2q_Moriond/ReminiAOD/resolved_btagged.txt");
	  //cout << "Reading resolved b-tagged fit file" << endl;
	}
	else if (nt == 11) {
	  parInput.open("/afs/cern.ch/user/t/tocheng/public/ZZ2l2q_Moriond/ReminiAOD/resolved_vbftagged.txt");
	  //cout << "Reading resolved vbf-tagged fit file" << endl;
	}
	else if (nt == 6) {
	  parInput.open("/afs/cern.ch/user/t/tocheng/public/ZZ2l2q_Moriond/ReminiAOD/merged_btagged.txt");
	  //cout << "Reading merged b-tagged fit file" << endl;
	}
	else if (nt == 10) {
	  parInput.open("/afs/cern.ch/user/t/tocheng/public/ZZ2l2q_Moriond/ReminiAOD/merged_vbftagged.txt");
	  //cout << "Reading merged vbf-tagged fit file" << endl;
	}
	for(int iline=0; iline<4; iline++){
	  parInput >> parName >> parValue[iline];
	  //cout << parValue[iline] << endl;
	}
	for(int iline=0; iline<4; iline++){
	  parInput >> parValue[4*(iline+1)] >> parValue[1+4*(iline+1)] >> parValue[2+4*(iline+1)] >> parValue[3+4*(iline+1)];
          //cout << parValue[4*(iline+1)] << " " << parValue[1+4*(iline+1)] << " "<< parValue[2+4*(iline+1)] << " " << parValue[3+4*(iline+1)] << endl;
	}

	if (nt == 3 || nt == 7 || nt == 11) {
	  
	  ffit[nt] = new TF1("pippo",myfunction,xMin[nt],xMax[nt],4);
	  ffitup[nt] = new TF1("pippoup",myfunctionErrUp,xMin[nt],xMax[nt],20);
	  ffitdown[nt] = new TF1("pippodown",myfunctionErrDown,xMin[nt],xMax[nt],20);
	} else {
	  ffit[nt] = new TF1("pippo",GaussExp,xMin[nt],xMax[nt],4);
	  ffitup[nt] = new TF1("pippoup",GaussExpErrUp,xMin[nt],xMax[nt],20);
	  ffitdown[nt] = new TF1("pippodown",GaussExpErrDown,xMin[nt],xMax[nt],20);
	}
	for (int i=0; i<20; i++) {
	  ffitup[nt]->SetParameter(i,parValue[i]);
	  ffitdown[nt]->SetParameter(i,parValue[i]);
	  if (i<4) ffit[nt]->SetParameter(i,parValue[i]);
	}
      }  else if (nt == 2) {
	TFile f1("/afs/cern.ch/user/t/tocheng/public/ZZ2l2q_Moriond/ReminiAOD/Template1D_spin0_merged_all.root");
	TH1F* htemp = (TH1F*)f1.Get("merged_DY"); 
	// rebin before density
	int ipast = 1;  float summa = 0.;  float summaErr = 0.;
	for(int ibin=1; ibin<=htemp->GetNbinsX(); ibin++){
          if (htemp->GetXaxis()->GetBinCenter(ibin) < xMin[nt]) continue;
	  if (htemp->GetXaxis()->GetBinCenter(ibin) > xMax[nt]) break;
	  if (htemp->GetXaxis()->GetBinCenter(ibin) < binsbtag [ipast]) {
	    summa += htemp->GetBinContent(ibin); 
	    summaErr += htemp->GetBinError(ibin)*htemp->GetBinError(ibin);
	    if (htemp->GetXaxis()->GetBinCenter(ibin+1) > binsbtag [ipast]) {
	      hbkg[nt]->SetBinContent(ipast,summa);
	      hbkg_up[nt]->SetBinContent(ipast,summa+sqrt(summaErr));
	      hbkg_down[nt]->SetBinContent(ipast,summa-sqrt(summaErr));
	      ipast++;         summa = 0.;      summaErr = 0.;
	    }
	  }
	}
	densityHist(hbkg[nt]);
	densityHist(hbkg_up[nt]);
	densityHist(hbkg_down[nt]);
      }    
    }

    for(int nt=0; nt<nType; nt++){
      if (nt == 3 || nt == 7 || nt == 11 || nt == 6 || nt == 10) {   // MODIFY HERE
	for(int ibin=1; ibin<=hbkg[nt]->GetNbinsX(); ibin++){
	  hbkg[nt]->SetBinContent(ibin,(ffit[nt]->Integral(hbkg[nt]->GetXaxis()->GetBinLowEdge(ibin),hbkg[nt]->GetXaxis()->GetBinUpEdge(ibin)))*50./hbkg[nt]->GetBinWidth(ibin));
	  hbkg_up[nt]->SetBinContent(ibin,(ffitup[nt]->Integral(hbkg[nt]->GetXaxis()->GetBinLowEdge(ibin),hbkg[nt]->GetXaxis()->GetBinUpEdge(ibin)))*50./hbkg[nt]->GetBinWidth(ibin));
	  hbkg_down[nt]->SetBinContent(ibin,(ffitdown[nt]->Integral(hbkg[nt]->GetXaxis()->GetBinLowEdge(ibin),hbkg[nt]->GetXaxis()->GetBinUpEdge(ibin)))*50./hbkg[nt]->GetBinWidth(ibin));
	}
      }
    }
  }
  // end background functions
  
  //---------- Will loop over all datasets
  for (int d=0; d<nDatasets; d++) {
    NGenEvt[d] = 0;      NEvtNarrow[d] = 0.;
    sumWeights[d] = 0.;  mass[d] = 0;
    inputTree[d] = new TChain("ZZTree/candTree");
  }

  ifstream list(theNtupleFile.c_str());
  char fileName[400];
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
          sumWeights[d] += hCounters[d]->GetBinContent(40);  // all weights
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

    if(debug_) cout<<"for "<<d<<"th dataset "<<Dsname[d]<<endl; 

    inputTree[d]->SetBranchAddress("Nvtx", &Nvtx);
    inputTree[d]->SetBranchAddress("RunNumber", &RunNumber);
    inputTree[d]->SetBranchAddress("EventNumber", &EventNumber);
    inputTree[d]->SetBranchAddress("LumiNumber", &LumiNumber);
    inputTree[d]->SetBranchAddress("genHEPMCweight", &genEventWeight);
    inputTree[d]->SetBranchAddress("overallEventWeight", &overallEventWeight);
    inputTree[d]->SetBranchAddress("ZZsel", &ZZsel);
    inputTree[d]->SetBranchAddress("ZZMass", &ZZMass);
    inputTree[d]->SetBranchAddress("ZZMass_up", &ZZMass_up);
    inputTree[d]->SetBranchAddress("ZZMass_dn", &ZZMass_dn);
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
    inputTree[d]->SetBranchAddress("LepPhi", &LepPhi);
    inputTree[d]->SetBranchAddress("LepEta", &LepEta);
    inputTree[d]->SetBranchAddress("JetPt", &JetPt);
    inputTree[d]->SetBranchAddress("JetEta", &JetEta);
    inputTree[d]->SetBranchAddress("JetPhi", &JetPhi);
    inputTree[d]->SetBranchAddress("JetIsInZZCand", &JetIsInZZCand);
    inputTree[d]->SetBranchAddress("JetBTagger", &JetBTagger);
    inputTree[d]->SetBranchAddress("JetIsBtaggedWithSF", &JetIsBtaggedWithSF);
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
    inputTree[d]->SetBranchAddress("p2bplus_VAJHU", &p2bplus_VAJHU );
    inputTree[d]->SetBranchAddress("pqqZJJ_VAMCFM", &pqqZJJ_VAMCFM );  
    inputTree[d]->SetBranchAddress("pvbf_VAJHU_highestPTJets", &pvbf_VAJHU_highestPTJets );
    inputTree[d]->SetBranchAddress("phjj_VAJHU_highestPTJets", &phjj_VAJHU_highestPTJets );
    inputTree[d]->SetBranchAddress("GenLep1Pt", &GenLep1Pt);
    inputTree[d]->SetBranchAddress("GenLep1Phi", &GenLep1Phi);
    inputTree[d]->SetBranchAddress("GenLep2Pt", &GenLep2Pt);
    inputTree[d]->SetBranchAddress("GenLep2Phi", &GenLep2Phi);
    inputTree[d]->SetBranchAddress("GenLep1Id", &GenLep1Id);
    inputTree[d]->SetBranchAddress("GenLep2Id", &GenLep2Id);


    //---------- Process tree

    Long64_t entries = inputTree[d]->GetEntries();

    // categorize datasets into process
    int process;
    if (d==0 || d==1 || d==2 || d==3 || d==4 || d==5) process=1;                  //DY
    else if (d==6 || d==7) process=2;                     //TT+WW
    else process=0;                                       //Data

    cout<<d<<"th dataet corresponds to process "<<process<<endl;

    // for synchronization
    ofstream myfile;
    int nPassMerged = 0 ;
    int nPassResol = 0 ;
    if (process == 0 && sync) myfile.open("synchronization.txt");
    /*
    if (process>0) {
    float eff = float(entries)/float(NGenEvt[d]);
    cout<<"Processing dataset "<<Dsname[d]<<" ("<<entries<<" entries of "<< NGenEvt[d] << " = " << eff*100. << "%) ..."<<endl;
    cout<<"process "<<process<<" d "<<d<<" "<<Dsname[d]<<endl;

    } else {
      cout<<"Processing dataset "<<d<<" ("<<entries<<" entries)"<<endl;
    }
    */

    if(debug_) cout<<"start looing over tree for"<<Dsname[d]<<endl;
    for (Long64_t z=0; z<entries; ++z){

      inputTree[d]->GetEntry(z);
      if(debug_) cout<<"tree entry "<<z;
      bool syncevents = process == 0 && ((RunNumber==273158 && EventNumber==402295935 && LumiNumber==256) ||
					(RunNumber==274200 && EventNumber==806347507 && LumiNumber==499) || 
					(RunNumber==278167 && EventNumber==1293600593  && LumiNumber==723));

      if(debug_) cout<<"ZZsel->size() "<<ZZsel->size()<<endl;
      if (ZZsel->size() < 1) continue;
      bool writeThis = (process==0 && sync); 

      Double_t eventWeight = 1. ;

      if (process > 0) eventWeight = ( lumin * 1000 * scaleF[process] * overallEventWeight  / sumWeights[d] ) * xsec ;
      if(debug_) cout<<"apply DY Zpt reweighting"<<endl;
      /*
      if(process==1){

        if(string(Dsname[d]).find("DYM50nlo") == std::string::npos) //DY but only for LO jet-binned
        {
           if(useDyNLO ) continue;
           double ZPt = sqrt(GenLep1Pt*GenLep1Pt + GenLep2Pt*GenLep2Pt + 2*GenLep1Pt*GenLep2Pt*cos(GenLep1Phi-GenLep2Phi));
           //cout<<"ZPT reweighting comparing GEN zpt "<<ZPt<<" and RECO zpt "<<double(Z2Pt->at(0))<<endl;
           double Kzpt = gKzpt_ratio->Eval(double(Z2Pt->at(0)));
	   eventWeight = eventWeight * Kzpt;
           //if(string(Dsname[d]).find("DYBFiltJet") == std::string::npos) continue;
        }
        else{ //DY for NLO
           if(!useDyNLO) continue;
           double ZPt = sqrt(GenLep1Pt*GenLep1Pt + GenLep2Pt*GenLep2Pt + 2*GenLep1Pt*GenLep2Pt*cos(GenLep1Phi-GenLep2Phi));
           double Kzpt = gdyzpt_dtmc_ratio->Eval(double(Z2Pt->at(0)));

           eventWeight = eventWeight*Kzpt/scaleF[3];
        }

      }
      */

      if(debug_) cout<<"Apply extra kFactors: LO to NNLO in QCD for TTbar xsec "<<Dsname[d]<<endl;
      if(string(Dsname[d]).find("TT") != std::string::npos) {
          eventWeight = eventWeight * 87.31/57.35;   
      }  
      if(debug_) cout<<"keep only events around nominal mass!!!"<<endl;  
      if (enforceNarrowWidth && mass[d] > 0 && fabs(genHMass-(float)mass[d]) > 0.05*mass[d]) continue;
      NEvtNarrow[d] = NEvtNarrow[d] + 1.;

      if (process > 0 && z == 0) 
         cout <<"sample "<<Dsname[d]<<"with cross-section = " << xsec << " pb; eventweight = " << eventWeight <<" with sumWeights[d]"<<sumWeights[d]<<endl;
      
      if(debug_) cout<<"find leading jets (notice this vector also includes subjets (identified by a btagger value of -1) which must be treated apart"<<endl;
      float pt1stJet = 0.0001;
      float pt2ndJet = 0.0001;
      float eta1stJet = -999;
      float eta2ndJet = -999;
      float phi1stJet = -999;
      float phi2ndJet = -999;

      float btag1stJet = -1.;
      float btag2ndJet = -1.;
      bool isB1stJet = false;
      bool isB2ndJet = false;
      float qglik1stJet = -1.;
      float qglik2ndJet = -1.;

      float pt1stSubjet = 0.0001;
      float pt2ndSubjet = 0.0001;
      float eta1stSubJet = -999;
      float eta2ndSubJet = -999;
      float phi1stSubJet = -999;
      float phi2ndSubJet = -999;

      float btag1stSubjet = -1.;
      float btag2ndSubjet = -1.;
 
      int nInJets = 0;
      int nExtraJets = 0;

      for (unsigned int nJet=0; nJet<JetPt->size(); nJet++) {
	if (JetQGLikelihood->at(nJet) > -800.) {  // ak4 jets
 	  if (JetIsInZZCand->at(nJet) ) {         
	    if (pt1stJet < JetPt->at(nJet)) { // find leading pT jet/subjet
	      pt2ndJet = pt1stJet;
              eta2ndJet = eta1stJet;
              phi2ndJet = phi1stJet;
              qglik2ndJet = qglik1stJet;
              btag2ndJet = btag1stJet;
              isB2ndJet = isB1stJet;

	      pt1stJet = JetPt->at(nJet);
              eta1stJet = JetEta->at(nJet);
              phi1stJet = JetPhi->at(nJet);
	      qglik1stJet = JetQGLikelihood->at(nJet);
              btag1stJet = JetBTagger->at(nJet);
              isB1stJet = (JetIsBtaggedWithSF->at(nJet) > 0);

	    } 
            else if (pt2ndJet < JetPt->at(nJet)) {
	      pt2ndJet = JetPt->at(nJet);
              eta2ndJet = JetEta->at(nJet);
              phi2ndJet = JetPhi->at(nJet);
	      qglik2ndJet = JetQGLikelihood->at(nJet);
              btag2ndJet = JetBTagger->at(nJet);
              isB2ndJet = (JetIsBtaggedWithSF->at(nJet) > 0);

	    }
	    nInJets++;
	  } else { 
            nExtraJets++;
	  }
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
        
      if(debug_) cout<<"find leading leptons"<<endl;
      float pt1stLep = 0.0001;
      float pt2ndLep = 0.0001;
      float eta1stLep = 0.0001;
      float eta2ndLep = 0.0001;
      float phi1stLep = 0.0001;
      float phi2ndLep = 0.0001;
      for (unsigned int nLep=0; nLep<LepPt->size(); nLep++) {
	if (pt1stLep < LepPt->at(nLep)) {
	  pt2ndLep = pt1stLep;
          eta2ndLep = eta1stLep;
          phi2ndLep = phi1stLep;
	  pt1stLep = LepPt->at(nLep);
          eta1stLep = LepEta->at(nLep);
          phi1stLep = LepPhi->at(nLep);
	} else if (pt2ndLep < LepPt->at(nLep)) {
	  pt2ndLep = LepPt->at(nLep);
          eta2ndLep = LepEta->at(nLep);
          phi2ndLep = LepPhi->at(nLep);
	}
      }  
   
      int fsstart,fsend;
      if (abs(Z2Flav->at(0))==121) {
        fsstart=0;
        fsend=2;
      } else {
        fsstart=1;
        fsend=3;
      }

      if(debug_) cout<<"set preferred type either resolved or merged"<<endl;
      int preferType = 0;    // choose merged or resolved
                             // and dump for synchronization 
      if (ZZMass->size() == 1 && abs(ZZCandType->at(0)) == 1) {
	if (ZZMass->at(0) > 400) {
	  nPassMerged++;
	}
	preferType = 1;
      }
      else if (ZZMass->size() == 1 && abs(ZZCandType->at(0)) == 2) {
	if (ZZMass->at(0) > 400) {
	  nPassResol++;
	}
        preferType = 2;
      }
      else if (ZZMass->size() == 2 && abs(ZZCandType->at(0)) == 1) {
	       nPassMerged++; nPassResol++;
	       // merged -> resolved (but ask for J quality)
               if ((Z2Pt->at(0) > 200. && Z1Pt->at(0) > 300. && Z1tau21->at(0) < 0.6)  || Z1Pt->at(1)<100. ) preferType = 1;
               else preferType = 2; 
      }
      // end dump for synchronization and choice

      if(debug_) cout<<"apply trigger weights"<<endl;
      bool weightMCtrig = true;
      if (weightMCtrig) {

	if (process > 0) {  //mc

          if ( abs(GenLep1Id)==11) {
            int bin1 = eff_mc_se->FindFixBin(eta1stLep,TMath::Min(pt1stLep,(float)199.0));
            int bin2 = eff_mc_sm->FindFixBin(eta2ndLep,TMath::Min(pt2ndLep,(float)199.0));

            double mceff_de = eff_mc_de1->GetEfficiency(bin1)*eff_mc_dm2->GetEfficiency(bin2);
            eventWeight *= mceff_de;
          } else {
            int bin1 = eff_mc_sm->FindFixBin(eta1stLep,TMath::Min(pt1stLep,(float)199.0));
            int bin2 = eff_mc_se->FindFixBin(eta2ndLep,TMath::Min(pt2ndLep,(float)199.0));

            double mceff_dm = eff_mc_dm1->GetEfficiency(bin1)*eff_mc_de2->GetEfficiency(bin2);
            eventWeight *= mceff_dm;
          }

        }

      }

      if(debug_) cout<<"fill histos"<<endl;
      for(int rs=fsstart; rs<fsend; rs++){
      
	for(unsigned int theCand=0; theCand<Z1Mass->size(); theCand++){
	
	  if (!CR) {   // ZZ2l2q events 

            if (abs(ZZCandType->at(theCand)) != preferType) continue;

            // random number to smear pruned jet mass
            //float Z1Mass_smear = Z1Mass->at(theCand);

            /*
            if(preferType==1 && process<0){
               TRandom3 rand;
               rand.SetSeed(abs(helcosthetaZ1->at(theCand)*1000000));
               float smear = rand.Gaus(0,1);
               float factor = 1.23;
               float sigma = sqrt(factor*factor-1.0)*8.01;
               Z1Mass_smear = sigma*smear+Z1Mass->at(theCand); 
               if(Z1Mass_smear<0) Z1Mass_smear = 0;
               if(debug_)
                 cout<<"Z1Mass_smear "<<Z1Mass_smear<<" Z1Mass "<<Z1Mass->at(theCand)<<endl;
            }
            */

            if(debug_)cout<<"mZZ cut"<<endl;

            if(preferType==1){ // merged
              if(ZZMass->at(theCand)<minMZZmerged) continue;
            }

            if(preferType==2){ // resolved
              if(isRefitSB && ZZMassRefit->at(theCand)<minMZZresolved)// refitting apply on sideband data as well
                continue;
              if(!isRefitSB && ZZMass->at(theCand)<minMZZresolved)// Not refitting apply on sideband data as well
                continue;
            }

            if(debug_) cout<<"blind higgs region with Z1Mass->at(theCand)"<<Z1Mass->at(theCand)<<" "<<Dsname[d]<<endl;
            if (Z1Mass->at(theCand) > 105. && Z1Mass->at(theCand) < 135.) continue;

	    if(debug_) cout<<"redefine type"<<endl;
	    int typ = -1;
	    if (abs(ZZCandType->at(theCand)) == 1) {
	      if (Z1Mass->at(theCand) > 60. && Z1Mass->at(theCand) < 105.) typ = 2;
	      else typ = 1;
	    } else {
	      if (Z1Mass->at(theCand) > 60. && Z1Mass->at(theCand) < 105.) typ = 3;
	      else typ = 0;
	    }

            // correct tau21
	    float t12weight = 1.;
            if (process == 3 && abs(ZZCandType->at(theCand)) == 1) {
               t12weight = 1.1;
              /*for (int itau = 0; itau < 24; itau++) {
                if (Z1tau21->at(theCand) > tau21bin[itau] && Z1tau21->at(theCand) < tau21bin[itau+1]) 
                   t12weight = 1.+tau21corr[itau];
              }*/
            }

            float mela = 1./(1.+getDZjjspin0Constant(ZZMass->at(theCand))*(pqqZJJ_VAMCFM->at(theCand)/p0plus_VAJHU->at(theCand)));
            float mela2 = 1./(1.+getDZjjspin2Constant(ZZMass->at(theCand))*(pqqZJJ_VAMCFM->at(theCand)/p2bplus_VAJHU->at(theCand)));
            float vbfmela = ((phjj_VAJHU_highestPTJets->at(theCand) > 0. && nExtraJets > 1) ? 1./(1.+getDVBF2jetsConstant(ZZMass->at(theCand))*(phjj_VAJHU_highestPTJets->at(theCand)/pvbf_VAJHU_highestPTJets->at(theCand))) : -1.); 
            // resovled
            if (typ==0 || typ==3) {  // resolved inclusive
              if (nExtraJets > 1 && vbfmela > 1.043-460./(ZZMass->at(theCand)+634.)) typ=typ+8; // vbf tag
              else if(isB1stJet && isB2ndJet) typ=typ+4; // b tag
              //else untagged
            }
            if(debug_) cout<<"typ==0 || typ==3 become resolved untagged (without vb tagged and b tagged)"<<endl;

            // merged
            if (typ==1 || typ==2) { // merged inclusive
               if (nExtraJets > 1 && vbfmela > 1.043-460./(ZZMass->at(theCand)+634.)) typ=typ+8; // vbf tag
               else if (btag1stSubjet > 0.5426 && btag2ndSubjet > 0.5426) typ=typ+4; // btag
               // else untagged 
            }
            if(debug_) cout<<"typ==1 || typ==2 become resolved untagged (without vb tagged and b tagged)"<<endl;

	    if ((typ==0 || typ==3) && nExtraJets > 1 && vbfmela > 1.043-460./(ZZMass->at(theCand)+634.)) typ=typ+8;
	    if ((typ==1 || typ==2) && nExtraJets > 1 && vbfmela > 1.043-460./(ZZMass->at(theCand)+634.)) typ=typ+8;             
	    if ((typ==0 || typ==3) && isB1stJet && isB2ndJet) typ=typ+4;
	    if ((typ==1 || typ==2) && btag1stSubjet > 0.5426 && btag2ndSubjet > 0.5426) typ=typ+4;

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
	    } 
            else {
	      tmvaJetPt1 = (float)pt1stSubjet;
	      tmvaJetPt2 = (float)pt2ndSubjet;
	    }              
	    
            if(debug_) cout<<"event filter final"<<endl;
            if(preferType==1)
              if(!passAdditionalCuts(tmvaZ2Mass,tmvaZ1Pt,tmvaZ2Pt,tmvaZ1tau21,true)) continue;
            if(preferType==2)
              if(!passAdditionalCuts(tmvaZ2Mass,tmvaZ1Pt,tmvaZ2Pt,tmvaZ1tau21,false)) continue;
	    
            //////////////////////
	    h1[0][process][rs][typ]->Fill(ZZMass->at(theCand),eventWeight*t12weight);
	    h1[1][process][rs][typ]->Fill(ZZPt->at(theCand),eventWeight*t12weight);
	    h1[2][process][rs][typ]->Fill(ZZMassRefit->at(theCand),eventWeight*t12weight);
	    h1[3][process][rs][typ]->Fill(Z1Mass->at(theCand),eventWeight*t12weight);
	    h1[4][process][rs][typ]->Fill(Z2Mass->at(theCand),eventWeight*t12weight);
	    h1[5][process][rs][typ]->Fill(Z1Pt->at(theCand),eventWeight*t12weight);
	    h1[6][process][rs][typ]->Fill(Z2Pt->at(theCand),eventWeight*t12weight);
	    h1[7][process][rs][typ]->Fill(Z2Flav->at(theCand),eventWeight*t12weight);
	    
	    if (typ==0 || typ==3 || typ==4 || typ==7 || typ==8 || typ==11) {   // only resolved
	      h1[8][process][rs][typ]->Fill(pt1stJet,eventWeight*t12weight);
	      h1[9][process][rs][typ]->Fill(pt2ndJet,eventWeight*t12weight);

              float dPhi = deltaPhi(phi1stJet,phi2ndJet);
              float dEta = eta1stJet-eta2ndJet;
              float dR = TMath::Sqrt(dPhi*dPhi+dEta*dEta);

              h2_mjjVSdR[process][rs][typ]->Fill(dR,Z1Mass->at(theCand),eventWeight*t12weight);

	    } 
            else {  // only merged
	      h1[8][process][rs][typ]->Fill(pt1stSubjet,eventWeight*t12weight);
	      h1[9][process][rs][typ]->Fill(pt2ndSubjet,eventWeight*t12weight);

              float dPhi = deltaPhi(phi1stSubJet,phi2ndSubJet);
              float dEta = eta1stSubJet-eta2ndSubJet;
              float dR = TMath::Sqrt(dPhi*dPhi+dEta*dEta);

              h2_mjjVSdR[process][rs][typ]->Fill(dR,Z1Mass->at(theCand),eventWeight*t12weight);

	    }
	    
	    if (typ==0 || typ==3 || typ==4 || typ==7  || typ==8 || typ==11) {   // only resolved
	      for (unsigned int nJet=0; nJet<JetPt->size(); nJet++) {
		if (JetIsInZZCand->at(nJet) && JetQGLikelihood->at(nJet) > -800. /*remove subjets of fat jet also included in this collection! */) {
		  h1[11][process][rs][typ]->Fill(JetBTagger->at(nJet),eventWeight*t12weight);
		  h1[10][process][rs][typ]->Fill(JetQGLikelihood->at(nJet),eventWeight*t12weight); 
		} 
	      }
	    } else {  
	      for (unsigned int nJet=0; nJet<JetPt->size(); nJet++) {
		if (JetIsInZZCand->at(nJet) && JetQGLikelihood->at(nJet) < -800. /*subjets of fat jet also included in this collection! */) {
		  h1[11][process][rs][typ]->Fill(JetBTagger->at(nJet),eventWeight*t12weight);
		} 
	      }
	    } 
	    
	    h1[12][process][rs][typ]->Fill(pt1stLep,eventWeight*t12weight);
	    h1[13][process][rs][typ]->Fill(pt2ndLep,eventWeight*t12weight);
	    
       	    h1[14][process][rs][typ]->Fill(abs(helcosthetaZ1->at(theCand)),eventWeight*t12weight);
	    h1[15][process][rs][typ]->Fill(helcosthetaZ2->at(theCand),eventWeight*t12weight);
	    h1[16][process][rs][typ]->Fill(costhetastar->at(theCand),eventWeight*t12weight);
	    h1[17][process][rs][typ]->Fill(phistarZ1->at(theCand),eventWeight*t12weight);
	    h1[18][process][rs][typ]->Fill(fabs(helphi->at(theCand)),eventWeight*t12weight);
	    
	    h1[19][process][rs][typ]->Fill(nExtraJets,eventWeight*t12weight);
	    h1[20][process][rs][typ]->Fill(Met,eventWeight*t12weight);

	    if (typ==1 || typ==2 || typ==5 || typ==6  || typ==9 || typ==10) 
	      h1[21][process][rs][typ]->Fill(Z1tau21->at(theCand),eventWeight*t12weight);   // only merged
	    else 
	      h1[22][process][rs][typ]->Fill(tmvaJetQGLikelihood1*tmvaJetQGLikelihood2,eventWeight*t12weight);     // only resolved

            h1[23][process][rs][typ]->Fill(mela,eventWeight*t12weight);
            h1[24][process][rs][typ]->Fill(vbfmela,eventWeight*t12weight);
            h1[26][process][rs][typ]->Fill(mela2,eventWeight*t12weight);
	    if (rs == onlyOneLep && (typ == 3 || typ == 7 || typ == 11) ) {
	      hmass[process][16]->Fill(mela,eventWeight*t12weight);
	      hmass[process][17]->Fill(mela2,eventWeight*t12weight);
	      hmass[process][18]->Fill(vbfmela,eventWeight*t12weight);
	    }
            if(debug_)cout<<"hmass for alpha method"<<endl;
            if (typ==0 || typ==3 || typ==4 || typ==7 || typ==8 || typ==11) {
              if (rs == onlyOneLep) {
                hmass[process][typ]->Fill(ZZMass->at(theCand),eventWeight*t12weight);
                hmassRefit[process][typ]->Fill(ZZMassRefit->at(theCand),eventWeight*t12weight);
                hmass_up[process][typ]->Fill(ZZMass_up->at(theCand),eventWeight*t12weight);
                hmass_down[process][typ]->Fill(ZZMass_dn->at(theCand),eventWeight*t12weight);
                if(debug_)cout<<"non-btagged cats "<<ZZMass->at(theCand)<<" "<<ZZMass_up->at(theCand)<<" "<<ZZMass_dn->at(theCand)<<endl;
                if ( typ == 0 || typ==8 ) {
                  hmass[process][12]->Fill(ZZMass->at(theCand),eventWeight*t12weight);
                  hmass_up[process][12]->Fill(ZZMass_up->at(theCand),eventWeight*t12weight);
                  hmass_down[process][12]->Fill(ZZMass_dn->at(theCand),eventWeight*t12weight);
                }
                if ( typ == 3 || typ==11 ) {
                  hmass[process][15]->Fill(ZZMass->at(theCand),eventWeight*t12weight);
                  hmass_up[process][15]->Fill(ZZMass_up->at(theCand),eventWeight*t12weight);
                  hmass_down[process][15]->Fill(ZZMass_dn->at(theCand),eventWeight*t12weight);
                }
              }
              if (mela > 0.5) h1[25][process][rs][typ]->Fill(ZZMassRefit->at(theCand),eventWeight*t12weight);
            } else {

              if (rs == onlyOneLep) {
                if(debug_)cout<<"merged"<<endl;
                hmass[process][typ]->Fill(ZZMass->at(theCand),eventWeight*t12weight);
                hmass_up[process][typ]->Fill(ZZMass_up->at(theCand),eventWeight*t12weight);
                hmass_down[process][typ]->Fill(ZZMass_dn->at(theCand),eventWeight*t12weight);
                if ( typ == 1 || typ== 9 ) {
                  hmass[process][13]->Fill(ZZMass->at(theCand),eventWeight*t12weight);
                  hmass_up[process][13]->Fill(ZZMass_up->at(theCand),eventWeight*t12weight);
                  hmass_down[process][13]->Fill(ZZMass_dn->at(theCand),eventWeight*t12weight);
                }
                if ( typ == 2 || typ== 10 ) {
                  hmass[process][14]->Fill(ZZMass->at(theCand),eventWeight*t12weight);
                  hmass_up[process][14]->Fill(ZZMass_up->at(theCand),eventWeight*t12weight);
                  hmass_down[process][14]->Fill(ZZMass_dn->at(theCand),eventWeight*t12weight);
                }
              }
              if (mela > 0.5) h1[25][process][rs][typ]->Fill(ZZMass->at(theCand),eventWeight*t12weight);
            }

            if(debug_) cout<<"end plotting"<<endl;
	    
	  } else { // control region for QG (only fille some variables)

	    for (unsigned int nJet=0; nJet<JetPt->size(); nJet++) {
	      if (JetIsInZZCand->at(nJet) && JetQGLikelihood->at(nJet) > -800. /*remove subjets of fat jet also included in this collection! */) {
		if (deltaPhi(phi1stLep, JetPhi->at(nJet)) < 2.2 || deltaPhi(phi2ndLep, JetPhi->at(nJet)) < 2.2) continue;
	      }
	    }
	    
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

	  } // if CR

	} // candidate type : resolved or merged

      } // fs 	 

    } // tree entry

    if (mass[d] > 0 && enforceNarrowWidth) NEvtNarrow[d] = NEvtNarrow[d] / entries;
    
    if (process==0 && sync) { 
      myfile.close(); 
      float eff = float(nPassMerged)/float(NGenEvt[d]);
      cout<<"Pass merged analysis = "<<nPassMerged<<"/"<< NGenEvt[d] << " (" << eff*100. << "%)"<<endl;
      eff = float(nPassResol)/float(NGenEvt[d]);
      cout<<"Pass resolved analysis = "<<nPassResol<<"/"<< NGenEvt[d] << " (" << eff*100. << "%) ..."<<endl;
      break;
    }  
  }

  cout<<"if draw"<<endl;

  if (draw) { 
    TCanvas c1;
    c1.cd();
    
    TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.35,0.95,0.97);
    pad1->Draw();
    TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.02,0.95,0.35);
    pad2->Draw();

    for(int rs=0; rs<nFS; rs++){   //ee, mumu, or all

      for(int v=0; v<nVariables; v++){ // variable

	for(int nt=0; nt<nType; nt++){ // type
	  
          // instead of THStack
          h1[v][1][rs][nt]->Add(h1[v][2][rs][nt]); // TT+WW+DY

	  h1[v][0][rs][nt]->GetXaxis()->SetTitle(varXLabel[v].c_str());
	  h1[v][0][rs][nt]->GetYaxis()->SetTitle(varYLabel[v].c_str());
	  h1[v][1][rs][nt]->SetFillStyle(1);
	  h1[v][1][rs][nt]->SetMinimum(0.1);
	  h1[v][1][rs][nt]->SetLineColor(kGreen+2);
	  h1[v][1][rs][nt]->SetFillColor(kGreen+2);
	  h1[v][2][rs][nt]->SetFillStyle(1);
	  h1[v][2][rs][nt]->SetMinimum(0.1);
	  h1[v][2][rs][nt]->SetLineColor(kYellow+2);
	  h1[v][2][rs][nt]->SetFillColor(kYellow+2);
	  h1[v][0][rs][nt]->SetMarkerStyle(20);
	  h1[v][0][rs][nt]->SetMinimum(0.1);
	  
	}
      }
    }

    for(int nt=0; nt<nType+7; nt++){
      
      if (nt<nType+4) {
	densityHist(hmass[0][nt]);
	densityHist(hmass[1][nt]);
	densityHist(hmass[2][nt]);
      } else {
        hmass[0][nt]->Scale(1./hmass[0][nt]->Integral());
	hmass[1][nt]->Scale(1./hmass[1][nt]->Integral());
	hmass[2][nt]->Scale(1./hmass[2][nt]->Integral());
      }
      
      hmass[0][nt]->SetMarkerStyle(20);
      
    }

    ofstream aa("unblind.txt");
    
    if (draw > 1) {

	for(int rs=0; rs<nFS; rs++){   // ee, mumu, or all

	  for(int v=0; v<nVariables; v++){

	    for(int nt=0; nt<nType; nt++){
	      
	      pad1->cd();  
	      
	      TLegend *legend = new TLegend(0.70,0.65,0.95,0.90,NULL,"brNDC");
	      legend->SetBorderSize(     0);
	      legend->SetFillColor (     0);
	      legend->SetTextAlign (    12);
	      legend->SetTextFont  (    42);
	      legend->SetTextSize  (0.03);
	      
	      legend->AddEntry(h1[v][0][rs][nt], processLabel[0].c_str() , "p");
	      for(int ipr=1; ipr<nProcesses; ipr++){ 
		string legType = "f";
		legend->AddEntry(h1[v][ipr][rs][nt], processLabel[ipr].c_str() , legType.c_str());
	      }
	      
	      if (nt==0 || nt==1 || nt==4 || nt==5 || nt==8 || nt==9 || CR || unblind) {
		
		 h1[v][0][rs][nt]->Draw("e");
		
		 h1[v][1][rs][nt]->Draw("histsame"); 
		 h1[v][2][rs][nt]->Draw("histsame");
	         h1[v][0][rs][nt]->Draw("esame");
		
	     } else {
		
		h1[v][1][rs][nt]->Draw("hist"); 
		h1[v][2][rs][nt]->Draw("histsame");

	      }
	      
	      gPad->SetLogy(varLogy[v]);
	      gPad->SetLogx(varLogx[v]);
	      
	      legend->Draw("same"); 
	      
	      pad2->cd();
	      
	      gPad->SetLogy(0);
	      TH1F* ratio = (TH1F*)h1[v][0][rs][nt]->Clone();//clone data
	      TH1F* mc = (TH1F*)h1[v][1][rs][nt]->Clone();//clone dy+tt/ww+vz
	      ratio->Divide(mc);
	      ratio->SetMinimum(0.5);
	      ratio->SetMaximum(2.0); 
	      ratio->GetYaxis()->SetTitle("Data/MC");
	      if (!(nt==0 || nt==1 || nt==4 || nt==5 || nt==8 || nt==9 || CR || unblind)) {
		ratio->SetLineColor(kWhite);
		ratio->SetMarkerColor(kWhite);
	      }
	      ratio->Draw("e");
	      TLine line(h1[v][0][rs][nt]->GetXaxis()->GetBinLowEdge(1),1.,
			 h1[v][0][rs][nt]->GetXaxis()->GetBinUpEdge(h1[v][0][rs][nt]->GetNbinsX()-1),1.);
	      line.SetLineColor(kRed);
	      line.SetLineStyle(kDashed);
	      line.Draw("same");  
	      cout<<varName[v].c_str()<<" "<<typeS[nt].c_str()<<" "<<mZZtype.c_str()<<endl;
	      c1.SaveAs(Form("~/www/HZZ/ZZ2l2q/%s/%s_%s_%s_mZZ%s.png",dirout.c_str(),varName[v].c_str(),sFS[rs].c_str(),typeS[nt].c_str(),mZZtype.c_str()));
              cout<<varName[v].c_str()<<" "<<typeS[nt].c_str()<<" "<<mZZtype.c_str()<<endl;

	  } // type

	} // variable

      }  // final state

    }  // if draw
    
    aa.close();

    cout << "NOW PAS-style" << endl;

       
  } else {
    
    if (!CR) {

      cout<<"writing output"<<endl;

      outputFile->cd();
      for(int nt=0; nt<nType+4; nt++){

        cout<<"mass"<<endl;
	hmass[0][nt]->Write();
	hmass[1][nt]->Write();
	hmass[2][nt]->Write();

        cout<<"mass refit"<<endl;
        hmassRefit[0][nt]->Write();
        hmassRefit[1][nt]->Write();
        hmassRefit[2][nt]->Write();

        cout<<"mass up"<<endl; 
	hmass_up[0][nt]->Write();
	hmass_up[1][nt]->Write();
	hmass_up[2][nt]->Write();

        cout<<"mass dn"<<endl;
	hmass_down[0][nt]->Write();
	hmass_down[1][nt]->Write();
	hmass_down[2][nt]->Write(); 

      }

      cout<<"close output file..."<<endl;
      outputFile->Close(); 
      cout<<"close output file"<<endl;
    }
  }
  mengch.close();
}
