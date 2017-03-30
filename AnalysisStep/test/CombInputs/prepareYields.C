///
/// usage:
/// -specify parameters (input/output location, luminosity, m4l window) at the end of this file
/// -run with:
///   root -q -l -b prepareYields.C++
/// -later, once yields have been stored, use parameter 0 to not rerun over all MC:
///   root -q -l -b prepareYields.C++(0)
///

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

#include "TCanvas.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPaveStats.h"
#include "TSpline.h"
#include "TSystem.h"
#include "TTree.h"

#include "../Plotter/tdrstyle.C"
#include "../Plotter/plotUtils.C"
#include "ZZAnalysis/AnalysisStep/src/kFactors.C"
#include "ZZAnalysis/AnalysisStep/src/cConstants.cc"
#include "ZZAnalysis/AnalysisStep/src/Discriminants.cc"
#include "ZZAnalysis/AnalysisStep/src/Category.cc"
#include "ZZAnalysis/AnalysisStep/src/bitops.cc"
#include "ZZAnalysis/AnalysisStep/interface/FinalStates.h"

using namespace std;

#define JUST125 0

#define DEBUG 0
#define UNBLINDSR 0
#define PRINTLATEXTABLES 0
#define SPLITVHHADLEP 1

#define APPLYKFACTORS 1
#define USEPUWEIGHT 1
#define USEVHMETTAG 1
#define USEBBH 0
#define MERGE2E2MU 1 // Don't change this unless you really know what you are doing.

#define DOZPLUSXFROMRUN2COMBINEDSHAPE 1
#define BUILDZPLUSXYIELDFROMSSCR 1 // Obtain Z+X yields in every category from rerunning over the SS CR. Overwrites hardcoded global arrays.
#define RESCALEZPLUSXTOCOMBINEDINCLUSIVE 1 // Renormalize all Z+X yields to the hardcoded SS/OS-combined inclusive numbers.

//Combined 0S+SS inclusive Z+X prediction
//Here: numbers for the full 2016 data set, computed by Roberto on March 3rd 2017
Float_t normCombFullRange4e    = 21.1;
Float_t normCombFullRange4mu   = 34.4;
Float_t normCombFullRange2e2mu = 59.9;

//Normalization of Z+X in final states and categories. 
//These numbers get overwritten if the BUILDZPLUSXYIELDFROMSSCR flag is on.
//Here: SS-based numbers from 170222 trees, for the full 2016 data set.
Float_t normSSFullRange4e[8] = {
  16.1702, //Untagged
  1.35714, //VBF1jTagged
  0.933429, //VBF2jTagged
  0.182454, //VHLeptTagged
  0.295943, //VHHadrTagged
  0.318589, //ttHTagged
  0.281382, //VHMETTagged    
  19.5391, //inclusive
};
Float_t normSSFullRange4mu[8] = {
  28.0966,
  2.43438,
  2.297,
  0.458604,
  0.929289,
  0.80248,
  0.945141,
  35.9635,
};
Float_t normSSFullRange2e2mu[8] = {
  44.4236,
  3.48285,
  3.12444,
  0.854668,
  1.28802,
  1.11163,
  1.17629,
  55.4615,
};



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////// global variables /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

enum Process {ggH=0, qqH=1, WH=2, WH_lep=3, WH_had=4, ZH=5, ZH_lep=6, ZH_had=7, ttH=8, bbH=9, qqZZ=10, ggZZ=11, ttZ=12, WWZ=13};
const int nProcesses = 14;
string sProcess[nProcesses] = {"ggH", "qqH", "WH", "WH_lep", "WH_had", "ZH", "ZH_lep", "ZH_had", "ttH", "bbH", "qqZZ", "ggZZ", "ttZ", "WWZ"};
bool isSignal[nProcesses] = {1,1,1,1,1,1,1,1,1,1,0,0,0,0,};
bool useProcess[nProcesses] = {1,1,!SPLITVHHADLEP,SPLITVHHADLEP,SPLITVHHADLEP,!SPLITVHHADLEP,SPLITVHHADLEP,SPLITVHHADLEP,1,USEBBH,1,1,0,0};
bool useProcessInMeas[nProcesses] = {1,1,!SPLITVHHADLEP,SPLITVHHADLEP,SPLITVHHADLEP,!SPLITVHHADLEP,SPLITVHHADLEP,SPLITVHHADLEP,1,0,1,1,0,0};

const int nMHPoints = 6;
string  sMHPoint[nMHPoints] = {"","120","124","125","126","130",};
Float_t fMHPoint[nMHPoints] = {0., 120., 124., 125., 126., 130.,};
int indexOf125 = 3;
Int_t nMHPointsProcess[nProcesses] = {5,5,5,5,5,5,5,5,5,0,0,0,0,0};
bool hasMHPoint[nProcesses][nMHPoints] = {
  {0,1,1,1,1,1,},//ggH
  {0,1,1,1,1,1,},//qqH
  {0,1,1,1,1,1,},//WH
  {0,1,1,1,1,1,},//WH_lep
  {0,1,1,1,1,1,},//WH_had
  {0,1,1,1,1,1,},//ZH
  {0,1,1,1,1,1,},//ZH_lep
  {0,1,1,1,1,1,},//ZH_had
  {0,1,1,1,1,1,},//ttH
  {0,0,0,1,0,0,},//bbH
  {1,0,0,0,0,0,},//qqZZ
  {1,0,0,0,0,0,},//ggZZ
  {1,0,0,0,0,0,},//ttZ
  {1,0,0,0,0,0,},//WWZ
};
const int orderOfPolynomial[nProcesses] = {2,2,2,2,2,2,2,2,2,0,0,0,0,0};

enum FinalState {fs4mu=0, fs4e=1, fs2e2mu=2, fs2mu2e=3};
const int nFinalStates = 4;
string sFinalState[nFinalStates+1] = {"4mu", "4e", "2e2mu", "2mu2e", "4l"};
Int_t fsMarkerStyle[nFinalStates+1] = {20,22,21,33,29};
int fs2[nFinalStates] = {1,0,2,4};
Float_t fsROSSS[nFinalStates] = { 1.22, 0.97, 1.30, 0.98 };
Float_t bkgdCompUncZpX[nFinalStates] = { 0.35, 0.32, 0.34, 0.34 }; // Numbers sent by Pedja on July 26th 2016

// Moriond 2017 categorization 
const int nCategories = 7;
string sCategory[nCategories+1] = {
  "UnTagged",
  "VBF1jTagged",
  "VBF2jTagged", 
  "VHLeptTagged",
  "VHHadrTagged",
  "ttHTagged",
  "VHMETTagged",
  "inclusive",
};
string sCategoryForSyst[nCategories] = {
  "UnTagged",
  "OneJetTagged",
  "VBFTagged", 
  "VHLeptTagged",
  "VHHadrTagged",
  "ttHTagged",
  "VHMETTagged",
};
string sCatLatex[nCategories] = {"Untagged","VBF-1j","VBF-2j","VH-lept.","VH-hadr.","ttH","VH-MET"};
string sCatLatexLong[nCategories+1] = {"Untagged","VBF-1jet-tagged","VBF-2jet-tagged","VH-leptonic-tagged","VH-hadronic-tagged","ttH-tagged","VH-MET-tagged","Inclusive"};
Color_t catColor[nCategories+1] = {kBlue-9, kCyan-6, kGreen-6, kRed-7, kOrange+6, kMagenta-6, kBlack, kGray };

enum ResonantStatus {resonant=0, nonresonant=1};
const int nRS = 2;
string sResonantStatus[nRS+1] = {"resonant", "nonresonant", "allres"};

enum SystShift {nominal=0, jecUp=1, jecDn=2, btagUp=3, btagDn=4};
const int nSystShifts = 5;
string sSystShift[nSystShifts] = {"SystNominal","JecUp","JecDn","bTagUp","bTagDn"};

Double_t deltaR(Double_t e1, Double_t p1, Double_t e2, Double_t p2) {
  Double_t deltaPhi = acos(cos(p1-p2));
  return TMath::Sqrt((e1-e2)*(e1-e2) + deltaPhi*deltaPhi);
}

Float_t fround(Float_t val, Int_t nDigits) { 
  Int_t f = pow(10,nDigits);
  return roundf(val * f) / f; 
}
string sround(Float_t val, Int_t nDigits) { 
  return string(Form(("%."+string(Form("%i",nDigits))+"f").c_str(),val)); 
}
string fixWidth(string str, unsigned n, bool atBeginning) {
  if(str.length()>n){
    cout<<"Error in function fixWidth('"<<str<<"',"<<n<<") : string '"<<str<<"' is too long."<<endl;
    return "";
  }else{
    string res = "";
    if(atBeginning) res += str;
    for(int i=0; i<(int)n-(int)(str.length()); i++) res += " ";
    if(!atBeginning) res += str;
    return res;
  }
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////// compute and save yields /////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void computeYields(string inputPathSignal, string inputPathqqZZ, string inputPathggZZ, double lumi, double sqrts, double m4lMin, double m4lMax)
{

  TH1::SetDefaultSumw2(true);

  //TFile* ggZZKFactorFile = TFile::Open("../../data/kfactors/Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
  //TSpline3* sp = (TSpline3*)ggZZKFactorFile->Get("sp_kfactor_Nominal");

  const int nDatasets = 16;
  string datasets[nDatasets] = {
    "ggH",
    "VBFH",
    "WplusH",
    "WminusH",
    "ZH",
    "ttH",
    "bbH",
    "ZZTo4l",
    "ggTo4e_Contin_MCFM701",//"ggZZ4e",
    "ggTo4mu_Contin_MCFM701",//"ggZZ4mu",
    "ggTo4tau_Contin_MCFM701",//"ggZZ4tau",
    "ggTo2e2mu_Contin_MCFM701",//"ggZZ2e2mu",
    "ggTo2e2tau_Contin_MCFM701",//"ggZZ2e2tau",
    "ggTo2mu2tau_Contin_MCFM701",//"ggZZ2mu2tau",
    "TTZJets_M10_MLM",
    "WWZ",
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
  Float_t PFMET;
  Float_t PFMET_jesUp;
  Float_t PFMET_jesDn;
  Short_t genExtInfo;
  Short_t Nvtx;
  Short_t NObsInt;
  Float_t NTrueInt;
  Float_t genHEPMCweight;
  Float_t PUWeight;
  Float_t dataMCWeight;
  Float_t overallEventWeight;
  Float_t KFactor_QCD_ggZZ_Nominal;
  Float_t KFactor_EW_qqZZ;
  Float_t KFactor_QCD_qqZZ_dPhi;
  Float_t KFactor_QCD_qqZZ_M;
  Float_t KFactor_QCD_qqZZ_Pt;
  Float_t xsec;
  Short_t ZZsel;
  Float_t ZZMass;
  Float_t ZZPt;
  Short_t Z1Flav;
  Short_t Z2Flav;
  Float_t p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal;
  Float_t p_JJVBF_SIG_ghv1_1_JHUGen_JECUp;
  Float_t p_JJVBF_SIG_ghv1_1_JHUGen_JECDn;
  Float_t p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal;
  Float_t p_JJQCD_SIG_ghg2_1_JHUGen_JECUp;
  Float_t p_JJQCD_SIG_ghg2_1_JHUGen_JECDn;
  Float_t p_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
  Float_t p_JVBF_SIG_ghv1_1_JHUGen_JECUp;
  Float_t p_JVBF_SIG_ghv1_1_JHUGen_JECDn;
  Float_t pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
  Float_t pAux_JVBF_SIG_ghv1_1_JHUGen_JECUp;
  Float_t pAux_JVBF_SIG_ghv1_1_JHUGen_JECDn;
  Float_t p_JQCD_SIG_ghg2_1_JHUGen_JECNominal;
  Float_t p_JQCD_SIG_ghg2_1_JHUGen_JECUp;
  Float_t p_JQCD_SIG_ghg2_1_JHUGen_JECDn;
  Float_t p_HadWH_SIG_ghw1_1_JHUGen_JECNominal;
  Float_t p_HadWH_SIG_ghw1_1_JHUGen_JECUp;
  Float_t p_HadWH_SIG_ghw1_1_JHUGen_JECDn;
  Float_t p_HadZH_SIG_ghz1_1_JHUGen_JECNominal;
  Float_t p_HadZH_SIG_ghz1_1_JHUGen_JECUp;
  Float_t p_HadZH_SIG_ghz1_1_JHUGen_JECDn;
  vector<Float_t> *LepEta = 0;
  vector<Float_t> *LepPhi = 0;
  Short_t nExtraLep;
  vector<Float_t> *ExtraLepEta = 0;
  vector<Float_t> *ExtraLepPhi = 0;
  Short_t nExtraZ;
  Short_t nCleanedJets;
  Short_t nCleanedJetsPt30;
  Short_t nCleanedJetsPt30_jecUp;
  Short_t nCleanedJetsPt30_jecDn;
  Short_t nCleanedJetsPt30BTagged_bTagSF;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jecUp;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jecDn;
  Short_t nCleanedJetsPt30BTagged_bTagSFUp;
  Short_t nCleanedJetsPt30BTagged_bTagSFDn;
  //vector<Float_t> *JetPt = 0;
  //vector<Float_t> *JetEta = 0;
  vector<Float_t> *JetPhi = 0;
  //vector<Float_t> *JetMass = 0;
  vector<Float_t> *JetQGLikelihood = 0;
  //Float_t jetPt[99];
  //Float_t jetEta[99];
  Float_t jetPhi[99];
  //Float_t jetMass[99];
  Float_t jetQGL[99]; 
  Float_t GenHMass;
  Float_t GenZ1Phi;
  Float_t GenZ2Phi;
  Float_t GenZ1Flav;
  Float_t GenZ2Flav;
  Float_t GenLep1Eta;
  Float_t GenLep1Phi;
  Short_t GenLep1Id;
  Float_t GenLep2Eta;
  Float_t GenLep2Phi;
  Short_t GenLep2Id;
  Float_t GenLep3Eta;
  Float_t GenLep3Phi;
  Short_t GenLep3Id;
  Float_t GenLep4Eta;
  Float_t GenLep4Phi;
  Short_t GenLep4Id;

  TH1F* hYield[nProcesses][nMHPoints][nFinalStates+1][nCategories+1][nRS+1][nSystShifts];
  for(int pr=0; pr<nProcesses; pr++)
    for(int mp=0; mp<nMHPoints; mp++)
      if(hasMHPoint[pr][mp])
	for(int fs=0; fs<nFinalStates+1; fs++)
	  for(int cat=0; cat<nCategories+1; cat++)
	    for(int rs=0; rs<nRS+1; rs++)
	      for(int sy=0; sy<nSystShifts; sy++){
		hYield[pr][mp][fs][cat][rs][sy] = new TH1F(Form("hYield_%s%s_%s_%s_%s_%s",sProcess[pr].c_str(),sMHPoint[mp].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str(),sResonantStatus[rs].c_str(),sSystShift[sy].c_str()),"",1,0,1);
		hYield[pr][mp][fs][cat][rs][sy]->Sumw2(true);
	      }
  
  int currentProcess;
  int currentFinalState;
  int currentCategory[nSystShifts];
  int currentResStatus;


  //---------- Will loop over all datasets

  for(int d=0; d<nDatasets; d++){

    //----- assign dataset to correct process
    currentProcess = -1;
    if(datasets[d]=="ggH") currentProcess = ggH;
    if(datasets[d]=="VBFH") currentProcess = qqH;
    if(datasets[d]=="WplusH") currentProcess = WH;
    if(datasets[d]=="WminusH") currentProcess = WH;
    if(datasets[d]=="ZH") currentProcess = ZH;
    if(datasets[d]=="ttH") currentProcess = ttH;
    if(datasets[d]=="bbH") currentProcess = bbH;
    if(datasets[d]=="ZZTo4l"||
       datasets[d]=="ZZTo4lamcatnlo") 
      currentProcess = qqZZ;
    if(datasets[d]=="ggZZ4e"||
       datasets[d]=="ggZZ4mu"||
       datasets[d]=="ggZZ4tau"||
       datasets[d]=="ggZZ2e2mu"||
       datasets[d]=="ggZZ2e2tau"||
       datasets[d]=="ggZZ2mu2tau"||
       datasets[d]=="ggTo4e_Contin_MCFM701"||
       datasets[d]=="ggTo4mu_Contin_MCFM701"||
       datasets[d]=="ggTo4tau_Contin_MCFM701"||
       datasets[d]=="ggTo2e2mu_Contin_MCFM701"||
       datasets[d]=="ggTo2e2tau_Contin_MCFM701"||
       datasets[d]=="ggTo2mu2tau_Contin_MCFM701") 
      currentProcess = ggZZ;
    if(datasets[d]=="TTZJets_M10_MLM") currentProcess = ttZ;
    if(datasets[d]=="WWZ") currentProcess = WWZ;

    int initialCurrentProcess = currentProcess;

    for(int mp=0; mp<nMHPoints; mp++){

      if(!hasMHPoint[currentProcess][mp] || (JUST125&&sMHPoint[mp]!="125"&&sMHPoint[mp]!="")) continue;
      
      string inputFileName = string(Form("%s%s%s/ZZ4lAnalysis.root",(currentProcess==ggZZ)?inputPathggZZ.c_str():(currentProcess==qqZZ)?inputPathqqZZ.c_str():inputPathSignal.c_str(),datasets[d].c_str(),sMHPoint[mp].c_str()));
      inputFile[d] = TFile::Open(inputFileName.c_str());
      
      hCounters[d] = (TH1F*)inputFile[d]->Get("ZZTree/Counters");
      NGenEvt[d] = (Long64_t)hCounters[d]->GetBinContent(1);
      gen_sumWeights[d] = (Long64_t)hCounters[d]->GetBinContent(USEPUWEIGHT?40:41);
      partialSampleWeight[d] = lumi * 1000 / gen_sumWeights[d] ;
      
      inputTree[d] = (TTree*)inputFile[d]->Get("ZZTree/candTree");
      inputTree[d]->SetBranchAddress("RunNumber", &nRun);
      inputTree[d]->SetBranchAddress("EventNumber", &nEvent);
      inputTree[d]->SetBranchAddress("LumiNumber", &nLumi);
      inputTree[d]->SetBranchAddress("PFMET", &PFMET);
      inputTree[d]->SetBranchAddress("PFMET_jesUp", &PFMET_jesUp);
      inputTree[d]->SetBranchAddress("PFMET_jesDn", &PFMET_jesDn);
      inputTree[d]->SetBranchAddress("genExtInfo", &genExtInfo);
      inputTree[d]->SetBranchAddress("Nvtx", &Nvtx);
      inputTree[d]->SetBranchAddress("NObsInt", &NObsInt);
      inputTree[d]->SetBranchAddress("NTrueInt", &NTrueInt);
      inputTree[d]->SetBranchAddress("genHEPMCweight",&genHEPMCweight);
      inputTree[d]->SetBranchAddress("PUWeight",&PUWeight);
      inputTree[d]->SetBranchAddress("dataMCWeight",&dataMCWeight);
      inputTree[d]->SetBranchAddress("overallEventWeight", &overallEventWeight);
      if(currentProcess==ggZZ){
	inputTree[d]->SetBranchAddress("KFactor_QCD_ggZZ_Nominal", &KFactor_QCD_ggZZ_Nominal);
      }
      if(currentProcess==qqZZ){
	inputTree[d]->SetBranchAddress("KFactor_EW_qqZZ", &KFactor_EW_qqZZ);
	inputTree[d]->SetBranchAddress("KFactor_QCD_qqZZ_dPhi", &KFactor_QCD_qqZZ_dPhi);
	inputTree[d]->SetBranchAddress("KFactor_QCD_qqZZ_M", &KFactor_QCD_qqZZ_M);
	inputTree[d]->SetBranchAddress("KFactor_QCD_qqZZ_Pt", &KFactor_QCD_qqZZ_Pt);
      }
      inputTree[d]->SetBranchAddress("xsec", &xsec);
      inputTree[d]->SetBranchAddress("ZZsel", &ZZsel);
      inputTree[d]->SetBranchAddress("ZZMass", &ZZMass);
      inputTree[d]->SetBranchAddress("ZZPt", &ZZPt);
      inputTree[d]->SetBranchAddress("Z1Flav", &Z1Flav);
      inputTree[d]->SetBranchAddress("Z2Flav", &Z2Flav);
      inputTree[d]->SetBranchAddress("LepEta", &LepEta);
      inputTree[d]->SetBranchAddress("LepPhi", &LepPhi);
      inputTree[d]->SetBranchAddress("nExtraLep", &nExtraLep);
      inputTree[d]->SetBranchAddress("ExtraLepEta", &ExtraLepEta);
      inputTree[d]->SetBranchAddress("ExtraLepPhi", &ExtraLepPhi);
      inputTree[d]->SetBranchAddress("nExtraZ", &nExtraZ);
      inputTree[d]->SetBranchAddress("nCleanedJets", &nCleanedJets);
      inputTree[d]->SetBranchAddress("nCleanedJetsPt30", &nCleanedJetsPt30);
      inputTree[d]->SetBranchAddress("nCleanedJetsPt30_jecUp", &nCleanedJetsPt30_jecUp);
      inputTree[d]->SetBranchAddress("nCleanedJetsPt30_jecDn", &nCleanedJetsPt30_jecDn);
      inputTree[d]->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSF", &nCleanedJetsPt30BTagged_bTagSF);
      inputTree[d]->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSF_jecUp", &nCleanedJetsPt30BTagged_bTagSF_jecUp);
      inputTree[d]->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSF_jecDn", &nCleanedJetsPt30BTagged_bTagSF_jecDn);
      inputTree[d]->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSFUp", &nCleanedJetsPt30BTagged_bTagSFUp);
      inputTree[d]->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSFDn", &nCleanedJetsPt30BTagged_bTagSFDn);
      //inputTree[d]->SetBranchAddress("JetPt", &JetPt);
      //inputTree[d]->SetBranchAddress("JetEta", &JetEta);
      inputTree[d]->SetBranchAddress("JetPhi", &JetPhi);
      //inputTree[d]->SetBranchAddress("JetMass", &JetMass);
      inputTree[d]->SetBranchAddress("JetQGLikelihood", &JetQGLikelihood);
      inputTree[d]->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal", &p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal);
      inputTree[d]->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JECUp", &p_JJVBF_SIG_ghv1_1_JHUGen_JECUp);
      inputTree[d]->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JECDn", &p_JJVBF_SIG_ghv1_1_JHUGen_JECDn);
      inputTree[d]->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal", &p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal);
      inputTree[d]->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JECUp", &p_JJQCD_SIG_ghg2_1_JHUGen_JECUp);
      inputTree[d]->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JECDn", &p_JJQCD_SIG_ghg2_1_JHUGen_JECDn);
      inputTree[d]->SetBranchAddress("p_JVBF_SIG_ghv1_1_JHUGen_JECNominal", &p_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
      inputTree[d]->SetBranchAddress("p_JVBF_SIG_ghv1_1_JHUGen_JECUp", &p_JVBF_SIG_ghv1_1_JHUGen_JECUp);
      inputTree[d]->SetBranchAddress("p_JVBF_SIG_ghv1_1_JHUGen_JECDn", &p_JVBF_SIG_ghv1_1_JHUGen_JECDn);
      inputTree[d]->SetBranchAddress("pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal", &pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
      inputTree[d]->SetBranchAddress("pAux_JVBF_SIG_ghv1_1_JHUGen_JECUp", &pAux_JVBF_SIG_ghv1_1_JHUGen_JECUp);
      inputTree[d]->SetBranchAddress("pAux_JVBF_SIG_ghv1_1_JHUGen_JECDn", &pAux_JVBF_SIG_ghv1_1_JHUGen_JECDn);
      inputTree[d]->SetBranchAddress("p_JQCD_SIG_ghg2_1_JHUGen_JECNominal", &p_JQCD_SIG_ghg2_1_JHUGen_JECNominal);
      inputTree[d]->SetBranchAddress("p_JQCD_SIG_ghg2_1_JHUGen_JECUp", &p_JQCD_SIG_ghg2_1_JHUGen_JECUp);
      inputTree[d]->SetBranchAddress("p_JQCD_SIG_ghg2_1_JHUGen_JECDn", &p_JQCD_SIG_ghg2_1_JHUGen_JECDn);
      inputTree[d]->SetBranchAddress("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal", &p_HadWH_SIG_ghw1_1_JHUGen_JECNominal);
      inputTree[d]->SetBranchAddress("p_HadWH_SIG_ghw1_1_JHUGen_JECUp", &p_HadWH_SIG_ghw1_1_JHUGen_JECUp);
      inputTree[d]->SetBranchAddress("p_HadWH_SIG_ghw1_1_JHUGen_JECDn", &p_HadWH_SIG_ghw1_1_JHUGen_JECDn);
      inputTree[d]->SetBranchAddress("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal", &p_HadZH_SIG_ghz1_1_JHUGen_JECNominal);
      inputTree[d]->SetBranchAddress("p_HadZH_SIG_ghz1_1_JHUGen_JECUp", &p_HadZH_SIG_ghz1_1_JHUGen_JECUp);
      inputTree[d]->SetBranchAddress("p_HadZH_SIG_ghz1_1_JHUGen_JECDn", &p_HadZH_SIG_ghz1_1_JHUGen_JECDn);
      inputTree[d]->SetBranchAddress("GenHMass", &GenHMass);
      inputTree[d]->SetBranchAddress("GenZ1Phi", &GenZ1Phi);
      inputTree[d]->SetBranchAddress("GenZ2Phi", &GenZ2Phi);
      inputTree[d]->SetBranchAddress("GenZ1Flav", &GenZ1Flav);
      inputTree[d]->SetBranchAddress("GenZ2Flav", &GenZ2Flav);
      inputTree[d]->SetBranchAddress("GenLep1Eta", &GenLep1Eta);
      inputTree[d]->SetBranchAddress("GenLep1Phi", &GenLep1Phi);
      inputTree[d]->SetBranchAddress("GenLep1Id", &GenLep1Id);
      inputTree[d]->SetBranchAddress("GenLep2Eta", &GenLep2Eta);
      inputTree[d]->SetBranchAddress("GenLep2Phi", &GenLep2Phi);
      inputTree[d]->SetBranchAddress("GenLep2Id", &GenLep2Id);
      inputTree[d]->SetBranchAddress("GenLep3Eta", &GenLep3Eta);
      inputTree[d]->SetBranchAddress("GenLep3Phi", &GenLep3Phi);
      inputTree[d]->SetBranchAddress("GenLep3Id", &GenLep3Id);
      inputTree[d]->SetBranchAddress("GenLep4Eta", &GenLep4Eta);
      inputTree[d]->SetBranchAddress("GenLep4Phi", &GenLep4Phi);
      inputTree[d]->SetBranchAddress("GenLep4Id", &GenLep4Id);


      //---------- Process tree

      Int_t nMCEvtInSR = 0; 
      Int_t nMCEvtInSRInMassWindow = 0; 
      float count01jet = 0.;
      float countdijet = 0.;

      Long64_t entries = inputTree[d]->GetEntries();
      cout<<"Processing dataset "<<datasets[d]<<sMHPoint[mp]<<" ("<<entries<<" entries) ..."<<endl;

      for (Long64_t z=0; z<entries; ++z){

	if(DEBUG && z>1000) break;

	inputTree[d]->GetEntry(z);
      
	if( !(ZZsel>=90) ) continue;
	nMCEvtInSR++;
	if(ZZMass<m4lMin || ZZMass>m4lMax) continue;
	nMCEvtInSRInMassWindow++;

	Float_t kfactor = 1.;
	if(APPLYKFACTORS){
	  if(currentProcess==qqZZ)      kfactor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M;
	  else if(currentProcess==ggZZ) kfactor = KFactor_QCD_ggZZ_Nominal;
	}

	Double_t eventWeight = partialSampleWeight[d] * xsec * kfactor * (USEPUWEIGHT ? overallEventWeight : genHEPMCweight*dataMCWeight) ;


	//----- find subprocess if needed

	if(SPLITVHHADLEP){
	  if(initialCurrentProcess==WH){
	    if(genExtInfo>10) currentProcess = WH_lep;
	    else currentProcess = WH_had;
	  }else if(initialCurrentProcess==ZH){
	    if(genExtInfo>10) currentProcess = ZH_lep;
	    else currentProcess = ZH_had;
	  }
	}


	//----- find final state

	currentFinalState = -1;
	if(Z1Flav==-121){
	  if(Z2Flav==-121)
	    currentFinalState = fs4e;
	  else if(Z2Flav==-169)
	    currentFinalState = fs2e2mu;
	  else
	    cerr<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", Z2Flav="<<Z2Flav<<endl;
	}else if(Z1Flav==-169){
	  if(Z2Flav==-121)
	    currentFinalState = fs2mu2e;
	  else if(Z2Flav==-169)
	    currentFinalState = fs4mu;
	  else
	    cerr<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", Z2Flav="<<Z2Flav<<endl;
	}else{
	  cerr<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", Z1Flav="<<Z1Flav<<endl;
	}
	if(MERGE2E2MU && currentFinalState==fs2mu2e) currentFinalState = fs2e2mu;


	//----- find category

	for(int j=0; j<nCleanedJets; j++){
	  //jetPt[j] = JetPt->at(j);
	  //jetEta[j] = JetEta->at(j);
	  jetPhi[j] = JetPhi->at(j);
	  //jetMass[j] = JetMass->at(j);
	  jetQGL[j] = JetQGLikelihood->at(j);
	}
	Short_t varied_nCleanedJetsPt30[nSystShifts] = {nCleanedJetsPt30,nCleanedJetsPt30_jecUp,nCleanedJetsPt30_jecDn,nCleanedJetsPt30,nCleanedJetsPt30};
	Short_t varied_nCleanedJetsPt30BTagged[nSystShifts] = {nCleanedJetsPt30BTagged_bTagSF,nCleanedJetsPt30BTagged_bTagSF_jecUp,nCleanedJetsPt30BTagged_bTagSF_jecDn,nCleanedJetsPt30BTagged_bTagSFUp,nCleanedJetsPt30BTagged_bTagSFDn};
	Float_t varied_p_JJQCD_SIG_ghg2_1_JHUGen  [nSystShifts] = {p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal  ,p_JJQCD_SIG_ghg2_1_JHUGen_JECUp  ,p_JJQCD_SIG_ghg2_1_JHUGen_JECDn  ,p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal  ,p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal  };
	Float_t varied_p_JQCD_SIG_ghg2_1_JHUGen   [nSystShifts] = {p_JQCD_SIG_ghg2_1_JHUGen_JECNominal   ,p_JQCD_SIG_ghg2_1_JHUGen_JECUp   ,p_JQCD_SIG_ghg2_1_JHUGen_JECDn   ,p_JQCD_SIG_ghg2_1_JHUGen_JECNominal   ,p_JQCD_SIG_ghg2_1_JHUGen_JECNominal   };
	Float_t varied_p_JJVBF_SIG_ghv1_1_JHUGen  [nSystShifts] = {p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal  ,p_JJVBF_SIG_ghv1_1_JHUGen_JECUp  ,p_JJVBF_SIG_ghv1_1_JHUGen_JECDn  ,p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal  ,p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal  };
	Float_t varied_p_JVBF_SIG_ghv1_1_JHUGen   [nSystShifts] = {p_JVBF_SIG_ghv1_1_JHUGen_JECNominal   ,p_JVBF_SIG_ghv1_1_JHUGen_JECUp   ,p_JVBF_SIG_ghv1_1_JHUGen_JECDn   ,p_JVBF_SIG_ghv1_1_JHUGen_JECNominal   ,p_JVBF_SIG_ghv1_1_JHUGen_JECNominal   };
	Float_t varied_pAux_JVBF_SIG_ghv1_1_JHUGen[nSystShifts] = {pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,pAux_JVBF_SIG_ghv1_1_JHUGen_JECUp,pAux_JVBF_SIG_ghv1_1_JHUGen_JECDn,pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal};
	Float_t varied_p_HadWH_SIG_ghw1_1_JHUGen  [nSystShifts] = {p_HadWH_SIG_ghw1_1_JHUGen_JECNominal  ,p_HadWH_SIG_ghw1_1_JHUGen_JECUp  ,p_HadWH_SIG_ghw1_1_JHUGen_JECDn  ,p_HadWH_SIG_ghw1_1_JHUGen_JECNominal  ,p_HadWH_SIG_ghw1_1_JHUGen_JECNominal  };
	Float_t varied_p_HadZH_SIG_ghz1_1_JHUGen  [nSystShifts] = {p_HadZH_SIG_ghz1_1_JHUGen_JECNominal  ,p_HadZH_SIG_ghz1_1_JHUGen_JECUp  ,p_HadZH_SIG_ghz1_1_JHUGen_JECDn  ,p_HadZH_SIG_ghz1_1_JHUGen_JECNominal  ,p_HadZH_SIG_ghz1_1_JHUGen_JECNominal  };
	Float_t varied_PFMET[nSystShifts] = {PFMET,PFMET_jesUp,PFMET_jesDn,PFMET,PFMET};
	for(int sy=0; sy<nSystShifts; sy++){
	  currentCategory[sy] = categoryMor17(
	    nExtraLep,
	    nExtraZ,
	    varied_nCleanedJetsPt30[sy],
	    varied_nCleanedJetsPt30BTagged[sy],
	    jetQGL,
	    varied_p_JJQCD_SIG_ghg2_1_JHUGen[sy],
	    varied_p_JQCD_SIG_ghg2_1_JHUGen[sy],
	    varied_p_JJVBF_SIG_ghv1_1_JHUGen[sy],
	    varied_p_JVBF_SIG_ghv1_1_JHUGen[sy],
	    varied_pAux_JVBF_SIG_ghv1_1_JHUGen[sy],
	    varied_p_HadWH_SIG_ghw1_1_JHUGen[sy],
	    varied_p_HadZH_SIG_ghz1_1_JHUGen[sy],
	    jetPhi,
	    ZZMass,
	    varied_PFMET[sy],
	    USEVHMETTAG,
	    false
	    );
	}
	if(nCleanedJetsPt30>=2) countdijet += eventWeight; else count01jet += eventWeight;


	/* //----- check if this is resonant signal (ie. H->4l with correct lepton choice)

	Short_t GenHLepId[4] = {GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id};
	Float_t GenHLepEta[4] = {GenLep1Eta,GenLep2Eta,GenLep3Eta,GenLep4Eta};
	Float_t GenHLepPhi[4] = {GenLep1Phi,GenLep2Phi,GenLep3Phi,GenLep4Phi};

	Int_t nGenHLep = 0;
	Int_t nRecoLepMatchedToGenHLep[4] = {0,0,0,0};
	Int_t nGenHLepMatchedToCandLep[4] = {0,0,0,0};
	for(Int_t iGenHLep=0; iGenHLep<4; iGenHLep++){
	  if(abs(GenHLepId[iGenHLep])==11 || abs(GenHLepId[iGenHLep])==13){
	    nGenHLep++;
	    for(Int_t iCandLep=0; iCandLep<4; iCandLep++){
	      if(deltaR(GenHLepEta[iGenHLep],GenHLepPhi[iGenHLep],LepEta->at(iCandLep),LepPhi->at(iCandLep)) < 0.1){
		nRecoLepMatchedToGenHLep[iGenHLep]++;
		nGenHLepMatchedToCandLep[iCandLep]++;
	      }
	    }
	    for(Int_t iExtraLep=0; iExtraLep<nExtraLep; iExtraLep++){
	      if(deltaR(GenHLepEta[iGenHLep],GenHLepPhi[iGenHLep],ExtraLepEta->at(iExtraLep),ExtraLepPhi->at(iExtraLep)) < 0.1){
		nRecoLepMatchedToGenHLep[iGenHLep]++;
	      }
	    }
	  }
	}
	Bool_t foundMatchingAmbiguity = false;
	for(Int_t iGenHLep=0; iGenHLep<4; iGenHLep++) if(nRecoLepMatchedToGenHLep[iGenHLep]>1){ foundMatchingAmbiguity = true; break; }
	for(Int_t iCandLep=0; iCandLep<4; iCandLep++) if(nGenHLepMatchedToCandLep[iCandLep]>1){ foundMatchingAmbiguity = true; break; }
	Int_t nMatches = 0;
	for(Int_t iGenHLep=0; iGenHLep<4; iGenHLep++) if(nRecoLepMatchedToGenHLep[iGenHLep]==1) nMatches++;

	if(nGenHLep==4 && !foundMatchingAmbiguity && nMatches==4) 
	  currentResStatus = resonant;
	else 
	  currentResStatus = nonresonant;
	//*/


	//----- fill yield histogram
	for(int sy=0; sy<nSystShifts; sy++)
	  hYield[currentProcess][mp][currentFinalState][currentCategory[sy]][nRS][sy]->Fill(0.5,eventWeight);
	//hYield[currentProcess][mp][currentFinalState][currentCategory[sy]][currentResStatus][sy]->Fill(0.5,eventWeight);


      } // end for entries

      //cout<<"nb of MC events in SR : "<<nMCEvtInSR<<endl;
      //cout<<"nb of MC events for requested m4l window : "<<nMCEvtInSRInMassWindow<<endl;
      //cout<<"expected yield with <2jets :  "<<count01jet<<endl;
      //cout<<"expected yield with >=2jets : "<<countdijet<<endl;
    
    } // end for mHPoints

  } // end for datasets


  //---------- Fill 'inclusive' counters
  for(int pr=0; pr<nProcesses; pr++)
    for(int mp=0; mp<nMHPoints; mp++)
      if(hasMHPoint[pr][mp])
	for(int fs=0; fs<nFinalStates; fs++)
	  for(int cat=0; cat<nCategories; cat++)
	    //for(int rs=0; rs<nRS; rs++)
	      for(int sy=0; sy<nSystShifts; sy++)
		hYield[pr][mp][nFinalStates][cat][nRS][sy]->Add(hYield[pr][mp][fs][cat][nRS][sy]);
  for(int pr=0; pr<nProcesses; pr++)
    for(int mp=0; mp<nMHPoints; mp++)
      if(hasMHPoint[pr][mp])
	for(int fs=0; fs<nFinalStates+1; fs++)
	  for(int cat=0; cat<nCategories; cat++)
	    //for(int rs=0; rs<nRS; rs++)
	      for(int sy=0; sy<nSystShifts; sy++)
		hYield[pr][mp][fs][nCategories][nRS][sy]->Add(hYield[pr][mp][fs][cat][nRS][sy]);

  //---------- Write yield arrays to a ROOT file
  TFile* fOutYields = new TFile(Form("yields_%iTeV_m4l%.1f-%.1f_%.3ffb-1%s.root",(int)sqrts,m4lMin,m4lMax,lumi,(MERGE2E2MU?"_m":"")),"recreate");
  fOutYields->cd();
  for(int pr=0; pr<nProcesses; pr++){
    for(int mp=0; mp<nMHPoints; mp++){
      if(!hasMHPoint[pr][mp] || (JUST125&&sMHPoint[mp]!="125"&&sMHPoint[mp]!="")) continue;
      for(int fs=0; fs<nFinalStates+1; fs++){
	if(MERGE2E2MU && fs==fs2mu2e) continue;
	for(int cat=0; cat<nCategories+1; cat++){
	  //for(int rs=0; rs<nRS+1; rs++)
	  for(int sy=0; sy<nSystShifts; sy++){
	    hYield[pr][mp][fs][cat][nRS][sy]->Write(hYield[pr][mp][fs][cat][nRS][sy]->GetName());
	    delete hYield[pr][mp][fs][cat][nRS][sy];
	  }
	}
      }
    }
  }
  fOutYields->Close();
  delete fOutYields; 

}


TGraphAsymmErrors* gr_FRmu_EB = 0;
TGraphAsymmErrors* gr_FRmu_EE = 0;
TGraphErrors* gr_FRel_EB = 0;
TGraphErrors* gr_FRel_EE = 0;

Float_t fakeRate13TeV(Float_t LepPt, Float_t LepEta, Int_t LepID, string syst="nominal") {
  Float_t myLepPt = LepPt>=80. ? 79. : LepPt;
  Int_t   myLepID = abs(LepID);

  int bin = 0;
  if(myLepPt > 5 && myLepPt<=7) bin = 0;
  else if(myLepPt > 7 && myLepPt<=10) bin = 1;
  else if(myLepPt > 10 && myLepPt<=20) bin = 2;
  else if(myLepPt > 20 && myLepPt<=30) bin = 3;
  else if(myLepPt > 30 && myLepPt<=40) bin = 4;
  else if(myLepPt > 40 && myLepPt<=50) bin = 5;
  else if(myLepPt > 50 && myLepPt<=80) bin = 6;
  if(fabs(myLepID)==11) bin = bin-1; // there is no [5, 7] bin in the electron fake rate      

  if(myLepID==11){
    if(fabs(LepEta)<1.479)
      return (gr_FRel_EB->GetY())[bin] + (syst=="up" ? (gr_FRel_EB->GetEY())[bin] : syst=="dn" ? -(gr_FRel_EB->GetEY())[bin] : 0.);
    else
      return (gr_FRel_EE->GetY())[bin] + (syst=="up" ? (gr_FRel_EE->GetEY())[bin] : syst=="dn" ? -(gr_FRel_EE->GetEY())[bin] : 0.);
  }else if(myLepID==13){
    if(fabs(LepEta)<1.2)
      return (gr_FRmu_EB->GetY())[bin] + (syst=="up" ? (gr_FRmu_EB->GetEYhigh())[bin] : syst=="dn" ? -(gr_FRmu_EB->GetEYlow())[bin] : 0.);
    else
      return (gr_FRmu_EE->GetY())[bin] + (syst=="up" ? (gr_FRmu_EE->GetEYhigh())[bin] : syst=="dn" ? -(gr_FRmu_EE->GetEYlow())[bin] : 0.);
  }else{
    cout<<"ERROR! wrong lepton ID : "<<myLepID<<endl;
    return 0.;
  }
}


void computeYieldsZPlusXSS(string inputPathData, string inputFileFakeRates, string outputDirectory, double lumi)
{

  cout<<"Preparing Z+X yields (SS method) from "<<inputPathData<<endl;

  TFile* fFakeRates = TFile::Open(inputFileFakeRates.c_str());
  gr_FRmu_EB = (TGraphAsymmErrors*)fFakeRates->Get("FR_SS_muon_EB");
  gr_FRmu_EE = (TGraphAsymmErrors*)fFakeRates->Get("FR_SS_muon_EE");
  gr_FRel_EB = (TGraphErrors*)fFakeRates->Get("FR_SS_electron_EB");
  gr_FRel_EE = (TGraphErrors*)fFakeRates->Get("FR_SS_electron_EE");

  TH1::SetDefaultSumw2(true);

  Int_t nRun;
  Long64_t nEvent;
  Int_t nLumi;
  Float_t PFMET;
  Int_t CRflag;
  Short_t ZZsel;
  Float_t ZZMass;
  Float_t ZZPt;
  Short_t Z1Flav;
  Short_t Z2Flav;
  Float_t p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal;
  Float_t p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal;
  Float_t p_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
  Float_t pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
  Float_t p_JQCD_SIG_ghg2_1_JHUGen_JECNominal;
  Float_t p_HadWH_SIG_ghw1_1_JHUGen_JECNominal;
  Float_t p_HadZH_SIG_ghz1_1_JHUGen_JECNominal;
  vector<Float_t> *LepPt = 0;
  vector<Float_t> *LepEta = 0;
  vector<Float_t> *LepLepId = 0;
  Short_t nExtraLep;
  Short_t nExtraZ;
  Short_t nCleanedJetsPt30;
  Short_t nCleanedJetsPt30BTagged;
  vector<Float_t> *JetPhi = 0;
  vector<Float_t> *JetQGLikelihood = 0;
  Float_t jetPhi[99];
  Float_t jetQGL[99]; 

  Float_t expectedYieldSR[nFinalStates+1][nCategories+1];
  Float_t expectedYieldSR_up[nFinalStates+1][nCategories+1];
  Float_t expectedYieldSR_dn[nFinalStates+1][nCategories+1];
  Int_t NumberOfEventsCR[nFinalStates+1][nCategories+1];
  for(int fs=0; fs<nFinalStates+1; fs++){
    for(int cat=0; cat<nCategories+1; cat++){
      expectedYieldSR[fs][cat] = 0.;
      expectedYieldSR_up[fs][cat] = 0.;
      expectedYieldSR_dn[fs][cat] = 0.;
      NumberOfEventsCR[fs][cat] = 0.;
    }
  }

  int currentFinalState;
  int currentCategory;

  TFile* dataFile = TFile::Open((inputPathData+"AllData/ZZ4lAnalysis.root").c_str());
  TTree* mytree = (TTree*)dataFile->Get("CRZLLTree/candTree");

  mytree->SetBranchAddress("RunNumber", &nRun);
  mytree->SetBranchAddress("EventNumber", &nEvent);
  mytree->SetBranchAddress("LumiNumber", &nLumi);
  mytree->SetBranchAddress("PFMET", &PFMET);
  mytree->SetBranchAddress("CRflag", &CRflag);
  mytree->SetBranchAddress("ZZsel", &ZZsel);
  mytree->SetBranchAddress("ZZMass", &ZZMass);
  mytree->SetBranchAddress("ZZPt", &ZZPt);
  mytree->SetBranchAddress("Z1Flav", &Z1Flav);
  mytree->SetBranchAddress("Z2Flav", &Z2Flav);
  mytree->SetBranchAddress("LepPt", &LepPt);
  mytree->SetBranchAddress("LepEta", &LepEta);
  mytree->SetBranchAddress("LepLepId", &LepLepId);
  mytree->SetBranchAddress("nExtraLep", &nExtraLep);
  mytree->SetBranchAddress("nExtraZ", &nExtraZ);
  mytree->SetBranchAddress("nCleanedJetsPt30", &nCleanedJetsPt30);
  mytree->SetBranchAddress("nCleanedJetsPt30BTagged", &nCleanedJetsPt30BTagged);
  mytree->SetBranchAddress("JetPhi", &JetPhi);
  mytree->SetBranchAddress("JetQGLikelihood", &JetQGLikelihood);
  mytree->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal", &p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal);
  mytree->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal", &p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal);
  mytree->SetBranchAddress("p_JVBF_SIG_ghv1_1_JHUGen_JECNominal", &p_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
  mytree->SetBranchAddress("pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal", &pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
  mytree->SetBranchAddress("p_JQCD_SIG_ghg2_1_JHUGen_JECNominal", &p_JQCD_SIG_ghg2_1_JHUGen_JECNominal);
  mytree->SetBranchAddress("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal", &p_HadWH_SIG_ghw1_1_JHUGen_JECNominal);
  mytree->SetBranchAddress("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal", &p_HadZH_SIG_ghz1_1_JHUGen_JECNominal);
  
  
  //---------- Process tree
  
  Long64_t entries = mytree->GetEntries();
  for (Long64_t z=0; z<entries; ++z){
    
    mytree->GetEntry(z);
    
    if(!CRflag) continue;
    if(!test_bit(CRflag,CRZLLss)) continue;
    
    //----- find final state
    currentFinalState = -1;
    if(Z1Flav==-121){
      if(Z2Flav==+121)
	currentFinalState = fs4e;
      else if(Z2Flav==+169)
	currentFinalState = fs2e2mu;
      else
	cerr<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", Z2Flav="<<Z2Flav<<endl;
    }else if(Z1Flav==-169){
      if(Z2Flav==+121)
	currentFinalState = fs2mu2e;
      else if(Z2Flav==+169)
	currentFinalState = fs4mu;
      else
	cerr<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", Z2Flav="<<Z2Flav<<endl;
    }else{
      cerr<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", Z1Flav="<<Z1Flav<<endl;
    }

    //----- find category
    for(int j=0; j<nCleanedJetsPt30; j++){
      jetPhi[j] = JetPhi->at(j);
      jetQGL[j] = JetQGLikelihood->at(j);
    }
    currentCategory = categoryMor17(
       nExtraLep,
       nExtraZ,
       nCleanedJetsPt30,
       nCleanedJetsPt30BTagged,
       jetQGL,
       p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
       p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
       p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
       p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
       pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
       p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
       p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
       jetPhi,
       ZZMass,
       PFMET,
       USEVHMETTAG,
       false
    );

    //----- update counters
    Float_t yieldSR = fsROSSS[currentFinalState] * fakeRate13TeV(LepPt->at(2),LepEta->at(2),LepLepId->at(2)) * fakeRate13TeV(LepPt->at(3),LepEta->at(3),LepLepId->at(3));
    Float_t yieldSR_up = fsROSSS[currentFinalState] * fakeRate13TeV(LepPt->at(2),LepEta->at(2),LepLepId->at(2),"up") * fakeRate13TeV(LepPt->at(3),LepEta->at(3),LepLepId->at(3),"up");
    Float_t yieldSR_dn = fsROSSS[currentFinalState] * fakeRate13TeV(LepPt->at(2),LepEta->at(2),LepLepId->at(2),"dn") * fakeRate13TeV(LepPt->at(3),LepEta->at(3),LepLepId->at(3),"dn");
    expectedYieldSR[currentFinalState][currentCategory] += yieldSR;
    expectedYieldSR_up[currentFinalState][currentCategory] += yieldSR_up;
    expectedYieldSR_dn[currentFinalState][currentCategory] += yieldSR_dn;
    NumberOfEventsCR[currentFinalState][currentCategory]++;

  }


  //---------- Fill 'inclusive' counters
  for(int fs=0; fs<nFinalStates; fs++){
    for(int cat=0; cat<nCategories; cat++){
      expectedYieldSR[nFinalStates][cat] += expectedYieldSR[fs][cat];
      expectedYieldSR_up[nFinalStates][cat] += expectedYieldSR_up[fs][cat];
      expectedYieldSR_dn[nFinalStates][cat] += expectedYieldSR_dn[fs][cat];
      NumberOfEventsCR[nFinalStates][cat] += NumberOfEventsCR[fs][cat];
    }
  }
  for(int fs=0; fs<nFinalStates+1; fs++){
    for(int cat=0; cat<nCategories; cat++){
      expectedYieldSR[fs][nCategories] += expectedYieldSR[fs][cat];
      expectedYieldSR_up[fs][nCategories] += expectedYieldSR_up[fs][cat];
      expectedYieldSR_dn[fs][nCategories] += expectedYieldSR_dn[fs][cat];
      NumberOfEventsCR[fs][nCategories] += NumberOfEventsCR[fs][cat];
    }
  }

  //---------- Refill global array 
  for(int cat=0; cat<nCategories+1; cat++){
    normSSFullRange4e[cat] = expectedYieldSR[fs4e][cat];
    normSSFullRange4mu[cat] = expectedYieldSR[fs4mu][cat];
    normSSFullRange2e2mu[cat] = expectedYieldSR[fs2e2mu][cat]+expectedYieldSR[fs2mu2e][cat];
  }

  //---------- Compute uncertainties
  Float_t StatUnc[nFinalStates+1][nCategories+1];
  Float_t SystUnc_up[nFinalStates+1][nCategories+1]; // from propagation of uncertainties of fake rates
  Float_t SystUnc_dn[nFinalStates+1][nCategories+1];
  Float_t BkgdUnc[nFinalStates+1][nCategories+1]; // from bkgd composition variation in OS method
  Float_t TotUnc_up[nFinalStates+1][nCategories+1]; // quadrature sum of the 3 components
  Float_t TotUnc_dn[nFinalStates+1][nCategories+1];
  for(int fs=0; fs<nFinalStates+1; fs++){
    for(int cat=0; cat<nCategories+1; cat++){
      StatUnc[fs][cat] = expectedYieldSR[fs][cat]/sqrt(NumberOfEventsCR[fs][cat]);
      SystUnc_up[fs][cat] = expectedYieldSR_up[fs][cat]-expectedYieldSR[fs][cat];
      SystUnc_dn[fs][cat] = expectedYieldSR[fs][cat]-expectedYieldSR_dn[fs][cat];
      BkgdUnc[fs][cat] = expectedYieldSR[fs][cat]*bkgdCompUncZpX[fs];
      TotUnc_up[fs][cat] = sqrt(StatUnc[fs][cat]*StatUnc[fs][cat] + SystUnc_up[fs][cat]*SystUnc_up[fs][cat] + BkgdUnc[fs][cat]*BkgdUnc[fs][cat]);
      TotUnc_dn[fs][cat] = sqrt(StatUnc[fs][cat]*StatUnc[fs][cat] + SystUnc_dn[fs][cat]*SystUnc_dn[fs][cat] + BkgdUnc[fs][cat]*BkgdUnc[fs][cat]);
    }
  }
  //properly combine 2e2mu and 2mu2e as TETM
  Float_t expectedYieldSR_TETM[nCategories+1];
  Float_t StatUnc_TETM[nCategories+1];
  Float_t SystUnc_up_TETM[nCategories+1];
  Float_t SystUnc_dn_TETM[nCategories+1];
  Float_t BkgdUnc_TETM[nCategories+1];
  Float_t TotUnc_up_TETM[nCategories+1];
  Float_t TotUnc_dn_TETM[nCategories+1];
  for(int cat=0; cat<nCategories+1; cat++){
    expectedYieldSR_TETM[cat] = expectedYieldSR[fs2e2mu][cat] + expectedYieldSR[fs2mu2e][cat];
    StatUnc_TETM[cat] = sqrt(StatUnc[fs2e2mu][cat]*StatUnc[fs2e2mu][cat] + StatUnc[fs2mu2e][cat]*StatUnc[fs2mu2e][cat]);
    SystUnc_up_TETM[cat] = sqrt(SystUnc_up[fs2e2mu][cat]*SystUnc_up[fs2e2mu][cat] + SystUnc_up[fs2mu2e][cat]*SystUnc_up[fs2mu2e][cat]);
    SystUnc_dn_TETM[cat] = sqrt(SystUnc_dn[fs2e2mu][cat]*SystUnc_dn[fs2e2mu][cat] + SystUnc_dn[fs2mu2e][cat]*SystUnc_dn[fs2mu2e][cat]);
    BkgdUnc_TETM[cat] = BkgdUnc[fs2e2mu][cat] + BkgdUnc[fs2mu2e][cat];
    TotUnc_up_TETM[cat] = sqrt(StatUnc_TETM[cat]*StatUnc_TETM[cat] + SystUnc_up_TETM[cat]*SystUnc_up_TETM[cat] + BkgdUnc_TETM[cat]*BkgdUnc_TETM[cat]);
    TotUnc_dn_TETM[cat] = sqrt(StatUnc_TETM[cat]*StatUnc_TETM[cat] + SystUnc_dn_TETM[cat]*SystUnc_dn_TETM[cat] + BkgdUnc_TETM[cat]*BkgdUnc_TETM[cat]);
  }

  /* //---------- Print Z+X expected yields
  for(int fs=0; fs<nFinalStates; fs++){
    cout<<sFinalState[fs]<<" : "
        <<expectedYieldSR[fs][nCategories]
        <<" +/- "<<expectedYieldSR[fs][nCategories]/sqrt(NumberOfEventsCR[fs][nCategories])
	<<" (stat., evt: "<<NumberOfEventsCR[fs][nCategories]<<")"
        <<endl;
  }
  cout<<"Total: "<<expectedYieldSR[nFinalStates][nCategories]<<endl;
  cout<<endl;
  //*/

  if(PRINTLATEXTABLES){
    //tables in LaTeX format for the Z+X chapter of the AN
    int w1 = 18;
    int w2 = 6;
    int nD = 3;
    for(int fs=0; fs<nFinalStates; fs++){
      cout<<"Z+X table for final state "<<sFinalState[fs]<<endl<<endl;
      cout<<"\\hline"<<endl<<"\\hline"<<endl;
      cout<<"Category & events in CR & exp. $N^{\\rm Z+X}$ in SR & (stat.) & (syst.) & (bkgd. compo.) & Total unc.\\\\"<<endl;
      cout<<"\\hline"<<endl<<"\\hline"<<endl;
      for(int cat=0; cat<=nCategories; cat++){
	if(cat==nCategories) cout<<"\\hline"<<endl;
	cout<<fixWidth(sCatLatexLong[cat],w1,1)<<" & "
	    <<fixWidth(string(Form("%i",NumberOfEventsCR[fs][cat])),w2,0)<<" & "
	    <<fixWidth(sround(expectedYieldSR[fs][cat],nD),w2,0)<<" & "
	    <<"$\\pm"<<fixWidth(sround(StatUnc[fs][cat],nD),w2,0)<<"$ & "
	    <<"$^{+"<<fixWidth(sround(SystUnc_up[fs][cat],nD),w2,0)<<"}_{-"<<fixWidth(sround(SystUnc_dn[fs][cat],nD),w2,0)<<"}$ & "
	    <<"$\\pm"<<fixWidth(sround(BkgdUnc[fs][cat],nD),w2,0)<<"$ & "
	    <<"$^{+"<<fixWidth(sround(TotUnc_up[fs][cat],nD),w2,0)<<"}_{-"<<fixWidth(sround(TotUnc_dn[fs][cat],nD),w2,0)<<"}$ \\\\"
	    <<endl<<"\\hline"<<endl;
      }
      cout<<"\\hline"<<endl;
      cout<<endl;
    }
    { // table combining 2e2mu and 2mu2e
      cout<<"Z+X table for sum of 2e2mu and 2mu2e"<<endl<<endl;
      cout<<"\\hline"<<endl<<"\\hline"<<endl;
      cout<<"Category & exp. $N^{\\rm Z+X}$ in SR & (stat.) & (syst.) & (bkgd. compo.) & Total unc.\\\\"<<endl;
      cout<<"\\hline"<<endl<<"\\hline"<<endl;
      for(int cat=0; cat<=nCategories; cat++){
	if(cat==nCategories) cout<<"\\hline"<<endl;
	cout<<fixWidth(sCatLatexLong[cat],w1,1)<<" & "
	    <<fixWidth(sround(expectedYieldSR_TETM[cat],nD),w2,0)<<" & "
	    <<"$\\pm"<<fixWidth(sround(StatUnc_TETM[cat],nD),w2,0)<<"$ & "
	    <<"$^{+"<<fixWidth(sround(SystUnc_up_TETM[cat],nD),w2,0)<<"}_{-"<<fixWidth(sround(SystUnc_dn_TETM[cat],nD),w2,0)<<"}$ & "
	    <<"$\\pm"<<fixWidth(sround(BkgdUnc_TETM[cat],nD),w2,0)<<"$ & "
	    <<"$^{+"<<fixWidth(sround(TotUnc_up_TETM[cat],nD),w2,0)<<"}_{-"<<fixWidth(sround(TotUnc_dn_TETM[cat],nD),w2,0)<<"}$ \\\\"
	    <<endl<<"\\hline"<<endl;
      }
      cout<<"\\hline"<<endl;
      cout<<endl;
    }
  }

  // Z+X systematic uncertainties for datacards
  {
    ofstream outFileSyst[nFinalStates];
    
    for(int fs=0; fs<nFinalStates; fs++){    
      if(fs==fs2mu2e) continue;
      
      outFileSyst[fs].open((outputDirectory+"/partof_systematics_13TeV_"+sFinalState[fs]+".yaml").c_str());

      outFileSyst[fs]<<"CMS_zjets_bkgdcompo:"<<endl
		     <<"    UnTagged: &zjets_bkgdcompo_untagged"<<endl
		     <<"        type : lnN"<<endl
		     <<"        zjets: "<<1+bkgdCompUncZpX[fs]<<endl;
      for(int cat=1; cat<nCategories; cat++)
	outFileSyst[fs]<<"    "<<sCategoryForSyst[cat]<<":"<<endl<<"        <<: *zjets_bkgdcompo_untagged"<<endl;
      outFileSyst[fs]<<endl;
      
      outFileSyst[fs]<<"CMS_zz"+sFinalState[fs]+"_zjets:"<<endl;
      for(int cat=0; cat<nCategories; cat++){
	outFileSyst[fs]<<"    "<<sCategoryForSyst[cat]<<":"<<endl
		       <<"        type : lnN"<<endl
		       <<"        zjets: ";
	if(fs!=fs2e2mu){
	  outFileSyst[fs]<<sround( 1. + sqrt( TotUnc_up[fs][cat]*TotUnc_up[fs][cat] - BkgdUnc[fs][cat]*BkgdUnc[fs][cat] ) / expectedYieldSR[fs][cat] ,3)<<"/"
			 <<sround( 1. - sqrt( TotUnc_dn[fs][cat]*TotUnc_dn[fs][cat] - BkgdUnc[fs][cat]*BkgdUnc[fs][cat] ) / expectedYieldSR[fs][cat] ,3)<<endl;
	}else{
	  outFileSyst[fs]<<sround( 1. + sqrt( TotUnc_up_TETM[cat]*TotUnc_up_TETM[cat] - BkgdUnc_TETM[cat]*BkgdUnc_TETM[cat] ) / expectedYieldSR_TETM[cat] ,3)<<"/"
			 <<sround( 1. - sqrt( TotUnc_dn_TETM[cat]*TotUnc_dn_TETM[cat] - BkgdUnc_TETM[cat]*BkgdUnc_TETM[cat] ) / expectedYieldSR_TETM[cat] ,3)<<endl;
	}
      }

      outFileSyst[fs].close();
    }
  }

}


void getZPlusXYields_Run2CombinedShape_InCateg(Float_t* yieldZPX4mu, Float_t* yieldZPX4e, Float_t* yieldZPX2e2mu, Float_t* yieldZPX4l, Int_t m4lMin, Int_t m4lMax) {

  /*// For Moriond'16: Z+X shapes sent by Pedja on March 1st 2016 (take the same for all categ)
  TF1 *f4eComb = new TF1("f4eComb", "landau(0)*(1 + exp( pol1(3))) + [5]*(TMath::Landau(x, [6], [7]))", 70, 3000);
  TF1 *f4muComb = new TF1("f4muComb","landau(0)",70,3000);
  TF1 *f2e2muComb = new TF1("f2e2muComb","landau(0)",70,3000);
  f4eComb->SetParameters(4.404e-05,151.2,36.6,7.06,-0.00497,0.01446,157.3,26.00);
  f4muComb->SetParameters(0.04276,134.6,24.4);
  f2e2muComb->SetParameters(0.04130,144.5,25.3);
  //*/

  //*// For ICHEP'16: Z+X shapes sent by Pedja on July 27th 2016 (take the same for all categ)
  TF1 *f4eComb    = new TF1("f4eComb"   ,"TMath::Landau(x, 141.9, 21.3)", 70, 3000);
  TF1 *f4muComb   = new TF1("f4muComb"  ,"TMath::Landau(x, 130.4, 15.6)", 70, 3000);
  TF1 *f2e2muComb = new TF1("f2e2muComb","0.45*TMath::Landau(x, 131.1, 18.1) + 0.55*TMath::Landau(x, 133.8, 18.9)", 70, 3000);
  //*/

  //----- compute normalization of the subrange of interest

  for(int cat=0; cat<nCategories+1; cat++){
    yieldZPX4mu  [cat] = normSSFullRange4mu  [cat] * (RESCALEZPLUSXTOCOMBINEDINCLUSIVE?(normCombFullRange4mu  /normSSFullRange4mu  [nCategories]):1.) * f4muComb  ->Integral(m4lMin,m4lMax) / f4muComb  ->Integral(70,3000);
    yieldZPX4e   [cat] = normSSFullRange4e   [cat] * (RESCALEZPLUSXTOCOMBINEDINCLUSIVE?(normCombFullRange4e   /normSSFullRange4e   [nCategories]):1.) * f4eComb   ->Integral(m4lMin,m4lMax) / f4eComb   ->Integral(70,3000);
    yieldZPX2e2mu[cat] = normSSFullRange2e2mu[cat] * (RESCALEZPLUSXTOCOMBINEDINCLUSIVE?(normCombFullRange2e2mu/normSSFullRange2e2mu[nCategories]):1.) * f2e2muComb->Integral(m4lMin,m4lMax) / f2e2muComb->Integral(70,3000);
    yieldZPX4l   [cat] = yieldZPX4mu[cat] + yieldZPX4e[cat] + yieldZPX2e2mu[cat];
    //cout<<sCategory[cat]<<"  "<<normSSFullRange4mu[cat]<<"  "<<normSSFullRange4e[cat]<<"  "<<normSSFullRange2e2mu[cat]<<endl;
    //cout<<sCategory[cat]<<"  "<<yieldZPX4mu[cat]<<"  "<<yieldZPX4e[cat]<<"  "<<yieldZPX2e2mu[cat]<<"  "<<yieldZPX4l[cat]<<endl;
  }

}


void dataEventCounts(string inputPathData, double m4lMin, double m4lMax)
{

  cout<<"Counting SR event in data in ["<<m4lMin<<","<<m4lMax<<"]GeV from "<<inputPathData<<endl;

  Int_t nRun;
  Long64_t nEvent;
  Int_t nLumi;
  Float_t PFMET;
  Short_t ZZsel;
  Float_t ZZMass;
  Float_t ZZPt;
  Short_t Z1Flav;
  Short_t Z2Flav;
  Float_t p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal;
  Float_t p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal;
  Float_t p_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
  Float_t pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
  Float_t p_JQCD_SIG_ghg2_1_JHUGen_JECNominal;
  Float_t p_HadWH_SIG_ghw1_1_JHUGen_JECNominal;
  Float_t p_HadZH_SIG_ghz1_1_JHUGen_JECNominal;
  Short_t nExtraLep;
  Short_t nExtraZ;
  Short_t nCleanedJetsPt30;
  Short_t nCleanedJetsPt30BTagged;
  vector<Float_t> *JetPhi = 0;
  vector<Float_t> *JetQGLikelihood = 0;
  Float_t jetPhi[99];
  Float_t jetQGL[99]; 

  Int_t NumberOfEventsSR[nFinalStates+1][nCategories+1];
  for(int fs=0; fs<nFinalStates+1; fs++)
    for(int cat=0; cat<nCategories+1; cat++)
      NumberOfEventsSR[fs][cat] = 0.;

  int currentFinalState;
  int currentCategory;

  TFile* dataFile = TFile::Open((inputPathData+"AllData/ZZ4lAnalysis.root").c_str());
  TTree* mytree = (TTree*)dataFile->Get("ZZTree/candTree");

  mytree->SetBranchAddress("RunNumber", &nRun);
  mytree->SetBranchAddress("EventNumber", &nEvent);
  mytree->SetBranchAddress("LumiNumber", &nLumi);
  mytree->SetBranchAddress("PFMET", &PFMET);
  mytree->SetBranchAddress("ZZsel", &ZZsel);
  mytree->SetBranchAddress("ZZMass", &ZZMass);
  mytree->SetBranchAddress("ZZPt", &ZZPt);
  mytree->SetBranchAddress("Z1Flav", &Z1Flav);
  mytree->SetBranchAddress("Z2Flav", &Z2Flav);
  mytree->SetBranchAddress("nExtraLep", &nExtraLep);
  mytree->SetBranchAddress("nExtraZ", &nExtraZ);
  mytree->SetBranchAddress("nCleanedJetsPt30", &nCleanedJetsPt30);
  mytree->SetBranchAddress("nCleanedJetsPt30BTagged", &nCleanedJetsPt30BTagged);
  mytree->SetBranchAddress("JetPhi", &JetPhi);
  mytree->SetBranchAddress("JetQGLikelihood", &JetQGLikelihood);
  mytree->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal", &p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal);
  mytree->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal", &p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal);
  mytree->SetBranchAddress("p_JVBF_SIG_ghv1_1_JHUGen_JECNominal", &p_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
  mytree->SetBranchAddress("pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal", &pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
  mytree->SetBranchAddress("p_JQCD_SIG_ghg2_1_JHUGen_JECNominal", &p_JQCD_SIG_ghg2_1_JHUGen_JECNominal);
  mytree->SetBranchAddress("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal", &p_HadWH_SIG_ghw1_1_JHUGen_JECNominal);
  mytree->SetBranchAddress("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal", &p_HadZH_SIG_ghz1_1_JHUGen_JECNominal);
  
  
  //---------- Process tree
  
  Long64_t entries = mytree->GetEntries();
  for (Long64_t z=0; z<entries; ++z){
    
    mytree->GetEntry(z);

    if( !(ZZsel>=90) ) continue;
    if(ZZMass<m4lMin || ZZMass>m4lMax) continue;
    
    //----- find final state
    currentFinalState = -1;
    if(Z1Flav==-121){
      if(Z2Flav==-121)
	currentFinalState = fs4e;
      else if(Z2Flav==-169)
	currentFinalState = fs2e2mu;
      else
	cerr<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", Z2Flav="<<Z2Flav<<endl;
    }else if(Z1Flav==-169){
      if(Z2Flav==-121)
	currentFinalState = fs2mu2e;
      else if(Z2Flav==-169)
	currentFinalState = fs4mu;
      else
	cerr<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", Z2Flav="<<Z2Flav<<endl;
    }else{
      cerr<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", Z1Flav="<<Z1Flav<<endl;
    }
    if(MERGE2E2MU && currentFinalState==fs2mu2e) currentFinalState=fs2e2mu;

    //----- find category
    for(int j=0; j<nCleanedJetsPt30; j++){
      jetPhi[j] = JetPhi->at(j);
      jetQGL[j] = JetQGLikelihood->at(j);
    }
    currentCategory = categoryMor17(
       nExtraLep,
       nExtraZ,
       nCleanedJetsPt30,
       nCleanedJetsPt30BTagged,
       jetQGL,
       p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
       p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
       p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
       p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
       pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
       p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
       p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
       jetPhi,
       ZZMass,
       PFMET,
       USEVHMETTAG,
       false
    );

    //----- update counter
    NumberOfEventsSR[currentFinalState][currentCategory]++;

  }

  //---------- Fill 'inclusive' counters
  for(int fs=0; fs<nFinalStates; fs++)
    for(int cat=0; cat<nCategories; cat++)
      NumberOfEventsSR[nFinalStates][cat] += NumberOfEventsSR[fs][cat];
  for(int fs=0; fs<nFinalStates+1; fs++)
    for(int cat=0; cat<nCategories; cat++)
      NumberOfEventsSR[fs][nCategories] += NumberOfEventsSR[fs][cat];

  //---------- Print a table
  int w1 = 12;
  int w2 = 5;
  cout<<fixWidth("",w1,1);
  for(int fs=0; fs<nFinalStates; fs++) cout<<"  "<<fixWidth(sFinalState[fs2[fs]],w2,0);
  cout<<endl;
  for(int cat=0; cat<nCategories+1; cat++){
    cout<<fixWidth(sCategory[cat],w1,1);
    for(int fs=0; fs<nFinalStates; fs++) cout<<"  "<<fixWidth(string(Form("%i",NumberOfEventsSR[fs2[fs]][cat])),w2,0); 
    cout<<endl;
  }
  cout<<endl;

}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////////////////////// fit signal yields vs. mH ///////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void DrawYieldFits(string outputDirectory, TGraphErrors* g[nProcesses][nFinalStates][nCategories+1], TF1* f[nProcesses][nFinalStates][nCategories+1], double lumi){

  setTDRStyle();

  string outputPath = string(Form("%s/CanvasesSignalFits/",outputDirectory.c_str()));
  gSystem->Exec(("mkdir -p "+outputPath).c_str());
  TCanvas* cYield[nProcesses][nCategories];
  TLegend* lgd[nProcesses][nCategories];
  bool first;

  for(int pr=0; pr<nProcesses; pr++){
    if(!useProcessInMeas[pr]) continue;
    if(!isSignal[pr]) continue;
    
    for(int cat=0; cat<nCategories; cat++){
      
      string canvasName = string(Form("cFits_%s_%s",sProcess[pr].c_str(),sCategory[cat].c_str()));
      cYield[pr][cat] = new TCanvas(canvasName.c_str(),canvasName.c_str(),500,500);
      lgd[pr][cat] = new TLegend(0.2,0.73,0.4,0.88);
      lgd[pr][cat]->SetFillStyle(0);
      //lgd[pr][cat]->SetBorderSize(0);
      first = true;
      
      for(int fs=0; fs<nFinalStates; fs++){
	if(MERGE2E2MU && fs==fs2mu2e) continue;

	g[pr][fs][cat]->GetXaxis()->SetTitle("generated m_{H}");
	g[pr][fs][cat]->GetYaxis()->SetTitle(Form("expected yield in %.1f fb^{-1}",lumi));
	g[pr][fs][cat]->GetXaxis()->SetTitleOffset(1.2);
	g[pr][fs][cat]->GetYaxis()->SetTitleOffset(1.6);
	g[pr][fs][cat]->GetXaxis()->SetLabelFont(42);
	g[pr][fs][cat]->GetYaxis()->SetLabelFont(42);
	g[pr][fs][cat]->GetXaxis()->SetLabelSize(0.04);
	g[pr][fs][cat]->GetYaxis()->SetLabelSize(0.04);
	g[pr][fs][cat]->GetXaxis()->SetTitleFont(42);
	g[pr][fs][cat]->GetYaxis()->SetTitleFont(42);
	g[pr][fs][cat]->GetXaxis()->SetTitleSize(0.04);
	g[pr][fs][cat]->GetYaxis()->SetTitleSize(0.04);
	g[pr][fs][cat]->SetMarkerStyle(fsMarkerStyle[fs]);
	g[pr][fs][cat]->SetMarkerColor(catColor[cat]);
	g[pr][fs][cat]->SetMinimum(0);
	g[pr][fs][cat]->GetXaxis()->SetLimits(fMHPoint[1]-5.,fMHPoint[nMHPoints-1]+5.);

	f[pr][fs][cat]->SetLineColor(catColor[cat]);
	f[pr][fs][cat]->SetLineWidth(1);

	lgd[pr][cat]->AddEntry(g[pr][fs][cat],sFinalState[fs].c_str(),"p");
      }

      for(int fs=nFinalStates-1; fs>=0; fs--){
	if(MERGE2E2MU && fs==fs2mu2e) continue;
	g[pr][fs][cat]->Draw(first?"AP":"P");
	f[pr][fs][cat]->Draw("SAME");
	first = false;
      }

      //lgd[pr][cat]->AddEntry(f[pr][0][cat],sCategory[cat].c_str(),"l");
      lgd[pr][cat]->Draw();

      TPaveText* pav = new TPaveText(0.13,0.93,0.8,1.,"brNDC");
      pav->SetFillStyle(0);
      pav->SetBorderSize(0);
      pav->SetTextAlign(11);
      pav->SetTextSize(0.04);
      pav->AddText((sProcess[pr]+", "+sCategory[cat]).c_str());
      pav->Draw();

      SaveCanvas(outputPath,cYield[pr][cat]);

    }
  }

}

void fitSignalYields(string outputDirectory, double lumi, double sqrts, double m4lMin, double m4lMax)
{

  //---------- Retrieve yields from the ROOT file

  TFile* fInYields = TFile::Open(Form("yields_%iTeV_m4l%.1f-%.1f_%.3ffb-1%s.root",(int)sqrts,m4lMin,m4lMax,lumi,(MERGE2E2MU?"_m":"")));
  TH1F* hTemp;
  Float_t yield[nProcesses][nMHPoints][nFinalStates][nCategories+1];
  Float_t yieldStatError[nProcesses][nMHPoints][nFinalStates][nCategories+1];
  for(int pr=0; pr<nProcesses; pr++){
    if(!useProcessInMeas[pr]) continue;
    if(!isSignal[pr]) continue;
    for(int mp=0; mp<nMHPoints; mp++){
      if(hasMHPoint[pr][mp] && isSignal[pr]){
	for(int fs=0; fs<nFinalStates; fs++){
	  if(MERGE2E2MU && fs==fs2mu2e) continue;
	  for(int cat=0; cat<nCategories+1; cat++){
	    hTemp = (TH1F*)fInYields->Get(Form("hYield_%s%s_%s_%s_%s_%s",sProcess[pr].c_str(),sMHPoint[mp].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str(),sResonantStatus[nRS].c_str(),sSystShift[nominal].c_str()));
	    yield[pr][mp][fs][cat] = hTemp->GetBinContent(1);
	    yieldStatError[pr][mp][fs][cat] = hTemp->GetBinError(1);//0.;//
	  }
	}
      }
    }
  }

  //---------- Fit yield vs. mH

  TFile* fOutYieldFunctions = new TFile(Form("yieldFunctions_%iTeV_m4l%.1f-%.1f_%.3ffb-1%s.root",(int)sqrts,m4lMin,m4lMax,lumi,(MERGE2E2MU?"_m":"")),"recreate");
  TGraphErrors* gYield[nProcesses][nFinalStates][nCategories+1];
  TF1* fYield[nProcesses][nFinalStates][nCategories+1];

  for(int pr=0; pr<nProcesses; pr++){
    if(!useProcessInMeas[pr]) continue;
    if(!isSignal[pr]) continue;

    for(int fs=0; fs<nFinalStates; fs++){
      if(MERGE2E2MU && fs==fs2mu2e) continue;

      for(int cat=0; cat<nCategories+1; cat++){

	gYield[pr][fs][cat] = new TGraphErrors(nMHPointsProcess[pr]);

	bool largeErrors = false;
	int iPoint = 0;
	for(int mp=0; mp<nMHPoints; mp++){
	  if(hasMHPoint[pr][mp]){
	    gYield[pr][fs][cat]->SetPoint(iPoint,fMHPoint[mp],yield[pr][mp][fs][cat]);
	    gYield[pr][fs][cat]->SetPointError(iPoint,0.,yieldStatError[pr][mp][fs][cat]);
	    iPoint++;
	    if(yieldStatError[pr][mp][fs][cat]>0.06*yield[pr][mp][fs][cat]) largeErrors = true;
	  }
	}
	int order = largeErrors ? 1 : orderOfPolynomial[pr] ; // if error bars are too large, take a linear function

	string fName = (string)Form("f_%s_%s_%s",sProcess[pr].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str());
	//fYield[pr][fs][cat] = new TF1(fName.c_str(),Form("pol%i",order),fMHPoint[1],fMHPoint[nMHPoints-1]); 
	fYield[pr][fs][cat] = new TF1(fName.c_str(),Form("pol%i",order),120,130); 
	
	cout<<endl<<sProcess[pr]<<", "<<sFinalState[fs]<<", "<<sCategory[cat]<<endl;
	gYield[pr][fs][cat]->Fit(fYield[pr][fs][cat],"N S");

	if(fYield[pr][fs][cat]->Eval(120.)<0. || fYield[pr][fs][cat]->Eval(130.)<fYield[pr][fs][cat]->Eval(120.)){ // if the fitted function decreases or goes negative, take a constant function
	  delete fYield[pr][fs][cat];
	  fYield[pr][fs][cat] = new TF1(fName.c_str(),Form("pol%i",0),120,130);
	  gYield[pr][fs][cat]->Fit(fYield[pr][fs][cat],"N S");
	}
	  
	//cout<<fYield[pr][fs][cat]->GetExpFormula("P")<<endl;
	fYield[pr][fs][cat]->Write();
	
      }
    }
  }

  DrawYieldFits(outputDirectory,gYield,fYield,lumi);

  fOutYieldFunctions->Close();
  delete fOutYieldFunctions; 

}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////// prepare fragments ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void generateFragments(string outputDirectory, double lumi, double sqrts, double m4lMin, double m4lMax, string mHoption)
{

  bool doSystStudy = (mHoption=="param"||mHoption=="125");
  for(int pr=0; pr<nProcesses; pr++) if(isSignal[pr] && useProcessInMeas[pr] && !hasMHPoint[pr][indexOf125]){ doSystStudy = false; break; }

  //---------- Retrieve yields and functions from the ROOT file

  TFile* fInYields = TFile::Open(Form("yields_%iTeV_m4l%.1f-%.1f_%.3ffb-1%s.root",(int)sqrts,m4lMin,m4lMax,lumi,(MERGE2E2MU?"_m":"")));
  TH1F* hTemp;
  Float_t yield[nProcesses][nFinalStates+1][nCategories+1];
  Float_t yieldStatError[nProcesses][nFinalStates+1][nCategories+1];
  Float_t yield_4l_forSyst[nProcesses][nCategories+1][nSystShifts];
  for(int pr=0; pr<nProcesses; pr++){
    if(isSignal[pr]&&(mHoption=="param"||mHoption=="125")&&!hasMHPoint[pr][indexOf125]) continue;
    for(int fs=0; fs<nFinalStates+1; fs++){
      if(MERGE2E2MU && fs==fs2mu2e) continue;
      for(int cat=0; cat<nCategories+1; cat++){
	hTemp = (TH1F*)fInYields->Get(Form("hYield_%s%s_%s_%s_%s_%s",sProcess[pr].c_str(),(isSignal[pr]?(mHoption!="param"?mHoption.c_str():"125"):""),sFinalState[fs].c_str(),sCategory[cat].c_str(),sResonantStatus[nRS].c_str(),sSystShift[nominal].c_str()));
	yield[pr][fs][cat] = hTemp->GetBinContent(1);
	yieldStatError[pr][fs][cat] = hTemp->GetBinError(1);
      }
    }
    if(doSystStudy){
      for(int cat=0; cat<nCategories+1; cat++){
	for(int sy=0; sy<nSystShifts; sy++){
	  hTemp = (TH1F*)fInYields->Get(Form("hYield_%s%s_%s_%s_%s_%s",sProcess[pr].c_str(),(isSignal[pr]?"125":""),sFinalState[nFinalStates].c_str(),sCategory[cat].c_str(),sResonantStatus[nRS].c_str(),sSystShift[sy].c_str()));
	  yield_4l_forSyst[pr][cat][sy] = hTemp->GetBinContent(1);
	}
      }
    }
  }

  TF1* fYield[nProcesses][nFinalStates][nCategories+1];
  TFile* fInYieldFunctions;
  if(mHoption=="param"){
    fInYieldFunctions = TFile::Open(Form("yieldFunctions_%iTeV_m4l%.1f-%.1f_%.3ffb-1%s.root",(int)sqrts,m4lMin,m4lMax,lumi,(MERGE2E2MU?"_m":"")));
    for(int pr=0; pr<nProcesses; pr++){
      if(!useProcessInMeas[pr]) continue;
      if(!isSignal[pr]) continue;
      for(int cat=0; cat<nCategories+1; cat++){
	for(int fs=0; fs<nFinalStates; fs++){
	  if(MERGE2E2MU && fs==fs2mu2e) continue;
	  fYield[pr][fs][cat] = (TF1*)fInYieldFunctions->Get(Form("f_%s_%s_%s",sProcess[pr].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str()));
	}	 
	if(isSignal[pr]&&(mHoption=="param"||mHoption=="125")&&!hasMHPoint[pr][indexOf125]){
	  yield[pr][nFinalStates][cat] = 0.;
	  for(int fs=0; fs<nFinalStates; fs++){
	    if(MERGE2E2MU && fs==fs2mu2e) continue;
	    yield[pr][fs][cat] = fYield[pr][fs][cat]->Eval(125.);
	    yield[pr][nFinalStates][cat] += yield[pr][fs][cat];
	  }
	}
      }
    }
  }

  //Z+X yield
  Float_t yieldZPlusX[nFinalStates+1][nCategories+1];
  if(DOZPLUSXFROMRUN2COMBINEDSHAPE && MERGE2E2MU){
    getZPlusXYields_Run2CombinedShape_InCateg(yieldZPlusX[fs4mu],yieldZPlusX[fs4e],yieldZPlusX[fs2e2mu],yieldZPlusX[nFinalStates],m4lMin,m4lMax);
  }


  //---------- Prepare the yaml fragments
  float totalYieldSgnl[nFinalStates+1][nCategories+1];
  float totalYieldBkgd[nFinalStates+1][nCategories+1];
  float totalYield[nFinalStates+1][nCategories+1];
  for(int fs=0; fs<nFinalStates+1; fs++)
    for(int cat=0; cat<nCategories+1; cat++){
      totalYieldSgnl[fs][cat] = 0.;
      totalYieldBkgd[fs][cat] = 0.;
      totalYield[fs][cat] = 0.;
    }

  string outputFileName[nFinalStates];
  ofstream outFile[nFinalStates];
  for(int fs=0; fs<nFinalStates; fs++){
    if(MERGE2E2MU && fs==fs2mu2e) continue;
    
    outputFileName[fs] = string(Form("%s/yields_per_tag_category_%iTeV_%s.yaml",outputDirectory.c_str(),(int)sqrts,sFinalState[fs].c_str())); // (not specifying MH and mass window in the name here)
    outFile[fs].open(outputFileName[fs]);
    
    outFile[fs]<<"---"<<endl;
    outFile[fs]<<"# sqrt(s) = "<<sqrts<<" TeV"<<endl;
    outFile[fs]<<"# integrated luminosity = "<<lumi<<" fb-1"<<endl;
    outFile[fs]<<endl;
    outFile[fs]<<"mass_range: '"<<m4lMin<<", "<<m4lMax<<"'"<<endl;
    outFile[fs]<<"kd_range: '0, 1'"<<endl;
    outFile[fs]<<endl;

    // outFile[fs]<<"# Category numbering convention:"<<endl;
    // for(int cat=0; cat<nCategories; cat++){
    //   outFile[fs]<<"# "<<cat<<" "<<sCategory[cat]<<endl;
    // }
    // outFile[fs]<<endl;

    for(int cat=0; cat<nCategories; cat++){
      if(!USEVHMETTAG && cat==6) continue;

      outFile[fs]<<sCategory[cat]<<": "<<endl;

      for(int pr=0; pr<nProcesses; pr++){
	if(!useProcess[pr]) continue;

	if(isSignal[pr]) totalYieldSgnl[fs][cat] += yield[pr][fs][cat];
	else totalYieldBkgd[fs][cat] += yield[pr][fs][cat];
	totalYield[fs][cat] += yield[pr][fs][cat];

	if(!useProcessInMeas[pr]) continue;
	if(isSignal[pr] && useProcessInMeas[pr] && mHoption=="param"){
	  outFile[fs]<<"    "<<sProcess[pr]<<": ";
	  for(int ord=0; ord<=orderOfPolynomial[pr]; ord++){
	    Float_t param = fYield[pr][fs][cat]->GetParameter(ord);
	    outFile[fs]<<"("<<(param!=param?0.:param); // if param is NaN, it means that a polynomial of lower order was used, so we put 0 for this coefficient
	    for(int ord2=0; ord2<=ord-1; ord2++) outFile[fs]<<"*@0";
	    outFile[fs]<<")";
	    if(ord<orderOfPolynomial[pr]) outFile[fs]<<"+";
	  }
	  outFile[fs]<<endl;
	}else{
	  outFile[fs]<<"    "<<sProcess[pr]<<": '"<<yield[pr][fs][cat]<<"'"<<endl;
	}

      }

      if(DOZPLUSXFROMRUN2COMBINEDSHAPE && MERGE2E2MU){
	outFile[fs]<<"    zjets: '"<<yieldZPlusX[fs][cat]<<"'"<<endl;
	totalYieldBkgd[fs][cat] += yieldZPlusX[fs][cat];
	totalYield[fs][cat] += yieldZPlusX[fs][cat];
      }

      outFile[fs]<<endl;
    }

    outFile[fs].close();

  }

  for(int fs=0; fs<nFinalStates; fs++)
    for(int cat=0; cat<nCategories; cat++){
      totalYieldSgnl[fs][nCategories] += totalYieldSgnl[fs][cat];
      totalYieldBkgd[fs][nCategories] += totalYieldBkgd[fs][cat];
      totalYield[fs][nCategories] += totalYield[fs][cat];
    }
  for(int fs=0; fs<nFinalStates; fs++)
    for(int cat=0; cat<nCategories+1; cat++){
      totalYieldSgnl[nFinalStates][cat] += totalYieldSgnl[fs][cat];
      totalYieldBkgd[nFinalStates][cat] += totalYieldBkgd[fs][cat];
      totalYield[nFinalStates][cat] += totalYield[fs][cat];
    }

  if(doSystStudy){
    float jecUpFactor, jecDnFactor, btagUpFactor, btagDnFactor;
    ofstream outFileSyst;
    
    outFileSyst.open((outputDirectory+"/partof_systematics_expt_13TeV.yaml").c_str());

    outFileSyst<<"JES:"<<endl;
    for(int cat=0; cat<nCategories; cat++){
      if(!USEVHMETTAG && cat==6) continue;
      outFileSyst<<"    "<<sCategory[cat]<<": "<<endl;
      outFileSyst<<"        type: lnN"<<endl;
      for(int pr=0; pr<nProcesses; pr++){
	if(!useProcessInMeas[pr]) continue;
	jecUpFactor = fround( yield_4l_forSyst[pr][cat][jecUp] / yield_4l_forSyst[pr][cat][nominal], 4);
	jecDnFactor = fround( yield_4l_forSyst[pr][cat][jecDn] / yield_4l_forSyst[pr][cat][nominal], 4);
	if(jecUpFactor==1 && jecDnFactor==1) continue;
	outFileSyst<<"        "<<sProcess[pr]<<": "<<Form("%.4f",jecUpFactor)<<"/"<<Form("%.4f",jecDnFactor)<<endl;
      }
    }
    outFileSyst<<endl;

    outFileSyst<<"bTagSF:"<<endl;
    for(int cat=0; cat<nCategories; cat++){
      if(!USEVHMETTAG && cat==6) continue;
      outFileSyst<<"    "<<sCategory[cat]<<": "<<endl;
      outFileSyst<<"        type: lnN"<<endl;
      for(int pr=0; pr<nProcesses; pr++){
	if(!useProcessInMeas[pr]) continue;
	btagUpFactor = fround( yield_4l_forSyst[pr][cat][btagUp] / yield_4l_forSyst[pr][cat][nominal], 4);
	btagDnFactor = fround( yield_4l_forSyst[pr][cat][btagDn] / yield_4l_forSyst[pr][cat][nominal], 4);
	if(btagUpFactor==1 && btagDnFactor==1) continue;
	outFileSyst<<"        "<<sProcess[pr]<<": "<<Form("%.4f",btagUpFactor)<<"/"<<Form("%.4f",btagDnFactor)<<endl;
      }
    }
    outFileSyst<<endl;

    outFileSyst.close();
  }


  //* //Control printout
  string sep = "----------------------------------------------------";
  cout<<sep<<endl<<"Expected yields for integrated luminosity "<<lumi<<"/fb, for m4l in ["<<m4lMin<<","<<m4lMax<<"]GeV"<<endl;
  cout<<sep<<endl<<"Process  Category  4e  4mu  2e2mu  4l"<<endl;
  cout<<sep<<endl;
  for(int pr=0; pr<nProcesses; pr++){
    for(int cat=0; cat<nCategories+1; cat++)
      cout<<sProcess[pr]<<"  "<<fixWidth(sCategory[cat],12,1)<<"  "<<yield[pr][1][cat]<<"  "<<yield[pr][0][cat]<<"  "<<yield[pr][2][cat]<<"  "<<yield[pr][4][cat]<<endl;
    cout<<sep<<endl;
  }
  if(DOZPLUSXFROMRUN2COMBINEDSHAPE && MERGE2E2MU)
    for(int cat=0; cat<nCategories+1; cat++)
      cout<<"Z+X"<<"  "<<fixWidth(sCategory[cat],12,1)<<"  "<<yieldZPlusX[1][cat]<<"  "<<yieldZPlusX[0][cat]<<"  "<<yieldZPlusX[2][cat]<<"  "<<yieldZPlusX[4][cat]<<endl;
  cout<<sep<<endl;
  for(int cat=0; cat<nCategories+1; cat++)
    cout<<"Total H125"<<"  "<<fixWidth(sCategory[cat],12,1)<<"  "<<totalYieldSgnl[1][cat]<<"  "<<totalYieldSgnl[0][cat]<<"  "<<totalYieldSgnl[2][cat]<<"  "<<totalYieldSgnl[4][cat]<<endl;
  cout<<sep<<endl;
  for(int cat=0; cat<nCategories+1; cat++)
    cout<<"Total bkgd"<<"  "<<fixWidth(sCategory[cat],12,1)<<"  "<<totalYieldBkgd[1][cat]<<"  "<<totalYieldBkgd[0][cat]<<"  "<<totalYieldBkgd[2][cat]<<"  "<<totalYieldBkgd[4][cat]<<endl;
  cout<<sep<<endl;
  for(int cat=0; cat<nCategories+1; cat++)
    cout<<"Total exp."<<"  "<<fixWidth(sCategory[cat],12,1)<<"  "<<totalYield[1][cat]<<"  "<<totalYield[0][cat]<<"  "<<totalYield[2][cat]<<"  "<<totalYield[4][cat]<<endl;
  cout<<sep<<endl;
  cout<<endl;
  //*/

  if(PRINTLATEXTABLES){

    //* //per-FS table in LaTeX format for the AN/PAS
    int w1 = 23;
    int w2 = 6;
    int ndec = 2; //1;
    cout<<fixWidth("Channel",w1,1);
    for(int fs=0; fs<nFinalStates; fs++) cout<<" & "<<fixWidth(sFinalState[fs2[fs]],w2+2,0)<<"                ";
    cout<<" \\\\"<<endl;
    cout<<"\\hline"<<endl;
    cout<<fixWidth("\\qqZZ",w1,1);
    for(int fs=0; fs<nFinalStates; fs++) cout<<" & $"<<fixWidth(sround(yield[qqZZ][fs2[fs]][nCategories],ndec),w2,0)<<"^{+ y.y}_{- z.z}$";
    cout<<" \\\\"<<endl;
    cout<<fixWidth("\\ggZZ",w1,1);
    for(int fs=0; fs<nFinalStates; fs++) cout<<" & $"<<fixWidth(sround(yield[ggZZ][fs2[fs]][nCategories],ndec),w2,0)<<"^{+ y.y}_{- z.z}$";
    cout<<" \\\\"<<endl;
    cout<<fixWidth("\\cPZ\\ + X",w1,1);
    for(int fs=0; fs<nFinalStates; fs++) cout<<" & $"<<fixWidth(sround(yieldZPlusX[fs2[fs]][nCategories],ndec),w2,0)<<"^{+ y.y}_{- z.z}$";
    cout<<" \\\\"<<endl;
    cout<<"\\hline"<<endl;
    cout<<fixWidth("Sum of backgrounds",w1,1);
    for(int fs=0; fs<nFinalStates; fs++) cout<<" & $"<<fixWidth(sround(totalYieldBkgd[fs2[fs]][nCategories],ndec),w2,0)<<"^{+ y.y}_{- z.z}$";
    cout<<" \\\\"<<endl;
    cout<<"\\hline"<<endl;
    cout<<fixWidth("Signal ($\\mH=125~\\GeV$)",w1,1);
    for(int fs=0; fs<nFinalStates; fs++) cout<<" & $"<<fixWidth(sround(totalYieldSgnl[fs2[fs]][nCategories],ndec),w2,0)<<"^{+ y.y}_{- z.z}$";
    cout<<" \\\\"<<endl;
    cout<<"\\hline"<<endl;
    cout<<fixWidth("Total expected",w1,1);
    for(int fs=0; fs<nFinalStates; fs++) cout<<" & $"<<fixWidth(sround(totalYield[fs2[fs]][nCategories],ndec),w2,0)<<"^{+ y.y}_{- z.z}$";
    cout<<" \\\\"<<endl;
    cout<<"\\hline"<<endl;
    cout<<fixWidth("Observed",w1,1);
    for(int fs=0; fs<nFinalStates; fs++) cout<<" & XXX";
    cout<<" \\\\"<<endl;
    
    cout<<"\\hline \\hline"<<endl;
    cout<<endl;
    //*/

    //* //per-category table in LaTeX format for the AN/PAS
    int w3 = 9;
    int w4 = 5;
    cout<<fixWidth("category",w3,1)<<" & ";
    for(int pr=0; pr<nProcesses; pr++) if(useProcess[pr]) cout<<fixWidth(sProcess[pr],w4+2,0)<<" & ";
    cout<<fixWidth("Z+X",w4+2,0)<<" & ";
    cout<<fixWidth("Signal",w4+2,0)<<" & ";
    cout<<fixWidth("Total",w4+2,0)<<" & Obs \\\\"<<endl;
    cout<<"\\hline"<<endl;
    for(int cat=0; cat<nCategories; cat++){
      cout<<fixWidth(sCatLatex[cat],w3,1)<<" & ";
      for(int pr=0; pr<nProcesses; pr++) if(useProcess[pr]) cout<<"$"<<fixWidth(sround(yield[pr][4][cat],2),w4,0)<<"$"<<" & ";
      cout<<"$"<<fixWidth(sround(yieldZPlusX[4][cat],2),w4,0)<<"$"<<" & ";
      cout<<"$"<<fixWidth(sround(totalYieldSgnl[4][cat],2),w4,0)<<"$"<<" & ";
      cout<<"$"<<fixWidth(sround(totalYield[4][cat],2),w4,0)<<"$"<<" & XXX \\\\"<<endl;
    }
    cout<<"\\hline"<<endl;
    cout<<fixWidth("Total",w3,1)<<" & ";
    for(int pr=0; pr<nProcesses; pr++) if(useProcess[pr]) cout<<"$"<<fixWidth(sround(yield[pr][4][nCategories],2),w4,0)<<"$"<<" & ";
    cout<<"$"<<fixWidth(sround(yieldZPlusX[4][nCategories],2),w4,0)<<"$"<<" & ";
    cout<<"$"<<fixWidth(sround(totalYieldSgnl[4][nCategories],2),w4,0)<<"$"<<" & ";
    cout<<"$"<<fixWidth(sround(totalYield[4][nCategories],2),w4,0)<<"$"<<" & XXX \\\\"<<endl;
    cout<<"\\hline \\hline"<<endl;
    cout<<endl;
    //*/
  }

}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// main function ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void prepareYields(bool recomputeYields = true) {


  // --------------- definitions ---------------

  // Specify input/output location

  string inputPathSignal = "";
  string inputFileData = "";

  string inputPathqqZZ = inputPathSignal;
  string inputPathggZZ = inputPathSignal;
  string inputFileFakeRates = "../../data/FakeRates/FakeRate_SS_Moriond368.root";

  string outputPath = "YieldFiles";
  string outputPathPlots = "$pl";//"YieldFiles";

  // Define c.o.m. energy (13TeV only for the moment)
  float sqrts = 13.;

  // Define the luminosity
  float lumi = 35.86706;

  // Choose an m4l window
  //* // for datacards
  float m4l_min = 105.;
  float m4l_max = 140.;
  //*/
  /* // for table in AN/PAS
  float m4l_min = 70.;
  float m4l_max = 3000.;
  //*/
  /* // for table in AN/PAS
  float m4l_min = 118.;
  float m4l_max = 130.;
  //*/
  /* // for sync
  float m4l_min = 110.;
  float m4l_max = 150.;
  //*/


  // --------------- processing ---------------

  outputPath = string(Form("%s_m4l_%.0f_%.0f_%.2ffbinv",outputPath.c_str(),m4l_min,m4l_max,lumi));
  gSystem->Exec(("mkdir -p "+outputPath).c_str());
  gSystem->Exec(("mkdir -p "+outputPathPlots).c_str());

  // Compute the yields for all available processes and mH values, and store them in a ROOT file 
  // (to be done only once, it can take a few minutes)
  if(recomputeYields)
    computeYields(inputPathSignal, inputPathqqZZ, inputPathggZZ, lumi, sqrts, m4l_min, m4l_max);

  // Prepare Z+X yields from SS CR, but in the full range, to be then normalized thanks to the shape
  if(BUILDZPLUSXYIELDFROMSSCR)
    computeYieldsZPlusXSS(inputPathData, inputFileFakeRates, outputPath, lumi);

  // /!\ Count observed events
  if(UNBLINDSR) dataEventCounts(inputPathData, m4l_min, m4l_max);

  // Parameterize signal yields as a function of mH
  if(!JUST125) 
    fitSignalYields(outputPathPlots, lumi, sqrts, m4l_min, m4l_max);

  // Prepare the yaml files
  if(JUST125)
    generateFragments(outputPath, lumi, sqrts, m4l_min, m4l_max, "125"); // This takes mH=125 for the signal yield.
  else
    generateFragments(outputPath, lumi, sqrts, m4l_min, m4l_max, "param"); // This puts the yield(mH) expression when possible.

  //generateFragments(outputPath, lumi, sqrts, m4l_min, m4l_max, "120"); // just for debugging

}

