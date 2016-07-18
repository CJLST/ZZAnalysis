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

#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "plotUtils.C"
#include "ZZAnalysis/AnalysisStep/src/kFactors.C"

#include <ZZAnalysis/AnalysisStep/src/Category.cc>
#include <ZZAnalysis/AnalysisStep/src/bitops.cc>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>

using namespace std;

#define DEBUG 0

#define DO1DPLOTS 1
#define DO2DPLOTS 1

#define FINALSTATE 4 // 0:4mu, 1:4e, 2:2e2mu, 3:2mu2e, 4:inclusive
#define MERGE2E2MU 1 // if activated, 2e2mu means "2e2mu and 2mu2e"

#define APPLYKFACTORS 1
#define RESCALETOSMPSIGNALSTRENGTH 0
#define SMPSIGNALSTRENGTH 0.99

#define USEDYANDTTBAR 0
#define REBINDYANDTTBAR 0

#define MASKH125FORLOWMASS 1
#define MASKH125FORHIGHMASS 1
#define MASKDATAFORHIGHMASS 0

#define USEZPLUSXRUN2COMBINEDSHAPE 0 // for variables where it is available (up to now, just m4l, but the current shape is from the 2015 dataset)
#define USEZPLUSXFULLRUN2SS 1
#define SMOOTHZPLUSXFULLRUN2SS 1
#define RENORMALIZEZPLUSXFULLRUN2SS 0
//Z+X normalization from combination of SS and OS method
Float_t normZPlusXFullSR4e    = 0.; // to be defined for the 2016 dataset
Float_t normZPlusXFullSR4mu   = 0.; // to be defined for the 2016 dataset
Float_t normZPlusXFullSR2e2mu = 0.; // to be defined for the 2016 dataset

#define STYLE1DPLOT 2 // 0:Legacy-like 1:Jamboree2015 2:Moriond2016
#define DRAWLINES (STYLE1DPLOT!=1)
#define LINEWIDTH (STYLE1DPLOT==0?2:STYLE1DPLOT==2?1:1)
#define DRAWLABELBYHAND 1
#define DRAWDATAMCRATIO 0

#define NSTYLES2DPLOT 8
#define STYLE2DPLOT 1 // 0:rainbow 1:gray 2:pink 3:orange 4:yellow 5:blue 6:teal 7:col+box 8:blue-yellow
#define MARKERSTYLE 0 // 0:full(legacy) 1:open




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


enum Blindings {fullyblind=0, blindabove110=1, blindbelow150=2, blind110150=3, blind110150andabove500=4, unblinded=5};
const int nBlindings = 6;
string sBlinding[nBlindings] = {"fullyblind", "M4l70To110", "M4l150ToInf", "blind110150", "blind110150andabove500", "unblinded"};
Bool_t plotThisBlinding[7][nBlindings] = {
  {1,1,1,1,1,0},
  {1,0,0,0,0,0},
  {0,1,0,0,0,0},
  {0,0,1,0,0,0},
  {0,0,0,1,0,0},
  {0,0,0,0,1,0},
  {0,1,1,1,0,0},
};
string blindingLabel[nBlindings] = {"", "70 < m_{4#font[12]{l}} < 110 GeV", "m_{4#font[12]{l}} > 150 GeV", "#splitline{m_{4#font[12]{l}} > 70 GeV}{m_{4#font[12]{l}} #notin [110, 150] GeV}", "#splitline{m_{4#font[12]{l}} #in [70, 500] GeV}{m_{4#font[12]{l}} #notin [110, 150] GeV}", ""};
Float_t xHistoBlindLow[nBlindings] = {  0.,  110.,   0., 110., 110.,  0. };
Float_t xHistoBlindUp [nBlindings] = { -1., 3000., 150., 150., 150., -1. };
Float_t xHistoBlind2Low[nBlindings] = {  0.,  0.,  0.,  0.,  500.,  0. };
Float_t xHistoBlind2Up [nBlindings] = { -1., -1., -1., -1., 3000., -1. };

const int nVariables = 35;
string varName[nVariables] = {
  "M4lV1",
  "M4lV1b",
  "M4lV2",
  "M4lV2Refit",
  "M4lV2b",
  "M4lV3",
  "M4lV4",
  "M4lV5",
  "M4l_70110B4",
  "M4l_70110B5",
  "M4l_110150",
  "M4l_105140",
  "M4l_70182",
  "M4l_above150",
  "MZ1V1",
  "MZ1V1Log",
  "MZ1V1_M4L118130",
  "MZ1V2",
  "MZ2V1",
  "MZ2V1Log",
  "MZ2V1_M4L118130",
  "MZ2V2",
  "KD",
  "KD_M4L118130",
  "DjetFisher",
  "DCombVbf2j",//"DjetVbfMela",
  "DCombVbf2j_M4L118130",//"DjetVbfMela_M4L118130",
  "DCombVbf2jLog_M4L118130",//"DjetVbfMelaLog_M4L118130",
  "Pt4l",
  "Eta4l",
  "NExtraLep",
  "NJets",
  "NJetsBTagged",
  "M4l_100180_HighKD",
  "M4l_110150_HighKD",
};
string varXLabel[nVariables] = {
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}}^{refit} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{Z1} (GeV)",
  "m_{Z1} (GeV)",
  "m_{Z1} (GeV)",
  "m_{Z1} (GeV)",
  "m_{Z2} (GeV)",
  "m_{Z2} (GeV)",
  "m_{Z2} (GeV)",
  "m_{Z2} (GeV)",
  "D_{bkg}^{kin}",
  "D_{bkg}^{kin}",
  "D_{jet}",
  "D_{2jet}^{comb.}",//"D_{jet}",
  "D_{2jet}^{comb.}",//"D_{jet}",
  "D_{2jet}^{comb.}",//"D_{jet}",
  "p_{T}^{4#font[12]{l}} (GeV)",
  "#eta^{4#font[12]{l}}",
  "number of additional leptons",
  "number of jets",
  "number of b-tagged jets",
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
};
string varYLabel[nVariables] = {
  "Events / 3 GeV",
  "Events / 3 GeV",
  "Events / 4 GeV",
  "Events / 4 GeV",
  "Events / 4 GeV",
  "Events / 5 GeV",
  "Events / 10 GeV",
  "Events / 10 GeV",
  "Events / 4 GeV",
  "Events / 5 GeV",
  "Events / 2 GeV",
  "Events / 1 GeV",
  "Events / 4 GeV",
  "Events / 20 GeV",
  "Events / 2 GeV",
  "Events / 2 GeV",
  "Events / 4 GeV",
  "Events / 5 GeV",
  "Events / 2 GeV",
  "Events / 2 GeV",
  "Events / 4 GeV",
  "Events / 5 GeV",
  "Events / 0.05",
  "Events / 0.1",
  "Events / 0.1",
  "Events / 0.05",
  "Events / 0.05",
  "Events / 0.05",
  "Events / 10 GeV",
  "Events / 0.5",
  "Events",
  "Events",
  "Events",
  "Events / 3 GeV",
  "Events / 2 GeV",
};
string varCutLabel[nVariables] = {
  "","","","","","","","","","","","","","","","","118 < m_{4#font[12]{l}} < 130 GeV","","","","118 < m_{4#font[12]{l}} < 130 GeV","","","118 < m_{4#font[12]{l}} < 130 GeV","","","118 < m_{4#font[12]{l}} < 130 GeV","118 < m_{4#font[12]{l}} < 130 GeV","","","","","","D_{bkg}^{kin} > 0.5","D_{bkg}^{kin} > 0.5",
};
Bool_t plotThisVar[7][nVariables] = {
  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,},
  {0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,},
  {0,0,1,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,},
  {0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,1,0,0,0,1,0,0,0,0,0,0,0,},//for AN
  {0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,},//for PAS
  {0,0,1,0,1,0,0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,1,0,0,1,0,0,0,0,0,0,0,1,},//unblinding
  {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,},
};
Int_t  varNbin[nVariables] = { 272,  333  , 204, 204, 204, 163,  82,  82,  10,   8,  20,  35,  28,  30,  40,  40,  20,/*  75,*/  30,  54,  54,  27,/*  75,*/  30, 20, 10, 20, 20, 20, 20,  40,  20, 6, 17, 8,  27,  20, };
Float_t varMin[nVariables] = {  70,    1.5,  70,  70,  71,  70,  70,  70,  70,  70, 110, 105,  70, 150,  40,  40,  40,/*   0,*/   0,  12,  12,  12,/*   0,*/   0,  0,  0,  0,  0,  0,  0,   0, -10, 0,  0, 0, 100, 110, };
Float_t varMax[nVariables] = { 886, 1000.5, 886, 886, 887, 885, 890, 890, 110, 110, 150, 140, 182, 750, 120, 120, 120,/* 150,*/ 150, 120, 120, 120,/* 150,*/ 150,  1,  1,  2,  1,  1,  1, 400,  10, 6, 17, 8, 181, 150, };
Bool_t varLogx[nVariables] = {1,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,};
Bool_t varLogy[nVariables] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,1,1,0,0,};
Int_t restrictCountVar[nVariables] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,4,2,0,0,};
Bool_t separateVbf[nVariables] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,};
Float_t varMinFactor[nVariables] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,2000.,0.,0.,0.,2000.,0.,0.,0.,0.,0.,0.,0.,70.,0.,0.,0.,20000.,2000.,100000.,0.,};
Int_t varCMSPos[nVariables] = {33,33, 0,33,33,33,33,33,11,11,11,11, 0,33, 0,33, 0,33, 0,33, 0,33, 0, 0,33,33,33, 0,33,11,11,11,11,11,11,};
Int_t varLegPos[nVariables] = {33,33,33,33,33,33,33,33,33,33,33,33,33,33,11,11,11,11,11,11,33,11,33,33,33,33,33,33,33,33,33,33,33,33,33,};
//Int_t rebinning[nVariables] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,2,1,1,2,1,1,1,1,1,1,1,1,1,1,1,};
Int_t rebinningDYTTbar[nVariables] = {16,1,1,1,1,1,8,1,2,2,1,1,1,3,5,5,5,2,5,5,5,2,4,4,4,4,4,2,2,2,1,1,1,1,1,};

Float_t varMaxCorrector[nBlindings][nVariables] = {
  { 1., 1. , 1.1, 1.1, 1., 1., 1., 1., 1.3, 1., 1., 1., 1., 1., 1., 1., 1.1, 1., 1., 1. , 1.9, 1. , 1.3, 3. , 1., 1., 1., 1., 1., 1.4, 1.,  1., 1., 1. , 1. , },
  { 1., 1. , 1. , 1. , 1., 1., 1., 1., 1. , 1., 1., 1., 1., 1., 1., 1., 1. , 1., 1., 1. , 1. , 1. , 1.1, 1. , 1., 1., 1., 1., 1., 1. , 1.,  1., 1., 1. , 1. , },
  { 1., 1.2, 1. , 1. , 1., 1., 1., 1., 1. , 1., 1., 1., 1., 1., 1., 1., 1. , 1., 1., 1. , 1. , 1. , 1.1, 1. , 1., 1., 1., 1., 1., 1. , 1., 10., 1., 1. , 1. , },
  { 1., 1. , 1. , 1. , 1., 1., 1., 1., 1. , 1., 1., 1., 1., 1., 1., 1., 1. , 1., 1., 1. , 1. , 1. , 1.1, 1. , 1., 1., 1., 1., 1., 1.5, 1., 10., 1., 1.6, 1. , },
  { 1., 1. , 1. , 1. , 1., 1., 1., 1., 1. , 1., 1., 1., 1., 1., 1., 1., 1. , 1., 1., 1. , 1. , 1. , 1.1, 1. , 1., 1., 1., 1., 1., 1.5, 1., 10., 1., 1.6, 1. , },
  { 1., 1. , 1. , 1. , 1., 1., 1., 1., 1. , 1., 1., 1., 1., 1., 1., 1., 1. , 1., 1., 1. , 1. , 1. , 1.1, 1. , 1., 1., 1., 1., 1., 1. , 1.,  5., 1., 1. , 1.2, },
};

const int nVarPairs = 12;
string varPairName[nVarPairs] = {
  "M4lVsKD",
  "M4lVsKD_M4L70110",
  "M4lVsKD_M4L100170",
  "M4lVsKD_M4L170780",
  "M4lVsKD_M4L150700",
  "M4lVsDjet_M4L100170",
  //  "M4lVsDjet_M4L114180",
  "MZ1VsMZ2V1",
  "MZ1VsMZ2V2",
  "MZ1VsMZ2V3",
  "MZ1VsMZ2V2_M4L118130",
  "M4lVsM4lRefit",
  "M4lVsM4lRefit_100180",
};
string varPairXLabel[nVarPairs] = {
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{Z1} (GeV)",
  "m_{Z1} (GeV)",
  "m_{Z1} (GeV)",
  "m_{Z1} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
};
string varPairYLabel[nVarPairs] = {
  "D_{bkg}^{kin}",
  "D_{bkg}^{kin}",
  "D_{bkg}^{kin}",
  "D_{bkg}^{kin}",
  "D_{bkg}^{kin}",
  "D_{jet}",
  "m_{Z2} (GeV)",
  "m_{Z2} (GeV)",
  "m_{Z2} (GeV)",
  "m_{Z2} (GeV)",
  "m_{4#font[12]{l}}^{refit} (GeV)",
  "m_{4#font[12]{l}}^{refit} (GeV)",
};
string varPairCutLabel[nVariables] = {
  "","","","","","","","","","118 < m_{4#font[12]{l}} < 130 GeV","","",
};
Int_t  varPairXNbin[nVarPairs] = { 262,  40,  35, 122, 110,  35,  40,  80,  60,  80, 302,  80, };
Float_t varPairXMin[nVarPairs] = { 100,  70, 100, 170, 150, 100,  40,  40,  75,  40,  70, 100, };
Float_t varPairXMax[nVarPairs] = { 886, 110, 170, 780, 700, 170, 120, 120, 105, 120, 886, 180, };
Int_t  varPairYNbin[nVarPairs] = { 30, 30, 30, 30, 30, 30,  54, 108,  60, 108, 302,  80, };
Float_t varPairYMin[nVarPairs] = {  0,  0,  0,  0,  0,  0,  12,  12,  75,  12,  70, 100, };
Float_t varPairYMax[nVarPairs] = {  1,  1,  1,  1,  1,  1, 120, 120, 105, 120, 886, 180, };
Bool_t varPairLogx[nVarPairs] = {1,0,0,0,0,0,0,0,0,0,0,0,};
Bool_t varPairLogy[nVarPairs] = {0,0,0,0,0,0,0,0,0,0,0,0,};
Int_t varPairLegPos[nVarPairs] = {33,33,33,33,33,33,11,11,11,33,11,};
Bool_t varPairLegIsWhite[nVarPairs] = {1,1,1,1,1,1,0,0,0,0,0,0,};
Bool_t varPairUseGrayStyle[nVarPairs] = {0,0,0,0,0,0,1,1,1,1,0,0,};
Bool_t plotThisVarPair[7][nVarPairs] = {
  {0,0,0,0,0,0,0,0,0,0,0,0,}, 
  {0,1,0,0,1,0,0,1,0,0,0,0,}, 
  {0,0,1,0,1,0,0,0,0,0,0,0,}, //for style tests
  {0,0,1,1,0,1,0,1,0,1,0,0,}, //for AN
  {0,0,1,0,0,1,0,0,0,1,0,0,}, //for PAS 
  {0,0,1,1,0,0,0,1,0,1,0,0,}, //for unblinding
  {1,1,1,1,1,1,1,1,1,1,1,1,}, 
};



enum Process {Data=0, H125=1, H125VBF=2, H125NONVBF=3, qqZZ=4, ggZZ=5, DY=6, ttbar=7};
const int nProcesses = 8;
string sProcess[nProcesses] = {"Data", "H125", "H125VBF", "H125NONVBF", "qqZZ", "ggZZ", "DY", "ttbar"};
string processLabel[nProcesses] = {" Data", " H(125)", " H(125), VBF", " H(125), other", " q#bar{q}#rightarrowZZ, Z#gamma*", " gg#rightarrowZZ, Z#gamma*", " Z + jets", " t#bar{t}"};
Int_t processFillColor[nProcesses] = {
  TColor::GetColor("#000000"), 
  TColor::GetColor((STYLE1DPLOT==0)?"#ffffff":(STYLE1DPLOT==1)?"#ff9090":"#ffb2b2"), //"#ffc18b"),
  TColor::GetColor((STYLE1DPLOT==0)?"#ffffff":(STYLE1DPLOT==1)?"#ff6868":"#ff9b9b"),
  TColor::GetColor((STYLE1DPLOT==0)?"#ff0000":(STYLE1DPLOT==1)?"#ff9090":"#ffdcdc"),
  TColor::GetColor((STYLE1DPLOT==0)?"#99ccff":(STYLE1DPLOT==1)?"#8bc5ff":"#99ccff"), 
  TColor::GetColor((STYLE1DPLOT==0)?"#3366ff":(STYLE1DPLOT==1)?"#2c5Df0":"#4b78ff"), 
  TColor::GetColor((STYLE1DPLOT==0)?"#669966":(STYLE1DPLOT==1)?"#6dae6d":"#669966"), 
  TColor::GetColor("#9a6666")
};
Int_t processLineColor[nProcesses] = {
  TColor::GetColor("#000000"),
  TColor::GetColor((STYLE1DPLOT==0)?"#ff0000":"#cc0000"),//"#770000"), 
  TColor::GetColor((STYLE1DPLOT==0)?"#ff0000":"#cc0000"),//"#770000"), 
  TColor::GetColor((STYLE1DPLOT==0)?"#ff0000":"#cc0000"),//"#770000"), 
  TColor::GetColor("#000099"),
  TColor::GetColor("#000099"),
  TColor::GetColor("#003300"),
  TColor::GetColor("#5f3f3f")
};
Bool_t useProcess[nProcesses] = {1,1,0,0,1,1,USEDYANDTTBAR,USEDYANDTTBAR,};

enum FinalState {fs4mu=0, fs4e=1, fs2e2mu=2, fs2mu2e=3};
const int nFinalStates = 4;
string sFinalState[nFinalStates+1] = {"4mu", "4e", "2e2mu", "2mu2e", "4l"};
string fsLabel[nFinalStates+1] = {" 4#mu", " 4e", " 2e2#mu", " 2#mu2e", " 4#font[12]{l}"};
Int_t fsMarkerStyleFull[nFinalStates+1] = {20,22,21,33,29};
Int_t fsMarkerStyleOpen[nFinalStates+1] = {5,2,4,25,29};
Int_t fsMarkerColor[nFinalStates+1] = {kRed+1,kGreen+2,kAzure,kViolet,kBlack};
//Int_t fsMarkerColor[nFinalStates+1] = {kAzure,kRed+1,kViolet+1,kMagenta+1,kBlack};
Float_t fsROSSS[nFinalStates] = { 1.22, 0.97, 1.30, 0.98 }; //FIXME: recompute this for Run II
string fsLabelForSS[nFinalStates] = {
  "Z1->mu+mu- + mumu(SS)",
  "Z1->e+e- + ee(SS)",
  "Z1->e+e- + mumu(SS)",
  "Z1->mu+mu- + ee(SS)",
};


/* ---------- early RunII categorization proposal 
const int nCategories = 6;
string sCategory[nCategories+1] = {
  "UnTagged",
  "OneJetTagged",
  "VBFTagged", 
  "VHLeptTagged",
  "VHHadrTagged",
  "ttHTagged",
  "inclusive",
};
//*/
/* ---------- Moriond 2016 categorization 
const int nCategories = 2;
string sCategory[nCategories+1] = {
  "UnTagged",
  "VBFTagged", 
  "inclusive", 
};
//*/
//* ---------- Ichep 2016 categorization 
const int nCategories = 6;
string sCategory[nCategories+1] = {
  "UnTagged",
  "VBF1jTagged",
  "VBF2jTagged", 
  "VHLeptTagged",
  "VHHadrTagged",
  "ttHTagged",
  "inclusive",
};
//*/

enum ResonantStatus {resonant=0, nonresonant=1};
const int nResStatuses = 2;
string sResonantStatus[nResStatuses+1] = {"resonant", "nonresonant", "allres"};





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void doHistograms(string inputFilePath_MC, string inputFilePath_Data, double lumi)
{

  //TFile* ggZZKFactorFile = TFile::Open("../../data/kfactors/Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
  //TSpline3* sp = (TSpline3*)ggZZKFactorFile->Get("sp_kfactor_Nominal");

  const int nDatasets = 16;
  string datasets[nDatasets] = {
    // "DoubleMu2015B",
    // "DoubleEG2015B",
    // "MuonEG2015B",
    // "SingleEle2015B",
    // "DoubleMu2015C",
    // "DoubleEG2015C",
    // "MuonEG2015C",
    // "SingleEle2015C",
    // "DoubleMu2015C_50ns",
    // "DoubleEG2015C_50ns",
    // "MuonEG2015C_50ns",
    // "SingleEle2015C_50ns",
    // "DoubleMu2015D",
    // "DoubleEG2015D",
    // "MuonEG2015D",
    // "SingleEle2015D",
    "AllData",
    "ggH125",
    "VBFH125",
    "WplusH125",
    "WminusH125",
    "ZH125",
    "ttH125",
    "ZZTo4l",//"ZZTo4lamcatnlo",//
    "ggTo4e_Contin_MCFM701",//"ggZZ4e",
    "ggTo4mu_Contin_MCFM701",//"ggZZ4mu",
    "ggTo4tau_Contin_MCFM701",//"ggZZ4tau",
    "ggZZ2e2mu",//"ggTo2e2mu_Contin_MCFM701",
    "ggZZ2e2tau",//"ggTo2e2tau_Contin_MCFM701",
    "ggTo2mu2tau_Contin_MCFM701",//"ggZZ2mu2tau",
    "DYJetsToLL_M50",
    "TTTo2L2Nu",//"TTJets",
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
  Float_t overallEventWeight;
  Float_t KFactor_QCD_ggZZ_Nominal;
  Float_t KFactor_EW_qqZZ;
  Float_t KFactor_QCD_qqZZ_dPhi;
  Float_t KFactor_QCD_qqZZ_M;
  Float_t KFactor_QCD_qqZZ_Pt;
  Float_t xsec;
  Short_t ZZsel;
  Float_t ZZMass;
  Float_t ZZMassErr;
  Float_t ZZMassErrCorr;
  Float_t ZZMassRefit;
  Float_t ZZMassRefitErr;
  Float_t ZZMassUnrefitErr;
  Float_t ZZPt;
  Float_t ZZEta;
  Float_t p0plus_VAJHU;
  Float_t bkg_VAMCFM;
  Float_t pvbf_VAJHU_highestPTJets;
  Float_t phjj_VAJHU_highestPTJets;
  Float_t pAux_vbf_VAJHU;
  Float_t phj_VAJHU;
  Float_t pwh_hadronic_VAJHU;
  Float_t pzh_hadronic_VAJHU;
  Float_t Z1Mass;
  Float_t Z2Mass;
  Short_t Z1Flav;
  Short_t Z2Flav;
  Float_t DiJetFisher;
  vector<Float_t> *LepEta = 0;
  vector<Float_t> *LepPhi = 0;
  Short_t nExtraLep;
  vector<Float_t> *ExtraLepEta = 0;
  vector<Float_t> *ExtraLepPhi = 0;
  Short_t nExtraZ;
  Short_t nJets;
  Short_t nJetsBTagged;
  vector<Float_t> *JetPt = 0;
  vector<Float_t> *JetEta = 0;
  vector<Float_t> *JetPhi = 0;
  vector<Float_t> *JetMass = 0;
  vector<Float_t> *JetQGLikelihood = 0;
  Float_t jetPt[99];
  Float_t jetEta[99];
  Float_t jetPhi[99];
  Float_t jetMass[99];
  Float_t jetQGL[99];
  Float_t jetPgOverPq[99];
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

  TH1F* h1[nVariables][nBlindings][nFinalStates+1][nCategories+1][nResStatuses+1][nProcesses];
  TH2F* h2[nVarPairs ][nBlindings][nFinalStates+1][nCategories+1][nResStatuses+1][nProcesses];
  for(int bl=0; bl<nBlindings; bl++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      for(int cat=0; cat<nCategories+1; cat++){
	for(int rs=0; rs<nResStatuses+1; rs++){
	  for(int pr=0; pr<nProcesses; pr++){
	    for(int v=0; v<nVariables; v++){
	      h1[v][bl][fs][cat][rs][pr] = new TH1F(
                 Form("h1_%s_%s_%s_%s_%s_%s",varName[v].c_str(),sBlinding[bl].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str(),sResonantStatus[rs].c_str(),sProcess[pr].c_str()),
		 Form(";%s;%s",varXLabel[v].c_str(),varYLabel[v].c_str()),
		 varNbin[v],varMin[v],varMax[v]);
	    }
	    for(int v2=0; v2<nVarPairs; v2++){
	      h2[v2][bl][fs][cat][rs][pr] = new TH2F(
                 Form("h2_%s_%s_%s_%s_%s_%s",varPairName[v2].c_str(),sBlinding[bl].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str(),sResonantStatus[rs].c_str(),sProcess[pr].c_str()),
		 Form(";%s;%s",varPairXLabel[v2].c_str(),varPairYLabel[v2].c_str()),
		 varPairXNbin[v2],varPairXMin[v2],varPairXMax[v2],
		 varPairYNbin[v2],varPairYMin[v2],varPairYMax[v2]);
	    }
	  }
	}
      }
    }
  }

  vector<Float_t> g2DataX[nVarPairs][nBlindings][nFinalStates+1][nCategories+1][nResStatuses+1];
  vector<Float_t> g2DataY[nVarPairs][nBlindings][nFinalStates+1][nCategories+1][nResStatuses+1];
  vector<Float_t> g2DataEX[nVarPairs][nBlindings][nFinalStates+1][nCategories+1][nResStatuses+1];
  vector<Float_t> g2DataEY[nVarPairs][nBlindings][nFinalStates+1][nCategories+1][nResStatuses+1];
  TGraphErrors* g2Data[nVarPairs][nBlindings][nFinalStates+1][nCategories+1][nResStatuses+1];
  
  int currentProcess;
  int currentFinalState;
  int currentCategory;
  int currentResStatus;


  //---------- Will loop over all datasets

  for(int d=0; d<nDatasets-2; d++){

    //----- assign dataset to correct process
    currentProcess = -1;
    if(datasets[d].find("2015")!=string::npos||
       datasets[d]=="AllData")
      currentProcess = Data;
    if(datasets[d]=="ggH125") currentProcess = H125NONVBF;
    if(datasets[d]=="VBFH125") currentProcess = H125VBF;
    if(datasets[d]=="WplusH125") currentProcess = H125NONVBF;
    if(datasets[d]=="WminusH125") currentProcess = H125NONVBF;
    if(datasets[d]=="ZH125") currentProcess = H125NONVBF;
    if(datasets[d]=="ttH125") currentProcess = H125NONVBF;
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
    if(datasets[d]=="DYJetsToLL_M50") currentProcess = DY;
    if(datasets[d]=="TTJets"||
       datasets[d]=="TTTo2L2Nu")
      currentProcess = ttbar;

    string inputFileName = string(Form("%s%s/ZZ4lAnalysis.root",(currentProcess==Data?inputFilePath_Data:inputFilePath_MC).c_str(),datasets[d].c_str()));
    inputFile[d] = TFile::Open(inputFileName.c_str());

    hCounters[d] = (TH1F*)inputFile[d]->Get("ZZTree/Counters");
    NGenEvt[d] = (Long64_t)hCounters[d]->GetBinContent(1);
    gen_sumWeights[d] = (Long64_t)hCounters[d]->GetBinContent(40);
    partialSampleWeight[d] = lumi * 1000 / gen_sumWeights[d] ;

    inputTree[d] = (TTree*)inputFile[d]->Get("ZZTree/candTree");
    inputTree[d]->SetBranchAddress("RunNumber", &nRun);
    inputTree[d]->SetBranchAddress("EventNumber", &nEvent);
    inputTree[d]->SetBranchAddress("LumiNumber", &nLumi);
    inputTree[d]->SetBranchAddress("Nvtx", &Nvtx);
    inputTree[d]->SetBranchAddress("NObsInt", &NObsInt);
    inputTree[d]->SetBranchAddress("NTrueInt", &NTrueInt);
    inputTree[d]->SetBranchAddress("overallEventWeight", &overallEventWeight);
    inputTree[d]->SetBranchAddress("KFactor_QCD_ggZZ_Nominal", &KFactor_QCD_ggZZ_Nominal);
    inputTree[d]->SetBranchAddress("KFactor_EW_qqZZ", &KFactor_EW_qqZZ);
    inputTree[d]->SetBranchAddress("KFactor_QCD_qqZZ_dPhi", &KFactor_QCD_qqZZ_dPhi);
    inputTree[d]->SetBranchAddress("KFactor_QCD_qqZZ_M", &KFactor_QCD_qqZZ_M);
    inputTree[d]->SetBranchAddress("KFactor_QCD_qqZZ_Pt", &KFactor_QCD_qqZZ_Pt);
    inputTree[d]->SetBranchAddress("xsec", &xsec);
    inputTree[d]->SetBranchAddress("ZZsel", &ZZsel);
    inputTree[d]->SetBranchAddress("ZZMass", &ZZMass);
    inputTree[d]->SetBranchAddress("ZZMassErr", &ZZMassErr);
    inputTree[d]->SetBranchAddress("ZZMassErrCorr", &ZZMassErrCorr);
    inputTree[d]->SetBranchAddress("ZZMassRefit", &ZZMassRefit);
    inputTree[d]->SetBranchAddress("ZZMassRefitErr", &ZZMassRefitErr);
    inputTree[d]->SetBranchAddress("ZZMassUnrefitErr", &ZZMassUnrefitErr);
    inputTree[d]->SetBranchAddress("ZZPt", &ZZPt);
    inputTree[d]->SetBranchAddress("ZZEta", &ZZEta);
    inputTree[d]->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU);
    inputTree[d]->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM);
    inputTree[d]->SetBranchAddress("pvbf_VAJHU_highestPTJets", &pvbf_VAJHU_highestPTJets);
    inputTree[d]->SetBranchAddress("phjj_VAJHU_highestPTJets", &phjj_VAJHU_highestPTJets);
    inputTree[d]->SetBranchAddress("pAux_vbf_VAJHU", &pAux_vbf_VAJHU);
    inputTree[d]->SetBranchAddress("phj_VAJHU", &phj_VAJHU);
    inputTree[d]->SetBranchAddress("pwh_hadronic_VAJHU", &pwh_hadronic_VAJHU);
    inputTree[d]->SetBranchAddress("pzh_hadronic_VAJHU", &pzh_hadronic_VAJHU);
    inputTree[d]->SetBranchAddress("Z1Mass", &Z1Mass);
    inputTree[d]->SetBranchAddress("Z2Mass", &Z2Mass);
    inputTree[d]->SetBranchAddress("Z1Flav", &Z1Flav);
    inputTree[d]->SetBranchAddress("Z2Flav", &Z2Flav);
    inputTree[d]->SetBranchAddress("LepEta", &LepEta);
    inputTree[d]->SetBranchAddress("LepPhi", &LepPhi);
    inputTree[d]->SetBranchAddress("nExtraLep", &nExtraLep);
    inputTree[d]->SetBranchAddress("ExtraLepEta", &ExtraLepEta);
    inputTree[d]->SetBranchAddress("ExtraLepPhi", &ExtraLepPhi);
    inputTree[d]->SetBranchAddress("nExtraZ", &nExtraZ);
    inputTree[d]->SetBranchAddress("nCleanedJetsPt30", &nJets);
    inputTree[d]->SetBranchAddress("nCleanedJetsPt30BTagged", &nJetsBTagged);
    inputTree[d]->SetBranchAddress("JetPt", &JetPt);
    inputTree[d]->SetBranchAddress("JetEta", &JetEta);
    inputTree[d]->SetBranchAddress("JetPhi", &JetPhi);
    inputTree[d]->SetBranchAddress("JetMass", &JetMass);
    inputTree[d]->SetBranchAddress("JetQGLikelihood", &JetQGLikelihood);
    inputTree[d]->SetBranchAddress("DiJetFisher", &DiJetFisher);
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

    Long64_t entries = inputTree[d]->GetEntries();
    cout<<"Processing dataset "<<datasets[d]<<" ("<<entries<<" entries) ..."<<endl;

    for (Long64_t z=0; z<entries; ++z){

      if(DEBUG && z>1000) break;

      inputTree[d]->GetEntry(z);
      
      if(LepEta->size()!=4){
        cout<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", stored "<<LepEta->size()<<" leptons instead of 4"<<endl;
        continue;
      }

      if( !(ZZsel>=90) ) continue;

      Float_t kfactor = 1.;
      if(APPLYKFACTORS){
	
	if(currentProcess==qqZZ){
	  
	  //kfactor = 1.065;
	  
	  //kfactor = (GenZ1Flav==GenZ2Flav) ? 1.09 : 1.11 ;
	  
	  //kfactor = 1.1; // Jamboree
	  
	  kfactor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M ; // as of Moriond2016 

	}else if(currentProcess==ggZZ){
	  
	  //kfactor = 2.;
	  
	  //kfactor = 1.7; // Jamboree
	  
	  kfactor = KFactor_QCD_ggZZ_Nominal; // as of Moriond2016

	  //kfactor = (float)sp->Eval(GenHMass); // (same as previous)
	  
	}
	
      }

      Double_t eventWeight = partialSampleWeight[d] * xsec * kfactor * overallEventWeight ;


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

      for(int j=0; j<nJets; j++){
	jetPt[j] = JetPt->at(j);
	jetEta[j] = JetEta->at(j);
	jetPhi[j] = JetPhi->at(j);
	jetMass[j] = JetMass->at(j);
	jetQGL[j] = JetQGLikelihood->at(j);
	jetPgOverPq[j] = 1./JetQGLikelihood->at(j) - 1.;
      }
      /* ---------- early RunII categorization proposal 
      currentCategory = category(
	 nExtraLep,
	 ZZPt,
	 ZZMass,
	 nJets, 
	 nJetsBTagged,
	 jetPt,
	 jetEta,
	 jetPhi,
	 jetMass,
	 DiJetFisher
	 );
      //*/
      /* ---------- Moriond 2016 categorization 
      currentCategory = categoryMor16(
	 nJets,
         pvbf_VAJHU_highestPTJets,
         phjj_VAJHU_highestPTJets
         );
      //*/
      //* ---------- Ichep 2016 categorization 
      currentCategory = categoryIchep16(
	 nExtraLep,
	 nExtraZ,
	 nJets,
	 nJetsBTagged,
	 jetQGL,
	 phjj_VAJHU_highestPTJets,
	 phj_VAJHU,
	 pvbf_VAJHU_highestPTJets,
	 pAux_vbf_VAJHU,
	 pwh_hadronic_VAJHU,
	 pzh_hadronic_VAJHU,
	 jetPhi,
	 ZZMass,
	 false
	 );
      //*/


      //----- here, define resonant signal as H->4l where l=e,mu (excluding decays to taus and 'wrong signal' from associated production)

      if(currentProcess==Data){
	currentResStatus = resonant;
      }else{	
	Short_t GenHLepId[4] = {GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id};
	Int_t nGenHLep = 0;
	for(Int_t iGenHLep=0; iGenHLep<4; iGenHLep++){
	  if(abs(GenHLepId[iGenHLep])==11 || abs(GenHLepId[iGenHLep])==13)
	    nGenHLep++;
	}
	currentResStatus = (nGenHLep==4) ? resonant : nonresonant ;
      }

      //----- fill histograms

      Float_t KD = p0plus_VAJHU / ( p0plus_VAJHU + bkg_VAMCFM ) ;
      Float_t vbfMela = pvbf_VAJHU_highestPTJets / ( phjj_VAJHU_highestPTJets + pvbf_VAJHU_highestPTJets );
      Float_t DCombVbf2j = (nJets>=2) ? 1/(1+ phjj_VAJHU_highestPTJets/pvbf_VAJHU_highestPTJets * TMath::Power(jetPgOverPq[0]*jetPgOverPq[1],1/3.) ) : -2 ;
      Float_t varVal[nVariables] = {
	ZZMass,
	ZZMass,
	ZZMass,
	ZZMassRefit,
	ZZMass,
	ZZMass,
	ZZMass,
	ZZMass,
	ZZMass,
	ZZMass,
	ZZMass,
	ZZMass,
	ZZMass,
	ZZMass,
	Z1Mass,
	Z1Mass,
	Z1Mass,
	Z1Mass,
	Z2Mass,
	Z2Mass,
	Z2Mass,
	Z2Mass,
	KD,
	KD,
	DiJetFisher,
	DCombVbf2j,//vbfMela,
	DCombVbf2j,//vbfMela,
	DCombVbf2j,//vbfMela,
	ZZPt,
	ZZEta,
	(Float_t)nExtraLep,
	(Float_t)nJets,
	(Float_t)nJetsBTagged,
	ZZMass,
	ZZMass,
      };
      Bool_t varPassCut[nVariables] = {
	1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,118<=ZZMass&&ZZMass<=130,1,1,1,118<=ZZMass&&ZZMass<=130,1,1,118<=ZZMass&&ZZMass<=130,nJets>=2,nJets>=2,nJets>=2&&118<=ZZMass&&ZZMass<=130,nJets>=2&&118<=ZZMass&&ZZMass<=130,1,1,1,1,1,KD>0.5,KD>0.5,
      };
      Float_t varPairVal[nVarPairs][2] = {
	{ ZZMass, KD },
	{ ZZMass, KD },
	{ ZZMass, KD },
	{ ZZMass, KD },
	{ ZZMass, KD },
	{ ZZMass, vbfMela },
	{ Z1Mass, Z2Mass },
	{ Z1Mass, Z2Mass },
	{ Z1Mass, Z2Mass },
	{ Z1Mass, Z2Mass },
	{ ZZMass, ZZMassRefit },
	{ ZZMass, ZZMassRefit },
      };
      Bool_t varPairPassCut[nVarPairs] = {
	1,1,1,1,1,nJets>=2,1,1,1,118<=ZZMass&&ZZMass<=130,1,1,
      };

      bool fillM4l[nBlindings] = {
	currentProcess!=Data,
	(currentProcess!=Data || ZZMass<110),
	(currentProcess!=Data || ZZMass>150),
	(currentProcess!=Data || ZZMass<110 || ZZMass>150),
	(currentProcess!=Data || ZZMass<110 || (ZZMass>150&&ZZMass<500)),
	true,
      };
      bool fillOtherThanM4l[nBlindings] = {
	currentProcess!=Data,
	ZZMass<110,
	ZZMass>150,
	ZZMass<110 || ZZMass>150,
	ZZMass<110 || (ZZMass>150&&ZZMass<500),
	true,
      };

      //float m4lerr = ZZMassErr;
      //float m4lerr = ZZMassErrCorr;
      float m4lerr = ZZMassUnrefitErr;

      for(int bl=0; bl<nBlindings; bl++){
	for(int v=0; v<nVariables; v++){
	  if( (varName[v].find("M4l")==0 && fillM4l         [bl]) ||
	      (varName[v].find("M4l")!=0 && fillOtherThanM4l[bl])    ){
	    if(varPassCut[v]){
	      h1[v][bl][currentFinalState][currentCategory][currentResStatus][currentProcess]->Fill(varVal[v],(currentProcess==Data)?1.:eventWeight);
	    }
	  }
	}
	for(int v2=0; v2<nVarPairs; v2++){
	  if( (varPairName[v2].find("M4l")==0 && fillM4l         [bl]) ||
	      (varPairName[v2].find("M4l")!=0 && fillOtherThanM4l[bl])    ){
	    if(varPairPassCut[v2]){
	      h2[v2][bl][currentFinalState][currentCategory][currentResStatus][currentProcess]->Fill(varPairVal[v2][0],varPairVal[v2][1],(currentProcess==Data)?1.:eventWeight);
	      if(currentProcess==Data){
		g2DataX[v2][bl][currentFinalState][currentCategory][currentResStatus].push_back(varPairVal[v2][0]);
		g2DataY[v2][bl][currentFinalState][currentCategory][currentResStatus].push_back(varPairVal[v2][1]);
		g2DataEX[v2][bl][currentFinalState][currentCategory][currentResStatus].push_back((varPairName[v2].find("M4l")==0)?m4lerr:0.);
		g2DataEY[v2][bl][currentFinalState][currentCategory][currentResStatus].push_back(0.);
	      }	      
	    }
	  }
	}
      }

    } // end for entries

  } // end for datasets


  //---------- Fill 'inclusive' histograms
  for(int bl=0; bl<nBlindings; bl++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      for(int cat=0; cat<nCategories+1; cat++){
	for(int rs=0; rs<nResStatuses+1; rs++){
	  for(int v=0; v<nVariables; v++){
	    h1[v][bl][fs][cat][rs][H125]->Add(h1[v][bl][fs][cat][rs][H125VBF]);
	    h1[v][bl][fs][cat][rs][H125]->Add(h1[v][bl][fs][cat][rs][H125NONVBF]);
	  }
	  for(int v2=0; v2<nVarPairs; v2++){
	    h2[v2][bl][fs][cat][rs][H125]->Add(h2[v2][bl][fs][cat][rs][H125VBF]);
	    h2[v2][bl][fs][cat][rs][H125]->Add(h2[v2][bl][fs][cat][rs][H125NONVBF]);
	  }
	}
      }
    }
    for(int pr=0; pr<nProcesses; pr++){
      for(int fs=0; fs<nFinalStates; fs++){
	for(int cat=0; cat<nCategories; cat++){
	  for(int rs=0; rs<nResStatuses; rs++){
	    for(int v=0; v<nVariables; v++){
	      h1[v][bl][nFinalStates][cat][rs][pr]->Add(h1[v][bl][fs][cat][rs][pr]);
	    }
	    for(int v2=0; v2<nVarPairs; v2++){
	      h2[v2][bl][nFinalStates][cat][rs][pr]->Add(h2[v2][bl][fs][cat][rs][pr]);
	      if(pr==Data){
		g2DataX[v2][bl][nFinalStates][cat][rs].insert(g2DataX[v2][bl][nFinalStates][cat][rs].end(), g2DataX[v2][bl][fs][cat][rs].begin(), g2DataX[v2][bl][fs][cat][rs].end());
		g2DataY[v2][bl][nFinalStates][cat][rs].insert(g2DataY[v2][bl][nFinalStates][cat][rs].end(), g2DataY[v2][bl][fs][cat][rs].begin(), g2DataY[v2][bl][fs][cat][rs].end());
		g2DataEX[v2][bl][nFinalStates][cat][rs].insert(g2DataEX[v2][bl][nFinalStates][cat][rs].end(), g2DataEX[v2][bl][fs][cat][rs].begin(), g2DataEX[v2][bl][fs][cat][rs].end());
		g2DataEY[v2][bl][nFinalStates][cat][rs].insert(g2DataEY[v2][bl][nFinalStates][cat][rs].end(), g2DataEY[v2][bl][fs][cat][rs].begin(), g2DataEY[v2][bl][fs][cat][rs].end());
	      }
	    }
	  }
	}
      }
      for(int fs=0; fs<nFinalStates+1; fs++){
	for(int cat=0; cat<nCategories; cat++){
	  for(int rs=0; rs<nResStatuses; rs++){
	    for(int v=0; v<nVariables; v++){
	      h1[v][bl][fs][nCategories][rs][pr]->Add(h1[v][bl][fs][cat][rs][pr]);
	    }
	    for(int v2=0; v2<nVarPairs; v2++){
	      h2[v2][bl][fs][nCategories][rs][pr]->Add(h2[v2][bl][fs][cat][rs][pr]);
	      if(pr==Data){
		g2DataX[v2][bl][fs][nCategories][rs].insert(g2DataX[v2][bl][fs][nCategories][rs].end(), g2DataX[v2][bl][fs][cat][rs].begin(), g2DataX[v2][bl][fs][cat][rs].end());
		g2DataY[v2][bl][fs][nCategories][rs].insert(g2DataY[v2][bl][fs][nCategories][rs].end(), g2DataY[v2][bl][fs][cat][rs].begin(), g2DataY[v2][bl][fs][cat][rs].end());
		g2DataEX[v2][bl][fs][nCategories][rs].insert(g2DataEX[v2][bl][fs][nCategories][rs].end(), g2DataEX[v2][bl][fs][cat][rs].begin(), g2DataEX[v2][bl][fs][cat][rs].end());
		g2DataEY[v2][bl][fs][nCategories][rs].insert(g2DataEY[v2][bl][fs][nCategories][rs].end(), g2DataEY[v2][bl][fs][cat][rs].begin(), g2DataEY[v2][bl][fs][cat][rs].end());
	      }
	    }
	  }
	}
	for(int cat=0; cat<nCategories+1; cat++){
	  for(int rs=0; rs<nResStatuses; rs++){
	    for(int v=0; v<nVariables; v++){
	      h1[v][bl][fs][cat][nResStatuses][pr]->Add(h1[v][bl][fs][cat][rs][pr]);
	    }
	    for(int v2=0; v2<nVarPairs; v2++){
	      h2[v2][bl][fs][cat][nResStatuses][pr]->Add(h2[v2][bl][fs][cat][rs][pr]);
	      if(pr==Data){
		g2DataX[v2][bl][fs][cat][nResStatuses].insert(g2DataX[v2][bl][fs][cat][nResStatuses].end(), g2DataX[v2][bl][fs][cat][rs].begin(), g2DataX[v2][bl][fs][cat][rs].end());
		g2DataY[v2][bl][fs][cat][nResStatuses].insert(g2DataY[v2][bl][fs][cat][nResStatuses].end(), g2DataY[v2][bl][fs][cat][rs].begin(), g2DataY[v2][bl][fs][cat][rs].end());
		g2DataEX[v2][bl][fs][cat][nResStatuses].insert(g2DataEX[v2][bl][fs][cat][nResStatuses].end(), g2DataEX[v2][bl][fs][cat][rs].begin(), g2DataEX[v2][bl][fs][cat][rs].end());
		g2DataEY[v2][bl][fs][cat][nResStatuses].insert(g2DataEY[v2][bl][fs][cat][nResStatuses].end(), g2DataEY[v2][bl][fs][cat][rs].begin(), g2DataEY[v2][bl][fs][cat][rs].end());
	      }
	    }
	  }
	}
      }
    }
  }
  

  //---------- Write histograms to a ROOT file
  string outFileName = string(Form("histos_plotDataVsMC_%.3ffb-1%s.root",lumi,(MERGE2E2MU?"_m":"")));
  cout<<"Writing MC histograms and data graphs into file "<<outFileName<<" ..."<<endl;
  TFile* fOutHistos = new TFile(outFileName.c_str(),"recreate");
  fOutHistos->cd();
  for(int bl=0; bl<nBlindings; bl++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      for(int cat=0; cat<nCategories+1; cat++){
	for(int rs=0; rs<nResStatuses+1; rs++){
	  for(int pr=0; pr<nProcesses; pr++){
	    if( (fs==nFinalStates || pr==DY) &&
		cat==nCategories &&
		rs==nResStatuses ){
	      for(int v=0; v<nVariables; v++){
		h1[v][bl][fs][cat][rs][pr]->Write(h1[v][bl][fs][cat][rs][pr]->GetName());
		delete h1[v][bl][fs][cat][rs][pr];
	      }
	      for(int v2=0; v2<nVarPairs; v2++){
		h2[v2][bl][fs][cat][rs][pr]->Write(h2[v2][bl][fs][cat][rs][pr]->GetName());
		delete h2[v2][bl][fs][cat][rs][pr];
	      }
	    }
	  }
	  if( cat==nCategories &&
	      rs==nResStatuses ){
	    for(int v2=0; v2<nVarPairs; v2++){
	      if(MERGE2E2MU && fs==fs2mu2e) continue;
	      g2Data[v2][bl][fs][cat][rs] = new TGraphErrors( g2DataX[v2][bl][fs][cat][rs].size(),
							      &(g2DataX[v2][bl][fs][cat][rs][0]),
							      &(g2DataY[v2][bl][fs][cat][rs][0]),
							      &(g2DataEX[v2][bl][fs][cat][rs][0]),
							      &(g2DataEY[v2][bl][fs][cat][rs][0]) );
	      g2Data[v2][bl][fs][cat][rs]->Write(Form("g2Data_%s_%s_%s_%s_%s",varPairName[v2].c_str(),sBlinding[bl].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str(),sResonantStatus[rs].c_str()));
	      delete g2Data[v2][bl][fs][cat][rs];
	    }
	  }
	}
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


TH1F* getM4lZPlusXHistogram_Run2CombinedShape(Int_t nbins, Int_t xmin, Int_t xmax, Int_t finalState, double lumi) {

  //----- take the Z+X shapes sent by Pedja on March 1st

  TF1 *f4eComb = new TF1("f4eComb", "landau(0)*(1 + exp( pol1(3))) + [5]*(TMath::Landau(x, [6], [7]))", 70, 1000);
  TF1 *f4muComb = new TF1("f4muComb","landau(0)",70,1000);
  TF1 *f2e2muComb = new TF1("f2e2muComb","landau(0)",70,1000);

  f4eComb->SetParameters(4.404e-05,151.2,36.6,7.06,-0.00497,0.01446,157.3,26.00);
  f4muComb->SetParameters(0.04276,134.6,24.4);
  f2e2muComb->SetParameters(0.04130,144.5,25.3);

  //----- compute normalization of the subrange of interest

  Int_t nentries = 1000000;

  TH1F* hFullRange4mu   = new TH1F("hFullRange4mu"  ,";;",1000,0,1000);
  TH1F* hFullRange4e    = new TH1F("hFullRange4e"   ,";;",1000,0,1000);
  TH1F* hFullRange2e2mu = new TH1F("hFullRange2e2mu",";;",1000,0,1000);
  hFullRange4mu  ->FillRandom("f4muComb"  ,nentries);
  hFullRange4e   ->FillRandom("f4eComb"   ,nentries);
  hFullRange2e2mu->FillRandom("f2e2muComb",nentries);
  Float_t norm4mu   = normZPlusXFullSR4mu   * hFullRange4mu  ->Integral(hFullRange4mu  ->FindBin(xmin), hFullRange4mu  ->FindBin(xmax)-1) / hFullRange4mu  ->Integral();
  Float_t norm4e    = normZPlusXFullSR4e    * hFullRange4e   ->Integral(hFullRange4e   ->FindBin(xmin), hFullRange4e   ->FindBin(xmax)-1) / hFullRange4e   ->Integral();
  Float_t norm2e2mu = normZPlusXFullSR2e2mu * hFullRange2e2mu->Integral(hFullRange2e2mu->FindBin(xmin), hFullRange2e2mu->FindBin(xmax)-1) / hFullRange2e2mu->Integral();
  delete hFullRange4mu;
  delete hFullRange4e;
  delete hFullRange2e2mu;

  //*
  cout<<"In function getM4lZPlusXHistogram_Run2CombinedShape, xmin="<<xmin<<", xmax="<<xmax<<","<<endl;
  cout<<"  yield in 4mu   = "<<norm4mu<<endl;
  cout<<"  yield in 4e    = "<<norm4e<<endl;
  cout<<"  yield in 2e2mu = "<<norm2e2mu<<endl;
  cout<<"  total yield in 4l = "<<norm4mu+norm4e+norm2e2mu<<endl;
  //*/

  //----- build final histograms

  TH1F* h4mu   = new TH1F("h4mu"  ,";;",nbins,xmin,xmax);
  TH1F* h4e    = new TH1F("h4e"   ,";;",nbins,xmin,xmax);
  TH1F* h2e2mu = new TH1F("h2e2mu",";;",nbins,xmin,xmax);
  h4mu  ->FillRandom("f4muComb"  ,nentries);
  h4e   ->FillRandom("f4eComb"   ,nentries);
  h2e2mu->FillRandom("f2e2muComb",nentries);
  h4mu  ->Scale(norm4mu   / h4mu  ->Integral());
  h4e   ->Scale(norm4e    / h4e   ->Integral());
  h2e2mu->Scale(norm2e2mu / h2e2mu->Integral());
  TH1F* h4l = new TH1F("h4l",";;",nbins,xmin,xmax);
  h4l->Add(h4mu  );
  h4l->Add(h4e   );
  h4l->Add(h2e2mu);

  //----- return

  if(finalState==fs4mu) return h4mu;
  else if(finalState==fs4e) return h4e;
  else if(finalState==fs2e2mu) return h2e2mu;
  else if(finalState==nFinalStates) return h4l;
  else{
    cout<<"ERROR computing Z+X histogram: wrong final state"<<endl;
    abort();
  } 

}

/* //for GianLuca's FR files
TH1F* h1D_FRmu_EB = 0;
TH1F* h1D_FRmu_EE = 0;
TH1F* h1D_FRel_EB = 0;
TH1F* h1D_FRel_EE = 0;

Float_t fakeRate13TeV(Float_t LepPt, Float_t LepEta, Int_t LepID) {
  Float_t myLepPt = LepPt>=80. ? 79. : LepPt;
  Int_t   myLepID = abs(LepID);
  if(myLepID==11){
    if(fabs(LepEta)<1.449)
      return h1D_FRel_EB->GetBinContent(h1D_FRel_EB->GetXaxis()->FindBin(myLepPt));
    else
      return h1D_FRel_EE->GetBinContent(h1D_FRel_EE->GetXaxis()->FindBin(myLepPt));
  }else if(myLepID==13){
    if(fabs(LepEta)<1.2)
      return h1D_FRmu_EB->GetBinContent(h1D_FRmu_EB->GetXaxis()->FindBin(myLepPt));
    else
      return h1D_FRmu_EE->GetBinContent(h1D_FRmu_EE->GetXaxis()->FindBin(myLepPt));
  }else{
    cout<<"ERROR! wrong lepton ID : "<<myLepID<<endl;
    return 0.;
  }  
}
//*/

//* //for Christophe's FR files
TGraph* gr_FRmu_EB = 0;
TGraph* gr_FRmu_EE = 0;
TGraph* gr_FRel_EB = 0;
TGraph* gr_FRel_EE = 0;

Float_t fakeRate13TeV(Float_t LepPt, Float_t LepEta, Int_t LepID) {
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
    if(fabs(LepEta)<1.449)
      return (gr_FRel_EB->GetY())[bin];
    else
      return (gr_FRel_EE->GetY())[bin];
  }else if(myLepID==13){
    if(fabs(LepEta)<1.2)
      return (gr_FRmu_EB->GetY())[bin];
    else
      return (gr_FRmu_EE->GetY())[bin];
  }else{
    cout<<"ERROR! wrong lepton ID : "<<myLepID<<endl;
    return 0.;
  }  
}
//*/

void doHistogramsZPlusXSS(string inputFileAllData, string inputFileFakeRates, double lumi)
{

  cout<<"Preparing Z+X histograms (SS method) from "<<inputFileAllData<<endl;

  TFile* fFakeRates = TFile::Open(inputFileFakeRates.c_str());

  /* //for GianLuca's FR files
  h1D_FRmu_EB = (TH1F*)fFakeRates->Get("NoWZ_h1D_FRmu_EB");
  h1D_FRmu_EE = (TH1F*)fFakeRates->Get("NoWZ_h1D_FRmu_EE");
  h1D_FRel_EB = (TH1F*)fFakeRates->Get("NoWZ_h1D_FRel_EB");
  h1D_FRel_EE = (TH1F*)fFakeRates->Get("NoWZ_h1D_FRel_EE");
  //*/

  //* //for Christophe's FR files
  gr_FRmu_EB = (TGraph*)fFakeRates->Get("FR_SS_muon_EB");
  gr_FRmu_EE = (TGraph*)fFakeRates->Get("FR_SS_muon_EE");
  gr_FRel_EB = (TGraph*)fFakeRates->Get("FR_SS_electron_EB");
  gr_FRel_EE = (TGraph*)fFakeRates->Get("FR_SS_electron_EE");
  //*/

  Int_t nRun;
  Long64_t nEvent;
  Int_t nLumi;
  Int_t CRflag;
  Short_t ZZsel;
  Float_t ZZMass;
  Float_t ZZMassRefit;
  Float_t ZZPt;
  Float_t ZZEta;
  Float_t p0plus_VAJHU;
  Float_t bkg_VAMCFM;
  Float_t pvbf_VAJHU_highestPTJets;
  Float_t phjj_VAJHU_highestPTJets;
  Float_t pAux_vbf_VAJHU;
  Float_t phj_VAJHU;
  Float_t pwh_hadronic_VAJHU;
  Float_t pzh_hadronic_VAJHU;
  Float_t Z1Mass;
  Float_t Z2Mass;
  Short_t Z1Flav;
  Short_t Z2Flav;
  Float_t DiJetFisher;
  vector<Float_t> *LepPt = 0;
  vector<Float_t> *LepEta = 0;
  vector<Float_t> *LepPhi = 0;
  vector<Short_t> *LepLepId = 0;
  Short_t nExtraLep;
  Short_t nExtraZ;
  Short_t nJets;
  Short_t nJetsBTagged;
  vector<Float_t> *JetQGLikelihood = 0;
  Float_t jetPgOverPq[99];

  TFile* dataFile = TFile::Open(inputFileAllData.c_str());
  TTree* mytree = (TTree*)dataFile->Get("CRZLLTree/candTree");

  mytree->SetBranchAddress("RunNumber", &nRun);
  mytree->SetBranchAddress("EventNumber", &nEvent);
  mytree->SetBranchAddress("LumiNumber", &nLumi);
  mytree->SetBranchAddress("CRflag", &CRflag);
  mytree->SetBranchAddress("ZZsel", &ZZsel);
  mytree->SetBranchAddress("ZZMass", &ZZMass);
  mytree->SetBranchAddress("ZZMassRefit", &ZZMassRefit);
  mytree->SetBranchAddress("ZZPt", &ZZPt);
  mytree->SetBranchAddress("ZZEta", &ZZEta);
  mytree->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU);
  mytree->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM);
  mytree->SetBranchAddress("pvbf_VAJHU_highestPTJets", &pvbf_VAJHU_highestPTJets);
  mytree->SetBranchAddress("phjj_VAJHU_highestPTJets", &phjj_VAJHU_highestPTJets);
  mytree->SetBranchAddress("pAux_vbf_VAJHU", &pAux_vbf_VAJHU);
  mytree->SetBranchAddress("phj_VAJHU", &phj_VAJHU);
  mytree->SetBranchAddress("pwh_hadronic_VAJHU", &pwh_hadronic_VAJHU);
  mytree->SetBranchAddress("pzh_hadronic_VAJHU", &pzh_hadronic_VAJHU);
  mytree->SetBranchAddress("Z1Mass", &Z1Mass);
  mytree->SetBranchAddress("Z2Mass", &Z2Mass);
  mytree->SetBranchAddress("Z1Flav", &Z1Flav);
  mytree->SetBranchAddress("Z2Flav", &Z2Flav);
  mytree->SetBranchAddress("DiJetFisher", &DiJetFisher);
  mytree->SetBranchAddress("LepPt", &LepPt);
  mytree->SetBranchAddress("LepEta", &LepEta);
  mytree->SetBranchAddress("LepPhi", &LepPhi);
  mytree->SetBranchAddress("LepLepId", &LepLepId);
  mytree->SetBranchAddress("nExtraLep", &nExtraLep);
  mytree->SetBranchAddress("nExtraZ", &nExtraZ);
  mytree->SetBranchAddress("nCleanedJetsPt30", &nJets);
  mytree->SetBranchAddress("nCleanedJetsPt30BTagged", &nJetsBTagged);
  mytree->SetBranchAddress("JetQGLikelihood", &JetQGLikelihood);

  TH1F* h1[nVariables][nBlindings][nFinalStates+1];
  for(int bl=0; bl<nBlindings; bl++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      for(int v=0; v<nVariables; v++){
	h1[v][bl][fs] = new TH1F(
           Form("h1_ZPlusXSS_%s_%s_%s",varName[v].c_str(),sBlinding[bl].c_str(),sFinalState[fs].c_str()),
	   Form(";%s;%s",varXLabel[v].c_str(),varYLabel[v].c_str()),
	   varNbin[v],varMin[v],varMax[v]);
      }
    }
  }

  Float_t expectedYieldSR[nFinalStates+1];
  Int_t NumberOfEventsCR[nFinalStates+1];
  for(int fs=0; fs<nFinalStates+1; fs++){
    expectedYieldSR[fs] = 0.;
    NumberOfEventsCR[fs] = 0.;
  }

  Int_t currentFinalState;

  //---------- Process tree

  for(Long64_t iEvt=0; iEvt<mytree->GetEntries(); ++iEvt){

    mytree->GetEntry(iEvt);

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

    //----- fill histograms, update counters

    for(int j=0; j<nJets; j++) jetPgOverPq[j] = 1./JetQGLikelihood->at(j) - 1.;
    Float_t KD = p0plus_VAJHU / ( p0plus_VAJHU + bkg_VAMCFM ) ;
    Float_t vbfMela = pvbf_VAJHU_highestPTJets / ( phjj_VAJHU_highestPTJets + pvbf_VAJHU_highestPTJets );
    Float_t DCombVbf2j = (nJets>=2) ? 1/(1+ phjj_VAJHU_highestPTJets/pvbf_VAJHU_highestPTJets * TMath::Power(jetPgOverPq[0]*jetPgOverPq[1],1/3.) ) : -2 ;
    Float_t varVal[nVariables] = {
      ZZMass,
      ZZMass,
      ZZMass,
      ZZMassRefit,
      ZZMass,
      ZZMass,
      ZZMass,
      ZZMass,
      ZZMass,
      ZZMass,
      ZZMass,
      ZZMass,
      ZZMass,
      ZZMass,
      Z1Mass,
      Z1Mass,
      Z1Mass,
      Z1Mass,
      Z2Mass,
      Z2Mass,
      Z2Mass,
      Z2Mass,
      KD,
      KD,
      DiJetFisher,
      DCombVbf2j,//vbfMela,
      DCombVbf2j,//vbfMela,
      DCombVbf2j,//vbfMela,
      ZZPt,
      ZZEta,
      (Float_t)nExtraLep,
      (Float_t)nJets,
      (Float_t)nJetsBTagged,
      ZZMass,
      ZZMass,
    };
    Bool_t varPassCut[nVariables] = {
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,118<=ZZMass&&ZZMass<=130,1,1,1,118<=ZZMass&&ZZMass<=130,1,1,118<=ZZMass&&ZZMass<=130,nJets>=2,nJets>=2,nJets>=2&&118<=ZZMass&&ZZMass<=130,nJets>=2&&118<=ZZMass&&ZZMass<=130,1,1,1,1,1,KD>0.5,KD>0.5,
    };
    
    bool varPassBlindingCut[nBlindings] = {
      true,
      ZZMass<110,
      ZZMass>150,
      ZZMass<110 || ZZMass>150,
      ZZMass<110 || (ZZMass>150&&ZZMass<500),
      true,
    };

    Float_t yieldSR = fsROSSS[currentFinalState] * fakeRate13TeV(LepPt->at(2),LepEta->at(2),LepLepId->at(2)) * fakeRate13TeV(LepPt->at(3),LepEta->at(3),LepLepId->at(3));

    for(int bl=0; bl<nBlindings; bl++){
      for(int v=0; v<nVariables; v++){
	if( varName[v].find("M4l")==0 || (varName[v].find("M4l")!=0 && varPassBlindingCut[bl]) ){
	  if(varPassCut[v]){
	    if(MERGE2E2MU && currentFinalState==fs2mu2e)
	      h1[v][bl][fs2e2mu]->Fill(varVal[v], yieldSR);
	    else
	      h1[v][bl][currentFinalState]->Fill(varVal[v], yieldSR);
	  }
	}
      }
    }

    expectedYieldSR[currentFinalState] += yieldSR;
    NumberOfEventsCR[currentFinalState]++;
		
  }

  //---------- Fill 'inclusive' histograms
  for(int bl=0; bl<nBlindings; bl++){
    for(int fs=0; fs<nFinalStates; fs++){
      for(int v=0; v<nVariables; v++){
	h1[v][bl][nFinalStates]->Add(h1[v][bl][fs]);
      }
    }
  }
  for(int fs=0; fs<nFinalStates; fs++){
    expectedYieldSR[nFinalStates] += expectedYieldSR[fs];
    NumberOfEventsCR[nFinalStates] += NumberOfEventsCR[fs];
  }

  //---------- Print Z+X expected yields
  for(int fs=0; fs<nFinalStates; fs++){
    cout<<fsLabelForSS[fs]<<" : "
	<<expectedYieldSR[fs]
	<<" +/- "<<expectedYieldSR[fs]/sqrt(NumberOfEventsCR[fs])<<" (stat., evt: "<<NumberOfEventsCR[fs]<<")" 
	<<" +/- "<<expectedYieldSR[fs]*0.50<< " (syst.)" 
	<<endl;
  }
  cout<<"Total: "<<expectedYieldSR[nFinalStates]<<endl;

  //---------- Smooth histogram if requested
  if(SMOOTHZPLUSXFULLRUN2SS){
    Float_t integral = 0;
    for(int bl=0; bl<nBlindings; bl++){
      for(int fs=0; fs<nFinalStates+1; fs++){
	for(int v=0; v<nVariables; v++){
	  integral = h1[v][bl][fs]->Integral();
	  h1[v][bl][fs]->Smooth(1);
	  h1[v][bl][fs]->Scale( integral / h1[v][bl][fs]->Integral() );
	}
      }
    }
  }

  //---------- Normalize to official Z+X yield (combination of OS and SS)
  if(RENORMALIZEZPLUSXFULLRUN2SS){
    if(MERGE2E2MU){
      Float_t normZPlusXFullSR[nFinalStates+1];
      normZPlusXFullSR[fs4e] = normZPlusXFullSR4e;
      normZPlusXFullSR[fs4mu] = normZPlusXFullSR4mu;
      normZPlusXFullSR[fs2e2mu] = normZPlusXFullSR2e2mu;
      normZPlusXFullSR[nFinalStates] = normZPlusXFullSR4e + normZPlusXFullSR4mu + normZPlusXFullSR2e2mu ;
      expectedYieldSR[fs2e2mu] += expectedYieldSR[fs2mu2e];//(imposed MERGE2E2MU)
      for(int bl=0; bl<nBlindings; bl++){
	for(int fs=0; fs<nFinalStates+1; fs++){
	  if(fs==fs2mu2e) continue;
	  for(int v=0; v<nVariables; v++){
	    h1[v][bl][fs]->Scale( normZPlusXFullSR[fs] / expectedYieldSR[fs] );
	  }
	}
      }
    }else{
      cout<<"WARNING: cannot renormalize Z+X histograms to official numbers when treating 2e2mu and 2mu2e separately"<<endl;
    }
  }

  //---------- Write histograms to a ROOT file
  string outFileName = string(Form("histos_plotDataVsMC_ZPlusXSS_%.3ffb-1%s.root",lumi,(MERGE2E2MU?"_m":"")));
  cout<<"Writing Z+X histograms into file "<<outFileName<<" ..."<<endl;
  TFile* fOutHistos = new TFile(outFileName.c_str(),"recreate");
  fOutHistos->cd();
  for(int bl=0; bl<nBlindings; bl++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      //if(fs==nFinalStates){
	for(int v=0; v<nVariables; v++){
	  h1[v][bl][fs]->Write(h1[v][bl][fs]->GetName());
	  delete h1[v][bl][fs];
	}
	//}
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


void DrawDataMC(TCanvas* c, TH1F** h, TH1F* hZPX, int v, int bl, double lumi, string lumiText, Bool_t large, Bool_t logX = false, Bool_t logY = false) {

  bool drawData = !( (MASKDATAFORHIGHMASS && bl==blindbelow150) || bl==fullyblind );
  bool maskH125 = 
    ( MASKH125FORLOWMASS && bl==blindabove110 && (varName[v].find("M4l")==string::npos || varName[v].find("M4l_70110")==0 || varName[v].find("M4l_70109")==0) ) ||
    ( MASKH125FORHIGHMASS && bl==blindbelow150 && (varName[v].find("M4l")==string::npos || varName[v].find("M4l_above150")==0) ) ;
  bool withRatioPlot = DRAWDATAMCRATIO && drawData;
  bool doBlindingHisto = xHistoBlindLow[bl]<xHistoBlindUp[bl] && varName[v].find("M4l")==0;
  bool doBlindingHisto2 = xHistoBlind2Low[bl]<xHistoBlind2Up[bl] && varName[v].find("M4l")==0;

  //----- prepare canvas
  TStyle* style = (TStyle*)gROOT->GetStyle("tdrStyle")->Clone();
  style->SetPadLeftMargin(0.14);
  style->SetPadRightMargin(0.04);
  style->cd();
  c->cd();
  c->UseCurrentStyle();
  if(logX) c->SetLogx();
  if(logY) c->SetLogy();

  //----- prepare MC histograms
  TH1F* hStacks[nProcesses];
  bool first = true;
  int previous = -1;
  for(int pr=nProcesses-1; pr>=1; pr--){
    if((useProcess[pr] && !((maskH125 || separateVbf[v]) && pr==H125)) || (separateVbf[v] && (pr==H125VBF || pr==H125NONVBF))){
      if(REBINDYANDTTBAR && (pr==DY || pr==ttbar)) h[pr] = Smooth(h[pr],rebinningDYTTbar[v]);
      hStacks[pr] = (TH1F*)h[pr]->Clone();
      if(!first) hStacks[pr]->Add(hStacks[previous]);
      hStacks[pr]->SetFillColor(processFillColor[pr]);
      hStacks[pr]->SetLineColor(DRAWLINES?processLineColor[pr]:processFillColor[pr]);
      hStacks[pr]->SetLineWidth(LINEWIDTH);
      if(!DRAWLINES) hStacks[pr]->SetLineColorAlpha(processFillColor[pr],0.);
      hStacks[pr]->GetXaxis()->SetTitleOffset(large?1.15:1.2);
      hStacks[pr]->GetYaxis()->SetTitleOffset(large?1.:1.4);
      if(large){
	hStacks[pr]->GetXaxis()->SetTitleSize(0.05);
	hStacks[pr]->GetYaxis()->SetTitleSize(0.05);
      }
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

  //----- prepare reducible background
  TH1F* hZPlusX = 0;
  Bool_t useZPlusX = USEZPLUSXFULLRUN2SS ||
    (USEZPLUSXRUN2COMBINEDSHAPE && (varName[v].find("M4l")==0 && varName[v].find("Refit")==string::npos && varCutLabel[v]=="") );
  useZPlusX = useZPlusX && hZPX->Integral()>0.;
  if(useZPlusX){
    hZPlusX = hZPX;
    if(logY) hZPlusX->SetMinimum(1e-20);
    hZPlusX->SetFillColor(processFillColor[DY]);
    hZPlusX->SetLineColor(DRAWLINES?processLineColor[DY]:processFillColor[DY]);
    hZPlusX->SetLineWidth(LINEWIDTH);
    for(int pr=nProcesses-1; pr>=1; pr--)
      if((useProcess[pr] && !((maskH125 || separateVbf[v]) && pr==H125)) || (separateVbf[v] && (pr==H125VBF || pr==H125NONVBF)))
	hStacks[pr]->Add(hZPlusX);
    cout<<"full expected yield for variable "<<varName[v]<<": "<<hStacks[idxSumMC]->Integral()<<endl;
    cout<<"Z+X yield for variable "<<varName[v]<<": "<<hZPlusX->Integral()<<endl;
  }

  //----- prepare data graph
  h[0]->SetMarkerStyle(20);
  h[0]->SetMarkerColor(kBlack);
  h[0]->SetMarkerSize(1.);
  TGraphAsymmErrors* gData = getDataGraph(h[0]);
  if(h[0]->GetNbinsX()>60)  gData->SetMarkerSize(0.8);
  if(h[0]->GetNbinsX()>100) gData->SetMarkerSize(0.7);

  //----- adjust Y axis
  const int npoints = gData->GetN();
  Float_t gDataErrorBarUp[npoints];
  for(int i=0; i<npoints; i++) gDataErrorBarUp[i] = gData->GetY()[i] + gData->GetEYhigh()[i] ;
  Float_t cmax;
  if(drawData) cmax = TMath::Max( (Float_t)hStacks[idxSumMC]->GetMaximum(), (Float_t)TMath::MaxElement(npoints,gDataErrorBarUp) );
  else cmax = (Float_t)hStacks[idxSumMC]->GetMaximum();
  cmax *= logY ? 30. : 1.1 ;
  cmax *= varMaxCorrector[bl][v];
  Float_t cminlog = hStacks[idxSumMC]->GetMaximum() / varMinFactor[v];
  hStacks[idxSumMC]->SetMaximum(cmax);
  hStacks[idxSumMC]->SetMinimum(logY?cminlog:0.);
  
  //----- prepare grey area/grid for blind region
  TH1F* hBlind = new TH1F(Form("hBlind_%s_%s",varName[v].c_str(),sBlinding[bl].c_str()),";;",1,xHistoBlindLow[bl],xHistoBlindUp[bl]);
  if(doBlindingHisto){
    hBlind->SetBinContent(1,cmax);
    //hBlind->SetFillColorAlpha(kBlack,0.2);
    //hBlind->SetLineColorAlpha(kBlack,0.2);
    hBlind->SetFillColor(kGray+3);
    hBlind->SetFillStyle(3013); //also tried 3001 and a few others, but they look bad in pdf format
    hBlind->SetLineColorAlpha(kWhite,0.);
  }
  TH1F* hBlind2 = new TH1F(Form("hBlind2_%s_%s",varName[v].c_str(),sBlinding[bl].c_str()),";;",1,xHistoBlind2Low[bl],xHistoBlind2Up[bl]);
  if(doBlindingHisto2){
    hBlind2->SetBinContent(1,cmax);
    hBlind2->SetFillColor(kGray+3);
    hBlind2->SetFillStyle(3013);
    hBlind2->SetLineColorAlpha(kWhite,0.);
  }

  //----- prepare legend
  bool useBlindingLabel = blindingLabel[bl]!="" && varName[v].find("M4l")==string::npos;
  float legTextSize = 0.034;
  int legPos = varLegPos[v];
  if(bl==blindabove110 && varName[v].find("MZ")!=string::npos) legPos = 33;
  float legLeft = (legPos==11) ? 0.20 : 0.68 ;
  float legWidth = 0.19;
  float legUp = 0.9;
  if(!withRatioPlot && legPos==varCMSPos[v] && !large) legUp -= 0.1;
  float legHeight = 0.16;
  if(!drawData) legHeight -= 0.04;
  if(maskH125) legHeight -= 0.04;
  if(useZPlusX) legHeight += 0.04;
  if(separateVbf[v]) legHeight += 0.04;
  if(withRatioPlot) legHeight /= 0.72;
  if(large){ legWidth *= 1.15; legHeight *= 1.15; legTextSize *= 1.15; }
  TLegend* lgd = new TLegend(legLeft,legUp-legHeight,legLeft+legWidth,legUp);
  lgd->SetFillStyle(0);
  lgd->SetBorderSize(0);
  lgd->SetTextSize(legTextSize);
  if(drawData)
    lgd->AddEntry(h[0],processLabel[0].c_str(),"p");
  for(int pr=1; pr<nProcesses; pr++)
    if((useProcess[pr] && !((maskH125 || separateVbf[v]) && pr==H125)) || (separateVbf[v] && (pr==H125VBF || pr==H125NONVBF)))
      lgd->AddEntry(hStacks[pr],processLabel[pr].c_str(),"f");
  if(useZPlusX)
    lgd->AddEntry(hZPlusX," Z+X","f");
    

  //----- prepare label for cut/blinding
  bool useLabel = useBlindingLabel || varCutLabel[v]!="";
  TPaveText* pav = 0;
  if(useLabel){
    pav = new TPaveText(legLeft-(legPos==11?0.02:0.02),legUp-legHeight-0.1,legLeft+legWidth+(legPos==11?0.06:0.06),legUp-legHeight+0.02,"brNDC");
    pav->SetFillStyle(0);
    pav->SetBorderSize(0);
    pav->SetTextAlign(legPos==11?13:33);
    pav->SetTextSize(0.037);
    pav->AddText((varCutLabel[v]!="")?varCutLabel[v].c_str():blindingLabel[bl].c_str());
  }

  //----- draw everything
  if(!withRatioPlot){

    //--- no Data/MC graph -> draw on main pad
    for(int pr=1; pr<nProcesses; pr++)
      if((useProcess[pr] && !((maskH125 || separateVbf[v]) && pr==H125)) || (separateVbf[v] && (pr==H125VBF || pr==H125NONVBF)))
	hStacks[pr]->Draw( pr==idxSumMC ? "HIST" : "HIST SAME" );
    if(useZPlusX) hZPlusX->Draw("HIST SAME");
    if(drawData) gData->Draw("P");
    if(doBlindingHisto) hBlind->Draw("HIST SAME");
    if(doBlindingHisto2) hBlind2->Draw("HIST SAME");
    lgd->Draw();
    if(useLabel) pav->Draw();
    gPad->RedrawAxis();

  }else{

    //--- with Data/MC graph -> need 2 pads
    TPad* pad1 = new TPad("pad1", "pad1", 0., 0.28, 1., 0.95);
    TPad* pad2 = new TPad("pad2", "pad2", 0., 0.13, 1., 0.28);
    if(logX){pad1->SetLogx(); pad2->SetLogx();}
    if(logY) pad1->SetLogy();
    pad1->SetMargin(0.14,0.04,0.03,0.);
    pad2->SetMargin(0.14,0.04,0.,0.);

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
      if((useProcess[pr] && !((maskH125 || separateVbf[v]) && pr==H125)) || (separateVbf[v] && (pr==H125VBF || pr==H125NONVBF))){
	hStacks[pr]->GetYaxis()->SetTitleOffset(1.);
	hStacks[pr]->GetYaxis()->SetTitleSize(0.06);
	hStacks[pr]->GetYaxis()->SetLabelSize(0.06);
	hStacks[pr]->Draw( pr==idxSumMC ? "HIST" : "HIST SAME" );
      }
    }
    if(useZPlusX) hZPlusX->Draw("HIST SAME");
    if(bl!=fullyblind) gData->Draw("P");
    hStacks[idxSumMC]->GetXaxis()->SetLabelSize(0.);
    if(doBlindingHisto) hBlind->Draw("HIST SAME");
    if(doBlindingHisto2) hBlind2->Draw("HIST SAME");
    lgd->Draw();
    if(useLabel) pav->Draw();
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
    if(doBlindingHisto){
      TH1F* hBlindDataMC = (TH1F*)hBlind->Clone();
      hBlindDataMC->SetBinContent(1,2.);
      hBlindDataMC->Draw("HIST SAME");
    }
    if(doBlindingHisto2){
      TH1F* hBlindDataMC2 = (TH1F*)hBlind2->Clone();
      hBlindDataMC2->SetBinContent(1,2.);
      hBlindDataMC2->Draw("HIST SAME");
    }
    pad2->RedrawAxis();

  }

  //----- customize m4l axis labels
  if(DRAWLABELBYHAND && logX && (varName[v]=="M4lV1"||varName[v]=="M4lV2"||varName[v]=="M4lV2Refit"||varName[v]=="M4lV2b"||varName[v]=="M4lV3"||varName[v]=="M4lV4")){
    c->cd();
    TText t;
    t.SetNDC();
    t.SetTextSize(0.04);
    t.SetTextFont(42);
    t.DrawText(large?0.811:0.805,0.093,"600");
    t.DrawText(large?0.904:0.898,0.093,"800");
  }

  //----- print official CMS labels and lumi
  writeExtraText = true;
  extraText  = "Preliminary";
  lumi_sqrtS = lumiText + " (13 TeV)";
  CMS_lumi( c, 0, (withRatioPlot || large) ? 0 : varCMSPos[v] );

  //----- print yields
  //if(varName[v]=="MZ1"){
  cout<<"Yields for blinding "<<bl<<" and variable "<<varName[v]<<endl;
  cout<<"  qqZZ             "<<h[qqZZ]->Integral()<<endl;
  cout<<"  ggZZ             "<<h[ggZZ]->Integral()<<endl;
  if(useZPlusX) cout<<"  Z+X              "<<hZPlusX->Integral()<<endl;
  cout<<"  all backgrounds  "<<hStacks[idxSumMC]->Integral() - h[H125]->Integral()<<endl;
  cout<<"  signal           "<<h[H125]->Integral()<<endl;
  cout<<"  total expected   "<<hStacks[idxSumMC]->Integral()<<endl;
  cout<<"  observed         "<<h[0]->Integral()<<endl;
    //}

}


void DrawDataMC2D(TCanvas* c, TH2F** h2, TGraphErrors** g2, int v2, int bl, string lumiText, int style2D, Bool_t logX = false, Bool_t logY = false) {

  bool maskH125 = 
    ( MASKH125FORLOWMASS  && bl==blindabove110 && (varPairName[v2].find("M4l")==string::npos || varPairName[v2].find("M4L70110")!=string::npos) ) ||
    ( MASKH125FORHIGHMASS && bl==blindbelow150 && (varPairName[v2].find("M4l")==string::npos || varPairName[v2].find("M4L180780")!=string::npos || varPairName[v2].find("M4L150700")!=string::npos) ) ;

  //----- prepare canvas
  TStyle* style = (TStyle*)gROOT->GetStyle("tdrStyle")->Clone();
  style->SetPadLeftMargin(0.12);
  style->SetPadRightMargin(0.14);
  style->SetEndErrorSize(0.);
  style->cd();

  bool LEGENDWHITE = 0;
  bool signalSeparately = (style2D==7);
  if(style2D==0 || varPairName[v2].find("M4lVsM4lRefit")==0){
    setColZGradient_Rainbow2();
    LEGENDWHITE = varPairLegIsWhite[v2];
  }else if(style2D==8){
    setColZGradient_TwoColors();
  }else{
    setColZGradient_OneColor(style2D,varPairName[v2].find("MZ1VsMZ2")==0);
  }
  bool withMassErrors = (varPairName[v2].find("M4l")==0 && MARKERSTYLE==0);

  c->cd();
  //c->SetFrameBorderMode(0);
  c->UseCurrentStyle();
  if(logX) c->SetLogx();
  if(logY) c->SetLogy();

  //----- prepare MC histogram
  TH2F* h2Stacked = 0;
  bool first = true;
  for(int pr=nProcesses-1; pr>=1; pr--){
    if(useProcess[pr] && !((maskH125 || signalSeparately) && pr==H125)){
      if(pr==DY || pr==ttbar) cout<<"WARNING in function DrawDataMC2D: including "<<sProcess[pr]<<endl;
      if(first) h2Stacked = (TH2F*)h2[pr]->Clone();
      else h2Stacked->Add(h2[pr]);
      first = false;
    }
  }
  h2Stacked->GetXaxis()->SetTitleOffset(1.2);
  h2Stacked->GetYaxis()->SetTitleOffset(1.4);
  h2Stacked->GetZaxis()->SetTitleOffset(1.35);
  h2Stacked->GetXaxis()->SetLabelFont(42);
  h2Stacked->GetYaxis()->SetLabelFont(42);
  h2Stacked->GetZaxis()->SetLabelFont(42);
  h2Stacked->GetXaxis()->SetTitleFont(42);
  h2Stacked->GetYaxis()->SetTitleFont(42);
  h2Stacked->GetZaxis()->SetTitleFont(42);
  h2Stacked->GetZaxis()->SetLabelSize(0.03);
  h2Stacked->GetZaxis()->SetTitle("Events / bin");
  if(logX) h2Stacked->GetXaxis()->SetMoreLogLabels();
  if(logX) h2Stacked->GetXaxis()->SetNoExponent();
  if(varPairName[v2].find("MZ1VsMZ2")==string::npos) h2Stacked->SetMinimum(-1e-20); // avoid white bins
  else h2Stacked->SetMinimum(+1e-20);//(+1e-5);//

  TH2F* h2Signal = 0;
  if(signalSeparately){
    h2Signal = (TH2F*)h2[H125]->Clone();
    h2Signal->SetFillColor(TColor::GetColor("#ff9090"));
  }

  //----- prepare data graphs
  for(int fs=0; fs<nFinalStates; fs++){
    if(MERGE2E2MU && fs==fs2mu2e) continue;
    g2[fs]->SetMarkerStyle(MARKERSTYLE==1?fsMarkerStyleOpen[fs]:fsMarkerStyleFull[fs]);
    g2[fs]->SetMarkerSize(MARKERSTYLE==1?1.:0.7);
    if(style2D==1){
      g2[fs]->SetMarkerColor(fsMarkerColor[fs]);
      g2[fs]->SetLineColor(fsMarkerColor[fs]);
    }else{
      g2[fs]->SetMarkerColor(kBlack);
      g2[fs]->SetLineColor(kBlack);
    }
  }
  
  //----- prepare grid for blind region
  bool doBlindingArea = xHistoBlindLow[bl]<xHistoBlindUp[bl] 
    && varPairName[v2].find("M4l")==0
    && xHistoBlindLow[bl]<h2Stacked->GetXaxis()->GetXmax()
    && xHistoBlindUp [bl]>h2Stacked->GetXaxis()->GetXmin(); 
  TBox* box = 0;
  if(doBlindingArea){
    box = new TBox(TMath::Max((Float_t)xHistoBlindLow[bl],(Float_t)h2Stacked->GetXaxis()->GetXmin()),h2Stacked->GetYaxis()->GetXmin(),TMath::Min((Float_t)xHistoBlindUp[bl],(Float_t)h2Stacked->GetXaxis()->GetXmax()),h2Stacked->GetYaxis()->GetXmax());
    box->SetFillColor(kBlack);
    box->SetFillStyle(3013);
  }
  bool doBlindingArea2 = xHistoBlind2Low[bl]<xHistoBlind2Up[bl] 
    && varPairName[v2].find("M4l")==0
    && xHistoBlind2Low[bl]<h2Stacked->GetXaxis()->GetXmax()
    && xHistoBlind2Up [bl]>h2Stacked->GetXaxis()->GetXmin(); 
  TBox* box2 = 0;
  if(doBlindingArea2){
    box2 = new TBox(TMath::Max((Float_t)xHistoBlind2Low[bl],(Float_t)h2Stacked->GetXaxis()->GetXmin()),h2Stacked->GetYaxis()->GetXmin(),TMath::Min((Float_t)xHistoBlind2Up[bl],(Float_t)h2Stacked->GetXaxis()->GetXmax()),h2Stacked->GetYaxis()->GetXmax());
    box2->SetFillColor(kBlack);
    box2->SetFillStyle(3013);
  }

  //----- prepare legend
  int legPos = varPairLegPos[v2];
  float legLeft = (legPos==11) ? 0.20 : 0.7 ;
  float legWidth = 0.12;
  float legUp = 0.92;
  float legHeight = 0.12;
  if(signalSeparately){
    legLeft = (legPos==11) ? 0.20 : 0.58 ;
    legWidth = 0.24;
    legHeight = 0.20;
  }
  TLegend* lgd = new TLegend(legLeft,legUp-legHeight,legLeft+legWidth,legUp);
  lgd->SetFillStyle(0);
  if(LEGENDWHITE){
    lgd->SetBorderSize(1);
    lgd->SetLineColor(kWhite);
    lgd->SetTextColor(kWhite);
  }else{
    lgd->SetBorderSize(1);
    lgd->SetLineColor(kBlack);
    lgd->SetTextColor(kBlack);
  }
  if(signalSeparately){
    lgd->AddEntry(h2Signal,processLabel[H125].c_str(),"f");
    h2Stacked->SetFillColor(TColor::GetColor("#74baff"));
    lgd->AddEntry(h2Stacked,"ZZ background","f");
  }
  TGraph* gLeg[nFinalStates];
  for(int fs=0; fs<nFinalStates; fs++){
    if(MERGE2E2MU && fs==fs2mu2e) continue;
    gLeg[fs] = (TGraph*)g2[fs]->Clone();
  }
  lgd->AddEntry(gLeg[fs4e   ],((signalSeparately?"Data, ":"")+fsLabel[fs4e   ]).c_str(),withMassErrors?"lp":"p");
  lgd->AddEntry(gLeg[fs4mu  ],((signalSeparately?"Data, ":"")+fsLabel[fs4mu  ]).c_str(),withMassErrors?"lp":"p");
  lgd->AddEntry(gLeg[fs2e2mu],((signalSeparately?"Data, ":"")+fsLabel[fs2e2mu]).c_str(),withMassErrors?"lp":"p");
  if(!MERGE2E2MU) lgd->AddEntry(gLeg[fs2mu2e],((signalSeparately?"Data, ":"")+fsLabel[fs2mu2e]).c_str(),withMassErrors?"lp":"p");

  //----- prepare label for cut/blinding
  bool useLabel = (blindingLabel[bl]!="" && varPairName[v2].find("M4l")==string::npos) || varPairCutLabel[v2]!="";
  TPaveText* pav = 0;
  if(useLabel){
    pav = new TPaveText(legLeft-(legPos==11?0.02:0.18),legUp-legHeight-0.1,legLeft+legWidth+(legPos==11?0.09:0.),legUp-legHeight+0.02,"brNDC");
    pav->SetFillStyle(0);
    pav->SetBorderSize(0);
    pav->SetTextAlign(13);
    pav->SetTextSize(0.034);
    pav->SetTextColor(LEGENDWHITE?kWhite:kBlack);
    pav->AddText((varPairCutLabel[v2]!="")?varPairCutLabel[v2].c_str():blindingLabel[bl].c_str());
  }

  //----- draw everything
  h2Stacked->Draw("COLZ");
  if(signalSeparately){
    h2Signal->Draw("BOX SAME");
    //readjust box size range:
    h2Signal->Scale(0.6);
    TH2F* hOne = (TH2F*)h2Signal->Clone();
    hOne->Reset();
    for(int i=1; i<=hOne->GetNbinsX()*(hOne->GetNbinsY()+4); i++) hOne->SetBinContent(i,-0.01*h2Signal->GetMaximum());
    h2Signal->Add(hOne);
  }
  if(bl!=fullyblind && varPairName[v2]!="M4lVsM4lRefit"){
    g2[fs2e2mu]->Draw(withMassErrors?"P":"XP");
    if(!MERGE2E2MU) g2[fs2mu2e]->Draw(withMassErrors?"P":"XP");
    g2[fs4mu]->Draw(withMassErrors?"P":"XP");
    g2[fs4e]->Draw(withMassErrors?"P":"XP");
  }
  if(doBlindingArea) box->Draw();
  if(doBlindingArea2) box2->Draw();
  if(bl!=fullyblind) lgd->Draw();
  if(useLabel) pav->Draw();
  gPad->RedrawAxis();

  //----- adjust color axis
  c->Update();
  TPaletteAxis* pal = (TPaletteAxis*)h2Stacked->GetListOfFunctions()->FindObject("palette");
  pal->SetX1NDC(0.87);
  pal->SetX2NDC(0.895);

  //----- customize m4l axis labels
  if(DRAWLABELBYHAND && logX && varPairName[v2]=="M4lVsKD"){
    c->cd();
    TText t;
    t.SetNDC();
    t.SetTextSize(0.04);
    t.SetTextFont(42);
    t.DrawText(0.122,0.093,"110");
    t.DrawText(0.720,0.093,"600");
    t.DrawText(0.820,0.093,"800");
  }

  //----- print official CMS labels and lumi
  writeExtraText = true;
  extraText  = "Preliminary";
  lumi_sqrtS = lumiText + " (13 TeV)";
  CMS_lumi( c, 0, 0 );

}


void doPlots(string outputDirectory, int variableList, int varPairList, int blindingList, double lumi, string lumiText)
{

  setTDRStyle();
  gStyle->SetGridColor(kGray);

  //---------- retrieve histograms from the main ROOT file
  TH1F* h1[nVariables][nBlindings][nProcesses];
  TH2F* h2[nVarPairs ][nBlindings][nProcesses];
  TGraphErrors* g2Data[nVarPairs][nBlindings][nFinalStates];
  string inFileName = string(Form("histos_plotDataVsMC_%.3ffb-1%s.root",lumi,(MERGE2E2MU?"_m":"")));
  cout<<"Retrieving MC histograms and data graphs from file "<<inFileName<<" ..."<<endl;
  TFile* fInHistos = TFile::Open(inFileName.c_str());
  for(int bl=0; bl<nBlindings; bl++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      for(int cat=0; cat<nCategories+1; cat++){
	for(int rs=0; rs<nResStatuses+1; rs++){
	  if( fs==FINALSTATE &&
	      cat==nCategories &&
	      rs==nResStatuses ){
	    for(int pr=0; pr<nProcesses; pr++){
	      for(int v=0; v<nVariables; v++){
		h1[v][bl][pr] = (TH1F*)fInHistos->Get( Form("h1_%s_%s_%s_%s_%s_%s",varName[v].c_str(),sBlinding[bl].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str(),sResonantStatus[rs].c_str(),sProcess[pr].c_str()) );
		h1[v][bl][pr]->GetXaxis()->SetTitle(ReplaceString(h1[v][bl][pr]->GetXaxis()->GetTitle(),"4#font[12]{l}",fsLabel[fs]).c_str());
		if(RESCALETOSMPSIGNALSTRENGTH && (pr==qqZZ||pr==ggZZ)) h1[v][bl][pr]->Scale(SMPSIGNALSTRENGTH); //FIXME: to be removed at some point ?
		if(restrictCountVar[v]) restrictXAxis(h1[v][bl][pr],restrictCountVar[v]);
	      }
	      for(int v2=0; v2<nVarPairs; v2++){
		h2[v2][bl][pr] = (TH2F*)fInHistos->Get( Form("h2_%s_%s_%s_%s_%s_%s",varPairName[v2].c_str(),sBlinding[bl].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str(),sResonantStatus[rs].c_str(),sProcess[pr].c_str()) );
		h2[v2][bl][pr]->GetXaxis()->SetTitle(ReplaceString(h2[v2][bl][pr]->GetXaxis()->GetTitle(),"4#font[12]{l}",fsLabel[fs]).c_str());
		h2[v2][bl][pr]->GetYaxis()->SetTitle(ReplaceString(h2[v2][bl][pr]->GetYaxis()->GetTitle(),"4#font[12]{l}",fsLabel[fs]).c_str());
		if(RESCALETOSMPSIGNALSTRENGTH && (pr==qqZZ||pr==ggZZ)) h2[v2][bl][pr]->Scale(SMPSIGNALSTRENGTH);  //FIXME: to be removed at some point ?
	      }
	    }
	  }
	  if( cat==nCategories &&
	      rs==nResStatuses ){
	    for(int v2=0; v2<nVarPairs; v2++){
	      if(MERGE2E2MU && fs==fs2mu2e) continue;
	      g2Data[v2][bl][fs] = (TGraphErrors*)fInHistos->Get( Form("g2Data_%s_%s_%s_%s_%s",varPairName[v2].c_str(),sBlinding[bl].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str(),sResonantStatus[rs].c_str()) );
	    }
	  }
	}
      }
    }
  }

  //---------- get Z+X histograms
  TH1F* h1_ZPlusX[nVariables][nBlindings];
  if(USEZPLUSXFULLRUN2SS){
    string inFileName_ZPlusXSS = string(Form("histos_plotDataVsMC_ZPlusXSS_%.3ffb-1%s.root",lumi,(MERGE2E2MU?"_m":"")));
    cout<<"Retrieving Z+X histograms from file "<<inFileName_ZPlusXSS<<" ..."<<endl;
    TFile* fInHistos_ZPlusXSS = TFile::Open(inFileName_ZPlusXSS.c_str());
    for(int bl=0; bl<nBlindings; bl++){
      for(int fs=0; fs<nFinalStates+1; fs++){
	if(fs==FINALSTATE){
	  for(int v=0; v<nVariables; v++){
	    h1_ZPlusX[v][bl] = (TH1F*)fInHistos_ZPlusXSS->Get( Form("h1_ZPlusXSS_%s_%s_%s",varName[v].c_str(),sBlinding[bl].c_str(),sFinalState[fs].c_str()) );
	    if(restrictCountVar[v]) restrictXAxis(h1_ZPlusX[v][bl],restrictCountVar[v]);
	  }
	}
      }
    }
  }
  if(USEZPLUSXRUN2COMBINEDSHAPE){
    if(!MERGE2E2MU){
      cout<<"ERROR: cannot get Z+X histograms from combined shape for 2e2mu and 2mu2e separately"<<endl;
      abort();
    }
    for(int bl=0; bl<nBlindings; bl++){
      if(bl==fullyblind || bl==unblinded){
	for(int v=0; v<nVariables; v++){
	  if(varName[v].find("M4l")==0 && varName[v].find("Refit")==string::npos && varCutLabel[v]==""){ 
	    h1_ZPlusX[v][bl] = getM4lZPlusXHistogram_Run2CombinedShape(varNbin[v],varMin[v],varMax[v],FINALSTATE,lumi); // overrides the histogram from SS method if any
	    if(restrictCountVar[v]) restrictXAxis(h1_ZPlusX[v][bl],restrictCountVar[v]);
	  }
	}
      }else{
	cout<<"WARNING: cannot get Z+X histogram from combined shape for requested blinding strategy"<<endl;
      }
    }
  }

  //---------- do the plots (1 canvas per variable and per blinding policy)

  string canvasName;
  TCanvas* c1[nVariables];
  if(DO1DPLOTS){
    cout<<"Doing 1D plots ..."<<endl;
    for(int bl=0; bl<nBlindings; bl++){
      if(!plotThisBlinding[blindingList][bl]) continue;
      for(int v=0; v<nVariables; v++){
	if(!plotThisVar[variableList][v]) continue;
	if(varName[v].find("M4L118130")!=string::npos && bl!=unblinded && bl!=fullyblind) continue;
	bool large = varName[v]=="M4lV2"||varName[v]=="M4lV2Refit"||varName[v]=="M4lV2b";//0;//
	gStyle->SetFrameLineWidth(large?2:1);
	canvasName = string(Form("c_%s_%s_%s",sBlinding[bl].c_str(),varName[v].c_str(),sFinalState[FINALSTATE].c_str()));
	c1[v] = new TCanvas(canvasName.c_str(),canvasName.c_str(),large?650:500,500);
	DrawDataMC(c1[v],h1[v][bl],h1_ZPlusX[v][bl],v,bl,lumi,lumiText,large,varLogx[v],varLogy[v]);
	SaveCanvas(outputDirectory,c1[v]);
      }
    }
  }


  //*
  TCanvas* c2[nVarPairs];
  if(DO2DPLOTS){
    cout<<"Doing 2D plots ..."<<endl;
    for(int bl=0; bl<nBlindings; bl++){
      if(!plotThisBlinding[blindingList][bl]) continue;
      for(int v2=0; v2<nVarPairs; v2++){
	if(!plotThisVarPair[varPairList][v2]) continue;
	if(varPairName[v2].find("M4L118130")!=string::npos && bl!=unblinded && bl!=fullyblind) continue;
	canvasName = string(Form("c_%s_2D_%s_%s",sBlinding[bl].c_str(),varPairName[v2].c_str(),sFinalState[FINALSTATE].c_str()));
	c2[v2] = new TCanvas(canvasName.c_str(),canvasName.c_str(),500,500);
	DrawDataMC2D(c2[v2],h2[v2][bl],g2Data[v2][bl],v2,bl,lumiText,STYLE2DPLOT,varPairLogx[v2],varPairLogy[v2]);
	SaveCanvas(outputDirectory,c2[v2]);
      }
    }
  }
  //*/
  /*
  TCanvas* c2[nVarPairs][NSTYLES2DPLOT];
  if(DO2DPLOTS){
    for(int st=0; st<NSTYLES2DPLOT; st++){
    for(int bl=0; bl<nBlindings; bl++){
      if(!plotThisBlinding[blindingList][bl]) continue;
      for(int v2=0; v2<nVarPairs; v2++){
	if(!plotThisVarPair[varPairList][v2]) continue;
	if(varPairName[v2].find("M4L118130")!=string::npos && bl!=unblinded) continue;
	canvasName = string(Form("c_%s_2D_%s_%s_st%i",sBlinding[bl].c_str(),varPairName[v2].c_str(),sFinalState[FINALSTATE].c_str(),st));
	c2[v2][st] = new TCanvas(canvasName.c_str(),canvasName.c_str(),500,500);
	DrawDataMC2D(c2[v2][st],h2[v2][bl],g2Data[v2][bl],v2,bl,lumiText,st,varPairLogx[v2],varPairLogy[v2]);
	SaveCanvas(outputDirectory,c2[v2][st]);
      }
    }
    }
  }
  //*/

}





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void plotDataVsMC(bool redoHistograms = true, string outputPath = "PlotsDataVsMC/") {

  // --------------- inputs ---------------

  // Define input/output location
  string inputPathMC   = "";
  string inputPathData = "";
  string inputFileDataForCR = "../DataTrees_160624/ZZ4lAnalysis.root";
  string inputFileFakeRates = "../../data/FakeRates/FakeRate_SS_2015.root"; //FIXME: to be replaced
  //string outputPath = "$pl/";

  // Define the luminosity
  float lumi = 2.6;
  string lumiText = "2.6 fb^{-1}";


  // Choose a list of 1D plots
  int variableList = 3;//4;//5;//6;//

  // Choose a list of 2D plots
  int varPairList = 3;//4;//5;//6;//

  // Choose a list of ways of blinding some m4l regions
  int blindingList = 5;//4;//1;//


  // --------------- processing ---------------

  gSystem->Exec(("mkdir -p "+outputPath).c_str());

  // Prepare the histograms and store them in a ROOT file (to be done only once, it can take a few minutes)
  if(redoHistograms)
    doHistograms(inputPathMC, inputPathData, lumi);

  // Prepare Z+X histograms
  if(USEZPLUSXFULLRUN2SS)
    doHistogramsZPlusXSS(inputFileDataForCR, inputFileFakeRates, lumi);

  // Do the plots (+- instantaneous)
  doPlots(outputPath, variableList, varPairList, blindingList, lumi, lumiText);

}

