/// 
/// usage: 
/// -specify parameters (input/output location, luminosity) at the end of this file
/// -run with:
///   root -l -b -q plotDataVsMC.C++
/// -later, once histograms have been stored, use parameter 0 to not rerun over all MC:
///   root -l -b -q "plotDataVsMC.C++(0)"
///

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

#include "Math/DistFunc.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TFrame.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
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
#include "ZZAnalysis/AnalysisStep/src/cConstants.cc"
#include "ZZAnalysis/AnalysisStep/src/Discriminants.cc"
#include "ZZAnalysis/AnalysisStep/src/Category.cc"
#include "ZZAnalysis/AnalysisStep/src/bitops.cc"
#include "ZZAnalysis/AnalysisStep/interface/FinalStates.h"

using namespace std;

#define DEBUG 0
#define VERBOSE2 0

#define DO1DPLOTS 1
#define DO2DPLOTS 1
#define PLOTLEVEL 2 // choose a set of plots: 0:test 1:PAS 2:AN 3:All
#define BLINDING 5 // 0:blind 5:unblinded, see below for intermediary options

#define APPLYKFACTORS 1
#define USEDYANDTTBAR 0
#define REBINDYANDTTBAR 0
#define MASKH125FORLOWMASS 1
#define MASKH125FORHIGHMASS 1
#define MASKDATAFORHIGHMASS 0
#define USEBBH 0
#define MERGE2E2MU 1 // If activated, 2e2mu means "2e2mu and 2mu2e". Don't change this unless you really know what you are doing.

#define USEZPLUSXANALYTICALSHAPE 1 // Obtain Z+X histograms from the analytical shapes when they are available (up to now, only for m4l).
#define BUILDZPLUSXHISTOFROMSSCR 1 // Obtain Z+X histograms by running over the SS CR. Overriden by the analytical shapes when relevant.
#define SMOOTHZPLUSXHISTOFROMSSCR 1 // Smooth the aforementioned SS-CR-based histos, preserving their normalization.
#define RENORMZPLUSXTOOFFICIALCATEG 0 // Normalize all Z+X histos to the hardcoded SS-based numbers in categories.
#define RENORMZPLUSXTOOFFICIALINCL 1 // (Re-)normalize all Z+X histos to the hardcoded SS/OS-combined inclusive numbers.

//Combined 0S+SS inclusive Z+X prediction
//Here: numbers for the full 2016 data set, computed by Roberto on March 3rd 2017
Float_t normZPlusXFullSR4e    = 21.1;
Float_t normZPlusXFullSR4mu   = 34.4;
Float_t normZPlusXFullSR2e2mu = 59.9;

//Normalization of Z+X in final states and categories.
//These numbers get overwritten if the RENORMZPLUSXTOOFFICIALCATEG flag is off.
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

#define STYLE1DPLOT 2 // 0:Legacy-like 1:Jamboree2015 2:Moriond2016
#define DRAWLINES (STYLE1DPLOT!=1)
#define LINEWIDTH (STYLE1DPLOT==0?2:STYLE1DPLOT==2?1:1)
#define CUSTOMXAXISLABELS 1
#define DRAWDATAMCRATIO 0
#define DRAWWP1D 1

#define NSTYLES2DPLOT 8
#define STYLE2DPLOT 1 // 0:rainbow 1:gray 2:pink 3:orange 4:yellow 5:blue 6:teal 7:col+box 8:blue-yellow
#define MARKERSTYLE 0 // 0:full(legacy) 1:open
#define LEGENDOUTOF2DFRAME 0
#define CATEGIN2D 1
#define DRAWWP2D 1

// These are the same 4 WP as in Category.cc
#define WP2J 0.437 // This is the value at 125GeV (0.431 at 118 GeV) of 1.043-460./(ZZMass+634.). The latter formula is also hardcoded in the definition of varPairExprWP.
#define WP1J 0.697
#define WPWH 0.951
#define WPZH 0.9937
// To change the c-constants for visualization purposes:
#define NEWWP2J 0.5
#define NEWWP1J 0.5
#define COMWPVH 0.5 //0.8
#define NEWWPWH COMWPVH
#define NEWWPZH COMWPVH
#define NEWWPVH COMWPVH // This one only makes sense if DWH and DZH have the same WP.
#define MINFACTVH 200. //80.


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

enum Process {Data=0, H125=1, H125VBF=2, H125VH=3, H125TTH=4, H125NONVBFVHTTH=5, qqZZ=6, ggZZ=7, DY=8, ttbar=9};
const int nProcesses = 10;
string sProcess[nProcesses] = {"Data", "H125", "H125VBF", "H125VH", "H125TTH", "H125NONVBFVHTTH", "qqZZ", "ggZZ", "DY", "ttbar"};
string processLabel[nProcesses] = {" Data", " H(125)", " H(125), VBF", " H(125), VH", " H(125), t#bar{t}H", " H(125), other", " q#bar{q}#rightarrowZZ, Z#gamma*", " gg#rightarrowZZ, Z#gamma*", " Z + jets", " t#bar{t}"};
Int_t processFillColor[nProcesses] = {
  TColor::GetColor("#000000"), 
  TColor::GetColor((STYLE1DPLOT==0)?"#ffffff":(STYLE1DPLOT==1)?"#ff9090":"#ffb2b2"), //"#ffc18b"),
  TColor::GetColor((STYLE1DPLOT==0)?"#ffffff":(STYLE1DPLOT==1)?"#ff6868":"#ff9b9b"),
  TColor::GetColor((STYLE1DPLOT==0)?"#ffffff":(STYLE1DPLOT==1)?"#ff6868":"#ff9b9b"),
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
  TColor::GetColor((STYLE1DPLOT==0)?"#ff0000":"#cc0000"),//"#770000"), 
  TColor::GetColor((STYLE1DPLOT==0)?"#ff0000":"#cc0000"),//"#770000"), 
  TColor::GetColor("#000099"),
  TColor::GetColor("#000099"),
  TColor::GetColor("#003300"),
  TColor::GetColor("#5f3f3f")
};
Bool_t useProcess[nProcesses] = {1,1,0,0,0,0,1,1,USEDYANDTTBAR,USEDYANDTTBAR,};


enum FinalState {fs4mu=0, fs4e=1, fs2e2mu=2, fs2mu2e=3};
const int nFinalStates = 4;
string sFinalState[nFinalStates+1] = {"4mu", "4e", "2e2mu", "2mu2e", "4l"};
string fsLabel[nFinalStates+1] = {"4#mu", "4e", "2e2#mu", "2#mu2e", "4#font[12]{l}"};
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


// Moriond 2017 categorization
const int nCat = 7;
string sCategory[nCat+1] = {
  "UnTagged",
  "VBF1jTagged",
  "VBF2jTagged",
  "VHLeptTagged",
  "VHHadrTagged",
  "ttHTagged",
  "VHMETTagged",
  "inclusive",
};
string categoryLabel[nCat+1] = {
  "untagged category",
  "VBF-1jet-tagged category",
  "VBF-2jet-tagged category", 
  "VH-leptonic-tagged category",
  "VH-hadronic-tagged category",
  "t#bar{t}H-tagged category",
  "VH-MET-tagged category",
  "all event categories",
};
string categoryLegLabel[nCat+1] = {
  "untagged",
  "VBF-1j tagged",
  "VBF-2j tagged", 
  "VH-lept. tagged",
  "VH-hadr. tagged",
  "t#bar{t}H tagged",
  "VH-MET tagged",
  "inclusive",
};
Int_t categMarkerStyle[nCat] = {20,26,32,28,27,30,25};
//Int_t categMarkerStyle[nCat] = {24,22,23,34,33,29,21};


enum ResonantStatus {resonant=0, nonresonant=1};
const int nRS = 2;
string sResonantStatus[nRS+1] = {"resonant", "nonresonant", "allres"};


enum Blindings {fullyblind=0, blindabove114=1, blindbelow130=2, blind114130=3, blind114130andabove300=4, unblinded=5};
const int nBlindings = 6;
string sBlinding[nBlindings] = {"fullyblind", "M4l70To114", "M4l130ToInf", "blind114130", "blind114130andabove300", "unblinded"};
string blindingLabel[nBlindings] = {"", "70 < m_{4#font[12]{l}} < 114 GeV", "m_{4#font[12]{l}} > 130 GeV", "#splitline{m_{4#font[12]{l}} > 70 GeV}{m_{4#font[12]{l}} #notin [114, 130] GeV}", "#splitline{m_{4#font[12]{l}} #in [70, 300] GeV}{m_{4#font[12]{l}} #notin [114, 130] GeV}", ""}; // corresponding cut have to be defined within the loops
Float_t xHistoBlindLow[nBlindings] = {  0.,  114.,   0., 114., 114.,  0. };
Float_t xHistoBlindUp [nBlindings] = { -1., 3000., 130., 130., 130., -1. };
Float_t xHistoBlind2Low[nBlindings] = {  0.,  0.,  0.,  0.,  300.,  0. };
Float_t xHistoBlind2Up [nBlindings] = { -1., -1., -1., -1., 3000., -1. };


struct Var1 {
  string Name;
  string XLabel;
  string YLabel;
  string CutLabel;
  Short_t plotLvl;
  Float_t Min;
  Float_t Max;
  Int_t Nbin;
  Bool_t isLogx;
  Bool_t isLogy;
  Bool_t Large;
  Int_t Categ;
  Bool_t InFS;
  Int_t restrXmax;
  Bool_t sepVbf;
  Bool_t sepVh;
  Bool_t sepTth;
  Float_t ValWP;
  Float_t MinFactor;
  Int_t CMSPos;
  Int_t LegPos;
  Int_t CutPos;
  Bool_t CXL;
  Int_t rebinDYTT;
  Float_t MaxCorrector[nBlindings];
  string prefix;
  Bool_t WithRatio;
  Bool_t SmoothZpX;
};

const int nVariables = 44;
Var1 myV1[nVariables] = {
//{ Name,                 XLabel,                           YLabel,           CutLabel,                                                        plotLvl, Min, Max, Nbin, isLogx, isLogy, Large, Categ, InFS, restrXmax, sepVbf, sepVh, sepTth, ValWP,   MinFactor, CMSPos, LegPos, CutPos, CXL,  rebinDYTT, MaxCorrector,              prefix, WithRatio, SmoothZpX},
  { "M4l_Full",           "m_{4#font[12]{l}} (GeV)",        "Events / 4 GeV", "",                                                              1,       70,  990, 230,  1,      0,      1,     nCat,  1,    0,         0,      0,     0,      0.,      0.,        0,      33,     0,      1,    1,         {1.1,1. ,1. ,1. ,1. ,1. }, "a",    0,         0,       },
  { "M4l_ForMergeRun1",   "m_{4#font[12]{l}} (GeV)",        "Events / 3 GeV", "",                                                              3,       1.5,1000.5,333, 0,      0,      0,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      33,     0,      0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "a",    0,         0,       },
//{ "M4l_Full_Bin5GeV",   "m_{4#font[12]{l}} (GeV)",        "Events / 5 GeV", "",                                                              3,       70,  1000,186,  1,      0,      1,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      33,     0,      0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "a",    0,         0,       },
  { "M4l_3p5TeV_Bin5GeV", "m_{4#font[12]{l}} (GeV)",        "Events / 5 GeV", "",                                                              3,       70,  3500,686,  1,      0,      1,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      33,     0,      0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "a",    0,         0,       },
  { "M4l_70110",          "m_{4#font[12]{l}} (GeV)",        "Events / 4 GeV", "",                                                              3,       70,  110, 10,   0,      0,      0,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      33,     0,      0,    1,         {1.3,1. ,1. ,1. ,1. ,1. }, "a",    0,         0,       },
  { "M4l_110150",         "m_{4#font[12]{l}} (GeV)",        "Events / 2 GeV", "",                                                              3,       110, 150, 20,   0,      0,      0,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      33,     0,      0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "a",    0,         0,       },
  { "M4l_70170",          "m_{4#font[12]{l}} (GeV)",        "Events / 2 GeV", "",                                                              1,       70,  170, 50,   0,      0,      0,     nCat,  1,    0,         0,      0,     0,      0.,      0.,        0,      33,     0,      0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "a",    0,         0,       },
  { "M4l_70170_0",        "m_{4#font[12]{l}} (GeV)",        "Events / 2 GeV", "",                                                              1,       70,  170, 50,   0,      0,      0,     0,     0,    0,         0,      0,     0,      0.,      0.,        0,      33,     0,      0,    1,         {1. ,1. ,1. ,1. ,1. ,1.1}, "g",    0,         0,       },
  { "M4l_70170_1",        "m_{4#font[12]{l}} (GeV)",        "Events / 4 GeV", "",                                                              1,       70,  170, 25,   0,      0,      0,     1,     0,    0,         1,      0,     0,      0.,      0.,        0,      33,     0,      0,    1,         {1. ,1. ,1. ,1. ,1. ,1.1}, "g",    0,         0,       },
  { "M4l_70170_2",        "m_{4#font[12]{l}} (GeV)",        "Events / 4 GeV", "",                                                              1,       70,  170, 25,   0,      0,      0,     2,     0,    0,         1,      0,     0,      0.,      0.,        0,      33,     0,      0,    1,         {1. ,1. ,1. ,1. ,1.,1.45}, "g",    0,         0,       },
  { "M4l_70170_3",        "m_{4#font[12]{l}} (GeV)",        "Events / 4 GeV", "",                                                              1,       70,  170, 25,   0,      0,      0,     3,     0,    0,         0,      1,     0,      0.,      0.,        0,      11,     0,      0,    1,         {1. ,1. ,1. ,1. ,1. ,1.5}, "g",    0,         0,       },
  { "M4l_70170_4",        "m_{4#font[12]{l}} (GeV)",        "Events / 4 GeV", "",                                                              1,       70,  170, 25,   0,      0,      0,     4,     0,    0,         0,      1,     0,      0.,      0.,        0,      33,     0,      0,    1,         {1. ,1. ,1. ,1. ,1. ,1.1}, "g",    0,         0,       },
  { "M4l_70170_5",        "m_{4#font[12]{l}} (GeV)",        "Events / 4 GeV", "",                                                              1,       70,  170, 25,   0,      0,      0,     5,     0,    0,         0,      0,     1,      0.,      0.,        0,      11,     0,      0,    1,         {1. ,1. ,1. ,1. ,1. ,1.1}, "g",    0,         0,       },
  { "M4l_70170_6",        "m_{4#font[12]{l}} (GeV)",        "Events / 4 GeV", "",                                                              1,       70,  170, 25,   0,      0,      0,     6,     0,    0,         0,      1,     0,      0.,      0.,        0,      11,     0,      0,    1,         {1. ,1. ,1. ,1. ,1. ,1.1}, "g",    0,         0,       },
  { "M4l_above170",       "m_{4#font[12]{l}} (GeV)",        "Events / 10 GeV","",                                                              1,       170, 1010,42,   0,      0,      0,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      33,     0,      0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "a",    0,         0,       },
  { "M4l_Full_HighKD",    "m_{4#font[12]{l}} (GeV)",        "Events / 4 GeV", "D_{bkg}^{kin} > 0.5            ",                               3,       70,  990, 230,  1,      0,      1,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      33,     0,      1,    1,         {1.1,1. ,1. ,1. ,1. ,1. }, "a",    0,         1,       },
  { "M4l_70170_HighKD",   "m_{4#font[12]{l}} (GeV)",        "Events / 2 GeV", "D_{bkg}^{kin} > 0.5         ",                                  2,       70,  170, 50,   0,      0,      0,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      33,     0,      0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "a",    0,         1,       },
  { "M4l_100170_HighKD",  "m_{4#font[12]{l}} (GeV)",        "Events / 2 GeV", "D_{bkg}^{kin} > 0.5",                                           3,       100, 170, 35,   0,      0,      0,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      33,     11,     0,    1,         {1. ,1. ,1. ,1.6,1.6,1. }, "a",    0,         1,       },
//{ "M4lrefit_Full",      "m_{4#font[12]{l}}^{refit} (GeV)","Events / 4 GeV", "",                                                              3,       70,  990, 230,  1,      0,      1,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      33,     0,      1,    1,         {1.1,1. ,1. ,1. ,1. ,1. }, "a",    0,         0,       },
  { "MZ1",                "m_{Z1} (GeV)",                   "Events / 2 GeV", "",                                                              2,       40,  120, 40,   0,      0,      0,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      11,     0,      0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "b",    0,         0,       },
  { "MZ1_M4L118130",      "m_{Z1} (GeV)",                   "Events / 4 GeV", "118 < m_{4#font[12]{l}} < 130 GeV",                             1,       40,  120, 20,   0,      0,      0,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      11,     0,      0,    1,         {1.1,1. ,1. ,1. ,1. ,1. }, "b",    0,         0,       },
  { "MZ2",                "m_{Z2} (GeV)",                   "Events / 2 GeV", "",                                                              2,       12,  120, 54,   0,      0,      0,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      11,     0,      0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "b",    0,         0,       },
  { "MZ2_M4L118130",      "m_{Z2} (GeV)",                   "Events / 4 GeV", "118 < m_{4#font[12]{l}} < 130 GeV",                             1,       12,  120, 27,   0,      0,      0,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      33,     0,      0,    1,         {1.9,1. ,1. ,1. ,1. ,1. }, "b",    0,         0,       },
  { "KD",                 "D_{bkg}^{kin}",                  "Events / 0.05",  "",                                                              2,       0,   1,   20,   0,      0,      0,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      33,     0,      0,    1,         {1.3,1.1,1.1,1.1,1.1,1.1}, "d",    0,         0,       },
//{ "KD_30bins",          "D_{bkg}^{kin}",                  "Events / bin",   "",                                                              3,       0,   1,   30,   0,      0,      0,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      33,     0,      0,    1,         {1.3,1.1,1.1,1.1,1.1,1.1}, "d",    0,         0,       },
  { "KD_30bins_M4Labove100","D_{bkg}^{kin}",                "Events / bin",   "m_{4#font[12]{l}} > 100 GeV",                                   3,       0,   1,   30,   0,      0,      0,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      33,     11,     0,    1,         {1.3,1.1,1.1,1.1,1.1,1.1}, "d",    0,         0,       },
  { "KD_M4L118130",       "D_{bkg}^{kin}",                  "Events / 0.1",   "118 < m_{4#font[12]{l}} < 130 GeV",                             1,       0,   1,   10,   0,      0,      0,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      11,     33,     0,    1,         {3. ,1. ,1. ,1. ,1. ,1.4}, "d",    0,         0,       },
  { "D2jet",              "D_{2jet}",                       "Events / 0.05",  "N(jets) #geq 2",                                                3,       0,   1,   20,   0,      0,      0,     nCat,  0,    0,         1,      0,     0,      0.,      0.,        0,      33,     11,     0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "e",    0,         0,       },
  { "D2jet_M4L118130",    "D_{2jet}",                       "Events / 0.1",   "#splitline{118 < m_{4#font[12]{l}} < 130 GeV}{N(jets) #geq 2}", 1,       0,   1,   10,   0,      0,      0,     nCat,  0,    0,         1,      0,     0,      NEWWP2J, 0.,        0,      33,     11,     0,    1,         {1. ,1. ,1. ,1. ,1. ,1.2}, "e",    0,         0,       },
  { "D1jet_M4L118130",    "D_{1jet}",                       "Events / 0.1",   "#splitline{118 < m_{4#font[12]{l}} < 130 GeV}{N(jets) = 1}",    1,       0,   1,   10,   0,      0,      0,     nCat,  0,    0,         1,      0,     0,      NEWWP1J, 0.,        0,      11,     33,     0,    1,         {1. ,1. ,1. ,1. ,1. ,1.3}, "e",    0,         0,       },
  { "DWH_M4L118130",      "D_{WH}",                         "Events / 0.1",   "#splitline{118 < m_{4#font[12]{l}} < 130 GeV}{N(jets) #geq 2}", 3,       0,   1,   10,   0,      0,      0,     nCat,  0,    0,         0,      1,     0,      NEWWPWH, 0.,        0,      33,     11,     0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "e",    0,         0,       },
  { "DZH_M4L118130",      "D_{ZH}",                         "Events / 0.1",   "#splitline{118 < m_{4#font[12]{l}} < 130 GeV}{N(jets) #geq 2}", 3,       0,   1,   10,   0,      0,      0,     nCat,  0,    0,         0,      1,     0,      NEWWPZH, 0.,        0,      33,     11,     0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "e",    0,         0,       },
  { "DVH_M4L118130",      "D_{VH}",                         "Events / 0.1",   "#splitline{118 < m_{4#font[12]{l}} < 130 GeV}{N(jets) #geq 2}", 1,       0,   1,   10,   0,      0,      0,     nCat,  0,    0,         0,      1,     0,      NEWWPVH, 0.,        0,      33,     11,     0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "e",    0,         0,       },
//{ "D2jet_M4L118130",    "D_{2jet}",                       "Events / 0.05",  "#splitline{118 < m_{4#font[12]{l}} < 130 GeV}{N(jets) #geq 2}", 1,       0,   1,   20,   0,      1,      0,     nCat,  0,    0,         1,      0,     0,      NEWWP2J, 40.,       0,      33,     11,     0,    1,         {1. ,1. ,1. ,1. ,1. ,3. }, "e",    0,         0,       },
//{ "D1jet_M4L118130",    "D_{1jet}",                       "Events / 0.05",  "#splitline{118 < m_{4#font[12]{l}} < 130 GeV}{N(jets) = 1}",    1,       0,   1,   20,   0,      0,      0,     nCat,  0,    0,         1,      0,     0,      NEWWP1J, 0.,        0,      11,     33,     0,    1,         {1. ,1. ,1. ,1. ,1. ,1.8}, "e",    0,         0,       },
//{ "DWH_M4L118130",      "D_{WH}",                         "Events / 0.05",  "#splitline{118 < m_{4#font[12]{l}} < 130 GeV}{N(jets) #geq 2}", 1,       0,   1,   20,   0,      1,      0,     nCat,  0,    0,         0,      1,     0,      NEWWPWH, MINFACTVH, 0,      33,     11,     0,    3,         {1. ,1. ,1. ,1. ,1. ,1. }, "e",    0,         0,       },
//{ "DZH_M4L118130",      "D_{ZH}",                         "Events / 0.05",  "#splitline{118 < m_{4#font[12]{l}} < 130 GeV}{N(jets) #geq 2}", 1,       0,   1,   20,   0,      1,      0,     nCat,  0,    0,         0,      1,     0,      NEWWPZH, MINFACTVH, 0,      33,     11,     0,    3,         {1. ,1. ,1. ,1. ,1. ,1. }, "e",    0,         0,       },
//{ "DVH_M4L118130",      "D_{VH}"/*"max(D_{WH},D_{ZH})"*/, "Events / 0.05",  "#splitline{118 < m_{4#font[12]{l}} < 130 GeV}{N(jets) #geq 2}", 1,       0,   1,   20,   0,      1,      0,     nCat,  0,    0,         0,      1,     0,      NEWWPVH, MINFACTVH, 0,      33,     11,     0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "e",    0,         0,       },
  { "Pt4l",               "p_{T}^{4#font[12]{l}} (GeV)",    "Events / 10 GeV","",                                                              3,       0,   400, 40,   0,      0,      0,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      33,     0,      0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "h",    0,         0,       },
  { "Pt4l_M4L118130",     "p_{T}^{4#font[12]{l}} (GeV)",    "Events / 10 GeV","118 < m_{4#font[12]{l}} < 130 GeV",                             3,       0,   400, 40,   0,      0,      0,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      33,     11,     0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "h",    0,         0,       },
  { "Eta4l",              "#eta^{4#font[12]{l}}",           "Events / 0.5",   "",                                                              3,       -8,  8,   32,   0,      0,      0,     nCat,  1,    0,         0,      0,     0,      0.,      0.,        0,      33,     0,      0,    1,         {1.4,1. ,1. ,1.5,1.5,1.2}, "h",    1,         0,       },
  { "Eta4l_M4L118130",    "#eta^{4#font[12]{l}}",           "Events / 0.8",   "118 < m_{4#font[12]{l}} < 130 GeV",                             3,       -8,  8,   20,   0,      0,      0,     nCat,  1,    0,         0,      0,     0,      0.,      0.,        0,      33,     11,     0,    1,         {1. ,1. ,1. ,1. ,1. ,1.2}, "h",    0,         0,       },
  { "Eta4l_M4L70110",     "#eta^{4#font[12]{l}}",           "Events / 0.8",   "70 < m_{4#font[12]{l}} < 110 GeV",                              3,       -8,  8,   20,   0,      0,      0,     nCat,  1,    0,         0,      0,     0,      0.,      0.,        0,      33,     11,     0,    1,         {1. ,1. ,1. ,1. ,1. ,1.2}, "h",    0,         0,       },
  { "Eta4l_M4Labove180",  "#eta^{4#font[12]{l}}",           "Events / 0.5",   "m_{4#font[12]{l}} > 180 GeV",                                   3,       -8,  8,   32,   0,      0,      0,     nCat,  1,    0,         0,      0,     0,      0.,      0.,        0,      33,     11,     0,    1,         {1. ,1. ,1. ,1. ,1. ,1.2}, "h",    0,         0,       },
  { "MET",                "E_{T}^{miss} (GeV)",             "Events / 10 GeV","",                                                              3,       0,   400, 40,   0,      0,      0,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      33,     0,      0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "h",    0,         0,       },
  { "MET_M4L118130",      "E_{T}^{miss} (GeV)",             "Events / 10 GeV","118 < m_{4#font[12]{l}} < 130 GeV",                             3,       0,   400, 40,   0,      0,      0,     nCat,  0,    0,         0,      0,     0,      0.,      0.,        0,      33,     11,     0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "h",    0,         0,       },
  { "NExtraLep",          "number of additional leptons",   "Events",         "",                                                              3,       0,   6,   6,    0,      1,      0,     nCat,  0,    2,         0,      0,     0,      0.,      20000.,    0,      33,     0,      0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "h",    0,         0,       },
  { "NExtraLep_M4L118130","number of additional leptons",   "Events",         "118 < m_{4#font[12]{l}} < 130 GeV",                             3,       0,   6,   6,    0,      1,      0,     nCat,  0,    2,         0,      0,     0,      0.,      20000.,    0,      33,     11,     0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "h",    0,         0,       },
  { "NJets",              "number of jets",                 "Events",         "",                                                              3,       0,   17,  17,   0,      1,      0,     nCat,  0,    4,         0,      0,     0,      0.,      200.,      0,      33,     0,      0,    1,         {1. ,1. ,10.,10.,10.,1. }, "h",    0,         0,       },
  { "NJets_M4L118130",    "number of jets",                 "Events",         "118 < m_{4#font[12]{l}} < 130 GeV",                             3,       0,   17,  17,   0,      1,      0,     nCat,  0,    4,         0,      0,     0,      0.,      200.,      0,      33,     11,     0,    1,         {1. ,1. ,10.,10.,10.,1. }, "h",    0,         0,       },
  { "NBtags",             "number of b-tagged jets",        "Events",         "",                                                              3,       0,   8,   8,    0,      1,      0,     nCat,  0,    2,         0,      0,     0,      0.,      5000.,     0,      33,     0,      0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "h",    0,         0,       },
  { "NBtags_M4L118130",   "number of b-tagged jets",        "Events",         "118 < m_{4#font[12]{l}} < 130 GeV",                             3,       0,   8,   8,    0,      1,      0,     nCat,  0,    2,         0,      0,     0,      0.,      5000.,     0,      33,     11,     0,    1,         {1. ,1. ,1. ,1. ,1. ,1. }, "h",    0,         0,       },
};


struct Var2 {
  string Name;
  string XLabel;
  string YLabel;
  string CutLabel;
  Short_t plotLvl;
  Float_t XMin;
  Float_t XMax;
  Int_t XNbin;
  Float_t YMin;
  Float_t YMax;
  Int_t YNbin;
  Bool_t isLogx;
  Bool_t isLogy;
  Bool_t WithCateg;
  string exprWP;
  Int_t LegPos;
  Bool_t LegIsWhite;
  string prefix;
};

const int nVarPairs = 9;
Var2 myV2[nVarPairs] = {
//{ Name,                   XLabel,                    YLabel,          CutLabel,                            plotLvl, XMin, XMax, XNbin, YMin, YMax, YNbin, isLogx, isLogy, WithCateg, exprWP,                       LegPos, LegIsWhite, prefix },
//{ "M4lVsKD",              "m_{4#font[12]{l}} (GeV)", "D_{bkg}^{kin}", "",                                  3,       100,  1000, 180,   0,    1,    30,    1,      0,      0,         "",                           33,     1,          "d"    },
//{ "M4lVsKD_M4L70110",     "m_{4#font[12]{l}} (GeV)", "D_{bkg}^{kin}", "",                                  3,       70,   110,  40,    0,    1,    30,    0,      0,      0,         "",                           33,     1,          "d"    },
  { "M4lVsKD_M4L100170",    "m_{4#font[12]{l}} (GeV)", "D_{bkg}^{kin}", "",                                  1,       100,  170,  35,    0,    1,    30,    0,      0,      1,         "",                           33,     1,          "d"    },
  { "M4lVsKD_M4L1701000",   "m_{4#font[12]{l}} (GeV)", "D_{bkg}^{kin}", "",                                  2,       170,  1000, 166,   0,    1,    30,    0,      0,      0,         "",                           33,     1,          "d"    },
//{ "M4lVsD2jet_M4L100170", "m_{4#font[12]{l}} (GeV)", "D_{2jet}",      "",                                  1,       100,  170,  35,    0,    1,    30,    0,      0,      1,         "1.043-460./(x+634.)",        33,     1,          "e"    },
  { "M4lVsD2jet_M4L100170", "m_{4#font[12]{l}} (GeV)", "D_{2jet}",      "",                                  2,       100,  170,  35,    0,    1,    30,    0,      0,      1,         string(Form("%.3f",NEWWP2J)), 33,     1,          "e"    },
  { "M4lVsD1jet_M4L100170", "m_{4#font[12]{l}} (GeV)", "D_{1jet}",      "",                                  2,       100,  170,  35,    0,    1,    30,    0,      0,      1,         string(Form("%.3f",NEWWP1J)), 33,     1,          "e"    },
  { "M4lVsDWH_M4L100170",   "m_{4#font[12]{l}} (GeV)", "D_{WH}",        "",                                  3,       100,  170,  35,    0,    1,    30,    0,      0,      1,         string(Form("%.3f",NEWWPWH)), 33,     1,          "e"    },
  { "M4lVsDZH_M4L100170",   "m_{4#font[12]{l}} (GeV)", "D_{ZH}",        "",                                  3,       100,  170,  35,    0,    1,    30,    0,      0,      1,         string(Form("%.3f",NEWWPZH)), 33,     1,          "e"    },
  { "M4lVsDVH_M4L100170",   "m_{4#font[12]{l}} (GeV)", "D_{VH}",        "",                                  2,       100,  170,  35,    0,    1,    30,    0,      0,      1,         string(Form("%.3f",NEWWPVH)), 33,     1,          "e"    },
  { "MZ1VsMZ2",             "m_{Z1} (GeV)",            "m_{Z2} (GeV)",  "",                                  2,       40,   120,  80,    12,   120,  108,   0,      0,      0,         "",                           11,     0,          "b"    },
  { "MZ1VsMZ2_M4L118130",   "m_{Z1} (GeV)",            "m_{Z2} (GeV)",  "118 < m_{4#font[12]{l}} < 130 GeV", 1,       40,   120,  80,    12,   120,  108,   0,      0,      0,         "",                           11,     0,          "b"    },
//{ "MZ1VsMZ2_alt1",        "m_{Z1} (GeV)",            "m_{Z2} (GeV)",  "",                                  3,       40,   120,  40,    12,   120,  54,    0,      0,      0,         "",                           11,     0,          "b"    },
//{ "MZ1VsMZ2_alt2",        "m_{Z1} (GeV)",            "m_{Z2} (GeV)",  "",                                  3,       75,   105,  60,    75,   105,  60,    0,      0,      0,         "",                           11,     0,          "b"    },
//{ "M4lVsM4lRefit",        "m_{4#font[12]{l}} (GeV)", "m_{4#font[12]{l}}^{refit} (GeV)", "",                3,       70,   1000, 310,   70,   1000, 310,   0,      0,      0,         "",                           33,     0,          "z"    },
//{ "M4lVsM4lRefit_100180", "m_{4#font[12]{l}} (GeV)", "m_{4#font[12]{l}}^{refit} (GeV)", "",                3,       100,  180,  80,    100,  180,  80,    0,      0,      0,         "",                           11,     0,          "z"    },
};




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void doHistograms(string inputPathMC, string inputPathData, double lumi)
{

  //TFile* ggZZKFactorFile = TFile::Open("../../data/kfactors/Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
  //TSpline3* sp = (TSpline3*)ggZZKFactorFile->Get("sp_kfactor_Nominal");

  const int nDatasets = 17;
  string datasets[nDatasets] = {
    "AllData",
    "ggH125",
    "VBFH125",
    "WplusH125",
    "WminusH125",
    "ZH125",
    "ttH125",
    "bbH125",
    "ZZTo4l",//"ZZTo4lamcatnlo",//
    "ggTo4e_Contin_MCFM701",//"ggZZ4e",
    "ggTo4mu_Contin_MCFM701",//"ggZZ4mu",
    "ggTo4tau_Contin_MCFM701",//"ggZZ4tau",
    "ggTo2e2mu_Contin_MCFM701",//"ggZZ2e2mu",
    "ggTo2e2tau_Contin_MCFM701",//"ggZZ2e2tau",
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
  Float_t PFMET;
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
  Float_t p_GG_SIG_ghg2_1_ghz1_1_JHUGen;
  Float_t p_QQB_BKG_MCFM;
  Float_t p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal;
  Float_t p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal;
  Float_t p_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
  Float_t pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
  Float_t p_JQCD_SIG_ghg2_1_JHUGen_JECNominal;
  Float_t p_HadWH_SIG_ghw1_1_JHUGen_JECNominal;
  Float_t p_HadZH_SIG_ghz1_1_JHUGen_JECNominal;
  Float_t p_HadWH_mavjj_JECNominal;
  Float_t p_HadZH_mavjj_JECNominal;
  Float_t Z1Mass;
  Float_t Z2Mass;
  Short_t Z1Flav;
  Short_t Z2Flav;
  Float_t DiJetFisher;
  vector<Float_t> *LepPt = 0;
  vector<Float_t> *LepEta = 0;
  vector<Float_t> *LepPhi = 0;
  vector<Float_t> *LepLepId = 0;
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

  TH1F* h1[nVariables][nBlindings][nFinalStates+1][nCat+1][nRS+1][nProcesses];
  TH2F* h2[nVarPairs ][nBlindings][nFinalStates+1][nCat+1][nRS+1][nProcesses];
  for(int bl=0; bl<nBlindings; bl++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      for(int cat=0; cat<nCat+1; cat++){
	for(int rs=0; rs<nRS+1; rs++){
	  for(int pr=0; pr<nProcesses; pr++){
	    for(int v1=0; v1<nVariables; v1++){
	      h1[v1][bl][fs][cat][rs][pr] = new TH1F(
                 Form("h1_%s_%s_%s_%s_%s_%s",myV1[v1].Name.c_str(),sBlinding[bl].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str(),sResonantStatus[rs].c_str(),sProcess[pr].c_str()),
		 Form(";%s;%s",myV1[v1].XLabel.c_str(),myV1[v1].YLabel.c_str()),
		 myV1[v1].Nbin,myV1[v1].Min,myV1[v1].Max);
	    }
	    for(int v2=0; v2<nVarPairs; v2++){
	      h2[v2][bl][fs][cat][rs][pr] = new TH2F(
                 Form("h2_%s_%s_%s_%s_%s_%s",myV2[v2].Name.c_str(),sBlinding[bl].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str(),sResonantStatus[rs].c_str(),sProcess[pr].c_str()),
		 Form(";%s;%s",myV2[v2].XLabel.c_str(),myV2[v2].YLabel.c_str()),
		 myV2[v2].XNbin,myV2[v2].XMin,myV2[v2].XMax,
		 myV2[v2].YNbin*2,myV2[v2].YMin,myV2[v2].YMax*2);
	    }
	  }
	}
      }
    }
  }

  vector<Float_t> g2DataX[nVarPairs][nBlindings][nFinalStates+1][nCat+1][nRS+1];
  vector<Float_t> g2DataY[nVarPairs][nBlindings][nFinalStates+1][nCat+1][nRS+1];
  vector<Float_t> g2DataEX[nVarPairs][nBlindings][nFinalStates+1][nCat+1][nRS+1];
  vector<Float_t> g2DataEY[nVarPairs][nBlindings][nFinalStates+1][nCat+1][nRS+1];
  TGraphErrors* g2Data[nVarPairs][nBlindings][nFinalStates+1][nCat+1][nRS+1];
  
  int currentProcess;
  int currentFinalState;
  int currentCategory;
  int currentResStatus;


  //---------- Will loop over all datasets

  for(int d=0; d<nDatasets-2; d++){

    if(datasets[d]=="bbH125" && !USEBBH) continue;

    //----- assign dataset to correct process
    currentProcess = -1;
    if(datasets[d]=="AllData") currentProcess = Data;
    if(datasets[d]=="ggH125") currentProcess = H125NONVBFVHTTH;
    if(datasets[d]=="VBFH125") currentProcess = H125VBF;
    if(datasets[d]=="WplusH125") currentProcess = H125VH;
    if(datasets[d]=="WminusH125") currentProcess = H125VH;
    if(datasets[d]=="ZH125") currentProcess = H125VH;
    if(datasets[d]=="ttH125") currentProcess = H125TTH;
    if(datasets[d]=="bbH125") currentProcess = H125NONVBFVHTTH;
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

    string inputFileName = string(Form("%s%s/ZZ4lAnalysis.root",(currentProcess==Data?inputPathData:inputPathMC).c_str(),datasets[d].c_str()));
    inputFile[d] = TFile::Open(inputFileName.c_str());

    hCounters[d] = (TH1F*)inputFile[d]->Get("ZZTree/Counters");
    NGenEvt[d] = (Long64_t)hCounters[d]->GetBinContent(1);
    gen_sumWeights[d] = (Long64_t)hCounters[d]->GetBinContent(40);
    partialSampleWeight[d] = lumi * 1000 / gen_sumWeights[d] ;

    inputTree[d] = (TTree*)inputFile[d]->Get("ZZTree/candTree");
    inputTree[d]->SetBranchAddress("RunNumber", &nRun);
    inputTree[d]->SetBranchAddress("EventNumber", &nEvent);
    inputTree[d]->SetBranchAddress("LumiNumber", &nLumi);
    inputTree[d]->SetBranchAddress("PFMET", &PFMET);
    inputTree[d]->SetBranchAddress("Nvtx", &Nvtx);
    inputTree[d]->SetBranchAddress("NObsInt", &NObsInt);
    inputTree[d]->SetBranchAddress("NTrueInt", &NTrueInt);
    if(currentProcess==ggZZ){
      inputTree[d]->SetBranchAddress("KFactor_QCD_ggZZ_Nominal", &KFactor_QCD_ggZZ_Nominal);
    }
    if(currentProcess==qqZZ){
      inputTree[d]->SetBranchAddress("KFactor_EW_qqZZ", &KFactor_EW_qqZZ);
      inputTree[d]->SetBranchAddress("KFactor_QCD_qqZZ_dPhi", &KFactor_QCD_qqZZ_dPhi);
      inputTree[d]->SetBranchAddress("KFactor_QCD_qqZZ_M", &KFactor_QCD_qqZZ_M);
      inputTree[d]->SetBranchAddress("KFactor_QCD_qqZZ_Pt", &KFactor_QCD_qqZZ_Pt);
    }
    inputTree[d]->SetBranchAddress("ZZsel", &ZZsel);
    inputTree[d]->SetBranchAddress("ZZMass", &ZZMass);
    inputTree[d]->SetBranchAddress("ZZMassErr", &ZZMassErr);
    inputTree[d]->SetBranchAddress("ZZMassErrCorr", &ZZMassErrCorr);
    inputTree[d]->SetBranchAddress("ZZMassRefit", &ZZMassRefit);
    inputTree[d]->SetBranchAddress("ZZMassRefitErr", &ZZMassRefitErr);
    inputTree[d]->SetBranchAddress("ZZMassUnrefitErr", &ZZMassUnrefitErr);
    inputTree[d]->SetBranchAddress("ZZPt", &ZZPt);
    inputTree[d]->SetBranchAddress("ZZEta", &ZZEta);
    inputTree[d]->SetBranchAddress("p_GG_SIG_ghg2_1_ghz1_1_JHUGen", &p_GG_SIG_ghg2_1_ghz1_1_JHUGen);
    inputTree[d]->SetBranchAddress("p_QQB_BKG_MCFM", &p_QQB_BKG_MCFM);
    inputTree[d]->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal", &p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal);
    inputTree[d]->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal", &p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal);
    inputTree[d]->SetBranchAddress("p_JVBF_SIG_ghv1_1_JHUGen_JECNominal", &p_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
    inputTree[d]->SetBranchAddress("pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal", &pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
    inputTree[d]->SetBranchAddress("p_JQCD_SIG_ghg2_1_JHUGen_JECNominal", &p_JQCD_SIG_ghg2_1_JHUGen_JECNominal);
    inputTree[d]->SetBranchAddress("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal", &p_HadWH_SIG_ghw1_1_JHUGen_JECNominal);
    inputTree[d]->SetBranchAddress("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal", &p_HadZH_SIG_ghz1_1_JHUGen_JECNominal);
    inputTree[d]->SetBranchAddress("p_HadWH_mavjj_JECNominal", &p_HadWH_mavjj_JECNominal);
    inputTree[d]->SetBranchAddress("p_HadZH_mavjj_JECNominal",&p_HadZH_mavjj_JECNominal);
    inputTree[d]->SetBranchAddress("Z1Mass", &Z1Mass);
    inputTree[d]->SetBranchAddress("Z2Mass", &Z2Mass);
    inputTree[d]->SetBranchAddress("Z1Flav", &Z1Flav);
    inputTree[d]->SetBranchAddress("Z2Flav", &Z2Flav);
    inputTree[d]->SetBranchAddress("LepPt", &LepPt);
    inputTree[d]->SetBranchAddress("LepEta", &LepEta);
    inputTree[d]->SetBranchAddress("LepPhi", &LepPhi);
    inputTree[d]->SetBranchAddress("LepLepId", &LepLepId);
    inputTree[d]->SetBranchAddress("nExtraLep", &nExtraLep);
    inputTree[d]->SetBranchAddress("ExtraLepEta", &ExtraLepEta);
    inputTree[d]->SetBranchAddress("ExtraLepPhi", &ExtraLepPhi);
    inputTree[d]->SetBranchAddress("nExtraZ", &nExtraZ);
    inputTree[d]->SetBranchAddress("nCleanedJetsPt30", &nJets);
    inputTree[d]->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSF", &nJetsBTagged);
    inputTree[d]->SetBranchAddress("JetPt", &JetPt);
    inputTree[d]->SetBranchAddress("JetEta", &JetEta);
    inputTree[d]->SetBranchAddress("JetPhi", &JetPhi);
    inputTree[d]->SetBranchAddress("JetMass", &JetMass);
    inputTree[d]->SetBranchAddress("JetQGLikelihood", &JetQGLikelihood);
    inputTree[d]->SetBranchAddress("DiJetFisher", &DiJetFisher);
    if(currentProcess!=Data){
      inputTree[d]->SetBranchAddress("overallEventWeight", &overallEventWeight);
      inputTree[d]->SetBranchAddress("xsec", &xsec);     
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
    }


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
	if(currentProcess==qqZZ)      kfactor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M;
	else if(currentProcess==ggZZ) kfactor = KFactor_QCD_ggZZ_Nominal;
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
      currentCategory = categoryMor17(
	 nExtraLep,
	 nExtraZ,
	 nJets,
	 nJetsBTagged,
	 jetQGL,
	 p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
	 p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
	 p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
	 p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
	 pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
	 p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
	 p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
    p_HadWH_mavjj_JECNominal,
    p_HadZH_mavjj_JECNominal,
	 jetPhi,
	 ZZMass,
	 PFMET,
	 true,
	 false
	 );


      /* //----- here, define resonant signal as H->4l where l=e,mu (excluding decays to taus and 'wrong signal' from associated production)
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
      //*/

      //----- fill histograms


      //FIXME: Switch to functions in Discriminants.h
      Float_t KD = p_GG_SIG_ghg2_1_ghz1_1_JHUGen / ( p_GG_SIG_ghg2_1_ghz1_1_JHUGen + p_QQB_BKG_MCFM*getDbkgkinConstant(Z1Flav*Z2Flav,ZZMass) );
      Float_t D2jet = (nJets>=2) ? DVBF2j_ME(p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass) : -2;
      Float_t D1jet = (nJets==1) ? DVBF1j_ME(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass) : -2;
      Float_t DWH = (nJets>=2) ? DWHh_ME(p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_HadWH_mavjj_JECNominal, ZZMass) : -2;
      Float_t DZH = (nJets>=2) ? DZHh_ME(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_HadZH_mavjj_JECNominal, ZZMass) : -2;
      Float_t oldCConstD2jet = getDVBF2jetsConstant(ZZMass);
      Float_t oldCConstD1jet = getDVBF1jetConstant(ZZMass);
      Float_t oldCConstDWH = getDWHhConstant(ZZMass);
      Float_t oldCConstDZH = getDZHhConstant(ZZMass);
      Float_t newCConstD2jet = getDVBF2jetsConstant_shiftWP(ZZMass,false,NEWWP2J);
      Float_t newCConstD1jet = getDVBF1jetConstant_shiftWP(ZZMass,false,NEWWP1J);
      Float_t newCConstDWH = getDWHhConstant_shiftWP(ZZMass,false,NEWWPWH);
      Float_t newCConstDZH = getDZHhConstant_shiftWP(ZZMass,false,NEWWPZH);
      //FIXME: What is happening here??? seems some leftover from a temporary hack.
      D2jet = 1/(newCConstD2jet/oldCConstD2jet*(1/D2jet-1)+1);
      D1jet = 1/(newCConstD1jet/oldCConstD1jet*(1/D1jet-1)+1);
      DWH = 1/(newCConstDWH/oldCConstDWH*(1/DWH-1)+1);
      DZH = 1/(newCConstDZH/oldCConstDZH*(1/DZH-1)+1);
      Float_t DVH = fmax(DWH,DZH);

      Float_t varVal[nVariables];
      Bool_t varPassCut[nVariables];
      for(int v1=0; v1<nVariables; v1++){
	string varStr = myV1[v1].Name.substr(0,myV1[v1].Name.find("_"));
	if     (varStr=="M4l") varVal[v1] = ZZMass;
	else if(varStr=="MZ1") varVal[v1] = Z1Mass;
	else if(varStr=="MZ2") varVal[v1] = Z2Mass;
	else if(varStr=="KD") varVal[v1] = KD;
	else if(varStr=="D2jet") varVal[v1] = D2jet;
	else if(varStr=="D1jet") varVal[v1] = D1jet;
	else if(varStr=="DWH") varVal[v1] = DWH;
	else if(varStr=="DZH") varVal[v1] = DZH;
	else if(varStr=="DVH") varVal[v1] = DVH;
	else if(varStr=="MET") varVal[v1] = PFMET;
	else if(varStr=="NExtraLep") varVal[v1] = (Float_t)nExtraLep;
	else if(varStr=="NJets") varVal[v1] = (Float_t)nJets;
	else if(varStr=="NBtags") varVal[v1] = (Float_t)nJetsBTagged;
	else if(varStr=="Pt4l") varVal[v1] = ZZPt;
	else if(varStr=="Eta4l") varVal[v1] = ZZEta;
	else if(varStr=="M4lrefit") varVal[v1] = ZZMassRefit;
	varPassCut[v1] = true;
	if(myV1[v1].Name.find("M4L118130")!=string::npos) varPassCut[v1] = varPassCut[v1] && 118<=ZZMass && ZZMass<=130; 
	if(myV1[v1].Name.find("M4L70110")!=string::npos) varPassCut[v1] = varPassCut[v1] && 70<=ZZMass && ZZMass<=110; 
	if(myV1[v1].Name.find("M4Labove100")!=string::npos) varPassCut[v1] = varPassCut[v1] && 100<=ZZMass; 
	if(myV1[v1].Name.find("M4Labove180")!=string::npos) varPassCut[v1] = varPassCut[v1] && 180<=ZZMass; 
	if(myV1[v1].Name.find("HighKD")!=string::npos) varPassCut[v1] = varPassCut[v1] && KD>0.5;
	if(varStr=="D2jet" || varStr=="DWH" || varStr=="DZH") varPassCut[v1] = varPassCut[v1] && nJets>=2; 
	else if(varStr=="D1jet") varPassCut[v1] = varPassCut[v1] && nJets==1; 
      }

      Float_t varPairVal[nVarPairs][2];
      Bool_t varPairPassCut[nVarPairs];
      for(int v2=0; v2<nVarPairs; v2++){
	string varPairStr = myV2[v2].Name.substr(0,myV2[v2].Name.find("_"));
	if     (varPairStr=="M4lVsKD") { varPairVal[v2][0] = ZZMass; varPairVal[v2][1] = KD; }
	else if(varPairStr=="M4lVsD2jet") { varPairVal[v2][0] = ZZMass; varPairVal[v2][1] = D2jet; }
	else if(varPairStr=="M4lVsD1jet") { varPairVal[v2][0] = ZZMass; varPairVal[v2][1] = D1jet; }
	else if(varPairStr=="M4lVsDWH") { varPairVal[v2][0] = ZZMass; varPairVal[v2][1] = DWH; }
	else if(varPairStr=="M4lVsDZH") { varPairVal[v2][0] = ZZMass; varPairVal[v2][1] = DZH; }
	else if(varPairStr=="M4lVsDVH") { varPairVal[v2][0] = ZZMass; varPairVal[v2][1] = DVH; }
	else if(varPairStr=="MZ1VsMZ2") { varPairVal[v2][0] = Z1Mass; varPairVal[v2][1] = Z2Mass; }
	else if(varPairStr=="M4lVsM4lRefit") { varPairVal[v2][0] = ZZMass; varPairVal[v2][1] = ZZMassRefit; }
	varPairPassCut[v2] = true;
	if(myV2[v2].Name.find("M4L118130")!=string::npos) varPairPassCut[v2] = varPairPassCut[v2] && 118<=ZZMass && ZZMass<=130; 
	if(varPairStr=="M4lVsD2jet" || varPairStr=="M4lVsDWH" || varPairStr=="M4lVsDZH") varPairPassCut[v2] = varPairPassCut[v2] && nJets>=2; 
	else if(varPairStr=="M4lVsD1jet") varPairPassCut[v2] = varPairPassCut[v2] && nJets==1; 
      }

      bool fillM4l[nBlindings] = {
	currentProcess!=Data,
	(currentProcess!=Data || ZZMass < xHistoBlindLow[1] ),
	(currentProcess!=Data || ZZMass > xHistoBlindUp[2]  ),
	(currentProcess!=Data || ZZMass < xHistoBlindLow[3] || ZZMass  > xHistoBlindUp[3] ),
	(currentProcess!=Data || ZZMass < xHistoBlindLow[4] || (ZZMass > xHistoBlindUp[4]  && ZZMass < xHistoBlind2Low[4])),
	true,
      };
      bool fillOtherThanM4l[nBlindings] = {
	currentProcess!=Data,
	ZZMass < xHistoBlindLow[1],
	ZZMass > xHistoBlindUp[2] ,
	ZZMass < xHistoBlindLow[3] || ZZMass  > xHistoBlindUp[3],
	ZZMass < xHistoBlindLow[4] || (ZZMass > xHistoBlindUp[4] && ZZMass < xHistoBlind2Low[4]),
	true,
      };

      float m4lerr = ZZMassErr;
      //float m4lerr = ZZMassErrCorr;
      //float m4lerr = ZZMassUnrefitErr;
      //cout<<"m4lerr"<<m4lerr<<endl;

      for(int bl=0; bl<nBlindings; bl++){
	for(int v1=0; v1<nVariables; v1++){
	  if( (myV1[v1].Name.find("M4l")==0 && fillM4l         [bl]) ||
	      (myV1[v1].Name.find("M4l")!=0 && fillOtherThanM4l[bl])    ){
	    if(varPassCut[v1]){
	      h1[v1][bl][currentFinalState][currentCategory][nRS][currentProcess]->Fill(varVal[v1],(currentProcess==Data)?1.:eventWeight);
	    }
	  }
	}
	for(int v2=0; v2<nVarPairs; v2++){
	  if( (myV2[v2].Name.find("M4l")==0 && fillM4l         [bl]) ||
	      (myV2[v2].Name.find("M4l")!=0 && fillOtherThanM4l[bl])    ){
	    if(varPairPassCut[v2]){
	      h2[v2][bl][currentFinalState][currentCategory][nRS][currentProcess]->Fill(varPairVal[v2][0],varPairVal[v2][1],(currentProcess==Data)?1.:eventWeight);
	      if(currentProcess==Data){
		g2DataX[v2][bl][currentFinalState][currentCategory][nRS].push_back(varPairVal[v2][0]);
		g2DataY[v2][bl][currentFinalState][currentCategory][nRS].push_back(varPairVal[v2][1]);
		g2DataEX[v2][bl][currentFinalState][currentCategory][nRS].push_back((myV2[v2].Name.find("M4l")==0)?m4lerr:0.);
		g2DataEY[v2][bl][currentFinalState][currentCategory][nRS].push_back(0.);
	      }      
	    }
	  }
	}
      }

    } // end for entries

  } // end for datasets


  //---------- Fill 'inclusive' histograms
  cout<<"Adding histograms ..."<<endl;
  for(int bl=0; bl<nBlindings; bl++){
    if(bl!=BLINDING) continue;
    for(int fs=0; fs<nFinalStates+1; fs++){
      for(int cat=0; cat<nCat+1; cat++){
	//for(int rs=0; rs<nRS+1; rs++){
	for(int v1=0; v1<nVariables; v1++){
	  h1[v1][bl][fs][cat][nRS][H125]->Add(h1[v1][bl][fs][cat][nRS][H125VBF]);
	  h1[v1][bl][fs][cat][nRS][H125]->Add(h1[v1][bl][fs][cat][nRS][H125VH]);
	  h1[v1][bl][fs][cat][nRS][H125]->Add(h1[v1][bl][fs][cat][nRS][H125TTH]);
	  h1[v1][bl][fs][cat][nRS][H125]->Add(h1[v1][bl][fs][cat][nRS][H125NONVBFVHTTH]);
	}
	for(int v2=0; v2<nVarPairs; v2++){
	  h2[v2][bl][fs][cat][nRS][H125]->Add(h2[v2][bl][fs][cat][nRS][H125VBF]);
	  h2[v2][bl][fs][cat][nRS][H125]->Add(h2[v2][bl][fs][cat][nRS][H125VH]);
	  h2[v2][bl][fs][cat][nRS][H125]->Add(h2[v2][bl][fs][cat][nRS][H125TTH]);
	  h2[v2][bl][fs][cat][nRS][H125]->Add(h2[v2][bl][fs][cat][nRS][H125NONVBFVHTTH]);
	}
	//}
      }
    }
    for(int pr=0; pr<nProcesses; pr++){
      for(int fs=0; fs<nFinalStates; fs++){
	for(int cat=0; cat<nCat; cat++){
	  //for(int rs=0; rs<nRS; rs++){
	  for(int v1=0; v1<nVariables; v1++){
	    h1[v1][bl][nFinalStates][cat][nRS][pr]->Add(h1[v1][bl][fs][cat][nRS][pr]);
	  }
	  for(int v2=0; v2<nVarPairs; v2++){
	    h2[v2][bl][nFinalStates][cat][nRS][pr]->Add(h2[v2][bl][fs][cat][nRS][pr]);
	    if(pr==Data){
	      g2DataX[v2][bl][nFinalStates][cat][nRS].insert(g2DataX[v2][bl][nFinalStates][cat][nRS].end(), g2DataX[v2][bl][fs][cat][nRS].begin(), g2DataX[v2][bl][fs][cat][nRS].end());
	      g2DataY[v2][bl][nFinalStates][cat][nRS].insert(g2DataY[v2][bl][nFinalStates][cat][nRS].end(), g2DataY[v2][bl][fs][cat][nRS].begin(), g2DataY[v2][bl][fs][cat][nRS].end());
	      g2DataEX[v2][bl][nFinalStates][cat][nRS].insert(g2DataEX[v2][bl][nFinalStates][cat][nRS].end(), g2DataEX[v2][bl][fs][cat][nRS].begin(), g2DataEX[v2][bl][fs][cat][nRS].end());
	      g2DataEY[v2][bl][nFinalStates][cat][nRS].insert(g2DataEY[v2][bl][nFinalStates][cat][nRS].end(), g2DataEY[v2][bl][fs][cat][nRS].begin(), g2DataEY[v2][bl][fs][cat][nRS].end());
	    }
	  }
	  //}
	}
      }
      for(int fs=0; fs<nFinalStates+1; fs++){
	for(int cat=0; cat<nCat; cat++){
	  //for(int rs=0; rs<nRS; rs++){
	  for(int v1=0; v1<nVariables; v1++){
	    h1[v1][bl][fs][nCat][nRS][pr]->Add(h1[v1][bl][fs][cat][nRS][pr]);
	  }
	  for(int v2=0; v2<nVarPairs; v2++){
	    h2[v2][bl][fs][nCat][nRS][pr]->Add(h2[v2][bl][fs][cat][nRS][pr]);
	    if(pr==Data){
	      g2DataX[v2][bl][fs][nCat][nRS].insert(g2DataX[v2][bl][fs][nCat][nRS].end(), g2DataX[v2][bl][fs][cat][nRS].begin(), g2DataX[v2][bl][fs][cat][nRS].end());
	      g2DataY[v2][bl][fs][nCat][nRS].insert(g2DataY[v2][bl][fs][nCat][nRS].end(), g2DataY[v2][bl][fs][cat][nRS].begin(), g2DataY[v2][bl][fs][cat][nRS].end());
	      g2DataEX[v2][bl][fs][nCat][nRS].insert(g2DataEX[v2][bl][fs][nCat][nRS].end(), g2DataEX[v2][bl][fs][cat][nRS].begin(), g2DataEX[v2][bl][fs][cat][nRS].end());
	      g2DataEY[v2][bl][fs][nCat][nRS].insert(g2DataEY[v2][bl][fs][nCat][nRS].end(), g2DataEY[v2][bl][fs][cat][nRS].begin(), g2DataEY[v2][bl][fs][cat][nRS].end());
	    }
	  }
	  //}
	}
	//for(int cat=0; cat<nCat+1; cat++){
	//  for(int rs=0; rs<nRS; rs++){
	//    for(int v1=0; v1<nVariables; v1++){
	//      h1[v1][bl][fs][cat][nRS][pr]->Add(h1[v1][bl][fs][cat][rs][pr]);
	//    }
	//    for(int v2=0; v2<nVarPairs; v2++){
	//      h2[v2][bl][fs][cat][nRS][pr]->Add(h2[v2][bl][fs][cat][rs][pr]);
	//      if(pr==Data){
	//	g2DataX[v2][bl][fs][cat][nRS].insert(g2DataX[v2][bl][fs][cat][nRS].end(), g2DataX[v2][bl][fs][cat][rs].begin(), g2DataX[v2][bl][fs][cat][rs].end());
	//	g2DataY[v2][bl][fs][cat][nRS].insert(g2DataY[v2][bl][fs][cat][nRS].end(), g2DataY[v2][bl][fs][cat][rs].begin(), g2DataY[v2][bl][fs][cat][rs].end());
	//	g2DataEX[v2][bl][fs][cat][nRS].insert(g2DataEX[v2][bl][fs][cat][nRS].end(), g2DataEX[v2][bl][fs][cat][rs].begin(), g2DataEX[v2][bl][fs][cat][rs].end());
	//	g2DataEY[v2][bl][fs][cat][nRS].insert(g2DataEY[v2][bl][fs][cat][nRS].end(), g2DataEY[v2][bl][fs][cat][rs].begin(), g2DataEY[v2][bl][fs][cat][rs].end());
	//      }
	//    }
	//  }
	//}
      }
    }
  }
  

  //---------- Write histograms to a ROOT file
  string outFileName = string(Form("histos_plotDataVsMC_%.3ffb-1%s.root",lumi,(MERGE2E2MU?"_m":"")));
  cout<<"Writing MC histograms and data graphs into file "<<outFileName<<" ..."<<endl;
  TFile* fOutHistos = new TFile(outFileName.c_str(),"recreate");
  fOutHistos->cd();
  for(int bl=0; bl<nBlindings; bl++){
    if(bl!=BLINDING) continue;
    for(int fs=0; fs<nFinalStates+1; fs++){
      //cout<<" Writing FS "<<sFinalState[fs]<<endl;
      for(int cat=0; cat<nCat+1; cat++){
	//cout<<" Writing cat. "<<sCategory[cat]<<endl;
	for(int rs=0; rs<nRS+1; rs++){
	  if(rs!=nRS) continue;
	  for(int pr=0; pr<nProcesses; pr++){
	    for(int v1=0; v1<nVariables; v1++){
	      if(myV1[v1].plotLvl>3) continue;
	      if( (fs==nFinalStates || myV1[v1].InFS) &&
		  (cat==myV1[v1].Categ) ){
		h1[v1][bl][fs][cat][rs][pr]->Write(h1[v1][bl][fs][cat][rs][pr]->GetName());
		delete h1[v1][bl][fs][cat][rs][pr];
	      }
	    }
	    for(int v2=0; v2<nVarPairs; v2++){
	      if(myV2[v2].plotLvl>3) continue;
	      if(fs==nFinalStates && cat==nCat){	      
		h2[v2][bl][fs][cat][rs][pr]->Write(h2[v2][bl][fs][cat][rs][pr]->GetName());
		delete h2[v2][bl][fs][cat][rs][pr];
	      }
	    }
	  }
	  for(int v2=0; v2<nVarPairs; v2++){
	    if(MERGE2E2MU && fs==fs2mu2e) continue;
	    g2Data[v2][bl][fs][cat][rs] = new TGraphErrors( g2DataX[v2][bl][fs][cat][rs].size(),
							    &(g2DataX[v2][bl][fs][cat][rs][0]),
							    &(g2DataY[v2][bl][fs][cat][rs][0]),
							    &(g2DataEX[v2][bl][fs][cat][rs][0]),
							    &(g2DataEY[v2][bl][fs][cat][rs][0]) );
	    g2Data[v2][bl][fs][cat][rs]->Write(Form("g2Data_%s_%s_%s_%s_%s",myV2[v2].Name.c_str(),sBlinding[bl].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str(),sResonantStatus[rs].c_str()));
	    delete g2Data[v2][bl][fs][cat][rs];
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
    if(fabs(LepEta)<1.479)
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

void doHistogramsZPlusXSS(string inputPathDataForCR, string inputFileFakeRates, double lumi)
{

  cout<<"Preparing Z+X histograms (SS method) from "<<inputPathDataForCR<<endl;

  TFile* fFakeRates = TFile::Open(inputFileFakeRates.c_str());
  gr_FRmu_EB = (TGraph*)fFakeRates->Get("FR_SS_muon_EB");
  gr_FRmu_EE = (TGraph*)fFakeRates->Get("FR_SS_muon_EE");
  gr_FRel_EB = (TGraph*)fFakeRates->Get("FR_SS_electron_EB");
  gr_FRel_EE = (TGraph*)fFakeRates->Get("FR_SS_electron_EE");

  Int_t nRun;
  Long64_t nEvent;
  Int_t nLumi;
  Float_t PFMET;
  Int_t CRflag;
  Short_t ZZsel;
  Float_t ZZMass;
  Float_t ZZMassRefit;
  Float_t ZZPt;
  Float_t ZZEta;
  Float_t p_GG_SIG_ghg2_1_ghz1_1_JHUGen;
  Float_t p_QQB_BKG_MCFM;
  Float_t p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal;
  Float_t p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal;
  Float_t p_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
  Float_t pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
  Float_t p_JQCD_SIG_ghg2_1_JHUGen_JECNominal;
  Float_t p_HadWH_SIG_ghw1_1_JHUGen_JECNominal;
  Float_t p_HadZH_SIG_ghz1_1_JHUGen_JECNominal;
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
  vector<Float_t> *JetPhi = 0;
  vector<Float_t> *JetQGLikelihood = 0;
  Float_t jetPhi[99];
  Float_t jetQGL[99];
  Float_t jetPgOverPq[99];

  TFile* dataFile = TFile::Open((inputPathDataForCR+"AllData/ZZ4lAnalysis.root").c_str());
  TTree* mytree = (TTree*)dataFile->Get("CRZLLTree/candTree");

  mytree->SetBranchAddress("RunNumber", &nRun);
  mytree->SetBranchAddress("EventNumber", &nEvent);
  mytree->SetBranchAddress("LumiNumber", &nLumi);
  mytree->SetBranchAddress("PFMET", &PFMET);
  mytree->SetBranchAddress("CRflag", &CRflag);
  mytree->SetBranchAddress("ZZsel", &ZZsel);
  mytree->SetBranchAddress("ZZMass", &ZZMass);
  mytree->SetBranchAddress("ZZMassRefit", &ZZMassRefit);
  mytree->SetBranchAddress("ZZPt", &ZZPt);
  mytree->SetBranchAddress("ZZEta", &ZZEta);
  mytree->SetBranchAddress("p_GG_SIG_ghg2_1_ghz1_1_JHUGen", &p_GG_SIG_ghg2_1_ghz1_1_JHUGen);
  mytree->SetBranchAddress("p_QQB_BKG_MCFM", &p_QQB_BKG_MCFM);
  mytree->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal", &p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal);
  mytree->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal", &p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal);
  mytree->SetBranchAddress("p_JVBF_SIG_ghv1_1_JHUGen_JECNominal", &p_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
  mytree->SetBranchAddress("pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal", &pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
  mytree->SetBranchAddress("p_JQCD_SIG_ghg2_1_JHUGen_JECNominal", &p_JQCD_SIG_ghg2_1_JHUGen_JECNominal);
  mytree->SetBranchAddress("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal", &p_HadWH_SIG_ghw1_1_JHUGen_JECNominal);
  mytree->SetBranchAddress("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal", &p_HadZH_SIG_ghz1_1_JHUGen_JECNominal);
  mytree->SetBranchAddress("p_HadWH_mavjj_JECNominal", &p_HadWH_mavjj_JECNominal);
  mytree->SetBranchAddress("p_HadZH_mavjj_JECNominal",&p_HadZH_mavjj_JECNominal);
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
  mytree->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSF", &nJetsBTagged);
  mytree->SetBranchAddress("JetPhi", &JetPhi);
  mytree->SetBranchAddress("JetQGLikelihood", &JetQGLikelihood);

  TH1F* h1[nVariables][nBlindings][nFinalStates+1][nCat+1];
  for(int bl=0; bl<nBlindings; bl++){
    for(int fs=0; fs<nFinalStates+1; fs++){
      for(int cat=0; cat<nCat+1; cat++){
	for(int v1=0; v1<nVariables; v1++){
	  h1[v1][bl][fs][cat] = new TH1F(
	     Form("h1_ZPlusXSS_%s_%s_%s_%s",myV1[v1].Name.c_str(),sBlinding[bl].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str()),
	     Form(";%s;%s",myV1[v1].XLabel.c_str(),myV1[v1].YLabel.c_str()),
	     myV1[v1].Nbin,myV1[v1].Min,myV1[v1].Max);
	}
      }
    }
  }

  Float_t expectedYieldSR[nFinalStates+1][nCat+1];
  Int_t NumberOfEventsCR[nFinalStates+1][nCat+1];
  for(int fs=0; fs<nFinalStates+1; fs++){
    for(int cat=0; cat<nCat+1; cat++){
      expectedYieldSR[fs][cat] = 0.;
      NumberOfEventsCR[fs][cat] = 0.;
    }
  }

  Int_t currentFinalState;
  Int_t currentCategory;

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


    //----- find category
    
    for(int j=0; j<nJets; j++){
      jetPhi[j] = JetPhi->at(j);
      jetQGL[j] = JetQGLikelihood->at(j);
      jetPgOverPq[j] = 1./JetQGLikelihood->at(j) - 1.;
    }
    currentCategory = categoryMor17(
       nExtraLep,
       nExtraZ,
       nJets,
       nJetsBTagged,
       jetQGL,
       p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
       p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
       p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
       p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
       pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
       p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
       p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
		 p_HadWH_mavjj_JECNominal,
		 p_HadZH_mavjj_JECNominal,
       jetPhi,
       ZZMass,
       PFMET,
       true,
       false
       );


    //----- fill histograms, update counters

    Float_t KD = p_GG_SIG_ghg2_1_ghz1_1_JHUGen / ( p_GG_SIG_ghg2_1_ghz1_1_JHUGen + p_QQB_BKG_MCFM*getDbkgkinConstant(Z1Flav*Z2Flav,ZZMass) );
    Float_t D2jet = (nJets>=2) ? DVBF2j_ME(p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass) : -2;
    Float_t D1jet = (nJets==1) ? DVBF1j_ME(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass) : -2;
    Float_t DWH = (nJets>=2) ? DWHh_ME(p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_HadWH_mavjj_JECNominal, ZZMass) : -2;
    Float_t DZH = (nJets>=2) ? DZHh_ME(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_HadZH_mavjj_JECNominal, ZZMass) : -2;
    Float_t oldCConstD2jet = getDVBF2jetsConstant(ZZMass);
    Float_t oldCConstD1jet = getDVBF1jetConstant(ZZMass);
    Float_t oldCConstDWH = getDWHhConstant(ZZMass);
    Float_t oldCConstDZH = getDZHhConstant(ZZMass);
    Float_t newCConstD2jet = getDVBF2jetsConstant_shiftWP(ZZMass,false,NEWWP2J);
    Float_t newCConstD1jet = getDVBF1jetConstant_shiftWP(ZZMass,false,NEWWP1J);
    Float_t newCConstDWH = getDWHhConstant_shiftWP(ZZMass,false,NEWWPWH);
    Float_t newCConstDZH = getDZHhConstant_shiftWP(ZZMass,false,NEWWPZH);
    D2jet = 1/(newCConstD2jet/oldCConstD2jet*(1/D2jet-1)+1);
    D1jet = 1/(newCConstD1jet/oldCConstD1jet*(1/D1jet-1)+1);
    DWH = 1/(newCConstDWH/oldCConstDWH*(1/DWH-1)+1);
    DZH = 1/(newCConstDZH/oldCConstDZH*(1/DZH-1)+1);
    Float_t DVH = fmax(DWH,DZH);

    Float_t varVal[nVariables];
    Bool_t varPassCut[nVariables];
    for(int v1=0; v1<nVariables; v1++){
      string varStr = myV1[v1].Name.substr(0,myV1[v1].Name.find("_"));
      if     (varStr=="M4l") varVal[v1] = ZZMass;
      else if(varStr=="MZ1") varVal[v1] = Z1Mass;
      else if(varStr=="MZ2") varVal[v1] = Z2Mass;
      else if(varStr=="KD") varVal[v1] = KD;
      else if(varStr=="D2jet") varVal[v1] = D2jet;
      else if(varStr=="D1jet") varVal[v1] = D1jet;
      else if(varStr=="DWH") varVal[v1] = DWH;
      else if(varStr=="DZH") varVal[v1] = DZH;
      else if(varStr=="DVH") varVal[v1] = DVH;
      else if(varStr=="MET") varVal[v1] = PFMET;
      else if(varStr=="NExtraLep") varVal[v1] = (Float_t)nExtraLep;
      else if(varStr=="NJets") varVal[v1] = (Float_t)nJets;
      else if(varStr=="NBtags") varVal[v1] = (Float_t)nJetsBTagged;
      else if(varStr=="Pt4l") varVal[v1] = ZZPt;
      else if(varStr=="Eta4l") varVal[v1] = ZZEta;
      else if(varStr=="M4lrefit") varVal[v1] = ZZMassRefit;
      varPassCut[v1] = true;
      if(myV1[v1].Name.find("M4L118130")!=string::npos) varPassCut[v1] = varPassCut[v1] && 118<=ZZMass && ZZMass<=130; 
      if(myV1[v1].Name.find("M4L70110")!=string::npos) varPassCut[v1] = varPassCut[v1] && 70<=ZZMass && ZZMass<=110; 
      if(myV1[v1].Name.find("M4Labove100")!=string::npos) varPassCut[v1] = varPassCut[v1] && 100<=ZZMass; 
      if(myV1[v1].Name.find("M4Labove180")!=string::npos) varPassCut[v1] = varPassCut[v1] && 180<=ZZMass; 
      if(myV1[v1].Name.find("HighKD")!=string::npos) varPassCut[v1] = varPassCut[v1] && KD>0.5;
      if(varStr=="D2jet" || varStr=="DWH" || varStr=="DZH") varPassCut[v1] = varPassCut[v1] && nJets>=2; 
      else if(varStr=="D1jet") varPassCut[v1] = varPassCut[v1] && nJets==1; 
    }
    
    bool varPassBlindingCut[nBlindings] = {
      true,
      ZZMass < xHistoBlindLow[1],
      ZZMass > xHistoBlindUp[2],
      ZZMass < xHistoBlindLow[3] || ZZMass  > xHistoBlindUp[3],
      ZZMass < xHistoBlindLow[4] || (ZZMass > xHistoBlindUp[4] && ZZMass < xHistoBlind2Low[4]),
      true,
    };

    Float_t yieldSR = fsROSSS[currentFinalState] * fakeRate13TeV(LepPt->at(2),LepEta->at(2),LepLepId->at(2)) * fakeRate13TeV(LepPt->at(3),LepEta->at(3),LepLepId->at(3));

    for(int bl=0; bl<nBlindings; bl++){
      for(int v1=0; v1<nVariables; v1++){
	if( myV1[v1].Name.find("M4l")==0 || (myV1[v1].Name.find("M4l")!=0 && varPassBlindingCut[bl]) ){
	  if(varPassCut[v1]){
	    h1[v1][bl][(MERGE2E2MU&&currentFinalState==fs2mu2e)?fs2e2mu:currentFinalState][currentCategory]->Fill(varVal[v1], yieldSR);
	  }
	}
      }
    }

    expectedYieldSR[currentFinalState][currentCategory] += yieldSR;
    NumberOfEventsCR[currentFinalState][currentCategory]++;
		
  }

  //---------- Fill 'inclusive' histograms
  for(int v1=0; v1<nVariables; v1++){
    for(int bl=0; bl<nBlindings; bl++){
      for(int fs=0; fs<nFinalStates; fs++){
	for(int cat=0; cat<nCat; cat++){
	  h1[v1][bl][nFinalStates][cat]->Add(h1[v1][bl][fs][cat]);
	}
      }
      for(int fs=0; fs<nFinalStates+1; fs++){
	for(int cat=0; cat<nCat; cat++){
	  h1[v1][bl][fs][nCat]->Add(h1[v1][bl][fs][cat]);
	}
      }
    }
  }
  for(int fs=0; fs<nFinalStates; fs++){
    for(int cat=0; cat<nCat; cat++){
      expectedYieldSR[nFinalStates][cat] += expectedYieldSR[fs][cat];
      NumberOfEventsCR[nFinalStates][cat] += NumberOfEventsCR[fs][cat];
    }
  }
  for(int fs=0; fs<nFinalStates+1; fs++){
    for(int cat=0; cat<nCat; cat++){
      expectedYieldSR[fs][nCat] += expectedYieldSR[fs][cat];
      NumberOfEventsCR[fs][nCat] += NumberOfEventsCR[fs][cat];
    }
  }

  //---------- Print Z+X expected yields
  if(VERBOSE2){
    for(int fs=0; fs<nFinalStates; fs++){
      cout<<fsLabelForSS[fs]<<" : "
	  <<expectedYieldSR[fs][nCat]
	  <<" +/- "<<expectedYieldSR[fs][nCat]/sqrt(NumberOfEventsCR[fs][nCat])
	  <<" (stat., evt: "<<NumberOfEventsCR[fs][nCat]<<")" 
	  <<endl;
    }
    cout<<"Total: "<<expectedYieldSR[nFinalStates][nCat]<<endl;
  }

  //---------- Smooth histogram if requested
  if(SMOOTHZPLUSXHISTOFROMSSCR){
    Float_t integral = 0;
    for(int v1=0; v1<nVariables; v1++){
      if(myV1[v1].SmoothZpX){
	for(int bl=0; bl<nBlindings; bl++){
	  for(int fs=0; fs<nFinalStates+1; fs++){
	    for(int cat=0; cat<nCat+1; cat++){
	      integral = h1[v1][bl][fs][cat]->Integral();
	      h1[v1][bl][fs][cat]->Smooth(1);
	      h1[v1][bl][fs][cat]->Scale( integral / h1[v1][bl][fs][cat]->Integral() );
	    }
	  }
	}
      }
    }
  }
  
  if(MERGE2E2MU)
    for(int cat=0; cat<nCat+1; cat++)
      expectedYieldSR[fs2e2mu][cat] += expectedYieldSR[fs2mu2e][cat];

  if(! RENORMZPLUSXTOOFFICIALCATEG){
    // Refill global array
    for(int cat=0; cat<nCat+1; cat++){
      normSSFullRange4e[cat] = expectedYieldSR[fs4e][cat];
      normSSFullRange4mu[cat] = expectedYieldSR[fs4mu][cat];
      normSSFullRange2e2mu[cat] = expectedYieldSR[fs2e2mu][cat];
    }
  }
  Float_t normSSFullRange[nFinalStates+1][nCat+1];
  for(int cat=0; cat<nCat+1; cat++){
    normSSFullRange[fs4e][cat] = normSSFullRange4e[cat];
    normSSFullRange[fs4mu][cat] = normSSFullRange4mu[cat];
    normSSFullRange[fs2e2mu][cat] = normSSFullRange2e2mu[cat];
    normSSFullRange[nFinalStates][cat] = normSSFullRange4e[cat] + normSSFullRange4mu[cat] + normSSFullRange2e2mu[cat] ;
  }
  if(RENORMZPLUSXTOOFFICIALCATEG){
    // Normalize to externally provided Z+X yield in categories from global arrays
    if(MERGE2E2MU){
      for(int bl=0; bl<nBlindings; bl++){
	for(int fs=0; fs<nFinalStates+1; fs++){
	  if(fs==fs2mu2e) continue;
	  for(int cat=0; cat<nCat+1; cat++){
	    for(int v1=0; v1<nVariables; v1++){
	      h1[v1][bl][fs][cat]->Scale( normSSFullRange[fs][cat] / expectedYieldSR[fs][cat] );
	    }
	  }
	}
      }
    }else{
      cout<<"WARNING: cannot renormalize Z+X histograms to official numbers when treating 2e2mu and 2mu2e separately"<<endl;
    }
  }
  if(RENORMZPLUSXTOOFFICIALINCL){
    // Normalize to externally provided inclusive Z+X yield from global variables
    if(MERGE2E2MU){
      Float_t normZPlusXFullSR[nFinalStates+1];
      normZPlusXFullSR[fs4e] = normZPlusXFullSR4e;
      normZPlusXFullSR[fs4mu] = normZPlusXFullSR4mu;
      normZPlusXFullSR[fs2e2mu] = normZPlusXFullSR2e2mu;
      normZPlusXFullSR[nFinalStates] = normZPlusXFullSR4e + normZPlusXFullSR4mu + normZPlusXFullSR2e2mu ;
      for(int bl=0; bl<nBlindings; bl++){
	for(int fs=0; fs<nFinalStates+1; fs++){
	  if(fs==fs2mu2e) continue;
	  for(int cat=0; cat<nCat+1; cat++){
	    for(int v1=0; v1<nVariables; v1++){
	      h1[v1][bl][fs][cat]->Scale( normZPlusXFullSR[fs] / normSSFullRange[fs][nCat] );
	    }
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
      for(int cat=0; cat<nCat+1; cat++){
	//if(fs==nFinalStates){
	for(int v1=0; v1<nVariables; v1++){
	  h1[v1][bl][fs][cat]->Write(h1[v1][bl][fs][cat]->GetName());
	  delete h1[v1][bl][fs][cat];
	}
	//}
      }
    }
  }
  fOutHistos->Close();
  delete fOutHistos;

}



void getM4lZPlusXHistoFromAnalyticalShape_InCateg(TH1F** hZPX4mu, TH1F** hZPX4e, TH1F** hZPX2e2mu, TH1F** hZPX4l, Int_t nbins, Float_t xmin, Float_t xmax, Int_t v1) {

  /*// For Moriond'16: Z+X shapes sent by Pedja on March 1st 2016 (take the same for all categ) 
  TF1 *f4eComb = new TF1("f4eComb", "landau(0)*(1 + exp( pol1(3))) + [5]*(TMath::Landau(x, [6], [7]))", 70, 1000);
  TF1 *f4muComb = new TF1("f4muComb","landau(0)",70,1000);
  TF1 *f2e2muComb = new TF1("f2e2muComb","landau(0)",70,1000);
  f4eComb->SetParameters(4.404e-05,151.2,36.6,7.06,-0.00497,0.01446,157.3,26.00);
  f4muComb->SetParameters(0.04276,134.6,24.4);
  f2e2muComb->SetParameters(0.04130,144.5,25.3);
  //*/

  //*// For ICHEP'16: Z+X shapes sent by Pedja on July 27th 2016 (take the same for all categ)
  TF1 *f4eComb    = new TF1("f4eComb"   ,"TMath::Landau(x, 141.9, 21.3)", 70, 3000);
  TF1 *f4muComb   = new TF1("f4muComb"  ,"TMath::Landau(x, 130.4, 15.6)", 70, 3000);
  TF1 *f2e2muComb = new TF1("f2e2muComb","0.45*TMath::Landau(x, 131.1, 18.1) + 0.55*TMath::Landau(x, 133.8, 18.9)", 70, 3000);
  //*/ 

  Int_t nentries = 1000000; //tried more entries but it takes too long

  for(int cat=0; cat<nCat+1; cat++){

    //----- compute normalization of the subrange of interest
    Float_t norm4mu   = normSSFullRange4mu  [cat] * f4muComb  ->Integral(xmin,xmax) / f4muComb  ->Integral(70,3000);
    Float_t norm4e    = normSSFullRange4e   [cat] * f4eComb   ->Integral(xmin,xmax) / f4eComb   ->Integral(70,3000);
    Float_t norm2e2mu = normSSFullRange2e2mu[cat] * f2e2muComb->Integral(xmin,xmax) / f2e2muComb->Integral(70,3000);

    //---------- Normalize to SS/OS-combined inclusive Z+X yields from global variables
    if(RENORMZPLUSXTOOFFICIALINCL){
      if(MERGE2E2MU){
	norm4mu   *= normZPlusXFullSR4mu   / normSSFullRange4mu  [nCat] ;
	norm4e    *= normZPlusXFullSR4e    / normSSFullRange4e   [nCat] ;
	norm2e2mu *= normZPlusXFullSR2e2mu / normSSFullRange2e2mu[nCat] ;
      }else{
	cout<<"WARNING: cannot renormalize Z+X histograms to official numbers when treating 2e2mu and 2mu2e separately"<<endl;
      }
    }

    //----- build final histograms
    hZPX4mu  [cat] = new TH1F(Form("h4mu_%s_%s"  ,sCategory[cat].c_str(),myV1[v1].Name.c_str()),";;",nbins,xmin,xmax);
    hZPX4e   [cat] = new TH1F(Form("h4e_%s_%s"   ,sCategory[cat].c_str(),myV1[v1].Name.c_str()),";;",nbins,xmin,xmax);
    hZPX2e2mu[cat] = new TH1F(Form("h2e2mu_%s_%s",sCategory[cat].c_str(),myV1[v1].Name.c_str()),";;",nbins,xmin,xmax);
    hZPX4mu  [cat]->FillRandom("f4muComb"  ,nentries);
    hZPX4e   [cat]->FillRandom("f4eComb"   ,nentries);
    hZPX2e2mu[cat]->FillRandom("f2e2muComb",nentries);
    hZPX4mu  [cat]->Scale(norm4mu   / hZPX4mu  [cat]->Integral());
    hZPX4e   [cat]->Scale(norm4e    / hZPX4e   [cat]->Integral());
    hZPX2e2mu[cat]->Scale(norm2e2mu / hZPX2e2mu[cat]->Integral());
    hZPX4l[cat] = new TH1F(Form("h4l_%s_%s"   ,sCategory[cat].c_str(),myV1[v1].Name.c_str()),";;",nbins,xmin,xmax);
    hZPX4l[cat]->Add(hZPX4mu  [cat]);
    hZPX4l[cat]->Add(hZPX4e   [cat]);
    hZPX4l[cat]->Add(hZPX2e2mu[cat]);

  }

}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void DrawDataMC1D(TCanvas* c, TH1F** h, TH1F* hZPX, int v1, int bl, double lumi, string lumiText, Int_t category, Bool_t logX = false, Bool_t logY = false) {

  bool drawData = !( (MASKDATAFORHIGHMASS && bl==blindbelow130) || bl==fullyblind );
  bool maskH125 = 
    ( MASKH125FORLOWMASS && bl==blindabove114 && (myV1[v1].Name.find("M4l")==string::npos || myV1[v1].Name.find("M4l_70110")==0 || myV1[v1].Name.find("M4l_70109")==0) ) ||
    ( MASKH125FORHIGHMASS && bl==blindbelow130 && (myV1[v1].Name.find("M4l")==string::npos || myV1[v1].Name.find("M4l_above150")==0) ) ;
  bool withRatioPlot = DRAWDATAMCRATIO && drawData && myV1[v1].WithRatio;
  bool doBlindingHisto = xHistoBlindLow[bl]<xHistoBlindUp[bl] && myV1[v1].Name.find("M4l")==0;
  bool doBlindingHisto2 = xHistoBlind2Low[bl]<xHistoBlind2Up[bl] && myV1[v1].Name.find("M4l")==0;
  bool withWP = DRAWWP1D && myV1[v1].ValWP!=0.;

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
    if((useProcess[pr] && !((maskH125 || myV1[v1].sepVbf || myV1[v1].sepVh || myV1[v1].sepTth) && pr==H125)) || (myV1[v1].sepVbf && (pr==H125VBF || pr==H125NONVBFVHTTH)) || (myV1[v1].sepVh && (pr==H125VH || pr==H125NONVBFVHTTH)) || (myV1[v1].sepTth && (pr==H125TTH || pr==H125NONVBFVHTTH))){
      if(REBINDYANDTTBAR && (pr==DY || pr==ttbar)) h[pr] = Smooth(h[pr],myV1[v1].rebinDYTT);
      hStacks[pr] = (TH1F*)h[pr]->Clone();
      if(myV1[v1].sepVbf && pr==H125NONVBFVHTTH){
	hStacks[pr]->Add((TH1F*)h[H125VH]->Clone());
	hStacks[pr]->Add((TH1F*)h[H125TTH]->Clone());
      }else if(myV1[v1].sepVh && pr==H125NONVBFVHTTH){
	hStacks[pr]->Add((TH1F*)h[H125VBF]->Clone());
	hStacks[pr]->Add((TH1F*)h[H125TTH]->Clone());
      }else if(myV1[v1].sepTth && pr==H125NONVBFVHTTH){
	hStacks[pr]->Add((TH1F*)h[H125VBF]->Clone());
	hStacks[pr]->Add((TH1F*)h[H125VH]->Clone());
      }
      if(!first) hStacks[pr]->Add(hStacks[previous]);
      hStacks[pr]->SetFillColor(processFillColor[pr]);
      hStacks[pr]->SetLineColor(DRAWLINES?processLineColor[pr]:processFillColor[pr]);
      hStacks[pr]->SetLineWidth(LINEWIDTH);
      if(!DRAWLINES) hStacks[pr]->SetLineColorAlpha(processFillColor[pr],0.);
      hStacks[pr]->GetXaxis()->SetTitleOffset(myV1[v1].Large?1.1:1.1);
      hStacks[pr]->GetYaxis()->SetTitleOffset(myV1[v1].Large?1.:1.2);
      hStacks[pr]->GetXaxis()->SetTitleSize(0.05);
      hStacks[pr]->GetYaxis()->SetTitleSize(0.05);
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
  Bool_t useZPlusX = BUILDZPLUSXHISTOFROMSSCR ||
    (USEZPLUSXANALYTICALSHAPE && (myV1[v1].Name.find("M4l")==0 /*&& myV1[v1].Name.find("Refit")==string::npos*/ && myV1[v1].CutLabel=="") );
  useZPlusX = useZPlusX && hZPX->Integral()>0.;
  if(useZPlusX){
    hZPlusX = hZPX;
    if(logY) hZPlusX->SetMinimum(1e-20);
    hZPlusX->SetFillColor(processFillColor[DY]);
    hZPlusX->SetLineColor(DRAWLINES?processLineColor[DY]:processFillColor[DY]);
    hZPlusX->SetLineWidth(LINEWIDTH);
    for(int pr=nProcesses-1; pr>=1; pr--)
      if((useProcess[pr] && !((maskH125 || myV1[v1].sepVbf || myV1[v1].sepVh || myV1[v1].sepTth) && pr==H125)) || (myV1[v1].sepVbf && (pr==H125VBF || pr==H125NONVBFVHTTH)) || (myV1[v1].sepVh && (pr==H125VH || pr==H125NONVBFVHTTH)) || (myV1[v1].sepTth && (pr==H125TTH || pr==H125NONVBFVHTTH)))
	hStacks[pr]->Add(hZPlusX);
    //cout<<"full expected yield for variable "<<myV1[v1].Name<<": "<<hStacks[idxSumMC]->Integral()<<endl;
    //cout<<"Z+X yield for variable "<<myV1[v1].Name<<": "<<hZPlusX->Integral()<<endl;
  }

  //----- prepare data graph
  h[0]->SetMarkerStyle(20);
  h[0]->SetMarkerColor(kBlack);
  h[0]->SetMarkerSize(1.);
  TGraphAsymmErrors* gData = getDataGraph(h[0]);
  if(h[0]->GetNbinsX()>39) gData->SetMarkerSize(0.8);
  if(h[0]->GetNbinsX()>99) gData->SetMarkerSize(0.6);

  //----- adjust Y axis
  const int npoints = gData->GetN();
  Float_t gDataErrorBarUp[npoints];
  for(int i=0; i<npoints; i++) gDataErrorBarUp[i] = gData->GetY()[i] + gData->GetEYhigh()[i] ;
  Float_t cmax;
  if(drawData && npoints>0) cmax = TMath::Max( (Float_t)hStacks[idxSumMC]->GetMaximum(), (Float_t)TMath::MaxElement(npoints,gDataErrorBarUp) );
  else cmax = (Float_t)hStacks[idxSumMC]->GetMaximum();
  cmax *= logY ? 30. : 1.1 ;
  cmax *= myV1[v1].MaxCorrector[bl];
  Float_t cminlog = hStacks[idxSumMC]->GetMaximum() / myV1[v1].MinFactor;
  hStacks[idxSumMC]->SetMaximum(cmax);
  hStacks[idxSumMC]->SetMinimum(logY?cminlog:0.);
  
  //----- prepare grey area/grid for blind region
  TH1F* hBlind = 0;
  if(doBlindingHisto){
    hBlind = new TH1F(Form("hBlind_%s_%s",myV1[v1].Name.c_str(),sBlinding[bl].c_str()),";;",1,xHistoBlindLow[bl],xHistoBlindUp[bl]);
    hBlind->SetBinContent(1,cmax);
    //hBlind->SetFillColorAlpha(kBlack,0.2);
    //hBlind->SetLineColorAlpha(kBlack,0.2);
    hBlind->SetFillColor(kGray+3);
    hBlind->SetFillStyle(3013); //also tried 3001 and a few others, but they look bad in pdf format
    hBlind->SetLineColorAlpha(kWhite,0.);
  }
  TH1F* hBlind2 = 0;
  if(doBlindingHisto2){
    hBlind2 = new TH1F(Form("hBlind2_%s_%s",myV1[v1].Name.c_str(),sBlinding[bl].c_str()),";;",1,xHistoBlind2Low[bl],xHistoBlind2Up[bl]);
    hBlind2->SetBinContent(1,cmax);
    hBlind2->SetFillColor(kGray+3);
    hBlind2->SetFillStyle(3013);
    hBlind2->SetLineColorAlpha(kWhite,0.);
  }

  //----- prepare legend
  bool useBlindingLabel = blindingLabel[bl]!="" && myV1[v1].Name.find("M4l")==string::npos;
  float legTextSize = 0.034;
  int legPos = myV1[v1].LegPos;
  if(bl==blindabove114 && myV1[v1].Name.find("MZ")!=string::npos) legPos = 33;
  float legLeft = (legPos==11) ? 0.20 : 0.68 ;
  float legWidth = 0.19;
  float legUp = 0.9;
  if(!withRatioPlot && legPos==myV1[v1].CMSPos && !myV1[v1].Large) legUp -= 0.1;
  float legHeight = 0.16;
  if(!drawData) legHeight -= 0.04;
  if(maskH125) legHeight -= 0.04;
  if(useZPlusX) legHeight += 0.04;
  if(myV1[v1].sepVbf||myV1[v1].sepVh||myV1[v1].sepTth) legHeight += 0.04;
  if(withRatioPlot){ legHeight *= 1.4; legUp = 0.95; }
  if(myV1[v1].Large){ legWidth *= 1.15; legHeight *= 1.15; legTextSize *= 1.15; }
  TLegend* lgd = new TLegend(legLeft,legUp-legHeight,legLeft+legWidth,legUp);
  lgd->SetFillStyle(0);
  lgd->SetBorderSize(0);
  lgd->SetTextSize(legTextSize);
  if(drawData)
    lgd->AddEntry(h[0],processLabel[0].c_str(),"p");
  for(int pr=1; pr<nProcesses; pr++)
    if((useProcess[pr] && !((maskH125 || myV1[v1].sepVbf || myV1[v1].sepVh || myV1[v1].sepTth) && pr==H125)) || (myV1[v1].sepVbf && (pr==H125VBF || pr==H125NONVBFVHTTH)) || (myV1[v1].sepVh && (pr==H125VH || pr==H125NONVBFVHTTH)) || (myV1[v1].sepTth && (pr==H125TTH || pr==H125NONVBFVHTTH)))
      lgd->AddEntry(hStacks[pr],processLabel[pr].c_str(),"f");
  if(useZPlusX)
    lgd->AddEntry(hZPlusX," Z+X","f");
    

  //----- prepare label for cut/blinding
  bool useLabel = useBlindingLabel || myV1[v1].CutLabel!="";
  TPaveText* pav = 0;
  if(useLabel){
    int pavPos = myV1[v1].CutPos;
    float pavLeft = (pavPos==11) ? 0.20 : (pavPos==33) ? 0.64 : legLeft-0.02 ;
    float pavWidth = 0.27;
    float pavUp = pavPos!=0 ? 0.93 : legUp-legHeight+0.05;
    float pavHeight = 0.06;
    if(myV1[v1].Name.find("DVH")!=string::npos) pavLeft += 0.04;
    //pav = new TPaveText(legLeft-(legPos==11?0.02:0.02),legUp-legHeight-0.1,legLeft+legWidth+(legPos==11?0.06:0.06),legUp-legHeight+0.02,"brNDC");
    pav = new TPaveText(pavLeft,pavUp-legHeight,pavLeft+pavWidth,pavUp,"brNDC");
    pav->SetFillStyle(0);
    pav->SetBorderSize(0);
    pav->SetTextAlign(pavPos==11?11:pavPos==33?31:legPos==11?13:33);
    pav->SetTextSize(0.037);
    if(myV1[v1].Large) pav->SetTextSize(0.04);
    pav->SetTextFont(42);
    pav->AddText((myV1[v1].CutLabel!="")?myV1[v1].CutLabel.c_str():blindingLabel[bl].c_str());
  }

  //----- prepare label for category
  bool useLabelCat = category!=nCat ; 
  TPaveText* pavCat = 0;
  if(useLabelCat){
    float pavCatLeft = (legPos==11) ? 0.66 : 0.18 ;
    float pavCatWidth = 0.27;
    float pavCatUp = legUp+0.02-0.05*useLabel;
    float pavCatHeight = 0.12;
    pavCat = new TPaveText(pavCatLeft,pavCatUp-pavCatHeight,pavCatLeft+pavCatWidth,pavCatUp,"brNDC");
    pavCat->SetFillStyle(0);
    pavCat->SetBorderSize(0);
    pavCat->SetTextAlign(legPos==11?31:11);
    pavCat->SetTextSize(0.037);
    pavCat->SetTextFont(42);
    pavCat->AddText(categoryLabel[category].c_str());
  }

  //----- line for WP
  TLine *lineWP = 0;
  if(withWP){
    lineWP = new TLine(myV1[v1].ValWP,logY?cminlog:0.,myV1[v1].ValWP,logY?cminlog*exp(0.6*log(cmax/cminlog)):0.6*cmax);
    lineWP->SetLineStyle(9);
    lineWP->SetLineWidth(8);
    lineWP->SetLineColor(13);
  }

  //----- draw everything
  if(!withRatioPlot){

    //--- no Data/MC graph -> draw on main pad
    for(int pr=1; pr<nProcesses; pr++)
      if((useProcess[pr] && !((maskH125 || myV1[v1].sepVbf || myV1[v1].sepVh || myV1[v1].sepTth) && pr==H125)) || (myV1[v1].sepVbf && (pr==H125VBF || pr==H125NONVBFVHTTH)) || (myV1[v1].sepVh && (pr==H125VH || pr==H125NONVBFVHTTH)) || (myV1[v1].sepTth && (pr==H125TTH || pr==H125NONVBFVHTTH)))
	hStacks[pr]->Draw( pr==idxSumMC ? "HIST" : "HIST SAME" );
    if(useZPlusX) hZPlusX->Draw("HIST SAME");
    if(withWP) lineWP->Draw();
    if(drawData) gData->Draw("P");
    if(doBlindingHisto) hBlind->Draw("HIST SAME");
    if(doBlindingHisto2) hBlind2->Draw("HIST SAME");
    lgd->Draw();
    if(useLabel) pav->Draw();
    if(useLabelCat) pavCat->Draw();
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
      if((useProcess[pr] && !((maskH125 || myV1[v1].sepVbf || myV1[v1].sepVh || myV1[v1].sepTth) && pr==H125)) || (myV1[v1].sepVbf && (pr==H125VBF || pr==H125NONVBFVHTTH)) || (myV1[v1].sepVh && (pr==H125VH || pr==H125NONVBFVHTTH)) || (myV1[v1].sepTth && (pr==H125TTH || pr==H125NONVBFVHTTH))){
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
    if(useLabelCat) pavCat->Draw();
    pad1->RedrawAxis();
    
    //--- Data/MC pad
    pad2->cd();
    pad2->SetGridy();
    TH1F* hOne = (TH1F*)h[0]->Clone();
    hOne->Reset();
    for(int i=1; i<=hOne->GetNbinsX(); i++) hOne->SetBinContent(i,1.);
    hOne->SetLineColor(kGray);
    hOne->Draw();
    //hOne->GetYaxis()->SetRangeUser(0.45,1.55);
    hOne->GetYaxis()->SetRangeUser(0.4,1.85);
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
  if(CUSTOMXAXISLABELS && myV1[v1].CXL){
    c->cd();
    TText t;
    t.SetNDC();
    t.SetTextSize(0.04);
    t.SetTextFont(42);
    t.DrawText(0.830,0.093,"700");
    t.DrawText(0.907,0.093,"900");
  }

  //----- print official CMS labels and lumi
  writeExtraText = true;
  extraText  = "Preliminary";
  lumi_sqrtS = lumiText + " (13 TeV)";
  CMS_lumi( c, 0, (withRatioPlot || myV1[v1].Large) ? 0 : myV1[v1].CMSPos );

  //----- print yields
  if(VERBOSE2){
    cout<<"Yields for variable "<<myV1[v1].Name<<":"<<endl;
    cout<<"  qqZZ             "<<h[qqZZ]->Integral()<<endl;
    cout<<"  ggZZ             "<<h[ggZZ]->Integral()<<endl;
    if(useZPlusX) cout<<"  Z+X              "<<hZPlusX->Integral()<<endl;
    cout<<"  all backgrounds  "<<hStacks[idxSumMC]->Integral() - h[H125]->Integral()<<endl;
    cout<<"  signal           "<<h[H125]->Integral()<<endl;
    cout<<"  total expected   "<<hStacks[idxSumMC]->Integral()<<endl;
    cout<<"  observed         "<<h[0]->Integral()<<endl;
  }

}


void DrawDataMC2D(TCanvas* c, TH2F** h2, TGraphErrors* g2[nFinalStates+1][nCat+1], int v2, int bl, string lumiText, int style2D, Bool_t logX = false, Bool_t logY = false) {

  bool maskH125 = 
    ( MASKH125FORLOWMASS  && bl==blindabove114 && (myV2[v2].Name.find("M4l")==string::npos || myV2[v2].Name.find("M4L70110")!=string::npos) ) ||
    ( MASKH125FORHIGHMASS && bl==blindbelow130 && (myV2[v2].Name.find("M4l")==string::npos || myV2[v2].Name.find("M4L180780")!=string::npos || myV2[v2].Name.find("M4L150700")!=string::npos) ) ;
  bool categMode = CATEGIN2D && myV2[v2].WithCateg;
  float factorHeight = 1.2;
  bool withWP = DRAWWP2D && myV2[v2].exprWP!="";
  bool legendOutMode = LEGENDOUTOF2DFRAME && !categMode;

  //----- prepare canvas
  TStyle* style = (TStyle*)gROOT->GetStyle("tdrStyle")->Clone();
  style->SetPadLeftMargin(0.12);
  style->SetPadRightMargin(0.14);
  style->SetEndErrorSize(0.);
  style->cd();

  bool LEGENDWHITE = 0;
  bool signalSeparately = (style2D==7);
  if(style2D==0 || myV2[v2].Name.find("M4lVsM4lRefit")==0){
    setColZGradient_Rainbow2();
    LEGENDWHITE = myV2[v2].LegIsWhite;
  }else if(style2D==8){
    setColZGradient_TwoColors();
  }else{
    setColZGradient_OneColor(style2D,myV2[v2].Name.find("MZ1VsMZ2")==0);
  }
  bool withMassErrors = (myV2[v2].Name.find("M4l")==0 && MARKERSTYLE==0);
  //cout<<"withMassErrors = "<<withMassErrors<<endl;

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
  h2Stacked->GetYaxis()->SetRangeUser(myV2[v2].YMin,(categMode?factorHeight-0.0001:1.)*myV2[v2].YMax);
  h2Stacked->GetXaxis()->SetTitleSize(0.05);
  h2Stacked->GetYaxis()->SetTitleSize(0.05);
  h2Stacked->GetXaxis()->SetTitleOffset(1.1);
  h2Stacked->GetYaxis()->SetTitleOffset(1.1);//1.17);
  h2Stacked->GetZaxis()->SetTitleOffset(1.);//1.15);
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
  if(myV2[v2].Name.find("MZ1VsMZ2")==string::npos) h2Stacked->SetMinimum(-1e-20); // avoid white bins
  else h2Stacked->SetMinimum(+1e-20);//(+1e-5);//

  TH2F* h2Signal = 0;
  if(signalSeparately){
    h2Signal = (TH2F*)h2[H125]->Clone();
    h2Signal->SetFillColor(TColor::GetColor("#ff9090"));
  }

  //----- prepare data graphs
  if(categMode){
    for(int fs=0; fs<nFinalStates; fs++){
      if(MERGE2E2MU && fs==fs2mu2e) continue;
      for(int cat=0; cat<nCat; cat++){
	g2[fs][cat]->SetMarkerStyle(categMarkerStyle[cat]);
	g2[fs][cat]->SetMarkerSize(0.9);
	g2[fs][cat]->SetMarkerColor(fsMarkerColor[fs]);
	g2[fs][cat]->SetLineColor(fsMarkerColor[fs]);
      }
    }
  }else{
    for(int fs=0; fs<nFinalStates; fs++){
      if(MERGE2E2MU && fs==fs2mu2e) continue;
      g2[fs][nCat]->SetMarkerStyle(MARKERSTYLE==1?fsMarkerStyleOpen[fs]:fsMarkerStyleFull[fs]);
      g2[fs][nCat]->SetMarkerSize(MARKERSTYLE==1?1.:0.7);
      if(style2D==1){
	g2[fs][nCat]->SetMarkerColor(fsMarkerColor[fs]);
	g2[fs][nCat]->SetLineColor(fsMarkerColor[fs]);
      }else{
	g2[fs][nCat]->SetMarkerColor(kBlack);
	g2[fs][nCat]->SetLineColor(kBlack);
      }
    }
  }
  
  //----- prepare grid for blind region
  bool doBlindingArea = xHistoBlindLow[bl]<xHistoBlindUp[bl] 
    && myV2[v2].Name.find("M4l")==0
    && xHistoBlindLow[bl]<h2Stacked->GetXaxis()->GetXmax()
    && xHistoBlindUp [bl]>h2Stacked->GetXaxis()->GetXmin(); 
  TBox* box = 0;
  if(doBlindingArea){
    box = new TBox(TMath::Max((Float_t)xHistoBlindLow[bl],(Float_t)h2Stacked->GetXaxis()->GetXmin()),h2Stacked->GetYaxis()->GetXmin(),TMath::Min((Float_t)xHistoBlindUp[bl],(Float_t)h2Stacked->GetXaxis()->GetXmax()),h2Stacked->GetYaxis()->GetXmax());
    box->SetFillColor(kBlack);
    box->SetFillStyle(3013);
  }
  bool doBlindingArea2 = xHistoBlind2Low[bl]<xHistoBlind2Up[bl] 
    && myV2[v2].Name.find("M4l")==0
    && xHistoBlind2Low[bl]<h2Stacked->GetXaxis()->GetXmax()
    && xHistoBlind2Up [bl]>h2Stacked->GetXaxis()->GetXmin(); 
  TBox* box2 = 0;
  if(doBlindingArea2){
    box2 = new TBox(TMath::Max((Float_t)xHistoBlind2Low[bl],(Float_t)h2Stacked->GetXaxis()->GetXmin()),h2Stacked->GetYaxis()->GetXmin(),TMath::Min((Float_t)xHistoBlind2Up[bl],(Float_t)h2Stacked->GetXaxis()->GetXmax()),h2Stacked->GetYaxis()->GetXmax());
    box2->SetFillColor(kBlack);
    box2->SetFillStyle(3013);
  }

  //----- prepare legend
  int legPos = myV2[v2].LegPos;
  float legLeft = legendOutMode ? 0.87 : (legPos==11) ? 0.17 : 0.69 ;
  float legWidth = 0.12;
  float legUp = legendOutMode ? 0.95 : 0.90;
  float legHeight = 0.12;
  if(signalSeparately){
    legLeft = (legPos==11) ? 0.20 : 0.58 ;
    legWidth = 0.24;
    legHeight = 0.20;
  }
  TLegend* lgd = 0;
  if(!categMode){
    lgd = new TLegend(legLeft,legUp-legHeight,legLeft+legWidth,legUp);
    lgd->SetFillStyle(0);
    if(LEGENDWHITE && !legendOutMode){
      lgd->SetBorderSize(1);
      lgd->SetLineColor(kWhite);
      lgd->SetTextColor(kWhite);
    }else{
      lgd->SetBorderSize(1);
      lgd->SetLineColor(kBlack);
      lgd->SetTextColor(kBlack);
    }
    //if(legendOutMode) lgd->SetBorderSize(0);
    if(signalSeparately){
      lgd->AddEntry(h2Signal,processLabel[H125].c_str(),"f");
      h2Stacked->SetFillColor(TColor::GetColor("#74baff"));
      lgd->AddEntry(h2Stacked,"ZZ background","f");
    }
    TGraph* gLeg[nFinalStates];
    for(int fs=0; fs<nFinalStates; fs++){
      if(MERGE2E2MU && fs==fs2mu2e) continue;
      gLeg[fs] = (TGraph*)g2[fs][nCat]->Clone();
    }
    lgd->AddEntry(gLeg[fs4e   ],(string(signalSeparately?"Data, ":"")+" "+fsLabel[fs4e   ]).c_str(),withMassErrors?"lp":"p");
    lgd->AddEntry(gLeg[fs4mu  ],(string(signalSeparately?"Data, ":"")+" "+fsLabel[fs4mu  ]).c_str(),withMassErrors?"lp":"p");
    lgd->AddEntry(gLeg[fs2e2mu],(string(signalSeparately?"Data, ":"")+" "+fsLabel[fs2e2mu]).c_str(),withMassErrors?"lp":"p");
    if(!MERGE2E2MU) lgd->AddEntry(gLeg[fs2mu2e],(string(signalSeparately?"Data, ":"")+" "+fsLabel[fs2mu2e]).c_str(),withMassErrors?"lp":"p");
  }

  //----- prepare label for cut/blinding
  bool useLabel = (blindingLabel[bl]!="" && myV2[v2].Name.find("M4l")==string::npos) || myV2[v2].CutLabel!="";
  TPaveText* pav = 0;
  if(useLabel){
    float pavLeft = legendOutMode ? 0.52 : 0.52; //legLeft-(legPos==11?0.02:0.18) ;
    float pavWidth = legWidth + 0.16 ;
    float pavUp = legendOutMode ? 0.92 : 0.92; //legUp-legHeight-0.1 ;
    float pavHeight = 0.03 ;
    pav = new TPaveText(pavLeft,pavUp,pavLeft+pavWidth,pavUp-pavHeight,"brNDC");
    pav->SetFillStyle(0);
    pav->SetBorderSize(0);
    pav->SetTextAlign(13);
    pav->SetTextSize(0.034);
    pav->SetTextFont(42);
    pav->SetTextColor(LEGENDWHITE?kWhite:kBlack);
    pav->AddText((myV2[v2].CutLabel!="")?myV2[v2].CutLabel.c_str():blindingLabel[bl].c_str());
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
  if(withWP){
    TF1* fWP = new TF1("fWP",myV2[v2].exprWP.c_str(),myV2[v2].XMin,myV2[v2].XMax);
    fWP->SetLineStyle(7);
    fWP->SetLineWidth(4);
    fWP->SetLineColor(13);
    fWP->Draw("same");
  }
  if(bl!=fullyblind && myV2[v2].Name!="M4lVsM4lRefit"){
    if(categMode){
      TBox* lgdBox = new TBox(myV2[v2].XMin,myV2[v2].YMax,myV2[v2].XMax,factorHeight*myV2[v2].YMax);
      lgdBox->SetFillColor(kWhite);
      lgdBox->SetLineWidth(1);
      lgdBox->SetLineColor(kBlack);
      lgdBox->Draw("l");
      float legTextSize = 0.032;
      TLegend* lgd1 = new TLegend(0.14,0.82,0.27,0.94);
      TLegend* lgd2 = new TLegend(0.28,0.82,0.55,0.94);
      TLegend* lgd3 = new TLegend(0.55,0.815,0.82,0.945);
      lgd1->SetFillStyle(0);
      lgd2->SetFillStyle(0);
      lgd3->SetFillStyle(0);
      lgd1->SetBorderSize(0);
      lgd2->SetBorderSize(0);
      lgd3->SetBorderSize(0);
      lgd1->SetTextFont(42);
      lgd2->SetTextFont(42);
      lgd3->SetTextFont(42);
      lgd1->SetTextSize(legTextSize);
      lgd2->SetTextSize(legTextSize);
      lgd3->SetTextSize(legTextSize);
      lgd1->AddEntry(g2[fs4e   ][0],(" "+fsLabel[fs4e   ]).c_str(),withMassErrors?"lp":"p");
      lgd1->AddEntry(g2[fs4mu  ][0],(" "+fsLabel[fs4mu  ]).c_str(),withMassErrors?"lp":"p");
      lgd1->AddEntry(g2[fs2e2mu][0],(" "+fsLabel[fs2e2mu]).c_str(),withMassErrors?"lp":"p");
      if(!MERGE2E2MU) lgd1->AddEntry(g2[fs2mu2e][0],(" "+fsLabel[fs2mu2e]).c_str(),withMassErrors?"lp":"p");
      TGraph* gLegCat[nCat];
      for(int cat=0; cat<nCat; cat++){ 
	gLegCat[cat] = (TGraph*)g2[fs4e][cat]->Clone();
	gLegCat[cat]->SetMarkerColor(kBlack);
	gLegCat[cat]->SetLineColor(kBlack);
      }
      lgd2->AddEntry(gLegCat[0],categoryLegLabel[0].c_str(),withMassErrors?"lp":"p");
      lgd2->AddEntry(gLegCat[1],categoryLegLabel[1].c_str(),withMassErrors?"lp":"p");
      lgd2->AddEntry(gLegCat[2],categoryLegLabel[2].c_str(),withMassErrors?"lp":"p");
      lgd3->AddEntry(gLegCat[4],categoryLegLabel[4].c_str(),withMassErrors?"lp":"p");
      lgd3->AddEntry(gLegCat[3],categoryLegLabel[3].c_str(),withMassErrors?"lp":"p");
      lgd3->AddEntry(gLegCat[6],categoryLegLabel[6].c_str(),withMassErrors?"lp":"p");
      lgd3->AddEntry(gLegCat[5],categoryLegLabel[5].c_str(),withMassErrors?"lp":"p");
      lgd1->Draw();
      lgd2->Draw();
      lgd3->Draw();
      for(int cat=0; cat<nCat; cat++){
	if(g2[fs2e2mu][cat]->GetN()>0) g2[fs2e2mu][cat]->Draw(withMassErrors?"P":"XP");
	if(!MERGE2E2MU) if(g2[fs2mu2e][cat]->GetN()>0) g2[fs2mu2e][cat]->Draw(withMassErrors?"P":"XP");
	if(g2[fs4mu][cat]->GetN()>0) g2[fs4mu][cat]->Draw(withMassErrors?"P":"XP");
	if(g2[fs4e][cat]->GetN()>0) g2[fs4e][cat]->Draw(withMassErrors?"P":"XP");
      }
      //h2Stacked->GetYaxis()->SetLimits(0.,1.199);
      h2Stacked->GetYaxis()->SetTitle(Form("%s          ",h2Stacked->GetYaxis()->GetTitle()));
      h2Stacked->GetYaxis()->SetTitleOffset(1.);//1.2);
    }else{
      if(g2[fs2e2mu][nCat]->GetN()>0) g2[fs2e2mu][nCat]->Draw(withMassErrors?"P":"XP");
      if(!MERGE2E2MU) if(g2[fs2mu2e][nCat]->GetN()>0) g2[fs2mu2e][nCat]->Draw(withMassErrors?"P":"XP");
      if(g2[fs4mu][nCat]->GetN()>0) g2[fs4mu][nCat]->Draw(withMassErrors?"P":"XP");
      if(g2[fs4e][nCat]->GetN()>0) g2[fs4e][nCat]->Draw(withMassErrors?"P":"XP");
    }
  }
  if(doBlindingArea) box->Draw();
  if(doBlindingArea2) box2->Draw();
  if(!categMode && bl!=fullyblind && myV2[v2].Name!="M4lVsM4lRefit") lgd->Draw();
  if(useLabel) pav->Draw();
  if(!categMode) gPad->RedrawAxis();

  //----- adjust color axis
  c->Update();
  TPaletteAxis* pal = (TPaletteAxis*)h2Stacked->GetListOfFunctions()->FindObject("palette");
  pal->SetX1NDC(0.875);
  pal->SetX2NDC(0.90);
  if(legendOutMode){
    pal->SetY1NDC(0.13);
    pal->SetY2NDC(0.8);
  }

  //----- customize m4l axis labels
  if(CUSTOMXAXISLABELS && logX && myV2[v2].Name=="M4lVsKD"){
    c->cd();
    TText t;
    t.SetNDC();
    t.SetTextSize(0.04);
    t.SetTextFont(42);
    t.DrawText(0.122,0.093,"110");
    t.DrawText(0.760,0.093,"700");
    t.DrawText(0.830,0.093,"900");
  }

  //----- print official CMS labels and lumi
  writeExtraText = true;
  extraText  = "Preliminary";
  lumi_sqrtS = lumiText + " (13 TeV)";
  CMS_lumi( c, 0, 0 );

  //----- mask last y axis label (should find a less bad way of doing this without ruining the whole layout)
  if(categMode){
    gPad->Update();
    TPad *p = new TPad("p","p",0.,0.,1.,1.); p->SetFillStyle(0); p->Draw(); p->cd();
    TBox* maskBox = new TBox(0.,0.9,0.115,1.);
    maskBox->SetFillColor(kWhite);
    maskBox->Draw();
  }

}


void doPlots(string outputDirectory, double lumi, string lumiText, Int_t FINALSTATE, Int_t CATEGORY)
{

  setTDRStyle();
  gStyle->SetGridColor(kGray);

  //---------- retrieve histograms from the main ROOT file
  TH1F* h1[nVariables][nFinalStates+1][nCat+1][nProcesses];
  TH2F* h2[nVarPairs ][nProcesses];
  TGraphErrors* g2Data[nVarPairs][nFinalStates+1][nCat+1];
  string inFileName = string(Form("histos_plotDataVsMC_%.3ffb-1%s.root",lumi,(MERGE2E2MU?"_m":"")));
  cout<<"Retrieving MC histograms and data graphs from file "<<inFileName<<" ..."<<endl;
  TFile* fInHistos = TFile::Open(inFileName.c_str());
  for(int fs=0; fs<nFinalStates+1; fs++){
    if(MERGE2E2MU && fs==fs2mu2e) continue;
    for(int cat=0; cat<nCat+1; cat++){
      for(int pr=0; pr<nProcesses; pr++){
	for(int v1=0; v1<nVariables; v1++){
	  if(myV1[v1].plotLvl>PLOTLEVEL) continue;
	  if(!( (fs==FINALSTATE || (FINALSTATE==nFinalStates && myV1[v1].InFS && cat==CATEGORY)) &&
		(cat==myV1[v1].Categ || (cat==CATEGORY && CATEGORY!=nCat && myV1[v1].Categ==nCat && fs==FINALSTATE)) )) continue;		
	  h1[v1][fs][cat][pr] = (TH1F*)fInHistos->Get( Form("h1_%s_%s_%s_%s_%s_%s",myV1[v1].Name.c_str(),sBlinding[BLINDING].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str(),sResonantStatus[nRS].c_str(),sProcess[pr].c_str()) );
	  h1[v1][fs][cat][pr]->GetXaxis()->SetTitle(ReplaceString(h1[v1][fs][cat][pr]->GetXaxis()->GetTitle(),"4#font[12]{l}",fsLabel[fs]).c_str());
	  if(myV1[v1].restrXmax) restrictXAxis(h1[v1][fs][cat][pr],myV1[v1].restrXmax);
	}
	for(int v2=0; v2<nVarPairs; v2++){
	  if(myV2[v2].plotLvl>PLOTLEVEL) continue;
	  if(!( fs==nFinalStates && cat==nCat )) continue;	      
	  h2[v2][pr] = (TH2F*)fInHistos->Get( Form("h2_%s_%s_%s_%s_%s_%s",myV2[v2].Name.c_str(),sBlinding[BLINDING].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str(),sResonantStatus[nRS].c_str(),sProcess[pr].c_str()) );
	  h2[v2][pr]->GetXaxis()->SetTitle(ReplaceString(h2[v2][pr]->GetXaxis()->GetTitle(),"4#font[12]{l}",fsLabel[fs]).c_str());
	  h2[v2][pr]->GetYaxis()->SetTitle(ReplaceString(h2[v2][pr]->GetYaxis()->GetTitle(),"4#font[12]{l}",fsLabel[fs]).c_str());
	}
      }
      for(int v2=0; v2<nVarPairs; v2++){
	if(MERGE2E2MU && fs==fs2mu2e) continue;
	g2Data[v2][fs][cat] = (TGraphErrors*)fInHistos->Get( Form("g2Data_%s_%s_%s_%s_%s",myV2[v2].Name.c_str(),sBlinding[BLINDING].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str(),sResonantStatus[nRS].c_str()) );
      }
    }
  }

  //---------- get Z+X histograms
  TH1F* h1_ZPlusX[nVariables][nFinalStates+1][nCat+1];
  if(BUILDZPLUSXHISTOFROMSSCR){
    string inFileName_ZPlusXSS = string(Form("histos_plotDataVsMC_ZPlusXSS_%.3ffb-1%s.root",lumi,(MERGE2E2MU?"_m":"")));
    cout<<"Retrieving Z+X histograms from file "<<inFileName_ZPlusXSS<<" ..."<<endl;
    TFile*fInHistos_ZPlusXSS = TFile::Open(inFileName_ZPlusXSS.c_str());
    for(int fs=0; fs<nFinalStates+1; fs++){
      if(MERGE2E2MU && fs==fs2mu2e) continue;
      for(int cat=0; cat<nCat+1; cat++){
	for(int v1=0; v1<nVariables; v1++){
	  if(myV1[v1].plotLvl>PLOTLEVEL) continue;
	  if(!( (fs==FINALSTATE || (FINALSTATE==nFinalStates && myV1[v1].InFS && cat==CATEGORY)) &&
		(cat==myV1[v1].Categ || (cat==CATEGORY && CATEGORY!=nCat && myV1[v1].Categ==nCat && fs==FINALSTATE)) )) continue;		
	  h1_ZPlusX[v1][fs][cat] = (TH1F*)fInHistos_ZPlusXSS->Get( Form("h1_ZPlusXSS_%s_%s_%s_%s",myV1[v1].Name.c_str(),sBlinding[BLINDING].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str()) );
	  if(myV1[v1].restrXmax) restrictXAxis(h1_ZPlusX[v1][fs][cat],myV1[v1].restrXmax);
	}
      }
    }
  }
  if(USEZPLUSXANALYTICALSHAPE){
    if(!MERGE2E2MU){
      cout<<"ERROR: cannot get Z+X histograms from combined shape for 2e2mu and 2mu2e separately"<<endl;
      abort();
    }else{ // override the histograms from SS method if any 
      for(int v1=0; v1<nVariables; v1++){
	if(myV1[v1].plotLvl>PLOTLEVEL) continue;
	if(myV1[v1].Name.find("M4l")==0 && myV1[v1].Name.find("refit")==string::npos && myV1[v1].CutLabel==""){
	  getM4lZPlusXHistoFromAnalyticalShape_InCateg(h1_ZPlusX[v1][fs4mu],h1_ZPlusX[v1][fs4e],h1_ZPlusX[v1][fs2e2mu],h1_ZPlusX[v1][nFinalStates],myV1[v1].Nbin,myV1[v1].Min,myV1[v1].Max,v1); 
	}
      }
    }
  }

  //---------- do the plots
  string canvasName;
  TCanvas* c1[nVariables];
  if(DO1DPLOTS){
    cout<<"Doing 1D plots ..."<<endl;
    for(int fs=0; fs<nFinalStates+1; fs++){
      if(MERGE2E2MU && fs==fs2mu2e) continue;
      for(int cat=0; cat<nCat+1; cat++){
	for(int v1=0; v1<nVariables; v1++){
	  if(myV1[v1].plotLvl>PLOTLEVEL) continue;
	  if(myV1[v1].Name.find("M4L118130")!=string::npos && BLINDING!=unblinded && BLINDING!=fullyblind) continue;
	  if(!( (fs==FINALSTATE || (FINALSTATE==nFinalStates && myV1[v1].InFS && cat==CATEGORY)) &&
		(cat==myV1[v1].Categ || (cat==CATEGORY && CATEGORY!=nCat && myV1[v1].Categ==nCat && fs==FINALSTATE)) )) continue;
	  //gStyle->SetFrameLineWidth(myV1[v1].Large?2:1);
	  canvasName = string(Form("%s%s_%s_%s_%s",fs!=nFinalStates?"f":cat!=nCat?"g":myV1[v1].prefix.c_str(),(BLINDING==unblinded?"":"_"+sBlinding[BLINDING]).c_str(),myV1[v1].Name.c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str()));
	  c1[v1] = new TCanvas(canvasName.c_str(),canvasName.c_str(),myV1[v1].Large?650:500,500);
	  DrawDataMC1D(c1[v1],h1[v1][fs][cat],h1_ZPlusX[v1][fs][cat],v1,BLINDING,lumi,lumiText,cat,myV1[v1].isLogx,myV1[v1].isLogy);
	  SaveCanvas(outputDirectory,c1[v1]);
	}
      }
    }
  }
  TCanvas* c2[nVarPairs];
  //TCanvas* c2[nVarPairs][NSTYLES2DPLOT]; //for style tests
  if(DO2DPLOTS && FINALSTATE==nFinalStates && CATEGORY==nCat){
    cout<<"Doing 2D plots ..."<<endl;
    for(int v2=0; v2<nVarPairs; v2++){
      if(myV2[v2].plotLvl>PLOTLEVEL) continue;
      if(myV2[v2].Name.find("M4L118130")!=string::npos && BLINDING!=unblinded && BLINDING!=fullyblind) continue;
      canvasName = string(Form("%s%s_2D_%s_%s",myV2[v2].prefix.c_str(),(BLINDING==unblinded?"":"_"+sBlinding[BLINDING]).c_str(),myV2[v2].Name.c_str(),sFinalState[FINALSTATE].c_str()));
      c2[v2] = new TCanvas(canvasName.c_str(),canvasName.c_str(),500,500);
      DrawDataMC2D(c2[v2],h2[v2],g2Data[v2],v2,BLINDING,lumiText,STYLE2DPLOT,myV2[v2].isLogx,myV2[v2].isLogy);
      SaveCanvas(outputDirectory,c2[v2]);
      /* //for style tests
      for(int st=0; st<NSTYLES2DPLOT; st++){
	canvasName = string(Form("z%s_2D_%s_%s_st%i",(BLINDING==unblinded?"":"_"+sBlinding[BLINDING]).c_str(),myV2[v2].Name.c_str(),sFinalState[FINALSTATE].c_str(),st));
	c2[v2][st] = new TCanvas(canvasName.c_str(),canvasName.c_str(),500,500);
	DrawDataMC2D(c2[v2][st],h2[v2],g2Data[v2],v2,BLINDING,lumiText,st,myV2[v2].isLogx,myV2[v2].isLogy);
	SaveCanvas(outputDirectory,c2[v2][st]);
      }
      //*/
    }
  }


}





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void plotDataVsMC(bool redoHistograms = true, Int_t FINALSTATE = nFinalStates, Int_t CATEGORY = nCat) {


  // --------------- definitions ---------------

  // Specify input/output location

  string inputPathMC   = "";
  string inputPathData = ""; // NB. Having at least the data tree copied locally really speeds things up.

  string inputPathDataForCR = inputPathData;
  string inputFileFakeRates = "../../data/FakeRates/FakeRate_SS_Moriond368.root";

  string outputPath = "";

  // Define the luminosity
  float lumi = ;
  string lumiText = "";
  // ***  full 2016 dataset  ***
  // float lumi = 35.86706;
  // string lumiText = "35.9 fb^{-1}";


  // --------------- processing ---------------

  gSystem->Exec(("mkdir -p "+outputPath).c_str());

  // Prepare the histograms and store them in a ROOT file 
  // (to be done only once, it can take a few minutes)
  if(redoHistograms)
    doHistograms(inputPathMC, inputPathData, lumi);

  // Prepare Z+X histograms
  if(BUILDZPLUSXHISTOFROMSSCR)
    doHistogramsZPlusXSS(inputPathDataForCR, inputFileFakeRates, lumi);

  // Do the plots (+- instantaneous)
  doPlots(outputPath, lumi, lumiText, FINALSTATE, CATEGORY);

}

