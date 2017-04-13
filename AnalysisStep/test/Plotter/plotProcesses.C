///
/// usage:
/// -specify parameters (input/output location, luminosity, m4l window) at the end of this file
/// -run with:
///   root -q -l -b plotProcesses.C++
///

#include <iostream>
#include <fstream>
#include <utility>
#include <map>
#include <vector>
#include <cmath>

#include "TCanvas.h"
#include "TColor.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TPaveText.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TSpline.h"
#include "TString.h"
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

using namespace std;

/// version of the categorization
#define BASKETLIST 12
   //  9 : ICHEP16 subcategories
   // 10 : ICHEP16 categories
   // 11 : Moriond17 subcategories
   // 12 : Moriond17 categories

/// choose whether to use q/g tagging on top of MELA 
#define USEQGTAGGING 0

/// working points for MELA-only categorization
#define WPVBF2JMELA  0.437  // ggH efficiency 0.25
#define USEMASSDEPWPVBF2JMELA 1
Float_t mdWPVBF2JMELA(Float_t m4l) { return 1.043-460./(m4l+634.); }
#define WPVBF1JMELA  0.697  // VBF efficiency 0.8
#define WPWHHADRMELA 0.951  // WH efficiency  0.495
#define WPZHHADRMELA 0.9937 // ZH efficiency  0.495

/// working points for MELA+q/g categorization
#define WPVBF2JMELAQG  0.363  // ggH efficiency 0.25
#define USEMASSDEPWPVBF2JMELAQG 0
Float_t mdWPVBF2JMELAQG(Float_t m4l) { return 1.04-239./(m4l+502.); } // not updated for Mor17
#define WPVBF1JMELAQG  0.716  // VBF efficiency 0.8
#define WPWHHADRMELAQG 0.965  // WH efficiency  0.495
#define WPZHHADRMELAQG 0.9952 // ZH efficiency  0.495

#define MASSZ 91.1876
#define CSVv2M 0.8484
#define CSVv2L 0.5426

#define DEBUG 0
#define PROCESSQQZZ 0
#define PROCESSGGZZ 0
#define USEBBH 0
#define CHECKFROMCATEGORYCC 1 // check that all events end up in the same category as provided by Category.cc

#define FINALSTATE 3 // 0:4mu, 1:4e, 2:2e2mu, 3:4mu+4e+2e2mu, 4: none('emptyEvents')
#define BTAGGINGSF 1 // 0:none, 1:central, 2:up, 3:down
#define APPLYKFACTORS 1
#define USEPUWEIGHT 0
#define OFFICIALQGTAGGER 1

#define VARLIST 3
#define VARPAIRLIST 2

/// plotting-related flags
#define SKIPPROCESSES 1
#define NORMPERSLICE 0 // 0:no 1:X 2:Y
#define PRINTYIELDS 0 // print yield numbers on basket efficiency plots
#define DRAWLINEALL 0 // draw inclusive line at the top of the category composition plot
#define DRAWWPMARKER 1
#define PRINTPLOTINFO 0
string INFO = "13TeV samples";
string INFO2 = "#sqrt{s} = 13 TeV";

/// what to do with non-resonant signal events (ZH,H->2l2X and ttH,H->2l2X):
#define excludeH2l2X     1 // completely exclude these events from the study
#define treatH2l2XAsBkgd 1 // treat these events as background in categorization-related plots 

/// choose which plots are made
#define doProdComp       1 // compare distribution of variables across processes (the variable list is controlled by VARLIST)
#define doProdCompMatch4 0
#define doMatch4OrNot    0
#define doMatchHLeps     0
#define doMatchAllLeps   0
#define doMatchWHZHttH   0 
#define do2DPlots        0 
#define doROCs           1 // ROCs comparing discriminants
#define doROCzooms       0
#define doEVWs           0
#define doBaskets        1 // categorization plots (efficiency, purity, etc.)
#define doYieldStudy     0
#define doPlotsByGenChan 0

/// a few flags that can restrict the study to a subset of the events (not very relevant anymore)
#define requireExactly4GoodLeps   0 
#define requireAtLeast5GoodLeps   0 
#define requireExactly5GoodLeps   0 
#define requireExactly6GoodLeps   0 
#define requireHLepsAreInEtaPtAcc 0 // impose that all 4 gen leptons from the Higgs are in the acceptance
#define requireHLepsAreGood       0 // impose that all 4 gen leptons from the Higgs are matched to good leptons (from the candidate or extra leptons)





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////// global variables /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


#define nGenChannels 11
#define nRecoChannels 7

#define nMatchHLepsStatuses 6
#define nMatchAllLepsStatuses 5
#define nMatchWHStatuses 5
#define nMatchZHStatuses 7
#define nMatchttHStatuses 9

#define nAssocWDecays 2
#define nAssocZDecays 3
#define nAssocttDecays 4
#define nAssocDecays (nAssocWDecays + nAssocZDecays + nAssocttDecays)



enum Process {ggH=0, qqH=1, WH=2, ZH=3, ttH=4, bbH=5, qqZZ=6, ggZZ=7};
const int nProcesses = 8;
const int nSignalProcesses = 6;
//string sProcess[nProcesses] = {"ggH", "qqH", "WH", "ZH", "ttH", "bbH", "qqZZ", "ggZZ"};
string sProcess[nProcesses] = {"ggH", "VBF", "WH", "ZH", "ttH", "bbH", "qqZZ", "ggZZ"};
string processLabel[nProcesses] = {"ggH", "VBF", "WH", "ZH", "t#bar{t}H", "b#bar{b}H", "q#bar{q}#rightarrowZZ", "gg#rightarrowZZ"};
bool isSignal[nProcesses] = {1,1,1,1,1,1,0,0,};
bool isProcessed[nProcesses] = {1,1,1,1,1,USEBBH,PROCESSQQZZ,PROCESSGGZZ,};
Int_t processColor[nProcesses] = { kBlue, kGreen+2, kRed, kOrange+1, kMagenta-7, kYellow+1, kBlack, kGray };



enum Variable {
  M4l,
  M4l2,
  MZ1,
  MZ2,
  Dkinbkg,
  DjetFisher,
  Pvbf,
  Phjj,
  Pvbf1j,
  Phj,
  Pwhhadr,
  Pzhhadr,
  Pwhlept,
  Pzhlept,
  D2jVbfHjj,
  D1jVbfHj,
  D2jWHHadrHjj,
  D2jZHHadrHjj,
  Pqj1,
  Pgj1,
  Pqj1Pqj2,
  Pgj1Pgj2,
  D2jqg,
  Dqgj1Dqgj2,
  Pqj1VbfTopo,
  Pgj1VbfTopo,
  Pqj1Pqj2VbfTopo,
  Pgj1Pgj2VbfTopo,
  D2jqgVbfTopo,
  Pq,
  Pg,
  D1jqg,
  D2jMelaQGVbfHjj,
  D2jMelaD2jQGVbfHjj,
  D1jMelaQGVbfHj,
  D1jMelaD1jQGVbfHj,
  D2jMelaQGWHHadrHjj,
  D2jMelaD2jQGWHHadrHjj,
  D2jMelaQGZHHadrHjj,
  D2jMelaD2jQGZHHadrHjj,
  RatioPvbfPhjj,
  RatioPqj1Pqj2Pgj1Pgj2,
  RatioPvbfPqj1Pqj2PhjjPgj1Pgj2,
  D2jMelaExpQGVbfHjj,
  D2jMelaSqQGVbfHjj,
  D2jMelaSqrtQGVbfHjj,
  D2jMelaCbrtQGVbfHjj,
  D2jMelaQrrtQGVbfHjj,
  D2jMelaQnrtQGVbfHjj,
  D1jMelaSqrtQGVbfHj,
  D1jMelaCbrtQGVbfHj,
  D1jMelaQrrtQGVbfHj,
  D1jMelaQnrtQGVbfHj,
  D2jMelaSqrtQGWHHadrHjj,
  D2jMelaCbrtQGWHHadrHjj,
  D2jMelaQrrtQGWHHadrHjj,
  D2jMelaQnrtQGWHHadrHjj,
  D2jMelaSqrtQGZHHadrHjj,
  D2jMelaCbrtQGZHHadrHjj,
  D2jMelaQrrtQGZHHadrHjj,
  D2jMelaQnrtQGZHHadrHjj,
  Pt4l,
  NGenLep,
  NGenLepInEtaPtAcc,
  NGenLepNotInEtaPtAcc,
  NGenHLepNotInEtaPtAcc,
  NGenAssocLepNotInEtaPtAcc,
  NGenLepMinusNGoodLep,
  NGenLepInEtaPtAccMinusNGoodLep,
  NExtraLep,
  NExtraZ,
  NJets,
  NBtaggedJets,
  MET,
};
const int nVariables = 74;
string varName[nVariables] = {
  "M4l",
  "M4l2",
  "MZ1",
  "MZ2",
  "Dkinbkg",
  "DjetFisher",
  "Pvbf",
  "Phjj",
  "Pvbf1j",
  "Phj",
  "Pwhhadr",
  "Pzhhadr",
  "Pwhlept",
  "Pzhlept",
  "D2jVbfHjj",
  "D1jVbfHj",
  "D2jWHHadrHjj",
  "D2jZHHadrHjj",
  "Pqj1",
  "Pgj1",
  "Pqj1Pqj2",
  "Pgj1Pgj2",
  "D2jqg",
  "Dqgj1Dqgj2",
  "Pqj1VbfTopo",
  "Pgj1VbfTopo",
  "Pqj1Pqj2VbfTopo",
  "Pgj1Pgj2VbfTopo",
  "D2jqgVbfTopo",
  "Pq",
  "Pg",
  "D1jqg",
  "D2jMelaQGVbfHjj",
  "D2jMelaD2jQGVbfHjj",
  "D1jMelaQGVbfHj",
  "D1jMelaD1jQGVbfHj",
  "D2jMelaQGWHHadrHjj",
  "D2jMelaD2jQGWHHadrHjj",
  "D2jMelaQGZHHadrHjj",
  "D2jMelaD2jQGZHHadrHjj",
  "RatioPvbfPhjj",
  "RatioPqj1Pqj2Pgj1Pgj2",
  "RatioPvbfPqj1Pqj2PhjjPgj1Pgj2",
  "D2jMelaExpQGVbfHjj",
  "D2jMelaSqQGVbfHjj",
  "D2jMelaSqrtQGVbfHjj",
  "D2jMelaCbrtQGVbfHjj",
  "D2jMelaQrrtQGVbfHjj",
  "D2jMelaQnrtQGVbfHjj",
  "D1jMelaSqrtQGVbfHj",
  "D1jMelaCbrtQGVbfHj",
  "D1jMelaQrrtQGVbfHj",
  "D1jMelaQnrtQGVbfHj",
  "D2jMelaSqrtQGWHHadrHjj",
  "D2jMelaCbrtQGWHHadrHjj",
  "D2jMelaQrrtQGWHHadrHjj",
  "D2jMelaQnrtQGWHHadrHjj",
  "D2jMelaSqrtQGZHHadrHjj",
  "D2jMelaCbrtQGZHHadrHjj",
  "D2jMelaQrrtQGZHHadrHjj",
  "D2jMelaQnrtQGZHHadrHjj",
  "Pt4l",
  "NGenLep",
  "NGenLepInEtaPtAcc",
  "NGenLepNotInEtaPtAcc",
  "NGenHLepNotInEtaPtAcc",
  "NGenAssocLepNotInEtaPtAcc",
  "NGenLepMinusNGoodLep",
  "NGenLepInEtaPtAccMinusNGoodLep",
  "NExtraLep",
  "NExtraZ",
  "NJets",
  "NBtaggedJets",
  "MET",
};
string varLabel[nVariables] = {
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{Z_{1}} (GeV)",
  "m_{Z_{2}} (GeV)",
  "D_{bkg}^{kin}",
  "D^{Fisher}",
  "P_{VBF}^{MELA}",
  "P_{ggH+2j}^{MELA}",
  "P_{VBF 1j}^{MELA}",
  "P_{ggH+1j}^{MELA}",
  "P_{WH-h}^{MELA}",
  "P_{ZH-h}^{MELA}",
  "P_{WH-l}^{MELA}",
  "P_{ZH-l}^{MELA}",
  "D_{VBF-2j}^{ME}",
  "D_{VBF-1j}^{ME}",
  "D_{WH-hadr.}^{ME}",
  "D_{ZH-hadr.}^{ME}",
  "P_{q}(j_{1})",
  "P_{g}(j_{1})",
  "P_{q}(j_{1})*P_{q}(j_{2})",
  "P_{g}(j_{1})*P_{g}(j_{2})",
  "D_{2jets}^{q/g}",
  "D_{1jet}^{q/g}(j_{1})*D_{1jet}^{q/g}(j_{2})",
  "P_{q}(j_{1})",
  "P_{g}(j_{1})",
  "P_{q}(j_{1})*P_{q}(j_{2})",
  "P_{g}(j_{1})*P_{g}(j_{2})",
  "D_{2jets}^{q/g}",
  "P_{q}",
  "P_{g}",
  "D_{1jet}^{q/g}",
  "D_{VBF-2j}^{ME, q/g}",
  "D_{VBF-2j}^{ME}*D_{2jets}^{q/g}",
  "D_{VBF-1j}^{ME, q/g}",
  "D_{VBF-1j}^{ME}*D_{1jet}^{q/g}",
  "D_{WH-hadr.}^{ME, q/g}",
  "D_{WH-hadr.}^{ME}*D_{2jets}^{q/g}",
  "D_{ZH-hadr.}^{ME, q/g}",
  "D_{ZH-hadr.}^{ME}*D_{2jets}^{q/g}",
  "P_{VBF}^{MELA} / P_{ggH+2j}^{MELA}",
  "[P_{q}(j_{1})*P_{q}(j_{2})] / [P_{g}(j_{1})*P_{g}(j_{2})]",
  "[P_{VBF}^{MELA}*P_{q}(j_{1})*P_{q}(j_{2})] / [P_{ggH+2j}^{MELA}*P_{g}(j_{1})*P_{g}(j_{2})]",
  "D_{VBF-2j}^{ME, exp(q/g)}",
  "D_{VBF-2j}^{ME, (q/g)^2}",
  "D_{VBF-2j}^{ME, (q/g)^(1/2)}",
  "D_{VBF-2j}^{ME, (q/g)^(1/3)}",
  "D_{VBF-2j}^{ME, (q/g)^(1/4)}",
  "D_{VBF-2j}^{ME, (q/g)^(1/5)}",
  "D_{VBF-1j}^{ME, (q/g)^(1/2)}",
  "D_{VBF-1j}^{ME, (q/g)^(1/3)}",
  "D_{VBF-1j}^{ME, (q/g)^(1/4)}",
  "D_{VBF-1j}^{ME, (q/g)^(1/5)}",
  "D_{WH-hadr.}^{ME, (q/g)^(1/2)}",
  "D_{WH-hadr.}^{ME, (q/g)^(1/3)}",
  "D_{WH-hadr.}^{ME, (q/g)^(1/4)}",
  "D_{WH-hadr.}^{ME, (q/g)^(1/5)}",
  "D_{ZH-hadr.}^{ME, (q/g)^(1/2)}",
  "D_{ZH-hadr.}^{ME, (q/g)^(1/3)}",
  "D_{ZH-hadr.}^{ME, (q/g)^(1/4)}",
  "D_{ZH-hadr.}^{ME, (q/g)^(1/5)}",
  "p_{T}^{4l} (GeV)",
  "# gen leptons",
  "# gen leptons in acceptance",
  "# gen leptons not in acceptance",
  "# gen leptons from H not in acceptance",
  "# gen associated leptons not in acceptance",
  "# gen leptons - # good leptons",
  "# gen leptons in acceptance - # good leptons",
  "number of additional leptons",//"# extra leptons",
  "number of additional #font[12]{l}^{+}#font[12]{l}^{-} pairs",//"# extra #font[12]{l}^{+}#font[12]{l}^{-} pairs",
  "number of selected jets",//"# jets",
  "number of selected b-tagged jets",//"# b-tagged jets",
  "E_{T}^{miss.}",
};
string varCutLabel[nVariables] = {
  "","","","","","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} = 1","N_{jets} = 1","N_{jets} #geq 2","N_{jets} #geq 2","N_{leptons} #geq 5","N_{extra Zs} #geq 1","N_{jets} #geq 2","N_{jets} = 1","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets}#geq2, D_{VBF/ggH+2j}^{MELA}>0.5","N_{jets}#geq2, D_{VBF/ggH+2j}^{MELA}>0.5","N_{jets}#geq2, D_{VBF/ggH+2j}^{MELA}>0.5","N_{jets}#geq2, D_{VBF/ggH+2j}^{MELA}>0.5","N_{jets}#geq2, D_{VBF/ggH+2j}^{MELA}>0.5","N_{jets} = 1","N_{jets} = 1","N_{jets} = 1","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} = 1","N_{jets} = 1","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} = 1","N_{jets} = 1","N_{jets} = 1","N_{jets} = 1","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","","","","","","","","","","","","","",
};
Bool_t plotThisVar[4][nVariables] = {
  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,},
  {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,},
  {1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,},
  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,},
};
Int_t  varNbin[nVariables] =    { 100,  70,  75,  75,  20,  50,    50,   50,   50,   50,     50,     50,   50,   50,   51,   51,   51,   51,   25,   25,   25,   25,   51,   51,   25,   25,   25,   25,   26,   25,   25,   51,   51,   51,   51,   51,   51,   51,   51,   51, 50, 50, 50,   51,   51,   51,   51,   51,   51,   51,   51,   51,   51,   51,   51,   51,   51,   51,   51,   51,   51,  50,  7,  7,  7,  5,  3,  6,  6,  6,  4, 15, 6, 100, };
Float_t varMin[nVariables] =    {  50, 105,   0,   0,   0,   0,    -2,   -5,  -10,  -10,    -10,    -10,   -1,   -1,    0,    0,    0,    0,-0.04,-0.04,-0.04,-0.04,    0,    0,-0.04,-0.04,-0.04,-0.04,    0,-0.04,-0.04,    0,    0,    0,    0,    0,    0,    0,    0,    0,  0,  0,  0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,   0,  0,  0,  0,  0,  0, -3, -3,  0,  0,  0, 0,   0, };
Float_t varMax[nVariables] =    { 850, 140, 150, 150,   1,   2,    98,  245,  490,  490,    490,    490,   49,   49, 1.02, 1.02, 1.02, 1.02, 0.96, 0.96, 0.96, 0.96, 1.02, 1.02, 0.96, 0.96, 0.96, 0.96, 1.04, 0.96, 0.96, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02,  5,  5,  5, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 500,  7,  7,  7,  5,  3,  3,  3,  6,  4, 15, 6, 200, };
Float_t varMaxRoc[nVariables] = { 850, 140, 150, 150,   1,  20, 10000,10000,10000,10000,1000000,1000000, 1000, 1000, 1.02, 1.02, 1.02, 1.02,  100,  100,  100,  100, 1.02, 1.02,  100,  100,  100,  100, 1.04,  100,  100, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02,100,100,100, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 500,  7,  7,  7,  5,  3,  3,  3,  6,  4, 15, 6, 200, };
Bool_t varLogy[nVariables] = {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,};
Float_t varMinFactor[nVariables] = {1000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,600,600,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000,10000,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,50000000,0,0,0,};

const Int_t nVarPairs = 3;
Int_t varPairRef[nVarPairs][3] = {
  { M4l, Dkinbkg, 1 },
  { MZ2, Dkinbkg, 1 },
  { D2jVbfHjj, D2jqg, 0 },
};
string varPairCutLabel[nVarPairs] = {
  "","","N_{jets} #geq 2",
};
Bool_t plotThisVarPair[3][nVarPairs] = {
  {0,0,0,},
  {1,1,1,},
  {0,0,1,},
};

const Int_t nROCs = 78;
Int_t rocRef[nROCs][4] = {
  { DjetFisher, qqH, ggH, 1 },
  { DjetFisher, qqH, qqZZ, 1 },
  { Pvbf, qqH, ggH, 1 },
  { Pvbf, qqH, qqZZ, 1 },
  { Phjj, qqH, ggH, 0 },
  { Phjj, qqH, qqZZ, 0 },
  { Phjj, WH, ggH, 0 },
  { Phjj, ZH, ggH, 0 },
  { Pvbf1j, qqH, ggH, 1 },
  { Phj, qqH, ggH, 0 },
  { Pwhhadr, WH, ggH, 1 },
  { Pzhhadr, ZH, ggH, 1 },
  { D2jVbfHjj, qqH, ggH, 1 },
  { D2jVbfHjj, qqH, qqZZ, 1 },
  { D1jVbfHj, qqH, ggH, 1 },
  { D2jWHHadrHjj, WH, ggH, 1 },
  { D2jZHHadrHjj, ZH, ggH, 1 },
  { Pqj1Pqj2, qqH, ggH, 1 },
  { Pgj1Pgj2, qqH, ggH, 0 },
  { D2jqg, qqH, ggH, 1 },
  { Dqgj1Dqgj2, qqH, ggH, 1 },
  { Dqgj1Dqgj2, qqH, ggH, 1 },
  { Pqj1Pqj2, qqH, qqZZ, 1 },
  { Pgj1Pgj2, qqH, qqZZ, 0 },
  { D2jqg, qqH, qqZZ, 1 },
  { Dqgj1Dqgj2, qqH, qqZZ, 1 },
  { Dqgj1Dqgj2, qqH, qqZZ, 1 },
  { Pq, qqH, ggH, 1 },
  { Pg, qqH, ggH, 0 },
  { D1jqg, qqH, ggH, 1 },
  { D1jqg, qqH, ggH, 1 },
  { Pqj1Pqj2, WH, ggH, 1 },
  { Pgj1Pgj2, WH, ggH, 0 },
  { D2jqg, WH, ggH, 1 },
  { Dqgj1Dqgj2, WH, ggH, 1 },
  { Dqgj1Dqgj2, WH, ggH, 1 },
  { Pqj1Pqj2, ZH, ggH, 1 },
  { Pgj1Pgj2, ZH, ggH, 0 },
  { D2jqg, ZH, ggH, 1 },
  { Dqgj1Dqgj2, ZH, ggH, 1 },
  { Dqgj1Dqgj2, ZH, ggH, 1 },
  { D2jMelaQGVbfHjj, qqH, ggH, 1 },
  { D2jMelaD2jQGVbfHjj, qqH, ggH, 1 },
  { D2jMelaQGVbfHjj, qqH, qqZZ, 1 },
  { D2jMelaD2jQGVbfHjj, qqH, qqZZ, 1 },
  { D1jMelaQGVbfHj, qqH, ggH, 1 },
  { D1jMelaD1jQGVbfHj, qqH, ggH, 1 },
  { D2jMelaQGWHHadrHjj, WH, ggH, 1 },
  { D2jMelaD2jQGWHHadrHjj, WH, ggH, 1 },
  { D2jMelaQGZHHadrHjj, ZH, ggH, 1 },
  { D2jMelaD2jQGZHHadrHjj, ZH, ggH, 1 },
  { RatioPvbfPhjj, qqH, ggH, 1 },
  { RatioPqj1Pqj2Pgj1Pgj2, qqH, ggH, 1 },
  { RatioPvbfPqj1Pqj2PhjjPgj1Pgj2, qqH, ggH, 1 },
  { D2jMelaExpQGVbfHjj, qqH, ggH, 1 },
  { D2jMelaSqQGVbfHjj, qqH, ggH, 1 },
  { D2jMelaSqrtQGVbfHjj, qqH, ggH, 1 },
  { D2jMelaCbrtQGVbfHjj, qqH, ggH, 1 },
  { D2jMelaQrrtQGVbfHjj, qqH, ggH, 1 },
  { D2jMelaQnrtQGVbfHjj, qqH, ggH, 1 },
  { D2jMelaExpQGVbfHjj, qqH, qqZZ, 1 },
  { D2jMelaSqQGVbfHjj, qqH, qqZZ, 1 },
  { D2jMelaSqrtQGVbfHjj, qqH, qqZZ, 1 },
  { D2jMelaCbrtQGVbfHjj, qqH, qqZZ, 1 },
  { D2jMelaQrrtQGVbfHjj, qqH, qqZZ, 1 },
  { D2jMelaQnrtQGVbfHjj, qqH, qqZZ, 1 },
  { D1jMelaSqrtQGVbfHj, qqH, ggH, 1 },
  { D1jMelaCbrtQGVbfHj, qqH, ggH, 1 },
  { D1jMelaQrrtQGVbfHj, qqH, ggH, 1 },
  { D1jMelaQnrtQGVbfHj, qqH, ggH, 1 },
  { D2jMelaSqrtQGWHHadrHjj, WH, ggH, 1 },
  { D2jMelaCbrtQGWHHadrHjj, WH, ggH, 1 },
  { D2jMelaQrrtQGWHHadrHjj, WH, ggH, 1 },
  { D2jMelaQnrtQGWHHadrHjj, WH, ggH, 1 },
  { D2jMelaSqrtQGZHHadrHjj, ZH, ggH, 1 },
  { D2jMelaCbrtQGZHHadrHjj, ZH, ggH, 1 },
  { D2jMelaQrrtQGZHHadrHjj, ZH, ggH, 1 },
  { D2jMelaQnrtQGZHHadrHjj, ZH, ggH, 1 },
};
//Float_t rocWP[nROCs] = {0,0,0,0,0,0,0,0,0,0,0,0,0.5,0.5,0.5,0.5,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//Float_t rocWP[nROCs] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//Float_t rocWP[nROCs] = {0.5,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,WPWHHADRMELAQG,0,WPZHHADRMELAQG,0,0,0,0,0,0,0,WPVBF2JMELAQG,0,0,0,0,0,WPVBF2JMELAQG,0,0,0,WPVBF1JMELAQG,0,0,0,0,0,0,0,0,0,0};
//Float_t rocWP[nROCs] = {0.5,0.5,0,0,0,0,0,0,0,0,0,0,WPVBF2JMELA,0,WPVBF1JMELA,WPWHHADRMELA,WPZHHADRMELA,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,WPVBF1JMELAQG,0,0,0,0,0,0,0,0,0,0,0,WPVBF2JMELAQG,0,0,0,0,0,WPVBF2JMELAQG,0,0,0,0,0,0,0,WPWHHADRMELAQG,0,0,0,WPZHHADRMELAQG,0,0};
Float_t rocWP[nROCs] = {0.5,0.5,0,0,0,0,0,0,0,0,0,0,WPVBF2JMELA,0,WPVBF1JMELA,WPWHHADRMELA,WPZHHADRMELA,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,WPVBF2JMELAQG,0,0,0,0,0,WPVBF2JMELAQG,0,0,0,WPVBF1JMELAQG,0,0,0,WPWHHADRMELAQG,0,0,0,WPZHHADRMELAQG,0,0};
string rocCutLabel[nROCs] = {
  "N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} = 1","N_{jets} = 1","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} = 1","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} = 1","N_{jets} = 1","N_{jets} = 1","N_{jets} = 1","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} = 1","N_{jets} = 1","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} = 1","N_{jets} = 1","N_{jets} = 1","N_{jets} = 1","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2",
};


const Int_t nEVWs = 12;
Int_t evwRef[nEVWs][3] = {
  { D2jVbfHjj, 1 },
  { D1jVbfHj, 1 },
  { D2jWHHadrHjj, 1 },
  { D2jZHHadrHjj, 1 },
  { D2jMelaCbrtQGVbfHjj, 1 },
  { D1jMelaCbrtQGVbfHj, 1 },
  //{ D1jMelaQGVbfHj, 1 },
  { D2jMelaCbrtQGWHHadrHjj, 1 },
  { D2jMelaCbrtQGZHHadrHjj, 1 },
  { D2jMelaQGWHHadrHjj, 1 },
  { D2jMelaQGZHHadrHjj, 1 },
  { Pwhlept, 1 },
  { Pzhlept, 1 },
};
string evwCutLabel[nEVWs] = {
  "N_{jets} #geq 2","N_{jets} = 1","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} = 1","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{jets} #geq 2","N_{leptons} #geq 5","N_{extra Zs} #geq 1",
};


const int nBasketLists = 13;
int basketListSize[nBasketLists] = { 7, 9, 24, 6, 38, 7, 10, 38, 7, 62, 7, 65, 8};
bool isSubcat[nBasketLists] = {0,0,1,0,1,0,0,1,0,1,0,1,0,};
Float_t longDim[nBasketLists] = {500,500,1000,500,1000,500,500,1000,500,1500,500,1500,500};
const int nBaskets = basketListSize[BASKETLIST];


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////////////////////// MISCELLANEOUS FUNCTIONS ////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


double rtd(double value, int digits)
{
  if (value == 0.0) // otherwise it will return 'nan' due to the log10() of zero
    return 0.0;

  double factor = pow(10.0, digits - ceil(log10(fabs(value))));
  return round(value * factor) / factor;   
}

bool test_bit_16( short mask, unsigned int iBit ) { return (mask >> iBit) & 1; }

Double_t deltaR(Double_t e1, Double_t p1, Double_t e2, Double_t p2) {
  Double_t deltaPhi = acos(cos(p1-p2));
  return TMath::Sqrt((e1-e2)*(e1-e2) + deltaPhi*deltaPhi);
}

string eventID(Int_t nRun, Int_t nEvent, Int_t nLumi) {
  return string(Form("run%i_lumi%i_event%i",nRun,nEvent,nLumi));
}

void printStatus(Int_t n, Int_t interval, Int_t total, string label) {
  if(n==0)              cout<<"  "<<1  <<" / "<<total<<" "<<label<<endl;
  if((n+1)%interval==0) cout<<"  "<<n+1<<" / "<<total<<" "<<label<<endl;
  if(n==total-1)        cout<<"  "<<n+1<<" / "<<total<<" "<<label<<endl;
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

string repeat(string pattern, int n) {
  string res = "";
  for(int i=0; i<n; i++) res += pattern;
  return res;
}

string rounding2(float x) {
  return string(Form("%.2f",x));
}
string rounding3(float x) {
  return string(Form("%.3f",x));
}

string percentage(float frac) {
  return string(Form("%.1f",100*frac));
}

void prepareBasketlabels(TH1F* h, vector<string> basketLabel){
  for(int i=1; i<=h->GetNbinsX(); i++)  h->GetXaxis()->SetBinLabel(i,basketLabel[i-1].c_str());
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetXaxis()->LabelsOption("h");
  h->GetXaxis()->SetNdivisions(-h->GetNbinsX());
}

void prepareBasketlabelsHorizontal(TH2F* h2, vector<string> basketLabel){
  int nbins = h2->GetNbinsY();
  for(int i=1; i<=nbins; i++) h2->GetYaxis()->SetBinLabel(nbins+1-i,basketLabel[i-1].c_str());
  h2->GetYaxis()->SetLabelSize(longDim[BASKETLIST]>500?0.03:0.06);
  h2->GetYaxis()->LabelsOption("h");
  h2->GetYaxis()->SetNdivisions(-h2->GetNbinsY());
}

void TriBulle(float* array, int length, bool ascendingOrder = true){
  int i = length;
  bool swap = true;
  while(i>0 && swap){
    swap = false;
    for(int j=0; j<i-1; j++){
      if( (array[j]>array[j+1] && ascendingOrder) || (array[j]<array[j+1] && !ascendingOrder) ){
	float tmp = array[j];
	array[j] = array[j+1];
	array[j+1] = tmp;
	swap = true;
      }
    }
    i--;
  }
}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////// PLOTTING FUNCTIONS ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void DrawByProdmodes(TCanvas* c, TH1F* hSource[nProcesses][nRecoChannels], int v) {

  gStyle->SetOptTitle(0);

  c->cd();
  Bool_t logY = false;
  if(v>=0) logY = varLogy[v];
  if(logY) c->SetLogy();
  c->SetTicks(0,0);

  TH1F* h[nProcesses];
  Bool_t skipProcess[nProcesses];
  Int_t nHist = 0;
  for(int pr=0; pr<nProcesses; pr++){
    if(!isProcessed[pr]) continue;
    skipProcess[pr] = //false;
      (v==Pvbf   && !(pr==ggH||pr==qqH)) ||
      (v==Phjj   && !(pr==ggH||pr==qqH||pr==WH||pr==ZH)) ||
      (v==Pvbf1j && !(pr==ggH||pr==qqH)) ||
      (v==Phj    && !(pr==ggH||pr==qqH)) ||
      (v==Pwhhadr&& !(pr==ggH||pr==WH )) ||
      (v==Pzhhadr&& !(pr==ggH||pr==ZH )) ||
      ((v==Pqj1||v==Pgj1||v==Pqj1Pqj2||v==Pgj1Pgj2||v==D2jqg||v==Pqj1VbfTopo||v==Pgj1VbfTopo||v==Pqj1Pqj2VbfTopo||v==Pgj1Pgj2VbfTopo||v==D2jqgVbfTopo) && !(pr==ggH||pr==qqH||pr==qqZZ)) ||
      //((v==D2jVbfHjj||v==D2jqg||v==D2jMelaQGVbfHjj||v==RatioPvbfPhjj||v==RatioPqj1Pqj2Pgj1Pgj2||v==RatioPvbfPqj1Pqj2PhjjPgj1Pgj2) && !(pr==ggH||pr==qqH||pr==qqZZ)) ||
      (v==D2jWHHadrHjj&& !(pr==ggH||pr==qqH||pr==WH||pr==ZH )) ||
      (v==D2jZHHadrHjj&& !(pr==ggH||pr==qqH||pr==WH||pr==ZH )) ;
    //*/
    if(SKIPPROCESSES && skipProcess[pr]) continue;
    nHist++;
    h[pr] = (TH1F*)hSource[pr][FINALSTATE]->Clone(); // includes the 3 decay channels !!
  }

  float lgdLow = 0.92-nHist*0.05;
  TLegend* lgd = new TLegend(0.72,lgdLow,0.92,0.92);
  lgd->SetFillStyle(0);
  lgd->SetBorderSize(0);
  lgd->SetTextSize(0.04);

  Float_t max = 0.;
  for(int pr=0; pr<nProcesses; pr++){
    if(!isProcessed[pr]) continue;
    if(SKIPPROCESSES && skipProcess[pr]) continue;

    h[pr]->Scale(1/h[pr]->Integral(0,h[pr]->GetNbinsX()+1));
    Float_t maxtemp = h[pr]->GetMaximum();
    if(maxtemp>max) max = maxtemp;

  }
  Float_t cmax = logY ? 10.*max : (v==Dkinbkg?2.:v==D1jVbfHj?1.3:1.1)*max ;
  Float_t cminlog = max / varMinFactor[v];

  for(int pr=0; pr<nProcesses; pr++){
    if(!isProcessed[pr]) continue;
    if(SKIPPROCESSES && skipProcess[pr]) continue;

    h[pr]->SetLineColor(processColor[pr]);
    h[pr]->SetLineWidth(2);
    h[pr]->SetFillStyle(0);
    h[pr]->SetStats(0);
    //if(!logY) h[pr]->SetMinimum(0);

    if(pr==0){
      h[pr]->SetMaximum(cmax);
      h[pr]->SetMinimum(logY?cminlog:0.);
      h[pr]->GetYaxis()->SetTitle("normalized to 1");
      h[pr]->Draw("hist"); 
    }else{
      h[pr]->Draw("histsame");
    }

    lgd->AddEntry( h[pr], processLabel[pr].c_str(), "l" );
  }

  lgd->Draw();

  if(v>=0 && varCutLabel[v]!="")
    printInfo(varCutLabel[v],0.72,lgdLow-0.09,1.,lgdLow-0.04);

  gPad->RedrawAxis();

  writeExtraText = true;
  extraText  = "Simulation";
  lumi_sqrtS = "13 TeV";
  CMS_lumi( c, 0, 0);

}

void DrawMatch4OrNot(TCanvas* c, TH1F* h, TH1F* hMatch, string title, Color_t color, Color_t colorNoMatch, Bool_t logY = false) {

  gStyle->SetOptTitle(1);

  c->cd();
  if(logY) c->SetLogy();
  c->SetTicks(0,0);

//   Int_t scale = 1/h->Integral();
//   h->Scale(scale);
//   hMatch->Scale(scale);

  h->SetLineColor(colorNoMatch);
  h->SetFillColor(colorNoMatch);
  hMatch->SetLineColor(color);
  hMatch->SetFillColor(color);
  h->SetTitle(title.c_str()); 
  hMatch->SetTitle(title.c_str());
  h->SetStats(0);
  hMatch->SetStats(0);
  //if(!logY) h->SetMinimum(0);
  
  h->Draw("hist"); 
  hMatch->Draw("histsames");
  
  TLegend* lgd = new TLegend(0.78,0.70,0.98,0.82);
  lgd->SetFillColor(0);
  lgd->AddEntry(hMatch,"match 4 gen","f");
  lgd->AddEntry(h,"other cases","f");
  lgd->Draw();

  gPad->RedrawAxis();

}

void DrawMatchHLeps(TCanvas* c, TH1F* h, TH1F** hMatch, string title, Color_t* color, string* matchHLepsKeys, Bool_t logY = false) {

  gStyle->SetOptTitle(1);

  c->cd();
  if(logY) c->SetLogy();
  c->SetTicks(0,0);

  h->SetLineColor(1);
  h->SetFillColor(1);
  h->SetTitle(title.c_str()); 
  h->SetStats(0);
  //if(!logY) h->SetMinimum(0);
  h->Draw("hist"); 

  TH1F* hStacks[nMatchHLepsStatuses];
  for(int m=0; m<nMatchHLepsStatuses; m++){
    hStacks[m] = (TH1F*)hMatch[m]->Clone();
    if(m>0) hStacks[m]->Add(hStacks[m-1]);
    hStacks[m]->SetLineColor(color[m]);
    hStacks[m]->SetFillColor(color[m]);
    hStacks[m]->SetTitle(title.c_str()); 
    hStacks[m]->SetStats(0);
  }
  for(int m=nMatchHLepsStatuses-1; m>=0; m--){
    hStacks[m]->Draw("histsame");
  }
  
  TLegend* lgd = new TLegend(0.78,0.45,0.98,0.82);
  lgd->SetFillColor(0);
  for(int m=0; m<nMatchHLepsStatuses; m++){
    lgd->AddEntry(hStacks[m],matchHLepsKeys[m].c_str(),"f");
  }
  lgd->Draw();

  gPad->RedrawAxis();

}

void DrawMatchAllLeps(TCanvas* c, TH1F* h, TH1F** hMatch, string title, Color_t* color, string* matchAllLepsKeys, Bool_t logY = false) {

  gStyle->SetOptTitle(1);

  c->cd();
  if(logY) c->SetLogy();
  c->SetTicks(0,0);

  h->SetLineColor(1);
  h->SetFillColor(1);
  h->SetTitle(title.c_str()); 
  h->SetStats(0);
  //if(!logY) h->SetMinimum(0);
  h->Draw("hist"); 

  TH1F* hStacks[nMatchAllLepsStatuses];
  for(int m=0; m<nMatchAllLepsStatuses; m++){
    hStacks[m] = (TH1F*)hMatch[m]->Clone();
    if(m>0) hStacks[m]->Add(hStacks[m-1]);
    hStacks[m]->SetLineColor(color[m]);
    hStacks[m]->SetFillColor(color[m]);
    hStacks[m]->SetTitle(title.c_str()); 
    hStacks[m]->SetStats(0);
  }
  for(int m=nMatchAllLepsStatuses-1; m>=0; m--){
    hStacks[m]->Draw("histsame");
  }
  
  TLegend* lgd = new TLegend(0.73,0.50,0.98,0.82);
  lgd->SetFillStyle(0);
  for(int m=0; m<nMatchAllLepsStatuses; m++){
    lgd->AddEntry(hStacks[m],matchAllLepsKeys[m].c_str(),"f");
  }
  lgd->Draw();

  gPad->RedrawAxis();

}

void DrawMatchCustom(TCanvas* c, TH1F* h, TH1F** hMatch, string title, Color_t* color, string* matchKeys, Int_t nStatuses, Float_t legLeft, Float_t legBottom, Bool_t logY = false) {

  gStyle->SetOptTitle(1);

  c->cd();
  if(logY) c->SetLogy();
  c->SetTicks(0,0);

  h->SetLineColor(1);
  h->SetFillColor(1);
  h->SetTitle(title.c_str()); 
  h->SetStats(0);
  //if(!logY) h->SetMinimum(0);
  h->Draw("hist"); 

  Float_t denom = h->Integral();
  string percentages[20];
  for(int m=0; m<nStatuses; m++){
    if(excludeH2l2X && matchKeys[m].find("2l2X")!=string::npos) continue;
    percentages[m] = percentage(hMatch[m]->Integral()/denom);
  }

  TH1F* hStacks[20];
  for(int m=0; m<nStatuses; m++){
    hStacks[m] = (TH1F*)hMatch[m]->Clone();
    if(m>0) hStacks[m]->Add(hStacks[m-1]);
    hStacks[m]->SetLineColor(color[m]);
    hStacks[m]->SetFillColor(color[m]);
    hStacks[m]->SetTitle(title.c_str()); 
    hStacks[m]->SetStats(0);
  }
  for(int m=nStatuses-1; m>=0; m--){
    hStacks[m]->Draw("histsame");
  }
  
  TLegend* lgd = new TLegend(legLeft,legBottom,legLeft+0.55,0.89);
  lgd->SetFillStyle(0);
  lgd->SetBorderSize(0);
  for(int m=0; m<nStatuses; m++){
    if(excludeH2l2X && matchKeys[m].find("2l2X")!=string::npos) continue;
    lgd->AddEntry(hStacks[m],matchKeys[m].c_str(),"f");
  }
  lgd->Draw();
  
  TPaveText* pav = new TPaveText(legLeft+0.03,legBottom,legLeft+0.11,0.89,"brNDC");
  pav->SetFillStyle(0);
  pav->SetBorderSize(0);
  pav->SetTextColor(0);
  for(int m=0; m<nStatuses; m++){
    if(excludeH2l2X && matchKeys[m].find("2l2X")!=string::npos) continue;
    pav->AddText((percentages[m]+" %").c_str());
  }
  pav->Draw();

  gPad->RedrawAxis();

}

void Draw2D(TCanvas* c, TH2F* h, int v2, string title) {
  TStyle* style = (TStyle*)gROOT->GetStyle("tdrStyle")->Clone();
  style->SetPadLeftMargin(0.14);
  style->SetPadRightMargin(0.12);
  style->cd();
  setColZGradient_Rainbow2();
  c->cd();
  c->UseCurrentStyle();
  //c->SetLogz();

  if(NORMPERSLICE){
    bool inX = NORMPERSLICE==1;
    NormalizePerSlice(h,inX);
    c->SetName(Form("%s_NormSlice%s",c->GetName(),inX?"X":"Y"));
  }

  h->GetXaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetZaxis()->SetTitleOffset(1.35);
  h->GetZaxis()->SetLabelSize(0.03);
  h->SetStats(0);
  h->Draw("COLZ");

  c->Update();
  TPaletteAxis* pal = (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette");
  pal->SetX1NDC(0.89);
  pal->SetX2NDC(0.915);

  string label = title;
  if(v2>=0 && varPairCutLabel[v2]!="") label += ", " + varPairCutLabel[v2];
  printInfo(label,0.13,0.95,0.5,0.995);
  if(PRINTPLOTINFO) printInfo(INFO2,0.1,0.9,0.75,0.94);
}

void DrawRoc(TGraph** gRoc, TGraph** gWp, string canvasname, vector<int> rocrefs, Int_t* colors, Int_t* widths, string outDir, string tagOut, bool zoom = false) {
  setTDRStyle();
  gStyle->SetOptTitle(0);
  gStyle->SetGridColor(kGray);
  TCanvas* c1 = new TCanvas(canvasname.c_str(),canvasname.c_str(),500,500);
  c1->cd();
  c1->SetGrid();
  Int_t nrocs = rocrefs.size();
  TLegend* lgd = new TLegend(0.65,0.17,0.93,0.17+nrocs*0.06);
  bool init = false;
  TGraph* myRoc[nrocs];
  TGraph* myWp[nrocs];
  TFile* fOutRocGraphs = new TFile("fROCs.root","update");
  fOutRocGraphs->cd();
  for(int r=0; r<nrocs; r++){
    myRoc[r] = (TGraph*)gRoc[rocrefs[r]]->Clone();
    if(zoom){
      myRoc[r]->GetXaxis()->SetRangeUser(0.,0.5);
      myRoc[r]->GetYaxis()->SetRangeUser(0.5,1.);
    }
    myRoc[r]->SetLineColor(colors[r]);
    myRoc[r]->SetMarkerColor(colors[r]);
    myRoc[r]->SetLineWidth(widths[r]);
    myRoc[r]->GetXaxis()->SetLabelSize(0.04);
    myRoc[r]->GetYaxis()->SetLabelSize(0.04);
    myRoc[r]->GetXaxis()->SetTitleSize(0.05);
    myRoc[r]->GetYaxis()->SetTitleSize(0.05);
    myRoc[r]->GetXaxis()->SetTitleOffset(1.1);
    //myRoc[r]->GetYaxis()->SetTitleOffset(1.3);
    //myRoc[r]->GetYaxis()->SetRangeUser(minY,1.);
    //myRoc[r]->GetXaxis()->SetRangeUser(minX,1.);
    myRoc[r]->Draw(!init?"AL":"L");
    printInfo(rocCutLabel[rocrefs[r]],0.65,0.17+nrocs*0.06+0.02,0.93,0.17+nrocs*0.06+0.07);
    lgd->AddEntry(myRoc[r],varLabel[rocRef[rocrefs[r]][0]].c_str(),"l");
    if(DRAWWPMARKER && gWp[rocrefs[r]]){
      myWp[r] = (TGraph*)gWp[rocrefs[r]]->Clone();
      myWp[r]->SetMarkerColor(colors[r]);
      myWp[r]->SetMarkerStyle(34);
      myWp[r]->SetMarkerSize(2);
      myWp[r]->Draw("LP");
    }
    init = true;
    myRoc[r]->Write(("gRoc_"+varName[rocRef[rocrefs[r]][0]]+"_"+sProcess[rocRef[rocrefs[r]][1]]+"_"+sProcess[rocRef[rocrefs[r]][2]]+"_"+tagOut).c_str());
  }
  lgd->SetFillStyle(0);
  lgd->SetBorderSize(0);
  lgd->SetTextSize(0.037);
  lgd->Draw();
  SaveCanvas(outDir,c1,tagOut);
  fOutRocGraphs->Close();
  delete fOutRocGraphs;
}

void DrawEvw(TCanvas* c1, TGraph** g, Int_t e) {

  setTDRStyle();
  gStyle->SetOptTitle(0);
  gStyle->SetGridColor(kGray);
  c1->cd();
  c1->SetGrid();

  TLegend* lgd = new TLegend(0.72,0.62,0.92,0.92);
  lgd->SetFillStyle(0);
  lgd->SetBorderSize(0);
  lgd->SetTextSize(0.04);

  bool init = true;
  for(int pr=0; pr<nProcesses; pr++){
    if(!isProcessed[pr]) continue;
    g[pr]->SetLineColor(processColor[pr]);
    g[pr]->SetMarkerColor(processColor[pr]);
    g[pr]->SetLineWidth(1);
    g[pr]->GetXaxis()->SetLabelSize(0.04);
    g[pr]->GetYaxis()->SetLabelSize(0.04);
    g[pr]->GetXaxis()->SetTitleSize(0.05);
    g[pr]->GetYaxis()->SetTitleSize(0.05);
    g[pr]->GetXaxis()->SetTitleOffset(1.1);
    //g[pr]->GetYaxis()->SetTitleOffset(1.3);
    g[pr]->Draw(init?"AL":"L");
    lgd->AddEntry( g[pr], processLabel[pr].c_str(), "l" );
    init = false;
  }

  lgd->Draw();

  printInfo(evwCutLabel[e],0.65,0.2,0.93,0.25);
}

void DrawHorizontal(TH1 *h, vector<string> basketLabel, Bool_t same = false, Int_t nDiv = 0, Bool_t without1stBin = false) {
  Double_t ymin = 0.; //h->GetMinimum();
  Double_t ymax = h->GetMaximum(); // 1.1*h->GetMaximum(); 
  TAxis *axis   = h->GetXaxis();
  Double_t xmin = axis->GetXmin();
  Double_t xmax = axis->GetXmax();
  Int_t nbins   = axis->GetNbins();
  TH2F *h2 = new TH2F(Form("h2_%i",rand()),Form("%s;%s;%s",h->GetTitle(),h->GetYaxis()->GetTitle(),""),10,ymin,ymax,nbins,xmin,xmax);
  h2->SetBit(kCanDelete);
  //h2->SetDirectory(0);
  h2->SetTickLength(0.01);
  h2->GetXaxis()->SetLabelSize(0.035);
  h2->GetXaxis()->SetTitleSize(0.045);
  h2->GetXaxis()->SetTitleOffset(1.1);
  if(without1stBin) h2->GetYaxis()->SetRangeUser(0,nBaskets-1);
  h2->SetStats(0);
  prepareBasketlabelsHorizontal(h2,basketLabel);
  h2->Draw(same?"same":"");
  if(nDiv!=0) h2->GetXaxis()->SetNdivisions(nDiv);
  TBox box;
  Int_t color = h->GetFillColor();
  if (color == 0) color = 1;
  box.SetFillColor(color);
  Double_t dy,x1,y1,x2,y2;
  for (Int_t i=1;i<=nbins;i++) {
    if(without1stBin && i==1) continue;
    dy = axis->GetBinWidth(nbins+1-i);
    x1 = 0.;
    y1 = axis->GetBinCenter(nbins+1-i)-0.4*dy;
    x2 = h->GetBinContent(i);
    y2 = axis->GetBinCenter(nbins+1-i)+0.4*dy;
    box.DrawBox(x1,y1,x2,y2);
  }
}

void DrawBasketEfficiencies(TCanvas* c, TH1F* h1, int pr, TH1F** hAssocDecay, string* assocDecayName, Bool_t H2l2XAsBkgd, vector<string> basketLabel, string lumiText, Bool_t logY = false) {

  bool without1stBin = !DRAWLINEALL;

  gStyle->SetOptTitle(0);
  //gStyle->SetTitleX(0.25);

  c->cd();
  if(logY) c->SetLogy();
  //c->SetTicks(0,0);
  c->SetLeftMargin(0.25);

  bool doAssoc = (pr==WH||pr==ZH||pr==ttH);
  Int_t nDecays = pr==WH ? nAssocWDecays : pr==ZH ? nAssocZDecays : pr==ttH ? nAssocttDecays : 0;
  Int_t startAt = pr==WH ? 0 : pr==ZH ? nAssocWDecays : pr==ttH ? nAssocWDecays+nAssocZDecays : 0;

  TH1F* h = (TH1F*)h1->Clone();
  if(H2l2XAsBkgd){
    if(pr==ZH) h->Add(hAssocDecay[nAssocWDecays+nAssocZDecays-1],-1);
    if(pr==ttH) h->Add(hAssocDecay[nAssocWDecays+nAssocZDecays+nAssocttDecays-1],-1);
  }
  if(doAssoc) h->SetLineColor(1);
  if(doAssoc) h->SetFillColor(1);
  h->SetTitle(processLabel[pr].c_str()); 
  h->SetStats(0);
  //if(!logY) h->SetMinimum(0);
  //h->GetXaxis()->SetRangeUser(1,nBaskets-1);
  //h->Draw(); 
  DrawHorizontal(h,basketLabel,false,508,without1stBin);

  TH1F* hStacks[nDecays];
  for(int a=0; a<nDecays; a++){
    if(H2l2XAsBkgd && (pr==ZH||pr==ttH) && a==nDecays-1) continue;
    hStacks[a] = (TH1F*)hAssocDecay[startAt+a]->Clone();
    if(a>0) hStacks[a]->Add(hStacks[a-1]);
    hStacks[a]->SetStats(0);
  }
  bool first = true;
  for(int a=nDecays-1; a>=0; a--){
    if(H2l2XAsBkgd && (pr==ZH||pr==ttH) && a==nDecays-1) continue;
    //hStacks[a]->Draw("same");
    DrawHorizontal(hStacks[a],basketLabel,!first,0,without1stBin);
    first = false;
  }
  
  if(doAssoc){
    TLegend* lgd = new TLegend((H2l2XAsBkgd?0.78:0.63),0.16,0.94,((pr==WH||(pr==ZH&&H2l2XAsBkgd))? 0.27 : (pr==ZH||(pr==ttH&&H2l2XAsBkgd))? 0.33 : pr==ttH? 0.39 : 0.));
    //lgd->SetFillColor(0);
    lgd->SetFillStyle(0);
    lgd->SetBorderSize(0);
    lgd->SetTextFont(42);
    for(int a=0; a<nDecays; a++){
      if(H2l2XAsBkgd && (pr==ZH||pr==ttH) && a==nDecays-1) continue;
      lgd->AddEntry(hStacks[a],assocDecayName[startAt+a].c_str(),"f");
    }
    lgd->Draw();
  }

  float prLabelY = (!doAssoc) ? 0.15 : (pr==WH||(pr==ZH&&H2l2XAsBkgd))? 0.28 : (pr==ZH||(pr==ttH&&H2l2XAsBkgd))? 0.34 : pr==ttH? 0.40 : 0. ;
  TPaveText* prLabel = printInfo(processLabel[pr],doAssoc?(H2l2XAsBkgd?0.78:0.63):0.7,prLabelY,0.94,prLabelY+0.07);
  if(!doAssoc) prLabel->SetTextAlign(32);

  if(PRINTYIELDS){
    TPaveText* pavNExp = new TPaveText(0.905,0.1,0.995,0.9,"brNDC");
    pavNExp->SetFillStyle(0);
    pavNExp->SetBorderSize(0);
    for(int b=2; b<=h->GetNbinsX(); b++){
      pavNExp->AddText(Form("%.4f",h->GetBinContent(b)));
    }
    pavNExp->Draw();
  }

  if(PRINTPLOTINFO) printInfo(INFO,0.23,0.95,0.75,0.99);

  gPad->RedrawAxis();

  //----- print official CMS labels and lumi
  writeExtraText = true;
  extraText  = "Preliminary";
  lumi_sqrtS = string(Form("%s (13 TeV)",lumiText.c_str()));
  CMS_lumi( c, 0, 0 );

}

vector<TH1F*> StackForPurity(vector<TH1F*> histos) {

  vector<TH1F*> result;

  int n = histos.size();
  for(int i=0; i<n; i++){
    TH1F* hTemp = (TH1F*)histos[i]->Clone();
    if(i>0) hTemp->Add(result[i-1]);
    result.push_back(hTemp);
  }

  int nBins = histos[0]->GetNbinsX();
  for(int b=1; b<=nBins; b++){
    Float_t norm = result[n-1]->GetBinContent(b);
    for(int i=0; i<n; i++){
      result[i]->SetBinContent(b,result[i]->GetBinContent(b)/norm);
    }
  }

  return result;
}

void DrawBasketPurities(TCanvas* c, TH1F** h, TH1F** hAssocDecay, string* assocDecayName, Bool_t withAssocDecays, Bool_t H2l2XAsBkgd, vector<string> basketLabel, string* labelMerge, string lumiText) {

  bool without1stBin = !DRAWLINEALL;

  float offset = 0.1;

  gStyle->SetOptTitle(0);
  //gStyle->SetFrameLineWidth(2);

  c->cd();
  c->SetLeftMargin(0.15+offset);
  c->SetRightMargin(0.35-offset);

  vector<TH1F*> histos;
  
  TLegend* lgd = new TLegend(0.67+offset,withAssocDecays?(H2l2XAsBkgd?0.5:0.4):0.7,((0.95+offset<0.99)?0.95+offset:0.99),0.95);
  lgd->SetFillColor(0);
  lgd->SetBorderSize(0);
  lgd->SetTextFont(42);
  
  TPaveText* pavNExp = new TPaveText(0.15+offset,0.13,0.4+offset,0.95,"brNDC");
  pavNExp->SetFillStyle(0);
  pavNExp->SetBorderSize(0);
  pavNExp->SetTextAlign(12);
  pavNExp->SetTextColor(0);
  pavNExp->SetTextFont(42);
  pavNExp->SetTextSize(longDim[BASKETLIST]>500?0.02:0.04);
  TH1F* hNExpectedEvts = (TH1F*)h[0]->Clone();

  TPaveText* pavLM = new TPaveText(0.65+offset,0.13,0.95,0.95,"brNDC");
  pavLM->SetFillStyle(0);
  pavLM->SetBorderSize(0);
  pavLM->SetTextAlign(12);
  pavLM->SetTextFont(42);
  pavLM->SetTextSize(0.02);

  for(int pr=0; pr<nSignalProcesses; pr++){
    if(!isProcessed[pr]) continue;
    bool doAssoc = (pr==WH||pr==ZH||pr==ttH);
    Int_t nDecays = pr==WH ? nAssocWDecays : pr==ZH ? nAssocZDecays : pr==ttH ? nAssocttDecays : 0;
    Int_t startAt = pr==WH ? 0 : pr==ZH ? nAssocWDecays : pr==ttH ? nAssocWDecays+nAssocZDecays : 0;

    TH1F* hTemp = (TH1F*)h[pr]->Clone();
    if(H2l2XAsBkgd){
      if(pr==ZH) hTemp->Add(hAssocDecay[nAssocWDecays+nAssocZDecays-1],-1);
      if(pr==ttH) hTemp->Add(hAssocDecay[nAssocWDecays+nAssocZDecays+nAssocttDecays-1],-1);
    }

    if(withAssocDecays && doAssoc){
      for(int a=0; a<nDecays; a++){
	bool isH2l2X = ((pr==ZH||pr==ttH) && a==nDecays-1); 
	if(!(H2l2XAsBkgd && isH2l2X)){
	  histos.push_back(hAssocDecay[startAt+a]);
	  lgd->AddEntry(hAssocDecay[startAt+a],assocDecayName[startAt+a].c_str(),"f");
	}
      }
    }else{
      histos.push_back(hTemp);
      lgd->AddEntry(hTemp,processLabel[pr].c_str(),"f");
    }

    if(pr!=0) hNExpectedEvts->Add(hTemp);
  }

  vector<TH1F*> stackedHistos = StackForPurity(histos); 
  int n = stackedHistos.size();
  for(int i=0; i<n; i++){
    stackedHistos[n-1-i]->SetStats(0);
    if(i==0){
      stackedHistos[n-1-i]->GetYaxis()->SetRangeUser(0.,1.);
      stackedHistos[n-1-i]->GetYaxis()->SetTitle("signal fraction");
      //stackedHistos[n-1-i]->Draw(); 
      DrawHorizontal(stackedHistos[n-1-i],basketLabel,false,-210,without1stBin);
    }else{
      //stackedHistos[n-1-i]->Draw("same");
      DrawHorizontal(stackedHistos[n-1-i],basketLabel,true,0,without1stBin);
    }
  }

  if(!isSubcat[BASKETLIST]) lgd->Draw();

  for(int b=(without1stBin?2:1); b<=hNExpectedEvts->GetNbinsX(); b++){
    //pavNExp->AddText(Form("%.3f exp. events in %s",hNExpectedEvts->GetBinContent(b),lumiText.c_str()));
    //pavNExp->AddText(Form(longDim[BASKETLIST]>500?"%.6f exp. events":"%.3f exp. events",hNExpectedEvts->GetBinContent(b)));
    pavNExp->AddText(Form(longDim[BASKETLIST]>500?"%.6f exp. events":"%.2f exp. events",hNExpectedEvts->GetBinContent(b)));
    if(isSubcat[BASKETLIST]) pavLM->AddText(labelMerge[b].c_str());
  }
  pavNExp->Draw();
  if(isSubcat[BASKETLIST]) pavLM->Draw();

  if(PRINTPLOTINFO) printInfo(INFO,0.23,0.95,0.75,0.99);
  
  gPad->RedrawAxis();

  //----- print official CMS labels and lumi
  writeExtraText = true;
  extraText  = "Preliminary";
  lumi_sqrtS = string(Form("%s (13 TeV)",lumiText.c_str()));
  CMS_lumi( c, 0, 0 );

}

void DrawBasketSOSPB(TCanvas* c, TH1F* h, vector<string> basketLabel) {

  gStyle->SetOptTitle(0);

  c->cd();
  c->SetLeftMargin(0.25);

  h->SetFillColor(kBlue+2);
  h->GetYaxis()->SetRangeUser(0.,1.);
  h->GetYaxis()->SetTitle("S/(S+B)");
  DrawHorizontal(h,basketLabel,false,-210);

  if(PRINTPLOTINFO) printInfo(INFO,0.23,0.95,0.75,0.99);

}

void DrawSorted(TCanvas* c, TH1F** hGen, TH1F** hReco, string pn, string channel, Float_t maxY = -1.) {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  c->cd();

  Color_t colorGen [4] = { kBlue, kRed+1, kGreen+2, kGray+1 };
  Color_t colorReco[4] = { kBlue, kRed+1, kGreen+2, kGray+1 };

  Float_t SF = 1/hGen[0]->Integral();

  for(int sl=0; sl<4; sl++){

    hGen[sl]->Scale(SF);
    hGen[sl]->SetLineColor(colorReco[sl]);
    hGen[sl]->SetLineWidth(2);
    hGen[sl]->SetLineStyle(2);
    hGen[sl]->SetFillColor(0);
    hGen[sl]->SetFillStyle(0);
    hGen[sl]->GetYaxis()->SetTitle("a.u.");
    //hGen[sl]->GetYaxis()->SetTitleOffset(0.5);
    hGen[sl]->GetYaxis()->SetLabelSize(0.03);
    if(maxY>0.) hGen[sl]->SetMaximum(maxY);

    hReco[sl]->Scale(SF);
    hReco[sl]->SetLineColor(colorGen[sl]);
    hReco[sl]->SetLineWidth(2);
    hReco[sl]->SetLineStyle(2);
    hReco[sl]->SetFillColor(colorGen[sl]);
    hReco[sl]->SetFillStyle(3002);
    hReco[sl]->GetYaxis()->SetTitle("a.u.");
    //hReco[sl]->GetYaxis()->SetTitleOffset(0.5);
    hReco[sl]->GetYaxis()->SetLabelSize(0.03);

  }

  hGen[3]->Draw();
  for(int sl=2;sl>=0;sl--) hGen [sl]->Draw("SAME");
  for(int sl=3;sl>=0;sl--) hReco[sl]->Draw("SAME");

  TLegend *leg = new TLegend(0.45,0.6,0.88,0.88);
  leg->SetLineColor(0);
  leg->SetFillColor(0);  
  TString legtitle; legtitle.Form("H#rightarrow ZZ* #rightarrow %s",channel.c_str());
  //TString massstring; massstring.Form("m_{H} = %d GeV",mass);
  leg->AddEntry((TObject*)0, legtitle.Data(), "");
  //leg->AddEntry((TObject*)0, massstring.Data(), "");
  TLegendEntry *before = leg->AddEntry("Before Analysis Selection","Before Analysis Selection","f");
  TLegendEntry *after  = leg->AddEntry("After Analysis Selection","After Analysis Selection","f");
  after->SetLineStyle(2);
  after->SetLineWidth(2);
  after->SetFillColor(kBlack);
  after->SetLineColor(kBlack);
  after->SetFillStyle(3002);
  before->SetLineStyle(2);
  before->SetLineWidth(2);
  before->SetFillColor(0);
  before->SetLineColor(kBlack);

  leg->Draw("SAME");

  gPad->RedrawAxis();

  if(PRINTPLOTINFO) printInfo(INFO2,0.1,0.9,0.75,0.94);

}

void Draw1D(TCanvas* c, TH1F* h) {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  c->cd();

  //Float_t SF = 1/h->Integral();
  //h->Scale(SF);
  h->SetLineColor(kBlue);
  h->SetLineWidth(2);
  h->SetLineStyle(1);
  h->SetFillColor(0);
  h->SetFillStyle(0);
  //h->GetYaxis()->SetTitle("a.u.");
  //h->GetYaxis()->SetTitleOffset(0.5);
  //h->GetYaxis()->SetLabelSize(0.03);

  h->Draw();

  if(PRINTPLOTINFO) printInfo(INFO2,0.1,0.9,0.75,0.94);

}





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// MAIN MACRO //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void run(	
	 string inputFilePath,	   
	 string outDir,
	 float lumi,
	 string lumiText,
	 double m4lMin, 
	 double m4lMax,
	 string tagOut = ""
	 )
{

  srand(time(0));

  ofstream txtOut;
  TString txtOutName = (TString)outDir.c_str()+"/matchingInfo.txt";
  txtOut.open(txtOutName);

  //TFile* ggZZKFactorFile = TFile::Open("../../data/kfactors/Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
  //TSpline3* sp = (TSpline3*)ggZZKFactorFile->Get("sp_kfactor_Nominal");


  // ------------------------------------------------------------
  // ---------------------- Definitions -------------------------
  // ------------------------------------------------------------

  const int nDatasets = 14;
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
  };

  Color_t allColorsMP[nMatchHLepsStatuses][nSignalProcesses] = {
    { kBlue+2 , kGreen+3, kRed+2 , kOrange+4, kMagenta+2 , kYellow+3 },
    { kBlue-4 , kGreen+2, kRed-4 , kOrange+9, kMagenta   , kYellow+2 },
    { kBlue-7 , kGreen+1, kRed-7 , kOrange+7, kMagenta-7 , kYellow+1 },
    { kBlue-9 , kGreen  , kRed-9 , kOrange+1, kMagenta-9 , kYellow-6 },
    { kBlue-10, kGreen-9, kRed-10, kOrange-9, kMagenta-10, kYellow-8 },
    { kGray+1 , kGray+1 , kGray+1, kGray+1  , kGray+1    , kGray+1   },
  };
  Color_t allColorsPM1[nSignalProcesses][nMatchHLepsStatuses];
  for(int pr=0; pr<nSignalProcesses; pr++){
    for(int m=0; m<nMatchHLepsStatuses; m++){
      allColorsPM1[pr][m] = allColorsMP[m][pr];
    }
  }
  Color_t allColorsPM2[nSignalProcesses][nMatchHLepsStatuses];
  for(int pr=0; pr<nSignalProcesses; pr++){
    for(int m=0; m<nMatchHLepsStatuses; m++){
      if(m<3)
	allColorsPM2[pr][m] = allColorsMP[m][pr];
      else
	allColorsPM2[pr][m] = allColorsMP[m+1][pr];
    }
  }

  string genChannels[nGenChannels] = {
    "4mu",
    "4e",
    "2e2mu",
    "4tau",
    "2e2tau",
    "2mu2tau",
    "2L2X1",
    "4l",
    "2L2tau",
    "2L2X2",
    "all",
  };
  string recoChannels[nRecoChannels] = {
    "4mu",
    "4e",
    "2e2mu",
    "4l",
    "not4l1",
    "not4l2",
  };

  string nbZDaughtersFromH[4] = { "2", "1", "0", "ambig." };

  string WHdecays[nAssocWDecays] = {
    "H->ZZ->4l, W->X     ",
    "H->ZZ->4l, W->lnu   ",
  };
  string ZHdecays[nAssocZDecays] = {
    "H->ZZ->4l, Z->X     ",
    "H->ZZ->4l, Z->2l    ",
    "H->ZZ->2l2X, Z->2l  ",
  };
  string ttHdecays[nAssocttDecays] = {
    "H->ZZ->4l, tt->X    ",
    "H->ZZ->4l, tt->lX   ",
    "H->ZZ->4l, tt->2lX  ",
    "H->ZZ->2l2X, tt->2lX",
  };
  string assocDecayName1[nAssocDecays] = {
    "H#rightarrowZZ#rightarrow4#font[12]{l}, W#rightarrowX",
    "H#rightarrowZZ#rightarrow4#font[12]{l}, W#rightarrow#font[12]{l}#nu",
    "H#rightarrowZZ#rightarrow4#font[12]{l}, Z#rightarrowX",
    "H#rightarrowZZ#rightarrow4#font[12]{l}, Z#rightarrow2#font[12]{l}",
    "H#rightarrowZZ#rightarrow2#font[12]{l}2X, Z#rightarrow2#font[12]{l}",
    "H#rightarrowZZ#rightarrow4#font[12]{l}, t#bar{t}#rightarrow0#font[12]{l}+X",
    "H#rightarrowZZ#rightarrow4#font[12]{l}, t#bar{t}#rightarrow1#font[12]{l}+X",
    "H#rightarrowZZ#rightarrow4#font[12]{l}, t#bar{t}#rightarrow2#font[12]{l}+X",
    "H#rightarrowZZ#rightarrow2#font[12]{l}2X, t#bar{t}#rightarrow2#font[12]{l}+X",
  };
  string assocDecayName2[nAssocDecays] = {
    " W#rightarrowX",
    " W#rightarrow#font[12]{l}#nu",
    " Z#rightarrowX",
    " Z#rightarrow2#font[12]{l}",
    " ERROR",
    " t#bar{t}#rightarrow0#font[12]{l}+X",
    " t#bar{t}#rightarrow1#font[12]{l}+X",
    " t#bar{t}#rightarrow2#font[12]{l}+X",
    " ERROR",
  };
  string assocDecayName3[nAssocDecays] = {
    "WH, H#rightarrow4#font[12]{l}, W#rightarrowX",
    "WH, H#rightarrow4#font[12]{l}, W#rightarrow#font[12]{l}#nu",
    "ZH, H#rightarrow4#font[12]{l}, Z#rightarrowX",
    "ZH, H#rightarrow4#font[12]{l}, Z#rightarrow2#font[12]{l}",
    "ZH, H#rightarrow2#font[12]{l}2X, Z#rightarrow2#font[12]{l}",
    "t#bar{t}H, H#rightarrow4#font[12]{l}, t#bar{t}#rightarrow0#font[12]{l}+X",
    "t#bar{t}H, H#rightarrow4#font[12]{l}, t#bar{t}#rightarrow1#font[12]{l}+X",
    "t#bar{t}H, H#rightarrow4#font[12]{l}, t#bar{t}#rightarrow2#font[12]{l}+X",
    "t#bar{t}H, H#rightarrow2#font[12]{l}2X, t#bar{t}#rightarrow2#font[12]{l}+X",
  };
  string assocDecayName4[nAssocDecays] = {
    "WH, W#rightarrowX",
    "WH, W#rightarrow#font[12]{l}#nu",
    "ZH, Z#rightarrowX",
    "ZH, Z#rightarrow2#font[12]{l}",
    "ERROR",
    "t#bar{t}H, t#bar{t}#rightarrow0#font[12]{l}+X",
    "t#bar{t}H, t#bar{t}#rightarrow1#font[12]{l}+X",
    "t#bar{t}H, t#bar{t}#rightarrow2#font[12]{l}+X",
    "ERROR",
  };
  Color_t assocDecayColor[nAssocDecays] = {
    kRed,
    kRed+1,
    kOrange+1, 
    kOrange+7,
    kOrange-8,
//     kViolet-4,
//     kViolet-5,
//     kViolet-6,
//     kViolet-8, 
    kMagenta-7, 
    kMagenta-3, 
    kMagenta+2, 
    kMagenta-8,
  };

  const Int_t nDecays = 5;
  string decayLabel[nDecays] = { ", all", ", H#rightarrow4l", ", H#rightarrow2l2X", ", H#rightarrow4l, 4 from H", ", H#rightarrow4l, <4 from H" };
  string decayInfix[nDecays] = { "HtoAny", "Hto4l", "Hto2l2X", "Hto4l4fromH", "Hto4lnot4fromH" };

  string matchHLepsKeys[nMatchHLepsStatuses] = {
    "4 matches",
    "3 matches",
    "2 matches",
    "1 match",
    "0 match",
    "ambiguous",
  };
  string matchHLepsInfix[nMatchHLepsStatuses] = {
    "4",
    "3",
    "2",
    "1",
    "0",
    "Ambig",
  };

  string matchAllLepsKeys[nMatchAllLepsStatuses] = {
    "4 from H",
    "3 from H, 1 assoc.",
    "2 from H, 2 assoc.",
    "< 4 matches",
    "ambiguous",
  };
  string matchAllLepsInfix[nMatchAllLepsStatuses] = {
    "4h0a",
    "3h1a",
    "2h2a",
    "Other",
    "Ambig",
  };

  string matchWH[nMatchWHStatuses] = {
    "H#rightarrowZZ#rightarrow4l, W#rightarrowX ; 4 from H",
    "H#rightarrowZZ#rightarrow4l, W#rightarrowl#nu ; 4 from H",
    "H#rightarrowZZ#rightarrow4l, W#rightarrowl#nu ; 3 from H, 1 from W",
    "< 4 matches",
    "ambiguous",
  };
  string matchZH[nMatchZHStatuses] = {
    "H#rightarrowZZ#rightarrow4l, Z#rightarrowX ; 4 from H",
    "H#rightarrowZZ#rightarrow4l, Z#rightarrow2l ; 4 from H",
    "H#rightarrowZZ#rightarrow4l, Z#rightarrow2l ; 3 from H, 1 from Z",
    "H#rightarrowZZ#rightarrow4l, Z#rightarrow2l ; 2 from H, 2 from Z",
    "H#rightarrowZZ#rightarrow2l2X, Z#rightarrow2l ; 2 from H, 2 from Z",
    "< 4 matches",
    "ambiguous",
  };
  string matchttH[nMatchttHStatuses] = {
    "H#rightarrowZZ#rightarrow4l, t#bar{t}#rightarrowX ; 4 from H",
    "H#rightarrowZZ#rightarrow4l, t#bar{t}#rightarrow1l+X ; 4 from H",
    "H#rightarrowZZ#rightarrow4l, t#bar{t}#rightarrow1l+X ; 3 from H, 1 from t#bar{t}",
    "H#rightarrowZZ#rightarrow4l, t#bar{t}#rightarrow2l+X ; 4 from H",
    "H#rightarrowZZ#rightarrow4l, t#bar{t}#rightarrow2l+X ; 3 from H, 1 from t#bar{t}",
    "H#rightarrowZZ#rightarrow4l, t#bar{t}#rightarrow2l+X ; 2 from H, 2 from t#bar{t}",
    "H#rightarrowZZ#rightarrow2l2X, t#bar{t}#rightarrow2l+X ; 2 from H, 2 from t#bar{t}",
    "< 4 matches",
    "ambiguous",
  };

  Color_t colorsMatchWH[nMatchWHStatuses] = {
    kGreen+2,
    kAzure,
    kAzure-4,
    kGray,
    kGray+1,
  };
  Color_t colorsMatchZH[nMatchZHStatuses] = {
    kGreen+2,
    kViolet+2,
    kViolet+1,
    kViolet-9,
    kOrange+7,
    kGray,
    kGray+1,
  };
  Color_t colorsMatchttH[nMatchttHStatuses] = {
    kGreen+2,
    kAzure,
    kAzure-4,
    kViolet+2,
    kViolet+1,
    kViolet-9,
    kOrange+7,
    kGray,
    kGray+1,
  };


  string labelMerge[100];
  vector<string> basketLabel[nBasketLists];
  basketLabel[9].push_back("all");
  basketLabel[9].push_back("0j 0l");
  basketLabel[9].push_back("0j 1l");
  basketLabel[9].push_back("0j #geq1l^{+}l^{-}");
  basketLabel[9].push_back("0j #geq2l else");
  basketLabel[9].push_back("1j 0b 0l VBF-1j");
  basketLabel[9].push_back("1j 0b 0l else");
  basketLabel[9].push_back("1j 0b 1l");
  basketLabel[9].push_back("1j 0b #geq1l^{+}l^{-}");
  basketLabel[9].push_back("1j 0b #geq2l else");
  basketLabel[9].push_back("1j 1b 0l VBF-1j");
  basketLabel[9].push_back("1j 1b 0l else");
  basketLabel[9].push_back("1j 1b 1l");
  basketLabel[9].push_back("1j 1b #geq1l^{+}l^{-}");
  basketLabel[9].push_back("1j 1b #geq2l else");
  basketLabel[9].push_back("2j 0b 0l VBF-2j");
  basketLabel[9].push_back("2j 0b 0l VH-had");
  basketLabel[9].push_back("2j 0b 0l else");
  basketLabel[9].push_back("2j 0b 1l");
  basketLabel[9].push_back("2j 0b #geq1l^{+}l^{-}");
  basketLabel[9].push_back("2j 0b #geq2l else");
  basketLabel[9].push_back("2j 1b 0l VBF-2j");
  basketLabel[9].push_back("2j 1b 0l VH-had");
  basketLabel[9].push_back("2j 1b 0l else");
  basketLabel[9].push_back("2j 1b 1l");
  basketLabel[9].push_back("2j 1b #geq1l^{+}l^{-}");
  basketLabel[9].push_back("2j 1b #geq2l else");
  basketLabel[9].push_back("2j 2b 0l VH-had");
  basketLabel[9].push_back("2j 2b 0l else");
  basketLabel[9].push_back("2j 2b 1l");
  basketLabel[9].push_back("2j 2b #geq1l^{+}l^{-}");
  basketLabel[9].push_back("2j 2b #geq2l else");
  basketLabel[9].push_back("3j 0b 0l VBF-2j");
  basketLabel[9].push_back("3j 0b 0l VH-had");
  basketLabel[9].push_back("3j 0b 0l else");
  basketLabel[9].push_back("3j 0b 1l");
  basketLabel[9].push_back("3j 0b #geq1l^{+}l^{-}");
  basketLabel[9].push_back("3j 0b #geq2l else");
  basketLabel[9].push_back("3j 1b 0l VBF-2j");
  basketLabel[9].push_back("3j 1b 0l VH-had");
  basketLabel[9].push_back("3j 1b 0l else");
  basketLabel[9].push_back("3j 1b 1l");
  basketLabel[9].push_back("3j 1b #geq1l^{+}l^{-}");
  basketLabel[9].push_back("3j 1b #geq2l else");
  basketLabel[9].push_back("3j #geq2b 0l VH-had");
  basketLabel[9].push_back("3j #geq2b 0l else");
  basketLabel[9].push_back("3j #geq2b 1l");
  basketLabel[9].push_back("3j #geq2b #geq1l^{+}l^{-}");
  basketLabel[9].push_back("3j #geq2b #geq2l else");
  basketLabel[9].push_back("#geq4j 0b 0l VBF-2j");
  basketLabel[9].push_back("#geq4j 0b 0l VH-had");
  basketLabel[9].push_back("#geq4j 0b 0l else");
  basketLabel[9].push_back("#geq4j 0b 1l");
  basketLabel[9].push_back("#geq4j 0b #geq1l^{+}l^{-}");
  basketLabel[9].push_back("#geq4j 0b #geq2l else");
  basketLabel[9].push_back("#geq4j 1b 0l VBF-2j");
  basketLabel[9].push_back("#geq4j 1b 0l VH-had");
  basketLabel[9].push_back("#geq4j 1b 0l else");
  basketLabel[9].push_back("#geq4j 1b #geq1l");
  basketLabel[9].push_back("#geq4j #geq2b 0l VH-had");
  basketLabel[9].push_back("#geq4j #geq2b 0l else");
  basketLabel[9].push_back("#geq4j #geq2b #geq1l");
  basketLabel[10].push_back("all");
  basketLabel[10].push_back("Untagged  ");
  //basketLabel[10].push_back("VBF-1j tagged");
  basketLabel[10].push_back("#splitline{VBF-1jet   }{  tagged}");
  //basketLabel[10].push_back("VBF-2j tagged");
  basketLabel[10].push_back("#splitline{VBF-2jet   }{  tagged}");
  basketLabel[10].push_back("#splitline{VH-hadronic}{     tagged}");
  basketLabel[10].push_back("#splitline{VH-leptonic}{    tagged}");
  basketLabel[10].push_back("t#bar{t}H tagged ");
  basketLabel[11].push_back("all");
  basketLabel[11].push_back("0j 0l");
  basketLabel[11].push_back("0j 1l");
  basketLabel[11].push_back("0j #geq1l^{+}l^{-}");
  basketLabel[11].push_back("0j #geq2l else");
  basketLabel[11].push_back("1j 0b 0l VBF-1j");
  basketLabel[11].push_back("1j 0b 0l else");
  basketLabel[11].push_back("1j 0b 1l");
  basketLabel[11].push_back("1j 0b #geq1l^{+}l^{-}");
  basketLabel[11].push_back("1j 0b #geq2l else");
  basketLabel[11].push_back("1j 1b 0l VBF-1j");
  basketLabel[11].push_back("1j 1b 0l else");
  basketLabel[11].push_back("1j 1b 1l");
  basketLabel[11].push_back("1j 1b #geq1l^{+}l^{-}");
  basketLabel[11].push_back("1j 1b #geq2l else");
  basketLabel[11].push_back("2j 0b 0l VBF-2j");
  basketLabel[11].push_back("2j 0b 0l VH-had");
  basketLabel[11].push_back("2j 0b 0l else");
  basketLabel[11].push_back("2j 0b 1l");
  basketLabel[11].push_back("2j 0b #geq1l^{+}l^{-}");
  basketLabel[11].push_back("2j 0b #geq2l else");
  basketLabel[11].push_back("2j 1b 0l VBF-2j");
  basketLabel[11].push_back("2j 1b 0l VH-had");
  basketLabel[11].push_back("2j 1b 0l else");
  basketLabel[11].push_back("2j 1b 1l");
  basketLabel[11].push_back("2j 1b #geq1l^{+}l^{-}");
  basketLabel[11].push_back("2j 1b #geq2l else");
  basketLabel[11].push_back("2j 2b 0l VH-had");
  basketLabel[11].push_back("2j 2b 0l else");
  basketLabel[11].push_back("2j 2b 1l");
  basketLabel[11].push_back("2j 2b #geq1l^{+}l^{-}");
  basketLabel[11].push_back("2j 2b #geq2l else");
  basketLabel[11].push_back("3j 0b 0l VBF-2j");
  basketLabel[11].push_back("3j 0b 0l VH-had");
  basketLabel[11].push_back("3j 0b 0l else");
  basketLabel[11].push_back("3j 0b 1l");
  basketLabel[11].push_back("3j 0b #geq1l^{+}l^{-}");
  basketLabel[11].push_back("3j 0b #geq2l else");
  basketLabel[11].push_back("3j 1b 0l VBF-2j");
  basketLabel[11].push_back("3j 1b 0l VH-had");
  basketLabel[11].push_back("3j 1b 0l else");
  basketLabel[11].push_back("3j 1b 1l");
  basketLabel[11].push_back("3j 1b #geq1l^{+}l^{-}");
  basketLabel[11].push_back("3j 1b #geq2l else");
  basketLabel[11].push_back("3j #geq2b 0l VH-had");
  basketLabel[11].push_back("3j #geq2b 0l else");
  basketLabel[11].push_back("3j #geq2b 1l");
  basketLabel[11].push_back("3j #geq2b #geq1l^{+}l^{-}");
  basketLabel[11].push_back("3j #geq2b #geq2l else");
  basketLabel[11].push_back("#geq4j 0b 0l VBF-2j");
  basketLabel[11].push_back("#geq4j 0b 0l VH-had");
  basketLabel[11].push_back("#geq4j 0b 0l else");
  basketLabel[11].push_back("#geq4j 0b 1l");
  basketLabel[11].push_back("#geq4j 0b #geq1l^{+}l^{-}");
  basketLabel[11].push_back("#geq4j 0b #geq2l else");
  basketLabel[11].push_back("#geq4j 1b 0l VBF-2j");
  basketLabel[11].push_back("#geq4j 1b 0l VH-had");
  basketLabel[11].push_back("#geq4j 1b 0l else");
  basketLabel[11].push_back("#geq4j 1b #geq1l");
  basketLabel[11].push_back("#geq4j #geq2b 0l VH-had");
  basketLabel[11].push_back("#geq4j #geq2b 0l else");
  basketLabel[11].push_back("#geq4j #geq2b #geq1l");
  basketLabel[11].push_back("0j 0l MET");
  basketLabel[11].push_back("1j 0b 0l MET");
  basketLabel[11].push_back("1j 1b 0l MET");
  basketLabel[12].push_back("all");
  basketLabel[12].push_back("Untagged  ");
  //basketLabel[12].push_back("VBF-1j tagged");
  basketLabel[12].push_back("#splitline{VBF-1jet   }{  tagged}");
  //basketLabel[12].push_back("VBF-2j tagged");
  basketLabel[12].push_back("#splitline{VBF-2jet   }{  tagged}");
  basketLabel[12].push_back("#splitline{VH-hadronic}{     tagged}");
  basketLabel[12].push_back("#splitline{VH-leptonic}{    tagged}");
  //basketLabel[12].push_back("VH-MET tagged");
  basketLabel[12].push_back("#splitline{VH-MET  }{  tagged}");
  basketLabel[12].push_back("t#bar{t}H tagged ");

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
  Short_t NRecoMu;
  Short_t NRecoEle;
  Short_t Nvtx;
  Short_t NObsInt;
  Float_t NTrueInt;
  Short_t trigWord;
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
  Float_t DiJetFisher;
  Float_t p_GG_SIG_ghg2_1_ghz1_1_JHUGen;
  Float_t p_QQB_BKG_MCFM;
  Float_t p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal;
  Float_t p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal;
  Float_t p_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
  Float_t pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
  Float_t p_JQCD_SIG_ghg2_1_JHUGen_JECNominal;
  Float_t p_HadWH_SIG_ghw1_1_JHUGen_JECNominal;
  Float_t p_LepWH_SIG_ghw1_1_JHUGen;
  Float_t p_HadZH_SIG_ghz1_1_JHUGen_JECNominal;
  Float_t p_LepZH_SIG_ghz1_1_JHUGen;
  Float_t Z1Mass;
  Float_t Z2Mass;
  Short_t Z1Flav;
  Short_t Z2Flav;
  vector<Float_t> *CandLepPt = 0;
  vector<Float_t> *CandLepEta = 0;
  vector<Float_t> *CandLepPhi = 0;
  vector<Short_t> *CandLepId = 0;
  Short_t nExtraLep;
  vector<Float_t> *ExtraLepPt = 0;
  vector<Float_t> *ExtraLepEta = 0;
  vector<Float_t> *ExtraLepPhi = 0;
  vector<Float_t> *ExtraLepId = 0;
  Short_t nExtraZ;
  Short_t nJets;
  Short_t nJetsBTagged;
  vector<Float_t> *JetPt = 0;
  vector<Float_t> *JetEta = 0;
  vector<Float_t> *JetPhi = 0;
  vector<Float_t> *JetMass = 0;
  vector<Float_t> *JetBTagger = 0;
  vector<Float_t> *JetIsBtagged = 0;
  vector<Float_t> *JetQGLikelihood = 0;
  vector<Float_t> *JetPQuark = 0;
  vector<Float_t> *JetPGluon = 0;
  Float_t jetPt[99];
  Float_t jetEta[99];
  Float_t jetPhi[99];
  Float_t jetMass[99];
  Float_t jetQGLikelihood[99];
  Float_t jetQGLikelihoodRaw[99];
  Float_t jetPQuark[99];
  Float_t jetPGluon[99];
  Float_t jetPgOverPq[99];
  Float_t DiJetMass;
  Float_t GenHMass;
  Float_t GenHPt;
  Float_t GenHRapidity;
  Float_t GenLep1Pt;
  Float_t GenLep1Eta;
  Float_t GenLep1Phi;
  Short_t GenLep1Id;
  Float_t GenLep2Pt;
  Float_t GenLep2Eta;
  Float_t GenLep2Phi;
  Short_t GenLep2Id;
  Float_t GenLep3Pt;
  Float_t GenLep3Eta;
  Float_t GenLep3Phi;
  Short_t GenLep3Id;
  Float_t GenLep4Pt;
  Float_t GenLep4Eta;
  Float_t GenLep4Phi;
  Short_t GenLep4Id;
  Float_t GenAssocLep1Pt;
  Float_t GenAssocLep1Eta;
  Float_t GenAssocLep1Phi;
  Short_t GenAssocLep1Id;
  Float_t GenAssocLep2Pt;
  Float_t GenAssocLep2Eta;
  Float_t GenAssocLep2Phi;
  Short_t GenAssocLep2Id;

  Int_t nbStored[nProcesses][nGenChannels];
  Int_t nbHLepsAreInEtaAcc[nProcesses][nGenChannels];
  Int_t nbHLepsAreInEtaPtAcc[nProcesses][nGenChannels];
  Int_t nb4GenLeps[nProcesses][nGenChannels];
  Int_t nb4RecoLeps[nProcesses][nGenChannels];
  Int_t nbWithBCInSRG[nProcesses][nGenChannels];
  Int_t nbPassTriggerG[nProcesses][nGenChannels];
  Int_t nbPassTriggerGNo1E[nProcesses][nGenChannels];
  Int_t nbPassTriggerGWithBCInSRG[nProcesses][nGenChannels];
  Int_t nbPassTriggerGNo1EWithBCInSRG[nProcesses][nGenChannels];
  Int_t nbHLepsAreInEtaPtAcc4RecoLeps[nProcesses][nGenChannels];
  Int_t nbHLepsAreInEtaPtAccWithBCInSRG[nProcesses][nGenChannels];
  Int_t nbHLepsAreInEtaPtAccPassTriggerG[nProcesses][nGenChannels];
  Int_t nbHLepsAreInEtaPtAccPassTriggerGNo1E[nProcesses][nGenChannels];
  Int_t nbHLepsAreInEtaPtAccPassTriggerGWithBCInSRG[nProcesses][nGenChannels];
  Int_t nbHLepsAreInEtaPtAccPassTriggerGNo1EWithBCInSRG[nProcesses][nGenChannels];
  Float_t preselExp[nProcesses];
  Float_t yieldStored[nProcesses][nGenChannels];
  Float_t yieldHLepsAreInEtaAcc[nProcesses][nGenChannels];
  Float_t yieldHLepsAreInEtaPtAcc[nProcesses][nGenChannels];
  Float_t yield4GenLeps[nProcesses][nGenChannels];
  Float_t yield4RecoLeps[nProcesses][nGenChannels];
  Float_t yieldWithBCInSRG[nProcesses][nGenChannels];
  Float_t yieldPassTriggerG[nProcesses][nGenChannels];
  Float_t yieldPassTriggerGNo1E[nProcesses][nGenChannels];
  Float_t yieldPassTriggerGWithBCInSRG[nProcesses][nGenChannels];
  Float_t yieldPassTriggerGNo1EWithBCInSRG[nProcesses][nGenChannels];
  Float_t yieldHLepsAreInEtaPtAcc4RecoLeps[nProcesses][nGenChannels];
  Float_t yieldHLepsAreInEtaPtAccWithBCInSRG[nProcesses][nGenChannels];
  Float_t yieldHLepsAreInEtaPtAccPassTriggerG[nProcesses][nGenChannels];
  Float_t yieldHLepsAreInEtaPtAccPassTriggerGNo1E[nProcesses][nGenChannels];
  Float_t yieldHLepsAreInEtaPtAccPassTriggerGWithBCInSRG[nProcesses][nGenChannels];
  Float_t yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCInSRG[nProcesses][nGenChannels];
  TH1F* h1GenptHLepsAreInEtaAcc[nProcesses][nGenChannels][4];
  TH1F* h1RecoptWithBCInSRGPassTriggerG[nProcesses][nGenChannels][4];
  TH1F* h1GenetaHLepsAreInPtAcc[nProcesses][nGenChannels][4];
  TH1F* h1RecoetaWithBCInSRGPassTriggerG[nProcesses][nGenChannels][4];
  TH1F* h1GenHPt[nProcesses];
  TH1F* h1GenHRapidity[nProcesses];
  TH2F* h2GenHRapidityVsPt[nProcesses];
  TH1F* h1NRecoLepHLepsAreInEtaPtAcc[nProcesses][nGenChannels];
  TH1F* h1NRecoMuoHLepsAreInEtaPtAcc[nProcesses][nGenChannels];
  TH1F* h1NRecoEleHLepsAreInEtaPtAcc[nProcesses][nGenChannels];

  Int_t nbWithBC[nProcesses][nRecoChannels];
  Int_t nbWithBCInSR[nProcesses][nRecoChannels];
  Int_t nbWithBCInSRExactly4GoodLeps[nProcesses][nRecoChannels];
  Int_t nbWithBCInSRAtLeast5GoodLeps[nProcesses][nRecoChannels];
  Int_t nbWithBCInSRExactly5GoodLeps[nProcesses][nRecoChannels];
  Int_t nbWithBCInSRExactly6GoodLeps[nProcesses][nRecoChannels];
  Int_t nbWithBCInSRHLepsAreInEtaPtAcc[nProcesses][nRecoChannels];
  Int_t nbWithBCInSRHLepsAreGood[nProcesses][nRecoChannels];
  Int_t nbWithBCInSRAll4LepRight[nProcesses][nRecoChannels];
  Int_t nbWithBCInSRMatchHLeps[nMatchHLepsStatuses][nProcesses][nRecoChannels];
  Int_t nbWithBCInSRMatchAllLeps[nMatchAllLepsStatuses][nProcesses][nRecoChannels];
  Int_t nbWithBCInSRMatchWH[nMatchWHStatuses][nRecoChannels];
  Int_t nbWithBCInSRMatchZH[nMatchZHStatuses][nRecoChannels];
  Int_t nbWithBCInSRMatchttH[nMatchttHStatuses][nRecoChannels];
  Float_t yieldWithBC[nProcesses][nRecoChannels];
  Float_t yieldWithBCInSR[nProcesses][nRecoChannels];

  TH1F* hBCInSR[nVariables][nProcesses][nRecoChannels];
  TH1F* hBCInSRMatchHLeps[nVariables][nMatchHLepsStatuses][nProcesses][nRecoChannels];
  TH1F* hBCInSRMatchAllLeps[nVariables][nMatchAllLepsStatuses][nProcesses][nRecoChannels];
  TH1F* hBCInSRMatchWH[nVariables][nMatchWHStatuses][nRecoChannels];
  TH1F* hBCInSRMatchZH[nVariables][nMatchZHStatuses][nRecoChannels];
  TH1F* hBCInSRMatchttH[nVariables][nMatchttHStatuses][nRecoChannels];

  TH2F* h2DBCInSR[nVarPairs][nProcesses][nRecoChannels];
  TH2F* h2DBCInSRDecays[nVarPairs][nProcesses][nDecays][nRecoChannels];

  TH1F* hBCInSRBaskets[nProcesses][nRecoChannels];
  TH1F* hBCInSRBasketsAssocDecays[nAssocDecays][nRecoChannels];

  for(int pr=0; pr<nProcesses; pr++){

    preselExp[pr] = 0.;

    for(int gc=0; gc<nGenChannels; gc++){

      string suffix = sProcess[pr]+"_"+genChannels[gc];
      string title = sProcess[pr]+", gen. channel "+genChannels[gc];

      nbStored[pr][gc] = 0;
      nbHLepsAreInEtaAcc[pr][gc] = 0;
      nbHLepsAreInEtaPtAcc[pr][gc] = 0;
      nb4GenLeps[pr][gc] = 0;
      nb4RecoLeps[pr][gc] = 0;
      nbWithBCInSRG[pr][gc] = 0;
      nbPassTriggerG[pr][gc] = 0;
      nbPassTriggerGNo1E[pr][gc] = 0;
      nbPassTriggerGWithBCInSRG[pr][gc] = 0;
      nbPassTriggerGNo1EWithBCInSRG[pr][gc] = 0;
      nbHLepsAreInEtaPtAcc4RecoLeps[pr][gc] = 0;
      nbHLepsAreInEtaPtAccWithBCInSRG[pr][gc] = 0;
      nbHLepsAreInEtaPtAccPassTriggerG[pr][gc] = 0;
      nbHLepsAreInEtaPtAccPassTriggerGNo1E[pr][gc] = 0;
      nbHLepsAreInEtaPtAccPassTriggerGWithBCInSRG[pr][gc] = 0;
      nbHLepsAreInEtaPtAccPassTriggerGNo1EWithBCInSRG[pr][gc] = 0;
      yieldStored[pr][gc] = 0.;
      yieldHLepsAreInEtaAcc[pr][gc] = 0.;
      yieldHLepsAreInEtaPtAcc[pr][gc] = 0.;
      yield4GenLeps[pr][gc] = 0.;
      yield4RecoLeps[pr][gc] = 0.;
      yieldWithBCInSRG[pr][gc] = 0.;
      yieldPassTriggerG[pr][gc] = 0.;
      yieldPassTriggerGNo1E[pr][gc] = 0.;
      yieldPassTriggerGWithBCInSRG[pr][gc] = 0.;
      yieldPassTriggerGNo1EWithBCInSRG[pr][gc] = 0.;
      yieldHLepsAreInEtaPtAcc4RecoLeps[pr][gc] = 0.;
      yieldHLepsAreInEtaPtAccWithBCInSRG[pr][gc] = 0.;
      yieldHLepsAreInEtaPtAccPassTriggerG[pr][gc] = 0.;
      yieldHLepsAreInEtaPtAccPassTriggerGNo1E[pr][gc] = 0.;
      yieldHLepsAreInEtaPtAccPassTriggerGWithBCInSRG[pr][gc] = 0.;
      yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCInSRG[pr][gc] = 0.;

      for(int sl=0; sl<4; sl++){
	h1GenptHLepsAreInEtaAcc[pr][gc][sl] = new TH1F(Form("h1GenptHLepsAreInEtaAcc_%s_%i",suffix.c_str(),sl),(title+";p_{T};").c_str(),80,0,80);
	h1RecoptWithBCInSRGPassTriggerG[pr][gc][sl] = new TH1F(Form("h1RecoptWithBCInSRGPassTriggerG_%s_%i",suffix.c_str(),sl),(title+";p_{T};").c_str(),100,0,100);
	h1GenetaHLepsAreInPtAcc[pr][gc][sl] = new TH1F(Form("h1GenetaHLepsAreInPtAcc_%s_%i",suffix.c_str(),sl),(title+";|#eta|;").c_str(),100,0,5);
	h1RecoetaWithBCInSRGPassTriggerG[pr][gc][sl] = new TH1F(Form("h1RecoetaWithBCInSRGPassTriggerG_%s_%i",suffix.c_str(),sl),(title+";|#eta|;").c_str(),100,0,5);
      }

      h1NRecoLepHLepsAreInEtaPtAcc[pr][gc] = new TH1F(Form("h1NRecoLepHLepsAreInEtaPtAcc_%s",suffix.c_str()),(title+";nb. of reco. leptons;").c_str(),10,0,10);
      h1NRecoMuoHLepsAreInEtaPtAcc[pr][gc] = new TH1F(Form("h1NRecoMuoHLepsAreInEtaPtAcc_%s",suffix.c_str()),(title+";nb. of reco. muons;").c_str(),10,0,10);
      h1NRecoEleHLepsAreInEtaPtAcc[pr][gc] = new TH1F(Form("h1NRecoEleHLepsAreInEtaPtAcc_%s",suffix.c_str()),(title+";nb. of reco. electrons;").c_str(),10,0,10);

    }

    h1GenHPt[pr] = new TH1F(Form("h1GenHPt_%s",sProcess[pr].c_str()),";p_{T}^{H};",100,0,100);
    h1GenHRapidity[pr] = new TH1F(Form("h1GenHRapidity_%s",sProcess[pr].c_str()),";y^{H};",100,0,5);
    h2GenHRapidityVsPt[pr] = new TH2F(Form("h2GenHRapidityVsPt_%s",sProcess[pr].c_str()),";p_{T}^{H};y^{H}",20,0,100,20,0,5);

    for(int rc=0; rc<nRecoChannels; rc++){
      
      string suffix = "_"+sProcess[pr]+"_"+recoChannels[rc];
      string title = sProcess[pr]+", channel "+recoChannels[rc];

      nbWithBC[pr][rc] = 0;
      nbWithBCInSR[pr][rc] = 0;
      nbWithBCInSRExactly4GoodLeps[pr][rc] = 0;
      nbWithBCInSRAtLeast5GoodLeps[pr][rc] = 0;
      nbWithBCInSRExactly5GoodLeps[pr][rc] = 0;
      nbWithBCInSRExactly6GoodLeps[pr][rc] = 0;
      nbWithBCInSRHLepsAreInEtaPtAcc[pr][rc] = 0;
      nbWithBCInSRHLepsAreGood[pr][rc] = 0;
      nbWithBCInSRAll4LepRight[pr][rc] = 0;
      for(int m=0; m<nMatchHLepsStatuses; m++) nbWithBCInSRMatchHLeps[m][pr][rc] = 0;
      for(int m=0; m<nMatchAllLepsStatuses; m++) nbWithBCInSRMatchAllLeps[m][pr][rc] = 0;
      if(pr==WH) for(int m=0; m<nMatchWHStatuses; m++) nbWithBCInSRMatchWH[m][rc] = 0;
      if(pr==ZH) for(int m=0; m<nMatchZHStatuses; m++) nbWithBCInSRMatchZH[m][rc] = 0;
      if(pr==ttH) for(int m=0; m<nMatchttHStatuses; m++) nbWithBCInSRMatchttH[m][rc] = 0;

      yieldWithBC[pr][rc] = 0.;
      yieldWithBCInSR[pr][rc] = 0.;

      for(int v=0; v<nVariables; v++){
	hBCInSR[v][pr][rc] = new TH1F(("hBCInSR_"+varName[v]+suffix).c_str(),Form("%s;%s;# exp. events in %s",title.c_str(),varLabel[v].c_str(),lumiText.c_str()),varNbin[v],varMin[v],varMax[v]);
	for(int m=0; m<nMatchHLepsStatuses; m++)
	  hBCInSRMatchHLeps[v][m][pr][rc] = new TH1F(("hBCInSRMatchHLeps_"+varName[v]+"_"+matchHLepsInfix[m]+suffix).c_str(),Form("%s;%s;# exp. events in %s",title.c_str(),varLabel[v].c_str(),lumiText.c_str()),varNbin[v],varMin[v],varMax[v]);
	for(int m=0; m<nMatchAllLepsStatuses; m++)
	  hBCInSRMatchAllLeps[v][m][pr][rc] = new TH1F(("hBCInSRMatchAllLeps_"+varName[v]+"_"+matchAllLepsInfix[m]+suffix).c_str(),Form("%s;%s;# exp. events in %s",title.c_str(),varLabel[v].c_str(),lumiText.c_str()),varNbin[v],varMin[v],varMax[v]);
	if(pr==WH)
	  for(int m=0; m<nMatchWHStatuses; m++)
	    hBCInSRMatchWH[v][m][rc] = new TH1F(Form("hBCInSRMatchWH_%s_%i_%s",varName[v].c_str(),m,recoChannels[rc].c_str()),Form("%s;%s;# exp. events in %s",title.c_str(),varLabel[v].c_str(),lumiText.c_str()),varNbin[v],varMin[v],varMax[v]);
	if(pr==ZH)
	  for(int m=0; m<nMatchZHStatuses; m++)
	    hBCInSRMatchZH[v][m][rc] = new TH1F(Form("hBCInSRMatchZH_%s_%i_%s",varName[v].c_str(),m,recoChannels[rc].c_str()),Form("%s;%s;# exp. events in %s",title.c_str(),varLabel[v].c_str(),lumiText.c_str()),varNbin[v],varMin[v],varMax[v]);
	if(pr==ttH)
	  for(int m=0; m<nMatchttHStatuses; m++)
	    hBCInSRMatchttH[v][m][rc] = new TH1F(Form("hBCInSRMatchttH_%s_%i_%s",varName[v].c_str(),m,recoChannels[rc].c_str()),Form("%s;%s;# exp. events in %s",title.c_str(),varLabel[v].c_str(),lumiText.c_str()),varNbin[v],varMin[v],varMax[v]);
      }

      for(int v2=0; v2<nVarPairs; v2++){
	h2DBCInSR[v2][pr][rc] = new TH2F(("h2DBCInSR_"+varName[varPairRef[v2][0]]+"_"+varName[varPairRef[v2][1]]+suffix).c_str(),(title+";"+varLabel[varPairRef[v2][0]]+";"+varLabel[varPairRef[v2][1]]).c_str(),varNbin[varPairRef[v2][0]],varMin[varPairRef[v2][0]],varMax[varPairRef[v2][0]],varNbin[varPairRef[v2][1]],varMin[varPairRef[v2][1]],varMax[varPairRef[v2][1]]);
	if(!varPairRef[v2][2]) continue;
	for(int d=0; d<nDecays; d++){
	  h2DBCInSRDecays[v2][pr][d][rc] = new TH2F(("h2DBCInSR_"+decayInfix[d]+"_"+varName[varPairRef[v2][0]]+"_"+varName[varPairRef[v2][1]]+suffix).c_str(),(title+";"+varLabel[varPairRef[v2][0]]+";"+varLabel[varPairRef[v2][1]]).c_str(),varNbin[varPairRef[v2][0]],varMin[varPairRef[v2][0]],varMax[varPairRef[v2][0]],varNbin[varPairRef[v2][1]],varMin[varPairRef[v2][1]],varMax[varPairRef[v2][1]]);
	}
      }
      
      hBCInSRBaskets[pr][rc] = new TH1F(("hBCInSR_baskets_"+suffix).c_str(),Form("%s;;# exp. events in %s",title.c_str(),lumiText.c_str()),nBaskets,0,nBaskets);
      prepareBasketlabels(hBCInSRBaskets[pr][rc], basketLabel[BASKETLIST]);
      hBCInSRBaskets[pr][rc]->SetLineColor(processColor[pr]);
      hBCInSRBaskets[pr][rc]->SetFillColor(processColor[pr]);
      
    }
    
  }

  for(int a=0; a<nAssocDecays; a++){
    for(int rc=0; rc<nRecoChannels; rc++){
      hBCInSRBasketsAssocDecays[a][rc] = new TH1F(Form("hBCInSR_baskets_assocDecays_%i_%s",a,recoChannels[rc].c_str()),Form(";;# exp. events in %s",lumiText.c_str()),nBaskets,0,nBaskets);
      prepareBasketlabels(hBCInSRBasketsAssocDecays[a][rc], basketLabel[BASKETLIST]);
      hBCInSRBasketsAssocDecays[a][rc]->SetLineColor(assocDecayColor[a]);
      hBCInSRBasketsAssocDecays[a][rc]->SetFillColor(assocDecayColor[a]);
    }
  }

  TH1F* hRocSgnl[nROCs];
  TH1F* hRocBkgd[nROCs];
  TGraph* gRoc[nROCs];
  TGraph* gWp[nROCs];
  for(int ro=0; ro<nROCs; ro++){
    hRocSgnl[ro] = new TH1F(Form("hRocSgnl_%s_%s_%s",varName[rocRef[ro][0]].c_str(),sProcess[rocRef[ro][1]].c_str(),sProcess[rocRef[ro][2]].c_str()),Form(";%s;",varLabel[rocRef[ro][0]].c_str()),10000,varMin[rocRef[ro][0]],varMaxRoc[rocRef[ro][0]]);
    hRocBkgd[ro] = new TH1F(Form("hRocBkgd_%s_%s_%s",varName[rocRef[ro][0]].c_str(),sProcess[rocRef[ro][1]].c_str(),sProcess[rocRef[ro][2]].c_str()),Form(";%s;",varLabel[rocRef[ro][0]].c_str()),10000,varMin[rocRef[ro][0]],varMaxRoc[rocRef[ro][0]]);
  }
  TH1F* hEvw[nEVWs][nProcesses];
  TGraph* gEvw[nEVWs][nProcesses];
  for(int e=0; e<nEVWs; e++){
    for(int pr=0; pr<nProcesses; pr++){
      hEvw[e][pr] = new TH1F(Form("hEvw_%s_%s",varName[evwRef[e][0]].c_str(),sProcess[pr].c_str()),Form(";%s;",varLabel[evwRef[e][0]].c_str()),10000,varMin[evwRef[e][0]],varMaxRoc[evwRef[e][0]]);
    }
  }

  Int_t nbTotalWH[nAssocWDecays] = {0,0};
  Int_t nbHLepsAreInEtaPtAccWH[nAssocWDecays] = {0,0};
  Int_t nbHLepsAreGoodWH[nAssocWDecays] = {0,0};
  Int_t nbAll4LepRightWH[nAssocWDecays] = {0,0};
  Int_t nbZ1DaughtersFromHWH[nAssocWDecays][4]; for(int i=0; i<nAssocWDecays; i++) for(int j=0; j<4; j++) nbZ1DaughtersFromHWH[i][j] = 0;
  Int_t nbZ2DaughtersFromHWH[nAssocWDecays][4]; for(int i=0; i<nAssocWDecays; i++) for(int j=0; j<4; j++) nbZ2DaughtersFromHWH[i][j] = 0;

  Int_t nbTotalZH[nAssocZDecays] = {0,0,0};
  Int_t nbHLepsAreInEtaPtAccZH[nAssocZDecays] = {0,0,0};
  Int_t nbHLepsAreGoodZH[nAssocZDecays] = {0,0,0};
  Int_t nbAll4LepRightZH[nAssocZDecays] = {0,0,0};
  Int_t nbZ1DaughtersFromHZH[nAssocZDecays][4]; for(int i=0; i<nAssocZDecays; i++) for(int j=0; j<4; j++) nbZ1DaughtersFromHZH[i][j] = 0;
  Int_t nbZ2DaughtersFromHZH[nAssocZDecays][4]; for(int i=0; i<nAssocZDecays; i++) for(int j=0; j<4; j++) nbZ2DaughtersFromHZH[i][j] = 0;

  Int_t nbTotalttH[nAssocttDecays] = {0,0,0,0};
  Int_t nbHLepsAreInEtaPtAccttH[nAssocttDecays] = {0,0,0,0};
  Int_t nbHLepsAreGoodttH[nAssocttDecays] = {0,0,0,0};
  Int_t nbAll4LepRightttH[nAssocttDecays] = {0,0,0,0};
  Int_t nbZ1DaughtersFromHttH[nAssocttDecays][4]; for(int i=0; i<nAssocttDecays; i++) for(int j=0; j<4; j++) nbZ1DaughtersFromHttH[i][j] = 0;
  Int_t nbZ2DaughtersFromHttH[nAssocttDecays][4]; for(int i=0; i<nAssocttDecays; i++) for(int j=0; j<4; j++) nbZ2DaughtersFromHttH[i][j] = 0;



  // ------------------------------------------------------------
  // ---------------------- Processing --------------------------
  // ------------------------------------------------------------

  int currentProcess;

  for(int d=0; d<nDatasets; d++){
    //cout<<"debug flag dataset "<<d<<endl;

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

    if(!isProcessed[currentProcess]) continue;

    cout<<datasets[d]<<endl;
    txtOut<<datasets[d]<<endl;

    string inputFileName = string(Form("%s%s%s/ZZ4lAnalysis.root",inputFilePath.c_str(),datasets[d].c_str(),isSignal[currentProcess]?"125":""));
    inputFile[d] = TFile::Open(inputFileName.c_str());

    hCounters[d] = (TH1F*)inputFile[d]->Get("ZZTree/Counters");
    NGenEvt[d] = (Long64_t)hCounters[d]->GetBinContent(1);
    cout<<"  number of generated events: "<<NGenEvt[d]<<endl;
    gen_sumWeights[d] = (Long64_t)hCounters[d]->GetBinContent(USEPUWEIGHT?40:41);
    partialSampleWeight[d] = lumi * 1000 / gen_sumWeights[d] ;

    inputTree[d] = (TTree*)inputFile[d]->Get("ZZTree/candTree");
    inputTree[d]->SetBranchAddress("RunNumber", &nRun);
    inputTree[d]->SetBranchAddress("EventNumber", &nEvent);
    inputTree[d]->SetBranchAddress("LumiNumber", &nLumi);
    inputTree[d]->SetBranchAddress("PFMET", &PFMET);
    inputTree[d]->SetBranchAddress("NRecoMu", &NRecoMu);
    inputTree[d]->SetBranchAddress("NRecoEle", &NRecoEle);
    inputTree[d]->SetBranchAddress("Nvtx", &Nvtx);
    inputTree[d]->SetBranchAddress("NObsInt", &NObsInt);
    inputTree[d]->SetBranchAddress("NTrueInt", &NTrueInt);
    inputTree[d]->SetBranchAddress("trigWord", &trigWord);
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
    inputTree[d]->SetBranchAddress("p_GG_SIG_ghg2_1_ghz1_1_JHUGen", &p_GG_SIG_ghg2_1_ghz1_1_JHUGen);
    inputTree[d]->SetBranchAddress("p_QQB_BKG_MCFM", &p_QQB_BKG_MCFM);
    inputTree[d]->SetBranchAddress("Z1Mass", &Z1Mass);
    inputTree[d]->SetBranchAddress("Z2Mass", &Z2Mass);
    inputTree[d]->SetBranchAddress("Z1Flav", &Z1Flav);
    inputTree[d]->SetBranchAddress("Z2Flav", &Z2Flav);
    inputTree[d]->SetBranchAddress("LepPt", &CandLepPt);
    inputTree[d]->SetBranchAddress("LepEta", &CandLepEta);
    inputTree[d]->SetBranchAddress("LepPhi", &CandLepPhi);
    inputTree[d]->SetBranchAddress("LepLepId", &CandLepId);
    inputTree[d]->SetBranchAddress("nExtraLep", &nExtraLep);
    inputTree[d]->SetBranchAddress("ExtraLepPt", &ExtraLepPt);
    inputTree[d]->SetBranchAddress("ExtraLepEta", &ExtraLepEta);
    inputTree[d]->SetBranchAddress("ExtraLepPhi", &ExtraLepPhi);
    inputTree[d]->SetBranchAddress("ExtraLepLepId", &ExtraLepId);
    inputTree[d]->SetBranchAddress("nExtraZ", &nExtraZ);
    inputTree[d]->SetBranchAddress("nCleanedJetsPt30", &nJets);
    inputTree[d]->SetBranchAddress( BTAGGINGSF==0 ? "nCleanedJetsPt30BTagged" :
				    BTAGGINGSF==1 ? "nCleanedJetsPt30BTagged_bTagSF" :
				    BTAGGINGSF==2 ? "nCleanedJetsPt30BTagged_bTagSFUp" :
				    BTAGGINGSF==3 ? "nCleanedJetsPt30BTagged_bTagSFDn" : ""
				    , &nJetsBTagged);
    inputTree[d]->SetBranchAddress("JetPt", &JetPt);
    inputTree[d]->SetBranchAddress("JetEta", &JetEta);
    inputTree[d]->SetBranchAddress("JetPhi", &JetPhi);
    inputTree[d]->SetBranchAddress("JetMass", &JetMass);
    inputTree[d]->SetBranchAddress("JetBTagger", &JetBTagger);
    inputTree[d]->SetBranchAddress("JetIsBtagged", &JetIsBtagged);
    inputTree[d]->SetBranchAddress("JetQGLikelihood", &JetQGLikelihood);
    //inputTree[d]->SetBranchAddress("JetPQuark", &JetPQuark);
    //inputTree[d]->SetBranchAddress("JetPGluon", &JetPGluon);
    inputTree[d]->SetBranchAddress("DiJetMass", &DiJetMass);
    inputTree[d]->SetBranchAddress("DiJetFisher", &DiJetFisher);
    inputTree[d]->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal", &p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal);
    inputTree[d]->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal", &p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal);
    inputTree[d]->SetBranchAddress("p_JVBF_SIG_ghv1_1_JHUGen_JECNominal", &p_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
    inputTree[d]->SetBranchAddress("pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal", &pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
    inputTree[d]->SetBranchAddress("p_JQCD_SIG_ghg2_1_JHUGen_JECNominal", &p_JQCD_SIG_ghg2_1_JHUGen_JECNominal);
    inputTree[d]->SetBranchAddress("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal", &p_HadWH_SIG_ghw1_1_JHUGen_JECNominal);
    inputTree[d]->SetBranchAddress("p_LepWH_SIG_ghw1_1_JHUGen", &p_LepWH_SIG_ghw1_1_JHUGen);
    inputTree[d]->SetBranchAddress("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal", &p_HadZH_SIG_ghz1_1_JHUGen_JECNominal);
    inputTree[d]->SetBranchAddress("p_LepZH_SIG_ghz1_1_JHUGen", &p_LepZH_SIG_ghz1_1_JHUGen);
    inputTree[d]->SetBranchAddress("GenHMass", &GenHMass);
    inputTree[d]->SetBranchAddress("GenHPt", &GenHPt);
    inputTree[d]->SetBranchAddress("GenHRapidity", &GenHRapidity);
    inputTree[d]->SetBranchAddress("GenLep1Pt", &GenLep1Pt);
    inputTree[d]->SetBranchAddress("GenLep1Eta", &GenLep1Eta);
    inputTree[d]->SetBranchAddress("GenLep1Phi", &GenLep1Phi);
    inputTree[d]->SetBranchAddress("GenLep1Id", &GenLep1Id);
    inputTree[d]->SetBranchAddress("GenLep2Pt", &GenLep2Pt);
    inputTree[d]->SetBranchAddress("GenLep2Eta", &GenLep2Eta);
    inputTree[d]->SetBranchAddress("GenLep2Phi", &GenLep2Phi);
    inputTree[d]->SetBranchAddress("GenLep2Id", &GenLep2Id);
    inputTree[d]->SetBranchAddress("GenLep3Pt", &GenLep3Pt);
    inputTree[d]->SetBranchAddress("GenLep3Eta", &GenLep3Eta);
    inputTree[d]->SetBranchAddress("GenLep3Phi", &GenLep3Phi);
    inputTree[d]->SetBranchAddress("GenLep3Id", &GenLep3Id);
    inputTree[d]->SetBranchAddress("GenLep4Pt", &GenLep4Pt);
    inputTree[d]->SetBranchAddress("GenLep4Eta", &GenLep4Eta);
    inputTree[d]->SetBranchAddress("GenLep4Phi", &GenLep4Phi);
    inputTree[d]->SetBranchAddress("GenLep4Id", &GenLep4Id);
    inputTree[d]->SetBranchAddress("GenAssocLep1Pt", &GenAssocLep1Pt);
    inputTree[d]->SetBranchAddress("GenAssocLep1Eta", &GenAssocLep1Eta);
    inputTree[d]->SetBranchAddress("GenAssocLep1Phi", &GenAssocLep1Phi);
    inputTree[d]->SetBranchAddress("GenAssocLep1Id", &GenAssocLep1Id);
    inputTree[d]->SetBranchAddress("GenAssocLep2Pt", &GenAssocLep2Pt);
    inputTree[d]->SetBranchAddress("GenAssocLep2Eta", &GenAssocLep2Eta);
    inputTree[d]->SetBranchAddress("GenAssocLep2Phi", &GenAssocLep2Phi);
    inputTree[d]->SetBranchAddress("GenAssocLep2Id", &GenAssocLep2Id);

    preselExp[currentProcess] += lumi * 1000 * xsec ;

    Long64_t entries = inputTree[d]->GetEntries();

    Float_t qgIsNormal = 0.;
    Float_t qgIsDefault = 0.;

    Float_t Vbf2jets = 0.;
    Float_t VbfLostJet = 0.;

    Float_t jetBtag = 0.;
    Float_t jetNonBtag = 0.;

    for (Long64_t z=0; z<entries; ++z){

      if(DEBUG && z>1000) break;

      printStatus(z,20000,entries,"entries");

      inputTree[d]->GetEntry(z);
      //cout<<"debug flag entry "<<z<<endl;

      Float_t kfactor = 1.;
      if(APPLYKFACTORS){
	if(currentProcess==qqZZ)      kfactor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M;
	else if(currentProcess==ggZZ) kfactor = KFactor_QCD_ggZZ_Nominal;
      }

      Double_t eventWeight = partialSampleWeight[d] * xsec * kfactor * (USEPUWEIGHT ? overallEventWeight : genHEPMCweight*dataMCWeight) ;

      //Bool_t FullSel70 = ZZsel>=90;
      //Bool_t FullSel100 = ZZsel>=100;
      Bool_t SR = ZZsel>=90;
      Bool_t passTrigger = test_bit_16(trigWord,0);
      Bool_t passTriggerNo1E = test_bit_16(trigWord,8);

      Short_t GenHLepId[4] = {GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id};
      Float_t GenHLepPt[4] = {GenLep1Pt,GenLep2Pt,GenLep3Pt,GenLep4Pt};
      Float_t GenHLepEta[4] = {GenLep1Eta,GenLep2Eta,GenLep3Eta,GenLep4Eta};
      Float_t GenHLepPhi[4] = {GenLep1Phi,GenLep2Phi,GenLep3Phi,GenLep4Phi};
      Short_t GenAssocLepId[2] = {GenAssocLep1Id,GenAssocLep2Id};
      Float_t GenAssocLepPt[2] = {GenAssocLep1Pt,GenAssocLep2Pt};
      Float_t GenAssocLepEta[2] = {GenAssocLep1Eta,GenAssocLep2Eta};
      Float_t GenAssocLepPhi[2] = {GenAssocLep1Phi,GenAssocLep2Phi};


      // ---------------------- reco decay channel --------------------------

      Int_t nCandEle = 0;
      Int_t nCandMu = 0;
      for(Int_t iCandLep=0; iCandLep<4; iCandLep++){
	if(abs(CandLepId->at(iCandLep))==11) nCandEle++;
	if(abs(CandLepId->at(iCandLep))==13) nCandMu++;
      }
      Int_t rc = 4;
      if     (nCandEle==0 && nCandMu==4) rc = 0;
      else if(nCandEle==4 && nCandMu==0) rc = 1;
      else if(nCandEle==2 && nCandMu==2) rc = 2;
      //if(rc==4) cout<<"error in channel definition"<<endl;
      Int_t rc2 = 5;
      if(rc==0||rc==1||rc==2) rc2 = 3;
      Int_t rc3 = 6;


      // ---------------------- Gen lepton counting --------------------------
 
      Int_t nGenHLep = 0;
      Int_t nGenHLepInEtaAcc   = 0;
      Int_t nGenHLepInPtAcc    = 0;
      Int_t nGenHLepInEtaPtAcc = 0;
      Int_t nGenHEle = 0;
      Int_t nGenHMu  = 0;
      Int_t nGenHTau = 0;
      Int_t nGenHLEP = 0;
      Int_t nGenAssocLep = 0;
      Int_t nGenAssocLepInEtaAcc   = 0;
      Int_t nGenAssocLepInPtAcc    = 0;
      Int_t nGenAssocLepInEtaPtAcc = 0;
      Int_t nGenAssocEle = 0;
      Int_t nGenAssocMu  = 0;
      Int_t nGenAssocTau = 0;
      Int_t nGenAssocLEP = 0;
      Int_t nGenLEPPlus  = 0;
      Int_t nGenLEPMinus = 0;
      Bool_t  GenHLepIsInEtaAcc  [4];
      Bool_t  GenHLepIsInPtAcc   [4];
      Bool_t  GenHLepIsInEtaPtAcc[4];
      Bool_t  GenAssocLepIsInEtaAcc  [4];
      Bool_t  GenAssocLepIsInPtAcc   [4];
      Bool_t  GenAssocLepIsInEtaPtAcc[4];
      for(Int_t iGenHLep=0; iGenHLep<4; iGenHLep++){
	if(abs(GenHLepId[iGenHLep])==11){
	  nGenHEle++; nGenHLep++; nGenHLEP++;
	  if(GenHLepId[iGenHLep]>0) nGenLEPPlus++; else nGenLEPMinus++;
	  GenHLepIsInEtaAcc[iGenHLep] = fabs(GenHLepEta[iGenHLep])<2.5 ;
	  GenHLepIsInPtAcc [iGenHLep] = GenHLepPt[iGenHLep]>7. ;
	  GenHLepIsInEtaPtAcc[iGenHLep] = GenHLepIsInEtaAcc[iGenHLep] && GenHLepIsInPtAcc[iGenHLep] ;
	  if(GenHLepIsInEtaAcc  [iGenHLep]) nGenHLepInEtaAcc  ++;
	  if(GenHLepIsInPtAcc   [iGenHLep]) nGenHLepInPtAcc   ++;
	  if(GenHLepIsInEtaPtAcc[iGenHLep]) nGenHLepInEtaPtAcc++;
	}else if(abs(GenHLepId[iGenHLep])==13){
	  nGenHMu ++; nGenHLep++; nGenHLEP++;
	  if(GenHLepId[iGenHLep]>0) nGenLEPPlus++; else nGenLEPMinus++;
	  GenHLepIsInEtaAcc[iGenHLep] = fabs(GenHLepEta[iGenHLep])<2.4 ;
	  GenHLepIsInPtAcc [iGenHLep] = GenHLepPt[iGenHLep]>5. ;
	  GenHLepIsInEtaPtAcc[iGenHLep] = GenHLepIsInEtaAcc[iGenHLep] && GenHLepIsInPtAcc[iGenHLep] ;
	  if(GenHLepIsInEtaAcc  [iGenHLep]) nGenHLepInEtaAcc  ++;
	  if(GenHLepIsInPtAcc   [iGenHLep]) nGenHLepInPtAcc   ++;
	  if(GenHLepIsInEtaPtAcc[iGenHLep]) nGenHLepInEtaPtAcc++;
	}else if(abs(GenHLepId[iGenHLep])==15){
	  nGenHTau++; nGenHLEP++;
	  if(GenHLepId[iGenHLep]>0) nGenLEPPlus++; else nGenLEPMinus++;
	}
      }
      for(Int_t iGenAssocLep=0; iGenAssocLep<2; iGenAssocLep++){
	if(abs(GenAssocLepId[iGenAssocLep])==11){
	  nGenAssocEle++; nGenAssocLep++; nGenAssocLEP++;
	  if(GenAssocLepId[iGenAssocLep]>0) nGenLEPPlus++; else nGenLEPMinus++;
	  GenAssocLepIsInEtaAcc[iGenAssocLep] = fabs(GenAssocLepEta[iGenAssocLep])<2.5 ;
	  GenAssocLepIsInPtAcc [iGenAssocLep] = GenAssocLepPt[iGenAssocLep]>7. ;
	  GenAssocLepIsInEtaPtAcc[iGenAssocLep] = GenAssocLepIsInEtaAcc[iGenAssocLep] && GenAssocLepIsInPtAcc[iGenAssocLep] ;
	  if(GenAssocLepIsInEtaAcc  [iGenAssocLep]) nGenAssocLepInEtaAcc  ++;
	  if(GenAssocLepIsInPtAcc   [iGenAssocLep]) nGenAssocLepInPtAcc   ++;
	  if(GenAssocLepIsInEtaPtAcc[iGenAssocLep]) nGenAssocLepInEtaPtAcc++;
	}else if(abs(GenAssocLepId[iGenAssocLep])==13){
	  nGenAssocMu ++; nGenAssocLep++; nGenAssocLEP++;
	  if(GenAssocLepId[iGenAssocLep]>0) nGenLEPPlus++; else nGenLEPMinus++;
	  GenAssocLepIsInEtaAcc[iGenAssocLep] = fabs(GenAssocLepEta[iGenAssocLep])<2.4 ;
	  GenAssocLepIsInPtAcc [iGenAssocLep] = GenAssocLepPt[iGenAssocLep]>5. ;
	  GenAssocLepIsInEtaPtAcc[iGenAssocLep] = GenAssocLepIsInEtaAcc[iGenAssocLep] && GenAssocLepIsInPtAcc[iGenAssocLep] ;
	  if(GenAssocLepIsInEtaAcc  [iGenAssocLep]) nGenAssocLepInEtaAcc  ++;
	  if(GenAssocLepIsInPtAcc   [iGenAssocLep]) nGenAssocLepInPtAcc   ++;
	  if(GenAssocLepIsInEtaPtAcc[iGenAssocLep]) nGenAssocLepInEtaPtAcc++;
	}else if(abs(GenAssocLepId[iGenAssocLep])==15){
	  nGenAssocTau++; nGenAssocLEP++;
	  if(GenAssocLepId[iGenAssocLep]>0) nGenLEPPlus++; else nGenLEPMinus++;
	}
      }
      Int_t nGenLep = nGenHLep + nGenAssocLep;
      Int_t nGenLepInEtaAcc   = nGenHLepInEtaAcc   + nGenAssocLepInEtaAcc  ;
      Int_t nGenLepInPtAcc    = nGenHLepInPtAcc    + nGenAssocLepInPtAcc   ;
      Int_t nGenLepInEtaPtAcc = nGenHLepInEtaPtAcc + nGenAssocLepInEtaPtAcc;
      Int_t nGenEle = nGenHEle + nGenAssocEle;
      Int_t nGenMu  = nGenHMu  + nGenAssocMu ;
      Int_t nGenTau = nGenHTau + nGenAssocTau;
      Int_t nGenLEP = nGenHLEP + nGenAssocLEP;
 
      Int_t currentGenHDecay = -1;
      if     (nGenHMu ==4               ) currentGenHDecay = 0;
      else if(nGenHEle==4               ) currentGenHDecay = 1;
      else if(nGenHEle==2 && nGenHMu ==2) currentGenHDecay = 2;
      else if(nGenHTau==4               ) currentGenHDecay = 3;
      else if(nGenHEle==2 && nGenHTau==2) currentGenHDecay = 4;
      else if(nGenHMu ==2 && nGenHTau==2) currentGenHDecay = 5;
      else                                currentGenHDecay = 6;
      Int_t gc1 = currentGenHDecay;
      Int_t gc2 = 9;
      if(0<=currentGenHDecay && currentGenHDecay<=2) gc2 = 7;
      if(3<=currentGenHDecay && currentGenHDecay<=5) gc2 = 8;
      Int_t gc3 = 10;


      // ---------------------- Control counters --------------------------

      nbStored[currentProcess][gc1]++;
      nbStored[currentProcess][gc2]++;
      nbStored[currentProcess][gc3]++;
      yieldStored[currentProcess][gc1] += eventWeight;
      yieldStored[currentProcess][gc2] += eventWeight;
      yieldStored[currentProcess][gc3] += eventWeight;
      if(nGenHLepInEtaAcc==4){
	nbHLepsAreInEtaAcc[currentProcess][gc1]++;
	nbHLepsAreInEtaAcc[currentProcess][gc2]++;
	nbHLepsAreInEtaAcc[currentProcess][gc3]++;
	yieldHLepsAreInEtaAcc[currentProcess][gc1] += eventWeight;
	yieldHLepsAreInEtaAcc[currentProcess][gc2] += eventWeight;
	yieldHLepsAreInEtaAcc[currentProcess][gc3] += eventWeight;
	Float_t SortedGenHLepPt[4];
	for(int sl=0; sl<4; sl++) SortedGenHLepPt[sl] = GenHLepPt[sl];
	TriBulle(SortedGenHLepPt,4,false);
	for(int sl=0; sl<4; sl++){
	  h1GenptHLepsAreInEtaAcc[currentProcess][gc1][sl]->Fill(SortedGenHLepPt[sl],eventWeight);
	  h1GenptHLepsAreInEtaAcc[currentProcess][gc2][sl]->Fill(SortedGenHLepPt[sl],eventWeight);
	  h1GenptHLepsAreInEtaAcc[currentProcess][gc3][sl]->Fill(SortedGenHLepPt[sl],eventWeight);
	}
      }
      if(nGenHLepInPtAcc==4){
	Float_t SortedGenHLepAbseta[4];
	for(int sl=0; sl<4; sl++) SortedGenHLepAbseta[sl] = fabs(GenHLepEta[sl]);
	TriBulle(SortedGenHLepAbseta,4,false);
	for(int sl=0; sl<4; sl++){
	  h1GenetaHLepsAreInPtAcc[currentProcess][gc1][sl]->Fill(SortedGenHLepAbseta[sl],eventWeight);
	  h1GenetaHLepsAreInPtAcc[currentProcess][gc2][sl]->Fill(SortedGenHLepAbseta[sl],eventWeight);
	  h1GenetaHLepsAreInPtAcc[currentProcess][gc3][sl]->Fill(SortedGenHLepAbseta[sl],eventWeight);
	}
      }
      if(nGenHLepInEtaPtAcc==4){
	nbHLepsAreInEtaPtAcc[currentProcess][gc1]++;
	nbHLepsAreInEtaPtAcc[currentProcess][gc2]++;
	nbHLepsAreInEtaPtAcc[currentProcess][gc3]++;
	yieldHLepsAreInEtaPtAcc[currentProcess][gc1] += eventWeight;
	yieldHLepsAreInEtaPtAcc[currentProcess][gc2] += eventWeight;
	yieldHLepsAreInEtaPtAcc[currentProcess][gc3] += eventWeight;
	h1NRecoLepHLepsAreInEtaPtAcc[currentProcess][gc1]->Fill(NRecoMu+NRecoEle,eventWeight);
	h1NRecoLepHLepsAreInEtaPtAcc[currentProcess][gc2]->Fill(NRecoMu+NRecoEle,eventWeight);
	h1NRecoLepHLepsAreInEtaPtAcc[currentProcess][gc3]->Fill(NRecoMu+NRecoEle,eventWeight);
	h1NRecoMuoHLepsAreInEtaPtAcc[currentProcess][gc1]->Fill(NRecoMu,eventWeight);
	h1NRecoMuoHLepsAreInEtaPtAcc[currentProcess][gc2]->Fill(NRecoMu,eventWeight);
	h1NRecoMuoHLepsAreInEtaPtAcc[currentProcess][gc3]->Fill(NRecoMu,eventWeight);
	h1NRecoEleHLepsAreInEtaPtAcc[currentProcess][gc1]->Fill(NRecoEle,eventWeight);
	h1NRecoEleHLepsAreInEtaPtAcc[currentProcess][gc2]->Fill(NRecoEle,eventWeight);
	h1NRecoEleHLepsAreInEtaPtAcc[currentProcess][gc3]->Fill(NRecoEle,eventWeight);
      }
      if(nGenLep>=4){
	nb4GenLeps[currentProcess][gc1]++;
	nb4GenLeps[currentProcess][gc2]++;
	nb4GenLeps[currentProcess][gc3]++;
	yield4GenLeps[currentProcess][gc1] += eventWeight;
	yield4GenLeps[currentProcess][gc2] += eventWeight;
	yield4GenLeps[currentProcess][gc3] += eventWeight;
      }
      if(NRecoMu+NRecoEle>=4){
	nb4RecoLeps[currentProcess][gc1]++;
	nb4RecoLeps[currentProcess][gc2]++;
	nb4RecoLeps[currentProcess][gc3]++;
	yield4RecoLeps[currentProcess][gc1] += eventWeight;
	yield4RecoLeps[currentProcess][gc2] += eventWeight;
	yield4RecoLeps[currentProcess][gc3] += eventWeight;
      }
      if(SR){
	nbWithBCInSRG[currentProcess][gc1]++;
	nbWithBCInSRG[currentProcess][gc2]++;
	nbWithBCInSRG[currentProcess][gc3]++;
	yieldWithBCInSRG[currentProcess][gc1] += eventWeight;
	yieldWithBCInSRG[currentProcess][gc2] += eventWeight;
	yieldWithBCInSRG[currentProcess][gc3] += eventWeight;
      }
      if(passTrigger){
	nbPassTriggerG[currentProcess][gc1]++;
	nbPassTriggerG[currentProcess][gc2]++;
	nbPassTriggerG[currentProcess][gc3]++;
	yieldPassTriggerG[currentProcess][gc1] += eventWeight;
	yieldPassTriggerG[currentProcess][gc2] += eventWeight;
	yieldPassTriggerG[currentProcess][gc3] += eventWeight;
      }
      if(passTriggerNo1E){
	nbPassTriggerGNo1E[currentProcess][gc1]++;
	nbPassTriggerGNo1E[currentProcess][gc2]++;
	nbPassTriggerGNo1E[currentProcess][gc3]++;
	yieldPassTriggerGNo1E[currentProcess][gc1] += eventWeight;
	yieldPassTriggerGNo1E[currentProcess][gc2] += eventWeight;
	yieldPassTriggerGNo1E[currentProcess][gc3] += eventWeight;
      }
      if(SR && passTrigger){
	nbPassTriggerGWithBCInSRG[currentProcess][gc1]++;
	nbPassTriggerGWithBCInSRG[currentProcess][gc2]++;
	nbPassTriggerGWithBCInSRG[currentProcess][gc3]++;
	yieldPassTriggerGWithBCInSRG[currentProcess][gc1] += eventWeight;
	yieldPassTriggerGWithBCInSRG[currentProcess][gc2] += eventWeight;
	yieldPassTriggerGWithBCInSRG[currentProcess][gc3] += eventWeight;
	Float_t SortedCandLepPt[4];
	for(int sl=0; sl<4; sl++) SortedCandLepPt[sl] = CandLepPt->at(sl);
	TriBulle(SortedCandLepPt,4,false);
	for(int sl=0; sl<4; sl++){
	  h1RecoptWithBCInSRGPassTriggerG[currentProcess][gc1][sl]->Fill(SortedCandLepPt[sl],eventWeight);
	  h1RecoptWithBCInSRGPassTriggerG[currentProcess][gc2][sl]->Fill(SortedCandLepPt[sl],eventWeight);
	  h1RecoptWithBCInSRGPassTriggerG[currentProcess][gc3][sl]->Fill(SortedCandLepPt[sl],eventWeight);
	}
	Float_t SortedCandLepAbseta[4];
	for(int sl=0; sl<4; sl++) SortedCandLepAbseta[sl] = fabs(CandLepEta->at(sl));
	TriBulle(SortedCandLepAbseta,4,false);
	for(int sl=0; sl<4; sl++){
	  h1RecoetaWithBCInSRGPassTriggerG[currentProcess][gc1][sl]->Fill(SortedCandLepAbseta[sl],eventWeight);
	  h1RecoetaWithBCInSRGPassTriggerG[currentProcess][gc2][sl]->Fill(SortedCandLepAbseta[sl],eventWeight);
	  h1RecoetaWithBCInSRGPassTriggerG[currentProcess][gc3][sl]->Fill(SortedCandLepAbseta[sl],eventWeight);
	}
      }
      if(SR && passTriggerNo1E){
	nbPassTriggerGNo1EWithBCInSRG[currentProcess][gc1]++;
	nbPassTriggerGNo1EWithBCInSRG[currentProcess][gc2]++;
	nbPassTriggerGNo1EWithBCInSRG[currentProcess][gc3]++;
	yieldPassTriggerGNo1EWithBCInSRG[currentProcess][gc1] += eventWeight;
	yieldPassTriggerGNo1EWithBCInSRG[currentProcess][gc2] += eventWeight;
	yieldPassTriggerGNo1EWithBCInSRG[currentProcess][gc3] += eventWeight;
      }
      if(nGenHLepInEtaPtAcc==4 && NRecoMu+NRecoEle>=4){
	nbHLepsAreInEtaPtAcc4RecoLeps[currentProcess][gc1]++;
	nbHLepsAreInEtaPtAcc4RecoLeps[currentProcess][gc2]++;
	nbHLepsAreInEtaPtAcc4RecoLeps[currentProcess][gc3]++;
	yieldHLepsAreInEtaPtAcc4RecoLeps[currentProcess][gc1] += eventWeight;
	yieldHLepsAreInEtaPtAcc4RecoLeps[currentProcess][gc2] += eventWeight;
	yieldHLepsAreInEtaPtAcc4RecoLeps[currentProcess][gc3] += eventWeight;
      }
      if(nGenHLepInEtaPtAcc==4 && SR){
	nbHLepsAreInEtaPtAccWithBCInSRG[currentProcess][gc1]++;
	nbHLepsAreInEtaPtAccWithBCInSRG[currentProcess][gc2]++;
	nbHLepsAreInEtaPtAccWithBCInSRG[currentProcess][gc3]++;
	yieldHLepsAreInEtaPtAccWithBCInSRG[currentProcess][gc1] += eventWeight;
	yieldHLepsAreInEtaPtAccWithBCInSRG[currentProcess][gc2] += eventWeight;
	yieldHLepsAreInEtaPtAccWithBCInSRG[currentProcess][gc3] += eventWeight;
      }
      if(nGenHLepInEtaPtAcc==4 && passTrigger){
	nbHLepsAreInEtaPtAccPassTriggerG[currentProcess][gc1]++;
	nbHLepsAreInEtaPtAccPassTriggerG[currentProcess][gc2]++;
	nbHLepsAreInEtaPtAccPassTriggerG[currentProcess][gc3]++;
	yieldHLepsAreInEtaPtAccPassTriggerG[currentProcess][gc1] += eventWeight;
	yieldHLepsAreInEtaPtAccPassTriggerG[currentProcess][gc2] += eventWeight;
	yieldHLepsAreInEtaPtAccPassTriggerG[currentProcess][gc3] += eventWeight;
      }
      if(nGenHLepInEtaPtAcc==4 && passTriggerNo1E){
	nbHLepsAreInEtaPtAccPassTriggerGNo1E[currentProcess][gc1]++;
	nbHLepsAreInEtaPtAccPassTriggerGNo1E[currentProcess][gc2]++;
	nbHLepsAreInEtaPtAccPassTriggerGNo1E[currentProcess][gc3]++;
	yieldHLepsAreInEtaPtAccPassTriggerGNo1E[currentProcess][gc1] += eventWeight;
	yieldHLepsAreInEtaPtAccPassTriggerGNo1E[currentProcess][gc2] += eventWeight;
	yieldHLepsAreInEtaPtAccPassTriggerGNo1E[currentProcess][gc3] += eventWeight;
      }
      if(nGenHLepInEtaPtAcc==4 && SR && passTrigger){
	nbHLepsAreInEtaPtAccPassTriggerGWithBCInSRG[currentProcess][gc1]++;
	nbHLepsAreInEtaPtAccPassTriggerGWithBCInSRG[currentProcess][gc2]++;
	nbHLepsAreInEtaPtAccPassTriggerGWithBCInSRG[currentProcess][gc3]++;
	yieldHLepsAreInEtaPtAccPassTriggerGWithBCInSRG[currentProcess][gc1] += eventWeight;
	yieldHLepsAreInEtaPtAccPassTriggerGWithBCInSRG[currentProcess][gc2] += eventWeight;
	yieldHLepsAreInEtaPtAccPassTriggerGWithBCInSRG[currentProcess][gc3] += eventWeight;
      }
      if(nGenHLepInEtaPtAcc==4 && SR && passTriggerNo1E){
	nbHLepsAreInEtaPtAccPassTriggerGNo1EWithBCInSRG[currentProcess][gc1]++;
	nbHLepsAreInEtaPtAccPassTriggerGNo1EWithBCInSRG[currentProcess][gc2]++;
	nbHLepsAreInEtaPtAccPassTriggerGNo1EWithBCInSRG[currentProcess][gc3]++;
	yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCInSRG[currentProcess][gc1] += eventWeight;
	yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCInSRG[currentProcess][gc2] += eventWeight;
	yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCInSRG[currentProcess][gc3] += eventWeight;
      }
      h1GenHPt[currentProcess]->Fill(GenHPt,eventWeight);
      h1GenHRapidity[currentProcess]->Fill(GenHRapidity,eventWeight);
      h2GenHRapidityVsPt[currentProcess]->Fill(GenHPt,GenHRapidity,eventWeight);
      //       nbStored[currentProcess][rc]++;
      //       nbStored[currentProcess][rc2]++;
      //       nbStored[currentProcess][rc3]++;
      //       yieldStored[currentProcess][rc] += eventWeight;
      //       yieldStored[currentProcess][rc2] += eventWeight;
      //       yieldStored[currentProcess][rc3] += eventWeight;


      // ---------------------- Cuts related to gen leptons --------------------------

      /*
      if(nGenLep<4){ // gen. 2L2tau events
	if( isSignal[currentProcess]) continue;
	if(!isSignal[currentProcess]) continue;
      } 

      if(
	 (!isSignal[currentProcess] && !( nGenHLep==4 && nGenAssocLep==0 ) ) ||
	 ( isSignal[currentProcess] && !( nGenHLep==4 || (nGenHLep==2&&nGenAssocLep==2&&(currentProcess==ZH||currentProcess==ttH)) ) )
	 ){
	cout<<"ERROR with nb. of gen leptons (dataset "<<datasets[d]<<", event "<<z<<", nGenHLep="<<nGenHLep<<", nGenAssocLep="<<nGenAssocLep<<")"<<endl;
	continue;
      }
      //*/

      ///////////////////////////
      if(excludeH2l2X)
	if(nGenHLep!=4 && isSignal[currentProcess]) continue;
      ///////////////////////////


      // ---------------------- Gen associated decay --------------------------

      Int_t currentAssocDecay = -1;
      if(currentProcess==WH){
	if(nGenHLep==4 && nGenAssocLep==0) currentAssocDecay = 0;
	else if(nGenHLep==4 && nGenAssocLep==1) currentAssocDecay = 1;
      }else if(currentProcess==ZH){
	if(nGenHLep==4 && nGenAssocLep==0) currentAssocDecay = 2;
	else if(nGenHLep==4 && nGenAssocLep==2) currentAssocDecay = 3;
	else if(nGenHLep==2 && nGenAssocLep==2) currentAssocDecay = 4;
      }else if(currentProcess==ttH){
	if(nGenHLep==4 && nGenAssocLep==0) currentAssocDecay = 5;
	else if(nGenHLep==4 && nGenAssocLep==1) currentAssocDecay = 6;
	else if(nGenHLep==4 && nGenAssocLep==2) currentAssocDecay = 7;
	else if(nGenHLep==2 && nGenAssocLep==2) currentAssocDecay = 8;
      }
      if((currentProcess==WH || currentProcess==ZH || currentProcess==ttH) && currentAssocDecay == -1)
	cout<<"ERROR with assoc decays: dataset="<<datasets[d]<<", nGenHLep="<<nGenHLep<<", nGenAssocLep="<<nGenAssocLep<<endl;


      // ---------------------- Successive selection steps --------------------------

      nbWithBC[currentProcess][rc]++;
      nbWithBC[currentProcess][rc2]++;
      yieldWithBC[currentProcess][rc] += eventWeight;
      yieldWithBC[currentProcess][rc2] += eventWeight;

      if( !(ZZsel>=90) ) continue;

      if(ZZMass<m4lMin || ZZMass>m4lMax) continue;

      nbWithBCInSR[currentProcess][rc]++;
      nbWithBCInSR[currentProcess][rc2]++;
      yieldWithBCInSR[currentProcess][rc] += eventWeight;
      yieldWithBCInSR[currentProcess][rc2] += eventWeight;


      // ---------------------- Gen to Good lepton matching --------------------------

      Int_t nRecoLepMatchedToGenHLep[4] = {0,0,0,0};
      Int_t nCandLepMatchedToGenHLep[4] = {0,0,0,0};
      Int_t nRecoLepMatchedToGenAssocLep[2] = {0,0};
      Int_t nCandLepMatchedToGenAssocLep[2] = {0,0};
      Int_t nGenLepMatchedToCandLep[4] = {0,0,0,0};
      Int_t nGenLepMatchedToRecoLep[7] = {0,0,0,0,0,0,0};
      Int_t nGenHLepMatchedToZ1Lep[4] = {0,0};
      Int_t nGenHLepMatchedToZ2Lep[4] = {0,0};
      for(Int_t iGenHLep=0; iGenHLep<4; iGenHLep++){
	if(abs(GenHLepId[iGenHLep])==11 || abs(GenHLepId[iGenHLep])==13){
	  for(Int_t iCandLep=0; iCandLep<4; iCandLep++){
	    if(deltaR(GenHLepEta[iGenHLep],GenHLepPhi[iGenHLep],CandLepEta->at(iCandLep),CandLepPhi->at(iCandLep)) < 0.1){
	      nRecoLepMatchedToGenHLep[iGenHLep]++;
	      nCandLepMatchedToGenHLep[iGenHLep]++;
	      nGenLepMatchedToRecoLep[iCandLep]++;
	      nGenLepMatchedToCandLep[iCandLep]++;
	      if(iCandLep<2){
		nGenHLepMatchedToZ1Lep[iCandLep]++;
	      }else{
		nGenHLepMatchedToZ2Lep[iCandLep-2]++;
	      }
	    }
	  }
	  for(Int_t iExtraLep=0; iExtraLep<nExtraLep; iExtraLep++){
	    if(deltaR(GenHLepEta[iGenHLep],GenHLepPhi[iGenHLep],ExtraLepEta->at(iExtraLep),ExtraLepPhi->at(iExtraLep)) < 0.1){
	      nRecoLepMatchedToGenHLep[iGenHLep]++;
	      nGenLepMatchedToRecoLep[4+iExtraLep]++;
	    }
	  }
	}
      }
      for(Int_t iGenAssocLep=0; iGenAssocLep<2; iGenAssocLep++){
	if(abs(GenAssocLepId[iGenAssocLep])==11 || abs(GenAssocLepId[iGenAssocLep])==13){
	  for(Int_t iCandLep=0; iCandLep<4; iCandLep++){
	    if(deltaR(GenAssocLepEta[iGenAssocLep],GenAssocLepPhi[iGenAssocLep],CandLepEta->at(iCandLep),CandLepPhi->at(iCandLep)) < 0.1){
	      nRecoLepMatchedToGenAssocLep[iGenAssocLep]++;
	      nCandLepMatchedToGenAssocLep[iGenAssocLep]++;
	      nGenLepMatchedToRecoLep[iCandLep]++;
	      nGenLepMatchedToCandLep[iCandLep]++;
	    }
	  }
	  for(Int_t iExtraLep=0; iExtraLep<nExtraLep; iExtraLep++){
	    if(deltaR(GenAssocLepEta[iGenAssocLep],GenAssocLepPhi[iGenAssocLep],ExtraLepEta->at(iExtraLep),ExtraLepPhi->at(iExtraLep)) < 0.1){
	      nRecoLepMatchedToGenAssocLep[iGenAssocLep]++;
	      nGenLepMatchedToRecoLep[4+iExtraLep]++;
	    }
	  }
	}
      }


      // ---------------------- Exploit matching information --------------------------
	    
      Bool_t foundMatchingAmbiguity = false;
      for(Int_t iGenHLep=0; iGenHLep<4; iGenHLep++) if(nRecoLepMatchedToGenHLep[iGenHLep]>1){ foundMatchingAmbiguity = true; break; }
      for(Int_t iGenAssocLep=0; iGenAssocLep<2; iGenAssocLep++) if(nRecoLepMatchedToGenAssocLep[iGenAssocLep]>1){ foundMatchingAmbiguity = true; break; }
      for(Int_t iRecoLep=0; iRecoLep<4+nExtraLep; iRecoLep++) if(nGenLepMatchedToRecoLep[iRecoLep]>1){ foundMatchingAmbiguity = true; break; }

      Int_t nOnes = 0;
      for(Int_t iGenHLep=0; iGenHLep<4; iGenHLep++) if(nRecoLepMatchedToGenHLep[iGenHLep]==1) nOnes++;

      Int_t nOnesHLeps = 0;
      for(Int_t iGenHLep=0; iGenHLep<4; iGenHLep++) if(nCandLepMatchedToGenHLep[iGenHLep]==1) nOnesHLeps++;
      Int_t nOnesAssocLeps = 0;
      for(Int_t iGenAssocLep=0; iGenAssocLep<2; iGenAssocLep++) if(nCandLepMatchedToGenAssocLep[iGenAssocLep]==1) nOnesAssocLeps++;
      Int_t currentMatchHLepsStatus = -1;
      if(foundMatchingAmbiguity){
	currentMatchHLepsStatus = 5;
      }else{
	if(nOnesHLeps==4) currentMatchHLepsStatus = 0;
	if(nOnesHLeps==3) currentMatchHLepsStatus = 1;
	if(nOnesHLeps==2) currentMatchHLepsStatus = 2;
	if(nOnesHLeps==1) currentMatchHLepsStatus = 3;
	if(nOnesHLeps==0) currentMatchHLepsStatus = 4;
      }
      Int_t currentMatchAllLepsStatus = -1;
      if(foundMatchingAmbiguity){
	currentMatchAllLepsStatus = 4;
      }else{
	if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchAllLepsStatus = 0;
	if(nOnesHLeps==3 && nOnesAssocLeps==1) currentMatchAllLepsStatus = 1;
	if(nOnesHLeps==2 && nOnesAssocLeps==2) currentMatchAllLepsStatus = 2;
	if(nOnesHLeps+nOnesAssocLeps<4) currentMatchAllLepsStatus = 3;
      }
      Int_t currentMatchWHStatus = -1;
      Int_t currentMatchZHStatus = -1;
      Int_t currentMatchttHStatus = -1;
      if(currentProcess==WH){
	if(foundMatchingAmbiguity){
	  currentMatchWHStatus = 4;
	}else{
	  if(nOnesHLeps+nOnesAssocLeps<4){
	    currentMatchWHStatus = 3;
	  }else{
	    if(nGenHLep==4 && nGenAssocLep==0){
	      if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchWHStatus = 0;
	      else cout<<"error nOnes"<<endl;
	    }else if(nGenHLep==4 && nGenAssocLep==1){
	      if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchWHStatus = 1;
	      else if(nOnesHLeps==3 && nOnesAssocLeps==1) currentMatchWHStatus = 2;
	      else cout<<"error nOnes"<<endl;
	    }else{
	      cout<<"error nGen"<<endl;
	    }
	  }
	} 
      }
      if(currentProcess==ZH){
	if(foundMatchingAmbiguity){
	  currentMatchZHStatus = 6;
	}else{
	  if(nOnesHLeps+nOnesAssocLeps<4){
	    currentMatchZHStatus = 5;
	  }else{
	    if(nGenHLep==4 && nGenAssocLep==0){
	      if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchZHStatus = 0;
	      else cout<<"error nOnes"<<endl;
	    }else if(nGenHLep==4 && nGenAssocLep==2){
	      if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchZHStatus = 1;
	      else if(nOnesHLeps==3 && nOnesAssocLeps==1) currentMatchZHStatus = 2;
	      else if(nOnesHLeps==2 && nOnesAssocLeps==2) currentMatchZHStatus = 3;
	      else cout<<"error nOnes"<<endl;
	    }else if(nGenHLep==2 && nGenAssocLep==2){
	      if(nOnesHLeps==2 && nOnesAssocLeps==2) currentMatchZHStatus = 4;
	      else cout<<"error nOnes"<<endl;
	    }else{
	      cout<<"error nGen"<<endl;
	    }
	  }
	} 
      }
      if(currentProcess==ttH){
	if(foundMatchingAmbiguity){
	  currentMatchttHStatus = 8;
	}else{
	  if(nOnesHLeps+nOnesAssocLeps<4){
	    currentMatchttHStatus = 7;
	  }else{
	    if(nGenHLep==4 && nGenAssocLep==0){
	      if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchttHStatus = 0;
	      else cout<<"error nOnes"<<endl;
	    }else if(nGenHLep==4 && nGenAssocLep==1){
	      if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchttHStatus = 1;
	      else if(nOnesHLeps==3 && nOnesAssocLeps==1) currentMatchttHStatus = 2;
	      else cout<<"error nOnes"<<endl;
	    }else if(nGenHLep==4 && nGenAssocLep==2){
	      if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchttHStatus = 3;
	      else if(nOnesHLeps==3 && nOnesAssocLeps==1) currentMatchttHStatus = 4;
	      else if(nOnesHLeps==2 && nOnesAssocLeps==2) currentMatchttHStatus = 5;
	      else cout<<"error nOnes"<<endl;
	    }else if(nGenHLep==2 && nGenAssocLep==2){
	      if(nOnesHLeps==2 && nOnesAssocLeps==2) currentMatchttHStatus = 6;
	      else cout<<"error nOnes"<<endl;
	    }else{
	      cout<<"error nGen"<<endl;
	    }
	  }
	} 
      } 

      Int_t currentZ1MatchStatus = -1;
      if(nGenHLepMatchedToZ1Lep[0]>1 || nGenHLepMatchedToZ1Lep[1]>1) currentZ1MatchStatus = 3;
      else currentZ1MatchStatus = 2 - (nGenHLepMatchedToZ1Lep[0] + nGenHLepMatchedToZ1Lep[1]);
      Int_t currentZ2MatchStatus = -1;
      if(nGenHLepMatchedToZ2Lep[0]>1 || nGenHLepMatchedToZ2Lep[1]>1) currentZ2MatchStatus = 3;
      else currentZ2MatchStatus = 2 - (nGenHLepMatchedToZ2Lep[0] + nGenHLepMatchedToZ2Lep[1]);


      // ---------------------- Cuts impacting counters/histograms --------------------------

      Bool_t Exactly4GoodLeps = nExtraLep==0 ;
      Bool_t AtLeast5GoodLeps = nExtraLep>=1 ;
      Bool_t Exactly5GoodLeps = nExtraLep==1 ;
      Bool_t Exactly6GoodLeps = nExtraLep==2 ;
      Bool_t HLepsAreInEtaPtAcc    = nGenHLepInEtaPtAcc==4 ;
      Bool_t HLepsAreGood     = !foundMatchingAmbiguity && nOnes==4 ;

      if(requireExactly4GoodLeps && !Exactly4GoodLeps) continue;
      if(requireAtLeast5GoodLeps && !AtLeast5GoodLeps) continue;
      if(requireExactly5GoodLeps && !Exactly5GoodLeps) continue;
      if(requireExactly6GoodLeps && !Exactly6GoodLeps) continue;
      if(requireHLepsAreInEtaPtAcc    && !HLepsAreInEtaPtAcc   ) continue;
      if(requireHLepsAreGood     && !HLepsAreGood    ) continue;


      // ---------------------- Fill histograms and increment counters --------------------------

      if(Exactly4GoodLeps){ nbWithBCInSRExactly4GoodLeps[currentProcess][rc]++; nbWithBCInSRExactly4GoodLeps[currentProcess][rc2]++; }
      if(AtLeast5GoodLeps){ nbWithBCInSRAtLeast5GoodLeps[currentProcess][rc]++; nbWithBCInSRAtLeast5GoodLeps[currentProcess][rc2]++; }
      if(Exactly5GoodLeps){ nbWithBCInSRExactly5GoodLeps[currentProcess][rc]++; nbWithBCInSRExactly5GoodLeps[currentProcess][rc2]++; }
      if(Exactly6GoodLeps){ nbWithBCInSRExactly6GoodLeps[currentProcess][rc]++; nbWithBCInSRExactly6GoodLeps[currentProcess][rc2]++; }
      if(HLepsAreInEtaPtAcc   ){ nbWithBCInSRHLepsAreInEtaPtAcc   [currentProcess][rc]++; nbWithBCInSRHLepsAreInEtaPtAcc   [currentProcess][rc2]++; }
      if(HLepsAreGood    ){ nbWithBCInSRHLepsAreGood    [currentProcess][rc]++; nbWithBCInSRHLepsAreGood    [currentProcess][rc2]++; }
      Float_t cQGUnoff = 1./1000.;
      Int_t nJetsBTagged_Redone = 0;
      Int_t nJetsBTaggedLoose_Redone = 0;
      for(int j=0; j<nJets; j++){
	if(JetBTagger->at(j)>CSVv2M) nJetsBTagged_Redone++;
	if(JetBTagger->at(j)>CSVv2L) nJetsBTaggedLoose_Redone++;
	jetPt[j] = JetPt->at(j);
	jetEta[j] = JetEta->at(j);
	jetPhi[j] = JetPhi->at(j);
	jetMass[j] = JetMass->at(j);
	if(j<2) { if(JetIsBtagged->at(j)) jetBtag += eventWeight; else jetNonBtag += eventWeight; }
	jetQGLikelihoodRaw[j] = JetQGLikelihood->at(j);
	jetQGLikelihood[j] = JetQGLikelihood->at(j);
	if(jetQGLikelihood[j]<0. && j<2){
	  //cout<<" jetQGLikelihood was "<<jetQGLikelihood[j];
	  TRandom3 rand;
	  rand.SetSeed(abs(static_cast<int>(sin(jetPhi[j])*100000)));
	  jetQGLikelihood[j] = rand.Uniform(); 
	  //cout<<", now "<<jetQGLikelihood[j]<<endl;
	  qgIsDefault += eventWeight;
	}else{
	  qgIsNormal += eventWeight;
	}
	jetPQuark[j] = 0.;//JetPQuark->at(j);//1000000. * JetPQuark->at(j);
	jetPGluon[j] = 0.;//JetPGluon->at(j);//1000. * JetPGluon->at(j);
	jetPgOverPq[j] = OFFICIALQGTAGGER ? (1./jetQGLikelihood[j] - 1.) : (cQGUnoff*jetPGluon[j]/jetPQuark[j]) ;
      }
      if(!BTAGGINGSF && nJetsBTagged!=nJetsBTagged_Redone) cout<<"ERROR : inconsistency in number of b-tagged jets"<<endl;
      if(nJets>=2) Vbf2jets += eventWeight; else VbfLostJet += eventWeight; 
      Float_t VbfLostJet = 0.;
      Float_t cVBF2j = getDVBF2jetsConstant(ZZMass);
      Float_t cVBF1j = getDVBF1jetConstant(ZZMass);
      Float_t cWH = getDWHhConstant(ZZMass);
      Float_t cZH = getDZHhConstant(ZZMass);
      Float_t cWHlept = 100000000000.;
      Float_t cZHlept = 100000000.;
      Float_t pwhlept = p_LepWH_SIG_ghw1_1_JHUGen/cWHlept;
      Float_t pzhlept = p_LepZH_SIG_ghz1_1_JHUGen/cZHlept;
      Float_t KD = 1./(1.+ getDbkgkinConstant(Z1Flav*Z2Flav,ZZMass)*p_QQB_BKG_MCFM/p_GG_SIG_ghg2_1_ghz1_1_JHUGen );//1./(1.+ p_QQB_BKG_MCFM/p_GG_SIG_ghg2_1_ghz1_1_JHUGen );
      Float_t d2jVbfHjj = 1./(1.+ (cVBF2j*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal) / p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal );
      Float_t d1jVbfHj = 1./(1.+ (cVBF1j*p_JQCD_SIG_ghg2_1_JHUGen_JECNominal) / (p_JVBF_SIG_ghv1_1_JHUGen_JECNominal*pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal) );
      Float_t d2jWHHadrHjj = 1./(1.+ (cWH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal) / p_HadWH_SIG_ghw1_1_JHUGen_JECNominal );
      Float_t d2jZHHadrHjj = 1./(1.+ (cZH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal) / p_HadZH_SIG_ghz1_1_JHUGen_JECNominal );
      Float_t pqj1Pqj2 = (nJets>=2) ? 1000*1000*jetPQuark[0]*jetPQuark[1]/(cQGUnoff*cQGUnoff) : -2 ;
      Float_t pgj1Pgj2 = (nJets>=2) ? 1000*1000*jetPGluon[0]*jetPGluon[1] : -2 ;
      Float_t d2jqg = (nJets>=2) ? 1./(1.+ jetPgOverPq[0]*jetPgOverPq[1] ) : -2 ;
      Float_t dqgj1Dqgj2 = (nJets>=2) ? 1./(1.+jetPgOverPq[0]) * 1./(1.+jetPgOverPq[1]) : -2 ;
      Float_t pq = (nJets>=1) ? 1000*jetPQuark[0]/cQGUnoff : -2 ;
      Float_t pg = (nJets>=1) ? 1000*jetPGluon[0] : -2 ;
      Float_t d1jqg = (nJets>=1) ? 1./(1.+jetPgOverPq[0]) : -2 ;
      Float_t d2jMelaQGVbfHjj = (nJets>=2) ? 1./(1.+ cVBF2j*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal * jetPgOverPq[0]*jetPgOverPq[1] ) : -2 ;
      Float_t d2jMelaD2jQGVbfHjj = d2jVbfHjj*d2jqg;
      Float_t d1jMelaQGVbfHj = (nJets>=1) ? 1./(1.+ (cVBF1j*p_JQCD_SIG_ghg2_1_JHUGen_JECNominal)/(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal*pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal) * jetPgOverPq[0] ) : -2 ;
      Float_t d1jMelaD1jQGVbfHj = d1jVbfHj*d1jqg;
      Float_t d2jMelaQGWHHadrHjj = (nJets>=2) ? 1./(1.+ cWH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadWH_SIG_ghw1_1_JHUGen_JECNominal * jetPgOverPq[0]*jetPgOverPq[1] ) : -2 ;
      Float_t d2jMelaD2jQGWHHadrHjj = d2jWHHadrHjj*d2jqg;
      Float_t d2jMelaQGZHHadrHjj = (nJets>=2) ? 1./(1.+ cZH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadZH_SIG_ghz1_1_JHUGen_JECNominal * jetPgOverPq[0]*jetPgOverPq[1] ) : -2 ;
      Float_t d2jMelaD2jQGZHHadrHjj = d2jZHHadrHjj*d2jqg;
      Float_t d2jMelaExpQGVbfHjj = (nJets>=2) ? 1./(1.+ cVBF2j*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal * TMath::Exp(jetPgOverPq[0]*jetPgOverPq[1]) ) : -2 ;
      Float_t d2jMelaSqQGVbfHjj = (nJets>=2) ? 1./(1.+ cVBF2j*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal * TMath::Power(jetPgOverPq[0]*jetPgOverPq[1],2.) ) : -2 ;
      Float_t d2jMelaSqrtQGVbfHjj = (nJets>=2) ? 1./(1.+ cVBF2j*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal * TMath::Power(jetPgOverPq[0]*jetPgOverPq[1],1/2.) ) : -2 ;
      Float_t d2jMelaCbrtQGVbfHjj = (nJets>=2) ? 1./(1.+ cVBF2j*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal * TMath::Power(jetPgOverPq[0]*jetPgOverPq[1],1/3.) ) : -2 ;
      Float_t d2jMelaQrrtQGVbfHjj = (nJets>=2) ? 1./(1.+ cVBF2j*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal * TMath::Power(jetPgOverPq[0]*jetPgOverPq[1],1/4.) ) : -2 ;
      Float_t d2jMelaQnrtQGVbfHjj = (nJets>=2) ? 1./(1.+ cVBF2j*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal * TMath::Power(jetPgOverPq[0]*jetPgOverPq[1],1/5.) ) : -2 ;
      Float_t d1jMelaSqrtQGVbfHj = (nJets>=1) ? 1./(1.+ (cVBF1j*p_JQCD_SIG_ghg2_1_JHUGen_JECNominal)/(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal*pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal) * TMath::Power(jetPgOverPq[0],1/2.) ) : -2 ;
      Float_t d1jMelaCbrtQGVbfHj = (nJets>=1) ? 1./(1.+ (cVBF1j*p_JQCD_SIG_ghg2_1_JHUGen_JECNominal)/(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal*pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal) * TMath::Power(jetPgOverPq[0],1/3.) ) : -2 ;
      Float_t d1jMelaQrrtQGVbfHj = (nJets>=1) ? 1./(1.+ (cVBF1j*p_JQCD_SIG_ghg2_1_JHUGen_JECNominal)/(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal*pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal) * TMath::Power(jetPgOverPq[0],1/4.) ) : -2 ;
      Float_t d1jMelaQnrtQGVbfHj = (nJets>=1) ? 1./(1.+ (cVBF1j*p_JQCD_SIG_ghg2_1_JHUGen_JECNominal)/(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal*pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal) * TMath::Power(jetPgOverPq[0],1/5.) ) : -2 ;
      Float_t d2jMelaSqrtQGWHHadrHjj = (nJets>=2) ? 1./(1.+ cWH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadWH_SIG_ghw1_1_JHUGen_JECNominal * TMath::Power(jetPgOverPq[0]*jetPgOverPq[1],1/2.) ) : -2 ;
      Float_t d2jMelaCbrtQGWHHadrHjj = (nJets>=2) ? 1./(1.+ cWH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadWH_SIG_ghw1_1_JHUGen_JECNominal * TMath::Power(jetPgOverPq[0]*jetPgOverPq[1],1/3.) ) : -2 ;
      Float_t d2jMelaQrrtQGWHHadrHjj = (nJets>=2) ? 1./(1.+ cWH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadWH_SIG_ghw1_1_JHUGen_JECNominal * TMath::Power(jetPgOverPq[0]*jetPgOverPq[1],1/4.) ) : -2 ;
      Float_t d2jMelaQnrtQGWHHadrHjj = (nJets>=2) ? 1./(1.+ cWH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadWH_SIG_ghw1_1_JHUGen_JECNominal * TMath::Power(jetPgOverPq[0]*jetPgOverPq[1],1/5.) ) : -2 ;
      Float_t d2jMelaSqrtQGZHHadrHjj = (nJets>=2) ? 1./(1.+ cZH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadZH_SIG_ghz1_1_JHUGen_JECNominal * TMath::Power(jetPgOverPq[0]*jetPgOverPq[1],1/2.) ) : -2 ;
      Float_t d2jMelaCbrtQGZHHadrHjj = (nJets>=2) ? 1./(1.+ cZH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadZH_SIG_ghz1_1_JHUGen_JECNominal * TMath::Power(jetPgOverPq[0]*jetPgOverPq[1],1/3.) ) : -2 ;
      Float_t d2jMelaQrrtQGZHHadrHjj = (nJets>=2) ? 1./(1.+ cZH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadZH_SIG_ghz1_1_JHUGen_JECNominal * TMath::Power(jetPgOverPq[0]*jetPgOverPq[1],1/4.) ) : -2 ;
      Float_t d2jMelaQnrtQGZHHadrHjj = (nJets>=2) ? 1./(1.+ cZH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadZH_SIG_ghz1_1_JHUGen_JECNominal * TMath::Power(jetPgOverPq[0]*jetPgOverPq[1],1/5.) ) : -2 ;
      Float_t varVal[nVariables] = {
	ZZMass,
	ZZMass,
	Z1Mass,
	Z2Mass,
	KD,
	DiJetFisher,
	p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal/cVBF2j,
	p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
	p_JVBF_SIG_ghv1_1_JHUGen_JECNominal*pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal/cVBF1j,
	p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
	p_HadWH_SIG_ghw1_1_JHUGen_JECNominal/cWH,
	p_HadZH_SIG_ghz1_1_JHUGen_JECNominal/cZH,
	pwhlept,
	pzhlept,
	d2jVbfHjj,
	d1jVbfHj,
	d2jWHHadrHjj,
	d2jZHHadrHjj,
	pq,
	pg,
	pqj1Pqj2,
	pgj1Pgj2,
	d2jqg,
	dqgj1Dqgj2,
	pq,
	pg,
	pqj1Pqj2,
	pgj1Pgj2,
	d2jqg,
	pq,
	pg,
	d1jqg,
	d2jMelaQGVbfHjj,
	d2jMelaD2jQGVbfHjj,
	d1jMelaQGVbfHj,
	d1jMelaD1jQGVbfHj,
	d2jMelaQGWHHadrHjj,
	d2jMelaD2jQGWHHadrHjj,
	d2jMelaQGZHHadrHjj,
	d2jMelaD2jQGZHHadrHjj,
	p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal / p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
	jetPgOverPq[0]*jetPgOverPq[1],
	p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal / p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal * jetPgOverPq[0]*jetPgOverPq[1],
	d2jMelaExpQGVbfHjj,
	d2jMelaSqQGVbfHjj,
	d2jMelaSqrtQGVbfHjj,
	d2jMelaCbrtQGVbfHjj,
	d2jMelaQrrtQGVbfHjj,
	d2jMelaQnrtQGVbfHjj,
	d1jMelaSqrtQGVbfHj,
	d1jMelaCbrtQGVbfHj,
	d1jMelaQrrtQGVbfHj,
	d1jMelaQnrtQGVbfHj,
	d2jMelaSqrtQGWHHadrHjj,
	d2jMelaCbrtQGWHHadrHjj,
	d2jMelaQrrtQGWHHadrHjj,
	d2jMelaQnrtQGWHHadrHjj,
	d2jMelaSqrtQGZHHadrHjj,
	d2jMelaCbrtQGZHHadrHjj,
	d2jMelaQrrtQGZHHadrHjj,
	d2jMelaQnrtQGZHHadrHjj,
	ZZPt,
	(Float_t)nGenLep,
	(Float_t)nGenLepInEtaPtAcc,
	(Float_t)(nGenLep-nGenLepInEtaPtAcc),
	(Float_t)(nGenHLep-nGenHLepInEtaPtAcc),
	(Float_t)(nGenAssocLep-nGenAssocLepInEtaPtAcc),
	(Float_t)(nGenLep-(4+nExtraLep)),
	(Float_t)(nGenLepInEtaPtAcc-(4+nExtraLep)),
	(Float_t)nExtraLep,
	(Float_t)nExtraZ,
	(Float_t)nJets,
	(Float_t)nJetsBTagged,
	PFMET,
      };
      Bool_t varPassCut[nVariables] = {
	1,1,1,1,1,nJets>=2,nJets>=2,nJets>=2,nJets==1,nJets==1,nJets>=2,nJets>=2,nExtraLep>=1,nExtraZ>=1,nJets>=2,nJets==1,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2&&d2jVbfHjj>0.5,nJets>=2&&d2jVbfHjj>0.5,nJets>=2&&d2jVbfHjj>0.5,nJets>=2&&d2jVbfHjj>0.5,nJets>=2&&d2jVbfHjj>0.5,nJets==1,nJets==1,nJets==1,nJets>=2,nJets>=2,nJets==1,nJets==1,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets==1,nJets==1,nJets==1,nJets==1,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,1,1,1,1,1,1,1,1,1,1,1,1,1,
      };
      for(int v=0; v<nVariables; v++){
	if(!varPassCut[v]) continue;
	hBCInSR[v][currentProcess][rc]->Fill(varVal[v],eventWeight);
	hBCInSR[v][currentProcess][rc2]->Fill(varVal[v],eventWeight);
      }

      Bool_t varPairPassCut[nVariables] = {
	1,1,nJets>=2,
      };
      for(int v2=0; v2<nVarPairs; v2++){
	if(!varPairPassCut[v2]) continue;
	h2DBCInSR[v2][currentProcess][rc]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
	h2DBCInSR[v2][currentProcess][rc2]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
	if(!varPairRef[v2][2]) continue;
	h2DBCInSRDecays[v2][currentProcess][0][rc]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
	h2DBCInSRDecays[v2][currentProcess][0][rc2]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
	if(nGenHLep==4){
	  h2DBCInSRDecays[v2][currentProcess][1][rc]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
	  h2DBCInSRDecays[v2][currentProcess][1][rc2]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
	  if(currentMatchHLepsStatus==0){
	    h2DBCInSRDecays[v2][currentProcess][3][rc]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
	    h2DBCInSRDecays[v2][currentProcess][3][rc2]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
	  }else if(currentMatchHLepsStatus<5){
	    h2DBCInSRDecays[v2][currentProcess][4][rc]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
	    h2DBCInSRDecays[v2][currentProcess][4][rc2]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
	  }
	}else{
	  h2DBCInSRDecays[v2][currentProcess][2][rc]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
	  h2DBCInSRDecays[v2][currentProcess][2][rc2]->Fill(varVal[varPairRef[v2][0]],varVal[varPairRef[v2][1]],eventWeight);
	}
      }

      nbWithBCInSRMatchHLeps[currentMatchHLepsStatus][currentProcess][rc]++;
      nbWithBCInSRMatchHLeps[currentMatchHLepsStatus][currentProcess][rc2]++;
      for(int v=0; v<nVariables; v++){
	hBCInSRMatchHLeps[v][currentMatchHLepsStatus][currentProcess][rc]->Fill(varVal[v],eventWeight);
	hBCInSRMatchHLeps[v][currentMatchHLepsStatus][currentProcess][rc2]->Fill(varVal[v],eventWeight);
      }

      nbWithBCInSRMatchAllLeps[currentMatchAllLepsStatus][currentProcess][rc]++;
      nbWithBCInSRMatchAllLeps[currentMatchAllLepsStatus][currentProcess][rc2]++;
      for(int v=0; v<nVariables; v++){
	hBCInSRMatchAllLeps[v][currentMatchAllLepsStatus][currentProcess][rc]->Fill(varVal[v],eventWeight);
	hBCInSRMatchAllLeps[v][currentMatchAllLepsStatus][currentProcess][rc2]->Fill(varVal[v],eventWeight);
      }
	    
      if(currentProcess==WH){
	nbWithBCInSRMatchWH[currentMatchWHStatus][rc]++;
	nbWithBCInSRMatchWH[currentMatchWHStatus][rc2]++;
	for(int v=0; v<nVariables; v++){
	  hBCInSRMatchWH[v][currentMatchWHStatus][rc]->Fill(varVal[v],eventWeight);
	  hBCInSRMatchWH[v][currentMatchWHStatus][rc2]->Fill(varVal[v],eventWeight);
	}
      }
      if(currentProcess==ZH){
	nbWithBCInSRMatchZH[currentMatchZHStatus][rc]++;
	nbWithBCInSRMatchZH[currentMatchZHStatus][rc2]++;
	for(int v=0; v<nVariables; v++){
	  hBCInSRMatchZH[v][currentMatchZHStatus][rc]->Fill(varVal[v],eventWeight);
	  hBCInSRMatchZH[v][currentMatchZHStatus][rc2]->Fill(varVal[v],eventWeight);
	}
      }
      if(currentProcess==ttH){
	nbWithBCInSRMatchttH[currentMatchttHStatus][rc]++;
	nbWithBCInSRMatchttH[currentMatchttHStatus][rc2]++;
	for(int v=0; v<nVariables; v++){
	  hBCInSRMatchttH[v][currentMatchttHStatus][rc]->Fill(varVal[v],eventWeight);
	  hBCInSRMatchttH[v][currentMatchttHStatus][rc2]->Fill(varVal[v],eventWeight);
	}
      }

      if(currentMatchAllLepsStatus==0){
	nbWithBCInSRAll4LepRight[currentProcess][rc]++;
	nbWithBCInSRAll4LepRight[currentProcess][rc2]++;
      }

      if(currentProcess==WH){
	Int_t currentWDecay = -1;
	if(nGenHLep==4 && nGenAssocLep==0) currentWDecay = 0;
	else if(nGenHLep==4 && nGenAssocLep==1) currentWDecay = 1;
	else cout<<"error"<<endl;
	nbTotalWH[currentWDecay]++;
	if(nGenHLepInEtaPtAcc==4) nbHLepsAreInEtaPtAccWH[currentWDecay]++;
	if(HLepsAreGood) nbHLepsAreGoodWH[currentWDecay]++;
	if(currentMatchAllLepsStatus==0) nbAll4LepRightWH[currentWDecay]++;
	nbZ1DaughtersFromHWH[currentWDecay][currentZ1MatchStatus]++;
	nbZ2DaughtersFromHWH[currentWDecay][currentZ2MatchStatus]++;
      }
      if(currentProcess==ZH){
	Int_t currentZDecay = -1;
	if(nGenHLep==4 && nGenAssocLep==0) currentZDecay = 0;
	else if(nGenHLep==4 && nGenAssocLep==2) currentZDecay = 1;
	else if(nGenHLep==2 && nGenAssocLep==2) currentZDecay = 2;
	else cout<<"error"<<endl;
	nbTotalZH[currentZDecay]++;
	if(nGenHLepInEtaPtAcc==4) nbHLepsAreInEtaPtAccZH[currentZDecay]++;
	if(HLepsAreGood) nbHLepsAreGoodZH[currentZDecay]++;
	if(currentMatchAllLepsStatus==0) nbAll4LepRightZH[currentZDecay]++;
	nbZ1DaughtersFromHZH[currentZDecay][currentZ1MatchStatus]++;
	nbZ2DaughtersFromHZH[currentZDecay][currentZ2MatchStatus]++;
      }
      if(currentProcess==ttH){
	Int_t currentttDecay = -1;
	if(nGenHLep==4 && nGenAssocLep==0) currentttDecay = 0;
	else if(nGenHLep==4 && nGenAssocLep==1) currentttDecay = 1;
	else if(nGenHLep==4 && nGenAssocLep==2) currentttDecay = 2;
	else if(nGenHLep==2 && nGenAssocLep==2) currentttDecay = 3;
	else cout<<"error"<<endl;
	nbTotalttH[currentttDecay]++;
	if(nGenHLepInEtaPtAcc==4) nbHLepsAreInEtaPtAccttH[currentttDecay]++;
	if(HLepsAreGood) nbHLepsAreGoodttH[currentttDecay]++;
	if(currentMatchAllLepsStatus==0) nbAll4LepRightttH[currentttDecay]++;
	nbZ1DaughtersFromHttH[currentttDecay][currentZ1MatchStatus]++;
	nbZ2DaughtersFromHttH[currentttDecay][currentZ2MatchStatus]++;
      }

      Bool_t rocPassCut[nROCs] = {
	nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets==1,nJets==1,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets==1,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets==1,nJets==1,nJets==1,nJets==1,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets==1,nJets==1,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets==1,nJets==1,nJets==1,nJets==1,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nJets>=2,
      };
      for(int ro=0; ro<nROCs; ro++){
	if(!rocPassCut[ro]) continue;
	if(currentProcess==rocRef[ro][1])
	  hRocSgnl[ro]->Fill(varVal[rocRef[ro][0]],eventWeight);
	else if(currentProcess==rocRef[ro][2]) 
	  hRocBkgd[ro]->Fill(varVal[rocRef[ro][0]],eventWeight);
      }
      Bool_t evwPassCut[nEVWs] = {
	nJets>=2,nJets==1,nJets>=2,nJets>=2,nJets>=2,nJets==1,nJets>=2,nJets>=2,nJets>=2,nJets>=2,nExtraLep>=1,nExtraZ>=1,
      };
      for(int e=0; e<nEVWs; e++){
	if(!evwPassCut[e]) continue;
	  hEvw[e][currentProcess]->Fill(varVal[evwRef[e][0]],eventWeight);
      }


      // ---------------------- Categorization stuff --------------------------

      /* // sanity check for conformity of discriminant values (OK)
      Float_t d2jVbfHjj_ext = DVBF2j_ME(p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
      Float_t d1jVbfHj_ext = DVBF1j_ME(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);;
      Float_t d2jWHHadrHjj_ext = DWHh_ME(p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
      Float_t d2jZHHadrHjj_ext = DZHh_ME(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
      Float_t d2jMelaCbrtQGVbfHjj_ext = DVBF2j_ME_QG(p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
      Float_t d1jMelaCbrtQGVbfHj_ext = DVBF1j_ME_QG(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
      Float_t d2jMelaCbrtQGWHHadrHjj_ext = DWHh_ME_QG(p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
      Float_t d2jMelaCbrtQGZHHadrHjj_ext = DZHh_ME_QG(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
      int nDigits = 5; //number of significant digits up to which we check equality
      if(nJets>=2 && rtd(d2jVbfHjj,nDigits)!=rtd(d2jVbfHjj_ext,nDigits)) cout<<"error: d2jVbfHjj is "<<d2jVbfHjj<<" from plotProcesses.C and "<<d2jVbfHjj_ext<<" from Discriminants.cc"<<endl;
      if(nJets==1 && rtd(d1jVbfHj,nDigits)!=rtd(d1jVbfHj_ext,nDigits)) cout<<"error: d1jVbfHj is "<<d1jVbfHj<<" from plotProcesses.C and "<<d1jVbfHj_ext<<" from Discriminants.cc"<<endl;
      if(nJets>=2 && rtd(d2jWHHadrHjj,nDigits)!=rtd(d2jWHHadrHjj_ext,nDigits)) cout<<"error: d2jWHHadrHjj is "<<d2jWHHadrHjj<<" from plotProcesses.C and "<<d2jWHHadrHjj_ext<<" from Discriminants.cc"<<endl;
      if(nJets>=2 && rtd(d2jZHHadrHjj,nDigits)!=rtd(d2jZHHadrHjj_ext,nDigits)) cout<<"error: d2jZHHadrHjj is "<<d2jZHHadrHjj<<" from plotProcesses.C and "<<d2jZHHadrHjj_ext<<" from Discriminants.cc"<<endl;
      if(nJets>=2 && rtd(d2jMelaCbrtQGVbfHjj,nDigits)!=rtd(d2jMelaCbrtQGVbfHjj_ext,nDigits)) cout<<"error: d2jMelaCbrtQGVbfHjj is "<<d2jMelaCbrtQGVbfHjj<<" from plotProcesses.C and "<<d2jMelaCbrtQGVbfHjj_ext<<" from Discriminants.cc"<<endl;
      if(nJets==1 && rtd(d1jMelaCbrtQGVbfHj,nDigits)!=rtd(d1jMelaCbrtQGVbfHj_ext,nDigits)) cout<<"error: d1jMelaCbrtQGVbfHj is "<<d1jMelaCbrtQGVbfHj<<" from plotProcesses.C and "<<d1jMelaCbrtQGVbfHj_ext<<" from Discriminants.cc"<<endl;
      if(nJets>=2 && rtd(d2jMelaCbrtQGWHHadrHjj,nDigits)!=rtd(d2jMelaCbrtQGWHHadrHjj_ext,nDigits)) cout<<"error: d2jMelaCbrtQGWHHadrHjj is "<<d2jMelaCbrtQGWHHadrHjj<<" from plotProcesses.C and "<<d2jMelaCbrtQGWHHadrHjj_ext<<" from Discriminants.cc"<<endl;
      if(nJets>=2 && rtd(d2jMelaCbrtQGZHHadrHjj,nDigits)!=rtd(d2jMelaCbrtQGZHHadrHjj_ext,nDigits)) cout<<"error: d2jMelaCbrtQGZHHadrHjj is "<<d2jMelaCbrtQGZHHadrHjj<<" from plotProcesses.C and "<<d2jMelaCbrtQGZHHadrHjj_ext<<" from Discriminants.cc"<<endl;
      //*/

      Bool_t Ichep16TagVBF2j  = 0.;
      Bool_t Ichep16TagVBF1j  = 0.;
      Bool_t Ichep16TagWHHadr = 0.;
      Bool_t Ichep16TagZHHadr = 0.;
      if(USEQGTAGGING){
	Ichep16TagVBF2j  = d2jMelaCbrtQGVbfHjj > 0.391;
	Ichep16TagVBF1j  = d1jMelaCbrtQGVbfHj > 0.72;
	Ichep16TagWHHadr = d2jMelaCbrtQGWHHadrHjj > 0.973;
	Ichep16TagZHHadr = d2jMelaCbrtQGZHHadrHjj > 0.996;
      }else{
	Ichep16TagVBF2j  = d2jVbfHjj > 1.043-460./(ZZMass+634.);
        Ichep16TagVBF1j  = d1jVbfHj > 0.699;
        Ichep16TagWHHadr = d2jWHHadrHjj > 0.959;
        Ichep16TagZHHadr = d2jZHHadrHjj > 0.9946;
      }

      Bool_t Mor17TagVBF2j  = 0.;
      Bool_t Mor17TagVBF1j  = 0.;
      Bool_t Mor17TagWHHadr = 0.;
      Bool_t Mor17TagZHHadr = 0.;
      if(USEQGTAGGING){
	Mor17TagVBF2j  = d2jMelaCbrtQGVbfHjj > (USEMASSDEPWPVBF2JMELAQG ? mdWPVBF2JMELAQG(ZZMass) : WPVBF2JMELAQG);
	Mor17TagVBF1j  = d1jMelaCbrtQGVbfHj > WPVBF1JMELAQG;
	Mor17TagWHHadr = d2jMelaCbrtQGWHHadrHjj > WPWHHADRMELAQG;
	Mor17TagZHHadr = d2jMelaCbrtQGZHHadrHjj > WPZHHADRMELAQG;
      }else{
	Mor17TagVBF2j  = d2jVbfHjj > (USEMASSDEPWPVBF2JMELA ? mdWPVBF2JMELA(ZZMass) : WPVBF2JMELA);
        Mor17TagVBF1j  = d1jVbfHj > WPVBF1JMELA;
        Mor17TagWHHadr = d2jWHHadrHjj > WPWHHADRMELA;
        Mor17TagZHHadr = d2jZHHadrHjj > WPZHHADRMELA;
      }

      Int_t currentBasket = -1;

      if(BASKETLIST==9 ||BASKETLIST==10){

	int cs = -1;
	int cm = -1;

	int un = 1;
	int vbf1j = 2;
	int vbf2j = 3;
	int vhl = 5;
	int vhh = 4;
	int tth = 6;
	string lab[7] = {"","Unt.","VBF-1j","VBF-2j","VH-had","VH-lep","ttH"};

	if(nJets==0){
	  if(nExtraLep==0)      { cs = 1; cm = un; }
	  else if(nExtraLep==1) { cs = 2; cm = vhl; }
	  else if(nExtraZ>=1)   { cs = 3; cm = vhl; }
	  else if(nExtraLep>=2) { cs = 4; cm = vhl; }
	  else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	}else if(nJets==1){
	  if(nJetsBTagged==0){
	    if(nExtraLep==0 && Ichep16TagVBF1j) 
	                          { cs = 5; cm = vbf1j; }
	    else if(nExtraLep==0) { cs = 6; cm = un; }
	    else if(nExtraLep==1) { cs = 7; cm = vhl; }
	    else if(nExtraZ>=1)   { cs = 8; cm = vhl; }
	    else if(nExtraLep>=2) { cs = 9; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else if(nJetsBTagged==1){
	    if(nExtraLep==0 && Ichep16TagVBF1j)
	                          { cs = 10; cm = vbf1j; }
	    else if(nExtraLep==0) { cs = 11; cm = un; }
	    else if(nExtraLep==1) { cs = 12; cm = tth; }
	    else if(nExtraZ>=1)   { cs = 13; cm = tth; }
	    else if(nExtraLep>=2) { cs = 14; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else cout<<"WARNING : inconsistent nJetsBTagged"<<endl;
	}else if(nJets==2){
	  if(nJetsBTagged==0){
	    if(nExtraLep==0 && Ichep16TagVBF2j)
                                  { cs = 15; cm = vbf2j; }
	    else if(nExtraLep==0 && (Ichep16TagWHHadr||Ichep16TagZHHadr))
	                          { cs = 16; cm = vhh; }
	    else if(nExtraLep==0) { cs = 17; cm = un; }
	    else if(nExtraLep==1) { cs = 18; cm = vhl; }
	    else if(nExtraZ>=1)   { cs = 19; cm = vhl; }
	    else if(nExtraLep>=2) { cs = 20; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else if(nJetsBTagged==1){
	    if(nExtraLep==0 && Ichep16TagVBF2j)
	                          { cs = 21; cm = vbf2j; }
	    else if(nExtraLep==0 && (Ichep16TagWHHadr||Ichep16TagZHHadr))
	                          { cs = 22; cm = vhh; }
	    else if(nExtraLep==0) { cs = 23; cm = un; }
	    else if(nExtraLep==1) { cs = 24; cm = tth; }
	    else if(nExtraZ>=1)   { cs = 25; cm = tth; }
	    else if(nExtraLep>=2) { cs = 26; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else if(nJetsBTagged==2){
	    if(nExtraLep==0 && (Ichep16TagWHHadr||Ichep16TagZHHadr))
                                  { cs = 27; cm = vhh; }
	    else if(nExtraLep==0) { cs = 28; cm = vhh; }
	    else if(nExtraLep==1) { cs = 29; cm = tth; }
	    else if(nExtraZ>=1)   { cs = 30; cm = tth; }
	    else if(nExtraLep>=2) { cs = 31; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else cout<<"WARNING : inconsistent nJetsBTagged"<<endl;
	}else if(nJets==3){
	  if(nJetsBTagged==0){
	    if(nExtraLep==0 && Ichep16TagVBF2j)
	                          { cs = 32; cm = vbf2j; }
	    else if(nExtraLep==0 && (Ichep16TagWHHadr||Ichep16TagZHHadr))
                                  { cs = 33; cm = vhh; }
	    else if(nExtraLep==0) { cs = 34; cm = un; }
	    else if(nExtraLep==1) { cs = 35; cm = vhl; }
	    else if(nExtraZ>=1)   { cs = 36; cm = vhl; }
	    else if(nExtraLep>=2) { cs = 37; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else if(nJetsBTagged==1){
	    if(nExtraLep==0 && Ichep16TagVBF2j)
	                          { cs = 38; cm = vbf2j; }
	    else if(nExtraLep==0 && (Ichep16TagWHHadr||Ichep16TagZHHadr))
	                          { cs = 39; cm = vhh; }
	    else if(nExtraLep==0) { cs = 40; cm = un; }
	    else if(nExtraLep==1) { cs = 41; cm = tth; }
	    else if(nExtraZ>=1)   { cs = 42; cm = tth; }
	    else if(nExtraLep>=2) { cs = 43; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else if(nJetsBTagged>=2){
	    if(nExtraLep==0 && (Ichep16TagWHHadr||Ichep16TagZHHadr))
	                          { cs = 44; cm = vhh; }
	    else if(nExtraLep==0) { cs = 45; cm = vhh; }
	    else if(nExtraLep==1) { cs = 46; cm = tth; }
	    else if(nExtraZ>=1)   { cs = 47; cm = tth; }
	    else if(nExtraLep>=2) { cs = 48; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else cout<<"WARNING : inconsistent nJetsBTagged"<<endl;
	}else if(nJets>=4){
	  if(nJetsBTagged==0){
	    if(nExtraLep==0 && Ichep16TagVBF2j)
	                          { cs = 49; cm = vbf2j; }
	    else if(nExtraLep==0 && (Ichep16TagWHHadr||Ichep16TagZHHadr))
	                          { cs = 50; cm = vhh; }
	    else if(nExtraLep==0) { cs = 51; cm = un; }
	    else if(nExtraLep==1) { cs = 52; cm = tth; }
	    else if(nExtraZ>=1)   { cs = 53; cm = tth; }
	    else if(nExtraLep>=2) { cs = 54; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else if(nJetsBTagged==1){
	    if(nExtraLep==0 && Ichep16TagVBF2j)
	                          { cs = 55; cm = tth; }
	    else if(nExtraLep==0 && (Ichep16TagWHHadr||Ichep16TagZHHadr))
	                          { cs = 56; cm = tth; }
	    else if(nExtraLep==0) { cs = 57; cm = tth; }
	    else if(nExtraLep>=1) { cs = 58; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else if(nJetsBTagged>=2){
	    if(nExtraLep==0 && (Ichep16TagWHHadr||Ichep16TagZHHadr))
	                          { cs = 59; cm = tth; }
	    else if(nExtraLep==0) { cs = 60; cm = tth; }
	    else if(nExtraLep>=1) { cs = 61; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else cout<<"WARNING : inconsistent nJetsBTagged"<<endl;
	}else{
	  cout<<"WARNING : inconsistent nJets"<<endl;
	}

	if(BASKETLIST==9) currentBasket = cs;
	if(BASKETLIST==10) currentBasket = cm;
	if(BASKETLIST==9) labelMerge[cs+1] = lab[cm];

	if(BASKETLIST==10 && CHECKFROMCATEGORYCC){
	  float currentBasket2 = categoryIchep16(
						 nExtraLep,
						 nExtraZ,
						 nJets, 
						 nJetsBTagged,
						 jetQGLikelihoodRaw,
						 p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
						 p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
						 p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
						 p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
						 pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
						 p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
						 p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
						 jetPhi,
						 ZZMass,
						 USEQGTAGGING
						 )
	    +1 ;
	  float newc2=-1;
	  if(currentBasket2==4) newc2=vhl;
	  if(currentBasket2==5) newc2=vhh;
	  if(newc2>-1) currentBasket2=newc2; 
	  if(currentBasket!=currentBasket2){
	    cout<<"error: categ from plotProcesses.C = "<<currentBasket<<", from Category.cc = "<<currentBasket2<<endl;
	  }
	  //currentBasket = currentBasket2;
	}

      } else if(BASKETLIST==11 || BASKETLIST==12){

	int cs = -1;
	int cm = -1;

	int un = 1;
	int vbf1j = 2;
	int vbf2j = 3;
	int vhl = 5;
	int vhh = 4;
	int tth = 7;
	int met = 6;//5;
	string lab[8] = {"","Unt.","VBF-1j","VBF-2j","VH-had","VH-lep","MET","ttH"};

	float METCUT = 100.;

	if(nJets==0){
	  if(nExtraLep==0 && PFMET>METCUT) { cs = 62; cm = met; }
	  else if(nExtraLep==0) { cs = 1; cm = un; }
	  else if(nExtraLep==1) { cs = 2; cm = vhl; }
	  else if(nExtraZ>=1)   { cs = 3; cm = vhl; }
	  else if(nExtraLep>=2) { cs = 4; cm = vhl; }
	  else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	}else if(nJets==1){
	  if(nJetsBTagged==0){
	    if(nExtraLep==0 && PFMET>METCUT) { cs = 63; cm = met; }
	    else if(nExtraLep==0 && Mor17TagVBF1j) 
	                          { cs = 5; cm = vbf1j; }
	    else if(nExtraLep==0) { cs = 6; cm = un; }
	    else if(nExtraLep==1) { cs = 7; cm = vhl; }
	    else if(nExtraZ>=1)   { cs = 8; cm = vhl; }
	    else if(nExtraLep>=2) { cs = 9; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else if(nJetsBTagged==1){
	    if(nExtraLep==0 && PFMET>METCUT) { cs = 64; cm = met; }
	    else if(nExtraLep==0 && Mor17TagVBF1j)
	                          { cs = 10; cm = vbf1j; }
	    else if(nExtraLep==0) { cs = 11; cm = un; }
	    else if(nExtraLep==1) { cs = 12; cm = tth; }
	    else if(nExtraZ>=1)   { cs = 13; cm = tth; }
	    else if(nExtraLep>=2) { cs = 14; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else cout<<"WARNING : inconsistent nJetsBTagged"<<endl;
	}else if(nJets==2){
	  if(nJetsBTagged==0){
	    if(nExtraLep==0 && Mor17TagVBF2j)
                                  { cs = 15; cm = vbf2j; }
	    else if(nExtraLep==0 && (Mor17TagWHHadr||Mor17TagZHHadr))
	                          { cs = 16; cm = vhh; }
	    else if(nExtraLep==0) { cs = 17; cm = un; }
	    else if(nExtraLep==1) { cs = 18; cm = vhl; }
	    else if(nExtraZ>=1)   { cs = 19; cm = vhl; }
	    else if(nExtraLep>=2) { cs = 20; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else if(nJetsBTagged==1){
	    if(nExtraLep==0 && Mor17TagVBF2j)
	                          { cs = 21; cm = vbf2j; }
	    else if(nExtraLep==0 && (Mor17TagWHHadr||Mor17TagZHHadr))
	                          { cs = 22; cm = vhh; }
	    else if(nExtraLep==0) { cs = 23; cm = un; }
	    else if(nExtraLep==1) { cs = 24; cm = tth; }
	    else if(nExtraZ>=1)   { cs = 25; cm = tth; }
	    else if(nExtraLep>=2) { cs = 26; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else if(nJetsBTagged==2){
	    if(nExtraLep==0 && (Mor17TagWHHadr||Mor17TagZHHadr))
                                  { cs = 27; cm = vhh; }
	    else if(nExtraLep==0) { cs = 28; cm = un; } // changed wrt Mor17
	    else if(nExtraLep==1) { cs = 29; cm = tth; }
	    else if(nExtraZ>=1)   { cs = 30; cm = tth; }
	    else if(nExtraLep>=2) { cs = 31; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else cout<<"WARNING : inconsistent nJetsBTagged"<<endl;
	}else if(nJets==3){
	  if(nJetsBTagged==0){
	    if(nExtraLep==0 && Mor17TagVBF2j)
	                          { cs = 32; cm = vbf2j; }
	    else if(nExtraLep==0 && (Mor17TagWHHadr||Mor17TagZHHadr))
                                  { cs = 33; cm = vhh; }
	    else if(nExtraLep==0) { cs = 34; cm = un; }
	    else if(nExtraLep==1) { cs = 35; cm = vhl; }
	    else if(nExtraZ>=1)   { cs = 36; cm = vhl; }
	    else if(nExtraLep>=2) { cs = 37; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else if(nJetsBTagged==1){
	    if(nExtraLep==0 && Mor17TagVBF2j)
	                          { cs = 38; cm = vbf2j; }
	    else if(nExtraLep==0 && (Mor17TagWHHadr||Mor17TagZHHadr))
	                          { cs = 39; cm = vhh; }
	    else if(nExtraLep==0) { cs = 40; cm = un; }
	    else if(nExtraLep==1) { cs = 41; cm = tth; }
	    else if(nExtraZ>=1)   { cs = 42; cm = tth; }
	    else if(nExtraLep>=2) { cs = 43; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else if(nJetsBTagged>=2){
	    if(nExtraLep==0 && (Mor17TagWHHadr||Mor17TagZHHadr))
	                          { cs = 44; cm = vhh; }
	    else if(nExtraLep==0) { cs = 45; cm = un; } // changed wrt Mor17
	    else if(nExtraLep==1) { cs = 46; cm = tth; }
	    else if(nExtraZ>=1)   { cs = 47; cm = tth; }
	    else if(nExtraLep>=2) { cs = 48; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else cout<<"WARNING : inconsistent nJetsBTagged"<<endl;
	}else if(nJets>=4){
	  if(nJetsBTagged==0){
	    if(nExtraLep==0 && Mor17TagVBF2j)
	                          { cs = 49; cm = vbf2j; }
	    else if(nExtraLep==0 && (Mor17TagWHHadr||Mor17TagZHHadr))
	                          { cs = 50; cm = vhh; }
	    else if(nExtraLep==0) { cs = 51; cm = un; }
	    else if(nExtraLep==1) { cs = 52; cm = tth; }
	    else if(nExtraZ>=1)   { cs = 53; cm = tth; }
	    else if(nExtraLep>=2) { cs = 54; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else if(nJetsBTagged==1){
	    if(nExtraLep==0 && Mor17TagVBF2j)
	                          { cs = 55; cm = tth; }
	    else if(nExtraLep==0 && (Mor17TagWHHadr||Mor17TagZHHadr))
	                          { cs = 56; cm = tth; }
	    else if(nExtraLep==0) { cs = 57; cm = tth; }
	    else if(nExtraLep>=1) { cs = 58; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else if(nJetsBTagged>=2){
	    if(nExtraLep==0 && (Mor17TagWHHadr||Mor17TagZHHadr))
	                          { cs = 59; cm = tth; }
	    else if(nExtraLep==0) { cs = 60; cm = tth; }
	    else if(nExtraLep>=1) { cs = 61; cm = tth; }
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else cout<<"WARNING : inconsistent nJetsBTagged"<<endl;
	}else{
	  cout<<"WARNING : inconsistent nJets"<<endl;
	}

	if(BASKETLIST==11) currentBasket = cs;
	if(BASKETLIST==12) currentBasket = cm;
	if(BASKETLIST==11) labelMerge[cs+1] = lab[cm];

	if(BASKETLIST==12 && CHECKFROMCATEGORYCC){
	  float currentBasket2 = categoryMor17(
					       nExtraLep,
					       nExtraZ,
					       nJets, 
					       nJetsBTagged,
					       jetQGLikelihoodRaw,
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
					       1,//whether to use VHMETTagged
					       USEQGTAGGING
					       )
	    +1 ;
	  float newc2=-1; 
	  if(currentBasket2==4) newc2=vhl; 
	  if(currentBasket2==5) newc2=vhh;
	  if(currentBasket2==6) newc2=tth;
	  if(currentBasket2==7) newc2=met;
	  if(newc2>-1) currentBasket2=newc2; 
	  if(currentBasket!=currentBasket2){
	    cout<<"error: categ from plotProcesses.C = "<<currentBasket<<", from Category.cc = "<<currentBasket2<<endl;
	  }
	  //currentBasket = currentBasket2;
	}

      }

      hBCInSRBaskets[currentProcess][rc]->Fill(currentBasket,eventWeight);
      hBCInSRBaskets[currentProcess][rc]->Fill(0.,eventWeight);
      hBCInSRBaskets[currentProcess][rc2]->Fill(currentBasket,eventWeight);
      hBCInSRBaskets[currentProcess][rc2]->Fill(0.,eventWeight);
      if(currentAssocDecay>=0){
	hBCInSRBasketsAssocDecays[currentAssocDecay][rc]->Fill(currentBasket,eventWeight);
	hBCInSRBasketsAssocDecays[currentAssocDecay][rc]->Fill(0.,eventWeight);
	hBCInSRBasketsAssocDecays[currentAssocDecay][rc2]->Fill(currentBasket,eventWeight);
	hBCInSRBasketsAssocDecays[currentAssocDecay][rc2]->Fill(0.,eventWeight);	      
      }


    } // end for entries

    //cout<<"fraction of 1st or 2nd jets with default value = "<<qgIsDefault/(qgIsDefault+qgIsNormal)<<endl;
    //cout<<"fraction of b-tagged jets = "<<jetBtag/(jetBtag+jetNonBtag)<<endl;
    //cout<<"fraction of events with less than 2 jets = "<<VbfLostJet/(VbfLostJet+Vbf2jets)<<endl;

  } // end for datasets
  



  // ------------------------------------------------------------
  // ------------------------ Printing --------------------------
  // ------------------------------------------------------------


  Int_t widthColumn2 = 7;

  for(int pr=0; pr<nProcesses; pr++){
    if(!isProcessed[pr]) continue;

    for(int rc=0; rc<nRecoChannels; rc++){
      txtOut<<" "<<recoChannels[rc]<<endl;
      txtOut<<"  has BC :     "<<fixWidth(Form("%i",nbWithBC[pr][rc]),widthColumn2,false)<<endl;
      txtOut<<"  in SR :      "<<fixWidth(Form("%i",nbWithBCInSR[pr][rc]),widthColumn2,false)<<endl;
      if(requireExactly4GoodLeps) txtOut<<"  [ From this point, require that there is exactly 4 good leptons ]"<<endl;
      if(requireAtLeast5GoodLeps) txtOut<<"  [ From this point, require that there is at least 5 good leptons ]"<<endl;
      if(requireExactly5GoodLeps) txtOut<<"  [ From this point, require that there is exactly 5 good leptons ]"<<endl;
      if(requireExactly6GoodLeps) txtOut<<"  [ From this point, require that there is exactly 6 good leptons ]"<<endl;
      if(requireHLepsAreInEtaPtAcc   ) txtOut<<"  [ From this point, require that the 4 gen-leptons from the H are in the acceptance ]"<<endl;
      if(requireHLepsAreGood    ) txtOut<<"  [ From this point, require that the 4 gen-leptons from the H are reconstructed as good leptons ]"<<endl;
      txtOut<<"  in SR, there is ==4 / >=5 / ==5 / ==6 good leptons : "<<nbWithBCInSRExactly4GoodLeps[pr][rc]<<" / "<<nbWithBCInSRAtLeast5GoodLeps[pr][rc]<<" / "<<nbWithBCInSRExactly5GoodLeps[pr][rc]<<" / "<<nbWithBCInSRExactly6GoodLeps[pr][rc]<<endl;
      txtOut<<"  in SR, the 4 gen-leptons from the H are in the acceptance : "<<nbWithBCInSRHLepsAreInEtaPtAcc[pr][rc]<<endl;
      txtOut<<"  in SR, the 4 gen-leptons from the H are reconstructed as good leptons : "<<nbWithBCInSRHLepsAreGood[pr][rc]<<endl;
      txtOut<<"  in SR, the 4 gen-leptons from the H are the 4 good leptons of the best candidate : "<<nbWithBCInSRAll4LepRight[pr][rc]<<endl;
    }

    int widthColumn1 = 22;
    int widthOtherColumns = 7;
    string separator = repeat("-",widthColumn1+4*(2+widthOtherColumns));
    if(pr==WH){
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<fixWidth("# of Z1 leptons from H",widthColumn1,true);
      for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(nbZDaughtersFromH[j],widthOtherColumns,false);
      txtOut<<endl;
      for(int i=0; i<nAssocWDecays; i++){
	txtOut<<" "<<fixWidth(WHdecays[i],widthColumn1,true);
	for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(percentage((float)nbZ1DaughtersFromHWH[i][j]/nbTotalWH[i])+" %",widthOtherColumns,false);
	txtOut<<endl;
      }
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<fixWidth("# of Z2 leptons from H",widthColumn1,true);
      for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(nbZDaughtersFromH[j],widthOtherColumns,false);
      txtOut<<endl;
      for(int i=0; i<nAssocWDecays; i++){
	txtOut<<" "<<fixWidth(WHdecays[i],widthColumn1,true);
	for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(percentage((float)nbZ2DaughtersFromHWH[i][j]/nbTotalWH[i])+" %",widthOtherColumns,false);
	txtOut<<endl;
      }
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<"the 4 gen-leptons from the H are in the acceptance :"<<endl;
      for(int i=0; i<nAssocWDecays; i++) txtOut<<" "<<fixWidth(WHdecays[i],widthColumn1,true)<<"   "<<percentage((float)nbHLepsAreInEtaPtAccWH[i]/nbTotalWH[i])+" %"<<endl;
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<"the 4 gen-leptons from the H are reconstructed as good leptons :"<<endl;
      for(int i=0; i<nAssocWDecays; i++) txtOut<<" "<<fixWidth(WHdecays[i],widthColumn1,true)<<"   "<<percentage((float)nbHLepsAreGoodWH[i]/nbTotalWH[i])+" %"<<endl;
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<"the 4 gen-leptons from the H are the 4 good leptons of the best candidate :"<<endl;
      for(int i=0; i<nAssocWDecays; i++) txtOut<<" "<<fixWidth(WHdecays[i],widthColumn1,true)<<"   "<<percentage((float)nbAll4LepRightWH[i]/nbTotalWH[i])+" %"<<endl;
    }
    if(pr==ZH){
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<fixWidth("# of Z1 leptons from H",widthColumn1,true);
      for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(nbZDaughtersFromH[j],widthOtherColumns,false);
      txtOut<<endl;
      for(int i=0; i<nAssocZDecays; i++){
	txtOut<<" "<<fixWidth(ZHdecays[i],widthColumn1,true);
	for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(percentage((float)nbZ1DaughtersFromHZH[i][j]/nbTotalZH[i])+" %",widthOtherColumns,false);
	txtOut<<endl;
      }
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<fixWidth("# of Z2 leptons from H",widthColumn1,true);
      for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(nbZDaughtersFromH[j],widthOtherColumns,false);
      txtOut<<endl;
      for(int i=0; i<nAssocZDecays; i++){
	txtOut<<" "<<fixWidth(ZHdecays[i],widthColumn1,true);
	for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(percentage((float)nbZ2DaughtersFromHZH[i][j]/nbTotalZH[i])+" %",widthOtherColumns,false);
	txtOut<<endl;
      }
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<"the 4 gen-leptons from the H are in the acceptance :"<<endl;
      for(int i=0; i<nAssocZDecays; i++) txtOut<<" "<<fixWidth(ZHdecays[i],widthColumn1,true)<<"   "<<percentage((float)nbHLepsAreInEtaPtAccZH[i]/nbTotalZH[i])+" %"<<endl;
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<"the 4 gen-leptons from the H are reconstructed as good leptons :"<<endl;
      for(int i=0; i<nAssocZDecays; i++) txtOut<<" "<<fixWidth(ZHdecays[i],widthColumn1,true)<<"   "<<percentage((float)nbHLepsAreGoodZH[i]/nbTotalZH[i])+" %"<<endl;
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<"the 4 gen-leptons from the H are the 4 good leptons of the best candidate :"<<endl;
      for(int i=0; i<nAssocZDecays; i++) txtOut<<" "<<fixWidth(ZHdecays[i],widthColumn1,true)<<"   "<<percentage((float)nbAll4LepRightZH[i]/nbTotalZH[i])+" %"<<endl;
    }
    if(pr==ttH){
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<fixWidth("# of Z1 leptons from H",widthColumn1,true);
      for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(nbZDaughtersFromH[j],widthOtherColumns,false);
      txtOut<<endl;
      for(int i=0; i<nAssocttDecays; i++){
	txtOut<<" "<<fixWidth(ttHdecays[i],widthColumn1,true);
	for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(percentage((float)nbZ1DaughtersFromHttH[i][j]/nbTotalttH[i])+" %",widthOtherColumns,false);
	txtOut<<endl;
      }
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<fixWidth("# of Z2 leptons from H",widthColumn1,true);
      for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(nbZDaughtersFromH[j],widthOtherColumns,false);
      txtOut<<endl;
      for(int i=0; i<nAssocttDecays; i++){
	txtOut<<" "<<fixWidth(ttHdecays[i],widthColumn1,true);
	for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(percentage((float)nbZ2DaughtersFromHttH[i][j]/nbTotalttH[i])+" %",widthOtherColumns,false);
	txtOut<<endl;
      }
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<"the 4 gen-leptons from the H are in the acceptance :"<<endl;
      for(int i=0; i<nAssocttDecays; i++) txtOut<<" "<<fixWidth(ttHdecays[i],widthColumn1,true)<<"   "<<percentage((float)nbHLepsAreInEtaPtAccttH[i]/nbTotalttH[i])+" %"<<endl;
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<"the 4 gen-leptons from the H are reconstructed as good leptons :"<<endl;
      for(int i=0; i<nAssocttDecays; i++) txtOut<<" "<<fixWidth(ttHdecays[i],widthColumn1,true)<<"   "<<percentage((float)nbHLepsAreGoodttH[i]/nbTotalttH[i])+" %"<<endl;
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<"the 4 gen-leptons from the H are the 4 good leptons of the best candidate :"<<endl;
      for(int i=0; i<nAssocttDecays; i++) txtOut<<" "<<fixWidth(ttHdecays[i],widthColumn1,true)<<"   "<<percentage((float)nbAll4LepRightttH[i]/nbTotalttH[i])+" %"<<endl;
    }
    
    txtOut<<endl;

  }

  txtOut.close();


  for(int pr=0; pr<nProcesses; pr++){
    if(!isProcessed[pr]) continue;

    if(doYieldStudy){

      cout<<"For process "<<sProcess[pr]<<": "<<endl;
      cout<<" Nb generated / total yield (L*xsec)   :     "<<NGenEvt[pr]<<"  /  "<<preselExp[pr]<<endl;
      cout<<" Nb / Yield of stored events           :     "<<nbStored[pr][10]<<"  /  "<<yieldStored[pr][10]<<endl;
      cout<<" Nb / Yield of stored events (H->4l)   :     "<<nbStored[pr][7]<<"  /  "<<yieldStored[pr][7]<<endl;
      cout<<"   same, 4 gen l in eta acceptance     :     "<<nbHLepsAreInEtaAcc[pr][7]<<"  /  "<<yieldHLepsAreInEtaAcc[pr][7]<<endl;
      cout<<"   same, 4 gen l in eta+pt acceptance  :     "<<nbHLepsAreInEtaPtAcc[pr][7]<<"  /  "<<yieldHLepsAreInEtaPtAcc[pr][7]<<endl;
      cout<<"   same, >=4 reconstructed leptons     :     "<<nb4RecoLeps[pr][7]<<"  /  "<<yield4RecoLeps[pr][7]<<endl;
      cout<<"   same, 4 gen in eta+pt acc. + 4 reco :     "<<nbHLepsAreInEtaPtAcc4RecoLeps[pr][7]<<"  /  "<<yieldHLepsAreInEtaPtAcc4RecoLeps[pr][7]<<endl;
      //cout<<" Nb / Yield of stored events (H->2L2t) :     "<<nbStored[pr][8]<<"  /  "<<yieldStored[pr][8]<<endl;
      //cout<<" Nb / Yield of stored events (H->else) :     "<<nbStored[pr][9]<<"  /  "<<yieldStored[pr][9]<<endl;
      //cout<<" Nb / Yield of stored events (HX->4lX) :     "<<nb4GenLeps[pr][10]<<"  /  "<<yield4GenLeps[pr][10]<<endl;
      //cout<<" Nb / Yield of events with BC          :     "<<nbWithBC[pr][3]<<"  /  "<<yieldWithBC[pr][3]<<endl;
      cout<<" Nb / Yield of events in SR            :     "<<nbWithBCInSR[pr][3]<<"  /  "<<yieldWithBCInSR[pr][3]<<endl;

      Int_t wi = 6;
      cout<<" Generated final state             ";
      cout<<fixWidth("   4mu",2*wi+1,1)<<"   ";
      cout<<fixWidth("   4e",2*wi+1,1)<<"   ";
      cout<<fixWidth("  2e2mu",2*wi+1,1)<<"   ";
      cout<<fixWidth("   all",2*wi+1,1)<<"   ";
      cout<<endl;
      cout<<" "<<repeat("-",92)<<endl;
      cout<<"                                   ";
      cout<<fixWidth("yield",wi,1)<<"  "<<fixWidth(" eff.",wi,1)<<"   ";
      cout<<fixWidth("yield",wi,1)<<"  "<<fixWidth(" eff.",wi,1)<<"   ";
      cout<<fixWidth("yield",wi,1)<<"  "<<fixWidth(" eff.",wi,1)<<"   ";
      cout<<fixWidth("yield",wi,1)<<"  "<<fixWidth(" eff.",wi,1)<<"   ";
      cout<<endl;
      cout<<" All generated events              ";
      cout<<fixWidth(rounding2(yieldStored[pr][0]),wi,0)<<"  "<<repeat(" ",wi)<<"   ";
      cout<<fixWidth(rounding2(yieldStored[pr][1]),wi,0)<<"  "<<repeat(" ",wi)<<"   ";
      cout<<fixWidth(rounding2(yieldStored[pr][2]),wi,0)<<"  "<<repeat(" ",wi)<<"   ";
      cout<<fixWidth(rounding2(yieldStored[pr][7]),wi,0)<<"  "<<repeat(" ",wi)<<"   ";
      cout<<endl;
      cout<<" 4 gen l in eta acceptance         ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaAcc[pr][0]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaAcc[pr][0]/yieldStored[pr][0])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaAcc[pr][1]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaAcc[pr][1]/yieldStored[pr][1])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaAcc[pr][2]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaAcc[pr][2]/yieldStored[pr][2])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaAcc[pr][7]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaAcc[pr][7]/yieldStored[pr][7])+"%",wi,0)<<"   ";
      cout<<endl;
      cout<<" 4 gen l in eta+pt acceptance      ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc[pr][0]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc[pr][0]/yieldHLepsAreInEtaAcc[pr][0])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc[pr][1]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc[pr][1]/yieldHLepsAreInEtaAcc[pr][1])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc[pr][2]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc[pr][2]/yieldHLepsAreInEtaAcc[pr][2])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc[pr][7]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc[pr][7]/yieldHLepsAreInEtaAcc[pr][7])+"%",wi,0)<<"   ";
      cout<<endl;
      cout<<"  (same but eff wrt. all)          ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc[pr][0]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc[pr][0]/yieldStored[pr][0])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc[pr][1]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc[pr][1]/yieldStored[pr][1])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc[pr][2]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc[pr][2]/yieldStored[pr][2])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc[pr][7]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc[pr][7]/yieldStored[pr][7])+"%",wi,0)<<"   ";
      cout<<endl;
      cout<<" eta+pt acc + >=4 reco leptons     ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc4RecoLeps[pr][0]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc4RecoLeps[pr][0]/yieldHLepsAreInEtaPtAcc[pr][0])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc4RecoLeps[pr][1]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc4RecoLeps[pr][1]/yieldHLepsAreInEtaPtAcc[pr][1])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc4RecoLeps[pr][2]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc4RecoLeps[pr][2]/yieldHLepsAreInEtaPtAcc[pr][2])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc4RecoLeps[pr][7]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc4RecoLeps[pr][7]/yieldHLepsAreInEtaPtAcc[pr][7])+"%",wi,0)<<"   ";
      cout<<endl;
      cout<<" eta+pt acc + in SR                ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccWithBCInSRG[pr][0]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccWithBCInSRG[pr][0]/yieldHLepsAreInEtaPtAcc4RecoLeps[pr][0])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccWithBCInSRG[pr][1]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccWithBCInSRG[pr][1]/yieldHLepsAreInEtaPtAcc4RecoLeps[pr][1])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccWithBCInSRG[pr][2]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccWithBCInSRG[pr][2]/yieldHLepsAreInEtaPtAcc4RecoLeps[pr][2])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccWithBCInSRG[pr][7]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccWithBCInSRG[pr][7]/yieldHLepsAreInEtaPtAcc4RecoLeps[pr][7])+"%",wi,0)<<"   ";
      cout<<endl;
      cout<<" eta+pt acc + in SR + trigger      ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGWithBCInSRG[pr][0]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGWithBCInSRG[pr][0]/yieldHLepsAreInEtaPtAccWithBCInSRG[pr][0])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGWithBCInSRG[pr][1]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGWithBCInSRG[pr][1]/yieldHLepsAreInEtaPtAccWithBCInSRG[pr][1])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGWithBCInSRG[pr][2]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGWithBCInSRG[pr][2]/yieldHLepsAreInEtaPtAccWithBCInSRG[pr][2])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGWithBCInSRG[pr][7]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGWithBCInSRG[pr][7]/yieldHLepsAreInEtaPtAccWithBCInSRG[pr][7])+"%",wi,0)<<"   ";
      cout<<endl;
      cout<<"  (same without single ele path)   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCInSRG[pr][0]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCInSRG[pr][0]/yieldHLepsAreInEtaPtAccWithBCInSRG[pr][0])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCInSRG[pr][1]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCInSRG[pr][1]/yieldHLepsAreInEtaPtAccWithBCInSRG[pr][1])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCInSRG[pr][2]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCInSRG[pr][2]/yieldHLepsAreInEtaPtAccWithBCInSRG[pr][2])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCInSRG[pr][7]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCInSRG[pr][7]/yieldHLepsAreInEtaPtAccWithBCInSRG[pr][7])+"%",wi,0)<<"   ";
      cout<<endl;
      cout<<" eta+pt acc + trigger (% wrt. acc) ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerG[pr][0]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerG[pr][0]/yieldHLepsAreInEtaPtAcc[pr][0])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerG[pr][1]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerG[pr][1]/yieldHLepsAreInEtaPtAcc[pr][1])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerG[pr][2]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerG[pr][2]/yieldHLepsAreInEtaPtAcc[pr][2])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerG[pr][7]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerG[pr][7]/yieldHLepsAreInEtaPtAcc[pr][7])+"%",wi,0)<<"   ";
      cout<<endl;
      cout<<"  (same without single ele path)   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGNo1E[pr][0]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGNo1E[pr][0]/yieldHLepsAreInEtaPtAcc[pr][0])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGNo1E[pr][1]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGNo1E[pr][1]/yieldHLepsAreInEtaPtAcc[pr][1])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGNo1E[pr][2]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGNo1E[pr][2]/yieldHLepsAreInEtaPtAcc[pr][2])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGNo1E[pr][7]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGNo1E[pr][7]/yieldHLepsAreInEtaPtAcc[pr][7])+"%",wi,0)<<"   ";
      cout<<endl;
      cout<<" "<<repeat("-",92)<<endl;

    }

  }
  

  // ------------------------------------------------------------
  // -------------------------- Plots ---------------------------
  // ------------------------------------------------------------

  if(doProdComp){
    TCanvas* cBCInSR[nVariables];
    for(int v=0; v<nVariables; v++){
      if(!plotThisVar[VARLIST][v]) continue;
      cBCInSR[v] = new TCanvas(Form("cBCInSR_%s",varName[v].c_str()),Form("cBCInSR_%s",varName[v].c_str()),500,500);
      DrawByProdmodes(cBCInSR[v],hBCInSR[v],v);
      SaveCanvas(outDir,cBCInSR[v],tagOut);
    }
  }

  if(doProdCompMatch4){
    TCanvas* cBCInSRMatch4[nVariables];
    for(int v=0; v<nVariables; v++){
      if(!plotThisVar[VARLIST][v]) continue;
      cBCInSRMatch4[v] = new TCanvas(Form("cBCInSRMatch4_%s",varName[v].c_str()),Form("cBCInSRMatch4_%s",varName[v].c_str()),500,500);
      DrawByProdmodes(cBCInSRMatch4[v],hBCInSRMatchHLeps[v][0],v);
      SaveCanvas(outDir,cBCInSRMatch4[v],tagOut);
    }
  }

  if(doMatch4OrNot){
    TCanvas* cBCInSRMatch4OrNot[nVariables][nSignalProcesses];
    for(int pr=0; pr<nSignalProcesses; pr++){
      if(!isProcessed[pr]) continue;
      for(int v=0; v<nVariables; v++){
	if(!plotThisVar[VARLIST][v]) continue;
	cBCInSRMatch4OrNot[v][pr] = new TCanvas(Form("cBCInSRMatch4OrNot_%s_%s",varName[v].c_str(),sProcess[pr].c_str()),Form("cBCInSRMatch4OrNot_%s_%s",varName[v].c_str(),sProcess[pr].c_str()),500,500);
	DrawMatch4OrNot(cBCInSRMatch4OrNot[v][pr],hBCInSR[v][pr][FINALSTATE],hBCInSRMatchHLeps[v][0][pr][FINALSTATE],sProcess[pr],allColorsMP[1][pr],allColorsMP[4][pr],v==0);
	SaveCanvas(outDir,cBCInSRMatch4OrNot[v][pr],tagOut);
      }
    }
  }

  if(doMatchHLeps){
    TCanvas* cBCInSRMatchHLeps[nVariables][nSignalProcesses];
    for(int pr=0; pr<nSignalProcesses; pr++){
      if(!isProcessed[pr]) continue;
      for(int v=0; v<nVariables; v++){
	if(!plotThisVar[VARLIST][v]) continue;
	TH1F* h[nMatchHLepsStatuses]; for(int m=0; m<nMatchHLepsStatuses; m++) h[m] = (TH1F*)hBCInSRMatchHLeps[v][m][pr][FINALSTATE];
	cBCInSRMatchHLeps[v][pr] = new TCanvas(Form("cBCInSRMatchHLeps_%s_%s",varName[v].c_str(),sProcess[pr].c_str()),Form("cBCInSRMatchHLeps_%s_%s",varName[v].c_str(),sProcess[pr].c_str()),500,500);
	DrawMatchHLeps(cBCInSRMatchHLeps[v][pr],hBCInSR[v][pr][FINALSTATE],h,sProcess[pr],allColorsPM1[pr],matchHLepsKeys,v==0);
	SaveCanvas(outDir,cBCInSRMatchHLeps[v][pr],tagOut);
      }
    }
  }

  if(doMatchAllLeps){
    TCanvas* cBCInSRMatchAllLeps[nVariables][nSignalProcesses];
    for(int pr=0; pr<nSignalProcesses; pr++){
      if(!isProcessed[pr]) continue;
      if(pr==ggH||pr==qqH) continue;
      for(int v=0; v<nVariables; v++){
	if(!plotThisVar[VARLIST][v]) continue;
	TH1F* h[nMatchAllLepsStatuses]; for(int m=0; m<nMatchAllLepsStatuses; m++) h[m] = (TH1F*)hBCInSRMatchAllLeps[v][m][pr][FINALSTATE];
	cBCInSRMatchAllLeps[v][pr] = new TCanvas(Form("cBCInSRMatchAllLeps_%s_%s",varName[v].c_str(),sProcess[pr].c_str()),Form("cBCInSRMatchAllLeps_%s_%s",varName[v].c_str(),sProcess[pr].c_str()),500,500);
	DrawMatchAllLeps(cBCInSRMatchAllLeps[v][pr],hBCInSR[v][pr][FINALSTATE],h,sProcess[pr],allColorsPM2[pr],matchAllLepsKeys,v==0);
	SaveCanvas(outDir,cBCInSRMatchAllLeps[v][pr],tagOut);
      }
    }
  }
  
  if(doMatchWHZHttH){

    if(isProcessed[WH]){
      TCanvas* cBCInSRMatchWH[nVariables];
      for(int v=0; v<nVariables; v++){
	if(!plotThisVar[VARLIST][v]) continue;
	TH1F* hMatchWH[nMatchWHStatuses]; for(int m=0; m<nMatchWHStatuses; m++) hMatchWH[m] = (TH1F*)hBCInSRMatchWH[v][m][FINALSTATE];
	cBCInSRMatchWH[v] = new TCanvas(Form("cBCInSRMatchWH_%s",varName[v].c_str()),Form("cBCInSRMatchWH_%s",varName[v].c_str()),500,500);
	if(v==1 || v==3)
	  DrawMatchCustom(cBCInSRMatchWH[v],hBCInSR[v][2][FINALSTATE],hMatchWH,sProcess[WH],colorsMatchWH,matchWH,nMatchWHStatuses,0.11,0.6,v==0);
	else
	  DrawMatchCustom(cBCInSRMatchWH[v],hBCInSR[v][2][FINALSTATE],hMatchWH,sProcess[WH],colorsMatchWH,matchWH,nMatchWHStatuses,0.34,0.6,v==0);
	SaveCanvas(outDir,cBCInSRMatchWH[v],tagOut);
      }
    }

    if(isProcessed[ZH]){
      TCanvas* cBCInSRMatchZH[nVariables];
      for(int v=0; v<nVariables; v++){
	if(!plotThisVar[VARLIST][v]) continue;
	TH1F* hMatchZH[nMatchZHStatuses]; for(int m=0; m<nMatchZHStatuses; m++) hMatchZH[m] = (TH1F*)hBCInSRMatchZH[v][m][FINALSTATE];
	cBCInSRMatchZH[v] = new TCanvas(Form("cBCInSRMatchZH_%s",varName[v].c_str()),Form("cBCInSRMatchZH_%s",varName[v].c_str()),500,500);
	if(v==1 || v==3)
	  DrawMatchCustom(cBCInSRMatchZH[v],hBCInSR[v][3][FINALSTATE],hMatchZH,sProcess[ZH],colorsMatchZH,matchZH,nMatchZHStatuses,0.11,0.5,v==0);
	else
	  DrawMatchCustom(cBCInSRMatchZH[v],hBCInSR[v][3][FINALSTATE],hMatchZH,sProcess[ZH],colorsMatchZH,matchZH,nMatchZHStatuses,0.34,0.5,v==0);
	SaveCanvas(outDir,cBCInSRMatchZH[v],tagOut);
      }
    }

    if(isProcessed[ttH]){
      TCanvas* cBCInSRMatchttH[nVariables];
      for(int v=0; v<nVariables; v++){
	if(!plotThisVar[VARLIST][v]) continue;
	TH1F* hMatchttH[nMatchttHStatuses]; for(int m=0; m<nMatchttHStatuses; m++) hMatchttH[m] = (TH1F*)hBCInSRMatchttH[v][m][FINALSTATE];
	cBCInSRMatchttH[v] = new TCanvas(Form("cBCInSRMatchttH_%s",varName[v].c_str()),Form("cBCInSRMatchttH_%s",varName[v].c_str()),500,500);
	if(v==1 || v==3)
	  DrawMatchCustom(cBCInSRMatchttH[v],hBCInSR[v][4][FINALSTATE],hMatchttH,sProcess[ttH],colorsMatchttH,matchttH,nMatchttHStatuses,0.11,0.4,v==0);
	else
	  DrawMatchCustom(cBCInSRMatchttH[v],hBCInSR[v][4][FINALSTATE],hMatchttH,sProcess[ttH],colorsMatchttH,matchttH,nMatchttHStatuses,0.34,0.4,v==0);
	SaveCanvas(outDir,cBCInSRMatchttH[v],tagOut);
      }
    }

  }
  
  if(do2DPlots){

    TCanvas* c2DBCInSR[nVarPairs][nProcesses];
    for(int v2=0; v2<nVarPairs; v2++){
      if(!plotThisVarPair[VARPAIRLIST][v2]) continue;
      for(int pr=0; pr<nProcesses; pr++){
	if(!isProcessed[pr]) continue;
	string name = "c2DBCInSR_"+varName[varPairRef[v2][0]]+"_"+varName[varPairRef[v2][1]]+"_"+sProcess[pr];
	c2DBCInSR[v2][pr] = new TCanvas(name.c_str(),name.c_str(),500,500);
	Draw2D(c2DBCInSR[v2][pr],h2DBCInSR[v2][pr][FINALSTATE],v2,sProcess[pr]);
	SaveCanvas(outDir,c2DBCInSR[v2][pr],tagOut);
      }
    }

    TCanvas* c2DBCInSRDecays[nVarPairs][nProcesses][nDecays];
    for(int v2=0; v2<nVarPairs; v2++){
      if(!plotThisVarPair[VARPAIRLIST][v2]) continue;
      if(!varPairRef[v2][2]) continue;
      for(int pr=0; pr<nProcesses; pr++){
	if(!isProcessed[pr]) continue;
	if(pr!=WH&&pr!=ZH&&pr!=ttH) continue;
	for(int d=0; d<nDecays; d++){
	  if(pr==WH&&(d==0||d==2)) continue;
	  string name = "c2DBCInSR_"+decayInfix[d]+"_"+varName[varPairRef[v2][0]]+"_"+varName[varPairRef[v2][1]]+"_"+sProcess[pr];
	  c2DBCInSRDecays[v2][pr][d] = new TCanvas(name.c_str(),name.c_str(),500,500);
	  Draw2D(c2DBCInSRDecays[v2][pr][d],h2DBCInSRDecays[v2][pr][d][FINALSTATE],v2,sProcess[pr]+decayLabel[d]);
	  SaveCanvas(outDir,c2DBCInSRDecays[v2][pr][d],tagOut);
	}
      }
    }

  }

  if(doROCs){
    for(int ro=0; ro<nROCs; ro++){
      gRoc[ro] = doROC(hRocSgnl[ro],hRocBkgd[ro],(bool)rocRef[ro][3],1,1,"",Form("%s efficiency",sProcess[rocRef[ro][1]].c_str()),Form("%s efficiency",sProcess[rocRef[ro][2]].c_str()));
      gWp[ro] = (!rocWP[ro]) ? 0 : doWP(hRocSgnl[ro],hRocBkgd[ro],rocWP[ro],(bool)rocRef[ro][3],1,1);
    }

    //Int_t rocColor[9] = {kRed+1,kAzure,kViolet,kGreen+2,kYellow+1,kOrange-3,kOrange+6,kOrange-8,kGray};
    //Int_t rocColor[9] = {kRed+1,kRed-7,kRed,kBlue+1,kBlue-9,kBlue-4,kGreen+1,kYellow+1,kGray+1};
    Int_t rocColor1[8] = {kRed+1,kRed-7,kRed,kBlue+1,kBlue-9,kBlue-4,kGreen+1,kGray+1};
    Int_t rocWidth1[9] = {1,1,2,1,1,2,2,2,2};
    Int_t rocColor2[8] = {kRed,kBlue-4,kGreen+1,kYellow+2,kOrange+1,kMagenta-9,kMagenta,kMagenta+3};//,kMagenta-9,kMagenta-7,kMagenta+1,kMagenta+3};
    //Int_t rocWidth2[9] = {2,2,2,2,2,2,2,2,2};
    Int_t rocWidth2[9] = {1,1,1,1,1,1,1,1,1};
    Int_t rocColor3[6] = {kRed,kBlue-4,kGreen+1,kMagenta-9,kMagenta,kMagenta+3};//,kMagenta-9,kMagenta-7,kMagenta+1,kMagenta+3};

    /*
    //DrawRoc(gRoc,gWp,"cCompareRocs_VBF2j",{2,4,12,17,18,19,20,21,0},rocColor1,rocWidth1,outDir,tagOut);
    DrawRoc(gRoc,gWp,"cCompareRocs_VBF2j",{2,4,12,17,18,19,41,0},rocColor1,rocWidth1,outDir,tagOut);
    //DrawRoc(gRoc,gWp,"cCompareRocs_VBF2j",{2,4,51,17,18,52,53,0},rocColor1,rocWidth1,outDir,tagOut); // This gives the same ROCs, as expected
    DrawRoc(gRoc,gWp,"cCompareRocs_VBF2j_qqZZ",{3,5,13,22,23,24,43,1},rocColor1,rocWidth1,outDir,tagOut);
    //DrawRoc(gRoc,gWp,"cCompareRocs_VBF1j",{8,9,14,27,28,29,30},rocColor1,rocWidth1,outDir,tagOut);
    DrawRoc(gRoc,gWp,"cCompareRocs_VBF1j",{8,9,14,27,28,29,45},rocColor1,rocWidth1,outDir,tagOut);
    DrawRoc(gRoc,gWp,"cCompareRocs_WHhadr",{10,6,15,31,32,33,47},rocColor1,rocWidth1,outDir,tagOut);
    DrawRoc(gRoc,gWp,"cCompareRocs_ZHhadr",{11,7,16,36,37,38,49},rocColor1,rocWidth1,outDir,tagOut);
    //*/

    // DrawRoc(gRoc,gWp,"cCompareRocs_VBF2j_transform",{12,19,41,54,55,56,57,58/*,59*/},rocColor2,rocWidth2,outDir,tagOut);
    // DrawRoc(gRoc,gWp,"cCompareRocs_VBF2j_qqZZ_transform",{13,24,43,60,61,62,63,64/*,65*/},rocColor2,rocWidth2,outDir,tagOut);
    // DrawRoc(gRoc,gWp,"cCompareRocs_VBF1j_transform",{14,29,45,66,67,68/*,69*/},rocColor3,rocWidth2,outDir,tagOut);
    // DrawRoc(gRoc,gWp,"cCompareRocs_WHhadr_transform",{15,33,47,70,71,72/*,73*/},rocColor3,rocWidth2,outDir,tagOut);
    // DrawRoc(gRoc,gWp,"cCompareRocs_ZHhadr_transform",{16,38,49,74,75,76/*,77*/},rocColor3,rocWidth2,outDir,tagOut);


    Int_t rocWidthOnes[9] = {1,1,1,1,1,1,1,1,1};
    Int_t rocColS[6] = {kGreen+1,kBlue-4,kRed,kMagenta-9,kMagenta,kMagenta+3};
    Int_t rocColS2[6] = {kGreen+1,kBlue-4,kRed,kGray+1};
    //* //for slides
    //DrawRoc(gRoc,gWp,"cCompareRocs_VBF2j_regular",{12,19,41,0},rocColS2,rocWidthOnes,outDir,tagOut);
    DrawRoc(gRoc,gWp,"cCompareRocs_VBF2j_transform",{12,19,41,56,57,58},rocColS,rocWidthOnes,outDir,tagOut);
    //DrawRoc(gRoc,gWp,"cCompareRocs_VBF2j_qqZZ_transform",{13,24,43,62,63,64},rocColS,rocWidthOnes,outDir,tagOut);
    DrawRoc(gRoc,gWp,"cCompareRocs_VBF1j_transform",{14,29,45,66,67,68/*,69*/},rocColS,rocWidthOnes,outDir,tagOut);
    DrawRoc(gRoc,gWp,"cCompareRocs_WHhadr_transform",{15,33,47,70,71,72/*,73*/},rocColS,rocWidthOnes,outDir,tagOut);
    DrawRoc(gRoc,gWp,"cCompareRocs_ZHhadr_transform",{16,38,49,74,75,76/*,77*/},rocColS,rocWidthOnes,outDir,tagOut);
    //*/

    if(doROCzooms){
      DrawRoc(gRoc,gWp,"cCompareRocs_VBF2j_transform_zoom",{12,19,41,54,55,56,57,58/*,59*/},rocColor2,rocWidth2,outDir,tagOut,true);
      //DrawRoc(gRoc,gWp,"cCompareRocs_VBF2j_qqZZ_transform_zoom",{13,24,43,60,61,62,63,64/*,65*/},rocColor2,rocWidth2,outDir,tagOut,true);
      DrawRoc(gRoc,gWp,"cCompareRocs_VBF1j_transform_zoom",{14,29,45,66,67,68/*,69*/},rocColor3,rocWidth2,outDir,tagOut,true);
      DrawRoc(gRoc,gWp,"cCompareRocs_WHhadr_transform_zoom",{15,33,47,70,71,72/*,73*/},rocColor3,rocWidth2,outDir,tagOut,true);
      DrawRoc(gRoc,gWp,"cCompareRocs_ZHhadr_transform_zoom",{16,38,49,74,75,76/*,77*/},rocColor3,rocWidth2,outDir,tagOut,true);
    }
  }

  if(doEVWs){
    TCanvas* cEvw[nEVWs];
    for(int e=0; e<8; e++){//nEVWs; e++){
      for(int pr=0; pr<nProcesses; pr++){
	if(!isProcessed[pr]) continue;
	gEvw[e][pr] = doEFF(hEvw[e][pr],(bool)evwRef[e][1],"process efficiency",1,1);
      }
      cEvw[e] = new TCanvas(Form("cEvw_%s",varName[evwRef[e][0]].c_str()),Form("cEvw_%s",varName[evwRef[e][0]].c_str()),500,500);
      DrawEvw(cEvw[e],gEvw[e],e);
      SaveCanvas(outDir,cEvw[e],tagOut);
    }
  }
  
  if(doBaskets){

    TCanvas* cBCInSRBasketsAll = new TCanvas("cBCInSR_baskets","cBCInSR_baskets",longDim[BASKETLIST],500);
    DrawByProdmodes(cBCInSRBasketsAll,hBCInSRBaskets,-1);
    SaveCanvas(outDir,cBCInSRBasketsAll,tagOut);
    
    TCanvas* cBCInSRBasketEfficiency[nProcesses];
    TH1F* hProdModes[nProcesses]; for(int pr=0; pr<nProcesses; pr++) hProdModes[pr] = (TH1F*)hBCInSRBaskets[pr][FINALSTATE];
    TH1F* hAssocDecays[nAssocDecays]; for(int a=0; a<nAssocDecays; a++) hAssocDecays[a] = (TH1F*)hBCInSRBasketsAssocDecays[a][FINALSTATE];
    for(int pr=0; pr<nProcesses; pr++){
      if(!isProcessed[pr]) continue;
      cBCInSRBasketEfficiency[pr] = new TCanvas(Form("cBCInSR_basketEfficiency_%s",sProcess[pr].c_str()),Form("cBCInSR_basketEfficiency_%s",sProcess[pr].c_str()),500,longDim[BASKETLIST]);
      DrawBasketEfficiencies(cBCInSRBasketEfficiency[pr],hProdModes[pr],pr,hAssocDecays,treatH2l2XAsBkgd?assocDecayName2:assocDecayName1,treatH2l2XAsBkgd,basketLabel[BASKETLIST],lumiText,false);
      SaveCanvas(outDir,cBCInSRBasketEfficiency[pr],tagOut);
    }

    TCanvas* cBCInSRBasketPurity = new TCanvas("cBCInSR_basketPurity","cBCInSR_basketPurity",800,longDim[BASKETLIST]);
    DrawBasketPurities(cBCInSRBasketPurity,hProdModes,hAssocDecays,treatH2l2XAsBkgd?assocDecayName4:assocDecayName3,false,treatH2l2XAsBkgd,basketLabel[BASKETLIST],labelMerge,lumiText);
    SaveCanvas(outDir,cBCInSRBasketPurity,tagOut);
    TCanvas* cBCInSRBasketPuritySplit = new TCanvas("cBCInSR_basketPuritySplit","cBCInSR_basketPuritySplit",800,longDim[BASKETLIST]);
    DrawBasketPurities(cBCInSRBasketPuritySplit,hProdModes,hAssocDecays,treatH2l2XAsBkgd?assocDecayName4:assocDecayName3,true,treatH2l2XAsBkgd,basketLabel[BASKETLIST],labelMerge,lumiText);
    SaveCanvas(outDir,cBCInSRBasketPuritySplit,tagOut);

    if(isProcessed[qqZZ]){

      TH1F* hSumSgnl = (TH1F*)hBCInSRBaskets[ggH][FINALSTATE]->Clone();
      for(int pr=ggH; pr<nSignalProcesses; pr++){
	if(!isProcessed[pr]) continue;
	hSumSgnl->Add(hBCInSRBaskets[pr][FINALSTATE]);
      }
      TH1F* hSumBkgd = (TH1F*)hBCInSRBaskets[qqZZ][FINALSTATE]->Clone();
      for(int pr=qqZZ; pr<nProcesses; pr++){
	if(!isProcessed[pr]) continue;
	hSumBkgd->Add(hBCInSRBaskets[pr][FINALSTATE]);
      }
      if(treatH2l2XAsBkgd){
	hSumSgnl->Add(hBCInSRBasketsAssocDecays[nAssocWDecays+nAssocZDecays-1][FINALSTATE],-1);
	hSumSgnl->Add(hBCInSRBasketsAssocDecays[nAssocWDecays+nAssocZDecays+nAssocttDecays-1][FINALSTATE],-1);
	hSumBkgd->Add(hBCInSRBasketsAssocDecays[nAssocWDecays+nAssocZDecays-1][FINALSTATE]);
	hSumBkgd->Add(hBCInSRBasketsAssocDecays[nAssocWDecays+nAssocZDecays+nAssocttDecays-1][FINALSTATE]);
      }

      TH1F* hDenom = (TH1F*)hSumSgnl->Clone();
      hDenom->Add(hSumBkgd);
      TH1F* hSOSPB = (TH1F*)hSumSgnl->Clone();
      hSOSPB->Divide(hDenom);
      TCanvas* cBCInSRBasketSOSPB = new TCanvas("cBCInSR_basketSOSPB","cBCInSR_basketSOSPB",500,longDim[BASKETLIST]);
      DrawBasketSOSPB(cBCInSRBasketSOSPB,hSOSPB,basketLabel[BASKETLIST]);
      SaveCanvas(outDir,cBCInSRBasketSOSPB,tagOut);

    }

    for(int pr=0; pr<nProcesses; pr++) 
      if(isProcessed[pr]) cout<<"process "<<sProcess[pr]<<": "<<0.5*hBCInSRBaskets[pr][FINALSTATE]->Integral()<<endl;

    //*
    for(int pr=0; pr<nProcesses; pr++) 
      for(int ba=0; ba<nBaskets; ba++) 
	if(isProcessed[pr]) cout<<" basket #"<<ba<<", process "<<sProcess[pr]<<": "<<hBCInSRBaskets[pr][FINALSTATE]->GetBinContent(ba+1)<<endl;
    //*/

  }

  if(doPlotsByGenChan){

    TCanvas* cSortedPt [nSignalProcesses][nGenChannels];
    TCanvas* cSortedEta[nSignalProcesses][nGenChannels];
    TCanvas* cNRecoLep[nSignalProcesses][nGenChannels];
    TCanvas* cNRecoMuo[nSignalProcesses][nGenChannels];
    TCanvas* cNRecoEle[nSignalProcesses][nGenChannels];
    TCanvas* cGenHPt[nSignalProcesses];
    TCanvas* cGenHRapidity[nSignalProcesses];
    TCanvas* cGenHRapidityVsPt[nSignalProcesses];

    for(int pr=0; pr<nSignalProcesses; pr++){
      if(!isProcessed[pr]) continue;
      if(!(pr==ggH)) continue;

      for(int gc=0; gc<3; gc++){
	string suffix = sProcess[pr]+"_"+genChannels[gc];

	cSortedPt[pr][gc] = new TCanvas(("cSortedPt_"+suffix).c_str(),("cSortedPt_"+suffix).c_str(),500,500);
	DrawSorted(cSortedPt[pr][gc],h1GenptHLepsAreInEtaAcc[pr][gc],h1RecoptWithBCInSRGPassTriggerG[pr][gc],sProcess[pr],genChannels[gc],0.1);
	SaveCanvas(outDir,cSortedPt[pr][gc],tagOut);
	cSortedEta[pr][gc] = new TCanvas(("cSortedEta_"+suffix).c_str(),("cSortedEta_"+suffix).c_str(),500,500);
	DrawSorted(cSortedEta[pr][gc],h1GenetaHLepsAreInPtAcc[pr][gc],h1RecoetaWithBCInSRGPassTriggerG[pr][gc],sProcess[pr],genChannels[gc],0.1);
	SaveCanvas(outDir,cSortedEta[pr][gc],tagOut);

	cNRecoLep[pr][gc] = new TCanvas(("cNRecoLep_"+suffix).c_str(),("cNRecoLep_"+suffix).c_str(),500,500);
	Draw1D(cNRecoLep[pr][gc],h1NRecoLepHLepsAreInEtaPtAcc[pr][gc]);
	SaveCanvas(outDir,cNRecoLep[pr][gc],tagOut);
	cNRecoMuo[pr][gc] = new TCanvas(("cNRecoMuo_"+suffix).c_str(),("cNRecoMuo_"+suffix).c_str(),500,500);
	Draw1D(cNRecoMuo[pr][gc],h1NRecoMuoHLepsAreInEtaPtAcc[pr][gc]);
	SaveCanvas(outDir,cNRecoMuo[pr][gc],tagOut);
	cNRecoEle[pr][gc] = new TCanvas(("cNRecoEle_"+suffix).c_str(),("cNRecoEle_"+suffix).c_str(),500,500);
	Draw1D(cNRecoEle[pr][gc],h1NRecoEleHLepsAreInEtaPtAcc[pr][gc]);
	SaveCanvas(outDir,cNRecoEle[pr][gc],tagOut);

      }

      cGenHPt[pr] = new TCanvas(("cGenHPt_"+sProcess[pr]).c_str(),("cGenHPt_"+sProcess[pr]).c_str(),500,500);
      Draw1D(cGenHPt[pr],h1GenHPt[pr]);
      SaveCanvas(outDir,cGenHPt[pr],tagOut);
      cGenHRapidity[pr] = new TCanvas(("cGenHRapidity_"+sProcess[pr]).c_str(),("cGenHRapidity_"+sProcess[pr]).c_str(),500,500);
      Draw1D(cGenHRapidity[pr],h1GenHRapidity[pr]);
      SaveCanvas(outDir,cGenHRapidity[pr],tagOut);
      cGenHRapidityVsPt[pr] = new TCanvas(("cGenHRapidityVsPt_"+sProcess[pr]).c_str(),("cGenHRapidityVsPt_"+sProcess[pr]).c_str(),500,500);
      Draw2D(cGenHRapidityVsPt[pr],h2GenHRapidityVsPt[pr],-1,"");
      SaveCanvas(outDir,cGenHRapidityVsPt[pr],tagOut);

    }

  }

}





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// MAIN MACRO //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void plotProcesses() {


  // --------------- definitions ---------------

  string inputPath = "";
  string treeTag = "";

  string outputPath = ""; 

  // Define the luminosity
  float lumi = 35.86706; 
  string lumiText = "35.9 fb^{-1}";

  // Choose an m4l window
  //* // window used in plots of AN/PAS
  float m4l_min = 118.;
  float m4l_max = 130.;
  //*/
  /* // full signal region 
  float m4l_min = 70.;
  float m4l_max = 3000.;
  //*/

  string mWindowTag = string(Form("m4l%ito%i",(int)m4l_min,(int)m4l_max));
  string basketListTag = string(Form("bList%i",BASKETLIST));

  // --------------- processing ---------------

  setTDRStyle();

  gSystem->Exec(("mkdir -p "+outputPath).c_str());	
	  
  run( 
      inputPath,
      outputPath,
      lumi,
      lumiText,
      m4l_min,
      m4l_max,
      mWindowTag+"_"+basketListTag+"_"+treeTag
       );
  
}
