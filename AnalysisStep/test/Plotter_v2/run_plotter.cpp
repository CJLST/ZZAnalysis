// C++
#include <iostream>
#include <fstream>
#include <string>

// ROOT
#include "TApplication.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"

// My own files
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/include/Plotter.h>
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/include/Settings.h>
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/src/setTDRStyle.cpp>

using namespace std;

int main( int argc, char *argv[] )
{
   setTDRStyle();

   TString eos_path  = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/";

   TString path_2018 = eos_path + "MC_2018/";
   TString path_2017 = eos_path + "MC_2017/";
   TString path_2016 = eos_path + "MC_2016_CorrectBTag/";
   TString file_name = "/ZZ4lAnalysis.root";

   TString eos_data_path  = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200430_LegacyRun2/";
   TString Data_2018    = eos_data_path + "Data_2018/AllData"        + file_name;
   TString Data_2017    = eos_data_path + "Data_2017/AllData"        + file_name;
   TString Data_2016    = eos_data_path + "Data_2016/AllData"        + file_name;

   TString FR_2018      = "../../data/FakeRates/newData_FakeRates_SS_2018.root";
   TString FR_2017      = "../../data/FakeRates/newData_FakeRates_SS_2017.root";
   TString FR_2016      = "../../data/FakeRates/newData_FakeRates_SS_2016.root";

   TString ggH125_2018      = path_2018 + "ggH125"     + file_name;
   TString VBFH125_2018     = path_2018 + "VBFH125"    + file_name;
   TString WpH125_2018      = path_2018 + "WplusH125"  + file_name;
   TString WmH125_2018      = path_2018 + "WminusH125" + file_name;
   TString ZH125_2018       = path_2018 + "ZH125"      + file_name;
   TString ggZH125_2018     = path_2018 + "ggZH125"    + file_name;
   TString ttH125_2018      = path_2018 + "ttH125"     + file_name;
   TString bbH125_2018      = path_2018 + "bbH125"     + file_name;
   TString tqH125_2018      = path_2018 + "tqH125"     + file_name;

   TString ZZTo4l_2018      = path_2018 + "ZZTo4lext1"                 + file_name;
   TString ggZZ4e_2018      = path_2018 + "ggTo4e_Contin_MCFM701"      + file_name;
   TString ggZZ4mu_2018     = path_2018 + "ggTo4mu_Contin_MCFM701"     + file_name;
   TString ggZZ4tau_2018    = path_2018 + "ggTo4tau_Contin_MCFM701"    + file_name;
   TString ggZZ2e2mu_2018   = path_2018 + "ggTo2e2mu_Contin_MCFM701"   + file_name;
   TString ggZZ2e2tau_2018  = path_2018 + "ggTo2e2tau_Contin_MCFM701"  + file_name;
   TString ggZZ2mu2tau_2018 = path_2018 + "ggTo2mu2tau_Contin_MCFM701" + file_name;

   // Triboson and TT-triboson like samples
   TString TTZZ_2018        = path_2018 + "TTZZ"                       + file_name;
   TString TTWW_2018        = path_2018 + "TTWW"                       + file_name;
   TString WWZ_2018         = path_2018 + "WWZ"                        + file_name;
   TString WZZ_2018         = path_2018 + "WZZ"                        + file_name;
   TString ZZZ_2018         = path_2018 + "ZZZ"                        + file_name;
   TString TTZJets_2018     = path_2018 + "TTZJets_M10_MLMext1"            + file_name;
   TString TTLLNuNu_2018    = path_2018 + "TTZToLLNuNu_M10ext1"            + file_name;
   TString TTLL_2018        = path_2018 + "TTZToLL_M1to10_MLM"         + file_name;
   // VBF off-shell sample
   TString VBFoff_2018      = path_2018 + "VBFToContinToZZ4l"          + file_name;

   TString ggH125_2017      = path_2017 + "ggH125"     + file_name;
   TString VBFH125_2017     = path_2017 + "VBFH125"    + file_name;
   TString WpH125_2017      = path_2017 + "WplusH125"  + file_name;
   TString WmH125_2017      = path_2017 + "WminusH125" + file_name;
   TString ZH125_2017       = path_2017 + "ZH125"      + file_name;
   TString ggZH125_2017     = path_2017 + "ggZH125"    + file_name;
   TString ttH125_2017      = path_2017 + "ttH125"     + file_name;
   TString bbH125_2017      = path_2017 + "bbH125"     + file_name;
   TString tqH125_2017      = path_2017 + "tqH125"     + file_name;

   TString ZZTo4l_2017      = path_2017 + "ZZTo4l"                     + file_name;
   TString ggZZ4e_2017      = path_2017 + "ggTo4e_Contin_MCFM701"      + file_name;
   TString ggZZ4mu_2017     = path_2017 + "ggTo4mu_Contin_MCFM701"     + file_name;
   TString ggZZ4tau_2017    = path_2017 + "ggTo4tau_Contin_MCFM701"    + file_name;
   TString ggZZ2e2mu_2017   = path_2017 + "ggTo2e2mu_Contin_MCFM701"   + file_name;
   TString ggZZ2e2tau_2017  = path_2017 + "ggTo2e2tau_Contin_MCFM701"  + file_name;
   TString ggZZ2mu2tau_2017 = path_2017 + "ggTo2mu2tau_Contin_MCFM701" + file_name;

   // Triboson and TT-triboson like samples
   TString TTZZ_2017        = path_2017 + "TTZZ"                       + file_name;
   TString TTWW_2017        = path_2017 + "TTWW"                       + file_name;
   TString WWZ_2017         = path_2017 + "WWZ"                        + file_name;
   TString WZZ_2017         = path_2017 + "WZZ"                        + file_name;
   TString ZZZ_2017         = path_2017 + "ZZZ"                        + file_name;
   TString TTZJets_2017     = path_2017 + "TTZJets_M10_MLM"            + file_name;
   TString TTLLNuNu_2017    = path_2017 + "TTZToLLNuNu_M10"            + file_name;
   TString TTLL_2017        = path_2017 + "TTZToLL_M1to10_MLM"         + file_name;
   // VBF off-shell sample
   TString VBFoff_2017      = path_2017 + "VBFToContinToZZ4l"          + file_name;

   TString ggH125_2016      = path_2016 + "ggH125"     + file_name;
   TString VBFH125_2016     = path_2016 + "VBFH125"    + file_name;
   TString WpH125_2016      = path_2016 + "WplusH125"  + file_name;
   TString WmH125_2016      = path_2016 + "WminusH125" + file_name;
   TString ZH125_2016       = path_2016 + "ZH125"      + file_name;
   TString ggZH125_2016     = path_2016 + "ggZH125"    + file_name;
   TString ttH125_2016      = path_2016 + "ttH125"     + file_name;
   TString bbH125_2016      = path_2016 + "bbH125"     + file_name;
   TString tqH125_2016      = path_2016 + "tqH125"     + file_name;

   TString ZZTo4l_2016      = path_2016 + "ZZTo4l"                     + file_name;
   TString ggZZ4e_2016      = path_2016 + "ggTo4e_Contin_MCFM701"      + file_name;
   TString ggZZ4mu_2016     = path_2016 + "ggTo4mu_Contin_MCFM701"     + file_name;
   TString ggZZ4tau_2016    = path_2016 + "ggTo4tau_Contin_MCFM701"    + file_name;
   TString ggZZ2e2mu_2016   = path_2016 + "ggTo2e2mu_Contin_MCFM701"   + file_name;
   TString ggZZ2e2tau_2016  = path_2016 + "ggTo2e2tau_Contin_MCFM701"  + file_name;
   TString ggZZ2mu2tau_2016 = path_2016 + "ggTo2mu2tau_Contin_MCFM701" + file_name;

   // Triboson and TT-triboson like samples
   TString TTZZ_2016        = path_2016 + "TTZZ"                       + file_name;
   TString TTWW_2016        = path_2016 + "TTWW"                       + file_name;
   TString WWZ_2016         = path_2016 + "WWZ"                        + file_name;
   TString WZZ_2016         = path_2016 + "WZZ"                        + file_name;
   TString ZZZ_2016         = path_2016 + "ZZZ"                        + file_name;
   TString TTZJets_2016     = path_2016 + "TTZJets_M10_MLM"            + file_name;
   TString TTLLNuNu_2016    = path_2016 + "TTZToLLNuNu_M10"            + file_name;
   TString TTLL_2016        = path_2016 + "TTZToLL_M1to10_MLM"         + file_name;
   // VBF off-shell sample
   TString VBFoff_2016      = path_2016 + "VBFToContinToZZ4l"          + file_name;
 

   Plotter *plotter = new Plotter( );

   plotter->SetBlinding(110, 138, 300, 1200);


////// 2016
// // Data
   plotter->MakeHistograms(Data_2016,2016);

// // Signal(s)
   plotter->MakeHistograms(ggH125_2016,2016);
   plotter->MakeHistograms(VBFH125_2016,2016);
   plotter->MakeHistograms(ZH125_2016,2016);
   plotter->MakeHistograms(ggZH125_2016, 2016);
   plotter->MakeHistograms(ttH125_2016,2016);
   plotter->MakeHistograms(bbH125_2016,2016);
   plotter->MakeHistograms(tqH125_2016,2016);
   plotter->MakeHistograms(WpH125_2016,2016);
   plotter->MakeHistograms(WmH125_2016,2016);

// // MC Backgrounds
   plotter->MakeHistograms(ZZTo4l_2016,2016);
   plotter->MakeHistograms(ggZZ4e_2016,2016);
   plotter->MakeHistograms(ggZZ4mu_2016,2016);
   plotter->MakeHistograms(ggZZ4tau_2016,2016);
   plotter->MakeHistograms(ggZZ2e2mu_2016,2016);
   plotter->MakeHistograms(ggZZ2e2tau_2016,2016);
   plotter->MakeHistograms(ggZZ2mu2tau_2016,2016);
   plotter->MakeHistograms(TTZZ_2016,2016);
   plotter->MakeHistograms(TTWW_2016,2016);
   plotter->MakeHistograms(WWZ_2016,2016);
   plotter->MakeHistograms(WZZ_2016,2016);
   plotter->MakeHistograms(ZZZ_2016,2016);
   plotter->MakeHistograms(TTZJets_2016, 2016);
   plotter->MakeHistograms(TTLLNuNu_2016, 2016);
   plotter->MakeHistograms(TTLL_2016, 2016);
   plotter->MakeHistograms(VBFoff_2016, 2016);

// // ZX Background
   plotter->MakeHistogramsZX(Data_2016, FR_2016, 2016);


////// 2017
// // Data
   plotter->MakeHistograms(Data_2017,2017);

// // Signal(s)
   plotter->MakeHistograms(ggH125_2017,2017);
   plotter->MakeHistograms(VBFH125_2017,2017);
   plotter->MakeHistograms(ZH125_2017,2017);
   plotter->MakeHistograms(ggZH125_2017, 2017);
   plotter->MakeHistograms(ttH125_2017,2017);
   plotter->MakeHistograms(bbH125_2017,2017);
   plotter->MakeHistograms(tqH125_2017,2017);
   plotter->MakeHistograms(WpH125_2017,2017);
   plotter->MakeHistograms(WmH125_2017,2017);

// // MC Backgrounds
   plotter->MakeHistograms(ZZTo4l_2017,2017);
   plotter->MakeHistograms(ggZZ4e_2017,2017);
   plotter->MakeHistograms(ggZZ4mu_2017,2017);
   plotter->MakeHistograms(ggZZ4tau_2017,2017);
   plotter->MakeHistograms(ggZZ2e2mu_2017,2017);
   plotter->MakeHistograms(ggZZ2e2tau_2017,2017);
   plotter->MakeHistograms(ggZZ2mu2tau_2017,2017);
   plotter->MakeHistograms(TTZZ_2017,2017);
   plotter->MakeHistograms(TTWW_2017,2017);
   plotter->MakeHistograms(WWZ_2017,2017);
   plotter->MakeHistograms(WZZ_2017,2017);
   plotter->MakeHistograms(ZZZ_2017,2017);
   plotter->MakeHistograms(TTZJets_2017, 2017);
   plotter->MakeHistograms(TTLLNuNu_2017, 2017);
   plotter->MakeHistograms(TTLL_2017, 2017);
   plotter->MakeHistograms(VBFoff_2017, 2017);

// // ZX Background
   plotter->MakeHistogramsZX(Data_2017, FR_2017, 2017);



////// 2018
// // Data
   plotter->MakeHistograms(Data_2018,2018);

// // Signal(s)
   plotter->MakeHistograms(ggH125_2018,2018);
   plotter->MakeHistograms(VBFH125_2018,2018);
   plotter->MakeHistograms(ZH125_2018,2018);
   plotter->MakeHistograms(ggZH125_2018, 2018);
   plotter->MakeHistograms(ttH125_2018,2018);
   plotter->MakeHistograms(bbH125_2018,2018);
   plotter->MakeHistograms(tqH125_2018,2018);
   plotter->MakeHistograms(WpH125_2018,2018);
   plotter->MakeHistograms(WmH125_2018,2018);

// // MC Backgrounds
   plotter->MakeHistograms(ZZTo4l_2018,2018);
   plotter->MakeHistograms(ggZZ4e_2018,2018);
   plotter->MakeHistograms(ggZZ4mu_2018,2018);
   plotter->MakeHistograms(ggZZ4tau_2018,2018);
   plotter->MakeHistograms(ggZZ2e2mu_2018,2018);
   plotter->MakeHistograms(ggZZ2e2tau_2018,2018);
   plotter->MakeHistograms(ggZZ2mu2tau_2018,2018);
   plotter->MakeHistograms(TTZZ_2018,2018);
   plotter->MakeHistograms(TTWW_2018,2018);
   plotter->MakeHistograms(WWZ_2018,2018);
   plotter->MakeHistograms(WZZ_2018,2018);
   plotter->MakeHistograms(ZZZ_2018,2018);
   plotter->MakeHistograms(TTZJets_2018, 2018);
   plotter->MakeHistograms(TTLLNuNu_2018, 2018);
   plotter->MakeHistograms(TTLL_2018, 2018);
   plotter->MakeHistograms(VBFoff_2018, 2018);

// // ZX Background
   plotter->MakeHistogramsZX(Data_2018, FR_2018, 2018);



   plotter->MakeM4lZX();

   plotter->FillInclusive();

   plotter->Save();


//===========================
// Plotting of blinded plots
//===========================
   // setTDRStyle();
   // plotter->GetHistos( "Blinded" );

   // plotter->plot_1D_all_cat("Blinded", "M4lMain",       "Plots/Blinded");
   // plotter->plot_1D_all_cat("Blinded", "M4lMainZoomed", "Plots/Blinded");

   // plotter->plot_1D_all_fs("Blinded", "M4lMain",       "Plots/Blinded");
   // plotter->plot_1D_all_fs("Blinded", "M4lMainZoomed", "Plots/Blinded");

   // plotter->plot_1D_single("Blinded", "M4lMainHighMass", "Plots/Blinded", Settings::fs4l, Settings::inclusive);
   // plotter->plot_1D_single("Blinded", "MZ1",             "Plots/Blinded", Settings::fs4l, Settings::inclusive);
   // plotter->plot_1D_single("Blinded", "MZ2",             "Plots/Blinded", Settings::fs4l, Settings::inclusive);
   // plotter->plot_1D_single("Blinded", "MZ1",             "Plots/Blinded", Settings::fs4e, Settings::inclusive);
   // plotter->plot_1D_single("Blinded", "MZ2",             "Plots/Blinded", Settings::fs4e, Settings::inclusive);
   // plotter->plot_1D_single("Blinded", "MZ1",             "Plots/Blinded", Settings::fs4mu, Settings::inclusive);
   // plotter->plot_1D_single("Blinded", "MZ2",             "Plots/Blinded", Settings::fs4mu, Settings::inclusive);
   // plotter->plot_1D_single("Blinded", "MZ1",             "Plots/Blinded", Settings::fs2e2mu, Settings::inclusive);
   // plotter->plot_1D_single("Blinded", "MZ2",             "Plots/Blinded", Settings::fs2e2mu, Settings::inclusive);

   // plotter->plot_1D_single("Blinded", "KD", "Plots/Blinded", Settings::fs4l, Settings::inclusive);
   // plotter->plot_1D_single("Blinded", "KD", "Plots/Blinded", Settings::fs4e, Settings::inclusive);
   // plotter->plot_1D_single("Blinded", "KD", "Plots/Blinded", Settings::fs4mu, Settings::inclusive);
   // plotter->plot_1D_single("Blinded", "KD", "Plots/Blinded", Settings::fs2e2mu, Settings::inclusive);

   // plotter->plot_1D_single("Blinded", "DVBFDEC", "Plots/Blinded", Settings::fs4l, Settings::VBF_2j_tagged);
   // plotter->plot_1D_single("Blinded", "DVHDEC",  "Plots/Blinded", Settings::fs4l, Settings::VH_hadron_tagged);

   // plotter->plot_1D_single("Blinded", "D1jet", "Plots/Blinded", Settings::fs4l, Settings::inclusive);
   // plotter->plot_1D_single("Blinded", "D2jet", "Plots/Blinded", Settings::fs4l, Settings::inclusive);
   // plotter->plot_1D_single("Blinded", "DVH",   "Plots/Blinded", Settings::fs4l, Settings::inclusive);

   // plotter->plot_2D_single("Blinded", "MZ1vsMZ2", "Plots/Blinded", Settings::inclusive);

   // plotter->plot_2D_error_single("Blinded", "KDvsM4l",          "Plots/Blinded", Settings::inclusive);
   // plotter->plot_2D_error_single("Blinded", "KDvsM4lZoomed",    "Plots/Blinded", Settings::inclusive);
   // plotter->plot_2D_error_single("Blinded", "KDvsM4lHighMass",  "Plots/Blinded", Settings::inclusive);
   // plotter->plot_2D_error_single("Blinded", "D1jetvsM4lZoomed", "Plots/Blinded", Settings::inclusive);
   // plotter->plot_2D_error_single("Blinded", "D2jetvsM4lZoomed", "Plots/Blinded", Settings::inclusive);
   // plotter->plot_2D_error_single("Blinded", "DWHvsM4lZoomed",   "Plots/Blinded", Settings::inclusive);
   // plotter->plot_2D_error_single("Blinded", "DZHvsM4lZoomed",   "Plots/Blinded", Settings::inclusive);
   // plotter->plot_2D_error_single("Blinded", "DVHvsM4lZoomed",   "Plots/Blinded", Settings::inclusive);

   // plotter->plot_2D_error_all_cat("Blinded", "KDvsM4lZoomed",    "Plots/Blinded");
   // plotter->plot_2D_error_all_cat("Blinded", "DVBFDECvsM4lZoomed",    "Plots/Blinded");
   // plotter->plot_2D_error_all_cat("Blinded", "DVHDECvsM4lZoomed",    "Plots/Blinded");


//=============================
// Plotting of unblinded plots
//=============================
/*
   setTDRStyle(); // Needed to reset margins set by 2D histograms

   plotter->GetHistos("Unblinded");

   plotter->plot_Purity("Unblinded", "Plots/Unblinded");setTDRStyle();
   //
   plotter->plot_STXS("Unblinded", "Plots/Unblinded");setTDRStyle();
   //
   plotter->plot_1D_all_cat("Unblinded", "M4lMain",       "Plots/Unblinded");
   plotter->plot_1D_all_cat("Unblinded", "M4lMainZoomed", "Plots/Unblinded");

   plotter->plot_1D_all_fs("Unblinded", "M4lMain",       "Plots/Unblinded");
   plotter->plot_1D_all_fs("Unblinded", "M4lMainZoomed", "Plots/Unblinded");
   
   plotter->plot_1D_single("Unblinded", "M4lMainHighMass", "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "MZ1",             "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "MZ2",             "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "MZ1_M4L118130",   "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "MZ2_M4L118130",   "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_2D_single("Unblinded", "MZ1vsMZ2",           "Plots/Unblinded", Settings::inclusive);
   plotter->plot_2D_single("Unblinded", "MZ1vsMZ2_M4L118130", "Plots/Unblinded", Settings::inclusive);setTDRStyle();
   
   plotter->plot_1D_single("Unblinded", "KD",               "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "DVBFDEC",          "Plots/Unblinded", Settings::fs4l, Settings::VBF_2j_tagged);
   plotter->plot_1D_single("Unblinded", "DVHDEC",           "Plots/Unblinded", Settings::fs4l, Settings::VH_hadron_tagged);
   plotter->plot_1D_single("Unblinded", "KD_M4L118130",     "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "DVBFDEC_M4L118130","Plots/Unblinded", Settings::fs4l, Settings::VBF_2j_tagged);
   plotter->plot_1D_single("Unblinded", "DVHDEC_M4L118130", "Plots/Unblinded", Settings::fs4l, Settings::VH_hadron_tagged);
   
   plotter->plot_1D_single("Unblinded", "PFMET",              "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "Pt4l",               "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "Eta4l",              "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "Pt_leading",         "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "Pt_trailing",        "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "Eta_leading",        "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "Eta_trailing",       "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "SIP_leading",        "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "SIP_trailing",       "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "ISO_leading",        "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "ISO_trailing",       "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "NExtraLep",          "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "NJets",              "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "NJetsBTagged",       "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "M4l_110150_HighKD",  "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   
   plotter->plot_1D_single("Unblinded", "D1jet",  "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "D2jet",  "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "DWH",    "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "DZH",    "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "DVH",    "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "D1jet_M4L118130",  "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "D2jet_M4L118130",  "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "DWH_M4L118130",    "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "DZH_M4L118130",    "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "DVH_M4L118130",    "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   
   plotter->plot_2D_error_single("Unblinded", "KDvsM4l",          "Plots/Unblinded", Settings::inclusive);
   plotter->plot_2D_error_single("Unblinded", "KDvsM4lZoomed",    "Plots/Unblinded", Settings::inclusive);
   plotter->plot_2D_error_single("Unblinded", "KDvsM4lHighMass",  "Plots/Unblinded", Settings::inclusive);
   plotter->plot_2D_error_single("Unblinded", "D1jetvsM4lZoomed", "Plots/Unblinded", Settings::inclusive);
   plotter->plot_2D_error_single("Unblinded", "D2jetvsM4lZoomed", "Plots/Unblinded", Settings::inclusive);
   plotter->plot_2D_error_single("Unblinded", "DWHvsM4lZoomed",   "Plots/Unblinded", Settings::inclusive);
   plotter->plot_2D_error_single("Unblinded", "DZHvsM4lZoomed",   "Plots/Unblinded", Settings::inclusive);
   plotter->plot_2D_error_single("Unblinded", "DVHvsM4lZoomed",   "Plots/Unblinded", Settings::inclusive);
   
   plotter->plot_2D_error_all_cat("Unblinded", "KDvsM4lZoomed",    "Plots/Unblinded");
   plotter->plot_2D_error_all_cat("Unblinded", "DVBFDECvsM4lZoomed",    "Plots/Unblinded");
   plotter->plot_2D_error_all_cat("Unblinded", "DVHDECvsM4lZoomed",    "Plots/Unblinded");
*/

//=============================
// PAS HIG-19-001 figures
//=============================
//
   setTDRStyle(); // Needed to reset margins set by 2D histograms

   plotter->GetHistos("Unblinded");

   plotter->plot_Purity("Unblinded", "Plots/Unblinded");setTDRStyle();
   //
   plotter->plot_STXS("Unblinded", "Plots/Unblinded");setTDRStyle();
   //
   plotter->plot_1D_all_fs("Unblinded", "M4lMain",       "Plots/Unblinded");
   plotter->plot_1D_all_fs("Unblinded", "M4lMainZoomed", "Plots/Unblinded");
   plotter->plot_1D_single("Unblinded", "M4lMain", "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "M4lMainZoomed", "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_all_cat("Unblinded", "M4lMainZoomed", "Plots/Unblinded");
   plotter->plot_1D_single("Unblinded", "MZ1_M4L118130",   "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "MZ2_M4L118130",   "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "KD_M4L118130",     "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "DVBFDEC_M4L118130","Plots/Unblinded", Settings::fs4l, Settings::VBF_2j_tagged);
   plotter->plot_1D_single("Unblinded", "DVHDEC_M4L118130", "Plots/Unblinded", Settings::fs4l, Settings::VH_hadron_tagged);
   plotter->plot_1D_single("Unblinded", "D1jet_M4L118130",  "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "D2jet_M4L118130",  "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "DVH_M4L118130",    "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_2D_single("Unblinded", "MZ1vsMZ2_M4L118130", "Plots/Unblinded", Settings::inclusive);
   plotter->plot_2D_error_all_cat("Unblinded", "KDvsM4lZoomed",    "Plots/Unblinded");
   plotter->plot_2D_error_all_cat("Unblinded", "DVBFDECvsM4lZoomed",    "Plots/Unblinded");
   plotter->plot_2D_error_all_cat("Unblinded", "DVHDECvsM4lZoomed",    "Plots/Unblinded");

   delete plotter;
}
