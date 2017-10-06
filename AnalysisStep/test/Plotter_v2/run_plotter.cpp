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
   
   TString path = "";
   TString file_name = "/ZZ4lAnalysis.root";
   TString file_name_FR = "/FakeRate_SS_Moriond368.root";
   
   TString Data        = path + "Data" + file_name;
   TString ggH125      = path + "ggH125" + file_name;
   TString VBFH125     = path + "VBFH125" + file_name;
   TString WpH125      = path + "WplusH125" + file_name;
   TString WmH125      = path + "WminusH125" + file_name;
   TString ZH125       = path + "ZH125" + file_name;
   TString ttH125      = path + "ttH125" + file_name;
   TString ZZTo4l      = path + "ZZTo4l" + file_name;
   TString ggZZ4e      = path + "ggTo4e" + file_name;
   TString ggZZ4mu     = path + "ggTo4mu" + file_name;
   TString ggZZ4tau    = path + "ggTo4tau" + file_name;
   TString ggZZ2e2mu   = path + "ggTo2e2mu" + file_name;
   TString ggZZ2e2tau  = path + "ggTo2e2tau" + file_name;
   TString ggZZ2mu2tau = path + "ggTo2mu2tau" + file_name;
   TString FakeRates   = "../../data/FakeRates" + file_name_FR;

   Plotter *plotter = new Plotter( 35.86706 ); 

   plotter->SetBlinding(110, 150, 300, 1000);

   plotter->MakeHistograms(Data);
   plotter->MakeHistograms(ggH125);
   plotter->MakeHistograms(VBFH125);
   plotter->MakeHistograms(ZH125);
   plotter->MakeHistograms(ttH125);
   plotter->MakeHistograms(WpH125);
   plotter->MakeHistograms(WmH125);
   plotter->MakeHistograms(ZZTo4l);
   plotter->MakeHistograms(ggZZ4e);
   plotter->MakeHistograms(ggZZ4mu);
   plotter->MakeHistograms(ggZZ4tau);
   plotter->MakeHistograms(ggZZ2e2mu);
   plotter->MakeHistograms(ggZZ2e2tau);
   plotter->MakeHistograms(ggZZ2mu2tau);


   plotter->MakeHistogramsZX(Data, FakeRates);
   plotter->MakeM4lZX();

   plotter->FillInclusive();

   plotter->Save();


//===========================
// Plotting of blinded plots
//===========================
   plotter->GetHistos( "Blinded" );
   
   plotter->plot_1D_all_cat("Blinded", "M4lMain",       "Plots/Blinded");
   plotter->plot_1D_all_cat("Blinded", "M4lMainZoomed", "Plots/Blinded");
   
   plotter->plot_1D_all_fs("Blinded", "M4lMain",       "Plots/Blinded");
   plotter->plot_1D_all_fs("Blinded", "M4lMainZoomed", "Plots/Blinded");
   
   plotter->plot_1D_single("Blinded", "M4lMainHighMass", "Plots/Blinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Blinded", "MZ1",             "Plots/Blinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Blinded", "MZ2",             "Plots/Blinded", Settings::fs4l, Settings::inclusive);
   
   plotter->plot_1D_single("Blinded", "KD", "Plots/Blinded", Settings::fs4l, Settings::inclusive);
   
   
   plotter->plot_2D_single("Blinded", "MZ1vsMZ2", "Plots/Blinded", Settings::inclusive);
   
   plotter->plot_2D_error_single("Blinded", "KDvsM4l",          "Plots/Blinded", Settings::inclusive);
   plotter->plot_2D_error_single("Blinded", "KDvsM4lZoomed",    "Plots/Blinded", Settings::inclusive);
   plotter->plot_2D_error_single("Blinded", "KDvsM4lHighMass",  "Plots/Blinded", Settings::inclusive);
   plotter->plot_2D_error_single("Blinded", "D1jetvsM4lZoomed", "Plots/Blinded", Settings::inclusive);
   plotter->plot_2D_error_single("Blinded", "D2jetvsM4lZoomed", "Plots/Blinded", Settings::inclusive);
   plotter->plot_2D_error_single("Blinded", "DWHvsM4lZoomed",   "Plots/Blinded", Settings::inclusive);
   plotter->plot_2D_error_single("Blinded", "DZHvsM4lZoomed",   "Plots/Blinded", Settings::inclusive);
   plotter->plot_2D_error_single("Blinded", "DVHvsM4lZoomed",   "Plots/Blinded", Settings::inclusive);

   plotter->plot_2D_error_all_cat("Blinded", "KDvsM4lZoomed",    "Plots/Blinded");
   plotter->plot_2D_error_all_cat("Blinded", "D1jetvsM4lZoomed", "Plots/Blinded");
   plotter->plot_2D_error_all_cat("Blinded", "D2jetvsM4lZoomed", "Plots/Blinded");
   plotter->plot_2D_error_all_cat("Blinded", "DWHvsM4lZoomed",   "Plots/Blinded");
   plotter->plot_2D_error_all_cat("Blinded", "DZHvsM4lZoomed",   "Plots/Blinded");
   plotter->plot_2D_error_all_cat("Blinded", "DVHvsM4lZoomed",   "Plots/Blinded");
   
   
   
//=============================
// Plotting of unblinded plots
//=============================

   setTDRStyle(); // Needed to reset margins set by 2D histograms

   plotter->GetHistos("Unblinded");
   
   plotter->plot_1D_all_cat("Unblinded", "M4lMain",       "Plots/Unblinded");
   plotter->plot_1D_all_cat("Unblinded", "M4lMainZoomed", "Plots/Unblinded");
      
   plotter->plot_1D_all_fs("Unblinded", "M4lMain",       "Plots/Unblinded");
   plotter->plot_1D_all_fs("Unblinded", "M4lMainZoomed", "Plots/Unblinded");

   plotter->plot_1D_single("Unblinded", "M4lMainHighMass", "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "MZ1",             "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "MZ2",             "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "MZ1_M4L118130",   "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "MZ2_M4L118130",   "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
  
   plotter->plot_1D_single("Unblinded", "KD",              "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "KD_M4L118130",    "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "D1jet_M4L118130", "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "D2jet_M4L118130", "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "DWH_M4L118130",   "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "DZH_M4L118130",   "Plots/Unblinded", Settings::fs4l, Settings::inclusive);
   plotter->plot_1D_single("Unblinded", "DVH_M4L118130",   "Plots/Unblinded", Settings::fs4l, Settings::inclusive);

   
   plotter->plot_2D_single("Unblinded", "MZ1vsMZ2",           "Plots/Unblinded", Settings::inclusive);
   plotter->plot_2D_single("Unblinded", "MZ1vsMZ2_M4L118130", "Plots/Unblinded", Settings::inclusive);
   
   plotter->plot_2D_error_single("Unblinded", "KDvsM4l",          "Plots/Unblinded", Settings::inclusive);
   plotter->plot_2D_error_single("Unblinded", "KDvsM4lZoomed",    "Plots/Unblinded", Settings::inclusive);
   plotter->plot_2D_error_single("Unblinded", "KDvsM4lHighMass",  "Plots/Unblinded", Settings::inclusive);
   plotter->plot_2D_error_single("Unblinded", "D1jetvsM4lZoomed", "Plots/Unblinded", Settings::inclusive);
   plotter->plot_2D_error_single("Unblinded", "D2jetvsM4lZoomed", "Plots/Unblinded", Settings::inclusive);
   plotter->plot_2D_error_single("Unblinded", "DWHvsM4lZoomed",   "Plots/Unblinded", Settings::inclusive);
   plotter->plot_2D_error_single("Unblinded", "DZHvsM4lZoomed",   "Plots/Unblinded", Settings::inclusive);
   plotter->plot_2D_error_single("Unblinded", "DVHvsM4lZoomed",   "Plots/Unblinded", Settings::inclusive);
   
   plotter->plot_2D_error_all_cat("Unblinded", "KDvsM4lZoomed",    "Plots/Unblinded");
   plotter->plot_2D_error_all_cat("Unblinded", "D1jetvsM4lZoomed", "Plots/Unblinded");
   plotter->plot_2D_error_all_cat("Unblinded", "D2jetvsM4lZoomed", "Plots/Unblinded");
   plotter->plot_2D_error_all_cat("Unblinded", "DWHvsM4lZoomed",   "Plots/Unblinded");
   plotter->plot_2D_error_all_cat("Unblinded", "DZHvsM4lZoomed",   "Plots/Unblinded");
   plotter->plot_2D_error_all_cat("Unblinded", "DVHvsM4lZoomed",   "Plots/Unblinded");
   
   delete plotter;
}
