// C++
#include <iostream>
#include <fstream>
#include <string>

// ROOT
#include "TApplication.h"
#include <TROOT.h>
#include "TStyle.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"

// My own files
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/include/Yields.h>
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/include/Variables.h>
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/src/setTDRStyle.cpp>

using namespace std;

int main( int argc, char *argv[] )
{
   setTDRStyle();
   
   TString path = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/";
   TString file_name = "/ZZ4lAnalysis.root";
   TString file_name_FR = "/FakeRates_SS_2016.root";
   
   //Data
   TString Data    = path + "Data_2016/AllData"        + file_name;
   TString FakeRates   = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/FRfiles" + file_name_FR;
   
   // Signal
   // TString ggH120      = path + "ggH120" + file_name;
   // TString ggH124      = path + "ggH124" + file_name;
   TString ggH125      = path + "MC_2016/"+ "ggH125" + file_name;
   // TString ggH126      = path + "MC_2016/"+ "ggH126" + file_name;
   // TString ggH130      = path + "MC_2016/"+ "ggH130" + file_name;
               
   // TString VBFH120     = path + "MC_2016/"+ "VBFH120" + file_name;
   // TString VBFH124     = path + "MC_2016/"+ "VBFH124" + file_name;
   TString VBFH125     = path + "MC_2016/"+ "VBFH125" + file_name;
   // TString VBFH126     = path + "MC_2016/"+ "VBFH126" + file_name;
   // TString VBFH130     = path + "MC_2016/"+ "VBFH130" + file_name;
   
   // TString WpH120      = path + "MC_2016/"+ "WplusH120" + file_name;
   // TString WpH124      = path + "MC_2016/"+ "WplusH124" + file_name;
   TString WpH125      = path + "MC_2016/"+ "WplusH125" + file_name;
   // TString WpH126      = path + "MC_2016/"+ "WplusH126" + file_name;
   // TString WpH130      = path + "MC_2016/"+ "WplusH130" + file_name;
    
   // TString WmH120      = path + "MC_2016/"+ "WminusH120" + file_name;
   // TString WmH124      = path + "MC_2016/"+ "WminusH124" + file_name;
   TString WmH125      = path + "MC_2016/"+ "WminusH125" + file_name;
   // TString WmH126      = path + "MC_2016/"+ "WminusH126" + file_name;
   // TString WmH130      = path + "MC_2016/"+ "WminusH130" + file_name;
      
   // TString ZH120       = path + "MC_2016/"+ "ZH120" + file_name;
   // TString ZH124       = path + "MC_2016/"+ "ZH124" + file_name;
   TString ZH125       = path + "MC_2016/"+ "ZH125" + file_name;
   // TString ZH126       = path + "MC_2016/"+ "ZH126" + file_name;
   // TString ZH130       = path + "MC_2016/"+ "ZH130" + file_name;
   
   // TString ttH120      = path + "MC_2016/"+ "ttH120" + file_name;
   // TString ttH124      = path + "MC_2016/"+ "ttH124" + file_name;
   TString ttH125      = path + "MC_2016/"+ "ttH125" + file_name;
   // TString ttH126      = path + "MC_2016/"+ "ttH126" + file_name;
   // TString ttH130      = path + "MC_2016/"+ "ttH130" + file_name;
   
   // TString bbH120      = path + "MC_2016/"+ "bbH120" + file_name;
 //   TString bbH124      = path + "MC_2016/"+ "bbH124" + file_name;
   TString bbH125      = path + "MC_2016/"+ "bbH125" + file_name;
   // TString bbH126      = path + "MC_2016/"+ "bbH126" + file_name;
   // TString bbH130      = path + "MC_2016/"+ "bbH130" + file_name;
   
   TString tqH125      = path + "MC_2016/"+ "tqH125" + file_name;

   // Backgrounds
   TString ZZTo4l      = path + "MC_2016/"+ "ZZTo4l" + file_name;
   TString ggZZ4e      = path + "MC_2016/"+ "ggTo4e_Contin_MCFM701" + file_name;
   TString ggZZ4mu     = path + "MC_2016/"+ "ggTo4mu_Contin_MCFM701" + file_name;
   TString ggZZ4tau    = path + "MC_2016/"+ "ggTo4tau_Contin_MCFM701" + file_name;
   TString ggZZ2e2mu   = path + "MC_2016/"+ "ggTo2e2mu_Contin_MCFM701" + file_name;
   TString ggZZ2e2tau  = path + "MC_2016/"+ "ggTo2e2tau_Contin_MCFM701" + file_name;
   TString ggZZ2mu2tau = path + "MC_2016/"+ "ggTo2mu2tau_Contin_MCFM701" + file_name;
   
   double lumi = 35.9;
   int year = 2016;
   Yields *yields = new Yields(lumi);
   
//===============
// Produce plots 
//===============
   
   yields->MakeHistograms(Data, year);

   // yields->MakeHistograms(ggH120, year);
   // yields->MakeHistograms(ggH124, year);
   yields->MakeHistograms(ggH125, year);
   // yields->MakeHistograms(ggH126, year);
   // yields->MakeHistograms(ggH130, year);

   // yields->MakeHistograms(VBFH120, year);
   // yields->MakeHistograms(VBFH124, year);
   yields->MakeHistograms(VBFH125, year);
   // yields->MakeHistograms(VBFH126, year);
   // yields->MakeHistograms(VBFH130, year);

   // yields->MakeHistograms(ZH120, year);
   // yields->MakeHistograms(ZH124, year);
   yields->MakeHistograms(ZH125, year);
   // yields->MakeHistograms(ZH126, year);
   // yields->MakeHistograms(ZH130, year);

   // yields->MakeHistograms(WpH120, year);
   // yields->MakeHistograms(WpH124, year);
   yields->MakeHistograms(WpH125, year);
   // yields->MakeHistograms(WpH126, year);
   // yields->MakeHistograms(WpH130, year);

   // yields->MakeHistograms(WmH120, year);
   // yields->MakeHistograms(WmH124, year);
   yields->MakeHistograms(WmH125, year);
   // yields->MakeHistograms(WmH126, year);
   // yields->MakeHistograms(WmH130, year);

   // yields->MakeHistograms(ttH120, year);
   // yields->MakeHistograms(ttH124, year);
   yields->MakeHistograms(ttH125, year);
   // yields->MakeHistograms(ttH126, year);
   // yields->MakeHistograms(ttH130, year);

// yields->MakeHistograms(bbH120, year);
   // yields->MakeHistograms(bbH124, year);
   yields->MakeHistograms(bbH125, year);
   // yields->MakeHistograms(bbH126, year);
   // yields->MakeHistograms(bbH130, year);

   yields->MakeHistograms(tqH125, year);

   yields->MakeHistograms(ZZTo4l, year);
   yields->MakeHistograms(ggZZ4e, year);
   yields->MakeHistograms(ggZZ4mu, year);
   yields->MakeHistograms(ggZZ4tau, year);
   yields->MakeHistograms(ggZZ2e2mu, year);
   yields->MakeHistograms(ggZZ2e2tau, year);
   yields->MakeHistograms(ggZZ2mu2tau, year);

   yields->FillInclusive();

   yields->Save();
   
//==============
// Print Yields
//==============
   
   yields->GetHistos("Yields");
   yields->Calculate_SS_ZX_Yields( Data, FakeRates);
   yields->Print("Yields");

   yields->Print("Yields", 118., 130.);
   yields->Print("Yields", 105., 140.);

   yields->PrintLatexTables("Yields", 118., 130.);
   yields->FillGraphs("Yields", 105., 140., "Q");
   yields->PrepareYamlFiles("Yields", "13", 105., 140.);
   
//==========================================
// Produce data ROOT files for datacard maker
//==========================================
   yields->ProduceDataROOTFiles( Data, "DataROOTFiles" );
   
   
   delete yields;
}
