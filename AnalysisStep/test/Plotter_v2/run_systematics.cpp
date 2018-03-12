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

// My own files
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/include/Systematics.h>
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/include/Variables.h>
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/src/setTDRStyle.cpp>

using namespace std;

int main( int argc, char *argv[] )
{
   setTDRStyle();
	
	TString path = "";
   TString file_name = "/ZZ4lAnalysis.root";
	
   // Signal
   TString ggH125         = path + "ggH125ext"       + file_name;
   TString ggH125_TU      = path + "ggH125_tuneup"   + file_name;
   TString ggH125_TD      = path + "ggH125_tunedown" + file_name;
	
   TString VBFH125        = path + "VBFH125ext"       + file_name;
   TString VBFH125_TU     = path + "VBFH125_tuneup"   + file_name;
   TString VBFH125_TD     = path + "VBFH125_tunedown" + file_name;
	
   TString WpH125         = path + "WplusH125ext" + file_name;
   TString WpH125_TU      = path + "WplusH125_tuneup" + file_name;
   TString WpH125_TD      = path + "WplusH125_tunedown" + file_name;
	
   TString WmH125         = path + "WminusH125ext" + file_name;
   TString WmH125_TU      = path + "WminusH125_tuneup" + file_name;
   TString WmH125_TD      = path + "WminusH125_tunedown" + file_name;
	
   TString ZH125          = path + "ZH125ext" + file_name;
   TString ZH125_TU       = path + "ZH125_tuneup" + file_name;
   TString ZH125_TD       = path + "ZH125_tunedown" + file_name;
	
   TString ttH125         = path + "ttH125ext" + file_name;
   TString ttH125_TU      = path + "ttH125_tuneup" + file_name;
   TString ttH125_TD      = path + "ttH125_tunedown" + file_name;
	
   TString bbH125         = path + "bbH125ext" + file_name;
   TString bbH125_TU      = path + "bbH125_tuneup" + file_name;
   TString bbH125_TD      = path + "bbH125_tunedown" + file_name;
	
   TString tqH125         = path + "tqH125ext" + file_name;
   TString tqH125_TU      = path + "tqH125_tuneup" + file_name;
   TString tqH125_TD      = path + "tqH125_tunedown" + file_name;
	
   Systematics *systematics = new Systematics();
	
	
//====================
// Print Systematics
//====================

	systematics->PrintSystematics_JEC(ggH125);
	
   delete systematics;
}
