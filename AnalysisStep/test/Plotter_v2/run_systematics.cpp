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
   TString ggH125         = path + "ggH125"          + file_name;
   TString ggH125ext      = path + "ggH125ext"       + file_name;
   TString ggH125_TU      = path + "ggH125_tuneup"   + file_name;
   TString ggH125_TD      = path + "ggH125_tunedown" + file_name;
	
	TString VBFH125        = path + "VBFH125"          + file_name;
   TString VBFH125ext     = path + "VBFH125ext"       + file_name;
   TString VBFH125_TU     = path + "VBFH125_tuneup"   + file_name;
   TString VBFH125_TD     = path + "VBFH125_tunedown" + file_name;
	
	TString WpH125         = path + "WplusH125"          + file_name;
   TString WpH125ext      = path + "WplusH125ext"       + file_name;
   TString WpH125_TU      = path + "WplusH125_tuneup"   + file_name;
   TString WpH125_TD      = path + "WplusH125_tunedown" + file_name;
	
	TString WmH125         = path + "WminusH125ext"       + file_name;
   TString WmH125ext      = path + "WminusH125ext"       + file_name;
   TString WmH125_TU      = path + "WminusH125_tuneup"   + file_name;
   TString WmH125_TD      = path + "WminusH125_tunedown" + file_name;
	
	TString ZH125          = path + "ZH125"          + file_name;
   TString ZH125ext       = path + "ZH125ext"       + file_name;
   TString ZH125_TU       = path + "ZH125_tuneup"   + file_name;
   TString ZH125_TD       = path + "ZH125_tunedown" + file_name;
	
	TString ttH125         = path + "ttH125"          + file_name;
   TString ttH125ext      = path + "ttH125ext"       + file_name;
   TString ttH125_TU      = path + "ttH125_tuneup"   + file_name;
   TString ttH125_TD      = path + "ttH125_tunedown" + file_name;
	
   TString bbH125         = path + "bbH125" + file_name;

   TString tqH125         = path + "tqH125" + file_name;

	
   // Backgrounds
   TString ZZTo4l      = path + "ZZTo4l"                     + file_name;
   TString ZZTo4lext   = path + "ZZTo4lext"                  + file_name;
   TString ggZZ4e      = path + "ggTo4e_Contin_MCFM701"      + file_name;
   TString ggZZ4mu     = path + "ggTo4mu_Contin_MCFM701"     + file_name;
   TString ggZZ4tau    = path + "ggTo4tau_Contin_MCFM701"    + file_name;
   TString ggZZ2e2mu   = path + "ggTo2e2mu_Contin_MCFM701"   + file_name;
   TString ggZZ2e2tau  = path + "ggTo2e2tau_Contin_MCFM701"  + file_name;
   TString ggZZ2mu2tau = path + "ggTo2mu2tau_Contin_MCFM701" + file_name;
	
   Systematics *systematics = new Systematics();
	
	
//====================
// Print Systematics
//====================

	systematics->FillSystematics(ggH125ext);
	systematics->FillSystematics(VBFH125ext);
	systematics->FillSystematics(WpH125ext);
	systematics->FillSystematics(WmH125ext);
	systematics->FillSystematics(ZH125ext);
	systematics->FillSystematics(ttH125ext);
	systematics->FillSystematics(bbH125);
	systematics->FillSystematics(tqH125);
	systematics->FillSystematics(ZZTo4l);
	systematics->FillSystematics(ggZZ4e);
	systematics->FillSystematics(ggZZ4mu);
	systematics->FillSystematics(ggZZ4tau);
	systematics->FillSystematics(ggZZ2e2mu);
	systematics->FillSystematics(ggZZ2e2tau);
	systematics->FillSystematics(ggZZ2mu2tau);

	systematics->FillSystematics_tuneUpDn(ggH125_TU);
	systematics->FillSystematics_tuneUpDn(ggH125_TD);
	systematics->FillSystematics_tuneUpDn(VBFH125_TU);
	systematics->FillSystematics_tuneUpDn(VBFH125_TD);
	systematics->FillSystematics_tuneUpDn(WpH125_TU);
	systematics->FillSystematics_tuneUpDn(WpH125_TD);
	systematics->FillSystematics_tuneUpDn(WmH125_TU);
	systematics->FillSystematics_tuneUpDn(WmH125_TD);
//	systematics->FillSystematics_tuneUpDn(ZH125_TU);
//	systematics->FillSystematics_tuneUpDn(ZH125_TD);
	systematics->FillSystematics_tuneUpDn(ttH125_TU);
	systematics->FillSystematics_tuneUpDn(ttH125_TD);

	systematics->PrintSystematics_PU();
	systematics->PrintSystematics_JEC();
	systematics->PrintSystematics_BTag();
	systematics->PrintSystematics_muRmuFScale();
	systematics->PrintSystematics_PythiaScale();
	systematics->PrintSystematics_PythiaTune();
	systematics->PrintSystematics_QCDScale();
	systematics->PrintSystematics_PDFScale();
	systematics->PrintSystematics_EWCorr();
	systematics->PrintSystematics_THU_ggH();
	
	
   delete systematics;
}
