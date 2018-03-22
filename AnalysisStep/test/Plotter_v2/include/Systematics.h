#ifndef Systematics_h
#define Systematics_h

// C++
#include <iostream>
#include <fstream>
#include <iomanip> // For setprecision
#include <vector>
#include <map>
#include <cstdlib>

// ROOT
#include "TApplication.h"
#include "TStyle.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"

// Include classes
#include "Tree.h"
#include "Histograms.h"
#include <ZZAnalysis/AnalysisStep/interface/Category.h>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>
#include <ZZAnalysis/AnalysisStep/interface/bitops.h>
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/include/FakeRates.h>
#include <ZZAnalysis/AnalysisStep/interface/cConstants.h>
#include <ZZAnalysis/AnalysisStep/interface/Discriminants.h>

using namespace std;

class Systematics: public Tree
{

public:
	
	Systematics();
	~Systematics();
	
	void FillSystematics( TString );
	void FillSystematics_tuneUpDn( TString );
	void PrintSystematics_PU( );
	void PrintSystematics_JEC( );
	void PrintSystematics_BTag( );
	void PrintSystematics_muRmuFScale( );
	void PrintSystematics_PythiaScale( );
	void PrintSystematics_PythiaTune( );
	void PrintSystematics_QCDScale( );
	void PrintSystematics_PDFScale( );
	void PrintSystematics_EWCorr( );
	void PrintSystematics_THU_ggH( );

	
private:

	float calculate_K_factor( TString );
   int FindFinalState();
   int find_current_production_mode( TString , int, int);
   int CountAssociatedLeptons();
	
   TFile *input_file, *input_file_data;
   TTree *input_tree, *input_tree_data;
	
   TH1F* hCounters;
	
   Long64_t n_gen_events;
	
   float jetPt[99];
   float jetEta[99];
   float jetPhi[99];
   float jetMass[99];
   float jetQGL[99];
   float jetPgOverPq[99];
	
   int _current_process, _current_production_mode , _current_final_state;
   int _current_category;
   int _n_gen_assoc_lep;
   float _k_factor, _yield_SR;
   double gen_sum_weights, _event_weight, _event_weight_UP, _event_weight_DN;
	
	vector<TString> _s_category, _s_production_mode;
	vector<int>   gen_assoc_lep_id_;
	
	int _current_category_JEC_UP, _current_category_JEC_DN;
   int _current_category_BTag_UP, _current_category_BTag_DN;
   float _k_factor_EWCorr_UP, _k_factor_EWCorr_DN;
	
	vector< vector <float> > _expected_yield_PU, _expected_yield_PU_UP, _expected_yield_PU_DN;
   vector< vector <float> > _expected_yield_JEC, _expected_yield_JEC_UP, _expected_yield_JEC_DN;
   vector< vector <float> > _expected_yield_BTag, _expected_yield_BTag_UP, _expected_yield_BTag_DN;
   vector< vector <float> > _expected_yield_muR, _expected_yield_muR_UP, _expected_yield_muR_DN;
   vector< vector <float> > _expected_yield_muF, _expected_yield_muF_UP, _expected_yield_muF_DN;
   vector< vector <float> > _expected_yield_As, _expected_yield_As_UP, _expected_yield_As_DN;
   vector< vector <float> > _expected_yield_PDF, _expected_yield_PDF_UP, _expected_yield_PDF_DN;
   vector< vector <float> > _expected_yield_EWCorr, _expected_yield_EWCorr_UP, _expected_yield_EWCorr_DN;
   vector< vector <float> > _expected_yield_PythiaScale, _expected_yield_PythiaScale_UP, _expected_yield_PythiaScale_DN;
   vector< vector <float> > _expected_yield_PythiaTune, _expected_yield_PythiaTune_UP, _expected_yield_PythiaTune_DN;
	
   vector< vector <float> > _expected_yield_THU_ggH_Mu, _expected_yield_THU_ggH_Mu_UP, _expected_yield_THU_ggH_Mu_DN;
   vector< vector <float> > _expected_yield_THU_ggH_Res, _expected_yield_THU_ggH_Res_UP, _expected_yield_THU_ggH_Res_DN;
   vector< vector <float> > _expected_yield_THU_ggH_Mig01, _expected_yield_THU_ggH_Mig01_UP, _expected_yield_THU_ggH_Mig01_DN;
   vector< vector <float> > _expected_yield_THU_ggH_Mig12, _expected_yield_THU_ggH_Mig12_UP, _expected_yield_THU_ggH_Mig12_DN;
   vector< vector <float> > _expected_yield_THU_ggH_VBF2j, _expected_yield_THU_ggH_VBF2j_UP, _expected_yield_THU_ggH_VBF2j_DN;
   vector< vector <float> > _expected_yield_THU_ggH_VBF3j, _expected_yield_THU_ggH_VBF3j_UP, _expected_yield_THU_ggH_VBF3j_DN;
   vector< vector <float> > _expected_yield_THU_ggH_PT60, _expected_yield_THU_ggH_PT60_UP, _expected_yield_THU_ggH_PT60_DN;
   vector< vector <float> > _expected_yield_THU_ggH_PT120, _expected_yield_THU_ggH_PT120_UP, _expected_yield_THU_ggH_PT120_DN;
   vector< vector <float> > _expected_yield_THU_ggH_qmtop, _expected_yield_THU_ggH_qmtop_UP, _expected_yield_THU_ggH_qmtop_DN;
	
};
#endif

