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
	
	void PrintSystematics_JEC( TString );
	
private:

	float calculate_K_factor( TString );
   int FindFinalState();
   int find_current_production_mode( TString , int, int);
   void SumInclusive();
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
   int _current_category, _current_category_UP, _current_category_DN;
   int _n_gen_assoc_lep;
   float _k_factor, _yield_SR;
   double gen_sum_weights, _event_weight, _event_weight_UP, _event_weight_DN;
	
	vector<TString> _s_category, _s_production_mode;
   vector< vector <float> > _expected_yield, _expected_yield_UP, _expected_yield_DN;
   vector<int>   gen_assoc_lep_id_;
};
#endif

