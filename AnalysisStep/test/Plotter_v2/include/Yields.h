#ifndef Plotter_h
#define Plotter_h

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
#include "TPaveText.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TColor.h"

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

class Yields: public Tree
{

public:

   Yields( double );
   ~Yields();
   
   void MakeHistograms( TString, int );
   void Calculate_SS_ZX_Yields( TString, TString );
   float calculate_K_factor( TString );
   int FindFinalState();
   int FindFinalStateZX();
   int find_current_process( TString, int, int);
   int CountAssociatedLeptons();
   void FillInclusive();
   void Save();
   void Delete();
   void GetHistos( TString );
   void FillGraphs( TString , float, float, TString );
   void PrepareYamlFiles( TString, TString , float , float );
   void Print( TString );
   void Print( TString, float, float );
   void PrintLatexTables( TString, float, float );
   void Split_2e2mu();
   void ProduceDataROOTFiles( TString, TString );
   
private:

   TFile *input_file, *input_file_data;
   TTree *input_tree, *input_tree_data;

   Histograms *yields_histos;
   
   map<TString, Histograms*> histo_map;
   
   TColor *_tclr;

   TH1F* hCounters;
   
   Long64_t n_gen_events;
   
   bool _merge_2e2mu;
   
   float jetPt[99];
   float jetEta[99];
   float jetPhi[99];
   float jetMass[99];
   float jetQGL[99];
   float jetPgOverPq[99];
   
   Float_t KD;
   Float_t D2jet;
   Float_t D1jet;
   Float_t DWH;
   Float_t DZH;

   int _current_process, _current_final_state, _current_category, _n_gen_assoc_lep;
   float _lumi, _k_factor, _SMP_signal_strength, _yield_SR, partial_sample_weight;
   double gen_sum_weights, _event_weight;
   
   vector< vector <float> > _expected_yield_SR, _number_of_events_CR;
   vector<TString> _s_category, _s_final_state;
   
   vector<int>   gen_assoc_lep_id_;
   vector<float> _fs_ROS_SS;
   vector<float> cb_SS;
   
   TFile *data_root_file;
   TTree *data_obs;
   Double_t mass4l;
   Double_t kd;
   vector<float> _mass[num_of_final_states][num_of_categories];
   vector<float> _kd[num_of_final_states][num_of_categories];
   
};
#endif
