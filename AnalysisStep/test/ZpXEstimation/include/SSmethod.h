#ifndef SSmethod_h
#define SSmethod_h

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
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TSystem.h"
#include "THStack.h"

// Include classes
#include <ZZAnalysis/AnalysisStep/test/ZpXEstimation/include/Tree.h>
#include <ZZAnalysis/AnalysisStep/test/ZpXEstimation/include/Settings.h>
#include <ZZAnalysis/AnalysisStep/test/ZpXEstimation/include/Plots.h>
#include <ZZAnalysis/AnalysisStep/interface/Category.h>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>
#include <ZZAnalysis/AnalysisStep/test/ZpXEstimation/include/FakeRates.h>
#include <ZZAnalysis/AnalysisStep/interface/bitops.h>
#include <ZZAnalysis/AnalysisStep/test/ZpXEstimation/include/CMS_lumi.h>

using namespace std;

const int num_of_processes         = Settings::num_of_processes;
const int num_of_flavours          = Settings::num_of_flavours;
const int num_of_final_states      = Settings::num_of_final_states;
const int num_of_categories        = Settings::num_of_categories;
const int num_of_regions_ss        = Settings::num_of_regions_ss;
const int num_of_eta_bins          = Settings::num_of_eta_bins;
const int num_of_fake_rates        = Settings::num_of_fake_rates;
const int num_of_z_mass_windows    = Settings::num_of_z_mass_windows;

class SSmethod: public Tree
{

public:
	
	SSmethod();
	~SSmethod();
   
   void FillFRHistos( TString );
   void FillDataMCPlots( TString );
   void MakeHistogramsZX( TString, TString );
   void SaveFRHistos( TString, bool , bool);
   void SaveDataMCHistos( TString );
   void SaveZXHistos( TString );
   void GetFRHistos( TString );
   void GetDataMCHistos( TString );
   void GetZXHistos( TString );
   void ProduceFakeRates( TString , TString input_file_data_name = "DONT_CORRECT");
   void PlotDataMC( TString , TString );
   void PlotZX( TString , TString );
   void FitZX( TString , TString );
   void Calculate_SSOS_Ratio( TString, TString, bool );
   void Set_pT_binning( int, float* );
   void SetLumi( float );
   
private:
   
   void DeclareFRHistos();
   void DeclareDataMCHistos();
   void DeclareZXHistos();
   void RemoveNegativeBins1D( TH1F* );
   void RemoveNegativeBins2D( TH2F* );
   void FillDataMCInclusive();
   void FillZXInclusive();
   void SubtractWZ();
   int find_current_process( TString );
   int FindFinalState();
   float calculate_K_factor( TString );
   bool GetVarLogX( TString );
   bool GetVarLogY( TString );
   void SavePlots( TCanvas*, TString );
   void PlotFR();
	void CorrectElectronFakeRate( TString );
	void Calculate_FR_nMissingHits( TString input_file_data_name, TGraphErrors *FR_MissingHits_graph[99][99]);
	void Fit_FRnMH_graphs( TGraphErrors *FR_MissingHits_graph[99][99] );
	void Correct_Final_FR( TString );
	int Find_Ele_pT_bin( Float_t );
	int Find_Ele_eta_bin( Float_t );
	
   TLegend* CreateLegend_FR( string , TGraphErrors*, TGraphErrors*,TGraphErrors*,TGraphErrors* );
   TLegend* CreateLegend_ZXcontr( string , TH1F*, TH1F*,TH1F*,TH1F*,TH1F* );
   TLegend* CreateLegend_ZLL( string , TH1F*, TH1F*,TH1F*,TH1F*,TH1F* );

   TFile *input_file, *input_file_data, *input_file_MC;
   TFile *fOutHistos;
   TTree *input_tree, *input_tree_data, *input_tree_MC;

   TH1F *hCounters;
   
   Long64_t n_gen_events;
   
   vector<string> _s_process, _s_flavour, _s_final_state, _s_category, _s_region;
   vector<float> _fs_ROS_SS;
   vector< vector <float> > _expected_yield_SR,_expected_yield_SR_up,_expected_yield_SR_dn, _number_of_events_CR;
   
   TString _histo_name, _histo_labels;
   
   float jetPt[99];
   float jetEta[99];
   float jetPhi[99];
   float jetMass[99];
   float jetQGL[99];
   float jetPgOverPq[99];
   
   float _pT_bins[99];
   
   float _N_SS_events[num_of_final_states][num_of_categories];
   float _N_OS_events[num_of_final_states][num_of_categories];
   
   int _current_process, _current_final_state, _current_category, _n_pT_bins, _current_pT_bin, _current_eta_bin;
   float _lumi, _yield_SR, _yield_SR_up, _yield_SR_dn , _k_factor;
   double gen_sum_weights, _event_weight, _f3, _f4;

   TH1F *histos_1D[num_of_regions_ss][num_of_processes][num_of_final_states][num_of_categories];
   
   TH1F *histos_ZX[num_of_regions_ss][num_of_processes][num_of_final_states][num_of_categories];
   
   TH2F *passing[num_of_processes][num_of_flavours], *failing[num_of_processes][num_of_flavours];
   
   TGraphErrors *FR_SS_electron_EB, *FR_SS_electron_EE, *FR_SS_muon_EB, *FR_SS_muon_EE;
   TGraphErrors *FR_SS_electron_EB_unc, *FR_SS_electron_EE_unc, *FR_SS_muon_EB_unc, *FR_SS_muon_EE_unc;
   TMultiGraph *mg_electrons, *mg_muons;
   
   vector<Float_t> vector_X[num_of_fake_rates][num_of_eta_bins][num_of_flavours];
   vector<Float_t> vector_Y[num_of_fake_rates][num_of_eta_bins][num_of_flavours];
   vector<Float_t> vector_EX[num_of_fake_rates][num_of_eta_bins][num_of_flavours];
   vector<Float_t> vector_EY[num_of_fake_rates][num_of_eta_bins][num_of_flavours];
	
	float _N_Passing[num_of_z_mass_windows][num_of_eta_bins][99];
	float _N_Failling[num_of_z_mass_windows][num_of_eta_bins][99];
	float _N_MissingHits[num_of_z_mass_windows][num_of_eta_bins][99];
	
	TF1 *Ele_FR_correction_function[num_of_eta_bins][99];

	
   
};
#endif
