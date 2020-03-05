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
#include "TLorentzVector.h"
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
#include <ZZAnalysis/AnalysisStep/interface/LeptonSFHelper.h>

// BOOLS
#define SMOOTH_ZX_FULL_RUN2_SS 0

#define NEWWP2J 0.5
#define NEWWP1J 0.7 //0.5
#define COMWPVH 0.5
#define NEWWPWH COMWPVH
#define NEWWPZH COMWPVH
#define NEWWPVH COMWPVH

using namespace std;

class Plotter: public Tree
{

public:
	
	Plotter();
	~Plotter();
   
   void MakeHistograms( TString , int );
   void MakeHistogramsZX( TString, TString , int);
   void MakeM4lZX( int );
   float calculate_K_factor( TString );
   int FindFinalState();
   int FindFinalStateZX();
   float ScaleToOSSSComb( int );
   int FindSTXSBin();
   int find_current_process( TString , int, int);
   int CountAssociatedLeptons();
   bool blind( float );
   void FillInclusive();
   void Save();
   void Delete();
   void GetHistos( TString );
   void PlotM4l( int );
   void PlotSTXS();
   void PlotPurity();
   void plot_STXS ( TString, TString);
   void plot_Purity ( TString, TString);
   void plot_1D_single( TString, TString, TString, int, int );
   void plot_1D_all_cat( TString, TString, TString);
   void plot_1D_all_fs( TString, TString, TString);
   void plot_2D_single( TString, TString, TString, int );
   void plot_2D_error_single( TString, TString, TString, int );
   void plot_2D_error_all_cat( TString , TString , TString );
   void SetBlinding( float, float);
   void SetBlinding( float, float, float, float);
   void Split_2e2mu();
   
private:

   TFile *input_file, *input_file_data;
   TTree *input_tree, *input_tree_data;
   
   LeptonSFHelper *lepSFHelper;
   float _updatedSF;

   Histograms *blinded_histos, *unblinded_histos, *combination_histos;
   
   map<TString, Histograms*> histo_map;
   
   TH1F* hCounters;
   
   Long64_t n_gen_events;
   
   float jetPt[99];
   float jetEta[99];
   float jetPhi[99];
   float jetMass[99];
   float jetQGL[99];
   float jetPgOverPq[99];
    
   float _ZZjjPt;
	
   float Pt_leading , Pt_trailing;
   float Eta_leading , Eta_trailing;
   float SIP_leading, SIP_trailing;
   float ISO_leading, ISO_trailing;
   
   Float_t KD;
   Float_t DVBFDEC;
   Float_t DVHDEC;
   Float_t D2jet;
   Float_t D1jet;
   Float_t DWH;
   Float_t DZH;
   Float_t CUSTOMCCONST, NEWWP, OLDWP;
	
   TH1F *VBF_nom;
   TH1F *VBF_tup;
   TH1F *VBF_tdn;
   TH1F *VBF_sup;
   TH1F *VBF_sdn; 
	
   
   float _blinding_lower[2], _blinding_upper[2];
   bool _merge_2e2mu;
   
   int _current_process, _current_final_state, _current_category, _current_category_stxs, _STXS_bin, _n_gen_assoc_lep;
   float _lumi, _k_factor, _SMP_signal_strength, _yield_SR, partial_sample_weight;
   double gen_sum_weights, _event_weight, _scaleToCombFactor;
   
   vector< vector <float> > _expected_yield_SR, _number_of_events_CR;
   vector<int>   gen_assoc_lep_id_;
   vector<float> _fs_ROS_SS;
};
#endif
