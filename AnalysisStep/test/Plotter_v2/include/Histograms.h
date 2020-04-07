#ifndef Histograms_h
#define Histograms_h

// C++
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

// ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TPaletteAxis.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TIterator.h"
#include "TROOT.h"

// Include classes
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/include/Settings.h>
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/include/M4lZX.h>
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/include/CMS_lumi.h>
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/include/Variables.h>
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/include/Cosmetics.h>


using namespace std;

const int num_of_production_modes    = Settings::num_of_production_modes;
const int num_of_processes           = Settings::num_of_processes;
const int num_of_processes_yields    = Settings::num_of_processes_yields;
const int num_of_final_states        = Settings::num_of_final_states;
const int num_of_categories          = Settings::num_of_categories;
const int num_of_STXS_categories     = Settings::num_of_STXS_categories;
const int num_of_STXS_bins           = Settings::num_of_STXS_bins;
const int num_of_1D_plot_names       = Settings::num_of_1D_plot_names;
const int num_of_2D_plot_names       = Settings::num_of_2D_plot_names;
const int num_of_2D_error_plot_names = Settings::num_of_2D_error_plot_names;


class Histograms
{
   
public:
   Histograms( double );
   Histograms( string );
   ~Histograms();
   
   void FillM4l( float, float, int, int, int );
   void FillM4lZX( float, float, int, int );
    
   void FillSTXS( float, float, int, int );
   void FillSTXSZX( float, float, int );
	
   void FillSTXSPurity( float, float, int, int );
	
   void FillMZ1( float, float, float, int, int, int );
   void FillMZ1ZX( float, float, float, int, int );
   
   void FillMZ2( float, float, float, int, int, int );
   void FillMZ2ZX( float, float, float, int, int );
   
   void FillKD( float, float, float, int, int, int );
   void FillKDZX( float, float, float, int, int );
	
   void FillDVBFDEC( float, float, float, int, int, int );
   void FillDVBFDECZX( float, float, float, int, int );
	
   void FillDVHDEC( float, float, float, int, int, int );
   void FillDVHDECZX( float, float, float, int, int );
   
   void FillD1jet( float, float, float, int, int, int );
   void FillD1jetZX( float, float, float, int, int );
   
   void FillD2jet( float, float, float, int, int, int );
   void FillD2jetZX( float, float, float, int, int );
   
   void FillDWH( float, float, float, int, int, int );
   void FillDWHZX( float, float, float, int, int );
   
   void FillDZH( float, float, float, int, int, int );
   void FillDZHZX( float, float, float, int, int );
   
   void FillDVH( float, float, float, int, int, int );
   void FillDVHZX( float, float, float, int, int );
   
   void FillMZ1vsMZ2( float, float, float, float, int, int, int );
   
   void FillVectors( float, float, float, float, float, int, float, float, float, float, float, int, int) ;
   void FillDvsM4l( float, float, float, float, int, float, float, float, float, float, float, int, int, int );
   void FillDvsM4l_ZX( float, float, float, float, int, float, int, int );
   
   void FillYields( float, float, int, int, int );
	
   void FillOthers( float, float, float, float, float, float, float, float, float, float, float, float, int, int, int, float, float, float, float, int, int, int );
   void FillOthersZX( float, float, float, float, float, float, float, float, float, float, float, float, int, int, int, float, float, float, float, int, int );
   
   void SaveHistos( string );
   void SaveYieldHistos( string );
   
   void DeleteHistos();
   void DeleteYieldsHistos();
   
   void FillInclusive();
   void FillInclusiveYields();
   
   void SmoothHistograms();
   void RenormalizeZX( int );
   
   void GetHistos( TString );
   void GetYieldsHistos( TString );
	
   void plot_Combination( TString, bool, int );
   void plot_STXS( TString );
   void plot_Purity( TString );
   void plot_1D_single( TString, TString, TString, int, int );
   void plot_1D_all_cat( TString, TString, TString );
   void plot_1D_all_fs( TString, TString, TString );
   void plot_2D_single( TString, TString, TString, int );
   void plot_2D_error_single( TString, TString, TString, int );
   void plot_2D_error_all_cat( TString , TString , TString );
   
   void FillYieldGraphs( float, float, TString);
   void PrepareYamlFiles( TString , float , float, vector<vector<float>> );
   void PrintYields( vector<vector<float>> );
   void PrintYields( float, float, vector<vector<float>> );
   void PrintLatexTables( float, float, vector<vector<float>> );
   void setColZGradient_OneColor( int , bool );
   void MakeZXShape( int, vector<vector<float>> );
   void DrawLogX( TCanvas *, int, int );
   void MakeCOLZGrey( bool );
   void SavePlots( TCanvas *, TString , TString);
   void Rebin( THStack * );
   void ChangeYaxisTitle( THStack * );
   int SetProcess( int, int );
   int SetPlotName( TString );
   float SetMassPoint( int );
   
   bool GetVarLogX( TString );
   bool GetVarLogY( TString );
	
   void Split_2e2mu();
	
   bool _merge_2e2mu;
   
   TLegend *CreateLegend( string, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F* );
   TLegend *CreateLegendVBF( string, TH1F*, TH1F*, TH1F*, TH1F*, TH1F* ,TH1F*, TH1F*, bool );
   TLegend *CreateLegendVH( string, TH1F*, TH1F*, TH1F*, TH1F*, TH1F* ,TH1F*, TH1F*, bool);
   TLegend *CreateLegendttH( string, TH1F*, TH1F*, TH1F*, TH1F*, TH1F* ,TH1F*, TH1F* );
   TLegend *CreateLegendSTXS( string , TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, bool);
   TLegend *Create2DLegend( string, TH2F*, TH2F*, TH2F* );
   TLegend *Create2DErrorLegend( string, TGraphErrors*, TGraphErrors*, TGraphErrors* );
   TLegend *Create2DLegendAllCat_KD( string, TGraphErrors*, TGraphErrors*, TGraphErrors*, TGraphErrors*, TGraphErrors*, TGraphErrors*, TGraphErrors*, TGraphErrors* );
   TLegend *Create2DLegendAllCat_DVBFDEC( string, TGraphErrors*, TGraphErrors*, TGraphErrors*, TGraphErrors*);
   TLegend *Create2DLegendAllCat_DVHDEC( string, TGraphErrors*, TGraphErrors*, TGraphErrors*, TGraphErrors*);
   
   TPaveText *CreateCutText( string, TString);
   TPaveText *CreateCatText( string, TString);
   
   TLine *CreateDashedLine(float, float, float, float);
   

private:
      
   float _lumi, _y_max;
   bool is_VBF_tagged_, is_VH_tagged_, is_ttH_tagged_, is_Djet_, is_DVH_;

   vector<double> mass_points;

   vector<TString> _s_category, _s_STXS_category, _s_STXS_bins, _s_STXS_bins_histoName, _s_category_label, _s_final_state, _s_process, _s_production_mode;
   string _histo_name, _histo_labels, _blinding;
   
   TString _graph_name, _fit_funct_name, _fs_label, _out_file_name;
   
//==========
// 1D plots
//==========
   TH1F *histos_1D[num_of_1D_plot_names][num_of_final_states][num_of_categories][num_of_processes_yields];
    
   TH1F *STXS_Categories[num_of_processes];
   TH1F *STXS_Yields[num_of_STXS_categories][num_of_processes];
   TH1F *STXS_Yields_ZX[num_of_STXS_categories];
    
   TH1F *Purity_Categories[num_of_STXS_bins];
   TH1F *Purity_Yields[num_of_STXS_categories][num_of_STXS_bins];
   float Purity_Normalization[num_of_STXS_categories];
   
   // Z+X
   TH1F *histos_1D_ZX[num_of_1D_plot_names][num_of_final_states][num_of_categories];
   TH1F *histos_1D_ZX_shape[num_of_1D_plot_names][num_of_final_states][num_of_categories];
   
//==========
// 2D plots
//==========
   TH2F *histos_2D[num_of_2D_plot_names][num_of_final_states][num_of_categories][num_of_processes];
   
//==========================
// 2D plots with mass error
//==========================
   
   // Data
   TGraphErrors *histos_2DError_data[num_of_2D_error_plot_names][num_of_final_states][num_of_categories];
   
   vector<Float_t> vector_X[num_of_2D_error_plot_names][num_of_final_states][num_of_categories];
   vector<Float_t> vector_Y[num_of_2D_error_plot_names][num_of_final_states][num_of_categories];
   vector<Float_t> vector_EX[num_of_2D_error_plot_names][num_of_final_states][num_of_categories];
   vector<Float_t> vector_EY[num_of_2D_error_plot_names][num_of_final_states][num_of_categories];
   
   // MC
   TH2F *histos_2DError[num_of_2D_error_plot_names][num_of_final_states][num_of_categories][num_of_processes];
   TH2F *histos_2DError_ZX[num_of_2D_error_plot_names][num_of_final_states][num_of_categories];
   
//==========================
// Graphs for yields vs mH
//==========================
   TGraphErrors *yields_graph[num_of_final_states][num_of_categories][num_of_production_modes];
};
#endif
