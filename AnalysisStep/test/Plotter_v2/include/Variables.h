#ifndef Variables_h
#define Variables_h

// C++
#include <iostream>
#include <vector>

// ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"


using namespace std;

class Variables
{
   
public:
   Variables ();
   ~Variables();

   
//=====
// M4l
//=====
   
   struct M4lMain
   {
      TString var_X_label       = "m_{4#font[12]{l}} (GeV)";
      TString var_X_label_4e    = "m_{4#font[12]{e}} (GeV)";
      TString var_X_label_4mu   = "m_{4#font[12]{#mu}} (GeV)";
      TString var_X_label_2e2mu = "m_{2#font[12]{e}2#font[12]{#mu}} (GeV)";
      TString var_Y_label       = "Events / 4 GeV";
      TString var_cut_label     = "";
      Int_t var_N_bin = 233;
      Float_t var_min = 70;
      Float_t var_max = 1002;
      Bool_t var_log_x = 1;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 1;
   };
   
   struct M4lYields
   {
      TString var_X_label = "";
      TString var_Y_label = "";
      TString var_cut_label = "";
      Int_t var_N_bin = 930;
      Float_t var_min = 70;
      Float_t var_max = 1000;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 1;
   };


   struct M4lMainZoomed
   {
      TString var_X_label       = "m_{4#font[12]{l}} (GeV)";
      TString var_X_label_4e    = "m_{4#font[12]{e}} (GeV)";
      TString var_X_label_4mu   = "m_{4#font[12]{#mu}} (GeV)";
      TString var_X_label_2e2mu = "m_{2#font[12]{e}2#font[12]{#mu}} (GeV)";
      TString var_Y_label       = "Events / 2 GeV";
      TString var_cut_label     = "";
      Int_t var_N_bin = 50;
      Float_t var_min = 70;
      Float_t var_max = 170;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 1;
   };
   
   struct M4lMainHighMass
   {
      TString var_X_label = "m_{4#font[12]{l}} (GeV)";
      TString var_Y_label = "Events / 10 GeV";
      TString var_cut_label = "";
      Int_t var_N_bin = 83;
      Float_t var_min = 170;
      Float_t var_max = 1000;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 1;
   };
   
   
//=====
// MZ1
//=====
   
   struct MZ1
   {
      TString var_X_label = "m_{Z_{1}} (GeV)";
      TString var_Y_label = "Events / 2 GeV";
      TString var_cut_label = "";
      Int_t var_N_bin = 40;
      Float_t var_min = 40;
      Float_t var_max = 120;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 11;
      Int_t rebinningDYTTbar = 5;
   };
   
   
   struct MZ1_M4L118130
   {
      TString var_X_label = "m_{Z_{1}} (GeV)";
      TString var_Y_label = "Events / 2 GeV";
      TString var_cut_label = "118 < m_{4#font[12]{l}} < 130 GeV";
      Int_t var_N_bin = 40;
      Float_t var_min = 40;
      Float_t var_max = 120;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 11;
      Int_t rebinningDYTTbar = 2;
      
      Float_t cut_d = 118.;
      Float_t cut_u = 130.;
   };
   
   
//=====
// MZ2
//=====
   
   struct MZ2
   {
      TString var_X_label = "m_{Z_{2}} (GeV)";
      TString var_Y_label = "Events / 2 GeV";
      TString var_cut_label = "";
      Int_t var_N_bin = 54;
      Float_t var_min = 12;
      Float_t var_max = 120;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 11;
      Int_t rebinningDYTTbar = 5;
   };
   
   
   struct MZ2_M4L118130
   {
      TString var_X_label = "m_{Z_{2}} (GeV)";
      TString var_Y_label = "Events / 2 GeV";
      TString var_cut_label = "118 < m_{4#font[12]{l}} < 130 GeV";
      Int_t var_N_bin = 54;
      Float_t var_min = 12;
      Float_t var_max = 120;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 5;
      
      Float_t cut_d = 118.;
      Float_t cut_u = 130.;
   };
   
   
//====
// KD 
//====
   
   struct KD
   {
      TString var_X_label = "D_{bkg}^{kin}";
      TString var_Y_label = "Events / 0.05";
      TString var_cut_label = "";
      Int_t var_N_bin = 20;
      Float_t var_min = 0;
      Float_t var_max = 1;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 4;
   };
   
   struct KD_M4L118130
   {
      TString var_X_label = "D_{bkg}^{kin}";
      TString var_Y_label = "Events / 0.1";
      TString var_cut_label = "118 < m_{4#font[12]{l}} < 130 GeV";
      Int_t var_N_bin = 20;
      Float_t var_min = 0;
      Float_t var_max = 1;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 4;
      
      Float_t cut_d = 118.;
      Float_t cut_u = 130.;
   };
   
   struct D1jet
   {
      TString var_X_label = "D_{1jet}";
      TString var_Y_label = "Events / 0.05";
      TString var_cut_label = "";
      Int_t var_N_bin = 20;
      Float_t var_min = 0;
      Float_t var_max = 1;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 4;
   };
   
   struct D1jet_M4L118130
   {
      TString var_X_label = "D_{1jet}";
      TString var_Y_label = "Events / 0.1";
      TString var_cut_label = "118 < m_{4#font[12]{l}} < 130 GeV";
      Int_t var_N_bin = 10;
      Float_t var_min = 0;
      Float_t var_max = 1;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 4;
      
      Float_t cut_d = 118.;
      Float_t cut_u = 130.;
   };
   
   struct D2jet
   {
      TString var_X_label = "D_{2jet}";
      TString var_Y_label = "Events / 0.05";
      TString var_cut_label = "";
      Int_t var_N_bin = 20;
      Float_t var_min = 0;
      Float_t var_max = 1;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 4;
   };
   
   struct D2jet_M4L118130
   {
      TString var_X_label = "D_{2jet}";
      TString var_Y_label = "Events / 0.1";
      TString var_cut_label = "118 < m_{4#font[12]{l}} < 130 GeV";
      Int_t var_N_bin = 10;
      Float_t var_min = 0;
      Float_t var_max = 1;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 4;
      
      Float_t cut_d = 118.;
      Float_t cut_u = 130.;
   };
   
   struct DWH
   {
      TString var_X_label = "D_{WH}";
      TString var_Y_label = "Events / 0.05";
      TString var_cut_label = "";
      Int_t var_N_bin = 20;
      Float_t var_min = 0;
      Float_t var_max = 1;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 4;
   };
   
   struct DWH_M4L118130
   {
      TString var_X_label = "D_{WH}";
      TString var_Y_label = "Events / 0.1";
      TString var_cut_label = "118 < m_{4#font[12]{l}} < 130 GeV";
      Int_t var_N_bin = 10;
      Float_t var_min = 0;
      Float_t var_max = 1;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 4;
      
      Float_t cut_d = 118.;
      Float_t cut_u = 130.;
   };
   
   struct DZH
   {
      TString var_X_label = "D_{ZH}";
      TString var_Y_label = "Events / 0.05";
      TString var_cut_label = "";
      Int_t var_N_bin = 20;
      Float_t var_min = 0;
      Float_t var_max = 1;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 4;
   };
   
   struct DZH_M4L118130
   {
      TString var_X_label = "D_{ZH}";
      TString var_Y_label = "Events / 0.1";
      TString var_cut_label = "118 < m_{4#font[12]{l}} < 130 GeV";
      Int_t var_N_bin = 10;
      Float_t var_min = 0;
      Float_t var_max = 1;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 4;
      
      Float_t cut_d = 118.;
      Float_t cut_u = 130.;
   };
   
   struct DVH
   {
      TString var_X_label = "D_{VH}";
      TString var_Y_label = "Events / 0.05";
      TString var_cut_label = "";
      Int_t var_N_bin = 20;
      Float_t var_min = 0;
      Float_t var_max = 1;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 4;
   };
   
   struct DVH_M4L118130
   {
      TString var_X_label = "D_{VH}";
      TString var_Y_label = "Events / 0.1";
      TString var_cut_label = "118 < m_{4#font[12]{l}} < 130 GeV";
      Int_t var_N_bin = 10;
      Float_t var_min = 0;
      Float_t var_max = 1;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 4;
      
      Float_t cut_d = 118.;
      Float_t cut_u = 130.;
   };
     
   
//=============
// MZ1 vs MZ2
//=============   
 
   struct MZ1vsMZ2
   {
      TString var_X_label = "m_{Z_{1}} (GeV)";
      TString var_Y_label = "m_{Z_{2}} (GeV)";
      TString var_cut_label = "";
      Int_t var_X_N_bin = 200;
      Float_t var_X_min = 40;
      Float_t var_X_max = 120;
      Int_t var_Y_N_bin = 220;
      Float_t var_Y_min = 12;
      Float_t var_Y_max = 120;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
   };
   
   struct MZ1vsMZ2_M4L118130
   {
      TString var_X_label = "m_{Z_{1}} (GeV)";
      TString var_Y_label = "m_{Z_{2}} (GeV)";
      TString var_cut_label = "";
      Int_t var_X_N_bin = 200;
      Float_t var_X_min = 40;
      Float_t var_X_max = 120;
      Int_t var_Y_N_bin = 220;
      Float_t var_Y_min = 12;
      Float_t var_Y_max = 120;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
      
      Float_t cut_d = 118.;
      Float_t cut_u = 130.;
   };
   

//===========
// KD vs M4l
//===========
   
   struct KDvsM4l
   {
      TString var_X_label = "m_{4#font[12]{l}} (GeV)";
      TString var_Y_label = "D_{bkg}^{kin}";
      TString var_cut_label = "";
      Int_t var_X_N_bin = 450;
      Float_t var_X_min = 100;
      Float_t var_X_max = 1000;
      Int_t var_Y_N_bin = 20;
      Float_t var_Y_min = 0;
      Float_t var_Y_max = 1;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
   };
   
   struct KDvsM4lZoomed
   {
      TString var_X_label = "m_{4#font[12]{l}} (GeV)";
      TString var_Y_label = "D_{bkg}^{kin}";
      TString var_cut_label = "";
      Int_t var_X_N_bin = 35;
      Float_t var_X_min = 100;
      Float_t var_X_max = 170;
      Int_t var_Y_N_bin = 20;
      Float_t var_Y_min = 0;
      Float_t var_Y_max = 1;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
   };
   
   struct KDvsM4lHighMass
   {
      TString var_X_label = "m_{4#font[12]{l}} (GeV)";
      TString var_Y_label = "D_{bkg}^{kin}";
      TString var_cut_label = "";
      Int_t var_X_N_bin = 415;
      Float_t var_X_min = 170;
      Float_t var_X_max = 1000;
      Int_t var_Y_N_bin = 20;
      Float_t var_Y_min = 0;
      Float_t var_Y_max = 1;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
   };
   
   struct D1jetvsM4lZoomed
   {
      TString var_X_label = "m_{4#font[12]{l}} (GeV)";
      TString var_Y_label = "D_{1jet}";
      TString var_cut_label = "";
      Int_t var_X_N_bin = 35;
      Float_t var_X_min = 100;
      Float_t var_X_max = 170;
      Int_t var_Y_N_bin = 20;
      Float_t var_Y_min = 0;
      Float_t var_Y_max = 1;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
   };
   
   struct D2jetvsM4lZoomed
   {
      TString var_X_label = "m_{4#font[12]{l}} (GeV)";
      TString var_Y_label = "D_{2jet}";
      TString var_cut_label = "";
      Int_t var_X_N_bin = 35;
      Float_t var_X_min = 100;
      Float_t var_X_max = 170;
      Int_t var_Y_N_bin = 20;
      Float_t var_Y_min = 0;
      Float_t var_Y_max = 1;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
   };
   
   struct DWHvsM4lZoomed
   {
      TString var_X_label = "m_{4#font[12]{l}} (GeV)";
      TString var_Y_label = "D_{WH}";
      TString var_cut_label = "";
      Int_t var_X_N_bin = 35;
      Float_t var_X_min = 100;
      Float_t var_X_max = 170;
      Int_t var_Y_N_bin = 20;
      Float_t var_Y_min = 0;
      Float_t var_Y_max = 1;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
   };
   
   struct DZHvsM4lZoomed
   {
      TString var_X_label = "m_{4#font[12]{l}} (GeV)";
      TString var_Y_label = "D_{ZH}";
      TString var_cut_label = "";
      Int_t var_X_N_bin = 35;
      Float_t var_X_min = 100;
      Float_t var_X_max = 170;
      Int_t var_Y_N_bin = 20;
      Float_t var_Y_min = 0;
      Float_t var_Y_max = 1;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
   };
   
   struct DVHvsM4lZoomed
   {
      TString var_X_label = "m_{4#font[12]{l}} (GeV)";
      TString var_Y_label = "D_{VH}";
      TString var_cut_label = "";
      Int_t var_X_N_bin = 35;
      Float_t var_X_min = 100;
      Float_t var_X_max = 170;
      Int_t var_Y_N_bin = 20;
      Float_t var_Y_min = 0;
      Float_t var_Y_max = 1;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 0;
      Int_t varLegPos = 33;
   };
   
   struct Pt4l
   {
      TString var_X_label = "p_{T}^{4#font[12]{l}} (GeV)";
      TString var_Y_label = "Events / 10";
      TString var_cut_label = "";
      Int_t var_N_bin = 40;
      Float_t var_min = 0;
      Float_t var_max = 400;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 33;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 2;
   };
   
   struct Eta4l
   {
      TString var_X_label = "#eta^{4#font[12]{l}}";
      TString var_Y_label = "Events / 0.5";
      TString var_cut_label = "";
      Int_t var_N_bin = 20;
      Float_t var_min = -10;
      Float_t var_max = 10;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 11;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 2;
   };
   
   struct NExtraLep
   {
      TString var_X_label = "number of additional leptons";
      TString var_Y_label = "Events";
      TString var_cut_label = "";
      Int_t var_N_bin = 6;
      Float_t var_min = 0;
      Float_t var_max = 6;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 1;
      Int_t restrict_count_var = 2;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 11;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 1;
   };
   
   struct NJets
   {
      TString var_X_label = "number of jets";
      TString var_Y_label = "Events";
      TString var_cut_label = "";
      Int_t var_N_bin = 17;
      Float_t var_min = 0;
      Float_t var_max = 17;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 1;
      Int_t restrict_count_var = 4;
      Float_t var_min_factor = 20000;
      Int_t var_CMS_pos = 11;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 1;
   };
   
   struct NJetsBTagged
   {
      TString var_X_label = "number of tagged b jets";
      TString var_Y_label = "Events";
      TString var_cut_label = "";
      Int_t var_N_bin = 8;
      Float_t var_min = 0;
      Float_t var_max = 8;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 1;
      Int_t restrict_count_var = 2;
      Float_t var_min_factor = 2000;
      Int_t var_CMS_pos = 11;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 1;
   };
   
   struct M4l_100180_HighKD
   {
      TString var_X_label = "m_{4#font[12]{l}} (GeV)";
      TString var_Y_label = "Events / 3 GeV";
      TString var_cut_label = "D_{bkg}^{kin} > 0.5";
      Int_t var_N_bin = 27;
      Float_t var_min = 100;
      Float_t var_max = 181;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 100000;
      Int_t var_CMS_pos = 11;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 1;
   };
   
   struct M4l_110150_HighKD
   {
      TString var_X_label = "m_{4#font[12]{l}} (GeV)";
      TString var_Y_label = "Events / 2 GeV";
      TString var_cut_label = "D_{bkg}^{kin} > 0.5";
      Int_t var_N_bin = 20;
      Float_t var_min = 110;
      Float_t var_max = 150;
      Bool_t var_log_x = 0;
      Bool_t var_log_y = 0;
      Int_t restrict_count_var = 0;
      Float_t var_min_factor = 0;
      Int_t var_CMS_pos = 11;
      Int_t varLegPos = 33;
      Int_t rebinningDYTTbar = 1;
   };
   
};
#endif
