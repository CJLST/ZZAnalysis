// Include classes
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/include/M4lZX.h>
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/include/ZXVariables.h>

using namespace std;

// Constructor
//============
M4lZX::M4lZX()
{
   _n_entries = 1000000; // Tried more entries, but it takes way too long to run
   _bin_down  = 70.;
   _bin_up    = 3000.; // Define full mass range
   
   f_4e_comb    = new TF1("f_4e_comb", "TMath::Landau(x, [0], [1])", _bin_down, _bin_up);
   f_4mu_comb   = new TF1("f_4mu_comb","TMath::Landau(x, [0], [1])", _bin_down, _bin_up);
   f_2e2mu_comb = new TF1("f_2e2mu_comb","[0]*TMath::Landau(x, [1], [2]) + [3]*TMath::Landau(x, [4], [5])", _bin_down, _bin_up);

   f_4e_comb->SetParameters(ZXVariables::ZX4e().par0, ZXVariables::ZX4e().par1);
   f_4mu_comb->SetParameters(ZXVariables::ZX4mu().par0, ZXVariables::ZX4mu().par1);
   f_2e2mu_comb->SetParameters(ZXVariables::ZX2e2mu().par0, ZXVariables::ZX2e2mu().par1, ZXVariables::ZX2e2mu().par2,
                               ZXVariables::ZX2e2mu().par3, ZXVariables::ZX2e2mu().par4, ZXVariables::ZX2e2mu().par5);

}
//=====================


//================
M4lZX::~M4lZX()
{
}
//================



//===================================================================================
void M4lZX::GetM4lZX(int n_bins, int x_min, int x_max, int category, vector< vector <float> > _norm_ZX_SS_SR, TH1F* h_4e, TH1F* h_4mu, TH1F* h_2e2mu, TH1F* h_4l )
{
   float ratio_4mu   = f_4mu_comb->Integral(x_min, x_max)/f_4mu_comb->Integral(_bin_down, _bin_up);
   float ratio_4e    = f_4e_comb->Integral(x_min, x_max)/f_4e_comb->Integral(_bin_down, _bin_up);
   float ratio_2e2mu = f_2e2mu_comb->Integral(x_min, x_max)/f_2e2mu_comb->Integral(_bin_down, _bin_up);

   _norm_4e    = (_norm_ZX_SS_SR[Settings::fs4e][category]    * ratio_4e    * ZXVariables::ZX4e().norm_Comb);
   _norm_4mu   = (_norm_ZX_SS_SR[Settings::fs4mu][category]   * ratio_4mu   * ZXVariables::ZX4mu().norm_Comb);
   _norm_2e2mu = (_norm_ZX_SS_SR[Settings::fs2e2mu][category] * ratio_2e2mu * ZXVariables::ZX2e2mu().norm_Comb);
   
   h_4mu  ->FillRandom("f_4mu_comb"  , _n_entries);
   h_4e   ->FillRandom("f_4e_comb"   , _n_entries);
   h_2e2mu->FillRandom("f_2e2mu_comb", _n_entries);
   
   h_4mu  ->Scale(_norm_4mu   / h_4mu->Integral());
   h_4e   ->Scale(_norm_4e    / h_4e->Integral());
   h_2e2mu->Scale(_norm_2e2mu / h_2e2mu->Integral());
  
   h_4l->Add(h_4mu);
   h_4l->Add(h_4e);
   h_4l->Add(h_2e2mu);
}
//===================================================================================

//===================================================================================
void M4lZX::GetM4lZXCombination( TH1F* h_4e, TH1F* h_4mu, TH1F* h_2e2mu, TH1F* h_4l)
{
   h_4mu  ->FillRandom("f_4mu_comb"  , _n_entries);
   h_4e   ->FillRandom("f_4e_comb"   , _n_entries);
   h_2e2mu->FillRandom("f_2e2mu_comb", _n_entries);
	
   h_4mu  ->Scale((ZXVariables::ZX4mu().yield_Comb_4mu_2016   + ZXVariables::ZX4mu().yield_Comb_4mu_2017   + ZXVariables::ZX4mu().yield_Comb_4mu_2018) / h_4mu->Integral());
   h_4e   ->Scale((ZXVariables::ZX4e().yield_Comb_4e_2016    + ZXVariables::ZX4e().yield_Comb_4e_2017    + ZXVariables::ZX4e().yield_Comb_4e_2018)    / h_4e->Integral());
   h_2e2mu->Scale((ZXVariables::ZX2e2mu().yield_Comb_2e2mu_2016 + ZXVariables::ZX2e2mu().yield_Comb_2e2mu_2017 + ZXVariables::ZX2e2mu().yield_Comb_2e2mu_2018) / h_2e2mu->Integral());
	
   h_4l->Add(h_4mu);
   h_4l->Add(h_4e);
   h_4l->Add(h_2e2mu);
}
//===================================================================================



//================================================================================
double M4lZX::GetM4lZX_Yields(vector< vector <float> > _norm_ZX_SS_SR, int x_min, int x_max, int final_state, int category)// [FIXME] Does not work correctly for the combination of 3 years
{
   float ratio_4mu   = f_4mu_comb->Integral(x_min, x_max)/f_4mu_comb->Integral(_bin_down, _bin_up);
   float ratio_4e    = f_4e_comb->Integral(x_min, x_max)/f_4e_comb->Integral(_bin_down, _bin_up);
   float ratio_2e2mu = f_2e2mu_comb->Integral(x_min, x_max)/f_2e2mu_comb->Integral(_bin_down, _bin_up);

   _norm_4mu   = (_norm_ZX_SS_SR[Settings::fs4mu][category]*ZXVariables::ZX4mu().norm_2016*ratio_4mu);
   _norm_4e    = (_norm_ZX_SS_SR[Settings::fs4e][category]*ZXVariables::ZX4e().yield_Comb_4e_2016*ratio_4e);
   _norm_2e2mu = (_norm_ZX_SS_SR[Settings::fs2e2mu][category]*ZXVariables::ZX2e2mu().yield_Comb_2e2mu_2016*ratio_2e2mu);
   
   if      ( final_state == Settings::fs4mu )   return _norm_4mu;
   else if ( final_state == Settings::fs4e )    return _norm_4e;
   else if ( final_state == Settings::fs2e2mu ) return _norm_2e2mu;
   else if ( final_state == Settings::fs2mu2e ) return _norm_2e2mu;
   else if ( final_state == Settings::fs4l )    return _norm_4e + _norm_4mu + _norm_2e2mu;
   else
   {
      cout << "[ERROR] Computing Z+X histogram: wrong final state: " << final_state << endl;
      abort();
   }       
}
//================================================================================



//==========================================================================
void M4lZX::RenormalizeZX( int category, int year, TH1F* histo4e_ZX, TH1F* histo4mu_ZX, TH1F* histo2e2mu_ZX)
{
   if (year == 2016 )
   {
      _norm_4e    = ZXVariables::ZX4e().norm_2016;
      _norm_4mu   = ZXVariables::ZX4mu().norm_2016;
      _norm_2e2mu = ZXVariables::ZX2e2mu().norm_2016;
   }
   else if (year == 2017 )
   {
      _norm_4e    = ZXVariables::ZX4e().norm_2017;
      _norm_4mu   = ZXVariables::ZX4mu().norm_2017;
      _norm_2e2mu = ZXVariables::ZX2e2mu().norm_2017;
   }
   else if (year == 2018 )
   {
      _norm_4e    = ZXVariables::ZX4e().norm_2018;
      _norm_4mu   = ZXVariables::ZX4mu().norm_2018;
      _norm_2e2mu = ZXVariables::ZX2e2mu().norm_2018;
   }
   
   histo4e_ZX   ->Scale( _norm_4mu );
   histo4mu_ZX  ->Scale( _norm_4e );
   histo2e2mu_ZX->Scale( _norm_2e2mu );

}
//==========================================================================

