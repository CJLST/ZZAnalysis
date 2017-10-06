// Include classes
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/include/Constants.h>

using namespace std;

// Constructor
//=======================
Constants::Constants()
{
   TString path = "../../data/cconstants/";
   
   f_4e_    = new TFile(path + "SmoothKDConstant_m4l_Dbkgkin_4e13TeV.root");
   f_4mu_   = new TFile(path + "SmoothKDConstant_m4l_Dbkgkin_4mu13TeV.root");
   f_2e2mu_ = new TFile(path + "SmoothKDConstant_m4l_Dbkgkin_2e2mu13TeV.root");

   spline_4e_    = (TSpline3*)f_4e_->Get("sp_gr_varTrue_Constant_Smooth");
   spline_4mu_   = (TSpline3*)f_4mu_->Get("sp_gr_varTrue_Constant_Smooth");
   spline_2e2mu_ = (TSpline3*)f_2e2mu_->Get("sp_gr_varTrue_Constant_Smooth");

}
//=======================


//========================
Constants::~Constants() {}
//========================


//====================================================
double Constants::getConstant(int ZZflav, float ZZMass)
{
   
   if ( abs(ZZflav) == 121*121 ) //4e
   {
      constant = spline_4e_->Eval(ZZMass);
   }
   else if ( abs(ZZflav) == 169*169 ) //4mu
   {
      constant = spline_4mu_->Eval(ZZMass);
   }
   else if ( abs(ZZflav) == 121*169 ) //2e2mu
   {
      constant = spline_2e2mu_->Eval(ZZMass);
   }
   else
   {
      cerr << "[ERROR] unknown final state." << endl;
   }
   
   return constant;
}
//====================================================
