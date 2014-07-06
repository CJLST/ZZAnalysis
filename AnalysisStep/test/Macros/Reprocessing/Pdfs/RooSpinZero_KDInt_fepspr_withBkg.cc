#include "Riostream.h" 
#include "RooSpinZero_KDInt_fepspr_withBkg.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <cmath>
#include "TMath.h"
#include "TH2F.h"

using namespace TMath;

ClassImp(RooSpinZero_KDInt_fepspr_withBkg) 

  RooSpinZero_KDInt_fepspr_withBkg::RooSpinZero_KDInt_fepspr_withBkg(const char *name, const char *title, 
				       RooAbsReal& _kd,
				       RooAbsReal& _kdint,
				       RooAbsReal& _fepspr,
				       RooAbsReal& _mu,
					   vector<TH2F*>& _histos,
					   vector<TH2F*>& _histos_bkg
				       ):
   RooAbsPdf(name,title), 
   kd("kd","kd",this,_kd),
   kdint("kdint","kdint",this,_kdint),
   fepspr("fepspr","fepspr",this,_fepspr),
   mu("mu","mu",this,_mu),
   histos(_histos),
   histos_bkg(_histos_bkg)
 { 
  if (histos.size()!=3){
    coutE(InputArguments) << "RooSpinZero_KDInt_fepspr_withBkg::RooSpinZero_KDInt_fepspr_withBkg(" << GetName() 
			  << ") number of histograms must be 3" << endl ;
    assert(0);
  };
   int nbinsx = histos[0]->GetXaxis()->GetNbins();
   int nbinsy = histos[0]->GetYaxis()->GetNbins();
   int binx_min=1,binx_max=nbinsx,biny_min=1,biny_max=nbinsy;
   for(int th=0;th<3;th++) IntegralT[th]=histos[th]->Integral(binx_min, binx_max, biny_min, biny_max);
   for(int th=0;th<1;th++) IntegralT_bkg[th]=histos_bkg[th]->Integral(binx_min, binx_max, biny_min, biny_max);

/*   cout << _kd.GetName() << '\t' << _kdint.GetName() << endl;
   for(int th=0;th<3;th++) cout << IntegralT[th] << '\t';
   for(int th=0;th<1;th++) cout << IntegralT_bkg[th] << '\t';
   cout << endl;
*/ }


 RooSpinZero_KDInt_fepspr_withBkg::RooSpinZero_KDInt_fepspr_withBkg(const RooSpinZero_KDInt_fepspr_withBkg& other, const char* name) :  
   RooAbsPdf(other,name), 
   kd("kd",this,other.kd),
   kdint("kdint",this,other.kdint),
   fepspr("fepspr",this,other.fepspr),
   mu("mu",this,other.mu),
   histos(other.histos),
   histos_bkg(other.histos_bkg)
 { 
   for(int th=0;th<3;th++) IntegralT[th]=other.IntegralT[th];
   for(int th=0;th<1;th++) IntegralT_bkg[th]=other.IntegralT_bkg[th];
 } 



 Double_t RooSpinZero_KDInt_fepspr_withBkg::evaluate() const 
 { 

   double intval = 0;
   double value = 0;

   int binx =  histos[0]->GetXaxis()->FindBin(kd);
   int biny =  histos[0]->GetYaxis()->FindBin(kdint);
   
   double T[3];
   for(int th=0;th<3;th++) T[th]=histos[th]->GetBinContent(binx,biny);
   double T_bkg[3];
   for(int th=0;th<1;th++) T_bkg[th]=histos_bkg[th]->GetBinContent(binx,biny);
   
   double f2 = abs(fepspr);
   double f1 = 1.0 - f2;
   double sgnf2=fepspr;
   if(fepspr!=0) sgnf2 /= f2;
   double f3 = -sgnf2*f1*f2;
   f1 *= f1;
   f2 *= f2;

   double xsecval = (f1*IntegralT[0] + f2*IntegralT[1] + f3*IntegralT[2])/IntegralT[0];
   if(xsecval!=0) value = ((f1*T[0] + f2*T[1] + f3*T[2])*mu)/xsecval + T_bkg[0];
   else value = T_bkg[0];
   if(xsecval!=0) intval = ((f1*IntegralT[0] + f2*IntegralT[1] + f3*IntegralT[2])/xsecval)*mu + IntegralT_bkg[0];
   else intval = IntegralT_bkg[0];

   if (value==0) value=1.0e-27;
   if (value<0 && intval<0) value=-value;
   if (intval==0) value=1.0e-20;
   if (value<0 && intval>0) value=1.0e-20;
   if (value>0 && intval<0) value=1.0e-20;

   return value ; 
   
 } 

Int_t RooSpinZero_KDInt_fepspr_withBkg::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{

//  cout << "no 4" << endl;
  if (matchArgs(allVars,analVars,kd,kdint)) return 3 ;
//  cout << "no 3" << endl;
  if (matchArgs(allVars,analVars,kd)) return 1 ;
//  cout << "no 2" << endl;
  if (matchArgs(allVars,analVars,kdint)) return 2 ;
//  cout << "no 1" << endl;

  return 0 ;
}

Double_t RooSpinZero_KDInt_fepspr_withBkg::analyticalIntegral(Int_t code, const char* /*rangeName*/) const
{

  int nbinsx = histos[0]->GetXaxis()->GetNbins();
  double xMin = histos[0]->GetXaxis()->GetBinLowEdge(1);
  double xMax = histos[0]->GetXaxis()->GetBinUpEdge(nbinsx);
  double dx = (xMax - xMin) / nbinsx; 

  
  int nbinsy = histos[0]->GetYaxis()->GetNbins();
  double yMin = histos[0]->GetYaxis()->GetBinLowEdge(1);
  double yMax = histos[0]->GetYaxis()->GetBinUpEdge(nbinsy);
  double dy = (yMax - yMin) / nbinsy;
  int binx_min=1,binx_max=nbinsx,biny_min=1,biny_max=nbinsy;
  double T[3]={0};
  double T_bkg[3]={0};
  double value = 0;
 
       // integrate out kd, depend on kdint
     if(code == 1){

	 int biny = histos[0]->GetYaxis()->FindBin(kdint);
	 biny_min=biny;
	 biny_max=biny;
     for(int th=0;th<3;th++) T[th]=histos[th]->Integral(binx_min, binx_max, biny_min, biny_max);
     for(int th=0;th<1;th++) T_bkg[th]=histos_bkg[th]->Integral(binx_min, binx_max, biny_min, biny_max);
       } // integrate out  kdint, depend on kd
     else if(code == 2){

	 int binx = histos[0]->GetXaxis()->FindBin(kd);
	 binx_min=binx;
	 binx_max=binx;
     for(int th=0;th<3;th++) T[th]=histos[th]->Integral(binx_min, binx_max, biny_min, biny_max);
     for(int th=0;th<1;th++) T_bkg[th]=histos_bkg[th]->Integral(binx_min, binx_max, biny_min, biny_max);
	   }
	 else if(code == 3){

		 for(int th=0;th<3;th++){T[th]=IntegralT[th];/* cout << T[th] << '\t' << IntegralT[th] << endl;*/};
		 for(int th=0;th<1;th++){T_bkg[th]=IntegralT_bkg[th];/* cout << T[th] << '\t' << IntegralT_bkg[th] << endl;*/};
	   }
	 else{
		   return 0;
	   };
   
   double f2 = abs(fepspr);
   double f1 = 1.0 - f2;
   double sgnf2=fepspr;
   if(fepspr!=0) sgnf2 /= f2;
   double f3 = -sgnf2*f1*f2;
   f1 *= f1;
   f2 *= f2;

   double xsecval = (f1*IntegralT[0] + f2*IntegralT[1] + f3*IntegralT[2])/IntegralT[0];
   if(xsecval!=0) value = ((f1*T[0] + f2*T[1] + f3*T[2])*mu)/xsecval + T_bkg[0];
   else value = T_bkg[0];

//  for(int o=0;o<6;o++) cout << T[o] << '\t';
//  cout << value << endl;
  if(value==0) cout << "WARNING: INTEGRAL 0" << endl;
  if(value<0) cout << "WARNING: INTEGRAL LT 0" << endl;
  if(value==0) value=1;
  if(value<0) value=-value;

  return value;
}

