#include "Riostream.h" 
#include "RooSpinZero_KD_withBkg.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <cmath>
#include "TMath.h"
#include "TH1F.h"

using namespace TMath;

ClassImp(RooSpinZero_KD_withBkg) 

RooSpinZero_KD_withBkg::RooSpinZero_KD_withBkg(const char *name, const char *title, 
				       RooAbsReal& _kd,
				       RooAbsReal& _fepspr,
				       RooAbsReal& _mu,
					   vector<TH1F*>& _histos,
					   vector<TH1F*>& _histos_bkg
				       ):
   RooAbsPdf(name,title), 
   kd("kd","kd",this,_kd),
   fepspr("fepspr","fepspr",this,_fepspr),
   mu("mu","mu",this,_mu),
   histos(_histos),
   histos_bkg(_histos_bkg)
 { 
  if (histos.size()!=3){
    coutE(InputArguments) << "RooSpinZero_KD_withBkg::RooSpinZero_KD_withBkg(" << GetName() 
			  << ") number of histograms must be 3" << endl ;
    assert(0);
  };
   int nbinsx = histos[0]->GetXaxis()->GetNbins();
   int binx_min=1,binx_max=nbinsx;
   for(int th=0;th<3;th++) IntegralT[th]=histos[th]->Integral(binx_min, binx_max);
   for(int th=0;th<1;th++) IntegralT_bkg[th]=histos_bkg[th]->Integral(binx_min, binx_max);
 }
  


 RooSpinZero_KD_withBkg::RooSpinZero_KD_withBkg(const RooSpinZero_KD_withBkg& other, const char* name) :  
   RooAbsPdf(other,name), 
   kd("kd",this,other.kd),
   fepspr("fepspr",this,other.fepspr),
   mu("mu",this,other.mu),
   histos(other.histos),
   histos_bkg(other.histos_bkg)
 { 
   for(int th=0;th<3;th++) IntegralT[th]=other.IntegralT[th];
   for(int th=0;th<1;th++) IntegralT_bkg[th]=other.IntegralT_bkg[th];
 } 



double RooSpinZero_KD_withBkg::evaluate() const 
 { 

   double intval = 0;
   double value = 0;

   int binx =  histos[0]->GetXaxis()->FindBin(kd);
   
   double T[3];
   for(int th=0;th<3;th++) T[th]=histos[th]->GetBinContent(binx);
   double T_bkg[3];
   for(int th=0;th<1;th++) T_bkg[th]=histos_bkg[th]->GetBinContent(binx);
   
   double f2 = abs(fepspr);
   double f1 = 1.0 - f2;
   double sgnf2=fepspr;
   if(fepspr!=0) sgnf2 /= f2;
   double f3 = sgnf2*sqrt(f1*f2);

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

int RooSpinZero_KD_withBkg::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{

  if (matchArgs(allVars,analVars,RooArgSet(*kd.absArg()))) return 1;

  return 0 ;
}


double RooSpinZero_KD_withBkg::analyticalIntegral(Int_t code, const char* rangeName) const
{
  double value = 0.;
  int nbinsx = histos[0]->GetNbinsX();
  double xMin = histos[0]->GetXaxis()->GetBinLowEdge(1);
  double xMax = histos[0]->GetXaxis()->GetBinUpEdge(nbinsx);
  double dx = (xMax - xMin) / ((double) nbinsx); 

   double T[3];
   for(int th=0;th<3;th++){T[th]=IntegralT[th];/* cout << T[th] << '\t' << IntegralT[th] << endl;*/};
   double T_bkg[3];
   for(int th=0;th<1;th++){T_bkg[th]=IntegralT_bkg[th];/* cout << T[th] << '\t' << IntegralT_bkg[th] << endl;*/};
   
   double f2 = abs(fepspr);
   double f1 = 1.0 - f2;
   double sgnf2=fepspr;
   if(fepspr!=0) sgnf2 /= f2;
   double f3 = sgnf2*sqrt(f1*f2);

   double xsecval = (f1*IntegralT[0] + f2*IntegralT[1] + f3*IntegralT[2])/IntegralT[0];
   if(xsecval!=0) value = ((f1*T[0] + f2*T[1] + f3*T[2])*mu)/xsecval + T_bkg[0];
   else value = T_bkg[0];

  if(value<0) value=-value;
  if(value==0) value=1;
  switch(code)
     {
     case 1: 
       {
	 return value;
       }
	 default:
       {
		   assert(0);
		   return 0;
	   }       
     }
}
