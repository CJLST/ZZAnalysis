#ifndef CMS_lumi_h
#define CMS_lumi_h

#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"

using namespace std;

class CMS_lumi
{

public:
      
   CMS_lumi();
   ~CMS_lumi();
   void set_lumi(TPad* pad, float lumi, int iPosX = 10 );
      
private:
      
   TString lumiText;
      
   TString cmsText   = "CMS";
   float cmsTextFont = 61; // default is helvetic-bold

   bool writeExtraText = true;
   TString extraText   = "Preliminary";
   float extraTextFont = 52;  // default is helvetica-italics

   float lumiTextSize   = 0.6;
   float lumiTextOffset = 0.2;
   float cmsTextSize    = 0.75;
   float cmsTextOffset  = 0.1;  // only used in outOfFrame version
   
   float relPosX    = 0.045;
   float relPosY    = 0.035;
   float relExtraDY = 1.2;
      
   // ratio of "CMS" and extra text size
   float extraOverCmsTextSize  = 0.76;
   
   TString lumi_sqrtS = " fb^{-1} (13 TeV)";
      
   bool drawLogo = false;    
};

#endif
