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
   void set_lumi( TPad* pad, float lumi );
      
private:
      
   TString CMS_text, lumi_text;
   TString lumi_sqrt = " fb^{-1} (13 TeV)";
   
   int CMS_text_font = 62;
   int extra_text_font = 52;

   float lumi_text_size   = 0.6;
   float lumi_text_offset = 0.2;
   float cms_text_size    = 0.75;
   float cms_text_offset  = 0.2;
         
   // ratio of CMS and extra text size
   float extra_over_CMS_text_size  = 0.76;
   
};
#endif