// Include classes
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/include/CMS_lumi.h>
#include <iostream>


// Constructor
//====================
CMS_lumi::CMS_lumi(){}
//====================


// Destructor
//=====================
CMS_lumi::~CMS_lumi(){}
//=====================


void CMS_lumi::set_lumi( TPad* pad, float lumi )
{           

   pad->cd();

//   int align_ = 11;
//
//   float H = pad->GetWh();
//   float W = pad->GetWw();
   float l = pad->GetLeftMargin();
   float t = pad->GetTopMargin();
   float r = pad->GetRightMargin();
//   float b = pad->GetBottomMargin();
//   float histo_width = 1-l-r;
//   float histo_height = 1-t-b;

   
//======
// LUMI
//======
   
   lumi_text = Form("%.1f", lumi) + lumi_sqrt;
   
   TLatex latex;
   latex.SetNDC();
   latex.SetTextAngle(0);
   latex.SetTextColor(kBlack);    

   latex.SetTextFont(42);
   latex.SetTextAlign(31); 
   latex.SetTextSize(lumi_text_size*t);    
   latex.DrawLatex(1-r, 1-t+lumi_text_offset*t, lumi_text);


//=====
// CMS
//=====
   
   latex.SetTextAlign(11); 
   latex.SetTextSize(cms_text_size*t);
   CMS_text = Form("#font[%i]{CMS} #scale[%.2f]{#font[%i]{Preliminary}}", CMS_text_font, extra_over_CMS_text_size, extra_text_font);
   latex.DrawLatex(l, 1-t+cms_text_offset*t, CMS_text);
}

void CMS_lumi::set_lumi_combination( TPad* pad )
{

   pad->cd();

//   int align_ = 11;
//
//   float H = pad->GetWh();
//   float W = pad->GetWw();
   float l = pad->GetLeftMargin();
   float t = pad->GetTopMargin();
   float r = pad->GetRightMargin();
//   float b = pad->GetBottomMargin();
//   float histo_width = 1-l-r;
//   float histo_height = 1-t-b;

	
//======
// LUMI
//======
	
   lumi_text = "137" + lumi_sqrt;
	
   TLatex latex;
   latex.SetNDC();
   latex.SetTextAngle(0);
   latex.SetTextColor(kBlack);

   latex.SetTextFont(42);
   latex.SetTextAlign(31);
   latex.SetTextSize(lumi_text_size*t);
   latex.DrawLatex(1-r, 1-t+lumi_text_offset*t, lumi_text);


//=====
// CMS
//=====
	
   latex.SetTextAlign(11);
   latex.SetTextSize(cms_text_size*t);
   CMS_text = Form("#font[%i]{CMS} #scale[%.2f]{#font[%i]{Preliminary 2016 + 2017 + 2018}}", CMS_text_font, extra_over_CMS_text_size, extra_text_font);
   latex.DrawLatex(l, 1-t+cms_text_offset*t, CMS_text);
}
