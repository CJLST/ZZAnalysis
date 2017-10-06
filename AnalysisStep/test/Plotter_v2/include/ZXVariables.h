#ifndef ZXVariables_h
#define ZXVariables_h

// C++
#include <iostream>

using namespace std;


static float ratio_4e    = 21.1/19.5391;
static float ratio_4mu   = 34.4/35.9635;
static float ratio_2e2mu = 59.9/55.4615;


class ZXVariables
{
   
public:
   ZXVariables ();
   ~ZXVariables();
   
   float yield_SS_4e    = 19.5391;
   float yield_SS_4mu   = 35.9635;
   float yield_SS_2e2mu = 55.4615;
   
   struct ZX4e
   {            
      float norm_inclusive        = ratio_4e;
      float norm_untagged         = ratio_4e;
      float norm_VBF_1j_tagged    = ratio_4e;
      float norm_VBF_2j_tagged    = ratio_4e;
      float norm_VH_lepton_tagged = ratio_4e;
      float norm_VH_hadron_tagged = ratio_4e;
      float norm_ttH_tagged       = ratio_4e;
      float norm_VH_MET_tagged    = ratio_4e;
       
      float par0 = 141.9;
      float par1 = 21.3;
   };
    
    struct ZX4mu
    {
      float norm_inclusive        = ratio_4mu;
      float norm_untagged         = ratio_4mu;
      float norm_VBF_1j_tagged    = ratio_4mu;
      float norm_VBF_2j_tagged    = ratio_4mu;
      float norm_VH_lepton_tagged = ratio_4mu;
      float norm_VH_hadron_tagged = ratio_4mu;
      float norm_ttH_tagged       = ratio_4mu;
      float norm_VH_MET_tagged    = ratio_4mu;
       
      float par0 = 130.4;
      float par1 = 15.6;
   };
   
    struct ZX2e2mu
    {
      float norm_inclusive        = ratio_2e2mu;
      float norm_untagged         = ratio_2e2mu;
      float norm_VBF_1j_tagged    = ratio_2e2mu;
      float norm_VBF_2j_tagged    = ratio_2e2mu;
      float norm_VH_lepton_tagged = ratio_2e2mu;
      float norm_VH_hadron_tagged = ratio_2e2mu;
      float norm_ttH_tagged       = ratio_2e2mu;
      float norm_VH_MET_tagged    = ratio_2e2mu;
        
      float par0 = 0.45;
      float par1 = 131.1;
      float par2 = 18.1;
      float par3 = 0.55;
      float par4 = 133.8;
      float par5 = 18.9;
   };
};
#endif