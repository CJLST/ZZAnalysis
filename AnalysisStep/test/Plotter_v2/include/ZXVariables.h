#ifndef ZXVariables_h
#define ZXVariables_h

// C++
#include <iostream>

using namespace std;

class ZXVariables
{
   
public:
   ZXVariables ();
   ~ZXVariables();
   
   struct ZX4e
   {
      // Updated 14/02/2020 with new yields from Elisa
      
      // 2016 Z+X Yields
      float yield_SS_4e_2016      = 13.12;
      float yield_Comb_4e_2016    = 16.22;
      // 2017 Z+X Yields
      float yield_SS_4e_2017      = 10.91;
      float yield_Comb_4e_2017    = 13.02;
      // 2018 Z+X Yields
      float yield_SS_4e_2018      = 16.05;
      float yield_Comb_4e_2018    = 19.4;
      
      float norm_2016         = yield_Comb_4e_2016/yield_SS_4e_2016;
      float norm_2017         = yield_Comb_4e_2017/yield_SS_4e_2017;
      float norm_2018         = yield_Comb_4e_2018/yield_SS_4e_2018;
      
      float norm_Comb = (yield_Comb_4e_2016 + yield_Comb_4e_2017 + yield_Comb_4e_2018) / (yield_SS_4e_2016 + yield_SS_4e_2017 + yield_SS_4e_2018);
       
      float par0 = 141.9;
      float par1 = 21.3;
   };
    
    struct ZX4mu
    {
      // 2016 Z+X Yields
      float yield_SS_4mu_2016     = 29.56;
      float yield_Comb_4mu_2016   = 28.21;
      // 2017 Z+X Yields
      float yield_SS_4mu_2017     = 33.31;
      float yield_Comb_4mu_2017   = 33.20;
      // 2018 Z+X Yields
      float yield_SS_4mu_2018     = 51.98;
      float yield_Comb_4mu_2018   = 51.33;
       
      float norm_2016         = yield_Comb_4mu_2016/yield_SS_4mu_2016;
      float norm_2017         = yield_Comb_4mu_2017/yield_SS_4mu_2017;
      float norm_2018         = yield_Comb_4mu_2018/yield_SS_4mu_2018;
       
      float norm_Comb = (yield_Comb_4mu_2016 + yield_Comb_4mu_2017 + yield_Comb_4mu_2018) / (yield_SS_4mu_2016 + yield_SS_4mu_2017 + yield_SS_4mu_2018);
       
      float par0 = 130.4;
      float par1 = 15.6;
   };
   
    struct ZX2e2mu
    {
      // 2016 Z+X Yields
      float yield_SS_2e2mu_2016   = 24.73;
      float yield_Comb_2e2mu_2016 = 44.48;
      // 2017 Z+X Yields
      float yield_SS_2e2mu_2017   = 26.11;
      float yield_Comb_2e2mu_2017 = 43.07;
      // 2018 Z+X Yields
      float yield_SS_2e2mu_2018   = 37.44;
      float yield_Comb_2e2mu_2018 = 64.06;
      
      float norm_2016         = yield_Comb_2e2mu_2016/yield_SS_2e2mu_2016;
      float norm_2017         = yield_Comb_2e2mu_2017/yield_SS_2e2mu_2017;
      float norm_2018         = yield_Comb_2e2mu_2018/yield_SS_2e2mu_2018;
       
      float norm_Comb = (yield_Comb_2e2mu_2016 + yield_Comb_2e2mu_2017 + yield_Comb_2e2mu_2018) / (yield_SS_2e2mu_2016 + yield_SS_2e2mu_2017 + yield_SS_2e2mu_2018);
        
      float par0 = 0.45;
      float par1 = 131.1;
      float par2 = 18.1;
      float par3 = 0.55;
      float par4 = 133.8;
      float par5 = 18.9;
   };
};
#endif
