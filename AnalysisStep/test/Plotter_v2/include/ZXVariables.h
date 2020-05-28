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
      // Updated 28/05/2020 with new yields from Elisa
      
      // 2016 Z+X Yields
      float yield_SS_4e_2016      = 13.03;
      float yield_Comb_4e_2016    = 16.13;
      // 2017 Z+X Yields
      float yield_SS_4e_2017      = 10.91;
      float yield_Comb_4e_2017    = 12.95;
      // 2018 Z+X Yields
      float yield_SS_4e_2018      = 15.99;
      float yield_Comb_4e_2018    = 19.42;
      
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
      float yield_SS_4mu_2016     = 29.66;
      float yield_Comb_4mu_2016   = 28.19;
      // 2017 Z+X Yields
      float yield_SS_4mu_2017     = 33.55;
      float yield_Comb_4mu_2017   = 33.13;
      // 2018 Z+X Yields
      float yield_SS_4mu_2018     = 52.17;
      float yield_Comb_4mu_2018   = 50.72;
       
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
      float yield_SS_2e2mu_2016   = 41.45;
      float yield_Comb_2e2mu_2016 = 44.39;
      // 2017 Z+X Yields
      float yield_SS_2e2mu_2017   = 40.96;
      float yield_Comb_2e2mu_2017 = 43.05;
      // 2018 Z+X Yields
      float yield_SS_2e2mu_2018   = 60.78;
      float yield_Comb_2e2mu_2018 = 63.87;
      
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
