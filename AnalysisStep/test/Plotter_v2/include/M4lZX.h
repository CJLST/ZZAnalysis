#ifndef M4lZX_h
#define M4lZX_h

// ROOT
#include "TF1.h"
#include "TH1F.h"

#include "Settings.h"

using namespace std;

class M4lZX
{

public:
   
   M4lZX();
   ~M4lZX();
   void GetM4lZX( int, int, int, int, vector< vector <float> >, TH1F*, TH1F*, TH1F*, TH1F*);
   double GetM4lZX_Yields(vector< vector <float> >, int , int , int , int );
   void RenormalizeZX( int, vector< vector <float> >, TH1F* , TH1F*, TH1F* );
   void SetNormalization ( int );

private:
   
   TF1 *f_4e_comb, *f_4mu_comb, *f_2e2mu_comb;
   
   int _n_entries;
   float _bin_down, _bin_up;
   float _norm_4mu, _norm_4e, _norm_2e2mu, _norm_ZX_OS_SR_4e, _norm_ZX_OS_SR_4mu, _norm_ZX_OS_SR_2e2mu;
};
#endif
