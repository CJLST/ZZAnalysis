#ifndef CONSTANTS_H
#define CONSTANTS_H

// C++
#include <iostream>
#include <fstream>

// ROOT
#include "TFile.h"
#include "TSpline.h"
#include "TString.h"

using namespace std;

class Constants
{
   
public:
   
   Constants();
   ~Constants();
   double getConstant( int, float );
   
private:
   
   TFile *f_4e_, *f_4mu_, *f_2e2mu_;
   TSpline3 *spline_4e_, *spline_4mu_, *spline_2e2mu_;
   
   double constant = -999;
};
#endif
