#ifndef WEIGHT_UTILS_H
#define WEIGHT_UTILS_H


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <iomanip>

#include <TSystem.h>
#include <TROOT.h>
#include "TH1F.h"
#include "TFile.h"

using namespace std;


class WeightUtils
{

 public:
  
  WeightUtils(const TString& weightFile)
    {
      fIn = TFile::Open(weightFile, "READ");
      h_weight_ggHToTot_7TeV = (TH1F*)fIn->Get("h_weight_ggHToTot_7TeV");
      h_weight_ggHToTot_8TeV = (TH1F*)fIn->Get("h_weight_ggHToTot_8TeV");
      h_weight_totToGGH_7TeV = (TH1F*)fIn->Get("h_weight_totToGGH_7TeV");
      h_weight_totToGGH_8TeV = (TH1F*)fIn->Get("h_weight_totToGGH_8TeV");
    }
  
  ~WeightUtils() {}
  
  Double_t scaleToTotXS(Double_t mH, Double_t sqrts)
  {
    if (sqrts == 7) return h_weight_ggHToTot_7TeV->GetBinContent(h_weight_ggHToTot_7TeV->FindBin(mH));
    else if (sqrts == 8) return h_weight_ggHToTot_8TeV->GetBinContent(h_weight_ggHToTot_8TeV->FindBin(mH));
    else {
      cout << "%% scaleToTotXS: sqrts must be either 7 or 8, returning 0" << endl;
      return 0;
    }
    
  }// end scaleToTotXS


  Double_t scaleToGGHXS(Double_t mH, Double_t sqrts)
  {
    if (sqrts == 7) return h_weight_totToGGH_7TeV->GetBinContent(h_weight_totToGGH_7TeV->FindBin(mH));
    else if (sqrts == 8) return h_weight_totToGGH_8TeV->GetBinContent(h_weight_totToGGH_8TeV->FindBin(mH));
    else {
      cout << "%% scaleToGGHXS: sqrts must be either 7 or 8, returning 0" << endl;
      return 0;
    }
    
  }// end scaleToGGHXS


 private:
  
  TFile* fIn;
  TH1F* h_weight_ggHToTot_7TeV;
  TH1F* h_weight_ggHToTot_8TeV;
  TH1F* h_weight_totToGGH_7TeV;
  TH1F* h_weight_totToGGH_8TeV;
  
  
}; // class weightUtils

#endif
