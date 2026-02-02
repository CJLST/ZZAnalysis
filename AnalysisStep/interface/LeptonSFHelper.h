#ifndef LEPTONSFHELPER_H
#define LEPTONSFHELPER_H

#include <string>
#include <iostream>
#include <vector>
#include <utility>

#include <cmath>
#include "TString.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH2D.h"

class LeptonSFHelper
{

 public:

  LeptonSFHelper(int year, std::string const &data_tag);
  ~LeptonSFHelper();

  /// return pair<SF, SFError>
  std::pair<float, float> getSF (int flav, float pt, float eta, float SCeta, float phi, bool isCrack) const;
   
 private:
  int theYear;
  std::string theDataTag;
  TH2F *h_Ele_ID;
  TH2F *h_Ele_ID_HoleBPix;
  TH2F *h_Ele_ID_Gap;
  TH2F *h_Ele_Reco_lowPt;
  TH2F *h_Ele_Reco_midPt;
  TH2F *h_Ele_Reco_highPt;
  TH2D *h_Mu_SF;
  TH2D *h_Mu_Unc;
  bool isPostBPix_;
};

#endif
