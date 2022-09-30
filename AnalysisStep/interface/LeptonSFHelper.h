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

enum SFsyst {central = 0, up = 1, down = 2};

class LeptonSFHelper
{

 public:

  LeptonSFHelper(bool preVFP);
  ~LeptonSFHelper();
  
  float getSF (int year, int flav, float pt, float eta, float SCeta, bool isCrack) const;
  float getSFError (int year, int flav, float pt, float eta, float SCeta, bool isCrack) const;
   
 private:
   TFile *root_file;
   
   // Electron SF map histograms
   TH2F *h_Ele_notCracks_2016, *h_Ele_notCracks_2017, *h_Ele_notCracks_2018;
   TH2F *h_Ele_Cracks_2016, *h_Ele_Cracks_2017, *h_Ele_Cracks_2018;
   
   TH2F *h_Ele_Reco_highPT_2016, *h_Ele_Reco_highPT_2017, *h_Ele_Reco_highPT_2018;
   TH2F *h_Ele_Reco_lowPT_2016, *h_Ele_Reco_lowPT_2017, *h_Ele_Reco_lowPT_2018;
   
   // Muons SF map histograms
   TH2D *h_Mu_SF_2016, *h_Mu_SF_2017, *h_Mu_SF_2018;
   TH2D *h_Mu_Unc_2016, *h_Mu_Unc_2017, *h_Mu_Unc_2018;

};


// Add bindings to call from python
extern "C" {
  void* get_LeptonSFHelper(void);
  void del_LeptonSFHelper(void* ptr);
  float LeptonSFHelper_getSF(void* ptr, int year, int flav, float pt, float eta, float SCeta, bool isCrack);
  float LeptonSFHelper_getSFError(void* ptr, int year, int flav, float pt, float eta, float SCeta, bool isCrack);
}

#endif
