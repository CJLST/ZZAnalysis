#ifndef BTAGGINGSFHELPER_H
#define BTAGGINGSFHELPER_H

//#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
//#include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include "BTagCalibrationStandalone.h"
#include <string>
#include <vector>
#include <utility>

#include "TFile.h"
#include "TH1F.h"

enum SFsyst {central = 0, up = 1, down = 2};

class BTaggingSFHelper
{

 public:

  BTaggingSFHelper(std::string SFfilename, std::string effFileName);
  ~BTaggingSFHelper();
  float getSF (SFsyst syst, int jetFlavor, float pt, float eta);
  float getEff (int jetFlavor, float pt, float eta);
   
 private:

  // related to scale factors
  BTagCalibration* m_calib;
  std::vector<BTagCalibrationReader*> m_readers; // Loads all of [central, up, down](b, c, udsg)

  // related to b tag efficiency
  TFile* m_fileEff;
  TH1F* m_hEff [3]; // [0: b, 1: c 2: udsg]

};

#endif
