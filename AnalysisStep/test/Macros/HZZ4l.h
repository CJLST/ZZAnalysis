#ifndef HZZ4l_h
#define HZZ4l_h

#include "HZZ4lBase.h"

#include <TString.h>
#include <string>
#include <vector>

class TH1F;

class HZZ4l : public HZZ4lBase{
public:
  
  HZZ4l(TChain *tree, TString sampleName);
  virtual ~HZZ4l() {};
  void Loop(const Int_t channelType, const TString outputName);

  //Methods to set the analyzer
  void set8TeV(bool myis8TeV){is8TeV = myis8TeV;}
  void setBSM(int myBSM_flag){BSM_flag = myBSM_flag;}
  void setKappa(float mykappa){kappa = mykappa;}
  void setCR(bool myisCR){isCR = myisCR;}

private:

  Float_t getMCWeight(const int year, const Int_t CandIndex) const;
  Float_t getAllWeight(const Float_t LepPt, const Float_t LepEta, const Int_t year, Int_t LepID) const;

  Float_t getNormalizedWeight(const Int_t channelType) const;

  std::string getSampleName(const std::string filestring) const;
  Int_t findBestCRCand(int) const;

  void getWeightFromFile(const TString& mPOLE, const Bool_t isVBF, const Bool_t isNewHighmass);
  Float_t getPwhgWeight(double mH, const Int_t WhichSide = 0) const;
  Float_t getHqTWeight(double mH, double genPt, TFile* f) const;

  Float_t getFakeWeight(const Float_t LepPt, const Float_t LepEta, const Int_t year, Int_t LepID, Int_t LepZ1ID);
  Float_t getZXfake_weight(const int year, const Int_t CandIndex);

  TString theSample;

  bool is8TeV;
  bool isCR;
  bool isSignal;
  int BSM_flag;
  float kappa;
  
  TH1F* nEventComplete;

  std::vector<double> pwhg_bincenters;
  std::vector<double> pwhg_weight;
  std::vector<double> pwhg_weightCPSP;
  std::vector<double> pwhg_weightCPSM;
  std::vector<double> pwhg_weightIntP;
  std::vector<double> pwhg_weightIntM;

  TFile* ZXWeightTables[4];

};
#endif
