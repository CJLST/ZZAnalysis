#ifndef HZZ4l_h
#define HZZ4l_h

#include "HZZ4lBase.h"

#include <string>
#include <vector>

class HZZ4l : public HZZ4lBase{
public:
  
  HZZ4l(TChain *tree=0);
  virtual ~HZZ4l() {};
  void Loop(const Int_t Nevt_gen, const Int_t channelType, const std::string outputName, const bool is8TeV, const Int_t whichVH);

private:

  Float_t getMCWeight(const int year, const Int_t CandIndex) const;
  Float_t getAllWeight(const Float_t LepPt, const Float_t LepEta, const Int_t year, Int_t LepID) const;

  Float_t getNormalizedWeight(const std::string theSample, const Int_t channelType, const Bool_t is8TeV, const Int_t whichVH) const;

  std::string getSampleName(const std::string filestring) const;
  Int_t findBestCRCand() const;

  void getWeightFromFile(const std::string mPOLE, const Bool_t is8TeV, const Bool_t isVBF);
  Float_t getPwhgWeight(double mH, const Int_t WhichSide = 0) const;

  std::vector<double> pwhg_bincenters;
  std::vector<double> pwhg_weight;
  std::vector<double> pwhg_weightP;
  std::vector<double> pwhg_weightM;

};
#endif
