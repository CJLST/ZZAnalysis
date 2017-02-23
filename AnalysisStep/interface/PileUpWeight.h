#ifndef PileUpWeight_h
#define PileUpWeight_h

/** \class PileUpWeight
 *
 *  Implement PU reweighting.
 *
 */

#include <TH1F.h>
#include <TFile.h>
#include <string>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

class PileUpWeight {
public:
  enum class PUvar : int {NOMINAL=0, VARUP=1, VARDOWN=2};

  PileUpWeight(int MC, int target); 

  float weight(float input, PUvar var = PUvar::NOMINAL);

  std::unique_ptr<TH1> h_nominal;
  std::unique_ptr<TH1> h_up;
  std::unique_ptr<TH1> h_down;
};
#endif
