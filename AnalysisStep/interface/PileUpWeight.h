#ifndef PileUpWeight_h
#define PileUpWeight_h

/** \class PileUpWeight
 *
 *  Implement PU reweighting.
 *  To prepare PU profile histograms, cf. utils/make_PU_weight_hist.py
 *
 */

#include <TH1F.h>
#include <TFile.h>
#include <string>

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
