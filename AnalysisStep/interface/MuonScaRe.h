/*  
 * Run3 Muon correction module adapted after https://gitlab.cern.ch/cms-muonPOG/muonscarekit/ .
 * See .cc file for techincal notes.
 */

#include <correction.h>
#include <TRandom3.h>
#include <string>

class MuonScaRe {
public:
  MuonScaRe(std::string json);

  double pt_resol(double pt, double eta, float nL, std::string var);
  
  double pt_scale(bool is_data, double pt, double eta, double phi, int charge, std::string var);

  // A per-muon seed can be optionally set to achieve deterministic random
  // smearing. 
  void setSeed(ULong_t seed){
    rng.SetSeed(seed);
  }

private:
  double get_k(double eta, std::string var);
  double get_std(double pt, double eta, float nL);
  double get_rndm(double eta, float nL);

  std::unique_ptr<correction::CorrectionSet> cset;
  TRandom3 rng;
};
