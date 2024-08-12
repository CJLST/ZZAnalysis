#include <correction.h>

class MuonScaleRe {
public:
  MuonScaleRe(std::string json);

  double pt_resol(double pt, double eta, float nL, std::string var);
  
  double pt_scale(bool is_data, double pt, double eta, double phi, int charge, std::string var);

private:
  double get_k(double eta, std::string var);
  double get_std(double pt, double eta, float nL);
  double get_rndm(double eta, float nL);

  std::unique_ptr<correction::CorrectionSet> cset;

};
