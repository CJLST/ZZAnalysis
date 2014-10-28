#include <ZZAnalysis/AnalysisStep/interface/utils.h>


namespace zzanalysis {
  TLorentzVector tlv(const math::XYZTLorentzVector& v){
    return TLorentzVector(v.x(),v.y(),v.z(),v.t());
  }
}

