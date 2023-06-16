#ifndef utils_h
#define utils_h

#include <TLorentzVector.h>
#include <DataFormats/Math/interface/LorentzVector.h>

namespace zzanalysis{
  TLorentzVector tlv(const math::XYZTLorentzVector& v);
}
double SetupToSqrts(int setup);



#endif

