#include <ZZAnalysis/AnalysisStep/interface/cConstants.h>

#include <iostream>
#include <cmath>
#include "TMath.h"


extern "C" float getDVBF2jetsConstant(float ZZMass){
  float par[9]={
    1.876,
    -55.488,
    403.32,
    0.3906,
    80.8,
    27.7,
    -0.06,
    54.97,
    309.96
  };
  float kappa =
    pow(1.-atan((ZZMass-par[1])/par[2])*2./TMath::Pi(), par[0])
    + par[3]*exp(-pow((ZZMass-par[4])/par[5], 2))
    + par[6]*exp(-pow((ZZMass-par[7])/par[8], 2));
  float constant = kappa/(1.-kappa);
  return constant;
}
extern "C" float getDVBF1jetConstant(float ZZMass){
  float par[8]={
    0.395,
    -0.07,
    85.,
    30.,
    -0.691,
    -5659.47,
    5734.37,
    0.75
  };
  float kappa =
    par[0]
    + par[1]*exp(-pow((ZZMass-par[2])/par[3], 2))
    + par[4]*pow(log((ZZMass-par[5])/par[6]), par[7])*(ZZMass>=(par[5]+par[6]));
  float constant = kappa/(1.-kappa);
  return constant;
}
extern "C" float getDbkgkinConstant(int ZZflav, float ZZMass){ // ZZflav==id1*id2*id3*id4
  float par[14]={
    0.775,
    -0.565,
    70.,
    5.90,
    -0.235,
    130.1,
    13.25,
    -0.33,
    191.04,
    16.05,
    187.47,
    -0.21,
    1700.,
    400.
  };
  if (abs(ZZflav)==121*121 || abs(ZZflav)==121*242 || abs(ZZflav)==242*242) par[11]=-0.42; // 4e
  float kappa =
    par[0]
    +par[1]*exp(-pow(((ZZMass-par[2])/par[3]), 2))
    +par[4]*exp(-pow(((ZZMass-par[5])/par[6]), 2))
    +par[7]*(
	     exp(-pow(((ZZMass-par[8])/par[9]), 2))*(ZZMass<par[8])
	     + exp(-pow(((ZZMass-par[8])/par[10]), 2))*(ZZMass>=par[8])
	     )
    + par[11]*exp(-pow(((ZZMass-par[12])/par[13]), 2));

  float constant = kappa/(1.-kappa);
  return constant;
}
extern "C" float getDbkgConstant(int ZZflav, float ZZMass){
  return getDbkgkinConstant(ZZflav, ZZMass);
}
