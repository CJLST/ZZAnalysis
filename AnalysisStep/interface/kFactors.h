#ifndef KFACTORS_H
#define KFACTORS_H
// We wrap kFactor functions in a struct so that it is possible to generate a
// dictionary for it

struct KFactors {
  static float kfactor_qqZZ_qcd_dPhi(float GENabsdPhiZZ, int finalState);  
  static float xsec_qqZZ_qcd_M(float GenMassZZ, int finalState, int order);
  static float kfactor_qqZZ_qcd_M(float GENmassZZ, int finalState, int order);
  static float kfactor_qqZZ_qcd_Pt(float GENpTZZ, int finalState);
};
#endif

