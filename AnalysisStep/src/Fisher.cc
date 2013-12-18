#include <ZZAnalysis/AnalysisStep/interface/Fisher.h>
#include <cmath>

float fisher(float mjj, float detajj) {
  return 0.18*fabs(detajj) + 1.92e-04*mjj;
}

