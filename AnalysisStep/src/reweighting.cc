#include <assert.h>

#include "TTree.h"

#include <ZZMatrixElement/MELA/interface/TVar.hh>

const int nReweightingSamples = 7;
//when MELA is updated, can remove these and replace
//all occurences with mela->selfDHvvcoupl
double selfDHvvcoupl[SIZE_HVV][2];
double selfDHwwcoupl[SIZE_HVV][2];

namespace {
  const int g1index = 0;
  const int g2index = 1;
  const int g4index = 3;
  const int g1prime2index = 11;

  const double g2mix = 1.65684;
  const double g4mix = 2.55052;
  const double g1prime2mix = -12100.42;
}

void setcouplings(double couplings[SIZE_HVV][2]) {
  for (int i = 0; i < SIZE_HVV; i++)
    for (int j = 0; j < 2; j++)
      selfDHvvcoupl[i][j] = couplings[i][j];
}

void setWWcouplings(double couplings[SIZE_HVV][2]) {
  for (int i = 0; i < SIZE_HVV; i++)
    for (int j = 0; j < 2; j++)
      selfDHwwcoupl[i][j] = couplings[i][j];
}

void setcouplings(int reweightinghypothesis) {
  double couplings[SIZE_HVV][2] = {{0}};
  switch (reweightinghypothesis) {
    case 0: couplings[g1index][0] = 1; break;                                            //0+m
    case 1: couplings[g2index][0] = 1; break;                                            //0+h
    case 2: couplings[g4index][0] = 1; break;                                            //0-
    case 3: couplings[g1prime2index][0] = 1; break;                                      //L1
    case 4: couplings[g1index][0] = 1; couplings[g2index][0] = g2mix; break;             //fa2=0.5
    case 5: couplings[g1index][0] = 1; couplings[g4index][0] = g4mix; break;             //fa3=0.5
    case 6: couplings[g1index][0] = 1; couplings[g1prime2index][0] = g1prime2mix; break; //fL1=0.5
    default: assert(false);
  }
  setcouplings(couplings);
}

TTree *fillcouplingstree(TTree *t) {
  double g1, g2, g4, g1prime2;
  t->Branch("g1", &g1, "g1/D");
  t->Branch("g2", &g2, "g2/D");
  t->Branch("g4", &g4, "g4/D");
  t->Branch("g1prime2", &g1prime2, "g1prime2/D");
  for (int i = 0; i < nReweightingSamples; i++) {
    setcouplings(i);
    g1 = selfDHvvcoupl[g1index][0];
    g2 = selfDHvvcoupl[g2index][0];
    g4 = selfDHvvcoupl[g4index][0];
    g1prime2 = selfDHvvcoupl[g1prime2index][0];
    t->Fill();
  }
  return t;
}
