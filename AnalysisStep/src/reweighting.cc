#include <assert.h>

#include "TTree.h"

#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>

#include <ZZMatrixElement/MELA/interface/Mela.h>
#include <ZZMatrixElement/MELA/interface/TVar.hh>

enum ReweightingType {
  NoReweighting = 0,
  HVV_spin0     = 1,
  HVV_spin012   = 2
};

class Reweighting {

private:

  //when MELA is updated, can remove these and replace
  //all occurences with mela->selfD???coupl
  double selfDHvvcoupl[SIZE_HVV][2];
  double selfDHwwcoupl[SIZE_HVV][2];
  //double selfDZqqcoupl[SIZE_ZQQ][2];
  double selfDZvvcoupl[SIZE_ZVV][2];
  //double selfDGqqcoupl[SIZE_GQQ][2];
  double selfDGggcoupl[SIZE_GGG][2];
  double selfDGvvcoupl[SIZE_GVV][2];

  double myHvvcoupl[SIZE_HVV][2];
  double myHwwcoupl[SIZE_HVV][2];
  //double myZqqcoupl[SIZE_ZQQ][2];
  double myZvvcoupl[SIZE_ZVV][2];
  //double myGqqcoupl[SIZE_GQQ][2];
  double myGggcoupl[SIZE_GGG][2];
  double myGvvcoupl[SIZE_GVV][2];

  const int g1index = 0;
  const int g2index = 1;
  const int g4index = 3;
  const int g1prime2index = 11;

  const double g2mix = 1.65684;
  const double g4mix = 2.55052;
  const double g1prime2mix = -12100.42;

  const int a1index = 0;
  const int b5index = 4;

  int spin;
  int myspin;
  Mela &myMela;

  void setHVVcouplings(double couplings[SIZE_HVV][2]) {
    for (int i = 0; i < SIZE_HVV; i++)
      for (int j = 0; j < 2; j++)
        selfDHvvcoupl[i][j] = couplings[i][j];
  }
  void setHWWcouplings(double couplings[SIZE_HVV][2]) {
    for (int i = 0; i < SIZE_HVV; i++)
      for (int j = 0; j < 2; j++)
        selfDHwwcoupl[i][j] = couplings[i][j];
  }
  void setZVVcouplings(double couplings[SIZE_ZVV][2]) {
    for (int i = 0; i < SIZE_ZVV; i++)
      for (int j = 0; j < 2; j++)
        selfDZvvcoupl[i][j] = couplings[i][j];
  }
  void setGVVcouplings(double couplings[SIZE_GVV][2]) {
    for (int i = 0; i < SIZE_GVV; i++)
      for (int j = 0; j < 2; j++)
        selfDGvvcoupl[i][j] = couplings[i][j];
  }
  void setGggcouplings(double couplings[SIZE_GGG][2]) {
    for (int i = 0; i < SIZE_GGG; i++)
      for (int j = 0; j < 2; j++)
        selfDGggcoupl[i][j] = couplings[i][j];
  }

  ReweightingType reweightingtypefromstring(std::string reweightingtypestring) {
    if (reweightingtypestring == "none") return NoReweighting;
    else if (reweightingtypestring == "HVV_spin0") return HVV_spin0;
    else if (reweightingtypestring == "HVV_spin012") return HVV_spin012;
    assert(false);
    return NoReweighting;
  }

  template <unsigned int size> void couplingsfromvectors(double (&couplings)[size][2], vector<double> real, vector<double> imaginary) {
    if (real.size() != size || imaginary.size() != size) {
      std::cout << "couplings should have size " << size << " but have size " << real.size() << " " << imaginary.size() << std::endl;
      assert(false);
    }
    for (unsigned int i = 0; i < size; i++) {
      couplings[i][0] = real[i];
      couplings[i][1] = imaginary[i];
    }
  }

public:
  int nReweightingSamples;
  ReweightingType reweightingtype;

  Reweighting(Mela &mela, std::string reweightingtypestring,
              int inputspin,
              std::vector<double> HVVcouplings_real,
              std::vector<double> HVVcouplings_imag,
              std::vector<double> HWWcouplings_real,
              std::vector<double> HWWcouplings_imag,
              std::vector<double> ZVVcouplings_real,
              std::vector<double> ZVVcouplings_imag,
              std::vector<double> Gggcouplings_real,
              std::vector<double> Gggcouplings_imag,
              std::vector<double> GVVcouplings_real,
              std::vector<double> GVVcouplings_imag
             ) : myspin(inputspin),
                 myMela(mela),
                 reweightingtype(reweightingtypefromstring(reweightingtypestring))
  {
    couplingsfromvectors(myHvvcoupl, HVVcouplings_real, HVVcouplings_imag);
    couplingsfromvectors(myHwwcoupl, HWWcouplings_real, HWWcouplings_imag);
    couplingsfromvectors(myZvvcoupl, ZVVcouplings_real, ZVVcouplings_imag);
    couplingsfromvectors(myGggcoupl, Gggcouplings_real, Gggcouplings_imag);
    couplingsfromvectors(myGvvcoupl, GVVcouplings_real, GVVcouplings_imag);
    switch (reweightingtype) {
      case NoReweighting: nReweightingSamples = 0; break;
      case HVV_spin0: nReweightingSamples = 7; break;
      default: assert(false);
    }
  }

  void setmycouplings() {
    spin = myspin;
    switch (myspin) {
      case 0:
        setHVVcouplings(myHvvcoupl);
        setHWWcouplings(myHwwcoupl);
        break;
      case 1:
        setZVVcouplings(myZvvcoupl);
        break;
      case 2:
        setGggcouplings(myGggcoupl);
        setGVVcouplings(myGvvcoupl);
        break;
      default:
        assert(false);
    }
  }

  void setcouplings(int reweightinghypothesis) {
    if (reweightingtype == HVV_spin0) {
      myMela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
      spin = 0;
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
      setHVVcouplings(couplings);
    }
  }

  bool canreweight(unsigned int nleptons, short genFinalState) {
    switch (reweightingtype) {
      case NoReweighting:
        return true;
      case HVV_spin0:
      case HVV_spin012:
        return nleptons == 4 && genFinalState != BUGGY;
      default:
        return false;
    }
  }

  TTree *fillcouplingstree(TTree *t) {
    vector<double> g1, g2, g4, g1prime2, a1, b5;
    vector<int> spin;
    t->Branch("spin", &spin);
    if (reweightingtype == HVV_spin0 || reweightingtype == HVV_spin012) {
      t->Branch("g1", &g1);
      t->Branch("g2", &g2);
      t->Branch("g4", &g4);
      t->Branch("g1prime2", &g1prime2);
    }
    if (reweightingtype == HVV_spin012) {
      t->Branch("a1", &a1);
      t->Branch("b5", &b5);
    }

    for (int i = 0; i < nReweightingSamples; i++) {
      setcouplings(i);
      spin.push_back(myspin);
      if (myspin == 0) {
        g1.push_back(selfDHvvcoupl[g1index][0]);
        g2.push_back(selfDHvvcoupl[g2index][0]);
        g4.push_back(selfDHvvcoupl[g4index][0]);
        g1prime2.push_back(selfDHvvcoupl[g1prime2index][0]);
      }
      else {
        g1.push_back(0);
        g2.push_back(0);
        g4.push_back(0);
        g1prime2.push_back(0);
      }

      if (myspin == 2) {
        a1.push_back(selfDGggcoupl[a1index][0]);
        b5.push_back(selfDGvvcoupl[b5index][0]);
      }
      else {
        a1.push_back(0);
        b5.push_back(0);
      }
    }

    t->Fill();
    return t;
  }

  float computeP(float mZZ, float mZ1, float mZ2,
                 float costhetastar, float costheta1, float costheta2, float phi, float phi1,
                 int flavor
                ) {
    float result;
    switch (spin) {
      case 0:
        myMela.computeP(mZZ, mZ1, mZ2, costhetastar, costheta1, costheta2, phi, phi1, flavor, selfDHvvcoupl, result);
        break;
      case 1:
        myMela.computeP_selfDspin1(mZZ, mZ1, mZ2, costhetastar, costheta1, costheta2, phi, phi1, flavor, selfDZvvcoupl, result);
        break;
      case 2:
        myMela.computeP_selfDspin2(mZZ, mZ1, mZ2, costhetastar, costheta1, costheta2, phi, phi1, flavor, selfDGggcoupl, selfDGvvcoupl, result);
        break;
    }
    return result;
  }
};
