#ifndef REWEIGHTING_H
#define REWEIGHTING_H


#include <assert.h>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>
#include <ZZMatrixElement/MELA/interface/Mela.h>
#include <ZZMatrixElement/MELA/interface/TVar.hh>
#include "TTree.h"

enum ReweightingType {
  NoReweighting = 0,
  HVV_spin0     = 1,
  HVV_spin012   = 2,
  VBFHZZ_spin0  = 3,
};

class Reweighting {

protected:

  double myHvvcoupl[SIZE_HVV][2];
  double myZvvcoupl[SIZE_ZVV][2];
  double myGggcoupl[SIZE_GGG][2];
  double myGvvcoupl[SIZE_GVV][2];

  static const int ghz1_index = 0;
  static const int ghz2_index = 1;
  static const int ghz4_index = 3;
  static const int ghz1_prime2_index = 11;

  static const constexpr double ghz2mix_decay = 1.663195;
  static const constexpr double ghz4mix_decay = 2.55502;
  static const constexpr double ghz1_prime2mix_decay = -12110.20;
  static const constexpr double ghz2mix_VBF = 0.271965;
  static const constexpr double ghz4mix_VBF = 0.297979;
  static const constexpr double ghz1_prime2mix_VBF = -2158.21;

  static const int a1_index = 0;
  static const int b5_index = 4;

  int spin;
  int myspin;
  typedef tuple<TVar::Process, TVar::MatrixElement, TVar::Production> MelaProcess;
  MelaProcess decayprocess;
  MelaProcess productionprocess;
  vector<double> cutoffs;

  Mela& mela;

  double& ghz1;
  double& ghz2;
  double& ghz4;
  double& ghz1_prime2;
  double& a1;
  double& b5;

  ReweightingType reweightingtypefromstring(const std::string& reweightingtypestring);

  template <unsigned int size> void couplingsfromvectors(double(&couplings)[size][2], vector<double> real, vector<double> imaginary) { // Could just have had a pair -- U. Sarica
    if (real.size() != size || imaginary.size() != size) {
      std::cout << "couplings should have size " << size << " but has size " << real.size() << " " << imaginary.size() << std::endl;
      assert(false);
    }
    for (unsigned int i = 0; i < size; i++) {
      couplings[i][0] = real[i];
      couplings[i][1] = imaginary[i];
    }
  }

  void fillreweightingweights(vector<float> &reweightingweights);

  bool canreweight();
  void setProcess(MelaProcess melaprocess) {
    mela.setProcess(get<0>(melaprocess), get<1>(melaprocess), get<2>(melaprocess));
  }
  void computeP(int reweightinghypothesis, float& prob, bool useConstant = true);

  void reset_SelfDCouplings() {
    mela.setInputEvent(0, 0, 0, true);
    float dummy;
    mela.computeP(dummy, false);
  }

  int nReweightingSamplesFromType(ReweightingType reweightingtype_);

public:

  ReweightingType reweightingtype;
  int nReweightingSamples;

  Reweighting(
    Mela& mela_,
    std::string reweightingtypestring,
    int inputspin,
    std::vector<double> HVVcouplings_real,
    std::vector<double> HVVcouplings_imag,
    std::vector<double> ZVVcouplings_real,
    std::vector<double> ZVVcouplings_imag,
    std::vector<double> Gggcouplings_real,
    std::vector<double> Gggcouplings_imag,
    std::vector<double> GVVcouplings_real,
    std::vector<double> GVVcouplings_imag,
    std::vector<double> cutoffs
    );
  virtual ~Reweighting(){}

  void setmycouplings();

  void setcouplings(int reweightinghypothesis);

  void fillcouplingstree(TTree* t);

  void fillreweightingweights(
    vector<float> &reweightingweights,
    SimpleParticleCollection_t* pDaughters,
    SimpleParticleCollection_t* pAssociated,
    SimpleParticleCollection_t* pMothers,
    bool isGen
  );

  void fillreweightingweights(
    vector<float> &reweightingweights,
    MELACandidate* candidate
  );
};


#endif
