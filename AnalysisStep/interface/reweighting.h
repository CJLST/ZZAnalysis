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
  HVV_spin012   = 2
};

class Reweighting {

protected:

  double myHvvcoupl[SIZE_HVV][2];
  double myZvvcoupl[SIZE_ZVV][2];
  double myGggcoupl[SIZE_GGG][2];
  double myGvvcoupl[SIZE_GVV][2];

  const int ghz1_index;
  const int ghz2_index;
  const int ghz4_index;
  const int ghz1_prime2_index;

  const double ghz2mix;
  const double ghz4mix;
  const double ghz1_prime2mix;

  const int a1_index;
  const int b5_index;

  int spin;
  int myspin;
  Mela& mela;

  ReweightingType reweightingtypefromstring(const std::string& reweightingtypestring);

  template <unsigned int size> void couplingsfromvectors(double(&couplings)[size][2], vector<double> real, vector<double> imaginary) { // Could just have had a pair -- U. Sarica
    if (real.size() != size || imaginary.size() != size) {
      std::cout << "couplings should have size " << size << " but have size " << real.size() << " " << imaginary.size() << std::endl;
      assert(false);
    }
    for (unsigned int i = 0; i < size; i++) {
      couplings[i][0] = real[i];
      couplings[i][1] = imaginary[i];
    }
  }

  void fillreweightingweights(vector<float> &reweightingweights);

  bool canreweight();

public:

  int nReweightingSamples;
  ReweightingType reweightingtype;

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
    std::vector<double> GVVcouplings_imag
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

