#ifndef REWEIGHTING_CC
#define REWEIGHTING_CC


#include <ZZAnalysis/AnalysisStep/interface/reweighting.h>


ReweightingType Reweighting::reweightingtypefromstring(const std::string& reweightingtypestring) {
  if (reweightingtypestring == "none") return NoReweighting;
  else if (reweightingtypestring == "HVV_spin0") return HVV_spin0;
  else if (reweightingtypestring == "HVV_spin012") return HVV_spin012;
  assert(false);
  return NoReweighting;
}

Reweighting::Reweighting(
  Mela* mela_,
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
  ) :
  ghz1_index(0),
  ghz2_index(1),
  ghz4_index(3),
  ghz1_prime2_index(11),
  ghz2mix(1.65684),
  ghz4mix(2.55052),
  ghz1_prime2mix(-12100.42),
  a1_index(0),
  b5_index(4),
  myspin(inputspin),
  mela(mela_),
  reweightingtype(reweightingtypefromstring(reweightingtypestring))
{
  couplingsfromvectors(myHvvcoupl, HVVcouplings_real, HVVcouplings_imag);
  couplingsfromvectors(myZvvcoupl, ZVVcouplings_real, ZVVcouplings_imag);
  couplingsfromvectors(myGggcoupl, Gggcouplings_real, Gggcouplings_imag);
  couplingsfromvectors(myGvvcoupl, GVVcouplings_real, GVVcouplings_imag);
  switch (reweightingtype) {
  case NoReweighting: nReweightingSamples = 0; break;
  case HVV_spin0: nReweightingSamples = 7; break;
  case HVV_spin012: nReweightingSamples = 8; break;
  default: assert(false);
  }
}

void Reweighting::setmycouplings() {
  spin = myspin;
  switch (myspin) {
  case 0: mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG); break;
  case 1: mela->setProcess(TVar::SelfDefine_spin1, TVar::JHUGen, TVar::ZZGG); break;
  case 2: mela->setProcess(TVar::SelfDefine_spin2, TVar::JHUGen, TVar::ZZGG); break;
  default: assert(false);
  }

  mela->selfDHggcoupl[0][0] = 1;
  for (unsigned int ic=0; ic<SIZE_HVV; ic++){ for (unsigned int im=0; im<2; im++) mela->selfDHzzcoupl[0][ic][im] = myHvvcoupl[ic][im]; }
  for (unsigned int ic=0; ic<SIZE_GGG; ic++){ for (unsigned int im=0; im<2; im++) mela->selfDGggcoupl[ic][im] = myGggcoupl[ic][im]; }
  for (unsigned int ic=0; ic<SIZE_GVV; ic++){ for (unsigned int im=0; im<2; im++) mela->selfDGvvcoupl[ic][im] = myGvvcoupl[ic][im]; }
}

void Reweighting::setcouplings(int reweightinghypothesis) {
  if (reweightingtype == NoReweighting) return;
  if (reweightingtype == HVV_spin0) {
    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    spin = 0;
    mela->selfDHggcoupl[0][0] = 1;
    switch (reweightinghypothesis) {
    case 0: mela->selfDHzzcoupl[0][ghz1_index][0] = 1; break;                                            //0+m
    case 1: mela->selfDHzzcoupl[0][ghz2_index][0] = 1; break;                                            //0+h
    case 2: mela->selfDHzzcoupl[0][ghz4_index][0] = 1; break;                                            //0-
    case 3: mela->selfDHzzcoupl[0][ghz1_prime2_index][0] = 1; break;                                      //L1
    case 4: mela->selfDHzzcoupl[0][ghz1_index][0] = 1; mela->selfDHzzcoupl[0][ghz2_index][0] = ghz2mix; break;             //fa2=0.5
    case 5: mela->selfDHzzcoupl[0][ghz1_index][0] = 1; mela->selfDHzzcoupl[0][ghz4_index][0] = ghz4mix; break;             //fa3=0.5
    case 6: mela->selfDHzzcoupl[0][ghz1_index][0] = 1; mela->selfDHzzcoupl[0][ghz1_prime2_index][0] = ghz1_prime2mix; break; //fL1=0.5
    default: assert(false);
    }
    return;
  }
  else if (reweightingtype == HVV_spin012) {
    mela->selfDHggcoupl[0][0] = 1;
    switch (reweightinghypothesis) {
    case 0: mela->selfDHzzcoupl[0][ghz1_index][0] = 1; spin = 0; break;                                               //0+m
    case 1: mela->selfDHzzcoupl[0][ghz2_index][0] = 1; spin = 0; break;                                               //0+h
    case 2: mela->selfDHzzcoupl[0][ghz4_index][0] = 1; spin = 0; break;                                               //0-
    case 3: mela->selfDHzzcoupl[0][ghz1_prime2_index][0] = 1; spin = 0; break;                                         //L1
    case 4: mela->selfDHzzcoupl[0][ghz1_index][0] = 1; mela->selfDHzzcoupl[0][ghz2_index][0] = ghz2mix; spin = 0; break;             //fa2=0.5
    case 5: mela->selfDHzzcoupl[0][ghz1_index][0] = 1; mela->selfDHzzcoupl[0][ghz4_index][0] = ghz4mix; spin = 0; break;             //fa3=0.5
    case 6: mela->selfDHzzcoupl[0][ghz1_index][0] = 1; mela->selfDHzzcoupl[0][ghz1_prime2_index][0] = ghz1_prime2mix; spin = 0; break; //fL1=0.5
    case 7: mela->selfDGggcoupl[a1_index][0] = 1; mela->selfDGvvcoupl[b5_index][0] = 1; spin = 2; break;                 //2b+
    default: assert(false);
    }
    switch (spin) {
    case 0: mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG); break;
    case 1: mela->setProcess(TVar::SelfDefine_spin1, TVar::JHUGen, TVar::ZZGG); break;
    case 2: mela->setProcess(TVar::SelfDefine_spin2, TVar::JHUGen, TVar::ZZGG); break;
    default: cout << spin << endl; assert(false);
    }
    return;
  }
  assert(false);
}

bool Reweighting::canreweight(unsigned int nleptons, short genFinalState) {
  switch (reweightingtype) {
  case NoReweighting:
    return true;
  case HVV_spin0:
  case HVV_spin012:
    return genFinalState != BUGGY;
  default:
    return false;
  }
}

void Reweighting::fillcouplingstree(TTree* t) {
  vector<double> ghz1, ghz2, ghz4, ghz1_prime2, a1, b5;
  vector<int> spin_v;
  t->Branch("spin", &spin_v);
  if (reweightingtype == HVV_spin0 || reweightingtype == HVV_spin012) {
    t->Branch("ghz1Re", &ghz1);
    t->Branch("ghz2Re", &ghz2);
    t->Branch("ghz4Re", &ghz4);
    t->Branch("ghz1_prime2Re", &ghz1_prime2);
  }
  if (reweightingtype == HVV_spin012) {
    t->Branch("a1Re", &a1);
    t->Branch("b5Re", &b5);
  }

  for (int i = 0; i < nReweightingSamples; i++) {
    setcouplings(i);
    spin_v.push_back(spin);
    if (spin == 0) {
      ghz1.push_back(mela->selfDHzzcoupl[0][ghz1_index][0]);
      ghz2.push_back(mela->selfDHzzcoupl[0][ghz2_index][0]);
      ghz4.push_back(mela->selfDHzzcoupl[0][ghz4_index][0]);
      ghz1_prime2.push_back(mela->selfDHzzcoupl[0][ghz1_prime2_index][0]);
    }
    else {
      ghz1.push_back(0);
      ghz2.push_back(0);
      ghz4.push_back(0);
      ghz1_prime2.push_back(0);
    }

    if (spin == 2) {
      a1.push_back(mela->selfDGggcoupl[a1_index][0]);
      b5.push_back(mela->selfDGvvcoupl[b5_index][0]);
    }
    else {
      a1.push_back(0);
      b5.push_back(0);
    }
  }
  t->Fill();
}

float Reweighting::computeP(
  float mzz, float m1, float m2,
  float hs, float h1, float h2, float phi, float phi1,
  int flavor // Need to replace the flavor and actuially this entire stuff, extremely inefficient
  ){
  float result;

  int idOrdered[4];
  if (flavor == 2){
    idOrdered[0]=13;
    idOrdered[1]=-13;
    idOrdered[2]=11;
    idOrdered[3]=-11;
  }
  else if (flavor == 1){
    idOrdered[0]=11;
    idOrdered[1]=-11;
    idOrdered[2]=11;
    idOrdered[3]=-11;
  }
  else if (flavor == 0){
    idOrdered[0]=13;
    idOrdered[1]=-13;
    idOrdered[2]=13;
    idOrdered[3]=-13;
  }
  else return 0;
  TLorentzVector pOrdered[4];
  std::vector<TLorentzVector> daus = mela->calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
  for (int ip=0; ip<min(4, (int)daus.size()); ip++) pOrdered[ip]=daus.at(ip);
  SimpleParticleCollection_t daughters;
  for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(idOrdered[idau], pOrdered[idau]));

  mela->setInputEvent(&daughters, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);
  mela->computeP(result, false);
  mela->resetInputEvent(); // Poor efficiency, result of passing vectors for each ME

  return result;
}


#endif

