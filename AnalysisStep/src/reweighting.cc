#ifndef REWEIGHTING_CC
#define REWEIGHTING_CC


#include <ZZAnalysis/AnalysisStep/interface/reweighting.h>


ReweightingType Reweighting::reweightingtypefromstring(const std::string& reweightingtypestring) {
  if (reweightingtypestring == "none") return NoReweighting;
  else if (reweightingtypestring == "HVV_spin0") return HVV_spin0;
  else if (reweightingtypestring == "HVV_spin012") return HVV_spin012;
  else if (reweightingtypestring == "VBFHZZ_spin0") return VBFHZZ_spin0;
  assert(false);
  return NoReweighting;
}

int Reweighting::nReweightingSamplesFromType(ReweightingType reweightingtype_) {
  switch (reweightingtype_) {
    case NoReweighting: return 0;
    case HVV_spin0: return 7;
    case HVV_spin012: return 8;
    case VBFHZZ_spin0: return 13;
    default: assert(false);
  }
}

Reweighting::Reweighting(
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
  std::vector<double> cutoffs_
  ) :
  myspin(inputspin),
  cutoffs(cutoffs_),
  mela(mela_),
  reweightingtype(reweightingtypefromstring(reweightingtypestring)),
  nReweightingSamples(nReweightingSamplesFromType(reweightingtype))
{
  if (cutoffs.size() == 0)
    for (int i = 0; i < nReweightingSamples; i++)
      cutoffs.push_back(-1);
  if (cutoffs.size() != (unsigned)nReweightingSamples) {
    cout << "Wrong number of cutoffs (" << cutoffs.size() << ", should be " << nReweightingSamples << ")!" << endl;
    assert(false);
  }
  couplingsfromvectors(myHvvcoupl, HVVcouplings_real, HVVcouplings_imag);
  couplingsfromvectors(myZvvcoupl, ZVVcouplings_real, ZVVcouplings_imag);
  couplingsfromvectors(myGggcoupl, Gggcouplings_real, Gggcouplings_imag);
  couplingsfromvectors(myGvvcoupl, GVVcouplings_real, GVVcouplings_imag);
}

void Reweighting::setmycouplings() {
  if (reweightingtype == NoReweighting) {
    return;
  } else if (reweightingtype == HVV_spin0 || reweightingtype == HVV_spin012) {
    productionprocess = make_tuple(TVar::Null, TVar::JHUGen, TVar::ZZGG);
    switch (myspin) {
    case 0: decayprocess = make_tuple(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG); break;
    case 1: decayprocess = make_tuple(TVar::SelfDefine_spin1, TVar::JHUGen, TVar::ZZGG); break;
    case 2: decayprocess = make_tuple(TVar::SelfDefine_spin2, TVar::JHUGen, TVar::ZZGG); break;
    default: assert(false);
    }
  } else if (reweightingtype == VBFHZZ_spin0) {
    assert(myspin == 0);
    productionprocess = make_tuple(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJVBF);
    decayprocess = make_tuple(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZINDEPENDENT);
  } else {
    assert(false);
  }

  spin = myspin;

  mela.selfDHggcoupl[0][0] = 1;
  for (unsigned int ic=0; ic<SIZE_HVV; ic++){ for (unsigned int im=0; im<2; im++) mela.selfDHzzcoupl[0][ic][im] = myHvvcoupl[ic][im]; }
  for (unsigned int ic=0; ic<SIZE_GGG; ic++){ for (unsigned int im=0; im<2; im++) mela.selfDGggcoupl[ic][im] = myGggcoupl[ic][im]; }
  for (unsigned int ic=0; ic<SIZE_GVV; ic++){ for (unsigned int im=0; im<2; im++) mela.selfDGvvcoupl[ic][im] = myGvvcoupl[ic][im]; }
}

void Reweighting::setcouplings(int reweightinghypothesis) {
  if (reweightingtype == NoReweighting) {
    return;
  } else if (reweightinghypothesis == -1) {
    return setmycouplings();
  } else if (reweightingtype == HVV_spin0) {
    decayprocess = make_tuple(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    productionprocess = make_tuple(TVar::Null, TVar::JHUGen, TVar::ZZGG);
    spin = 0;
    mela.selfDHggcoupl[0][0] = 1;
    switch (reweightinghypothesis) {
    case 0: mela.selfDHzzcoupl[0][ghz1_index][0] = 1; break;                                            //0+m
    case 1: mela.selfDHzzcoupl[0][ghz2_index][0] = 1; break;                                            //0+h
    case 2: mela.selfDHzzcoupl[0][ghz4_index][0] = 1; break;                                            //0-
    case 3: mela.selfDHzzcoupl[0][ghz1_prime2_index][0] = 1; break;                                      //L1
    case 4: mela.selfDHzzcoupl[0][ghz1_index][0] = 1; mela.selfDHzzcoupl[0][ghz2_index][0] = ghz2mix_decay; break;             //fa2=0.5
    case 5: mela.selfDHzzcoupl[0][ghz1_index][0] = 1; mela.selfDHzzcoupl[0][ghz4_index][0] = ghz4mix_decay; break;             //fa3=0.5
    case 6: mela.selfDHzzcoupl[0][ghz1_index][0] = 1; mela.selfDHzzcoupl[0][ghz1_prime2_index][0] = ghz1_prime2mix_decay; break; //fL1=0.5
    default: assert(false);
    }
    return;
  } else if (reweightingtype == HVV_spin012) {
    productionprocess = make_tuple(TVar::Null, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][0] = 1;
    switch (reweightinghypothesis) {
    case 0: mela.selfDHzzcoupl[0][ghz1_index][0] = 1; spin = 0; break;                                               //0+m
    case 1: mela.selfDHzzcoupl[0][ghz2_index][0] = 1; spin = 0; break;                                               //0+h
    case 2: mela.selfDHzzcoupl[0][ghz4_index][0] = 1; spin = 0; break;                                               //0-
    case 3: mela.selfDHzzcoupl[0][ghz1_prime2_index][0] = 1; spin = 0; break;                                         //L1
    case 4: mela.selfDHzzcoupl[0][ghz1_index][0] = 1; mela.selfDHzzcoupl[0][ghz2_index][0] = ghz2mix_decay; spin = 0; break;             //fa2=0.5
    case 5: mela.selfDHzzcoupl[0][ghz1_index][0] = 1; mela.selfDHzzcoupl[0][ghz4_index][0] = ghz4mix_decay; spin = 0; break;             //fa3=0.5
    case 6: mela.selfDHzzcoupl[0][ghz1_index][0] = 1; mela.selfDHzzcoupl[0][ghz1_prime2_index][0] = ghz1_prime2mix_decay; spin = 0; break; //fL1=0.5
    case 7: mela.selfDGggcoupl[a1_index][0] = 1; mela.selfDGvvcoupl[b5_index][0] = 1; spin = 2; break;                 //2b+
    default: assert(false);
    }
    switch (spin) {
    case 0: decayprocess = make_tuple(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG); break;
    case 1: decayprocess = make_tuple(TVar::SelfDefine_spin1, TVar::JHUGen, TVar::ZZGG); break;
    case 2: decayprocess = make_tuple(TVar::SelfDefine_spin2, TVar::JHUGen, TVar::ZZGG); break;
    default: cout << spin << endl; assert(false);
    }
    return;
  } else if (reweightingtype == VBFHZZ_spin0) {
    decayprocess = make_tuple(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZINDEPENDENT);
    productionprocess = make_tuple(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJVBF);
    spin = 0;
    mela.selfDHggcoupl[0][0] = 1;
    switch (reweightinghypothesis) {
    case 0: mela.selfDHzzcoupl[0][ghz1_index][0] = 1; break;                                            //0+m
    case 1: mela.selfDHzzcoupl[0][ghz2_index][0] = 1; break;                                            //0+h
    case 2: mela.selfDHzzcoupl[0][ghz4_index][0] = 1; break;                                            //0-
    case 3: mela.selfDHzzcoupl[0][ghz1_prime2_index][0] = 1; break;                                     //L1
    case 4: mela.selfDHzzcoupl[0][ghz1_index][0] = 1; mela.selfDHzzcoupl[0][ghz2_index][0] = ghz2mix_decay; break;               //fa2_decay=0.5
    case 5: mela.selfDHzzcoupl[0][ghz1_index][0] = 1; mela.selfDHzzcoupl[0][ghz4_index][0] = ghz4mix_decay; break;               //fa3_decay=0.5
    case 6: mela.selfDHzzcoupl[0][ghz1_index][0] = 1; mela.selfDHzzcoupl[0][ghz1_prime2_index][0] = ghz1_prime2mix_decay; break; //fL1_decay=0.5
    case 7: mela.selfDHzzcoupl[0][ghz1_index][0] = 1; mela.selfDHzzcoupl[0][ghz2_index][0] = ghz2mix_VBF; break;                 //fa2_VBF=0.5
    case 8: mela.selfDHzzcoupl[0][ghz1_index][0] = 1; mela.selfDHzzcoupl[0][ghz4_index][0] = ghz4mix_VBF; break;                 //fa3_VBF=0.5
    case 9: mela.selfDHzzcoupl[0][ghz1_index][0] = 1; mela.selfDHzzcoupl[0][ghz1_prime2_index][0] = ghz1_prime2mix_VBF; break;   //fL1_VBF=0.5
    case 10: mela.selfDHzzcoupl[0][ghz1_index][0] = 1; mela.selfDHzzcoupl[0][ghz2_index][0] = -sqrt(ghz2mix_VBF*ghz2mix_decay); break; //fa2_VBFdecay=-0.5
    case 11: mela.selfDHzzcoupl[0][ghz1_index][0] = 1; mela.selfDHzzcoupl[0][ghz4_index][0] = -sqrt(ghz4mix_VBF*ghz4mix_decay); break; //fa3_VBFdecay=-0.5
    case 12: mela.selfDHzzcoupl[0][ghz1_index][0] = 1; mela.selfDHzzcoupl[0][ghz1_prime2_index][0] = -sqrt(ghz1_prime2mix_VBF*ghz1_prime2mix_decay); break; //fL1_VBFdecay=-0.5
    default: assert(false);
    }
    return;
  }
  assert(false);
}

bool Reweighting::canreweight() {
  //probably can just return true or delete this function
  switch (reweightingtype) {
  case NoReweighting:
    return true;
  case HVV_spin0:
  case HVV_spin012:
  case VBFHZZ_spin0:
    return true; //gives something sensible if not enough daughters
  default:
    return false;
  }
}

void Reweighting::fillcouplingstree(TTree* t) {
  vector<double> ghz1, ghz2, ghz4, ghz1_prime2, a1, b5;
  vector<int> spin_v;
  t->Branch("spin", &spin_v);
  if (reweightingtype == HVV_spin0 || reweightingtype == HVV_spin012 || reweightingtype == VBFHZZ_spin0) {
    t->Branch("ghz1Re", &ghz1);
    t->Branch("ghz2Re", &ghz2);
    t->Branch("ghz4Re", &ghz4);
    t->Branch("ghz1_prime2Re", &ghz1_prime2);
  }
  if (reweightingtype == HVV_spin012) {
    t->Branch("a1Re", &a1);
    t->Branch("b5Re", &b5);
  }

  cerr << "There are about to be " << nReweightingSamples << " warnings about no daughters.  Please ignore them." << endl;
  for (int i = 0; i < nReweightingSamples; i++) {
    reset_SelfDCouplings();
    setcouplings(i);
    spin_v.push_back(spin);
    if (spin == 0) {
      ghz1.push_back(mela.selfDHzzcoupl[0][ghz1_index][0]);
      ghz2.push_back(mela.selfDHzzcoupl[0][ghz2_index][0]);
      ghz4.push_back(mela.selfDHzzcoupl[0][ghz4_index][0]);
      ghz1_prime2.push_back(mela.selfDHzzcoupl[0][ghz1_prime2_index][0]);
    }
    else {
      ghz1.push_back(0);
      ghz2.push_back(0);
      ghz4.push_back(0);
      ghz1_prime2.push_back(0);
    }

    if (spin == 2) {
      a1.push_back(mela.selfDGggcoupl[a1_index][0]);
      b5.push_back(mela.selfDGvvcoupl[b5_index][0]);
    }
    else {
      a1.push_back(0);
      b5.push_back(0);
    }
  }
  t->Fill();
}

void Reweighting::fillreweightingweights(
  vector<float>& reweightingweights,
  SimpleParticleCollection_t* pDaughters,
  SimpleParticleCollection_t* pAssociated,
  SimpleParticleCollection_t* pMothers,
  bool isGen
) {
  if (nReweightingSamples == 0) return;
  mela.setInputEvent(pDaughters, pAssociated, pMothers, isGen);
  fillreweightingweights(reweightingweights);
  mela.resetInputEvent();
}


void Reweighting::fillreweightingweights(
  vector<float>& reweightingweights,
  MELACandidate* candidate
) {
  if (nReweightingSamples == 0) return;

  //mela.setCurrentCandidate(candidate);  //doesn't work, does not append the candidate to TEvtProb::candList and so it doesn't know that it exists
  /////////////////////////////////////////////////////////////
  //Hopefully Ulascan will tell me the easier way to do this...
  SimpleParticleCollection_t pDaughters, pAssociated, pMothers;
  for (int i = 0; i < candidate->getNDaughters(); i++) {
    auto daughter = candidate->getSortedDaughter(i);
    pDaughters.emplace_back(daughter->id, daughter->p4);
  }
  for (int i = 0; i < candidate->getNAssociatedLeptons(); i++) {
    auto associated = candidate->getAssociatedLepton(i);
    pAssociated.emplace_back(associated->id, associated->p4);
  }
  //skip neutrinos, https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement/blob/71f2f5458464e45ad85a9386be2c1247332476da/MELA/src/MELACandidate.cc#L448
  for (int i = 0; i < candidate->getNAssociatedPhotons(); i++) {
    auto associated = candidate->getAssociatedPhoton(i);
    pAssociated.emplace_back(associated->id, associated->p4);
  }
  for (int i = 0; i < candidate->getNAssociatedJets(); i++) {
    auto associated = candidate->getAssociatedJet(i);
    pAssociated.emplace_back(associated->id, associated->p4);
  }
  for (int i = 0; i < candidate->getNMothers(); i++) {
    auto mother = candidate->getMother(i);
    pMothers.emplace_back(mother->id, mother->p4);
  }
  mela.setInputEvent(&pDaughters, &pAssociated, &pMothers, true);
  /////////////////////////////////////////////////////////////
  fillreweightingweights(reweightingweights);
  mela.resetInputEvent();
}

void Reweighting::computeP(int reweightinghypothesis, float& prob, bool useConstant) {
  prob = 1;
  float tmp = 0;
  setcouplings(reweightinghypothesis);
  if (get<0>(decayprocess) != TVar::Null) {
    setProcess(decayprocess);
    mela.computeP(tmp, useConstant);
    prob *= tmp;
  }
  setcouplings(reweightinghypothesis);
  if (get<0>(productionprocess) != TVar::Null) {
    setProcess(productionprocess);
    mela.computeProdP(tmp, useConstant);
    prob *= tmp;
  }
}

//this version requires the input event or current candidate to already be set in mela.
//it's protected and is called by the other ones.
void Reweighting::fillreweightingweights(vector<float>& reweightingweights) {

  reweightingweights.clear();
  assert(canreweight());

  float myprobability, probability;
  computeP(-1, myprobability, false);

  if (myprobability == 0) {
    for (int i = 0; i < nReweightingSamples; i++) {
      reweightingweights.push_back(0);
    }
    return;
  }

  for (int i = 0; i < nReweightingSamples; i++) {
    computeP(i, probability, false);
    float ratio = probability / myprobability;
    if (cutoffs[i] > 0 && ratio > cutoffs[i]) ratio = cutoffs[i] / (ratio*ratio);
    reweightingweights.push_back(ratio);
  }
}


#endif

