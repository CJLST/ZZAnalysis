#ifndef REWEIGHTING_CC
#define REWEIGHTING_CC


#include <ZZAnalysis/AnalysisStep/interface/reweighting.h>


ReweightingType Reweighting::reweightingtypefromstring(const std::string& reweightingtypestring) {
  if (reweightingtypestring == "none") return NoReweighting;
  else if (reweightingtypestring == "HVV_spin0") return HVV_spin0;
  else if (reweightingtypestring == "HVV_spin012") return HVV_spin012;
  else if (reweightingtypestring == "VBFHZZ_spin0") return VBFHZZ_spin0;
  else if (reweightingtypestring == "ZHZZ_spin0") return ZHZZ_spin0;
  else if (reweightingtypestring == "WHZZ_spin0") return WHZZ_spin0;
  assert(false);
  return NoReweighting;
}

int Reweighting::nReweightingSamplesFromType(ReweightingType reweightingtype_) {
  switch (reweightingtype_) {
    case NoReweighting: return 0;
    case HVV_spin0: return 7;
    case HVV_spin012: return 8;
    case VBFHZZ_spin0:
    case ZHZZ_spin0:
    case WHZZ_spin0:
      return 13;
    default: assert(false);
  }
}

Reweighting::Reweighting(
  Mela& mela_,
  std::string reweightingtypestring,
  int inputspin,
  std::vector<double> Hzzcouplings_real,
  std::vector<double> Hzzcouplings_imag,
  std::vector<double> Zvvcouplings_real,
  std::vector<double> Zvvcouplings_imag,
  std::vector<double> Gggcouplings_real,
  std::vector<double> Gggcouplings_imag,
  std::vector<double> Gvvcouplings_real,
  std::vector<double> Gvvcouplings_imag,
  std::vector<double> cutoffs_
  ) :
  myspin(inputspin),
  cutoffs(cutoffs_),
  fillingcouplingstree(false),
  mela(mela_),
  ghg2(mela.selfDHggcoupl[ghg2_index][0]),
  ghz1(mela.selfDHzzcoupl[0][ghz1_index][0]),
  ghz2(mela.selfDHzzcoupl[0][ghz2_index][0]),
  ghz4(mela.selfDHzzcoupl[0][ghz4_index][0]),
  ghz1_prime2(mela.selfDHzzcoupl[0][ghz1_prime2_index][0]),
  a1(mela.selfDGggcoupl[a1_index][0]),
  b5(mela.selfDGvvcoupl[b5_index][0]),
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
  couplingsfromvectors(myHvvcoupl, Hzzcouplings_real, Hzzcouplings_imag);
  couplingsfromvectors(myZvvcoupl, Zvvcouplings_real, Zvvcouplings_imag);
  couplingsfromvectors(myGggcoupl, Gggcouplings_real, Gggcouplings_imag);
  couplingsfromvectors(myGvvcoupl, Gvvcouplings_real, Gvvcouplings_imag);
}

void Reweighting::setmycouplings() {
  if (reweightingtype == NoReweighting) {
    return;
  } else if (reweightingtype == HVV_spin0 || reweightingtype == HVV_spin012) {
    productionprocess = make_tuple(TVar::nProcesses, TVar::JHUGen, TVar::ZZGG);
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
  } else if (reweightingtype == WHZZ_spin0 || reweightingtype == ZHZZ_spin0) {
    assert(myspin == 0);
    productionprocess = make_tuple(TVar::SelfDefine_spin0, TVar::JHUGen, getVHtype());
    decayprocess = make_tuple(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZINDEPENDENT);
  } else {
    assert(false);
  }

  spin = myspin;

  ghg2 = 1;
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
    productionprocess = make_tuple(TVar::nProcesses, TVar::JHUGen, TVar::ZZGG);
    spin = 0;
    ghg2 = 1;
    switch (reweightinghypothesis) {
    case 0: ghz1 = 1; break;                                     //0+m
    case 1: ghz2 = 1; break;                                     //0+h
    case 2: ghz4 = 1; break;                                     //0-
    case 3: ghz1_prime2 = 1; break;                              //L1
    case 4: ghz1 = 1; ghz2 = ghz2mix_decay; break;               //fa2=0.5
    case 5: ghz1 = 1; ghz4 = ghz4mix_decay; break;               //fa3=0.5
    case 6: ghz1 = 1; ghz1_prime2 = ghz1_prime2mix_decay; break; //fL1=0.5
    default: assert(false);
    }
    return;
  } else if (reweightingtype == HVV_spin012) {
    productionprocess = make_tuple(TVar::nProcesses, TVar::JHUGen, TVar::ZZGG);
    ghg2 = 1;
    switch (reweightinghypothesis) {
    case 0: ghz1 = 1; spin = 0; break;                                     //0+m
    case 1: ghz2 = 1; spin = 0; break;                                     //0+h
    case 2: ghz4 = 1; spin = 0; break;                                     //0-
    case 3: ghz1_prime2 = 1; spin = 0; break;                              //L1
    case 4: ghz1 = 1; ghz2 = ghz2mix_decay; spin = 0; break;               //fa2=0.5
    case 5: ghz1 = 1; ghz4 = ghz4mix_decay; spin = 0; break;               //fa3=0.5
    case 6: ghz1 = 1; ghz1_prime2 = ghz1_prime2mix_decay; spin = 0; break; //fL1=0.5
    case 7: a1 = 1; b5 = 1; spin = 2; break;                               //2b+
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
    switch (reweightinghypothesis) {
    case 0: ghz1 = 1; break;                                                                //0+m
    case 1: ghz2 = 1; break;                                                                //0+h
    case 2: ghz4 = 1; break;                                                                //0-
    case 3: ghz1_prime2 = 1; break;                                                         //L1
    case 4: ghz1 = 1; ghz2 = ghz2mix_decay; break;                                          //fa2_decay=0.5
    case 5: ghz1 = 1; ghz4 = ghz4mix_decay; break;                                          //fa3_decay=0.5
    case 6: ghz1 = 1; ghz1_prime2 = ghz1_prime2mix_decay; break;                            //fL1_decay=0.5
    case 7: ghz1 = 1; ghz2 = ghz2mix_VBF; break;                                            //fa2_VBF=0.5
    case 8: ghz1 = 1; ghz4 = ghz4mix_VBF; break;                                            //fa3_VBF=0.5
    case 9: ghz1 = 1; ghz1_prime2 = ghz1_prime2mix_VBF; break;                              //fL1_VBF=0.5
    case 10: ghz1 = 1; ghz2 = -sqrt(ghz2mix_VBF*ghz2mix_decay); break;                      //fa2_VBFdecay=-0.5
    case 11: ghz1 = 1; ghz4 = -sqrt(ghz4mix_VBF*ghz4mix_decay); break;                      //fa3_VBFdecay=-0.5
    case 12: ghz1 = 1; ghz1_prime2 = -sqrt(ghz1_prime2mix_VBF*ghz1_prime2mix_decay); break; //fL1_VBFdecay=-0.5
    default: assert(false);
    }
    return;
  } else if (reweightingtype == WHZZ_spin0) {
    decayprocess = make_tuple(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZINDEPENDENT);
    productionprocess = make_tuple(TVar::SelfDefine_spin0, TVar::JHUGen, getVHtype());
    spin = 0;
    switch (reweightinghypothesis) {
    case 0: ghz1 = 1; break;                                                                //0+m
    case 1: ghz2 = 1; break;                                                                //0+h
    case 2: ghz4 = 1; break;                                                                //0-
    case 3: ghz1_prime2 = 1; break;                                                         //L1
    case 4: ghz1 = 1; ghz2 = ghz2mix_decay; break;                                          //fa2_decay=0.5
    case 5: ghz1 = 1; ghz4 = ghz4mix_decay; break;                                          //fa3_decay=0.5
    case 6: ghz1 = 1; ghz1_prime2 = ghz1_prime2mix_decay; break;                            //fL1_decay=0.5
    case 7: ghz1 = 1; ghz2 = ghz2mix_WH; break;                                             //fa2_WH=0.5
    case 8: ghz1 = 1; ghz4 = ghz4mix_WH; break;                                             //fa3_WH=0.5
    case 9: ghz1 = 1; ghz1_prime2 = ghz1_prime2mix_WH; break;                               //fL1_WH=0.5
    case 10: ghz1 = 1; ghz2 = -sqrt(ghz2mix_WH*ghz2mix_decay); break;                       //fa2_WHdecay=-0.5
    case 11: ghz1 = 1; ghz4 = -sqrt(ghz4mix_WH*ghz4mix_decay); break;                       //fa3_WHdecay=-0.5
    case 12: ghz1 = 1; ghz1_prime2 = -sqrt(ghz1_prime2mix_WH*ghz1_prime2mix_decay); break;  //fL1_WHdecay=-0.5
    default: assert(false);
    }
    return;
  } else if (reweightingtype == ZHZZ_spin0) {
    decayprocess = make_tuple(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZINDEPENDENT);
    productionprocess = make_tuple(TVar::SelfDefine_spin0, TVar::JHUGen, getVHtype());
    spin = 0;
    switch (reweightinghypothesis) {
    case 0: ghz1 = 1; break;                                                                //0+m
    case 1: ghz2 = 1; break;                                                                //0+h
    case 2: ghz4 = 1; break;                                                                //0-
    case 3: ghz1_prime2 = 1; break;                                                         //L1
    case 4: ghz1 = 1; ghz2 = ghz2mix_decay; break;                                          //fa2_decay=0.5
    case 5: ghz1 = 1; ghz4 = ghz4mix_decay; break;                                          //fa3_decay=0.5
    case 6: ghz1 = 1; ghz1_prime2 = ghz1_prime2mix_decay; break;                            //fL1_decay=0.5
    case 7: ghz1 = 1; ghz2 = ghz2mix_ZH; break;                                             //fa2_ZH=0.5
    case 8: ghz1 = 1; ghz4 = ghz4mix_ZH; break;                                             //fa3_ZH=0.5
    case 9: ghz1 = 1; ghz1_prime2 = ghz1_prime2mix_ZH; break;                               //fL1_ZH=0.5
    case 10: ghz1 = 1; ghz2 = -sqrt(ghz2mix_ZH*ghz2mix_decay); break;                       //fa2_ZHdecay=-0.5
    case 11: ghz1 = 1; ghz4 = -sqrt(ghz4mix_ZH*ghz4mix_decay); break;                       //fa3_ZHdecay=-0.5
    case 12: ghz1 = 1; ghz1_prime2 = -sqrt(ghz1_prime2mix_ZH*ghz1_prime2mix_decay); break;  //fL1_ZHdecay=-0.5
    default: assert(false);
    }
    return;
  }
  assert(false);
}

void Reweighting::fillcouplingstree(TTree* t) {
  vector<double> ghz1_, ghz2_, ghz4_, ghz1_prime2_, a1_, b5_;
  vector<int> spin_v;
  t->Branch("spin", &spin_v);
  if (reweightingtype == HVV_spin0 || reweightingtype == HVV_spin012 || reweightingtype == VBFHZZ_spin0 || reweightingtype == ZHZZ_spin0 || reweightingtype == WHZZ_spin0) {
    t->Branch("ghz1Re", &ghz1_);
    t->Branch("ghz2Re", &ghz2_);
    t->Branch("ghz4Re", &ghz4_);
    t->Branch("ghz1_prime2Re", &ghz1_prime2_);
  }
  if (reweightingtype == HVV_spin012) {
    t->Branch("a1Re", &a1_);
    t->Branch("b5Re", &b5_);
  }

  fillingcouplingstree = true;
  cerr << "There are about to be " << nReweightingSamples << " warnings about no daughters.  Please ignore them." << endl;
  for (int i = 0; i < nReweightingSamples; i++) {
    reset_SelfDCouplings();
    setcouplings(i);
    spin_v.push_back(spin);
    if (spin == 0) {
      ghz1_.push_back(ghz1);
      ghz2_.push_back(ghz2);
      ghz4_.push_back(ghz4);
      ghz1_prime2_.push_back(ghz1_prime2);
    }
    else {
      ghz1_.push_back(0);
      ghz2_.push_back(0);
      ghz4_.push_back(0);
      ghz1_prime2_.push_back(0);
    }

    if (spin == 2) {
      a1_.push_back(a1);
      b5_.push_back(b5);
    }
    else {
      a1_.push_back(0);
      b5_.push_back(0);
    }
  }
  t->Fill();
  fillingcouplingstree = false;
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

  mela.setCurrentCandidate(candidate);
  fillreweightingweights(reweightingweights);
  mela.resetInputEvent();
}

void Reweighting::computeP(int reweightinghypothesis, float& prob, bool useConstant) {
  prob = 1;
  float tmp = 0;
  if (get<0>(decayprocess) != TVar::nProcesses) {
    setcouplings(reweightinghypothesis);
    setProcess(decayprocess);
    mela.computeP(tmp, useConstant);
    prob *= tmp;
  }
  if (get<0>(productionprocess) != TVar::nProcesses) {
    setcouplings(reweightinghypothesis);
    setProcess(productionprocess);
    mela.computeProdP(tmp, useConstant);
    prob *= tmp;
  }
}

//this version requires the input event or current candidate to already be set in mela.
//it's protected and is called by the other ones.
void Reweighting::fillreweightingweights(vector<float>& reweightingweights) {

  reweightingweights.clear();

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
    if (cutoffs[i] > 0 && ratio > cutoffs[i]) ratio = cutoffs[i]*cutoffs[i] / ratio;
    reweightingweights.push_back(ratio);
  }
}

TVar::Production Reweighting::getVHtype() {
  //Very rudimentary check.  If there are too many jets and leptons
  // it will just give up.  Also doesn't check that the ids make sense
  // (e.g. e+ mu- --> VHLep even though this is nonsense, and even if set to WH)
  if (reweightingtype == WHZZ_spin0 || reweightingtype == ZHZZ_spin0) {
    //find if it's leptonic or hadronic
    auto candidate = mela.getCurrentCandidate();
    if (!candidate) {
      if (fillingcouplingstree) { //so that we can call setcouplings() to get the values while filling the tree
        return TVar::Had_WH;      //doesn't actually matter what we return here
      } else {
        cerr << "!candidate" << endl;
        assert(false);
      }
    }
    int jets = candidate->getNAssociatedJets();
    int leptons = candidate->getNAssociatedLeptons();
    int neutrinos = candidate->getNAssociatedNeutrinos();
    if (jets >= 2 && leptons+neutrinos <= 1) {
      if (reweightingtype == WHZZ_spin0) 
        return TVar::Had_WH;
      else
        return TVar::Had_ZH;
    }
    if (leptons+neutrinos >= 2 && jets <= 1) {
      if (reweightingtype == WHZZ_spin0) 
        return TVar::Lep_WH;
      else
        return TVar::Lep_ZH;
    }
    cerr << "Warning: ambiguous VH event, could be hadronic or leptonic." << endl
         << "It's possible that a more precise check of the ids would give an unambiguous result." << endl
         << "Also you could try getting this information right from the LHE event with mother information" << endl
         << "Associated ids:" << endl;
    for (int id : candidate->getAssociatedParticleIds()) cerr << id << " ";
    cerr << endl;
  } else {
    cerr << "Should not call getVHtype for reweightingtype = " << reweightingtype << endl;
  }
  assert(false);
  return TVar::ZZINDEPENDENT;
}



#endif
