#include <ZZAnalysis/AnalysisStep/interface/MELAHypothesis.h>
#include <iostream>


#define CALL_CLASSMEMBER_OBJ(object,ptrToMember)  ((object).*(ptrToMember))
#define CALL_CLASSMEMBER_REF(object,ptrToMember)  ((object)->*(ptrToMember))

using namespace std;


MELAHypothesis::MELAHypothesis(
  Mela* mela_,
  MELAOptionParser* opt_
  ) :
  mela(mela_),
  opt(opt_),
  optIsOwned(false),
  hasMaximizationClients(false)
{ reset(); }
MELAHypothesis::MELAHypothesis(
  Mela* mela_,
  string stropt
  ) :
  mela(mela_),
  optIsOwned(true),
  hasMaximizationClients(false)
{
  opt = new MELAOptionParser(stropt);
  reset();
}

void MELAHypothesis::reset(){
  isUpdated=false;
  pME = (opt ? opt->getDefaultME() : 0);
  pAux=1.;
  cMEAvg=1.;
}

void MELAHypothesis::computeP(MELACandidate* cand){
  if (cand!=0) mela->setCurrentCandidate(cand);
  computeP();
}
void MELAHypothesis::computeP(unsigned int index){
  mela->setCurrentCandidateFromIndex(index);
  computeP();
}
void MELAHypothesis::computeP(){
  if (opt->usePM4L()){ computePM4l(); return; }
  else if (opt->usePMaVJJ()){ computePMAVJJ(); return; }
  else if (opt->usePropagator()){ computePropagator(); return; }

  if (isUpdated && !hasMaximizationClients) return; // Avoid further computations if there are no clients
  reset(); // Note: Sets isUpdated=false.

  bool isGen = opt->isGen();
  MELACandidate* melaCand = mela->getCurrentCandidate();
  if (melaCand!=0){
    // Set the couplings
    // Comment: We have more couplings than we will ever need in the next 50 years!
    mela->differentiate_HWW_HZZ = opt->coupl_H.separateWWZZcouplings;
    for (unsigned int im=0; im<2; im++){
      //****Spin-0****//
      // First resonance parameters
      for (int ic=0; ic<(int)SIZE_HQQ; ic++){
        mela->selfDHqqcoupl[0][ic][im] = opt->coupl_H.Hqqcoupl[ic][im];
        mela->selfDHbbcoupl[0][ic][im] = opt->coupl_H.Hbbcoupl[ic][im];
        mela->selfDHttcoupl[0][ic][im] = opt->coupl_H.Httcoupl[ic][im];
        mela->selfDHb4b4coupl[0][ic][im] = opt->coupl_H.Hb4b4coupl[ic][im];
        mela->selfDHt4t4coupl[0][ic][im] = opt->coupl_H.Ht4t4coupl[ic][im];
      }
      for (int ic=0; ic<(int)SIZE_HGG; ic++){
        mela->selfDHggcoupl[0][ic][im] = opt->coupl_H.Hggcoupl[ic][im];
        mela->selfDHg4g4coupl[0][ic][im] = opt->coupl_H.Hg4g4coupl[ic][im];
      }
      for (int ic=0; ic<(int)SIZE_HVV; ic++){
        mela->selfDHzzcoupl[0][ic][im] = opt->coupl_H.Hzzcoupl[ic][im];
        mela->selfDHwwcoupl[0][ic][im] = opt->coupl_H.Hwwcoupl[ic][im];
      }
      if (im==0){ // Only real numbers
        for (int iq=0; iq<(int)SIZE_HVV_CQSQ; iq++){
          for (int ic=0; ic<(int)SIZE_HVV_LAMBDAQSQ; ic++){
            mela->selfDHzzLambda_qsq[0][ic][iq] = opt->coupl_H.HzzLambda_qsq[ic][iq];
            mela->selfDHwwLambda_qsq[0][ic][iq] = opt->coupl_H.HwwLambda_qsq[ic][iq];
          }
          mela->selfDHzzCLambda_qsq[0][iq] = opt->coupl_H.HzzCLambda_qsq[iq];
          mela->selfDHwwCLambda_qsq[0][iq] = opt->coupl_H.HwwCLambda_qsq[iq];
        }
      }

      for (int ic=0; ic<(int)SIZE_HQQ; ic++){
        mela->selfDHqqcoupl[1][ic][im] = opt->coupl_H.H2qqcoupl[ic][im];
        mela->selfDHbbcoupl[1][ic][im] = opt->coupl_H.H2bbcoupl[ic][im];
        mela->selfDHttcoupl[1][ic][im] = opt->coupl_H.H2ttcoupl[ic][im];
        mela->selfDHb4b4coupl[1][ic][im] = opt->coupl_H.H2b4b4coupl[ic][im];
        mela->selfDHt4t4coupl[1][ic][im] = opt->coupl_H.H2t4t4coupl[ic][im];
      }
      for (int ic=0; ic<(int)SIZE_HGG; ic++){
        mela->selfDHggcoupl[1][ic][im] = opt->coupl_H.H2ggcoupl[ic][im];
        mela->selfDHg4g4coupl[1][ic][im] = opt->coupl_H.H2g4g4coupl[ic][im];
      }
      for (int ic=0; ic<(int)SIZE_HVV; ic++){
        mela->selfDHzzcoupl[1][ic][im] = opt->coupl_H.H2zzcoupl[ic][im];
        mela->selfDHwwcoupl[1][ic][im] = opt->coupl_H.H2wwcoupl[ic][im];
      }
      if (im==0){ // Only real numbers
        for (int iq=0; iq<(int)SIZE_HVV_CQSQ; iq++){
          for (int ic=0; ic<(int)SIZE_HVV_LAMBDAQSQ; ic++){
            mela->selfDHzzLambda_qsq[1][ic][iq] = opt->coupl_H.H2zzLambda_qsq[ic][iq];
            mela->selfDHwwLambda_qsq[1][ic][iq] = opt->coupl_H.H2wwLambda_qsq[ic][iq];
          }
          mela->selfDHzzCLambda_qsq[1][iq] = opt->coupl_H.H2zzCLambda_qsq[iq];
          mela->selfDHwwCLambda_qsq[1][iq] = opt->coupl_H.H2wwCLambda_qsq[iq];
        }
      }
      //****Spin-1****//
      for (int ic=0; ic<(int)SIZE_ZQQ; ic++) mela->selfDZqqcoupl[ic][im] = opt->coupl_Zp.Zqqcoupl[ic][im];
      for (int ic=0; ic<(int)SIZE_ZVV; ic++) mela->selfDZvvcoupl[ic][im] = opt->coupl_Zp.Zvvcoupl[ic][im];
      //****Spin-2****//
      for (int ic=0; ic<(int)SIZE_GQQ; ic++) mela->selfDGqqcoupl[ic][im] = opt->coupl_X.Gqqcoupl[ic][im];
      for (int ic=0; ic<(int)SIZE_GGG; ic++) mela->selfDGggcoupl[ic][im] = opt->coupl_X.Gggcoupl[ic][im];
      for (int ic=0; ic<(int)SIZE_GVV; ic++) mela->selfDGvvcoupl[ic][im] = opt->coupl_X.Gvvcoupl[ic][im];
    }
    // That is a lot of them!

    // Set the masses
    // Notice that >=-1 is used. -1 is the value to disable a resonance in the options
    if (opt->hmass>=-1.) mela->setMelaHiggsMassWidth(opt->hmass, opt->hwidth, 0);
    if (opt->h2mass>=-1.) mela->setMelaHiggsMassWidth(opt->h2mass, opt->h2width, 1);

    // Set the process elemenets
    TVar::Process theProc = opt->proc;
    TVar::Production theProd = opt->prod;
    TVar::MatrixElement theME = opt->ME;
    // In case the VH ME has a mismatch with the event, switch the production.
    unsigned int nactivejets=0;
    unsigned int nactiveleps=0;
    for (auto const& part : melaCand->getAssociatedLeptons()){ if (part->passSelection) nactiveleps++; }
    for (auto const& part : melaCand->getAssociatedJets()){ if (part->passSelection) nactivejets++; }
    if (
      isGen && nactivejets<2 && nactiveleps>=2
      ){
      if (theProd==TVar::Had_WH) theProd=TVar::Lep_WH;
      else if (theProd==TVar::Had_ZH) theProd=TVar::Lep_ZH;
      else if (theProd==TVar::Had_WH_S) theProd=TVar::Lep_WH_S;
      else if (theProd==TVar::Had_ZH_S) theProd=TVar::Lep_ZH_S;
      else if (theProd==TVar::Had_WH_TU) theProd=TVar::Lep_WH_TU;
      else if (theProd==TVar::Had_ZH_TU) theProd=TVar::Lep_ZH_TU;
      else if (theProd==TVar::JJEW || theProd==TVar::JJEW_S || theProd==TVar::JJEW_TU){
        MELAParticle* aV=nullptr;
        for (MELAParticle* tmp:melaCand->getAssociatedSortedVs()){
          if (tmp!=0 && tmp->passSelection && (PDGHelpers::isAZBoson(tmp->id) || PDGHelpers::isAWBoson(tmp->id))){
            if (tmp->getNDaughters()==2 &&
              tmp->getDaughter(0)->passSelection && PDGHelpers::isALepton(tmp->getDaughter(0)->id)
              &&
              tmp->getDaughter(1)->passSelection && PDGHelpers::isALepton(tmp->getDaughter(1)->id)
              ){
              aV=tmp;
              break;
            }
          }
        }
        if (aV){
          if (PDGHelpers::isAZBoson(aV->id)){
            if (theProd==TVar::JJEW) theProd=TVar::Lep_ZH;
            else if (theProd==TVar::JJEW_S) theProd=TVar::Lep_ZH_S;
            else if (theProd==TVar::JJEW_TU) theProd=TVar::Lep_ZH_TU;
          }
          else if (PDGHelpers::isAWBoson(aV->id)){
            if (theProd==TVar::JJEW) theProd=TVar::Lep_WH;
            else if (theProd==TVar::JJEW_S) theProd=TVar::Lep_WH_S;
            else if (theProd==TVar::JJEW_TU) theProd=TVar::Lep_WH_TU;
          }
        }
      }
    }
    if (
      isGen && nactiveleps<2 && nactivejets>=2
      ){
      if (theProd==TVar::Lep_WH) theProd=TVar::Had_WH;
      else if (theProd==TVar::Lep_ZH) theProd=TVar::Had_ZH;
      else if (theProd==TVar::Lep_WH_S) theProd=TVar::Had_WH_S;
      else if (theProd==TVar::Lep_ZH_S) theProd=TVar::Had_ZH_S;
      else if (theProd==TVar::Lep_WH_TU) theProd=TVar::Had_WH_TU;
      else if (theProd==TVar::Lep_ZH_TU) theProd=TVar::Had_ZH_TU;
    }
    mela->setProcess(theProc, theME, theProd);
    typedef void (Mela::*MelaComputePFcn)(float&, bool);
    MelaComputePFcn computePFcn=nullptr;
    if (
      theProd==TVar::Lep_WH || theProd==TVar::Had_WH || theProd==TVar::Lep_ZH || theProd==TVar::Had_ZH || theProd==TVar::JJVBF || theProd==TVar::JJEW || theProd==TVar::JJQCD
      || theProd==TVar::Lep_WH_S || theProd==TVar::Had_WH_S || theProd==TVar::Lep_ZH_S || theProd==TVar::Had_ZH_S || theProd==TVar::JJVBF_S || theProd==TVar::JJEW_S || theProd==TVar::JJQCD_S
      || theProd==TVar::Lep_WH_TU || theProd==TVar::Had_WH_TU || theProd==TVar::Lep_ZH_TU || theProd==TVar::Had_ZH_TU || theProd==TVar::JJVBF_TU || theProd==TVar::JJEW_TU || theProd==TVar::JJQCD_TU
      || theProd == TVar::GammaH
      || theProd==TVar::JQCD
      || theProd==TVar::ttH || theProd==TVar::bbH
      ){
      if (theME != TVar::MCFM) computePFcn = &Mela::computeProdP;
      else computePFcn = &Mela::computeProdDecP;
    }
    else computePFcn = &Mela::computeP;
    float pMEprior=pME;
    CALL_CLASSMEMBER_REF(mela, computePFcn)(pME, !isGen);
    if (pME<0. && pME!=pMEprior){
      TVar::VerbosityLevel bkpverbosity = mela->getVerbosity();
      mela->setVerbosity(TVar::DEBUG_VERBOSE);
      CALL_CLASSMEMBER_REF(mela, computePFcn)(pME, !isGen);
      TUtil::PrintCandidateSummary(mela->getCurrentCandidate());
      mela->setVerbosity(bkpverbosity);
    }
    if (!isGen){
      mela->getPAux(pAux);
      mela->getConstant(cMEAvg);
    }

    isUpdated = true;
  }
}

void MELAHypothesis::computePM4l(MELACandidate* cand){
  if (cand!=0) mela->setCurrentCandidate(cand);
  computePM4l();
}
void MELAHypothesis::computePM4l(unsigned int index){
  mela->setCurrentCandidateFromIndex(index);
  computePM4l();
}
void MELAHypothesis::computePM4l(){
  if (isUpdated && !hasMaximizationClients) return; // Avoid further computations if there are no clients
  reset(); // Note: Sets isUpdated=false.
  if (mela->getCurrentCandidate()!=0){
    // Override the ME and the production
    mela->setProcess(opt->proc, TVar::JHUGen, TVar::ZZGG);
    mela->computePM4l(opt->superSyst, pME);
    isUpdated = true;
  }
}

void MELAHypothesis::computePMAVJJ(MELACandidate* cand){
  if (cand!=0) mela->setCurrentCandidate(cand);
  computePMAVJJ();
}
void MELAHypothesis::computePMAVJJ(unsigned int index){
  mela->setCurrentCandidateFromIndex(index);
  computePMAVJJ();
}
void MELAHypothesis::computePMAVJJ(){
  if (isUpdated && !hasMaximizationClients) return; // Avoid further computations if there are no clients
  reset(); // Note: Sets isUpdated=false.
  if (mela->getCurrentCandidate()!=0){
    // Override the ME and the production
    mela->setProcess(opt->proc, opt->ME, opt->prod);
    mela->computeDijetConvBW(pME, opt->usePMaVJJTrue());
    isUpdated = true;
  }
}

void MELAHypothesis::computePropagator(MELACandidate* cand){
  if (cand!=0) mela->setCurrentCandidate(cand);
  computePropagator();
}
void MELAHypothesis::computePropagator(unsigned int index){
  mela->setCurrentCandidateFromIndex(index);
  computePropagator();
}
void MELAHypothesis::computePropagator(){
  if (isUpdated && !hasMaximizationClients) return; // Avoid further computations if there are no clients
  reset(); // Note: Sets isUpdated=false.
  MELACandidate* melaCand = mela->getCurrentCandidate();
  if (melaCand!=0){
    if (opt->hmass>=-1.) mela->setMelaHiggsMassWidth(opt->hmass, opt->hwidth, 0);
    mela->getXPropagator(opt->propScheme, pME);
    isUpdated = true;
  }
}


Float_t MELAHypothesis::getVal(METype valtype) const{
  if (valtype==UseME) return pME;
  else if (valtype==UsePAux) return pAux;
  else if (valtype==UsePConstant) return cMEAvg;
  else return -1;
}


