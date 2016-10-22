#include <ZZAnalysis/AnalysisStep/interface/MELAHypothesis.h>
#include <iostream>

using namespace std;


MELAHypothesis::MELAHypothesis(
  Mela* mela_,
  MELAOptionParser* opt_
  ) :
  mela(mela_),
  opt(opt_),
  optIsOwned(false)
{ reset(); }
MELAHypothesis::MELAHypothesis(
  Mela* mela_,
  string stropt
  ) :
  mela(mela_),
  optIsOwned(true)
{ opt = new MELAOptionParser(stropt); reset(); }

void MELAHypothesis::reset(){
  pME=0.; if (opt!=0) pME = opt->getDefaultME();
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

  reset();

  bool isGen = opt->isGen();
  MELACandidate* melaCand = mela->getCurrentCandidate();
  if (melaCand!=0){
    // Set the couplings
    // Comment: We have more couplings than we will ever need in the next 50 years!
    mela->differentiate_HWW_HZZ = opt->coupl_H.separateWWZZcouplings;
    for (unsigned int im=0; im<2; im++){
      //****Spin-0****//
      for (int ic=0; ic<(int)SIZE_HQQ; ic++) mela->selfDHqqcoupl[ic][im] = opt->coupl_H.Hqqcoupl[ic][im];
      for (int ic=0; ic<(int)SIZE_HGG; ic++) mela->selfDHggcoupl[ic][im] = opt->coupl_H.Hggcoupl[ic][im];
      // The first dimension (of size [nSupportedHiggses=2]) supports a second resonance present in MCFM
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
    mela->setMelaHiggsMassWidth(opt->hmass, opt->hwidth, 0);
    mela->setMelaHiggsMassWidth(opt->h2mass, opt->h2width, 1);

    // Set the process elemenets
    TVar::Process theProc = opt->proc;
    TVar::Production theProd = opt->prod;
    TVar::MatrixElement theME = opt->ME;
    if (
      isGen && theME==TVar::JHUGen && melaCand->getNAssociatedJets()<2
      ){
      if (theProd==TVar::Had_WH) theProd=TVar::Lep_WH;
      else if (theProd==TVar::Had_ZH) theProd=TVar::Lep_ZH;
    }
    if (
      isGen && theME==TVar::JHUGen && melaCand->getNAssociatedLeptons()<2
      ){
      if (theProd==TVar::Lep_WH) theProd=TVar::Had_WH;
      else if (theProd==TVar::Lep_ZH) theProd=TVar::Had_ZH;
    }
    mela->setProcess(theProc, theME, theProd);
    if (
      theProd==TVar::Lep_WH || theProd==TVar::Had_WH || theProd==TVar::Lep_ZH || theProd==TVar::Had_ZH || theProd == TVar::GammaH
      ||
      theProd==TVar::JJVBF || theProd==TVar::JJQCD || theProd==TVar::JQCD
      ||
      theProd==TVar::ttH || theProd==TVar::bbH
      ){
      if (theME != TVar::MCFM) mela->computeProdP(pME, !isGen);
      else mela->computeProdDecP(pME, !isGen);
    }
    else mela->computeP(pME, !isGen);

    if (!isGen){
      mela->getPAux(pAux);
    }
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
  reset();
  if (mela->getCurrentCandidate()!=0) mela->computePM4l(opt->superSyst, pME);
}

Float_t MELAHypothesis::getVal(METype valtype) const{
  if (valtype==UseME) return pME;
  else if (valtype==UsePAux) return pAux;
  else if (valtype==UsePConstant) return cMEAvg;
  else return -1;
}


