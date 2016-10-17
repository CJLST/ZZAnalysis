#include <ZZAnalysis/AnalysisStep/interface/MELAHypothesis.h>
#include <iostream>


MELAHypothesis::MELAHypothesis(
  Mela* mela_,
  METype calc_type_,

  SpinZeroCouplings* coupl_H_,
  SpinOneCouplings* coupl_Zp_,
  SpinTwoCouplings* coupl_X_,

  TVar::Process proc_,
  TVar::Production prod_,
  TVar::MatrixElement me_,

  float hmass_,
  float hwidth_,
  float h2mass_,
  float h2width_
  ) :
  mela(mela_),
  calc_type(calc_type_),

  coupl_H(coupl_H_),
  coupl_Zp(coupl_Zp_),
  coupl_X(coupl_X_),

  proc(proc_),
  prod(prod_),
  me(me_),

  hmass(hmass_),
  hwidth(hwidth_),
  h2mass(h2mass_),
  h2width(h2width_),

  pME(-1.),
  pAux(1.),
  cMEAvg(1.),

  value(0.)
{}

void MELAHypothesis::reset(){
  pME=-1.;
  pAux=1.;
  cMEAvg=1.;
  value=0.;
}

void MELAHypothesis::computeP(MELACandidate* cand, bool isGen){
  reset();
  if (cand!=0) mela->setCurrentCandidate(cand);
  MELACandidate* melaCand = mela->getCurrentCandidate();

  if (melaCand!=0){
    // Set the couplings
    // Comment: We have more couplings than we will ever need in the next 50 years!
    mela->differentiate_HWW_HZZ = coupl_H->separateWWZZcouplings;
    for (unsigned int im=0; im<2; im++){
      //****Spin-0****//
      mela->selfDHvvcoupl_freenorm[im] = coupl_H->Hvvcoupl_freenorm[im];
      for (int ic=0; ic<(int)SIZE_HQQ; ic++) mela->selfDHqqcoupl[ic][im] = coupl_H->Hqqcoupl[ic][im];
      for (int ic=0; ic<(int)SIZE_HGG; ic++) mela->selfDHggcoupl[ic][im] = coupl_H->Hggcoupl[ic][im];
      // The first dimension (of size [nSupportedHiggses=2]) supports a second resonance present in MCFM
      for (int ic=0; ic<(int)SIZE_HVV; ic++){
        mela->selfDHzzcoupl[0][ic][im] = coupl_H->Hzzcoupl[ic][im];
        mela->selfDHwwcoupl[0][ic][im] = coupl_H->Hwwcoupl[ic][im];
      }
      if (im==0){ // Only real numbers
        for (int iq=0; iq<(int)SIZE_HVV_CQSQ; iq++){
          for (int ic=0; ic<(int)SIZE_HVV_LAMBDAQSQ; ic++){
            mela->selfDHzzLambda_qsq[0][ic][iq] = coupl_H->HzzLambda_qsq[ic][iq];
            mela->selfDHwwLambda_qsq[0][ic][iq] = coupl_H->HwwLambda_qsq[ic][iq];
          }
          mela->selfDHzzCLambda_qsq[0][iq] = coupl_H->HzzCLambda_qsq[iq];
          mela->selfDHwwCLambda_qsq[0][iq] = coupl_H->HwwCLambda_qsq[iq];
        }
      }
      for (int ic=0; ic<(int)SIZE_HVV; ic++){
        mela->selfDHzzcoupl[1][ic][im] = coupl_H->H2zzcoupl[ic][im];
        mela->selfDHwwcoupl[1][ic][im] = coupl_H->H2wwcoupl[ic][im];
      }
      if (im==0){ // Only real numbers
        for (int iq=0; iq<(int)SIZE_HVV_CQSQ; iq++){
          for (int ic=0; ic<(int)SIZE_HVV_LAMBDAQSQ; ic++){
            mela->selfDHzzLambda_qsq[1][ic][iq] = coupl_H->H2zzLambda_qsq[ic][iq];
            mela->selfDHwwLambda_qsq[1][ic][iq] = coupl_H->H2wwLambda_qsq[ic][iq];
          }
          mela->selfDHzzCLambda_qsq[1][iq] = coupl_H->H2zzCLambda_qsq[iq];
          mela->selfDHwwCLambda_qsq[1][iq] = coupl_H->H2wwCLambda_qsq[iq];
        }
      }
      //****Spin-1****//
      for (int ic=0; ic<(int)SIZE_ZQQ; ic++) mela->selfDZqqcoupl[ic][im] = coupl_Zp->Zqqcoupl[ic][im];
      for (int ic=0; ic<(int)SIZE_ZVV; ic++) mela->selfDZvvcoupl[ic][im] = coupl_Zp->Zvvcoupl[ic][im];
      //****Spin-2****//
      for (int ic=0; ic<(int)SIZE_GQQ; ic++) mela->selfDGqqcoupl[ic][im] = coupl_X->Gqqcoupl[ic][im];
      for (int ic=0; ic<(int)SIZE_GGG; ic++) mela->selfDGggcoupl[ic][im] = coupl_X->Gggcoupl[ic][im];
      for (int ic=0; ic<(int)SIZE_GVV; ic++) mela->selfDGvvcoupl[ic][im] = coupl_X->Gvvcoupl[ic][im];
    }
    // That is a lot of them!

    // Set the masses
    mela->setMelaHiggsMassWidth(hmass, hwidth, 0);
    mela->setMelaHiggsMassWidth(h2mass, h2width, 1);

    // Set the process elemenets
    mela->setProcess(proc, me, prod);
    if (
      prod==TVar::Lep_WH || prod==TVar::Had_WH || prod==TVar::Lep_ZH || prod==TVar::Had_ZH || prod == TVar::GammaH
      ||
      prod==TVar::JJVBF || prod==TVar::JJQCD || prod==TVar::JQCD
      ||
      prod==TVar::ttH || prod==TVar::bbH
      ){
      if (me != TVar::MCFM) mela->computeProdP(pME, !isGen);
      else mela->computeProdDecP(pME, !isGen);
    }
    else mela->computeP(pME, !isGen);

    if (!isGen){
      mela->getPAux(pAux);
    }
  }

  setVal();
}

void MELAHypothesis::computePM4l(MELACandidate* cand, TVar::SuperMelaSyst syst){
  reset();
  if (cand!=0) mela->setCurrentCandidate(cand);
  MELACandidate* melaCand = mela->getCurrentCandidate();

  if (melaCand!=0) mela->computePM4l(syst, pME);
  setVal();
}

void MELAHypothesis::setVal(){
  if (calc_type==MELAHypothesis::UseME) value = pME;
  else if (calc_type==MELAHypothesis::UsePAux) value = pAux;
  else if (calc_type==MELAHypothesis::UsePConstant) value = cMEAvg;
  else{
    std::cerr << "MELAHypothesis:: : Could not determine the type (current type=" << calc_type << ") of the requested value. Returning -1!" << std::endl;
    value = -1;
  }
}


