/** \class MELABranch
*
*
*  \author N. Amapane - Torino
*  \author U. Sarica - JHU
*/
#ifndef MELAHYPOTHESIS_H
#define MELAHYPOTHESIS_H

#include <ZZMatrixElement/MELA/interface/Mela.h>


class MELAHypothesis{

protected:

  Mela* mela;

  SpinZeroCouplings* coupl_H;
  SpinOneCouplings* coupl_Zp;
  SpinTwoCouplings* coupl_X;

  TVar::Process proc;
  TVar::Production prod;
  TVar::MatrixElement me;

  float hmass;
  float hwidth;
  float h2mass;
  float h2width;

public:

  float pME;
  float pAux;
  float cMEAvg;

  MELAHypothesis(
    Mela* mela_,

    SpinZeroCouplings* coupl_H_,
    SpinOneCouplings* coupl_Zp_,
    SpinTwoCouplings* coupl_X_,

    TVar::Process proc_,
    TVar::Production prod_,
    TVar::MatrixElement me_,

    float hmass_=-1.,
    float hwidth_=0.,
    float h2mass_=-1.,
    float h2width_=0.
    );
  virtual ~MELAHypothesis(){}

  float compute(bool isGen, MELACandidate* cand=0);

};


#endif
