/** \class MELABranch
*
*
*  \author N. Amapane - Torino
*  \author U. Sarica - JHU
*/
#ifndef MELABRANCH_H
#define MELABRANCH_H

#include <ZZMatrixElement/MELA/interface/Mela.h>
#include "ExtendedBranch.h"
#include "MELAHypothesis.h"


class MELABranch : public ExtendedBranch<Float_t>{

public:

  MELAHypothesis* targetP;
  MELAHypothesis* originP;

  MELABranch(
    TTree* theTree_, TString bname_, Float_t defVal_,
    MELAHypothesis* targetP_, MELAHypothesis* originP_ = 0
    ) :
    ExtendedBranch(theTree_, bname_, defVal_),
    targetP(targetP_), originP(originP_)
  {}
  virtual ~MELABranch(){}

  void compute(bool isGen, MELACandidate* cand=0){ value=1.; if (targetP!=0) value *= targetP->compute(isGen, cand); if (originP!=0) value /= originP->compute(isGen, cand); }

};


#endif
