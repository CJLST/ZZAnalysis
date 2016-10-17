/** \class MELABranch
*
*
*  \author N. Amapane - Torino
*  \author U. Sarica - JHU
*/
#ifndef MELABRANCH_H
#define MELABRANCH_H

#include "ExtendedBranch.h"
#include "MELAHypothesis.h"

namespace BranchHelpers{

  class MELABranch : public ExtendedBranch<Float_t>{

  public:

    MELAHypothesis* targetP;
    MELAHypothesis* originP;

    MELABranch(
      TTree* theTree_, TString bname_, Float_t defVal_,
      MELAHypothesis* targetP_, MELAHypothesis* originP_ = 0
      ) :
      ExtendedBranch(theTree_, bname_, defVal_, false),
      targetP(targetP_), originP(originP_)
    {}
    virtual ~MELABranch(){}

    void setVal(){ Float_t tmp=1.; if (targetP!=0) tmp *= targetP->value; if (originP!=0) tmp /= originP->value; setVal(tmp); }

  };

}

#endif
