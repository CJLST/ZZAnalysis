/** \class MELABranch
*
*
*  \author N. Amapane - Torino
*  \author U. Sarica - JHU
*/
#ifndef MELABRANCH_H
#define MELABRANCH_H

#include "ExtendedBranch.h"
#include "MELAComputation.h"


namespace BranchHelpers{

  class MELABranch : public ExtendedBranch<Float_t>{

  public:

    MELAComputation* computer;
    MELAHypothesis::METype valtype;

    MELABranch(
      TTree* theTree_, TString bname_, Float_t defVal_,
      MELAComputation* computer_, TTree* failedTree_=0
      );
    virtual ~MELABranch(){}

    // This function should be run after MELAComputation classes are all "update"d.
    void setVal();

    virtual void Print();

  };

}

#endif
