/** \class ExtendedBranch
 *
 *
 *  \author N. Amapane - Torino
 *  \author U. Sarica - JHU
 */
#ifndef EXTENDEDBRANCH_H
#define EXTENDEDBRANCH_H

#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <cassert>
#include "TString.h"
#include "TTree.h"


typedef std::vector<Bool_t> vectorBool_t;
typedef std::vector<Char_t> vectorChar_t;
typedef std::vector<Short_t> vectorShort_t;
typedef std::vector<Int_t> vectorInt_t;
typedef std::vector<Float_t> vectorFloat_t;
typedef std::vector<Double_t> vectorDouble_t;

namespace BranchHelpers{
  enum BranchTypes{
    bBool, bChar, bShort,
    bInt, bLong, bFloat, bDouble,

    bVectorBool, bVectorChar, bVectorShort,
    bVectorInt, bVectorFloat, bVectorDouble,

    nBranchTypes
  };


  template<typename varType> class ExtendedBranch{

  public:

    void reset(){
      if (
        btype == BranchHelpers::bVectorBool
        ||
        btype == BranchHelpers::bVectorChar
        ||
        btype == BranchHelpers::bVectorShort
        ||
        btype == BranchHelpers::bVectorInt
        ||
        btype == BranchHelpers::bVectorFloat
        ||
        btype == BranchHelpers::bVectorDouble
        ) value.clear();
      else if(defVal!=0) value = *defVal;
    }

    template<typename inType> void setVal(inType inVal){
      if (
        btype == BranchHelpers::bVectorBool
        ||
        btype == BranchHelpers::bVectorChar
        ||
        btype == BranchHelpers::bVectorShort
        ||
        btype == BranchHelpers::bVectorInt
        ||
        btype == BranchHelpers::bVectorFloat
        ||
        btype == BranchHelpers::bVectorDouble
        ) value.push_back(inVal);
      else value = inVal;
    }
    varType getVal(){ return value; }
    varType& getRef(){ return value; }
    varType* getPtr(){ return &value; }
    TBranch* getBranch(){ return theTree->GetBranch(bname); }
    BranchHelpers::BranchTypes getType(){ return btype; }

    ExtendedBranch(TTree* theTree_, TString bname_, varType defVal_) : theTree(theTree_), bname(bname_), defVal(0)
    {
      btype = getBranchType();
      createBranch(&defVal);
    }
    ExtendedBranch(TTree* theTree_, TString bname_) : theTree(theTree_), bname(bname_), defVal(0)
    {
      btype = getBranchType();
      createBranch((varType*)0);
    }
    ExtendedBranch(){ /* Do nothing */}
    virtual ~ExtendedBranch(){
      delete defVal;
    }


  protected:

    TTree* theTree;
    TString bname;
    varType* defVal;
    varType value;
    BranchHelpers::BranchTypes btype;

    BranchHelpers::BranchTypes getBranchType(){
      if (dynamic_cast<Bool_t*>(&value)) btype = BranchHelpers::bBool;
      else if (dynamic_cast<Char_t*>(&value)) btype = BranchHelpers::bChar;
      else if (dynamic_cast<Short_t*>(&value)) btype = BranchHelpers::bShort;

      else if (dynamic_cast<Int_t*>(&value)) btype = BranchHelpers::bInt;
      else if (dynamic_cast<Long64_t*>(&value)) btype = BranchHelpers::bLong;
      else if (dynamic_cast<Float_t*>(&value)) btype = BranchHelpers::bFloat;
      else if (dynamic_cast<Double_t*>(&value)) btype = BranchHelpers::bDouble;

      else if (dynamic_cast<vectorBool_t*>(&value)) btype = BranchHelpers::bVectorBool;
      else if (dynamic_cast<vectorChar_t*>(&value)) btype = BranchHelpers::bVectorChar;
      else if (dynamic_cast<vectorShort_t*>(&value)) btype = BranchHelpers::bVectorShort;

      else if (dynamic_cast<vectorInt_t*>(&value)) btype = BranchHelpers::bVectorInt;
      else if (dynamic_cast<vectorFloat_t*>(&value)) btype = BranchHelpers::bVectorFloat;
      else if (dynamic_cast<vectorDouble_t*>(&value)) btype = BranchHelpers::bVectorDouble;

      else{
        std::cerr << "ExtendedBranch::getBranchType: Could not determine the branch type for branch " << bname << std::endl;
        assert(0);
      }
    }
    void createBranch(varType* defVal_){
      if (defVal_!=0){
        if (defVal!=0) delete defVal; defVal = new varType;
        *defVal = *defVal_;
      }
      reset();

      TString leaftypename = "";
      if (btype == BranchHelpers::bBool) leaftypename = "O";
      else if (btype == BranchHelpers::bChar) leaftypename = "B";
      else if (btype == BranchHelpers::bShort) leaftypename = "S";
      else if (btype == BranchHelpers::bInt) leaftypename = "I";
      else if (btype == BranchHelpers::bLong) leaftypename = "L";
      else if (btype == BranchHelpers::bFloat) leaftypename = "F";
      else if (btype == BranchHelpers::bDouble) leaftypename = "D";

      if (leaftypename=="") theTree->Branch(bname, &value);
      else theTree->Branch(bname, &value, (bname + "/" + leaftypename));
    }

  };

}

#endif
