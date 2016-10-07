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
    bBool, bChar,
    bShort, bInt, bLong, bFloat, bDouble,

    bVectorBool, bVectorChar, bVectorShort,
    bVectorInt, bVectorFloat, bVectorDouble,

    nBranchTypes
  };


  template<typename varType> class ExtendedBranch{

  public:

    void setTree(TTree* theTree_){ theTree = theTree_; }
    BranchHelpers::BranchTypes getBranchType(){
      if (dynamic_cast<Bool_t>(value)) return BranchHelpers::bBool;
      else if (dynamic_cast<Char_t>(value)) return BranchHelpers::bChar;
      else if (dynamic_cast<Short_t>(value)) return BranchHelpers::bShort;

      else if (dynamic_cast<Int_t>(value)) return BranchHelpers::bInt;
      else if (dynamic_cast<Long64_t>(value)) return BranchHelpers::bLong;
      else if (dynamic_cast<Float_t>(value)) return BranchHelpers::bFloat;
      else if (dynamic_cast<Double_t>(value)) return BranchHelpers::bDouble;

      else if (dynamic_cast<vectorBool_t>(value)) return BranchHelpers::bVectorBool;
      else if (dynamic_cast<vectorChar_t>(value)) return BranchHelpers::bVectorChar;
      else if (dynamic_cast<vectorShort_t>(value)) return BranchHelpers::bVectorShort;

      else if (dynamic_cast<vectorInt_t>(value)) return BranchHelpers::bVectorInt;
      else if (dynamic_cast<vectorFloat_t>(value)) return BranchHelpers::bVectorFloat;
      else if (dynamic_cast<vectorDouble_t>(value)) return BranchHelpers::bVectorDouble;

      else return BranchHelpers::nBranchTypes;
    }
    void reset(){
      BranchHelpers::BranchTypes btype = getBranchType();
      assert(btype!=BranchHelpers::nBranchTypes);
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
      else value = defVal;
    }
    void createBranch(TString bname_, varType defVal_){
      bname = bname_;
      defVal = defVal_;
      reset();
      BranchHelpers::BranchTypes btype = getBranchType(); // Assert is already called in reset

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
    template<typename inType> void setVal(inType inVal){
      BranchHelpers::BranchTypes btype = getBranchType();
      assert(btype!=BranchHelpers::nBranchTypes);
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

    ExtendedBranch(TTree* theTree_, TString bname_, varType defVal_){
      setTree(theTree_);
      createBranch(bname_, defVal);
    }
    virtual ~ExtendedBranch(){}


  protected:

    TTree* theTree;
    TString bname;
    varType value;
    varType defVal;

  };

}

#endif
