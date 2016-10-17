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


namespace BranchHelpers{
  enum BranchTypes{
    bBool, bChar, bShort,
    bInt, bLong, bFloat, bDouble,
    nBranchTypes
  };


  template<typename varType> class ExtendedBranch{

  protected:

    TTree* theTree;
    varType* defVal;

  public:

    TString bname;
    varType value;
    std::vector<varType> valueArray;
    BranchHelpers::BranchTypes btype;
    Bool_t isVector;

  public:

    void reset(){
      valueArray.clear();
      if (defVal!=0) value = *defVal;
    }

    void setVal(varType inVal){
      if (isVector) value.push_back(inVal);
      else value = inVal;
    }
    varType getDefVal() const{
      if (!isVector && defVal!=0) return *defVal;
      else return 0;
    }
    varType getVal() const{
      if (!isVector) return value;
      else if (valueArray.size()>0) return valueArray.at(valueArray.size()-1);
      else return 0;
    }
    std::vector<varType> getArray() const{ return valueArray; }

    TBranch* getBranch(){ return theTree->GetBranch(bname); }

    ExtendedBranch(TTree* theTree_, TString bname_, varType defVal_, Bool_t isVector_=false) : theTree(theTree_), defVal(0), bname(bname_), isVector(isVector_)
    {
      getBranchType();
      createBranch(&defVal_);
    }
    virtual ~ExtendedBranch(){
      delete defVal;
    }


  protected:

    void getBranchType(){}
    /*
    void getBranchType(){
      if (dynamic_cast<Bool_t*>(&value)) btype = BranchHelpers::bBool;
      else if (dynamic_cast<Char_t*>(&value)) btype = BranchHelpers::bChar;
      else if (dynamic_cast<Short_t*>(&value)) btype = BranchHelpers::bShort;

      else if (dynamic_cast<Int_t*>(&value)) btype = BranchHelpers::bInt;
      else if (dynamic_cast<Long64_t*>(&value)) btype = BranchHelpers::bLong;
      else if (dynamic_cast<Float_t*>(&value)) btype = BranchHelpers::bFloat;
      else if (dynamic_cast<Double_t*>(&value)) btype = BranchHelpers::bDouble;

      else{
        std::cerr << "ExtendedBranch::getBranchType: Could not determine the branch type for branch " << bname << std::endl;
        assert(0);
      }
    }
    */
    void createBranch(varType* defVal_){
      if (!isVector && defVal_!=0){
        if (defVal!=0) delete defVal; defVal = new varType;
        *defVal = *defVal_;
      }
      reset();

      TString leaftypename = "";
      if (!isVector){
        if (btype == BranchHelpers::bBool) leaftypename = "O";
        else if (btype == BranchHelpers::bChar) leaftypename = "B";
        else if (btype == BranchHelpers::bShort) leaftypename = "S";
        else if (btype == BranchHelpers::bInt) leaftypename = "I";
        else if (btype == BranchHelpers::bLong) leaftypename = "L";
        else if (btype == BranchHelpers::bFloat) leaftypename = "F";
        else if (btype == BranchHelpers::bDouble) leaftypename = "D";
      }
      if (theTree!=0 && bname!=""){
        if (leaftypename=="") theTree->Branch(bname, &valueArray);
        else theTree->Branch(bname, &value, (bname + "/" + leaftypename));
      }
    }

  };

}

#endif
