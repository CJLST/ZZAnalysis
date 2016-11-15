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

  // Type checking, adaptation of suggestion by Michael Aaron Safyan at http://stackoverflow.com/questions/2728078/is-there-an-easy-way-to-check-a-fundamental-type
  template<BranchTypes VAL> struct constant_type_value{
    static const BranchTypes type = VAL;
  };
  typedef constant_type_value<nBranchTypes> no_type;
  typedef constant_type_value<bBool> bool_type;
  typedef constant_type_value<bChar> char_type;
  typedef constant_type_value<bShort> short_type;
  typedef constant_type_value<bInt> int_type;
  typedef constant_type_value<bLong> long_type;
  typedef constant_type_value<bFloat> float_type;
  typedef constant_type_value<bDouble> double_type;

  template<typename T> struct hasType : public no_type{};
  template<> struct hasType<Bool_t> : public bool_type{};
  template<> struct hasType<Char_t> : public char_type{};
  template<> struct hasType<Short_t> : public short_type{};
  template<> struct hasType<Int_t> : public int_type{};
  template<> struct hasType<Long64_t> : public long_type{};
  template<> struct hasType<Float_t> : public float_type{};
  template<> struct hasType<Double_t> : public double_type{};

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

    void setValue(varType inVal){
      if (isVector) valueArray.push_back(inVal);
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

    virtual void Print(){
      std::cout << "**********\n";
      std::cout << "ExtendedBranch summary for " << bname << ":\n";
      std::cout << "\tType=" << btype << "\n";
      if (isVector){
        std::cout << "\tValues = ";
        for (unsigned int iv=0; iv<valueArray.size(); iv++) std::cout << valueArray.at(iv) << " ";
        std::cout << "\n";
      }
      else std::cout << "\tValue = " << value << " (default = " << (defVal!=0 ? *defVal : 0) << ")\n";
      std::cout << "\tTree address: " << theTree << "\n";
      std::cout << "**********";
      std::cout << std::endl;
    }

    ExtendedBranch(TTree* theTree_, TString bname_, varType defVal_, Bool_t isVector_=false) : theTree(theTree_), defVal(0), bname(bname_), isVector(isVector_)
    {
      getBranchType();
      createBranch(&defVal_);
    }
    virtual ~ExtendedBranch(){
      delete defVal;
    }


  protected:

    void getBranchType(){
      btype = hasType<varType>::type;
      if (btype == nBranchTypes){
        std::cerr << "ExtendedBranch::getBranchType: Could not determine the branch type for branch " << bname << std::endl;
        assert(0);
      }
    }
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
