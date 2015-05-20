#include <cassert>
#include <iostream>

#include "HZZ4lNtupleFactory.h"

/*
bool addKinRefit = false;
bool addVtxFit = false;
bool writePhotons = false;  // Write photons in the tree. Must be set also in HZZ4lNtupleMaker.cc
*/
using namespace std;

HZZ4lNtupleFactory::HZZ4lNtupleFactory(TTree* outTree_input)
{
  //---- create output tree ----
  _outTree = outTree_input;
  /*
  cout<<"Factory!"<<endl;
  for(int i=0;i<99;i++){
     intVector[i]=0;
  shortVector[i]=0;
   boolVector[i]=0;
   longVector[i]=0;
   charVector[i]=0;
   floatVector[i]=0;
  //vectorVector[i]=0;
  }
  //std::vector<float> defaultVector[6];
  //std::vector<TString> nameVector[7];

  cout<<"end factory"<<endl;
  
  for(int i=0;i<7;i++)nBranches[i]=0;
*/
  //InitializeVariables();
  
  _firstZStored = false;
  _LeptonIndex = 1;
  _LeptonIsoIndex = 1;

}

///--- Destructor ---
HZZ4lNtupleFactory::~HZZ4lNtupleFactory()
{
}

///---- Write an event to TTree ----
void HZZ4lNtupleFactory::FillEvent()
{
  _outTree->Fill();
  InitializeVariables(); // Reset all values and clean vectors
}

///---- Write to a text file branches declaration ----
void HZZ4lNtupleFactory::DumpBranches(TString filename) const
{
  //----- symply use MakeClass
  _outTree->MakeClass(filename);
  return;
}

void HZZ4lNtupleFactory::Book(TString name, Float_t &variable){
  TString leafname=name.Data();
  leafname.Append("/F");
  defaultsFloat[&variable] = variable; // defaultT is a map<T*, value>
  _outTree->Branch(name.Data(), &variable, leafname.Data());
}
void HZZ4lNtupleFactory::Book(TString name, Int_t &value){
  TString leafname=name.Data();
  leafname.Append("/I");
  defaultsInt[&value] = value; // defaultT is a map<T*, value>
  _outTree->Branch(name.Data(), &value, leafname.Data());
}
void HZZ4lNtupleFactory::Book(TString name, Bool_t &value){
  TString leafname=name.Data();
  leafname.Append("/O");
  defaultsBool[&value] = value; // defaultT is a map<T*, value>
  _outTree->Branch(name.Data(), &value, leafname.Data());
}
void HZZ4lNtupleFactory::Book(TString name, Short_t &value){
  TString leafname=name.Data();
  leafname.Append("/S");
  defaultsShort[&value] = value; // defaultT is a map<T*, value>
  _outTree->Branch(name.Data(), &value, leafname.Data());
}
void HZZ4lNtupleFactory::Book(TString name, Long64_t &value){
  TString leafname=name.Data();
  leafname.Append("/L");
  defaultsLong[&value] = value; // defaultT is a map<T*, value>
  _outTree->Branch(name.Data(), &value, leafname.Data());
}
void HZZ4lNtupleFactory::Book(TString name, Char_t &value){
  TString leafname=name.Data();
  leafname.Append("/B");
  defaultsChar[&value] = value; // defaultT is a map<T*, value>
  _outTree->Branch(name.Data(), &value, leafname.Data());
}
void HZZ4lNtupleFactory::Book(TString name, std::vector<float> &value){
  defaultsVector[&value] = value; // defaultT is a map<T*, value>
  _outTree->Branch(name.Data(), &value);
}


void HZZ4lNtupleFactory::InitializeVariables()
{
 for (auto it = defaultsFloat.begin(); it != defaultsFloat.end(); ++it ) *(it->first) = it->second;
 for (auto it = defaultsInt.begin(); it != defaultsInt.end(); ++it ) *(it->first) = it->second;
 for (auto it = defaultsBool.begin(); it != defaultsBool.end(); ++it ) *(it->first) = it->second;
 for (auto it = defaultsShort.begin(); it != defaultsShort.end(); ++it ) *(it->first) = it->second;
 for (auto it = defaultsLong.begin(); it != defaultsLong.end(); ++it ) *(it->first) = it->second;
 for (auto it = defaultsChar.begin(); it != defaultsChar.end(); ++it ) *(it->first) = it->second;
 for (auto it = defaultsVector.begin(); it != defaultsVector.end(); ++it ) it->first->clear();//*(it->first)->clear();// = it->second;

/*
  for(int i=0;i<nBranches[kBool];i++){boolVector[i]=defaultVector[kBool].at(i);}
  for(int i=0;i<nBranches[kShort];i++){shortVector[i]=defaultVector[kShort].at(i);}
  for(int i=0;i<nBranches[kInt];i++){intVector[i]=defaultVector[kInt].at(i);}
  for(int i=0;i<nBranches[kChar];i++){charVector[i]=defaultVector[kChar].at(i);}
  for(int i=0;i<nBranches[kLong];i++){longVector[i]=defaultVector[kLong].at(i);}
  for(int i=0;i<nBranches[kFloat];i++){floatVector[i]=defaultVector[kFloat].at(i);}
  for(int i=0;i<nBranches[kVectorFloat];i++){vectorVector[i].clear();}
*/
}


void HZZ4lNtupleFactory::createNewCandidate()
{
  _firstZStored = false;
  _LeptonIndex = 1;
  _LeptonIsoIndex = 1;

  return;
}



