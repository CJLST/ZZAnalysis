#ifndef HZZ4lNtupleFactory_h
#define HZZ4lNtupleFactory_h

#include <vector>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>

class HZZ4lNtupleFactory{
  
 
 protected:
  
 public:
  HZZ4lNtupleFactory(TTree* outTree_input);
  ~HZZ4lNtupleFactory();

  enum varTypes {kBool,kShort,kInt,kChar,kLong,kFloat,kVectorFloat}; 

  //Fill the tree and initialize the variables
  void FillEvent();
  
  //Fill the tree without resetting the variables
  void FillCurrentTree(){_outTree->Fill();}
  void DumpBranches(TString filename) const;
  
  void Book(TString branchName, float defaultValue=0,int varType=kFloat);
  void Book(TString *names, int type, int nVarsToFill, float *defaultValues);
  int SetVariable(TString varName, double value);
  int SetVariableLong(TString varName, Long64_t value);
  int SetVariables(TString *varName, double *value, int nVars);
  //int SetVariable(TString varName, float value){return SetVariable(varName, (double)value);}
  //int SetVariables(TString *varName, float *value, int nVars);
  
  void createNewCandidate();
  void InitializeVariables();

  //default fillers (old way of filling, only valid if branch names are not changed)
  void FillHGenInfo(const math::XYZTLorentzVector Hp, float w);
  void FillZGenInfo(const math::XYZTLorentzVector pZ1, const math::XYZTLorentzVector pZ2);
  void FillLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, Short_t Lep3Id, Short_t Lep4Id, 
    const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2, const math::XYZTLorentzVector Lep3, 
    const math::XYZTLorentzVector Lep4, float weight);
  void FillAssocLepGenInfo(std::vector<const reco::Candidate *>& AssocLeps);
  void FillZInfo(Float_t ZMass, Float_t ZPt, short ZFlav, Float_t Z1MassRefit);
  void FillExtraLepInfo(int extraLeptonIndex, bool extraLeptonExists, const reco::CandidatePtr ExtraLep);

 private:

  TTree* _outTree;
  
  int intVector[99];
  Short_t shortVector[99];
  Bool_t boolVector[99];
  Long64_t longVector[99];
  char charVector[99];
  Float_t floatVector[299];
  std::vector<float> vectorVector[99];
  
  std::vector<float> defaultVector[6];
  std::vector<TString> nameVector[7];
  Int_t nBranches[7];
  //std::vector<int> typeVector;
  
  bool _firstZStored;
  int _LeptonIndex;
  int _LeptonIsoIndex;

  void InitializeBranches();

};

#endif
