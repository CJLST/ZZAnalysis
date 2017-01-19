#ifndef HZZ4lNtupleFactory_h
#define HZZ4lNtupleFactory_h

#include <vector>
#include <unordered_map>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>
#include <ZZAnalysis/AnalysisStep/interface/MELABranch.h>
#include <ZZAnalysis/AnalysisStep/interface/MELACluster.h>

using namespace BranchHelpers;

class HZZ4lNtupleFactory{
  
 
 protected:
  
 public:
  HZZ4lNtupleFactory(TTree* outTree_input, TTree* failedTree_input=0);
  ~HZZ4lNtupleFactory();

  enum varTypes {kBool,kShort,kInt,kChar,kLong,kFloat,kVectorFloat}; 

  //Fill the tree and initialize the variables
  void FillEvent(bool passed=true);
  
  //Fill the tree without resetting the variables
  void FillCurrentTree(bool passed=true);
  void DumpBranches(TString filename) const;
  
  void Book(TString branchName, Float_t &value, bool putinfailedtree=false);
  void Book(TString branchName, Char_t &value, bool putinfailedtree=false);
  void Book(TString branchName, Int_t &value, bool putinfailedtree=false);
  void Book(TString branchName, Bool_t &value, bool putinfailedtree=false);
  void Book(TString branchName, Long64_t &value, bool putinfailedtree=false);
  void Book(TString branchName, std::vector<float> &value, bool putinfailedtree=false);
  void Book(TString branchName, Short_t &value, bool putinfailedtree=false);
  void Book(TString branchName, std::vector<short> &value, bool putinfailedtree=false);
  void Book(TString branchName, std::vector<char> &value, bool putinfailedtree=false);
  void Book(TString branchName, std::vector<bool> &value, bool putinfailedtree=false);
  void BookMELABranches(MELAOptionParser* me_opt, bool isGen, MELAComputation* computer_);
  //void Book(TString branchName, float defaultValue=0,int varType=kFloat);

  //void Book(TString *names, int type, int nVarsToFill, float *defaultValues);
  //int SetVariable(TString varName, double value){return 0;}
  //int SetVariableLong(TString varName, Long64_t value);
  //int SetVariables(TString *varName, double *value, int nVars){return 0;}
  //int SetVariable(TString varName, float value){return SetVariable(varName, (double)value);}
  //int SetVariables(TString *varName, float *value, int nVars){return 0;}
  
  //void createNewCandidate();
  void InitializeVariables();

  //default fillers (old way of filling, only valid if branch names are not changed)
  //void FillHGenInfo(const math::XYZTLorentzVector Hp, float w);
  //void FillZGenInfo(const math::XYZTLorentzVector pZ1, const math::XYZTLorentzVector pZ2);
  //void FillLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, Short_t Lep3Id, Short_t Lep4Id, 
  //  const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2, const math::XYZTLorentzVector Lep3, 
  //  const math::XYZTLorentzVector Lep4, float weight);
  //void FillAssocLepGenInfo(std::vector<const reco::Candidate *>& AssocLeps);
  //void FillZInfo(Float_t ZMass, Float_t ZPt, short ZFlav);
  //void FillExtraLepInfo(int extraLeptonIndex, bool extraLeptonExists, const reco::CandidatePtr ExtraLep);

  std::vector<MELABranch*>* getRecoMELABranches();
  std::vector<MELABranch*>* getLHEMELABranches();


 private:

  TTree* _outTree;
  TTree* _failedTree;
  // MELABranches: Only branches are owned. Every other related object is owned by HZZ4lNtupleMaker.
  std::vector<MELABranch*> recome_branches;
  std::vector<MELABranch*> lheme_branches;

  /*
  int intVector[99];
  Short_t shortVector[99];
  Bool_t boolVector[99];
  Long64_t longVector[99];
  char charVector[99];
  Float_t floatVector[299];
  std::vector<float> vectorVector[99];
  */
  /*
  std::vector<float> defaultVector[6];
  std::vector<TString> nameVector[7];
    Int_t nBranches[7];
  std::vector<TString> NotBookedBranches;
  //std::vector<int> typeVector;
  */
  
  std::unordered_map<Float_t*, Float_t> defaultsFloat;
  std::unordered_map<Int_t*, Int_t> defaultsInt;
  std::unordered_map<Short_t*, Short_t> defaultsShort;
  std::unordered_map<Long64_t*, Long64_t> defaultsLong;
  std::unordered_map<Bool_t*, Bool_t> defaultsBool;
  std::unordered_map<Char_t*, Char_t> defaultsChar;
  std::unordered_map<std::vector<float>*, std::vector<float>> defaultsVectorFloat;
  std::unordered_map<std::vector<short>*, std::vector<short>> defaultsVectorShort;
  std::unordered_map<std::vector<char>*, std::vector<char>> defaultsVectorChar;
  std::unordered_map<std::vector<bool>*, std::vector<bool>> defaultsVectorBool;
  
  //bool _firstZStored;
  //int _LeptonIndex;
  //int _LeptonIsoIndex;

  //void InitializeBranches();

};

#endif
