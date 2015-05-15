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
  */
  for(int i=0;i<7;i++)nBranches[i]=0;

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

void HZZ4lNtupleFactory::Book(TString *names, int type, int nBranchessToFill, float *defaultValues){
  for(int i=0;i<nBranchessToFill;i++)Book(names[i],defaultValues[i],type);
}

void HZZ4lNtupleFactory::Book(TString branchName, float defaultValue, int varType){
  //should add a protection for existing branches
  //kBool,kShort,kInt,kFloat,kVectorFloat
  TString leafName=branchName.Data();
  nameVector[varType].push_back(branchName);

  if(varType!=kVectorFloat)defaultVector[varType].push_back(defaultValue);
  if(varType==kBool){
    leafName.Append("/O");
    boolVector[nBranches[kBool]]=(bool)defaultValue;
    _outTree->Branch(branchName.Data(),&boolVector[nBranches[varType]],leafName.Data());
  }
  else if(varType==kShort){
    leafName.Append("/S");
    shortVector[nBranches[kShort]]=(Short_t)defaultValue;
    _outTree->Branch(branchName.Data(),&shortVector[nBranches[varType]],leafName.Data());
  }
  else if(varType==kInt){
    leafName.Append("/I");
    intVector[nBranches[kInt]]=(Int_t)defaultValue;
    _outTree->Branch(branchName.Data(),&(intVector[nBranches[varType]]),leafName.Data());
  }
  else if(varType==kChar){
    leafName.Append("/B");
    charVector[nBranches[kChar]]=(char)defaultValue;
    _outTree->Branch(branchName.Data(),&charVector[nBranches[varType]],leafName.Data());
  }
    else if(varType==kLong){
    leafName.Append("/L");
    longVector[nBranches[kLong]]=(Long64_t)defaultValue;
    _outTree->Branch(branchName.Data(),&longVector[nBranches[varType]],leafName.Data());
  }

  else if(varType==kFloat){
    leafName.Append("/F");
    floatVector[nBranches[kFloat]]=defaultValue;
    _outTree->Branch(branchName.Data(),&floatVector[nBranches[varType]],leafName.Data());
  }
  else if(varType==kVectorFloat){
    std::vector<float> tempVector;
    vectorVector[nBranches[kVectorFloat]]=tempVector;
    _outTree->Branch(branchName.Data(),&vectorVector[nBranches[varType]]);
  }
  nBranches[varType]++;
}

int HZZ4lNtupleFactory::SetVariables(TString *branchName, double *value, int nBranchessToSet){
  int setVars=0;
  for(int i=0;i<nBranchessToSet;i++)setVars+=SetVariable(branchName[i],value[i]);
  return setVars;
}

//int HZZ4lNtupleFactory::SetVariables(TString *branchName, float *value, int nBranchessToSet){
//  int setVars=0;
//  for(int i=0;i<nBranchessToSet;i++)setVars+=SetVariable(branchName[i],value[i]);
//  return setVars;
//}

int HZZ4lNtupleFactory::SetVariableLong(TString branchName, Long64_t value){
int iname=999;
    for(int in=0;in<(int)nameVector[kLong].size();in++){
      if(branchName.CompareTo(nameVector[kLong].at(in))==0){
        iname=in;
        break;
      }
    }
 if(iname<(int)nameVector[kLong].size())longVector[iname]=value;
 return 0;

}

int HZZ4lNtupleFactory::SetVariable(TString branchName, double value){
  //cout<<"Setting "<<branchName.Data()<<endl;
  bool found=false;
  int iname=0,itype=0;
  for(int it=0;it<=kVectorFloat && !found;it++){
    for(int in=0;in<(int)nameVector[it].size();in++){
      if(branchName.CompareTo(nameVector[it].at(in))==0){
        found=true;
        iname=in;
        itype=it;
        break;
      }
    }
  }
  if(!found){
    bool notBooked=true;
    for(int i=0;i<(int)NotBookedBranches.size();i++){
      if(branchName.CompareTo(NotBookedBranches.at(i).Data())==0){
        notBooked=false;
        break;
      }
    }
    if(notBooked){
      cout<<"Warning!!! Variable "<<branchName.Data()<<" is not booked in NtupleFactory"<<endl;
      NotBookedBranches.push_back(branchName.Data());
    }
    return 0;
  }
  if(itype==kBool) boolVector[iname]=(bool)value;
  else if(itype==kShort) shortVector[iname]=(Short_t)value;
  else if(itype==kInt) intVector[iname]=(int)value;
  else if(itype==kChar) charVector[iname]=(char)value;
  else if(itype==kLong) longVector[iname]=(Long64_t)value;
  else if(itype==kFloat) floatVector[iname]=value;
  else if(itype==kVectorFloat){
    vectorVector[iname].push_back((float)value);
  }else return -1; //should be impossible to get here...
  return 1;
}

void HZZ4lNtupleFactory::InitializeVariables()
{

  for(int i=0;i<nBranches[kBool];i++){boolVector[i]=defaultVector[kBool].at(i);}
  for(int i=0;i<nBranches[kShort];i++){shortVector[i]=defaultVector[kShort].at(i);}
  for(int i=0;i<nBranches[kInt];i++){intVector[i]=defaultVector[kInt].at(i);}
  for(int i=0;i<nBranches[kChar];i++){charVector[i]=defaultVector[kChar].at(i);}
  for(int i=0;i<nBranches[kLong];i++){longVector[i]=defaultVector[kLong].at(i);}
  for(int i=0;i<nBranches[kFloat];i++){floatVector[i]=defaultVector[kFloat].at(i);}
  for(int i=0;i<nBranches[kVectorFloat];i++){vectorVector[i].clear();}
}


void HZZ4lNtupleFactory::createNewCandidate()
{
  _firstZStored = false;
  _LeptonIndex = 1;
  _LeptonIsoIndex = 1;

  return;
}

void HZZ4lNtupleFactory::FillZInfo(Float_t ZMass, Float_t ZPt, short ZFlav)
{
  if(!_firstZStored){
    SetVariable("Z1Mass",ZMass);
    SetVariable("Z1Pt",ZPt);
    SetVariable("Z1Flav",ZFlav);
    _firstZStored = true;
  }
  else{
    SetVariable("Z2Mass",ZMass);
    SetVariable("Z2Pt",ZPt);
    SetVariable("Z2Flav",ZFlav);
  }

  return;
}


void HZZ4lNtupleFactory::FillHGenInfo(const math::XYZTLorentzVector pH, float w)
{
  SetVariable("GenHMass",pH.M());
  SetVariable("GenHPt",pH.Pt());

  SetVariable("HqTMCweight",w);

  return;
}

void HZZ4lNtupleFactory::FillZGenInfo(const math::XYZTLorentzVector pZ1, const math::XYZTLorentzVector pZ2)
{
  SetVariable("GenZ1Mass", pZ1.M());
  SetVariable("GenZ1Pt", pZ1.Pt());

  SetVariable("GenZ2Mass", pZ2.M());
  SetVariable("GenZ2Pt", pZ2.Pt());

  return;
}

void HZZ4lNtupleFactory::FillLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, Short_t Lep3Id, Short_t Lep4Id, 
					const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2, 
					const math::XYZTLorentzVector Lep3, const math::XYZTLorentzVector Lep4,float weight)
{
  SetVariable("GenLep1Pt",Lep1.Pt());
  SetVariable("GenLep1Eta",Lep1.Eta());
  SetVariable("GenLep1Phi",Lep1.Phi());
  SetVariable("GenLep1Id",Lep1Id);

  SetVariable("GenLep2Pt",Lep2.Pt());
  SetVariable("GenLep2Eta",Lep2.Eta());
  SetVariable("GenLep2Phi",Lep2.Phi());
  SetVariable("GenLep2Id",Lep2Id);

  SetVariable("GenLep3Pt",Lep3.Pt());
  SetVariable("GenLep3Eta",Lep3.Eta());
  SetVariable("GenLep3Phi",Lep3.Phi());
  SetVariable("GenLep3Id",Lep3Id);

  SetVariable("GenLep4Pt",Lep4.Pt());
  SetVariable("GenLep4Eta",Lep4.Eta());
  SetVariable("GenLep4Phi",Lep4.Phi());
  SetVariable("GenLep4Id",Lep4Id);

  SetVariable("dataMCWeight",weight);
  
  return;
}

void HZZ4lNtupleFactory::FillAssocLepGenInfo(std::vector<const reco::Candidate *>& AssocLeps)
{
  if (AssocLeps.size() >= 1) {
    SetVariable("genAssocLep1Pt ",AssocLeps.at(0)->p4().Pt());
    SetVariable("genAssocLep1Eta",AssocLeps.at(0)->p4().Eta());
    SetVariable("genAssocLep1Phi",AssocLeps.at(0)->p4().Phi());
    SetVariable("genAssocLep1Id ",AssocLeps.at(0)->pdgId());
  }
  if (AssocLeps.size() >= 2) {
    SetVariable("genAssocLep2Pt ",AssocLeps.at(1)->p4().Pt());
    SetVariable("genAssocLep2Eta",AssocLeps.at(1)->p4().Eta());
    SetVariable("genAssocLep2Phi",AssocLeps.at(1)->p4().Phi());
    SetVariable("genAssocLep2Id ",AssocLeps.at(1)->pdgId());
  }

  return;
}

void HZZ4lNtupleFactory::FillExtraLepInfo(int extraLeptonIndex, bool extraLeptonExists, const reco::CandidatePtr ExtraLep)
{
  Float_t Pt         = extraLeptonExists ? ExtraLep->pt()    : -9999. ;
  Float_t Eta        = extraLeptonExists ? ExtraLep->eta()   : -9999. ;
  Float_t Phi        = extraLeptonExists ? ExtraLep->phi()   : -9999. ;
  Int_t   LepId      = extraLeptonExists ? ExtraLep->pdgId() :     0  ;
//   Float_t SIP        = extraLeptonExists ? userdatahelpers::getUserFloat(&*ExtraLep,"SIP")              : -9999. ;
//   Bool_t  isID       = extraLeptonExists ? userdatahelpers::getUserFloat(&*ExtraLep,"isID")             :     0  ;
//   Float_t BDT        = extraLeptonExists ? userdatahelpers::getUserFloat(&*ExtraLep,"BDT")              : -9999. ;
//   Char_t  missingHit = extraLeptonExists ? (char)userdatahelpers::getUserFloat(&*ExtraLep,"missingHit") :     0  ;
//   Float_t chargedHadIso = extraLeptonExists ? userdatahelpers::getUserFloat(&*ExtraLep,"PFChargedHadIso") : -9999. ;
//   Float_t neutralHadIso = extraLeptonExists ? userdatahelpers::getUserFloat(&*ExtraLep,"PFNeutralHadIso") : -9999. ;
//   Float_t photonIso     = extraLeptonExists ? userdatahelpers::getUserFloat(&*ExtraLep,"PFPhotonIso")     : -9999. ;
//   Float_t combRelIsoPF  = extraLeptonExists ? userdatahelpers::getUserFloat(&*ExtraLep,"combRelIsoPF")    : -9999. ;

  switch(extraLeptonIndex){

  case 1:
    SetVariable("ExtraLep1Pt"        ,Pt);
    SetVariable("ExtraLep1Eta"       ,Eta);
    SetVariable("ExtraLep1Phi"       ,Phi);
    SetVariable("ExtraLep1LepId"     ,LepId);
//     SetVariable("ExtraLep1SIP       ,SIP);
//     SetVariable("ExtraLep1isID      ,isID);
//     SetVariable("ExtraLep1BDT       ,BDT);
//     SetVariable("ExtraLep1missingHit,missingHit);
//     SetVariable("ExtraLep1chargedHadIso,chargedHadIso);
//     SetVariable("ExtraLep1neutralHadIso,neutralHadIso);
//     SetVariable("ExtraLep1photonIso    ,photonIso);
//     SetVariable("ExtraLep1combRelIsoPF ,combRelIsoPF);
    break;

  case 2:
    SetVariable("ExtraLep2Pt"        ,Pt);
    SetVariable("ExtraLep2Eta"       ,Eta);
    SetVariable("ExtraLep2Phi"       ,Phi);
    SetVariable("ExtraLep2LepId"     ,LepId);
//     SetVariable("ExtraLep2SIP       ,SIP);
//     SetVariable("ExtraLep2isID      ,isID);
//     SetVariable("ExtraLep2BDT       ,BDT);
//     SetVariable("ExtraLep2missingHit,missingHit);
//     SetVariable("ExtraLep2chargedHadIso,chargedHadIso);
//     SetVariable("ExtraLep2neutralHadIso,neutralHadIso);
//     SetVariable("ExtraLep2photonIso    ,photonIso);
//     SetVariable("ExtraLep2combRelIsoPF ,combRelIsoPF);
    break;

  case 3:
    SetVariable("ExtraLep3Pt"        ,Pt);
    SetVariable("ExtraLep3Eta"       ,Eta);
    SetVariable("ExtraLep3Phi"       ,Phi);
    SetVariable("ExtraLep3LepId"     ,LepId);
//     SetVariable("ExtraLep3SIP       ,SIP;
//     SetVariable("ExtraLep3isID      ,isID;
//     SetVariable("ExtraLep3BDT       ,BDT;
//     SetVariable("ExtraLep3missingHit,missingHit;
//     SetVariable("ExtraLep3chargedHadIso,chargedHadIso;
//     SetVariable("ExtraLep3neutralHadIso,neutralHadIso;
//     SetVariable("ExtraLep3photonIso    ,photonIso;
//     SetVariable("ExtraLep3combRelIsoPF ,combRelIsoPF;
    break;

  default:
    std::cout << "Error in indexing the extra leptons ! Will abort..." << std::endl;
    assert(0);
  }
  //Fill the tree once for each extra lepton (no vectors anymore)
  //SetVariable("outTree->Fill();
  return;
}


