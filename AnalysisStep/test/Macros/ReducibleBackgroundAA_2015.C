#include <string>
#include <iostream>

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

#include <ZZAnalysis/AnalysisStep/src/bitops.cc>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>

using namespace std;

enum FinalState {fs4mu=0, fs4e=1, fs2e2mu=2, fs2mu2e=3};
const int nFinalStates = 4;
Float_t fsROSSS[nFinalStates] = { 1.22, 0.97, 1.30, 0.98 };
string fsLabelForSS[nFinalStates] = {
  "Z1->mu+mu- + mumu(SS)",
  "Z1->e+e- + ee(SS)",
  "Z1->e+e- + mumu(SS)",
  "Z1->mu+mu- + ee(SS)",
};

TH1F* h1D_FRmu_EB = 0;
TH1F* h1D_FRmu_EE = 0;
TH1F* h1D_FRel_EB = 0;
TH1F* h1D_FRel_EE = 0;

Float_t fakeRate13TeV(Float_t LepPt, Float_t LepEta, Int_t LepID) {
  Float_t myLepPt = LepPt>=80. ? 79. : LepPt;
  Int_t   myLepID = abs(LepID);
  if(myLepID==11){
    if(fabs(LepEta)<1.449)
      return h1D_FRel_EB->GetBinContent(h1D_FRel_EB->GetXaxis()->FindBin(myLepPt));
    else
      return h1D_FRel_EE->GetBinContent(h1D_FRel_EE->GetXaxis()->FindBin(myLepPt));
  }else if(myLepID==13){
    if(fabs(LepEta)<1.2)
      return h1D_FRmu_EB->GetBinContent(h1D_FRmu_EB->GetXaxis()->FindBin(myLepPt));
    else
      return h1D_FRmu_EE->GetBinContent(h1D_FRmu_EE->GetXaxis()->FindBin(myLepPt));
  }else{
    cout<<"ERROR! wrong lepton ID : "<<myLepID<<endl;
    return 0.;
  }  
}

void ReducibleBackgroundAA_2015()
{

  TFile* fFakeRates = TFile::Open("fakeRates_20151202.root");
  /*
  h1D_FRmu_EB = (TH1F*)fFakeRates->Get("h1D_FRmu_EB");
  h1D_FRmu_EE = (TH1F*)fFakeRates->Get("h1D_FRmu_EE");
  h1D_FRel_EB = (TH1F*)fFakeRates->Get("h1D_FRel_EB");
  h1D_FRel_EE = (TH1F*)fFakeRates->Get("h1D_FRel_EE");
  //*/
  //*
  h1D_FRmu_EB = (TH1F*)fFakeRates->Get("NoWZ_h1D_FRmu_EB");
  h1D_FRmu_EE = (TH1F*)fFakeRates->Get("NoWZ_h1D_FRmu_EE");
  h1D_FRel_EB = (TH1F*)fFakeRates->Get("NoWZ_h1D_FRel_EB");
  h1D_FRel_EE = (TH1F*)fFakeRates->Get("NoWZ_h1D_FRel_EE");
  //*/

  TFile* dataFile = TFile::Open("../DataTrees_151202/ZZ4lAnalysis_allData.root") ;
  TTree* mytree = (TTree*)dataFile->Get("CRZLLTree/candTree") ;


  Int_t nRun;
  Long64_t nEvent;
  Int_t nLumi;
  Int_t CRflag;
  vector<Float_t> *LepPt = 0;
  vector<Float_t> *LepEta = 0;
  vector<Short_t> *LepLepId = 0;
  Short_t Z1Flav;
  Short_t Z2Flav;
  mytree->SetBranchAddress("RunNumber", &nRun);
  mytree->SetBranchAddress("EventNumber", &nEvent);
  mytree->SetBranchAddress("LumiNumber", &nLumi);
  mytree->SetBranchAddress("CRflag", &CRflag);
  mytree->SetBranchAddress("LepPt", &LepPt);
  mytree->SetBranchAddress("LepEta", &LepEta);
  mytree->SetBranchAddress("LepLepId", &LepLepId);
  mytree->SetBranchAddress("Z1Flav", &Z1Flav);
  mytree->SetBranchAddress("Z2Flav", &Z2Flav);

  Int_t currentFinalState;

  Float_t expectedNumberOfEvents[nFinalStates];
  Float_t NumberOfEvents[nFinalStates];
  for(int fs=0; fs<nFinalStates; fs++){
    expectedNumberOfEvents[fs] = 0.;
    NumberOfEvents[fs] = 0.;
  }

  for(Long64_t iEvt=0; iEvt<mytree->GetEntries(); ++iEvt){

    mytree->GetEntry(iEvt);

    if(!CRflag) continue;
    if(!test_bit(CRflag,CRZLLss)) continue;

    currentFinalState = -1;
    if(Z1Flav==-121){
      if(Z2Flav==+121)
	currentFinalState = fs4e;
      else if(Z2Flav==+169)
	currentFinalState = fs2e2mu;
      else
	cerr<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", Z2Flav="<<Z2Flav<<endl;
    }else if(Z1Flav==-169){
      if(Z2Flav==+121)
	currentFinalState = fs2mu2e;
      else if(Z2Flav==+169)
	currentFinalState = fs4mu;
      else
	cerr<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", Z2Flav="<<Z2Flav<<endl;
    }else{
      cerr<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", Z1Flav="<<Z1Flav<<endl;
    }

    expectedNumberOfEvents[currentFinalState] += fsROSSS[currentFinalState] * fakeRate13TeV(LepPt->at(2),LepEta->at(2),LepLepId->at(2)) * fakeRate13TeV(LepPt->at(3),LepEta->at(3),LepLepId->at(3));
    NumberOfEvents[currentFinalState]++;
		
  }

  for(int fs=0; fs<nFinalStates; fs++){
    cout<<fsLabelForSS[fs]<<" : "
	<<expectedNumberOfEvents[fs]
	<<" +/- "<<expectedNumberOfEvents[fs]/sqrt(NumberOfEvents[fs])<<" (stat., evt: "<<NumberOfEvents[fs]<<")" 
	<<" +/- "<<expectedNumberOfEvents[fs]*0.50<< " (syst.)" 
	<<endl ;
  }

  cout<<"Total: "<<expectedNumberOfEvents[0]+expectedNumberOfEvents[1]+expectedNumberOfEvents[2]+expectedNumberOfEvents[3]<<endl;
	
}


