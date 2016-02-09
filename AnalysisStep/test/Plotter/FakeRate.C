///
/// Use the interpreter and not the compiler 
/// root -l
/// .L FakeRate.C
/// FakeRate() if in lxcms03 or
/// FakeRate(inputFilePath_Data,inputFilePath_MC)
///

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Math/DistFunc.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "plotUtils.C"

#include <ZZAnalysis/AnalysisStep/src/Category.cc>
#include <ZZAnalysis/AnalysisStep/src/bitops.cc>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>

using namespace std;


void FakeRate(Float_t lumi,  string inputFilePath_Data = "/data3/Higgs/160203/", 
		string inputFilePath_MC = "/data3/Higgs/160203/"){

  const int nDatasets = 20;

  string datasets[nDatasets] = {
    //    "AllData", //Change name folder to use this. need 2015
    "DoubleMu2015B",
    "DoubleEG2015B",
    "MuonEG2015B",
    "SingleEle2015B",
    "DoubleMu2015C",
    "DoubleEG2015C",
    "MuonEG2015C",
    "SingleEle2015C",
    "DoubleMu2015C_50ns",
    "DoubleEG2015C_50ns",
    "MuonEG2015C_50ns",
    "SingleEle2015C_50ns",
    "DoubleMu2015D",
    "DoubleEG2015D",
    "MuonEG2015D",
    "SingleEle2015D",
    "WZTo3LNu",
    "DYJetsToLL_M50",
    "TTTo2L2Nu",
    "TTJets"
  };
    
  string finalState[2] = {
    "el",
    "mu",
  };
  
  
  string variables[8] = {
    "LepCombRelIsoPF",
    "Z1Mass",
    "LepSIP",
    "LepPt",
    "FakeRate_num_barrel_pt",        
    "FakeRate_denom_barrel_pt",       
    "FakeRate_num_endcap_pt",        
    "FakeRate_denom_endcap_pt"         
  };
  
  enum var {
    LepCombRelIsoPF_,
    Z1Mass_,
    LepSIP_,    
    LepPt_,
    LepPtPassBar_,
    LepPtTotBar_,
    LepPtPassEnd_,
    LepPtTotEnd_,    	 
  };
  
  TTree* inputTree[nDatasets];
  TFile* inputFile[nDatasets];

  TH1F* hCounters[nDatasets];

  Long64_t NGenEvt[nDatasets];
  Double_t gen_sumWeights[nDatasets];
  Float_t partialSampleWeight[nDatasets];

  
  Long64_t nEvent;
  Float_t Z1Mass;
  Float_t PFMET;
  Float_t PFMETNoHF;
  Float_t overallEventWeight;
  Float_t xsec;
  
  vector<float> *LepPt=0;
  vector<float> *LepEta=0;
  vector<short> *LepLepId=0;
  vector<float> *LepSIP=0;
  vector<float> *LepCombRelIsoPF = 0; 
  vector<bool>  *LepisID = 0;

  bool isData[nDatasets];  
  TH1F* h1[nDatasets][3][5];

  gSystem->mkdir("fakeRateResults");
  string bashCommand = "hadd -f fakeRateResults/data.root "  ;

  for(int d=0; d<nDatasets; d++){
    
    cout<<"Processing dataset "<<datasets[d]<<endl;
    
    if(datasets[d].find("2015")<100) {
      isData[d]=kTRUE;
      string addData = "fakeRateResults/" +datasets[d]+".root ";
      bashCommand += addData;
    }





    string inputFileName = string(Form("%s%s/ZZ4lAnalysis.root",(isData[d]?inputFilePath_Data:inputFilePath_MC).c_str(),datasets[d].c_str()));
    
    inputFile[d] = new TFile(inputFileName.c_str());
    

    if(inputFile[d]->IsZombie()){
      std::cout<<"file "<<inputFileName<<" not found or zombie"<<std::endl;  
      continue;
    }


    hCounters[d] = (TH1F*)inputFile[d]->Get("CRZLTree/Counters");
    NGenEvt[d] = (Long64_t)hCounters[d]->GetBinContent(1);
    gen_sumWeights[d] = (Long64_t)hCounters[d]->GetBinContent(40);

    inputTree[d] = (TTree*)inputFile[d]->Get("CRZLTree/candTree");
    inputTree[d]->SetBranchAddress("overallEventWeight", &overallEventWeight);
    inputTree[d]->SetBranchAddress("EventNumber", &nEvent);
    inputTree[d]->SetBranchAddress("xsec", &xsec);
    inputTree[d]->SetBranchAddress("Z1Mass", &Z1Mass);
    inputTree[d]->SetBranchAddress("LepPt",&LepPt);
    inputTree[d]->SetBranchAddress("LepEta",&LepEta);
    inputTree[d]->SetBranchAddress("LepLepId",&LepLepId);
    inputTree[d]->SetBranchAddress("LepSIP",&LepSIP);
    inputTree[d]->SetBranchAddress("LepCombRelIsoPF",&LepCombRelIsoPF);
    inputTree[d]->SetBranchAddress("PFMET",&PFMET);
    inputTree[d]->SetBranchAddress("PFMETNoHF",&PFMETNoHF);
    inputTree[d]->SetBranchAddress("LepisID",&LepisID);


    partialSampleWeight[d] = lumi * 1000 / gen_sumWeights[d] ;

    Float_t xbins[8] = {5,7,10,20,30,40,50,80};
    //    Float_t xbins[4] = {0,30,60,200};   

    Long64_t entries = inputTree[d]->GetEntries();
    
    for(int v=0; v<=7; v++){
      for(int f=0;f<=1;f++){
	if(LepPtPassBar_    || v==var:: LepPtTotBar_  || v==var::  LepPtPassEnd_ || v==var::  LepPtTotEnd_ ) h1[d][f][v] = new TH1F(("h_"+variables[v]+"_"+finalState[f]).c_str(),";;",7,xbins);
	else if(v==var::LepCombRelIsoPF_)   h1[d][f][v] = new TH1F(("h_"+variables[v]+"_"+finalState[f]).c_str(),";;",10,0.,2.);
	else  h1[d][f][v] = new TH1F(("h_"+variables[v]+"_"+finalState[f]).c_str(),";;",200,0.,100.);
      }
    }
    



    for (Long64_t z=0; z<entries; ++z){
      
      inputTree[d]->GetEntry(z);   
      

      
      bool isTrigger =  ( ((LepPt->at(0)>20.) && (LepPt->at(1)>10.) ) || ( (LepPt->at(1)>20.) && (LepPt->at(0)>10.) ) );

      if(  (abs(Z1Mass-91.19)>7) || (LepCombRelIsoPF->at(0)>0.35) || (LepCombRelIsoPF->at(1)>0.35) || (LepSIP->at(0)>4.) || (LepSIP->at(1)>4.) || (LepSIP->at(2)>4.) ||  !isTrigger  || PFMETNoHF>25. ) continue;

      Double_t eventWeight = partialSampleWeight[d] * xsec ;// * overallEventWeight ; //FIXME no weight for CRZL 
    

      bool FullSell  = ((LepCombRelIsoPF->at(2)<0.35) && (LepisID->at(2)));
      
      if(abs(LepLepId->at(2))==11){ 
	h1[d][0][var::LepPt_]->Fill(LepPt->at(2),(isData[d]==kTRUE)?1.:eventWeight);
	h1[d][0][var::LepCombRelIsoPF_]->Fill(LepCombRelIsoPF->at(1),(isData[d]==kTRUE)?1.:eventWeight);
	
	if(abs(LepEta->at(2))<1.45){

	  h1[d][0][var::LepPtTotBar_]->Fill(LepPt->at(2),(isData[d]==kTRUE)?1.:eventWeight);

	  if(FullSell) h1[d][0][var::LepPtPassBar_]->Fill(LepPt->at(2),(isData[d]==kTRUE)?1.:eventWeight);
	} 
	else{

	  h1[d][0][var::LepPtTotEnd_]->Fill(LepPt->at(2),(isData[d]==kTRUE)?1.:eventWeight);
	  if(FullSell) h1[d][0][var::LepPtPassEnd_]->Fill(LepPt->at(2),(isData[d]==kTRUE)?1.:eventWeight);
	}
      }

      else  if(abs(LepLepId->at(2))==13){
	h1[d][1][var::LepPt_]->Fill(LepPt->at(2),(isData[d]==kTRUE)?1.:eventWeight);
	h1[d][1][var::LepCombRelIsoPF_]->Fill(LepCombRelIsoPF->at(1),(isData[d]==kTRUE)?1.:eventWeight);
	
	if(abs(LepEta->at(2))<1.45){

	  h1[d][1][var::LepPtTotBar_]->Fill(LepPt->at(2),(isData[d]==kTRUE)?1.:eventWeight);
	  if(FullSell) h1[d][1][var::LepPtPassBar_]->Fill(LepPt->at(2),(isData[d]==kTRUE)?1.:eventWeight);
	} 
	else{
	  h1[d][1][var::LepPtTotEnd_]->Fill(LepPt->at(2),(isData[d]==kTRUE)?1.:eventWeight);
	  if(FullSell) h1[d][1][var::LepPtPassEnd_]->Fill(LepPt->at(2),(isData[d]==kTRUE)?1.:eventWeight);
	}
      }
    }


    TFile* fOutHistos = new TFile(string("fakeRateResults/"+datasets[d]+".root").c_str(),"recreate");
       
    for(int v=0;v<=7;v++){
      for(int f=0;f<=1;f++){
	h1[d][f][v]->Write();
      }
    }
    fOutHistos->Close();

  } //End Samples loop

  std::cout<<bashCommand<<"\n";
  gROOT->ProcessLine(Form(".!%s",bashCommand.c_str()));
}

