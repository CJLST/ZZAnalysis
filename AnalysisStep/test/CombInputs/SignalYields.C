/* 
 * usage: 
 * -set all input variables in Config.h
 * -run with:
 * root -q -l -b SignalYields.C++
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TMath.h"

#include <ZZAnalysis/AnalysisStep/src/Category.cc>
#include "Config.h"

using namespace std;


enum Process {ggH=0, qqH=1, WH=2, ZH=3, ttH=4};
const int nProcesses = 5;
string sProcess[nProcesses] = {"ggH", "qqH", "WH", "ZH", "ttH"};

enum FinalState {fs4mu=0, fs4e=1, fs2e2mu=2, fs2mu2e=3};
const int nFinalStates = 4;
string sFinalState[nFinalStates] = {"4mu", "4e", "2e2mu", "2mu2e"};

const int nCategories = 6;
string sCategory[nCategories] = {
  "Untagged",
  "OneJetTagged",
  "VBFTagged", 
  "VHLeptTagged",
  "VHHadrTagged",
  "ttHTagged",
};

enum ResonantStatus {resonant=0, nonresonant=1};
const int nResStatuses = 2;
TString sResonantStatus[nResStatuses] = {"resonant", "nonresonant"};

Double_t deltaR(Double_t e1, Double_t p1, Double_t e2, Double_t p2) {
  Double_t deltaPhi = acos(cos(p1-p2));
  return TMath::Sqrt((e1-e2)*(e1-e2) + deltaPhi*deltaPhi);
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void computeSignalYields(string outputDirectory, double lumi, double sqrts, double mHPoint, double m4lMin, double m4lMax)
{

  const int nDatasets = 6;
  string datasets[nDatasets] = {
    "ggH",
    "VBFH",
    "WplusH",
    "WminusH",
    "ZH",
    "ttH",
  };

  TFile* inputFile[nDatasets];
  TTree* inputTree[nDatasets];
  TH1F* hCounters[nDatasets];
  Long64_t NGenEvt[nDatasets];
  Double_t gen_sumWeights[nDatasets];
  Float_t partialSampleWeight[nDatasets];

  Int_t nRun;
  Long64_t nEvent;
  Int_t nLumi;
  Short_t Nvtx;
  Short_t NObsInt;
  Float_t NTrueInt;
  Float_t overallEventWeight;
  Float_t xsec;
  Short_t ZZsel;
  Float_t ZZMass;
  Float_t ZZPt;
  Short_t Z1Flav;
  Short_t Z2Flav;
  Float_t DiJetFisher;
  vector<Float_t> *CandLepEta = 0;
  vector<Float_t> *CandLepPhi = 0;
  Short_t nExtraLep;
  vector<Float_t> *ExtraLepEta = 0;
  vector<Float_t> *ExtraLepPhi = 0;
  Short_t nJets;
  Short_t nJetsBTagged;
  vector<Float_t> *JetPt = 0;
  vector<Float_t> *JetEta = 0;
  vector<Float_t> *JetPhi = 0;
  vector<Float_t> *JetMass = 0;
  Float_t jetPt[99];
  Float_t jetEta[99];
  Float_t jetPhi[99];
  Float_t jetMass[99];
  Float_t GenLep1Eta;
  Float_t GenLep1Phi;
  Short_t GenLep1Id;
  Float_t GenLep2Eta;
  Float_t GenLep2Phi;
  Short_t GenLep2Id;
  Float_t GenLep3Eta;
  Float_t GenLep3Phi;
  Short_t GenLep3Id;
  Float_t GenLep4Eta;
  Float_t GenLep4Phi;
  Short_t GenLep4Id;

  Float_t yield[nProcesses][nFinalStates][nCategories+1][nResStatuses+1];
  for(int pr=0; pr<nProcesses; pr++)
    for(int fs=0; fs<nFinalStates; fs++)
      for(int cat=0; cat<nCategories+1; cat++)
	for(int rs=0; rs<nResStatuses+1; rs++)
	  yield[pr][fs][cat][rs] = 0.;
  
  int currentProcess;
  int currentFinalState;
  int currentCategory;
  int currentResStatus;


  //---------- Will loop over all datasets

  for(int d=0; d<nDatasets; d++){

    if(datasets[d]!="ggH" && mHPoint!=125.) continue; //FIXME

    cout<<"Processing dataset "<<datasets[d]<<mHPoint<<"..."<<endl;

    string inputFileName = string(Form("%s%s%i/ZZ4lAnalysis.root",inputFilePath.c_str(),datasets[d].c_str(),(int)mHPoint));
    inputFile[d] = TFile::Open(inputFileName.c_str());

    hCounters[d] = (TH1F*)inputFile[d]->Get("ZZTree/Counters");
    NGenEvt[d] = (Long64_t)hCounters[d]->GetBinContent(1);
    gen_sumWeights[d] = (Long64_t)hCounters[d]->GetBinContent(40);
    partialSampleWeight[d] = lumi * 1000 / gen_sumWeights[d] ;

    inputTree[d] = (TTree*)inputFile[d]->Get("ZZTree/candTree");
    inputTree[d]->SetBranchAddress("RunNumber", &nRun);
    inputTree[d]->SetBranchAddress("EventNumber", &nEvent);
    inputTree[d]->SetBranchAddress("LumiNumber", &nLumi);
    inputTree[d]->SetBranchAddress("Nvtx", &Nvtx);
    inputTree[d]->SetBranchAddress("NObsInt", &NObsInt);
    inputTree[d]->SetBranchAddress("NTrueInt", &NTrueInt);
    inputTree[d]->SetBranchAddress("overallEventWeight", &overallEventWeight);
    inputTree[d]->SetBranchAddress("xsec", &xsec);
    inputTree[d]->SetBranchAddress("ZZsel", &ZZsel);
    inputTree[d]->SetBranchAddress("ZZMass", &ZZMass);
    inputTree[d]->SetBranchAddress("ZZPt", &ZZPt);
    inputTree[d]->SetBranchAddress("Z1Flav", &Z1Flav);
    inputTree[d]->SetBranchAddress("Z2Flav", &Z2Flav);
    inputTree[d]->SetBranchAddress("LepEta", &CandLepEta);
    inputTree[d]->SetBranchAddress("LepPhi", &CandLepPhi);
    inputTree[d]->SetBranchAddress("nExtraLep", &nExtraLep);
    inputTree[d]->SetBranchAddress("ExtraLepEta", &ExtraLepEta);
    inputTree[d]->SetBranchAddress("ExtraLepPhi", &ExtraLepPhi);
    inputTree[d]->SetBranchAddress("nCleanedJetsPt30", &nJets);
    inputTree[d]->SetBranchAddress("nCleanedJetsPt30BTagged", &nJetsBTagged);
    inputTree[d]->SetBranchAddress("JetPt", &JetPt);
    inputTree[d]->SetBranchAddress("JetEta", &JetEta);
    inputTree[d]->SetBranchAddress("JetPhi", &JetPhi);
    inputTree[d]->SetBranchAddress("JetMass", &JetMass);
    inputTree[d]->SetBranchAddress("DiJetFisher", &DiJetFisher);
    inputTree[d]->SetBranchAddress("GenLep1Eta", &GenLep1Eta);
    inputTree[d]->SetBranchAddress("GenLep1Phi", &GenLep1Phi);
    inputTree[d]->SetBranchAddress("GenLep1Id", &GenLep1Id);
    inputTree[d]->SetBranchAddress("GenLep2Eta", &GenLep2Eta);
    inputTree[d]->SetBranchAddress("GenLep2Phi", &GenLep2Phi);
    inputTree[d]->SetBranchAddress("GenLep2Id", &GenLep2Id);
    inputTree[d]->SetBranchAddress("GenLep3Eta", &GenLep3Eta);
    inputTree[d]->SetBranchAddress("GenLep3Phi", &GenLep3Phi);
    inputTree[d]->SetBranchAddress("GenLep3Id", &GenLep3Id);
    inputTree[d]->SetBranchAddress("GenLep4Eta", &GenLep4Eta);
    inputTree[d]->SetBranchAddress("GenLep4Phi", &GenLep4Phi);
    inputTree[d]->SetBranchAddress("GenLep4Id", &GenLep4Id);


    //----- assign dataset to correct process
    
    currentProcess = -1;
    if(datasets[d]=="ggH") currentProcess = ggH;
    if(datasets[d]=="VBFH") currentProcess = qqH;
    if(datasets[d]=="WplusH") currentProcess = WH;
    if(datasets[d]=="WminusH") currentProcess = WH;
    if(datasets[d]=="ZH") currentProcess = ZH;
    if(datasets[d]=="ttH") currentProcess = ttH;


    //---------- Process tree

    Long64_t entries = inputTree[d]->GetEntries();

    for (Long64_t z=0; z<entries; ++z){

      inputTree[d]->GetEntry(z);
      
      if( !(ZZsel>=90) ) continue;
      if(ZZMass<m4lMin || ZZMass>m4lMax) continue;

      Double_t eventWeight = partialSampleWeight[d] * xsec * overallEventWeight ;


      //----- find final state

      currentFinalState = -1;
      if(Z1Flav==-121){
	if(Z2Flav==-121)
	  currentFinalState = fs4e;
	else if(Z2Flav==-169)
	  currentFinalState = fs2e2mu;
	else
	  cerr<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", Z2Flav="<<Z2Flav<<endl;
      }else if(Z1Flav==-169){
	if(Z2Flav==-121)
	  currentFinalState = fs2mu2e;
	else if(Z2Flav==-169)
	  currentFinalState = fs4mu;
	else
	  cerr<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", Z2Flav="<<Z2Flav<<endl;
      }else{
	cerr<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", Z1Flav="<<Z1Flav<<endl;
      }


      //----- find category

      for(int j=0; j<nJets; j++){
	jetPt[j] = JetPt->at(j);
	jetEta[j] = JetEta->at(j);
	jetPhi[j] = JetPhi->at(j);
	jetMass[j] = JetMass->at(j);
      }
      currentCategory = category(
	     nExtraLep,
	     ZZPt,
	     ZZMass,
	     nJets, 
	     nJetsBTagged,
	     jetPt,
	     jetEta,
	     jetPhi,
	     jetMass,
	     DiJetFisher
	     );


      //----- check if this is resonant signal (ie. H->4l with correct lepton choice)

      Short_t GenHLepId[4] = {GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id};
      Float_t GenHLepEta[4] = {GenLep1Eta,GenLep2Eta,GenLep3Eta,GenLep4Eta};
      Float_t GenHLepPhi[4] = {GenLep1Phi,GenLep2Phi,GenLep3Phi,GenLep4Phi};

      Int_t nGenHLep = 0;
      Int_t nRecoLepMatchedToGenHLep[4] = {0,0,0,0};
      Int_t nGenHLepMatchedToCandLep[4] = {0,0,0,0};
      for(Int_t iGenHLep=0; iGenHLep<4; iGenHLep++){
	if(abs(GenHLepId[iGenHLep])==11 || abs(GenHLepId[iGenHLep])==13){
	  nGenHLep++;
	  for(Int_t iCandLep=0; iCandLep<4; iCandLep++){
	    if(deltaR(GenHLepEta[iGenHLep],GenHLepPhi[iGenHLep],CandLepEta->at(iCandLep),CandLepPhi->at(iCandLep)) < 0.1){
	      nRecoLepMatchedToGenHLep[iGenHLep]++;
	      nGenHLepMatchedToCandLep[iCandLep]++;
	    }
	  }
	  for(Int_t iExtraLep=0; iExtraLep<nExtraLep; iExtraLep++){
	    if(deltaR(GenHLepEta[iGenHLep],GenHLepPhi[iGenHLep],ExtraLepEta->at(iExtraLep),ExtraLepPhi->at(iExtraLep)) < 0.1){
	      nRecoLepMatchedToGenHLep[iGenHLep]++;
	    }
	  }
	}
      }
      Bool_t foundMatchingAmbiguity = false;
      for(Int_t iGenHLep=0; iGenHLep<4; iGenHLep++) if(nRecoLepMatchedToGenHLep[iGenHLep]>1){ foundMatchingAmbiguity = true; break; }
      for(Int_t iCandLep=0; iCandLep<4; iCandLep++) if(nGenHLepMatchedToCandLep[iCandLep]>1){ foundMatchingAmbiguity = true; break; }
      Int_t nMatches = 0;
      for(Int_t iGenHLep=0; iGenHLep<4; iGenHLep++) if(nRecoLepMatchedToGenHLep[iGenHLep]==1) nMatches++;

      if(nGenHLep==4 && !foundMatchingAmbiguity && nMatches==4) 
	currentResStatus = resonant;
      else 
	currentResStatus = nonresonant;


      //----- fill counter

      yield[currentProcess][currentFinalState][currentCategory][currentResStatus] += eventWeight;


    } // end for entries

  } // end for datasets


  //---------- Fill 'inclusive' counters
  for(int pr=0; pr<nProcesses; pr++)
    for(int fs=0; fs<nFinalStates; fs++)
      for(int cat=0; cat<nCategories; cat++)
	for(int rs=0; rs<nResStatuses; rs++)
	  yield[pr][fs][nCategories][rs] += yield[pr][fs][cat][rs];
  for(int pr=0; pr<nProcesses; pr++)
    for(int fs=0; fs<nFinalStates; fs++)
      for(int cat=0; cat<nCategories+1; cat++)
	for(int rs=0; rs<nResStatuses; rs++)
	  yield[pr][fs][cat][nResStatuses] += yield[pr][fs][cat][rs];


  //---------- Write yields to txt file

  string outputFileName = string(Form("%s/signalYields_M%i_%iTeV_window%i-%i.txt",outputDirectory.c_str(),(int)mHPoint,(int)sqrts,(int)m4lMin,(int)m4lMax));
  ofstream outFile;
  outFile.open(outputFileName);

  outFile<<"############### SIGNAL YIELDS ###############"<<endl;
  outFile<<"mass window: "<<m4lMin<<" <= m4l <= "<<m4lMax<<endl;
  outFile<<"sqrt(s) = "<<sqrts<<" TeV"<<endl;
  outFile<<"integrated lumi. = "<<lumi<<" fb-1"<<endl;
  outFile<<endl;
  
  for(int pr=0; pr<nProcesses; pr++){
    if(pr>0 && mHPoint!=125.) continue; //FIXME
    for(int fs=0; fs<nFinalStates; fs++){
      for(int cat=0; cat<nCategories; cat++){
	outFile<<sProcess[pr]<<" "
	       <<sFinalState[fs]<<" "
	       <<sCategory[cat]<<" "
	       <<yield[pr][fs][cat][nResStatuses]
	       <<endl;
      }
    }
    outFile<<endl;
  }

  outFile.close();

}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void SignalYields() {

  string outDir = "YieldFiles";
  gSystem->Exec(("mkdir -p "+outDir).c_str());

  for(int mp=0; mp<nMassPoints; mp++){
    computeSignalYields(outDir, lumi13TeV, 13., massPoints[mp], 105., 140.);
  }

}

