/* 
 * usage: 
 * -set all input variables in Config.h
 * -run with:
 * root -q -l -b IrrBkgdYields.C++
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


enum Process {qqZZ=0, ggZZ=1};
const int nProcesses = 2;
string sProcess[nProcesses] = {"qqZZ", "ggZZ"};

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



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void computeIrrBkgdYields(string outputDirectory, double lumi, double sqrts, double m4lMin, double m4lMax)
{

  const int nDatasets = 7;
  string datasets[nDatasets] = {
    "ZZTo4l",
    "ggZZ4e",
    "ggZZ4mu",
    "ggZZ4tau",
    "ggZZ2e2mu",
    "ggZZ2e2tau",
    "ggZZ2mu2tau",
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
  Short_t nExtraLep;
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

  Float_t yield[nProcesses][nFinalStates][nCategories+1];
  for(int pr=0; pr<nProcesses; pr++)
    for(int fs=0; fs<nFinalStates; fs++)
      for(int cat=0; cat<nCategories+1; cat++)
	yield[pr][fs][cat] = 0.;
  
  int currentProcess;
  int currentFinalState;
  int currentCategory;


  //---------- Will loop over all datasets

  for(int d=0; d<nDatasets; d++){

    cout<<"Processing dataset "<<datasets[d]<<"..."<<endl;

    string inputFileName = string(Form("%s%s/ZZ4lAnalysis.root",inputFilePath.c_str(),datasets[d].c_str()));
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
    inputTree[d]->SetBranchAddress("nExtraLep", &nExtraLep);
    inputTree[d]->SetBranchAddress("nCleanedJetsPt30", &nJets);
    inputTree[d]->SetBranchAddress("nCleanedJetsPt30BTagged", &nJetsBTagged);
    inputTree[d]->SetBranchAddress("JetPt", &JetPt);
    inputTree[d]->SetBranchAddress("JetEta", &JetEta);
    inputTree[d]->SetBranchAddress("JetPhi", &JetPhi);
    inputTree[d]->SetBranchAddress("JetMass", &JetMass);
    inputTree[d]->SetBranchAddress("DiJetFisher", &DiJetFisher);


    //----- assign dataset to correct process
    
    currentProcess = -1;
    if(datasets[d]=="ZZTo4l") currentProcess = qqZZ;
    if(datasets[d]=="ggZZ4e"||
       datasets[d]=="ggZZ4mu"||
       datasets[d]=="ggZZ4tau"||
       datasets[d]=="ggZZ2e2mu"||
       datasets[d]=="ggZZ2e2tau"||
       datasets[d]=="ggZZ2mu2tau") 
      currentProcess = ggZZ;


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


      //----- fill counter

      yield[currentProcess][currentFinalState][currentCategory] += eventWeight;
      yield[currentProcess][currentFinalState][nCategories] += eventWeight;


    } // end for entries

  } // end for datasets


  //---------- Write yields to txt file

  string outputFileName = string(Form("%s/irrBkgdYields_%iTeV_window%i-%i.txt",outputDirectory.c_str(),(int)sqrts,(int)m4lMin,(int)m4lMax));
  ofstream outFile;
  outFile.open(outputFileName);

  outFile<<"############### IRREDUCIBLE BACKGROUND YIELDS ###############"<<endl;
  outFile<<"mass window: "<<m4lMin<<" <= m4l <= "<<m4lMax<<endl;
  outFile<<"sqrt(s) = "<<sqrts<<" TeV"<<endl;
  outFile<<"integrated lumi. = "<<lumi<<" fb-1"<<endl;
  outFile<<endl;
  
  for(int pr=0; pr<nProcesses; pr++){
    for(int fs=0; fs<nFinalStates; fs++){
      for(int cat=0; cat<nCategories; cat++){
	outFile<<sProcess[pr]<<" "
	       <<sFinalState[fs]<<" "
	       <<sCategory[cat]<<" "
	       <<yield[pr][fs][cat]
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


void IrrBkgdYields() {

  string outDir = "YieldFiles";
  gSystem->Exec(("mkdir -p "+outDir).c_str());

  computeIrrBkgdYields(outDir, lumi13TeV, 13., 105.,  140.);
  computeIrrBkgdYields(outDir, lumi13TeV, 13.,  70., 1000.);

}

