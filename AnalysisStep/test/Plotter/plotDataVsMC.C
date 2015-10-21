/* 
 * usage: 
 * -specify 'inputPath' at the end of this file
 * -run with:
 *   root -l plotDataVsMC.C++
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TMath.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "Math/DistFunc.h"

#include "tdrstyle.C"
#include "CMS_lumi.C"

#include <ZZAnalysis/AnalysisStep/src/Category.cc>

using namespace std;

#define DEBUG 0
#define DRAWLINES 1
#define DRAWLABELBYHAND 1
#define APPLYKFACTORS 1
#define REBINDYTTBAR 0





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

enum Blindings {fullyblind=0, blindabove110=1, blindbelow150=2, blind110150=3, unblinded=4};
const int nBlindings = 5;
string sBlinding[nBlindings] = {"fullyblind", "blindabove110", "blindbelow150", "blind110150", "unblinded"};
Bool_t plotThisBlinding[5][nBlindings] = {
  {1,1,1,1,0,},
  {1,0,0,0,0,},
  {0,1,0,0,0,},
  {0,0,1,0,0,},
  {0,0,0,1,0,},
};
string blindingLabel[nBlindings] = {"", "m_{4#font[12]{l}} < 110 GeV", "m_{4#font[12]{l}} > 150 GeV", "m_{4#font[12]{l}} #notin [110, 150] GeV", ""};
Float_t xHistoBlindLow[nBlindings] = {  0.,  110.,   0., 110.,  0. };
Float_t xHistoBlindUp [nBlindings] = { -1., 1000., 150., 150., -1. };

const int nVariables = 12;
string varName[nVariables] = {
  "M4l",
  "M4l2",
  "M4lLow",
  "M4lHigh",
  "MZ1",
  "MZ2",
  "KD",
  "Djet",
  "Pt4l",
  "NExtraLep",
  "NJets",
  "NJetsBTagged",
};
string varXLabel[nVariables] = {
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{4#font[12]{l}} (GeV)",
  "m_{Z_{1}} (GeV)",
  "m_{Z_{2}} (GeV)",
  "D_{bkg}^{kin}",
  "D_{jet}",
  "p_{T}^{4l} (GeV)",
  "number of additional leptons",
  "number of jets",
  "number of b-tagged jets",
};
string varYLabel[nVariables] = {
  "Events / 3 GeV",
  "Events / 6 GeV",//"Events / 12 GeV",//
  "Events / 4 GeV",
  "Events / 20 GeV",
  "Events / 2 GeV",//"Events / 5 GeV",//
  "Events / 2 GeV",//"Events / 5 GeV",//
  "Events / 0.05",
  "Events / 0.1",
  "Events / 10 GeV",
  "Events",
  "Events",
  "Events",
};
Bool_t plotThisVar[3][nVariables] = {
  {1,0,0,0,1,1,1,0,0,0,1,0,},
  {1,1,1,1,1,1,1,0,0,0,1,0,},
  {1,1,1,1,1,1,1,1,1,1,1,1,},
};
//*
Int_t  varNbin[nVariables] = { 272, 136,  10,  30,  75,  75, 20, 20,  40, 6, 17, 8, };
Float_t varMin[nVariables] = {  70,  70,  70, 150,   0,   0,  0,  0,   0, 0,  0, 0, };
Float_t varMax[nVariables] = { 886, 886, 110, 750, 150, 150,  1,  2, 400, 6, 17, 8, };
//*/
/*
Int_t  varNbin[nVariables] = { 272,  72,  10,  30,  30,  30, 20, 20,  40, 6, 17, 8, };
Float_t varMin[nVariables] = {  70,  70,  70, 150,   0,   0,  0,  0,   0, 0,  0, 0, };
Float_t varMax[nVariables] = { 886, 886, 110, 750, 150, 150,  1,  2, 400, 6, 17, 8, };
//*/
Bool_t varLogx[nVariables] = {1,1,0,0,0,0,0,0,0,0,0,0,};
Bool_t varLogy[nVariables] = {0,0,0,0,0,0,0,0,0,1,1,1,};
Int_t restrictCountVar[nVariables] = {0,0,0,0,0,0,0,0,0,2,4,2,};
Float_t varMinFactor[nVariables] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,1000000.,5000.,100000.,};
Int_t varCMSPos[nVariables] = {33,33,11,33,11,11,11,33,33,33,33,33};
Int_t varLegPos[nVariables] = {3,3,1,3,1,1,3,3,3,3,3,3,};
Int_t rebinning[nVariables] = {16,8,2,3,2,2,4,4,2,1,1,1,};

Float_t varMaxCorrector[nBlindings][nVariables] = {
  { 1. , 1., 1., 1., 1., 1. , 1.3, 1., 1., 1.,  1., 1., },
  { 1. , 1., 1., 1., 1., 1.1, 1.1, 1., 1., 1.,  1., 1., },
  { 1.2, 1., 1., 1., 1., 1. , 1.1, 1., 1., 1., 10., 1., },
  { 1. , 1., 1., 1., 1., 1. , 1.1, 1., 1., 1., 10., 1., },
  { 1. , 1., 1., 1., 1., 1. , 1.1, 1., 1., 1.,  5., 1., },
};



enum Process {Data=0, H125=1, qqZZ=2, ggZZ=3, DY=4, ttbar=5};
const int nProcesses = 6;
string sProcess[nProcesses] = {"Data", "H125", "qqZZ", "ggZZ", "DY", "ttbar"};
string processLabel[nProcesses] = {"Data", "m_{H} = 125 GeV", "q#bar{q}#rightarrowZZ", "gg#rightarrowZZ", "Z + jets", "t#bar{t}"};
Int_t processFillColor[nProcesses] = {TColor::GetColor("#000000"), TColor::GetColor("#ffafaf"), TColor::GetColor("#99ccff"), TColor::GetColor("#3366ff"), TColor::GetColor("#669966"), TColor::GetColor("#9a6666")};
Int_t processLineColor[nProcesses] = {TColor::GetColor("#000000"), TColor::GetColor("#cc0000"), TColor::GetColor("#000099"), TColor::GetColor("#000099"), TColor::GetColor("#003300"), TColor::GetColor("#5f3f3f")};
Float_t kFactor[nProcesses] = {1.,1.,1.065,2.,1.,1.,}; //FIXME: these may change in the future
Bool_t useProcess[nProcesses] = {1,1,1,1,0,0,};

enum FinalState {fs4mu=0, fs4e=1, fs2e2mu=2, fs2mu2e=3};
const int nFinalStates = 4;
string sFinalState[nFinalStates+1] = {"4mu", "4e", "2e2mu", "2mu2e", "4l"};

const int nCategories = 6;
string sCategory[nCategories+1] = {
  "Untagged",
  "OneJetTagged",
  "VBFTagged", 
  "VHLeptTagged",
  "VHHadrTagged",
  "ttHTagged",
  "inclusive",
};

enum ResonantStatus {resonant=0, nonresonant=1};
const int nResStatuses = 2;
string sResonantStatus[nResStatuses+1] = {"resonant", "nonresonant", "allres"};

Double_t deltaR(Double_t e1, Double_t p1, Double_t e2, Double_t p2) {
  Double_t deltaPhi = acos(cos(p1-p2));
  return TMath::Sqrt((e1-e2)*(e1-e2) + deltaPhi*deltaPhi);
}





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void doHistograms(string inputFilePath, double lumi)
{

  const int nDatasets = 27;
  string datasets[nDatasets] = {
    "DoubleMu2015B",
    "DoubleEG2015B",
    "MuonEG2015B",
    "DoubleMu2015C_50ns",
    "DoubleEG2015C_50ns",
    "MuonEG2015C_50ns",
    "DoubleMu2015D",
    "DoubleEG2015D",
    "MuonEG2015D",
    "DoubleMu2015Dv4",
    "DoubleEG2015Dv4",
    "MuonEG2015Dv4",
    "ggH125",
    "VBFH125",
    "WplusH125",
    "WminusH125",
    "ZH125",
    "ttH125",
    "ZZTo4l",
    "ggZZ4e",
    "ggZZ4mu",
    "ggZZ4tau",
    "ggZZ2e2mu",
    "ggZZ2e2tau",
    "ggZZ2mu2tau",
    "DYJetsToLL_M50",
    "TTTo2L2nu",//"TTJets",
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
  Float_t p0plus_VAJHU;
  Float_t bkg_VAMCFM;
  Float_t Z1Mass;
  Float_t Z2Mass;
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

  TH1F* h1[nVariables][nBlindings][nFinalStates+1][nCategories+1][nResStatuses+1][nProcesses];
  for(int v=0; v<nVariables; v++){
    for(int bl=0; bl<nBlindings; bl++){
      for(int fs=0; fs<nFinalStates+1; fs++){
	for(int cat=0; cat<nCategories+1; cat++){
	  for(int rs=0; rs<nResStatuses+1; rs++){
	    for(int pr=0; pr<nProcesses; pr++){
	      h1[v][bl][fs][cat][rs][pr] = new TH1F(
                 Form("h1_%s_%s_%s_%s_%s_%s",varName[v].c_str(),sBlinding[bl].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str(),sResonantStatus[rs].c_str(),sProcess[pr].c_str()),
		 Form(";%s;%s",varXLabel[v].c_str(),varYLabel[v].c_str()),
		 varNbin[v],varMin[v],varMax[v]);
	    }
	  }
	}
      }
    }
  }
  
  int currentProcess;
  int currentFinalState;
  int currentCategory;
  int currentResStatus;


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
    inputTree[d]->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU);
    inputTree[d]->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM);
    inputTree[d]->SetBranchAddress("Z1Mass", &Z1Mass);
    inputTree[d]->SetBranchAddress("Z2Mass", &Z2Mass);
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
    if(datasets[d]=="DoubleMu2015B"||
       datasets[d]=="DoubleEG2015B"||
       datasets[d]=="MuonEG2015B"||
       datasets[d]=="DoubleMu2015C_50ns"||
       datasets[d]=="DoubleEG2015C_50ns"||
       datasets[d]=="MuonEG2015C_50ns"||
       datasets[d]=="DoubleMu2015D"||
       datasets[d]=="DoubleEG2015D"||
       datasets[d]=="MuonEG2015D"||
       datasets[d]=="DoubleMu2015Dv4"||
       datasets[d]=="DoubleEG2015Dv4"||
       datasets[d]=="MuonEG2015Dv4") 
      currentProcess = Data;
    if(datasets[d]=="ggH125") currentProcess = H125;
    if(datasets[d]=="VBFH125") currentProcess = H125;
    if(datasets[d]=="WplusH125") currentProcess = H125;
    if(datasets[d]=="WminusH125") currentProcess = H125;
    if(datasets[d]=="ZH125") currentProcess = H125;
    if(datasets[d]=="ttH125") currentProcess = H125;
    if(datasets[d]=="ZZTo4l") currentProcess = qqZZ;
    if(datasets[d]=="ggZZ4e"||
       datasets[d]=="ggZZ4mu"||
       datasets[d]=="ggZZ4tau"||
       datasets[d]=="ggZZ2e2mu"||
       datasets[d]=="ggZZ2e2tau"||
       datasets[d]=="ggZZ2mu2tau") 
      currentProcess = ggZZ;
    if(datasets[d]=="DYJetsToLL_M50") currentProcess = DY;
    if(datasets[d]=="TTJets"||
       datasets[d]=="TTTo2L2nu")
      currentProcess = ttbar;


    //---------- Process tree

    Long64_t entries = inputTree[d]->GetEntries();

    for (Long64_t z=0; z<entries; ++z){

      if(DEBUG && z>1000) break;

      inputTree[d]->GetEntry(z);
      
      if( !(ZZsel>=90) ) continue;

      Double_t eventWeight = partialSampleWeight[d] * xsec * (APPLYKFACTORS?kFactor[currentProcess]:1.) * overallEventWeight ;


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


      //----- here, define resonant signal as H->4l where l=e,mu (excluding decays to taus and 'wrong signal' from associated production)

      if(currentProcess==Data){
	currentResStatus = resonant;
      }else{	
	Short_t GenHLepId[4] = {GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id};
	Int_t nGenHLep = 0;
	for(Int_t iGenHLep=0; iGenHLep<4; iGenHLep++){
	  if(abs(GenHLepId[iGenHLep])==11 || abs(GenHLepId[iGenHLep])==13)
	    nGenHLep++;
	}
	currentResStatus = (nGenHLep==4) ? resonant : nonresonant ;
      }

      //----- fill histograms

      Float_t varVal[nVariables] = {
	ZZMass,
	ZZMass,
	ZZMass,
	ZZMass,
	Z1Mass,
	Z2Mass,
	p0plus_VAJHU / ( p0plus_VAJHU + bkg_VAMCFM ),
	DiJetFisher,
	ZZPt,
	(Float_t)nExtraLep,
	(Float_t)nJets,
	(Float_t)nJetsBTagged,
      };

      bool fillM4l[nBlindings] = {
	currentProcess!=Data,
	(currentProcess!=Data || ZZMass<110),
	(currentProcess!=Data || ZZMass>150),
	(currentProcess!=Data || ZZMass<110 || ZZMass>150),
	true,
      };
      bool fillOtherThanM4l[nBlindings] = {
	currentProcess!=Data,
	ZZMass<110,
	ZZMass>150,
	ZZMass<110 || ZZMass>150,
	true,
      };

      for(int v=0; v<nVariables; v++){
	for(int bl=0; bl<nBlindings; bl++){
	  if( (varName[v].find("M4l")!=string::npos && fillM4l         [bl]) ||
	      (varName[v].find("M4l")==string::npos && fillOtherThanM4l[bl])    ){
	    h1[v][bl][currentFinalState][currentCategory][currentResStatus][currentProcess]->Fill(varVal[v],(currentProcess==Data)?1.:eventWeight);
	  }
	}
      }

    } // end for entries

  } // end for datasets


  //---------- Fill 'inclusive' counters
  for(int v=0; v<nVariables; v++){
    for(int bl=0; bl<nBlindings; bl++){
      for(int pr=0; pr<nProcesses; pr++){
	for(int fs=0; fs<nFinalStates; fs++){
	  for(int cat=0; cat<nCategories; cat++){
	    for(int rs=0; rs<nResStatuses; rs++){
	      h1[v][bl][nFinalStates][cat][rs][pr]->Add(h1[v][bl][fs][cat][rs][pr]);
	    }
	  }
	}
	for(int fs=0; fs<nFinalStates+1; fs++){
	  for(int cat=0; cat<nCategories; cat++){
	    for(int rs=0; rs<nResStatuses; rs++){
	      h1[v][bl][fs][nCategories][rs][pr]->Add(h1[v][bl][fs][cat][rs][pr]);
	    }
	  }
	  for(int cat=0; cat<nCategories+1; cat++){
	    for(int rs=0; rs<nResStatuses; rs++){
	      h1[v][bl][fs][cat][nResStatuses][pr]->Add(h1[v][bl][fs][cat][rs][pr]);
	    }
	  }
	}
      }
    }
  }


  //---------- Write histograms to a ROOT file
  TFile* fOutHistos = new TFile("histos_plotRunII.root","recreate");
  fOutHistos->cd();
  for(int v=0; v<nVariables; v++){
    for(int bl=0; bl<nBlindings; bl++){
      for(int fs=0; fs<nFinalStates+1; fs++){
	for(int cat=0; cat<nCategories+1; cat++){
	  for(int rs=0; rs<nResStatuses+1; rs++){
	    for(int pr=0; pr<nProcesses; pr++){
	      if( fs==nFinalStates &&
		  cat==nCategories &&
		  rs==nResStatuses )
		h1[v][bl][fs][cat][rs][pr]->Write(h1[v][bl][fs][cat][rs][pr]->GetName());
	      delete h1[v][bl][fs][cat][rs][pr];
	    }
	  }
	}
      }
    }
  }  
  fOutHistos->Close();
  delete fOutHistos; 

}





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void printInfo(string info, Double_t x1, Double_t y1, Double_t x2, Double_t y2){
  TPaveText* pav = new TPaveText(x1,y1,x2,y2,"brNDC");
  pav->SetFillStyle(0);
  pav->SetBorderSize(0);
  pav->SetTextAlign(12);
  pav->AddText(info.c_str());
  pav->Draw();
}

void SaveCanvas(string directory, TCanvas* c, string tag = "") {
  c->SaveAs(Form("%s%s_%s.root",directory.c_str(),c->GetName(),tag.c_str()));
  c->SaveAs(Form("%s%s_%s.pdf" ,directory.c_str(),c->GetName(),tag.c_str()));  
  c->SaveAs(Form("%s%s_%s.png" ,directory.c_str(),c->GetName(),tag.c_str()));
}

TGraphAsymmErrors* getDataGraph(TH1F* h, bool drawZeroBins=false) {

  float fX[1000];
  float fY[1000];
  float fEXlow[1000];
  float fEXhigh[1000];
  float fEYlow[1000];
  float fEYhigh[1000];  

  TAxis *xaxis = ((TH1*)h)->GetXaxis();

  double q=(1-0.6827)/2.;

  int ibin=0;
  for (Int_t i=0; i<h->GetNbinsX(); ++i) {
    float yy = h->GetBinContent(i+1);
    if (drawZeroBins || yy > 0){
      float xx = xaxis->GetBinCenter(i+1);
      fX[ibin] = xx;
      fY[ibin] = yy;
      fEXlow[ibin]  = 0.;
      fEXhigh[ibin] = 0.;
      double N = yy;
      fEYlow[ibin]  = (N==0)?0:(N-ROOT::Math::chisquared_quantile_c(1-q,2*N)/2.);
      fEYhigh[ibin] = ROOT::Math::chisquared_quantile_c(q,2*(N+1))/2.-N;
      //      cout << fEYlow[ibin] << " " << N << " " << fEYhigh[ibin] << endl;
// Old formula from http://www-cdf.fnal.gov/physics/statistics/notes/pois_eb.txt 
//       fEYlow[ibin]  = -0.5+sqrt(yy+0.25); 
//       fEYhigh[ibin] =  0.5+sqrt(yy+0.25);
      ++ibin;
    }
  }

  TGraphAsymmErrors* g = new TGraphAsymmErrors(ibin, fX, fY, fEXlow, fEXhigh, fEYlow, fEYhigh);
  h->TAttLine::Copy(*g);
  h->TAttFill::Copy(*g);
  h->TAttMarker::Copy(*g);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.9);

//   g->SetName(h->GetName());
//   g->SetTitle(h->GetTitle());

  return g;
}

void restrictXAxis(TH1F* h, Int_t countmax){
  for(Int_t i=1; i<=h->GetNbinsX(); i++)
    h->GetXaxis()->SetBinLabel(i,Form("%i",(Int_t)(h->GetXaxis()->GetBinCenter(i))));
  h->GetXaxis()->SetBinLabel(countmax+1,Form("#geq%i",(Int_t)(h->GetXaxis()->GetBinCenter(countmax+1))));
  for(Int_t i=countmax+2; i<=h->GetNbinsX()+1; i++)
    h->SetBinContent(countmax+1,h->GetBinContent(countmax+1)+h->GetBinContent(i));
  h->GetXaxis()->SetLabelSize(0.06);
  h->GetXaxis()->SetNdivisions(-h->GetNbinsX());
  h->GetXaxis()->SetRangeUser(0., (float)(countmax+1));
}

TH1F* Smooth(TH1F* hIn, Int_t factor){
  TH1F* hOut = (TH1F*)hIn->Clone();
  TH1F* hInRebinned = (TH1F*)hIn->Rebin(factor,"hInRebinned");
  hInRebinned->Scale(1./(Float_t)factor);
  for(int ibin=1; ibin<=hIn->GetNbinsX(); ibin++)
    hOut->SetBinContent(ibin, hInRebinned->GetBinContent(1+(ibin-1)/factor));
  return hOut;
}

void DrawDataMC(TCanvas* c, TH1F** h, int v, int bl, string lumiText, Bool_t logX = false, Bool_t logY = false) {

  //----- prepare canvas
  c->cd();
  if(logX) c->SetLogx();
  if(logY) c->SetLogy();
  bool first;

  //----- prepare MC histograms
  TH1F* hStacks[nProcesses];
  first = true;
  int previous = -1;;
  for(int pr=nProcesses-1; pr>=1; pr--){
    if(useProcess[pr]){
      if(REBINDYTTBAR && (pr==DY || pr==ttbar)) h[pr] = Smooth(h[pr],rebinning[v]);
      hStacks[pr] = (TH1F*)h[pr]->Clone();
      if(!first) hStacks[pr]->Add(hStacks[previous]);
      hStacks[pr]->SetFillColor(processFillColor[pr]);
      hStacks[pr]->SetLineColor(DRAWLINES?processLineColor[pr]:processFillColor[pr]);
      hStacks[pr]->GetXaxis()->SetTitleOffset(1.2);
      hStacks[pr]->GetYaxis()->SetTitleOffset(1.4);
      hStacks[pr]->GetXaxis()->SetLabelFont(42);
      hStacks[pr]->GetYaxis()->SetLabelFont(42);
      hStacks[pr]->GetXaxis()->SetTitleFont(42);
      hStacks[pr]->GetYaxis()->SetTitleFont(42);
      if(logX) hStacks[pr]->GetXaxis()->SetMoreLogLabels();
      if(logX) hStacks[pr]->GetXaxis()->SetNoExponent();
      first = false;
      previous = pr;
    }
  }

  //----- prepare data graph
  h[0]->SetMarkerStyle(20);
  h[0]->SetMarkerColor(kBlack);
  h[0]->SetMarkerSize(1.);
  //h[0]->Draw("P0 E1 SAME"); //replaced by:
  TGraphAsymmErrors* gData = getDataGraph(h[0]);

  //----- adjust Y axis and draw MC and data
  const int npoints = gData->GetN();
  Float_t gDataErrorBarUp[npoints];
  for(int i=0; i<npoints; i++) gDataErrorBarUp[i] = gData->GetY()[i] + gData->GetEYhigh()[i] ;
  Float_t cmax = TMath::Max( (Float_t)hStacks[1]->GetMaximum(), (Float_t)TMath::MaxElement(npoints,gDataErrorBarUp) );
  cmax *= logY ? 30. : 1.1 ;
  cmax *= varMaxCorrector[bl][v];
  Float_t cminlog = hStacks[1]->GetMaximum() / varMinFactor[v];
  first = true;
  for(int pr=1; pr<nProcesses; pr++){
    if(useProcess[pr]){
      if(first){
	hStacks[pr]->SetMaximum(cmax);
	if(logY) hStacks[pr]->SetMinimum(cminlog);
	hStacks[pr]->Draw("HIST");
	first = false;
      }else{
	hStacks[pr]->Draw("HIST SAME");
      }
    }
  }
  if( !(bl==fullyblind) ) gData->Draw("P");
  
  //----- draw grey area or grid over blind region
  bool doBlindingHisto = xHistoBlindLow[bl]<xHistoBlindUp[bl] && varName[v].find("M4l")!=string::npos;
  if(doBlindingHisto){
    TH1F* hBlind = new TH1F("hBlind",";;",1,xHistoBlindLow[bl],xHistoBlindUp[bl]);
    hBlind->SetBinContent(1,cmax);
    //hBlind->SetFillColorAlpha(kBlack,0.2);
    //hBlind->SetLineColorAlpha(kBlack,0.2);
    hBlind->SetFillColor(kGray+3);
    hBlind->SetFillStyle(3013); //also tried 3001 and a few others, but they look bad in pdf format
    hBlind->SetLineColorAlpha(kWhite,0.);
    hBlind->Draw("HIST SAME");
  }

  //----- draw legend
  bool doBlindingLabel = blindingLabel[bl]!="" && varName[v].find("M4l")==string::npos;
  float legLeft  = (varLegPos[v]==1 && !(bl==blindabove110 && (varName[v]=="MZ1" || varName[v]=="MZ2")) ) ? 0.2 : 0.65 ;
  float legRight = (varLegPos[v]==1 && !(bl==blindabove110 && (varName[v]=="MZ1" || varName[v]=="MZ2")) ) ? 0.5 : 0.95 ;
  if(doBlindingLabel && bl==blind110150) legRight += 0.05;
  float legUp   = 0.8;
  float legDown = doBlindingLabel ? 0.55 : (bl==fullyblind) ? 0.65 : 0.6;
  TLegend* lgd = new TLegend(legLeft,legDown,legRight,legUp);
  if(doBlindingLabel) lgd->SetHeader(blindingLabel[bl].c_str());
  lgd->SetFillStyle(0);
  lgd->SetBorderSize(0);
  if( !(bl==fullyblind) )
    lgd->AddEntry(h[0],processLabel[0].c_str(),"p");
  for(int pr=1; pr<nProcesses; pr++)
    if(useProcess[pr])
      lgd->AddEntry(hStacks[pr],processLabel[pr].c_str(),"f");
  lgd->Draw();

  //----- customize m4l axis labels
  if(DRAWLABELBYHAND && logX && (varName[v]=="M4l"||varName[v]=="M4l2")){
    TText t;
    t.SetNDC();
    t.SetTextSize(0.04);
    t.SetTextFont(42);
    t.DrawText(0.825,0.093,"600");
    t.DrawText(0.918,0.093,"800");
  }

  //----- print official CMS labels and lumi
  writeExtraText = true;
  extraText  = "Preliminary";
  lumi_sqrtS = lumiText + " (13 TeV)";
  CMS_lumi( c, 0, varCMSPos[v] );

  gPad->RedrawAxis();

  //----- print yields
  if(varName[v]=="MZ1"){
    cout<<"Yields for blinding "<<bl<<":"<<endl;
    cout<<"  expected: "<<hStacks[1]->Integral()<<endl;
    cout<<"  observed: "<<h[0]->Integral()<<endl;
  }

}

void doPlots(string outputDirectory, int variableList, int blindingList, string lumiText)
{

  setTDRStyle();

  //---------- retrieve histograms from the ROOT file
  TH1F* h1[nVariables][nBlindings][nProcesses];
  TFile* fInHistos = TFile::Open("histos_plotRunII.root");
  for(int v=0; v<nVariables; v++)
    for(int bl=0; bl<nBlindings; bl++)
      for(int fs=0; fs<nFinalStates+1; fs++)
	for(int cat=0; cat<nCategories+1; cat++)
	  for(int rs=0; rs<nResStatuses+1; rs++)
	    for(int pr=0; pr<nProcesses; pr++)
	      if( fs==nFinalStates &&
		  cat==nCategories &&
		  rs==nResStatuses ){
		h1[v][bl][pr] = (TH1F*)fInHistos->Get( Form("h1_%s_%s_%s_%s_%s_%s",varName[v].c_str(),sBlinding[bl].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str(),sResonantStatus[rs].c_str(),sProcess[pr].c_str()) );
		if(restrictCountVar[v]) restrictXAxis(h1[v][bl][pr],restrictCountVar[v]);
	      }

  //---------- do the plots (1 canvas per variable and per blinding policy)
  TCanvas* c1[nVariables];
  for(int v=0; v<nVariables; v++){
    if(!plotThisVar[variableList][v]) continue;
    for(int bl=0; bl<nBlindings; bl++){
      if(!plotThisBlinding[blindingList][bl]) continue;
      c1[v] = new TCanvas(Form("c1_%s_%s",varName[v].c_str(),sBlinding[bl].c_str()),Form("c1_%s_%s",varName[v].c_str(),sBlinding[bl].c_str()),500,500);
      DrawDataMC(c1[v],h1[v][bl],v,bl,lumiText,varLogx[v],varLogy[v]);
      SaveCanvas(outputDirectory,c1[v]);
    }
  }

}





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void plotDataVsMC() {

  // --------------- inputs ---------------

  // Define input/output location
  string inputPath = "";
  string outputPath = "PlotsDataVsMC/";
  gSystem->Exec(("mkdir -p "+outputPath).c_str());

  // Define the luminosity
  // 2015D (not v4) : 578.3/pb
  // 50ns (2015B + 2015C_50ns) : 71.52/pb
  // 2015D + 2015Dv4 : 832.31/pb
  float lumi = 0.90383;
  string lumiText = "903.8 pb^{-1}"; //"0.904 fb^{-1}"

  // Choose the list of variables that you want to plot
  int variableList = 0;

  // Choose a list of ways of blinding some m4l regions
  int blindingList = 1; 


  // --------------- processing ---------------

  // Prepare the histograms and store them in a ROOT file (to be done only once, it can take a few minutes)
  doHistograms(inputPath, lumi);

  // Do the plots (+- instantaneous)
  doPlots(outputPath, variableList, blindingList, lumiText);

}

