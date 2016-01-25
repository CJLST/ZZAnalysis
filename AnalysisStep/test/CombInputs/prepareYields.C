/* 
 * usage: 
 * -at the end of the present file, specify input trees, luminosity and m4l window
 * -run with:
 *   root -q -l -b prepareYields.C++
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFitResult.h"
#include "TPaveStats.h"
#include "TMath.h"

#include "../Plotter/tdrstyle.C"
#include "../Plotter/plotUtils.C"
#include <ZZAnalysis/AnalysisStep/src/Category.cc>

using namespace std;

#define DEBUG 0
#define MERGE2E2MU 1

#define APPLYKFACTORS 1
#define RESCALETOSMPSIGNALSTRENGTH 0
#define SMPSIGNALSTRENGTH 0.99



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////// global variables /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

enum Process {ggH=0, qqH=1, WH=2, ZH=3, ttH=4, qqZZ=5, ggZZ=6};
const int nProcesses = 7;
string sProcess[nProcesses] = {"ggH", "qqH", "WH", "ZH", "ttH", "qqZZ", "ggZZ"};

const int nMHPoints = 4;
string sMHPoint[nMHPoints] = {"", "124", "125", "126"};
Float_t fMHPoint[nMHPoints] = {0., 124., 125., 126.};
bool isSignal[nProcesses] = {1,1,1,1,1,0,0,};
Int_t nMHPointsProcess[nProcesses] = {3,1,1,1,1,0,0};
bool hasMHPoint[nProcesses][nMHPoints] = {
  {0,1,1,1},
  {0,0,1,0},
  {0,0,1,0},
  {0,0,1,0},
  {0,0,1,0},
  {1,0,0,0},
  {1,0,0,0},
};

enum FinalState {fs4mu=0, fs4e=1, fs2e2mu=2, fs2mu2e=3};
const int nFinalStates = 4;
string sFinalState[nFinalStates+1] = {"4mu", "4e", "2e2mu", "2mu2e", "4l"};
Int_t fsMarkerStyle[nFinalStates+1] = {20,22,21,33,29};

const int nCategories = 6;
string sCategory[nCategories+1] = {
  "UnTagged",
  "OneJetTagged",
  "VBFTagged", 
  "VHLeptTagged",
  "VHHadrTagged",
  "ttHTagged",
  "inclusive",
};
Color_t catColor[nCategories+1] = {kBlue-9, kCyan-6, kGreen-6, kRed-7, kOrange+6, kMagenta-6, kBlack, };

enum ResonantStatus {resonant=0, nonresonant=1};
const int nResStatuses = 2;
string sResonantStatus[nResStatuses+1] = {"resonant", "nonresonant", "allres"};

Double_t deltaR(Double_t e1, Double_t p1, Double_t e2, Double_t p2) {
  Double_t deltaPhi = acos(cos(p1-p2));
  return TMath::Sqrt((e1-e2)*(e1-e2) + deltaPhi*deltaPhi);
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////// compute and save yields /////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void computeYields(string inputFilePath, double lumi, double sqrts, double m4lMin, double m4lMax)
{

  const int nDatasets = 13;
  string datasets[nDatasets] = {
    "ggH",
    "VBFH",
    "WplusH",
    "WminusH",
    "ZH",
    "ttH",
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

  Float_t yield[nProcesses][nMHPoints][nFinalStates+1][nCategories+1][nResStatuses+1];
  for(int pr=0; pr<nProcesses; pr++)
    for(int mp=0; mp<nMHPoints; mp++)
      if(hasMHPoint[pr][mp])
	for(int fs=0; fs<nFinalStates+1; fs++)
	  for(int cat=0; cat<nCategories+1; cat++)
	    for(int rs=0; rs<nResStatuses+1; rs++)
	      yield[pr][mp][fs][cat][rs] = 0.;
  
  int currentProcess;
  int currentFinalState;
  int currentCategory;
  int currentResStatus;


  //---------- Will loop over all datasets

  for(int d=0; d<nDatasets; d++){

    //----- assign dataset to correct process
    currentProcess = -1;
    if(datasets[d]=="ggH") currentProcess = ggH;
    if(datasets[d]=="VBFH") currentProcess = qqH;
    if(datasets[d]=="WplusH") currentProcess = WH;
    if(datasets[d]=="WminusH") currentProcess = WH;
    if(datasets[d]=="ZH") currentProcess = ZH;
    if(datasets[d]=="ttH") currentProcess = ttH;
    if(datasets[d]=="ZZTo4l"||
       datasets[d]=="ZZTo4lamcatnlo") 
      currentProcess = qqZZ;
    if(datasets[d]=="ggZZ4e"||
       datasets[d]=="ggZZ4mu"||
       datasets[d]=="ggZZ4tau"||
       datasets[d]=="ggZZ2e2mu"||
       datasets[d]=="ggZZ2e2tau"||
       datasets[d]=="ggZZ2mu2tau") 
      currentProcess = ggZZ;
    
    for(int mp=0; mp<nMHPoints; mp++){

      if(!hasMHPoint[currentProcess][mp]) continue;
      
      string inputFileName = string(Form("%s%s%s/ZZ4lAnalysis.root",inputFilePath.c_str(),datasets[d].c_str(),sMHPoint[mp].c_str()));
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


      //---------- Process tree

      Long64_t entries = inputTree[d]->GetEntries();
      cout<<"Processing dataset "<<datasets[d]<<sMHPoint[mp]<<" ("<<entries<<" entries) ..."<<endl;

      for (Long64_t z=0; z<entries; ++z){

	if(DEBUG && z>1000) break;

	inputTree[d]->GetEntry(z);
      
	if( !(ZZsel>=90) ) continue;
	if(ZZMass<m4lMin || ZZMass>m4lMax) continue;

	Float_t kfactor = 1.;
	if(APPLYKFACTORS){
	  
	  if(currentProcess==qqZZ){
	    
	    //kfactor = 1.065;
	    
	    // if(GenZ1Flav==GenZ2Flav)
	    //   kfactor = 1.09;
	    // else
	    //   kfactor = 1.11;
	    
	    kfactor = 1.1;	  
	    
	  }else if(currentProcess==ggZZ){
	    //kfactor = 2.;
	    kfactor = 1.7;
	  }
	  
	}

	Double_t eventWeight = partialSampleWeight[d] * xsec * kfactor * overallEventWeight ;


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

	yield[currentProcess][mp][currentFinalState][currentCategory][currentResStatus] += eventWeight;


      } // end for entries

    } // end for mHPoints

  } // end for datasets


  //---------- Fill 'inclusive' counters
  for(int pr=0; pr<nProcesses; pr++)
    for(int mp=0; mp<nMHPoints; mp++)
      if(hasMHPoint[pr][mp])
	for(int fs=0; fs<nFinalStates; fs++)
	  for(int cat=0; cat<nCategories; cat++)
	    for(int rs=0; rs<nResStatuses; rs++)
	      yield[pr][mp][nFinalStates][cat][rs] += yield[pr][mp][fs][cat][rs];
  for(int pr=0; pr<nProcesses; pr++)
    for(int mp=0; mp<nMHPoints; mp++)
      if(hasMHPoint[pr][mp])
	for(int fs=0; fs<nFinalStates+1; fs++)
	  for(int cat=0; cat<nCategories; cat++)
	    for(int rs=0; rs<nResStatuses; rs++)
	      yield[pr][mp][fs][nCategories][rs] += yield[pr][mp][fs][cat][rs];
  for(int pr=0; pr<nProcesses; pr++)
    for(int mp=0; mp<nMHPoints; mp++)
      if(hasMHPoint[pr][mp])
	for(int fs=0; fs<nFinalStates+1; fs++)
	  for(int cat=0; cat<nCategories+1; cat++)
	    for(int rs=0; rs<nResStatuses; rs++)
	      yield[pr][mp][fs][cat][nResStatuses] += yield[pr][mp][fs][cat][rs];
  

  //---------- Write yield arrays to a ROOT file

  TFile* fOutYields = new TFile(Form("yields_%iTeV_m4l%.1f-%.1f_%.3ffb-1.root",(int)sqrts,m4lMin,m4lMax,lumi),"recreate");
  fOutYields->cd();
  for(int pr=0; pr<nProcesses; pr++)
    for(int mp=0; mp<nMHPoints; mp++)
      if(hasMHPoint[pr][mp])
	for(int fs=0; fs<nFinalStates+1; fs++)
	  for(int cat=0; cat<nCategories+1; cat++)
	    for(int rs=0; rs<nResStatuses+1; rs++){
	      TH1F* hTemp = new TH1F(Form("h1_%s%s_%s_%s_%s",sProcess[pr].c_str(),sMHPoint[mp].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str(),sResonantStatus[rs].c_str()),"",1,0,1);
	      hTemp->SetBinContent(1,yield[pr][mp][fs][cat][rs]);
	      hTemp->Write(hTemp->GetName());
	      delete hTemp;
	    }
  fOutYields->Close();
  delete fOutYields; 

}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////////////////////// fit signal yields vs. mH ///////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void DrawYieldFits(string outputDirectory, TGraph* g[nProcesses][nFinalStates][nCategories+1], TF1* f[nProcesses][nFinalStates][nCategories+1]){

  setTDRStyle();

  string outputPath = string(Form("%s/CanvasesSignalFits/",outputDirectory.c_str()));
  gSystem->Exec(("mkdir -p "+outputPath).c_str());
  TCanvas* cYield[nProcesses];

  TLegend* lgd = new TLegend(0.65,0.4,0.93,0.8);
  lgd->SetFillStyle(0);
  lgd->SetBorderSize(0);

  bool first = true;   

  for(int pr=0; pr<nProcesses; pr++){
    if(!isSignal[pr]) continue;    

    if(pr!=ggH) continue; //FIXME: for the moment, only ggH has more than one mass point

    string canvasName = string(Form("cFits_%s",sProcess[pr].c_str()));
    cYield[pr] = new TCanvas(canvasName.c_str(),canvasName.c_str(),500,500);

    for(int fs=0; fs<nFinalStates; fs++){
      if(MERGE2E2MU && fs==fs2mu2e) continue;
      for(int cat=0; cat<nCategories+1; cat++){

	g[pr][fs][cat]->GetXaxis()->SetTitle("generated m_{H}");
	g[pr][fs][cat]->GetYaxis()->SetTitle("expected yield");
	g[pr][fs][cat]->GetXaxis()->SetTitleOffset(1.2);
	g[pr][fs][cat]->GetYaxis()->SetTitleOffset(1.6);
	g[pr][fs][cat]->GetXaxis()->SetLabelFont(42);
	g[pr][fs][cat]->GetYaxis()->SetLabelFont(42);
	g[pr][fs][cat]->GetXaxis()->SetLabelSize(0.04);
	g[pr][fs][cat]->GetYaxis()->SetLabelSize(0.04);
	g[pr][fs][cat]->GetXaxis()->SetTitleFont(42);
	g[pr][fs][cat]->GetYaxis()->SetTitleFont(42);
	g[pr][fs][cat]->GetXaxis()->SetTitleSize(0.04);
	g[pr][fs][cat]->GetYaxis()->SetTitleSize(0.04);
	g[pr][fs][cat]->SetMarkerStyle(fsMarkerStyle[fs]);
	g[pr][fs][cat]->SetMarkerColor(catColor[cat]);
	g[pr][fs][cat]->SetMinimum(0);
	g[pr][fs][cat]->GetXaxis()->SetLimits(123.,129.2);

	if(cat==nCategories) continue;

	g[pr][fs][cat]->Draw(first?"AP":"P");

	f[pr][fs][cat]->SetLineColor(catColor[cat]);
	f[pr][fs][cat]->SetLineWidth(1);
	f[pr][fs][cat]->Draw("SAME");
	if(fs==0) lgd->AddEntry(f[pr][fs][cat],sCategory[cat].c_str(),"l");

	first = false;
      }
    }

    for(int fs=0; fs<nFinalStates; fs++)
      if(!(MERGE2E2MU && fs==fs2mu2e))
	lgd->AddEntry(g[pr][fs][nCategories],sFinalState[fs].c_str(),"p");

    lgd->Draw();
    SaveCanvas(outputPath,cYield[pr]);
  }

}

void fitSignalYields(string outputDirectory, double lumi, double sqrts, double m4lMin, double m4lMax)
{

  //---------- Retrieve yields from the ROOT file

  TFile* fInYields = TFile::Open(Form("yields_%iTeV_m4l%.1f-%.1f_%.3ffb-1.root",(int)sqrts,m4lMin,m4lMax,lumi));
  Float_t yield[nProcesses][nMHPoints][nFinalStates][nCategories+1];
  for(int pr=0; pr<nProcesses; pr++)
    for(int mp=0; mp<nMHPoints; mp++)
      if(hasMHPoint[pr][mp] && isSignal[pr])
	for(int fs=0; fs<nFinalStates; fs++)
	  for(int cat=0; cat<nCategories+1; cat++){
	    yield[pr][mp][fs][cat] = ((TH1F*)fInYields->Get(Form("h1_%s%s_%s_%s_%s",sProcess[pr].c_str(),sMHPoint[mp].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str(),sResonantStatus[nResStatuses].c_str())))->GetBinContent(1);
	    if(MERGE2E2MU && fs==fs2e2mu) 
	      yield[pr][mp][fs][cat] += yield[pr][mp][fs2mu2e][cat];
	  }

  //---------- Fit yield vs. mH

  TFile* fOutYieldFunctions = new TFile(Form("yieldFunctions_%iTeV_m4l%.1f-%.1f_%.3ffb-1.root",(int)sqrts,m4lMin,m4lMax,lumi),"recreate");
  TGraph* gYield[nProcesses][nFinalStates][nCategories+1];
  TF1* fYield[nProcesses][nFinalStates][nCategories+1];

  const int nParameters = 3; // Let's take a 2nd order polynomial for now.
  Float_t fitParameters[nParameters][nProcesses][nFinalStates][nCategories+1];

  for(int pr=0; pr<nProcesses; pr++){
    if(!isSignal[pr]) continue;

    if(pr!=ggH) continue; //FIXME: for the moment, only ggH has more than one mass point

    for(int fs=0; fs<nFinalStates; fs++){
      if(MERGE2E2MU && fs==fs2mu2e) continue;

      for(int cat=0; cat<nCategories+1; cat++){

	gYield[pr][fs][cat] = new TGraph(nMHPointsProcess[pr]);

	int iPoint = 0;
	for(int mp=0; mp<nMHPoints; mp++){
	  if(hasMHPoint[pr][mp]){
	    gYield[pr][fs][cat]->SetPoint(iPoint,fMHPoint[mp],yield[pr][mp][fs][cat]);
	    iPoint++;
	  }
	}

	string fName = (string)Form("f_%s_%s_%s",sProcess[pr].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str());
	fYield[pr][fs][cat] = new TF1(fName.c_str(),"pol2",124,126); 
	//fYield[pr][fs][cat] = new TF1(fName.c_str(),"[0]+[1]*x+[2]*x*x",124,126); // the fit doesn't work with this
	
	cout<<endl<<sProcess[pr]<<", "<<sFinalState[fs]<<", "<<sCategory[cat]<<endl;
	gYield[pr][fs][cat]->Fit(fYield[pr][fs][cat],"N S");

	//cout<<fYield[pr][fs][cat]->GetExpFormula("P")<<endl;
	fYield[pr][fs][cat]->Write();
	
      }
    }
  }

  DrawYieldFits(outputDirectory,gYield,fYield);

  fOutYieldFunctions->Close();
  delete fOutYieldFunctions; 

}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////// prepare fragments ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void generateFragments(string outputDirectory, double lumi, double sqrts, double m4lMin, double m4lMax, string mHoption)
{

  //---------- Retrieve yields and functions from the ROOT file

  Float_t yield[nProcesses][nFinalStates+1][nCategories+1];
  TFile* fInYields = TFile::Open(Form("yields_%iTeV_m4l%.1f-%.1f_%.3ffb-1.root",(int)sqrts,m4lMin,m4lMax,lumi));
  for(int pr=0; pr<nProcesses; pr++)
    for(int fs=0; fs<nFinalStates+1; fs++)
      for(int cat=0; cat<nCategories+1; cat++){
	yield[pr][fs][cat] = ((TH1F*)fInYields->Get(Form("h1_%s%s_%s_%s_%s",sProcess[pr].c_str(),(isSignal[pr]?(mHoption!="param"?mHoption.c_str():"125"):""),sFinalState[fs].c_str(),sCategory[cat].c_str(),sResonantStatus[nResStatuses].c_str())))->GetBinContent(1);
	if(RESCALETOSMPSIGNALSTRENGTH && (pr==qqZZ||pr==ggZZ)) yield[pr][fs][cat] *= SMPSIGNALSTRENGTH; //FIXME: to be removed at some point ?
      }

  TF1* fYield[nProcesses][nFinalStates][nCategories+1];
  if(mHoption=="param"){
    TFile* fInYieldFunctions = TFile::Open(Form("yieldFunctions_%iTeV_m4l%.1f-%.1f_%.3ffb-1.root",(int)sqrts,m4lMin,m4lMax,lumi));
    for(int pr=0; pr<nProcesses; pr++){
      if(!isSignal[pr]) continue;
      if(pr!=ggH) continue; //FIXME: for the moment, only ggH has more than one mass point
      for(int fs=0; fs<nFinalStates; fs++){
	if(MERGE2E2MU && fs==fs2mu2e) continue;
	for(int cat=0; cat<nCategories+1; cat++)
	  fYield[pr][fs][cat] = (TF1*)fInYieldFunctions->Get(Form("f_%s_%s_%s",sProcess[pr].c_str(),sFinalState[fs].c_str(),sCategory[cat].c_str()));	  
      }
    }
  }


  //---------- Prepare the yaml fragments

  float totalYield = 0.;
  string outputFileName[nFinalStates];
  ofstream outFile[nFinalStates];
  for(int fs=0; fs<nFinalStates; fs++){
    if(MERGE2E2MU && fs==fs2mu2e) continue;
    
    outputFileName[fs] = string(Form("%s/yields_%iTeV_%s.yaml",outputDirectory.c_str(),(int)sqrts,sFinalState[fs].c_str())); // (not specifying MH and mass window in the name here)
    outFile[fs].open(outputFileName[fs]);
    
    outFile[fs]<<"---"<<endl;
    outFile[fs]<<"# mass window: "<<m4lMin<<" <= m4l <= "<<m4lMax<<endl;
    outFile[fs]<<"# sqrt(s) = "<<sqrts<<" TeV"<<endl;
    outFile[fs]<<"# integrated luminosity = "<<lumi<<" fb-1"<<endl;
    outFile[fs]<<endl;

    // outFile[fs]<<"# Category numbering convention:"<<endl;
    // for(int cat=0; cat<nCategories; cat++){
    //   outFile[fs]<<"# "<<cat<<" "<<sCategory[cat]<<endl;
    // }
    // outFile[fs]<<endl;
    
    for(int cat=0; cat<nCategories; cat++){
      outFile[fs]<<sCategory[cat]<<": "<<endl;

      for(int pr=0; pr<nProcesses; pr++){
	float y = yield[pr][fs][cat];
	if(MERGE2E2MU && fs==fs2e2mu) y += yield[pr][fs2mu2e][cat];

	if(pr==ggH && mHoption=="param"){ //FIXME: for the moment, only ggH has more than one mass point
	  outFile[fs]<<"    "<<sProcess[pr]<<": "
		     <<"'("<<fYield[pr][fs][cat]->GetParameter(0)<<")+("
		     <<fYield[pr][fs][cat]->GetParameter(1)<<"*MH)+("
		     <<fYield[pr][fs][cat]->GetParameter(2)<<"*MH*MH)'"<<endl; //FIXME: expression is hardcoded here
	}else{
	  outFile[fs]<<"    "<<sProcess[pr]<<": "<<y<<endl;
	}

	totalYield += y;
      }
      outFile[fs]<<endl;
    }

    outFile[fs].close();
  }

  cout<<"total signal yield in mass window ["<<m4lMin<<","<<m4lMax<<"] : "<<totalYield<<endl;
  for(int pr=0; pr<nProcesses; pr++)
    cout<<" "<<sProcess[pr]<<": "<<yield[pr][nFinalStates][nCategories]<<endl;

}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// main function ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void prepareYields(bool recomputeYields = true) {


  // --------------- definitions ---------------

  // Define input/output location
  string inputPath = "";
  string outputPath = "YieldFiles";

  // Define c.o.m. energy (13TeV only for the moment)
  float sqrts = 13.;

  // Define the luminosity

  // all 50ns (2015B + 2015C_50ns) : 71.52/pb
  // 25ns 2015D (no v4) : 578.3/pb
  // 25ns 2015D + 2015Dv4 (Oct. 17th JSON) : 832.31/pb
  // all 25ns (2015C + 2015D + 2015Dv4) (Nov. 13th Silver JSON) : 2.46/fb

  //float lumi = 0.90383;
  float lumi = 2.6;

  // m4l window
  float m4l_min = 105.;
  float m4l_max = 140.;


  // --------------- processing ---------------

  gSystem->Exec(("mkdir -p "+outputPath).c_str());

  // Compute the yields for all available processes and mH values, and store them in a ROOT file 
  // (to be done only once, it can take a few minutes)
  if(recomputeYields)
    computeYields(inputPath, lumi, sqrts, m4l_min, m4l_max);

  // Parameterize signal yields as a function of mH
  fitSignalYields(outputPath, lumi, sqrts, m4l_min, m4l_max);

  // Prepare the yaml files
  //generateFragments(outputPath, lumi, sqrts, m4l_min, m4l_max, "125"); // This takes mH=125 for the signal yield.
  generateFragments(outputPath, lumi, sqrts, m4l_min, m4l_max, "param"); // This puts the yield(mH) expression when possible.


}

