#define HZZ4l_cxx

#include "HZZ4l.h"

//Root includes
#include <TH1F.h>
#include <TH2D.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TString.h>
#include <TRegexp.h>
#include <TCanvas.h>

//RooFit includes
#include <RooRandom.h>

//Std includes
#include <boost/regex.hpp>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "../Plotter/root_lib/XSecReader.h"
#include "../../interface/PUReweight.h"
#include "../../interface/FinalStates.h"
#include "../../src/bitops.cc"
#include "../../src/Fisher.cc"

using namespace std;

namespace {
  const bool doMuEffCorr    = true;
  const bool doEleEffCorr   = true;
  const bool doHighmassCorr = true; 
  const bool doHqTCorr      = false;
  const bool saveJets       = true;
  const bool applySystMu    = false;
  const bool applySystEle   = false;

  const Int_t pwhg_flag = 0;    // 0 means standard high mass weights
                                // 1 means CPS+ + Interference
                                // 2 means CPS- + Interference
                                // 3 means CPS  + Interference+
                                // 4 means CPS  + Interference-

  const Float_t Run2011AFraction = 0.465;

  XSecReader xsecRead7TeV("../Plotter/Xsection7TeV_YR3.txt","../Plotter/Luminosity.txt");
  XSecReader xsecRead8TeV("../Plotter/Xsection8TeV_YR3.txt","../Plotter/Luminosity.txt");

  PUReweight PUWeighter(PUReweight::LEGACY);

  //This is used to always give the same shifts when doing systematics
  int Seeder = 0;
}

HZZ4l::HZZ4l(TChain *tree, TString sampleName) : HZZ4lBase(tree), theSample(sampleName)
{
  isCR      = false;
  isSignal      = false;
  if (theSample.BeginsWith("ZZ4lAnalysis_")) theSample.Remove(0,13);
  for (int i=0; i<4; ++i) ZXWeightTables[i]=0;
}

void HZZ4l::Loop(Int_t channelType, const TString outputName)
{
  if (fChain == 0) return;

  TCanvas dummy; // Hack o avoid "no dictionary for class TPaletteAxis is available" warning
  
  //Set the channel name
  string channelname;
  if(channelType == 0) channelname = "4mu";
  else if(channelType == 1) channelname = "4e";
  else if(channelType == 2) channelname = "2e2mu";
  else abort();

  TString chainName(fChain->GetName());
  static TDirectory *filedir = gDirectory;
  if (chainName.Contains("CR")) nEventComplete = (TH1F*)filedir->Get("CRZLLTree/Counters");
  else nEventComplete = (TH1F*)filedir->Get(TString("ZZ" + channelname + "Tree/Counters"));

  Float_t Nevt_Gen = nEventComplete->GetBinContent(1);            // 
  Float_t Nevt_Gen_weighted = nEventComplete->GetBinContent(0); // Weighted #of events, can be used only when no filter was applied

  // Sample types 
  bool isData        = false;
  bool isVBF         = false;
  bool isNewHighmass = false; // for POWHEG15
  bool isMinlo       = false;
  Float_t mPOLE      = 0.;    // nominal H mass, for signals

  // Normalization variables
  Float_t MC_weight_initial        = 1.;
  Float_t MC_weight_norm_initial   = 1.;
  Float_t MC_weight_noxsec_initial = 1.;
  Float_t MC_weight_interm         = 1.;
  Float_t MC_weight_norm_interm    = 1.;
  Float_t MC_weight_noxsec_interm  = 1.;
  Float_t MC_weight                = 1.;
  Float_t MC_weight_norm           = 1.;
  Float_t MC_weight_noxsec         = 1.;
  Float_t PUWeight                 = 1.;
  Float_t powheg_weight            = 1.;
  Float_t dataMCWeight             = 1.;
  Float_t HqTWeight                = 1.;
  Float_t ZXfake_weight            = 1.;

  cout << "Sample name :                          " << theSample << std::endl;
  cout << "Name of the output file :              " << outputName << std::endl;
  cout << "Number of total entries :              " << fChain->GetEntries() << std::endl ;
  cout << "Number of generated events :           " << int(Nevt_Gen) << std::endl; 


  if(theSample.Contains("DoubleMu") || theSample.Contains("DoubleEle") || theSample.Contains("DoubleOr") || theSample.Contains("MuEG")) isData = true;
  if(!isData){
    if(theSample.Contains("H")) isSignal = true;
    if(theSample.Contains("VBF")) isVBF = true;
    if(theSample.Contains("powheg15")) isNewHighmass = true;
    if(theSample.Contains("minlo")) isMinlo = true;

    //Reweighting for high mass
    if(isSignal && doHighmassCorr){
      TString massString = theSample(TRegexp("[^H]*$")); // ASSUME THAT SIGNAL SAMPLE NAMES ENDS WITH "H[mass]" eg "H125.6"
      if (massString(TRegexp("[a-zA-Z]"))!="") cout << "ERROR: malformed signal sample name " << theSample << " " << massString << endl<<endl<<endl;
      mPOLE = massString.Atof();
      if(mPOLE > 399.) getWeightFromFile(massString, isVBF, isNewHighmass);
    }

    cout << "Number of generated events, weighted : " << Nevt_Gen_weighted << std::endl;
    if (isSignal) cout << "mPOLE: " << mPOLE << " VBF: " << isVBF << " powheg15: " << isNewHighmass << " minlo: " << isMinlo << endl;

    if (!isMinlo && Nevt_Gen_weighted!=0 && Nevt_Gen_weighted != Nevt_Gen)  cout << "WARNING: Nevt_Gen_weighted!= Nevt_Gen" << endl;

    float Nevt_norm = Nevt_Gen;

    // normalize to the sum of weights, instead of generated events (only for MINLO at the moment; all other samples have genMCWeight=1
    // (although it would make sense to always normalize to sum of weights that also includes the PU weight; but filtered samples have to be
    // handled appropriately)
    if (isMinlo && Nevt_Gen_weighted!=0) { 
      Nevt_norm = Nevt_Gen_weighted;
    }

    //Initial values of MC_weight, MC_weight_norm, MC_weight_noxsec
    if(!is8TeV) MC_weight_initial = xsecRead7TeV.getWeight(theSample, "1fb-1","all", true)/Nevt_norm;
    else MC_weight_initial = xsecRead8TeV.getWeight(theSample, "1fb-1","all", true)/Nevt_norm;

    MC_weight_norm_initial = getNormalizedWeight(channelType);
    if (MC_weight_norm_initial<0){
      cout << "ERROR: MC normalization negative " << MC_weight_norm_initial  << endl;
      abort();
    }
    MC_weight_noxsec_initial = 1./Nevt_norm;

  }


  TFile* fHqT;
  if (doHqTCorr){
    // Get the file for HqT reweighting
    fHqT = TFile::Open("HqTWeights.root","READ");
  }

  const Long64_t nentries = fChain->GetEntriesFast();

  //counters
  Int_t nTOTEv      = 0;
  Float_t nBestCand = 0.;
  Int_t nSel        = 0;

  //Variables necessary for systematics
  Int_t N_syst_Iter = 500;
  Float_t nBestSyst[N_syst_Iter];


  //Event variables
  Long64_t myEventNumber  = 0;
  Int_t    myRunNumber    = 0;
  Int_t    myLumiNumber    = 0;
  Short_t  mygenProcessId = 0;
  Float_t  mygenHEPMCweight = 1.;
  Float_t  mygenhpt       = 0;
  Float_t  mygenhmass     = 0;
  Int_t    myCRflag       = 0;

  //ZZ variables
  Float_t myZZMass        = 0.;
  Float_t myZZPt          = 0.;
  Float_t myZZEta         = 0.;
  Float_t myZZPhi         = 0.;
  //  Float_t myZZRapidity    = 0.;
  Float_t myZZMassErr     = 0.;
  Float_t myZZMassErrCorr = 0.;
  Float_t myZZFisher      = 0.;

  //Z variables
  Float_t myZ1Mass        = 0.;
  Float_t myZ1Pt          = 0.;
  Short_t myZ1ids         = 0.;
  Float_t myZ2Mass        = 0.;
  Float_t myZ2Pt          = 0.;
  Short_t myZ2ids         = 0.;
  Float_t mycosthetastar  = 0.;
  Float_t myhelphi        = 0.;
  Float_t myhelcosthetaZ1 = 0.;
  Float_t myhelcosthetaZ2 = 0.;
  Float_t myphistarZ1     = 0.;
  Float_t myphistarZ2     = 0.;
  Float_t myxi            = 0.;
  Float_t myxistar        = 0.;

  //Lepton variables
  Float_t myLep1Pt        = 0.;
  Float_t myLep1Eta       = 0.;
  Int_t   myLep1ID        = 0.;
  Float_t myLep2Pt        = 0.;
  Float_t myLep2Eta       = 0.;
  Int_t   myLep2ID        = 0.;
  Float_t myLep3Pt        = 0.;
  Float_t myLep3Eta       = 0.;
  Int_t   myLep3ID        = 0.;
  Float_t myLep4Pt        = 0.;
  Float_t myLep4Eta       = 0.;
  Int_t   myLep4ID        = 0.;

  // Lepton scale and resolution systematics on probabilities
  Float_t myp0plus_m4l_ScaleUp   = 0.;
  Float_t myp0plus_m4l_ScaleDown = 0.;
  Float_t myp0plus_m4l_ResUp     = 0.;
  Float_t myp0plus_m4l_ResDown   = 0.;
  Float_t mybkg_m4l_ScaleUp      = 0.;
  Float_t mybkg_m4l_ScaleDown    = 0.;
  Float_t mybkg_m4l_ResUp        = 0.;
  Float_t mybkg_m4l_ResDown      = 0.;

  //KD variables
  Float_t myp0plus_VAJHU  = 0.;
  Float_t myp0hplus_VAJHU = 0.;
  Float_t myp0minus_VAJHU = 0.;
  Float_t myp0plus_VAMCFM   = 0.;
  Float_t myp1_VAJHU      = 0.;
  Float_t myp1_prodIndep_VAJHU      = 0.;
  Float_t myp1plus_VAJHU  = 0.;
  Float_t myp1plus_prodIndep_VAJHU      = 0.;
  Float_t myp2_prodIndep_VAJHU      = 0.;
  Float_t myp2hplus_VAJHU  = 0.;
  Float_t myp2hminus_VAJHU  = 0.;
  Float_t myp2bplus_VAJHU  = 0.;
  Float_t myp2_VAJHU      = 0.;
  Float_t myp2qqb_VAJHU   = 0.;

  Float_t myp2hplus_qqb_VAJHU = 0;
  Float_t myp2hplus_prodIndep_VAJHU = 0;
  Float_t myp2hminus_qqb_VAJHU = 0;
  Float_t myp2hminus_prodIndep_VAJHU = 0;
  Float_t myp2bplus_qqb_VAJHU = 0;
  Float_t myp2bplus_prodIndep_VAJHU = 0;
  Float_t myp2h2plus_gg_VAJHU = 0;
  Float_t myp2h2plus_qqbar_VAJHU = 0;
  Float_t myp2h2plus_prodIndep_VAJHU = 0;
  Float_t myp2h3plus_gg_VAJHU = 0;
  Float_t myp2h3plus_qqbar_VAJHU = 0;
  Float_t myp2h3plus_prodIndep_VAJHU = 0;
  Float_t myp2h6plus_gg_VAJHU = 0;
  Float_t myp2h6plus_qqbar_VAJHU = 0;
  Float_t myp2h6plus_prodIndep_VAJHU = 0;
  Float_t myp2h7plus_gg_VAJHU = 0;
  Float_t myp2h7plus_qqbar_VAJHU = 0;
  Float_t myp2h7plus_prodIndep_VAJHU = 0;
  Float_t myp2h9minus_gg_VAJHU = 0;
  Float_t myp2h9minus_qqbar_VAJHU = 0;
  Float_t myp2h9minus_prodIndep_VAJHU = 0;
  Float_t myp2h10minus_gg_VAJHU = 0;
  Float_t myp2h10minus_qqbar_VAJHU = 0;
  Float_t myp2h10minus_prodIndep_VAJHU = 0;
  Float_t mybkg_VAMCFM    = 0.;
  Float_t mybkg_prodIndep_VAMCFM= 0.;
  Float_t myggzz_VAMCFM   = 0.;
  Float_t myggzz_p0plus_VAMCFM   = 0.;
  Float_t myggzz_c1_VAMCFM   = 0.;
  Float_t myggzz_c5_VAMCFM   = 0.;
  Float_t myggzz_ci_VAMCFM   = 0.;
  Float_t myphjj_VAJHU_old = 0.;
  Float_t mypvbf_VAJHU_old = 0.;
  Float_t myphjj_VAJHU_old_up = 0.;
  Float_t mypvbf_VAJHU_old_up = 0.;
  Float_t myphjj_VAJHU_old_dn = 0.;
  Float_t mypvbf_VAJHU_old_dn = 0.;
  Float_t myphjj_VAJHU_new = 0.;
  Float_t mypvbf_VAJHU_new = 0.;
  Float_t myphjj_VAJHU_new_up = 0.;
  Float_t mypvbf_VAJHU_new_up = 0.;
  Float_t myphjj_VAJHU_new_dn = 0.;
  Float_t mypvbf_VAJHU_new_dn = 0.;
  Float_t myp0_g1prime2_VAJHU = 0;
  Float_t mypg1g1prime2_VAJHU = 0;
  Float_t myDgg10_VAMCFM= 0;
  Float_t myp0plus_m4l    = 0.;
  Float_t mybkg_m4l       = 0.;
  Float_t mypg1g4_mela = 0.;
  Float_t mypg1g4_VAJHU= 0.;
  Float_t mypg1g4_pi2_VAJHU= 0.;
  Float_t mypg1g2_pi2_VAJHU= 0.;
  Float_t mypg1g2_mela= 0.;
  Float_t mypg1g2_VAJHU= 0.; 

  Float_t mypzzzg_VAJHU=0.;
  Float_t mypzzgg_VAJHU=0.;
  Float_t mypzzzg_PS_VAJHU=0.;
  Float_t mypzzgg_PS_VAJHU=0.;
  Float_t myp0Zgs_VAJHU=0.;
  Float_t myp0gsgs_VAJHU=0.;
  Float_t myp0Zgs_PS_VAJHU=0.;
  Float_t myp0gsgs_PS_VAJHU=0.;
  //Jet variables
  Float_t myDiJetMass      = -99.;
  Float_t myDiJetMassPlus  = -99.;
  Float_t myDiJetMassMinus = -99.;
  Float_t myDiJetDEta      = -99.;
  Short_t myNJets30        = -99.; 
  vector<double> myJetPt;
  vector<double> myJetSigma;
  vector<double> myJetEta;
  vector<double> myJetPhi;
  vector<double> myJetMass;
  vector<double> myJetBTag;

  //Firstly create the File before the tree, because ROOT is designed my fucking monkeys high on crack
  TFile fOut(outputName,"RECREATE");

  //Validation plots
  TH1F hLepGenPt("hLepGenPt","hLepGenPt",1000,0,1000);

  //Weight plots
  TH1F hNvtxNoWeight("hNvtxNoWeight","hNvtxNoWeight",100,0,100);
  TH1F hNvtxWeight("hNvtxWeight","hNvtxWeight",100,0,100);
  TH1F hPUWeight("hPUWeight","hPUWeight",250,0,25);
  TH1F hTagProbeWeight("hTagProbeWeight","hTAGProbeWeight",300,0,10);

  //ZZMass plot with and without high mass weight
  TH1F hZZmassNoPoweg("hZZmassNoPowheg","hZZmassNoPowheg",700,0.,1800.);
  TH1F hZZmassPoweg("hZZmassPowheg","hZZmassPowheg",700,0.,1800.);
  TH1F hZZmassPowegReco("hZZmassPowhegReco","hZZmassPowhegReco",700,0.,1800.);

  TTree SelTree("SelectedTree","The Selection Tree");
  SelTree.Branch("EventNumber",&myEventNumber,"EventNumber/L");
  SelTree.Branch("LumiNumber",&myLumiNumber,"LumiNumber/I");
  SelTree.Branch("RunNumber",&myRunNumber,"RunNumber/I");
  SelTree.Branch("ZZMass",&myZZMass,"ZZMass/F");
  SelTree.Branch("ZZMassErr",&myZZMassErr,"ZZMassErr/F");
  SelTree.Branch("ZZMassErrCorr",&myZZMassErrCorr,"ZZMassErrCorr/F");
  SelTree.Branch("ZZPt",&myZZPt,"ZZPt/F");
  SelTree.Branch("ZZEta",&myZZEta,"ZZEta/F");
  SelTree.Branch("ZZPhi",&myZZPhi,"ZZPhi/F");
  SelTree.Branch("ZZFisher",&myZZFisher,"ZZFisher/F");
  SelTree.Branch("p0plus_VAJHU",&myp0plus_VAJHU,"p0plus_VAJHU/F");
  SelTree.Branch("p0hplus_VAJHU",&myp0hplus_VAJHU,"p0hplus_VAJHU/F");
  SelTree.Branch("p0minus_VAJHU",&myp0minus_VAJHU,"p0minus_VAJHU/F");
  SelTree.Branch("p0plus_VAMCFM",&myp0plus_VAMCFM,"p0plus_VAMCFM/F");
  SelTree.Branch("p1_VAJHU",&myp1_VAJHU,"p1_VAJHU/F");
  SelTree.Branch("p1_prodIndep_VAJHU",&myp1_prodIndep_VAJHU,"p1_prodIndep_VAJHU/F");
  SelTree.Branch("p1plus_VAJHU",&myp1plus_VAJHU,"p1plus_VAJHU/F");
  SelTree.Branch("p1plus_prodIndep_VAJHU",&myp1plus_prodIndep_VAJHU,"p1plus_prodIndep_VAJHU/F");
  SelTree.Branch("p2_prodIndep_VAJHU",&myp2_prodIndep_VAJHU,"p2_prodIndep_VAJHU/F");
  SelTree.Branch("p2hplus_VAJHU",&myp2hplus_VAJHU,"p2hplus_VAJHU/F");
  SelTree.Branch("p2hminus_VAJHU",&myp2hminus_VAJHU,"p2hminus_VAJHU/F");
  SelTree.Branch("p2bplus_VAJHU",&myp2bplus_VAJHU,"p2bplus_VAJHU/F");
  SelTree.Branch("p2_VAJHU",&myp2_VAJHU,"p2_VAJHU/F");
  SelTree.Branch("p2qqb_VAJHU",&myp2qqb_VAJHU,"p2qqb_VAJHU/F");

  SelTree.Branch("p2hplus_qqb_VAJHU", &myp2hplus_qqb_VAJHU, "p2hplus_qqb_VAJHU/F");
  SelTree.Branch("p2hplus_prodIndep_VAJHU", &myp2hplus_prodIndep_VAJHU, "p2hplus_prodIndep_VAJHU/F");
  SelTree.Branch("p2hminus_qqb_VAJHU", &myp2hminus_qqb_VAJHU, "p2hminus_qqb_VAJHU/F");
  SelTree.Branch("p2hminus_prodIndep_VAJHU", &myp2hminus_prodIndep_VAJHU, "p2hminus_prodIndep_VAJHU/F");
  SelTree.Branch("p2bplus_qqb_VAJHU", &myp2bplus_qqb_VAJHU, "p2bplus_qqb_VAJHU/F");
  SelTree.Branch("p2bplus_prodIndep_VAJHU", &myp2bplus_prodIndep_VAJHU, "p2bplus_prodIndep_VAJHU/F");
  SelTree.Branch("p2h2plus_gg_VAJHU", &myp2h2plus_gg_VAJHU, "p2h2plus_gg_VAJHU/F");
  SelTree.Branch("p2h2plus_qqbar_VAJHU", &myp2h2plus_qqbar_VAJHU, "p2h2plus_qqbar_VAJHU/F");
  SelTree.Branch("p2h2plus_prodIndep_VAJHU", &myp2h2plus_prodIndep_VAJHU, "p2h2plus_prodIndep_VAJHU/F");
  SelTree.Branch("p2h3plus_gg_VAJHU", &myp2h3plus_gg_VAJHU, "p2h3plus_gg_VAJHU/F");
  SelTree.Branch("p2h3plus_qqbar_VAJHU", &myp2h3plus_qqbar_VAJHU, "p2h3plus_qqbar_VAJHU/F");
  SelTree.Branch("p2h3plus_prodIndep_VAJHU", &myp2h3plus_prodIndep_VAJHU, "p2h3plus_prodIndep_VAJHU/F");
  SelTree.Branch("p2h6plus_gg_VAJHU", &myp2h6plus_gg_VAJHU, "p2h6plus_gg_VAJHU/F");
  SelTree.Branch("p2h6plus_qqbar_VAJHU", &myp2h6plus_qqbar_VAJHU, "p2h6plus_qqbar_VAJHU/F");
  SelTree.Branch("p2h6plus_prodIndep_VAJHU", &myp2h6plus_prodIndep_VAJHU, "p2h6plus_prodIndep_VAJHU/F");
  SelTree.Branch("p2h7plus_gg_VAJHU", &myp2h7plus_gg_VAJHU, "p2h7plus_gg_VAJHU/F");
  SelTree.Branch("p2h7plus_qqbar_VAJHU", &myp2h7plus_qqbar_VAJHU, "p2h7plus_qqbar_VAJHU/F");
  SelTree.Branch("p2h7plus_prodIndep_VAJHU", &myp2h7plus_prodIndep_VAJHU, "p2h7plus_prodIndep_VAJHU/F");
  SelTree.Branch("p2h9minus_gg_VAJHU", &myp2h9minus_gg_VAJHU, "p2h9minus_gg_VAJHU/F");
  SelTree.Branch("p2h9minus_qqbar_VAJHU", &myp2h9minus_qqbar_VAJHU, "p2h9minus_qqbar_VAJHU/F");
  SelTree.Branch("p2h9minus_prodIndep_VAJHU", &myp2h9minus_prodIndep_VAJHU, "p2h9minus_prodIndep_VAJHU/F");
  SelTree.Branch("p2h10minus_gg_VAJHU", &myp2h10minus_gg_VAJHU, "p2h10minus_gg_VAJHU/F");
  SelTree.Branch("p2h10minus_qqbar_VAJHU", &myp2h10minus_qqbar_VAJHU, "p2h10minus_qqbar_VAJHU/F");
  SelTree.Branch("p2h10minus_prodIndep_VAJHU", &myp2h10minus_prodIndep_VAJHU, "p2h10minus_prodIndep_VAJHU/F");
  SelTree.Branch("bkg_VAMCFM",&mybkg_VAMCFM,"bkg_VAMCFM/F");
  SelTree.Branch("bkg_prodIndep_VAMCFM",&mybkg_prodIndep_VAMCFM,"bkg_prodIndep_VAMCFM/F");
  SelTree.Branch("ggzz_VAMCFM",&myggzz_VAMCFM,"ggzz_VAMCFM/F");
  SelTree.Branch("ggzz_p0plus_VAMCFM",&myggzz_p0plus_VAMCFM,"ggzz_p0plus_VAMCFM/F");
  SelTree.Branch("ggzz_c1_VAMCFM",&myggzz_c1_VAMCFM,"ggzz_c1_VAMCFM/F");
  SelTree.Branch("ggzz_c5_VAMCFM",&myggzz_c5_VAMCFM,"ggzz_c5_VAMCFM/F");
  SelTree.Branch("ggzz_ci_VAMCFM",&myggzz_ci_VAMCFM,"ggzz_ci_VAMCFM/F");
  SelTree.Branch("phjj_VAJHU_old",&myphjj_VAJHU_old,"phjj_VAJHU_old/F");
  SelTree.Branch("pvbf_VAJHU_old",&mypvbf_VAJHU_old,"pvbf_VAJHU_old/F");
  SelTree.Branch("phjj_VAJHU_old_up",&myphjj_VAJHU_old_up,"phjj_VAJHU_old_up/F");
  SelTree.Branch("pvbf_VAJHU_old_up",&mypvbf_VAJHU_old_up,"pvbf_VAJHU_old_up/F");
  SelTree.Branch("phjj_VAJHU_old_dn",&myphjj_VAJHU_old_dn,"phjj_VAJHU_old_dn/F");
  SelTree.Branch("pvbf_VAJHU_old_dn",&mypvbf_VAJHU_old_dn,"pvbf_VAJHU_old_dn/F");
  SelTree.Branch("phjj_VAJHU_new",&myphjj_VAJHU_new,"phjj_VAJHU_new/F");
  SelTree.Branch("pvbf_VAJHU_new",&mypvbf_VAJHU_new,"pvbf_VAJHU_new/F");
  SelTree.Branch("phjj_VAJHU_new_up",&myphjj_VAJHU_new_up,"phjj_VAJHU_new_up/F");
  SelTree.Branch("pvbf_VAJHU_new_up",&mypvbf_VAJHU_new_up,"pvbf_VAJHU_new_up/F");
  SelTree.Branch("phjj_VAJHU_new_dn",&myphjj_VAJHU_new_dn,"phjj_VAJHU_new_dn/F");
  SelTree.Branch("pvbf_VAJHU_new_dn",&mypvbf_VAJHU_new_dn,"pvbf_VAJHU_new_dn/F");
  SelTree.Branch("p0_g1prime2_VAJHU",&myp0_g1prime2_VAJHU,"p0_g1prime2_VAJHU/F");
  SelTree.Branch("pg1g1prime2_VAJHU",&mypg1g1prime2_VAJHU,"pg1g1prime2_VAJHU/F");
  SelTree.Branch("Dgg10_VAMCFM",&myDgg10_VAMCFM,"Dgg10_VAMCFM/F");

  SelTree.Branch("p0plus_m4l",&myp0plus_m4l,"p0plus_m4l/F");
  SelTree.Branch("bkg_m4l",&mybkg_m4l,"bkg_m4l/F");
  SelTree.Branch("Z1Mass",&myZ1Mass,"Z1Mass/F");
  SelTree.Branch("Z1Pt",&myZ1Pt,"Z1Pt/F");
  SelTree.Branch("Z1ids",&myZ1ids,"Z1ids/S");
  SelTree.Branch("Z2Mass",&myZ2Mass,"Z2Mass/F");
  SelTree.Branch("Z2Pt",&myZ2Pt,"Z2Pt/F");
  SelTree.Branch("Z2ids",&myZ2ids,"Z2ids/S");
  SelTree.Branch("MC_weight",&MC_weight,"MC_weight/F");
  SelTree.Branch("MC_weight_norm",&MC_weight_norm,"MC_weight_norm/F");
  SelTree.Branch("MC_weight_noxsec",&MC_weight_noxsec,"MC_weight_noxsec/F");
  SelTree.Branch("MC_weight_PUWeight",&PUWeight,"MC_weight_PUWeight/F");
  SelTree.Branch("MC_weight_powhegWeight",&powheg_weight,"MC_weight_powhegWeight/F");
  SelTree.Branch("MC_weight_dataMC",&dataMCWeight,"MC_weight_dataMC/F");
  SelTree.Branch("MC_weight_HqT",&HqTWeight,"MC_weight_HqT/F");
  SelTree.Branch("costhetastar",&mycosthetastar,"costhetastar/F");
  SelTree.Branch("helphi",&myhelphi,"helphi/F");
  SelTree.Branch("helcosthetaZ1",&myhelcosthetaZ1,"helcosthetaZ1/F");
  SelTree.Branch("helcosthetaZ2",&myhelcosthetaZ2,"helcosthetaZ2/F");
  SelTree.Branch("phistarZ1",&myphistarZ1,"phistarZ1/F");
  SelTree.Branch("phistarZ2",&myphistarZ2,"phistarZ2/F");
  SelTree.Branch("xi",&myxi,"xi/F");
  SelTree.Branch("xistar",&myxistar,"xistar/F");
  SelTree.Branch("pg1g4_mela",&mypg1g4_mela,"pg1g4_mela/F");
  SelTree.Branch("pg1g4_VAJHU",&mypg1g4_VAJHU,"pg1g4_VAJHU/F");
  SelTree.Branch("pg1g4_pi2_VAJHU",&mypg1g4_pi2_VAJHU,"pg1g4_pi2_VAJHU/F");
  SelTree.Branch("pg1g2_pi2_VAJHU",&mypg1g2_pi2_VAJHU,"pg1g2_pi2_VAJHU/F");
  SelTree.Branch("pg1g2_mela",&mypg1g2_mela,"pg1g2_mela/F");
  SelTree.Branch("pg1g2_VAJHU",&mypg1g2_VAJHU,"pg1g2_VAJHU/F");

  SelTree.Branch("pzzzg_VAJHU",    &mypzzzg_VAJHU,    "pzzzg_VAJHU/F");
  SelTree.Branch("pzzgg_VAJHU",    &mypzzgg_VAJHU,    "pzzgg_VAJHU/F");
  SelTree.Branch("pzzzg_PS_VAJHU", &mypzzzg_PS_VAJHU, "pzzzg_PS_VAJHU/F");
  SelTree.Branch("pzzgg_PS_VAJHU", &mypzzgg_PS_VAJHU, "pzzgg_PS_VAJHU/F");
  SelTree.Branch("p0Zgs_VAJHU",    &myp0Zgs_VAJHU,     "p0Zgs_VAJHU/F");
  SelTree.Branch("p0gsgs_VAJHU",   &myp0gsgs_VAJHU,   "p0gsgs_VAJHU/F");
  SelTree.Branch("p0Zgs_PS_VAJHU", &myp0Zgs_PS_VAJHU, "p0Zgs_PS_VAJHU/F");
  SelTree.Branch("p0gsgs_PS_VAJHU",&myp0gsgs_PS_VAJHU,"p0gsgs_PS_VAJHU/F");
  SelTree.Branch("genProcessId",&mygenProcessId,"genProcessId/S");
  SelTree.Branch("genHEPMCweight",&mygenHEPMCweight,"genHEPMCweight/F");
  SelTree.Branch("GenHPt",&mygenhpt,"GenHPt/F");
  SelTree.Branch("GenHMass",&mygenhmass,"GenHMass/F");
  SelTree.Branch("p0plus_m4l_ScaleUp",&myp0plus_m4l_ScaleUp,"p0plus_m4l_ScaleUp/F");
  SelTree.Branch("p0plus_m4l_ScaleDown",&myp0plus_m4l_ScaleDown,"p0plus_m4l_ScaleDown/F");
  SelTree.Branch("p0plus_m4l_ResUp",&myp0plus_m4l_ResUp,"p0plus_m4l_ResUp/F");
  SelTree.Branch("p0plus_m4l_ResDown",&myp0plus_m4l_ResDown,"p0plus_m4l_ResDown/F");
  SelTree.Branch("bkg_m4l_ScaleUp",&mybkg_m4l_ScaleUp,"bkg_m4l_ScaleUp/F");
  SelTree.Branch("bkg_m4l_ScaleDown",&mybkg_m4l_ScaleDown,"bkg_m4l_ScaleDown/F");
  SelTree.Branch("bkg_m4l_ResUp",&mybkg_m4l_ResUp,"bkg_m4l_ResUp/F");
  SelTree.Branch("bkg_m4l_ResDown",&mybkg_m4l_ResDown,"bkg_m4l_ResDown/F");

  if(saveJets){
    SelTree.Branch("DiJetMass",&myDiJetMass,"DiJetMass/F");
    SelTree.Branch("DiJetMassPlus",&myDiJetMassPlus,"DiJetMassPlus/F");
    SelTree.Branch("DiJetMassMinus",&myDiJetMassMinus,"DiJetMassMinus/F");
    SelTree.Branch("DiJetDEta",&myDiJetDEta,"DiJetDEta/F");
    SelTree.Branch("NJets30",&myNJets30,"NJets30/S");
    SelTree.Branch("JetPt",&myJetPt);
    SelTree.Branch("JetSigma",&myJetSigma);
    SelTree.Branch("JetEta",&myJetEta);
    SelTree.Branch("JetPhi",&myJetPhi);
    SelTree.Branch("JetMass",&myJetMass);
    SelTree.Branch("JetBTag",&myJetBTag);
  }  

  if(isCR){
    // 	SelTree.Branch("Lep1Pt",&myLep1Pt,"Lep1Pt/F");
    // 	SelTree.Branch("Lep1Eta",&myLep1Eta,"Lep1Eta/F");
    // 	SelTree.Branch("Lep1ID",&myLep1ID,"Lep1ID/I");
    // 	SelTree.Branch("Lep2Pt",&myLep2Pt,"Lep2Pt/F");
    // 	SelTree.Branch("Lep2Eta",&myLep2Eta,"Lep2Eta/F");
    // 	SelTree.Branch("Lep2ID",&myLep2ID,"Lep2ID/I");
    // 	SelTree.Branch("Lep3Pt",&myLep3Pt,"Lep3Pt/F");
    // 	SelTree.Branch("Lep3Eta",&myLep3Eta,"Lep3Eta/F");
    // 	SelTree.Branch("Lep3ID",&myLep3ID,"Lep3ID/I");
    // 	SelTree.Branch("Lep4Pt",&myLep4Pt,"Lep4Pt/F");
    // 	SelTree.Branch("Lep4Eta",&myLep4Eta,"Lep4Eta/F");
    // 	SelTree.Branch("Lep4ID",&myLep4ID,"Lep4ID/I");
    SelTree.Branch("CRflag",&myCRflag,"CRflag/I");
    SelTree.Branch("ZXfake_weight",&ZXfake_weight,"ZXfake_weight/F");
  }

  //Start looping over the events
  Long64_t nb = 0;
  for (Long64_t jentry=0;jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    //    static const string filestring = fChain->GetCurrentFile()->GetName();
    nTOTEv++;
    // MC weighting, including data/MC corrections and PU
    if(!isData){
      MC_weight_interm = MC_weight_initial;
      MC_weight_norm_interm = MC_weight_norm_initial;
      MC_weight_noxsec_interm = MC_weight_noxsec_initial;

      // FIXME: at the moment we use the gen weight only for minlo, cf. normalization to Nevt_Gen_weighted above.
      // Also, only for MC_weight and MC_weight_noxsec, as we do not have the counters for weighted events per final state to be used for MC_weight_norm.
      if (isMinlo) {
	MC_weight_interm*=genHEPMCweight;
	MC_weight_noxsec_interm*=genHEPMCweight;
      }

      if(!is8TeV) PUWeight = PUWeighter.weight(2011,2011,NTrueInt);
      else if(is8TeV) PUWeight = PUWeighter.weight(2012,2012,NTrueInt);
      else abort();
      
      if (fabs(PUWeight-PUWeight12)>0.001){
	cout << "ERROR: PUweights: " << PUWeight << " " << PUWeight12 << " " << is8TeV << " " << NTrueInt << endl;
	//abort();
      }

      MC_weight_interm *= PUWeight;
      MC_weight_norm_interm *= PUWeight;
      MC_weight_noxsec_interm *= PUWeight;
      
      //A couple of validation plots
      if (genFinalState==channelType) { // only for "proper" events
	hNvtxNoWeight.Fill(Nvtx);
	hNvtxWeight.Fill(Nvtx,PUWeight);
	hPUWeight.Fill(PUWeight);
      }
      mygenProcessId = genProcessId;
      mygenHEPMCweight = genHEPMCweight;
      mygenhpt = GenHPt;
      mygenhmass = GenHMass;
    }

    //Number of Higgs candidates in the event
    const Int_t NHiggs = ZZMass->size();

    bool foundBestCand(false);

    //change the lineshape at high mass (Giampi's corrections)
    if(isSignal) {
      // High-mass weights, only above 400 GeV
      if(mPOLE > 399. && doHighmassCorr && !isMinlo){
	powheg_weight = getPwhgWeight(GenHMass, pwhg_flag);
	    
	//cout<<"giampi weight "<<powheg_weight <<endl;

	MC_weight_interm *= powheg_weight;
	MC_weight_norm_interm *= powheg_weight;
	MC_weight_noxsec_interm *= powheg_weight;

	if (doHqTCorr){
	  // HqT weight
	  HqTWeight = getHqTWeight(GenHMass,GenHPt,fHqT);
	  //cout<<"HqT weight = "<<HqTWeight<<endl;
	    
	  MC_weight_interm *= HqTWeight;
	  MC_weight_norm_interm *= HqTWeight;
	  MC_weight_noxsec_interm *= HqTWeight;
	}
      }
      
      if (genFinalState==channelType) { // only for "proper" events
	hZZmassNoPoweg.Fill(GenHMass);  // Gen H mass, not reweighted	
	hZZmassPoweg.Fill(GenHMass,powheg_weight*HqTWeight);
      }
    }

    //fill the plot with the highest generated pt for signal
    if(isSignal && genFinalState==channelType){ // only for "proper" events
      Float_t muptmax = -1.;
      if(GenLep1Pt > muptmax) muptmax = GenLep1Pt;
      if(GenLep2Pt > muptmax) muptmax = GenLep2Pt;
      if(GenLep3Pt > muptmax) muptmax = GenLep3Pt;
      if(GenLep4Pt > muptmax) muptmax = GenLep4Pt;
      hLepGenPt.Fill(muptmax);
    }
      // Jet information: at present the list of jets is relative to the best candidate only
      // It is empty in the CR, where only "ZZFisher" is saved!
      // Note that the list includes jets below 30, for JES syst studies
      if(saveJets){
	//Clear the Jet collection to be saved
	myJetPt.clear();
	myJetSigma.clear();
	myJetEta.clear();
	myJetPhi.clear();
	myJetMass.clear();
	myJetBTag.clear();
	int countJet = 0;
	for(int i=0;i<JetPt->size();i++){
	  myJetPt.push_back(JetPt->at(i));
	  myJetSigma.push_back(JetSigma->at(i));
	  myJetEta.push_back(JetEta->at(i));
	  myJetPhi.push_back(JetPhi->at(i));
	  myJetMass.push_back(JetMass->at(i));
	  myJetBTag.push_back(JetBTag->at(i));
	  if (JetPt->at(i)>30.) countJet++; // We want Fisher only for 2 jets > 30; jets below 30 are kept for JES syst studies
	}
	myDiJetMass = DiJetMass; //if <2 jets, DiJetMass=DiJetDEta=-99;
	myDiJetMassPlus = DiJetMassPlus;
	myDiJetMassMinus = DiJetMassMinus;
	myDiJetDEta = DiJetDEta;
	myNJets30 = countJet;
      }

    //Loop over all the H-> ZZ candidates in the event
    for(int nH=0; nH<NHiggs;nH++){

      //Do the final selection and best candidate selection. Depends on signal region or control region
      if(isCR && !CRflag->at(nH)) continue; // Belongs to at least 1 CR
//      if(!isCR && (ZZsel->at(nH) < 100. || nH != iBC) ) continue; // Only save events passing the FULL selection 
      if(!isCR && (ZZsel->at(nH) < 90 || nH != iBC) ) continue;  // Only 

      Int_t RunFraction = -1;
      Float_t whatPeriod = RooRandom::randomGenerator()->Uniform();
      if(!is8TeV && whatPeriod < Run2011AFraction) RunFraction = 0;
      else if(!is8TeV && whatPeriod > Run2011AFraction) RunFraction = 1;
      else if(is8TeV) RunFraction = 2;
      else abort();

      if(!isData){
	dataMCWeight = getMCWeight(RunFraction, nH);
        MC_weight = MC_weight_interm * dataMCWeight;
        MC_weight_norm = MC_weight_norm_interm * dataMCWeight;
        MC_weight_noxsec = MC_weight_noxsec_interm * dataMCWeight;
	hTagProbeWeight.Fill(dataMCWeight);
      }
      
      // Z+X Fake Weight
      if(isCR) ZXfake_weight = getZXfake_weight(is8TeV, nH); 
      //cout << "ZXfakeweight = " << ZXfake_weight << endl;
      //cout << RunNumber << " " << LumiNumber << " " << EventNumber << " " << ZXfake_weight << endl;


      nBestCand += MC_weight;
      nSel++;
      foundBestCand = true;

      //cout << getMCWeight(RunFraction, nH) << endl;

      hZZmassPowegReco.Fill(ZZMass->at(nH),MC_weight); 

      if(isSignal && !isVBF && !isData && !isCR && ( (applySystMu && channelType != 1) || (applySystEle && channelType != 0) ) ){
	
	for(int i=0; i<N_syst_Iter; i++){
	  //No need of other normalizations, as they wipe in the ratio
	  Seeder = i;
	  nBestSyst[i] += getMCWeight(RunFraction, nH);
	}
      }
      
      TLorentzVector m[4];
      // FIXME: these are without FSR, do not use to build the ZZ kinematics!!!
//       m[0].SetPtEtaPhiM(Lep1Pt->at(nH),Lep1Eta->at(nH),Lep1Phi->at(nH),0.1506583);
//       m[1].SetPtEtaPhiM(Lep2Pt->at(nH),Lep2Eta->at(nH),Lep2Phi->at(nH),0.1506583);
//       m[2].SetPtEtaPhiM(Lep3Pt->at(nH),Lep3Eta->at(nH),Lep3Phi->at(nH),0.1506583);
//       m[3].SetPtEtaPhiM(Lep4Pt->at(nH),Lep4Eta->at(nH),Lep4Phi->at(nH),0.1506583);      
//      TLorentzVector momH = m[0] + m[1] + m[2] + m[3]; // Withouth FSR, do not use

      //Fill the variables that go into the final tree
      myZZMass = ZZMass->at(nH);
      myZZMassErr = ZZMassErr->at(nH);
      myZZMassErrCorr = ZZMassErrCorr->at(nH);
      myCRflag = (int)CRflag->at(nH);
      myZZPt = ZZPt->at(nH);
      myZZEta = ZZEta->at(nH);
      myZZPhi = ZZPhi->at(nH);
      //      myZZRapidity = momH.Rapidity();
      myZZFisher = ZZFisher->at(nH);
      myp0plus_VAJHU = p0plus_VAJHU->at(nH);
      myp0hplus_VAJHU = p0hplus_VAJHU->at(nH);
      myp0minus_VAJHU = p0minus_VAJHU->at(nH);
      myp1_VAJHU = p1_VAJHU->at(nH);
      myp1_prodIndep_VAJHU = p1_prodIndep_VAJHU->at(nH);
      myp1plus_VAJHU = p1plus_VAJHU->at(nH);
      myp1plus_prodIndep_VAJHU = p1plus_prodIndep_VAJHU->at(nH);
      myp2_prodIndep_VAJHU = p2_prodIndep_VAJHU->at(nH);
      myp2_VAJHU = p2_VAJHU->at(nH);
      myp2qqb_VAJHU = p2qqb_VAJHU->at(nH);
      myp2hplus_VAJHU = p2hplus_VAJHU->at(nH);
      myp2hminus_VAJHU = p2hminus_VAJHU->at(nH);
      myp2bplus_VAJHU = p2bplus_VAJHU->at(nH);

 			myp2hplus_qqb_VAJHU =          p2hplus_qqb_VAJHU ->at(nH);           
 			myp2hplus_prodIndep_VAJHU =    p2hplus_prodIndep_VAJHU ->at(nH);     
 			myp2hminus_qqb_VAJHU =         p2hminus_qqb_VAJHU ->at(nH);          
 			myp2hminus_prodIndep_VAJHU =   p2hminus_prodIndep_VAJHU ->at(nH);    
 			myp2bplus_qqb_VAJHU =          p2bplus_qqb_VAJHU ->at(nH);           
 			myp2bplus_prodIndep_VAJHU =    p2bplus_prodIndep_VAJHU ->at(nH);     
 			myp2h2plus_gg_VAJHU =          p2h2plus_gg_VAJHU ->at(nH);           
 			myp2h2plus_qqbar_VAJHU =       p2h2plus_qqbar_VAJHU ->at(nH);        
 			myp2h2plus_prodIndep_VAJHU =   p2h2plus_prodIndep_VAJHU ->at(nH);    
 			myp2h3plus_gg_VAJHU =          p2h3plus_gg_VAJHU ->at(nH);           
 			myp2h3plus_qqbar_VAJHU =       p2h3plus_qqbar_VAJHU ->at(nH);        
 			myp2h3plus_prodIndep_VAJHU =   p2h3plus_prodIndep_VAJHU ->at(nH);    
 			myp2h6plus_gg_VAJHU =          p2h6plus_gg_VAJHU ->at(nH);           
 			myp2h6plus_qqbar_VAJHU =       p2h6plus_qqbar_VAJHU ->at(nH);        
 			myp2h6plus_prodIndep_VAJHU =   p2h6plus_prodIndep_VAJHU ->at(nH);    
 			myp2h7plus_gg_VAJHU =          p2h7plus_gg_VAJHU ->at(nH);           
 			myp2h7plus_qqbar_VAJHU =       p2h7plus_qqbar_VAJHU ->at(nH);        
 			myp2h7plus_prodIndep_VAJHU =   p2h7plus_prodIndep_VAJHU ->at(nH);    
 			myp2h9minus_gg_VAJHU =         p2h9minus_gg_VAJHU ->at(nH);          
 			myp2h9minus_qqbar_VAJHU =      p2h9minus_qqbar_VAJHU ->at(nH);       
 			myp2h9minus_prodIndep_VAJHU =  p2h9minus_prodIndep_VAJHU ->at(nH);   
 			myp2h10minus_gg_VAJHU =        p2h10minus_gg_VAJHU ->at(nH);         
 			myp2h10minus_qqbar_VAJHU =     p2h10minus_qqbar_VAJHU ->at(nH);      
 			myp2h10minus_prodIndep_VAJHU = p2h10minus_prodIndep_VAJHU ->at(nH);  
      mybkg_VAMCFM = bkg_VAMCFM->at(nH);
      mybkg_prodIndep_VAMCFM = bkg_prodIndep_VAMCFM->at(nH);
      myggzz_VAMCFM = ggzz_VAMCFM->at(nH);
      myggzz_p0plus_VAMCFM = ggzz_p0plus_VAMCFM->at(nH);
      myp0plus_VAMCFM = p0plus_VAMCFM->at(nH);
      myggzz_c1_VAMCFM = ggzz_c1_VAMCFM->at(nH);
      myggzz_c5_VAMCFM = ggzz_c5_VAMCFM->at(nH);
      myggzz_ci_VAMCFM = ggzz_ci_VAMCFM->at(nH);
      myphjj_VAJHU_old = phjj_VAJHU_old->at(nH);
      mypvbf_VAJHU_old = pvbf_VAJHU_old->at(nH);
      myphjj_VAJHU_old_up = phjj_VAJHU_old_up->at(nH);
      mypvbf_VAJHU_old_up = pvbf_VAJHU_old_up->at(nH);
      myphjj_VAJHU_old_dn = phjj_VAJHU_old_dn->at(nH);
      mypvbf_VAJHU_old_dn = pvbf_VAJHU_old_dn->at(nH);
      myphjj_VAJHU_new = phjj_VAJHU_new->at(nH);
      mypvbf_VAJHU_new = pvbf_VAJHU_new->at(nH);
      myphjj_VAJHU_new_up = phjj_VAJHU_new_up->at(nH);
      mypvbf_VAJHU_new_up = pvbf_VAJHU_new_up->at(nH);
      myphjj_VAJHU_new_dn = phjj_VAJHU_new_dn->at(nH);
      mypvbf_VAJHU_new_dn = pvbf_VAJHU_new_dn->at(nH);
	    myp0_g1prime2_VAJHU = p0_g1prime2_VAJHU->at(nH);
	    mypg1g1prime2_VAJHU = pg1g1prime2_VAJHU->at(nH);
	    myDgg10_VAMCFM = Dgg10_VAMCFM->at(nH);
      myp0plus_m4l = p0plus_m4l->at(nH);
      mybkg_m4l = bkg_m4l->at(nH);
      myZ1Mass = Z1Mass->at(nH);
      myZ1Pt = Z1Pt->at(nH);
      myZ2Mass = Z2Mass->at(nH);
      myZ2Pt = Z2Pt->at(nH);
      mycosthetastar = costhetastar->at(nH);
      myhelphi = helphi->at(nH);
      myhelcosthetaZ1 = helcosthetaZ1->at(nH);
      myhelcosthetaZ2 = helcosthetaZ2->at(nH);
      myphistarZ1 = phistarZ1->at(nH);
      myphistarZ2 = phistarZ2->at(nH);
      myxi = xi->at(nH);
      myxistar = xistar->at(nH);
      mypg1g4_mela =     pg1g4_mela->at(nH);      
      mypg1g4_VAJHU=     pg1g4_VAJHU->at(nH);     
      mypg1g4_pi2_VAJHU= pg1g4_pi2_VAJHU->at(nH); 
      mypg1g2_pi2_VAJHU= pg1g2_pi2_VAJHU->at(nH); 
      mypg1g2_mela=      pg1g2_mela->at(nH);      
      mypg1g2_VAJHU=     pg1g2_VAJHU->at(nH);     
  		mypzzzg_VAJHU=    pzzzg_VAJHU->at(nH);
  		mypzzgg_VAJHU=    pzzgg_VAJHU->at(nH);
  		mypzzzg_PS_VAJHU= pzzzg_PS_VAJHU->at(nH);
  		mypzzgg_PS_VAJHU= pzzgg_PS_VAJHU->at(nH);
  		myp0Zgs_VAJHU=    p0Zgs_VAJHU->at(nH);
  		myp0gsgs_VAJHU=   p0gsgs_VAJHU->at(nH);
  		myp0Zgs_PS_VAJHU= p0Zgs_PS_VAJHU->at(nH);
  		myp0gsgs_PS_VAJHU=p0gsgs_PS_VAJHU->at(nH);

      myp0plus_m4l_ScaleUp = p0plus_m4l_ScaleUp->at(nH);
      myp0plus_m4l_ScaleDown = p0plus_m4l_ScaleDown->at(nH);
      myp0plus_m4l_ResUp = p0plus_m4l_ResUp->at(nH);
      myp0plus_m4l_ResDown = p0plus_m4l_ResDown->at(nH);
      mybkg_m4l_ScaleUp = bkg_m4l_ScaleUp->at(nH);
      mybkg_m4l_ScaleDown = bkg_m4l_ScaleDown->at(nH);
      mybkg_m4l_ResUp = bkg_m4l_ResUp->at(nH);
      mybkg_m4l_ResDown = bkg_m4l_ResDown->at(nH);

      if(isCR){
	// for CR, the jet list is empty; we "fake" the NJets30 variable for consistency
	if (myZZFisher>=0) myNJets30=2;
	else myNJets30=0;
      }
      myLep1Pt  = Lep1Pt->at(nH);
      myLep1Eta = Lep1Eta->at(nH);
      myLep1ID = Lep1LepId->at(nH);
      myLep2Pt  = Lep2Pt->at(nH);
      myLep2Eta = Lep2Eta->at(nH);
      myLep2ID = Lep2LepId->at(nH);
      myLep3Pt  = Lep3Pt->at(nH);
      myLep3Eta = Lep3Eta->at(nH);
      myLep3ID = Lep3LepId->at(nH);
      myLep4Pt  = Lep4Pt->at(nH);
      myLep4Eta = Lep4Eta->at(nH);
      myLep4ID = Lep4LepId->at(nH);
      myZ1ids = myLep1ID*myLep2ID;
      myZ2ids = myLep3ID*myLep4ID;

      myEventNumber = EventNumber; 
      myRunNumber = RunNumber;       
      myLumiNumber = LumiNumber;

      if(isCR) SelTree.Fill();
      else break; //found best cand, no need to continue

    }//end of loop over the candidates

    if(foundBestCand && !isCR) SelTree.Fill();

  }//end of loop over the events

  if (fHqT!=0) fHqT->Close();


  // Check the effect of the weights on the normalization
  float normDiffPUWeight = 0;
  if (hNvtxNoWeight.Integral()>0) {
    normDiffPUWeight = (hNvtxWeight.Integral()-hNvtxNoWeight.Integral())/hNvtxNoWeight.Integral();
    if ( fabs(normDiffPUWeight)>0.01 ) { //1%
      cout<<"++++++++"<<endl;
      cout<<"++++ WARNING: The normalization after PU weights changes by "<<setprecision(2)<<normDiffPUWeight*100<<" %"<<endl;
    }
  }
  

  float normDiffPowhegWeight = 0;
  if (hZZmassNoPoweg.Integral()>0) {
    normDiffPowhegWeight = (hZZmassPoweg.Integral()-hZZmassNoPoweg.Integral())/hZZmassNoPoweg.Integral();
    if ( fabs(normDiffPowhegWeight)>0.008 ) { //0.8%
      cout<<"++++++++"<<endl;
      cout<<"++++ WARNING: The normalization after POWHEG weight changes by "<<setprecision(2)<<normDiffPowhegWeight*100<<" %"<<endl;
    }
  }

  //Save the final tree
  fOut.cd();
  SelTree.Write();
  hNvtxNoWeight.Write();
  hNvtxWeight.Write();
  hLepGenPt.Write();
  hPUWeight.Write();
  hTagProbeWeight.Write();
  hZZmassNoPoweg.Write();
  hZZmassPoweg.Write();
  hZZmassPowegReco.Write();
  nEventComplete->Write("hCounters");
  fOut.Close();

  //Save the number of events for the systematic on normalization
  if(isSignal && !isVBF && !isData && !isCR && (applySystMu || applySystEle)){

    TH1F h_aux_syst("h_aux_syst","h_aux_syst",1000,0.,10.);
    for(int i=0; i<N_syst_Iter; i++) h_aux_syst.Fill(nBestSyst[i]/nSel);
    Float_t corr_syst = h_aux_syst.GetRMS();

    string outfileName;
    if(channelType == 0) outfileName = "Nevt_4mu.txt";
    if(channelType == 1) outfileName = "Nevt_4e.txt";
    if(channelType == 2) outfileName = "Nevt_2mu2e.txt";
    ofstream ofsN(outfileName.c_str(),fstream::out | fstream::app);
    ofsN << corr_syst << endl;
  }

  cout << "Total number of events analyzed :       " << nTOTEv << endl;
  cout << "Total number of events selected :       " << nSel << endl;
  cout << "Total number of best candidates found : " << nBestCand << endl;
  cout << "################################################" << endl;

  return;
}


Int_t HZZ4l::findBestCRCand(int CR) const
{
  Int_t bestCand   = -1;
  Float_t bestMMPt = -1.;

  for(int nH=0; nH<CRflag->size();nH++){
    if (test_bit(CRflag->at(nH),CR)) {
      bestCand=nH;
      break;
    }
  }

  return bestCand;
}

Float_t HZZ4l::getAllWeight(const Float_t LepPt, const Float_t LepEta, const Int_t year, Int_t LepID) const
{
  Float_t weight  = 1.; 
  Float_t errCorr = 0.;
  Float_t errCorrSyst = 0.;

  Float_t myLepPt = LepPt;
  Float_t myLepEta = LepEta;
  Int_t   myLepID = abs(LepID);
  
  //avoid to go out of the TH boundary
  if(myLepID == 13 && myLepPt > 99.) myLepPt = 99.;
  if(myLepID == 11 && myLepPt > 199.) myLepPt = 199.;
  if(myLepID == 11) myLepEta = fabs(myLepEta);

  //Scale factors for data/MC efficiency
  static TFile fMuWeight("scale_factors_muons2011.root");
  static TFile fMuWeight12("scale_factors_muons2012.root");
  static TFile fElWeight("scale_factors_ele2011.root");
  static TFile fElWeight12("scale_factors_ele2012.root");
 
  static TH2D *hTH2D_Mu_All_2011A = (TH2D*)fMuWeight.Get("TH2D_ALL_2011A"); 
  static TH2D *hTH2D_Mu_All_2011B = (TH2D*)fMuWeight.Get("TH2D_ALL_2011B"); 
  static TH2D *hTH2D_Mu_All_2012  = (TH2D*)fMuWeight12.Get("TH2D_ALL_2012"); 
  
  TString eleSFname="h_electronScaleFactor_RecoIdIsoSip";
  static TH2D *hTH2D_El_All_2011A = (TH2D*)fElWeight.Get(eleSFname.Data());
  static TH2D *hTH2D_El_All_2011B = (TH2D*)fElWeight.Get(eleSFname.Data()); 
  static TH2D *hTH2D_El_All_2012  = (TH2D*)fElWeight12.Get(eleSFname.Data());

  if(year == 0){
    if(myLepID == 13){                                               
      weight  = hTH2D_Mu_All_2011A->GetBinContent(hTH2D_Mu_All_2011A->GetXaxis()->FindBin(myLepPt),hTH2D_Mu_All_2011A->GetYaxis()->FindBin(LepEta));
      errCorr = hTH2D_Mu_All_2011A->GetBinError(hTH2D_Mu_All_2011A->GetXaxis()->FindBin(myLepPt),hTH2D_Mu_All_2011A->GetYaxis()->FindBin(LepEta));
      
    }
    
    else if(myLepID == 11){   

      weight  = hTH2D_El_All_2011A->GetBinContent(hTH2D_El_All_2011A->GetXaxis()->FindBin(myLepPt),hTH2D_El_All_2011A->GetYaxis()->FindBin(myLepEta));
      errCorr = hTH2D_El_All_2011A->GetBinError(hTH2D_El_All_2011A->GetXaxis()->FindBin(myLepPt),hTH2D_El_All_2011A->GetYaxis()->FindBin(myLepEta));   
    }
    else {
      abort();
    }    
  }
  else if(year == 1){
    if( myLepID == 13){                                               
      weight  = hTH2D_Mu_All_2011B->GetBinContent(hTH2D_Mu_All_2011B->GetXaxis()->FindBin(myLepPt),hTH2D_Mu_All_2011B->GetYaxis()->FindBin(LepEta));
      errCorr = hTH2D_Mu_All_2011B->GetBinError(hTH2D_Mu_All_2011B->GetXaxis()->FindBin(myLepPt),hTH2D_Mu_All_2011B->GetYaxis()->FindBin(LepEta));
      
    }
    
    else  if(myLepID == 11){   

      weight  = hTH2D_El_All_2011B->GetBinContent(hTH2D_El_All_2011B->GetXaxis()->FindBin(myLepPt),hTH2D_El_All_2011B->GetYaxis()->FindBin(myLepEta));
      errCorr = hTH2D_El_All_2011B->GetBinError(hTH2D_El_All_2011B->GetXaxis()->FindBin(myLepPt),hTH2D_El_All_2011B->GetYaxis()->FindBin(myLepEta));
    }
    
    else {
      abort();
    }    
  }
  else if(year == 2){
    if( myLepID == 13){
      weight  = hTH2D_Mu_All_2012->GetBinContent(hTH2D_Mu_All_2012->GetXaxis()->FindBin(myLepPt),hTH2D_Mu_All_2012->GetYaxis()->FindBin(LepEta));
      errCorr = hTH2D_Mu_All_2012->GetBinError(hTH2D_Mu_All_2012->GetXaxis()->FindBin(myLepPt),hTH2D_Mu_All_2012->GetYaxis()->FindBin(LepEta));
    }
    
    else if(myLepID == 11){
      weight  = hTH2D_El_All_2012->GetBinContent(hTH2D_El_All_2012->GetXaxis()->FindBin(myLepPt),hTH2D_El_All_2012->GetYaxis()->FindBin(myLepEta));
      errCorr = hTH2D_El_All_2012->GetBinError(hTH2D_El_All_2012->GetXaxis()->FindBin(myLepPt),hTH2D_El_All_2012->GetYaxis()->FindBin(myLepEta));
    }
    else {
      abort();  
    }
  }
  else abort();
  
  //add the systematics on T&P corrections (for muons only, electrons have them already included)
  if(myLepID == 13){
    if(myLepPt >= 15.) errCorrSyst = 0.005;
    else errCorrSyst = 0.015;
  }

  //FIXME
  if(myLepPt < 5. && myLepID == 13) weight = 1.;

  if(weight < 0.001 || weight > 10.){
    cout << "myLepPt = " << myLepPt << " myLepEta = " << myLepEta << " weight = " << weight << endl;
    abort();  //no correction should be zero, if you find one, stop
  }

  static TRandom3 randomToss;
  if( (applySystMu && myLepID == 13) || (applySystEle && myLepID == 11) ){

    //apply correlation matrix by assigning the proper seed
    Int_t CorrSeeder = Seeder;
    if(myLepID == 13){
      if(myLepPt < 20. && fabs(myLepEta) < 1.2) CorrSeeder += 100001;
      else if(myLepPt < 20. && fabs(myLepEta) >= 1.2) CorrSeeder += 100002;
    }

    randomToss.SetSeed(CorrSeeder);
    weight = randomToss.Gaus(weight,errCorr);

    //apply systematic (totally correlated in eta) for muons
    if(myLepID == 13){
      randomToss.SetSeed(Seeder);
      weight = randomToss.Gaus(weight,errCorrSyst);
    }
  }

  return weight;
}

Float_t HZZ4l::getMCWeight(const int year, const Int_t CandIndex) const
{
  Float_t eff_weight = 1.;
  
  TLorentzVector m[4];
  m[0].SetPtEtaPhiM(Lep1Pt->at(CandIndex),Lep1Eta->at(CandIndex),Lep1Phi->at(CandIndex),0.1506583);
  m[1].SetPtEtaPhiM(Lep2Pt->at(CandIndex),Lep2Eta->at(CandIndex),Lep2Phi->at(CandIndex),0.1506583);
  m[2].SetPtEtaPhiM(Lep3Pt->at(CandIndex),Lep3Eta->at(CandIndex),Lep3Phi->at(CandIndex),0.1506583);
  m[3].SetPtEtaPhiM(Lep4Pt->at(CandIndex),Lep4Eta->at(CandIndex),Lep4Phi->at(CandIndex),0.1506583);

  Int_t LepId[4];
  LepId[0] = Lep1LepId->at(CandIndex);
  LepId[1] = Lep2LepId->at(CandIndex);
  LepId[2] = Lep3LepId->at(CandIndex);
  LepId[3] = Lep4LepId->at(CandIndex);

  for(int nLep=0;nLep<4;nLep++){

    if(fabs(LepId[nLep]) == 13 && doMuEffCorr) eff_weight *= getAllWeight(m[nLep].Pt(), m[nLep].Eta(), year, LepId[nLep]);

    if(fabs(LepId[nLep]) == 11 && doEleEffCorr) eff_weight *= getAllWeight(m[nLep].Pt(), m[nLep].Eta(), year, LepId[nLep]);

  }

  return eff_weight;
}

//Get the normalization for the specific final state (4mu, 4e, 2e2mu), that is used for MC_weight_norm (efficiency weights).
Float_t HZZ4l::getNormalizedWeight(const Int_t channelType) const
{
  // For signals, take directly the values from the Counters histogram.
  if (isSignal || theSample=="ZZJetsTo4L"){
    Int_t gen_4l_H = nEventComplete->GetBinContent(channelType + 2);
    return 1./gen_4l_H;
  } 

  /* The mechanism for efficiency weight for ZZ samples needs to be reviewed.
     In particular, the 4tau sample is not included in the table.
     These values are not used anyhow.

  else if (theSample.Contains("ZZ")) { // For ZZ samples, use cached values to account for tau contribution properly
    //Set the channel name
    string channelname;
    if(channelType == 0) channelname = "4mu";
    else if(channelType == 1) channelname = "4e";
    else if(channelType == 2) channelname = "2e2mu";
    else abort();

    TString tagName = channelname + "_" + theSample;

    ifstream inputFile;
    if(is8TeV) inputFile.open("saveNormValues_8TeV.txt");
    else inputFile.open("saveNormValues_7TeV.txt");
    
    map<string,Float_t> NormMap;
    char normName[50];
    Float_t normValue = 0.;

    while(inputFile >> normName >> normValue) NormMap[normName] = normValue;

    if (NormMap.find(tagName.Data())!=NormMap.end()){
      return (NormMap.find(tagName.Data()))->second;
    } else {
      cout << "Sample " << theSample << " not found in table of normalization values" << endl;
      return 0.;
    }
  }
  */
  // Not relevant for other samples
  return 0.;     
}


void HZZ4l::getWeightFromFile(const TString& mPOLE, const Bool_t isVBF, const Bool_t isNewHighmass)
{
  TString Energy = is8TeV ? "8TeV" : "7TeV";

  TString filename;
  if(!isVBF) filename = "WeightCards/mZZ_Higgs" + mPOLE + "_" + Energy + "_Lineshape+Interference.txt";
  else filename = "WeightCards/VBF/VBF_ratio_" + Energy + "_" + mPOLE + ".txt";

  if(BSM_flag == 1){
    int kappa_bin = kappa/0.2;
    switch(kappa_bin){

    case 1:
      filename = "WeightCards/EWKsinglet/kappa_02/mZZ_Higgs" + mPOLE + "_" + Energy + "_Lineshape+Interference.txt";

    case 2:
      filename = "WeightCards/EWKsinglet/kappa_04/mZZ_Higgs" + mPOLE + "_" + Energy + "_Lineshape+Interference.txt";

    case 3:
      filename = "WeightCards/EWKsinglet/kappa_06/mZZ_Higgs" + mPOLE + "_" + Energy + "_Lineshape+Interference.txt";

    case 4:
      filename = "WeightCards/EWKsinglet/kappa_08/mZZ_Higgs" + mPOLE + "_" + Energy + "_Lineshape+Interference.txt";

    case 5:
      filename = "WeightCards/EWKsinglet/kappa_10/mZZ_Higgs" + mPOLE + "_" + Energy + "_Lineshape+Interference.txt";
    }
  }

  cout << "File used to extract high mass weights: " << filename << endl;

  double bincenter, initial, cps, cps_int, cpsp_int, cpsm_int, cps_intp, cps_intm;

  vector<double> Vinitial;
  vector<double> Vcps_int;
  vector<double> Vcps;
  vector<double> Vcpsp_int;
  vector<double> Vcpsm_int;
  vector<double> Vcps_intp;
  vector<double> Vcps_intm;

  std::ifstream ifs(filename.Data());
  while( ifs.good() ) {
    ifs >> bincenter >> initial >> cps >> cps_int >> cpsp_int >> cpsm_int >> cps_intp >> cps_intm;

    pwhg_bincenters.push_back(bincenter);
    Vinitial.push_back(initial);
    Vcps.push_back(cps);
    Vcps_int.push_back(cps_int);
    Vcpsp_int.push_back(cpsp_int);
    Vcpsm_int.push_back(cpsm_int);
    Vcps_intp.push_back(cps_intp);
    Vcps_intm.push_back(cps_intm);
  }
 
  double Iinitial  = 0.;
  double Icps      = 0.;
  double Icps_int  = 0.;
  double Icpsp_int = 0.;
  double Icpsm_int = 0.;
  double Icps_intp = 0.;
  double Icps_intm = 0.;

  for(int i=0;i<Vinitial.size();i++){
    Iinitial  += Vinitial.at(i);
    Icps      += Vcps.at(i);
    Icps_int  += Vcps_int.at(i);
    Icpsp_int += Vcpsp_int.at(i);
    Icpsm_int += Vcpsm_int.at(i);
    Icps_intp += Vcps_intp.at(i);
    Icps_intm += Vcps_intm.at(i);
  }

  cout << "Integral of the high mass distributions = " << Iinitial << " / " << Icps << " / " << Icps_int << endl;

  const double cps_int_cps_corr = Icps/Icps_int;
  const double cpsp_int_cps_corr = Icps/Icpsp_int;
  const double cpsm_int_cps_corr = Icps/Icpsm_int;
  const double cps_intp_cps_corr = Icps/Icps_intp;
  const double cps_intm_cps_corr = Icps/Icps_intm;

  const double cps_int_initial_corr = Iinitial/Icps_int;
  const double cpsp_int_initial_corr = Iinitial/Icpsp_int;
  const double cpsm_int_initial_corr = Iinitial/Icpsm_int;
  const double cps_intp_initial_corr = Iinitial/Icps_intp;
  const double cps_intm_initial_corr = Iinitial/Icps_intm;

  for(int i=0;i<Vinitial.size();i++){

    // Old MC (not powheg15), both gg and VBF
    if(Vinitial.at(i) > 0. && !isNewHighmass){
      pwhg_weight.push_back(max(0.,cps_int_initial_corr*Vcps_int.at(i)/Vinitial.at(i)));
      pwhg_weightCPSP.push_back(max(0.,cpsp_int_initial_corr*Vcpsp_int.at(i)/Vinitial.at(i)));
      pwhg_weightCPSM.push_back(max(0.,cpsm_int_initial_corr*Vcpsm_int.at(i)/Vinitial.at(i)));
      pwhg_weightIntP.push_back(max(0.,cps_intp_initial_corr*Vcps_intp.at(i)/Vinitial.at(i)));
      pwhg_weightIntM.push_back(max(0.,cps_intm_initial_corr*Vcps_intm.at(i)/Vinitial.at(i)));
    }
    // powheg15, !VBF
    else if(Vcps.at(i) > 0. && !isVBF && isNewHighmass){
      pwhg_weight.push_back(max(0.,cps_int_cps_corr*Vcps_int.at(i)/Vcps.at(i)));
      pwhg_weightCPSP.push_back(max(0.,cpsp_int_cps_corr*Vcpsp_int.at(i)/Vcps.at(i)));
      pwhg_weightCPSM.push_back(max(0.,cpsm_int_cps_corr*Vcpsm_int.at(i)/Vcps.at(i)));
      pwhg_weightIntP.push_back(max(0.,cps_intp_cps_corr*Vcps_intp.at(i)/Vcps.at(i)));
      pwhg_weightIntM.push_back(max(0.,cps_intm_cps_corr*Vcps_intm.at(i)/Vcps.at(i)));
    } 
    // powheg15, VBF
    else if(isVBF && isNewHighmass){
      pwhg_weight.push_back(1.);
      pwhg_weightCPSP.push_back(1.);
      pwhg_weightCPSM.push_back(1.);
      pwhg_weightIntP.push_back(1.);
      pwhg_weightIntM.push_back(1.);      
    }
    
    else{//weights are not defined if initial distribution is 0 => set weight to 0
      pwhg_weight.push_back(0.);
      pwhg_weightCPSP.push_back(0.);
      pwhg_weightCPSM.push_back(0.);
      pwhg_weightIntP.push_back(0.);
      pwhg_weightIntM.push_back(0.);
    }

  }

  return;
}

Float_t HZZ4l::getPwhgWeight(double mH, const Int_t WhichSide) const
{
  //Weight to be returned and +- error
  Float_t weight     = 1.;
  Float_t weightCPSp = 1.;
  Float_t weightCPSm = 1.;
  Float_t weightIntp = 1.;
  Float_t weightIntm = 1.;

  if(mH < pwhg_bincenters.front() || mH >  pwhg_bincenters.back()) return 0.; // set weights to 0 if out of range

  vector<double>::const_iterator low = lower_bound(pwhg_bincenters.begin(), pwhg_bincenters.end(),mH); 
  int lowindex = (low -  pwhg_bincenters.begin());

  if(mH == *low){//exact match
    weight     = pwhg_weight[lowindex];
    weightCPSp = pwhg_weightCPSP[lowindex];
    weightCPSm = pwhg_weightCPSM[lowindex];
    weightIntp = pwhg_weightIntP[lowindex];
    weightIntm = pwhg_weightIntM[lowindex];
  }else{//linear interpolation
    lowindex--; // lower_bound finds the first element not smaller than X
    weight     = pwhg_weight[lowindex]     + (mH - pwhg_bincenters[lowindex]) * (pwhg_weight[lowindex+1] - pwhg_weight[lowindex])         / (pwhg_bincenters[lowindex+1] - pwhg_bincenters[lowindex]);
    weightCPSp = pwhg_weightCPSP[lowindex] + (mH - pwhg_bincenters[lowindex]) * (pwhg_weightCPSP[lowindex+1] - pwhg_weightCPSP[lowindex]) / (pwhg_bincenters[lowindex+1] - pwhg_bincenters[lowindex]);
    weightCPSm = pwhg_weightCPSM[lowindex] + (mH - pwhg_bincenters[lowindex]) * (pwhg_weightCPSM[lowindex+1] - pwhg_weightCPSM[lowindex]) / (pwhg_bincenters[lowindex+1] - pwhg_bincenters[lowindex]);
    weightIntp = pwhg_weightIntP[lowindex] + (mH - pwhg_bincenters[lowindex]) * (pwhg_weightIntP[lowindex+1] - pwhg_weightIntP[lowindex]) / (pwhg_bincenters[lowindex+1] - pwhg_bincenters[lowindex]);
    weightIntm = pwhg_weightIntM[lowindex] + (mH - pwhg_bincenters[lowindex]) * (pwhg_weightIntM[lowindex+1] - pwhg_weightIntM[lowindex]) / (pwhg_bincenters[lowindex+1] - pwhg_bincenters[lowindex]);
  }

  if(weight > 20. ){
    cout << "Removing event with very sick weight!!! " << mH << " " << *low << " " << lowindex << " " << pwhg_weight[lowindex] << " " << weight << endl;
    weight = 0.;
  }

  switch(WhichSide){

  case 0:
    return weight;

  case 1:
    return weightCPSm;

  case 2:
    return weightCPSp;

  case 3:
    return weightIntm;

  case 4:
    return weightIntp;

  default:
    abort();

  }

  return -999.;
}


Float_t HZZ4l::getHqTWeight(double mH, double genPt, TFile* f) const
{

  //cout<<"mH = "<<mH<<", genPt = "<<genPt<<endl;
  if (mH<400 || genPt>250) return 1.;
  
  double weight = 1.;
  
  const int masses[4] = {400,600,800,1000};
  double massDiff = 1000;
  int iMass = -1;
  for (int i=0; i<4; ++i){
    double massDiffTmp = std::fabs(mH-masses[i]);
    if (massDiffTmp<massDiff){
      massDiff = massDiffTmp;
      iMass = i;
    }
  }
  
  if (iMass<0) return weight;
  else{
    TH1D* h_weight = (TH1D*)f->Get(Form("wH_%d",masses[iMass]));
    weight = h_weight->GetBinContent(h_weight->FindBin(genPt));
    if (weight!=0) return weight;
  }
}



// Added by CO
Float_t HZZ4l::getFakeWeight(const Float_t LepPt, const Float_t LepEta, const Int_t year, Int_t LepID, Int_t LepZ1ID)
{
  // year 0 = 2011
  // year 1 = 2012

  Float_t weight  = 1.; 
  
  Float_t myLepPt   = LepPt;
  Float_t myLepEta  = fabs(LepEta);
  Int_t   myLepID   = abs(LepID);
  Int_t   myZ1LepID = abs(LepZ1ID);

  //cout << " pt = " << myLepPt << " eta = " << myLepEta << " ZID = " << myZ1LepID << " LepID = " << myLepID << endl;

  //avoid to go out of the TH boundary
  if(myLepPt > 79.) myLepPt = 79.;

  int n_file = myLepID-11 + year;

  TFile*  FileZXWeight = ZXWeightTables[n_file];
  if (FileZXWeight ==0) { // init
    // Fake Rate File 
    TString file_name[4] = {
      "FR2_2011_AA_electron.root",       // 2011 e
      "FR2_AA_ControlSample_ABCD.root",  // 2012 e
      "FR2_2011_AA_muon.root",           // 2011 mu
      "FR2_AA_muon.root"                 // 2012 mu
    };

    for (int i=0; i<4; ++i) {
      ZXWeightTables[i] = new TFile(file_name[i]);
      if (ZXWeightTables[i]==0) ;  
    }
    FileZXWeight = ZXWeightTables[n_file];
  }
  
  //cout << " File ? " << file_name[n_file] << " n_file = " << n_file << endl; 

  // Fake Rate Histo
  // TString histo_name [4] = {
  //     "eff_Z1ee_plus_electron",
  //     "eff_Z1ee_plus_muon",
  //     "eff_Z1mumu_plus_electron",
  //     "eff_Z1mumu_plus_muon"
  //   };
  
  TString Z1flavor = "Z1ee";     if(myZ1LepID==13) Z1flavor = "Z1mumu";
  TString Z2flavor = "electron"; if(myLepID==13)   Z2flavor = "muon";
  TString histo_name = "eff_"+Z1flavor+"_plus_"+Z2flavor;
  
  //cout << " histo = " << histo_name << endl;
  TH2D *h2_fake = (TH2D*)FileZXWeight->Get(histo_name);

  weight = h2_fake->GetBinContent(h2_fake->GetXaxis()->FindBin(myLepPt), h2_fake->GetYaxis()->FindBin(myLepEta));
  // cout << " binx     = " << h2_fake->GetXaxis()->FindBin(myLepPt) << " biny = " << h2_fake->GetYaxis()->FindBin(myLepEta) << endl;
  //   cout << " binXedge = " << h2_fake->GetXaxis()->GetBinLowEdge(h2_fake->GetXaxis()->FindBin(myLepPt))
  //        << " binYedge = " << h2_fake->GetYaxis()->GetBinLowEdge(h2_fake->GetYaxis()->FindBin(myLepEta)) << endl;
  //cout << " weight = " << weight << endl;

  return weight;

} // end of getFakeWeight


// Added by CO
Float_t HZZ4l::getZXfake_weight(const int year, const Int_t CandIndex)
{
  Float_t zx_weight = 1.;

  TLorentzVector m[2]; // only the Z2 legs
  m[0].SetPtEtaPhiM(Lep3Pt->at(CandIndex),Lep3Eta->at(CandIndex),Lep3Phi->at(CandIndex),0.1506583);
  m[1].SetPtEtaPhiM(Lep4Pt->at(CandIndex),Lep4Eta->at(CandIndex),Lep4Phi->at(CandIndex),0.1506583);

  int Z1_id = Lep1LepId->at(CandIndex); // assuming the Z1 pair is made from SF leptons...
  Int_t LepId[2];
  LepId[0] = Lep3LepId->at(CandIndex);
  LepId[1] = Lep4LepId->at(CandIndex);

  for(int ilep=0;ilep<2;ilep++) {
    zx_weight *= getFakeWeight(m[ilep].Pt(), m[ilep].Eta(), year, LepId[ilep], Z1_id);
  } // for loop on Z2 legs

  return zx_weight;
}
