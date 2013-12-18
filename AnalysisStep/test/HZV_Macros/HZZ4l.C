#define HZZ4l_cxx

#include "HZZ4l.h"

//Root includes
#include <TH1F.h>
#include <TH2D.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TString.h>

//RooFit includes
#include <RooRandom.h>

//Std includes
#include <boost/regex.hpp>
#include <fstream>
#include <iostream>

#include "../Plotter/root_lib/XSecReader.h"
#include "../../interface/PUReweight.h"

using namespace std;

namespace {
  const bool doMuEffCorr    = true;
  const bool doEleEffCorr   = true;
  const bool doHighmassCorr = true; 
  const bool saveJets       = true;
  const bool applySystMu    = false;
  const bool applySystEle   = false;

  const Int_t pwhg_flag = 0;    // 0 means standard high mass weights, -1 means weights - syst and +1 means weights + syst

  const Float_t Run2011AFraction = 0.465;

  XSecReader xsecRead7TeV("../Plotter/Xsection_v1.txt","../Plotter/Luminosity.txt");
  XSecReader xsecRead8TeV("../Plotter/Xsection8TeV_v2.txt","../Plotter/Luminosity.txt");

  static PUReweight PUWeighter(PUReweight::LEGACY);

  //This is used to always give the same shifts when doing systematics
  int Seeder = 0;
}

HZZ4l::HZZ4l(TChain *tree) : HZZ4lBase(tree) 
{
}

void HZZ4l::Loop(const Int_t Nevt_Gen, Int_t channelType, const string outputName, const bool is8TeV, const Int_t whichVH)
{
  if (fChain == 0) return;

  const Long64_t nentries = fChain->GetEntriesFast();

  //counters
  Int_t nTOTEv      = 0;
  Float_t nBestCand = 0.;
  Int_t nSel        = 0;

  //Variables necessary for systematics
  Int_t N_syst_Iter = 500;
  Float_t nBestSyst[N_syst_Iter];

  //bools
  bool isData(false);
  bool isCR(false);
  bool isSignal(false);
  bool isVBF(false);
  bool isVH(false);

  //normalization variables
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
  Float_t mPOLE                    = 0.;
  Float_t widthH                   = 0.; 
  Bool_t isMCRescaled(false);

  //Final tree for the UL
  Float_t myZZMass        = 0.;
  Float_t myZZPt          = 0.;
  Float_t myZZRapidity    = 0.;
  Float_t myZZMassErr     = 0.;
  Float_t myZZMassErrCorr = 0.;
  Float_t myZ2MassErr     = 0.;
  Float_t myZZLD          = 0.;
  Float_t myZZLDPSig      = 0.;
  Float_t myZZLDPBkg      = 0.;
  Float_t myZZpseudoLD    = 0.;
  Float_t myZZgravLD      = 0.;
  Float_t myZZVAKD        = 0.;
  Float_t myZZMEKDLD      = 0.;
  Float_t myZZMEKDpseudoLD= 0.;
  Float_t myZZMEKDgravLD  = 0.;
  Float_t myp0plus_melaNorm = 0.;
  Float_t myp0plus_mela   = 0.;
  Float_t myp0hplus_mela  = 0.;
  Float_t myp0minus_mela  = 0.;
  Float_t myp0plus_VAJHU  = 0.;
  Float_t myp0hplus_VAJHU = 0.;
  Float_t myp0minus_VAJHU = 0.;
  Float_t myp0plus_VAMCFM = 0.;
  Float_t myp1_mela       = 0.;
  Float_t myp1plus_mela   = 0.;
  Float_t myp1_VAJHU      = 0.;
  Float_t myp1plus_VAJHU  = 0.;
  Float_t myp2_mela       = 0.;
  Float_t myp2qqb_mela    = 0.;
  Float_t myp2_VAJHU      = 0.;
  Float_t myp2qqb_VAJHU   = 0.;
  Float_t mybkg_mela      = 0.;
  Float_t mybkg_VAMCFM    = 0.;
  Float_t mybkg_VAMCFMNorm= 0.;
  Float_t myggzz_VAMCFM   = 0.;
  Float_t myp0_pt         = 0.;
  Float_t myp0_y          = 0.;
  Float_t mybkg_pt        = 0.;
  Float_t mybkg_y         = 0.;
  Float_t myp0plus_m4l    = 0.;
  Float_t mybkg_m4l       = 0.;
  Float_t myZ1Mass        = 0.;
  Float_t myZ2Mass        = 0.;
  Float_t mycosthetastar  = 0.;
  Float_t myhelphi        = 0.;
  Float_t myhelcosthetaZ1 = 0.;
  Float_t myhelcosthetaZ2 = 0.;
  Float_t myphistarZ1     = 0.;
  Float_t myphistarZ2     = 0.;
  Float_t myLIWeight     = 0.;
  Long64_t myEventNumber  = 0;
  Int_t   myRunNumber     = 0;
  Int_t  mygenProcessId   = 0;
  Float_t mygenhpt        = 0;

  // Lepton scale and resolution systematics on probabilities
  Float_t myp0plus_m4l_ScaleUp   = 0;
  Float_t myp0plus_m4l_ScaleDown = 0;
  Float_t myp0plus_m4l_ResUp     = 0;
  Float_t myp0plus_m4l_ResDown   = 0;
  Float_t mybkg_m4l_ScaleUp      = 0;
  Float_t mybkg_m4l_ScaleDown    = 0;
  Float_t mybkg_m4l_ResUp        = 0;
  Float_t mybkg_m4l_ResDown      = 0;
  Int_t myCRflag                 = 0;


  Float_t myLep1Pt        = 0.;
  Float_t myLep1Eta       = 0.;
  Float_t myLep1ID        = 0.;
  Float_t myLep2Pt        = 0.;
  Float_t myLep2Eta       = 0.;
  Float_t myLep2ID        = 0.;
  Float_t myLep3Pt        = 0.;
  Float_t myLep3Eta       = 0.;
  Float_t myLep3ID        = 0.;
  Float_t myLep4Pt        = 0.;
  Float_t myLep4Eta       = 0.;
  Float_t myLep4ID        = 0.;

  Float_t myDiJetMass     = -99;
  Float_t myDiJetMassPlus = -99;
  Float_t myDiJetMassMinus= -99;
  Float_t myDiJetDEta     = -99;
  Float_t myFisher        = -99;
  Int_t myNJets           = 0; 
  std::vector<Double_t> myJetPt;
  std::vector<Double_t> myJetSigma;
  std::vector<Double_t> myJetEta;
  std::vector<Double_t> myJetPhi;
  std::vector<Double_t> myJetMass;
  std::vector<Double_t> myJetBTag;

  //Firstly create the File before the tree, because ROOT is designed my fucking monkeys high on crack
  TFile fOut(outputName.c_str(),"RECREATE");

  //Validation plots
  TH1F *hNvtxNoWeight = new TH1F("hNvtxNoWeight","hNvtxNoWeight",100,0,100);
  TH1F *hNvtxWeight = new TH1F("hNvtxWeight","hNvtxWeight",100,0,100);

  TH1F *hLepGenPt = new TH1F("hLepGenPt","hLepGenPt",1000,0,1000);

  //Weight plots
  TH1F *hPUWeight = new TH1F("hPUWeight","hPUWeight",250,0,25);
  TH1F *hTagProbeWeight = new TH1F("hTagPropeWeight","hTAGProbeWeight",300,0,10);

  //ZZMass plot with and without Giampi weight
  TH1F *hZZmassNoPoweg = new TH1F("hZZmassNoPoweg","hZZmassNoPoweg",700,0.,1800.);
  TH1F *hZZmassPoweg = new TH1F("hZZmassPoweg","hZZmassPoweg",700,0.,1800.);
  TH1F *hZZmassPowegReco = new TH1F("hZZmassPowegReco","hZZmassPowegReco",700,0.,1800.);

  TTree *SelTree = new TTree("SelectedTree","The Selection Tree");
  SelTree->Branch("ZZMass",&myZZMass,"ZZMass/F");
  SelTree->Branch("ZZMassErr",&myZZMassErr,"ZZMassErr/F");
  SelTree->Branch("ZZMassErrCorr",&myZZMassErrCorr,"ZZMassErrCorr/F");
  SelTree->Branch("Z2MassErr",&myZ2MassErr,"Z2MassErr/F");
  SelTree->Branch("ZZPt",&myZZPt,"ZZPt/F");
  SelTree->Branch("ZZRapidity",&myZZRapidity,"ZZRapidity/F");
  SelTree->Branch("ZZLD",&myZZLD,"ZZLD/F");
  SelTree->Branch("ZZgravLD",&myZZgravLD,"ZZgravLD/F");
  SelTree->Branch("ZZLD_PSig",&myZZLDPSig,"ZZLD_PSig/F");
  SelTree->Branch("ZZLD_PBkg",&myZZLDPBkg,"ZZLD_PBkg/F");
  SelTree->Branch("ZZpseudoLD",&myZZpseudoLD,"ZZpseudoLD/F");
  SelTree->Branch("ZZVAKD",&myZZVAKD,"ZZVAKD/F");
  SelTree->Branch("ZZMEKDLD",&myZZMEKDLD,"ZZMEKDLD/F");
  SelTree->Branch("ZZMEKDpseudoLD",&myZZMEKDpseudoLD,"ZZMEKDpseudoLD/F");
  SelTree->Branch("ZZMEKDgravLD",&myZZMEKDgravLD,"ZZMEKDgravLD/F");
  SelTree->Branch("p0plus_melaNorm",&myp0plus_melaNorm,"p0plus_melaNorm/F");
  SelTree->Branch("p0plus_mela",&myp0plus_mela,"p0plus_mela/F");
  SelTree->Branch("p0hplus_mela",&myp0hplus_mela,"p0hplus_mela/F");
  SelTree->Branch("p0minus_mela",&myp0minus_mela,"p0minus_mela/F");
  SelTree->Branch("p0plus_VAJHU",&myp0plus_VAJHU,"p0plus_VAJHU/F");
  SelTree->Branch("p0hplus_VAJHU",&myp0hplus_VAJHU,"p0hplus_VAJHU/F");
  SelTree->Branch("p0minus_VAJHU",&myp0minus_VAJHU,"p0minus_VAJHU/F");
  SelTree->Branch("p0plus_VAMCFM",&myp0plus_VAMCFM,"p0plus_VAMCFM/F");
  SelTree->Branch("p1_mela",&myp1_mela,"p1_mela/F");
  SelTree->Branch("p1plus_mela",&myp1plus_mela,"p1plus_mela/F");
  SelTree->Branch("p1_VAJHU",&myp1_VAJHU,"p1_VAJHU/F");
  SelTree->Branch("p1plus_VAJHU",&myp1plus_VAJHU,"p1plus_VAJHU/F");
  SelTree->Branch("p2_mela",&myp2_mela,"p2_mela/F");
  SelTree->Branch("p2qqb_mela",&myp2qqb_mela,"p2qqb_mela/F");
  SelTree->Branch("p2_VAJHU",&myp2_VAJHU,"p2_VAJHU/F");
  SelTree->Branch("p2qqb_VAJHU",&myp2qqb_VAJHU,"p2qqb_VAJHU/F");
  SelTree->Branch("bkg_mela",&mybkg_mela,"bkg_mela/F");
  SelTree->Branch("bkg_VAMCFM",&mybkg_VAMCFM,"bkg_VAMCFM/F");
  SelTree->Branch("bkg_VAMCFMNorm",&mybkg_VAMCFMNorm,"bkg_VAMCFMNorm/F");
  SelTree->Branch("ggzz_VAMCFM",&myggzz_VAMCFM,"ggzz_VAMCFM/F");
  SelTree->Branch("p0_pt",&myp0_pt,"p0_pt/F");
  SelTree->Branch("p0_y",&myp0_y,"p0_y/F");
  SelTree->Branch("bkg_pt",&mybkg_pt,"bkg_pt/F");
  SelTree->Branch("bkg_y",&mybkg_y,"bkg_y/F");
  SelTree->Branch("p0plus_m4l",&myp0plus_m4l,"p0plus_m4l/F");
  SelTree->Branch("bkg_m4l",&mybkg_m4l,"bkg_m4l/F");
  SelTree->Branch("Z1Mass",&myZ1Mass,"Z1Mass/F");
  SelTree->Branch("Z2Mass",&myZ2Mass,"Z2Mass/F");
  SelTree->Branch("MC_weight",&MC_weight,"MC_weight/F");
  SelTree->Branch("MC_weight_norm",&MC_weight_norm,"MC_weight_norm/F");
  SelTree->Branch("MC_weight_noxsec",&MC_weight_noxsec,"MC_weight_noxsec/F");
  SelTree->Branch("MC_weight_PUWeight",&PUWeight,"MC_weight_PUWeight/F");
  SelTree->Branch("MC_weight_powhegWeight",&powheg_weight,"MC_weight_powhegWeight/F");
  SelTree->Branch("MC_weight_dataMC",&dataMCWeight,"MC_weight_dataMC/F");
  SelTree->Branch("costhetastar",&mycosthetastar,"costhetastar/F");
  SelTree->Branch("helphi",&myhelphi,"helphi/F");
  SelTree->Branch("helcosthetaZ1",&myhelcosthetaZ1,"helcosthetaZ1/F");
  SelTree->Branch("helcosthetaZ2",&myhelcosthetaZ2,"helcosthetaZ2/F");
  SelTree->Branch("phistarZ1",&myphistarZ1,"phistarZ1/F");
  SelTree->Branch("phistarZ2",&myphistarZ2,"phistarZ2/F");
  SelTree->Branch("LIWeight",&myLIWeight,"LIWeight/F");
  SelTree->Branch("genProcessId",&mygenProcessId,"genProcessId/I");
  SelTree->Branch("GenHPt",&mygenhpt,"GenHPt/F");
  SelTree->Branch("p0plus_m4l_ScaleUp",&myp0plus_m4l_ScaleUp,"p0plus_m4l_ScaleUp/F");
  SelTree->Branch("p0plus_m4l_ScaleDown",&myp0plus_m4l_ScaleDown,"p0plus_m4l_ScaleDown/F");
  SelTree->Branch("p0plus_m4l_ResUp",&myp0plus_m4l_ResUp,"p0plus_m4l_ResUp/F");
  SelTree->Branch("p0plus_m4l_ResDown",&myp0plus_m4l_ResDown,"p0plus_m4l_ResDown/F");
  SelTree->Branch("bkg_m4l_ScaleUp",&mybkg_m4l_ScaleUp,"bkg_m4l_ScaleUp/F");
  SelTree->Branch("bkg_m4l_ScaleDown",&mybkg_m4l_ScaleDown,"bkg_m4l_ScaleDown/F");
  SelTree->Branch("bkg_m4l_ResUp",&mybkg_m4l_ResUp,"bkg_m4l_ResUp/F");
  SelTree->Branch("bkg_m4l_ResDown",&mybkg_m4l_ResDown,"bkg_m4l_ResDown/F");

  if(saveJets){
    SelTree->Branch("DiJetMass",&myDiJetMass,"DiJetMass/F");
    SelTree->Branch("DiJetMassPlus",&myDiJetMassPlus,"DiJetMassPlus/F");
    SelTree->Branch("DiJetMassMinus",&myDiJetMassMinus,"DiJetMassMinus/F");
    SelTree->Branch("DiJetDEta",&myDiJetDEta,"DiJetDEta/F");
    SelTree->Branch("Fisher",&myFisher,"Fisher/F");
    SelTree->Branch("NJets",&myNJets,"NJets/I");
    SelTree->Branch("JetPt",&myJetPt);
    SelTree->Branch("JetSigma",&myJetSigma);
    SelTree->Branch("JetEta",&myJetEta);
    SelTree->Branch("JetPhi",&myJetPhi);
    SelTree->Branch("JetMass",&myJetMass);
    SelTree->Branch("JetBTag",&myJetBTag);
  }

  //Start looping over the events
  Long64_t nb = 0;
  for (Long64_t jentry=0;jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);

    static const string filestring = fChain->GetCurrentFile()->GetName();
    string theSample = getSampleName(filestring);

    //We have to deal with VH splitting 
    Float_t MC_weight_initial_VH[3];    
    Float_t MC_weight_norm_initial_VH[3];    

    if(jentry == 0){
      //Check if we are dealing with CR
      TString TreeName(fChain->GetName());
      if(TreeName.Contains("CR")) isCR = true;
      if(isCR){
	SelTree->Branch("Lep1Pt",&myLep1Pt,"Lep1Pt/F");
	SelTree->Branch("Lep1Eta",&myLep1Eta,"Lep1Eta/F");
	SelTree->Branch("Lep1ID",&myLep1ID,"Lep1ID/F");
	SelTree->Branch("Lep2Pt",&myLep2Pt,"Lep2Pt/F");
	SelTree->Branch("Lep2Eta",&myLep2Eta,"Lep2Eta/F");
	SelTree->Branch("Lep2ID",&myLep2ID,"Lep2ID/F");
	SelTree->Branch("Lep3Pt",&myLep3Pt,"Lep3Pt/F");
	SelTree->Branch("Lep3Eta",&myLep3Eta,"Lep3Eta/F");
	SelTree->Branch("Lep3ID",&myLep3ID,"Lep3ID/F");
	SelTree->Branch("Lep4Pt",&myLep4Pt,"Lep4Pt/F");
	SelTree->Branch("Lep4Eta",&myLep4Eta,"Lep4Eta/F");
	SelTree->Branch("Lep4ID",&myLep4ID,"Lep4ID/F");
	SelTree->Branch("CRflag",&myCRflag,"CRflag/I");
      }

      //check wether this is data or MC, and which MC sample
      if(TString(filestring).Contains("Double")|| TString(filestring).Contains("MuEG")) isData = true;

      //If we are dealing with data
      if(isData){
	SelTree->Branch("EventNumber",&myEventNumber,"EventNumber/L");
	SelTree->Branch("RunNumber",&myRunNumber,"RunNumber/I");
      }
      //Now the MC samples
      else{
	if(TString(filestring).Contains("H")) isSignal = true;
	if(TString(filestring).Contains("VBF")) isVBF = true;
 	if(TString(filestring).Contains("VH")) {
	  isVH=true;
	  if (whichVH==0){
	    cout<<"ERROR: whichVH is not set for a VH sample!"<<endl;
	    return;
	  } 
 	}
	
	string VH[3] = {"WH", "ZH", "ttH"};
	if (isVH){
	  std::string::size_type pos = theSample.find("VH");
	  for (UInt_t iVH=0; iVH<3; ++iVH){
	    theSample.replace(pos,2,VH[iVH]);
	    Float_t xsec(0.);
	    if (!is8TeV) xsec = xsecRead7TeV.getWeight(theSample, "1fb-1","all", true);
	    else xsec = xsecRead8TeV.getWeight(theSample, "1fb-1","all", true);
	    MC_weight_initial_VH[iVH] = xsec/Nevt_Gen;
	    MC_weight_norm_initial_VH[iVH] = getNormalizedWeight(filestring,channelType,is8TeV,iVH+1);
	  }
	}


	//Reweighting for high mass
	TString sampleString(theSample);
	if(isSignal){
	  sampleString.Remove(0,1);
	  if(isVBF) sampleString.Remove(0,3);
	  mPOLE = sampleString.Atof();
	  string dummy_sampleString(sampleString.Data());
          getWeightFromFile(dummy_sampleString, is8TeV, isVBF);
	}

	//Cross section normalization
	if(!is8TeV) MC_weight_initial = xsecRead7TeV.getWeight(theSample, "1fb-1","all", true)/Nevt_Gen;
	else MC_weight_initial = xsecRead8TeV.getWeight(theSample, "1fb-1","all", true)/Nevt_Gen;

	cout << MC_weight_initial << endl;

	MC_weight_norm_initial = getNormalizedWeight(filestring,channelType,is8TeV,whichVH);
	if (MC_weight_norm_initial<0){
	  cout<<"'whichVH' flag is not set correctly: choose 0 (no VH), 1 (WH), 2 (ZH) or 3 (ttH)"<<endl;
	  return;
	}
	MC_weight_noxsec_initial = 1./Nevt_Gen;
	
	//find which PU reweighting
	//if(TString(filestring).Contains("ggZZ4l") || TString(filestring).Contains("ggZZ2l2l") /*|| TString(filestring).Contains("TTTo2L2Nu2B") || TString(filestring).Contains("ZZTo4tau")*/) {
	//  isMCRescaled = true;
	//}

      }//End of MC if

    }

    // Split VH samples
    if ( (whichVH == 1 && !(genProcessId == 26)) ||
	 (whichVH == 2 && !(genProcessId == 24)) ||
	 (whichVH == 3 && !(genProcessId == 121 || genProcessId == 122)) ) continue;
    
    if (whichVH!=0){
      MC_weight_initial = MC_weight_initial_VH[whichVH-1];
      MC_weight_norm_initial = MC_weight_norm_initial_VH[whichVH-1];
    }
    
    
    
    nTOTEv++;
    
    // MC weighting, including data/MC corrections and PU
    if(!isData){
      MC_weight_interm = MC_weight_initial;
      MC_weight_norm_interm = MC_weight_norm_initial;
      MC_weight_noxsec_interm = MC_weight_noxsec_initial;

      if(!is8TeV) PUWeight = PUWeighter.weight(2011,2011,NTrueInt);
      else if(is8TeV && isMCRescaled) PUWeight = PUWeighter.weight(2011,2012,NTrueInt);
      else if(is8TeV && !isMCRescaled) PUWeight = PUWeighter.weight(2012,2012,NTrueInt);
      else abort();
      
      if (fabs(PUWeight-PUWeight12)>0.001){
	cout << "ERROR: PUweights: " << PUWeight << " " << PUWeight12 << " " << is8TeV << " " << isMCRescaled << " " << NTrueInt << endl;
	//abort();
      }

      MC_weight_interm *= PUWeight;
      MC_weight_norm_interm *= PUWeight;
      MC_weight_noxsec_interm *= PUWeight;
      
      //A couple of validation plots
      hNvtxNoWeight->Fill(Nvtx);
      hNvtxWeight->Fill(Nvtx,PUWeight);
      hPUWeight->Fill(PUWeight);

      mygenProcessId = genProcessId;
      mygenhpt = GenHPt;
    }

    //Number of Higgs candidates in the event
    const Int_t NHiggs = ZZMass->size();

    bool foundBestCand(false);

    //Find best candidate in case we are looking at the control regions
    //Int_t iCRBC = -1;
    //if(isCR) iCRBC = findBestCRCand();


    //change the lineshape at high mass (Giampi's corrections)
    if(isSignal){
      if(mPOLE > 399.){
	hZZmassNoPoweg->Fill(GenHMass,MC_weight_interm); 

	if(doHighmassCorr){
	  powheg_weight = getPwhgWeight(GenHMass, pwhg_flag);
	    
	  //cout<<"giampi weight "<<powheg_weight <<endl;

	  MC_weight_interm *= powheg_weight;
	  MC_weight_norm_interm *= powheg_weight;
	  MC_weight_noxsec_interm *= powheg_weight;
	}

	hZZmassPoweg->Fill(GenHMass,MC_weight_interm); 
      }
    }

    //fill the plot with the highest generated pt for signal
    if(isSignal){
      Float_t muptmax = -1.;
      if(GenLep1Pt > muptmax) muptmax = GenLep1Pt;
      if(GenLep2Pt > muptmax) muptmax = GenLep2Pt;
      if(GenLep3Pt > muptmax) muptmax = GenLep3Pt;
      if(GenLep4Pt > muptmax) muptmax = GenLep4Pt;

      hLepGenPt->Fill(muptmax);
    }

    //Consider only the jets with a pt cut
    if(saveJets){
      //Clear the Jet collection to be saved
      myJetPt.clear();
      myJetSigma.clear();
      myJetEta.clear();
      myJetPhi.clear();
      myJetMass.clear();
      myJetBTag.clear();
      int countJet=0;
      for(int i=0;i<JetPt->size();i++){
	//if(JetPt->at(i) < 30.) continue; //Redundant

	myJetPt.push_back(JetPt->at(i));
	myJetSigma.push_back(JetSigma->at(i));
	myJetEta.push_back(JetEta->at(i));
	myJetPhi.push_back(JetPhi->at(i));
	myJetMass.push_back(JetMass->at(i));
	myJetBTag.push_back(JetBTag->at(i));
	if (myJetPt>30.) countJet++;
      }
      myDiJetMass = DiJetMass; //if <2 jets, DiJetMass=DiJetDEta=-99;
      myDiJetMassPlus = DiJetMassPlus;
      myDiJetMassMinus = DiJetMassMinus;
      myDiJetDEta = DiJetDEta;
      myNJets = countJet;
      if(myDiJetMass>-99)
	myFisher = 0.18*fabs(myDiJetDEta) + 1.92e-04*myDiJetMass;
      else
	myFisher=-99;
    }

    //Loop over all the H-> ZZ candidates in the event
    for(int nH=0; nH<NHiggs;nH++){

      //Do the final selection and best candidate selection. Depends on signal region or control region
      if(isCR && !CRflag->at(nH)) continue;
      if(!isCR && (ZZsel->at(nH) < 100. || nH != iBC) ) continue;

      //if(mPOLE > 399. && doHighmassCorr){
      //hZZmassNoPoweg->Fill(GenHMass,MC_weight/getPwhgWeight(GenHMass)); 
      //hZZmassPoweg->Fill(GenHMass,MC_weight); 
      //}

 
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
	hTagProbeWeight->Fill(dataMCWeight);
	MC_weight *= 19.64; //FIXME
      }
      
      nBestCand += MC_weight;
      nSel++;
      foundBestCand = true;

      //cout << getMCWeight(RunFraction, nH) << endl;

      hZZmassPowegReco->Fill(ZZMass->at(nH),MC_weight); 

      if(isSignal && !isVBF && !isData && !isCR && ( (applySystMu && channelType != 1) || (applySystEle && channelType != 0) ) ){
	
	for(int i=0; i<N_syst_Iter; i++){
	  //No need of other normalizations, as they wipe in the ratio
	  Seeder = i;
	  nBestSyst[i] += getMCWeight(RunFraction, nH);
	}
      }
      
      TLorentzVector m[4];
      m[0].SetPtEtaPhiM(Lep1Pt->at(nH),Lep1Eta->at(nH),Lep1Phi->at(nH),0.1506583);
      m[1].SetPtEtaPhiM(Lep2Pt->at(nH),Lep2Eta->at(nH),Lep2Phi->at(nH),0.1506583);
      m[2].SetPtEtaPhiM(Lep3Pt->at(nH),Lep3Eta->at(nH),Lep3Phi->at(nH),0.1506583);
      m[3].SetPtEtaPhiM(Lep4Pt->at(nH),Lep4Eta->at(nH),Lep4Phi->at(nH),0.1506583);
      
      TLorentzVector momH = m[0] + m[1] + m[2] + m[3];

      //Fill the variables that go into the final tree
      myZZMass = ZZMass->at(nH);
      myZZMassErr = ZZMassErr->at(nH);
      myZZMassErrCorr = ZZMassErrCorr->at(nH);
      myZ2MassErr = Z2MassErr->at(nH);
      myCRflag = (int)CRflag->at(nH);
      myZZPt = ZZPt->at(nH);
      myZZRapidity = momH.Rapidity();
      myZZLD = ZZLD->at(nH);
      myZZLDPSig = ZZLD_PSig->at(nH);
      myZZLDPBkg = ZZLD_PBkg->at(nH);
      myZZpseudoLD = ZZpseudoLD->at(nH);
      myZZgravLD = ZZgravLD->at(nH);
      myZZVAKD = ZZVAKD->at(nH); 
      myZZMEKDLD = ZZMEKDLD->at(nH);
      myZZMEKDpseudoLD = ZZMEKDpseudoLD->at(nH);
      myZZMEKDgravLD = ZZMEKDgravLD->at(nH);
      myp0plus_melaNorm = p0plus_melaNorm->at(nH);
      myp0plus_mela = p0plus_mela->at(nH);
      myp0hplus_mela = p0hplus_mela->at(nH);
      myp0minus_mela = p0minus_mela->at(nH);
      myp0plus_VAJHU = p0plus_VAJHU->at(nH);
      myp0hplus_VAJHU = p0hplus_VAJHU->at(nH);
      myp0minus_VAJHU = p0minus_VAJHU->at(nH);
      myp0plus_VAMCFM = p0plus_VAMCFM->at(nH);
      myp1_mela = p1_mela->at(nH);
      myp1plus_mela = p1plus_mela->at(nH);
      myp1_VAJHU = p1_VAJHU->at(nH);
      myp1plus_VAJHU = p1plus_VAJHU->at(nH);
      myp2_mela = p2_mela->at(nH);
      myp2qqb_mela = p2qqb_mela->at(nH);
      myp2_VAJHU = p2_VAJHU->at(nH);
      myp2qqb_VAJHU = p2qqb_VAJHU->at(nH);
      mybkg_mela = bkg_mela->at(nH);
      mybkg_VAMCFM = bkg_VAMCFM->at(nH);
      mybkg_VAMCFMNorm = bkg_VAMCFMNorm->at(nH);
      myggzz_VAMCFM = ggzz_VAMCFM->at(nH);
      myp0_pt = p0_pt->at(nH);
      myp0_y = p0_y->at(nH);
      mybkg_pt = bkg_pt->at(nH);
      mybkg_y = bkg_y->at(nH);
      myp0plus_m4l = p0plus_m4l->at(nH);
      mybkg_m4l = bkg_m4l->at(nH);
      myZ1Mass = Z1Mass->at(nH);
      myZ2Mass = Z2Mass->at(nH);
      mycosthetastar = costhetastar->at(nH);
      myhelphi = helphi->at(nH);
      myhelcosthetaZ1 = helcosthetaZ1->at(nH);
      myhelcosthetaZ2 = helcosthetaZ2->at(nH);
      myphistarZ1 = phistarZ1->at(nH);
      myphistarZ2 = phistarZ2->at(nH);
      myLIWeight = LIWeight->at(nH);
      myp0plus_m4l_ScaleUp = p0plus_m4l_ScaleUp->at(nH);
      myp0plus_m4l_ScaleDown = p0plus_m4l_ScaleDown->at(nH);
      myp0plus_m4l_ResUp = p0plus_m4l_ResUp->at(nH);
      myp0plus_m4l_ResDown = p0plus_m4l_ResDown->at(nH);
      mybkg_m4l_ScaleUp = bkg_m4l_ScaleUp->at(nH);
      mybkg_m4l_ScaleDown = bkg_m4l_ScaleDown->at(nH);
      mybkg_m4l_ResUp = bkg_m4l_ResUp->at(nH);
      mybkg_m4l_ResDown = bkg_m4l_ResDown->at(nH);
      
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

      if(isData){
	myEventNumber = EventNumber; 
	myRunNumber = RunNumber;       
      }

      break; //found best cand, no need to continue

    }//end of loop over the candidates

    if(foundBestCand) SelTree->Fill();

  }//end of loop over the events

  //Save the final tree
  fOut.cd();
  SelTree->Write();
  hNvtxNoWeight->Write();
  hNvtxWeight->Write();
  hLepGenPt->Write();
  hPUWeight->Write();
  hTagProbeWeight->Write();
  hZZmassNoPoweg->Write();
  hZZmassPoweg->Write();
  hZZmassPowegReco->Write();

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

  cout << "Total number of events analyzed " << nTOTEv << endl;
  cout << "Total number of events selected " << nSel << endl;
  cout << "Total number of best candidates found " << nBestCand << endl;
  cout << "################################################" << endl;

  return;
}

Int_t HZZ4l::findBestCRCand() const
{
  Int_t bestCand   = -1;
  Float_t bestMMPt = -1.;

  for(int nH=0; nH<ZZMass->size();nH++){
    const Float_t MMPt = Lep3Pt->at(nH) + Lep4Pt->at(nH);
    if(MMPt > bestMMPt){
      bestMMPt = MMPt;
      bestCand = nH;
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
  Float_t myLepID = abs(LepID);
  
  //avoid to go out of the TH boundary
  if(abs(LepID) == 13 && myLepPt > 99.) myLepPt = 99.;
  if(abs(LepID) == 11 && myLepPt > 199.) myLepPt = 199.;
  myLepEta = fabs(myLepEta);

  //Scale factors for data/MC efficiency
  static TFile fMuWeight("../Macros/scale_factors_muons2011.root");
  static TFile fMuWeight12("../Macros/scale_factors_muons2012.root");
  static TFile fElWeight("../Macros/scale_factors_ele2011.root");
  static TFile fElWeight12("../Macros/scale_factors_ele2012.root");
 

  static TH2D *hTH2D_Mu_All_2011A = (TH2D*)fMuWeight.Get("TH2D_ALL_2011A"); 
  static TH2D *hTH2D_Mu_All_2011B = (TH2D*)fMuWeight.Get("TH2D_ALL_2011B"); 
  static TH2D *hTH2D_Mu_All_2012  = (TH2D*)fMuWeight12.Get("TH2D_ALL_2012"); 
  
  static TH2D *hTH2D_El_All_2011A = (TH2D*)fElWeight.Get("h_electron_scale_factor_RECO_ID_ISO_SIP");
  static TH2D *hTH2D_El_All_2011B = (TH2D*)fElWeight.Get("h_electron_scale_factor_RECO_ID_ISO_SIP"); 
  static TH2D *hTH2D_El_All_2012  = (TH2D*)fElWeight12.Get("h_electron_scale_factor_RECO_ID_ISO_SIP");
    
  if(year == 0){
    if( myLepID == 13){                                               
      weight  = hTH2D_Mu_All_2011A->GetBinContent(hTH2D_Mu_All_2011A->GetXaxis()->FindBin(myLepPt),hTH2D_Mu_All_2011A->GetYaxis()->FindBin(LepEta));
      errCorr = hTH2D_Mu_All_2011A->GetBinError(hTH2D_Mu_All_2011A->GetXaxis()->FindBin(myLepPt),hTH2D_Mu_All_2011A->GetYaxis()->FindBin(LepEta));
      
    }
    
    else if( myLepID == 11){   

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
    
    else  if( myLepID == 11){   

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
    
    else if( myLepID == 11){
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
  if(myLepPt < 5. && myLepID == 13) weight = 1;

  if(weight < 0.001 || weight > 10.){
    cout << "pt = " << myLepPt << " eta = " << myLepEta << " weight = " << weight << endl;
    abort();  //no correction should be zero, if you find one, stop  //FIXME
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

Float_t HZZ4l::getNormalizedWeight(const string filestring, const Int_t channelType, const Bool_t is8TeV, const Int_t whichVH) const
{
  Float_t MC_normweight;

  string theSample = getSampleName(filestring);

  //Set the channel name
  string channelname;
  if(channelType == 0) channelname = "4mu";
  else if(channelType == 1) channelname = "4e";
  else if(channelType == 2) channelname = "2e2mu";
  else abort();

  string tagName = channelname + "_" + theSample;

  ifstream inputFile;
  if(is8TeV) inputFile.open("../Macros/saveNormValues_8TeV.txt");
  else inputFile.open("../Macros/saveNormValues_7TeV.txt");

  map<string,Float_t> NormMap;
  char normName[50];
  Float_t normValue = 0.;

  while(inputFile >> normName >> normValue) NormMap[normName] = normValue;

  MC_normweight = (NormMap.find(tagName))->second;

  TString mysample(theSample.c_str());

  if(mysample.Contains("H")){
    //Get the normalization. It's the first bin of the histo
    TFile fIn(filestring.c_str());
    string histoName = "ZZ" + channelname + "Tree/Counters";
    TH1F *nEventComplete = (TH1F*)fIn.Get(histoName.c_str());

    Int_t gen_4l_H = -1;
    
    // Convention on "whichVH flag":
    // 0 == no VH
    // 1 == WH
    // 2 == ZH
    // 3 == ttH
    if (whichVH==0) gen_4l_H = nEventComplete->GetBinContent(channelType + 2);
    else if (whichVH==1) gen_4l_H = nEventComplete->GetBinContent(channelType + 21);
    else if (whichVH==2) gen_4l_H = nEventComplete->GetBinContent(channelType + 24);
    else if (whichVH==3) gen_4l_H = nEventComplete->GetBinContent(channelType + 27);

    //cout << "Check " << nEventComplete->GetBinContent(1) << "  " << nEventComplete->GetBinContent(2) << "  " << nEventComplete->GetBinContent(3) << "  " << nEventComplete->GetBinContent(4) << endl;

    MC_normweight = 1./gen_4l_H;
  }

  //cout << MC_normweight << endl;
  return MC_normweight;
}

string HZZ4l::getSampleName(const string filestring) const
{
  //cout << "filestring = " << filestring.c_str() << endl;
  string theSample = "";

  string tmpre =  "(.*)_(.*).root";
  boost::regex re;
  boost::cmatch matches;

  TString mysample(filestring.c_str());
  if(mysample.Contains("_tau")) tmpre =  "(.*)_(.*)_(.*).root";

  try{   //make sure the regex is correctly formatted
    re = tmpre;
  }
  catch (boost::regex_error& e){
    cout << tmpre << " is not a valid regular expression: \""
	 << e.what() << "\"" << endl;
  }

  if (boost::regex_match(filestring.c_str(), matches, re)){
    // matches[0] contains the original string  
    // matches[n] contains a sub_match object for each matching subexpression
    for (int i = 1; i < matches.size(); i++){
      // sub_match::first and sub_match::second refer to the first and one past the last chars of the matching subexpression
      string match(matches[i].first, matches[i].second);
      //cout << "\tmatches[" << i << "] = " << match << endl;
      theSample = match;
      if(mysample.Contains("_tau") && i == 2) break; 
    }      
  }
  else{
    cout << "The regexp \"" << re << "\" does not match \"" << filestring << "\"" << endl;
  }

  return theSample;
}

void HZZ4l::getWeightFromFile(const string mPOLE, const Bool_t is8TeV, const Bool_t isVBF)
{
  string Energy = is8TeV ? "8TeV" : "7TeV";

  string filename;
  if(!isVBF) filename = "../Macros/WeightCards/mZZ_Higgs" + mPOLE + "_" + Energy + "_Lineshape+Interference.txt";
  //FIXME else filename = "WeightCards/CPS_only/vbfHZZ4L" + mPOLE + "_" + Energy + "_run_Lineshape.txt";

  cout << "File used to extract high mass weights: " << filename << endl;

  double bincenter, initial, pow, powp, powm, out, outp, outm;

  std::ifstream ifs(filename.c_str());
  while( ifs.good() ) {
    if(!isVBF) ifs >> bincenter >> initial >> pow >> powp >>  powm >> out>> outp >> outm;
    else ifs >> bincenter >> initial >> pow >> powp >>  powm;

    pwhg_bincenters.push_back(bincenter);
    if(initial > 0. && !isVBF){
      pwhg_weight.push_back(max(0.,out/initial));
      pwhg_weightP.push_back(max(0.,outp/initial));
      pwhg_weightM.push_back(max(0.,outm/initial));
    }
    else if(initial > 0. && isVBF){
      pwhg_weight.push_back(max(0.,pow/initial));
      pwhg_weightP.push_back(max(0.,powp/initial));
      pwhg_weightM.push_back(max(0.,powm/initial));
    }
    else{//weights are not defined if initial distribution is 0 => set weight to 0
      pwhg_weight.push_back(0.);
      pwhg_weightP.push_back(0.);
      pwhg_weightM.push_back(0.);
    }

  }

  return;
}

Float_t HZZ4l::getPwhgWeight(double mH, const Int_t WhichSide) const
{
  //Weight to be returned and +- error
  Float_t weight  = 1.;
  Float_t weightp = 1.;
  Float_t weightm = 1.;
  if(mH < pwhg_bincenters.front() || mH >  pwhg_bincenters.back()) return 0.; // set weights to 0 if out of range

  vector<double>::const_iterator low = lower_bound(pwhg_bincenters.begin(), pwhg_bincenters.end(),mH); 
  int lowindex = (low -  pwhg_bincenters.begin());

  if(mH == *low){//exact match
    weight  = pwhg_weight[lowindex];
    weightp = pwhg_weightP[lowindex];
    weightm = pwhg_weightM[lowindex];
  }else{//linear interpolation
    lowindex--; // lower_bound finds the first element not smaller than X
    weight  = pwhg_weight[lowindex]  + (mH - pwhg_bincenters[lowindex]) * (pwhg_weight[lowindex+1] - pwhg_weight[lowindex]) / (pwhg_bincenters[lowindex+1] - pwhg_bincenters[lowindex]);
    weightp = pwhg_weightP[lowindex] + (mH - pwhg_bincenters[lowindex]) * (pwhg_weight[lowindex+1] - pwhg_weight[lowindex]) / (pwhg_bincenters[lowindex+1] - pwhg_bincenters[lowindex]);
    weightm = pwhg_weightM[lowindex] + (mH - pwhg_bincenters[lowindex]) * (pwhg_weight[lowindex+1] - pwhg_weight[lowindex]) / (pwhg_bincenters[lowindex+1] - pwhg_bincenters[lowindex]);
  }

  if(weight > 20. ) cout << mH << " " << *low << " " << lowindex << " " << pwhg_weight[lowindex] << " " << weight << endl;

  switch(WhichSide){

  case 0:
    return weight;

  case -1:
    return weightm;

  case 1:
    return weightp;

  default:
    abort();

  }

  return -999.;
}
