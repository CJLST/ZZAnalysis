// -*- C++ -*-
// 
// Fill a tree for selected candidates.
//


// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/LuminosityBlock.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <DataFormats/Common/interface/TriggerResults.h>
#include <FWCore/Common/interface/TriggerNames.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include <DataFormats/Common/interface/View.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <DataFormats/JetReco/interface/PFJet.h>
#include <DataFormats/JetReco/interface/PFJetCollection.h>
#include <DataFormats/Math/interface/LorentzVector.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <DataFormats/Common/interface/MergeableCounter.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>
#include <ZZAnalysis/AnalysisStep/interface/MCHistoryTools.h>
#include <ZZAnalysis/AnalysisStep/interface/PUReweight.h>
#include <ZZAnalysis/AnalysisStep/interface/bitops.h>
#include <ZZAnalysis/AnalysisStep/interface/Fisher.h>
#include "ZZ4lConfigHelper.h"
#include "HZZ4lNtupleFactory.h"

#include <TRandom3.h>
#include <TH2D.h>
#include "TLorentzVector.h"

#include "ZZAnalysis/AnalysisStep/interface/PUReweight.h"

#include <string>

namespace {
  bool writePhotons = false;  // Write photons in the tree. Note: must be set also in HZZ4lNtupleFactory.cc
  bool writeJets = true;     // Write jets in the tree. FIXME: make this configurable
  bool addKinRefit = false;
  bool addVtxFit = false;
}

using namespace std;
using namespace edm;
//
// class declaration
//
class HZZ4lNtupleMaker : public edm::EDAnalyzer {
public:
  explicit HZZ4lNtupleMaker(const edm::ParameterSet&);
  ~HZZ4lNtupleMaker();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  void BookAllBranches();  
  virtual void FillCandidate(const pat::CompositeCandidate& higgs, bool evtPass, const edm::Event&, const Int_t CRflag);
  virtual void FillPhoton(const pat::Photon& photon);
  virtual void FillJet(const pat::Jet& jet);
  virtual void endJob() ;

  Float_t getAllWeight(const Float_t LepPt, const Float_t LepEta, Int_t LepID) const;
  Float_t getHqTWeight(double mH, double genPt) const;
  Float_t getFakeWeight(const Float_t LepPt, const Float_t LepEta, Int_t LepID, Int_t LepZ1ID);

  // ----------member data ---------------------------
  ZZ4lConfigHelper myHelper;
  int theChannel;
  std::string theCandLabel;
  TString theFileName;

  HZZ4lNtupleFactory *myTree;
  TH1F *hCounter;

  Bool_t isMC;

  bool applyTrigger;    // Keep only events passing trigger
  bool applySkim;       //   "     "      "     skim
  bool skipEmptyEvents; // Skip events whith no selected candidate (otherwise, gen info is preserved for all events)
  Float_t xsec;
  int year;
  
  PUReweight reweight;

  //counters
  Float_t Nevt_Gen;
  Float_t Nevt_Gen_lumiBlock;

  Float_t gen_ZZ4mu;
  Float_t gen_ZZ4e;
  Float_t gen_ZZ2mu2e;
  Float_t gen_ZZ2l2tau;
  Float_t gen_ZZ4mu_EtaAcceptance;
  Float_t gen_ZZ4mu_LeptonAcceptance;
  Float_t gen_ZZ4e_EtaAcceptance;
  Float_t gen_ZZ4e_LeptonAcceptance;
  Float_t gen_ZZ2mu2e_EtaAcceptance; 
  Float_t gen_ZZ2mu2e_LeptonAcceptance; 
  Float_t gen_BUGGY;
  Float_t gen_Unknown;
  
  Float_t gen_sumPUWeight;
  Float_t gen_sumGenMCWeight;
  Float_t gen_sumWeights;

  string sampleName;
   
  TH2D *hTH2D_Mu_All;// = (TH2D*)fMuWeight.Get("TH2D_ALL_2011A"); 
  TH2D *hTH2D_El_All;//  = (TH2D*)fElWeight12.Get(eleSFname.Data());
  TH2D* h_weight; //HqT weights
  //TH2F *h_ZXWeightMuo;
  //TH2F *h_ZXWeightEle;
  TH2D* h_ZXWeight[4];

};

//
// constructors and destructor
//
HZZ4lNtupleMaker::HZZ4lNtupleMaker(const edm::ParameterSet& pset) :
  myHelper(pset),
  reweight(),
  hTH2D_Mu_All(0),
  hTH2D_El_All(0),
  h_weight(0)
  //h_ZXWeight(0)
{
  //cout<< "Beginning Constructor\n\n\n" <<endl;
  theCandLabel = pset.getUntrackedParameter<string>("CandCollection"); // Name of input ZZ collection
  theChannel = myHelper.channel(); // Valid options: ZZ, ZLL, ZL 
  theFileName = pset.getUntrackedParameter<string>("fileName"); 
  skipEmptyEvents = pset.getParameter<bool>("skipEmptyEvents"); // Do not store 
  sampleName = pset.getParameter<string>("sampleName");
  xsec = pset.getParameter<double>("xsec");
  year = pset.getParameter<int>("setup");

  if (skipEmptyEvents) {
    applyTrigger=true;
    applySkim=true;
  } else {
    applyTrigger=false;
    applySkim=false;    
  }

  isMC = myHelper.isMC();

  Nevt_Gen = 0;
  Nevt_Gen_lumiBlock = 0;

  //For Efficiency studies
  gen_ZZ4mu = 0;
  gen_ZZ4e = 0;
  gen_ZZ2mu2e = 0;
  gen_ZZ2l2tau = 0;
  gen_ZZ4mu_EtaAcceptance = 0;
  gen_ZZ4mu_LeptonAcceptance = 0;
  gen_ZZ4e_EtaAcceptance = 0;
  gen_ZZ4e_LeptonAcceptance = 0;
  gen_ZZ2mu2e_EtaAcceptance = 0;
  gen_ZZ2mu2e_LeptonAcceptance = 0;
  gen_BUGGY = 0;
  gen_Unknown = 0;

  gen_sumPUWeight = 0.f;
  gen_sumGenMCWeight = 0.f;
  gen_sumWeights =0.f;
  
  //Scale factors for data/MC efficiency
  //FIXME: to adjust for 13 TeV
  float Run2011AFraction=0.465;
  TString yearString;yearString.Form("TH2D_ALL_%d",year);
  if(year==2011){
    TRandom3 randomGenerator(0);
    Float_t whatPeriod = randomGenerator.Uniform();
    if(whatPeriod < Run2011AFraction) yearString.Append("A");
    else yearString.Append("B");
  }
  TString filename;filename.Form("ZZAnalysis/AnalysisStep/test/Macros/scale_factors_muons%d.root",year); 
  if(year==2015){
    filename.Form("ZZAnalysis/AnalysisStep/test/Macros/scale_factors_muons%d.root",2012); //FIXME
    yearString.Form("TH2D_ALL_%d",2012);
  }
  edm::FileInPath fip(filename.Data());
  std::string fipPath=fip.fullPath();
  TFile *fMuWeight = TFile::Open(fipPath.data(),"READ");
  hTH2D_Mu_All = (TH2D*)fMuWeight->Get(yearString.Data())->Clone(); 

  filename.Form("ZZAnalysis/AnalysisStep/test/Macros/scale_factors_ele%d.root",year); 
  if(year==2015)filename.Form("ZZAnalysis/AnalysisStep/test/Macros/scale_factors_ele%d.root",2012);//FIXME  
  edm::FileInPath fipEle(filename.Data());
  fipPath=fipEle.fullPath();
  TFile *fEleWeight = TFile::Open(fipPath.data(),"READ");
  hTH2D_El_All = (TH2D*)fEleWeight->Get("h_electronScaleFactor_RecoIdIsoSip")->Clone();
     
  //HqT weights
  edm::FileInPath HqTfip("ZZAnalysis/AnalysisStep/test/Macros/HqTWeights.root");
  fipPath=HqTfip.fullPath();
  TFile *fHqt = TFile::Open(fipPath.data(),"READ");
  h_weight = (TH2D*)fHqt->Get("wH_400")->Clone();//FIXME: Ask simon to provide the 2D histo

  //CR fake rate weight
  filename.Form("ZZAnalysis/AnalysisStep/test/Macros/FR2_2011_AA_electron.root");
  if(year==2015)filename.Form("ZZAnalysis/AnalysisStep/test/Macros/FR2_AA_ControlSample_ABCD.root");
  edm::FileInPath fipEleZX(filename.Data());
  fipPath=fipEleZX.fullPath();
  TFile *FileZXWeightEle = TFile::Open(fipPath.data(),"READ");
  
  filename.Form("ZZAnalysis/AnalysisStep/test/Macros/FR2_2011_AA_muon.root");
  if(year==2015)filename.Form("ZZAnalysis/AnalysisStep/test/Macros/FR2_AA_muon.root");
  edm::FileInPath fipMuZX(filename.Data());
  fipPath=fipMuZX.fullPath();  
  TFile *FileZXWeightMuo = TFile::Open(fipPath.data(),"READ");
  
  h_ZXWeight[0]=(TH2D*)FileZXWeightEle->Get("eff_Z1ee_plus_electron")->Clone();
  h_ZXWeight[1]=(TH2D*)FileZXWeightEle->Get("eff_Z1mumu_plus_electron")->Clone();
  h_ZXWeight[2]=(TH2D*)FileZXWeightMuo->Get("eff_Z1ee_plus_muon")->Clone();
  h_ZXWeight[3]=(TH2D*)FileZXWeightMuo->Get("eff_Z1mumu_plus_muon")->Clone();

  fMuWeight->Close();
  fEleWeight->Close();
  fHqt->Close();
  FileZXWeightEle->Close();
  FileZXWeightMuo->Close();
    
  //cout<<"Booking branches\n\n\n"<<endl;

}

HZZ4lNtupleMaker::~HZZ4lNtupleMaker()
{
}


// ------------ method called for each event  ------------
void HZZ4lNtupleMaker::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  // Primary vertices
  Handle<vector<reco::Vertex> >  vertexs;
  event.getByLabel("goodPrimaryVertices",vertexs);
  
  //----------------------------------------------------------------------
  // Analyze MC truth; collect MC weights and update counters (this is done for all generated events, 
  // including those that do not pass skim, trigger etc!)
  int nObsInt  = -1;
  float nTrueInt = -1.;
  Float_t weight2 = 1.;
  Int_t genFinalState = -1;
  Int_t genProcessId = -1;
  Float_t genHEPMCweight = 1.;
  Short_t genExtInfo = -1;

  bool gen_ZZ4lInEtaAcceptance = false;   // All 4 gen leptons in eta acceptance
  bool gen_ZZ4lInEtaPtAcceptance = false; // All 4 gen leptons in eta,pT acceptance

  const reco::Candidate * genH = 0;
  std::vector<const reco::Candidate *> genZLeps;
  std::vector<const reco::Candidate *> genAssocLeps;
  if (isMC) {
    // get PU weights
    Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    event.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);

    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      if(PVI->getBunchCrossing() == 0) { 
	nObsInt  = PVI->getPU_NumInteractions();
	nTrueInt = PVI->getTrueNumInteractions();
	break;
      } 
    }

    int source = myHelper.sampleType();
    int target = myHelper.setup();
    weight2 = reweight.weight(source,target,nTrueInt);

    MCHistoryTools mch(event, sampleName);
    genFinalState = mch.genFinalState();
    genProcessId = mch.getProcessID();
    genHEPMCweight = mch.gethepMCweight();
    genExtInfo = mch.genAssociatedFS();

    // keep track of sum of weights
    gen_sumPUWeight    += weight2;
    gen_sumGenMCWeight += genHEPMCweight;
    //weight = weight2*genHEPMCweight;
    gen_sumWeights     += weight2*genHEPMCweight;

    mch.genAcceptance(gen_ZZ4lInEtaAcceptance, gen_ZZ4lInEtaPtAcceptance);

    ++Nevt_Gen_lumiBlock;
    if (genFinalState == EEEE) {
      ++gen_ZZ4e;
      if (gen_ZZ4lInEtaAcceptance) ++gen_ZZ4e_EtaAcceptance;
      if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ4e_LeptonAcceptance;
    } else if (genFinalState == MMMM) {
      ++gen_ZZ4mu;
      if (gen_ZZ4lInEtaAcceptance) ++gen_ZZ4mu_EtaAcceptance;
      if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ4mu_LeptonAcceptance;
    } else if (genFinalState == EEMM) {
      ++gen_ZZ2mu2e;
      if (gen_ZZ4lInEtaAcceptance) ++gen_ZZ2mu2e_EtaAcceptance;
      if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ2mu2e_LeptonAcceptance;
    } else if (genFinalState == LLTT){
      ++gen_ZZ2l2tau;
    } else if (genFinalState == BUGGY){ // handle H->ddbar 2012 generator bug!!!
      ++gen_BUGGY;
      return; // BUGGY events are skipped
    } else {
      ++gen_Unknown;
    }
    

    //Information on generated candidates, will be used later
    genH = mch.genH();
    genZLeps     = mch.sortedGenZZLeps();
    genAssocLeps = mch.genAssociatedLeps();
  }
  // End of MC history analysis ------------------------------------------


  // Get candidate collection
  edm::Handle<edm::View<pat::CompositeCandidate> > candHandle;
  event.getByLabel(theCandLabel, candHandle);
  const edm::View<pat::CompositeCandidate>* cands = candHandle.product();
  myTree->InitializeVariables();

  if (skipEmptyEvents && cands->size() == 0) return; // Skip events with no candidate, unless skipEmptyEvents = false

  // For Z+L CRs, we want only events with exactly 1 Z+l candidate. FIXME: this has to be reviewed.
  if (theChannel==ZL && cands->size() != 1) return;


  // Apply MC filter (skip event)
  if (isMC && !(myHelper.passMCFilter(event))) return;

  // Apply skim
  Short_t trigWord=0;
  bool evtPassSkim = myHelper.passSkim(event, trigWord);
  if (applySkim && !evtPassSkim) return;

  // Apply trigger request (skip event)
  bool evtPassTrigger = myHelper.passTrigger(event, trigWord);
  if (applyTrigger && !evtPassTrigger) return;

  //Fill MC truth information
  if (isMC) {
    if(genH != 0){
      myTree->FillHGenInfo(genH->p4(),getHqTWeight(genH->p4().M(),genH->p4().Pt()));
    }
    else if(genZLeps.size()==4){ // for 4l events take the mass of the ZZ(4l) system
      myTree->FillHGenInfo((genZLeps.at(0)->p4()+genZLeps.at(1)->p4()+genZLeps.at(2)->p4()+genZLeps.at(3)->p4()),0);
    }

    if (genFinalState!=BUGGY) {
    
      float eff_weight =1.;
      for(int nLep=0;nLep<(int)genZLeps.size();nLep++){
        int lepid = genZLeps.at(nLep)->pdgId();
        if(fabs(lepid) ==13 || fabs(lepid) ==11)eff_weight *= getAllWeight(genZLeps.at(nLep)->pt(), genZLeps.at(nLep)->eta(),genZLeps.at(nLep)->pdgId());
      }
      if (genZLeps.size()==4) {
	
	// "generated Zs" defined with standard pairing applied on gen leptons (genZLeps is sorted by MCHistoryTools)
	myTree->FillZGenInfo(genZLeps.at(0)->p4()+genZLeps.at(1)->p4(),
			     genZLeps.at(2)->p4()+genZLeps.at(3)->p4());

	// Gen leptons
	myTree->FillLepGenInfo(genZLeps.at(0)->pdgId(), genZLeps.at(1)->pdgId(), genZLeps.at(2)->pdgId(), genZLeps.at(3)->pdgId(),
			       genZLeps.at(0)->p4(), genZLeps.at(1)->p4(), genZLeps.at(2)->p4(), genZLeps.at(3)->p4(),eff_weight);	
      }

      if (genZLeps.size()==3) {
	myTree->FillLepGenInfo(genZLeps.at(0)->pdgId(), genZLeps.at(1)->pdgId(), genZLeps.at(2)->pdgId(), 0,
			       genZLeps.at(0)->p4(), genZLeps.at(1)->p4(), genZLeps.at(2)->p4(), *(new math::XYZTLorentzVector),eff_weight);
      }
      if (genZLeps.size()==2) {
	myTree->FillLepGenInfo(genZLeps.at(0)->pdgId(), genZLeps.at(1)->pdgId(), 0, 0,
			       genZLeps.at(0)->p4(), genZLeps.at(1)->p4(), *(new math::XYZTLorentzVector), *(new math::XYZTLorentzVector),eff_weight);
      }

      if (genAssocLeps.size()==1 || genAssocLeps.size()==2) {
	myTree->FillAssocLepGenInfo(genAssocLeps);
      }

    }
  }

  // Photons (store them only for events with at least 1 candidate)
  // FIXME: should rather write used FSR photons, with info on matching lepton.
  if (writePhotons && cands->size()!=0) {
    edm::Handle<vector<pat::Photon> > photons;
    if (event.getByLabel("cmgPhotonSel",photons)) {
      for( vector<pat::Photon>::const_iterator ph = photons->begin(); ph != photons->end(); ++ph) {
	if(ph->pt() > 2. && fabs(ph->eta()) < 2.8) FillPhoton(*ph);
      }
    }
  }


  // Jets
  Handle<edm::View<pat::Jet> > CleanedJets;
  event.getByLabel("cleanJets", CleanedJets);
  vector<const pat::Jet*> cleanedJets;
  vector<const pat::Jet*> cleanedJetsPt30;
  int nCleanedJetsPt30BTagged = 0;
  for(edm::View<pat::Jet>::const_iterator jet = CleanedJets->begin(); jet != CleanedJets->end(); ++jet){
    cleanedJets.push_back(&*jet);
    if(jet->pt()>30){
      cleanedJetsPt30.push_back(&*jet);
      if(jet->userFloat("isBtagged")) nCleanedJetsPt30BTagged++;
    }
  }

  float detajj = -99.f;
  float Mjj    = -99.f;
  float Fisher = -99.f;
  if (cleanedJetsPt30.size()>=2) {
    detajj = cleanedJetsPt30[0]->eta()-cleanedJetsPt30[1]->eta();
    Mjj = (cleanedJetsPt30[0]->p4()+cleanedJetsPt30[1]->p4()).M();
    Fisher = fisher(Mjj,detajj);      
  }

  if (writeJets){

    if(theChannel!=ZL){
      for (unsigned int i=0; i<cleanedJets.size(); i++) {
	FillJet(*(cleanedJets.at(i)));
      }

      // Note that jets variables are filled for jets above 20 GeV, to allow JES studies.
      // detajj, Mjj and ZZFisher are filled only for true dijet events (jets above 30 GeV)
      if(cleanedJets.size()>1){ 
	const pat::Jet& myjet1 = *(cleanedJets.at(0)); 
	const pat::Jet& myjet2 = *(cleanedJets.at(1));
	math::XYZTLorentzVector jet1 = myjet1.p4();
	math::XYZTLorentzVector jet2 = myjet2.p4();
	Float_t jesUnc1 = 0.;//myjet1.uncOnFourVectorScale();
	Float_t jesUnc2 = 0.;//myjet2.uncOnFourVectorScale();
	math::XYZTLorentzVector jetScalePlus1 = jet1*(1+jesUnc1);
	math::XYZTLorentzVector jetScaleMinus1 = jet1*(1-jesUnc1);
	math::XYZTLorentzVector jetScalePlus2 = jet2*(1+jesUnc2);
	math::XYZTLorentzVector jetScaleMinus2 = jet2*(1-jesUnc2);
	Float_t MjjPlus = (jetScalePlus1+jetScalePlus2).M();
	Float_t MjjMinus = (jetScaleMinus1+jetScaleMinus2).M();
	//myTree->FillDiJetInfo(Mjj,MjjPlus,MjjMinus,detajj,Fisher);
	TString names[]={"DiJetMass","DiJetMassPlus","DiJetMassMinus","DiJetDEta","DiJetFisher"};
	double dijetinfo[]={Mjj,MjjPlus,MjjMinus,detajj,Fisher};
	myTree->SetVariables((TString *)names,(double *)dijetinfo,5);
      }
    }
  }
  
  Handle<vector<reco::MET> > pfmetcoll;
  event.getByLabel("slimmedMETs", pfmetcoll);
  float pfmet = -1;
  if(pfmetcoll.isValid()){
    pfmet = pfmetcoll->front().pt();
  }
  
    //Save general event info in the tree
  Int_t NbestCand = -1; //FIXME now store only 1 candidate in the SR, but we still have to save iBC correctly in the SR
  if (theChannel==ZZ) NbestCand=0;
  //myTree->FillEventInfo(event.id().run(), event.id().event(), event.luminosityBlock(), NbestCand, vertexs->size(), nObsInt, nTrueInt, weight2, pfmet, pfjetscoll->size(), cleanedJets.size(), cleanedJetsPt30.size(), nCleanedJetsPt30BTagged, genFinalState, genProcessId, genHEPMCweight, trigWord, genExtInfo,xsec);
  double eventinfo[17]={(double)event.id().run(),(double)event.luminosityBlock(),(double)NbestCand, (double)vertexs->size(), 
    (double)nObsInt, (double)nTrueInt, (double)weight2, (double)pfmet, (double)cleanedJets.size(), 
    (double)cleanedJetsPt30.size(), (double)nCleanedJetsPt30BTagged, (double)genFinalState, (double)genProcessId, 
    (double)genHEPMCweight, (double)trigWord, (double)genExtInfo,(double)xsec};
  TString eventnames[17]={
    "RunNumber",
    "LumiNumber",
    "iBC",
    "Nvtx",
    "NObsInt",
    "NTrueInt",
    "PUWeight12",
    "PFMET",
    //"nJets",
    "nCleanedJets",
    "nCleanedJetsPt30",
    "nCleanedJetsPt30BTagged",
    "genFinalState",
    "genProcessId",
    "genHEPMCweight",
    "trigWord",
    "genExtInfo",
    "xsec",
  };
  myTree->SetVariables((TString *)eventnames,(double *)eventinfo,17);
  myTree->SetVariableLong("EventNumber",event.id().event());
  
  //Loop on the candidates
  int nFilled=0;
  for( edm::View<pat::CompositeCandidate>::const_iterator cand = cands->begin(); cand != cands->end(); ++cand) {
    Int_t CRFLAG=0;
    bool candIsBest = cand->userFloat("isBestCand");

    //    int candChannel = cand->userFloat("candChannel"); // This is currently the product of pdgId of leptons (eg 14641, 28561, 20449)
    
    if (theChannel==ZLL) {
      // AA CRs
      if(cand->userFloat("isBestCRZLLss")&&cand->userFloat("CRZLLss"))set_bit(CRFLAG,CRZLLss);      

      // A CRs
      if(cand->userFloat("isBestCRZLLos_2P2F")&&cand->userFloat("CRZLLos_2P2F"))set_bit(CRFLAG,CRZLLos_2P2F);      
      if(cand->userFloat("isBestCRZLLos_3P1F")&&cand->userFloat("CRZLLos_3P1F"))set_bit(CRFLAG,CRZLLos_3P1F);
    }
    
    //    if(theChannel==ZL){} Nothing special in this case
 
    if (!(candIsBest||CRFLAG)) continue; // Skip events other than the best cand (or CR candidates in the CR)
    
    //For the SR, also fold information about acceptance in CRflag 
    if (isMC && (theChannel==EEEE||theChannel==MMMM||theChannel==EEMM)) {
      if (gen_ZZ4lInEtaAcceptance)   set_bit(CRFLAG,28);
      if (gen_ZZ4lInEtaPtAcceptance) set_bit(CRFLAG,29);
    }
    FillCandidate(*cand, evtPassTrigger&&evtPassSkim, event, CRFLAG);
    if(CRFLAG && !candIsBest && cleanedJets.size()==0 && Fisher>=0){
      // for CR, the jet list is empty; we "fake" the NJets30 variable for consistency
      myTree->SetVariable("nCleanedJets",2);
    }
    myTree->FillCurrentTree();//Do not reinitialize the event variables
    ++nFilled;
    /*
    //Assume only one isBestCand, and that for SR CRflag=0
    if(nFilled==0){
      FillCandidate(*cand, evtPassTrigger&&evtPassSkim, event, CRFLAG);
    }else {
      myTree->FillCurrentTree();
      isAlreadyFilled=true;
    }
    ++nFilled;
    */
  }//end loop over candidates
  //cout<<"nFilled "<<nFilled<<endl;
  // Events with no ZZ candidate have already been skipped at the beginning, but in case of CRs there could be no cands with CRflag!=0, so none was filled. 
  //if (skipEmptyEvents && nFilled==0 && theChannel==ZLL) return;

  if (skipEmptyEvents==false && nFilled==0)myTree->FillCurrentTree();
 
  // Final call to save the tree entry, and reset tree variables
  /*
  if(!isAlreadyFilled){
    myTree->FillEvent();
  }else{
    myTree->InitializeVariables();
  }
  */
  return;
}

void HZZ4lNtupleMaker::FillPhoton(const pat::Photon& photon)
{
  const Float_t photPt  = photon.pt();
  const Float_t photEta = photon.eta();
  const Float_t photPhi = photon.phi();

  TString names[3]={"PhotPt","PhotEta","PhotPhi"};
  double vars[3]={photPt, photEta, photPhi};
  myTree->SetVariables((TString *)names,(double *)vars,3);
  return;
}


void HZZ4lNtupleMaker::FillJet(const pat::Jet& jet)
{
  const Float_t jetPt  = jet.pt();
  const Float_t jetEta = jet.eta();
  const Float_t jetPhi = jet.phi();
  const Float_t jetMass = jet.p4().M();
  const Float_t jetBTagger = jet.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
  const Float_t jetIsBtagged = jet.userFloat("isBtagged");
  const Float_t jetQGLikelihood = jet.userFloat("qgLikelihood");
  const Float_t jesUnc = 0.;//jet.uncOnFourVectorScale();

  //myTree->FillJetInfo( );
  TString names[8]={"JetPt","JetEta","JetPhi","JetMass","JetBTagger","JetIsBtagged","JetQGLikelihood","JetSigma"};
  double vars[8]={jetPt, jetEta, jetPhi, jetMass, jetBTagger, jetIsBtagged, jetQGLikelihood, jesUnc};
  myTree->SetVariables((TString *)names,(double *)vars,8);
   
  return;
}


void HZZ4lNtupleMaker::FillCandidate(const pat::CompositeCandidate& cand, bool evtPass, const edm::Event& event, Int_t CRflag)
{
  //Initialize a new candidate into the tree
  myTree->createNewCandidate();

  //Find out if this is a proper 4l candidate or a Z+l candidate

  //  cout << "AAA " << cand.daughter(1)->numberOfDaughters() <<endl;  

  //FIXME: Dobbiamo salvare info su qualita' del fit a H?
  //Chi2 e chi2constr

  //Fill the info on the Higgs candidate
  const Float_t ZZMass = cand.p4().mass();
  const Float_t ZZMassErr = cand.userFloat("massError");
  const Float_t ZZMassErrCorr = cand.userFloat("massErrorCorr");
  const Float_t ZZMassPreFSR = cand.userFloat("m4l");
  const Float_t ZZMassRefit = cand.userFloat("m4lRef");
  const Float_t Chi2KinFit = cand.userFloat("chi2Fit");
  const Float_t ZZMassCFit = cand.userFloat("CFitM");
  const Float_t Chi2CFit = cand.userFloat("CFitChi2");
  
  const Float_t ZZPt  = cand.p4().pt();
  const Float_t ZZEta = cand.p4().eta();
  const Float_t ZZPhi = cand.p4().phi();
  //const Float_t ZZLD = cand.userFloat("LD");
  //const Float_t ZZLDPSig = cand.userFloat("PSig");
  //const Float_t ZZLDPBkg = cand.userFloat("PBkg");
  //const Float_t ZZpseudoLD = cand.userFloat("pseudoLD");
  //const Float_t ZZgravLD = cand.userFloat("spin2PMLD");
  //const Float_t ZZMEKDLD = cand.userFloat("MEKD_LD");
  //const Float_t ZZMEKDpseudoLD = cand.userFloat("MEKD_PseudoLD");
  //const Float_t ZZMEKDgravLD = cand.userFloat("MEKD_GravLD");

  //const Float_t p0plus_melaNorm = cand.userFloat("p0plus_melaNorm");
  //const Float_t p0plus_mela = cand.userFloat("p0plus_mela");
  //const Float_t p0minus_mela = cand.userFloat("p0minus_mela");
  //const Float_t p0hplus_mela = cand.userFloat("p0hplus_mela"); // 0h+, analytic distribution
  const Float_t p0plus_VAJHU = cand.userFloat("p0plus_VAJHU");
  const Float_t p0minus_VAJHU = cand.userFloat("p0minus_VAJHU");
  const Float_t p0plus_VAMCFM = cand.userFloat("p0plus_VAMCFM");
  const Float_t p0hplus_VAJHU = cand.userFloat("p0hplus_VAJHU"); // 0h+ (high dimensional operator), vector algebra, JHUgen
  //const Float_t p1_mela = cand.userFloat("p1_mela");
  //const Float_t p1_prodIndep_mela = cand.userFloat("p1_prodIndep_mela");
  //const Float_t p1plus_mela = cand.userFloat("p1plus_mela"); // 1+, analytic distribution 
  //const Float_t p1plus_prodIndep_mela = cand.userFloat("p1plus_prodIndep_mela"); // 1+, analytic distribution 
  const Float_t p1_VAJHU = cand.userFloat("p1_VAJHU");
  const Float_t p1_prodIndep_VAJHU = cand.userFloat("p1_prodIndep_VAJHU");
  const Float_t p1plus_VAJHU = cand.userFloat("p1plus_VAJHU"); // 1+ (axial vector), vector algebra, JHUgen,
  const Float_t p1plus_prodIndep_VAJHU = cand.userFloat("p1plus_prodIndep_VAJHU"); // 1+ (axial vector), vector algebra, JHUgen,
  //const Float_t p2_mela  = cand.userFloat("p2_mela");
  //const Float_t p2_prodIndep_mela  = cand.userFloat("p2_prodIndep_mela");
  //const Float_t p2qqb_mela = cand.userFloat("p2qqb_mela"); // graviton produced by qqbar vector algebra, analytical,
  //const Float_t p2hplus_mela = cand.userFloat("p2hplus_mela"); // graviton produced by qqbar vector algebra, analytical,
  //const Float_t p2hminus_mela = cand.userFloat("p2hminus_mela"); // graviton produced by qqbar vector algebra, analytical,
  //const Float_t p2bplus_mela = cand.userFloat("p2bplus_mela"); // graviton produced by qqbar vector algebra, analytical,
  const Float_t p2_VAJHU = cand.userFloat("p2_VAJHU");
  const Float_t p2_prodIndep_VAJHU = cand.userFloat("p2_prodIndep_VAJHU");
  const Float_t p2qqb_VAJHU = cand.userFloat("p2qqb_VAJHU");
  const Float_t p2hplus_VAJHU = cand.userFloat("p2hplus_VAJHU");
  const Float_t p2hminus_VAJHU = cand.userFloat("p2hminus_VAJHU");
  const Float_t p2bplus_VAJHU = cand.userFloat("p2bplus_VAJHU");
  const Float_t p2hplus_qqb_VAJHU= cand.userFloat(					"p2hplus_qqb_VAJHU");					
  const Float_t p2hplus_prodIndep_VAJHU= cand.userFloat(		"p2hplus_prodIndep_VAJHU");	
  const Float_t p2hminus_qqb_VAJHU= cand.userFloat(				"p2hminus_qqb_VAJHU");				
  const Float_t p2hminus_prodIndep_VAJHU= cand.userFloat(	"p2hminus_prodIndep_VAJHU");	
  const Float_t p2bplus_qqb_VAJHU= cand.userFloat(					"p2bplus_qqb_VAJHU"				);	
  const Float_t p2bplus_prodIndep_VAJHU= cand.userFloat(		"p2bplus_prodIndep_VAJHU"		);
  const Float_t p2h2plus_gg_VAJHU= cand.userFloat(      		"p2h2plus_gg_VAJHU"      		);
  const Float_t p2h2plus_qqbar_VAJHU= cand.userFloat(   		"p2h2plus_qqbar_VAJHU"   		);
  const Float_t p2h2plus_prodIndep_VAJHU= cand.userFloat(	"p2h2plus_prodIndep_VAJHU");	
  const Float_t p2h3plus_gg_VAJHU= cand.userFloat(       	"p2h3plus_gg_VAJHU"       );	
  const Float_t p2h3plus_qqbar_VAJHU= cand.userFloat(    	"p2h3plus_qqbar_VAJHU"    );	
  const Float_t p2h3plus_prodIndep_VAJHU= cand.userFloat(	"p2h3plus_prodIndep_VAJHU");	
  const Float_t p2h6plus_gg_VAJHU= cand.userFloat(       	"p2h6plus_gg_VAJHU"       );	
  const Float_t p2h6plus_qqbar_VAJHU= cand.userFloat(    	"p2h6plus_qqbar_VAJHU"    );	
  const Float_t p2h6plus_prodIndep_VAJHU= cand.userFloat(	"p2h6plus_prodIndep_VAJHU");	
  const Float_t p2h7plus_gg_VAJHU= cand.userFloat(       	"p2h7plus_gg_VAJHU"       );	
  const Float_t p2h7plus_qqbar_VAJHU= cand.userFloat(    	"p2h7plus_qqbar_VAJHU"    );	
  const Float_t p2h7plus_prodIndep_VAJHU= cand.userFloat(	"p2h7plus_prodIndep_VAJHU");	
  const Float_t p2h9minus_gg_VAJHU= cand.userFloat(       	"p2h9minus_gg_VAJHU"       	);
  const Float_t p2h9minus_qqbar_VAJHU= cand.userFloat(    	"p2h9minus_qqbar_VAJHU"    	);
  const Float_t p2h9minus_prodIndep_VAJHU= cand.userFloat(	"p2h9minus_prodIndep_VAJHU"	);
  const Float_t p2h10minus_gg_VAJHU= cand.userFloat(       "p2h10minus_gg_VAJHU"      ); 
  const Float_t p2h10minus_qqbar_VAJHU= cand.userFloat(    "p2h10minus_qqbar_VAJHU"     ); 
  const Float_t p2h10minus_prodIndep_VAJHU= cand.userFloat("p2h10minus_prodIndep_VAJHU" ); 
//const Float_t bkg_mela = cand.userFloat("bkg_mela");
	const Float_t bkg_VAMCFM = cand.userFloat("bkg_VAMCFM");
  const Float_t bkg_prodIndep_VAMCFM = cand.userFloat("bkg_prodIndep_VAMCFM");
  const Float_t ggzz_VAMCFM = cand.userFloat("ggzz_VAMCFM");
  const Float_t ggzz_p0plus_VAMCFM = cand.userFloat("ggzz_p0plus_VAMCFM");
  const Float_t ggzz_c1_VAMCFM = cand.userFloat("ggzz_c1_VAMCFM");
  const Float_t ggzz_c5_VAMCFM = cand.userFloat("ggzz_c5_VAMCFM");
  const Float_t ggzz_ci_VAMCFM = cand.userFloat("ggzz_ci_VAMCFM");
  //const Float_t bkg_VAMCFMNorm = cand.userFloat("bkg_VAMCFMNorm");
  //const Float_t p0_pt = cand.userFloat("p0_pt");
  //const Float_t p0_y = cand.userFloat("p0_y");
  //const Float_t bkg_pt = cand.userFloat("bkg_pt");
  //const Float_t bkg_y = cand.userFloat("bkg_y");

  const Float_t p0plus_m4l = cand.userFloat("p0plus_m4l");
  const Float_t bkg_m4l = cand.userFloat("bkg_m4l");
  const Float_t pg1g4_mela = cand.userFloat("pg1g4_mela");
  const Float_t pg1g4_VAJHU = cand.userFloat("pg1g4_VAJHU");
  const Float_t pg1g4_pi2_VAJHU = cand.userFloat("pg1g4_pi2_VAJHU");
  const Float_t pg1g2_pi2_VAJHU = cand.userFloat("pg1g2_pi2_VAJHU");
  const Float_t pg1g2_mela = cand.userFloat("pg1g2_mela");
  const Float_t pg1g2_VAJHU = cand.userFloat("pg1g2_VAJHU");
  const Float_t p0plus_m4l_ScaleUp = cand.userFloat("p0plus_m4l_ScaleUp");// signal m4l probability for systematics
  const Float_t bkg_m4l_ScaleUp = cand.userFloat("bkg_m4l_ScaleUp");// backgroun m4l probability for systematics
  const Float_t p0plus_m4l_ScaleDown = cand.userFloat("p0plus_m4l_ScaleDown");// signal m4l probability for systematics
  const Float_t bkg_m4l_ScaleDown = cand.userFloat("bkg_m4l_ScaleDown");// backgroun m4l probability for systematics
  const Float_t p0plus_m4l_ResUp = cand.userFloat("p0plus_m4l_ResUp");// signal m4l probability for systematics
  const Float_t bkg_m4l_ResUp = cand.userFloat("bkg_m4l_ResUp");// backgroun m4l probability for systematics
  const Float_t p0plus_m4l_ResDown = cand.userFloat("p0plus_m4l_ResDown");// signal m4l probability for systematics
  const Float_t bkg_m4l_ResDown = cand.userFloat("bkg_m4l_ResDown");// backgroun m4l probability for systematics

  const Float_t phjj_VAJHU_old = cand.userFloat("phjj_VAJHU_old");
  const Float_t pvbf_VAJHU_old = cand.userFloat("pvbf_VAJHU_old");
  const Float_t phjj_VAJHU_old_up = cand.userFloat("phjj_VAJHU_old_up");
  const Float_t pvbf_VAJHU_old_up = cand.userFloat("pvbf_VAJHU_old_up");
  const Float_t phjj_VAJHU_old_dn = cand.userFloat("phjj_VAJHU_old_dn");
  const Float_t pvbf_VAJHU_old_dn = cand.userFloat("pvbf_VAJHU_old_dn");
  const Float_t phjj_VAJHU_new = cand.userFloat("phjj_VAJHU_new");
  const Float_t pvbf_VAJHU_new = cand.userFloat("pvbf_VAJHU_new");
  const Float_t phjj_VAJHU_new_up = cand.userFloat("phjj_VAJHU_new_up");
  const Float_t pvbf_VAJHU_new_up = cand.userFloat("pvbf_VAJHU_new_up");
  const Float_t phjj_VAJHU_new_dn = cand.userFloat("phjj_VAJHU_new_dn");
  const Float_t pvbf_VAJHU_new_dn = cand.userFloat("pvbf_VAJHU_new_dn");
  const Float_t p0_g1prime2_VAJHU= cand.userFloat("p0_g1prime2_VAJHU");
  const Float_t pg1g1prime2_VAJHU= cand.userFloat("pg1g1prime2_VAJHU");
  const Float_t Dgg10_VAMCFM= cand.userFloat("Dgg10_VAMCFM");
  const Float_t pzzzg_VAJHU= cand.userFloat(    "pzzzg_VAJHU");
  const Float_t pzzgg_VAJHU= cand.userFloat(    "pzzgg_VAJHU");
  const Float_t pzzzg_PS_VAJHU= cand.userFloat( "pzzzg_PS_VAJHU");
  const Float_t pzzgg_PS_VAJHU= cand.userFloat( "pzzgg_PS_VAJHU");
  const Float_t p0Zgs_VAJHU= cand.userFloat(    "p0Zgs_VAJHU");
  const Float_t p0gsgs_VAJHU= cand.userFloat(   "p0gsgs_VAJHU");
  const Float_t p0Zgs_PS_VAJHU= cand.userFloat( "p0Zgs_PS_VAJHU");
  const Float_t p0gsgs_PS_VAJHU= cand.userFloat("p0gsgs_PS_VAJHU");

  //Int_t isSignal = -1;
  //Int_t isRightPair = -1;
  //if(isMC){
  //  isSignal = cand.userFloat("MC_isRight");
  //  isRightPair = cand.userFloat("MC_isRightPair");
  //}

  //const Float_t mZa = cand.userFloat("mZa");
  //const Float_t mZb = cand.userFloat("mZb");
  //const Float_t mLL4 = cand.userFloat("mLL4");
  //const Float_t mLL6 = cand.userFloat("mLL6");
  //const Float_t SIP4 = cand.userFloat("SIP4");
  //const Float_t iso34 = cand.userFloat("iso34");

  //Z1 and Z2 variables
  const reco::Candidate* Z1;
  const reco::Candidate* Z2;
  vector<const reco::Candidate*> leptons;
  vector<string> labels;

  if (theChannel!=ZL) { // Regular 4l candidates
    Z1   = cand.daughter("Z1");
    Z2   = cand.daughter("Z2");
    userdatahelpers::getSortedLeptons(cand, leptons, labels);
  } else {              // Special handling of Z+l candidates 
    Z1   = cand.daughter(0); // the Z
    Z2   = cand.daughter(1); // This is actually the additional lepton!
    userdatahelpers::getSortedLeptons(cand, leptons, labels, false); // note: we get just 3 leptons in this case.
  }

  Float_t Z1Mass = Z1->mass();
  Float_t Z1Pt =   Z1->pt();
  Float_t Z1MassRefit = cand.userFloat("mZ1Ref");
  short Z1Flav = Z1->daughter(0)->pdgId()*Z1->daughter(1)->pdgId();
  
  Float_t Z2Mass = Z2->mass();
  Float_t Z2Pt =   Z2->pt();
  short Z2Flav = Z2->daughter(0)->pdgId()*Z2->daughter(1)->pdgId();

  // Precomputed selections
  bool candPass70Z2Loose   = cand.userFloat("Z2Mass") && 
                             cand.userFloat("MAllComb") &&
                             cand.userFloat("pt1")>20 && cand.userFloat("pt2")>10. &&
                             ZZMass>70.;
  bool candPassFullSel70   = cand.userFloat("FullSel70");
  bool candPassFullSel   = cand.userFloat("FullSel");
  bool candIsBest = cand.userFloat("isBestCand");
  bool passMz_zz = (Z1Mass>60. && Z1Mass<120. && Z2Mass>60. && Z2Mass<120.);   //FIXME hardcoded cut

  Int_t sel = 0;

//   if(candPassFullSel){ // FIXME: this step (fullSel, without candIsBest) is not part of the normal selection flow.
//     sel = 50;
//   }
  
  if (candIsBest) {
    //    sel = 10; //FIXME see above
    if (candPass70Z2Loose) sel=70;
    if (candPassFullSel70){ // includes MZ2 > 12
      sel = 90;
      if (candPassFullSel){
	sel=100;
	if (passMz_zz) sel = 120;
      }
    }
  }
  if (!(evtPass)) {sel = -sel;} // avoid confusion when we write events which do not pass trigger/skim

//cout<<"KD= "<<p0plus_VAJHU+p0minus_VAJHU<<endl;
  //myTree->FillProbability(
  double probabilities[84]={ p0plus_VAJHU,
                          p0minus_VAJHU,
                          p0plus_VAMCFM,
                          p0hplus_VAJHU, // 0h+ (high dimensional operator), vector algebra, JHUgen
                          p1_VAJHU,
                          p1_prodIndep_VAJHU,
                          p1plus_VAJHU, // 1+ (axial vector), vector algebra, JHUgen,
                          p1plus_prodIndep_VAJHU, // 1+ (axial vector), vector algebra, JHUgen,
                          p2_VAJHU ,
                          p2_prodIndep_VAJHU ,
                          p2qqb_VAJHU,
                          p2hplus_VAJHU,
                          p2hminus_VAJHU,
                          p2bplus_VAJHU,
                          p2hplus_qqb_VAJHU,
                          p2hplus_prodIndep_VAJHU,
                          p2hminus_qqb_VAJHU,
                          p2hminus_prodIndep_VAJHU,
                          p2bplus_qqb_VAJHU,
                          p2bplus_prodIndep_VAJHU,
                          p2h2plus_gg_VAJHU,
                          p2h2plus_qqbar_VAJHU,
                          p2h2plus_prodIndep_VAJHU,
                          p2h3plus_gg_VAJHU,
                          p2h3plus_qqbar_VAJHU,
                          p2h3plus_prodIndep_VAJHU,
                          p2h6plus_gg_VAJHU,
                          p2h6plus_qqbar_VAJHU,
                          p2h6plus_prodIndep_VAJHU,
                          p2h7plus_gg_VAJHU,
                          p2h7plus_qqbar_VAJHU,
                          p2h7plus_prodIndep_VAJHU,
                          p2h9minus_gg_VAJHU,
                          p2h9minus_qqbar_VAJHU,
                          p2h9minus_prodIndep_VAJHU,
                          p2h10minus_gg_VAJHU,
                          p2h10minus_qqbar_VAJHU,
                          p2h10minus_prodIndep_VAJHU,
                          bkg_VAMCFM,
                          bkg_prodIndep_VAMCFM,
                          ggzz_VAMCFM,
                          ggzz_p0plus_VAMCFM,
                          ggzz_c1_VAMCFM,
                          ggzz_c5_VAMCFM,
                          ggzz_ci_VAMCFM,
                          phjj_VAJHU_old,
                          pvbf_VAJHU_old,
                          phjj_VAJHU_old_up,
                          pvbf_VAJHU_old_up,
                          phjj_VAJHU_old_dn,
                          pvbf_VAJHU_old_dn,
                          phjj_VAJHU_new,
                          pvbf_VAJHU_new,
                          phjj_VAJHU_new_up,
                          pvbf_VAJHU_new_up,
                          phjj_VAJHU_new_dn,
                          pvbf_VAJHU_new_dn,
                          p0_g1prime2_VAJHU,
                          pg1g1prime2_VAJHU,
                          Dgg10_VAMCFM,
                          pg1g4_mela,
                          pg1g4_VAJHU,
                          pg1g4_pi2_VAJHU,
                          pg1g2_pi2_VAJHU,
                          pg1g2_mela,
                          pg1g2_VAJHU,
                          pzzzg_VAJHU,
                          pzzgg_VAJHU,
                          pzzzg_PS_VAJHU,
                          pzzgg_PS_VAJHU,
                          p0Zgs_VAJHU,
                          p0gsgs_VAJHU,
                          p0Zgs_PS_VAJHU,
                          p0gsgs_PS_VAJHU,
                          			 bkg_m4l,   // backgroun m4l probability as in datacards
			 p0plus_m4l_ScaleUp,  // signal m4l probability for systematics
			 bkg_m4l_ScaleUp,     // backgroun m4l probability for systematics
			 p0plus_m4l_ScaleDown,  // signal m4l probability for systematics
			 bkg_m4l_ScaleDown,     // backgroun m4l probability for systematics
			 p0plus_m4l_ResUp,  // signal m4l probability for systematics
			 bkg_m4l_ResUp,     // backgroun m4l probability for systematics
			 p0plus_m4l_ResDown,  // signal m4l probability for systematics
			 bkg_m4l_ResDown,     // backgroun m4l probability for systematics
			 p0plus_m4l
                          };
    TString probNames[84]={ "p0plus_VAJHU",
                          "p0minus_VAJHU",
                          "p0plus_VAMCFM",
                          "p0hplus_VAJHU", // 0h+ (high dimensional operator)", vector algebra", JHUgen
                          "p1_VAJHU",
                          "p1_prodIndep_VAJHU",
                          "p1plus_VAJHU", // 1+ (axial vector)", vector algebra", JHUgen",
                          "p1plus_prodIndep_VAJHU", // 1+ (axial vector)", vector algebra", JHUgen",
                          "p2_VAJHU",
                          "p2_prodIndep_VAJHU",
                          "p2qqb_VAJHU",
                          "p2hplus_VAJHU",
                          "p2hminus_VAJHU",
                          "p2bplus_VAJHU",
                          "p2hplus_qqb_VAJHU",
                          "p2hplus_prodIndep_VAJHU",
                          "p2hminus_qqb_VAJHU",
                          "p2hminus_prodIndep_VAJHU",
                          "p2bplus_qqb_VAJHU",
                          "p2bplus_prodIndep_VAJHU",
                          "p2h2plus_gg_VAJHU",
                          "p2h2plus_qqbar_VAJHU",
                          "p2h2plus_prodIndep_VAJHU",
                          "p2h3plus_gg_VAJHU",
                          "p2h3plus_qqbar_VAJHU",
                          "p2h3plus_prodIndep_VAJHU",
                          "p2h6plus_gg_VAJHU",
                          "p2h6plus_qqbar_VAJHU",
                          "p2h6plus_prodIndep_VAJHU",
                          "p2h7plus_gg_VAJHU",
                          "p2h7plus_qqbar_VAJHU",
                          "p2h7plus_prodIndep_VAJHU",
                          "p2h9minus_gg_VAJHU",
                          "p2h9minus_qqbar_VAJHU",
                          "p2h9minus_prodIndep_VAJHU",
                          "p2h10minus_gg_VAJHU",
                          "p2h10minus_qqbar_VAJHU",
                          "p2h10minus_prodIndep_VAJHU",
                          "bkg_VAMCFM",
                          "bkg_prodIndep_VAMCFM",
                          "ggzz_VAMCFM",
                          "ggzz_p0plus_VAMCFM",
                          "ggzz_c1_VAMCFM",
                          "ggzz_c5_VAMCFM",
                          "ggzz_ci_VAMCFM",
                          "phjj_VAJHU_old",
                          "pvbf_VAJHU_old",
                          "phjj_VAJHU_old_up",
                          "pvbf_VAJHU_old_up",
                          "phjj_VAJHU_old_dn",
                          "pvbf_VAJHU_old_dn",
                          "phjj_VAJHU_new",
                          "pvbf_VAJHU_new",
                          "phjj_VAJHU_new_up",
                          "pvbf_VAJHU_new_up",
                          "phjj_VAJHU_new_dn",
                          "pvbf_VAJHU_new_dn",
                          "p0_g1prime2_VAJHU",
                          "pg1g1prime2_VAJHU",
                          "Dgg10_VAMCFM",
                          "pg1g4_mela",
                          "pg1g4_VAJHU",
                          "pg1g4_pi2_VAJHU",
                          "pg1g2_pi2_VAJHU",
                          "pg1g2_mela",
                          "pg1g2_VAJHU",
                          "pzzzg_VAJHU",
                          "pzzgg_VAJHU",
                          "pzzzg_PS_VAJHU",
                          "pzzgg_PS_VAJHU",
                          "p0Zgs_VAJHU",
                          "p0gsgs_VAJHU",
                          "p0Zgs_PS_VAJHU",
                          "p0gsgs_PS_VAJHU",
                          "bkg_m4l",   // backgroun m4l probability as in datacards
			                    "p0plus_m4l_ScaleUp",  // signal m4l probability for systematics
			                    "bkg_m4l_ScaleUp",     // backgroun m4l probability for systematics
			                    "p0plus_m4l_ScaleDown",  // signal m4l probability for systematics
			                    "bkg_m4l_ScaleDown",     // backgroun m4l probability for systematics
			                    "p0plus_m4l_ResUp",  // signal m4l probability for systematics
			                    "bkg_m4l_ResUp",     // backgroun m4l probability for systematics
			                    "p0plus_m4l_ResDown",  // signal m4l probability for systematics
			                    "bkg_m4l_ResDown",     // backgroun m4l probability for systematics
			                    "p0plus_m4l" // signal m4l probability as in datacards
                          };
  myTree->SetVariables((TString *)probNames,(double *)probabilities,84);
	
	//FIXME? Not filled in the miniAOD branch, is it right?
	/*
	float Hadditional[6]={mZa, mZb, mLL4, mLL6, SIP4, iso34};		    
  TString addName[6]={};
  myTree->SetVariables(addName,Hadditional,6);
  */
  myTree->FillZInfo(Z1Mass, Z1Pt, Z1Flav, Z1MassRefit);
  myTree->FillZInfo(Z2Mass, Z2Pt, Z2Flav, Z2Mass);

  //Fill the angular variables
  Float_t costheta1 = cand.userFloat("costheta1");
  Float_t costheta2 = cand.userFloat("costheta2");
  Float_t phi       = cand.userFloat("phi");
  Float_t costhetastar = cand.userFloat("costhetastar");
  Float_t phistar1      = cand.userFloat("phistar1");
  Float_t phistar2      = cand.userFloat("phistar2");
  Float_t xi            = cand.userFloat("xi");
  Float_t xistar        = cand.userFloat("xistar");
  double angularinfo[8]={costhetastar, phi, costheta1, costheta2, phistar1, phistar2,xi,xistar};
  TString angleName[8]={
    "costhetastar",
    "helphi",
    "helcosthetaZ1",
    "helcosthetaZ2",
    "phistarZ1",
    "phistarZ2",
    "xi",
    "xistar"};
  myTree->SetVariables((TString *)angleName,(double *)angularinfo,8);

  // Retrieve the userFloat of the leptons in vectors ordered in the same way.
  vector<float> SIP(4);
  vector<float> PFChargedHadIso(4);
  vector<float> PFNeutralHadIso(4);
  vector<float> PFPhotonIso(4);
  vector<float> combRelIsoPF(4);
  vector<bool>  isID(4);
  
  for (unsigned int i=0; i<leptons.size(); ++i){
    SIP[i]             = userdatahelpers::getUserFloat(leptons[i],"SIP");
    PFChargedHadIso[i] = userdatahelpers::getUserFloat(leptons[i],"PFChargedHadIso");
    PFNeutralHadIso[i] = userdatahelpers::getUserFloat(leptons[i],"PFNeutralHadIso");
    PFPhotonIso[i]     = userdatahelpers::getUserFloat(leptons[i],"PFPhotonIso");
    isID[i]            = userdatahelpers::getUserFloat(leptons[i],"ID");

    if (theChannel==ZL) {
      combRelIsoPF[i]    = userdatahelpers::getUserFloat(leptons[i],"combRelIsoPF");
      //FIXME cannot take labels[i]+"SIP", that info only attached to the Z!!
    } else {
      combRelIsoPF[i]    = cand.userFloat(labels[i]+"combRelIsoPFFSRCorr"); // Note: the FSR-corrected iso is attached to the Z, not to the lepton!
      // Check that I don't mess up with labels[] and leptons[]
      //cout<<"SIP: "<<SIP[i]<<"  labels "<<labels[i]<<cand.userFloat(labels[i]+"SIP")<<endl;
      SIP[i] = cand.userFloat(labels[i]+"SIP");
    }

    //Fill the info on the lepton candidates  
    TString lepnames[12];
    TString outlepnames[12]={"Pt","Eta","Phi", "LepId", "SIP", "isID", "BDT","missingHit","chargedHadIso","neutralHadIso","photonIso","combRelIsoPF"};
    for(int ilep=0;ilep<12;ilep++)lepnames[ilep].Form("Lep%d%s",i+1,outlepnames[ilep].Data());
    double lepVariables[]={(double)leptons[i]->pt(),
			(double)leptons[i]->eta(),
			(double)leptons[i]->phi(),
			(double)leptons[i]->pdgId(),
			(double)SIP[i],
			(double)isID[i],
			(double)userdatahelpers::getUserFloat(leptons[i],"BDT"),
			//userdatahelpers::getUserFloat(leptons[i],"MCParentCode"), // This is dummy as of now.
			(double)userdatahelpers::getUserFloat(leptons[i],"missingHit"),
			(double)PFChargedHadIso[i],
			(double)PFNeutralHadIso[i],
			(double)PFPhotonIso[i], 
			(double)combRelIsoPF[i]
			};
    myTree->SetVariables((TString *)lepnames,(double *)lepVariables,8);
  }
  
  //convention: 0 -> 4mu   1 -> 4e   2 -> 2mu2e
  float zxw = 1.0;
  if(CRflag){
    for(int izx=0;izx<2;izx++)
      zxw *= getFakeWeight(Z2->daughter(izx)->pt(),Z2->daughter(izx)->eta(),Z2->daughter(izx)->pdgId(),Z1->daughter(0)->pdgId());
  }
  
  double hinfo[14]={ZZMass, ZZMassErr, ZZMassErrCorr, ZZMassPreFSR, ZZMassRefit, (double)sel, ZZPt, ZZEta, ZZPhi, (double)CRflag,
    Chi2KinFit, ZZMassCFit, Chi2CFit, zxw};
	TString hinfonames[14]={
	"ZZMass",
"ZZMassErr",
"ZZMassErrCorr",
"ZZMassPreFSR",
"ZZMassRefit",
"ZZsel",
"ZZPt",
"ZZEta",
"ZZPhi",
"CRflag",
"ZZChi2KinFit",
"ZZMassCFit",
"ZZChi2CFit",
"ZXFakeweight"
	};
  myTree->SetVariables((TString *)hinfonames,(double *)hinfo,14);

  //Fill the info on categorization
  myTree->SetVariable("nExtraLep",cand.userFloat("nExtraLep"));
  myTree->SetVariable("nExtraZ",cand.userFloat("nExtraZ"));

  //Fill the info on the extra leptons
  for(int iExtraLep=1;iExtraLep<4;iExtraLep++){
    TString extraString;extraString.Form("ExtraLep%d",iExtraLep);
    if(cand.hasUserCand(extraString.Data())){
      TString varExtra[4]={"Pt","Eta","Phi","LepId"};
      for(int iextra=0;iextra<4;iextra++)varExtra[iextra].Prepend(extraString.Data());
      reco::CandidatePtr candPtr=cand.userCand(extraString.Data());
      double valuesExtra[4]={candPtr->pt(),candPtr->eta(),candPtr->phi(),(double)candPtr->pdgId()};
      myTree->SetVariables((TString *)varExtra,(double *)valuesExtra,4);
    }
  }
  //myTree->FillExtraLepInfo( 1, cand.hasUserCand("ExtraLep1"), (cand.hasUserCand("ExtraLep1") ? cand.userCand("ExtraLep1") : *(new reco::CandidatePtr)) );
  //myTree->FillExtraLepInfo( 2, cand.hasUserCand("ExtraLep2"), (cand.hasUserCand("ExtraLep2") ? cand.userCand("ExtraLep2") : *(new reco::CandidatePtr)) );
  //myTree->FillExtraLepInfo( 3, cand.hasUserCand("ExtraLep3"), (cand.hasUserCand("ExtraLep3") ? cand.userCand("ExtraLep3") : *(new reco::CandidatePtr)) );

}



// ------------ method called once each job just before starting event loop  ------------
void HZZ4lNtupleMaker::beginJob()
{
  edm::Service<TFileService> fs;
  myTree = new HZZ4lNtupleFactory( fs->make<TTree>(theFileName,"Event Summary"));
  hCounter = fs->make<TH1F>("Counters", "Counters", 45, 0., 45.);
  BookAllBranches();
}

// ------------ method called once each job just after ending the event loop  ------------
void HZZ4lNtupleMaker::endJob()
{
  hCounter->SetBinContent(0 ,gen_sumWeights); // also stored in bin 40
  hCounter->SetBinContent(1 ,Nevt_Gen-gen_BUGGY);
  hCounter->SetBinContent(2 ,gen_ZZ4mu);
  hCounter->SetBinContent(3 ,gen_ZZ4e);
  hCounter->SetBinContent(4 ,gen_ZZ2mu2e);
  hCounter->SetBinContent(5 ,gen_ZZ2l2tau);
  hCounter->SetBinContent(6 ,gen_ZZ4mu_EtaAcceptance);
  hCounter->SetBinContent(7 ,gen_ZZ4mu_LeptonAcceptance);
  hCounter->SetBinContent(10,gen_ZZ4e_EtaAcceptance);
  hCounter->SetBinContent(11,gen_ZZ4e_LeptonAcceptance);
  hCounter->SetBinContent(14,gen_ZZ2mu2e_EtaAcceptance);
  hCounter->SetBinContent(15,gen_ZZ2mu2e_LeptonAcceptance);
  hCounter->SetBinContent(19,gen_BUGGY);
  hCounter->SetBinContent(20,gen_Unknown);

  hCounter->SetBinContent(40,gen_sumWeights); // Also stored in underflow bin; added here for convenience
  hCounter->SetBinContent(41,gen_sumGenMCWeight);
  hCounter->SetBinContent(42,gen_sumPUWeight);


  hCounter->GetXaxis()->SetBinLabel(1 ,"Nevt_Gen");
  hCounter->GetXaxis()->SetBinLabel(2 ,"gen_ZZ4mu");
  hCounter->GetXaxis()->SetBinLabel(3 ,"gen_ZZ4e");
  hCounter->GetXaxis()->SetBinLabel(4 ,"gen_ZZ2mu2e");
  hCounter->GetXaxis()->SetBinLabel(5 ,"gen_ZZ2l2tau");
  hCounter->GetXaxis()->SetBinLabel(6 ,"gen_ZZ4mu_EtaAcceptance");
  hCounter->GetXaxis()->SetBinLabel(7 ,"gen_ZZ4mu_LeptonAcceptance");
  hCounter->GetXaxis()->SetBinLabel(10,"gen_ZZ4e_EtaAcceptance");
  hCounter->GetXaxis()->SetBinLabel(11,"gen_ZZ4e_LeptonAcceptance");
  hCounter->GetXaxis()->SetBinLabel(14,"gen_ZZ2mu2e_EtaAcceptance");
  hCounter->GetXaxis()->SetBinLabel(15,"gen_ZZ2mu2e_LeptonAcceptance");
  hCounter->GetXaxis()->SetBinLabel(19,"gen_BUGGY");
  hCounter->GetXaxis()->SetBinLabel(20,"gen_Unknown");

  hCounter->GetXaxis()->SetBinLabel(40,"gen_sumWeights");
  hCounter->GetXaxis()->SetBinLabel(41,"gen_sumGenMCWeight");
  hCounter->GetXaxis()->SetBinLabel(42,"gen_sumPUWeight");

  return;
}

// ------------ method called when starting to processes a run  ------------
void HZZ4lNtupleMaker::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void HZZ4lNtupleMaker::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void HZZ4lNtupleMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
  Nevt_Gen_lumiBlock = 0;
}

// ------------ method called when ending the processing of a luminosity block  ------------
void HZZ4lNtupleMaker::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup)
{
  Float_t Nevt_preskim = -1.;
  edm::Handle<edm::MergeableCounter> preSkimCounter;
  if (iLumi.getByLabel("preSkimCounter", preSkimCounter)) { // Counter before skim. Does not exist for non-skimmed samples.
    Nevt_preskim = preSkimCounter->value;
  }  

  // Nevt_gen: this is the number before any skim
  if (Nevt_preskim>=0.) {
    Nevt_Gen += Nevt_preskim; 
  } else {
    Nevt_Gen += Nevt_Gen_lumiBlock ;
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HZZ4lNtupleMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


Float_t HZZ4lNtupleMaker::getAllWeight(const Float_t LepPt, const Float_t LepEta, Int_t LepID) const
{
  Float_t weight  = 1.; 
  //Float_t errCorr = 0.;
  //Float_t errCorrSyst = 0.;

  Float_t myLepPt = LepPt;
  Float_t myLepEta = LepEta;
  Int_t   myLepID = abs(LepID);
  
  //avoid to go out of the TH boundary
  if(myLepID == 13 && myLepPt > 99.) myLepPt = 99.;
  if(myLepID == 11 && myLepPt > 199.) myLepPt = 199.;
  if(myLepID == 11) myLepEta = fabs(myLepEta);

  if(myLepID == 13){                                               
      weight  = hTH2D_Mu_All->GetBinContent(hTH2D_Mu_All->GetXaxis()->FindBin(myLepPt),hTH2D_Mu_All->GetYaxis()->FindBin(LepEta));
    //  errCorr = hTH2D_Mu_All_2011A->GetBinError(hTH2D_Mu_All_2011A->GetXaxis()->FindBin(myLepPt),hTH2D_Mu_All_2011A->GetYaxis()->FindBin(LepEta));
      
  }
  else if(myLepID == 11){   

    weight  = hTH2D_El_All->GetBinContent(hTH2D_El_All->GetXaxis()->FindBin(myLepPt),hTH2D_El_All->GetYaxis()->FindBin(myLepEta));
      //errCorr = hTH2D_El_All_2011A->GetBinError(hTH2D_El_All_2011A->GetXaxis()->FindBin(myLepPt),hTH2D_El_All_2011A->GetYaxis()->FindBin(myLepEta));   
  }else {
    cout<<"ERROR! wrong lepton ID "<<myLepID<<endl;
    //abort();
    weight=0;
  }   

  
  //add the systematics on T&P corrections (for muons only, electrons have them already included)
  //if(myLepID == 13){
  //  if(myLepPt >= 15.) errCorrSyst = 0.005;
  //  else errCorrSyst = 0.015;
  //}

  //FIXME
  if(myLepPt < 5. && myLepID == 13) weight = 1.;

  if(weight < 0.001 || weight > 10.){
    cout << "ERROR! LEP out of range! myLepPt = " << myLepPt << " myLepEta = " << myLepEta <<" myLepID "<<myLepID<< " weight = " << weight << endl;
    //abort();  //no correction should be zero, if you find one, stop
  }
/*
//FIXME: in HZZ4l was working because seeder was the event number (1=first event, 10=10th event)
//Here I don't have this info. Random seed?
  static TRandom3 randomToss;
  if( ( myLepID == 13) || (myLepID == 11) ){

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
*/
  return weight;
}

Float_t HZZ4lNtupleMaker::getHqTWeight(double mH, double genPt) const
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
  
  if (iMass>=0) {
    weight = h_weight->GetBinContent(h_weight->FindBin(genPt));
  }
  return weight;
}


// Added by CO
Float_t HZZ4lNtupleMaker::getFakeWeight(const Float_t LepPt, const Float_t LepEta, Int_t LepID, Int_t LepZ1ID)
{
  // year 0 = 2011
  // year 1 = 2012

  Float_t weight  = 1.; 
  
  Int_t nHisto=0;
  
  Float_t myLepPt   = LepPt;
  Float_t myLepEta  = fabs(LepEta);
  Int_t   myLepID   = abs(LepID);
  Int_t   myZ1LepID = abs(LepZ1ID);

  //cout << " pt = " << myLepPt << " eta = " << myLepEta << " ZID = " << myZ1LepID << " LepID = " << myLepID << endl;

  //avoid to go out of the TH boundary
  if(myLepPt > 79.) myLepPt = 79.;
  if(myLepID==13)nHisto+=2;
  if(myZ1LepID==13)nHisto+=1;
 
 
  TString Z1flavor = "Z1ee";     if(myZ1LepID==13) Z1flavor = "Z1mumu";
  TString Z2flavor = "electron"; if(myLepID==13)   Z2flavor = "muon";
  TString histo_name = "eff_"+Z1flavor+"_plus_"+Z2flavor;
  
  //cout << " histo = " << histo_name << endl;

  weight = h_ZXWeight[nHisto]->GetBinContent(h_ZXWeight[nHisto]->GetXaxis()->FindBin(myLepPt), h_ZXWeight[nHisto]->GetYaxis()->FindBin(myLepEta));
  /*
  h_ZXWeight[0]=FileZXWeightEle->Get("eff_Z1ee_plus_electron")->Clone();
  h_ZXWeight[2]=FileZXWeightEle->Get("eff_Z1mumu_plus_electron")->Clone();
  
  h_ZXWeight[5]=FileZXWeightMuo->Get("eff_Z1ee_plus_muon")->Clone();
  h_ZXWeight[7]=FileZXWeightMuo->Get("eff_Z1mumu_plus_muon")->Clone();
*/
  return weight;

} // end of getFakeWeight

void HZZ4lNtupleMaker::BookAllBranches(){
   //Event variables
  myTree->Book("RunNumber",0,HZZ4lNtupleFactory::kInt);
  myTree->Book("EventNumber",0,HZZ4lNtupleFactory::kLong);
  myTree->Book("LumiNumber",0,HZZ4lNtupleFactory::kInt);
  myTree->Book("iBC",0,HZZ4lNtupleFactory::kInt);
  myTree->Book("Nvtx",0,HZZ4lNtupleFactory::kInt);
  myTree->Book("NObsInt",0,HZZ4lNtupleFactory::kInt);
  myTree->Book("NTrueInt",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("PUWeight12",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("PFMET",-99,HZZ4lNtupleFactory::kFloat);
  //myTree->Book("nJets",0,HZZ4lNtupleFactory::kInt);
  myTree->Book("nCleanedJets",0,HZZ4lNtupleFactory::kInt);
  myTree->Book("nCleanedJetsPt30",0,HZZ4lNtupleFactory::kInt);
  myTree->Book("nCleanedJetsPt30BTagged",0,HZZ4lNtupleFactory::kInt);
  myTree->Book("trigWord",0,HZZ4lNtupleFactory::kShort);
  
  //H variables
  myTree->Book("ZZMass",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("ZZMassErr",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("ZZMassErrCorr",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("ZZMassPreFSR",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("ZZsel",0,HZZ4lNtupleFactory::kShort);
  myTree->Book("ZZPt",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("ZZEta",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("ZZPhi",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("CRflag",0,HZZ4lNtupleFactory::kInt);


  //Z1 variables
  myTree->Book("Z1Mass",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Z1Pt",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Z1Flav",0,HZZ4lNtupleFactory::kShort);

  //Kin refitted info
  if (addKinRefit) {
    myTree->Book("ZZMassRefit",0,HZZ4lNtupleFactory::kFloat);
    myTree->Book("ZZChi2KinFit",0,HZZ4lNtupleFactory::kFloat);
    myTree->Book("Z1MassRefit",0,HZZ4lNtupleFactory::kFloat);
  }

  if (addVtxFit){
    myTree->Book("ZZMassCFit",0,HZZ4lNtupleFactory::kFloat);
    myTree->Book("ZZChi2CFit",0,HZZ4lNtupleFactory::kFloat);

  }

  //Z2 variables
  myTree->Book("Z2Mass",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Z2Pt",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Z2Flav",0,HZZ4lNtupleFactory::kShort);

  //Angular variables
  myTree->Book("costhetastar",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("helphi",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("helphiZ1",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("helphiZ2",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("helcosthetaZ1",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("helcosthetaZ2",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("phistarZ1",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("phistarZ2",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("xi",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("xistar",0,HZZ4lNtupleFactory::kFloat);

  //Lepton variables
  myTree->Book("Lep1Pt",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep1Eta",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep1Phi",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep1LepId",0,HZZ4lNtupleFactory::kInt);
  myTree->Book("Lep1SIP",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep1isID",0,HZZ4lNtupleFactory::kBool);
  myTree->Book("Lep1BDT",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep1missingHit",0,HZZ4lNtupleFactory::kChar);
  //--> Commented as these parentIDs are not computed at the moment
  //  myTree->Book("Lep1ParentId",0,HZZ4lNtupleFactory::kFloat);

  myTree->Book("Lep2Pt",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep2Eta",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep2Phi",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep2LepId",0,HZZ4lNtupleFactory::kInt);
  myTree->Book("Lep2SIP",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep2isID",0,HZZ4lNtupleFactory::kBool);
  myTree->Book("Lep2BDT",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep2missingHit",0,HZZ4lNtupleFactory::kChar);
  //  myTree->Book("Lep2ParentId",0,HZZ4lNtupleFactory::kFloat);

  myTree->Book("Lep3Pt",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep3Eta",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep3Phi",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep3LepId",0,HZZ4lNtupleFactory::kInt);
  myTree->Book("Lep3SIP",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep3isID",0,HZZ4lNtupleFactory::kBool);
  myTree->Book("Lep3BDT",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep3missingHit",0,HZZ4lNtupleFactory::kChar);
  //  myTree->Book("Lep3ParentId",0,HZZ4lNtupleFactory::kFloat);

  myTree->Book("Lep4Pt",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep4Eta",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep4Phi",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep4LepId",0,HZZ4lNtupleFactory::kInt);
  myTree->Book("Lep4SIP",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep4isID",0,HZZ4lNtupleFactory::kBool);
  myTree->Book("Lep4BDT",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep4missingHit",0,HZZ4lNtupleFactory::kChar);
  //  myTree->Book("Lep4ParentId",0,HZZ4lNtupleFactory::kFloat);

  //Lepton isolation variables
  myTree->Book("Lep1chargedHadIso",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep1neutralHadIso",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep1photonIso",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep1combRelIsoPF",0,HZZ4lNtupleFactory::kFloat);

  myTree->Book("Lep2chargedHadIso",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep2neutralHadIso",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep2photonIso",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep2combRelIsoPF",0,HZZ4lNtupleFactory::kFloat);

  myTree->Book("Lep3chargedHadIso",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep3neutralHadIso",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep3photonIso",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep3combRelIsoPF",0,HZZ4lNtupleFactory::kFloat);

  myTree->Book("Lep4chargedHadIso",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep4neutralHadIso",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep4photonIso",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Lep4combRelIsoPF",0,HZZ4lNtupleFactory::kFloat);

  //Photon variables
  if (writePhotons) {
    myTree->Book("PhotPt",0,HZZ4lNtupleFactory::kFloat);
    myTree->Book("PhotEta",0,HZZ4lNtupleFactory::kFloat);
    myTree->Book("PhotPhi",0,HZZ4lNtupleFactory::kFloat);
  }

  //Discriminants
  myTree->Book("p0plus_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p0minus_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p0plus_VAMCFM",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p0hplus_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p1_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p1_prodIndep_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p1plus_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p1plus_prodIndep_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2_prodIndep_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2qqb_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2hplus_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2hminus_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2bplus_VAJHU",0,HZZ4lNtupleFactory::kFloat);

  myTree->Book("p2hplus_qqb_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2hplus_prodIndep_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2hminus_qqb_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2hminus_prodIndep_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2bplus_qqb_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2bplus_prodIndep_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2h2plus_gg_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2h2plus_qqbar_VAJHU"               ,               0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2h2plus_prodIndep_VAJHU"   ,       0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2h3plus_gg_VAJHU"          ,       0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2h3plus_qqbar_VAJHU"       ,       0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2h3plus_prodIndep_VAJHU"   ,       0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2h6plus_gg_VAJHU"          ,       0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2h6plus_qqbar_VAJHU"       ,       0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2h6plus_prodIndep_VAJHU"   ,       0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2h7plus_gg_VAJHU"          ,       0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2h7plus_qqbar_VAJHU"       ,       0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2h7plus_prodIndep_VAJHU"   ,       0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2h9minus_gg_VAJHU"         ,               0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2h9minus_qqbar_VAJHU"      ,               0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2h9minus_prodIndep_VAJHU"  ,               0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2h10minus_gg_VAJHU"       ,                0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2h10minus_qqbar_VAJHU"    ,                0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p2h10minus_prodIndep_VAJHU",                0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("bkg_VAMCFM",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("bkg_prodIndep_VAMCFM",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("ggzz_VAMCFM",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("ggzz_p0plus_VAMCFM",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("ggzz_c1_VAMCFM",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("ggzz_c5_VAMCFM",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("ggzz_ci_VAMCFM",0,HZZ4lNtupleFactory::kFloat);
  
  myTree->Book("p0plus_m4l",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("bkg_m4l",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("pg1g4_mela",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("pg1g4_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("pg1g4_pi2_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("pg1g2_pi2_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("pg1g2_mela",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("pg1g2_VAJHU",0,HZZ4lNtupleFactory::kFloat);

  myTree->Book("p0plus_m4l_ScaleUp",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("bkg_m4l_ScaleUp",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p0plus_m4l_ScaleDown",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("bkg_m4l_ScaleDown",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p0plus_m4l_ResUp",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("bkg_m4l_ResUp",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p0plus_m4l_ResDown",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("bkg_m4l_ResDown",0,HZZ4lNtupleFactory::kFloat);

  //Production MELA
  myTree->Book("phjj_VAJHU_old",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("pvbf_VAJHU_old",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("phjj_VAJHU_old_up",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("pvbf_VAJHU_old_up",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("phjj_VAJHU_old_dn",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("pvbf_VAJHU_old_dn",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("phjj_VAJHU_new",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("pvbf_VAJHU_new",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("phjj_VAJHU_new_up",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("pvbf_VAJHU_new_up",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("phjj_VAJHU_new_dn",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("pvbf_VAJHU_new_dn",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p0_g1prime2_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("pg1g1prime2_VAJHU",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("Dgg10_VAMCFM",0,HZZ4lNtupleFactory::kFloat);

  myTree->Book("pzzzg_VAJHU",      0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("pzzgg_VAJHU",      0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("pzzzg_PS_VAJHU",   0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("pzzgg_PS_VAJHU",   0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p0Zgs_VAJHU",      0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p0gsgs_VAJHU",     0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p0Zgs_PS_VAJHU",   0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("p0gsgs_PS_VAJHU",  0,HZZ4lNtupleFactory::kFloat);
  
  //Jet variables
  myTree->Book("JetPt",0,HZZ4lNtupleFactory::kVectorFloat);
  myTree->Book("JetEta",0,HZZ4lNtupleFactory::kVectorFloat);
  myTree->Book("JetPhi",0,HZZ4lNtupleFactory::kVectorFloat);
  myTree->Book("JetMass",0,HZZ4lNtupleFactory::kVectorFloat);
  myTree->Book("JetBTagger",0,HZZ4lNtupleFactory::kVectorFloat);
  myTree->Book("JetIsBtagged",0,HZZ4lNtupleFactory::kVectorFloat);
  myTree->Book("JetQGLikelihood",0,HZZ4lNtupleFactory::kVectorFloat);
  myTree->Book("JetSigma",0,HZZ4lNtupleFactory::kVectorFloat);
  myTree->Book("DiJetMass",-99,HZZ4lNtupleFactory::kFloat);
  myTree->Book("DiJetMassPlus",-99,HZZ4lNtupleFactory::kFloat);
  myTree->Book("DiJetMassMinus",-99,HZZ4lNtupleFactory::kFloat);
  myTree->Book("DiJetDEta",-99,HZZ4lNtupleFactory::kFloat);
  myTree->Book("DiJetFisher",-99,HZZ4lNtupleFactory::kFloat);

  //Categorization-related variables
  myTree->Book("nExtraLep",0,HZZ4lNtupleFactory::kInt);
  myTree->Book("nExtraZ",0,HZZ4lNtupleFactory::kInt);

  //Variables of extra leptons
  myTree->Book("ExtraLepPt",0,HZZ4lNtupleFactory::kVectorFloat);
  myTree->Book("ExtraLepEta",0,HZZ4lNtupleFactory::kVectorFloat);
  myTree->Book("ExtraLepPhi",0,HZZ4lNtupleFactory::kVectorFloat);
  myTree->Book("ExtraLepLepId",0,HZZ4lNtupleFactory::kVectorFloat);


// Extended information on extra leptons.
// Currently not used, skip them for the time being
//   myTree->Book("ExtraLep1SIP",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("ExtraLep1isID",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("ExtraLep1BDT",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("ExtraLep1missingHit",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("ExtraLep1chargedHadIso",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("ExtraLep1neutralHadIso",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("ExtraLep1photonIso",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("ExtraLep1combRelIsoPF",0,HZZ4lNtupleFactory::kFloat);
//
//   myTree->Book("ExtraLep2SIP",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("ExtraLep2isID",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("ExtraLep2BDT",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("ExtraLep2missingHit",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("ExtraLep2chargedHadIso",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("ExtraLep2neutralHadIso",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("ExtraLep2photonIso",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("ExtraLep2combRelIsoPF",0,HZZ4lNtupleFactory::kFloat);
//
//   myTree->Book("ExtraLep3SIP",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("ExtraLep3isID",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("ExtraLep3BDT",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("ExtraLep3missingHit",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("ExtraLep3chargedHadIso",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("ExtraLep3neutralHadIso",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("ExtraLep3photonIso",0,HZZ4lNtupleFactory::kFloat);
//   myTree->Book("ExtraLep3combRelIsoPF",0,HZZ4lNtupleFactory::kFloat);

  // Gen variables
  // FIXME: don't book these fo data so we save some disk space...
  myTree->Book("genFinalState",0,HZZ4lNtupleFactory::kInt);
  myTree->Book("genProcessId",0,HZZ4lNtupleFactory::kInt);
  myTree->Book("genHEPMCweight",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("genExtInfo",0,HZZ4lNtupleFactory::kShort);
  myTree->Book("xsec",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("dataMCWeight",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("HqTMCweight",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("ZXFakeweight",0,HZZ4lNtupleFactory::kFloat);

  myTree->Book("GenHMass",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("GenHPt",0,HZZ4lNtupleFactory::kFloat);

  myTree->Book("GenZ1Mass",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("GenZ1Pt",0,HZZ4lNtupleFactory::kFloat);

  myTree->Book("GenZ2Mass",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("GenZ2Pt",0,HZZ4lNtupleFactory::kFloat);

  myTree->Book("GenLep1Pt",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("GenLep1Eta",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("GenLep1Phi",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("GenLep1Id",0,HZZ4lNtupleFactory::kShort);

  myTree->Book("GenLep2Pt",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("GenLep2Eta",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("GenLep2Phi",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("GenLep2Id",0,HZZ4lNtupleFactory::kShort);

  myTree->Book("GenLep3Pt",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("GenLep3Eta",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("GenLep3Phi",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("GenLep3Id",0,HZZ4lNtupleFactory::kShort);

  myTree->Book("GenLep4Pt",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("GenLep4Eta",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("GenLep4Phi",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("GenLep4Id",0,HZZ4lNtupleFactory::kShort);

  myTree->Book("GenAssocLep1Pt",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("GenAssocLep1Eta",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("GenAssocLep1Phi",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("GenAssocLep1Id",0,HZZ4lNtupleFactory::kShort);

  myTree->Book("GenAssocLep2Pt",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("GenAssocLep2Eta",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("GenAssocLep2Phi",0,HZZ4lNtupleFactory::kFloat);
  myTree->Book("GenAssocLep2Id",0,HZZ4lNtupleFactory::kShort);


}

//define this as a plug-in
DEFINE_FWK_MODULE(HZZ4lNtupleMaker);
