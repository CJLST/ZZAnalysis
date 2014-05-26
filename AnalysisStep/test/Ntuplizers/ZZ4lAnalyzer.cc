/** \class ZZ4lAnalyzer  
 *
 *  For the time being: retrieving all relevant information from candidates and filling the 2011 "flagship" plots
 *
 *  $Date: 2013/12/05 19:26:27 $
 *  $Revision: 1.111 $
 *   C. Botta - Torino 
 *   N. Amapane - Torino
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/LuminosityBlock.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <DataFormats/Common/interface/View.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/Math/interface/LorentzVector.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/JetReco/interface/PFJet.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <DataFormats/JetReco/interface/PFJetCollection.h>
#include <AnalysisDataFormats/CMGTools/interface/BaseMET.h>
#include <AnalysisDataFormats/CMGTools/interface/PFJet.h>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include <FWCore/ParameterSet/interface/FileInPath.h>


#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <DataFormats/Common/interface/MergeableCounter.h>

#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>
#include <ZZAnalysis/AnalysisStep/interface/MCHistoryTools.h>
#include <ZZAnalysis/AnalysisStep/interface/PUReweight.h>
#include <ZZAnalysis/AnalysisStep/interface/VBFCandidateJetSelector.h>
#include <ZZAnalysis/AnalysisStep/interface/Fisher.h>
#include "ZZ4lConfigHelper.h"
#include <boost/lexical_cast.hpp>

#include "knownCandidates.h"

#include <iostream>
#include <iterator>
#include <string>

#include <Math/VectorUtil.h>

#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"

#define HCMSSW
#include "ZZAnalysis/AnalysisStep/interface/Histograms.h"
#include <algorithm>

using namespace std;
using namespace edm;
using namespace reco;




float sqr(float x) {return x*x;}

class ZZ4lAnalyzer: public edm::EDAnalyzer {
public:

  ZZ4lAnalyzer(const edm::ParameterSet& pset);
  
  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
  void endLuminosityBlock(const edm::LuminosityBlock &, const edm::EventSetup&);
  virtual void beginJob();
  virtual void endJob();  
  
  void printParticle(const reco::Candidate* c=0, string idx="", int pdgId=0, 
		     float iso=-1, float SIP=-1);
  
private:
  ZZ4lConfigHelper myHelper;
  std::string theCandLabel;
  Channel theChannel;
  bool isMC;
  bool dumpMC;
  PUReweight reweight;

  TH2F* corrSigmaMu;
  TH2F* corrSigmaEle;

  // Selection counters
  double Nevt_Gen;
  double gen_ZZ4e;
  double gen_ZZ4mu;
  double gen_ZZ2mu2e;
  double gen_ZZtau;
  double gen_BUGGY;      //Events to be filtered out.
  double Nevt_PAT;
  double Nevt_Skim;
  double Nevt_TriggerBit;
  double Nevt_4l;         // Step 3 (4l, AFAS)
  double Nevt_Z1;
  double Nevt_Z1Mass;
  double Nevt_Z1Plus1Lep;
  double Nevt_ZZCand;
  double Nevt_ZZBestCand;
  double Nevt_ZZCandpT;
  double Nevt_ZZCandMZ1;
  double Nevt_ZZCandMZ2;
  double Nevt_ZZCandMAllComb;
  double Nevt_ZZCandM70;
  double Nevt_ZZCandM100;
  double Nevt_ZZCandM100_1Jet;
  double Nevt_ZZCandM100_2Jets;
  double Nevt_ZZCandM100_VBF;
  double Nevt_ZZCand_zz;
  double Nevt_MELA;
  double Nevt_psMELA;
  double Nevt_grMELA;
  double Nevt_WithFSR;
  double Nevt_With1FSR;
  double Nevt_With2FSR;
  double Nevt_WithFSRImprove;
  double weight;

  //histograms
  TH1F* nEventComplete;
  TH1F* hZ1Mass;
  TH1F* hZ1MassEBEB;
  TH1F* hZ1MassEEEE;
  TH1F* hZ1Mass_w;
  TH1F* hZ1MassEBEB_w;
  TH1F* hZ1MassEEEE_w;

//   HCand* hCandZZBestCand;
//   HCand* hCandZZCandMZ2;
//   HCand* hCandZZCandpT;
//   HCand* hCandZZCandMAllComb;
//   HCand* hCandZZCandM70;
//   HCand* hCandZZCandM100;
//   HCand* hCandZZCand;
//   HCand* hCandZZCandM70_w;
//   HCand* hCandZZCandM100_w;
//   HCand* hCandZZCand_w;

  TH1F* h_dR_all;
  TH1F* h_dR_pure;
  TH1F* h_ZZMass_fsr;
  TH1F* h_ZZMassPreFsr_fsr;
  TH1F* hNvtxNoWeight;
  TH1F* hNvtxWeight;

  string sampleName;
};

// Constructor
ZZ4lAnalyzer::ZZ4lAnalyzer(const ParameterSet& pset) :
  myHelper(pset),
  theCandLabel(pset.getUntrackedParameter<string>("candCollection")),
  dumpMC(pset.getUntrackedParameter<bool>("dumpMC",false)),
  reweight(PUReweight::LEGACY),
  Nevt_Gen(0),
  gen_ZZ4e(0),
  gen_ZZ4mu(0),
  gen_ZZ2mu2e(0),
  gen_ZZtau(0),
  gen_BUGGY(0),
  Nevt_PAT(0),
  Nevt_Skim(0),
  Nevt_TriggerBit(0),
  Nevt_4l(0),
  Nevt_Z1(0),
  Nevt_Z1Mass(0),
  Nevt_Z1Plus1Lep(0),
  Nevt_ZZCand(0),
  Nevt_ZZBestCand(0),
  Nevt_ZZCandpT(0),
  Nevt_ZZCandMZ1(0), 
  Nevt_ZZCandMZ2(0), 
  Nevt_ZZCandMAllComb(0), 
  Nevt_ZZCandM70(0),
  Nevt_ZZCandM100(0),
  Nevt_ZZCandM100_1Jet(0),
  Nevt_ZZCandM100_2Jets(0),
  Nevt_ZZCandM100_VBF(0),
  Nevt_ZZCand_zz(0), 
  Nevt_MELA(0), 
  Nevt_psMELA(0), 
  Nevt_grMELA(0), 
  Nevt_WithFSR(0),
  Nevt_With1FSR(0),
  Nevt_With2FSR(0),
  Nevt_WithFSRImprove(0),
  weight(1.),
  sampleName(pset.getParameter<string>("sampleName"))
{

  isMC = myHelper.isMC();

  theChannel = myHelper.channel();
  if (theChannel>3) {
    cout << "ERROR: ZZ4lAnalyzer: channel "<< theChannel << " is not valid" <<endl;
    abort();
  }
}


void ZZ4lAnalyzer::beginJob(){
  // Book histograms
  edm::Service<TFileService> fileService;

  // Counting Histograms. We do not want sumw2 on these!
  TH1F::SetDefaultSumw2(kFALSE);
  nEventComplete = fileService->make<TH1F>("nEventComplete", "nEventComplete", 12, 0., 12.);
  hZ1Mass = fileService->make<TH1F>("hZ1Mass", "hZ1Mass", 400, 0., 200.);
  hZ1Mass_w = fileService->make<TH1F>("hZ1Mass_w", "hZ1Mass_w", 400, 0., 200.);
  
  if (theChannel==EEEE) {
    hZ1MassEBEB = fileService->make<TH1F>("hZ1MassEBEB", "hZ1MassEBEB", 400, 0., 200.);
    hZ1MassEEEE = fileService->make<TH1F>("hZ1MassEEEE", "hZ1MassEEEE", 400, 0., 200.);
    hZ1MassEBEB_w = fileService->make<TH1F>("hZ1MassEBEB_w", "hZ1MassEBEB_w", 400, 0., 200.);
    hZ1MassEEEE_w = fileService->make<TH1F>("hZ1MassEEEE_w", "hZ1MassEEEE_w", 400, 0., 200.);
  }

  h_dR_all = fileService->make<TH1F>("h_dR_all","dR(#gamma radiating lepton) for all selected fsr photons",50.,0.,0.5);
  h_dR_pure = fileService->make<TH1F>("h_dR_pure","dR(#gamma radiating lepton) for pure selected fsr photons",50.,0.,0.5);
  h_ZZMass_fsr = fileService->make<TH1F>("h_ZZMass_fsr","ZZ mass after fsr recovery for events affected by fsr algorithm",100,80.,160.);
  h_ZZMassPreFsr_fsr = fileService->make<TH1F>("h_ZZMassPreFsr_fsr","ZZ mass before fsr recovery for events affected by fsr algorithm",100,80.,160.);

  hNvtxNoWeight = fileService->make<TH1F>("hNvtxNoWeight","hNvtxNoWeight",100,0,100);
  hNvtxWeight = fileService->make<TH1F>("hNvtxWeight","hNvtxWeight",100,0,100);


  // Distributions
  TH1F::SetDefaultSumw2(kTRUE);
  // hCand Histogram for step plots >= 6
//   hCandZZBestCand = new HCand("hCandZZBestCand ");
//   hCandZZCandMZ2 = new HCand("hCandZZCandMZ2");
//   hCandZZCandpT = new HCand("hCandZZCandpT");
//   hCandZZCandMAllComb = new HCand("hCandZZCandMAllComb");
//   hCandZZCandM70 = new HCand("hCandZZCandM70");
//   hCandZZCandM100 = new HCand("hCandZZCandM100");
//   hCandZZCand = new HCand("hCandZZCand");

//   hCandZZCandM70_w = new HCand("hCandZZCandM70_w");
//   hCandZZCandM100_w = new HCand("hCandZZCandM100_w");
//   hCandZZCand_w = new HCand("hCandZZCand_w");
   
}


void ZZ4lAnalyzer::endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup) 
{

  double Nevt_preskim = -1.;
  edm::Handle<edm::MergeableCounter> preSkimCounter;
  if (iLumi.getByLabel("preSkimCounter", preSkimCounter)) { // Counter before skim. Does not exist for non-skimmed samples.
    Nevt_preskim = preSkimCounter->value;
  }  
  
  edm::Handle<edm::MergeableCounter> prePathCounter;
  iLumi.getByLabel("prePathCounter", prePathCounter);       // Counter of input events in the input pattuple

  // Nevt_gen: this is the number before any skim
  if (Nevt_preskim>=0.) {
    Nevt_Gen = Nevt_Gen + Nevt_preskim; 
  } else {
    Nevt_Gen = Nevt_Gen + prePathCounter->value;    
  }

  // Nevt_PAT: number of events in the pattuple
  Nevt_PAT = Nevt_PAT + prePathCounter->value;


  // cout << "Nevt_Gen_lumi: " << Nevt_Gen << endl;
  // cout << "Nevt_afterSkim_lumi: " << Nevt_afterSkim << endl;
  
}



void ZZ4lAnalyzer::endJob(){

  string sFS=finalStateNiceName(theChannel);
  sFS = sFS+"_";
  cout << endl;
  cout << "*************************" <<endl;
  cout << "Final state: " << sFS << endl;
  cout << "*************************" <<endl;
  cout << sFS+"Nevt_Gen:         " << Nevt_Gen << endl;
  if (isMC) {
    cout << sFS+"gen_BUGGY:        " << gen_BUGGY << endl;
    cout << sFS+"gen_ZZ4mu:        " << gen_ZZ4mu << endl;
    cout << sFS+ "gen_ZZ4e:         " << gen_ZZ4e << endl;
    cout << sFS+ "gen_ZZ2mu2e:      " << gen_ZZ2mu2e << endl;
    cout << sFS+ "gen_ZZtau:        " << gen_ZZtau << endl;
  }
  cout << sFS+ "Nevt_PAT:         " << Nevt_PAT << endl;
  cout << sFS+ "Nevt_Skim:        " << Nevt_Skim << endl;
  cout << sFS+ "Nevt_TriggerBit:  " << Nevt_TriggerBit << endl;


  { // 2012 selection
    cout << sFS+ "Nevt_Z1:               " << Nevt_Z1 << endl;  
    cout << sFS+ "Nevt_Z1Mass:           " << Nevt_Z1Mass << endl;  
    cout << sFS+ "Nevt_Z1Plus1Lep:       " << Nevt_Z1Plus1Lep << endl;  
    cout << sFS+ "Nevt_ZZBestCand:       " << Nevt_ZZBestCand << endl;
    cout << sFS+ "Nevt_ZZCandMZ2:        " << Nevt_ZZCandMZ2 << endl;
    cout << sFS+ "Nevt_ZZCandpT:         " << Nevt_ZZCandpT << endl;
    cout << sFS+ "Nevt_ZZCandMAllComb:   " << Nevt_ZZCandMAllComb << endl;
    cout << sFS+ "Nevt_ZZCandM70:        " << Nevt_ZZCandM70 << endl;
	cout << sFS+ "Nevt_ZZCandM100:       " << Nevt_ZZCandM100 << endl;
	cout << sFS+ "Nevt_ZZCandM100_1Jet:  " << Nevt_ZZCandM100_1Jet << endl;
	cout << sFS+ "Nevt_ZZCandM100_2Jets: " << Nevt_ZZCandM100_2Jets << endl;
	cout << sFS+ "Nevt_ZZCandM100_VBF:   " << Nevt_ZZCandM100_VBF << endl;
    cout << sFS+ "Nevt_MELA:             " << Nevt_MELA    << endl;
    cout << sFS+ "Nevt_pseudoMELA:       " << Nevt_psMELA    << endl;
    cout << sFS+ "Nevt_graviMELA:        " << Nevt_grMELA    << endl;
    cout << sFS+ "Nevt_WithFSR:          " << Nevt_WithFSR    << endl;
    cout << sFS+ "Nevt_With1FSR:         " << Nevt_With1FSR   << endl;
    cout << sFS+ "Nevt_With2FSR:         " << Nevt_With2FSR    << endl;
    cout << sFS+ "Nevt_WithFSRImprove:   " << Nevt_WithFSRImprove    << endl;



    //----------------------------------------------------------------------
    // nEvent Histos for Step Plot: "Nominal Signal Selection"
    //----------------------------------------------------------------------
    
    nEventComplete->GetXaxis()->SetBinLabel(1,"Skim");             // HZZ Skim
    nEventComplete->GetXaxis()->SetBinLabel(2,"HLT");              // TriggerBit 
    nEventComplete->GetXaxis()->SetBinLabel(3,"Z1");               // At least one Z1
    nEventComplete->GetXaxis()->SetBinLabel(4,"m_{Z1}");           // At least one Z1Mass
    nEventComplete->GetXaxis()->SetBinLabel(5,"Z1+l");             // At least one Z1+l (good lepton)
    nEventComplete->GetXaxis()->SetBinLabel(6,"ZZ");               // At least one LLLL candidate with Z1  
    nEventComplete->GetXaxis()->SetBinLabel(7,"m_{Z2}");           // At least one LLLL candidate with 4 < Z2M < 12
    nEventComplete->GetXaxis()->SetBinLabel(8,"pT>20,10");         // At least one LLLL candidate with pT 20, 10
    nEventComplete->GetXaxis()->SetBinLabel(9,"m_{ll}>4");         // At least one LLLL candidate with 6/6 mll > 4 
    nEventComplete->GetXaxis()->SetBinLabel(10,"m_{4l}>70");       // At least one LLLL candidate with mlll > 70
    nEventComplete->GetXaxis()->SetBinLabel(11,"m_{4l}>100");      // At least one LLLL candidate with mlll > 100
    nEventComplete->GetXaxis()->SetBinLabel(12,"ZZ sel.");              // High-mass (ZZ selection) 

    nEventComplete->SetBinContent(1,Nevt_Skim); 
    nEventComplete->SetBinContent(2,Nevt_TriggerBit); 
    nEventComplete->SetBinContent(3,Nevt_Z1);
    nEventComplete->SetBinContent(4,Nevt_Z1Mass);
    nEventComplete->SetBinContent(5,Nevt_Z1Plus1Lep);
    nEventComplete->SetBinContent(6,Nevt_ZZBestCand);
    nEventComplete->SetBinContent(7,Nevt_ZZCandMZ2);
    nEventComplete->SetBinContent(8,Nevt_ZZCandpT);
    nEventComplete->SetBinContent(9,Nevt_ZZCandMAllComb);
    nEventComplete->SetBinContent(10,Nevt_ZZCandM70);
    nEventComplete->SetBinContent(11,Nevt_ZZCandM100);
    nEventComplete->SetBinContent(12,Nevt_ZZCand_zz);
  
  }
  
  
  if (isMC){
    nEventComplete->SetBinContent(0,Nevt_Gen-gen_BUGGY); // Save normalization factor in underflow bin (correcting for buggy MC events)
  } else {
    nEventComplete->SetBinContent(0,1.);
  }

}



void ZZ4lAnalyzer::analyze(const Event & event, const EventSetup& eventSetup){  
  int irun=event.id().run();
  long long int ievt=event.id().event(); 
  int ils =event.luminosityBlock();

  //  int nObsInt  = -1;
  float nTrueInt = -1.;
  float PUweight = 1.;

  Handle<vector<Vertex> >  vertexs;
  event.getByLabel("goodPrimaryVertices",vertexs);
  int Nvtx = vertexs->size();

  const reco::Candidate * genH = 0;
//   std::vector<const reco::Candidate *> genZs;
//   std::vector<const reco::Candidate *> genZLeps;

  int genFinalState = NONE;
  if (isMC) {
    MCHistoryTools mch(event,sampleName);
    genFinalState = mch.genFinalState();

    if (genFinalState == EEEE) {
      ++gen_ZZ4e;
//       if (gen_ZZ4lInEtaAcceptance) ++gen_ZZ4e_EtaAcceptance;
//       if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ4e_LeptonAcceptance;
//       if (gen_ZZInAcceptance) {
// 	++gen_ZZ4e_MassAcceptance;
// 	if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ4e_MassPtAcceptance;
//      }
    } else if (genFinalState == MMMM) {
      ++gen_ZZ4mu;
//       if (gen_ZZ4lInEtaAcceptance) ++gen_ZZ4mu_EtaAcceptance;
//       if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ4mu_LeptonAcceptance;
//       if (gen_ZZInAcceptance) {
// 	++gen_ZZ4mu_MassAcceptance;
// 	if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ4mu_MassPtAcceptance;
//       }
    } else if (genFinalState == EEMM) {
      ++gen_ZZ2mu2e;
//       if (gen_ZZ4lInEtaAcceptance) ++gen_ZZ2mu2e_EtaAcceptance;
//       if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ2mu2e_LeptonAcceptance;
//       if (gen_ZZInAcceptance) {
// 	++gen_ZZ2mu2e_MassAcceptance;
// 	if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ2mu2e_MassPtAcceptance;
//       }
    } else if (genFinalState == LLTT){
      ++gen_ZZtau;
    } else if (genFinalState == BUGGY){ // handle H->ddbar gen bug!!!
      ++gen_BUGGY;
      return;


      genH = mch.genH();
//       genZs = mch.genZs();
//       genZLeps = mch.genZLeps();

    }


    Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    event.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);

    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      if(PVI->getBunchCrossing() == 0) { 
	//	nObsInt  = PVI->getPU_NumInteractions();
	nTrueInt = PVI->getTrueNumInteractions();
	break;
      } 
    }

    int source = myHelper.sampleType();
    int target = myHelper.setup();
    PUweight = reweight.weight(source,target,nTrueInt);
    hNvtxNoWeight->Fill(Nvtx);    
    hNvtxWeight->Fill(Nvtx,PUweight);

  }

  // theChannel Candidates
  Handle<View<pat::CompositeCandidate> > candHandle;
  event.getByLabel(theCandLabel, candHandle);
  const View<pat::CompositeCandidate>* cands = candHandle.product();

  // Muons 
  Handle<pat::MuonCollection> muons;
  event.getByLabel("softMuons", muons);

  // Electrons
  Handle<pat::ElectronCollection> electrons;
  event.getByLabel("cleanSoftElectrons", electrons);

  // Zs (for step plots)
  Handle<View<pat::CompositeCandidate> > Z;
  event.getByLabel("ZCand", Z);
 

  // Z+l (for step plots)
  Handle<View<reco::CompositeCandidate> > Zl;
  event.getByLabel("ZlCand", Zl);
  const View<reco::CompositeCandidate>* trileps = Zl.product();
  
  Handle<double> rhoHandle;
  event.getByLabel(InputTag("kt6PFJetsForIso","rho"), rhoHandle);

  float pfmet = -1;  
//   Handle<reco::PFMETCollection> pfmetcoll;
//   event.getByLabel("patMETs", pfmetcoll);      // This the type 1 corrected MET
//   //  event.getByLabel("patMETsRaw", pfmetcoll);   // This is the raw PFMET
//   float pfmet = -1;
//   if(pfmetcoll.isValid()) { 
//     const PFMETCollection *pfmetcol = pfmetcoll.product();
//     const PFMET *pfmetObj = &(pfmetcol->front());
//     pfmet = pfmetObj->pt();
//   }
  
  //MET
  Handle<vector<cmg::BaseMET> > pfmetcoll;
  event.getByLabel("cmgPFMET", pfmetcoll);
  if(pfmetcoll.isValid()){
    pfmet = pfmetcoll->front().pt();
  }

  // Jet collection (preselected with pT>10)
  Handle<edm::View<cmg::PFJet> > pfjetscoll;
  event.getByLabel("cmgPFJetSel", pfjetscoll);

  // Apply MC filter (skip event)
  if (isMC && !(myHelper.passMCFilter(event))) return;

  // Skim
  bool evtPassSkim = myHelper.passSkim(event);
  
  // Trigger requests
  bool evtPassTrigger = myHelper.passTrigger(event);
  

  // "Z1" and "Z1+mass" steps
  bool evtPassZ1 = false;
  bool evtPassZ1Mass = false;
  float ZCandMass = 0;

  // print a status 
  if (!isMC&&theChannel==EEMM) cout << irun << ":" << ils << ":" << ievt << " " << evtPassSkim << " " << evtPassTrigger << endl;

  if(!evtPassSkim) return;
  Nevt_Skim+= weight;
  if(!evtPassTrigger) return;
  Nevt_TriggerBit+= weight;
  
  for( View<pat::CompositeCandidate>::const_iterator cand = Z->begin(); cand != Z->end(); ++cand) {
    
    
    //     cout << "Z: " << cand->mass()
    // 	 << " " << cand->userFloat("GoodLeptons")
    // 	 << " " << cand->userFloat("isBestZee")
    // 	 << " " << cand->userFloat("isBestZmm")
    // 	 << " " << cand->userFloat("isBestZ") <<endl;
    if (cand->userFloat("GoodLeptons")) {
      ZCandMass = cand->mass();      
      if ((theChannel==EEEE && cand->userFloat("isBestZee")) ||
	  (theChannel==MMMM && cand->userFloat("isBestZmm"))  ||
	  (((theChannel==EEMM||theChannel==ZZ)) && cand->userFloat("isBestZ"))) {
	evtPassZ1 = true;
	if (ZCandMass > 40 && ZCandMass <120) { 
	  evtPassZ1Mass = true;
	  if ((theChannel==EEEE && cand->userFloat("isBestZee")) ||
	      (theChannel==MMMM && cand->userFloat("isBestZmm")) ||
	      (theChannel==EEMM && cand->userFloat("isBestZ"))) {
	    
	    hZ1Mass->Fill(ZCandMass, 1.);
	    hZ1Mass_w->Fill(ZCandMass, PUweight);
	    if (theChannel==EEEE) {
	      const pat::Electron* e1 = dynamic_cast<const pat::Electron*>(cand->daughter(0)->masterClone().get());
	      const pat::Electron* e2 = dynamic_cast<const pat::Electron*>(cand->daughter(1)->masterClone().get());

	      if (e1->isEB() && e2->isEB()){
		hZ1MassEBEB->Fill(ZCandMass, 1.);
		hZ1MassEBEB_w->Fill(ZCandMass, PUweight);
	      }
	      if (e1->isEE() && e2->isEE()){
		hZ1MassEEEE->Fill(ZCandMass, 1.);
		hZ1MassEEEE_w->Fill(ZCandMass, PUweight);
	      }	      
	    }
	  }
	}
      }     
    }  
  }
  
  // Trilepton step. We just ask a third good lepton on top of evtPassZ1
  bool evtPassTriGoodLep = false;
  for(View<reco::CompositeCandidate>::const_iterator trilep = trileps->begin(); trilep != trileps->end(); ++trilep) {
    if (userdatahelpers::getUserFloat(trilep->daughter(0),"GoodLeptons") && 
	userdatahelpers::getUserFloat(trilep->daughter(1),"isGood") && userdatahelpers::getUserFloat(trilep->daughter(1),"combRelIsoPF")<0.4) { //FIXME this is ISO without FSR!!!!
      evtPassTriGoodLep=true;
    }
  }  

  //----------------------------------------------------------------------
  // Loop on ZZ candidates
  //----------------------------------------------------------------------

  // Event Counters for step plot (PRL selection)

  // Event Counters for step plot (new selection)
  bool evtPassZZ   = false;
  bool evtPassIsBestZ1 = false;
  bool evtPass20_10 = false;
  //  bool evtPassZ1Mass_B = false;
  bool evtPassZ2Mass = false;
  bool evtPassMLLallComb = false;
  bool evtPassM70 = false;
  bool evtPassM100 = false;
  bool evtPass1Jet = false;
  bool evtPass2Jets = false;
  bool evtPassVBF = false;
  bool evtPassM100_zz = false;
  bool evtPassMELA = false;
  bool evtPasspsMELA = false;
  bool evtPassgrMELA = false;
  for( View<pat::CompositeCandidate>::const_iterator cand = cands->begin(); cand != cands->end(); ++cand) {
    //--- Extract all relevant candidate information and increment counters
      
    // Pointers to Z and leptons
    const Candidate* Z1   = cand->daughter("Z1");
    const Candidate* Z2   = cand->daughter("Z2");
    vector<const Candidate*> leptons(4);
    vector<string> labels(4);
    userdatahelpers::getSortedLeptons(*cand, leptons, labels);
    const Candidate* Z1Lp = leptons[0];
    const Candidate* Z1Ln = leptons[1];
    const Candidate* Z2Lp = leptons[2];
    const Candidate* Z2Ln = leptons[3];

    int Z1dauFSR = userdatahelpers::getUserFloat(Z1,"dauWithFSR");
    int Z2dauFSR = userdatahelpers::getUserFloat(Z2,"dauWithFSR");
    const reco::Candidate* fsrPhotonZ1 = 0;
    const reco::Candidate* fsrPhotonZ2 = 0;
    if (Z1dauFSR >= 0) fsrPhotonZ1 = cand->daughter("Z1")->daughter(2);
    if (Z2dauFSR >= 0) fsrPhotonZ2 = cand->daughter("Z2")->daughter(2);
    
    // Retrieve the userFloat of the leptons in vectors ordered in the same way.
    vector<float> SIP(4);
//     vector<float> looseIso(4);
//     vector<float> iso(4);
    vector<float> isoPF(4);
    vector<float> pt(4);
    for (int i=0; i<4; ++i){
      SIP[i]      = userdatahelpers::getUserFloat(leptons[i],"SIP");
//       looseIso[i] = userdatahelpers::getUserFloat(leptons[i],"looseIso");
//       iso[i]      = userdatahelpers::getUserFloat(leptons[i],"combRelIso");
      isoPF[i]    = cand->userFloat(labels[i]+"combRelIsoPFFSRCorr"); // Note: the FSR-corrected iso is attached to the Z, not to the lepton!
      pt[i]       = leptons[i]->pt();

      // Check that I don't mess up with labels[] and leptons[]
      assert(SIP[i] == cand->userFloat(labels[i]+"SIP"));
    }
      

    // Sort lepton variables
    vector< float> ptS(pt);
//     vector<float> isoS(iso);
    vector<float> isoPFS(isoPF);
    vector<float> SIPS(SIP);
    sort(ptS.begin(),ptS.end());
//    sort(isoS.begin(),isoS.end());
    sort(isoPFS.begin(),isoPFS.end());
    sort(SIPS.begin(),SIPS.end());



    // Masses
    float ZZMass = cand->p4().mass();
    float ZZPt = cand->pt();
    float ZZMassErr = cand->userFloat("massError");
    float ZZMassErrCorr = cand->userFloat("massErrorCorr");
    float Z1Mass = Z1->mass();
    float Z2Mass = Z2->mass();	
    float Z1OtherCombMass = cand->userFloat("mZa");
    float Z2OtherCombMass = cand->userFloat("mZb");
      
    // Cut variables
    //    float relIso_sum2least = cand->userFloat("iso34");
    //    float SIP4      = cand->userFloat("SIP4");
      
    // Decay Angles and MELA discriminant
    float costheta1 = cand->userFloat("costheta1");
    float costheta2 = cand->userFloat("costheta2");
    float Phi       = cand->userFloat("phi");
    float costhetastar = cand->userFloat("costhetastar");
    float Phi1      = cand->userFloat("phistar1");
	//    bool  swapZ1Z2  = false;
    //    float LD        = cand->userFloat("LD");
    //    float pseudoLD  = cand->userFloat("pseudoLD");
    //    float graviLD   = cand->userFloat("spin2PMLD");
    //    float psig        = cand->userFloat("PSig");
    //    float pbkg        = cand->userFloat("PBkg");
    
	 //here all the KDs
     //probabilities from the tree
    float p0plus_VAJHU = cand->userFloat("p0plus_VAJHU");
    float bkg_VAMCFM = cand->userFloat("bkg_VAMCFM");
    float p0minus_VAJHU = cand->userFloat("p0minus_VAJHU");
    float p0hplus_VAJHU = cand->userFloat("p0hplus_VAJHU");
    float p1plus_VAJHU = cand->userFloat("p1plus_VAJHU");
    float p1_VAJHU = cand->userFloat("p1_VAJHU");
    float p2_VAJHU = cand->userFloat("p2_VAJHU");
    float p2qqb_VAJHU = cand->userFloat("p2qqb_VAJHU");
    float p2hplus_VAJHU = cand->userFloat("p2hplus_VAJHU");
    float p2hminus_VAJHU = cand->userFloat("p2hminus_VAJHU");      
    float p2bplus_VAJHU = cand->userFloat("p2bplus_VAJHU");       
	  
    //KDs computation
    float KD = p0plus_VAJHU/(p0plus_VAJHU + bkg_VAMCFM) ;
    float KD_pseudo = p0plus_VAJHU/(p0plus_VAJHU + p0minus_VAJHU) ;
    float KD_highdim = p0plus_VAJHU/(p0plus_VAJHU + p0hplus_VAJHU) ;
    float KD_vec = p0plus_VAJHU/(p0plus_VAJHU + p1plus_VAJHU) ;
    float KD_psvec = p0plus_VAJHU/(p0plus_VAJHU + p1_VAJHU) ;
    float KD_gggrav = p0plus_VAJHU/(p0plus_VAJHU + p2_VAJHU) ;
    float KD_qqgrav = p0plus_VAJHU/(p0plus_VAJHU + p2qqb_VAJHU) ;
    float KD_2hp = p0plus_VAJHU/(p0plus_VAJHU + p2hplus_VAJHU) ;
    float KD_2hm = p0plus_VAJHU/(p0plus_VAJHU + p2hminus_VAJHU) ;
    float KD_2bp = p0plus_VAJHU/(p0plus_VAJHU + p2bplus_VAJHU) ;
	  
	  
    // Final selection as specified in the .py file. The selected candidate is
    // the one with (pass && isBest).
    bool candPass = cand->userFloat("FullSel");
    bool candPass70 = cand->userFloat("FullSel70");
    bool candIsBest = cand->userFloat("isBestCand");

    float VD   = cand->userFloat("VD");
    float deta = cand->userFloat("detajj");
    float mjj  = cand->userFloat("mjj");
    float j1pT =-1.;
    float j2pT =-1.;
    int nJet30 = -1;
    VBFCandidateJetSelector myVBFCandidateJetSelector;

    // Pick jets, since they are not associated to the candidate yet
    std::vector<const cmg::PFJet*> cleanedJets;
    std::vector<const cmg::PFJet*> cleanedJetsPt30;
    if ( candIsBest && candPass70 ) {
      cleanedJets = myVBFCandidateJetSelector.cleanJets(*cand,pfjetscoll,myHelper.setup()); 
    
      // Create a collection implementing the official jet selection (pt>30)
      for (unsigned int i=0; i < cleanedJets.size(); ++i){
	const cmg::PFJet& myjet = *(cleanedJets.at(i));  
	if (myjet.pt()>30) cleanedJetsPt30.push_back(&myjet);
      }

      // Also x-check VBF variables. FIXME: remove this once it's checked
      float VD_chk   =-99.;
      float deta_chk =-99.;
      float mjj_chk  =-99.;
      nJet30 =cleanedJetsPt30.size();
      if (nJet30>0) j1pT=cleanedJetsPt30[0]->pt();
      if (nJet30>1) {
	j2pT=cleanedJetsPt30[1]->pt();      
	deta_chk = cleanedJetsPt30[0]->eta()-cleanedJetsPt30[1]->eta();
	mjj_chk = (cleanedJetsPt30[0]->p4()+cleanedJetsPt30[1]->p4()).M();
	VD_chk = fisher(mjj,deta);
      }
      if (VD_chk!=VD || deta_chk != deta || mjj_chk != mjj) {
	cout << "ERROR: ZZ4lAnalyzer: " << VD_chk << " " << VD << " " <<  deta_chk << " " <<  deta << " " << mjj_chk << " " << mjj << endl;
      }
    }
	  
    //----------------------------------------------------------------------
    // count events at different steps (step plots)
    //----------------------------------------------------------------------
    int sel = 0; // Selection level, per candidate

//     bool passMzLoose = (Z1Mass>50. && Z1Mass<120. && Z2Mass>12. && Z2Mass<120.); //FIXME hardcoded cut
    bool passMz_zz = (Z1Mass>60. && Z1Mass<120. && Z2Mass>60. && Z2Mass<120.);   //FIXME hardcoded cut      
    { // New 2012 Nominal Signal Selection; split in all syncronization steps
      
      if (cand->userFloat("GoodLeptons")) { // One ZZ (not necessarily including the best Z)
	evtPassZZ = true; 
	vector<float> ptS(pt);
	sort(ptS.begin(),ptS.end());
	sel = 10;
	
	if (userdatahelpers::getUserFloat(Z1,"isBestZ")) {
	  evtPassIsBestZ1 = true;
	  sel = 20;
	  
	  if (Z1Mass > 40 && Z1Mass <120) { // 4c?
	    sel = 30;
	    
	    if (candIsBest) { // step 4c:  Pair #2 built (SF/OS, highest pT leptons) 

	      if (cand->userFloat("Z2Mass")) { // Step 4d. Apply mll cuts (4<mZ2<120)
		evtPassZ2Mass = true;
		sel = 40;
	      
		if (ptS[3] > 20 && ptS[2] > 10) { // Step  5: pt > 20/10 in ZZ
		  evtPass20_10 = true;
		  sel = 50;

		  if (cand->userFloat("MAllComb")) { // Step 6 : QCD suppression (mll>4 GeV cut on OS lepton pairs) 
		    evtPassMLLallComb = true;
		    sel = 60;
		  
		    if (ZZMass>70){  // Step 7: m(4l) > 70 (Z->4l phasespace)
		      sel=70;

		      if (Z2Mass>12){
			evtPassM70 = true;
			sel = 90;

			if (ZZMass>100){  // Step 8: m(4l) > 100
			  evtPassM100 = true;
			  sel =100; 
		      
			  //here include the VBF stuff
			  if (cleanedJetsPt30.size() > 0)
			  {
				 evtPass1Jet = true;
				  if (cleanedJetsPt30.size() == 2)
				  {
					  evtPass2Jets = true;
// 					  cmg::PFJet myjet1 = *(cleanedJetsPt30.at(0));
// 					  cmg::PFJet myjet2 = *(cleanedJetsPt30.at(1));
// 					  math::XYZTLorentzVector jet1 = myjet1.p4();
// 					  math::XYZTLorentzVector jet2 = myjet2.p4();
//					  double VD  = 0.18*fabs(jet1.Eta()-jet2.Eta()) + 1.92E-4*(jet1+jet2).M();
					  if(VD>0.4)
					    evtPassVBF = true;
				  }  
			  }
				
			  if (passMz_zz) {
			    evtPassM100_zz = true;
			  }
		      
			  if(KD > 0.1){
			    evtPassMELA = true;
			  }
			
			  if (KD_pseudo > 0.3) {
			    evtPasspsMELA = true;
			  }
			  if (KD_gggrav > 0.15) {
			    evtPassgrMELA = true;
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      if ((sel>=100) != (candPass&&candIsBest)) {
	cout << "ERROR2: " << theCandLabel << " " << cand->mass() << " " << candPass << " sel: " << sel << endl;
	abort();
      }
      
      //      if (candPass&&candIsBest) sel = 100;
      if (candPass&&candIsBest&&passMz_zz) sel=120;
    }

//     if(evtPassM70 && candIsBest) {
//       // Plot at step 7 (for single resonant)
//       hCandZZCandM70->Fill(ZZMass, Z1Mass, Z2Mass, Z1OtherCombMass, Z2OtherCombMass, ptS, isoPFS, SIPS, 1.);      
//       hCandZZCandM70_w->Fill(ZZMass, Z1Mass, Z2Mass, Z1OtherCombMass, Z2OtherCombMass, ptS, isoPFS, SIPS, PUweight);
//       if (Z2Mass>12) { 
// 	//Final sel plot, but with mass>70 instead of 100
// 	hCandZZCandM100->Fill(ZZMass, Z1Mass, Z2Mass, Z1OtherCombMass, Z2OtherCombMass, ptS, isoPFS, SIPS, 1.);
// 	hCandZZCandM100_w->Fill(ZZMass, Z1Mass, Z2Mass, Z1OtherCombMass, Z2OtherCombMass, ptS, isoPFS, SIPS, PUweight);
//       }
//     }
//     if (candPass&&candIsBest) { //Full selection
//       hCandZZCand->Fill(ZZMass, Z1Mass, Z2Mass, Z1OtherCombMass, Z2OtherCombMass, ptS, isoPFS, SIPS, 1.);
//       hCandZZCand_w->Fill(ZZMass, Z1Mass, Z2Mass, Z1OtherCombMass, Z2OtherCombMass, ptS, isoPFS, SIPS, PUweight); 
//    }
    


    float ZZMassPreFsr = (Z1Ln->p4()+Z1Lp->p4()+Z2Ln->p4()+Z2Lp->p4()).mass();
    if(sel>=100 && evtPassMELA)  { // Add Mela cut to compare with sync wiki
      int nFSR = 0;
      if (Z1dauFSR>=0) ++nFSR;
      if (Z2dauFSR>=0) ++nFSR;

      if (nFSR>0) {
	Nevt_WithFSR+= weight;
	float nominalMass = 125; //FIXME hardcoded
	if (genH) nominalMass =	genH->mass(); //use genH mass in MC
	//	bool isPure = false;
	h_ZZMass_fsr->Fill(ZZMass);
	h_ZZMassPreFsr_fsr->Fill(ZZMassPreFsr);
	float dR1(1000.), dR2(1000.);
	const reco::Candidate* dauWithFSRZ1 = 0;
	const reco::Candidate* dauWithFSRZ2 = 0;
	if (Z1dauFSR >= 0) {
	  dauWithFSRZ1 = cand->daughter("Z1")->daughter(Z1dauFSR);
	  dR1 = ROOT::Math::VectorUtil::DeltaR(dauWithFSRZ1->momentum(),fsrPhotonZ1->momentum());
	}
	if (Z2dauFSR >= 0) {
	  dauWithFSRZ2 = cand->daughter("Z2")->daughter(Z2dauFSR);
	  dR2 = ROOT::Math::VectorUtil::DeltaR(dauWithFSRZ2->momentum(),fsrPhotonZ2->momentum());
	}
//	if (dR1 == 0 || dR2 == 0) cout<<"dR1 = "<<dR1<<", dR2 = "<<dR2<<endl;

	float dRmin = min(dR1, dR2);
	//	float dRmax = max(dR1, dR2);
	h_dR_all->Fill(dRmin); 
	// if (dRmax<1.) ... // plot for second FSR
	if (abs(ZZMass-nominalMass)<abs(ZZMassPreFsr-nominalMass)){
	  //	  isPure=true;
	  Nevt_WithFSRImprove+= weight;
	  h_dR_pure->Fill(dRmin);
	// if (dRmax<1.) ...
	}

	// cout << "FSR " << ievt << " " << ZZMassPreFsr << " " << ZZMass << " " << endl;
	if(nFSR>1){
	  Nevt_With2FSR+= weight;
	}else {
	  Nevt_With1FSR+= weight;
	}

      }
    }
    
      
    //--- Print standard text summary
    cout.setf(ios::fixed);
    cout << setprecision(3);

    //FIXME: to avoid double counting we should check the correct PD
    if (candIsBest && 
	(dumpMC || (!isMC && 
		    (myHelper.PD!="MuEG" || theChannel==EEMM))) // Do not take 4mu, 4e from MuEG
	){

      //      bool unblinded = ((ZZMass > 70 && ZZMass < 110 ) || (ZZMass > 140 && ZZMass < 300 ));
      bool unblinded = (ZZMass > 140 && ZZMass < 300 );
      string sblabel ="b";
      if (unblinded) sblabel = "";

      string label = findLabel(irun, ils, ievt);
      string grp;
      string srun = boost::lexical_cast<string>( irun );

      if (candPass) {
	grp = "%% " + srun + " %100" + sblabel + "% ";
      } else if (candPass70 && Z2Mass>12) {
	grp = "%% " + srun + " %70" + sblabel + "% ";
      } else { // loose sel
        grp = "** " + srun + " % ";
      }
      
      cout << grp << "Candidate in: " << finalStateNiceName(theChannel) << " sel= " << sel;
      if (label!="???") cout << " (" << label <<")";
      cout << endl;

      string withFSR = "";
      if (Z1dauFSR >= 0 || Z2dauFSR >= 0) withFSR="+g";
      

      cout << grp << endl
// 	   << grp << endl
// 	   << grp << "####" << endl 
// 	   << grp << "# " << label << " #" << endl 
// 	   << grp << "####" << endl
	   << grp << endl;
		cout << grp << "run= " << irun << " evt= " << ievt << " ls= " << ils << " m" << finalStateNiceName(theChannel) << withFSR << "= " << setprecision(2) << ZZMass << " mZ1= " <<  Z1Mass << " mZ2= " << Z2Mass << " massError= " << ZZMassErr << " massErrorCorr= " << ZZMassErrCorr
		<< " KD= " << setprecision(3) << KD << setprecision(2) << " pt4l= " << ZZPt  
		     << setprecision(1) << " Njets= " << (float)nJet30 << setprecision(2) << " ptj1= " << setprecision(2) << j1pT << " ptj2= " << j2pT << setprecision(2) << " mjj= " << mjj << setprecision(3) << " deta= " << deta << " VD= " << VD 
		//here the additional KDs
		<< setprecision(3)
		<< " KD_pseudo= " <<  KD_pseudo 
		<< " KD_highdim= " <<  KD_highdim 
		<< " KD_vec= " <<  KD_vec 
		<< " KD_psvec= " <<  KD_psvec 
		<< " KD_gggrav= " <<  KD_gggrav 
		<< " KD_qqgrav= " <<  KD_qqgrav 
	    << " KD_2hp= " << KD_2hp 
		<< " KD_2hm= " << KD_2hm
		<< " KD_2bp= " << KD_2bp
		<< endl;
		
      cout << grp << "Other pairs: mZ'1= " << Z1OtherCombMass << " mZ'2= " << Z2OtherCombMass << endl;
      //      cout << grp << "pt_4l= " << cand->pt() << " eta_4l= " << cand->eta() << " y_4l " << cand->p4().Rapidity()
      //	   << " phi_4l= " << cand->phi() << " N_vtx " << vertexs->size() <<endl;
      cout << grp << "costheta1= " << costheta1 << " costheta2= " << costheta2 << " Phi= " << Phi << " costhetastar= " << costhetastar
	   << " Phi1= " << Phi1
	   << endl;
      cout << grp << "nVtx= " << Nvtx << " PFMET= " << pfmet << endl;
      cout << grp << "b% " << " VBF= " << evtPassVBF << endl;

      //      <<  " rho= " << *rhoHandle 	

      // Print formatted table
      cout << grp; printParticle();
      cout << grp; printParticle(&*cand, "ZZ", 0);
      cout << grp; printParticle(Z1,   "Z1", 0);
      cout << grp; printParticle(Z2,   "Z2", 0);
      const string roles[] = {"L11","L12","L21", "L22"};
      for (int i = 0; i<4; ++i) {
	cout << grp; printParticle(leptons[i], roles[i], 0, isoPF[i], SIP[i]);	
      }
      if (Z1dauFSR >= 0) cout << grp << "L1" << Z1dauFSR+1 << " FSR: pT= " << fsrPhotonZ1->pt() << " eta= " << fsrPhotonZ1->eta() << " phi= " << fsrPhotonZ1->phi() << " mll= " << userdatahelpers::getUserFloat(Z1,"mll") << endl;
      if (Z2dauFSR >= 0) cout << grp << "L2" << Z2dauFSR+1 << " FSR: pT= " << fsrPhotonZ2->pt() << " eta= " << fsrPhotonZ2->eta() << " phi= " << fsrPhotonZ2->phi() << " mll= " << userdatahelpers::getUserFloat(Z2,"mll") << endl;
      if ((Z1dauFSR >= 0)&&(Z2dauFSR >= 0)) cout << grp << "m4l (no FSR)= " << ZZMassPreFsr << endl;
    }

  } // End of loop on candidates


  //----------------------------------------------------------------------
  // Update per-event counters 
  //----------------------------------------------------------------------

  {
    if (evtPassZ1) {
      //	  cout <<"SEL3a " << ievt << endl;
      Nevt_Z1+= weight;	
      if (evtPassZ1Mass) {
	++Nevt_Z1Mass;
	if(evtPassTriGoodLep) Nevt_Z1Plus1Lep+= weight; // Do not make the following dependent on this request, since 
	                                                // at the moment the l in Z+l is isolated WITHOUT FSR. 
	{
	  if(evtPassZZ){ 
	    // FIXME does not work with FSR
// 	    if (evtPassZ1Mass_B!=evtPassZ1Mass) { // Check consistency between Z and ZZ collections
// 	      cout << "ERROR " << " " << evtPassZ1Mass_B << " " << evtPassZ1Mass << endl;
// 	      abort();
// 	    }
	    ++Nevt_ZZCand;
	    if(evtPassIsBestZ1) Nevt_ZZBestCand+= weight;
	    if(evtPassZ1Mass) Nevt_ZZCandMZ1+= weight; //FIXME REDUNDANT
	    if(evtPassZ2Mass) Nevt_ZZCandMZ2+= weight;
	    if(evtPass20_10) Nevt_ZZCandpT+= weight;
	    if(evtPassMLLallComb) Nevt_ZZCandMAllComb+= weight;
	    if(evtPassM70) Nevt_ZZCandM70+= weight;
	    if(evtPassM100) { 
	      Nevt_ZZCandM100+= weight;
//	      cout << "SEL " << irun << " " << ils << " " << ievt << endl;
	    }
		if(evtPass1Jet) Nevt_ZZCandM100_1Jet += weight;  
		if(evtPass2Jets) Nevt_ZZCandM100_2Jets += weight;  
		if(evtPassVBF) Nevt_ZZCandM100_VBF += weight;  
	    if(evtPassM100_zz) Nevt_ZZCand_zz+= weight;

	    if(evtPassMELA) {
	      Nevt_MELA+= weight;
	    }
	    if (evtPasspsMELA) {
	      Nevt_psMELA+= weight;
	    }
	    if (evtPassgrMELA) {
	      Nevt_grMELA+= weight;
	    }
	  }
	}
      }
    }
  }
}  



void ZZ4lAnalyzer::printParticle(const reco::Candidate* c, string idx, int pdgId, 
				 float iso, float SIP) {
	
  const int precision = 3;
  const int colw      = 9;
  int pprec = cout.precision();
  
  if (c==0) { //Header
    cout << setw (4) << "    " << setw (4) << "Id" 
	 << setw (colw) << "pt" << setw (colw) << "eta" << setw (colw) << "phi" 
	 << setw (colw) << "px" << setw (colw)<< "py" << setw (colw) << "pz" 
      //	 << setw (colw) << "E" 
	 << setw (colw) << "y" << setw (colw) << "iso"  << setw (colw) << "SIP"
	 << endl;

    
  } else {
    cout << setprecision (precision);
    if (pdgId==0) pdgId = c->pdgId();
    cout << setw (4) << idx << setw (4) << pdgId 
	 << setw (colw) << c->pt() << setw (colw) << c->eta() << setw (colw) << c->phi()
	 << setw (colw) << c->px() << setw (colw) << c->py() << setw (colw) << c->pz() 
//	 << setw (colw) << c->energy() 
	 << setw (colw) << c->p4().Rapidity();
    if (SIP>-0.5) { 
      cout << setw (colw) << iso << setw (colw) << SIP;
    } else {
      cout << setw (colw) << "-" << setw (colw) << "-";      
    }
    
    
    cout << endl;
    cout << setprecision (pprec);
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ZZ4lAnalyzer);

