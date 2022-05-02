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
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/Math/interface/LorentzVector.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/JetReco/interface/PFJet.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <DataFormats/JetReco/interface/PFJetCollection.h>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include <FWCore/ParameterSet/interface/FileInPath.h>

//ATjets Additional libraries for GenJet variables
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <JetMETCorrections/Objects/interface/JetCorrectionsRecord.h>
#include <JetMETCorrections/Modules/interface/JetResolution.h>


#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <DataFormats/Common/interface/MergeableCounter.h>

#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>
#include <ZZAnalysis/AnalysisStep/interface/MCHistoryTools.h>
#include <ZZAnalysis/AnalysisStep/interface/PileUpWeight.h>
#include <ZZAnalysis/AnalysisStep/interface/Fisher.h>
#include "ZZ4lConfigHelper.h"
#include <boost/lexical_cast.hpp>

#include <iostream>
#include <fstream>
#include <iterator>
#include <string>

#include <Math/VectorUtil.h>

#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"

#define HCMSSW
#include "ZZAnalysis/AnalysisStep/interface/Histograms.h"
#include <ZZAnalysis/AnalysisStep/interface/bitops.h>
#include <algorithm>

using namespace std;
using namespace edm;
using namespace reco;
using namespace userdatahelpers;

namespace {
//float sqr(float x) {return x*x;}

  bool lesspT(PhotonPtr a, PhotonPtr b) {
    return (a->pt()<b->pt());
  }
}


class ZZ4lAnalyzer: public edm::EDAnalyzer {
public:

  explicit ZZ4lAnalyzer(const edm::ParameterSet& pset);
  ~ZZ4lAnalyzer();

  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
  void endLuminosityBlock(const edm::LuminosityBlock &, const edm::EventSetup&);
  virtual void beginJob();
  virtual void endJob();

  void printParticle(const reco::Candidate* c=0, string idx="", int pdgId=0,
		     float iso=-1, float SIP=-1);

private:
  ZZ4lConfigHelper myHelper;
  std::string theCandLabel;
  bool isMC;
  bool dumpMC;
  Float_t GenZZMass;
  Float_t genHEPMCweight;
  PileUpWeight* pileUpReweight;

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
  double Nevt_passDiMu;
  double Nevt_passDiEle;
  double Nevt_passMuEle;
  double Nevt_passTriEle;
  double Nevt_passTriMu;
  double Nevt_passSingleEle;

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

  //histograms
  //  TH1F* nEventComplete;
//   TH1F* hZ1Mass;
//   TH1F* hZ1MassEBEB;
//   TH1F* hZ1MassEEEE;
//   TH1F* hZ1Mass_w;
//   TH1F* hZ1MassEBEB_w;
//   TH1F* hZ1MassEEEE_w;

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

  TH1F* hNvtxNoWeight;
  TH1F* hNvtxWeight;
  TH1F* hNTrueIntNoWeight;
  TH1F* hNTrueIntWeight;
  TH1F* hRhoNoWeight;
  TH1F* hRhoWeight;

  TH1F* hGenZZMass;
  TH1F* hGenZZMass_4e;
  TH1F* hGenZZMass_4mu;
  TH1F* hGenZZMass_2e2mu;
  TH1F* hGenZZMass_2l2tau;

  edm::EDGetTokenT<vector<Vertex> > vtxToken;
  edm::EDGetTokenT<double> rhoToken;
  edm::EDGetTokenT<edm::View<reco::Candidate> > genParticleToken;
  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken;
  edm::EDGetTokenT<edm::View<reco::GenJet> > genJetsToken; //ATjets
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedgenParticlesToken; //ATbbf
  //edm::EDGetTokenT<edm::View<pat::Jet> > jetToken;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultToken;
  edm::EDGetTokenT<edm::View<reco::Candidate> > softLeptonToken;

  edm::EDGetTokenT<edm::MergeableCounter> preSkimToken;

  string sampleName;
  bool dumpForSync;
  ofstream leptonSyncFile;

};

// Constructor
ZZ4lAnalyzer::ZZ4lAnalyzer(const ParameterSet& pset) :
  myHelper(pset),
  theCandLabel(pset.getUntrackedParameter<string>("candCollection")),
  dumpMC(pset.getUntrackedParameter<bool>("dumpMC",false)),
  GenZZMass(0.),
  genHEPMCweight(0.),
  pileUpReweight(nullptr),
  Nevt_Gen(0),
  gen_ZZ4e(0),
  gen_ZZ4mu(0),
  gen_ZZ2mu2e(0),
  gen_ZZtau(0),
  gen_BUGGY(0),
  Nevt_PAT(0),
  Nevt_Skim(0),
  Nevt_TriggerBit(0),
  Nevt_passDiMu(0),
  Nevt_passDiEle(0),
  Nevt_passMuEle(0),
  Nevt_passTriEle(0),
  Nevt_passTriMu(0),
  Nevt_passSingleEle(0),
  Nevt_4l(0),
  sampleName(pset.getParameter<string>("sampleName")),
  dumpForSync(pset.getUntrackedParameter<bool>("dumpForSync",false))
{

  vtxToken = consumes<vector<Vertex> >(edm::InputTag("goodPrimaryVertices"));
  rhoToken = consumes<double>(edm::InputTag("fixedGridRhoFastjetAll",""));
  genParticleToken = consumes<edm::View<reco::Candidate> >( edm::InputTag("prunedGenParticles"));
  genInfoToken = consumes<GenEventInfoProduct>( edm::InputTag("generator"));
  genJetsToken = consumes<edm::View<reco::GenJet> >(edm::InputTag("slimmedGenJets")); //AT jets (Word between "" not so sure, BBF puts "genJetsSrc")
  packedgenParticlesToken = consumes<edm::View<pat::PackedGenParticle> > (edm::InputTag("packedGenParticles")); //ATbbf
  consumesMany<std::vector< PileupSummaryInfo > >();
  //jetToken = consumes<edm::View<pat::Jet> >(edm::InputTag("cleanJets"));
  triggerResultToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults"));
  softLeptonToken = consumes<edm::View<reco::Candidate> >(edm::InputTag("softLeptons"));

  preSkimToken = consumes<edm::MergeableCounter,edm::InLumi>(edm::InputTag("preSkimCounter"));

  isMC = myHelper.isMC();
  if (isMC) pileUpReweight = new PileUpWeight(myHelper.sampleType(), myHelper.setup());
}

ZZ4lAnalyzer::~ZZ4lAnalyzer()
{
  delete pileUpReweight;
}

void ZZ4lAnalyzer::beginJob(){
  // Book some control histograms
  edm::Service<TFileService> fileService;

  hNvtxNoWeight = fileService->make<TH1F>("hNvtxNoWeight","hNvtxNoWeight",100,0,100);
  hNvtxWeight = fileService->make<TH1F>("hNvtxWeight","hNvtxWeight",100,0,100);
  hNTrueIntNoWeight = fileService->make<TH1F>("hNTrueIntNoWeight","hNTrueIntNoWeight",100,0,100);
  hNTrueIntWeight = fileService->make<TH1F>("hNTrueIntWeight","hNTrueIntWeight",100,0,100);
  hRhoNoWeight = fileService->make<TH1F>("hRhoNoWeight","hRhoNoWeight",100,0,50);
  hRhoWeight = fileService->make<TH1F>("hRhoWeight","hRhoWeight",100,0,50);

  hGenZZMass = fileService->make<TH1F>("hGenZZMass","hGenZZMass",13000,0,13000);
  hGenZZMass_4e     = fileService->make<TH1F>("hGenZZMass_4e"    ,"hGenZZMass_4e"    ,13000,0,13000);
  hGenZZMass_4mu    = fileService->make<TH1F>("hGenZZMass_4mu"   ,"hGenZZMass_4mu"   ,13000,0,13000);
  hGenZZMass_2e2mu  = fileService->make<TH1F>("hGenZZMass_2e2mu" ,"hGenZZMass_2e2mu" ,13000,0,13000);
  hGenZZMass_2l2tau = fileService->make<TH1F>("hGenZZMass_2l2tau","hGenZZMass_2l2tau",13000,0,13000);

  // Counting Histograms. We do not want sumw2 on these!
  TH1F::SetDefaultSumw2(kFALSE);
  //  nEventComplete = fileService->make<TH1F>("nEventComplete", "nEventComplete", 12, 0., 12.);
//   hZ1Mass = fileService->make<TH1F>("hZ1Mass", "hZ1Mass", 400, 0., 200.);
//   hZ1Mass_w = fileService->make<TH1F>("hZ1Mass_w", "hZ1Mass_w", 400, 0., 200.);

//   hZ1MassEBEB = fileService->make<TH1F>("hZ1MassEBEB", "hZ1MassEBEB", 400, 0., 200.);
//   hZ1MassEEEE = fileService->make<TH1F>("hZ1MassEEEE", "hZ1MassEEEE", 400, 0., 200.);
//   hZ1MassEBEB_w = fileService->make<TH1F>("hZ1MassEBEB_w", "hZ1MassEBEB_w", 400, 0., 200.);
//   hZ1MassEEEE_w = fileService->make<TH1F>("hZ1MassEEEE_w", "hZ1MassEEEE_w", 400, 0., 200.);



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


  if (dumpForSync) {
    leptonSyncFile.open("leptons.txt");
    leptonSyncFile.setf(ios::fixed);
  }

}


void ZZ4lAnalyzer::endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup)
{

  double Nevt_preskim = -1.;
  edm::Handle<edm::MergeableCounter> preSkimCounter;
  if (iLumi.getByToken(preSkimToken, preSkimCounter)) { // Counter before skim. Does not exist for non-skimmed samples.
    Nevt_preskim = preSkimCounter->value;
  }

  // Nevt_gen: this is the number before any skim
  if (Nevt_preskim>=0.) {
    Nevt_Gen += Nevt_preskim;
  } else {
    //    Nevt_Gen = Nevt_Gen + prePathCounter->value;
  }

  // Nevt_PAT: number of events in the pattuple
  //  Nevt_PAT = Nevt_PAT + prePathCounter->value;


  // cout << "Nevt_Gen_lumi: " << Nevt_Gen << endl;
  // cout << "Nevt_afterSkim_lumi: " << Nevt_afterSkim << endl;

}



void ZZ4lAnalyzer::endJob(){

  cout << endl;
  cout << "*************************" <<endl;
  cout << "Nevt_Gen:         " << Nevt_Gen << endl;
  if (isMC) {
    cout << "gen_BUGGY:        " << gen_BUGGY << endl;
    cout << "gen_ZZ4mu:        " << gen_ZZ4mu << endl;
    cout <<  "gen_ZZ4e:         " << gen_ZZ4e << endl;
    cout <<  "gen_ZZ2mu2e:      " << gen_ZZ2mu2e << endl;
    cout <<  "gen_ZZtau:        " << gen_ZZtau << endl;
  }
  cout <<  "Nevt_PAT:         " << Nevt_PAT << endl;
  cout <<  "Nevt_Skim:        " << Nevt_Skim << endl;
  cout <<  "Nevt_TriggerBit:  " << Nevt_TriggerBit
       << " ( " << Nevt_passDiEle << ":" << Nevt_passDiMu  << ":" << Nevt_passMuEle << ":" << Nevt_passTriEle << ":"
       << Nevt_passTriMu<< ":" << Nevt_passSingleEle << " )" << endl;


  if (false) { // 2012 selection
    //----------------------------------------------------------------------
    // nEvent Histos for Step Plot: "Nominal Signal Selection"
    //----------------------------------------------------------------------

//     nEventComplete->GetXaxis()->SetBinLabel(1,"Skim");             // HZZ Skim
//     nEventComplete->GetXaxis()->SetBinLabel(2,"HLT");              // TriggerBit
//     nEventComplete->GetXaxis()->SetBinLabel(3,"Z1");               // At least one Z1
//     nEventComplete->GetXaxis()->SetBinLabel(4,"m_{Z1}");           // At least one Z1Mass
//     nEventComplete->GetXaxis()->SetBinLabel(5,"Z1+l");             // At least one Z1+l (good lepton)
//     nEventComplete->GetXaxis()->SetBinLabel(6,"ZZ");               // At least one LLLL candidate with Z1
//     nEventComplete->GetXaxis()->SetBinLabel(7,"m_{Z2}");           // At least one LLLL candidate with 4 < Z2M < 12
//     nEventComplete->GetXaxis()->SetBinLabel(8,"pT>20,10");         // At least one LLLL candidate with pT 20, 10
//     nEventComplete->GetXaxis()->SetBinLabel(9,"m_{ll}>4");         // At least one LLLL candidate with 6/6 mll > 4
//     nEventComplete->GetXaxis()->SetBinLabel(10,"m_{4l}>70");       // At least one LLLL candidate with mlll > 70
//     nEventComplete->GetXaxis()->SetBinLabel(11,"m_{4l}>100");      // At least one LLLL candidate with mlll > 100
//     nEventComplete->GetXaxis()->SetBinLabel(12,"ZZ sel.");              // High-mass (ZZ selection)

//     nEventComplete->SetBinContent(1,Nevt_Skim);
//     nEventComplete->SetBinContent(2,Nevt_TriggerBit);
//     nEventComplete->SetBinContent(3,Nevt_Z1);
//     nEventComplete->SetBinContent(4,Nevt_Z1Mass);
//     nEventComplete->SetBinContent(5,Nevt_Z1Plus1Lep);
//     nEventComplete->SetBinContent(6,Nevt_ZZBestCand);
//     nEventComplete->SetBinContent(7,Nevt_ZZCandMZ2);
//     nEventComplete->SetBinContent(8,Nevt_ZZCandpT);
//     nEventComplete->SetBinContent(9,Nevt_ZZCandMAllComb);
//     nEventComplete->SetBinContent(10,Nevt_ZZCandM70);
//     nEventComplete->SetBinContent(11,Nevt_ZZCandM100);
//     nEventComplete->SetBinContent(12,Nevt_ZZCand_zz);

//     if (isMC){
//       nEventComplete->SetBinContent(0,Nevt_Gen-gen_BUGGY); // Save normalization factor in underflow bin (correcting for buggy MC events)
//     } else {
//       nEventComplete->SetBinContent(0,1.);
//     }
  }

  if (dumpForSync) {
    leptonSyncFile.close();
  }

}



void ZZ4lAnalyzer::analyze(const Event & event, const EventSetup& eventSetup){
  int irun=event.id().run();
  long long int ievt=event.id().event();
  int ils =event.luminosityBlock();

  // int nObsInt  = -1;
  float nTrueInt = -1.;
  float PUweight = 1.;

  Handle<vector<Vertex> > vertices;
  event.getByToken(vtxToken,vertices);
  int Nvtx = vertices->size();

  Handle<double> rhoHandle;
  event.getByToken(rhoToken, rhoHandle);

  int genFinalState = NONE;
  if (isMC) {

    edm::Handle<edm::View<reco::Candidate> > genParticles;
    event.getByToken(genParticleToken, genParticles);
    edm::Handle<GenEventInfoProduct> genInfo;
    event.getByToken(genInfoToken, genInfo);
    edm::Handle<edm::View<reco::GenJet> > genJets; //ATjets
    event.getByToken(genJetsToken, genJets); //ATjets
    edm::Handle<edm::View<pat::PackedGenParticle> > packedgenParticles; //ATbbf
    event.getByToken(packedgenParticlesToken, packedgenParticles); //ATbbf

    MCHistoryTools mch(event, sampleName, genParticles, genInfo, genJets, packedgenParticles);
    genFinalState = mch.genFinalState();

    const reco::Candidate* genH = mch.genH();
    std::vector<const reco::Candidate *> genZLeps = mch.sortedGenZZLeps();
    if(genH!=0){
      GenZZMass = genH->mass();
    }else if(genZLeps.size()==4){ // for bkgd 4l events, take the mass of the ZZ(4l) system
      GenZZMass = (genZLeps.at(0)->p4()+genZLeps.at(1)->p4()+genZLeps.at(2)->p4()+genZLeps.at(3)->p4()).M();
    }
    genHEPMCweight = mch.gethepMCweight();

    hGenZZMass->Fill( GenZZMass, genHEPMCweight*PUweight );
    if (genFinalState == EEEE) {
      ++gen_ZZ4e;
      hGenZZMass_4e->Fill( GenZZMass, genHEPMCweight*PUweight );
    } else if (genFinalState == MMMM) {
      ++gen_ZZ4mu;
      hGenZZMass_4mu->Fill( GenZZMass, genHEPMCweight*PUweight );
    } else if (genFinalState == EEMM) {
      ++gen_ZZ2mu2e;
      hGenZZMass_2e2mu->Fill( GenZZMass, genHEPMCweight*PUweight );
    } else if (genFinalState == LLTT || genFinalState == llTT || genFinalState == TTTT){
      ++gen_ZZtau;
      hGenZZMass_2l2tau->Fill( GenZZMass, genHEPMCweight*PUweight );
    } else if (genFinalState == BUGGY){
      ++gen_BUGGY;
      return;
    }

    vector<Handle<std::vector< PileupSummaryInfo > > >  PupInfos; //FIXME support for miniAOD v1/v2 where name changed; catch does not work...
    event.getManyByType(PupInfos);
    Handle<std::vector< PileupSummaryInfo > > PupInfo = PupInfos.front();
//     try {
//       cout << "TRY HZZ4lNtupleMaker" <<endl;
//       event.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
//     } catch (const cms::Exception& e){
//       cout << "FAIL HZZ4lNtupleMaker" <<endl;
//       event.getByLabel(edm::InputTag("slimmedAddPileupInfo"), PupInfo);
//     }

    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      if(PVI->getBunchCrossing() == 0) {
	//	nObsInt  = PVI->getPU_NumInteractions();
	nTrueInt = PVI->getTrueNumInteractions();
	break;
      }
    }

    PUweight = pileUpReweight->weight(nTrueInt);

    // These plots are intended to cross-check the PU reweighting procedure, 
    // by comparison of the PU-reweighted and original distributions.
    hNvtxWeight->Fill(Nvtx,PUweight);
    hNTrueIntWeight->Fill(nTrueInt,PUweight);
    hRhoWeight->Fill(*rhoHandle,PUweight);

  }
  hNvtxNoWeight->Fill(Nvtx);
  hNTrueIntNoWeight->Fill(nTrueInt);
  hRhoNoWeight->Fill(*rhoHandle);

  // theChannel Candidates
//   Handle<View<pat::CompositeCandidate> > candHandle;
//   event.getByLabel(theCandLabel, candHandle);
//   const View<pat::CompositeCandidate>* cands = candHandle.product();

  // Muons
  // Handle<pat::MuonCollection> muons;
  // event.getByLabel("softMuons", muons);

  // Electrons
  // Handle<pat::ElectronCollection> electrons;
  // event.getByLabel("cleanSoftElectrons", electrons);

  // Zs (for step plots)
  // Handle<View<pat::CompositeCandidate> > Z;
  // event.getByLabel("ZCand", Z);

  // Z+l (for step plots)
//   Handle<View<reco::CompositeCandidate> > Zl;
//   event.getByLabel("ZlCand", Zl);
//   const View<reco::CompositeCandidrate>* trileps = Zl.product();

  //  float pfmet = -1;
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
//   Handle<vector<reco::MET> > pfmetcoll;
//   event.getByLabel("slimmedMETs", pfmetcoll);
//   if(pfmetcoll.isValid()){
//     pfmet = pfmetcoll->front().pt();
//   }

  // Jet collection
  // Handle<edm::View<pat::Jet> > CleanJets;
  // event.getByToken(jetToken, CleanJets);
  // vector<const pat::Jet*> cleanedJetsPt30;
  // for(edm::View<pat::Jet>::const_iterator jet = CleanJets->begin(); jet != CleanJets->end(); ++jet){
  //   if(jet->pt()>30) cleanedJetsPt30.push_back(&*jet);
  // }

  // Trigger results
  Handle<edm::TriggerResults> triggerResults;
  event.getByToken(triggerResultToken, triggerResults);

  // Apply MC filter (skip event)
  if (isMC && !(myHelper.passMCFilter(event,triggerResults))) return;

  // Skim
  bool evtPassSkim = myHelper.passSkim(event,triggerResults);

  // Trigger requests
  short trigworld = 0;
  bool evtPassTrigger = myHelper.passTrigger(event,triggerResults,trigworld);
  if(test_bit_16(trigworld,1)) ++Nevt_passDiMu;
  if(test_bit_16(trigworld,2)) ++Nevt_passDiEle;
  if(test_bit_16(trigworld,3)) ++Nevt_passMuEle;
  if(test_bit_16(trigworld,4)) ++Nevt_passTriEle;
  if(test_bit_16(trigworld,5)) ++Nevt_passTriMu;
  if(test_bit_16(trigworld,6)) ++Nevt_passSingleEle;

  // print a status
  //if (!isMC) cout << irun << ":" << ils << ":" << ievt << " " << evtPassSkim << " " << evtPassTrigger << endl;

  if(!evtPassSkim) return;
  ++Nevt_Skim;
  if(!evtPassTrigger) return;
  ++Nevt_TriggerBit;

  // Fill Lepton sync dump
  if (dumpForSync) {
    Handle<View<reco::Candidate> > softleptoncoll;
    event.getByToken(softLeptonToken, softleptoncoll);

    // {run}:{lumi}:{event}:{pdgId:d}:{pT:.2f}:{eta:.2f}:{phi:.2f}{SIP:.2f}:{PFChargedHadIso:.2f}:{PFNeutralHadIso:.2f}:{PFPhotonIso:.2f}:{rho:.2f}:{combRelIsoPF:.3f}:{eleBDT:.3f}

    for (View<reco::Candidate>::const_iterator lep = softleptoncoll->begin(); lep != softleptoncoll->end(); ++ lep) {
      if((bool)getUserFloat(&*lep,"ID")){ // Tight leptons (no SIP or ISO)
	int id = lep->pdgId();
	leptonSyncFile << irun << ":" << ils << ":" << ievt << ":" << id << ":"
		       << setprecision(2) << lep->pt() << ":" << lep->eta() << ":" << lep->phi() << ":" << getUserFloat(&*lep,"SIP") << ":"
		       << getUserFloat(&*lep,"PFChargedHadIso") << ":" << getUserFloat(&*lep,"PFNeutralHadIso") << ":"
		       << getUserFloat(&*lep,"PFPhotonIso") << ":"
		       << (abs(id)==11?getUserFloat(&*lep,"rho"):getUserFloat(&*lep,"PFPUChargedHadIso")) << ":"
		       << setprecision(3) << getUserFloat(&*lep,"combRelIsoPF") << ":" << getUserFloat(&*lep,"BDT");


	const PhotonPtrVector* gammas = getUserPhotons(&*lep);
	if (gammas) {
	  const pat::PFParticle* gamma = (std::max_element(gammas->begin(), gammas->end(), lesspT))->get();
	  double dR = ROOT::Math::VectorUtil::DeltaR(gamma->momentum(),lep->momentum());
	  leptonSyncFile << ":" << setprecision(2) << gamma->pt() << ":" << dR; // << ":" << gRelIso; Must be recomputed
	}
	leptonSyncFile << endl;
      }
    }
  }

  return;

  //----------------------------------------------------------------------
  // Loop on ZZ candidates


//   for( View<pat::CompositeCandidate>::const_iterator cand = cands->begin(); cand != cands->end(); ++cand) {

//     //--- Print standard text summary
//     cout.setf(ios::fixed);
//     cout << setprecision(3);

//     //FIXME: to avoid double counting we should check the correct PD
//     if (candIsBest &&
// 	(dumpMC || (!isMC &&
// 		    (myHelper.PD!="MuEG" || theChannel==EEMM))) // Do not take 4mu, 4e from MuEG
// 	){

//       //      bool unblinded = ((ZZMass > 70 && ZZMass < 110 ) || (ZZMass > 140 && ZZMass < 300 ));
//       bool unblinded = (ZZMass > 140 && ZZMass < 300 );
//       string sblabel ="b";
//       if (unblinded) sblabel = "";

//       string grp;
//       string srun = boost::lexical_cast<string>( irun );

//       if (candPass) {
// 	grp = "%% " + srun + " %100" + sblabel + "% ";
//       } else if (candPass70 && Z2Mass>12) {
// 	grp = "%% " + srun + " %70" + sblabel + "% ";
//       } else { // loose sel
//         grp = "** " + srun + " % ";
//       }

//       cout << grp << "Candidate in: " << finalStateNiceName(theChannel) << " sel= " << sel;
//       cout << endl;

//       string withFSR = "";
//       if (Z1dauFSR >= 0 || Z2dauFSR >= 0) withFSR="+g";


//       cout << grp << endl
// // 	   << grp << endl
// // 	   << grp << "####" << endl
// // 	   << grp << "# " << label << " #" << endl
// // 	   << grp << "####" << endl
// 	   << grp << endl;
// 		cout << grp << "run= " << irun << " evt= " << ievt << " ls= " << ils << " m" << finalStateNiceName(theChannel) << withFSR << "= " << setprecision(2) << ZZMass << " mZ1= " <<  Z1Mass << " mZ2= " << Z2Mass << " massError= " << ZZMassErr << " massErrorCorr= " << ZZMassErrCorr
// 		<< " KD= " << setprecision(3) << KD << setprecision(2) << " pt4l= " << ZZPt
// 		     << setprecision(1) << " Njets= " << (float)nJet30 << setprecision(2) << " ptj1= " << setprecision(2) << j1pT << " ptj2= " << j2pT << setprecision(2) << " mjj= " << mjj << setprecision(3) << " deta= " << deta << " VD= " << VD
// 		//here the additional KDs
// 		<< setprecision(3)
// 		<< " KD_pseudo= " <<  KD_pseudo
// 		<< " KD_highdim= " <<  KD_highdim
// 		<< " KD_vec= " <<  KD_vec
// 		<< " KD_psvec= " <<  KD_psvec
// 		<< " KD_gggrav= " <<  KD_gggrav
// 		<< " KD_qqgrav= " <<  KD_qqgrav
// 	    << " KD_2hp= " << KD_2hp
// 		<< " KD_2hm= " << KD_2hm
// 		<< " KD_2bp= " << KD_2bp
// 		<< endl;

//       cout << grp << "Other pairs: mZ'1= " << Z1OtherCombMass << " mZ'2= " << Z2OtherCombMass << endl;
//       //      cout << grp << "pt_4l= " << cand->pt() << " eta_4l= " << cand->eta() << " y_4l " << cand->p4().Rapidity()
//       //	   << " phi_4l= " << cand->phi() << " N_vtx " << vertices->size() <<endl;
//       cout << grp << "costheta1= " << costheta1 << " costheta2= " << costheta2 << " Phi= " << Phi << " costhetastar= " << costhetastar
// 	   << " Phi1= " << Phi1
// 	   << endl;
//       cout << grp << "nVtx= " << Nvtx << " PFMET= " << pfmet << endl;
//       cout << grp << "b% " << " VBF= " << evtPassVBF << endl;

//       //      <<  " rho= " << *rhoHandle

//       // Print formatted table
//       cout << grp; printParticle();
//       cout << grp; printParticle(&*cand, "ZZ", 0);
//       cout << grp; printParticle(Z1,   "Z1", 0);
//       cout << grp; printParticle(Z2,   "Z2", 0);
//       const string roles[] = {"L11","L12","L21", "L22"};
//       for (int i = 0; i<4; ++i) {
// 	cout << grp; printParticle(leptons[i], roles[i], 0, isoPF[i], SIP[i]);
//       }
//       if (Z1dauFSR >= 0) cout << grp << "L1" << Z1dauFSR+1 << " FSR: pT= " << fsrPhotonZ1->pt() << " eta= " << fsrPhotonZ1->eta() << " phi= " << fsrPhotonZ1->phi() << " mll= " << userdatahelpers::getUserFloat(Z1,"mll") << endl;
//       if (Z2dauFSR >= 0) cout << grp << "L2" << Z2dauFSR+1 << " FSR: pT= " << fsrPhotonZ2->pt() << " eta= " << fsrPhotonZ2->eta() << " phi= " << fsrPhotonZ2->phi() << " mll= " << userdatahelpers::getUserFloat(Z2,"mll") << endl;
//       if ((Z1dauFSR >= 0)&&(Z2dauFSR >= 0)) cout << grp << "m4l (no FSR)= " << ZZMassPreFsr << endl;
//     }

//   } // End of loop on candidates
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
