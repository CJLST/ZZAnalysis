// system include files
#include <memory>
#include <cmath>
#include <string>
#include <cassert>

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
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <ZZAnalysis/AnalysisStep/interface/METCorrectionHandler.h>
#include <ZZAnalysis/AnalysisStep/interface/PhotonIDHelper.h>
#include <DataFormats/Math/interface/LorentzVector.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <DataFormats/Common/interface/MergeableCounter.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include <SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h>
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>

#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>
#include <ZZAnalysis/AnalysisStep/interface/PileUpWeight.h>
#include <ZZAnalysis/AnalysisStep/interface/bitops.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonSFHelper.h>
#include "ZZ4lConfigHelper.h"
#include "HZZ4lNtupleFactory.h"

#include "Math/VectorUtil.h"
#include "TH1.h"

using namespace std;
using namespace edm;

namespace {
  bool skipMuDataMCWeight = false; // skip computation of data/MC weight for mu
  bool skipEleDataMCWeight = false; // skip computation of data/MC weight for ele
  bool addJets = true;
  bool addPhotons = true;
  bool addQGLInputs = false;
  bool addAddLept = false;

  //List of variables with default values
  Int_t RunNumber = 0;
  Long64_t EventNumber = 0;
  Int_t LumiNumber = 0;
  Short_t Nvtx = 0;
  Short_t NObsInt = 0;
  Float_t NTrueInt = 0;
  Float_t PUWeight = 0;
  Short_t trigWord = 0;
  Short_t Zsel = 0;
  Float_t ZMass = 0;
  Float_t ZMassPreFSR = 0;
  Float_t ZPt = 0;
  Float_t ZEta = 0;
  Float_t ZPhi = 0;
  Short_t ZFlav = 0;
  Short_t evtPassMETTrigger = 0;
  std::vector<float> LepPt;
  std::vector<float> LepEta;
  std::vector<float> LepPhi;
  std::vector<float> LepSCEta;
  std::vector<short> LepLepId;
  std::vector<float> LepSIP;
  std::vector<float> Lepdxy;
  std::vector<float> Lepdz;
  std::vector<float> LepTime;
  std::vector<bool> LepisID;
  std::vector<float> LepBDT;
  std::vector<bool> LepisCrack;
  std::vector<char> LepMissingHit;
  std::vector<float> LepChargedHadIso;
  std::vector<float> LepNeutralHadIso;
  std::vector<float> LepPhotonIso;
  std::vector<float> LepCombRelIsoPF;
  std::vector<float> LepCombRelIsoPFPreFSR;
  std::vector<float> LepScale_Total_Up;
  std::vector<float> LepScale_Total_Dn;
  std::vector<float> LepScale_Stat_Up;
  std::vector<float> LepScale_Stat_Dn;
  std::vector<float> LepScale_Syst_Up;
  std::vector<float> LepScale_Syst_Dn;
  std::vector<float> LepScale_Gain_Up;
  std::vector<float> LepScale_Gain_Dn;
  std::vector<float> LepSigma_Total_Up;
  std::vector<float> LepSigma_Total_Dn;
  std::vector<float> LepSigma_Rho_Up;
  std::vector<float> LepSigma_Rho_Dn;
  std::vector<float> LepSigma_Phi_Up;
  std::vector<float> LepSigma_Phi_Dn;
  std::vector<float> fsrPt;
  std::vector<float> fsrEta;
  std::vector<float> fsrPhi;
  std::vector<float> fsrDR;
  std::vector<short> fsrLept;
  std::vector<short> fsrLeptID;
  Short_t nAddEle = 0;
  std::vector<float> AddElePt;
  std::vector<float> AddEleEta;
  std::vector<float> AddElePhi;
  std::vector<short> AddEleLepId;
  std::vector<float> AddEleSIP;
  std::vector<bool> AddEleisID;
  std::vector<float> AddEleBDT;
  std::vector<char> AddEleMissingHit;
  std::vector<float> AddEleChargedHadIso;
  std::vector<float> AddEleNeutralHadIso;
  std::vector<float> AddElePhotonIso;
  std::vector<float> AddEleCombRelIsoPF;
  Short_t nAddMu = 0;
  std::vector<float> AddMuPt;
  std::vector<float> AddMuEta;
  std::vector<float> AddMuPhi;
  std::vector<short> AddMuLepId;
  std::vector<float> AddMuSIP;
  std::vector<float> AddMuTime;
  std::vector<bool> AddMuisID;
  std::vector<float> AddMuChargedHadIso;
  std::vector<float> AddMuNeutralHadIso;
  std::vector<float> AddMuPhotonIso;
  std::vector<float> AddMuCombRelIsoPF;
  Short_t nCleanedJetsPt30 = 0;
  Short_t nCleanedJetsPt30_jesUp = 0;
  Short_t nCleanedJetsPt30_jesUp_Total = 0;
  Short_t nCleanedJetsPt30_jesUp_Abs = 0;
  Short_t nCleanedJetsPt30_jesUp_Abs_year = 0;
  Short_t nCleanedJetsPt30_jesUp_BBEC1 = 0;
  Short_t nCleanedJetsPt30_jesUp_BBEC1_year = 0;
  Short_t nCleanedJetsPt30_jesUp_EC2 = 0;
  Short_t nCleanedJetsPt30_jesUp_EC2_year = 0;
  Short_t nCleanedJetsPt30_jesUp_FlavQCD = 0;
  Short_t nCleanedJetsPt30_jesUp_HF = 0;
  Short_t nCleanedJetsPt30_jesUp_HF_year = 0;
  Short_t nCleanedJetsPt30_jesUp_RelBal = 0;
  Short_t nCleanedJetsPt30_jesUp_RelSample_year = 0;
  Short_t nCleanedJetsPt30_jesDn = 0;
  Short_t nCleanedJetsPt30_jesDn_Total = 0;
  Short_t nCleanedJetsPt30_jesDn_Abs = 0;
  Short_t nCleanedJetsPt30_jesDn_Abs_year = 0;
  Short_t nCleanedJetsPt30_jesDn_BBEC1 = 0;
  Short_t nCleanedJetsPt30_jesDn_BBEC1_year = 0;
  Short_t nCleanedJetsPt30_jesDn_EC2 = 0;
  Short_t nCleanedJetsPt30_jesDn_EC2_year = 0;
  Short_t nCleanedJetsPt30_jesDn_FlavQCD = 0;
  Short_t nCleanedJetsPt30_jesDn_HF = 0;
  Short_t nCleanedJetsPt30_jesDn_HF_year = 0;
  Short_t nCleanedJetsPt30_jesDn_RelBal = 0;
  Short_t nCleanedJetsPt30_jesDn_RelSample_year = 0;
  Short_t nCleanedJetsPt30_jerUp = 0;
  Short_t nCleanedJetsPt30_jerDn = 0;
  Short_t nCleanedJetsPt30BTagged  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_Total = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_Abs = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_Abs_year = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_BBEC1 = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_BBEC1_year = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_EC2 = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_EC2_year = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_FlavQCD = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_HF = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_HF_year = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_RelBal = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp_RelSample_year = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_Total = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_Abs = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_Abs_year = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_BBEC1 = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_BBEC1_year = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_EC2 = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_EC2_year = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_FlavQCD = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_HF = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_HF_year = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_RelBal = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn_RelSample_year = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jerUp  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jerDn  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSFUp  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSFDn  = 0;
  std::vector<float> JetPt;
  std::vector<float> JetEta;
  std::vector<float> JetPhi;
  std::vector<float> JetMass;
  std::vector<float> JetEnergy;
  std::vector<float> JetBTagger;
  std::vector<float> JetIsBtagged;
  std::vector<float> JetIsBtaggedWithSF;
  std::vector<float> JetIsBtaggedWithSFUp;
  std::vector<float> JetIsBtaggedWithSFDn;
  std::vector<float> JetQGLikelihood;
  std::vector<float> JetAxis2;
  std::vector<float> JetMult;
  std::vector<float> JetPtD;
  std::vector<float> JetSigma;
  std::vector<float> JetSigma_Total ;
  std::vector<float> JetSigma_Abs ;
  std::vector<float> JetSigma_Abs_year ;
  std::vector<float> JetSigma_BBEC1 ;
  std::vector<float> JetSigma_BBEC1_year ;
  std::vector<float> JetSigma_EC2 ;
  std::vector<float> JetSigma_EC2_year ;
  std::vector<float> JetSigma_FlavQCD ;
  std::vector<float> JetSigma_HF ;
  std::vector<float> JetSigma_HF_year ;
  std::vector<float> JetSigma_RelBal ;
  std::vector<float> JetSigma_RelSample_year ;
  std::vector<short> JetHadronFlavour;
  std::vector<short> JetPartonFlavour;
  std::vector<float> JetPUValue;
  std::vector<short> JetPUID;
  std::vector<short> JetID;
  std::vector<float> JetJESUp ;
  std::vector<float> JetJESUp_Total ;
  std::vector<float> JetJESUp_Abs ;
  std::vector<float> JetJESUp_Abs_year ;
  std::vector<float> JetJESUp_BBEC1 ;
  std::vector<float> JetJESUp_BBEC1_year ;
  std::vector<float> JetJESUp_EC2 ;
  std::vector<float> JetJESUp_EC2_year ;
  std::vector<float> JetJESUp_FlavQCD ;
  std::vector<float> JetJESUp_HF ;
  std::vector<float> JetJESUp_HF_year ;
  std::vector<float> JetJESUp_RelBal ;
  std::vector<float> JetJESUp_RelSample_year ;
  std::vector<float> JetJESDown ;
  std::vector<float> JetJESDown_Total ;
  std::vector<float> JetJESDown_Abs ;
  std::vector<float> JetJESDown_Abs_year ;
  std::vector<float> JetJESDown_BBEC1 ;
  std::vector<float> JetJESDown_BBEC1_year ;
  std::vector<float> JetJESDown_EC2 ;
  std::vector<float> JetJESDown_EC2_year ;
  std::vector<float> JetJESDown_FlavQCD ;
  std::vector<float> JetJESDown_HF ;
  std::vector<float> JetJESDown_HF_year ;
  std::vector<float> JetJESDown_RelBal ;
  std::vector<float> JetJESDown_RelSample_year ;
  std::vector<float> JetPt_JERUp;
  std::vector<float> JetPt_JERDown;
  std::vector<float> JetPtJEC_noJER;
  std::vector<float> JetRawPt;
   
  std::vector<float> PhotonPt;
  std::vector<float> PhotonEta;
  std::vector<float> PhotonPhi;
  std::vector<bool> PhotonIsCutBasedLooseID;

	
  // Generic MET object
  METObject metobj;
  METObject metobj_corrected;
  Float_t GenMET = -99;
  Float_t GenMETPhi = -99;

  Float_t genHEPMCweight = 0;
  Float_t xsection = 0;
  Float_t dataMCWeight = 0;
  Float_t overallEventWeight = 0;
   
  Float_t L1prefiringWeight = 0;
  Float_t L1prefiringWeightUp = 0;
  Float_t L1prefiringWeightDn = 0;
}

//
// class declaration
//
class ZNtupleMaker : public edm::EDAnalyzer {
public:
  explicit ZNtupleMaker(const edm::ParameterSet&);
  ~ZNtupleMaker();
	
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
private:
  virtual void beginJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  void BookAllBranches();
  virtual void FillCandidate(const pat::CompositeCandidate& higgs, bool evtPass, const edm::Event&);
  virtual void FillJet(const pat::Jet& jet);
  virtual void FillPhoton(int year, const pat::Photon& photon);
  virtual void endJob();

  Float_t getAllWeight(const reco::Candidate* Lep);

  // ----------member data ---------------------------
  ZZ4lConfigHelper myHelper;
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
  edm::InputTag metTag;
   
  METCorrectionHandler* metCorrHandler;

  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken;
  edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > candToken;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultToken;
  edm::EDGetTokenT<vector<reco::Vertex> > vtxToken;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetToken;
  edm::EDGetTokenT<pat::METCollection> metToken;
  edm::EDGetTokenT<pat::PhotonCollection> photonToken;

  edm::EDGetTokenT<edm::MergeableCounter> preSkimToken;

  edm::EDGetTokenT<vector<pat::Electron> > electronToken;
  const vector<pat::Electron>* softElectrons;
  edm::EDGetTokenT<vector<pat::Muon> > muonToken;
  const vector<pat::Muon>* softMuons;
   
  edm::EDGetTokenT< double > prefweight_token;
  edm::EDGetTokenT< double > prefweightup_token;
  edm::EDGetTokenT< double > prefweightdown_token;

  PileUpWeight pileUpReweight;

  //counters
  Float_t Nevt_Gen;
  Float_t Nevt_Gen_lumiBlock;
	
  Float_t gen_sumPUWeight;
  Float_t gen_sumGenMCWeight;
  Float_t gen_sumWeights;

  string sampleName;

  LeptonSFHelper *lepSFHelper;

};

//
// constructors and destructor
//
ZNtupleMaker::ZNtupleMaker(const edm::ParameterSet& pset) :
  myHelper(pset),
  pileUpReweight(myHelper.sampleType(), myHelper.setup())
{
  theCandLabel = pset.getUntrackedParameter<string>("CandCollection"); // Name of input Z collection
  theFileName = pset.getUntrackedParameter<string>("fileName");
  skipEmptyEvents = pset.getParameter<bool>("skipEmptyEvents"); // Do not store
  sampleName = pset.getParameter<string>("sampleName");
  xsec = pset.getParameter<double>("xsec");
  year = pset.getParameter<int>("setup");
  metTag = pset.getParameter<edm::InputTag>("metSrc");

  consumesMany<std::vector< PileupSummaryInfo > >();
  genInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  candToken = consumes<edm::View<pat::CompositeCandidate> >(edm::InputTag(theCandLabel));
  triggerResultToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults"));
  vtxToken = consumes<vector<reco::Vertex> >(edm::InputTag("goodPrimaryVertices"));
  jetToken = consumes<edm::View<pat::Jet> >(edm::InputTag("cleanJets"));
  
  metToken = consumes<pat::METCollection>(metTag);
  metCorrHandler = new METCorrectionHandler(Form("%i", year));
  photonToken = consumes<pat::PhotonCollection>(edm::InputTag("slimmedPhotons"));

  electronToken = consumes<vector<pat::Electron> >(edm::InputTag("softElectrons"));
  //softElectrons->clear();
  muonToken = consumes<vector<pat::Muon> >(edm::InputTag("softMuons"));
  //softMuons->clear();

  preSkimToken = consumes<edm::MergeableCounter,edm::InLumi>(edm::InputTag("preSkimCounter"));

  if (skipEmptyEvents) {
    applyTrigger=true;
    applySkim=true;
  } else {
    applyTrigger=false;
    applySkim=false;
  }

  isMC = myHelper.isMC();
   
  if( isMC && (year == 2016 || year == 2017))
  {
     prefweight_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProb"));
     prefweightup_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbUp"));
     prefweightdown_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbDown"));
  }

  Nevt_Gen = 0;
  Nevt_Gen_lumiBlock = 0;

  gen_sumPUWeight = 0.f;
  gen_sumGenMCWeight = 0.f;
  gen_sumWeights =0.f;

   //Scale factors for data/MC efficiency
   if (!skipEleDataMCWeight && isMC) { lepSFHelper = new LeptonSFHelper(); }
	
}


ZNtupleMaker::~ZNtupleMaker()
{
   delete metCorrHandler;
}


// ------------ method called for each event  ------------
void ZNtupleMaker::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  myTree->InitializeVariables();


  if(isMC){

    // get PU weights
    vector<Handle<std::vector< PileupSummaryInfo > > >  PupInfos; //FIXME support for miniAOD v1/v2 where name changed; catch does not work...
    event.getManyByType(PupInfos);
    Handle<std::vector< PileupSummaryInfo > > PupInfo = PupInfos.front();
	  
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      if(PVI->getBunchCrossing() == 0) {
        NObsInt  = PVI->getPU_NumInteractions();
        NTrueInt = PVI->getTrueNumInteractions();
        break;
      }
    }

    PUWeight = pileUpReweight.weight(NTrueInt);

    edm::Handle<GenEventInfoProduct> gen;
    event.getByToken(genInfoToken, gen);
    GenEventInfoProduct  genInfo = *(gen.product());
    genHEPMCweight = genInfo.weight();

    // keep track of sum of weights
    gen_sumPUWeight    += PUWeight;
    gen_sumGenMCWeight += genHEPMCweight;
    gen_sumWeights     += PUWeight*genHEPMCweight;
     
    // L1 prefiring weights
    if( year == 2016 || year == 2017 )
    {
       edm::Handle< double > theprefweight;
       event.getByToken(prefweight_token, theprefweight ) ;
       L1prefiringWeight =(*theprefweight);
        
       edm::Handle< double > theprefweightup;
       event.getByToken(prefweightup_token, theprefweightup ) ;
       L1prefiringWeightUp =(*theprefweightup);
        
       edm::Handle< double > theprefweightdown;
       event.getByToken(prefweightdown_token, theprefweightdown ) ;
       L1prefiringWeightDn =(*theprefweightdown);
    }
    else if ( year == 2018 )
    {
       L1prefiringWeight   = 1.;
       L1prefiringWeightUp = 1.;
       L1prefiringWeightDn = 1.;
    }

  }


  // Get candidate collection
  edm::Handle<edm::View<pat::CompositeCandidate> > candHandle;
  event.getByToken(candToken, candHandle);
  const edm::View<pat::CompositeCandidate>* cands = candHandle.product();

  if (skipEmptyEvents && cands->size() == 0) return; // Skip events with no candidate, unless skipEmptyEvents = false


  // Retrieve trigger results
  Handle<edm::TriggerResults> triggerResults;
  event.getByToken(triggerResultToken, triggerResults);

  // Apply skim
  bool evtPassSkim = myHelper.passSkim(event,triggerResults,trigWord);
  if (applySkim && !evtPassSkim) return;

  // Apply trigger request (skip event)
  bool evtPassTrigger = myHelper.passTrigger(event,triggerResults,trigWord);
  if (applyTrigger && !evtPassTrigger) return;
	
  // Apply MET trigger request (skip event)
  evtPassMETTrigger = myHelper.passMETTrigger(event,triggerResults);


  // General event information
  RunNumber = event.id().run();
  LumiNumber = event.luminosityBlock();
  EventNumber = event.id().event();
  xsection = xsec;

  // Primary vertices
  Handle<vector<reco::Vertex> > vertices;
  event.getByToken(vtxToken,vertices);
  Nvtx = vertices->size();

   // Jets (cleaned wrt all tight isolated leptons)
   Handle<edm::View<pat::Jet> > CleanedJets;
   event.getByToken(jetToken, CleanedJets);
   vector<const pat::Jet*> cleanedJets;
   for(edm::View<pat::Jet>::const_iterator jet = CleanedJets->begin(); jet != CleanedJets->end(); ++jet){
      cleanedJets.push_back(&*jet);
   }
   
   // Count and store jets, after additional cleaning for CRs...
   for (unsigned i=0; i<cleanedJets.size(); ++i) {
       
     // count jes up/down njets pt30
     float jes_unc = cleanedJets[i]->userFloat("jes_unc");
     float jes_unc_Total = cleanedJets[i]->userFloat("jes_unc_split_Total");
     float jes_unc_Abs = cleanedJets[i]->userFloat("jes_unc_split_Abs");
     float jes_unc_Abs_year = cleanedJets[i]->userFloat("jes_unc_split_Abs_year");
     float jes_unc_BBEC1 = cleanedJets[i]->userFloat("jes_unc_split_BBEC1");
     float jes_unc_BBEC1_year = cleanedJets[i]->userFloat("jes_unc_split_BBEC1_year");
     float jes_unc_EC2 = cleanedJets[i]->userFloat("jes_unc_split_EC2");
     float jes_unc_EC2_year = cleanedJets[i]->userFloat("jes_unc_split_EC2_year");
     float jes_unc_FlavQCD = cleanedJets[i]->userFloat("jes_unc_split_FlavQCD");
     float jes_unc_HF = cleanedJets[i]->userFloat("jes_unc_split_HF");
     float jes_unc_HF_year = cleanedJets[i]->userFloat("jes_unc_split_HF_year");
     float jes_unc_RelBal = cleanedJets[i]->userFloat("jes_unc_split_RelBal");
     float jes_unc_RelSample_year = cleanedJets[i]->userFloat("jes_unc_split_RelSample_year");

     float pt_nominal = cleanedJets[i]->pt();
     float pt_jes_up = pt_nominal * (1.0 + jes_unc);
     float pt_jes_up_Total = pt_nominal * (1.0 + jes_unc_Total);
     float pt_jes_up_Abs = pt_nominal * (1.0 + jes_unc_Abs);
     float pt_jes_up_Abs_year = pt_nominal * (1.0 + jes_unc_Abs_year);
     float pt_jes_up_BBEC1 = pt_nominal * (1.0 + jes_unc_BBEC1);
     float pt_jes_up_BBEC1_year = pt_nominal * (1.0 + jes_unc_BBEC1_year);
     float pt_jes_up_EC2 = pt_nominal * (1.0 + jes_unc_EC2);
     float pt_jes_up_EC2_year = pt_nominal * (1.0 + jes_unc_EC2_year);
     float pt_jes_up_FlavQCD = pt_nominal * (1.0 + jes_unc_FlavQCD);
     float pt_jes_up_HF = pt_nominal * (1.0 + jes_unc_HF);
     float pt_jes_up_HF_year = pt_nominal * (1.0 + jes_unc_HF_year);
     float pt_jes_up_RelBal = pt_nominal * (1.0 + jes_unc_RelBal);
     float pt_jes_up_RelSample_year = pt_nominal * (1.0 + jes_unc_RelSample_year);
     float pt_jes_dn = pt_nominal * (1.0 - jes_unc);
     float pt_jes_dn_Total = pt_nominal * (1.0 - jes_unc_Total);
     float pt_jes_dn_Abs = pt_nominal * (1.0 - jes_unc_Abs);
     float pt_jes_dn_Abs_year = pt_nominal * (1.0 - jes_unc_Abs_year);
     float pt_jes_dn_BBEC1 = pt_nominal * (1.0 - jes_unc_BBEC1);
     float pt_jes_dn_BBEC1_year = pt_nominal * (1.0 - jes_unc_BBEC1_year);
     float pt_jes_dn_EC2 = pt_nominal * (1.0 - jes_unc_EC2);
     float pt_jes_dn_EC2_year = pt_nominal * (1.0 - jes_unc_EC2_year);
     float pt_jes_dn_FlavQCD = pt_nominal * (1.0 - jes_unc_FlavQCD);
     float pt_jes_dn_HF = pt_nominal * (1.0 - jes_unc_HF);
     float pt_jes_dn_HF_year = pt_nominal * (1.0 - jes_unc_HF_year);
     float pt_jes_dn_RelBal = pt_nominal * (1.0 - jes_unc_RelBal);
     float pt_jes_dn_RelSample_year = pt_nominal * (1.0 - jes_unc_RelSample_year);

     if(pt_nominal>30){
       ++nCleanedJetsPt30;
       if(cleanedJets[i]->userFloat("isBtagged")) ++nCleanedJetsPt30BTagged;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF_Up")) ++nCleanedJetsPt30BTagged_bTagSFUp;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF_Dn")) ++nCleanedJetsPt30BTagged_bTagSFDn;
     }
     if(pt_jes_up>30){
       ++nCleanedJetsPt30_jesUp;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp;
     }
     if(pt_jes_up_Total>30){
       ++nCleanedJetsPt30_jesUp_Total;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_Total;
     }
     if(pt_jes_up_Abs>30){
       ++nCleanedJetsPt30_jesUp_Abs;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_Abs;
     }
     if(pt_jes_up_Abs_year>30){
       ++nCleanedJetsPt30_jesUp_Abs_year;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_Abs_year;
     }
     if(pt_jes_up_BBEC1>30){
       ++nCleanedJetsPt30_jesUp_BBEC1;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_BBEC1;
     }
     if(pt_jes_up_BBEC1_year>30){
       ++nCleanedJetsPt30_jesUp_BBEC1_year;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_BBEC1_year;
     }
     if(pt_jes_up_EC2>30){
       ++nCleanedJetsPt30_jesUp_EC2;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_EC2;
     }
     if(pt_jes_up_EC2_year>30){
       ++nCleanedJetsPt30_jesUp_EC2_year;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_EC2_year;
     }
     if(pt_jes_up_FlavQCD>30){
       ++nCleanedJetsPt30_jesUp_FlavQCD;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_FlavQCD;
     }
     if(pt_jes_up_HF>30){
       ++nCleanedJetsPt30_jesUp_HF;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_HF;
     }
     if(pt_jes_up_HF_year>30){
       ++nCleanedJetsPt30_jesUp_HF_year;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_HF_year;
     }
     if(pt_jes_up_RelBal>30){
       ++nCleanedJetsPt30_jesUp_RelBal;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_RelBal;
     }
     if(pt_jes_up_RelSample_year>30){
       ++nCleanedJetsPt30_jesUp_RelSample_year;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp_RelSample_year;
     }

     if(pt_jes_dn>30){
       ++nCleanedJetsPt30_jesDn;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn;
     }
     if(pt_jes_dn_Total>30){
       ++nCleanedJetsPt30_jesDn_Total;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_Total;
     }
     if(pt_jes_dn_Abs>30){
       ++nCleanedJetsPt30_jesDn_Abs;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_Abs;
     }
     if(pt_jes_dn_Abs_year>30){
       ++nCleanedJetsPt30_jesDn_Abs_year;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_Abs_year;
     }
     if(pt_jes_dn_BBEC1>30){
       ++nCleanedJetsPt30_jesDn_BBEC1;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_BBEC1;
     }
     if(pt_jes_dn_BBEC1_year>30){
       ++nCleanedJetsPt30_jesDn_BBEC1_year;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_BBEC1_year;
     }
     if(pt_jes_dn_EC2>30){
       ++nCleanedJetsPt30_jesDn_EC2;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_EC2;
     }
     if(pt_jes_dn_EC2_year>30){
       ++nCleanedJetsPt30_jesDn_EC2_year;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_EC2_year;
     }
     if(pt_jes_dn_FlavQCD>30){
       ++nCleanedJetsPt30_jesDn_FlavQCD;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_FlavQCD;
     }
     if(pt_jes_dn_HF>30){
       ++nCleanedJetsPt30_jesDn_HF;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_HF;
     }
     if(pt_jes_dn_HF_year>30){
       ++nCleanedJetsPt30_jesDn_HF_year;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_HF_year;
     }
     if(pt_jes_dn_RelBal>30){
       ++nCleanedJetsPt30_jesDn_RelBal;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_RelBal;
     }
     if(pt_jes_dn_RelSample_year>30){
       ++nCleanedJetsPt30_jesDn_RelSample_year;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn_RelSample_year;
     }

     
     // count jer up/down njets pt30
     float pt_jer_up = cleanedJets[i]->userFloat("pt_jerup");
     float pt_jer_dn = cleanedJets[i]->userFloat("pt_jerdn");
     
     if(pt_jer_up>30){
       ++nCleanedJetsPt30_jerUp;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jerUp;
     }
     if(pt_jer_dn>30){
       ++nCleanedJetsPt30_jerDn;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jerDn;
     }
     if (addJets) FillJet(*(cleanedJets.at(i)));
   }
   
  // MET
  Handle<pat::METCollection> metHandle;
  event.getByToken(metToken, metHandle);
  GenMET=GenMETPhi=-99;
  if (metHandle.isValid()){
     const pat::MET &met = metHandle->front();
     
     metobj.extras.met = metobj.extras.met_original = metobj.extras.met_raw
     = metobj.extras.met_METup = metobj.extras.met_METdn
     = metobj.extras.met_JERup = metobj.extras.met_JERdn
     = metobj.extras.met_PUup = metobj.extras.met_PUdn
     
     = metobj_corrected.extras.met = metobj_corrected.extras.met_original = metobj_corrected.extras.met_raw
     = metobj_corrected.extras.met_METup = metobj_corrected.extras.met_METdn
     = metobj_corrected.extras.met_JERup = metobj_corrected.extras.met_JERdn
     = metobj_corrected.extras.met_PUup = metobj_corrected.extras.met_PUdn
     
     = met.pt();
     metobj.extras.phi = metobj.extras.phi_original = metobj.extras.phi_raw
     = metobj.extras.phi_METup = metobj.extras.phi_METdn
     = metobj.extras.phi_JECup = metobj.extras.phi_JECdn
     = metobj.extras.phi_JERup = metobj.extras.phi_JERdn
     = metobj.extras.phi_PUup = metobj.extras.phi_PUdn
     
     = metobj_corrected.extras.phi = metobj_corrected.extras.phi_original = metobj_corrected.extras.phi_raw
     = metobj_corrected.extras.phi_METup = metobj_corrected.extras.phi_METdn
     = metobj_corrected.extras.phi_JECup = metobj_corrected.extras.phi_JECdn
     = metobj_corrected.extras.phi_JERup = metobj_corrected.extras.phi_JERdn
     = metobj_corrected.extras.phi_PUup = metobj_corrected.extras.phi_PUdn
     
     = met.phi();
     
     metobj.extras.met_JECup = metobj_corrected.extras.met_JECup = met.shiftedPt(pat::MET::JetEnUp);
     metobj.extras.met_JECdn = metobj_corrected.extras.met_JECdn = met.shiftedPt(pat::MET::JetEnDown);
     metobj.extras.phi_JECup = metobj_corrected.extras.phi_JECup = met.shiftedPhi(pat::MET::JetEnUp);
     metobj.extras.phi_JECdn = metobj_corrected.extras.phi_JECdn = met.shiftedPhi(pat::MET::JetEnDown);
     
     if (isMC && metCorrHandler && met.genMET()){
        GenMET = met.genMET()->pt();
        GenMETPhi = met.genMET()->phi();
        metCorrHandler->correctMET(GenMET, GenMETPhi, &metobj_corrected, false); // FIXME: Last argument should be for isFastSim, but we don't have it yet
     }
     else if (isMC){
        cms::Exception e("METCorrectionHandler");
        e << "Either no met.genMET or metCorrHandler!";
        throw e;
     }
  }
  else{
     metobj.extras.met = metobj.extras.met_original = metobj.extras.met_raw
     = metobj.extras.met_METup = metobj.extras.met_METdn
     = metobj.extras.met_JECup = metobj.extras.met_JECdn
     = metobj.extras.met_JERup = metobj.extras.met_JERdn
     = metobj.extras.met_PUup = metobj.extras.met_PUdn
     
     = metobj_corrected.extras.met = metobj_corrected.extras.met_original = metobj_corrected.extras.met_raw
     = metobj_corrected.extras.met_METup = metobj_corrected.extras.met_METdn
     = metobj_corrected.extras.met_JECup = metobj_corrected.extras.met_JECdn
     = metobj_corrected.extras.met_JERup = metobj_corrected.extras.met_JERdn
     = metobj_corrected.extras.met_PUup = metobj_corrected.extras.met_PUdn
     
     = metobj.extras.phi = metobj.extras.phi_original = metobj.extras.phi_raw
     = metobj.extras.phi_METup = metobj.extras.phi_METdn
     = metobj.extras.phi_JECup = metobj.extras.phi_JECdn
     = metobj.extras.phi_JERup = metobj.extras.phi_JERdn
     = metobj.extras.phi_PUup = metobj.extras.phi_PUdn
     
     = metobj_corrected.extras.phi = metobj_corrected.extras.phi_original = metobj_corrected.extras.phi_raw
     = metobj_corrected.extras.phi_METup = metobj_corrected.extras.phi_METdn
     = metobj_corrected.extras.phi_JECup = metobj_corrected.extras.phi_JECdn
     = metobj_corrected.extras.phi_JERup = metobj_corrected.extras.phi_JERdn
     = metobj_corrected.extras.phi_PUup = metobj_corrected.extras.phi_PUdn
     
     = -99;
  }
   
  // Photons
  Handle<pat::PhotonCollection> photonCands;
  event.getByToken(photonToken, photonCands);
  vector<const pat::Photon*> photons;
   
  for(unsigned int i = 0; i< photonCands->size(); ++i){
     const pat::Photon* photon = &((*photonCands)[i]);
     photons.push_back(&*photon);
  }
   
   
  if (addPhotons){
     for (unsigned i=0; i<photons.size(); ++i) {
        FillPhoton(year, *(photons.at(i)));
     }
  }
	
  // all soft leptons
  edm::Handle<vector<pat::Electron> > electronHandle;
  event.getByToken(electronToken, electronHandle);
  softElectrons = electronHandle.product();
  edm::Handle<vector<pat::Muon> > muonHandle;
  event.getByToken(muonToken, muonHandle);
  softMuons = muonHandle.product();

  // Loop on candidates
  int nFilled = 0;
  for( edm::View<pat::CompositeCandidate>::const_iterator cand = cands->begin(); cand != cands->end(); ++cand) {

    // let's take only the selected Z's now, since we don't plot other events anyway
    if(!((bool)(cand->userFloat("isBestZ")))) continue;

    FillCandidate(*cand, evtPassTrigger&&evtPassSkim, event);
    // Fill the candidate as one entry in the tree. Do not reinitialize the event variables, as there could be several candidates per event.
    myTree->FillCurrentTree();
    ++nFilled;
  }
	
  // If no candidate was filled but we still want to keep gen-level and weights, we need to fill one entry anyhow.
  if (skipEmptyEvents==false && nFilled==0) myTree->FillCurrentTree();
}


void ZNtupleMaker::FillCandidate(const pat::CompositeCandidate& cand, bool evtPass, const edm::Event& event)
{
  //Initialize a new candidate into the tree
  //myTree->createNewCandidate();

  //Fill the info on the Z candidate
  ZMass = cand.p4().mass();
  if (cand.numberOfDaughters()>2) ZMassPreFSR = cand.userFloat("mll");
  ZPt  = cand.p4().pt();
  ZEta = cand.p4().eta();
  ZPhi = cand.p4().phi();
  ZFlav = abs(cand.daughter(0)->pdgId()) * cand.daughter(0)->charge() * abs(cand.daughter(1)->pdgId()) * cand.daughter(1)->charge(); // FIXME: temporarily changed, waiting for a fix to the mismatch of charge() and pdgId() for muons with BTT=4

  // Precomputed selections
  bool candPassGoodLeptons = cand.userFloat("GoodLeptons");
  bool candPassZ1Presel = cand.userFloat("Z1Presel");
  bool candisBestZ = cand.userFloat("isBestZ");
  if(candPassGoodLeptons) Zsel = 10;
  if(!candPassZ1Presel &&  candisBestZ) Zsel = 30;
  if( candPassZ1Presel && !candisBestZ) Zsel = 40;
  if( candPassZ1Presel &&  candisBestZ) Zsel = 50;
  if(!evtPass) Zsel = -Zsel; // avoid confusion when we write events which do not pass trigger/skim

  // Retrieve the userFloats of the leptons in vectors ordered in the same way.
  vector<const reco::Candidate*> leptons;
  vector<const reco::Candidate*> fsrPhot;
  vector<short> fsrIndex;
  vector<string> labels;
  userdatahelpers::getSortedZLeptons(cand, leptons, labels, fsrPhot, fsrIndex);
	
  LepPt.clear();
  LepEta.clear();
  LepPhi.clear();
  LepSCEta.clear();
  LepLepId.clear();
  LepSIP.clear();
  Lepdxy.clear();
  Lepdz.clear();
  LepTime.clear();
  LepisID.clear();
  LepBDT.clear();
  LepisCrack.clear();
  LepScale_Total_Up.clear();
  LepScale_Total_Dn.clear();
  LepScale_Stat_Up.clear();
  LepScale_Stat_Dn.clear();
  LepScale_Syst_Up.clear();
  LepScale_Syst_Dn.clear();
  LepScale_Gain_Up.clear();
  LepScale_Gain_Dn.clear();
  LepSigma_Total_Up.clear();
  LepSigma_Total_Dn.clear();
  LepSigma_Rho_Up.clear();
  LepSigma_Rho_Dn.clear();
  LepSigma_Phi_Up.clear();
  LepSigma_Phi_Dn.clear();
  LepMissingHit.clear();
  LepChargedHadIso.clear();
  LepNeutralHadIso.clear();
  LepPhotonIso.clear();
  LepCombRelIsoPF.clear();
  LepCombRelIsoPFPreFSR.clear();

  for (unsigned int i=0; i<leptons.size(); ++i){
    short lepFlav = std::abs(leptons[i]->pdgId());

    // Check that I don't mess up with labels[] and leptons[]
    assert(userdatahelpers::getUserFloat(leptons[i],"SIP") == cand.userFloat(labels[i]+"SIP"));

    //Fill the info on the lepton candidates
    LepPt .push_back( leptons[i]->pt() );
    LepEta.push_back( leptons[i]->eta() );
    LepPhi.push_back( leptons[i]->phi() );
    LepSCEta.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"SCeta") : -99. );
    LepLepId.push_back( leptons[i]->pdgId() );
    LepSIP  .push_back( userdatahelpers::getUserFloat(leptons[i],"SIP") );
    Lepdxy  .push_back( userdatahelpers::getUserFloat(leptons[i],"dxy") );
    Lepdz   .push_back( userdatahelpers::getUserFloat(leptons[i],"dz") );
    LepTime .push_back( lepFlav==13 ? userdatahelpers::getUserFloat(leptons[i],"time") : 0. );
    LepisID .push_back( userdatahelpers::getUserFloat(leptons[i],"ID") );
    LepBDT  .push_back( userdatahelpers::getUserFloat(leptons[i],"BDT"));
    LepisCrack.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"isCrack") : 0 );
    LepMissingHit.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"missingHit") : 0 );
    LepScale_Total_Up.push_back( userdatahelpers::getUserFloat(leptons[i],"scale_total_up") );
    LepScale_Total_Dn.push_back( userdatahelpers::getUserFloat(leptons[i],"scale_total_dn") );
    LepScale_Stat_Up.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"scale_stat_up") : -99. );
    LepScale_Stat_Dn.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"scale_stat_dn") : -99. );
    LepScale_Syst_Up.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"scale_syst_up") : -99. );
    LepScale_Syst_Dn.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"scale_syst_dn") : -99. );
    LepScale_Gain_Up.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"scale_gain_up") : -99. );
    LepScale_Gain_Dn.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"scale_gain_dn") : -99. );
    LepSigma_Total_Up.push_back( userdatahelpers::getUserFloat(leptons[i],"sigma_total_up") );
    LepSigma_Total_Dn.push_back( userdatahelpers::getUserFloat(leptons[i],"sigma_total_dn") );
    LepSigma_Rho_Up.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"sigma_rho_up") : -99. );
    LepSigma_Rho_Dn.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"sigma_rho_dn") : -99. );
    LepSigma_Phi_Up.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"sigma_phi_up") : -99. );
    LepSigma_Phi_Dn.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"sigma_phi_dn") : -99. );
    LepChargedHadIso.push_back( userdatahelpers::getUserFloat(leptons[i],"PFChargedHadIso") );
    LepNeutralHadIso.push_back( userdatahelpers::getUserFloat(leptons[i],"PFNeutralHadIso") );
    LepPhotonIso    .push_back( userdatahelpers::getUserFloat(leptons[i],"PFPhotonIso") );
    LepCombRelIsoPF      .push_back( cand.userFloat(labels[i]+"combRelIsoPFFSRCorr") );
    LepCombRelIsoPFPreFSR.push_back( userdatahelpers::getUserFloat(leptons[i],"combRelIsoPF") );
  }

  // FSR
  fsrPt.clear();
  fsrEta.clear();
  fsrPhi.clear();
  fsrLept.clear();
  fsrLeptID.clear();
  fsrDR.clear();
  for (unsigned i=0; i<fsrPhot.size(); ++i) {
    math::XYZTLorentzVector fsr = fsrPhot[i]->p4();
    fsrPt.push_back(fsr.pt());
    fsrEta.push_back(fsr.eta());
    fsrPhi.push_back(fsr.phi());
    fsrLept.push_back(fsrIndex[i]+1);
    fsrLeptID.push_back(leptons[fsrIndex[i]]->pdgId());
    fsrDR.push_back(ROOT::Math::VectorUtil::DeltaR(leptons[fsrIndex[i]]->momentum(), fsrPhot[i]->momentum()));
  }

  //additional leptons
  AddElePt.clear();
  AddEleEta.clear();
  AddElePhi.clear();
  AddEleLepId.clear();
  AddEleSIP.clear();
  AddEleisID.clear();
  AddEleBDT.clear();
  AddEleMissingHit.clear();
  AddEleChargedHadIso.clear();
  AddEleNeutralHadIso.clear();
  AddElePhotonIso.clear();
  AddEleCombRelIsoPF.clear();
  for (unsigned i=0; i<softElectrons->size(); ++i) {
    pat::Electron lep = softElectrons->at(i);
    if(reco::deltaR(lep.p4(), leptons[0]->p4()) > 0.02 && reco::deltaR(lep.p4(), leptons[1]->p4()) > 0.02){
      nAddEle++;
      if (addAddLept) {
	AddElePt .push_back( lep.pt() );
	AddEleEta.push_back( lep.eta() );
	AddElePhi.push_back( lep.phi() );
	AddEleLepId.push_back( lep.pdgId() );
	AddEleSIP  .push_back( lep.userFloat("SIP") );
	AddEleisID .push_back( lep.userFloat("ID") );
	AddEleBDT  .push_back( lep.userFloat("BDT") );
	AddEleMissingHit.push_back( lep.userFloat("missingHit") );
	AddEleChargedHadIso.push_back( lep.userFloat("PFChargedHadIso") );
	AddEleNeutralHadIso.push_back( lep.userFloat("PFNeutralHadIso") );
	AddElePhotonIso    .push_back( lep.userFloat("PFPhotonIso") );
	AddEleCombRelIsoPF .push_back( lep.userFloat("combRelIsoPF") );
      }
    }
  }
  AddMuPt.clear();
  AddMuEta.clear();
  AddMuPhi.clear();
  AddMuLepId.clear();
  AddMuSIP.clear();
  AddMuTime.clear();
  AddMuisID.clear();
  AddMuChargedHadIso.clear();
  AddMuNeutralHadIso.clear();
  AddMuPhotonIso.clear();
  AddMuCombRelIsoPF.clear();
  for (unsigned i=0; i<softMuons->size(); ++i) {
    pat::Muon lep = softMuons->at(i);
    if(reco::deltaR(lep.p4(), leptons[0]->p4()) > 0.02 && reco::deltaR(lep.p4(), leptons[1]->p4()) > 0.02){
      nAddMu++;
      if (addAddLept) {
	AddMuPt .push_back( lep.pt() );
	AddMuEta.push_back( lep.eta() );
	AddMuPhi.push_back( lep.phi() );
	AddMuLepId.push_back( lep.pdgId() );
	AddMuSIP  .push_back( lep.userFloat("SIP") );
	AddMuTime .push_back( lep.userFloat("time") );
	AddMuisID .push_back( lep.userFloat("ID") );
	AddMuChargedHadIso.push_back( lep.userFloat("PFChargedHadIso") );
	AddMuNeutralHadIso.push_back( lep.userFloat("PFNeutralHadIso") );
	AddMuPhotonIso    .push_back( lep.userFloat("PFPhotonIso") );
	AddMuCombRelIsoPF .push_back( lep.userFloat("combRelIsoPF") );
      }
    }
  }

  //Compute the data/MC weight and overall event weight
  dataMCWeight = 1.;
  if(isMC){
    for(unsigned int i=0; i<leptons.size(); ++i){
      dataMCWeight *= getAllWeight(leptons[i]);
    }
  }
  overallEventWeight = PUWeight * genHEPMCweight * dataMCWeight;

}


void ZNtupleMaker::FillJet(const pat::Jet& jet)
{
  JetPt  .push_back( jet.pt());
  JetEta .push_back( jet.eta());
  JetPhi .push_back( jet.phi());
  JetMass .push_back( jet.p4().M());
  JetEnergy .push_back( jet.p4().energy());
  JetBTagger .push_back( jet.userFloat("bTagger"));
  JetIsBtagged .push_back( jet.userFloat("isBtagged"));
  JetIsBtaggedWithSF .push_back( jet.userFloat("isBtaggedWithSF"));
  JetIsBtaggedWithSFUp .push_back( jet.userFloat("isBtaggedWithSF_Up"));
  JetIsBtaggedWithSFDn .push_back( jet.userFloat("isBtaggedWithSF_Dn"));
  JetQGLikelihood .push_back( jet.userFloat("qgLikelihood"));
  if(addQGLInputs){
    JetAxis2 .push_back( jet.userFloat("axis2"));
    JetMult .push_back( jet.userFloat("mult"));
    JetPtD .push_back( jet.userFloat("ptD"));
  }
  JetSigma .push_back(jet.userFloat("jes_unc"));
  JetHadronFlavour .push_back(jet.hadronFlavour());
  JetPartonFlavour .push_back(jet.partonFlavour());
	
  JetJESUp .push_back(jet.userFloat("pt_jesup"));
  JetJESUp_Total .push_back(jet.userFloat("pt_jesup_split_Total"));
  JetJESUp_Abs .push_back(jet.userFloat("pt_jesup_split_Abs"));
  JetJESUp_Abs_year .push_back(jet.userFloat("pt_jesup_split_Abs_year"));
  JetJESUp_BBEC1 .push_back(jet.userFloat("pt_jesup_split_BBEC1"));
  JetJESUp_BBEC1_year .push_back(jet.userFloat("pt_jesup_split_BBEC1_year"));
  JetJESUp_EC2 .push_back(jet.userFloat("pt_jesup_split_EC2"));
  JetJESUp_EC2_year .push_back(jet.userFloat("pt_jesup_split_EC2_year"));
  JetJESUp_FlavQCD .push_back(jet.userFloat("pt_jesup_split_FlavQCD"));
  JetJESUp_HF .push_back(jet.userFloat("pt_jesup_split_HF"));
  JetJESUp_HF_year .push_back(jet.userFloat("pt_jesup_split_HF_year"));
  JetJESUp_RelBal .push_back(jet.userFloat("pt_jesup_split_RelBal"));
  JetJESUp_RelSample_year .push_back(jet.userFloat("pt_jesup_split_RelSample_year"));
  JetJESDown .push_back(jet.userFloat("pt_jesdn"));
  JetJESDown_Total .push_back(jet.userFloat("pt_jesdn_split_Total"));
  JetJESDown_Abs .push_back(jet.userFloat("pt_jesdn_split_Abs"));
  JetJESDown_Abs_year .push_back(jet.userFloat("pt_jesdn_split_Abs_year"));
  JetJESDown_BBEC1 .push_back(jet.userFloat("pt_jesdn_split_BBEC1"));
  JetJESDown_BBEC1_year .push_back(jet.userFloat("pt_jesdn_split_BBEC1_year"));
  JetJESDown_EC2 .push_back(jet.userFloat("pt_jesdn_split_EC2"));
  JetJESDown_EC2_year .push_back(jet.userFloat("pt_jesdn_split_EC2_year"));
  JetJESDown_FlavQCD .push_back(jet.userFloat("pt_jesdn_split_FlavQCD"));
  JetJESDown_HF .push_back(jet.userFloat("pt_jesdn_split_HF"));
  JetJESDown_HF_year .push_back(jet.userFloat("pt_jesdn_split_HF_year"));
  JetJESDown_RelBal .push_back(jet.userFloat("pt_jesdn_split_RelBal"));
  JetJESDown_RelSample_year .push_back(jet.userFloat("pt_jesdn_split_RelSample_year"));
  
  JetPt_JERUp .push_back(jet.userFloat("pt_jerup"));
  JetPt_JERDown .push_back(jet.userFloat("pt_jerdn"));
   
  JetRawPt  .push_back( jet.userFloat("RawPt"));
  JetPtJEC_noJER .push_back( jet.userFloat("pt_JEC_noJER"));
    
  JetID.push_back(jet.userFloat("JetID"));
  JetPUID.push_back(jet.userFloat("PUjetID"));
    
  if (jet.hasUserFloat("pileupJetIdUpdated:fullDiscriminant")) { // if JEC is reapplied, we set this
	  JetPUValue.push_back(jet.userFloat("pileupJetIdUpdated:fullDiscriminant"));
   } else {
     JetPUValue.push_back(jet.userFloat("pileupJetId:fullDiscriminant"));
   }
}

void ZNtupleMaker::FillPhoton(int year, const pat::Photon& photon)
{
   PhotonPt  .push_back( photon.pt());
   PhotonEta .push_back( photon.eta());
   PhotonPhi .push_back( photon.phi());
   
   PhotonIsCutBasedLooseID .push_back( PhotonIDHelper::isCutBasedID_Loose(year, photon) );
}

// ------------ method called once each job just before starting event loop  ------------
void ZNtupleMaker::beginJob()
{
  edm::Service<TFileService> fs;
  myTree = new HZZ4lNtupleFactory( fs->make<TTree>(theFileName,"Event Summary"));
  hCounter = fs->make<TH1F>("Counters", "Counters", 5, 0., 5.);
  BookAllBranches();
}

// ------------ method called once each job just after ending the event loop  ------------
void ZNtupleMaker::endJob()
{
  hCounter->SetBinContent(1,gen_sumWeights);
  hCounter->SetBinContent(2,gen_sumGenMCWeight);
  hCounter->SetBinContent(3,gen_sumPUWeight);

  hCounter->GetXaxis()->SetBinLabel(1,"gen_sumWeights");
  hCounter->GetXaxis()->SetBinLabel(2,"gen_sumGenMCWeight");
  hCounter->GetXaxis()->SetBinLabel(3,"gen_sumPUWeight");

  return;
}

// ------------ method called when starting to processes a run  ------------
void ZNtupleMaker::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void ZNtupleMaker::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void ZNtupleMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
  Nevt_Gen_lumiBlock = 0;
}

// ------------ method called when ending the processing of a luminosity block  ------------
void ZNtupleMaker::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup)
{
  Float_t Nevt_preskim = -1.;
  edm::Handle<edm::MergeableCounter> preSkimCounter;
  if (iLumi.getByToken(preSkimToken, preSkimCounter)) { // Counter before skim. Does not exist for non-skimmed samples.
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
void ZNtupleMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


Float_t ZNtupleMaker::getAllWeight(const reco::Candidate* Lep)
{
 Int_t   myLepID = abs(Lep->pdgId());
 if (skipMuDataMCWeight&& myLepID==13) return 1.;
 if (skipEleDataMCWeight&& myLepID==11) return 1.;

 Float_t myLepPt = Lep->pt();
 Float_t myLepEta = Lep->eta();
   
 float SF = 1.0;
 //float SF_Unc = 0.0;


 Float_t SCeta;
 if (myLepID == 11) SCeta = userdatahelpers::getUserFloat(Lep,"SCeta");
 else SCeta = myLepEta;
   
 Float_t mySCeta;
   
 // Deal with very rare cases when SCeta is out of 2.5 bonds
 if ( myLepEta <= 2.5 && SCeta >= 2.5) mySCeta = 2.49;
 else if ( myLepEta >= -2.5 && SCeta <= -2.5) mySCeta = -2.49;
 else mySCeta = SCeta;
   
 bool isCrack;
 if (myLepID == 11) isCrack = userdatahelpers::getUserFloat(Lep,"isCrack");
 else isCrack = false;
   
   
 SF = lepSFHelper->getSF(year,myLepID,myLepPt,myLepEta, mySCeta, isCrack);
 //SF_Unc = lepSFHelper->getSFError(year,myLepID,myLepPt,myLepEta, mySCeta, isCrack);
   
 return SF;
}


void ZNtupleMaker::BookAllBranches(){
  myTree->Book("RunNumber",RunNumber);
  myTree->Book("EventNumber",EventNumber);
  myTree->Book("LumiNumber",LumiNumber);
  myTree->Book("Nvtx",Nvtx);
  myTree->Book("NObsInt",NObsInt);
  myTree->Book("NTrueInt",NTrueInt);
  myTree->Book("PUWeight",PUWeight);
  myTree->Book("trigWord",trigWord);
  myTree->Book("Zsel",Zsel);
  myTree->Book("ZMass",ZMass);
  myTree->Book("evtPassMETFilter",evtPassMETTrigger);
  myTree->Book("ZMassPreFSR",ZMassPreFSR);
  myTree->Book("ZPt",ZPt);
  myTree->Book("ZEta",ZEta);
  myTree->Book("ZPhi",ZPhi);
  myTree->Book("ZFlav",ZFlav);
  myTree->Book("LepPt",LepPt);
  myTree->Book("LepEta",LepEta);
  myTree->Book("LepPhi",LepPhi);
  myTree->Book("LepSCEta",LepSCEta);
  myTree->Book("LepLepId",LepLepId);
  myTree->Book("LepSIP",LepSIP);
  myTree->Book("Lepdxy",Lepdxy);
  myTree->Book("Lepdz",Lepdz);
  myTree->Book("LepTime",LepTime);
  myTree->Book("LepisID",LepisID);
  myTree->Book("LepBDT",LepBDT);
  myTree->Book("LepisCrack",LepisCrack, false);
  myTree->Book("LepScale_Total_Up",LepScale_Total_Up, false);
  myTree->Book("LepScale_Total_Dn",LepScale_Total_Dn, false);
  myTree->Book("LepScale_Stat_Up",LepScale_Stat_Up, false);
  myTree->Book("LepScale_Stat_Dn",LepScale_Stat_Dn, false);
  myTree->Book("LepScale_Syst_Up",LepScale_Syst_Up, false);
  myTree->Book("LepScale_Syst_Dn",LepScale_Syst_Dn, false);
  myTree->Book("LepScale_Gain_Up",LepScale_Gain_Up, false);
  myTree->Book("LepScale_Gain_Dn",LepScale_Gain_Dn, false);
  myTree->Book("LepSigma_Total_Up",LepSigma_Total_Up, false);
  myTree->Book("LepSigma_Total_Dn",LepSigma_Total_Dn, false);
  myTree->Book("LepSigma_Rho_Up",LepSigma_Rho_Up, false);
  myTree->Book("LepSigma_Rho_Dn",LepSigma_Rho_Dn, false);
  myTree->Book("LepSigma_Phi_Up",LepSigma_Phi_Up, false);
  myTree->Book("LepSigma_Phi_Dn",LepSigma_Phi_Up, false);
  myTree->Book("LepMissingHit",LepMissingHit);
  myTree->Book("LepChargedHadIso",LepChargedHadIso);
  myTree->Book("LepNeutralHadIso",LepNeutralHadIso);
  myTree->Book("LepPhotonIso",LepPhotonIso);
  myTree->Book("LepCombRelIsoPF",LepCombRelIsoPF);
  myTree->Book("LepCombRelIsoPFPreFSR",LepCombRelIsoPFPreFSR);
  myTree->Book("fsrPt",fsrPt);
  myTree->Book("fsrEta",fsrEta);
  myTree->Book("fsrPhi",fsrPhi);
  myTree->Book("fsrLept",fsrLept);
  myTree->Book("fsrDR",fsrDR);
  myTree->Book("fsrLeptId",fsrLeptID);
  myTree->Book("nAddEle",nAddEle);
  myTree->Book("nAddMu",nAddMu);
  if (addAddLept) {
    myTree->Book("AddElePt",AddElePt);
    myTree->Book("AddEleEta",AddEleEta);
    myTree->Book("AddElePhi",AddElePhi);
    myTree->Book("AddEleLepId",AddEleLepId);
    myTree->Book("AddEleSIP",AddEleSIP);
    myTree->Book("AddEleisID",AddEleisID);
    myTree->Book("AddEleBDT",AddEleBDT);
    myTree->Book("AddEleMissingHit",AddEleMissingHit);
    myTree->Book("AddEleChargedHadIso",AddEleChargedHadIso);
    myTree->Book("AddEleNeutralHadIso",AddEleNeutralHadIso);
    myTree->Book("AddElePhotonIso",AddElePhotonIso);
    myTree->Book("AddEleCombRelIsoPF",AddEleCombRelIsoPF);
    myTree->Book("AddMuPt",AddMuPt);
    myTree->Book("AddMuEta",AddMuEta);
    myTree->Book("AddMuPhi",AddMuPhi);
    myTree->Book("AddMuLepId",AddMuLepId);
    myTree->Book("AddMuSIP",AddMuSIP);
    myTree->Book("AddMuTime",AddMuTime);
    myTree->Book("AddMuisID",AddMuisID);
    myTree->Book("AddMuChargedHadIso",AddMuChargedHadIso);
    myTree->Book("AddMuNeutralHadIso",AddMuNeutralHadIso);
    myTree->Book("AddMuPhotonIso",AddMuPhotonIso);
    myTree->Book("AddMuCombRelIsoPF",AddMuCombRelIsoPF);
  }
	
   myTree->Book("nCleanedJetsPt30",nCleanedJetsPt30);
   myTree->Book("nCleanedJetsPt30_jesUp",nCleanedJetsPt30_jesUp);
   myTree->Book("nCleanedJetsPt30_jesUp_Total",nCleanedJetsPt30_jesUp_Total);
   myTree->Book("nCleanedJetsPt30_jesUp_Abs",nCleanedJetsPt30_jesUp_Abs);
   myTree->Book("nCleanedJetsPt30_jesUp_Abs_year",nCleanedJetsPt30_jesUp_Abs_year);
   myTree->Book("nCleanedJetsPt30_jesUp_BBEC1",nCleanedJetsPt30_jesUp_BBEC1);
   myTree->Book("nCleanedJetsPt30_jesUp_BBEC1_year",nCleanedJetsPt30_jesUp_BBEC1_year);
   myTree->Book("nCleanedJetsPt30_jesUp_EC2",nCleanedJetsPt30_jesUp_EC2);
   myTree->Book("nCleanedJetsPt30_jesUp_EC2_year",nCleanedJetsPt30_jesUp_EC2_year);
   myTree->Book("nCleanedJetsPt30_jesUp_FlavQCD",nCleanedJetsPt30_jesUp_FlavQCD);
   myTree->Book("nCleanedJetsPt30_jesUp_HF",nCleanedJetsPt30_jesUp_HF);
   myTree->Book("nCleanedJetsPt30_jesUp_HF_year",nCleanedJetsPt30_jesUp_HF_year);
   myTree->Book("nCleanedJetsPt30_jesUp_RelBal",nCleanedJetsPt30_jesUp_RelBal);
   myTree->Book("nCleanedJetsPt30_jesUp_RelSample_year",nCleanedJetsPt30_jesUp_RelSample_year);
   myTree->Book("nCleanedJetsPt30_jesDn",nCleanedJetsPt30_jesDn);
   myTree->Book("nCleanedJetsPt30_jesDn_Total",nCleanedJetsPt30_jesDn_Total);
   myTree->Book("nCleanedJetsPt30_jesDn_Abs",nCleanedJetsPt30_jesDn_Abs);
   myTree->Book("nCleanedJetsPt30_jesDn_Abs_year",nCleanedJetsPt30_jesDn_Abs_year);
   myTree->Book("nCleanedJetsPt30_jesDn_BBEC1",nCleanedJetsPt30_jesDn_BBEC1);
   myTree->Book("nCleanedJetsPt30_jesDn_BBEC1_year",nCleanedJetsPt30_jesDn_BBEC1_year);
   myTree->Book("nCleanedJetsPt30_jesDn_EC2",nCleanedJetsPt30_jesDn_EC2);
   myTree->Book("nCleanedJetsPt30_jesDn_EC2_year",nCleanedJetsPt30_jesDn_EC2_year);
   myTree->Book("nCleanedJetsPt30_jesDn_FlavQCD",nCleanedJetsPt30_jesDn_FlavQCD);
   myTree->Book("nCleanedJetsPt30_jesDn_HF",nCleanedJetsPt30_jesDn_HF);
   myTree->Book("nCleanedJetsPt30_jesDn_HF_year",nCleanedJetsPt30_jesDn_HF_year);
   myTree->Book("nCleanedJetsPt30_jesDn_RelBal",nCleanedJetsPt30_jesDn_RelBal);
   myTree->Book("nCleanedJetsPt30_jesDn_RelSample_year",nCleanedJetsPt30_jesDn_RelSample_year);
   myTree->Book("nCleanedJetsPt30_jerUp",nCleanedJetsPt30_jerUp);
   myTree->Book("nCleanedJetsPt30_jerDn",nCleanedJetsPt30_jerDn);
   myTree->Book("nCleanedJetsPt30BTagged",nCleanedJetsPt30BTagged);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF",nCleanedJetsPt30BTagged_bTagSF);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp",nCleanedJetsPt30BTagged_bTagSF_jesUp);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_Total",nCleanedJetsPt30BTagged_bTagSF_jesUp_Total);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_Abs",nCleanedJetsPt30BTagged_bTagSF_jesUp_Abs);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_Abs_year",nCleanedJetsPt30BTagged_bTagSF_jesUp_Abs_year);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_BBEC1",nCleanedJetsPt30BTagged_bTagSF_jesUp_BBEC1);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_BBEC1_year",nCleanedJetsPt30BTagged_bTagSF_jesUp_BBEC1_year);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_EC2",nCleanedJetsPt30BTagged_bTagSF_jesUp_EC2);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_EC2_year",nCleanedJetsPt30BTagged_bTagSF_jesUp_EC2_year);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_FlavQCD",nCleanedJetsPt30BTagged_bTagSF_jesUp_FlavQCD);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_HF",nCleanedJetsPt30BTagged_bTagSF_jesUp_HF);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_HF_year",nCleanedJetsPt30BTagged_bTagSF_jesUp_HF_year);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_RelBal",nCleanedJetsPt30BTagged_bTagSF_jesUp_RelBal);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp_RelSample_year",nCleanedJetsPt30BTagged_bTagSF_jesUp_RelSample_year);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn",nCleanedJetsPt30BTagged_bTagSF_jesDn);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_Total",nCleanedJetsPt30BTagged_bTagSF_jesDn_Total);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_Abs",nCleanedJetsPt30BTagged_bTagSF_jesDn_Abs);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_Abs_year",nCleanedJetsPt30BTagged_bTagSF_jesDn_Abs_year);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_BBEC1",nCleanedJetsPt30BTagged_bTagSF_jesDn_BBEC1);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_BBEC1_year",nCleanedJetsPt30BTagged_bTagSF_jesDn_BBEC1_year);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_EC2",nCleanedJetsPt30BTagged_bTagSF_jesDn_EC2);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_EC2_year",nCleanedJetsPt30BTagged_bTagSF_jesDn_EC2_year);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_FlavQCD",nCleanedJetsPt30BTagged_bTagSF_jesDn_FlavQCD);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_HF",nCleanedJetsPt30BTagged_bTagSF_jesDn_HF);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_HF_year",nCleanedJetsPt30BTagged_bTagSF_jesDn_HF_year);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_RelBal",nCleanedJetsPt30BTagged_bTagSF_jesDn_RelBal);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn_RelSample_year",nCleanedJetsPt30BTagged_bTagSF_jesDn_RelSample_year);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jerUp",nCleanedJetsPt30BTagged_bTagSF_jerUp);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jerDn",nCleanedJetsPt30BTagged_bTagSF_jerDn);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSFUp",nCleanedJetsPt30BTagged_bTagSFUp);
   myTree->Book("nCleanedJetsPt30BTagged_bTagSFDn",nCleanedJetsPt30BTagged_bTagSFDn);

  if (addJets) {
    myTree->Book("JetPt",JetPt);
    myTree->Book("JetEta",JetEta);
    myTree->Book("JetPhi",JetPhi);
    myTree->Book("JetMass",JetMass);
    myTree->Book("JetEnergy",JetEnergy);
    myTree->Book("JetBTagger",JetBTagger);
    myTree->Book("JetIsBtagged",JetIsBtagged);
    myTree->Book("JetIsBtaggedWithSF",JetIsBtaggedWithSF);
    myTree->Book("JetIsBtaggedWithSFUp",JetIsBtaggedWithSFUp);
    myTree->Book("JetIsBtaggedWithSFDn",JetIsBtaggedWithSFDn);
    myTree->Book("JetQGLikelihood",JetQGLikelihood);
    myTree->Book("JetSigma",JetSigma);
    myTree->Book("JetHadronFlavour",JetHadronFlavour);
    myTree->Book("JetPartonFlavour",JetPartonFlavour);
    myTree->Book("JetPt_JESUp",JetJESUp);
    myTree->Book("JetPt_JESUp_Total",JetJESUp_Total);
    myTree->Book("JetPt_JESUp_Abs",JetJESUp_Abs);
    myTree->Book("JetPt_JESUp_Abs_year",JetJESUp_Abs_year);
    myTree->Book("JetPt_JESUp_BBEC1",JetJESUp_BBEC1);
    myTree->Book("JetPt_JESUp_BBEC1_year",JetJESUp_BBEC1_year);
    myTree->Book("JetPt_JESUp_EC2",JetJESUp_EC2);
    myTree->Book("JetPt_JESUp_EC2_year",JetJESUp_EC2_year);
    myTree->Book("JetPt_JESUp_FlavQCD",JetJESUp_FlavQCD);
    myTree->Book("JetPt_JESUp_HF",JetJESUp_HF);
    myTree->Book("JetPt_JESUp_HF_year",JetJESUp_HF_year);
    myTree->Book("JetPt_JESUp_RelBal",JetJESUp_RelBal);
    myTree->Book("JetPt_JESUp_RelSample_year",JetJESUp_RelSample_year);
    myTree->Book("JetPt_JESDown",JetJESDown);
    myTree->Book("JetPt_JESDown_Total",JetJESDown_Total);
    myTree->Book("JetPt_JESDown_Abs",JetJESDown_Abs);
    myTree->Book("JetPt_JESDown_Abs_year",JetJESDown_Abs_year);
    myTree->Book("JetPt_JESDown_BBEC1",JetJESDown_BBEC1);
    myTree->Book("JetPt_JESDown_BBEC1_year",JetJESDown_BBEC1_year);
    myTree->Book("JetPt_JESDown_EC2",JetJESDown_EC2);
    myTree->Book("JetPt_JESDown_EC2_year",JetJESDown_EC2_year);
    myTree->Book("JetPt_JESDown_FlavQCD",JetJESDown_FlavQCD);
    myTree->Book("JetPt_JESDown_HF",JetJESDown_HF);
    myTree->Book("JetPt_JESDown_HF_year",JetJESDown_HF_year);
    myTree->Book("JetPt_JESDown_RelBal",JetJESDown_RelBal);
    myTree->Book("JetPt_JESDown_RelSample_year",JetJESDown_RelSample_year);
    myTree->Book("JetPt_JERUp",JetPt_JERUp);
    myTree->Book("JetPt_JERDown",JetPt_JERDown);
    myTree->Book("JetRawPt",JetRawPt);
    myTree->Book("JetPtJEC_noJER",JetPtJEC_noJER);
    myTree->Book("JetPUID", JetPUID);
    myTree->Book("JetID", JetID);
    myTree->Book("JetPUValue", JetPUValue);
    myTree->Book("GenMET", GenMET);
    myTree->Book("GenMETPhi", GenMETPhi);
    myTree->Book("PFMET", metobj.extras.met);
    myTree->Book("PFMET_jesUp", metobj.extras.met_JECup);
    myTree->Book("PFMET_jesDn", metobj.extras.met_JECdn);
    myTree->Book("PFMETPhi", metobj.extras.phi);
    myTree->Book("PFMETPhi_jesUp", metobj.extras.phi_JECup);
    myTree->Book("PFMETPhi_jesDn", metobj.extras.phi_JECdn);
    myTree->Book("PFMET_corrected", metobj_corrected.extras.met);
    myTree->Book("PFMET_corrected_jesUp", metobj_corrected.extras.met_JECup);
    myTree->Book("PFMET_corrected_jesDn", metobj_corrected.extras.met_JECdn);
    myTree->Book("PFMET_corrected_jerUp", metobj_corrected.extras.met_JERup);
    myTree->Book("PFMET_corrected_jerDn", metobj_corrected.extras.met_JERdn);
    myTree->Book("PFMET_corrected_puUp", metobj_corrected.extras.met_PUup);
    myTree->Book("PFMET_corrected_puDn", metobj_corrected.extras.met_PUdn);
    myTree->Book("PFMET_corrected_metUp", metobj_corrected.extras.met_METup);
    myTree->Book("PFMET_corrected_metDn", metobj_corrected.extras.met_METdn);
    myTree->Book("PFMETPhi_corrected", metobj_corrected.extras.phi);
    myTree->Book("PFMETPhi_corrected_jesUp", metobj_corrected.extras.phi_JECup);
    myTree->Book("PFMETPhi_corrected_jesDn", metobj_corrected.extras.phi_JECdn);
    myTree->Book("PFMETPhi_corrected_jerUp", metobj_corrected.extras.phi_JERup);
    myTree->Book("PFMETPhi_corrected_jerDn", metobj_corrected.extras.phi_JERdn);
    myTree->Book("PFMETPhi_corrected_puUp", metobj_corrected.extras.phi_PUup);
    myTree->Book("PFMETPhi_corrected_puDn", metobj_corrected.extras.phi_PUdn);
    myTree->Book("PFMETPhi_corrected_metUp", metobj_corrected.extras.phi_METup);
    myTree->Book("PFMETPhi_corrected_metDn", metobj_corrected.extras.phi_METdn);
  }
   
  if(addPhotons){
    myTree->Book("PhotonPt",PhotonPt);
    myTree->Book("PhotonEta",PhotonEta);
    myTree->Book("PhotonPhi",PhotonPhi);
    myTree->Book("PhotonIsCutBasedLooseID",PhotonIsCutBasedLooseID);
  }
	
  if(addQGLInputs){
    myTree->Book("JetAxis2",JetAxis2);
    myTree->Book("JetMult",JetMult);
    myTree->Book("JetPtD",JetPtD);
  }
  if (isMC) {
    myTree->Book("genHEPMCweight",genHEPMCweight);
    myTree->Book("xsec",xsection);
    myTree->Book("dataMCWeight",dataMCWeight);
    myTree->Book("overallEventWeight",overallEventWeight);
    myTree->Book("L1prefiringWeight", L1prefiringWeight);
    myTree->Book("L1prefiringWeightUp", L1prefiringWeightUp);
    myTree->Book("L1prefiringWeightDn", L1prefiringWeightDn);
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZNtupleMaker);

