// system include files
#include <memory>
#include <cmath>
#include <string>

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
#include <DataFormats/Math/interface/LorentzVector.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <DataFormats/Common/interface/MergeableCounter.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include <SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h>
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>

#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>
#include <ZZAnalysis/AnalysisStep/interface/PUReweight.h>
#include <ZZAnalysis/AnalysisStep/interface/bitops.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>
#include <ZZAnalysis/AnalysisStep/interface/PUReweight.h>
#include "ZZ4lConfigHelper.h"
#include "HZZ4lNtupleFactory.h"

#include "Math/VectorUtil.h"
#include "TH1.h"

using namespace std;
using namespace edm;

namespace {
  //bool skipDataMCWeight = true; // skip computation of data/MC weight 

  //List of variables with default values
  Int_t RunNumber = 0;
  Long64_t EventNumber = 0;
  Int_t LumiNumber = 0;
  Short_t Nvtx = 0;
  Short_t NObsInt = 0;
  Float_t NTrueInt = 0;
  Float_t PUWeight = 0;
  Short_t trigWord = 0;
  Short_t nCleanedJetsPt30 = 0;
  Short_t Zsel = 0;
  Float_t ZMass = 0;
  Float_t ZPt = 0;
  Float_t ZEta = 0;
  Float_t ZPhi = 0;
  Short_t ZFlav = 0;
  std::vector<float> LepPt;
  std::vector<float> LepEta;
  std::vector<float> LepPhi;
  std::vector<short> LepLepId;
  std::vector<float> LepSIP;
  std::vector<float> LepTime;
  std::vector<bool> LepisID;
  std::vector<float> LepBDT;
  std::vector<char> LepMissingHit;
  std::vector<float> LepChargedHadIso;
  std::vector<float> LepNeutralHadIso;
  std::vector<float> LepPhotonIso;
  std::vector<float> LepCombRelIsoPF;
  std::vector<float> LepCombRelIsoPFPreFSR;
  std::vector<float> fsrPt;
  std::vector<float> fsrEta; 
  std::vector<float> fsrPhi;
  std::vector<float> fsrDR;
  std::vector<short> fsrLept;
  std::vector<short> fsrLeptID;
  Float_t genHEPMCweight = 0;
  Float_t xsection = 0;
  Float_t dataMCWeight = 0;
  Float_t overallEventWeight = 0;
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
  virtual void endJob();

  Float_t getAllWeight(Float_t LepPt, Float_t LepEta, Int_t LepID) const;

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

  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken;
  edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > candToken;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultToken;
  edm::EDGetTokenT<vector<reco::Vertex> > vtxToken;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetToken;

  edm::EDGetTokenT<edm::MergeableCounter> preSkimToken;
  
  PUReweight reweight;

  //counters
  Float_t Nevt_Gen;
  Float_t Nevt_Gen_lumiBlock;
  
  Float_t gen_sumPUWeight;
  Float_t gen_sumGenMCWeight;
  Float_t gen_sumWeights;

  string sampleName;
};

//
// constructors and destructor
//
ZNtupleMaker::ZNtupleMaker(const edm::ParameterSet& pset) :
  myHelper(pset),
  reweight()
{
  theCandLabel = pset.getUntrackedParameter<string>("CandCollection"); // Name of input Z collection
  theFileName = pset.getUntrackedParameter<string>("fileName"); 
  skipEmptyEvents = pset.getParameter<bool>("skipEmptyEvents"); // Do not store 
  sampleName = pset.getParameter<string>("sampleName");
  xsec = pset.getParameter<double>("xsec");
  year = pset.getParameter<int>("setup");

  consumesMany<std::vector< PileupSummaryInfo > >();
  genInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  candToken = consumes<edm::View<pat::CompositeCandidate> >(edm::InputTag(theCandLabel));
  triggerResultToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults"));
  vtxToken = consumes<vector<reco::Vertex> >(edm::InputTag("goodPrimaryVertices"));
  jetToken = consumes<edm::View<pat::Jet> >(edm::InputTag("cleanJets"));

  preSkimToken = consumes<edm::MergeableCounter,edm::InLumi>(edm::InputTag("preSkimCounter"));

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

  gen_sumPUWeight = 0.f;
  gen_sumGenMCWeight = 0.f;
  gen_sumWeights =0.f;
  
}

ZNtupleMaker::~ZNtupleMaker()
{
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

    PUWeight = reweight.weight(myHelper.sampleType(), myHelper.setup(), NTrueInt);

    edm::Handle<GenEventInfoProduct> gen;
    event.getByToken(genInfoToken, gen);
    GenEventInfoProduct  genInfo = *(gen.product());
    genHEPMCweight = genInfo.weight();

    // keep track of sum of weights
    gen_sumPUWeight    += PUWeight;
    gen_sumGenMCWeight += genHEPMCweight;
    gen_sumWeights     += PUWeight*genHEPMCweight;

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


  // General event information
  RunNumber = event.id().run();
  LumiNumber = event.luminosityBlock();
  EventNumber = event.id().event();
  xsection = xsec;

  // Primary vertices
  Handle<vector<reco::Vertex> > vertices;
  event.getByToken(vtxToken,vertices);
  Nvtx = vertices->size();

  // Jets
  Handle<edm::View<pat::Jet> > CleanedJets;
  event.getByToken(jetToken, CleanedJets);
  for(edm::View<pat::Jet>::const_iterator jet = CleanedJets->begin(); jet != CleanedJets->end(); ++jet)
    if(jet->pt()>30)
      nCleanedJetsPt30++;
  

  // Loop on candidates
  int nFilled = 0;
  for( edm::View<pat::CompositeCandidate>::const_iterator cand = cands->begin(); cand != cands->end(); ++cand) {
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
  LepLepId.clear();
  LepSIP.clear();
  LepTime.clear();
  LepisID.clear();
  LepBDT.clear();
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
    LepLepId.push_back( leptons[i]->pdgId() );
    LepSIP  .push_back( userdatahelpers::getUserFloat(leptons[i],"SIP") );
    LepTime .push_back( lepFlav==13 ? userdatahelpers::getUserFloat(leptons[i],"time") : 0. );
    LepisID .push_back( userdatahelpers::getUserFloat(leptons[i],"ID") );
    LepBDT  .push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"BDT") : 0. );
    LepMissingHit.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"missingHit") : 0 );
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

  //Compute the data/MC weight and overall event weight
  dataMCWeight = 1.;
  for(unsigned int i=0; i<leptons.size(); ++i){
    dataMCWeight *= getAllWeight(leptons[i]->pt(), leptons[i]->eta(), leptons[i]->pdgId());
  }
  overallEventWeight = PUWeight * genHEPMCweight * dataMCWeight;

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


Float_t ZNtupleMaker::getAllWeight(Float_t LepPt, Float_t LepEta, Int_t LepID) const
{
  // deactivated for the moment
  //if(skipDataMCWeight)
    return 1.;
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
  myTree->Book("nCleanedJetsPt30",nCleanedJetsPt30);
  myTree->Book("Zsel",Zsel);
  myTree->Book("ZMass",ZMass);
  myTree->Book("ZPt",ZPt);
  myTree->Book("ZEta",ZEta);
  myTree->Book("ZPhi",ZPhi);
  myTree->Book("ZFlav",ZFlav);
  myTree->Book("LepPt",LepPt);
  myTree->Book("LepEta",LepEta);
  myTree->Book("LepPhi",LepPhi);
  myTree->Book("LepLepId",LepLepId);
  myTree->Book("LepSIP",LepSIP);
  myTree->Book("LepTime",LepTime);
  myTree->Book("LepisID",LepisID);
  myTree->Book("LepBDT",LepBDT);
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
  myTree->Book("genHEPMCweight",genHEPMCweight);
  myTree->Book("xsec",xsection);
  myTree->Book("dataMCWeight",dataMCWeight);
  myTree->Book("overallEventWeight",overallEventWeight);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZNtupleMaker);
