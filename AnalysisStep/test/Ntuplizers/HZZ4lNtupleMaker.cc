// -*- C++ -*-
//
// Package:    HZZ4lNtupleMaker
// Class:      HZZ4lNtupleMaker
// 
/**\class HZZ4lNtupleMaker HZZ4lNtupleMaker.cc HZZ4lAnalysis/HZZ4lNtupleMaker/src/HZZ4lNtupleMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Stefano Casasso,,,
//         Created:  Tue Feb 28 14:33:03 CET 2012
// $Id: HZZ4lNtupleMaker.cc,v 1.101 2013/12/12 01:20:02 anderso Exp $
//
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
#include <ZZAnalysis/AnalysisStep/interface/VBFCandidateJetSelector.h>
#include <ZZAnalysis/AnalysisStep/interface/bitops.h>
#include <ZZAnalysis/AnalysisStep/interface/Fisher.h>
#include "ZZ4lConfigHelper.h"

// #include <ZZMatrixElement/MELA/interface/Mela.h>
// #include <ZZMatrixElement/MELA/src/computeAngles.h>

#include "HZZ4lNtupleFactory.h"

#include "TLorentzVector.h"

#include <string>

namespace {
  bool writePhotons = false;  // Write photons in the tree. Note: must be set also in HZZ4lNtupleFactory.cc
  bool writeJets = true;     // Write jets in the tree. FIXME: make this configurable
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
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void FillCandidate(const pat::CompositeCandidate& higgs, bool evtPass, const edm::Event&, const Int_t CRflag);
  virtual void FillPhoton(const pat::Photon& photon);
  virtual void FillJet(const pat::Jet& jet);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  // ----------member data ---------------------------
  ZZ4lConfigHelper myHelper;
  int theChannel;
  std::string theCandLabel;
  TString theFileName;

  HZZ4lNtupleFactory *myTree;
  TH1F *hCounter;

  Bool_t isMC;

  bool applyTrigger;    // Only events passing trigger
  bool applySkim;       //  "     "      "     skim
  bool skipEmptyEvents; // Skip events whith no candidate in the collection
  bool writeBestCandOnly;

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
  Float_t gen_ZZ4mu_MassAcceptance;
  Float_t gen_ZZ4mu_MassPtAcceptance;
  Float_t gen_ZZ4e_EtaAcceptance;
  Float_t gen_ZZ4e_LeptonAcceptance;
  Float_t gen_ZZ4e_MassAcceptance;
  Float_t gen_ZZ4e_MassPtAcceptance;
  Float_t gen_ZZ2mu2e_EtaAcceptance; 
  Float_t gen_ZZ2mu2e_LeptonAcceptance; 
  Float_t gen_ZZ2mu2e_MassAcceptance; 
  Float_t gen_ZZ2mu2e_MassPtAcceptance;
  Float_t gen_WH4mu;
  Float_t gen_WH4e;
  Float_t gen_WH2mu2e;
  Float_t gen_ZH4mu;
  Float_t gen_ZH4e;
  Float_t gen_ZH2mu2e;
  Float_t gen_ttH4mu;
  Float_t gen_ttH4e;
  Float_t gen_ttH2mu2e;
  Float_t gen_4mu_m180;
  Float_t gen_4e_m180;
  Float_t gen_2e2mu_m180;  
  Float_t gen_BUGGY;
  Float_t gen_Unknown;
  
  Float_t gen_sumPUWeight;
  Float_t gen_sumGenMCWeight;
  Float_t gen_sumWeights;

  string sampleName;
  //  Mela mela;

};

//
// constructors and destructor
//
HZZ4lNtupleMaker::HZZ4lNtupleMaker(const edm::ParameterSet& pset) :
  myHelper(pset),
  reweight()
  //  mela((pset.getParameter<int>("setup")==2011)?7:8,pset.getParameter<double>("superMelaMass")) //FIXME: need to handle cases where setup>2012
{
  theCandLabel = pset.getUntrackedParameter<string>("CandCollection");
  theChannel = myHelper.channel();
  theFileName = pset.getUntrackedParameter<string>("fileName");
  skipEmptyEvents = pset.getParameter<bool>("skipEmptyEvents");
  writeBestCandOnly = pset.getParameter<bool>("onlyBestCandidate");
  sampleName = pset.getParameter<string>("sampleName");
  
  //  mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);

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
  gen_ZZ4mu_MassAcceptance = 0;
  gen_ZZ4mu_MassPtAcceptance = 0;
  gen_ZZ4e_EtaAcceptance = 0;
  gen_ZZ4e_LeptonAcceptance = 0;
  gen_ZZ4e_MassAcceptance = 0;
  gen_ZZ4e_MassPtAcceptance = 0;
  gen_ZZ2mu2e_EtaAcceptance = 0;
  gen_ZZ2mu2e_LeptonAcceptance = 0;
  gen_ZZ2mu2e_MassAcceptance = 0;
  gen_ZZ2mu2e_MassPtAcceptance = 0;
  gen_WH4mu = 0;
  gen_WH4e = 0;
  gen_WH2mu2e = 0;
  gen_ZH4mu = 0;
  gen_ZH4e = 0;
  gen_ZH2mu2e = 0;
  gen_ttH4mu = 0;
  gen_ttH4e = 0;
  gen_ttH2mu2e = 0;
  gen_4mu_m180 = 0;
  gen_4e_m180 = 0;
  gen_2e2mu_m180 = 0;
  gen_BUGGY = 0;
  gen_Unknown = 0;

  gen_sumPUWeight = 0.f;
  gen_sumGenMCWeight = 0.f;
  gen_sumWeights =0.f;
}

HZZ4lNtupleMaker::~HZZ4lNtupleMaker()
{
}


//
// member functions
//

// ------------ method called for each event  ------------
void HZZ4lNtupleMaker::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  Handle<vector<reco::Vertex> >  vertexs;
  event.getByLabel("goodPrimaryVertices",vertexs);
  
  //----------------------------------------------------------------------
  // Analyze MC history. THIS HAS TO BE DONE BEFORE ANY RETURN STATEMENT
  // (eg skim or trigger), in order to update the gen counters correctly!!!
  int nObsInt  = -1;
  float nTrueInt = -1.;
  Float_t weight2 = 1.;
  Int_t genFinalState = -1;
  Int_t genProcessId = -1;
  Float_t genHEPMCweight = 1.;
  Short_t genExtInfo = -1;

  bool gen_ZZ4lInEtaAcceptance = false;   // All 4 gen leptons in eta acceptance
  bool gen_ZZ4lInEtaPtAcceptance = false; // All 4 gen leptons in eta,pT acceptance
  bool gen_m4l_180 = false;               // gen_m4l > 180


  const reco::Candidate * genH = 0;
  //  std::vector<const reco::Candidate *> genZs;
  std::vector<const reco::Candidate *> genZLeps;
  std::vector<const reco::Candidate *> genAssocLeps;
  if (isMC) {

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

    // keep value of sum of weights
    gen_sumPUWeight    += weight2;
    gen_sumGenMCWeight += genHEPMCweight;
    gen_sumWeights     += weight2*genHEPMCweight;

    bool gen_ZZInAcceptance = false; // Unused; old ZZ phase space

    mch.genAcceptance(gen_ZZInAcceptance, gen_ZZ4lInEtaAcceptance, gen_ZZ4lInEtaPtAcceptance, gen_m4l_180);

    ++Nevt_Gen_lumiBlock;
    if (genFinalState == EEEE) {
      ++gen_ZZ4e;
      if (gen_ZZ4lInEtaAcceptance) ++gen_ZZ4e_EtaAcceptance;
      if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ4e_LeptonAcceptance;
      if (gen_ZZInAcceptance) {
	++gen_ZZ4e_MassAcceptance;
	if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ4e_MassPtAcceptance;	
      }
      if (gen_m4l_180) ++gen_4e_m180;
      if (genProcessId==26) ++gen_WH4e;
      if (genProcessId==24) ++gen_ZH4e;
      if (genProcessId==121 || genProcessId==122) ++gen_ttH4e;
    } else if (genFinalState == MMMM) {
      ++gen_ZZ4mu;
      if (gen_ZZ4lInEtaAcceptance) ++gen_ZZ4mu_EtaAcceptance;
      if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ4mu_LeptonAcceptance;
      if (gen_ZZInAcceptance) {
	++gen_ZZ4mu_MassAcceptance;
	if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ4mu_MassPtAcceptance;
      }
      if (gen_m4l_180) ++gen_4mu_m180;
      if (genProcessId==26) ++gen_WH4mu;
      if (genProcessId==24) ++gen_ZH4mu;
      if (genProcessId==121 || genProcessId==122) ++gen_ttH4mu;
    } else if (genFinalState == EEMM) {
      ++gen_ZZ2mu2e;
      if (gen_ZZ4lInEtaAcceptance) ++gen_ZZ2mu2e_EtaAcceptance;
      if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ2mu2e_LeptonAcceptance;
      if (gen_ZZInAcceptance) {
	++gen_ZZ2mu2e_MassAcceptance;
	if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ2mu2e_MassPtAcceptance;
      }
      if (gen_m4l_180) ++gen_2e2mu_m180;
      if (genProcessId==26) ++gen_WH2mu2e;
      if (genProcessId==24) ++gen_ZH2mu2e;
      if (genProcessId==121 || genProcessId==122) ++gen_ttH2mu2e;
    } else if (genFinalState == LLTT){
      ++gen_ZZ2l2tau;
    } else if (genFinalState == BUGGY){ // handle H->ddbar gen bug!!!
      ++gen_BUGGY;
      return; // BUGGY events are skipped
    } else {
      ++gen_Unknown;
    }
    

    //Information on generated candidates, will be used later
    genH = mch.genH();
    //    genZs = mch.genZs(); // These are the pdgID=23 particles in the MC history, if present.
    //    genZLeps = mch.genZLeps();
    genZLeps     = mch.sortedGenZZLeps();
    genAssocLeps = mch.genAssociatedLeps();
  }
  //----------------------------------------------------------------------

  //Get candidate collection
  edm::Handle<edm::View<pat::CompositeCandidate> > candHandle;
  event.getByLabel(theCandLabel, candHandle);
  const edm::View<pat::CompositeCandidate>* cands = candHandle.product();

  myTree->InitializeVariables();

  if (skipEmptyEvents) {
    if (cands->size() == 0) return; // Skip events with no candidate
  } else { // Storing events with no candidate is useful to keep information on gen particles on signal. For these events the rest of the tree is irrelevant.
    if (cands->size() == 0 && genFinalState!=theChannel) return; // Skip empty events of the "wrong" gen final state in SR of signal samples
  }
  

  // For Z+L CRs, we want only events with exactly 1 Z+l candidate.
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


  // This is very important in order not to carry information from the previous event. 
  // NOTHING IN THE NTUPLE MUST BE FILLED BEFORE THIS LINE
  // because the cleaning of the ntuple objects is done at the end, in the FillEvent() method.

  //Counter to find the best candidate
  Int_t NbestCand = -1;
  Int_t CandCounter = 0;
  bool wasFilled = false;

  for( edm::View<pat::CompositeCandidate>::const_iterator cand = cands->begin(); cand != cands->end(); ++cand) {
    Int_t CRFLAG=0;
    bool candIsBest = cand->userFloat("isBestCand");

    int candChannel = cand->userFloat("candChannel");
    // In the signal region, select only candidates of the final state we are looking at.
    if ((theChannel==MMMM && candChannel != 28561) ||
	(theChannel==EEEE && candChannel != 14641) ||
	(theChannel==EEMM && candChannel != 20449)) continue;    
    
    if (theChannel==ZLL) {
      // AA CRs
      if(cand->userFloat("isBestCRZLLss")&&cand->userFloat("CRZLLss"))
	set_bit(CRFLAG,CRZLLss);      

      // Older Z2 ID/noSIP CRs
      if(cand->userFloat("isBestCRZLL")&&cand->userFloat("CRZLL"))
	set_bit(CRFLAG,CRZLL);
    }
    
    //    if(theChannel==ZL){} Nothing special in this case
 
    if (writeBestCandOnly && !(candIsBest||CRFLAG)) continue; // Skip events other than the best cand (or CR candidates in the CR)
    
    //For the SR, also fold information about acceptance in CRflag 
    if (isMC && (theChannel==EEEE||theChannel==MMMM||theChannel==EEMM)) {
      if (gen_ZZ4lInEtaAcceptance)   set_bit(CRFLAG,28);
      if (gen_ZZ4lInEtaPtAcceptance) set_bit(CRFLAG,29);
      if (gen_m4l_180)               set_bit(CRFLAG,30);
    }

    FillCandidate(*cand, evtPassTrigger&&evtPassSkim, event,CRFLAG);
    wasFilled = true;
    if(candIsBest) NbestCand = CandCounter;
    CandCounter++;
  }

  // Events with no candidate have already been skipped at the beginning, but in case of CRs there could be no cands with CRflag!=0, so none was filled. 
  if (skipEmptyEvents && (!wasFilled) && theChannel==ZLL) return;

  if (isMC) {

    if(genH != 0){ // H explicit
      myTree->FillHGenInfo(genH->p4());
    }
    else if(genZLeps.size()==4){ // take the ZZ(4l) system 
      myTree->FillHGenInfo((genZLeps.at(0)->p4()+genZLeps.at(1)->p4()+genZLeps.at(2)->p4()+genZLeps.at(3)->p4()));
    }

    //if (genFinalState!=BUGGY && genFinalState!=NONE) {
    if (genFinalState!=BUGGY) { // removed genFinalState!=NONE requirement, so as to include events where the ZZ system doesn't decay to 4 leptons, e.g. ZH with Z->ll, H->ZZ->llqq ... 

//       if (genZs.size()==2){
// 	myTree->FillZGenInfo(genZs.at(0)->p4(), genZs.at(1)->p4());
//       }
      if (genZLeps.size()==4) {
	
	myTree->FillZGenInfo(genZLeps.at(0)->p4()+genZLeps.at(1)->p4(),
			     genZLeps.at(2)->p4()+genZLeps.at(3)->p4());

	myTree->FillLepGenInfo(genZLeps.at(0)->pdgId(), genZLeps.at(1)->pdgId(), genZLeps.at(2)->pdgId(), genZLeps.at(3)->pdgId(),
			       genZLeps.at(0)->p4(), genZLeps.at(1)->p4(), genZLeps.at(2)->p4(), genZLeps.at(3)->p4());

	//FIXME
// 	math::XYZTLorentzVector p11 = genZLeps.at(0)->p4();
// 	math::XYZTLorentzVector p12 = genZLeps.at(1)->p4();
// 	math::XYZTLorentzVector p21 = genZLeps.at(2)->p4();
// 	math::XYZTLorentzVector p22 = genZLeps.at(3)->p4();
//      TLorentzVector pL11(p11.x(),p11.y(),p11.z(),p11.t());
//      TLorentzVector pL12(p12.x(),p12.y(),p12.z(),p12.t());
//      TLorentzVector pL21(p21.x(),p21.y(),p21.z(),p21.t());
//      TLorentzVector pL22(p22.x(),p22.y(),p22.z(),p22.t());
// 	int id11 = genZLeps.at(0)->pdgId();
// 	int id12 = genZLeps.at(1)->pdgId();
// 	int id21 = genZLeps.at(2)->pdgId();
// 	int id22 = genZLeps.at(3)->pdgId();
// 	float gencostheta1=0, gencostheta2=0, genphi=0, gencosthetastar=0, genphistar1=0;
//	mela::computeAngles(pL11,id11,pL12,id12,pL21,id21,pL22,id22,gencosthetastar,gencostheta1,gencostheta2,genphi,genphistar1);
	
      }

      if (genZLeps.size()==3) {
	myTree->FillLepGenInfo(genZLeps.at(0)->pdgId(), genZLeps.at(1)->pdgId(), genZLeps.at(2)->pdgId(), 0,
			       genZLeps.at(0)->p4(), genZLeps.at(1)->p4(), genZLeps.at(2)->p4(), *(new math::XYZTLorentzVector));
      }
      if (genZLeps.size()==2) {
	myTree->FillLepGenInfo(genZLeps.at(0)->pdgId(), genZLeps.at(1)->pdgId(), 0, 0,
			       genZLeps.at(0)->p4(), genZLeps.at(1)->p4(), *(new math::XYZTLorentzVector), *(new math::XYZTLorentzVector));
      }

      if (genAssocLeps.size()==1 || genAssocLeps.size()==2) {
	myTree->FillAssocLepGenInfo(genAssocLeps);
      }

    }

  }


  

  // cmgPhotons (store them only for events with at least 1 candidate)
  if (writePhotons && cands->size()!=0) {
    edm::Handle<vector<pat::Photon> > photons;
    if (event.getByLabel("cmgPhotonSel",photons)) {
      //cout << "Photons: " << photons->size() <<endl;
      for( vector<pat::Photon>::const_iterator ph = photons->begin(); ph != photons->end(); ++ph) {
	//cout << "photon: pt= " << ph->pt() << endl;
	if(ph->pt() > 2. && fabs(ph->eta()) < 2.8) FillPhoton(*ph);
      }
    }
  }


  // cmgJets (store them only for events with at least 1 candidate)
//   if (writeJets && cands->size()!=0) {
//     edm::Handle<vector<pat::Jet> > jets;
//     if (event.getByLabel("cmgPFJetSel",jets)) {
//       //cout << "Jets: " << jets->size() <<endl;
//       for( vector<pat::Jet>::const_iterator j = jets->begin(); j != jets->end(); ++j) {
// 	//cout << "jet: pt= " << j->pt() << endl;
// 	bool looseJetID =(j->getSelection("cuts_looseJetId")>0);
// 	if(looseJetID && j->passPuJetId("full", PileupJetIdentifier::kLoose )) FillJet(*j);
//       }
//     }
//   }

  // Jet collection (preselected with pT>10)
  Handle<edm::View<pat::Jet> > pfjetscoll;
  event.getByLabel("slimmedJets", pfjetscoll);

  // lepton collection
  Handle<View<reco::Candidate> > softleptoncoll;
  event.getByLabel("softLeptons", softleptoncoll);
  vector<const reco::Candidate*> goodisoleptons;
  for( View<reco::Candidate>::const_iterator lep = softleptoncoll->begin(); lep != softleptoncoll->end(); ++ lep ){ 
    if((bool)userdatahelpers::getUserFloat(&*lep,"isGood") && (bool)userdatahelpers::getUserFloat(&*lep,"isIsoFSRUncorr")){
      goodisoleptons.push_back(&*lep);
    }
  }

  std::vector<const pat::Jet*> cleanedJets;
  VBFCandidateJetSelector myVBFCandidateJetSelector;
  cleanedJets = myVBFCandidateJetSelector.cleanJets(goodisoleptons,pfjetscoll,myHelper.setup());
  vector<const pat::Jet*> cleanedJetsPt30;
  int nCleanedJetsPt30BTagged = 0;
  for (unsigned int i=0; i<cleanedJets.size(); ++i){
    const pat::Jet& myjet = *(cleanedJets.at(i));  
    if (myjet.pt()>30) {
      cleanedJetsPt30.push_back(&myjet);
      if(myjet.bDiscriminator("combinedSecondaryVertexBJetTags")>0.679) nCleanedJetsPt30BTagged++; // CSV Medium WP
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
    }

    // Note that jets variables are filled for jets above 20 GeV, to allow JES studies.
    // detajj, Mjj and ZZFisher are filled only for true dijet events (jets above 30 GeV)
    if(cleanedJets.size()>1 && theChannel!=ZL){ 
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
      myTree->FillDiJetInfo(Mjj,MjjPlus,MjjMinus,detajj,Fisher);
    }

  }
  

  //MET info
  //   Handle<reco::PFMETCollection> pfmetcoll;
  //   event.getByLabel("patMETs", pfmetcoll);
  //   float pfmet = -1;
  //   if(pfmetcoll.isValid()){
  //     const reco::PFMETCollection *pfmetcol = pfmetcoll.product();
  //     const reco::PFMET *pfmetObj = &(pfmetcol->front());
  //     pfmet = pfmetObj->pt();
  //     //cout << pfmet << endl;
  //   }

  Handle<vector<reco::MET> > pfmetcoll;
  event.getByLabel("slimmedMETs", pfmetcoll);
  float pfmet = -1;
  if(pfmetcoll.isValid()){
    pfmet = pfmetcoll->front().pt();
  }


  //Save general event info in the tree. This must be done after the loop on the candidates so that we know the best candidate position in the list
  myTree->FillEventInfo(event.id().run(), event.id().event(), event.luminosityBlock(), NbestCand, vertexs->size(), nObsInt, nTrueInt, weight2, pfmet, pfjetscoll->size(), cleanedJets.size(), cleanedJetsPt30.size(), nCleanedJetsPt30BTagged, genFinalState, genProcessId, genHEPMCweight, trigWord, genExtInfo);

  myTree->FillEvent();

  return;
}

void HZZ4lNtupleMaker::FillPhoton(const pat::Photon& photon)
{
  const Float_t photPt  = photon.pt();
  const Float_t photEta = photon.eta();
  const Float_t photPhi = photon.phi();

  myTree->FillPhotonInfo(photPt, photEta, photPhi);

  return;
}


void HZZ4lNtupleMaker::FillJet(const pat::Jet& jet)
{
  const Float_t jetPt  = jet.pt();
  const Float_t jetEta = jet.eta();
  const Float_t jetPhi = jet.phi();
  const Float_t jetMass = jet.p4().M();
  const Float_t jetBTag = jet.bDiscriminator("combinedSecondaryVertexBJetTags");
  const Float_t jesUnc = 0.;//jet.uncOnFourVectorScale();

  myTree->FillJetInfo(jetPt, jetEta, jetPhi, jetMass, jetBTag, jesUnc );

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

  Int_t isSignal = -1;
  Int_t isRightPair = -1;
  if(isMC){
    isSignal = cand.userFloat("MC_isRight");
    isRightPair = cand.userFloat("MC_isRightPair");
  }

  const Float_t mZa = cand.userFloat("mZa");
  const Float_t mZb = cand.userFloat("mZb");
  const Float_t mLL4 = cand.userFloat("mLL4");
  const Float_t mLL6 = cand.userFloat("mLL6");
  const Float_t SIP4 = cand.userFloat("SIP4");
  const Float_t iso34 = cand.userFloat("iso34");

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

  const Float_t Z1Mass = Z1->mass();
  const Float_t Z1Pt =   Z1->pt();
  const Float_t Z1MassRefit = cand.userFloat("mZ1Ref");

  const Float_t Z2Mass = Z2->mass();
  const Float_t Z2Pt =   Z2->pt();

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

  myTree->FillProbability(//p0plus_melaNorm,
			  //p0plus_mela,
			  //p0minus_mela,
			  //p0hplus_mela, // 0h+, analytic distribution
			  p0plus_VAJHU,
			  p0minus_VAJHU,
			  p0plus_VAMCFM,
			  p0hplus_VAJHU, // 0h+ (high dimensional operator), vector algebra, JHUgen
			  //p1_mela,
			  //p1_prodIndep_mela,
			  //p1plus_mela, // 1+, analytic distribution 
			  //p1plus_prodIndep_mela, // 1+, analytic distribution 
			  p1_VAJHU,
			  p1_prodIndep_VAJHU,
			  p1plus_VAJHU, // 1+ (axial vector), vector algebra, JHUgen,
			  p1plus_prodIndep_VAJHU, // 1+ (axial vector), vector algebra, JHUgen,
			  //p2_mela ,
			  //p2_prodIndep_mela ,
			  //p2qqb_mela, // graviton produced by qqbar vector algebra, analytical,
			  //p2hplus_mela, // graviton produced by qqbar vector algebra, analytical,
			  //p2hminus_mela, // graviton produced by qqbar vector algebra, analytical,
			  //p2bplus_mela, // graviton produced by qqbar vector algebra, analytical,
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
			  //bkg_mela,
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
			  p0gsgs_PS_VAJHU
			  //bkg_VAMCFMNorm,
			  //p0_pt,
			  //p0_y,
			  //bkg_pt,
			  //bkg_y
			  );
			  
  myTree->FillSuperMela( p0plus_m4l,  // signal m4l probability as in datacards
			 bkg_m4l,   // backgroun m4l probability as in datacards
			 p0plus_m4l_ScaleUp,  // signal m4l probability for systematics
			 bkg_m4l_ScaleUp,     // backgroun m4l probability for systematics
			 p0plus_m4l_ScaleDown,  // signal m4l probability for systematics
			 bkg_m4l_ScaleDown,     // backgroun m4l probability for systematics
			 p0plus_m4l_ResUp,  // signal m4l probability for systematics
			 bkg_m4l_ResUp,     // backgroun m4l probability for systematics
			 p0plus_m4l_ResDown,  // signal m4l probability for systematics
			 bkg_m4l_ResDown);     // backgroun m4l probability for systematics
  
  myTree->FillHAdditionalInfo(mZa, mZb, mLL4, mLL6, SIP4, iso34);

  myTree->FillZInfo(Z1Mass, Z1Pt, Z1MassRefit);
  myTree->FillZInfo(Z2Mass, Z2Pt, Z2Mass);

  //Fill the angular variables
  const Float_t costheta1 = cand.userFloat("costheta1");
  const Float_t costheta2 = cand.userFloat("costheta2");
  const Float_t phi       = cand.userFloat("phi");
  const Float_t costhetastar = cand.userFloat("costhetastar");
  const Float_t phistar1      = cand.userFloat("phistar1");
  const Float_t phistar2      = cand.userFloat("phistar2");
  const Float_t xi            = cand.userFloat("xi");
  const Float_t xistar        = cand.userFloat("xistar");
//   const Float_t phi1       = cand.userFloat("helphiZl1"); // unused
//   const Float_t phi2       = cand.userFloat("helphiZl2");
  myTree->FillAngularInfo(costhetastar, phi, costheta1, costheta2, phistar1, phistar2,xi,xistar);

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
    myTree->FillLepInfo(leptons[i]->pt(),
			leptons[i]->eta(),
			leptons[i]->phi(),
			leptons[i]->pdgId(),
			SIP[i],
			isID[i],
			userdatahelpers::getUserFloat(leptons[i],"BDT"),
			userdatahelpers::getUserFloat(leptons[i],"MCParentCode"), // This is dummy as of now.
			userdatahelpers::getUserFloat(leptons[i],"missingHit"));

    //Isolation variables
    myTree->FillLepIsolInfo(PFChargedHadIso[i],
			    PFNeutralHadIso[i],
			    PFPhotonIso[i], 
			    combRelIsoPF[i]);
  }
  
  //convention: 0 -> 4mu   1 -> 4e   2 -> 2mu2e
  myTree->FillHInfo(ZZMass, ZZMassErr, ZZMassErrCorr, ZZMassPreFSR, ZZMassRefit, Chi2KinFit, ZZMassCFit, Chi2CFit,  sel, ZZPt, ZZEta, ZZPhi,
		    isSignal, isRightPair, CRflag);

  //Fill the info on categorization
  const Int_t nExtraLep = cand.userFloat("nExtraLep");
  const Int_t nExtraZ = cand.userFloat("nExtraZ");
  myTree->FillCategorizationInfo(nExtraLep, nExtraZ);

  //Fill the info on the extra leptons
  myTree->FillExtraLepInfo( 1, cand.hasUserCand("ExtraLep1"), (cand.hasUserCand("ExtraLep1") ? cand.userCand("ExtraLep1") : *(new reco::CandidatePtr)) );
  myTree->FillExtraLepInfo( 2, cand.hasUserCand("ExtraLep2"), (cand.hasUserCand("ExtraLep2") ? cand.userCand("ExtraLep2") : *(new reco::CandidatePtr)) );
  myTree->FillExtraLepInfo( 3, cand.hasUserCand("ExtraLep3"), (cand.hasUserCand("ExtraLep3") ? cand.userCand("ExtraLep3") : *(new reco::CandidatePtr)) );

}



// ------------ method called once each job just before starting event loop  ------------
void HZZ4lNtupleMaker::beginJob()
{
  edm::Service<TFileService> fs;
  myTree = new HZZ4lNtupleFactory( fs->make<TTree>(theFileName,"Event Summary"));
  hCounter = fs->make<TH1F>("Counters", "Counters", 45, 0., 45.);
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
  hCounter->SetBinContent(8 ,gen_ZZ4mu_MassAcceptance);
  hCounter->SetBinContent(9 ,gen_ZZ4mu_MassPtAcceptance);
  hCounter->SetBinContent(10,gen_ZZ4e_EtaAcceptance);
  hCounter->SetBinContent(11,gen_ZZ4e_LeptonAcceptance);
  hCounter->SetBinContent(12,gen_ZZ4e_MassAcceptance);
  hCounter->SetBinContent(13,gen_ZZ4e_MassPtAcceptance);
  hCounter->SetBinContent(14,gen_ZZ2mu2e_EtaAcceptance);
  hCounter->SetBinContent(15,gen_ZZ2mu2e_LeptonAcceptance);
  hCounter->SetBinContent(16,gen_ZZ2mu2e_MassAcceptance);
  hCounter->SetBinContent(17,gen_ZZ2mu2e_MassPtAcceptance);
  hCounter->SetBinContent(19,gen_BUGGY);
  hCounter->SetBinContent(20,gen_Unknown);
  hCounter->SetBinContent(21,gen_WH4mu);
  hCounter->SetBinContent(22,gen_WH4e);
  hCounter->SetBinContent(23,gen_WH2mu2e);
  hCounter->SetBinContent(24,gen_ZH4mu);
  hCounter->SetBinContent(25,gen_ZH4e);
  hCounter->SetBinContent(26,gen_ZH2mu2e);
  hCounter->SetBinContent(27,gen_ttH4mu);
  hCounter->SetBinContent(28,gen_ttH4e);
  hCounter->SetBinContent(29,gen_ttH2mu2e);
  hCounter->SetBinContent(30,gen_4mu_m180);
  hCounter->SetBinContent(31,gen_4e_m180);
  hCounter->SetBinContent(32,gen_2e2mu_m180);

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
  hCounter->GetXaxis()->SetBinLabel(8 ,"gen_ZZ4mu_MassAcceptance");
  hCounter->GetXaxis()->SetBinLabel(9 ,"gen_ZZ4mu_MassPtAcceptance");
  hCounter->GetXaxis()->SetBinLabel(10,"gen_ZZ4e_EtaAcceptance");
  hCounter->GetXaxis()->SetBinLabel(11,"gen_ZZ4e_LeptonAcceptance");
  hCounter->GetXaxis()->SetBinLabel(12,"gen_ZZ4e_MassAcceptance");
  hCounter->GetXaxis()->SetBinLabel(13,"gen_ZZ4e_MassPtAcceptance");
  hCounter->GetXaxis()->SetBinLabel(14,"gen_ZZ2mu2e_EtaAcceptance");
  hCounter->GetXaxis()->SetBinLabel(15,"gen_ZZ2mu2e_LeptonAcceptance");
  hCounter->GetXaxis()->SetBinLabel(16,"gen_ZZ2mu2e_MassAcceptance");
  hCounter->GetXaxis()->SetBinLabel(17,"gen_ZZ2mu2e_MassPtAcceptance");
  hCounter->GetXaxis()->SetBinLabel(19,"gen_BUGGY");
  hCounter->GetXaxis()->SetBinLabel(20,"gen_Unknown");
  hCounter->GetXaxis()->SetBinLabel(21,"gen_WH4mu");
  hCounter->GetXaxis()->SetBinLabel(22,"gen_WH4e");
  hCounter->GetXaxis()->SetBinLabel(23,"gen_WH2mu2e");
  hCounter->GetXaxis()->SetBinLabel(24,"gen_ZH4mu");
  hCounter->GetXaxis()->SetBinLabel(25,"gen_ZH4e");
  hCounter->GetXaxis()->SetBinLabel(26,"gen_ZH2mu2e");
  hCounter->GetXaxis()->SetBinLabel(27,"gen_ttH4mu");
  hCounter->GetXaxis()->SetBinLabel(28,"gen_ttH4e");
  hCounter->GetXaxis()->SetBinLabel(29,"gen_ttH2mu2e");
  hCounter->GetXaxis()->SetBinLabel(30,"gen_4mu_m180");
  hCounter->GetXaxis()->SetBinLabel(31,"gen_4e_m180");
  hCounter->GetXaxis()->SetBinLabel(32,"gen_2e2mu_m180");

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

//define this as a plug-in
DEFINE_FWK_MODULE(HZZ4lNtupleMaker);
