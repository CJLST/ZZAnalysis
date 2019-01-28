/** \class ZZ4lAnalyzerCR  
 *
 *  For the time being: retrieving all relevant information from CR collections and filling the 2011 "flagship" plots
 *
 *  $Date: 2013/07/01 17:29:01 $
 *  $Revision: 1.14 $
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
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include "FWCore/ServiceRegistry/interface/Service.h"


#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <DataFormats/Common/interface/MergeableCounter.h>

#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>
#include <ZZAnalysis/AnalysisStep/interface/PileUpWeight.h>
#include "ZZ4lConfigHelper.h"


#include <iostream>
#include <iterator>
#include <string>

#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"

#define HCMSSW
#include "ZZAnalysis/AnalysisStep/interface/Histograms.h"
#include <algorithm>

using namespace std;
using namespace edm;
using namespace reco;


class ZZ4lAnalyzerCR: public edm::EDAnalyzer {
public:

  explicit ZZ4lAnalyzerCR(const edm::ParameterSet& pset);
  ~ZZ4lAnalyzerCR();
  
  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
  virtual void beginJob();
  virtual void endJob();  
  
  void printParticle(const reco::Candidate* c=0, int idx=0, int pdgId=0, 
		     float looseIso=-1, float iso=-1, float SIP=-1);
  
private:
  ZZ4lConfigHelper myHelper;
  Channel theChannel;
  bool isMC;
  PileUpWeight* pileUpReweight;

  double weight;

  edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > candToken;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultToken;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > PupInfoToken;

  //histograms
  TH1F* nEventComplete;
  HCand* hCandCR;
  HCand* hCandCR_w;

};

// Constructor
ZZ4lAnalyzerCR::ZZ4lAnalyzerCR(const ParameterSet& pset) :
  myHelper(pset),
  pileUpReweight(nullptr),
  weight(1.),
  candToken(consumes<edm::View<pat::CompositeCandidate> >(pset.getUntrackedParameter<edm::InputTag>("candCollection")))
{
  isMC = myHelper.isMC();
  theChannel = myHelper.channel();
//   if (theChannel!=18) {
//     cout << "ERROR: ZZ4lAnalyzerCR: channel "<< theChannel << " is not valid" <<endl;
//     abort();
//   }

  triggerResultToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults"));
  PupInfoToken = consumes<std::vector<PileupSummaryInfo> >(edm::InputTag("addPileupInfo"));
  if (isMC) pileUpReweight = new PileUpWeight(myHelper.sampleType(), myHelper.setup());
}

ZZ4lAnalyzerCR::~ZZ4lAnalyzerCR()
{
  delete pileUpReweight;
}


void ZZ4lAnalyzerCR::beginJob(){
  // Book histograms
  edm::Service<TFileService> fileService;

  // Distributions
  TH1F::SetDefaultSumw2(kTRUE);
  hCandCR = new HCand("hCandCR");
  hCandCR_w = new HCand("hCandCR_w");
   
}

void ZZ4lAnalyzerCR::endJob(){
   
}


void ZZ4lAnalyzerCR::analyze(const Event & event, const EventSetup& eventSetup){
    

  // Trigger results
  Handle<edm::TriggerResults> triggerResults;
  event.getByToken(triggerResultToken, triggerResults);

  // Apply MC filter (skip event)
  if (isMC && !(myHelper.passMCFilter(event,triggerResults))) return;
 
  // Skim
  bool evtPassSkim = myHelper.passSkim(event,triggerResults);

  // Trigger requests
  bool evtPassTrigger = myHelper.passTrigger(event,triggerResults); 
  if (!(evtPassTrigger && evtPassSkim))return;

 
  // PU reweight
  float PUweight = 1.;
  if (isMC) {
    float nTrueInt = -1.;
    Handle<std::vector<PileupSummaryInfo> > PupInfo;
    event.getByToken(PupInfoToken, PupInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      if(PVI->getBunchCrossing() == 0) { 
	//	nObsInt  = PVI->getPU_NumInteractions();
	nTrueInt = PVI->getTrueNumInteractions();
	break;
      } 
    }
    PUweight = pileUpReweight->weight(nTrueInt);
  }
  

   // CR Candidates
  Handle<View<pat::CompositeCandidate> > candHandle;
  event.getByToken(candToken, candHandle);
  const View<pat::CompositeCandidate>* cands = candHandle.product();

 
  //----------------------------------------------------------------------
  // Loop on CR candidates
  //----------------------------------------------------------------------

  View<pat::CompositeCandidate>::const_iterator bestCand = cands->end();
  float maxPtSum = -1.;
  for( View<pat::CompositeCandidate>::const_iterator cand = cands->begin(); cand != cands->end(); ++cand) {
    //----------------------------------------------------------------------
    // Choosing the best CR Cand
    //----------------------------------------------------------------------
    
    // Look for the daughter with largesr-pT leptons
    float ptSum = cand->daughter(1)->daughter(0)->pt()+ cand->daughter(1)->daughter(1)->pt();
    if (ptSum > maxPtSum){
      maxPtSum = ptSum;
      bestCand = cand;
            
    } 
    
  } // End of loop on candidates

  if (bestCand == cands->end()) return;
  
  //--- Extract all relevant best candidate information 
  
  // Pointers to Z and leptons, prefixes for UserFloats
  const Candidate* Z1   = bestCand->daughter(0);
  const Candidate* Z2   = bestCand->daughter(1);
  vector<const Candidate*> leptons(4);
  vector<const Candidate*> fsr;
  vector<short> fsrIndex;
  vector<string> labels(4);
  userdatahelpers::getSortedLeptons(*bestCand, leptons, labels, fsr, fsrIndex);
    
  // Retrieve the userFloat of the leptons in vectors ordered in the same way.
  vector<float> SIP(4);
//   vector<float> looseIso(4);
//   vector<float> iso(4);
  vector<float> isoPF(4);
  vector<float> pt(4);
  for (int i=0; i<4; ++i){
    SIP[i]      = userdatahelpers::getUserFloat(leptons[i],"SIP");
//     looseIso[i] = userdatahelpers::getUserFloat(leptons[i],"looseIso");
//     iso[i]      = userdatahelpers::getUserFloat(leptons[i],"combRelIso");
    isoPF[i]    = bestCand->userFloat(labels[i]+"combRelIsoPFFSRCorr"); // Note: the FSR-corrected iso is attached to the Z, not to the lepton!

    pt[i]       = leptons[i]->pt();
  }
  
  
  // Sort lepton variables
  vector<float> ptS(pt);
  //  vector<float> isoS(iso);
  vector<float> isoPFS(isoPF);
  vector<float> SIPS(SIP);
  sort(ptS.begin(),ptS.end());
  //  sort(isoS.begin(),isoS.end());
  sort(isoPFS.begin(),isoPFS.end());
  sort(SIPS.begin(),SIPS.end());
  
  
  // Masses
  float ZZMass = bestCand->p4().mass();
//  float ZZMassErr = bestCand->userFloat("massError");
  float Z1Mass = Z1->mass();
  float Z2Mass = Z2->mass();	
  float Z1OtherCombMass = bestCand->userFloat("mZa");
  float Z2OtherCombMass = bestCand->userFloat("mZb");
  
  // Cut variables
// float relIso_sum2least = bestCand->userFloat("iso34");
// float SIP4      = bestCand->userFloat("SIP4");  
  
  hCandCR->Fill(ZZMass, Z1Mass, Z2Mass, Z1OtherCombMass, Z2OtherCombMass, ptS, isoPFS, SIPS, 1.);
  hCandCR_w->Fill(ZZMass, Z1Mass, Z2Mass, Z1OtherCombMass, Z2OtherCombMass, ptS, isoPFS, SIPS, PUweight);
  
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ZZ4lAnalyzerCR);

