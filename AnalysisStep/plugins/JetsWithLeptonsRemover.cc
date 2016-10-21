/** \class to remove jets overlapping with leptons
 *  Producers to purge jets made with isolated leptons.
 *
 *  $Date:  $
 *  $Revision: $
 *  \author R. Covarelli
 */


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include <ZZAnalysis/AnalysisStep/interface/PhotonFwd.h>
#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/Vector3D.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

class  JetsWithLeptonsRemover: public edm::EDProducer {
public:
  
  enum MatchingType{byConstituents, byDeltaR};

  explicit JetsWithLeptonsRemover(const edm::ParameterSet & iConfig);
  virtual ~JetsWithLeptonsRemover() { }

  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);
  bool checkLeptonJet(const edm::Event & event, const pat::Jet& jet);
  bool isMatchingWithZZLeptons(const edm::Event & event, const pat::Jet& jet);
  template <typename LEP>
  bool isMatchingWith(const edm::EDGetTokenT<edm::View<LEP> >& token, const StringCutObjectSelector<LEP>& presel, const edm::Event & event, const pat::Jet& jet);

private:
  /// Labels for input collections
  MatchingType matchingType_;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetToken_;
  edm::EDGetTokenT<edm::View<pat::Muon> > muonToken_;
  edm::EDGetTokenT<edm::View<pat::Electron> > electronToken_;
  edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > diBosonToken_;
  bool cleaningFromDiboson_;
 
  /// Preselection cut
  StringCutObjectSelector<pat::Jet> preselectionJ_;
  StringCutObjectSelector<pat::Muon> preselectionMu_;
  StringCutObjectSelector<pat::Electron> preselectionEle_;
  StringCutObjectSelector<pat::CompositeCandidate> preselectionVV_;

  bool cleanFSRFromLeptons_;

  // Some istograms for monitoring
  bool  activateDebugPrintOuts_;
  bool  doDebugPlots_;
  double theDeltaRCut_;

  TH1F *hNLeptonJets;
  TH1F *hDeltaPt_jet_lepton;
  TH1F *hDeltaPt_jetcomp_lepton;
  TH1F *hDeltaPt_jet_fsr;
  TH1F *hDeltaPt_jetcomp_fsr;
  TH1F *hDeltaPhi_jet_lepton;
  TH1F *hDeltaPhi_jet_fsr;
  TH1F *hDeltaEta_jet_lepton;
  TH1F *hDeltaEta_jetcomp_lepton;
  TH1F *hDeltaEta_jet_fsr;
  TH1F *hDeltaEta_jetcomp_fsr;
};


JetsWithLeptonsRemover::JetsWithLeptonsRemover(const edm::ParameterSet & iConfig)
  : jetToken_     (consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("Jets")))
  , muonToken_    (consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("Muons")))
  , electronToken_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("Electrons")))
  , diBosonToken_ (consumes<edm::View<pat::CompositeCandidate> >(iConfig.getParameter<edm::InputTag>("Diboson")))
  , preselectionJ_    (iConfig.getParameter<std::string>("JetPreselection"))
  , preselectionMu_   (iConfig.getParameter<std::string>("MuonPreselection"))
  , preselectionEle_  (iConfig.getParameter<std::string>("ElectronPreselection"))
  , preselectionVV_   (iConfig.getParameter<std::string>("DiBosonPreselection"))
  , cleanFSRFromLeptons_   (iConfig.getParameter<bool>("cleanFSRFromLeptons"))
  , activateDebugPrintOuts_ (iConfig.getUntrackedParameter<bool>("DebugPrintOuts",false))   
  , doDebugPlots_           (iConfig.getUntrackedParameter<bool>("DebugPlots",false))
  , theDeltaRCut_  (iConfig.getUntrackedParameter<double>("DeltaRCut",0.4))
{

  std::string matchingType = iConfig.getParameter<std::string>("MatchingType"); 
  if(matchingType == "byConstituents") matchingType_ = JetsWithLeptonsRemover::byConstituents;
  else if(matchingType == "byDeltaR")  matchingType_ = JetsWithLeptonsRemover::byDeltaR;
  else std::cout << "Not making any matching, the matching you choose is not foreseen: " << matchingType << std::endl; 
  
  produces<std::vector<pat::Jet> >(); 

  edm::InputTag diBosonSrc = iConfig.getParameter<edm::InputTag>("Diboson");
  if(diBosonSrc.label() != "") cleaningFromDiboson_ = true;
  else cleaningFromDiboson_ = false;

  if(doDebugPlots_){
    edm::Service<TFileService> fs;
    hNLeptonJets              = fs->make<TH1F>("hNLeptonJets"            , "Number of lepton-jets found",  10,   0, 10);
    hDeltaPt_jet_lepton       = fs->make<TH1F>("hDeltaPt_jet_lepton"     , "#Delta p_T (jet, l)"        , 100, -50, 50);
    hDeltaPt_jetcomp_lepton   = fs->make<TH1F>("hDeltaPt_jetcomp_lepton" , "#Delta p_T (jetcomp, l)"    , 100, -50, 50);
    hDeltaPt_jet_fsr          = fs->make<TH1F>("hDeltaPt_jet_fsr"        , "#Delta p_T (jet, fsr)"      , 100, -50, 50);
    hDeltaPt_jetcomp_fsr      = fs->make<TH1F>("hDeltaPt_jetcomp_fsr"    , "#Delta p_T (jetcomp, fsr)"  , 100, -50, 50);
    hDeltaPhi_jet_lepton      = fs->make<TH1F>("hDeltaPhi_jet_lepton"    , "#Delta #phi (jet, l)"       , 100,  -4,  4);
    hDeltaPhi_jet_fsr	      = fs->make<TH1F>("hDeltaPhi_jet_fsr"       , "#Delta #phi (jet, fsr)"     , 100,  -4,  4);
    hDeltaEta_jet_lepton      = fs->make<TH1F>("hDeltaEta_jet_lepton"    , "#Delta #eta (jet, l)"       , 100,  -5,  5);
    hDeltaEta_jetcomp_lepton  = fs->make<TH1F>("hDeltaEta_jetcomp_lepton", "#Delta #eta (jetcomp, l)"   , 100,  -5,  5);
    hDeltaEta_jet_fsr	      = fs->make<TH1F>("hDeltaEta_jet_fsr"       , "#Delta #eta (jet, fsr)"     , 100,  -5,  5);
    hDeltaEta_jetcomp_fsr     = fs->make<TH1F>("hDeltaEta_jetcomp_fsr"   , "#Delta #eta (jetcomp, fsr)" , 100,  -5,  5);
  }
}

void JetsWithLeptonsRemover::produce(edm::Event & event, const edm::EventSetup & iSetup) {
  using namespace edm;
  using namespace std;
  
  Handle<edm::View<pat::Jet> > jets;
  event.getByToken(jetToken_, jets);


  //std::cout<<"----------- Muon -----------"<<std::endl;
  //foreach(const pat::Muon& muon, *muons)
  //  std::cout<<"pt: " << muon.pt() << " eta: " << muon.eta() << " phi: " << muon.phi() << " p: " << muon.p() <<std::endl;
  

  if(activateDebugPrintOuts_) std::cout << "\n\n----------- NEW EVENT ----------- number of jets: " << jets->size() << std::endl;
  int passPresel = 0;
  int numLepJets = 0;
  auto_ptr<vector<pat::Jet> > out(new vector<pat::Jet>());
  foreach(const pat::Jet& jet, *jets){

    if(!preselectionJ_(jet)) continue;
    ++passPresel;
      
    if(activateDebugPrintOuts_) std::cout<<"\n+++++ Jet +++++ pt: " << jet.pt() << " eta: " << jet.eta() << " phi: " << jet.phi() << std::endl;
    
    if(checkLeptonJet(event, jet))
      ++numLepJets;
    else
      out->push_back(jet);
  }

  if(doDebugPlots_) hNLeptonJets->Fill(numLepJets);
  
  if(activateDebugPrintOuts_) std::cout << "Pass Presel: " << passPresel << " pass cleaning: " << out->size()
					<< "\n-------------------------------------------------------------------------" << std::endl;
  
  event.put(out);
}


bool JetsWithLeptonsRemover::isMatchingWithZZLeptons(const edm::Event & event, const pat::Jet& jet) {

  edm::Handle<edm::View<pat::CompositeCandidate> > VV; event.getByToken(diBosonToken_, VV);
    pat::CompositeCandidate bestVV; bool found = false;

    // Search for the best ZZ pair that satisfy the preselection requirements
    foreach(const pat::CompositeCandidate &vv, *VV)
      if (preselectionVV_(vv)){
	bestVV = vv; 
	found = true;
	break;
      }
      
    
    if(found) { 
      // loop over the Zs
      for(int i=0; i<2; ++i){
	
	const pat::CompositeCandidate *v =  dynamic_cast<const pat::CompositeCandidate*>(bestVV.daughter(i)->masterClone().get());
	
	// loop over the leptons of each Z
	for(int j=0; j<2; ++j){
	  bool checkingVariable = false;
	  if (matchingType_ == JetsWithLeptonsRemover::byDeltaR)
	    checkingVariable = reco::deltaR(*v->daughter(j), jet) < theDeltaRCut_;
	  else
	    std::cout << "Not making any matching, the matching you choose is not foreseen" << std::endl;
	    
	  
	  if(checkingVariable){
	    if(activateDebugPrintOuts_) std::cout << "\t\t !!! Found a matching lepton-jet (from VV candidate) !!!"<<std::endl;
	    if(doDebugPlots_){
	      hDeltaPt_jet_lepton     ->Fill(v->daughter(j)->pt()  - jet.pt());
	      hDeltaPhi_jet_lepton    ->Fill(v->daughter(j)->phi() - jet.phi());
	      hDeltaEta_jet_lepton    ->Fill(v->daughter(j)->eta() - jet.eta());  
	    }
	    return true;
	  }
	}
	  
	// Check if the jet matches FSR photons
	for (unsigned jfsr=2; jfsr<v->numberOfDaughters(); ++jfsr) {
	  const pat::PFParticle* fsr = static_cast<const pat::PFParticle*>(v->daughter(jfsr));
	  double photon_en_frac = v->daughter(jfsr)->energy()/jet.energy();
	  if(activateDebugPrintOuts_) {
	    int ilep = fsr->userFloat("leptIdx");
	    std::cout << "Sister of " << ilep << " (" <<  fsr->pdgId() << "), pt: "
		      << v->daughter(jfsr)->pt()   << " eta: " << v->daughter(jfsr)->eta() << " phi: " << fsr->phi() << " p: " << fsr->p()
		      << " Photon energy fraction in the jet: " <<  photon_en_frac 
		      << std::endl;
	  }

	  // If the FSR photon matches, reject the jet
	  if(jet.photonMultiplicity() > 0 && photon_en_frac > 0.5 && reco::deltaR(*fsr, jet) < 0.05){
	    if(activateDebugPrintOuts_) std::cout << "\t\t !!! Found a matching FSR lepton-jet (from VV candidate) !!!"<<std::endl;	  
	    if(doDebugPlots_){
	      hDeltaPt_jet_fsr     ->Fill(fsr->pt()  - jet.pt());
	      hDeltaPhi_jet_fsr    ->Fill(fsr->phi() - jet.phi());
	      hDeltaEta_jet_fsr    ->Fill(fsr->eta() - jet.eta());  
	    }
	    return true;
	  }
	}
      }
    }
    return false;
}


template <typename LEP>
bool JetsWithLeptonsRemover::isMatchingWith(const edm::EDGetTokenT<edm::View<LEP> >& token, const StringCutObjectSelector<LEP>& presel, const edm::Event & event, const pat::Jet& jet){

  // Check for muon-originated jets   
  edm::Handle<edm::View<LEP> > leptons; event.getByToken(token, leptons);
  
  foreach(const LEP& lepton, *leptons){

    if(!presel(lepton)) continue;
    
    if(activateDebugPrintOuts_) std::cout<<"Lepton pt: " << lepton.pt()   << " eta: " << lepton.eta()    << " phi: " << lepton.phi() << " p: " << lepton.p() << std::endl;  
      
    bool checkingVariable = false;
    if (matchingType_ == JetsWithLeptonsRemover::byDeltaR)
      checkingVariable = reco::deltaR(lepton, jet) < theDeltaRCut_;
    else
      std::cout << "Not making any matching, the matching you choose is not foreseen" << std::endl;
    
    if(checkingVariable){
      if(activateDebugPrintOuts_) std::cout << "\t\t !!! Found a matching lepton-jet !!!"<<std::endl;
      return true;
    }

    if (cleanFSRFromLeptons_) {
      const PhotonPtrVector* gammas = userdatahelpers::getUserPhotons(&lepton);
      if (gammas==0) continue;
      assert(gammas->size()<=1); // Must have already been preselected, so there should be at most 1 per l
      if (gammas->size()==1){
	const pat::PFParticle* fsr = gammas->begin()->get();
	if (matchingType_ == JetsWithLeptonsRemover::byDeltaR) {
	  checkingVariable = (reco::deltaR(*fsr, jet) < theDeltaRCut_);
	}
	if(checkingVariable){
	  if(activateDebugPrintOuts_) {
	    double photon_en_frac = fsr->energy()/jet.energy();
	    std::cout << "\t\t !!! Found a matching FSR-jet !!! " <<  fsr->energy() << " " << jet.energy() << " " << photon_en_frac<<  std::endl;
	  }
	  return true;
	}
      }
    }


  }
  return false;
}


bool JetsWithLeptonsRemover::checkLeptonJet(const edm::Event & event, const pat::Jet& jet){

  if(cleaningFromDiboson_ && isMatchingWithZZLeptons(event,jet)) return true;
  
  if(isMatchingWith<pat::Muon>    (muonToken_,     preselectionMu_,  event, jet) || 
     isMatchingWith<pat::Electron>(electronToken_, preselectionEle_, event, jet)) return true;
  
  return false;
}









#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(JetsWithLeptonsRemover);
