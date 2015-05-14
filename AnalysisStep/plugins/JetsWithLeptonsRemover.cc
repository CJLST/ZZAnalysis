/** \class to remove jets overlapping with leptons
 *  Producers to purge jets made with isolated leptons.
 *
 *  $Date:  $
 *  $Revision: $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
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
  bool isGood(const pat::Jet &jet) const;
  bool checkLeptonJet(const edm::Event & event, const pat::Jet& jet);
  bool isMatchingWithZZLeptons(const edm::Event & event, const pat::Jet& jet);
  template <typename LEP>
  bool isMatchingWith(const edm::InputTag& src, const StringCutObjectSelector<LEP>& presel, const edm::Event & event, const pat::Jet& jet);

private:
  /// Labels for input collections
  MatchingType matchingType_;
  edm::InputTag jetSrc_;
  edm::InputTag muonSrc_;
  edm::InputTag electronSrc_;
  edm::InputTag diBosonSrc_;
  bool cleaningFromDiboson_;
 
  /// Preselection cut
  StringCutObjectSelector<pat::Jet> preselectionJ_;
  StringCutObjectSelector<pat::Muon> preselectionMu_;
  StringCutObjectSelector<pat::Electron> preselectionEle_;
  StringCutObjectSelector<pat::CompositeCandidate> preselectionVV_;

  // Some istograms for monitoring
  bool  activateDebugPrintOuts_;
  bool  doDebugPlots_;
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
  : jetSrc_           (iConfig.getParameter<edm::InputTag>("Jets"))
  , muonSrc_          (iConfig.getParameter<edm::InputTag>("Muons"))
  , electronSrc_      (iConfig.getParameter<edm::InputTag>("Electrons"))
  , diBosonSrc_       (iConfig.getParameter<edm::InputTag>("Diboson"))  
  , preselectionJ_    (iConfig.getParameter<std::string>("JetPreselection"))
  , preselectionMu_   (iConfig.getParameter<std::string>("MuonPreselection"))
  , preselectionEle_  (iConfig.getParameter<std::string>("ElectronPreselection"))
  , preselectionVV_   (iConfig.getParameter<std::string>("DiBosonPreselection"))

  , activateDebugPrintOuts_ (iConfig.getUntrackedParameter<bool>("DebugPrintOuts",false))   
  , doDebugPlots_           (iConfig.getUntrackedParameter<bool>("DebugPlots",false)) 
{

  std::string matchingType = iConfig.getParameter<std::string>("MatchingType");
  if(matchingType == "byConstituents") matchingType_ = JetsWithLeptonsRemover::byConstituents;
  else if(matchingType == "byDeltaR")  matchingType_ = JetsWithLeptonsRemover::byDeltaR;
  else std::cout << "Not making any matching, the matching you choose is not foreseen: " << matchingType << std::endl; 
  
  produces<std::vector<pat::Jet> >(); 

  if(diBosonSrc_.label() != "") cleaningFromDiboson_ = true;
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
  event.getByLabel(jetSrc_, jets);


  //std::cout<<"----------- Muon -----------"<<std::endl;
  //foreach(const pat::Muon& muon, *muons)
  //  std::cout<<"pt: " << muon.pt() << " eta: " << muon.eta() << " phi: " << muon.phi() << " p: " << muon.p() <<std::endl;
  

  if(activateDebugPrintOuts_) std::cout << "\n\n----------- NEW EVENT ----------- number of jets: " << jets->size() << std::endl;
  int passPresel = 0;
  int numLepJets = 0;
  auto_ptr<vector<pat::Jet> > out(new vector<pat::Jet>());
  foreach(const pat::Jet& jet, *jets){

    if(!preselectionJ_(jet) || !isGood(jet)) continue;
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



bool JetsWithLeptonsRemover::isGood(const pat::Jet& jet) const {
  
  float jeta=fabs(jet.eta());

  // cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
  bool looseJetID = (jet.neutralHadronEnergyFraction() < 0.99 && 
		     jet.neutralEmEnergyFraction() < 0.99 &&
		     jet.numberOfDaughters() > 1 &&
		     jet.chargedMultiplicity()+jet.neutralMultiplicity() >0 &&
		     jet.muonEnergyFraction() < 0.8 &&
		     jet.chargedEmEnergyFraction() < 0.9 &&
		     ( jet.chargedHadronEnergyFraction() > 0 || jeta > 2.4 )  &&
		     ( jet.chargedMultiplicity() > 0 || jeta > 2.4 ) &&
		     ( jet.chargedEmEnergyFraction() < 0.99 || jeta > 2.4 ) );
  if(!looseJetID) return false;

  bool passPU = true;
  
  float jpt=jet.pt();
  float jpumva=jet.userFloat("pileupJetId:fullDiscriminant");

  if(jpt>20){
    if(jeta > 3.)       { if(jpumva <= -0.45) passPU=false;}
    else if(jeta > 2.75){ if(jpumva <= -0.55) passPU=false;}
    else if(jeta > 2.50){ if(jpumva <= -0.60) passPU=false;}
    else                { if(jpumva <= -0.63) passPU=false;}
  }
  else{
    if(jeta > 3.)       { if(jpumva <= -0.95) passPU=false;}
    else if(jeta > 2.75){ if(jpumva <= -0.94) passPU=false;}
    else if(jeta > 2.50){ if(jpumva <= -0.96) passPU=false;}
    else                { if(jpumva <= -0.95) passPU=false;}
  }
  

  return looseJetID && passPU;
}


bool JetsWithLeptonsRemover::isMatchingWithZZLeptons(const edm::Event & event, const pat::Jet& jet) {

  edm::Handle<edm::View<pat::CompositeCandidate> > VV   ; event.getByLabel(diBosonSrc_, VV);
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
	    checkingVariable = reco::deltaR(*v->daughter(j), jet) < 0.4;
	  else
	    std::cout << "Not making any matching, the matching you choose is not foreseen" << std::endl;
	    
	  
	  if(checkingVariable){
	    if(activateDebugPrintOuts_) std::cout << "\t\t !!! Found a matching lepton-jet !!!"<<std::endl;
	    if(doDebugPlots_){
	      hDeltaPt_jet_lepton     ->Fill(v->daughter(j)->pt()  - jet.pt());
	      hDeltaPhi_jet_lepton    ->Fill(v->daughter(j)->phi() - jet.phi());
	      hDeltaEta_jet_lepton    ->Fill(v->daughter(j)->eta() - jet.eta());  
	    }
	    return true;
	  }
	}
	  
	// Check if the jet matches FSR photons
	if(v->hasUserFloat("dauWithFSR") && v->userFloat("dauWithFSR") >= 0){
	    
	  double photon_en_frac = v->daughter(2)->energy()/jet.energy();
	  if(activateDebugPrintOuts_)
	    std::cout << "Sister of " << v->userFloat("dauWithFSR") << " (" <<  v->daughter(v->userFloat("dauWithFSR"))->pdgId() << "), pt: "
		      << v->daughter(2)->pt()   << " eta: " << v->daughter(2)->eta() << " phi: " << v->daughter(2)->phi() << " p: " << v->daughter(2)->p()
		      << " Photon energy fraction in the jet: " <<  photon_en_frac 
		      << std::endl;

	  // If the FSR photon matches, reject the jet
	  if(jet.photonMultiplicity() > 0 && photon_en_frac > 0.5 && reco::deltaR(*v->daughter(2), jet) < 0.05){
	    if(activateDebugPrintOuts_) std::cout << "\t\t !!! Found a matching FSR lepton-jet !!!"<<std::endl;	  
	    if(doDebugPlots_){
	      hDeltaPt_jet_fsr     ->Fill(v->daughter(2)->pt()  - jet.pt());
	      hDeltaPhi_jet_fsr    ->Fill(v->daughter(2)->phi() - jet.phi());
	      hDeltaEta_jet_fsr    ->Fill(v->daughter(2)->eta() - jet.eta());  
	    }
	    return true;
	  }
	}
      }
    }
    return false;
}


template <typename LEP>
bool JetsWithLeptonsRemover::isMatchingWith(const edm::InputTag& src, const StringCutObjectSelector<LEP>& presel, const edm::Event & event, const pat::Jet& jet){

  // Check for muon-originated jets   
  edm::Handle<std::vector<LEP> >  leptons; event.getByLabel(src, leptons);
  
  foreach(const LEP& lepton, *leptons){

    if(!presel(lepton)) continue;
    
    if(activateDebugPrintOuts_) std::cout<<"Lepton pt: " << lepton.pt()   << " eta: " << lepton.eta()    << " phi: " << lepton.phi() << " p: " << lepton.p() << std::endl;  
      
    bool checkingVariable = false;
    if (matchingType_ == JetsWithLeptonsRemover::byDeltaR)
      checkingVariable = reco::deltaR(lepton, jet) < 0.4;
    else
      std::cout << "Not making any matching, the matching you choose is not foreseen" << std::endl;
    
    if(checkingVariable){
      if(activateDebugPrintOuts_) std::cout << "\t\t !!! Found a matching lepton-jet (muon not coming from ZZ decay) !!!"<<std::endl;
      return true;
    }
  }
  return false;
}


bool JetsWithLeptonsRemover::checkLeptonJet(const edm::Event & event, const pat::Jet& jet){

  if(cleaningFromDiboson_ && isMatchingWithZZLeptons(event,jet)) return true;
  
  if(isMatchingWith<pat::Muon>    (muonSrc_,     preselectionMu_,  event, jet) || 
     isMatchingWith<pat::Electron>(electronSrc_, preselectionEle_, event, jet)) return true;
  
  return false;
}









#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(JetsWithLeptonsRemover);
