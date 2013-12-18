/** \class MuFiller
 *
 *  No description available.
 *
 *  $Date: 2013/05/14 10:08:19 $
 *  $Revision: 1.18 $
 *  \author N. Amapane (Torino)
 *  \author S. Bolognesi (JHU)
 *  \author C. Botta (CERN)
 *  \author S. Casasso (Torino)
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>

#include <DataFormats/PatCandidates/interface/Muon.h>
#include <Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h>
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>
#include <ZZAnalysis/AnalysisStep/interface/MCHistoryTools.h>
#include <ZZAnalysis/AnalysisStep/interface/SIPCalculator.h>


#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;


class MuFiller : public edm::EDProducer {
 public:
  /// Constructor
  explicit MuFiller(const edm::ParameterSet&);
    
  /// Destructor
  ~MuFiller(){};  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  const edm::InputTag theCandidateTag;
  int sampleType;
  int setup;
  const StringCutObjectSelector<pat::Muon, true> cut;
  const CutSet<pat::Muon> flags;
  SIPCalculator *sipCalculator_;

};


MuFiller::MuFiller(const edm::ParameterSet& iConfig) :
  theCandidateTag(iConfig.getParameter<InputTag>("src")),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  cut(iConfig.getParameter<std::string>("cut")),
  flags(iConfig.getParameter<edm::ParameterSet>("flags"))
{
  sipCalculator_ = new SIPCalculator();
  produces<pat::MuonCollection>();
}


void
MuFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //Initialize SIP calculator
  sipCalculator_->initialize(iSetup);

  //--- Get leptons and rho
  edm::Handle<pat::MuonRefVector> muonHandle;
  iEvent.getByLabel(theCandidateTag, muonHandle);

  InputTag theRhoTag = LeptonIsoHelper::getMuRhoTag(sampleType, setup);
  edm::Handle<double> rhoHandle;
  iEvent.getByLabel(theRhoTag, rhoHandle);
  double rho = *rhoHandle;

  edm::Handle<vector<Vertex> >  vertexs;
  iEvent.getByLabel("goodPrimaryVertices",vertexs);

  // Output collection
  auto_ptr<pat::MuonCollection> result( new pat::MuonCollection() );

  //FIXME: effective areas to be updated!
  const float AreaEcal[2]    = {0.074, 0.045};   //   barrel/endcap
  const float AreaHcal[2]    = {0.022, 0.030};   //   barrel/endcap

  for (unsigned int i = 0; i< muonHandle->size(); ++i){
    //---Clone the pat::Muon
    pat::Muon l(*((*muonHandle)[i].get()));

    //--- Rho-corrected isolation and loose iso
    float tkIso   = l.userIsolation(pat::User1Iso);
    float ecalIso = l.ecalIso();
    float hcalIso = l.hcalIso();
    float feta = fabs(l.eta());
    float pt = l.pt();

    Int_t ifid = (feta < 1.479) ? 0 : 1;
    ecalIso = ecalIso - AreaEcal[ifid]*rho;
    hcalIso = hcalIso - AreaHcal[ifid]*rho;
    
    float combRelIso = (ecalIso + hcalIso + tkIso)/pt;

    //--- PF ISO
    float PFChargedHadIso   = l.chargedHadronIso();
    float PFNeutralHadIso   = l.neutralHadronIso();
    float PFPhotonIso       = l.photonIso();
    
    float combRelIsoPF = LeptonIsoHelper::combRelIsoPF(sampleType, setup, rho, l);

    bool isGlb = l.isGlobalMuon();
    bool isTk = l.isTrackerMuon();
    
    float mvaIsoRings = l.userFloat("mvaIsoRings");
    bool isMvaIsoRings = false;
    // "reference WP" from https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateMuonSelection#Reference_Working_Point r15
    if ((isGlb && isTk && pt < 10 && feta < 1.5   && mvaIsoRings > -0.593) ||
	(isGlb && isTk && pt >= 10 && feta < 1.5  && mvaIsoRings > 0.337) ||
	(isGlb && isTk && pt < 10 && feta >= 1.5  && mvaIsoRings > -0.767) ||
	(isGlb && isTk && pt >= 10 && feta >= 1.5 && mvaIsoRings > 0.410) ||
	(!isGlb && isTk                           && mvaIsoRings > -0.989) ||
	(isGlb && !isTk                           && mvaIsoRings > -0.995)) {
      isMvaIsoRings = true;
    }
    


    //--- SIP, dxy, dz
    float IP      = fabs(l.dB(pat::Muon::PV3D));
    float IPError = l.edB(pat::Muon::PV3D);
    float SIP     = IP/IPError;

    float dxy = 999.;
    float dz  = 999.;
    const Vertex* vertex = 0;
    if (vertexs->size()>0) {
      vertex = &(vertexs->front());
      dxy = fabs(l.innerTrack()->dxy(vertex->position()));
      dz  = fabs(l.innerTrack()->dz(vertex->position()));
    }

    //--- Trigger matching
    bool HLTMatch = ((!l.triggerObjectMatchesByFilter("hltSingleMu13L3Filtered17").empty())||
		     ((!l.triggerObjectMatchesByFilter("hltSingleMu13L3Filtered13").empty()) && 
		      (l.triggerObjectMatchesByFilter("hltSingleMu13L3Filtered13").at(0).pt()>17)) || 
		     ((!l.triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").empty()) && 
		      (l.triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").at(0).pt()>17)) || 
		     ((!l.triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").empty()) && 
		      (l.triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").at(0).pt()>17)));
    //FIXME
    
    //--- Embed user variables
    l.addUserFloat("looseIso",tkIso/l.pt());
    l.addUserFloat("combRelIso",combRelIso);
    l.addUserFloat("PFChargedHadIso",PFChargedHadIso);
    l.addUserFloat("PFNeutralHadIso",PFNeutralHadIso);
    l.addUserFloat("PFPhotonIso",PFPhotonIso);
    l.addUserFloat("combRelIsoPF",combRelIsoPF);
    l.addUserFloat("rho",rho);
    l.addUserFloat("isMvaIsoRings", isMvaIsoRings);    

    if(vertexs->size()>0) {
      SIP=sipCalculator_->calculate(l,vertexs->front());
    }
    l.addUserFloat("SIP",SIP);

    l.addUserFloat("dxy",dxy);
    l.addUserFloat("dz",dz);
    l.addUserFloat("HLTMatch", HLTMatch);
    // l.addUserCand("MCMatch",genMatch); // FIXME

    //--- isPFMuon flag - in old samples, l.isPFMuon() is not functional, so this has to be filled
    //    beforehand with the module PATPFMuonEmbedder.
    if(!l.hasUserFloat("isPFMuon")) {
      l.addUserFloat("isPFMuon",l.isPFMuon());
    }
    
    //--- MC parent code 
    MCHistoryTools mch(iEvent);
    if (mch.isMC()) {
      int MCParentCode = 0;//FIXME: does not work on cmg mch.getParentCode((l.genParticleRef()).get());
      l.addUserFloat("MCParentCode",MCParentCode);
    }

    //--- Check selection cut. Being done here, flags are not available; but this way we 
    //    avoid wasting time on rejected leptons.
    if (!cut(l)) continue;

    //--- Embed flags (ie flags specified in the "flags" pset)
    for(CutSet<pat::Muon>::const_iterator flag = flags.begin(); flag != flags.end(); ++flag) {
      l.addUserFloat(flag->first,int((*(flag->second))(l)));
    }
    
    result->push_back(l);
  }
  iEvent.put(result);
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(MuFiller);

