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
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>

#include "MuonMVAReader/Reader/interface/MuonGBRForestReader.hpp"

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
   ~MuFiller();
   
   private:
   virtual void beginJob(){};
   virtual void produce(edm::Event&, const edm::EventSetup&);
   virtual void endJob(){};
   
   edm::EDGetTokenT<pat::MuonRefVector> muonToken;
   int sampleType;
   int setup;
   const StringCutObjectSelector<pat::Muon, true> cut;
   const edm::EDGetTokenT< edm::TriggerResults > triggerResultsToken_;
   const CutSet<pat::Muon> flags;
   edm::EDGetTokenT<double> rhoToken;
   edm::EDGetTokenT<vector<Vertex> > vtxToken;
   
   // MVA Reader
   MuonGBRForestReader *r;
};


MuFiller::MuFiller(const edm::ParameterSet& iConfig) :
muonToken(consumes<pat::MuonRefVector>(iConfig.getParameter<edm::InputTag>("src"))),
sampleType(iConfig.getParameter<int>("sampleType")),
setup(iConfig.getParameter<int>("setup")),
cut(iConfig.getParameter<std::string>("cut")),
triggerResultsToken_( consumes< edm::TriggerResults >( iConfig.getParameter< edm::InputTag >( "TriggerResults" ) ) ),
flags(iConfig.getParameter<edm::ParameterSet>("flags"))
{
   rhoToken = consumes<double>(LeptonIsoHelper::getMuRhoTag(sampleType, setup));
   vtxToken = consumes<vector<Vertex> >(edm::InputTag("goodPrimaryVertices"));
   produces<pat::MuonCollection>();
   
   // MVA Reader
   r = new MuonGBRForestReader(setup, 2);
   
}

MuFiller::~MuFiller(){}


void
MuFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   
   //--- Get leptons and rho
   edm::Handle<pat::MuonRefVector> muonHandle;
   iEvent.getByToken(muonToken, muonHandle);
   
   edm::Handle< edm::TriggerResults > triggerResults;
   iEvent.getByToken( triggerResultsToken_, triggerResults );
   
   edm::Handle<double> rhoHandle;
   iEvent.getByToken(rhoToken, rhoHandle);
   double rho = *rhoHandle;
   
   edm::Handle<vector<Vertex> > vertices;
   iEvent.getByToken(vtxToken,vertices);
   
   
   // Output collection
   auto result = std::make_unique<pat::MuonCollection>();
   
   for (unsigned int i = 0; i< muonHandle->size(); ++i){
      //---Clone the pat::Muon
      pat::Muon l(*((*muonHandle)[i].get()));
      
      //--- PF ISO
      // for cone size R=0.3 :
      float PFChargedHadIso   = l.pfIsolationR03().sumChargedHadronPt;
      float PFNeutralHadIso   = l.pfIsolationR03().sumNeutralHadronEt;
      float PFPhotonIso       = l.pfIsolationR03().sumPhotonEt;
      float PFPUChargedHadIso = l.pfIsolationR03().sumPUPt;
      
      float combRelIsoPF = LeptonIsoHelper::combRelIsoPF(sampleType, setup, rho, l);
      
      //--- SIP, dxy, dz
      float IP      = std::abs(l.dB(pat::Muon::PV3D));
      float IPError = l.edB(pat::Muon::PV3D);
      float SIP     = IP/IPError;
      
      float dxy = 999.;
      float dz  = 999.;
      const Vertex* vertex = 0;
      if (vertices->size()>0) {
         vertex = &(vertices->front());
         dxy = fabs(l.muonBestTrack()->dxy(vertex->position()));
         dz  = fabs(l.muonBestTrack()->dz(vertex->position()));
      }
      
      
//=================
// Begin MVA reader
//=================

      float glb_valid_mu_hits_, glb_chi2_, tk_valid_pixel_hits_, tk_valid_hits_;
//      float glb_pT_error_, tk_tracker_lay_, tk_pixel_lay_, tk_chi2_, tk_pT_error_;
      
//      bool is_global_mu_  = l.isGlobalMuon();
      
      // GlobalTrackQualityVariables
      if ( l.globalTrack().isNonnull() )
      {
         glb_valid_mu_hits_ = (l.globalTrack()->hitPattern().numberOfValidMuonHits());
         glb_chi2_          = (l.globalTrack()->normalizedChi2());
//         glb_pT_error_      = (l.globalTrack()->ptError()); 
      }
      else
      {
         glb_valid_mu_hits_ = -1;
         glb_chi2_          = -1;
//         glb_pT_error_      = -1;
      }
      
      
      // TrackQualityVariables
      reco::TrackRef track = l.innerTrack();
//      valid_KF = (myTrackRef.isAvailable());
//      valid_KF = (myTrackRef.isNonnull());


      if ( track.isNonnull() )
      {
         tk_valid_hits_       = l.innerTrack()->hitPattern().numberOfValidHits();
         tk_valid_pixel_hits_ = l.innerTrack()->hitPattern().numberOfValidPixelHits();
//         tk_tracker_lay_      = l.innerTrack()->hitPattern().trackerLayersWithMeasurement();
//         tk_pixel_lay_        = l.innerTrack()->hitPattern().pixelLayersWithMeasurement();
//         tk_chi2_             = l.innerTrack()->normalizedChi2();
//         tk_pT_error_         = l.innerTrack()->ptError();
      }
      else
      {
         tk_valid_hits_       = -1;
         tk_valid_pixel_hits_ = -1;
//         tk_tracker_lay_      = -1;
//         tk_pixel_lay_        = -1;
//         tk_chi2_             = -1;
//         tk_pT_error_         = -1;
      }

      
      float BDT = -99;
      float pt  = l.pt();
      float eta = l.eta();

      
      BDT = r->Get_MVA_value_two_eta_bins(pt, eta, glb_valid_mu_hits_, tk_valid_pixel_hits_, tk_valid_hits_, glb_chi2_, PFPhotonIso, PFChargedHadIso, PFNeutralHadIso, rho);
//       cout << BDT << endl;
      
      
      bool isBDT = false;
      // New WPs at 96-99.7% of efficiency in low-high pT respectively 
      if ( setup == 2016 )
      {  
         isBDT = ((pt < 10 && abs(eta) < 1.2 && BDT > -0.798) || (pt < 10 && abs(eta) >= 1.2 && BDT > 0.226)
              || (pt >= 10 && abs(eta) < 1.2 && BDT > -3.743) || (pt >= 10 && abs(eta) >= 1.2 && BDT > -3.522));
      }
      else if ( setup == 2017 )
      {
         isBDT = ((pt < 10 && abs(eta) < 1.2 && BDT > -0.578) || (pt < 10 && abs(eta) >= 1.2 && BDT > 0.252)
              || (pt >= 10 && abs(eta) < 1.2 && BDT > -3.688)  || (pt >= 10 && abs(eta) >= 1.2 && BDT > -3.970));
      }
      else if ( setup == 2018 )
      {
         isBDT = ((pt < 10 && abs(eta) < 1.2 && BDT > -0.643)  || (pt < 10 && abs(eta) >= 1.2 && BDT > 0.197)
              || (pt >= 10 && abs(eta) < 1.2 && BDT > -2.737) || (pt >= 10 && abs(eta) >= 1.2 && BDT > -3.804));
      }
      else
      {
         std::cerr << "[ERROR] MuFiller: no MVA setup for: " << setup << " year!" << std::endl;
      }

//=================
// End MVA reader
//================= 
      
      
      //--- Trigger matching
      bool HLTMatch = false;
      pat::TriggerObjectStandAloneCollection obj= l.triggerObjectMatches();
      for ( size_t iTrigObj = 0; iTrigObj < obj.size(); ++iTrigObj ) {
         obj.at( iTrigObj ).unpackFilterLabels(iEvent,*triggerResults );
      }
      for ( size_t i = 0; i < obj.size(); ++i ) {
         if ( obj.at( i ).hasFilterLabel( "hltSingleMu13L3Filtered17"))
         HLTMatch=true;
         if ( obj.at( i ).hasFilterLabel( "hltSingleMu13L3Filtered13") && obj.at(i).pt()>17)
         HLTMatch=true;
         if ( obj.at( i ).hasFilterLabel( "hltDiMuonL3PreFiltered5") && obj.at(i).pt()>17)
         HLTMatch=true;
         if ( obj.at( i ).hasFilterLabel( "hltDiMuonL3PreFiltered7") && obj.at(i).pt()>17)
         HLTMatch=true;
      }
      
      //--- Muon Timing
      float muontime = 0;
      if (l.time().nDof>4) muontime= l.time().timeAtIpInOut;
      
      
      //--- Embed user variables
      l.addUserFloat("PFChargedHadIso",PFChargedHadIso);
      l.addUserFloat("PFNeutralHadIso",PFNeutralHadIso);
      l.addUserFloat("PFPhotonIso",PFPhotonIso);
      l.addUserFloat("PFPUChargedHadIso",PFPUChargedHadIso);
      l.addUserFloat("combRelIsoPF",combRelIsoPF);
      l.addUserFloat("rho",rho);
      l.addUserFloat("SIP",SIP);
      l.addUserFloat("dxy",dxy);
      l.addUserFloat("dz",dz);
      l.addUserFloat("BDT",BDT);
      l.addUserFloat("isBDT",isBDT);
      l.addUserFloat("HLTMatch", HLTMatch);
      // l.addUserCand("MCMatch",genMatch); // FIXME
      l.addUserFloat("time",muontime);
      
      
      //--- isPFMuon flag - in old samples, l.isPFMuon() is not functional, so this has to be filled
      //    beforehand with the module PATPFMuonEmbedder.
      if(!l.hasUserFloat("isPFMuon")) {
         l.addUserFloat("isPFMuon",l.isPFMuon());
      }
      
      // isHighPtTrackerMuon - used in 2016 tight muon ID
      bool isTrackerHighPt = ( l.numberOfMatchedStations() > 1
                              && (l.muonBestTrack()->ptError()/l.muonBestTrack()->pt()) < 0.3
                              && dxy < 0.2
                              && dz < 0.5
                              && l.innerTrack()->hitPattern().numberOfValidPixelHits() > 0
                              && l.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 );
      l.addUserFloat("isTrackerHighPtMuon",isTrackerHighPt);
      
      
      if (!l.hasUserFloat("correctedPtError")) {
         l.addUserFloat("correctedPtError",l.muonBestTrack()->ptError()); //This is expected by the kin fitter
      }
      
      //--- MC parent code
      //     MCHistoryTools mch(iEvent);
      //     if (mch.isMC()) {
      //       int MCParentCode = 0;//FIXME: does not work on cmg mch.getParentCode((l.genParticleRef()).get());
      //       l.addUserFloat("MCParentCode",MCParentCode);
      //     }
      
      //--- Check selection cut. Being done here, flags are not available; but this way we
      //    avoid wasting time on rejected leptons.
      if (!cut(l)) continue;
      
      //--- Embed flags (ie flags specified in the "flags" pset)
      for(CutSet<pat::Muon>::const_iterator flag = flags.begin(); flag != flags.end(); ++flag) {
         l.addUserFloat(flag->first,int((*(flag->second))(l)));
      }
      
      result->push_back(l);
   }
   iEvent.put(std::move(result));
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(MuFiller);

