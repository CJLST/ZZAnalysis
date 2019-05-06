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

#include <MuonMVAReader/Reader/interface/MuonMVAReader.hpp>

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
      MuonMVAReader *r;
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
   if (setup == 2016)
   {
      r = new MuonMVAReader("../../../MuonMVAReader/Reader/data/MVAWeightFiles/Summer_16_ID_ISO_SIP/");
      r->Initialize();
   }
   else if (setup == 2017)
   {
      r = new MuonMVAReader("../../../MuonMVAReader/Reader/data/MVAWeightFiles/Summer_16_ID_ISO_SIP/");
      r->Initialize();
   }
   else if (setup == 2018)
   {
      r = new MuonMVAReader("../../../MuonMVAReader/Reader/data/MVAWeightFiles/Summer_16_ID_ISO_SIP/");
      r->Initialize();
   }
   else
   {
      throw cms::Exception("MuonMVA") << "Muon MVA is not defined for the given setup (" << setup << ")!";
   }

}
MuFiller::~MuFiller(){
}


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

      
      // MVA Reader begin
      float mu_N_hits_, mu_chi_square_, mu_N_pixel_hits_, mu_N_tracker_hits_;
      
      bool is_global_mu_  = l.isGlobalMuon();
      if ( is_global_mu_ )
      {
         // Number of muon chamber hits included in the the global muon track fit
         mu_N_hits_ = (l.globalTrack()->hitPattern().numberOfValidMuonHits());
         
         // Chi2 of the global track fit
         mu_chi_square_ = (l.globalTrack()->normalizedChi2());
      }
      else
      {
         mu_N_hits_     = -1;
         mu_chi_square_ = -1;
      }
      
      
      // Number of hits in the pixel detector
      bool valid_KF = false;
      reco::TrackRef myTrackRef = l.innerTrack();
      valid_KF = (myTrackRef.isAvailable());
      valid_KF = (myTrackRef.isNonnull());
      
      
      if ( valid_KF )
      {
         // Number of pixel hits
         mu_N_pixel_hits_ = l.innerTrack()->hitPattern().numberOfValidPixelHits();
         
         // Number of hits in the tracker layers
         mu_N_tracker_hits_ = l.innerTrack()->hitPattern().trackerLayersWithMeasurement();
      }
      else
      {
         mu_N_pixel_hits_ = -1;
         mu_N_tracker_hits_ = -1;
      }
      
      float BDT = -99;
      BDT = r->Get_MVA_value(l.pt(), l.eta(), mu_N_hits_, mu_N_pixel_hits_, mu_N_tracker_hits_, mu_chi_square_, PFPhotonIso, PFChargedHadIso, PFNeutralHadIso, rho, SIP, dxy, dz);
      
      cout << BDT << endl;
      
      bool isBDT = false;
      // MVA Reader end
      
      
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

