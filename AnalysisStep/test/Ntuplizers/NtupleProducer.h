#ifndef NtupleProducer_H
#define NtupleProducer_H

#include <memory>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <DataFormats/HepMCCandidate/interface/GenParticleFwd.h>
#include <FWCore/Framework/interface/ESHandle.h>

#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TVector3.h"

#include <DataFormats/PatCandidates/interface/Isolation.h>

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

class EgammaTowerIsolation ;
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include <vector>


/*
 namespace edm {
 class ParameterSet;
 class Event;
 class EventSetup;
 class LuminosityBlock;
 }
 */

using namespace std;

namespace pat {
	class Muon;
	class Electron;
}

namespace reco {
  class GenParticle;
}


class NtupleProducer: public edm::EDAnalyzer {
public:
	
	//typedef std::map<pat::IsolationKeys, float> IsolationMap;
	//typedef std::map<const pat::Muon*, IsolationMap> CandIsolationMap;
	typedef std::map<const pat::Muon*, float> CandIsolationMap;
	typedef math::XYZTLorentzVector LorentzVector ;
	
	/// Constructor
	NtupleProducer(const edm::ParameterSet&);
	
	/// Destructor
	virtual ~NtupleProducer();
	
	// Operations
	void endLuminosityBlock(const edm::LuminosityBlock &, const edm::EventSetup&);
	void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void beginJob() ;
	virtual void endJob() ;
	
	mutable CandIsolationMap candiso;
	mutable CandIsolationMap candisoPF;
    CandIsolationMap computeIso(const std::vector<pat::Muon>& muons, const vector<pat::Electron>& electrons);	
	
	void FillEvent (const edm::Event&, const edm::EventSetup&);
	void FillVertices (const edm::Event&, const edm::EventSetup&);
	void FillMuons (const edm::Event&, const edm::EventSetup&);
	void FillElectrons(const edm::Event&, const edm::EventSetup&);
	void FillSC(const edm::Event&, const edm::EventSetup&);
	void FillPhotons(const edm::Event&, const edm::EventSetup&);
	void FillJets(const edm::Event&, const edm::EventSetup&);
	void FillTruth(const edm::Event&, const edm::EventSetup&);
	void FillMET (const edm::Event& iEvent, const edm::EventSetup& iSetup);
	void setMomentum (TLorentzVector & myvector , const LorentzVector & mom) ;
	void Init();

	bool isHLTMatch(const pat::Muon* mu);
	bool isHLTMatch1(const pat::Muon* mu, bool isMC, int irun);
	bool isHLTMatch2(const pat::Muon* mu, bool isMC, int irun);
	bool isHLTMatch3(const pat::Muon* mu, bool isMC, int irun);

	bool isFSR(const reco::Candidate* genLep);
	const reco::Candidate* getParent(const reco::Candidate* genLep);
	

private:
	
	//inputTag
	edm::InputTag EleTag_;
	edm::InputTag MuonTag_;
	edm::InputTag JetTag_;
	edm::InputTag PhotonTag_;
	edm::InputTag VerticesTag_;
	// Trigger Stuff
	edm::InputTag HLTTag_; 
	edm::InputTag triggerEventTag_;
	edm::InputTag MCTag_ ;
	bool isMC_;	
	int lepton_setup;
	
	//edm::InputTag MuRhoCorrection_;
	//edm::InputTag EleRhoCorrection_;
	//edm::InputTag SigmaRhoCorrection_;
	edm::InputTag PileupSrc_;	
		
	std::vector<edm::InputTag > HLT_Filters_;
	edm::InputTag SCTag_;

	//tree
	TTree *mytree_;
	TLorentzVector myvector ;    
	
	//counters
	unsigned int Nevt_Gen;
	unsigned int Nevt_H4lFilter;
	unsigned int Nevt_afterCleaning;
	unsigned int Nevt_afterSkim;
	
	// For Custom HCAL
	edm::Handle<CaloTowerCollection> * m_calotowers; //towersH_ ;
	edm::InputTag hcalTowers_ ;
	EgammaTowerIsolation * hadDepth1Isolation03_; //towerIso1_ ;
	EgammaTowerIsolation * hadDepth2Isolation03_; //:towerIso2_ ;
	double hOverEConeSize_ ;
	double hOverEPtMin_ ;            

	//global variables
	int _nEvent, _nRun, _nLumi;
	//pile-up
	int _PU_N;
	double _PU_Murho,_PU_Elerho, _PU_sigma, _PU_beta ;  //corrections from FastJets
	
	//trigger fired names
	char trig_fired_names[10000];

	//vertices
	int _vtx_N;
	double _vtx_x[50], _vtx_y[50], _vtx_z[50];
	double _vtx_normalizedChi2[50], _vtx_ndof[50], _vtx_nTracks[50], _vtx_d0[50];
	GlobalPoint vertexPosition ;

	// MET
	double _met_calo_et,_met_calo_px, _met_calo_py, _met_calo_phi, _met_calo_set, _met_calo_sig; 
	double _met_calomu_et,_met_calomu_px, _met_calomu_py, _met_calomu_phi, _met_calomu_set, _met_calomu_sig; 
	double _met_tc_et,_met_tc_px, _met_tc_py, _met_tc_phi, _met_tc_set, _met_tc_sig; 
	double _met_pf_et,_met_pf_px, _met_pf_py, _met_pf_phi, _met_pf_set, _met_pf_sig; 

	//muons
	int _muons_N;
	TClonesArray * m_muons;
	int _muons_charge[50];
	int _muons_isPF[50];
	int _muons_istracker[50], _muons_isstandalone[50], _muons_isglobal[50];
	double _muons_dxy[50], _muons_dz[50], _muons_dxyPV[50], _muons_dzPV[50], _muons_normalizedChi2[50];
	int  _muons_NtrackerHits[50], _muons_NpixelHits[50], _muons_NmuonHits[50], _muons_Nmatches[50];
	int _muons_nTkIsoR03[50], _muons_nTkIsoR05[50];
	double _muons_tkIsoR03[50],_muons_tkIsoR05[50],_muons_emIsoR03[50],_muons_emIsoR05[50],_muons_hadIsoR03[50],_muons_hadIsoR05[50];
	double _muons_trkDxy[50], _muons_trkDxyError[50], _muons_trkDxyB[50],
	_muons_trkDz[50], _muons_trkDzError[50], _muons_trkDzB[50], _muons_trkChi2PerNdof[50], 
	_muons_trkCharge[50],_muons_trkNHits[50],_muons_trkNPixHits[50];
	double _muons_trkmuArbitration[50], _muons_trkmu2DCompatibilityLoose[50], 
	_muons_trkmu2DCompatibilityTight[50], _muons_trkmuOneStationLoose[50], 
	_muons_trkmuOneStationTight[50], _muons_trkmuLastStationLoose[50],
	_muons_trkmuLastStationTight[50], _muons_trkmuOneStationAngLoose[50], 
	_muons_trkmuOneStationAngTight[50], _muons_trkmuLastStationAngLoose[50],
	_muons_trkmuLastStationAngTight[50], _muons_trkmuLastStationOptimizedLowPtLoose[50], 
	_muons_trkmuLastStationOptimizedLowPtTight[50];
	double _muons_caloCompatibility[50], _muons_segmentCompatibility[50], _muons_glbmuPromptTight[50] ; 
	double _muons_HZZisoTk[50], _muons_HZZisoTk5[50],_muons_HZZisoEcal[50], _muons_HZZisoHcal[50], _muons_HZZisoComb[50] ;
	double _muons_FsrIsoEcalDr005[50], _muons_FsrIsoEcalDr007[50],_muons_FsrIsoEcalDr010[50];
	double _muons_isoecal[50], _muons_isohcal[50];
	double _muons_pfChargedHadIso[50], _muons_pfNeutralHadIso[50], _muons_pfPhotonIso[50], _muons_pfChargedHadPUIso[50],_muons_pfCombRelIso[50];
	double _muons_IP[50], _muons_IPError[50], _muons_SIP[50] ;
	int _muons_isHLTMatch[50], _muons_isHLTMatch8[50], _muons_isHLTMatch13[50], _muons_isHLTMatch17[50] ;
	
	double _muons_Err_00[50] ;
	double _muons_Err_01[50] ;
	double _muons_Err_02[50] ;
	double _muons_Err_10[50] ;
	double _muons_Err_11[50] ;
	double _muons_Err_12[50] ;
	double _muons_Err_20[50] ;
	double _muons_Err_21[50] ;
	double _muons_Err_22[50] ;
	//electrons
	int ele_N;
	TClonesArray * m_electrons;
	int ele_echarge[50];
	double ele_he[50] , ele_eseedpout[50] , ele_ep[50] , ele_eseedp[50] , ele_eelepout[50] ;       
	double ele_deltaetaseed[50] , ele_deltaetaele[50] , ele_deltaphiseed[50] , ele_deltaphiele[50] , ele_deltaetain[50] , ele_deltaphiin[50] ;
	double ele_sigmaietaieta[50] , ele_sigmaetaeta[50] , ele_e15[50] , ele_e25max[50] , ele_e55[50] , ele_e1[50] ; //, ele_e33[50] , ele_e2overe9[50] ;
	double ele_pin_mode[50] , ele_pout_mode[50] , ele_pTin_mode[50] , ele_pTout_mode[50] ; 
	double ele_fbrem[50],ele_SCfbrem[50],ele_pfSCfbrem[50], ele_mva[50] ;
	int ele_isbarrel[50] , ele_isendcap[50] , 
	ele_isEBetaGap[50] , ele_isEBphiGap[50] , ele_isEEdeeGap[50] , ele_isEEringGap[50] ,
	ele_isecalDriven[50] , ele_istrackerDriven[50] ,
	ele_eClass[50];
	int ele_missing_hits[50], ele_lost_hits[50]; 
	double ele_chi2_hits[50]; 
	double ele_dxyB[50], ele_dxy[50], ele_dzB[50], ele_dz[50], ele_dszB[50], ele_dsz[50];              
	double ele_tkSumPt_dr03[50] , ele_ecalRecHitSumEt_dr03[50] , ele_hcalDepth1TowerSumEt_dr03[50] , ele_hcalDepth2TowerSumEt_dr03[50] ,
	ele_tkSumPt_dr04[50] , ele_ecalRecHitSumEt_dr04[50] , ele_hcalDepth1TowerSumEt_dr04[50] , ele_hcalDepth2TowerSumEt_dr04[50] ;
	double ele_conv_dcot[50];
	double ele_conv_dist[50];
	int ele_expected_inner_hits[50];
	double ele_eidVeryLoose[50], ele_eidLoose[50], ele_eidMedium[50], ele_eidTight[50]; 
	double ele_eidHZZVeryLoose[50], ele_eidHZZLoose[50], ele_eidHZZMedium[50], ele_eidHZZTight[50], ele_eidHZZSuperTight[50]; 
	double ele_eidMVATrig[50], ele_eidMVANoTrig[50];
	double ele_HZZisoTk[50],ele_HZZisoTk5[50],  ele_HZZisoEcal[50], ele_HZZisoHcal[50],ele_HZZisoComb[50] ;
	double ele_pfChargedHadIso[50], ele_pfNeutralHadIso[50], ele_pfPhotonIso[50], ele_pfChargedHadPUIso[50],ele_pfCombRelIso[50] ;
	double ele_IP[50], ele_IPError[50], ele_SIP[50] ;
	double ele_sclE[50], ele_sclEt[50], ele_sclEta[50], ele_sclPhi[50], ele_sclRawE[50];
	int ele_sclNclus[50];
	double ele_ecalErr[50], ele_trackErr[50], ele_combErr[50], ele_PFcombErr[50];
	double ele_ecalRegressionEnergy[50],ele_ecalRegressionError[50], ele_ecalTrackRegressionEnergy[50],ele_ecalTrackRegressionError[50],ele_ecalScale[50],
	      ele_ecalSmear[50],ele_ecalRegressionScale[50],ele_ecalRegressionSmear[50],
	      ele_ecalTrackRegressionScale[50],ele_ecalTrackRegressionSmear[50];

	
	
	double ele_HCALFullConeSum[50];
	double ele_sclphiwidth[50];

	double ele_mvafbrem[50], ele_mvadetain[50], ele_mvadphiin[50], ele_mvasieie[50], ele_mvahoe[50], ele_mvaeop[50], 
	  ele_mvae1x5e5x5[50], ele_mvaeleopout[50], ele_mvakfchi2[50], ele_mvadist[50],  ele_mvadcot[50], ele_mvaeta[50],
	  ele_mvapt[50];
	int ele_mvakfhits[50], ele_mvamishits[50], ele_mvaecalseed	[50];

	//Photons
	int photon_N;
	TClonesArray * m_photons;
	int _pho_isFromMuon[50]; 
	int _pho_isEB[50], _pho_isEE[50];
	double _pho_sigmaietaieta [50], _pho_he [50], _pho_r9[50],  
	  _pho_TkIso03[50], _pho_HCTkIso03 [50], _pho_emIso03 [50],_pho_hadIso03[50];

	// SuperClusters
	int _sc_N; 
	double _sc_E[50], _sc_Et[50], _sc_Eta[50], _sc_Phi[50], _sc_TkIso[50]; 

	// L1 EM candidates
	int _trig_L1emIso_N; 
	TClonesArray * m_L1emIso;
	int _trig_L1emNonIso_N; 
	TClonesArray * m_L1emNonIso;

	// HLT	
	int _trig_HLT_N; 
	double _trig_HLT_eta[500], _trig_HLT_phi[500], _trig_HLT_energy[500], _trig_HLT_pt[500];
	std::vector<std::string> *  _trig_HLT_name;
	//TClonesArray * m_HLT;

	// Calo Jets 
	//TClonesArray * _m_jets_calo ;
	//int _jets_calo_N;
	//double _jets_calo_E[50], _jets_calo_pT[50], _jets_calo_px[50], _jets_calo_py[50], _jets_calo_pz[50], _jets_calo_eta[50], _jets_calo_phi[50];

	// PF Jets
	TClonesArray * _m_jets_pf;
	int _jets_pf_N;
	
	double jets_pf_chargedHadEFrac[50], jets_pf_chargedEmEFrac[50], jets_pf_chargedMuEFrac[50]; 
	double jets_pf_neutralHadEFrac[50], jets_pf_neutralEmEFrac[50], jets_pf_PhotonEFrac[50];
	  
	int jets_pf_chargedHadMultiplicity[50], jets_pf_neutralHadMultiplicity[50];
	int jets_pf_chargedMultiplicity[50], jets_pf_neutralMultiplicity[50];
	int jets_pf_nConstituents[50];

	//MC
	TClonesArray * _m_MC_gen_V;
	TClonesArray * _m_MC_gen_Higgs;
	TClonesArray * _m_MC_gen_photons;
	TClonesArray * _m_MC_gen_leptons;
	TClonesArray * _m_MC_gen_leptons_status1;
	TClonesArray * _m_MC_gen_leptons_status2;
	double _MC_gen_V_pdgid[10];
	double _MC_gen_leptons_pdgid[10];
	double _MC_gen_leptons_status1_pdgid[10];
	double _MC_gen_leptons_status2_pdgid[10];
	double _MC_pthat;
	int _MC_flavor[2];

	int _MC_gen_photons_isFSR[5000];

	//TClonesArray &MC_gen_leptons_status2 = *_m_MC_gen_leptons_status2;
	//TClonesArray &MC_gen_leptons_status1 = *_m_MC_gen_leptons_status1;
	
	bool fill_L1trigger;
	bool fill_SC;

	edm::ESHandle<MagneticField> magfield;
};
#endif


