
#include "NtupleProducer.h"

#include <FWCore/Framework/interface/MakerMacros.h>
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/LuminosityBlock.h>
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <DataFormats/Common/interface/MergeableCounter.h>
#include <DataFormats/Common/interface/View.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>
//
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
//"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
// SuperCluster
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"

// Trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

// MET
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

// PF Jets
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
// Calo Jets
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"


//MC
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/CartesianTrajectoryError.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianCurvilinearToCartesian.h"
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>

#include <cmath>

#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>



using namespace std;
using namespace edm;
using namespace reco;

#include <algorithm>

/// Constructor
NtupleProducer::NtupleProducer(const ParameterSet& iConfig) :
EleTag_ (iConfig.getParameter<edm::InputTag> ("EleTag")),
MuonTag_ (iConfig.getParameter<edm::InputTag> ("MuonTag")),
JetTag_ (iConfig.getParameter<edm::InputTag> ("JetTag")),
PhotonTag_ (iConfig.getParameter<edm::InputTag> ("PhotonTag")),
VerticesTag_(iConfig.getParameter<edm::InputTag> ("VerticesTag")),
HLTTag_(iConfig.getParameter<edm::InputTag> ("HLTTag")),
triggerEventTag_(iConfig.getParameter<edm::InputTag> ("TriggerEventTag")),
MCTag_(iConfig.getParameter<edm::InputTag> ("MCTag")),
isMC_ (iConfig.getParameter<bool>("isMC")),
lepton_setup(iConfig.getParameter<int>("lepton_setup")),

//RhoCorrection_("kt6PFJetsForIso:rho"),//2011
//MuRhoCorrection_("kt6PFJetsCentralNeutral:rho"),//2012
//EleRhoCorrection_("kt6PFJets:rho"),//2012
//SigmaRhoCorrection_("kt6PFJetsForIso:sigma"),
PileupSrc_ ("addPileupInfo")
{ 
	//cout << " filter ?" << endl;
	HLT_Filters_   = iConfig.getParameter<std::vector<edm::InputTag > >("HLTFilters");
	//cout << " end filter ?" << endl;
	fill_L1trigger=false;
	fill_SC=false;
	if (fill_SC) SCTag_ = iConfig.getParameter<edm::InputTag> ("SCTag");

}

/// Destructor
NtupleProducer::~NtupleProducer(){
	
        delete m_electrons ;
	if (fill_L1trigger) {
	  delete m_L1emIso;
	  delete m_L1emNonIso;
	}
	delete m_muons;
	delete _m_jets_pf;
	delete m_photons;
	
	if(isMC_ ) {
		delete _m_MC_gen_V;
		delete _m_MC_gen_photons;
		delete _m_MC_gen_leptons;
		delete _m_MC_gen_Higgs;
		delete _m_MC_gen_leptons_status1;
		delete _m_MC_gen_leptons_status2;
	} // if MC
	
}

void NtupleProducer::beginJob(){
	
	// Book histograms
	edm::Service<TFileService> fs ;
	mytree_  = fs->make <TTree>("simpleRootTree","simpleRootTree"); 
	
	//// Counters
	//mytree_->Branch("Nevt_Gen",&Nevt_Gen,"Nevt_Gen/I");
	//mytree_->Branch("Nevt_Skim",&Nevt_afterSkim,"Nevt_Skim/I");
	
	// Global
	mytree_->Branch("nEvent",&_nEvent,"nEvent/I");
	mytree_->Branch("nRun",&_nRun,"nRun/I");
	mytree_->Branch("nLumi",&_nLumi,"nLumi/I");
	
	// Pile UP
	mytree_->Branch("PU_N",&_PU_N,"PU_N/I");
	mytree_->Branch("PU_MurhoCorr",&_PU_Murho,"PU_MurhoCorr/D");
	mytree_->Branch("PU_ElerhoCorr",&_PU_Elerho,"PU_ElerhoCorr/D");
	mytree_->Branch("PU_sigmaCorr",&_PU_sigma,"PU_sigmaCorr/D");
	mytree_->Branch("PU_beta",&_PU_beta,"PU_beta/D");
	
	// Trigger
	mytree_->Branch("trig_fired_names",&trig_fired_names,"trig_fired_names[10000]/C");
	
	// Trigger: L1 EM objects
	if (fill_L1trigger) {
	  mytree_->Branch("trig_L1emIso_N",&_trig_L1emIso_N,"trig_L1emIso_N/I");
	  m_L1emIso = new TClonesArray ("TLorentzVector");
	  mytree_->Branch ("trig_L1emIso", "TClonesArray", &m_L1emIso, 256000,0);
	  mytree_->Branch("trig_L1emNonIso_N",&_trig_L1emNonIso_N,"trig_L1emNonIso_N/I");
	  m_L1emNonIso = new TClonesArray ("TLorentzVector");
	  mytree_->Branch ("trig_L1emsNonIso", "TClonesArray", &m_L1emNonIso, 256000,0);
	}
	// Triger: HLT objects
	//cout << "HLT" << endl;
	mytree_->Branch("trig_HLT_N",      &_trig_HLT_N,     "trig_HLT_N/I");
	//m_HLT = new TClonesArray ("TLorentzVector");
	//mytree_->Branch ("trig_HLT", "TClonesArray", &m_HLT, 256000,0);
	mytree_->Branch("trig_HLT_eta",    &_trig_HLT_eta,   "trig_HLT_eta[500]/D");
	mytree_->Branch("trig_HLT_phi",    &_trig_HLT_phi,   "trig_HLT_phi[500]/D");
	mytree_->Branch("trig_HLT_energy", &_trig_HLT_energy,"trig_HLT_energy[500]/D");
	mytree_->Branch("trig_HLT_pt",     &_trig_HLT_pt,    "trig_HLT_pt[500]/D");
	mytree_->Branch("trig_HLT_name",  "std::vector<std::string>",&_trig_HLT_name);

	//cout << "end HLT" << endl;
	
	// Vertices
	mytree_->Branch("vtx_N",&_vtx_N,"vtx_N/I");
	mytree_->Branch("vtx_normalizedChi2",&_vtx_normalizedChi2,"vtx_normalizedChi2[50]/D");
	mytree_->Branch("vtx_ndof",&_vtx_ndof,"vtx_ndof[50]/D");
	mytree_->Branch("vtx_nTracks",&_vtx_nTracks,"vtx_nTracks[50]/D");
	mytree_->Branch("vtx_d0",&_vtx_d0,"vtx_d0[50]/D");
	mytree_->Branch("vtx_x",&_vtx_x,"vtx_x[50]/D");
	mytree_->Branch("vtx_y",&_vtx_y,"vtx_y[50]/D");
	mytree_->Branch("vtx_z",&_vtx_z,"vtx_z[50]/D");
	
	// Muons
	mytree_->Branch("muons_N",&_muons_N,"muons_N/I");
	m_muons = new TClonesArray ("TLorentzVector");
	mytree_->Branch("muons", "TClonesArray", &m_muons, 256000,0);
	mytree_->Branch("muons_charge",&_muons_charge,"muons_charge[50]/I");
	mytree_->Branch("muons_isPF",&_muons_isPF,"muons_isPF[50]/I");
	mytree_->Branch("muons_istracker",&_muons_istracker,"muons_istracker[50]/I");
	mytree_->Branch("muons_isstandalone",&_muons_isstandalone,"muons_isstandalone[50]/I");
	mytree_->Branch("muons_isglobal",&_muons_isglobal,"muons_isglobal[50]/I");
	//
	mytree_->Branch("muons_dxy",&_muons_dxy,"muons_dxy[50]/D");
	mytree_->Branch("muons_dz",&_muons_dz,"muons_dz[50]/D");
	mytree_->Branch("muons_dxyPV",&_muons_dxyPV,"muons_dxyPV[50]/D");
	mytree_->Branch("muons_dzPV",&_muons_dzPV,"muons_dzPV[50]/D");
	mytree_->Branch("muons_normalizedChi2",&_muons_normalizedChi2,"muons_normalizedChi2[50]/D");
	mytree_->Branch("muons_NmuonHits",&_muons_NmuonHits,"muons_NmuonHits[50]/I");
	mytree_->Branch("muons_NtrackerHits",&_muons_NtrackerHits,"muons_NtrackerHits[50]/I");
	mytree_->Branch("muons_NpixelHits",&_muons_NpixelHits,"muons_NpixelHits[50]/I");
	mytree_->Branch("muons_Nmatches",&_muons_Nmatches,"muons_Nmatches[50]/I");
	//
	mytree_->Branch("muons_nTkIsoR03",&_muons_nTkIsoR03,"muons_nTkIsoR03[50]/I");
	mytree_->Branch("muons_nTkIsoR05",&_muons_nTkIsoR05,"muons_nTkIsoR05[50]/I");
	mytree_->Branch("muons_tkIsoR03",&_muons_tkIsoR03,"muons_tkIsoR03[50]/D");
	mytree_->Branch("muons_tkIsoR05",&_muons_tkIsoR05,"muons_tkIsoR05[50]/D");
	mytree_->Branch("muons_emIsoR03",&_muons_emIsoR03,"muons_emIsoR03[50]/D");
	mytree_->Branch("muons_emIsoR05",&_muons_emIsoR05,"muons_emIsoR05[50]/D");
	mytree_->Branch("muons_hadIsoR03",&_muons_hadIsoR03,"muons_hadIsoR03[50]/D");
	mytree_->Branch("muons_hadIsoR05",&_muons_hadIsoR05,"muons_hadIsoR05[50]/D");
	//
	mytree_->Branch("muons_trkDxy",&_muons_trkDxy,"muons_trkDxy[50]/D");
	mytree_->Branch("muons_trkDxyError",&_muons_trkDxyError,"muons_trkDxyError[50]/D");
	mytree_->Branch("muons_trkDxyB",&_muons_trkDxyB,"muons_trkDxyB[50]/D");
	mytree_->Branch("muons_trkDz",&_muons_trkDz,"muons_trkDz[50]/D");
	mytree_->Branch("muons_trkDzError",&_muons_trkDzError,"muons_trkDzError[50]/D");
	mytree_->Branch("muons_trkDzB",&_muons_trkDzB,"muons_trkDzB[50]/D"); 
	mytree_->Branch("muons_trkChi2PerNdof",&_muons_trkChi2PerNdof,"muons_trkChi2PerNdof[50]/D");
	mytree_->Branch("muons_trkCharge",&_muons_trkCharge,"muons_trkCharge[50]/D");
	mytree_->Branch("muons_trkNHits",&_muons_trkNHits,"muons_trkNHits[50]/D");
	mytree_->Branch("muons_trkNPixHits",&_muons_trkNPixHits,"muons_trkNPixHits[50]/D");
	mytree_->Branch("muons_trkmuArbitration",&_muons_trkmuArbitration,"muons_trkmuArbitration[50]/D");
	mytree_->Branch("muons_trkmu2DCompatibilityLoose",&_muons_trkmu2DCompatibilityLoose,"muons_trkmu2DCompatibilityLoose[50]/D");
	mytree_->Branch("muons_trkmu2DCompatibilityTight",&_muons_trkmu2DCompatibilityTight,"muons_trkmu2DCompatibilityTight[50]/D");
	mytree_->Branch("muons_trkmuOneStationLoose",&_muons_trkmuOneStationLoose,"muons_trkmuOneStationLoose[50]/D");
	mytree_->Branch("muons_trkmuOneStationTight",&_muons_trkmuOneStationTight,"muons_trkmuOneStationTight[50]/D");
	mytree_->Branch("muons_trkmuLastStationLoose",&_muons_trkmuLastStationLoose,"muons_trkmuLastStationLoose[50]/D");
	mytree_->Branch("muons_trkmuLastStationTight",&_muons_trkmuLastStationTight,"muons_trkmuLastStationTight[50]/D");
	mytree_->Branch("muons_trkmuOneStationAngLoose",&_muons_trkmuOneStationAngLoose,"muons_trkmuOneStationAngLoose[50]/D");
	mytree_->Branch("muons_trkmuOneStationAngTight",&_muons_trkmuOneStationAngTight,"muons_trkmuOneStationAngTight[50]/D");
	mytree_->Branch("muons_trkmuLastStationAngLoose",&_muons_trkmuLastStationAngLoose,"muons_trkmuLastStationAngLoose[50]/D");
	mytree_->Branch("muons_trkmuLastStationAngTight",&_muons_trkmuLastStationAngTight,"muons_trkmuLastStationAngTight[50]/D");
	mytree_->Branch("muons_trkmuLastStationOptimizedLowPtLoose",&_muons_trkmuLastStationOptimizedLowPtLoose,"muons_trkmuLastStationOptimizedLowPtLoose[50]/D");
	mytree_->Branch("muons_trkmuLastStationOptimizedLowPtTight",&_muons_trkmuLastStationOptimizedLowPtTight,"muons_trkmuLastStationOptimizedLowPtTight[50]/D");
	mytree_->Branch("muons_HZZisoTk",&_muons_HZZisoTk,"muons_HZZisoTk[50]/D");
	mytree_->Branch("muons_HZZisoTk5",&_muons_HZZisoTk5,"muons_HZZisoTk5[50]/D");
	mytree_->Branch("muons_HZZisoEcal",&_muons_HZZisoEcal,"muons_HZZisoEcal[50]/D");
	mytree_->Branch("muons_HZZisoHcal",&_muons_HZZisoHcal,"muons_HZZisoHcal[50]/D");
	mytree_->Branch("muons_HZZisoComb",&_muons_HZZisoComb,"muons_HZZisoComb[50]/D");
	mytree_->Branch("muons_FsrIsoEcalDr005",&_muons_FsrIsoEcalDr005,"muons_FsrIsoEcalDr005[50]/D");
	mytree_->Branch("muons_FsrIsoEcalDr007",&_muons_FsrIsoEcalDr007,"muons_FsrIsoEcalDr007[50]/D");
	mytree_->Branch("muons_FsrIsoEcalDr010",&_muons_FsrIsoEcalDr010,"muons_FsrIsoEcalDr010[50]/D");
	
	mytree_->Branch("muons_pfChargedHadIso",&_muons_pfChargedHadIso,"muons_pfChargedHadIso[50]/D");
	mytree_->Branch("muons_pfNeutralHadIso",&_muons_pfNeutralHadIso,"muons_pfNeutralHadIso[50]/D");
	mytree_->Branch("muons_pfPhotonIso",&_muons_pfPhotonIso,"muons_pfPhotonIso[50]/D");
	mytree_->Branch("muons_pfChargedHadPUIso",&_muons_pfChargedHadPUIso,"muons_pfChargedHadPUIso[50]/D");
	mytree_->Branch("muons_pfCombRelIso",&_muons_pfCombRelIso,"muons_pfCombRelIso[50]/D");
	// 
	mytree_->Branch("muons_IP",&_muons_IP,"muons_IP[50]/D");
	mytree_->Branch("muons_IPError",&_muons_IPError,"muons_IPError[50]/D");	
	mytree_->Branch("muons_SIP",&_muons_SIP,"muons_SIP[50]/D");
	//
	mytree_->Branch ("muons_isHLTMatch",  &_muons_isHLTMatch,  "muons_isHLTMatch[50]/I");
	mytree_->Branch ("muons_isHLTMatch8", &_muons_isHLTMatch8, "muons_isHLTMatch8[50]/I");
	mytree_->Branch ("muons_isHLTMatch13",&_muons_isHLTMatch13,"muons_isHLTMatch13[50]/I");
	mytree_->Branch ("muons_isHLTMatch17",&_muons_isHLTMatch17,"muons_isHLTMatch17[50]/I");
	
	// MET
	mytree_->Branch("met_calo_et",&_met_calo_et,"met_calo_et/D");
	mytree_->Branch("met_calo_px",&_met_calo_px,"met_calo_px/D");
	mytree_->Branch("met_calo_py",&_met_calo_py,"met_calo_py/D");
	mytree_->Branch("met_calo_phi",&_met_calo_phi,"met_calo_phi/D");
	mytree_->Branch("met_calo_set",&_met_calo_set,"met_calo_set/D");
	mytree_->Branch("met_calo_sig",&_met_calo_sig,"met_calo_sig/D");
	
	mytree_->Branch("met_calomu_et",&_met_calomu_et,"met_calomu_et/D");
	mytree_->Branch("met_calomu_px",&_met_calomu_px,"met_calomu_px/D");
	mytree_->Branch("met_calomu_py",&_met_calomu_py,"met_calomu_py/D");
	mytree_->Branch("met_calomu_phi",&_met_calomu_phi,"met_calomu_phi/D");
	mytree_->Branch("met_calomu_set",&_met_calomu_set,"met_calomu_set/D");
	mytree_->Branch("met_calomu_sig",&_met_calomu_sig,"met_calomu_sig/D");
	
	mytree_->Branch("met_tc_et",&_met_tc_et,"met_tc_et/D");
	mytree_->Branch("met_tc_px",&_met_tc_px,"met_tc_px/D");
	mytree_->Branch("met_tc_py",&_met_tc_py,"met_tc_py/D");
	mytree_->Branch("met_tc_phi",&_met_tc_phi,"met_tc_phi/D");
	mytree_->Branch("met_tc_set",&_met_tc_set,"met_tc_set/D");
	mytree_->Branch("met_tc_sig",&_met_tc_sig,"met_tc_sig/D");
	
	mytree_->Branch("met_pf_et",&_met_pf_et,"met_pf_et/D");
	mytree_->Branch("met_pf_px",&_met_pf_px,"met_pf_px/D");
	mytree_->Branch("met_pf_py",&_met_pf_py,"met_pf_py/D");
	mytree_->Branch("met_pf_phi",&_met_pf_phi,"met_pf_phi/D");
	mytree_->Branch("met_pf_set",&_met_pf_set,"met_pf_set/D");
	mytree_->Branch("met_pf_sig",&_met_pf_sig,"met_pf_sig/D");
	
	// Electrons
	mytree_->Branch("ele_N",&ele_N,"ele_N/I");
	m_electrons = new TClonesArray ("TLorentzVector");
	mytree_->Branch ("electrons", "TClonesArray", &m_electrons, 256000,0);
	mytree_->Branch("ele_echarge",ele_echarge,"ele_echarge[50]/I");
	//
	mytree_->Branch("ele_he",ele_he,"ele_he[50]/D");
	mytree_->Branch("ele_eseedpout",ele_eseedpout,"ele_eseedpout[50]/D");
	mytree_->Branch("ele_ep",ele_ep,"ele_ep[50]/D");
	mytree_->Branch("ele_eseedp",ele_eseedp,"ele_eseedp[50]/D");
	mytree_->Branch("ele_eelepout",ele_eelepout,"ele_eelepout[50]/D");
	//
	mytree_->Branch("ele_pin_mode",ele_pin_mode,"ele_pin_mode[50]/D");
	mytree_->Branch("ele_pout_mode",ele_pout_mode,"ele_pout_mode[50]/D");
	mytree_->Branch("ele_pTin_mode",ele_pTin_mode,"ele_pTin_mode[50]/D");
	mytree_->Branch("ele_pTout_mode",ele_pTout_mode,"ele_pTout_mode[50]/D");
	//
	mytree_->Branch("ele_deltaetaseed",ele_deltaetaseed,"ele_deltaetaseed[50]/D");
	mytree_->Branch("ele_deltaphiseed",ele_deltaphiseed,"ele_deltaphiseed[50]/D");
	mytree_->Branch("ele_deltaetaele",ele_deltaetaele,"ele_deltaetaele[50]/D");
	mytree_->Branch("ele_deltaphiele",ele_deltaphiele,"ele_deltaphiele[50]/D");
	mytree_->Branch("ele_deltaetain",ele_deltaetain,"ele_deltaetain[50]/D");
	mytree_->Branch("ele_deltaphiin",ele_deltaphiin,"ele_deltaphiin[50]/D");
	//
	mytree_->Branch("ele_sigmaietaieta",ele_sigmaietaieta,"ele_sigmaietaieta[50]/D");
	mytree_->Branch("ele_sigmaetaeta",ele_sigmaetaeta,"ele_sigmaetaeta[50]/D");
	mytree_->Branch("ele_e15",ele_e15,"ele_e15[50]/D");
	mytree_->Branch("ele_e25max",ele_e25max,"ele_e25max[50]/D");
	mytree_->Branch("ele_e55",ele_e55,"ele_e55[50]/D");
	mytree_->Branch("ele_e1",ele_e1,"ele_e1[50]/D");
	//mytree_->Branch("ele_e33",ele_e33,"ele_e33[50]/D");  ---> No ECAL Reduced Collection
	//mytree_->Branch("ele_e2overe9",ele_e2overe9,"ele_e2overe9[50]/D");  ---> No ECAL Reduced Collection
	// 
	mytree_->Branch("ele_fbrem",ele_fbrem,"ele_fbrem[50]/D");
	mytree_->Branch("ele_SCfbrem",ele_SCfbrem,"ele_SCfbrem[50]/D");
	mytree_->Branch("ele_pfSCfbrem",ele_pfSCfbrem,"ele_pfSCfbrem[50]/D");
	mytree_->Branch("ele_mva",ele_mva,"ele_mva[50]/D");
	//
	mytree_->Branch("ele_isbarrel",ele_isbarrel,"ele_isbarrel[50]/I");
	mytree_->Branch("ele_isendcap",ele_isendcap,"ele_isendcap[50]/I");
	mytree_->Branch("ele_isEBetaGap",ele_isEBetaGap,"ele_isEBetaGap[50]/I");
	mytree_->Branch("ele_isEBphiGap",ele_isEBphiGap,"ele_isEBphiGap[50]/I");
	mytree_->Branch("ele_isEEdeeGap",ele_isEEdeeGap,"ele_isEEdeeGap[50]/I");
	mytree_->Branch("ele_isEEringGap",ele_isEEringGap,"ele_isEEringGap[50]/I");
	mytree_->Branch("ele_isecalDriven",ele_isecalDriven,"ele_isecalDriven[50]/I");
	mytree_->Branch("ele_istrackerDriven",ele_istrackerDriven,"ele_istrackerDriven[50]/I");
	mytree_->Branch("ele_eClass",ele_eClass,"ele_eClass[50]/I");
	//
	mytree_->Branch("ele_missing_hits",ele_missing_hits,"ele_missing_hits[50]/I");
	mytree_->Branch("ele_lost_hits",ele_lost_hits,"ele_lost_hits[50]/I");
	mytree_->Branch("ele_chi2_hits",ele_chi2_hits,"ele_chi2_hits[50]/D");	
	//
	mytree_->Branch("ele_dxy",ele_dxy,"ele_dxy[50]/D");
	mytree_->Branch("ele_dxyB",ele_dxyB,"ele_dxyB[50]/D");
	mytree_->Branch("ele_dz",ele_dz,"ele_dz[50]/D");
	mytree_->Branch("ele_dzB",ele_dzB,"ele_dzB[50]/D");
	mytree_->Branch("ele_dsz",ele_dsz,"ele_dsz[50]/D");
	mytree_->Branch("ele_dszB",ele_dszB,"ele_dszB[50]/D");
	//
	mytree_->Branch("ele_tkSumPt_dr03",ele_tkSumPt_dr03,"ele_tkSumPt_dr03[50]/D"); 
	mytree_->Branch("ele_ecalRecHitSumEt_dr03",ele_ecalRecHitSumEt_dr03,"ele_ecalRecHitSumEt_dr03[50]/D"); 
	mytree_->Branch("ele_hcalDepth1TowerSumEt_dr03",ele_hcalDepth1TowerSumEt_dr03,"ele_hcalDepth1TowerSumEt_dr03[50]/D"); 
	mytree_->Branch("ele_hcalDepth2TowerSumEt_dr03",ele_hcalDepth2TowerSumEt_dr03,"ele_hcalDepth2TowerSumEt_dr03[50]/D"); 
	mytree_->Branch("ele_tkSumPt_dr04",ele_tkSumPt_dr04,"ele_tkSumPt_dr04[50]/D"); 
	mytree_->Branch("ele_ecalRecHitSumEt_dr04",ele_ecalRecHitSumEt_dr04,"ele_ecalRecHitSumEt_dr04[50]/D"); 
	mytree_->Branch("ele_hcalDepth1TowerSumEt_dr04",ele_hcalDepth1TowerSumEt_dr04,"ele_hcalDepth1TowerSumEt_dr04[50]/D"); 
	mytree_->Branch("ele_hcalDepth2TowerSumEt_dr04",ele_hcalDepth2TowerSumEt_dr04,"ele_hcalDepth2TowerSumEt_dr04[50]/D"); 
	//
	mytree_->Branch("ele_conv_dist",&ele_conv_dist,"ele_conv_dist[50]/D");
	mytree_->Branch("ele_conv_dcot",&ele_conv_dcot,"ele_conv_dcot[50]/D");
	mytree_->Branch("ele_expected_inner_hits",ele_expected_inner_hits,"ele_expected_inner_hits[50]/I");
	//
	mytree_->Branch("ele_eidVeryLoose",ele_eidVeryLoose,"ele_eidVeryLoose[50]/D");
	mytree_->Branch("ele_eidLoose",ele_eidLoose,"ele_eidLoose[50]/D");
	mytree_->Branch("ele_eidMedium",ele_eidMedium,"ele_eidMedium[50]/D");
	mytree_->Branch("ele_eidTight",ele_eidTight,"ele_eidTight[50]/D"); 
	mytree_->Branch("ele_eidHZZVeryLoose",ele_eidHZZVeryLoose,"ele_eidHZZVeryLoose[50]/D");
	mytree_->Branch("ele_eidHZZLoose",ele_eidHZZLoose,"ele_eidHZZLoose[50]/D");
	mytree_->Branch("ele_eidHZZMedium",ele_eidHZZMedium,"ele_eidHZZMedium[50]/D");
	mytree_->Branch("ele_eidHZZTight",ele_eidHZZTight,"ele_eidHZZTight[50]/D"); 
	mytree_->Branch("ele_eidHZZSuperTight",ele_eidHZZSuperTight,"ele_eidHZZSuperTight[50]/D"); 
	mytree_->Branch("ele_eidMVATrig",ele_eidMVATrig,"ele_eidMVATrig[50]/D"); 
	mytree_->Branch("ele_eidMVANoTrig",ele_eidMVANoTrig,"ele_eidMVANoTrig[50]/D"); 
 	//
	mytree_->Branch("ele_HZZisoTk",ele_HZZisoTk,"ele_HZZisoTk[50]/D");
	mytree_->Branch("ele_HZZisoTk5",ele_HZZisoTk5,"ele_HZZisoTk5[50]/D");
	mytree_->Branch("ele_HZZisoEcal",ele_HZZisoEcal,"ele_HZZisoEcal[50]/D");
	mytree_->Branch("ele_HZZisoHcal",ele_HZZisoHcal,"ele_HZZisoHcal[50]/D");
	mytree_->Branch("ele_HZZisoComb",ele_HZZisoComb,"ele_HZZisoComb[50]/D");
	mytree_->Branch("ele_pfChargedHadIso",&ele_pfChargedHadIso,"ele_pfChargedHadIso[50]/D");
	mytree_->Branch("ele_pfNeutralHadIso",&ele_pfNeutralHadIso,"ele_pfNeutralHadIso[50]/D");
	mytree_->Branch("ele_pfPhotonIso",&ele_pfPhotonIso,"ele_pfPhotonIso[50]/D");
	mytree_->Branch("ele_pfChargedHadPUIso",&ele_pfChargedHadPUIso,"ele_pfChargedHadPUIso[50]/D");
	mytree_->Branch("ele_pfCombRelIso",&ele_pfCombRelIso,"ele_pfCombRelIso[50]/D");
	// 
	mytree_->Branch("ele_IP",ele_IP,"ele_IP[50]/D");
	mytree_->Branch("ele_IPError",ele_IPError,"ele_IPError[50]/D");	
	mytree_->Branch("ele_SIP",ele_SIP,"ele_SIP[50]/D");
	//
	  mytree_->Branch("ele_sclRawE", ele_sclRawE, "ele_sclRawE[50]/D");
	  mytree_->Branch("ele_sclE",    ele_sclE,    "ele_sclE[50]/D");
	  mytree_->Branch("ele_sclEt",   ele_sclEt,   "ele_sclEt[50]/D");
	  mytree_->Branch("ele_sclEta",  ele_sclEta,  "ele_sclEta[50]/D");
	  mytree_->Branch("ele_sclPhi",  ele_sclPhi,  "ele_sclPhi[50]/D");
	  mytree_->Branch("ele_sclNclus",  ele_sclNclus,  "ele_sclNclus[50]/I");
	  mytree_->Branch("ele_sclphiwidth", ele_sclphiwidth, "ele_sclphiwidth[50]/F");
	  //
	mytree_->Branch("ele_ecalErr",   ele_ecalErr,   "ele_ecalErr[50]/D");
	mytree_->Branch("ele_trackErr",  ele_trackErr,  "ele_trackErr[50]/D");
	mytree_->Branch("ele_combErr",   ele_combErr,   "ele_combErr[50]/D");
	mytree_->Branch("ele_PFcombErr", ele_PFcombErr, "ele_PFcombErr[50]/D");
	
	mytree_->Branch("ele_ecalRegressionEnergy",   ele_ecalRegressionEnergy,   "ele_ecalRegressionEnergy[50]/D");
	mytree_->Branch("ele_ecalRegressionError", ele_ecalRegressionError, "ele_ecalRegressionError[50]/D");
	mytree_->Branch("ele_ecalTrackRegressionEnergy",ele_ecalTrackRegressionEnergy,"ele_ecalTrackRegressionEnergy[50]/D");
	mytree_->Branch("ele_ecalTrackRegressionError",ele_ecalTrackRegressionError,"ele_ecalTrackRegressionError[50]/D");
	mytree_->Branch("ele_ecalScale",ele_ecalScale,"ele_ecalScale[50]/D");
	mytree_->Branch("ele_ecalSmear",ele_ecalSmear,"ele_ecalSmear[50]/D");
	mytree_->Branch("ele_ecalRegressionScale",ele_ecalRegressionScale,"ele_ecalRegressionScale[50]/D");
	mytree_->Branch("ele_ecalRegressionSmear",ele_ecalRegressionSmear,"ele_ecalRegressionSmear[50]/D");
	mytree_->Branch("ele_ecalTrackRegressionScale",ele_ecalTrackRegressionScale,"ele_ecalTrackRegressionScale[50]/D");
	mytree_->Branch("ele_ecalTrackRegressionSmear",ele_ecalTrackRegressionSmear,"ele_ecalTrackRegressionSmear[50]/D");

		
		
		
	// Variables for the mva. Most of them are duplicated, but since they are corrected at analysis level, it could be dangerous
	mytree_->Branch("ele_mvafbrem", ele_mvafbrem,"ele_mvafbrem[50]/D");
	mytree_->Branch("ele_mvadetain", ele_mvadetain,"ele_mvadetain[50]/D");
	mytree_->Branch("ele_mvadphiin", ele_mvadphiin,"ele_mvadphiin[50]/D");
	mytree_->Branch("ele_mvasieie", ele_mvasieie,"ele_mvasiesie[50]/D");
	mytree_->Branch("ele_mvahoe", ele_mvahoe,"ele_mvahoe[50]/D");
	mytree_->Branch("ele_mvaeop", ele_mvaeop,"ele_mvaeop[50]/D");
	mytree_->Branch("ele_mvae1x5e5x5", ele_mvae1x5e5x5,"ele_mvae1x5e5x5[50]/D");
	mytree_->Branch("ele_mvaeleopout", ele_mvaeleopout,"ele_mvaeleopout[50]/D");
	mytree_->Branch("ele_mvakfchi2", ele_mvakfchi2,"ele_mvakfchi2[50]/D");
	mytree_->Branch("ele_mvakfhits", ele_mvakfhits,"ele_mvakfhits[50]/I");
	mytree_->Branch("ele_mvamishits", ele_mvamishits,"ele_mvamisthits[50]/I");
	mytree_->Branch("ele_mvadist", ele_mvadist,"ele_mvadist[50]/D");
	mytree_->Branch("ele_mvadcot", ele_mvadcot,"ele_mvadcot[50]/D");
	mytree_->Branch("ele_mvaeta", ele_mvaeta,"ele_mvaeta[50]/D");
	mytree_->Branch("ele_mvapt", ele_mvapt,"ele_mvapt[50]/D");
	mytree_->Branch("ele_mvaecalseed", ele_mvaecalseed,"ele_mvaecalseed[50]/I");
	
	//cout << "photon" << endl;
	// Photons
	mytree_->Branch("photon_N",&photon_N,"photon_N/I");
	m_photons = new TClonesArray ("TLorentzVector");
	mytree_->Branch ("photons", "TClonesArray", &m_photons, 256000,0);
	
	mytree_->Branch("pho_isFromMuon",_pho_isFromMuon, "pho_isFromMuon[50]/I");
	mytree_->Branch("pho_isEB",_pho_isEB, "pho_isEB[50]/I");
	mytree_->Branch("pho_isEE",_pho_isEE, "pho_isEE[50]/I");
	mytree_->Branch("pho_sigmaietaieta ", _pho_sigmaietaieta , "pho_sigmaietaieta[50]/D");
	mytree_->Branch("pho_he",             _pho_he , "pho_he[50]/D");
	mytree_->Branch("pho_r9",             _pho_r9 , "pho_r9[50]/D");
	mytree_->Branch("pho_TkIso03",        _pho_TkIso03 ,   "pho_TkIso03[50]/D");
	mytree_->Branch("pho_HCTkIso03",      _pho_HCTkIso03 , "pho_HCTkIso03[50]/D");
	mytree_->Branch("pho_emIso03",        _pho_emIso03,    "pho_emIso03[50]/D");
	mytree_->Branch("pho_hadIso03",       _pho_hadIso03 ,  "pho_hadIso03[50]/D");
	
	// SuperClusters
	if (fill_SC) {
	//cout << "SC" << endl;
	  mytree_->Branch("sc_N",   &_sc_N,   "sc_N/I");
	  mytree_->Branch("sc_E",   &_sc_E,   "sc_E[50]/D");
	  mytree_->Branch("sc_Et",  &_sc_Et,  "sc_Et[50]/D");
	  mytree_->Branch("sc_Eta", &_sc_Eta, "sc_Eta[50]/D");
	  mytree_->Branch("sc_Phi", &_sc_Phi, "sc_Phi[50]/D");
	  mytree_->Branch("sc_TkIso", &_sc_TkIso, "sc_TkIso[50]/D");
	  //cout << "end SC" << endl;
	}
	// PF jets
	_m_jets_pf   = new TClonesArray ("TLorentzVector");
	mytree_->Branch("jets_pf_N",  &_jets_pf_N,  "jets_pf_N/I");
	mytree_->Branch ("jets_pf",   "TClonesArray", &_m_jets_pf, 256000,0);
	
	mytree_->Branch ("jets_pf_chargedHadEFrac", &jets_pf_chargedHadEFrac,"jets_pf_chargedHadEFrac[50]/D]");
	mytree_->Branch ("jets_pf_chargedEmEFrac",  &jets_pf_chargedEmEFrac, "jets_pf_chargedEmEFrac[50]/D");
	mytree_->Branch ("jets_pf_chargedMuEFrac",  &jets_pf_chargedMuEFrac, "jets_pf_chargedMuEFrac[50]/D");
	
	mytree_->Branch ("jets_pf_neutralHadEFrac", &jets_pf_neutralHadEFrac, "jets_pf_neutralHadEFrac[50]/D");
	mytree_->Branch ("jets_pf_neutralEmEFrac",  &jets_pf_neutralEmEFrac,  "jets_pf_neutralEmEFrac[50]/D");
	mytree_->Branch ("jets_pf_PhotonEFrac",     &jets_pf_PhotonEFrac,     "jets_pf_PhotonEFrac[50]/D");
	
	mytree_->Branch ("jets_pf_chargedHadMultiplicity", &jets_pf_chargedHadMultiplicity, "jets_pf_chargedHadMultiplicity[50]/I");
	mytree_->Branch ("jets_pf_neutralHadMultiplicity", &jets_pf_neutralHadMultiplicity, "jets_pf_neutralHadMultiplicity[50]/I");
	
	mytree_->Branch ("jets_pf_chargedMultiplicity",    &jets_pf_chargedMultiplicity,    "jets_pf_chargedMultiplicity[50]/I");
	mytree_->Branch ("jets_pf_neutralMultiplicity",    &jets_pf_neutralMultiplicity,    "jets_pf_neutralMultiplicity[50]/I");
	
	mytree_->Branch ("jets_pf_nConstituents",          &jets_pf_nConstituents,          "jets_pf_nConstituents[50]/I");
	
	
	// Truth Leptons
	//cout << "truth leptons" << endl;
	_m_MC_gen_V = new TClonesArray ("TLorentzVector");
	mytree_->Branch ("MC_gen_V", "TClonesArray", &_m_MC_gen_V, 256000,0);
	mytree_->Branch ("MC_gen_V_pdgid",&_MC_gen_V_pdgid, "MC_gen_V_pdgid[10]/D");
	//
	_m_MC_gen_Higgs = new TClonesArray ("TLorentzVector");
	mytree_->Branch ("MC_gen_Higgs", "TClonesArray", &_m_MC_gen_Higgs, 256000,0);
	//mytree_->Branch ("MC_gen_Higgs_pdgid",&_MC_gen_Higgs_pdgid, "MC_gen_Higgs_pdgid[10]/D");
	//
	_m_MC_gen_leptons         = new TClonesArray ("TLorentzVector");
	_m_MC_gen_leptons_status1 = new TClonesArray ("TLorentzVector");
	_m_MC_gen_leptons_status2 = new TClonesArray ("TLorentzVector");
	mytree_->Branch ("MC_gen_leptons", "TClonesArray", &_m_MC_gen_leptons, 256000,0);
	mytree_->Branch ("MC_gen_leptons_status1", "TClonesArray", &_m_MC_gen_leptons_status1, 256000,0);
	mytree_->Branch ("MC_gen_leptons_status2", "TClonesArray", &_m_MC_gen_leptons_status2, 256000,0);
	mytree_->Branch ("MC_gen_leptons_pdgid",&_MC_gen_leptons_pdgid, "MC_gen_leptons_pdgid[10]/D");
	mytree_->Branch ("MC_gen_leptons_status1_pdgid",&_MC_gen_leptons_status1_pdgid, "MC_gen_leptons_status1_pdgid[10]/D");
	mytree_->Branch ("MC_gen_leptons_status2_pdgid",&_MC_gen_leptons_status2_pdgid, "MC_gen_leptons_status2_pdgid[10]/D");
	mytree_->Branch ("MC_pthat",&_MC_pthat,"MC_pthat/D");
	mytree_->Branch ("MC_flavor",&_MC_flavor,"MC_flavor[2]/I");
	_m_MC_gen_photons         = new TClonesArray ("TLorentzVector");
	mytree_->Branch ("MC_gen_photons", "TClonesArray", &_m_MC_gen_photons, 256000,0);
	mytree_->Branch ("MC_gen_photons_isFSR",&_MC_gen_photons_isFSR,"MC_gen_photons_isFSR[5000]/I");
	//cout << "end truth leptons" << endl;
	
	Nevt_Gen = 0;
	Nevt_H4lFilter = 0;
	Nevt_afterCleaning = 0;
	Nevt_afterSkim = 0;
       _trig_HLT_name= new std::vector<std::string>;
	
}

void NtupleProducer::endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup) 
{
	
	// 	edm::Handle<edm::MergeableCounter> Nevt_Gen_lumi;
	// 	iLumi.getByLabel("onTopCounter", Nevt_Gen_lumi);
	// 	Nevt_Gen = Nevt_Gen_lumi->value + Nevt_Gen;
	
	// 	edm::Handle<edm::MergeableCounter> Nevt_H4lFilter_lumi;
	// 	iLumi.getByLabel("firstCounter", Nevt_H4lFilter_lumi);
	// 	Nevt_H4lFilter = Nevt_H4lFilter_lumi->value + Nevt_H4lFilter; 
	
	// 	edm::Handle<edm::MergeableCounter> Nevt_afterCleaning_lumi;
	// 	iLumi.getByLabel("secondCounter", Nevt_afterCleaning_lumi);
	// 	Nevt_afterCleaning = Nevt_afterCleaning_lumi->value + Nevt_afterCleaning;
	
	// 	edm::Handle<edm::MergeableCounter> Nevt_afterSkim_lumi;
	// 	iLumi.getByLabel("thirdCounter",Nevt_afterSkim_lumi);
	// 	Nevt_afterSkim = Nevt_afterSkim_lumi->value + Nevt_afterSkim;
	
	// 	cout << "Nevt_Gen_lumi: " << Nevt_Gen << endl;
	// 	cout << "Nevt_H4lFilter_lumi: " << Nevt_H4lFilter << endl;
	// 	cout << "Nevt_afterCleaning_lumi: " << Nevt_afterCleaning << endl;
	// 	cout << "Nevt_afterSkim_lumi: " << Nevt_afterSkim<< endl;
	
}

void NtupleProducer::endJob(){
}

void NtupleProducer::analyze(const Event & iEvent, const EventSetup& iSetup){
	
	//inizialazed variables
  //	cout << "init" << endl;
	Init();
	
	//cout << "clear" << endl;
	//Fill Event Branches
	if (fill_L1trigger){
	  m_L1emIso->Clear();
	  m_L1emNonIso->Clear();
	}
	//	cout << "event" << endl;
	FillEvent (iEvent, iSetup);
	
	//   //Fill Vertices Branches
	//	cout << "vertices" << endl;
	FillVertices (iEvent, iSetup);
	
	//   //Fill Muons Branches
	//	cout << "muon" << endl;
	m_muons -> Clear() ;
	FillMuons (iEvent, iSetup);
	
	//   //Fill Electrons Branches
	//	cout << "electron clear" << endl;
	m_electrons -> Clear() ;
	//	cout << "electron" << endl;
	FillElectrons (iEvent, iSetup);
	
	//cout << "photon" << endl;
	//m_photons -> Clear() ;
	FillPhotons (iEvent, iSetup);
	if (fill_SC){
	//   cout << "SC" << endl;
	FillSC(iEvent, iSetup);
	}
	//	cout << "jets" << endl;
	FillJets(iEvent, iSetup);
	//	cout << "met" << endl;
	FillMET (iEvent, iSetup);
	
	if(isMC_ ) {
	  //	cout << "truth2" << endl;
		_m_MC_gen_V->Clear();
		_m_MC_gen_Higgs->Clear();
		_m_MC_gen_photons->Clear();
		_m_MC_gen_leptons->Clear();
		_m_MC_gen_leptons_status1->Clear();
		_m_MC_gen_leptons_status2->Clear();
		//cout << "truth" << endl;
		FillTruth(iEvent, iSetup);
	}
	
	//	cout << "tree" << endl;
	//Fill the tree
	mytree_->Fill();
	
} 


// ====================================================================================
void NtupleProducer::FillEvent (const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	_nEvent = iEvent.id().event();
	_nRun   = iEvent.id().run();
	_nLumi  = iEvent.luminosityBlock();
	
	// ----------------------------------------------------
	// Pile up
	// ----------------------------------------------------
	if(isMC_ ) {
		Handle<vector<PileupSummaryInfo> > PupInfo;
		iEvent.getByLabel(PileupSrc_, PupInfo);
		for (vector<PileupSummaryInfo>::const_iterator cand = PupInfo->begin();cand != PupInfo->end(); ++ cand) {
			_PU_N = cand->getPU_NumInteractions();
		} // loop on Pile up
	} // if MC
	
	
	InputTag theEleRhoTag = LeptonIsoHelper::getEleRhoTag(lepton_setup,lepton_setup);
	InputTag theMuRhoTag = LeptonIsoHelper::getMuRhoTag(lepton_setup,lepton_setup);
	//rh0/sigma correction
	Handle<double> MurhoHandle, sigmaHandle;
	iEvent.getByLabel(theMuRhoTag, MurhoHandle);
	//iEvent.getByLabel(SigmaRhoCorrection_, sigmaHandle);
	_PU_Murho   = *MurhoHandle;

	Handle<double> ElerhoHandle;
	iEvent.getByLabel(theEleRhoTag, ElerhoHandle);
	//iEvent.getByLabel(SigmaRhoCorrection_, sigmaHandle);
	_PU_Elerho   = *ElerhoHandle;
	//_PU_sigma = *sigmaHandle;
	
	
// 	//rh0/sigma correction
// 	Handle<double> MurhoHandle, sigmaHandle;
// 	iEvent.getByLabel(MuRhoCorrection_, MurhoHandle);
// 	//iEvent.getByLabel(SigmaRhoCorrection_, sigmaHandle);
// 	_PU_Murho   = *MurhoHandle;
// 
// 	Handle<double> ElerhoHandle;
// 	iEvent.getByLabel(EleRhoCorrection_, ElerhoHandle);
// 	//iEvent.getByLabel(SigmaRhoCorrection_, sigmaHandle);
// 	_PU_Elerho   = *ElerhoHandle;
// 	//_PU_sigma = *sigmaHandle;
	

	 
	
	
	// ----------------------------------------------------
	// Trigger
	// ----------------------------------------------------
	
	
	if(fill_L1trigger) {
		// %%%%%%%%%%%%%%
		// L1 Trigger
		// %%%%%%%%%%%%%%
		edm::Handle< l1extra::L1EmParticleCollection > emNonisolColl ;
		iEvent.getByLabel("l1extraParticles","NonIsolated", emNonisolColl ) ;
		edm::Handle< l1extra::L1EmParticleCollection > emIsolColl ;
		iEvent.getByLabel("l1extraParticles","Isolated", emIsolColl ) ;
		
		TClonesArray &L1emIso = *m_L1emIso;//electrons;
		
		// Isolated candidates
		int isocounter = 0;
		_trig_L1emIso_N = emIsolColl->size();
		
		for( l1extra::L1EmParticleCollection::const_iterator emItr = emIsolColl->begin(); emItr != emIsolColl->end() ;++emItr) {    
			setMomentum (myvector, emItr->p4());
			new (L1emIso[isocounter]) TLorentzVector (myvector);
			// From Trigger twiki
			// _trig_L1emIso_eta[isocounter]    = emItr->eta();
			//     _trig_L1emIso_phi[isocounter]    = emItr->phi();
			//     _trig_L1emIso_energy[isocounter] = emItr->energy();
			//     _trig_L1emIso_et[isocounter]     = emItr->et();
			isocounter++;
		}
		
		//cout << "L1 non iso" << endl;
		// Non Isolated candidates
		TClonesArray &L1emNonIso = *m_L1emNonIso;
		int nonisocounter = 0;
		_trig_L1emNonIso_N = emNonisolColl->size();
		
		for( l1extra::L1EmParticleCollection::const_iterator emItr = emNonisolColl->begin(); emItr != emNonisolColl->end() ;++emItr){  
			setMomentum (myvector, emItr->p4());
			new (L1emNonIso[nonisocounter]) TLorentzVector (myvector);
			
			nonisocounter++;
		} // for loop on Non Iso cand
		
	} // L1 trigger
	
	
	// %%%%%%%%%%%%%%
	// Fired Triggers
	// %%%%%%%%%%%%%%
	//cout << "fired" << endl;
	Handle<edm::TriggerResults> triggerResultsHandle;
	iEvent.getByLabel (HLTTag_,triggerResultsHandle);
	const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResultsHandle);
	
	//Get List of available Triggers
	//for (int in=0;in<(int)triggerNames.size();in++) {
	//cout << " Trigger Names " << in << " = " << triggerNames.triggerName(in) << endl;
	//} // for loop in triggernames
    
	// LOOP Over Trigger Results
	char trig_fired_names_local[10000];
	strcpy(trig_fired_names_local,"*");
	for (int iHLT = 0 ; 
	     iHLT<static_cast<int>(triggerResultsHandle->size()); 
	     ++iHLT) {	
	  
	  if (triggerResultsHandle->accept (iHLT)) {
	    if ( strlen(trig_fired_names_local) <= 9950) {
	      {
		const char* c_str();
		string hlt_string = triggerNames.triggerName(iHLT);
		strcat(trig_fired_names_local,hlt_string.c_str());
		strcat(trig_fired_names_local,"*");
	      }
	    }
	  } // if HLT
	}
	strcpy(trig_fired_names,trig_fired_names_local);

	// ----------------------
	//  get HLT candidates
	// ----------------------
	//cout << "get HLT" << endl;
	edm::Handle<trigger::TriggerEvent> trigEvent;
// 	iEvent.getByLabel(triggerEventTag_, trigEvent);
	
// 	const Int_t N_filter(trigEvent->sizeFilters());
// 	std::vector<Int_t> ID_filter; 
	
	// Print Official Filters
	//for(int ifi=0;ifi<N_filter;ifi++) {
	//cout << "filter tag " << ifi << " = " << trigEvent->filterTag(ifi) << endl;
	//} // for loop on filters
	
	int hlt_counter = 0;
	//TClonesArray &HLTobj = *m_HLT; //L1emIso;//electrons;
	
	// Loop on user's Filters
// 	for(int itrig=0;itrig< (int) HLT_Filters_.size();itrig++) {
// 		//cout << "itrig  = " << itrig << endl;
// 		//cout << "filter = " << HLT_Filters_[itrig] << endl;
// 		
// 		ID_filter.push_back(trigEvent->filterIndex(HLT_Filters_[itrig])); 
// 		
// 		const trigger::TriggerObjectCollection& TOC(trigEvent->getObjects());
// 		if( ID_filter[itrig] <  N_filter) { // !!! To be checked !!! trigEvent->size() ) {
// 			const trigger::Keys& keys( trigEvent->filterKeys(ID_filter[itrig])); 
// 			
// 			// Loop on HLT objects
// 			for ( int hlto = 0; hlto < (int) keys.size(); hlto++ ) {
// 				if(hlt_counter>499) continue;
// 				
// 				trigger::size_type hltf = keys[hlto];
// 				const trigger::TriggerObject& TrigObj(TOC[hltf]);
// 				
// 				_trig_HLT_eta[hlt_counter]    = TrigObj.eta();
// 				_trig_HLT_phi[hlt_counter]    = TrigObj.phi();
// 				_trig_HLT_energy[hlt_counter] = TrigObj.energy();
// 				_trig_HLT_pt[hlt_counter]     = TrigObj.pt();
// 				
// 				const std::string encodedFilterTag(HLT_Filters_[itrig].encode());
// 				//cout << "encoded = " << encodedFilterTag << endl;
// 				(*_trig_HLT_name)[hlt_counter]=encodedFilterTag;
// 					       
// 				hlt_counter++;
// 			} // for loop on HLT objects
// 		} // if idfilter<trigevent size
// 	} // for loop on filters
// 	
	_trig_HLT_N = hlt_counter;
	if(hlt_counter>499) { _trig_HLT_N = 500; cout << "Number of HLT Objects>500, trig_HLT_N set to 500" << endl;}
	
}	

// ====================================================================================
void NtupleProducer::FillVertices (const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	
  
// 	edm::Handle<reco::VertexCollection> recoPrimaryVertexCollection;
// 	iEvent.getByLabel(VerticesTag_,recoPrimaryVertexCollection);
	Handle<vector<reco::Vertex> >  recoPrimaryVertexCollection;
	iEvent.getByLabel("goodPrimaryVertices",recoPrimaryVertexCollection);
	
	const reco::VertexCollection & vertices = *recoPrimaryVertexCollection.product();

// 	edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
// 	iEvent.getByType(recoBeamSpotHandle);
// 	const reco::BeamSpot bs = *recoBeamSpotHandle;

	int vtx_counter=0;
	_vtx_N = recoPrimaryVertexCollection->size();
	
	// select the primary vertex as the one with higest sum of (pt)^2 of tracks  --> FIXME it should be now the default                                                                             
	//PrimaryVertexSorter PVSorter;
	//std::vector<reco::Vertex> sortedVertices = PVSorter.sortedList( *(recoPrimaryVertexCollection.product()) );

	if(_vtx_N > 0) {
		reco::VertexCollection::const_iterator firstPV = vertices.begin();
		GlobalPoint local_vertexPosition(firstPV->position().x(),
					  firstPV->position().y(),
					  firstPV->position().z());
		vertexPosition = local_vertexPosition;
	}
	else {
// 		GlobalPoint local_vertexPosition(bs.position().x(),
// 					  bs.position().y(),
// 					  bs.position().z());
		GlobalPoint local_vertexPosition(0,0,0); //TMP

		vertexPosition = local_vertexPosition;
	}
	
	for(reco::VertexCollection::const_iterator PV=vertices.begin() ; PV!=vertices.end() ; ++PV) { 
		if(vtx_counter > 49 ) continue;
		
		_vtx_normalizedChi2[vtx_counter] = PV->normalizedChi2(); 
		_vtx_ndof[vtx_counter] = PV->ndof();
		_vtx_nTracks[vtx_counter] = PV->tracksSize(); 
		_vtx_d0[vtx_counter] = PV->position().Rho();
		_vtx_x[vtx_counter] = PV->x();
		_vtx_y[vtx_counter] = PV->y();
		_vtx_z[vtx_counter] = PV->z();
		
		vtx_counter++;
	} 
    if(vtx_counter>49) { _vtx_N = 50; cout << "Number of primary vertices>9, vtx_N set to 100" << endl;}
	
}	

// ====================================================================================
void NtupleProducer::FillMuons(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	//Get the Beam Spot
// 	edm::Handle<reco::BeamSpot> recoBeamSpotHandle ;
// 	iEvent.getByType(recoBeamSpotHandle) ;
// 	const reco::BeamSpot bs = *recoBeamSpotHandle ;
	
	//Get the Muon collection
	Handle<pat::MuonCollection> muonsCol;
	iEvent.getByLabel(MuonTag_, muonsCol);
	
	//Get the Electron collection
	edm::Handle<pat::ElectronCollection> electronsCol;
	iEvent.getByLabel(EleTag_, electronsCol);
	
	//For computing the right isolation
	candiso =  computeIso(*muonsCol, *electronsCol);  
	
	iSetup.get<IdealMagneticFieldRecord>().get(magfield);  


	TClonesArray &muons = *m_muons;
	int mu_counter = 0;
	_muons_N = muonsCol->size() ;

	for (pat::MuonCollection::const_iterator imuons = muonsCol->begin();  
		 imuons != muonsCol->end(); 
		 ++imuons){
		
		if(mu_counter>49) continue;
		
		edm::Ref<pat::MuonCollection> muonsEdmRef(muonsCol, mu_counter);	
		
		//check if the muon is a PF one
		
		//int muPF_counter = 0;
		bool boolisPFMuon = false;
		
		//is a PF muon  
		if (imuons->isPFMuon()) boolisPFMuon =true;
		if ( boolisPFMuon)  _muons_isPF[mu_counter] = 1;
		else  _muons_isPF[mu_counter] = 0;
		
		// 4-vector
		setMomentum (myvector, imuons->p4());
		new (muons[mu_counter]) TLorentzVector (myvector);
		
		// charge
		_muons_charge[mu_counter] = imuons->charge(); 
		// provenance
		if(imuons->isTrackerMuon())    _muons_istracker[mu_counter]    = 1;
		if(imuons->isStandAloneMuon()) _muons_isstandalone[mu_counter] = 1;
		if(imuons->isGlobalMuon())     _muons_isglobal[mu_counter]     = 1;
		
		//vertex and hits
		reco::TrackRef gm = imuons->globalTrack();
		reco::TrackRef tk = imuons->innerTrack();
		if(imuons->isGlobalMuon()==1) {
// 			_muons_dxy[mu_counter]            = gm->dxy(bs.position()); 
// 			_muons_dz[mu_counter]             = gm->dz(bs.position()); 
			_muons_dxyPV[mu_counter]          = gm->dxy(math::XYZPoint(vertexPosition)); 
			_muons_dzPV[mu_counter]           = gm->dz(math::XYZPoint(vertexPosition)); 
			_muons_normalizedChi2[mu_counter] = gm->normalizedChi2(); 
			_muons_NmuonHits[mu_counter]      = gm->hitPattern().numberOfValidMuonHits(); 
		} 
		if(imuons->isGlobalMuon()==1 || imuons->isTrackerMuon()==1) {
			_muons_NtrackerHits[mu_counter] = tk->hitPattern().numberOfValidTrackerHits();
			_muons_NpixelHits[mu_counter]   = tk->hitPattern().numberOfValidPixelHits();
		} 
		_muons_Nmatches[mu_counter]             = imuons->numberOfMatches(); 
		_muons_caloCompatibility[mu_counter]    = imuons->caloCompatibility() ;
		_muons_segmentCompatibility[mu_counter] = ( muon::segmentCompatibility ( (*imuons) , reco::Muon::SegmentAndTrackArbitration) ) ;
		_muons_glbmuPromptTight[mu_counter]     = ( muon::isGoodMuon( (*imuons) , muon::GlobalMuonPromptTight) );
		
		if(imuons->innerTrack().isAvailable()){
			_muons_trkDxy[mu_counter]=imuons->innerTrack()->dxy();
			_muons_trkDxyError[mu_counter]=imuons->innerTrack()->dxyError();
// 			_muons_trkDxyB[mu_counter]=imuons->innerTrack()->dxy(bs.position()) ;
			_muons_trkDz[mu_counter]=imuons->innerTrack()->dz();
			_muons_trkDzError[mu_counter]=imuons->innerTrack()->dzError();
// 			_muons_trkDzB[mu_counter]=imuons->innerTrack()->dz(bs.position());
			_muons_trkChi2PerNdof[mu_counter]=imuons->innerTrack()->normalizedChi2();
			_muons_trkCharge[mu_counter]=imuons->innerTrack()->charge();
			_muons_trkNHits[mu_counter]=imuons->innerTrack()->numberOfValidHits();
			_muons_trkNPixHits[mu_counter]=imuons->innerTrack()->hitPattern().numberOfValidPixelHits();
			// Tracker muon properties
			_muons_trkmuArbitration[mu_counter]=(muon::segmentCompatibility( (*imuons),reco::Muon::SegmentAndTrackArbitration));
			_muons_trkmu2DCompatibilityLoose[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TM2DCompatibilityLoose));
			_muons_trkmu2DCompatibilityTight[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TM2DCompatibilityTight));
			_muons_trkmuOneStationLoose[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMOneStationLoose));
			_muons_trkmuOneStationTight[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMOneStationTight));
			_muons_trkmuLastStationLoose[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMLastStationLoose));
			_muons_trkmuLastStationTight[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMLastStationTight));
			_muons_trkmuOneStationAngLoose[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMOneStationAngLoose));
			_muons_trkmuOneStationAngTight[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMOneStationAngTight));
			_muons_trkmuLastStationAngLoose[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMLastStationAngLoose));
			_muons_trkmuLastStationAngTight[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMLastStationAngTight));
			_muons_trkmuLastStationOptimizedLowPtLoose[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMLastStationOptimizedLowPtLoose));
			_muons_trkmuLastStationOptimizedLowPtTight[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMLastStationOptimizedLowPtTight));
		}
	    
		// standard isolation
		_muons_nTkIsoR03[mu_counter] = imuons->isolationR03().nTracks; 
		_muons_nTkIsoR05[mu_counter] = imuons->isolationR05().nTracks;
		_muons_tkIsoR03[mu_counter]  = imuons->isolationR03().sumPt;
		_muons_tkIsoR05[mu_counter]  = imuons->isolationR05().sumPt;
		_muons_emIsoR03[mu_counter]  = imuons->isolationR03().emEt;
		_muons_emIsoR05[mu_counter]  = imuons->isolationR05().emEt;
		_muons_hadIsoR03[mu_counter] = imuons->isolationR03().hadEt,
		_muons_hadIsoR05[mu_counter] = imuons->isolationR05().hadEt;
		// HZZ isolation
		_muons_HZZisoTk[mu_counter]     = candiso[&(*imuons)];
		_muons_HZZisoTk5[mu_counter]    = candiso[&(*imuons)];

		//_muons_HZZisoTk[mu_counter]     = imuons->userIsolation(pat::User1Iso);
		//_muons_HZZisoTk5[mu_counter]    = imuons->userIsolation(pat::User1Iso);
		_muons_HZZisoEcal[mu_counter]   = imuons->ecalIso();
		_muons_HZZisoHcal[mu_counter]   = imuons->hcalIso();
		_muons_HZZisoComb[mu_counter]   = _muons_HZZisoTk[mu_counter] +_muons_HZZisoEcal[mu_counter] + _muons_HZZisoHcal[mu_counter] ;
		
		_muons_FsrIsoEcalDr005[mu_counter]= imuons->userIsolation(pat::User2Iso);
		_muons_FsrIsoEcalDr007[mu_counter]= imuons->userIsolation(pat::User3Iso);
		_muons_FsrIsoEcalDr010[mu_counter]= imuons->userIsolation(pat::User4Iso);
		
		//isolation for muons in the cone 0.4
		_muons_pfChargedHadIso[mu_counter]   = imuons->chargedHadronIso();
		_muons_pfNeutralHadIso[mu_counter]   = imuons->neutralHadronIso();
		_muons_pfPhotonIso[mu_counter]       = imuons->photonIso();
		_muons_pfCombRelIso[mu_counter]       =  LeptonIsoHelper::combRelIsoPF(lepton_setup, lepton_setup, _PU_Murho, *imuons);
		//FIXME these two for applying the beta corrections
		//_muons_pfChargedHadPUIso[mu_counter] = imuons->chargedAllIso();
		//_muons_pfChargedHadPUIso[mu_counter] = imuons->puChargedHadronIso();
		
		// SIP3D
		_muons_IP[mu_counter] = fabs(imuons->dB(pat::Muon::PV3D));
		_muons_IPError[mu_counter] = imuons->edB(pat::Muon::PV3D);	
		_muons_SIP[mu_counter] = _muons_IP[mu_counter]/_muons_IPError[mu_counter];
		
		_muons_isHLTMatch[mu_counter]   = isHLTMatch(&*imuons); // strange &*...
		_muons_isHLTMatch8[mu_counter]  = isHLTMatch1(&*imuons, isMC_, iEvent.id().run());
		_muons_isHLTMatch13[mu_counter] = isHLTMatch2(&*imuons, isMC_, iEvent.id().run());
		_muons_isHLTMatch17[mu_counter] = isHLTMatch3(&*imuons, isMC_, iEvent.id().run());
		
		bool fillerr=true;
		if (fillerr){
		  //cout<<"debug err 0"<<endl;
		  if (!imuons->track()) continue;
		  const Track* t = imuons->track().get();
		  //cout<<"debug err 1"<<endl;

		  //     TransientTrack tt = trackBuilder->build(t);
		  //     const CartesianTrajectoryError& cartErr = tt.initialFreeState().cartesianError();
		  const GlobalPoint mupoint(t->vx(), t->vy(),  t->vz());
		  //cout<<t->vx()<<" "<< t->vy()<< " "<<  t->vz()<<endl;
		  // cout<<mupoint<<endl;
		  // cout<<"debug err 1a1"<<endl;
		  const GlobalVector muvector(t->px(),t->py(),t->pz());
		  //cout<<muvector<<endl;
		  //cout<<"debug err 1a2"<<endl;
// 		  TrackCharge trcharge= t->charge();
		  // cout<<trcharge<<endl;
		  // cout<<"debug err 1a3"<<endl;
// 		  const MagneticField* B= magfield.product();
		  //cout<<B<<endl;

		  //cout<<"debug err 1a4"<<endl;
				  
		  GlobalTrajectoryParameters gp(GlobalPoint(t->vx(), t->vy(),  t->vz()),
						GlobalVector(t->px(),t->py(),t->pz()),
						t->charge(),
						magfield.product());
		  //cout<<"debug err 1a"<<endl;
		  JacobianCurvilinearToCartesian curv2cart(gp);
		  //cout<<"debug err 1b"<<endl;
		  CartesianTrajectoryError cartErr= ROOT::Math::Similarity(curv2cart.jacobian(), t->covariance());
		  //cout<<"debug err 1c"<<endl;
		  const AlgebraicSymMatrix66 & m = cartErr.matrix();
		  //cout<<"debug err 2"<<endl;

		  _muons_Err_00[mu_counter] = m[3][3];
		  _muons_Err_01[mu_counter] = m[3][4];
		  _muons_Err_02[mu_counter] = m[3][5];
		  _muons_Err_10[mu_counter] = m[4][3];
		  _muons_Err_11[mu_counter] = m[4][4];
		  _muons_Err_12[mu_counter] = m[4][5];
		  _muons_Err_20[mu_counter] = m[5][3];
		  _muons_Err_21[mu_counter] = m[5][4];
		  _muons_Err_22[mu_counter] = m[5][5];
		  //cout<<"debug err 3"<<endl;
		  
		}//fillerr
		++mu_counter;
	} // for loop on muons
	
	if(mu_counter>49) { _muons_N = 50; cout << "Number of muons>49, muons_N set to 50" << endl;}
	
	
	
	
	
	//	std::cout << muon->isGlobalMuon() << std::endl;
	//}
	
}	

// ====================================================================================
void NtupleProducer::FillElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	//Get the Beam Spot
// 	edm::Handle<reco::BeamSpot> recoBeamSpotHandle ;
// 	iEvent.getByType(recoBeamSpotHandle) ;
// 	const reco::BeamSpot bs = *recoBeamSpotHandle ;
	
	//Get the Electron collection
	edm::Handle<pat::ElectronCollection> electronsCol;
	iEvent.getByLabel(EleTag_, electronsCol);
	
	TClonesArray &electrons = *m_electrons;
	int counter = 0;
	ele_N = electronsCol->size() ;
	
	for (pat::ElectronCollection::const_iterator ielectrons = electronsCol->begin();  
		 ielectrons != electronsCol->end(); ++ielectrons)
	{
		
		if(counter>49) continue;
		

		// 4-vector
		setMomentum (myvector, ielectrons->p4());
		new (electrons[counter]) TLorentzVector (myvector);
		//
		ele_echarge[counter] = ielectrons->charge(); 
		//
		ele_he[counter]      = ielectrons->hadronicOverEm() ;
		ele_eseedpout[counter] = ielectrons->eSeedClusterOverPout();
		ele_ep[counter]        = ielectrons->eSuperClusterOverP() ;        
		ele_eseedp[counter]    = ielectrons->eSeedClusterOverP() ;         
		ele_eelepout[counter]  = ielectrons->eEleClusterOverPout() ;       
		//
		ele_pin_mode[counter]    = ielectrons->trackMomentumAtVtx().R() ; 
		ele_pout_mode[counter]   = ielectrons->trackMomentumOut().R() ; 
		ele_pTin_mode[counter]   = ielectrons->trackMomentumAtVtx().Rho() ; 
		ele_pTout_mode[counter]  = ielectrons->trackMomentumOut().Rho() ; 
		//
		ele_deltaetaseed[counter] = ielectrons->deltaEtaSeedClusterTrackAtCalo() ; 
		ele_deltaphiseed[counter] = ielectrons->deltaPhiSeedClusterTrackAtCalo() ;  
		ele_deltaetaele[counter]  = ielectrons->deltaEtaEleClusterTrackAtCalo() ;  
		ele_deltaphiele[counter]  = ielectrons->deltaPhiEleClusterTrackAtCalo() ; 
		ele_deltaetain[counter]   = ielectrons->deltaEtaSuperClusterTrackAtVtx();
		ele_deltaphiin[counter]   = ielectrons->deltaPhiSuperClusterTrackAtVtx();   
		//
		ele_sigmaietaieta[counter] = ielectrons->sigmaIetaIeta() ; 
		ele_sigmaetaeta[counter]   = ielectrons->sigmaEtaEta() ;
		ele_e15[counter]           = ielectrons->e1x5() ;
		ele_e25max[counter]        = ielectrons->e2x5Max() ;
		ele_e55[counter]           = ielectrons->e5x5() ;
		//ele_e1[counter]            = FIXME ;
		//ele_e33[counter]           = FIXME ;
		//
		// E/P combination
		ele_ecalErr[counter]   = ielectrons->ecalEnergyError();
		ele_trackErr[counter]  = ielectrons->trackMomentumError();
		ele_combErr[counter]   = ielectrons->p4Error(GsfElectron::P4_COMBINATION);
		ele_PFcombErr[counter] = ielectrons->p4Error(GsfElectron::P4_PFLOW_COMBINATION);
		//cout << "Errors (ecal/track/p4comb/PFprcomb) :" <<ele_ecalErr[counter] <<" " << ele_trackErr[counter]<<" "<< ele_combErr[counter]<<" "<< ele_PFcombErr[counter] <<endl;
		
		//regression
		ele_ecalRegressionEnergy[counter]  = ielectrons->ecalRegressionEnergy();
		ele_ecalRegressionError[counter] = ielectrons->ecalRegressionError();
		ele_ecalTrackRegressionEnergy[counter] = ielectrons->ecalTrackRegressionEnergy();
		ele_ecalTrackRegressionError[counter] = ielectrons->ecalTrackRegressionError();
		ele_ecalScale[counter] = ielectrons->ecalScale();               
		ele_ecalSmear[counter] = ielectrons->ecalSmear();                
		ele_ecalRegressionScale[counter] = ielectrons->ecalRegressionScale();     
		ele_ecalRegressionSmear[counter] = ielectrons->ecalRegressionSmear();     
		ele_ecalTrackRegressionScale[counter] = ielectrons->ecalTrackRegressionScale();
		ele_ecalTrackRegressionSmear[counter] = ielectrons->ecalTrackRegressionSmear();
		
		
// 		cout << "REGRESSION " << ele_ecalRegressionEnergy[counter]<< " " << ele_ecalRegressionError[counter] << endl;
		//
		ele_mva[counter]   = ielectrons->mva() ;
		//
		if (ielectrons->isEB()) ele_isbarrel[counter] = 1 ; 
		else  ele_isbarrel[counter] = 0 ;
		if (ielectrons->isEE()) ele_isendcap[counter] = 1 ; 
		else  ele_isendcap[counter] = 0 ;
		if (ielectrons->isEBEtaGap()) ele_isEBetaGap[counter] = 1 ;  
		if (ielectrons->isEBPhiGap()) ele_isEBphiGap[counter] = 1 ;  
		if (ielectrons->isEEDeeGap()) ele_isEEdeeGap[counter] = 1 ;  
		if (ielectrons->isEERingGap()) ele_isEEringGap[counter] = 1 ;
		if (ielectrons->ecalDrivenSeed()) ele_isecalDriven[counter] = 1 ;
		if (ielectrons->trackerDrivenSeed()) ele_istrackerDriven[counter] = 1 ;
		ele_eClass[counter]   = ielectrons->classification() ;
		//
		ele_missing_hits[counter] = ielectrons->gsfTrack()->numberOfLostHits();
		ele_lost_hits[counter]    = ielectrons->gsfTrack()->numberOfValidHits() ;
		ele_chi2_hits[counter]    = ielectrons->gsfTrack()->normalizedChi2() ;
		//
// 		ele_dxyB[counter] = ielectrons->gsfTrack()->dxy(bs.position()) ;
		ele_dxy[counter]  = ielectrons->gsfTrack()->dxy() ;
// 		ele_dzB[counter]  = ielectrons->gsfTrack()->dz(bs.position()) ;
		ele_dz[counter]   = ielectrons->gsfTrack()->dz() ;
// 		ele_dszB[counter] = ielectrons->gsfTrack()->dsz(bs.position()) ;
		ele_dsz[counter]  = ielectrons->gsfTrack()->dsz() ;
		//
		// Isolation variables
		ele_tkSumPt_dr03[counter]              = ielectrons->dr03TkSumPt() ;
		ele_ecalRecHitSumEt_dr03[counter]      = ielectrons->dr03EcalRecHitSumEt() ;
		ele_hcalDepth1TowerSumEt_dr03[counter] = ielectrons->dr03HcalDepth1TowerSumEt() ;
		ele_hcalDepth2TowerSumEt_dr03[counter] = ielectrons->dr03HcalDepth2TowerSumEt() ;
		ele_tkSumPt_dr04[counter]              = ielectrons->dr04TkSumPt() ;
		ele_ecalRecHitSumEt_dr04[counter]      = ielectrons->dr04EcalRecHitSumEt() ;
		ele_hcalDepth1TowerSumEt_dr04[counter] = ielectrons->dr04HcalDepth1TowerSumEt() ;
		ele_hcalDepth2TowerSumEt_dr04[counter] = ielectrons->dr04HcalDepth2TowerSumEt() ;
		//
		// Custom HCAL
		//  m_calotowers = new edm::Handle<CaloTowerCollection>() ;
		//     if (!iEvent.getByLabel("towerMaker",*m_calotowers)) //hcalTowers_
		//       { edm::LogError("ElectronHcalHelper::readEvent")<<"failed to get the hcal towers of label "<<hcalTowers_ ; }
		
		//     edm::Ref<pat::ElectronCollection> electronEdmRef(electronsCol, counter);
		//     double egHcalIsoConeSizeOutSmall = 0.3;
		//     double egHcalIsoConeSizeIn       = 0.0, egHcalIsoPtMin=0.0;
		//     int egHcalDepth1 = 1; 
		//     int egHcalDepth2 = 2;
		//hadDepth1Isolation03_  = new EgammaTowerIsolation(egHcalIsoConeSizeOutSmall,egHcalIsoConeSizeIn,egHcalIsoPtMin,egHcalDepth1,m_calotowers->product()) ;
		//hadDepth2Isolation03_  = new EgammaTowerIsolation(egHcalIsoConeSizeOutSmall,egHcalIsoConeSizeIn,egHcalIsoPtMin,egHcalDepth2,m_calotowers->product()) ;
		//double hcalDepth1TowerSumEt03 = hadDepth1Isolation03_->getTowerEtSum(&(*electronEdmRef)); //ielectrons)); //electronRef));
		//double hcalDepth2TowerSumEt03 = hadDepth2Isolation03_->getTowerEtSum(&(*electronEdmRef)); // electronEdmRefielectrons)); //electronRef));
		//ele_HCALFullConeSum[counter]  = hcalDepth1TowerSumEt03+hcalDepth2TowerSumEt03;
		//&((*EleHandle)[i])
		//
		ele_conv_dcot[counter] = ielectrons->convDist(); //userFloat("dcot");
		ele_conv_dist[counter] = ielectrons->convDcot(); //ielectrons->userFloat("dist");
		
		// FIXME: Always returns 0 :(
		//cout << " dcot = " << ielectrons->userFloat("dcot") << " dist = " << ielectrons->userFloat("dist") << endl;
		
		ele_expected_inner_hits[counter] = ielectrons->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
		//
		ele_eidVeryLoose[counter] = ielectrons->electronID("eidLoose"); //CiCVeryLoose"); 
		ele_eidLoose[counter] = ielectrons->electronID("eidLoose"); //CiCLoose");
		ele_eidMedium[counter] = ielectrons->electronID("eidMedium"); //CiCMedium");
		ele_eidTight[counter] = ielectrons->electronID("eidTight"); //eidCiCTight"); 
		ele_eidHZZVeryLoose[counter] = 0; //ielectrons->electronID("eidCiCHZZVeryLoose"); 
		ele_eidHZZLoose[counter] = 0; //ielectrons->electronID("eidCiCHZZLoose"); 
		ele_eidHZZMedium[counter] = 0; //ielectrons->electronID("eidCiCHZZMedium"); 
		ele_eidHZZTight[counter] = 0; //ielectrons->electronID("eidCiCHZZTight"); 
		ele_eidHZZSuperTight[counter] = 0; //ielectrons->electronID("eidCiCHZZSuperTight"); 
		ele_eidMVATrig[counter] = ielectrons->electronID("mvaTrigV0") ;
		ele_eidMVANoTrig[counter] = ielectrons->electronID("mvaNonTrigV0") ;
		
		// HZZ isolation
		ele_HZZisoTk[counter]   = ielectrons->userIsolation(pat::User1Iso);
		ele_HZZisoTk5[counter]   = ielectrons->userIsolation(pat::User1Iso);
		ele_HZZisoEcal[counter]   = ielectrons->dr03EcalRecHitSumEt();
		ele_HZZisoHcal[counter]   = ielectrons->dr03HcalTowerSumEt();
		ele_HZZisoComb[counter] = ele_HZZisoTk[counter] + ele_HZZisoEcal[counter] + ele_HZZisoHcal[counter];
		
		//isolation for electrons in the cone 0.4
		ele_pfChargedHadIso[counter]   = ielectrons->chargedHadronIso();
		ele_pfNeutralHadIso[counter]   = ielectrons->neutralHadronIso();
		ele_pfPhotonIso[counter]       = ielectrons->photonIso();
		ele_pfCombRelIso[counter]      = LeptonIsoHelper::combRelIsoPF(lepton_setup, lepton_setup, _PU_Elerho, *ielectrons);
		//FIXME these two for applying the beta corrections
		//ele_pfChargedHadPUIso[counter]  = ielectrons->chargedAllIso();
		//ele_pfChargedHadPUIso[counter]  = ielectrons->puChargedHadronIso();
		
		// SIP3D
		ele_IP[counter] = fabs(ielectrons->dB(pat::Electron::PV3D));
		ele_IPError[counter] = ielectrons->edB(pat::Electron::PV3D);	
		ele_SIP[counter] = ele_IP[counter]/ele_IPError[counter];
		
		// Get SuperCluster Informations
		//cout << " SuperCluster "<< endl;
		//	if(ielectrons->ecalDrivenSeed()) {
		reco::SuperClusterRef sclRef = ielectrons->superCluster();
		//cout << " SuperClusterRef" << endl;
		////math::XYZPoint sclPos        = ielectrons->superClusterPosition();
		//      cout << " pflow" << endl;
		
		//if (!ielectrons->ecalDrivenSeed() && ielectrons->trackerDrivenSeed()) 
		//sclRef = ielectrons->pflowSuperCluster();
			
			
		double R  = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y() +sclRef->z()*sclRef->z());
		double Rt = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y());
		ele_sclRawE[counter]  = sclRef->rawEnergy() ;

			ele_sclE[counter]     = sclRef->energy() ;
// 		ele_sclE[counter]     = ielectrons->correctedEcalEnergy();  //for 5XY
		ele_sclEt[counter]    = sclRef->energy()*(Rt/R) ;
		ele_sclEta[counter]   = sclRef->eta() ;
		ele_sclPhi[counter]   = sclRef->phi() ;
		ele_sclNclus[counter] = sclRef->clustersSize();
		
		ele_sclphiwidth[counter] = sclRef->phiWidth();
		//	} // if ECAL driven
		
		//cout << "sc Et = " << ele_sclEt[counter] << endl;
		ele_mvafbrem[counter] = ielectrons->fbrem();
		ele_mvadetain[counter] = ielectrons->deltaEtaSuperClusterTrackAtVtx();
		ele_mvadphiin[counter] = ielectrons->deltaPhiSuperClusterTrackAtVtx();
		ele_mvasieie[counter] = ielectrons->sigmaIetaIeta();
		ele_mvahoe[counter] = ielectrons->hcalOverEcal();
		ele_mvaeop[counter] = ielectrons->eSuperClusterOverP();
		ele_mvae1x5e5x5[counter] = (ielectrons->e5x5()) !=0. ? ielectrons->e1x5()/ielectrons->e5x5() : -1. ;
		ele_mvaeleopout[counter] = ielectrons->eEleClusterOverPout();
		bool validKF= (ielectrons->track().isNonnull());
		ele_mvakfchi2[counter] = validKF ? ielectrons->track()->normalizedChi2() : 0 ;
		ele_mvakfhits[counter] = validKF ? ielectrons->track()->hitPattern().trackerLayersWithMeasurement() : 0 ;
		ele_mvamishits[counter] = ielectrons->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
		ele_mvadist[counter] = ielectrons->convDist();
		ele_mvadcot[counter] = ielectrons->convDcot();
		ele_mvaeta[counter] = ielectrons->eta();
		ele_mvapt[counter] = ielectrons->pt();
		ele_mvaecalseed[counter] = ielectrons->ecalDrivenSeed();
		
		//fbrem
		ele_fbrem[counter] = ielectrons->fbrem();
// 		ele_fbrem[counter] = ielectrons->trackFbrem(); //works in 53X
// 		ele_SCfbrem[counter] = ielectrons->superClusterFbrem();//works in 53X
// 		ele_pfSCfbrem[counter] = ielectrons->pfSuperClusterFbrem();//works in 53X
		//cout<<"electron fbrem track= "<<ielectrons->trackFbrem()<<" fbrem eg calo= "<< ielectrons->superClusterFbrem()<<" fbrem pf calo= "<< ielectrons->pfSuperClusterFbrem()<<endl;


		++counter;
	}// for loop on electrons
	
	if(counter>49) { ele_N = 50; cout << "Number of electrons>49, electrons_N set to 50" << endl;}
	
	
}	


// ====================================================================================
void NtupleProducer::FillPhotons(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	//Get the Photon collection
  // edm::Handle<pat::PhotonCollection> photonsCol;
  // iEvent.getByLabel(PhotonTag_, photonsCol);
	//Get the Photon collection
  // edm::Handle<reco::PFCandidateCollection> photonsCol;
  edm::Handle<vector<pat::Photon> > photonsCol;
  iEvent.getByLabel(PhotonTag_, photonsCol);


	TClonesArray &photons = *m_photons;
	int counter = 0;
	photon_N = photonsCol->size() ;
	
	for( vector<pat::Photon>::const_iterator ipho = photonsCol->begin(); ipho != photonsCol->end(); ++ipho) {
	//	for (pat::PhotonCollection::const_iterator ipho = photonsCol->begin();  
	//	for (reco::PFCandidateCollection::const_iterator ipho = photonsCol->begin();  
	  //ipho != photonsCol->end(); 
	  //++ipho){

	  
	//	for( vector<cmg::Photon>::const_iterator ipho = photonsCol->begin(); ipho != photonsCol->end(); ++ipho) {

		
		if(counter>49) continue;
		
		// 4-vector
		//	cout << "photon: pt= " << ipho->pt() << endl;

		setMomentum (myvector, ipho->p4());
		new (photons[counter]) TLorentzVector (myvector);
		
		if (ipho->isFromMuon()){_pho_isFromMuon[counter]= 1;}//cout<<"isfrommuon =1 "<<endl;}
		else {_pho_isFromMuon[counter]= 0;}//cout<<"isfrommuon =0 "<<endl;}
		
		
	// 	_pho_isEB[counter] =  ipho->isEB();
// 		_pho_isEE[counter] =  ipho->isEE();
		
// 		_pho_sigmaietaieta[counter] = ipho->sigmaIetaIeta();
// 		_pho_he[counter]            = ipho->hadronicOverEm();
// 		_pho_r9[counter]            = ipho->r9();
		
// 		_pho_TkIso03[counter]         = ipho->trkSumPtSolidConeDR03();
// 		_pho_HCTkIso03[counter]       = ipho->trkSumPtHollowConeDR03();
// 		_pho_emIso03[counter]         = ipho->ecalRecHitSumEtConeDR03();
// 		_pho_hadIso03[counter]        = ipho->hcalTowerSumEtConeDR03();
		
		++counter;
	}// for loop on photons
	
	if(counter>49) { photon_N = 50; cout << "Number of photons>49, photon_N set to 50" << endl;}
  	
}	



// ====================================================================================
void NtupleProducer::FillSC(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	//Get the SC collection
	edm::Handle<reco::SuperClusterCollection> sc_coll;
	iEvent.getByLabel(SCTag_, sc_coll);
	
	_sc_N =  sc_coll->size();
	
	int index_sc = 0;
	
	edm::Handle<edm::ValueMap<double> > sc_tkiso_map;
	iEvent.getByLabel(edm::InputTag("SuperClusterIsolation"), sc_tkiso_map);
	
	// --------------------------
	//   Loop on SuperClusters 
	// --------------------------
	for( reco::SuperClusterCollection::const_iterator isc=sc_coll->begin(); isc!=sc_coll->end(); isc++) {
		if(index_sc>49) continue;
		double R  = TMath::Sqrt(isc->x()*isc->x() + isc->y()*isc->y() +isc->z()*isc->z());
		double Rt = TMath::Sqrt(isc->x()*isc->x() + isc->y()*isc->y());
		
		_sc_E[index_sc]   = isc->energy();
		_sc_Et[index_sc]  = isc->energy()*(Rt/R);
		_sc_Eta[index_sc] = isc->eta();
		_sc_Phi[index_sc] = isc->phi();
		
		edm::Ref<reco::SuperClusterCollection> SCRef(sc_coll, index_sc); //isc);
		//cout << "iso = " << (* sc_tkiso_map)[SCRef] << endl;
		
		_sc_TkIso[index_sc] = (* sc_tkiso_map)[SCRef];
		
		//edm::Ref<pat::ElectronCollection> electronEdmRef(electronsCol, counter);
		
		index_sc++;
	} // for loop on SC
	
	
	if(index_sc>49) { _sc_N = 50; cout << "Number of SC>49, SC_N set to 50" << endl;}
}	


// ====================================================================================
void NtupleProducer::FillMET (const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	
	// caloMET object (negative vector sum of calorimeter towers)
	//edm::Handle< edm::View<reco::CaloMET> > caloMEThandle;
	//iEvent.getByLabel("met", caloMEThandle);
	
	// MET object that corrects the basic calorimeter MET for muons
	// edm::Handle< edm::View<reco::CaloMET> > muCorrMEThandle;
	//   iEvent.getByLabel("corMetGlobalMuons", muCorrMEThandle);
	
	// MET object that corrects the basic calorimeter MET for muons and tracks
	//edm::Handle< edm::View<reco::MET> > tcMEThandle;
	//iEvent.getByLabel("tcMet", tcMEThandle);
	
	// MET object built as the (negative) vector sum of all particles (PFCandidates) reconstructed in the event
// 	edm::Handle< edm::View<pat::MET> > pfMEThandle;
// 	iEvent.getByLabel("patMETs", pfMEThandle);
	
	edm::Handle< edm::View<reco::MET> > pfMEThandle;
	iEvent.getByLabel("slimmedMETs", pfMEThandle);
	
	
	// CALO MET
	//  _met_calo_et  = (caloMEThandle->front() ).et();
	//   _met_calo_px  = (caloMEThandle->front() ).px();
	//   _met_calo_py  = (caloMEThandle->front() ).py();
	//   _met_calo_phi = (caloMEThandle->front() ).phi();
	//   _met_calo_set = (caloMEThandle->front() ).sumEt();
	//   _met_calo_sig = (caloMEThandle->front() ).mEtSig();
	
	// CALOMU MET
	//  _met_calomu_et  = (muCorrMEThandle->front() ).et();
	//   _met_calomu_px  = (muCorrMEThandle->front() ).px();
	//   _met_calomu_py  = (muCorrMEThandle->front() ).py();
	//   _met_calomu_phi = (muCorrMEThandle->front() ).phi();
	//   _met_calomu_set = (muCorrMEThandle->front() ).sumEt();
	//   _met_calomu_sig = (muCorrMEThandle->front() ).mEtSig();
	
	// TC MET
	// _met_tc_et  = (tcMEThandle->front() ).et();
	//   _met_tc_px  = (tcMEThandle->front() ).px();
	//   _met_tc_py  = (tcMEThandle->front() ).py();
	//   _met_tc_phi = (tcMEThandle->front() ).phi();
	//   _met_tc_set = (tcMEThandle->front() ).sumEt();
	//   _met_tc_sig = (tcMEThandle->front() ).mEtSig();
	
	// PFMET
	_met_pf_et  = (pfMEThandle->front() ).et();
	_met_pf_px  = (pfMEThandle->front() ).px();
	_met_pf_py  = (pfMEThandle->front() ).py();
	_met_pf_phi = (pfMEThandle->front() ).phi();
	_met_pf_set = (pfMEThandle->front() ).sumEt();
// 	_met_pf_sig = (pfMEThandle->front() ).mEtSig();
	
} // end of Fill MET


// ====================================================================================
void NtupleProducer::FillJets(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	
	// --------------------------------------------------
	// PF Jets
	// --------------------------------------------------
	//edm::Handle<pat::JetCollection>  pfjets;
	//  iEvent.getByLabel("cleanPatJets", pfjets);
// 	edm::Handle<edm::View<pat::Jet> >  pfjets;
// 	edm::Handle<std::vector<cmg::PFJet> >  pfjets;
	edm::Handle<edm::View<pat::Jet> >  pfjets;
	iEvent.getByLabel(JetTag_, pfjets);
	
	_jets_pf_N = pfjets->size();
	int index_pf_jets = 0;
	
	TClonesArray &jets_pf = *_m_jets_pf; 
	
	// Loop on Pf Jets
	//for ( pat::JetCollection::const_iterator ijets=pfjets->begin(); ijets!=pfjets->end(); ijets++) {  
// 	for(edm::View<pat::Jet>::const_iterator ijets=pfjets->begin(); ijets!=pfjets->end(); ijets++) {  
	for(edm::View<pat::Jet>::const_iterator ijets=pfjets->begin(); ijets!=pfjets->end(); ijets++) {  

// 	for(unsigned int ijets=0; ijets < pfjets->size(); ijets++) {  
		if (index_pf_jets>49) continue;
		
		setMomentum (myvector, ijets->p4());
		new (jets_pf[index_pf_jets]) TLorentzVector(myvector);
		
		jets_pf_chargedHadEFrac[index_pf_jets] = ijets->component(reco::PFCandidate::h).fraction();
		jets_pf_chargedEmEFrac[index_pf_jets]  = ijets->component(reco::PFCandidate::e).fraction();
		jets_pf_chargedMuEFrac[index_pf_jets]  = ijets->component(reco::PFCandidate::mu).fraction();
		
		jets_pf_neutralHadEFrac[index_pf_jets] = ijets->component(reco::PFCandidate::h0).fraction();
// 		jets_pf_neutralEmEFrac[index_pf_jets]  = ijets->neutralEmEnergyFraction ();
		jets_pf_PhotonEFrac[index_pf_jets]     = ijets->component(reco::PFCandidate::gamma).fraction();
		
		jets_pf_chargedHadMultiplicity[index_pf_jets] = ijets->component(reco::PFCandidate::h).number();
		jets_pf_neutralHadMultiplicity[index_pf_jets] = ijets->component(reco::PFCandidate::h0).number();
		
		jets_pf_chargedMultiplicity[index_pf_jets] = ijets->component(reco::PFCandidate::h).number()
						    +ijets->component(reco::PFCandidate::e).number()
						    +ijets->component(reco::PFCandidate::mu).number();
						    
		jets_pf_neutralMultiplicity[index_pf_jets] =  ijets->component(reco::PFCandidate::h0).number()
						    +ijets->component(reco::PFCandidate::gamma).number();
		
// 		jets_pf_chargedHadEFrac[index_pf_jets] = ijets->chargedHadronEnergyFraction ();
// 		jets_pf_chargedEmEFrac[index_pf_jets]  = ijets->chargedEmEnergyFraction ();
// 		jets_pf_chargedMuEFrac[index_pf_jets]  = ijets->chargedMuEnergyFraction ();
// 		
// 		jets_pf_neutralHadEFrac[index_pf_jets] = ijets->neutralHadronEnergyFraction ();
// 		jets_pf_neutralEmEFrac[index_pf_jets]  = ijets->neutralEmEnergyFraction ();
// 		jets_pf_PhotonEFrac[index_pf_jets]     = ijets->photonEnergyFraction();
// 		
// 		jets_pf_chargedHadMultiplicity[index_pf_jets] = ijets->chargedHadronMultiplicity ();
// 		jets_pf_neutralHadMultiplicity[index_pf_jets] = ijets->neutralHadronMultiplicity ();
// 		
// 		jets_pf_chargedMultiplicity[index_pf_jets] = ijets->chargedMultiplicity ();
// 		jets_pf_neutralMultiplicity[index_pf_jets] = ijets->neutralMultiplicity ();
// 		
		jets_pf_nConstituents[index_pf_jets]       = ijets->nConstituents();
		
		index_pf_jets++;
	} // for loop on pf jets
	
	if(index_pf_jets>49) { _jets_pf_N = 50; cout << "Number of pfjets>49, RECO_PFJETS_N set to 50" << endl;}
	
} // end of FillJets


// ====================================================================================
void NtupleProducer::FillTruth(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	
	//edm::Handle< GenEventInfoProduct > HepMCEvt;
	//iEvent.getByLabel(MCTag_, HepMCEvt);
	//if(HepMCEvt->hasBinningValues()) _MC_pthat = (HepMCEvt->binningValues())[0];
	//else  _MC_pthat = 0.0;
	
	edm::Handle<View<Candidate> > genCandidatesCollection;
	//iEvent.getByLabel("prunedGen", genCandidatesCollection);
	iEvent.getByLabel("genParticlesPruned", genCandidatesCollection);
	
	TClonesArray &MC_gen_V               = *_m_MC_gen_V;
	TClonesArray &MC_gen_Higgs           = *_m_MC_gen_Higgs;
	//cout << " photon" << endl;
	TClonesArray &MC_gen_photons         = *_m_MC_gen_photons;
	TClonesArray &MC_gen_leptons         = *_m_MC_gen_leptons;
	TClonesArray &MC_gen_leptons_status2 = *_m_MC_gen_leptons_status2;
	TClonesArray &MC_gen_leptons_status1 = *_m_MC_gen_leptons_status1;
	
	int counter             = 0;
	int counter_higgs       = 0;
	int counter_daughters   = 0;
	int counter_lep_status2 = 0;
	int counter_lep_status1 = 0;
	int counter_photon      = 0;
	
	// ----------------------------
	//      Loop on particles
	// ----------------------------
	for( View<Candidate>::const_iterator p = genCandidatesCollection->begin();p != genCandidatesCollection->end(); ++ p ) {
		
		// %%%%%%%%%%%%%%%%%%
		// If Higgs
		// %%%%%%%%%%%%%%%%%%
		if (p->pdgId() == 25 && p->status()==3) {
			setMomentum (myvector,p->p4());
			// 		  cout << "Higgs PdgId=" << p->pdgId() << " Higgs status=" << p->status() << " Mass=" << myvector.M() << endl;
			new (MC_gen_Higgs[counter_higgs]) TLorentzVector(myvector);
			counter_higgs++;
		} // if Higgs
		
		// %%%%%%%%%%%%%%%%%%
		// If Leptons from Z
		// %%%%%%%%%%%%%%%%%%
		if(fabs(p->pdgId())==11 || fabs(p->pdgId())==13 ||  fabs(p->pdgId())==15) {
			
			//   if(p->status()==1) {
			// 	cout << "Status1 pdgid = " << fabs(p->pdgId()) << endl;
			// 	cout << " Nmother = " << p->numberOfMothers() << endl;
			// 	for(unsigned int i=0;i<p->numberOfMothers();i++) {
			// 	  cout << " mother pdgid = " << p->mother(i)->pdgId() << " status = " << p->mother(i)->status() << endl;
			// 	} // for loop on mothers
			//       }
			
			if(p->status()==1) {
				if(p->numberOfMothers()>0) { // Need a mother...
					if(p->mother(0)->pdgId()== p->pdgId()) {
						setMomentum(myvector, p->p4());
						new (MC_gen_leptons_status1[counter_lep_status1]) TLorentzVector(myvector);
						_MC_gen_leptons_status1_pdgid[counter_lep_status1] = p->pdgId();
						counter_lep_status1++;
					}// if pdgid
				} // if mother
			} // if status 1
			
			if(p->status()==3) {
				if(p->numberOfMothers()>0) { // Need a mother...
					if(p->mother(0)->pdgId()==23) {  // If Mother is a Z 
						
						//cout << "number of daughters = " << p->numberOfDaughters() << " mother id = " << p->pdgId() << endl;
						
						if(p->numberOfDaughters()>0) { // Need a daughter...
							
							//cout << " status of daughter = " << p->daughter(0)->status() << " pdgid = " << p->daughter(0)->pdgId() << endl;
							
							// Status 2 Leptons
							if(p->daughter(0)->pdgId()==p->pdgId() && p->daughter(0)->status()==2) { // if daughter is lepton & status 2
								setMomentum(myvector, p->daughter(0)->p4());
								new (MC_gen_leptons_status2[counter_lep_status2]) TLorentzVector(myvector);
								_MC_gen_leptons_status2_pdgid[counter_lep_status2] = p->daughter(0)->pdgId();
								counter_lep_status2++;
								
								//cout << "dod = " << p->daughter(0)->daughter(0)->pdgId() << " status = " << p->daughter(0)->daughter(0)->status() << endl;
								
								// 		if(p->daughter(0)->daughter(0)->status()==2) {
								// 		  //cout << "Ndodod = " << p->daughter(0)->daughter(0)->numberOfDaughters() << endl;
								// 		  for(unsigned int i=0;i<p->numberOfDaughters();i++) {
								// 		    cout << " Dodod pdgid = " << 
								// 		  }
								// 		}
								
								// 	// Status 1 Leptons, from Status 2
								// 		if(p->daughter(0)->daughter(0)->pdgId()==p->pdgId() && p->daughter(0)->daughter(0)->status()==1) {
								// 		  setMomentum(myvector, p->daughter(0)->daughter(0)->p4());
								// 		  new (MC_gen_leptons_status1[counter_lep_status1]) TLorentzVector(myvector);
								// 		  _MC_gen_leptons_status1_pdgid[counter_lep_status1] = p->daughter(0)->daughter(0)->pdgId();
								// 		  counter_lep_status1++;
								// 		} // if status 1 from status 2 from status 3
								// 		else {
								// 		  if(fabs(p->pdgId())==11)
								// 		    //cout << "(status2) mother   id ? " << p->daughter(0)->pdgId() << " status ? " << p->daughter(0)->status() << endl;
								// 		    //cout << "(status2) daughter id ? " << p->daughter(0)->daughter(0)->pdgId() << " status ? " << p->daughter(0)->daughter(0)->status() << endl;
								// 		}
								
							} // if Daughter Status = 2
							
							//  // Status 1 Leptons, from Status 3
							// 	      if(p->daughter(0)->pdgId()==p->pdgId() && p->daughter(0)->status()==1) {
							// 		setMomentum(myvector, p->daughter(0)->p4());
							// 		new (MC_gen_leptons_status1[counter_lep_status1]) TLorentzVector(myvector);
							// 		_MC_gen_leptons_status1_pdgid[counter_lep_status1] = p->daughter(0)->pdgId();
							// 		counter_lep_status1++;
							// 	      } // if Daughters Status =1 (from status 3)
							// 	      else {
							// 		if(fabs(p->pdgId())==11)
							// 		  //cout << "(status3) daughter id ? " << p->daughter(0)->pdgId() << " status ? " << p->daughter(0)->status() << endl;
							// 	      }
							
							
						} // Need a daughter
					} // if MOther = Z
				} // if Nmother (status3) >0
			} // if hard scatter electron (status3)
		} // if leptons
		
		//  if(p->status()==2 && p->mother(i)->pdgId() ==  p->pdgId() && p->mother(i)->status()==3 
		// 	 &&  p->mother(i)->mother(0)==23) { 
		
		//       } // status 2 && mother, same pdgid and status 3
		
		
		//       if(p->status()==1 && p->mother(i)->pdgId() == p->pdgId() 
		// 	 && (p->mother(i)->status()==2 && p->mother(i)->mother(0)->pdgId()== p->pdgId())
		// || p->mother(i)->status()==3) {
		
		//       } // if 
		
		//     } // if leptons
		
		
		//       if(p->status()==1) { 
		// 	cout << "Ele Status 1, Nmother =  " << p->numberOfMothers() << " Pt = " << p->p4().Pt() << endl;
		// 	for(unsigned int i=0;i<p->numberOfMothers();i++) {
		// 	  cout << " Mother1 pdgid = " << p->mother(i)->pdgId() << " status = " << p->mother(i)->status() << endl;
		// 	} //mother
		//        } // status
		
		//       if(p->status()==2) { 
		// 	cout << "Ele Status 2, Nmother =  " << p->numberOfMothers() << " Pt = " << p->p4().Pt() << endl;
		// 	for(unsigned int i=0;i<p->numberOfMothers();i++) {
		// 	  cout << " Mother2 pdgid = " << p->mother(i)->pdgId() << " status = " << p->mother(i)->status() << endl;
		// 	} //mother
		//       } // status
		
		//       if(p->status()==3) { 
		// 	cout << "Ele Status 3, Nmother =  " << p->numberOfMothers() << " Pt = " << p->p4().Pt() << endl;
		// 	for(unsigned int i=0;i<p->numberOfMothers();i++) {
		// 	  cout << " Mother3 pdgid = " << p->mother(i)->pdgId() << " status = " << p->mother(i)->status() << endl;
		// 	} //mother
		//       } // status
		
		//     } // if electron
		
		// %%%%%%%%%%%%%%%%%%
		//     If W or Z
		// %%%%%%%%%%%%%%%%%%
		if (p->pdgId() == 23 || fabs(p->pdgId())==24) {
			
			//cout << "Z status= " << p->status() << " mass = " << p->mass() << endl;
			
			//  if(p->status()==2) {
			// 	//cout << "status2,NB daugthers ? " << p->numberOfDaughters() << endl;
			// 	for(unsigned int i=0;i<p->numberOfDaughters();i++) {
			// 	  //cout << "Status2 id daughters = " << p->daughter(i)->pdgId() << " Status = " << p->daughter(i)->status()  << " Pt = " << p->daughter(i)->p4().Pt() << endl;
			// 	} // loop on daughters
			//       }
			
			if(p->status()==3) {
				// Fill truth W,Z
				setMomentum (myvector,p->p4());
				new (MC_gen_V[counter]) TLorentzVector(myvector);
				_MC_gen_V_pdgid[counter] = p->pdgId();
				
				//size_t nfirstdau = p->numberOfDaughters();
				
				// Loop on daughters
				for(unsigned int i=0;i<p->numberOfDaughters();i++) {
					bool islep = false;
					if(fabs(p->daughter(i)->pdgId())==11) { _MC_flavor[counter] = 0; islep=true;} // electron
					if(fabs(p->daughter(i)->pdgId())==13) { _MC_flavor[counter] = 1; islep=true;} // muon
					if(fabs(p->daughter(i)->pdgId())==15) { _MC_flavor[counter] = 2; islep=true;} // taus
					
					//cout << "Status3 id daughters = " << p->daughter(i)->pdgId() << " Status = " << p->daughter(i)->status()  << " Pt = " << p->daughter(i)->p4().Pt() << endl;
					
					
					//  for(unsigned int j=0;j<p->daughter(i)->numberOfDaughters();j++) {
					// 	    //cout << "DoD status = " <<  p->daughter(i)->daughter(j)->status() << " Pt = " << p->daughter(i)->daughter(j)->p4().Pt() << endl;
					
					// 	    for(unsigned int k=0;k<p->daughter(i)->daughter(j)->numberOfDaughters();k++) {
					// 	      //cout << "DoDoD status = " <<  p->daughter(i)->daughter(j)->daughter(k)->status() << " Pt = " << p->daughter(i)->daughter(j)->daughter(k)->p4().Pt() << endl;
					// 	    } // loop on k
					
					// 	  } //
					
					if(islep) { // p->daughter(i)->status()==1) { ?!
						setMomentum(myvector, p->daughter(i)->p4());
						new (MC_gen_leptons[counter_daughters]) TLorentzVector(myvector);
						_MC_gen_leptons_pdgid[counter_daughters] = p->daughter(i)->pdgId();
						
						counter_daughters++;
					} // if is lepton
				} // for loop on daughters
				counter++;
			} // if status stable
		} // if W or Z
		
		// -------------
		//   Photons
		// -------------
		
		if(fabs(p->pdgId())==22) { 
			if(p->status()==1) {
				
				//cout << " if photon" << endl;
				
				setMomentum(myvector, p->p4());
				new (MC_gen_photons[counter_photon]) TLorentzVector(myvector);
				
				const Candidate * gen_photon = 0; 
				//const GenParticle * gen_photon2 = 0; 
				gen_photon  = &*p;
				//gen_photon2 =  &*p;
				//cout << "is fsr ?" << endl;
				//const reco::GenParticle* photon = getParent(gen_photon);
				bool fsr = isFSR(gen_photon); //const reco::GenParticle* genLep)
				
				//cout << "fsr = " << fsr << endl;
				_MC_gen_photons_isFSR[counter_photon] = fsr;
				//" pT = " <<p->p4().Pt()  << endl;
				
				//_MC_photon_isFSR[counter_photon] = isFSR(p, 
				
				counter_photon++;
			} // if status 1
		} // if photon
		
		
	} // for loop on particles
	
	
	
	
	
	// -------------------------
	// Only 4e events
	// -------------------------
	if(_MC_flavor[0]==0 && _MC_flavor[1]==0) {
		// ----------------------------
		//      Loop on particles
		// ----------------------------
		for( View<Candidate>::const_iterator p = genCandidatesCollection->begin();p != genCandidatesCollection->end(); ++ p ) {
			
			if(fabs(p->pdgId())==11) { 
				
				if(p->status()==1) {
					
					//cout << "N mother = " << p->numberOfMothers() << endl;
					for(unsigned int i=0;i<p->numberOfMothers();i++) {
						//cout << " mother pdgid = " << p->mother(i)->pdgId() << " status = " << p->mother(i)->status() << endl;
					} // for loop on mothers
				} // if status 1
				
				if(p->status()==2) {
					
					//cout << "N daughters (status2) = " << p->numberOfDaughters() << endl;
					for(unsigned int i=0;i<p->numberOfDaughters();i++) {
						//cout << " daughter pdgid = " << p->daughter(i)->pdgId() << " status = " << p->daughter(i)->status() << endl;
					} // for loop on mothers
				} // if status 2
				
				if(p->status()==3) {
					
					//cout << "N daughters (status3) = " << p->numberOfDaughters() << endl;
					for(unsigned int i=0;i<p->numberOfDaughters();i++) {
						//cout << " daughter pdgid = " << p->daughter(i)->pdgId() << " status = " << p->daughter(i)->status() << endl;
					} // for loop on mothers
				} // if status 3
				
				
				
			} // if electron
		} // for loop on gen particles
		
	} // if 4e event
	
	
	
	
	
	
	
	
} // end of FillTruth

// ====================================================================================
bool NtupleProducer::isFSR(const reco::Candidate* genLep) { //GenParticle* genLep){
	// ====================================================================================
	//if(genLep==0) return false;
	bool isFSR=false;
	//if(isPhoton(genLep)){
	
	//cout << "get parent" << endl;
	
	const reco::Candidate* parent=getParent(genLep);
	
	//cout << "parent" << endl;
	
	//cout << " pdgid = " << parent->pdgId() << endl;
	
	if(parent!=0) {
		if((fabs(parent->pdgId())==11
			|| fabs(parent->pdgId())==13
			|| fabs(parent->pdgId())==15
			) && parent->status()==2) {
			//cout << "parent2" << endl;
			const reco::Candidate* parent2=getParent(parent);
			//cout << "end parent2" << endl;
			if(parent2->pdgId()==23) isFSR = true;
		}
		//}
	} // parent!=0
	
	if(isFSR) return true;
	else return false;
}

// ====================================================================================
const reco::Candidate* NtupleProducer::getParent(const reco::Candidate* genLep) {
	// ====================================================================================
	if(genLep==0) return 0;
	
	//cout << "is flavor" << endl;
	
	// Find the actual lepton parent  
	int flavor = genLep->pdgId();
	
	//cout << " flavor " << flavor << endl;
	const reco::Candidate* m_mother = 0;
	
	//cout << "genLel mother = " << endl;
	//cout << genLep->mother() << endl;
	//cout << "end genLel mother = " << endl;
	
	while (genLep->mother()!=0 && genLep->mother()->pdgId() == flavor) {
		//cout  << " mother = " << genLep->mother()->pdgId() << endl;
		genLep = (const Candidate*) genLep->mother();
	}
	if(genLep->mother()!=0) return (const Candidate*) genLep->mother();
	else return m_mother;
}

// ====================================================================================
void NtupleProducer::setMomentum (TLorentzVector &myvector, const LorentzVector & mom)
// ====================================================================================
{
	myvector.SetPx (mom.Px());
	myvector.SetPy (mom.Py());
	myvector.SetPz (mom.Pz());
	myvector.SetE (mom.E());
}

// ====================================================================================
void NtupleProducer::Init()
// ====================================================================================
{
	_PU_beta = 0;
	_PU_Murho = 0;
	_PU_Elerho = 0;
	_PU_N =-999;
	//counters
	ele_N    = 0;
	_muons_N = 0;
	_vtx_N   = 0;
	photon_N = 0;
	_sc_N    = 0; 
	
	for (int i = 0 ; i < 50 ; ++i) {
		//muons
		_muons_charge[i] = 0 ;
		_muons_isPF[i] = 0 ;
		_muons_istracker[i] = 0 ; 
		_muons_isstandalone[i] = 0 ; 
		_muons_isglobal[i] = 0 ;
		_muons_dxy[i] = 0 ; 
		_muons_dz[i] = 0 ; 
		_muons_dxyPV[i] = 0 ; 
		_muons_dzPV[i] = 0 ; 
		_muons_normalizedChi2[i] = 0 ;
		_muons_NtrackerHits[i] = 0 ; 
		_muons_NpixelHits[i] = 0 ; 
		_muons_NmuonHits[i] = 0 ; 
		_muons_Nmatches[i] = 0 ;
		_muons_nTkIsoR03[i] = 0 ; 
		_muons_nTkIsoR05[i] = 0 ;
		_muons_tkIsoR03[i] = 0 ; 
		_muons_tkIsoR05[i] = 0 ; 
		_muons_emIsoR03[i] = 0 ;
		_muons_emIsoR05[i] = 0 ;
		_muons_hadIsoR03[i] = 0 ; 
		_muons_hadIsoR05[i] = 0 ;
		_muons_trkDxy[i] = 0 ; 
		_muons_trkDxyError[i] = 0 ; 
		_muons_trkDxyB[i] = 0 ;
		_muons_trkDz[i] = 0 ; 
		_muons_trkDzError[i] = 0 ; 
		_muons_trkDzB[i] = 0 ; 
		_muons_trkChi2PerNdof[i] = 0 ; 
		_muons_trkCharge[i] = 0 ; 
		_muons_trkNHits[i] = 0 ; 
		_muons_trkNPixHits[i] = 0 ;
		_muons_trkmuArbitration[i] = 0 ; 
		_muons_trkmu2DCompatibilityLoose[i] = 0 ; 
		_muons_trkmu2DCompatibilityTight[i] = 0 ; 
		_muons_trkmuOneStationLoose[i] = 0 ; 
		_muons_trkmuOneStationTight[i] = 0 ; 
		_muons_trkmuLastStationLoose[i] = 0 ;
		_muons_trkmuLastStationTight[i] = 0 ; 
		_muons_trkmuOneStationAngLoose[i] = 0 ; 
		_muons_trkmuOneStationAngTight[i] = 0 ; 
		_muons_trkmuLastStationAngLoose[i] = 0 ;
		_muons_trkmuLastStationAngTight[i] = 0 ; 
		_muons_trkmuLastStationOptimizedLowPtLoose[i] = 0 ; 
		_muons_trkmuLastStationOptimizedLowPtTight[i] = 0 ;
		_muons_caloCompatibility[i] = 0 ; 
		_muons_segmentCompatibility[i] = 0 ; 
		_muons_glbmuPromptTight[i] = 0 ; 
		_muons_HZZisoTk[i] = 0 ; 
		_muons_HZZisoTk5[i] = 0 ;
		_muons_HZZisoEcal[i]   = 0;
		_muons_HZZisoHcal[i]   = 0;
		_muons_HZZisoComb[i] = 0 ;
		_muons_FsrIsoEcalDr005[i]= 0;
		_muons_FsrIsoEcalDr007[i]= 0;
		_muons_FsrIsoEcalDr010[i]= 0;
		_muons_pfChargedHadIso[i]   = 0;
		_muons_pfNeutralHadIso[i]   = 0;
		_muons_pfPhotonIso[i]       = 0;
		_muons_pfChargedHadPUIso[i] = 0;
		_muons_pfCombRelIso[i]       = 0;
		_muons_IP[i] = 0 ; 
		_muons_IPError[i] = 0 ; 
		_muons_SIP[i] = 0 ;
		_muons_isHLTMatch[i] = 0;
		_muons_isHLTMatch8[i] = 0;
		_muons_isHLTMatch13[i] = 0;
		_muons_isHLTMatch17[i] = 0;
		_muons_Err_00[i] = 0;
		_muons_Err_01[i] = 0 ;
		_muons_Err_02[i] = 0 ;
		_muons_Err_10[i] = 0 ;
		_muons_Err_11[i] = 0 ;
		_muons_Err_12[i] = 0 ;
		_muons_Err_20[i] = 0 ;
		_muons_Err_21[i] = 0;
		_muons_Err_22[i] = 0 ;
		
		//vertices
		_vtx_x[i] = 0 ; 
		_vtx_y[i] = 0 ; 
		_vtx_z[i] = 0 ;
		_vtx_normalizedChi2[i] = 0 ; 
		_vtx_ndof[i] = 0 ; 
		_vtx_nTracks[i] = 0 ; 
		_vtx_d0[i] = 0 ;
		
		//electrons
		ele_echarge[i] = 0 ;
		ele_he[i] = 0 ;  
		ele_eseedpout[i] = 0 ;  
		ele_ep[i] = 0 ;  
		ele_eseedp[i] = 0 ;  
		ele_eelepout[i] = 0 ;        
		ele_deltaetaseed[i] = 0 ;  
		ele_deltaetaele[i] = 0 ;  
		ele_deltaphiseed[i] = 0 ;  
		ele_deltaphiele[i] = 0 ;  
		ele_deltaetain[i] = 0 ;  
		ele_deltaphiin[i] = 0 ; 
		ele_sigmaietaieta[i] = 0 ;  
		ele_sigmaetaeta[i] = 0 ;  
		ele_e15[i] = 0 ;  
		ele_e25max[i] = 0 ;  
		ele_e55[i] = 0 ;  
		ele_e1[i] = 0 ;  
		//ele_e33[i] = 0 ;  
		//ele_e2overe9[i] = 0 ; 
		ele_pin_mode[i] = 0 ;  
		ele_pout_mode[i] = 0 ;  
		ele_pTin_mode[i] = 0 ;  
		ele_pTout_mode[i] = 0 ;  
		ele_fbrem[i] = 0 ;  
		ele_SCfbrem[i] = 0 ;  
		ele_pfSCfbrem[i] = 0 ;  
		ele_mva[i] = 0 ; 
		ele_isbarrel[i] = 0 ;  
		ele_isendcap[i] = 0 ;  
		ele_isEBetaGap[i] = 0 ;  
		ele_isEBphiGap[i] = 0 ;  
		ele_isEEdeeGap[i] = 0 ;  
		ele_isEEringGap[i] = 0 ; 
		ele_isecalDriven[i] = 0 ;  
		ele_istrackerDriven[i] = 0 ; 
		ele_eClass[i] = 0 ;
		ele_missing_hits[i] = 0 ; 
		ele_lost_hits[i] = 0 ; 
		ele_chi2_hits[i] = 0 ; 
		ele_dxyB[i] = 0 ; 
		ele_dxy[i] = 0 ; 
		ele_dzB[i] = 0 ; 
		ele_dz[i] = 0 ; 
		ele_dszB[i] = 0 ; 
		ele_dsz[i] = 0 ;              
		ele_tkSumPt_dr03[i] = 0 ;  
		ele_ecalRecHitSumEt_dr03[i] = 0 ;  
		ele_hcalDepth1TowerSumEt_dr03[i] = 0 ;  
		ele_hcalDepth2TowerSumEt_dr03[i] = 0 ; 
		ele_tkSumPt_dr04[i] = 0 ;  
		ele_ecalRecHitSumEt_dr04[i] = 0 ;  
		ele_hcalDepth1TowerSumEt_dr04[i] = 0 ;  
		ele_hcalDepth2TowerSumEt_dr04[i] = 0 ; 
		ele_conv_dcot[i] = 0 ;
		ele_conv_dist[i] = 0 ;
		ele_expected_inner_hits[i] = 0;
		ele_eidVeryLoose[i] = 0 ; 
		ele_eidLoose[i] = 0 ; 
		ele_eidMedium[i] = 0 ; 
		ele_eidTight[i] = 0 ;
		ele_eidHZZVeryLoose[i] = 0 ; 
		ele_eidHZZLoose[i] = 0 ; 
		ele_eidHZZMedium[i] = 0 ; 
		ele_eidHZZTight[i] = 0 ; 
		ele_eidHZZSuperTight[i] = 0 ; 
		ele_eidMVANoTrig[i] = -999 ;
		ele_eidMVATrig[i] = -999 ;
		
		
		ele_HZZisoTk[i] = 0 ; 
		ele_HZZisoTk5[i] = 0 ; 
		ele_HZZisoEcal[i]   = 0;
		ele_HZZisoHcal[i]   = 0;
		ele_HZZisoComb[i] = 0 ; 
		ele_pfChargedHadIso[i]   = 0;
		ele_pfNeutralHadIso[i]   = 0;
		ele_pfPhotonIso[i]       = 0;
		ele_pfChargedHadPUIso[i] = 0;
		ele_pfCombRelIso[i] = 0;
		
		ele_IP[i] = 0 ; 
		ele_IPError[i] = 0 ; 
		ele_SIP[i] = 0 ; 
		//
		ele_sclRawE[i]=0;
		ele_sclE[i]=0;
		ele_sclEt[i]=0;
		ele_sclEta[i]=0;
		ele_sclPhi[i]=0;
		ele_sclNclus[i]=0;
		ele_sclphiwidth[i]=0;
		
		ele_ecalErr[i]=0;
		ele_trackErr[i]=0;
		ele_combErr[i]=0;
		ele_PFcombErr[i]=0;
		
		
		ele_ecalRegressionEnergy[i]  = 0;
		ele_ecalRegressionError[i] = 0;
		ele_ecalTrackRegressionEnergy[i]  = 0;
		ele_ecalTrackRegressionError[i]  = 0;
		ele_ecalScale[i]  = 0;
		ele_ecalSmear[i]  = 0;
		ele_ecalRegressionScale[i]  = 0;
		ele_ecalRegressionSmear[i]  = 0;
		ele_ecalTrackRegressionScale[i]  = 0;
		ele_ecalTrackRegressionSmear[i]  = 0;
		
		ele_mvafbrem[i]=0;
		ele_mvadetain[i]=0;
		ele_mvadphiin[i]=0;
		ele_mvasieie[i]=0;
		ele_mvahoe[i]=0;
		ele_mvaeop[i]=0;
		ele_mvae1x5e5x5[i]=0;
		ele_mvaeleopout[i]=0;
		ele_mvakfchi2[i]=0;
		ele_mvakfhits[i]=0;
		ele_mvamishits[i]=0;
		ele_mvadist[i]=0;
		ele_mvadcot[i]=0;
		ele_mvaeta[i]=0;
		ele_mvapt[i]=0;
		ele_mvaecalseed[i]=0;
		
		ele_sclphiwidth[i]= 0;
		
		// Photons
		_pho_isFromMuon[i]= -999;
		_pho_isEB [i]=0;
		_pho_isEE [i]=0;
		
		_pho_sigmaietaieta [i]=0;
		_pho_he [i]=0;
		_pho_r9 [i]=0; 
		
		_pho_TkIso03 [i]=0; 
		_pho_HCTkIso03 [i]=0;
		_pho_emIso03 [i]=0; 
		_pho_hadIso03[i]=0;
		
		// SuperCluster
		if (fill_SC){
		  _sc_E[i]   = 0.; 
		  _sc_Et[i]  = 0.; 
		  _sc_Eta[i] = 0.; 
		  _sc_Phi[i] = 0.; 
		  _sc_TkIso[i] = 0;
		}
		
	} // for loop on 100
	
	// MET
	_met_calo_et  = 0.;
	_met_calo_px  = 0.; 
	_met_calo_py  = 0.; 
	_met_calo_phi = 0.; 
	_met_calo_set = 0.; 
	_met_calo_sig = 0.; 
	
	_met_calomu_et  = 0.;
	_met_calomu_px  = 0.; 
	_met_calomu_py  = 0.;
	_met_calomu_phi = 0.; 
	_met_calomu_set = 0.; 
	_met_calomu_sig = 0; 
	
	_met_tc_et  = 0.;
	_met_tc_px  = 0.; 
	_met_tc_py  = 0.; 
	_met_tc_phi = 0.; 
	_met_tc_set = 0.; 
	_met_tc_sig = 0.; 
	
	_met_pf_et  = 0.;
	_met_pf_px  = 0.; 
	_met_pf_py  = 0.; 
	_met_pf_phi = 0.; 
	_met_pf_set = 0.; 
	_met_pf_sig = 0.; 
	
	
	//cout << "HLT init" << endl;
	
	// HLT
	_trig_HLT_N = 0;
	_trig_HLT_name->clear();
	_trig_HLT_name->reserve(500); 
	for(int ihlt=0;ihlt<500;ihlt++) {
		//cout << "i = " << ihlt << endl;
		_trig_HLT_eta[ihlt]    = 0.; 
		_trig_HLT_phi[ihlt]    = 0.; 
		_trig_HLT_energy[ihlt] = 0.; 
		_trig_HLT_pt[ihlt]     = 0.;
		_trig_HLT_name->push_back("*");
	} // for loop on hlt
	//cout << "HLT end init" << endl;
	
	//_trig_HLT_name[ihlt]   = "";
	//}
	
	
	for (int i = 0 ; i < 10 ; ++i) {
		_MC_gen_V_pdgid[i]               = 0;
		_MC_gen_leptons_pdgid[i]         = 0;
		_MC_gen_leptons_status1_pdgid[i] = 0;
		_MC_gen_leptons_status2_pdgid[i] = 0;
	}
	_MC_pthat = 0;
	_MC_flavor[0] = 0;
	_MC_flavor[1] = 0;
	
	for (int i = 0 ; i < 5000 ; ++i) {
		_MC_gen_photons_isFSR[i] = 0;
	}
	
	// PF Jets
	_jets_pf_N = 0;
	
	for(int ipfjet=0;ipfjet<50;ipfjet++) {
		jets_pf_chargedHadEFrac[ipfjet] = 0.;
		jets_pf_chargedEmEFrac[ipfjet]  = 0.;
		jets_pf_chargedMuEFrac[ipfjet]  = 0.;
		
		jets_pf_neutralHadEFrac[ipfjet] = 0.;
		jets_pf_neutralEmEFrac[ipfjet]  = 0.;
		jets_pf_PhotonEFrac[ipfjet]     = 0.;
		
		jets_pf_chargedHadMultiplicity[ipfjet] = 0;
		jets_pf_neutralHadMultiplicity[ipfjet] = 0;
		
		jets_pf_chargedMultiplicity[ipfjet] = 0;
		jets_pf_neutralMultiplicity[ipfjet] = 0;
		
		jets_pf_nConstituents[ipfjet]      = 0;
	} // for loop on PFjets
	
	
}


//=========================================================
//This is for the LOWER di-lepton threshold
bool NtupleProducer::isHLTMatch(const pat::Muon* mu) {
	//=========================================================
	if (mu==0) return false;
	
	//For Winter10 HighPU - SingleMu 
	//return !(mu->triggerObjectMatchesByPath("HLT_Mu15_v1").empty());
	//For Winter10 and Dec22 2010 Data
	//return (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu3").empty())||!(mu->triggerObjectMatchesByPath("HLT_DoubleMu3_v*").empty()));
	//For 2011 Data + Spring11 sample
	return (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu7_v*").empty()) || 
			(!(mu->triggerObjectMatchesByPath("HLT_DoubleMu5_v*").empty()) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").at(0).pt()>7)));
	
} // if

//==================================================================================================================
//This is for the HIGHER di-lepton threshold
bool NtupleProducer::isHLTMatch3(const pat::Muon* mu, bool isMC, int irun) {
	//==================================================================================================================
	//  if (mu==0) return false;
	
	//   // --> For 2011 Data 5E32 + Spring11 sample
	//   //   return (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu7_v*").empty()) || 
	//   // 	    (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu5_v*").empty()) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").at(0).pt()>7)));
	//   // <--
	
	//   if (!isMC && irun<=149442){ //Special for emulation of 1e33 menu on 2010 data
	//     //    return ((!(mu->triggerObjectMatchesByPath("HLT_DoubleMu3").empty()) || !(mu->triggerObjectMatchesByPath("HLT_DoubleMu3_v*").empty())) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered").at(0).pt()>13));
	//     return mu->pt()>17;
	//   } else {
	//     //For 2011 Data 1E33 + Spring11 sample
	//     return ((!(mu->triggerObjectMatchesByPath("HLT_Mu13_Mu8_v*").empty()) && (!mu->triggerObjectMatchesByFilter("hltSingleMu13L3Filtered13").empty()))|| (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu5_v*").empty()) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").at(0).pt()>17)) || (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu7_v*").empty()) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").at(0).pt()>17)));
	//   }
	
	if (mu==0) return false;
	
	
	// --> For 2011 Data 5E32 + Spring11 sample
	//   return (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu7_v*").empty()) || 
	// 	    (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu5_v*").empty()) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").at(0).pt()>7)));
	// <--
	
	if (!isMC && irun<=149442){ //Special for emulation of 1e33 menu on 2010 data
		//    return ((!(mu->triggerObjectMatchesByPath("HLT_DoubleMu3").empty()) || !(mu->triggerObjectMatchesByPath("HLT_DoubleMu3_v*").empty())) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered").at(0).pt()>13));
		return mu->pt()>13;
	} else {
		//For 2011 Data up to 3E33 + Spring11/Summer11 sample
		return ((!mu->triggerObjectMatchesByFilter("hltSingleMu13L3Filtered17").empty())||
				((!mu->triggerObjectMatchesByFilter("hltSingleMu13L3Filtered13").empty()) && (mu->triggerObjectMatchesByFilter("hltSingleMu13L3Filtered13").at(0).pt()>17)) || 
				((!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").at(0).pt()>17)) || 
				((!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").at(0).pt()>17)));
	} // else
	
}


//==================================================================================================================
//This is for the HIGHER di-lepton threshold
bool NtupleProducer::isHLTMatch2(const pat::Muon* mu, bool isMC, int irun) {
	//==================================================================================================================
	if (mu==0) return false;
	
	// --> For 2011 Data 5E32 + Spring11 sample
	//   return (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu7_v*").empty()) || 
	// 	    (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu5_v*").empty()) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").at(0).pt()>7)));
	// <--
	
	if (!isMC && irun<=149442){ //Special for emulation of 1e33 menu on 2010 data
		//    return ((!(mu->triggerObjectMatchesByPath("HLT_DoubleMu3").empty()) || !(mu->triggerObjectMatchesByPath("HLT_DoubleMu3_v*").empty())) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered").at(0).pt()>13));
		return mu->pt()>13;
	} else {
		////For 2011 Data 1E33 + Spring11 sample
		//return ((!(mu->triggerObjectMatchesByPath("HLT_Mu13_Mu8_v*").empty()) && (!mu->triggerObjectMatchesByFilter("hltSingleMu13L3Filtered13").empty()))|| (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu5_v*").empty()) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").at(0).pt()>13)) || (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu7_v*").empty()) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").at(0).pt()>13)));
		//}
		
		return ((!mu->triggerObjectMatchesByFilter("hltSingleMu13L3Filtered13").empty())||((!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").at(0).pt()>13)) || ((!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").at(0).pt()>13)));
	}// else
	
	
}


//==================================================================================================================
//This is for the LOWER di-lepton threshold
bool NtupleProducer::isHLTMatch1(const pat::Muon* mu, bool isMC, int irun) {
	//==================================================================================================================
	if (mu==0) return false;
	
	//---> For 2011 Data 5E32 + Spring11 sample
	//   return (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu7_v*").empty()) || 
	// 	  (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu5_v*").empty()) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").at(0).pt()>7.)));
	// <---
	
	
	if (!isMC && irun<=149442){ //Special for emulation of 1e33 menu on 2010 data
		//    return ((!(mu->triggerObjectMatchesByPath("HLT_DoubleMu3").empty()) || !(mu->triggerObjectMatchesByPath("HLT_DoubleMu3_v*").empty())) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered").at(0).pt()>8));
		return mu->pt()>8;
	} else {
		////For 2011 Data 1E33 + Spring11 sample
		//return ((!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered8").empty())|| (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu5_v*").empty()) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").at(0).pt()>8))|| (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu7_v*").empty()) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").at(0).pt()>8)));
		//}  
		//For 2011 Data 3E33 + Spring11/Summer11 sample
		return ((!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered8").empty())|| (!mu->triggerObjectMatchesByFilter("hltDiMuonL3p5PreFiltered8").empty()) || ((!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").at(0).pt()>8))|| ((!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").at(0).pt()>8)));
	} // else  
	
}



NtupleProducer::CandIsolationMap
NtupleProducer::computeIso(const vector<pat::Muon>& muons, const vector<pat::Electron>& electrons) {
	// Here we want to compute isolation after vetoing all, and only, the four leptons of the H candidate.
	
	CandIsolationMap result; 
	const double coneSize     = 0.3;
	const double vetoConeSize = 0.015;
	
	for (vector<pat::Muon>::const_iterator mu = muons.begin(); mu != muons.end(); ++mu) {
		const IsoDeposit* depTk   = mu->isoDeposit(pat::TrackIso); // Standard
        
		//Fill vetos for all other muons. There's no need to check that they are in the cone, IsoDeposit::depositWithin does that for us.
		reco::IsoDeposit::Vetos myVetos;    
		for (vector<pat::Muon>::const_iterator othermu = muons.begin(); othermu != muons.end(); ++othermu) {
			if(((othermu)->isGlobalMuon() || (othermu)->isTrackerMuon()) && (othermu)->pt()>5.) { //|| (othermu)->isTrackerMuon()) && (othermu)->pt()>3.)// for PF and if the IDminPT is moved at 3
				reco::IsoDeposit::Veto myVeto(othermu->isoDeposit(pat::TrackIso)->direction(), vetoConeSize);
				myVetos.push_back(myVeto);   
			}
		}
		
		for (vector<pat::Electron>::const_iterator othere = electrons.begin(); othere != electrons.end(); ++othere) {
			if((othere)->pt()>7){
				reco::IsoDeposit::Veto myVeto( reco::IsoDeposit::Direction(othere->gsfTrack()->eta(), othere->gsfTrack()->phi()), vetoConeSize);
				myVetos.push_back(myVeto);   
			}
		}
		result[&(*mu)] = depTk->depositWithin(coneSize, myVetos, true);
	}
	return result;
}




// //==================================================================================================================
// //This is for the LOWER di-lepton threshold
// bool NtupleProducer::isHLTMatch1(const pat::Muon* mu, bool isMC, int irun) {
//   //==================================================================================================================
//   if (mu==0) return false;

//   //For Winter10 HighPU - SingleMu 
//   //return !(mu->triggerObjectMatchesByPath("HLT_Mu15_v1").empty());
//   //For Winter10 and Dec22 2010 Data
//   //return (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu3").empty())||!(mu->triggerObjectMatchesByPath("HLT_DoubleMu3_v*").empty()));


//   //---> For 2011 Data 5E32 + Spring11 sample
//   //   return (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu7_v*").empty()) || 
//   // 	  (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu5_v*").empty()) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").at(0).pt()>7.)));
//   // <---

//   if (!isMC && irun<=149442){ //Special for emulation of 1e33 menu on 2010 data
//     return ((!(mu->triggerObjectMatchesByPath("HLT_DoubleMu3").empty()) || !(mu->triggerObjectMatchesByPath("HLT_DoubleMu3_v*").empty())) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered").at(0).pt()>13));
//   } else {
//     //For 2011 Data 1E33 + Spring11 sample
//     return ((!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered8").empty())|| (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu5_v*").empty()) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").at(0).pt()>8))|| (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu7_v*").empty()) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").at(0).pt()>8)));
//   }  
// }

// //==================================================================================================================
// //This is for the HIGHER di-lepton threshold
// bool NtupleProducer::isHLTMatch2(const pat::Muon* mu, bool isMC, int irun) {
// //==================================================================================================================
//   if (mu==0) return false;


//   // --> For 2011 Data 5E32 + Spring11 sample
//   //   return (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu7_v*").empty()) || 
//   // 	    (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu5_v*").empty()) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").at(0).pt()>7)));
//   // <--

//   if (!isMC && irun<=149442){ //Special for emulation of 1e33 menu on 2010 data
//     return ((!(mu->triggerObjectMatchesByPath("HLT_DoubleMu3").empty()) || !(mu->triggerObjectMatchesByPath("HLT_DoubleMu3_v*").empty())) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered").at(0).pt()>13));
//   } else {
//     //For 2011 Data 1E33 + Spring11 sample
//     return ( (!(mu->triggerObjectMatchesByPath("HLT_Mu13_Mu8_v*").empty()) && (!mu->triggerObjectMatchesByFilter("hltSingleMu13L3Filtered13").empty()))
// 	     || (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu5_v*").empty()) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered5").at(0).pt()>13)) 
// 	     || (!(mu->triggerObjectMatchesByPath("HLT_DoubleMu7_v*").empty()) && (!mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").empty()) && (mu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered7").at(0).pt()>13)));
//   }
// }







DEFINE_FWK_MODULE(NtupleProducer);


