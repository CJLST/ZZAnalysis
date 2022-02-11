// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "TH1.h"
#include "TTree.h"

#include <string>
#include <iostream>

using namespace std;
using namespace edm;

class LeptonNtupleMaker : public edm::EDAnalyzer {

public:
  explicit LeptonNtupleMaker(const edm::ParameterSet&);
  ~LeptonNtupleMaker();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------

  edm::EDGetTokenT<pat::ElectronCollection> electronToken;
  edm::EDGetTokenT<pat::MuonCollection> muonToken;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken;
  edm::EDGetTokenT<vector<reco::Vertex> > vtxToken;
  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken;
  TString theFileName;
  TTree* myTree;
  TH1F *hCounter;
  Float_t gen_sumGenMCWeight;

  //event variables
  Int_t nEvent;
  Int_t nRun;
  Int_t nLumi;
  Int_t Nvtx;
  Int_t NObsInt;
  Float_t NTrueInt;
  Float_t genHEPMCweight;

  //gen electrons
  Int_t ge_n;
  Float_t ge_pt[100];
  Float_t ge_eta[100];
  Float_t ge_phi[100];

  //gen muons
  Int_t gm_n;
  Float_t gm_pt[100];
  Float_t gm_eta[100];
  Float_t gm_phi[100];

  //reco electrons  
  Int_t nele;
  Int_t ele_pdgId[100];
  Float_t ele_pt[100];
  Float_t ele_eta[100];
  Float_t ele_sclEta[100];
  Int_t ele_isEB[100];
  Int_t ele_isEE[100];
  Float_t ele_phi[100];
  Float_t ele_dxy[100];
  Float_t ele_dz[100];
  Int_t ele_missingHit[100];
  Float_t ele_BDT[100];
  Int_t ele_isBDT[100];
  Float_t ele_SIP[100];
  Float_t ele_PFChargedHadIso[100];
  Float_t ele_PFNeutralHadIso[100];
  Float_t ele_PFPhotonIso[100];
  Float_t ele_combRelIsoPF[100];
  //Float_t ele_combRelIsoPFNon[100];
  //Float_t ele_combRelIsoPFRho[100];
  //Float_t ele_combRelIsoPFDbt[100];
  //Float_t ele_combRelIsoPFPfw[100];
  //Float_t ele_combRelIsoPFPup[100];
  Float_t ele_rho[100];
  Int_t ele_isLooseEle[100];
  Int_t ele_isID[100];
  Int_t ele_isSIP[100];
  Int_t ele_isIsoFSRUncorr[100];
  Int_t ele_isTightEle[100];
  Int_t ele_isGood[100];
  Int_t ele_isFullSel[100];

  Float_t fMVAVar_kfhits[100];
  Float_t fMVAVar_see[100];
  Float_t fMVAVar_spp[100];
  Float_t fMVAVar_OneMinusE1x5E5x5[100];
  Float_t fMVAVar_R9[100];
  Float_t fMVAVar_etawidth[100];
  Float_t fMVAVar_phiwidth[100];
  Float_t fMVAVar_HoE[100];
  Float_t fMVAVar_PreShowerOverRaw[100];
  Float_t fMVAVar_kfchi2[100];
  Float_t fMVAVar_gsfchi2[100];
  Float_t fMVAVar_fbrem[100];
  Float_t fMVAVar_EoP[100];
  Float_t fMVAVar_eleEoPout[100];
  Float_t fMVAVar_IoEmIoP[100];
  Float_t fMVAVar_deta[100];
  Float_t fMVAVar_dphi[100];
  Float_t fMVAVar_detacalo[100]; 

  //reco muons  
  Int_t nmu;
  Int_t mu_pdgId[100];
  Float_t mu_pt[100];
  Float_t mu_eta[100];
  Float_t mu_phi[100];
  Float_t mu_dxy[100];
  Float_t mu_dz[100];
  Int_t mu_isGlobalMuon[100];
  Int_t mu_isTrackerMuon[100];
  Int_t mu_numberOfMatches[100];
  Int_t mu_muonBestTrackType[100];
  Int_t mu_isPFMuon[100];
  Float_t mu_SIP[100];
  Float_t mu_PFChargedHadIso[100];
  Float_t mu_PFNeutralHadIso[100];
  Float_t mu_PFPhotonIso[100];
  Float_t mu_combRelIsoPF[100];
  //Float_t mu_combRelIsoPFNon[100];
  //Float_t mu_combRelIsoPFRho[100];
  //Float_t mu_combRelIsoPFDbt[100];
  //Float_t mu_combRelIsoPFPfw[100];
  //Float_t mu_combRelIsoPFPup[100];
  Float_t mu_rho[100];
  Int_t mu_isLooseMu[100];
  Int_t mu_isID[100];
  Int_t mu_isSIP[100];
  Int_t mu_isIsoFSRUncorr[100];
  Int_t mu_isTightMu[100];
  Int_t mu_isGood[100];
  Int_t mu_isFullSel[100];

};


LeptonNtupleMaker::LeptonNtupleMaker(const edm::ParameterSet& iConfig) :
  electronToken(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronSrc"))),
  muonToken(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
  genParticleToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleSrc")))
{
  theFileName = iConfig.getUntrackedParameter<string>("fileName"); 

  vtxToken = consumes<vector<reco::Vertex> >(edm::InputTag("goodPrimaryVertices"));
  genInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  consumesMany<std::vector< PileupSummaryInfo > >();

}


LeptonNtupleMaker::~LeptonNtupleMaker() {}


bool isGoodGenLepton(reco::GenParticleRef p, int leptonType){

  if(fabs(p->pdgId())==leptonType){

    if(p->status()==1){
      if(p->numberOfMothers()>0){
	const reco::Candidate *Mom = p->mother();
	while(Mom!=0){
	  if(fabs(Mom->pdgId()) == 23){ 
	    return true;
	  }else if(fabs(Mom->pdgId())==15 || fabs(Mom->pdgId())>50){ 
	    break;
	  }
	  Mom = Mom->mother();
	}
      }
    }

  }

  return false;

}


// ------------ method called for each event  ------------
void
LeptonNtupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  //----- Event-level information

  nEvent = iEvent.id().event();
  nRun   = iEvent.id().run();
  nLumi  = iEvent.luminosityBlock();

  vector<Handle<std::vector< PileupSummaryInfo > > > PupInfos;
  iEvent.getManyByType(PupInfos);
  Handle<std::vector< PileupSummaryInfo > > PupInfo = PupInfos.front();   
  std::vector<PileupSummaryInfo>::const_iterator PVI;
  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
    if(PVI->getBunchCrossing() == 0) { 
      NObsInt  = PVI->getPU_NumInteractions();
      NTrueInt = PVI->getTrueNumInteractions();
      break;
    } 
  }

  // Generator weigths
  edm::Handle<GenEventInfoProduct> gen;
  iEvent.getByToken(genInfoToken, gen);
  GenEventInfoProduct  genInfo = *(gen.product());
  genHEPMCweight = genInfo.weight();
  gen_sumGenMCWeight += genHEPMCweight;

  // Primary vertices
  Handle<vector<reco::Vertex> > vertices;
  iEvent.getByToken(vtxToken,vertices);
  Nvtx = vertices->size();


  //----- Gen leptons

  edm::Handle<reco::GenParticleCollection> gpH;
  iEvent.getByToken(genParticleToken, gpH);
  //gp_n = 0;
  ge_n = 0;
  gm_n = 0;
  for (size_t i=0; i<gpH->size(); ++i) {
    const reco::GenParticleRef gp(gpH, i);
    if(isGoodGenLepton(gp,11)){
      ge_pt [ge_n] = gp->pt ();
      ge_eta[ge_n] = gp->eta();
      ge_phi[ge_n] = gp->phi();
      ge_n++;
    }
    if(isGoodGenLepton(gp,13)){
      gm_pt [gm_n] = gp->pt ();
      gm_eta[gm_n] = gp->eta();
      gm_phi[gm_n] = gp->phi();
      gm_n++;
    }
  }
    

  //----- Electrons
  
  edm::Handle<pat::ElectronCollection> electronHandle;
  iEvent.getByToken(electronToken, electronHandle);

  nele = 0;
  //cout<<"in LeptonNtupleMaker.cc : electronHandle->size() = "<<electronHandle->size()<<endl;

  for(unsigned int iele=0; iele<electronHandle->size(); iele++){

    ele_pdgId[nele] = -9999;
    ele_pt[nele] = -9999.;
    ele_eta[nele] = -9999.;
    ele_sclEta[nele] = -9999.;
    ele_isEB[nele] = -9999;
    ele_isEE[nele] = -9999;
    ele_phi[nele] = -9999.;
    ele_dxy[nele] = -9999.;
    ele_dz[nele] = -9999.;
    ele_missingHit[nele] = -9999;
    ele_BDT[nele] = -9999.;
    ele_isBDT[nele] = -9999;
    ele_SIP[nele] = -9999.;
    ele_PFChargedHadIso[nele] = -9999.;
    ele_PFNeutralHadIso[nele] = -9999.;
    ele_PFPhotonIso[nele] = -9999.;
    ele_combRelIsoPF[nele] = -9999.;
    //ele_combRelIsoPFNon[nele] = -9999.;
    //ele_combRelIsoPFRho[nele] = -9999.;
    //ele_combRelIsoPFDbt[nele] = -9999.;
    //ele_combRelIsoPFPfw[nele] = -9999.;
    //ele_combRelIsoPFPup[nele] = -9999.;
    ele_rho[nele] = -9999.;
    ele_isLooseEle[nele] = -9999;
    ele_isID[nele] = -9999;
    ele_isSIP[nele] = -9999;
    ele_isIsoFSRUncorr[nele] = -9999;
    ele_isTightEle[nele] = -9999;
    ele_isGood[nele] = -9999;
    ele_isFullSel[nele] = -9999;

    fMVAVar_kfhits[nele] = -9999.;
    fMVAVar_see[nele] = -9999.;
    fMVAVar_spp[nele] = -9999.;
    fMVAVar_OneMinusE1x5E5x5[nele] = -9999.;
    fMVAVar_R9[nele] = -9999.;
    fMVAVar_etawidth[nele] = -9999.;
    fMVAVar_phiwidth[nele] = -9999.;
    fMVAVar_HoE[nele] = -9999.;
    fMVAVar_PreShowerOverRaw[nele] = -9999.;
    fMVAVar_kfchi2[nele] = -9999.;
    fMVAVar_gsfchi2[nele] = -9999.;
    fMVAVar_fbrem[nele] = -9999.;
    fMVAVar_EoP[nele] = -9999.;
    fMVAVar_eleEoPout[nele] = -9999.;
    fMVAVar_IoEmIoP[nele] = -9999.;
    fMVAVar_deta[nele] = -9999.;
    fMVAVar_dphi[nele] = -9999.;
    fMVAVar_detacalo[nele] = -9999.; 
      
    const pat::Electron* e = &((*electronHandle)[iele]);

    if(e->pt()<5.) continue;

    ele_pdgId[nele] = e->pdgId();
    ele_pt[nele] = e->pt(); // ->et() gives same result
    ele_eta[nele] = e->eta();
    ele_sclEta[nele] = e->superCluster()->eta();
    ele_isEB[nele] = e->isEB();
    ele_isEE[nele] = e->isEE();
    ele_phi[nele] = e->phi();
    ele_dxy[nele] = e->dB(pat::Electron::PV2D);
    ele_dz[nele] = e->dB(pat::Electron::PVDZ);
    ele_missingHit[nele] = e->userFloat("missingHit");
    ele_BDT[nele] = e->userFloat("BDT");
    ele_isBDT[nele] = e->userFloat("isBDT");
    ele_SIP[nele] = e->userFloat("SIP");
    ele_PFChargedHadIso[nele] = e->userFloat("PFChargedHadIso");
    ele_PFNeutralHadIso[nele] = e->userFloat("PFNeutralHadIso");
    ele_PFPhotonIso[nele] = e->userFloat("PFPhotonIso");
    ele_combRelIsoPF[nele] = e->userFloat("combRelIsoPF");
    //ele_combRelIsoPFNon[nele] = e->userFloat("combRelIsoPFNon");
    //ele_combRelIsoPFRho[nele] = e->userFloat("combRelIsoPFRho");
    //ele_combRelIsoPFDbt[nele] = e->userFloat("combRelIsoPFDbt");
    //ele_combRelIsoPFPfw[nele] = e->userFloat("combRelIsoPFPfw");
    //ele_combRelIsoPFPup[nele] = e->userFloat("combRelIsoPFPup");
    ele_rho[nele] = e->userFloat("rho");
    ele_isLooseEle[nele] = ele_pt[nele]>7. && fabs(ele_eta[nele])<2.5 && ele_dxy[nele]<0.5 && ele_dz[nele]<1. ;
    ele_isID[nele] = e->userFloat("ID");
    ele_isSIP[nele] = e->userFloat("isSIP");
    ele_isIsoFSRUncorr[nele] = e->userFloat("isIsoFSRUncorr");
    ele_isTightEle[nele] = ele_isLooseEle[nele] && ele_isID[nele];
    ele_isGood[nele] = ele_isTightEle[nele] && ele_isSIP[nele]; //e->userFloat("isGood")
    ele_isFullSel[nele] = ele_isGood[nele] && ele_isIsoFSRUncorr[nele];

    bool validKF= false;
    reco::TrackRef myTrackRef = e->closestCtfTrackRef();
    validKF = (myTrackRef.isAvailable());
    validKF = (myTrackRef.isNonnull());
    fMVAVar_kfhits[nele]    = (validKF) ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1. ;
    fMVAVar_see[nele]       = e->full5x5_sigmaIetaIeta();
    fMVAVar_spp[nele]       = e->full5x5_sigmaIphiIphi();
    fMVAVar_OneMinusE1x5E5x5[nele] = (e->full5x5_e5x5()) !=0. ? 1.-(e->full5x5_e1x5()/e->full5x5_e5x5()) : -1. ;
    fMVAVar_R9[nele]        = e->full5x5_r9();
    fMVAVar_etawidth[nele]  = e->superCluster()->etaWidth();
    fMVAVar_phiwidth[nele]  = e->superCluster()->phiWidth();
    fMVAVar_HoE[nele]       = e->hadronicOverEm();
    fMVAVar_PreShowerOverRaw[nele] = e->superCluster()->preshowerEnergy() / e->superCluster()->rawEnergy();
    fMVAVar_kfchi2[nele]    = (validKF) ? myTrackRef->normalizedChi2() : 0 ;
    fMVAVar_gsfchi2[nele]   = e->gsfTrack()->normalizedChi2();
    fMVAVar_fbrem[nele]     = e->fbrem();
    fMVAVar_EoP[nele]       = e->eSuperClusterOverP();;
    fMVAVar_eleEoPout[nele] = e->eEleClusterOverPout();
    fMVAVar_IoEmIoP[nele]   = (1.0/e->ecalEnergy()) - (1.0 / e->trackMomentumAtVtx().R());
    fMVAVar_deta[nele]      = e->deltaEtaSuperClusterTrackAtVtx();
    fMVAVar_dphi[nele]      = e->deltaPhiSuperClusterTrackAtVtx();
    fMVAVar_detacalo[nele]  = e->deltaEtaSeedClusterTrackAtCalo(); 

    nele++;

  }


  //----- Muons
  
  edm::Handle<pat::MuonCollection> muonHandle;
  iEvent.getByToken(muonToken, muonHandle);

  nmu = 0;
  //cout<<"in LeptonNtupleMaker.cc : muonHandle->size() = "<<muonHandle->size()<<endl;

  for(unsigned int imu=0; imu<muonHandle->size(); imu++){

    mu_pdgId[nmu] = -9999;
    mu_pt[nmu] = -9999.;
    mu_eta[nmu] = -9999.;
    mu_phi[nmu] = -9999.;
    mu_dxy[nmu] = -9999.;
    mu_dz[nmu] = -9999.;
    mu_isGlobalMuon[nmu] = -9999;
    mu_isTrackerMuon[nmu] = -9999;
    mu_numberOfMatches[nmu] = -9999;
    mu_muonBestTrackType[nmu] = -9999;
    mu_isPFMuon[nmu] = -9999;
    mu_SIP[nmu] = -9999.;
    mu_PFChargedHadIso[nmu] = -9999.;
    mu_PFNeutralHadIso[nmu] = -9999.;
    mu_PFPhotonIso[nmu] = -9999.;
    mu_combRelIsoPF[nmu] = -9999.;
    //mu_combRelIsoPFNon[nmu] = -9999.;
    //mu_combRelIsoPFRho[nmu] = -9999.;
    //mu_combRelIsoPFDbt[nmu] = -9999.;
    //mu_combRelIsoPFPfw[nmu] = -9999.;
    //mu_combRelIsoPFPup[nmu] = -9999.;
    mu_rho[nmu] = -9999.;
    mu_isLooseMu[nmu] = -9999;
    mu_isID[nmu] = -9999;
    mu_isSIP[nmu] = -9999;
    mu_isIsoFSRUncorr[nmu] = -9999;
    mu_isTightMu[nmu] = -9999;
    mu_isGood[nmu] = -9999;
    mu_isFullSel[nmu] = -9999;
      
    const pat::Muon* m = &((*muonHandle)[imu]);

    if(m->pt()<5.) continue;

    mu_pdgId[nmu] = m->pdgId();
    mu_pt[nmu] = m->pt();
    mu_eta[nmu] = m->eta();
    mu_phi[nmu] = m->phi();
    mu_dxy[nmu] = m->dB(pat::Muon::PV2D);
    mu_dz[nmu] = m->dB(pat::Muon::PVDZ);
    mu_isGlobalMuon[nmu] = m->isGlobalMuon();
    mu_isTrackerMuon[nmu] = m->isTrackerMuon();
    mu_numberOfMatches[nmu] = m->numberOfMatches();
    mu_muonBestTrackType[nmu] = m->muonBestTrackType();
    mu_isPFMuon[nmu] = m->isPFMuon();
    mu_SIP[nmu] = m->userFloat("SIP");
    mu_PFChargedHadIso[nmu] = m->userFloat("PFChargedHadIso");
    mu_PFNeutralHadIso[nmu] = m->userFloat("PFNeutralHadIso");
    mu_PFPhotonIso[nmu] = m->userFloat("PFPhotonIso");
    mu_combRelIsoPF[nmu] = m->userFloat("combRelIsoPF");
    //mu_combRelIsoPFNon[nmu] = m->userFloat("combRelIsoPFNon");
    //mu_combRelIsoPFRho[nmu] = m->userFloat("combRelIsoPFRho");
    //mu_combRelIsoPFDbt[nmu] = m->userFloat("combRelIsoPFDbt");
    //mu_combRelIsoPFPfw[nmu] = m->userFloat("combRelIsoPFPfw");
    //mu_combRelIsoPFPup[nmu] = m->userFloat("combRelIsoPFPup");
    mu_rho[nmu] = m->userFloat("rho");
    mu_isLooseMu[nmu] = mu_pt[nmu]>5. && fabs(mu_eta[nmu])<2.4 && mu_dxy[nmu]<0.5 && mu_dz[nmu]<1. && (mu_isGlobalMuon[nmu] || (mu_isTrackerMuon[nmu] && mu_numberOfMatches[nmu]>0)) && mu_muonBestTrackType[nmu]!=2;
    mu_isID[nmu] = m->userFloat("ID");
    mu_isSIP[nmu] = m->userFloat("isSIP");
    mu_isIsoFSRUncorr[nmu] = m->userFloat("isIsoFSRUncorr");
    mu_isTightMu[nmu] = mu_isLooseMu[nmu] && mu_isID[nmu];
    mu_isGood[nmu] = mu_isTightMu[nmu] && mu_isSIP[nmu]; //m->userFloat("isGood")
    mu_isFullSel[nmu] = mu_isGood[nmu] && mu_isIsoFSRUncorr[nmu];

    nmu++;

  }
  
  myTree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
LeptonNtupleMaker::beginJob(){

  edm::Service<TFileService> fs;

  myTree = fs->make<TTree>(theFileName,"The tree");
  
  myTree->Branch("nEvent", &nEvent, "nEvent/I");
  myTree->Branch("nRun", &nRun, "nRun/I");
  myTree->Branch("nLumi", &nLumi, "nLumi/I");
  myTree->Branch("Nvtx", &Nvtx, "Nvtx/I");
  myTree->Branch("NObsInt", &NObsInt, "NObsInt/I");
  myTree->Branch("NTrueInt", &NTrueInt, "NTrueInt/F");
  myTree->Branch("genHEPMCweight", &genHEPMCweight, "genHEPMCweight/F");

  myTree->Branch("gen", &ge_n, "gen/I");
  myTree->Branch("gept", &ge_pt, "gept[gen]/F");
  myTree->Branch("geeta", &ge_eta, "geeta[gen]/F");
  myTree->Branch("gephi", &ge_phi, "gephi[gen]/F");
  myTree->Branch("gmn", &gm_n, "gmn/I");
  myTree->Branch("gmpt", &gm_pt, "gmpt[gmn]/F");
  myTree->Branch("gmeta", &gm_eta, "gmeta[gmn]/F");
  myTree->Branch("gmphi", &gm_phi, "gmphi[gmn]/F");
  
  myTree->Branch("nele", &nele, "nele/I");
  myTree->Branch("ele_pdgId", &ele_pdgId, "ele_pdgId[nele]/I");
  myTree->Branch("ele_pt" , &ele_pt , "ele_pt[nele]/F" );
  myTree->Branch("ele_eta", &ele_eta, "ele_eta[nele]/F");
  myTree->Branch("ele_sclEta", &ele_sclEta, "ele_sclEta[nele]/F");
  myTree->Branch("ele_isEB", &ele_isEB, "ele_isEB[nele]/I");
  myTree->Branch("ele_isEE", &ele_isEE, "ele_isEE[nele]/I");
  myTree->Branch("ele_phi", &ele_phi, "ele_phi[nele]/F");
  myTree->Branch("ele_dxy", &ele_dxy, "ele_dxy[nele]/F");
  myTree->Branch("ele_dz", &ele_dz, "ele_dz[nele]/F");
  myTree->Branch("ele_missingHit", &ele_missingHit, "ele_missingHit[nele]/I");
  myTree->Branch("ele_BDT", &ele_BDT, "ele_BDT[nele]/F");
  myTree->Branch("ele_isBDT", &ele_isBDT, "ele_isBDT[nele]/I");
  myTree->Branch("ele_SIP", &ele_SIP, "ele_SIP[nele]/F");
  myTree->Branch("ele_PFChargedHadIso", &ele_PFChargedHadIso, "ele_PFChargedHadIso[nele]/F");
  myTree->Branch("ele_PFNeutralHadIso", &ele_PFNeutralHadIso, "ele_PFNeutralHadIso[nele]/F");
  myTree->Branch("ele_PFPhotonIso", &ele_PFPhotonIso, "ele_PFPhotonIso[nele]/F");
  myTree->Branch("ele_combRelIsoPF", &ele_combRelIsoPF, "ele_combRelIsoPF[nele]/F");
  //myTree->Branch("ele_combRelIsoPFNon", &ele_combRelIsoPFNon, "ele_combRelIsoPFNon[nele]/F");
  //myTree->Branch("ele_combRelIsoPFRho", &ele_combRelIsoPFRho, "ele_combRelIsoPFRho[nele]/F");
  //myTree->Branch("ele_combRelIsoPFDbt", &ele_combRelIsoPFDbt, "ele_combRelIsoPFDbt[nele]/F");
  //myTree->Branch("ele_combRelIsoPFPfw", &ele_combRelIsoPFPfw, "ele_combRelIsoPFPfw[nele]/F");
  //myTree->Branch("ele_combRelIsoPFPup", &ele_combRelIsoPFPup, "ele_combRelIsoPFPup[nele]/F");
  myTree->Branch("ele_rho", &ele_rho, "ele_rho[nele]/F");
  myTree->Branch("ele_isLooseEle", &ele_isLooseEle, "ele_isLooseEle[nele]/I");
  myTree->Branch("ele_isID", &ele_isID, "ele_isID[nele]/I");
  myTree->Branch("ele_isSIP", &ele_isSIP, "ele_isSIP[nele]/I");
  myTree->Branch("ele_isIsoFSRUncorr", &ele_isIsoFSRUncorr, "ele_isIsoFSRUncorr[nele]/I");
  myTree->Branch("ele_isTightEle", &ele_isTightEle, "ele_isTightEle[nele]/I");
  myTree->Branch("ele_isGood", &ele_isGood, "ele_isGood[nele]/I");
  myTree->Branch("ele_isFullSel", &ele_isFullSel, "ele_isFullSel[nele]/I");

  myTree->Branch("ele_kfhits", &fMVAVar_kfhits, "ele_kfhits[nele]/F");
  myTree->Branch("ele_oldsigmaietaieta", &fMVAVar_see, "ele_oldsigmaietaieta[nele]/F");
  myTree->Branch("ele_oldsigmaiphiiphi", &fMVAVar_spp, "ele_oldsigmaiphiiphi[nele]/F");
  myTree->Branch("ele_oldcircularity", &fMVAVar_OneMinusE1x5E5x5, "ele_oldcircularity[nele]/F");
  myTree->Branch("ele_oldr9", &fMVAVar_R9, "ele_oldr9[nele]/F");
  myTree->Branch("ele_scletawidth", &fMVAVar_etawidth, "ele_scletawidth[nele]/F");
  myTree->Branch("ele_sclphiwidth", &fMVAVar_phiwidth, "ele_sclphiwidth[nele]/F");
  myTree->Branch("ele_he", &fMVAVar_HoE, "ele_he[nele]/F");
  myTree->Branch("ele_psEoverEraw", &fMVAVar_PreShowerOverRaw, "ele_psEoverEraw[nele]/F");
  myTree->Branch("ele_kfchi2", &fMVAVar_kfchi2, "ele_kfchi2[nele]/F");
  myTree->Branch("ele_chi2_hits", &fMVAVar_gsfchi2, "ele_chi2_hits[nele]/F");
  myTree->Branch("ele_fbrem", &fMVAVar_fbrem, "ele_fbrem[nele]/F");
  myTree->Branch("ele_ep", &fMVAVar_EoP, "ele_ep[nele]/F");
  myTree->Branch("ele_eelepout", &fMVAVar_eleEoPout, "ele_eelepout[nele]/F");
  myTree->Branch("ele_IoEmIop", &fMVAVar_IoEmIoP, "ele_IoEmIop[nele]/F");
  myTree->Branch("ele_deltaetain", &fMVAVar_deta, "ele_deltaetain[nele]/F");
  myTree->Branch("ele_deltaphiin", &fMVAVar_dphi, "ele_deltaphiin[nele]/F");
  myTree->Branch("ele_deltaetaseed", &fMVAVar_detacalo, "ele_deltaetaseed[nele]/F");

  myTree->Branch("nmu", &nmu, "nmu/I");
  myTree->Branch("mu_pdgId", &mu_pdgId, "mu_pdgId[nmu]/I");
  myTree->Branch("mu_pt" , &mu_pt , "mu_pt[nmu]/F" );
  myTree->Branch("mu_eta", &mu_eta, "mu_eta[nmu]/F");
  myTree->Branch("mu_phi", &mu_phi, "mu_phi[nmu]/F");
  myTree->Branch("mu_dxy", &mu_dxy, "mu_dxy[nmu]/F");
  myTree->Branch("mu_dz", &mu_dz, "mu_dz[nmu]/F");
  myTree->Branch("mu_isGlobalMuon", &mu_isGlobalMuon, "mu_isGlobalMuon[nmu]/I");
  myTree->Branch("mu_isTrackerMuon", &mu_isTrackerMuon, "mu_isTrackerMuon[nmu]/I");
  myTree->Branch("mu_numberOfMatches", &mu_numberOfMatches, "mu_numberOfMatches[nmu]/I");
  myTree->Branch("mu_muonBestTrackType", &mu_muonBestTrackType, "mu_muonBestTrackType[nmu]/I");
  myTree->Branch("mu_isPFMuon", &mu_isPFMuon, "mu_isPFMuon[nmu]/I");
  myTree->Branch("mu_SIP", &mu_SIP, "mu_SIP[nmu]/F");
  myTree->Branch("mu_PFChargedHadIso", &mu_PFChargedHadIso, "mu_PFChargedHadIso[nmu]/F");
  myTree->Branch("mu_PFNeutralHadIso", &mu_PFNeutralHadIso, "mu_PFNeutralHadIso[nmu]/F");
  myTree->Branch("mu_PFPhotonIso", &mu_PFPhotonIso, "mu_PFPhotonIso[nmu]/F");
  myTree->Branch("mu_combRelIsoPF", &mu_combRelIsoPF, "mu_combRelIsoPF[nmu]/F");
  //myTree->Branch("mu_combRelIsoPFNon", &mu_combRelIsoPFNon, "mu_combRelIsoPFNon[nmu]/F");
  //myTree->Branch("mu_combRelIsoPFRho", &mu_combRelIsoPFRho, "mu_combRelIsoPFRho[nmu]/F");
  //myTree->Branch("mu_combRelIsoPFDbt", &mu_combRelIsoPFDbt, "mu_combRelIsoPFDbt[nmu]/F");
  //myTree->Branch("mu_combRelIsoPFPfw", &mu_combRelIsoPFPfw, "mu_combRelIsoPFPfw[nmu]/F");
  //myTree->Branch("mu_combRelIsoPFPup", &mu_combRelIsoPFPup, "mu_combRelIsoPFPup[nmu]/F");
  myTree->Branch("mu_rho", &mu_rho, "mu_rho[nmu]/F");
  myTree->Branch("mu_isLooseMu", &mu_isLooseMu, "mu_isLooseMu[nmu]/I");
  myTree->Branch("mu_isID", &mu_isID, "mu_isID[nmu]/I");
  myTree->Branch("mu_isSIP", &mu_isSIP, "mu_isSIP[nmu]/I");
  myTree->Branch("mu_isIsoFSRUncorr", &mu_isIsoFSRUncorr, "mu_isIsoFSRUncorr[nmu]/I");
  myTree->Branch("mu_isTightMu", &mu_isTightMu, "mu_isTightMu[nmu]/I");
  myTree->Branch("mu_isGood", &mu_isGood, "mu_isGood[nmu]/I");
  myTree->Branch("mu_isFullSel", &mu_isFullSel, "mu_isFullSel[nmu]/I");

  hCounter = fs->make<TH1F>("Counters", "Counters", 5, 0., 5.);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
LeptonNtupleMaker::endJob() {

  hCounter->SetBinContent(1,gen_sumGenMCWeight);
  hCounter->GetXaxis()->SetBinLabel(1,"gen_sumGenMCWeight");

  return;
}

// ------------ method called when starting to processes a run  ------------
/*
void 
LeptonNtupleMaker::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
LeptonNtupleMaker::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
LeptonNtupleMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
LeptonNtupleMaker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LeptonNtupleMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LeptonNtupleMaker);
