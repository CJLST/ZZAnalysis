// -*- C++ -*-
//
// Fill a tree for selected candidates.
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
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <DataFormats/JetReco/interface/PFJet.h>
#include <DataFormats/JetReco/interface/PFJetCollection.h>
#include <DataFormats/Math/interface/LorentzVector.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <DataFormats/Common/interface/MergeableCounter.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include <SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h>
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>

#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>
#include <ZZAnalysis/AnalysisStep/interface/MCHistoryTools.h>
#include <ZZAnalysis/AnalysisStep/interface/PUReweight.h>
#include "ZZAnalysis/AnalysisStep/interface/EwkCorrections.h"
#include "ZZAnalysis/AnalysisStep/interface/reweighting.h"
#include "ZZAnalysis/AnalysisStep/interface/LHEHandler.h"
#include "ZZAnalysis/AnalysisStep/src/kFactors.C"
#include <ZZAnalysis/AnalysisStep/interface/bitops.h>
#include <ZZAnalysis/AnalysisStep/interface/Fisher.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>
#include <ZZAnalysis/AnalysisStep/interface/JetCleaner.h>
#include <ZZAnalysis/AnalysisStep/interface/utils.h>

#include <ZZMatrixElement/MELA/interface/Mela.h>

#include "ZZ4lConfigHelper.h"
#include "HZZ4lNtupleFactory.h"

#include <TRandom3.h>
#include <TH2D.h>
#include "TLorentzVector.h"
#include "TSpline.h"

#include <string>

namespace {
  bool writeJets = true;     // Write jets in the tree. FIXME: make this configurable
  bool addKinRefit = true;
  bool addVtxFit = false;
  bool addFSRDetails = false;
  bool addQGLInputs = true;
  bool skipMuDataMCWeight = false; // skip computation of data/MC weight for mu
  bool skipEleDataMCWeight = false; // skip computation of data/MC weight for ele
  bool skipFakeWeight = true;   // skip computation of fake rate weight for CRs
  bool skipHqTWeight = true;    // skip computation of hQT weight

  //List of variables with default values
  Int_t RunNumber  = 0;
  Long64_t EventNumber  = 0;
  Int_t LumiNumber  = 0;
  Short_t NRecoMu  = 0;
  Short_t NRecoEle  = 0;
  Short_t Nvtx  = 0;
  Short_t NObsInt  = 0;
  Float_t NTrueInt  = 0;
  Float_t PUWeight  = 0;
  Float_t KFactor_QCD_ggZZ_Nominal = 0;
  Float_t KFactor_QCD_ggZZ_PDFScaleDn = 0;
  Float_t KFactor_QCD_ggZZ_PDFScaleUp = 0;
  Float_t KFactor_QCD_ggZZ_QCDScaleDn = 0;
  Float_t KFactor_QCD_ggZZ_QCDScaleUp = 0;
  Float_t KFactor_QCD_ggZZ_AsDn = 0;
  Float_t KFactor_QCD_ggZZ_AsUp = 0;
  Float_t KFactor_QCD_ggZZ_PDFReplicaDn = 0;
  Float_t KFactor_QCD_ggZZ_PDFReplicaUp = 0;
  Float_t KFactor_EW_qqZZ = 0;
  Float_t KFactor_EW_qqZZ_unc = 0;
  Float_t KFactor_QCD_qqZZ_dPhi = 0;
  Float_t KFactor_QCD_qqZZ_M = 0;
  Float_t KFactor_QCD_qqZZ_Pt = 0;
  Float_t PFMET  =  -99;
  Float_t PFMETPhi  =  -99;
  Float_t PFMETNoHF  =  -99;
  Float_t PFMETNoHFPhi  =  -99;
  Short_t nCleanedJets  =  0;
  Short_t nCleanedJetsPt30  = 0;
  Short_t nCleanedJetsPt30_jecUp  = 0;
  Short_t nCleanedJetsPt30_jecDn  = 0;
  Short_t nCleanedJetsPt30BTagged  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSFUp  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSFDn  = 0;
  Short_t trigWord  = 0;
  Float_t ZZMass  = 0;
  Float_t ZZMassErr  = 0;
  Float_t ZZMassErrCorr  = 0;
  Float_t ZZMassPreFSR  = 0;
  Short_t ZZsel  = 0;
  Float_t ZZPt  = 0;
  Float_t ZZEta  = 0;
  Float_t ZZPhi  = 0;
  Int_t CRflag  = 0;
  Float_t Z1Mass  = 0;
  Float_t Z1Pt  = 0;
  Short_t Z1Flav  = 0;
  Float_t ZZMassRefit  = 0;
  Float_t ZZMassRefitErr  = 0;
  Float_t ZZMassUnrefitErr  = 0;
  Float_t ZZMassCFit  = 0;
  Float_t ZZChi2CFit  = 0;
  Float_t Z2Mass  = 0;
  Float_t Z2Pt  = 0;
  Short_t Z2Flav  = 0;
  Float_t costhetastar  = 0;
  Float_t helphi  = 0;
  Float_t helcosthetaZ1  = 0;
  Float_t helcosthetaZ2  = 0;
  Float_t phistarZ1  = 0;
  Float_t phistarZ2  = 0;
  Float_t xi  = 0;
  Float_t xistar  = 0;
  Float_t TLE_dR_Z = -1; // Delta-R between a TLE and the Z it does not belong to.
  Float_t TLE_min_dR_3l = 999; // Minimum DR between a TLE and any of the other leptons

  std::vector<float> LepPt;
  std::vector<float> LepEta;
  std::vector<float> LepPhi;
  std::vector<short> LepLepId;
  std::vector<float> LepSIP;
  std::vector<float> LepTime;
  std::vector<bool> LepisID;
  std::vector<float> LepBDT;
  std::vector<char> LepMissingHit;
  //std::vector<float> LepChargedHadIso;
  //std::vector<float> LepNeutralHadIso;
  //std::vector<float> LepPhotonIso;
  std::vector<float> LepCombRelIsoPF;
  std::vector<short> LepisLoose;
  std::vector<float> fsrPt;
  std::vector<float> fsrEta;
  std::vector<float> fsrPhi;
  std::vector<float> fsrDR;
  std::vector<short> fsrLept;
  std::vector<short> fsrLeptID;
  std::vector<float> fsrGenPt;
  Bool_t passIsoPreFSR = 0;

  std::vector<float> JetPt ;
  std::vector<float> JetEta ;
  std::vector<float> JetPhi ;
  std::vector<float> JetMass ;
  std::vector<float> JetBTagger ;
  std::vector<float> JetIsBtagged;
  std::vector<float> JetIsBtaggedWithSF;
  std::vector<float> JetIsBtaggedWithSFUp;
  std::vector<float> JetIsBtaggedWithSFDn;
  std::vector<float> JetQGLikelihood;
  std::vector<float> JetAxis2;
  std::vector<float> JetMult;
  std::vector<float> JetPtD;
  std::vector<float> JetSigma ;
  std::vector<short> JetHadronFlavour;
  Float_t DiJetMass  = -99;
//   Float_t DiJetMassPlus  = -99;
//   Float_t DiJetMassMinus  = -99;
  Float_t DiJetDEta  = -99;
  Float_t DiJetFisher  = -99;
  Short_t nExtraLep  = 0;
  Short_t nExtraZ  = 0;
  std::vector<float> ExtraLepPt;
  std::vector<float> ExtraLepEta;
  std::vector<float> ExtraLepPhi ;
  std::vector<short> ExtraLepLepId;
  Short_t genFinalState  = 0;
  Int_t genProcessId  = 0;
  Float_t genHEPMCweight  = 0;
  std::vector<float> reweightingweights;

  std::vector<float> LHEMotherPz;
  std::vector<float> LHEMotherE;
  std::vector<short> LHEMotherId;
  std::vector<float> LHEDaughterPt;
  std::vector<float> LHEDaughterEta;
  std::vector<float> LHEDaughterPhi;
  std::vector<float> LHEDaughterMass;
  std::vector<short> LHEDaughterId;
  std::vector<float> LHEAssociatedParticlePt;
  std::vector<float> LHEAssociatedParticleEta;
  std::vector<float> LHEAssociatedParticlePhi;
  std::vector<float> LHEAssociatedParticleMass;
  std::vector<short> LHEAssociatedParticleId;

  Float_t LHEPDFScale = 0;
  Float_t LHEweight_QCDscale_muR1_muF1  = 0;
  Float_t LHEweight_QCDscale_muR1_muF2  = 0;
  Float_t LHEweight_QCDscale_muR1_muF0p5  = 0;
  Float_t LHEweight_QCDscale_muR2_muF1  = 0;
  Float_t LHEweight_QCDscale_muR2_muF2  = 0;
  Float_t LHEweight_QCDscale_muR2_muF0p5  = 0;
  Float_t LHEweight_QCDscale_muR0p5_muF1  = 0;
  Float_t LHEweight_QCDscale_muR0p5_muF2  = 0;
  Float_t LHEweight_QCDscale_muR0p5_muF0p5  = 0;
  Float_t LHEweight_PDFVariation_Up = 0;
  Float_t LHEweight_PDFVariation_Dn = 0;
  Float_t LHEweight_AsMZ_Up = 0;
  Float_t LHEweight_AsMZ_Dn = 0;

  Short_t genExtInfo  = 0;
  Float_t xsection  = 0;
  Float_t dataMCWeight  = 0;
  Float_t trigEffWeight  = 0;
  Float_t HqTMCweight  = 0;
  Float_t ZXFakeweight  = 0;
  Float_t overallEventWeight  = 0;
  Float_t GenHMass  = 0;
  Float_t GenHPt  = 0;
  Float_t GenHRapidity  = 0;
  Float_t GenZ1Mass  = 0;
  Float_t GenZ1Eta  = 0;
  Float_t GenZ1Pt  = 0;
  Float_t GenZ1Phi  = 0;
  Float_t GenZ1Flav  = 0;
  Float_t GenZ2Mass  = 0;
  Float_t GenZ2Eta  = 0;
  Float_t GenZ2Pt  = 0;
  Float_t GenZ2Phi  = 0;
  Float_t GenZ2Flav  = 0;
  Float_t GenLep1Pt  = 0;
  Float_t GenLep1Eta  = 0;
  Float_t GenLep1Phi  = 0;
  Short_t GenLep1Id  = 0;
  Float_t GenLep2Pt  = 0;
  Float_t GenLep2Eta  = 0;
  Float_t GenLep2Phi  = 0;
  Short_t GenLep2Id  = 0;
  Float_t GenLep3Pt  = 0;
  Float_t GenLep3Eta  = 0;
  Float_t GenLep3Phi  = 0;
  Short_t GenLep3Id  = 0;
  Float_t GenLep4Pt  = 0;
  Float_t GenLep4Eta  = 0;
  Float_t GenLep4Phi  = 0;
  Short_t GenLep4Id  = 0;
  Float_t GenAssocLep1Pt  = 0;
  Float_t GenAssocLep1Eta  = 0;
  Float_t GenAssocLep1Phi  = 0;
  Short_t GenAssocLep1Id  = 0;
  Float_t GenAssocLep2Pt  = 0;
  Float_t GenAssocLep2Eta  = 0;
  Float_t GenAssocLep2Phi  = 0;
  Short_t GenAssocLep2Id  = 0;


//FIXME: temporary fix to the mismatch of charge() and sign(pdgId()) for muons with BTT=4
  int getPdgId(const reco::Candidate* p) {
    int id = p->pdgId();
    if (id!=22 && //for TLEs
	signbit(id) && p->charge()<0) id*=-1; // negative pdgId must be positive charge
    return id;
  }

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
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  void BookAllBranches();
  virtual void FillKFactors(edm::Handle<GenEventInfoProduct>& genInfo, std::vector<const reco::Candidate *>& genZLeps);
  virtual void FillLHECandidate();
  virtual void FillCandidate(const pat::CompositeCandidate& higgs, bool evtPass, const edm::Event&, const Int_t CRflag);
  virtual void FillJet(const pat::Jet& jet);
  virtual void endJob() ;

  void FillHGenInfo(const math::XYZTLorentzVector Hp, float w);
  void FillZGenInfo(Short_t Z1Id, Short_t Z2Id,
                    const math::XYZTLorentzVector pZ1, const math::XYZTLorentzVector pZ2);
  void FillLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, Short_t Lep3Id, Short_t Lep4Id,
    const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2, const math::XYZTLorentzVector Lep3, const math::XYZTLorentzVector Lep4);
  void FillAssocLepGenInfo(std::vector<const reco::Candidate *>& AssocLeps);

  Float_t getAllWeight(const reco::Candidate* Lep) const;
  Float_t getHqTWeight(double mH, double genPt) const;
  Float_t getFakeWeight(Float_t LepPt, Float_t LepEta, Int_t LepID, Int_t LepZ1ID);

  void addweight(float &weight, vector<float> &weight_reweighted, float weighttoadd);

  void getCheckedUserFloat(const pat::CompositeCandidate& cand, const std::string& strval, Float_t& setval, Float_t defaultval=0);

  void buildMELABranches();
  void addToMELACluster(MELAComputation* me_computer, std::vector<MELACluster*>& me_clusters);
  void computeMELABranches(MELACandidate* cand);
  void updateMELAClusters_Common();
  void pushRecoMELABranches(const pat::CompositeCandidate& cand);
  void pushLHEMELABranches();
  void clearMELABranches();

  // ----------member data ---------------------------
  ZZ4lConfigHelper myHelper;
  int theChannel;
  std::string theCandLabel;
  TString theFileName;

  HZZ4lNtupleFactory *myTree;
  TH1F *hCounter;
  TH2F *hCounter_reweighted;
  TTree *couplingstree;

  bool isMC;
  bool is_loose_ele_selection; // Collection includes candidates with loose electrons/TLEs
  bool applyTrigger;    // Keep only events passing trigger (if skipEmptyEvents=true)
  bool applySkim;       //   "     "      "         skim (if skipEmptyEvents=true)
  bool skipEmptyEvents; // Skip events whith no selected candidate (otherwise, gen info is preserved for all events; candidates not passing trigger&&skim are flagged with negative ZZsel)
  bool applyTrigEffWeight;// apply trigger efficiency weight (concerns samples where trigger is not applied)
  Float_t xsec;
  int year;
  double sqrts;
  double Hmass;

  Mela mela;
  std::vector<std::string> recoMElist;
  std::vector<MELAOptionParser*> recome_originalopts;
  std::vector<MELAOptionParser*> recome_copyopts;
  std::vector<std::string> lheMElist;
  //std::vector<MELAOptionParser*> lheme_originalopts;
  std::vector<MELAOptionParser*> lheme_copyopts;
  std::vector<MELAHypothesis*> lheme_units;
  std::vector<MELAHypothesis*> lheme_aliased_units;
  std::vector<MELAComputation*> lheme_computers;
  std::vector<MELACluster*> lheme_clusters;

  Reweighting reweighting;
  int nReweightingSamples;
  bool doreweighting;
  ReweightingType reweightingtype;

  bool addLHEKinematics;
  LHEHandler* lheHandler;
  int apply_K_NNLOQCD_ZZGG; // 0: Do not; 1: NNLO/LO; 2: NNLO/NLO; 3: NLO/LO
  bool apply_K_NNLOQCD_ZZQQB;
  bool apply_K_NLOEW_ZZQQB;

  bool addProdAnomalousProbabilities;

  edm::EDGetTokenT<edm::View<reco::Candidate> > genParticleToken;
  edm::Handle<edm::View<reco::Candidate> > genParticles;
  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken;
  edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > candToken;
  edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > lhecandToken;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultToken;
  edm::EDGetTokenT<vector<reco::Vertex> > vtxToken;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetToken;
  edm::EDGetTokenT<pat::METCollection> metToken;
  edm::EDGetTokenT<pat::METCollection> metNoHFToken;
  edm::EDGetTokenT<pat::MuonCollection> muonToken;
  edm::EDGetTokenT<pat::ElectronCollection> electronToken;

  edm::EDGetTokenT<edm::MergeableCounter> preSkimToken;

  PUReweight reweight;

  //counters
  Float_t Nevt_Gen;
  Float_t Nevt_Gen_lumiBlock;

  Float_t gen_ZZ4mu;
  Float_t gen_ZZ4e;
  Float_t gen_ZZ2mu2e;
  Float_t gen_ZZ2l2tau;
  Float_t gen_ZZ2emu2tau;
  Float_t gen_ZZ4tau;
  Float_t gen_ZZ4mu_EtaAcceptance;
  Float_t gen_ZZ4mu_LeptonAcceptance;
  Float_t gen_ZZ4e_EtaAcceptance;
  Float_t gen_ZZ4e_LeptonAcceptance;
  Float_t gen_ZZ2mu2e_EtaAcceptance;
  Float_t gen_ZZ2mu2e_LeptonAcceptance;
  Float_t gen_BUGGY;
  Float_t gen_Unknown;

  Float_t gen_sumPUWeight;
  Float_t gen_sumGenMCWeight;
  Float_t gen_sumWeights;

  std::vector<Float_t> Nevt_Gen_reweighted;
  std::vector<Float_t> Nevt_Gen_lumiBlock_reweighted;

  std::vector<Float_t> gen_ZZ4mu_reweighted;
  std::vector<Float_t> gen_ZZ4e_reweighted;
  std::vector<Float_t> gen_ZZ2mu2e_reweighted;
  std::vector<Float_t> gen_ZZ2l2tau_reweighted;
  std::vector<Float_t> gen_ZZ2emu2tau_reweighted;
  std::vector<Float_t> gen_ZZ4tau_reweighted;
  std::vector<Float_t> gen_ZZ4mu_EtaAcceptance_reweighted;
  std::vector<Float_t> gen_ZZ4mu_LeptonAcceptance_reweighted;
  std::vector<Float_t> gen_ZZ4e_EtaAcceptance_reweighted;
  std::vector<Float_t> gen_ZZ4e_LeptonAcceptance_reweighted;
  std::vector<Float_t> gen_ZZ2mu2e_EtaAcceptance_reweighted;
  std::vector<Float_t> gen_ZZ2mu2e_LeptonAcceptance_reweighted;
  std::vector<Float_t> gen_BUGGY_reweighted;
  std::vector<Float_t> gen_Unknown_reweighted;

  std::vector<Float_t> gen_sumGenMCWeight_reweighted;
  std::vector<Float_t> gen_sumWeights_reweighted;

  string sampleName;

  std::vector<const reco::Candidate *> genFSR;

  std::vector<std::vector<float> > ewkTable;
  TSpline3* spkfactor_ggzz_nnlo[9]; // Nominal, PDFScaleDn, PDFScaleUp, QCDScaleDn, QCDScaleUp, AsDn, AsUp, PDFReplicaDn, PDFReplicaUp
  TSpline3* spkfactor_ggzz_nlo[9]; // Nominal, PDFScaleDn, PDFScaleUp, QCDScaleDn, QCDScaleUp, AsDn, AsUp, PDFReplicaDn, PDFReplicaUp

  TH2D *hTH2D_Mu_All;
  TH2F *hTH2F_El_Reco;
  TH1 *hTH2D_El_IdIsoSip_notCracks;
  TH1 *hTH2D_El_IdIsoSip_Cracks;
  TH2D* h_weight; //HqT weights
  //TH2F *h_ZXWeightMuo;
  //TH2F *h_ZXWeightEle;
  TH2D* h_ZXWeight[4];

};

//
// constructors and destructor
//
HZZ4lNtupleMaker::HZZ4lNtupleMaker(const edm::ParameterSet& pset) :
  myHelper(pset),
  theChannel(myHelper.channel()), // Valid options: ZZ, ZLL, ZL
  theCandLabel(pset.getUntrackedParameter<string>("CandCollection")), // Name of input ZZ collection
  theFileName(pset.getUntrackedParameter<string>("fileName")),
  skipEmptyEvents(pset.getParameter<bool>("skipEmptyEvents")), // Do not store
  applyTrigEffWeight(pset.getParameter<bool>("applyTrigEff")),
  xsec(pset.getParameter<double>("xsec")),
  year(pset.getParameter<int>("setup")),
  sqrts(SetupToSqrts(year)),
  Hmass(pset.getParameter<double>("superMelaMass")),
  mela(sqrts, Hmass, TVar::ERROR),
  recoMElist(pset.getParameter<std::vector<std::string>>("recoProbabilities")),
  lheMElist(pset.getParameter<std::vector<std::string>>("lheProbabilities")),
  reweighting(
              mela,
              pset.getParameter<std::string>("reweightingtype"),
              pset.getParameter<int>("spin"),
              pset.getParameter<std::vector<double> >("Hzzcouplings_real"),
              pset.getParameter<std::vector<double> >("Hzzcouplings_imag"),
              pset.getParameter<std::vector<double> >("Zvvcouplings_real"),
              pset.getParameter<std::vector<double> >("Zvvcouplings_imag"),
              pset.getParameter<std::vector<double> >("Gggcouplings_real"),
              pset.getParameter<std::vector<double> >("Gggcouplings_imag"),
              pset.getParameter<std::vector<double> >("Gvvcouplings_real"),
              pset.getParameter<std::vector<double> >("Gvvcouplings_imag"),
              pset.getParameter<std::vector<double> >("reweightingcutoffs")
             ),
  nReweightingSamples(reweighting.nReweightingSamples),
  doreweighting(nReweightingSamples != 0),
  reweightingtype(reweighting.reweightingtype),
  addLHEKinematics(pset.getParameter<bool>("AddLHEKinematics")),
  lheHandler(0),
  apply_K_NNLOQCD_ZZGG(pset.getParameter<int>("Apply_K_NNLOQCD_ZZGG")),
  apply_K_NNLOQCD_ZZQQB(pset.getParameter<bool>("Apply_K_NNLOQCD_ZZQQB")),
  apply_K_NLOEW_ZZQQB(pset.getParameter<bool>("Apply_K_NLOEW_ZZQQB")),
  addProdAnomalousProbabilities(pset.getParameter<bool>("addProdAnomalousProbabilities")),
  reweight(),
  sampleName(pset.getParameter<string>("sampleName")),
  hTH2D_Mu_All(0),
  hTH2F_El_Reco(0),
  hTH2D_El_IdIsoSip_notCracks(0),
  hTH2D_El_IdIsoSip_Cracks(0),
  h_weight(0)
{
  //cout<< "Beginning Constructor\n\n\n" <<endl;

  Nevt_Gen_reweighted.resize(nReweightingSamples);
  Nevt_Gen_lumiBlock_reweighted.resize(nReweightingSamples);
  gen_ZZ4mu_reweighted.resize(nReweightingSamples);
  gen_ZZ4e_reweighted.resize(nReweightingSamples);
  gen_ZZ2mu2e_reweighted.resize(nReweightingSamples);
  gen_ZZ2l2tau_reweighted.resize(nReweightingSamples);
  gen_ZZ2emu2tau_reweighted.resize(nReweightingSamples);
  gen_ZZ4tau_reweighted.resize(nReweightingSamples);
  gen_ZZ4mu_EtaAcceptance_reweighted.resize(nReweightingSamples);
  gen_ZZ4mu_LeptonAcceptance_reweighted.resize(nReweightingSamples);
  gen_ZZ4e_EtaAcceptance_reweighted.resize(nReweightingSamples);
  gen_ZZ4e_LeptonAcceptance_reweighted.resize(nReweightingSamples);
  gen_ZZ2mu2e_EtaAcceptance_reweighted.resize(nReweightingSamples);
  gen_ZZ2mu2e_LeptonAcceptance_reweighted.resize(nReweightingSamples);
  gen_BUGGY_reweighted.resize(nReweightingSamples);
  gen_Unknown_reweighted.resize(nReweightingSamples);
  gen_sumGenMCWeight_reweighted.resize(nReweightingSamples);
  gen_sumWeights_reweighted.resize(nReweightingSamples);

  consumesMany<std::vector< PileupSummaryInfo > >();
  genParticleToken = consumes<edm::View<reco::Candidate> >(edm::InputTag("prunedGenParticles"));
  genInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  consumesMany<LHEEventProduct>();
  candToken = consumes<edm::View<pat::CompositeCandidate> >(edm::InputTag(theCandLabel));

  is_loose_ele_selection = false;
  if(pset.exists("is_loose_ele_selection")) { 
    is_loose_ele_selection = pset.getParameter<bool>("is_loose_ele_selection");
  }
  triggerResultToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults"));
  vtxToken = consumes<vector<reco::Vertex> >(edm::InputTag("goodPrimaryVertices"));
  jetToken = consumes<edm::View<pat::Jet> >(edm::InputTag("cleanJets"));
  metToken     = consumes<pat::METCollection>(edm::InputTag("slimmedMETs"));
  metNoHFToken = consumes<pat::METCollection>(edm::InputTag("slimmedMETsNoHF"));
  muonToken = consumes<pat::MuonCollection>(edm::InputTag("slimmedMuons"));
  electronToken = consumes<pat::ElectronCollection>(edm::InputTag("slimmedElectrons"));

  preSkimToken = consumes<edm::MergeableCounter,edm::InLumi>(edm::InputTag("preSkimCounter"));

  if (skipEmptyEvents) {
    applyTrigger=true;
    applySkim=true;
  } else {
    applyTrigger=false;
    applySkim=false;
  }

  isMC = myHelper.isMC();
  if (isMC) lheHandler = new LHEHandler(pset.getParameter<int>("VVMode"), pset.getParameter<int>("VVDecayMode"), addLHEKinematics || doreweighting);

  Nevt_Gen = 0;
  Nevt_Gen_lumiBlock = 0;

  //For Efficiency studies
  gen_ZZ4mu = 0;
  gen_ZZ4e = 0;
  gen_ZZ2mu2e = 0;
  gen_ZZ2l2tau = 0;
  gen_ZZ2emu2tau = 0;
  gen_ZZ4tau = 0;
  gen_ZZ4mu_EtaAcceptance = 0;
  gen_ZZ4mu_LeptonAcceptance = 0;
  gen_ZZ4e_EtaAcceptance = 0;
  gen_ZZ4e_LeptonAcceptance = 0;
  gen_ZZ2mu2e_EtaAcceptance = 0;
  gen_ZZ2mu2e_LeptonAcceptance = 0;
  gen_BUGGY = 0;
  gen_Unknown = 0;

  gen_sumPUWeight = 0.f;
  gen_sumGenMCWeight = 0.f;
  gen_sumWeights =0.f;

  std::string fipPath;

  // Read EWK K-factor table from file
  edm::FileInPath ewkFIP("ZZAnalysis/AnalysisStep/data/kfactors/ZZ_EwkCorrections.dat");
  fipPath=ewkFIP.fullPath();
  ewkTable = EwkCorrections::readFile_and_loadEwkTable(fipPath.data());

  // Read the ggZZ k-factor shape from file
  TString strZZGGKFVar[9]={
    "Nominal", "PDFScaleDn", "PDFScaleUp", "QCDScaleDn", "QCDScaleUp", "AsDn", "AsUp", "PDFReplicaDn", "PDFReplicaUp"
  };
  edm::FileInPath ggzzFIP_NNLO("ZZAnalysis/AnalysisStep/data/kfactors/Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
  fipPath=ggzzFIP_NNLO.fullPath();
  TFile* ggZZKFactorFile = TFile::Open(fipPath.data());
  for (unsigned int ikf=0; ikf<9; ikf++) spkfactor_ggzz_nnlo[ikf] = (TSpline3*)ggZZKFactorFile->Get(Form("sp_kfactor_%s", strZZGGKFVar[ikf].Data()))->Clone(Form("sp_kfactor_%s_NNLO", strZZGGKFVar[ikf].Data()));
  ggZZKFactorFile->Close();
  edm::FileInPath ggzzFIP_NLO("ZZAnalysis/AnalysisStep/data/kfactors/Kfactor_Collected_ggHZZ_2l2l_NLO_NNPDF_NarrowWidth_13TeV.root");
  fipPath=ggzzFIP_NLO.fullPath();
  ggZZKFactorFile = TFile::Open(fipPath.data());
  for (unsigned int ikf=0; ikf<9; ikf++) spkfactor_ggzz_nlo[ikf] = (TSpline3*)ggZZKFactorFile->Get(Form("sp_kfactor_%s", strZZGGKFVar[ikf].Data()))->Clone(Form("sp_kfactor_%s_NLO", strZZGGKFVar[ikf].Data()));
  ggZZKFactorFile->Close();

  //Scale factors for data/MC efficiency
  if (!skipMuDataMCWeight) {
    TString filename;
    filename.Form("ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ScaleFactors_mu_%d.root",year);
    edm::FileInPath fipMu(filename.Data());
    fipPath = fipMu.fullPath();
    TFile *fMuWeight = TFile::Open(fipPath.data(),"READ");
    hTH2D_Mu_All = (TH2D*)fMuWeight->Get("FINAL")->Clone();
    fMuWeight->Close();
  }
  
  if (!skipEleDataMCWeight) {

    if(year>=2016) {
        TString filename("ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ScaleFactors_ele_2016_v3.root");
        //TString filename("ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ele_scale_factors_2016_v1.root");  
        edm::FileInPath fipEleNotCracks(filename.Data());
        fipPath = fipEleNotCracks.fullPath();
        TFile *root_file = TFile::Open(fipPath.data(),"READ");
        hTH2D_El_IdIsoSip_notCracks = (TH1*) root_file->Get("ele_scale_factors")->Clone();
        hTH2D_El_IdIsoSip_Cracks = (TH1*) root_file->Get("ele_scale_factors_gap")->Clone();
        root_file->Close();

        TString filenameEleReco("ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/egammaEffi.txt_SF2D.root");
        edm::FileInPath fipEleReco(filenameEleReco.Data());
        fipPath = fipEleReco.fullPath();
        TFile *root_file_reco = TFile::Open(fipPath.data(),"READ");
        hTH2F_El_Reco = (TH2F*) root_file_reco->Get("EGamma_SF2D")->Clone();
        root_file_reco->Close();

    } else {
        TString filename;
        filename.Form("ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ScaleFactors_ele_%d_IdIsoSip.root",year);
        edm::FileInPath fipEleNotCracks(filename.Data());
        fipPath = fipEleNotCracks.fullPath();
        TFile *fEleWeightNotCracks = TFile::Open(fipPath.data(),"READ");
        hTH2D_El_IdIsoSip_notCracks = (TH2D*)fEleWeightNotCracks->Get("hScaleFactors_IdIsoSip")->Clone();
        fEleWeightNotCracks->Close();


        filename.Form("ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ScaleFactors_ele_%d_IdIsoSip_Cracks.root",year);
        edm::FileInPath fipEleCracks(filename.Data());
        fipPath = fipEleCracks.fullPath();
        TFile *fEleWeightCracks = TFile::Open(fipPath.data(),"READ");
        hTH2D_El_IdIsoSip_Cracks = (TH2D*)fEleWeightCracks->Get("hScaleFactors_IdIsoSip_Cracks")->Clone();
        fEleWeightCracks->Close();
    }
  }

  if (!skipHqTWeight) {
    //HqT weights
    edm::FileInPath HqTfip("ZZAnalysis/AnalysisStep/test/Macros/HqTWeights.root");
    std::string fipPath=HqTfip.fullPath();
    TFile *fHqt = TFile::Open(fipPath.data(),"READ");
    h_weight = (TH2D*)fHqt->Get("wH_400")->Clone();//FIXME: Ask simon to provide the 2D histo
    fHqt->Close();
  }

  if (!skipFakeWeight) {
    //CR fake rate weight
    TString filename;filename.Form("ZZAnalysis/AnalysisStep/test/Macros/FR2_2011_AA_electron.root");
    if(year==2015)filename.Form("ZZAnalysis/AnalysisStep/test/Macros/FR2_AA_ControlSample_ABCD.root");
    edm::FileInPath fipEleZX(filename.Data());
    std::string fipPath=fipEleZX.fullPath();
    TFile *FileZXWeightEle = TFile::Open(fipPath.data(),"READ");

    filename.Form("ZZAnalysis/AnalysisStep/test/Macros/FR2_2011_AA_muon.root");
    if(year==2015)filename.Form("ZZAnalysis/AnalysisStep/test/Macros/FR2_AA_muon.root");
    edm::FileInPath fipMuZX(filename.Data());
    fipPath=fipMuZX.fullPath();
    TFile *FileZXWeightMuo = TFile::Open(fipPath.data(),"READ");

    h_ZXWeight[0]=(TH2D*)FileZXWeightEle->Get("eff_Z1ee_plus_electron")->Clone();
    h_ZXWeight[1]=(TH2D*)FileZXWeightEle->Get("eff_Z1mumu_plus_electron")->Clone();
    h_ZXWeight[2]=(TH2D*)FileZXWeightMuo->Get("eff_Z1ee_plus_muon")->Clone();
    h_ZXWeight[3]=(TH2D*)FileZXWeightMuo->Get("eff_Z1mumu_plus_muon")->Clone();

    FileZXWeightEle->Close();
    FileZXWeightMuo->Close();
  }
}

HZZ4lNtupleMaker::~HZZ4lNtupleMaker()
{
  clearMELABranches(); // Cleans LHE branches
  if (lheHandler!=0) delete lheHandler;
}


// ------------ method called for each event  ------------
void HZZ4lNtupleMaker::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  myTree->InitializeVariables();

  //----------------------------------------------------------------------
  // Analyze MC truth; collect MC weights and update counters (this is done for all generated events,
  // including those that do not pass skim, trigger etc!)

  bool gen_ZZ4lInEtaAcceptance = false;   // All 4 gen leptons in eta acceptance
  bool gen_ZZ4lInEtaPtAcceptance = false; // All 4 gen leptons in eta,pT acceptance

  const reco::Candidate * genH = 0;
  std::vector<const reco::Candidate *> genZLeps;
  std::vector<const reco::Candidate *> genAssocLeps;

  edm::Handle<GenEventInfoProduct> genInfo;

  if (isMC) {
    // get PU weights
    vector<Handle<std::vector< PileupSummaryInfo > > >  PupInfos; //FIXME support for miniAOD v1/v2 where name changed; catch does not work...
    event.getManyByType(PupInfos);
    Handle<std::vector< PileupSummaryInfo > > PupInfo = PupInfos.front();
//     try {
//       cout << "TRY HZZ4lNtupleMaker" <<endl;
//       event.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
//     } catch (const cms::Exception& e){
//       cout << "FAIL HZZ4lNtupleMaker" <<endl;
//       event.getByLabel(edm::InputTag("slimmedAddPileupInfo"), PupInfo);
//     }

    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      if(PVI->getBunchCrossing() == 0) {
        NObsInt  = PVI->getPU_NumInteractions();
        NTrueInt = PVI->getTrueNumInteractions();
        break;
      }
    }

    // get PU weight
    PUWeight = reweight.weight(myHelper.sampleType(), myHelper.setup(), NTrueInt);

    event.getByToken(genParticleToken, genParticles);
    event.getByToken(genInfoToken, genInfo);

    MCHistoryTools mch(event, sampleName, genParticles, genInfo);
    genFinalState = mch.genFinalState();
    genProcessId = mch.getProcessID();
    genHEPMCweight = mch.gethepMCweight();
    genExtInfo = mch.genAssociatedFS();

    //Information on generated candidates, will be used later
    genH = mch.genH();
    genZLeps     = mch.sortedGenZZLeps();
    genAssocLeps = mch.genAssociatedLeps();
    genFSR       = mch.genFSR();

    if(genH != 0){
      FillHGenInfo(genH->p4(),getHqTWeight(genH->p4().M(),genH->p4().Pt()));
    }
    else if(genZLeps.size()==4){ // for 4l events take the mass of the ZZ(4l) system
      FillHGenInfo((genZLeps.at(0)->p4()+genZLeps.at(1)->p4()+genZLeps.at(2)->p4()+genZLeps.at(3)->p4()),0);
    }

    if (genFinalState!=BUGGY) {

      if (genZLeps.size()==4) {

        // "generated Zs" defined with standard pairing applied on gen leptons (genZLeps is sorted by MCHistoryTools)
        FillZGenInfo(genZLeps.at(0)->pdgId()*genZLeps.at(1)->pdgId(),
                     genZLeps.at(2)->pdgId()*genZLeps.at(3)->pdgId(),
                     genZLeps.at(0)->p4()+genZLeps.at(1)->p4(),
                     genZLeps.at(2)->p4()+genZLeps.at(3)->p4());

        // Gen leptons
        FillLepGenInfo(genZLeps.at(0)->pdgId(), genZLeps.at(1)->pdgId(), genZLeps.at(2)->pdgId(), genZLeps.at(3)->pdgId(),
           genZLeps.at(0)->p4(), genZLeps.at(1)->p4(), genZLeps.at(2)->p4(), genZLeps.at(3)->p4());

      }

      if (genZLeps.size()==3) {
        FillLepGenInfo(genZLeps.at(0)->pdgId(), genZLeps.at(1)->pdgId(), genZLeps.at(2)->pdgId(), 0,
                       genZLeps.at(0)->p4(), genZLeps.at(1)->p4(), genZLeps.at(2)->p4(), *(new math::XYZTLorentzVector));
      }
      if (genZLeps.size()==2) {
        FillLepGenInfo(genZLeps.at(0)->pdgId(), genZLeps.at(1)->pdgId(), 0, 0,
                       genZLeps.at(0)->p4(), genZLeps.at(1)->p4(), *(new math::XYZTLorentzVector), *(new math::XYZTLorentzVector));
      }

      if (genAssocLeps.size()==1 || genAssocLeps.size()==2) {
        FillAssocLepGenInfo(genAssocLeps);
      }

    }

    // LHE information
    if (isMC){
      edm::Handle<LHEEventProduct> lhe_evt;
      vector<edm::Handle<LHEEventProduct> > lhe_handles;
      event.getManyByType(lhe_handles);
      if (lhe_handles.size()>0){
        lhe_evt = lhe_handles.front();
        lheHandler->setHandle(&lhe_evt);
        lheHandler->extract();
        FillLHECandidate();
        if (doreweighting) {
          reweighting.fillreweightingweights(reweightingweights, lheHandler->getBestCandidate());
        }
        lheHandler->clear();
      }
      //else cerr << "lhe_handles.size()==0" << endl;
    }
    //

    // keep track of sum of weights
    gen_sumPUWeight += PUWeight;
    addweight(gen_sumGenMCWeight, gen_sumGenMCWeight_reweighted, genHEPMCweight);
    addweight(gen_sumWeights, gen_sumWeights_reweighted, PUWeight*genHEPMCweight);

    mch.genAcceptance(gen_ZZ4lInEtaAcceptance, gen_ZZ4lInEtaPtAcceptance);

    addweight(Nevt_Gen_lumiBlock, Nevt_Gen_lumiBlock_reweighted, 1);
    if (genFinalState == EEEE) {
      addweight(gen_ZZ4e, gen_ZZ4e_reweighted, 1);
      if (gen_ZZ4lInEtaAcceptance) addweight(gen_ZZ4e_EtaAcceptance, gen_ZZ4e_EtaAcceptance_reweighted, 1);
      if (gen_ZZ4lInEtaPtAcceptance) addweight(gen_ZZ4e_LeptonAcceptance, gen_ZZ4e_LeptonAcceptance_reweighted, 1);;
    } else if (genFinalState == MMMM) {
      addweight(gen_ZZ4mu, gen_ZZ4mu_reweighted, 1);
      if (gen_ZZ4lInEtaAcceptance) addweight(gen_ZZ4mu_EtaAcceptance, gen_ZZ4mu_EtaAcceptance_reweighted, 1);
      if (gen_ZZ4lInEtaPtAcceptance) addweight(gen_ZZ4mu_LeptonAcceptance, gen_ZZ4mu_LeptonAcceptance_reweighted, 1);;
    } else if (genFinalState == EEMM) {
      addweight(gen_ZZ2mu2e, gen_ZZ2mu2e_reweighted, 1);
      if (gen_ZZ4lInEtaAcceptance) addweight(gen_ZZ2mu2e_EtaAcceptance, gen_ZZ2mu2e_EtaAcceptance_reweighted, 1);
      if (gen_ZZ4lInEtaPtAcceptance) addweight(gen_ZZ2mu2e_LeptonAcceptance, gen_ZZ2mu2e_LeptonAcceptance_reweighted, 1);;
    } else if (genFinalState == llTT){
      addweight(gen_ZZ2emu2tau, gen_ZZ2emu2tau_reweighted, 1);
      addweight(gen_ZZ2l2tau, gen_ZZ2l2tau_reweighted, 1);
    } else if (genFinalState == TTTT){
      addweight(gen_ZZ4tau, gen_ZZ4tau_reweighted, 1);
      addweight(gen_ZZ2l2tau, gen_ZZ2l2tau_reweighted, 1);
    } else if (genFinalState == BUGGY){ // handle H->ddbar 2012 generator bug!!!
      addweight(gen_BUGGY, gen_BUGGY_reweighted, 1);
      return; // BUGGY events are skipped
    } else {
      addweight(gen_Unknown, gen_Unknown_reweighted, 1);
    }

// End of MC history analysis ------------------------------------------
  } else {
    ++Nevt_Gen_lumiBlock; // keep track of # events for data as well
  }



  // Get candidate collection
  edm::Handle<edm::View<pat::CompositeCandidate> > candHandle;
  event.getByToken(candToken, candHandle);
  if(candHandle.failedToGet()) {
    if(is_loose_ele_selection) return; // The collection can be missing in this case since we have a filter to skip the module when a regular candidate is present.
    else edm::LogError("") << "ZZ collection not found in non-loose electron flow. This should never happen";
  }
  const edm::View<pat::CompositeCandidate>* cands = candHandle.product();

  if (skipEmptyEvents && cands->size() == 0) return; // Skip events with no candidate, unless skipEmptyEvents = false

  // For Z+L CRs, we want only events with exactly 1 Z+l candidate. FIXME: this has to be reviewed.
  if (theChannel==ZL && cands->size() != 1) return;


  // Retrieve trigger results
  Handle<edm::TriggerResults> triggerResults;
  event.getByToken(triggerResultToken, triggerResults);

  // Apply MC filter (skip event)
  if (isMC && !(myHelper.passMCFilter(event,triggerResults))) return;

  // Apply skim
  bool evtPassSkim = myHelper.passSkim(event,triggerResults,trigWord);
  if (applySkim && !evtPassSkim) return;

  // Apply trigger request (skip event)
  bool evtPassTrigger = myHelper.passTrigger(event,triggerResults,trigWord);
  if (applyTrigger && !evtPassTrigger) return;


  //Fill MC truth information
  if (isMC) FillKFactors(genInfo, genZLeps);

  // General event information
  RunNumber=event.id().run();
  LumiNumber=event.luminosityBlock();
  EventNumber=event.id().event();
  xsection=xsec;


  // Primary vertices
  Handle<vector<reco::Vertex> > vertices;
  event.getByToken(vtxToken,vertices);
  Nvtx=vertices->size();


  // Jets (cleaned wrt all tight isolated leptons)
  Handle<edm::View<pat::Jet> > CleanedJets;
  event.getByToken(jetToken, CleanedJets);
  vector<const pat::Jet*> cleanedJets;
  for(edm::View<pat::Jet>::const_iterator jet = CleanedJets->begin(); jet != CleanedJets->end(); ++jet){
    cleanedJets.push_back(&*jet);
  }

  // MET
  Handle<pat::METCollection> metHandle;
  Handle<pat::METCollection> metNoHFHandle;
  event.getByToken(metToken, metHandle);
  event.getByToken(metNoHFToken, metNoHFHandle);
  if(metHandle.isValid()){
    PFMET = metHandle->front().pt();
    PFMETPhi = metHandle->front().phi();
  }
  if(metNoHFHandle.isValid()){
    PFMETNoHF = metNoHFHandle->front().pt();
    PFMETNoHFPhi = metNoHFHandle->front().phi();
  }


  // number of reconstructed leptons
  edm::Handle<pat::MuonCollection> muonHandle;
  event.getByToken(muonToken, muonHandle);
  for(unsigned int i = 0; i< muonHandle->size(); ++i){
    const pat::Muon* m = &((*muonHandle)[i]);
    if(m->pt()>5 && m->isPFMuon()) // these cuts are implicit in miniAOD
      NRecoMu++;
  }
  edm::Handle<pat::ElectronCollection> electronHandle;
  event.getByToken(electronToken, electronHandle);
  for(unsigned int i = 0; i< electronHandle->size(); ++i){
    const pat::Electron* e = &((*electronHandle)[i]);
    if(e->pt()>5) // this cut is implicit in miniAOD
      NRecoEle++;
  }


  //Loop on the candidates
  vector<Int_t> CRFLAG(cands->size());
  for( edm::View<pat::CompositeCandidate>::const_iterator cand = cands->begin(); cand != cands->end(); ++cand) {
    size_t icand= cand-cands->begin();

    //    int candChannel = cand->userFloat("candChannel"); // This is currently the product of pdgId of leptons (eg 14641, 28561, 20449)

    if (theChannel==ZLL) {
      // Cross check region for Z + 1 loose electron + 1 loose TLE (defined only in loose_ele paths) 
      if (is_loose_ele_selection) {
        if (cand->userFloat("isBestCRZLL")&&cand->userFloat("CRZLL")) set_bit(CRFLAG[icand],ZLL);
      }

      // AA CRs
      if (cand->userFloat("isBestCRZLLss")&&cand->userFloat("CRZLLss")) set_bit(CRFLAG[icand],CRZLLss);

      // A CRs
      if (cand->userFloat("isBestCRZLLos_2P2F")&&cand->userFloat("CRZLLos_2P2F")) set_bit(CRFLAG[icand],CRZLLos_2P2F);
      if (cand->userFloat("isBestCRZLLos_3P1F")&&cand->userFloat("CRZLLos_3P1F")) set_bit(CRFLAG[icand],CRZLLos_3P1F);

      if (CRFLAG[icand]) { // This candidate belongs to one of the CRs: perform additional jet cleaning.
        // Note that this is (somewhat incorrectly) done per-event, so there could be some over-cleaning in events with >1 CR candidate.
        for (unsigned i=0; i<cleanedJets.size(); ++i) {
          if (cleanedJets[i]!=0  && (!jetCleaner::isGood(*cand, *(cleanedJets[i])))) {
            cleanedJets[i]=0;
          }
        }
      }
    }
  }

  // Count and store jets, after additional cleaning for CRs...
  for (unsigned i=0; i<cleanedJets.size(); ++i) {
    if (cleanedJets[i]==0) {
      continue; // Jet has been suppressed by additional cleaning
    }

    ++nCleanedJets;
    if(cleanedJets[i]->pt()>30){
      ++nCleanedJetsPt30;
      if(cleanedJets[i]->userFloat("isBtagged")) ++nCleanedJetsPt30BTagged;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF_Up")) ++nCleanedJetsPt30BTagged_bTagSFUp;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF_Dn")) ++nCleanedJetsPt30BTagged_bTagSFDn;
    }

    // count jec up/down njets pt30
    float jec_unc = cleanedJets[i]->userFloat("jec_unc");

    float pt_up = cleanedJets[i]->pt() * (1.0 + jec_unc);
    float pt_dn = cleanedJets[i]->pt() * (1.0 - jec_unc);

    if (pt_up>30) ++nCleanedJetsPt30_jecUp;
    if (pt_dn>30) ++nCleanedJetsPt30_jecDn;


    if (writeJets && theChannel!=ZL) FillJet(*(cleanedJets.at(i))); // No additional pT cut (for JEC studies)
  }

  // Now we can write the variables for candidates
  int nFilled=0;
  for( edm::View<pat::CompositeCandidate>::const_iterator cand = cands->begin(); cand != cands->end(); ++cand) {
    size_t icand= cand-cands->begin();

    if (!( theChannel==ZL || CRFLAG[icand] || (bool)(cand->userFloat("isBestCand")) )) continue; // Skip events other than the best cand (or CR candidates in the CR)

    //For the SR, also fold information about acceptance in CRflag.
    if (isMC && (theChannel==ZZ)) {
      if (gen_ZZ4lInEtaAcceptance)   set_bit(CRFLAG[icand],28);
      if (gen_ZZ4lInEtaPtAcceptance) set_bit(CRFLAG[icand],29);
    }
    FillCandidate(*cand, evtPassTrigger&&evtPassSkim, event, CRFLAG[icand]);

    // Fill the candidate as one entry in the tree. Do not reinitialize the event variables, as in CRs
    // there could be several candidates per event.
    myTree->FillCurrentTree();
    ++nFilled;
  }

  // If no candidate was filled but we still want to keep gen-level and weights, we need to fill one entry anyhow.
  if (skipEmptyEvents==false && nFilled==0) myTree->FillCurrentTree();
}


void HZZ4lNtupleMaker::FillJet(const pat::Jet& jet)
{
   JetPt  .push_back( jet.pt());
   JetEta .push_back( jet.eta());
   JetPhi .push_back( jet.phi());
   JetMass .push_back( jet.p4().M());
   JetBTagger .push_back( jet.userFloat("bTagger"));
   JetIsBtagged .push_back( jet.userFloat("isBtagged"));
   JetIsBtaggedWithSF .push_back( jet.userFloat("isBtaggedWithSF"));
   JetIsBtaggedWithSFUp .push_back( jet.userFloat("isBtaggedWithSF_Up"));
   JetIsBtaggedWithSFDn .push_back( jet.userFloat("isBtaggedWithSF_Dn"));
   JetQGLikelihood .push_back( jet.userFloat("qgLikelihood"));
   if(addQGLInputs){
     JetAxis2 .push_back( jet.userFloat("axis2"));
     JetMult .push_back( jet.userFloat("mult"));
     JetPtD .push_back( jet.userFloat("ptD"));
   }
   JetSigma .push_back(jet.userFloat("jec_unc"));
   JetHadronFlavour .push_back(jet.hadronFlavour());
}


void HZZ4lNtupleMaker::FillKFactors(edm::Handle<GenEventInfoProduct>& genInfo, std::vector<const reco::Candidate *>& genZLeps){
  KFactor_QCD_ggZZ_Nominal=1;
  KFactor_QCD_ggZZ_PDFScaleDn=1;
  KFactor_QCD_ggZZ_PDFScaleUp=1;
  KFactor_QCD_ggZZ_QCDScaleDn=1;
  KFactor_QCD_ggZZ_QCDScaleUp=1;
  KFactor_QCD_ggZZ_AsDn=1;
  KFactor_QCD_ggZZ_AsUp=1;
  KFactor_QCD_ggZZ_PDFReplicaDn=1;
  KFactor_QCD_ggZZ_PDFReplicaUp=1;
  KFactor_QCD_qqZZ_dPhi=1;
  KFactor_QCD_qqZZ_M=1;
  KFactor_QCD_qqZZ_Pt=1;
  KFactor_EW_qqZZ=1;
  KFactor_EW_qqZZ_unc=0;

  if (isMC){
    GenEventInfoProduct  genInfoP = *(genInfo.product());
    if (apply_K_NNLOQCD_ZZGG>0 && apply_K_NNLOQCD_ZZGG!=3){
      if (spkfactor_ggzz_nnlo[0]!=0) KFactor_QCD_ggZZ_Nominal = (float)spkfactor_ggzz_nnlo[0]->Eval(GenHMass);
      if (spkfactor_ggzz_nnlo[1]!=0) KFactor_QCD_ggZZ_PDFScaleDn = (float)spkfactor_ggzz_nnlo[1]->Eval(GenHMass);
      if (spkfactor_ggzz_nnlo[2]!=0) KFactor_QCD_ggZZ_PDFScaleUp = (float)spkfactor_ggzz_nnlo[2]->Eval(GenHMass);
      if (spkfactor_ggzz_nnlo[3]!=0) KFactor_QCD_ggZZ_QCDScaleDn = (float)spkfactor_ggzz_nnlo[3]->Eval(GenHMass);
      if (spkfactor_ggzz_nnlo[4]!=0) KFactor_QCD_ggZZ_QCDScaleUp = (float)spkfactor_ggzz_nnlo[4]->Eval(GenHMass);
      if (spkfactor_ggzz_nnlo[5]!=0) KFactor_QCD_ggZZ_AsDn = (float)spkfactor_ggzz_nnlo[5]->Eval(GenHMass);
      if (spkfactor_ggzz_nnlo[6]!=0) KFactor_QCD_ggZZ_AsUp = (float)spkfactor_ggzz_nnlo[6]->Eval(GenHMass);
      if (spkfactor_ggzz_nnlo[7]!=0) KFactor_QCD_ggZZ_PDFReplicaDn = (float)spkfactor_ggzz_nnlo[7]->Eval(GenHMass);
      if (spkfactor_ggzz_nnlo[8]!=0) KFactor_QCD_ggZZ_PDFReplicaUp = (float)spkfactor_ggzz_nnlo[8]->Eval(GenHMass);
      if (apply_K_NNLOQCD_ZZGG==2){
        if (spkfactor_ggzz_nlo[0]!=0){
          float divisor = (float)spkfactor_ggzz_nlo[0]->Eval(GenHMass);
          KFactor_QCD_ggZZ_Nominal /= divisor;
          KFactor_QCD_ggZZ_PDFScaleDn /= divisor;
          KFactor_QCD_ggZZ_PDFScaleUp /= divisor;
          KFactor_QCD_ggZZ_QCDScaleDn /= divisor;
          KFactor_QCD_ggZZ_QCDScaleUp /= divisor;
          KFactor_QCD_ggZZ_AsDn /= divisor;
          KFactor_QCD_ggZZ_AsUp /= divisor;
          KFactor_QCD_ggZZ_PDFReplicaDn /= divisor;
          KFactor_QCD_ggZZ_PDFReplicaUp /= divisor;
        }
        else{
          KFactor_QCD_ggZZ_Nominal=0;
          KFactor_QCD_ggZZ_PDFScaleDn=0;
          KFactor_QCD_ggZZ_PDFScaleUp=0;
          KFactor_QCD_ggZZ_QCDScaleDn=0;
          KFactor_QCD_ggZZ_QCDScaleUp=0;
          KFactor_QCD_ggZZ_AsDn=0;
          KFactor_QCD_ggZZ_AsUp=0;
          KFactor_QCD_ggZZ_PDFReplicaDn=0;
          KFactor_QCD_ggZZ_PDFReplicaUp=0;
        }
      }
    }
    else if (apply_K_NNLOQCD_ZZGG==3){
      if (spkfactor_ggzz_nlo[0]!=0) KFactor_QCD_ggZZ_Nominal = (float)spkfactor_ggzz_nlo[0]->Eval(GenHMass);
      if (spkfactor_ggzz_nlo[1]!=0) KFactor_QCD_ggZZ_PDFScaleDn = (float)spkfactor_ggzz_nlo[1]->Eval(GenHMass);
      if (spkfactor_ggzz_nlo[2]!=0) KFactor_QCD_ggZZ_PDFScaleUp = (float)spkfactor_ggzz_nlo[2]->Eval(GenHMass);
      if (spkfactor_ggzz_nlo[3]!=0) KFactor_QCD_ggZZ_QCDScaleDn = (float)spkfactor_ggzz_nlo[3]->Eval(GenHMass);
      if (spkfactor_ggzz_nlo[4]!=0) KFactor_QCD_ggZZ_QCDScaleUp = (float)spkfactor_ggzz_nlo[4]->Eval(GenHMass);
      if (spkfactor_ggzz_nlo[5]!=0) KFactor_QCD_ggZZ_AsDn = (float)spkfactor_ggzz_nlo[5]->Eval(GenHMass);
      if (spkfactor_ggzz_nlo[6]!=0) KFactor_QCD_ggZZ_AsUp = (float)spkfactor_ggzz_nlo[6]->Eval(GenHMass);
      if (spkfactor_ggzz_nlo[7]!=0) KFactor_QCD_ggZZ_PDFReplicaDn = (float)spkfactor_ggzz_nlo[7]->Eval(GenHMass);
      if (spkfactor_ggzz_nlo[8]!=0) KFactor_QCD_ggZZ_PDFReplicaUp = (float)spkfactor_ggzz_nlo[8]->Eval(GenHMass);
    }

    if (genFinalState!=BUGGY){
      if (genZLeps.size()==4) {
        // Calculate NNLO/NLO QCD K factors for qqZZ
        if (apply_K_NNLOQCD_ZZQQB){
          bool sameflavor=(genZLeps.at(0)->pdgId()*genZLeps.at(1)->pdgId() == genZLeps.at(2)->pdgId()*genZLeps.at(3)->pdgId());
          KFactor_QCD_qqZZ_dPhi = kfactor_qqZZ_qcd_dPhi(fabs(GenZ1Phi-GenZ2Phi), (sameflavor) ? 1 : 2);
          KFactor_QCD_qqZZ_M    = kfactor_qqZZ_qcd_M(GenHMass, (sameflavor) ? 1 : 2, 2) / kfactor_qqZZ_qcd_M(GenHMass, (sameflavor) ? 1 : 2, 1);
          KFactor_QCD_qqZZ_Pt   = kfactor_qqZZ_qcd_Pt(GenHPt, (sameflavor) ? 1 : 2);
        }
        // Calculate NLO EWK K factors for qqZZ
        if (apply_K_NLOEW_ZZQQB){
          TLorentzVector GENZ1Vec, GENZ2Vec, GENZZVec;
          GENZ1Vec.SetPtEtaPhiM(GenZ1Pt, GenZ1Eta, GenZ1Phi, GenZ1Mass);
          GENZ2Vec.SetPtEtaPhiM(GenZ2Pt, GenZ2Eta, GenZ2Phi, GenZ2Mass);
          GENZZVec = GENZ1Vec + GENZ2Vec;

          KFactor_EW_qqZZ = EwkCorrections::getEwkCorrections(genParticles, ewkTable, genInfoP, GENZ1Vec, GENZ2Vec);

          bool sameflavor=(genZLeps.at(0)->pdgId()*genZLeps.at(1)->pdgId() == genZLeps.at(2)->pdgId()*genZLeps.at(3)->pdgId());
          float K_NNLO_LO = kfactor_qqZZ_qcd_M(GenHMass, (sameflavor) ? 1 : 2, 2);
          float rho = GENZZVec.Pt()/(GenLep1Pt+GenLep2Pt+GenLep3Pt+GenLep4Pt);
          if (rho<0.3) KFactor_EW_qqZZ_unc = fabs((K_NNLO_LO-1.)*(1.-KFactor_EW_qqZZ));
          else KFactor_EW_qqZZ_unc = fabs(1.-KFactor_EW_qqZZ);
        }
      }
    }
  }

}


void HZZ4lNtupleMaker::FillLHECandidate(){
  LHEMotherPz.clear(); //FIXME
  LHEMotherE.clear();
  LHEMotherId.clear();
  LHEDaughterPt.clear();
  LHEDaughterEta.clear();
  LHEDaughterPhi.clear();
  LHEDaughterMass.clear();
  LHEDaughterId.clear();
  LHEAssociatedParticlePt.clear();
  LHEAssociatedParticleEta.clear();
  LHEAssociatedParticlePhi.clear();
  LHEAssociatedParticleMass.clear();
  LHEAssociatedParticleId.clear();

  LHEPDFScale = 0;
  LHEweight_QCDscale_muR1_muF1=0;
  LHEweight_QCDscale_muR1_muF2=0;
  LHEweight_QCDscale_muR1_muF0p5=0;
  LHEweight_QCDscale_muR2_muF1=0;
  LHEweight_QCDscale_muR2_muF2=0;
  LHEweight_QCDscale_muR2_muF0p5=0;
  LHEweight_QCDscale_muR0p5_muF1=0;
  LHEweight_QCDscale_muR0p5_muF2=0;
  LHEweight_QCDscale_muR0p5_muF0p5=0;

  MELACandidate* cand = lheHandler->getBestCandidate();
  if (cand!=0 && addLHEKinematics){
    for (int imot=0; imot<cand->getNMothers(); imot++){
      MELAParticle* apart = cand->getMother(imot);
      if (apart==0){ LHEMotherPz.clear(); LHEMotherE.clear(); LHEMotherId.clear(); break; } // Something went wrong
      LHEMotherPz.push_back(apart->z());
      LHEMotherE.push_back(apart->t());
      LHEMotherId.push_back((short)apart->id);
    }

    for (int iV=0; iV<min(2, cand->getNSortedVs()); iV++){
      MELAParticle* Vi = cand->getSortedV(iV);
      if (Vi!=0){
        for (int iVj=0; iVj<Vi->getNDaughters(); iVj++){
          MELAParticle* Vij = Vi->getDaughter(iVj);
          if (Vij!=0){
            LHEDaughterPt.push_back(Vij->pt());
            LHEDaughterEta.push_back(Vij->eta());
            LHEDaughterPhi.push_back(Vij->phi());
            LHEDaughterMass.push_back(Vij->m());
            LHEDaughterId.push_back((short)Vij->id);
          }
        }
      }
    }

    vector<MELAParticle*> AssociatedParticle;
    vector<MELAParticle*> tmpAssociatedParticle;
    for (int aa=0; aa<cand->getNAssociatedJets(); aa++){
      MELAParticle* apart = cand->getAssociatedJet(aa);
      tmpAssociatedParticle.push_back(apart);
    }
    for (int aa=0; aa<cand->getNAssociatedLeptons(); aa++){
      MELAParticle* apart = cand->getAssociatedLepton(aa);
      if (!PDGHelpers::isANeutrino(apart->id)) tmpAssociatedParticle.push_back(apart);
    }
    for (int aa=0; aa<cand->getNAssociatedNeutrinos(); aa++){
      MELAParticle* apart = cand->getAssociatedNeutrino(aa);
      tmpAssociatedParticle.push_back(apart);
    }
    while (tmpAssociatedParticle.size()>0){ // Re-sort all associated particles by leading pT (categories are individually sorted, but mixing categories loses this sorting)
      MELAParticle* tmpPart=0;
      int pos=0;
      for (unsigned int el=0; el<tmpAssociatedParticle.size(); el++){
        if (tmpPart==0){ tmpPart = tmpAssociatedParticle.at(el); pos=el; }
        else if (tmpPart->pt()<tmpAssociatedParticle.at(el)->pt()){ tmpPart = tmpAssociatedParticle.at(el); pos=el; } // Safer to do in two steps
      }
      AssociatedParticle.push_back(tmpPart);
      tmpAssociatedParticle.erase(tmpAssociatedParticle.begin()+pos);
    }
    for (unsigned int aa=0; aa<AssociatedParticle.size(); aa++){
      MELAParticle* apart = AssociatedParticle.at(aa);
      if (apart!=0){
        LHEAssociatedParticlePt.push_back(apart->pt());
        LHEAssociatedParticleEta.push_back(apart->eta());
        LHEAssociatedParticlePhi.push_back(apart->phi());
        LHEAssociatedParticleMass.push_back(apart->m());
        LHEAssociatedParticleId.push_back((short)apart->id);
      }
    }

    /*
    cout << "NEW EVENT:" << endl;
    cout << "Mothers:" << endl;
    for (unsigned int ipart=0; ipart<LHEMotherId.size(); ipart++) cout << "\t Mot" << ipart << " (pz, E, id) = " << LHEMotherPz.at(ipart) << " " << LHEMotherE.at(ipart) << " " << LHEMotherId.at(ipart) << endl;
    cout << "Daughters:" << endl;
    for (unsigned int ipart=0; ipart<LHEDaughterId.size(); ipart++) cout << "\t Dau" << ipart << " (pt, eta, phi, m, id) = " << LHEDaughterPt.at(ipart) << " " << LHEDaughterEta.at(ipart) << " " << LHEDaughterPhi.at(ipart) << " " << LHEDaughterMass.at(ipart) << " " << LHEDaughterId.at(ipart) << endl;
    cout << "Associated:" << endl;
    for (unsigned int ipart=0; ipart<LHEAssociatedParticleId.size(); ipart++) cout << "\t APart" << ipart << " (pt, eta, phi, m, id) = " << LHEAssociatedParticlePt.at(ipart) << " " << LHEAssociatedParticleEta.at(ipart) << " " << LHEAssociatedParticlePhi.at(ipart) << " " << LHEAssociatedParticleMass.at(ipart) << " " << LHEAssociatedParticleId.at(ipart) << endl;
    cout << endl;
    */

    computeMELABranches(cand);
    pushLHEMELABranches();
  }

  LHEPDFScale = lheHandler->getPDFScale();
  LHEweight_QCDscale_muR1_muF1 = lheHandler->getLHEWeight(0, 1.);
  LHEweight_QCDscale_muR1_muF2 = lheHandler->getLHEWeight(1, 1.);
  LHEweight_QCDscale_muR1_muF0p5 = lheHandler->getLHEWeight(2, 1.);
  LHEweight_QCDscale_muR2_muF1 = lheHandler->getLHEWeight(3, 1.);
  LHEweight_QCDscale_muR2_muF2 = lheHandler->getLHEWeight(4, 1.);
  LHEweight_QCDscale_muR2_muF0p5 = lheHandler->getLHEWeight(5, 1.);
  LHEweight_QCDscale_muR0p5_muF1 = lheHandler->getLHEWeight(6, 1.);
  LHEweight_QCDscale_muR0p5_muF2 = lheHandler->getLHEWeight(7, 1.);
  LHEweight_QCDscale_muR0p5_muF0p5 = lheHandler->getLHEWeight(8, 1.);
  LHEweight_PDFVariation_Up = lheHandler->getLHEWeight_PDFVariationUpDn(1, 1.);
  LHEweight_PDFVariation_Dn = lheHandler->getLHEWeight_PDFVariationUpDn(-1, 1.);
  LHEweight_AsMZ_Up = lheHandler->getLHEWeigh_AsMZUpDn(1, 1.);
  LHEweight_AsMZ_Dn = lheHandler->getLHEWeigh_AsMZUpDn(-1, 1.);
}


void HZZ4lNtupleMaker::FillCandidate(const pat::CompositeCandidate& cand, bool evtPass, const edm::Event& event, Int_t CRFLAG)
{
  //Initialize a new candidate into the tree
  //myTree->createNewCandidate(); // this doesn't do anything anymore

  //Reinitialize the per-candidate vectors (necessary because in CRs we can store more than 1 candidate per event)
  LepPt.clear();
  LepEta.clear();
  LepPhi.clear();
  LepLepId.clear();
  LepSIP.clear();
  LepTime.clear();
  LepisID.clear();
  LepBDT.clear();
  LepMissingHit.clear();
  //LepChargedHadIso.clear();
  //LepNeutralHadIso.clear();
  //LepPhotonIso.clear();
  LepCombRelIsoPF.clear();
  TLE_dR_Z = -1;
  fsrPt.clear();
  fsrEta.clear();
  fsrPhi.clear();
  fsrLept.clear();
  fsrLeptID.clear();
  fsrDR.clear();
  fsrGenPt.clear();
  ExtraLepPt.clear();
  ExtraLepEta.clear();
  ExtraLepPhi.clear();
  ExtraLepLepId.clear();

  CRflag = CRFLAG;

  if(theChannel!=ZL){
    //Fill the info on the Higgs candidate
    ZZMass = cand.p4().mass();
    ZZMassErr = cand.userFloat("massError");
    ZZMassErrCorr = cand.userFloat("massErrorCorr");
    ZZMassPreFSR = cand.userFloat("m4l");

    ZZPt  = cand.p4().pt();
    ZZEta = cand.p4().eta();
    ZZPhi = cand.p4().phi();

    if(addKinRefit){
      if (cand.hasUserFloat("ZZMassRefit")) {
  ZZMassRefit = cand.userFloat("ZZMassRefit");
  ZZMassRefitErr = cand.userFloat("ZZMassRefitErr");
  ZZMassUnrefitErr = cand.userFloat("ZZMassUnrefitErr");
      }
    }
    if(addVtxFit){
      ZZMassCFit = cand.userFloat("CFitM");
      ZZChi2CFit = cand.userFloat("CFitChi2");
    }

    DiJetMass  = cand.userFloat("DiJetMass");
    DiJetDEta  = cand.userFloat("DiJetDEta");
    DiJetFisher  = cand.userFloat("DiJetFisher");
    //    DiJetMassPlus  = cand.userFloat("DiJetMassPlus");
    //    DiJetMassMinus  = cand.userFloat("DiJetMassMinus");

    //Fill the angular variables
    helcosthetaZ1 = cand.userFloat("costheta1");
    helcosthetaZ2 = cand.userFloat("costheta2");
    helphi       = cand.userFloat("phi");
    costhetastar = cand.userFloat("costhetastar");
    phistarZ1      = cand.userFloat("phistar1");
    //phistarZ2      = cand.userFloat("phistar2");
    xi            = cand.userFloat("xi");
    xistar        = cand.userFloat("xistar");

    // Get MELA probabilities
    pushRecoMELABranches(cand);

  }

  //Z1 and Z2 variables
  const reco::Candidate* Z1;
  const reco::Candidate* Z2;
  vector<const reco::Candidate*> leptons;
  vector<const reco::Candidate*> fsrPhot;
  vector<short> fsrIndex;
  vector<string> labels;

  if (theChannel!=ZL) { // Regular 4l candidates
    Z1   = cand.daughter("Z1");
    Z2   = cand.daughter("Z2");
    userdatahelpers::getSortedLeptons(cand, leptons, labels, fsrPhot, fsrIndex);
  } else {              // Special handling of Z+l candidates
    Z1   = cand.daughter(0); // the Z
    Z2   = cand.daughter(1); // This is actually the additional lepton!
    userdatahelpers::getSortedLeptons(cand, leptons, labels, fsrPhot, fsrIndex, false); // note: we get just 3 leptons in this case.
  }

  Z1Mass = Z1->mass();
  Z1Pt =   Z1->pt();
  Z1Flav =  getPdgId(Z1->daughter(0)) * getPdgId(Z1->daughter(1));
 
  Z2Mass = Z2->mass();
  Z2Pt =   Z2->pt();
  Z2Flav = theChannel==ZL ? getPdgId(Z2) : getPdgId(Z2->daughter(0)) * getPdgId(Z2->daughter(1));

  const reco::Candidate* non_TLE_Z = nullptr;
  size_t TLE_index = 999;
  if(abs(Z1Flav) == 11*11 || abs(Z1Flav) == 13*13) non_TLE_Z = Z1;
  if(abs(Z2Flav) == 11*11 || abs(Z2Flav) == 13*13) non_TLE_Z = Z2;
  for (size_t i=0; i<leptons.size(); ++i){
    if(abs(leptons[i]->pdgId()) == 22) TLE_index = i;
  }
  if(TLE_index < 999 && non_TLE_Z != nullptr) {
    TLE_dR_Z = reco::deltaR(non_TLE_Z->p4(), leptons[TLE_index]->p4()); 
  }

  Int_t sel = 0;
  if(theChannel==ZZ){

    // Precomputed selections
    bool candPass70Z2Loose = cand.userFloat("Z2Mass") &&
                             cand.userFloat("MAllComb") &&
                             cand.userFloat("pt1")>20 && cand.userFloat("pt2")>10. &&
                             ZZMass>70.;
    bool candPassFullSel70 = cand.userFloat("SR");
    bool candPassFullSel   = cand.userFloat("FullSel");
    bool candIsBest = cand.userFloat("isBestCand");
    bool passMz_zz = (Z1Mass>60. && Z1Mass<120. && Z2Mass>60. && Z2Mass<120.);   //FIXME hardcoded cut

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

  }

  // Retrieve the userFloat of the leptons in vectors ordered in the same way.
  vector<float> SIP(4);
  vector<float> combRelIsoPF(4);
  passIsoPreFSR = true;
  for (unsigned int i=0; i<leptons.size(); ++i){
    float curr_dR = 999;
    if(i != TLE_index && TLE_index < 999)
      curr_dR = reco::deltaR(leptons[i]->p4(), leptons[TLE_index]->p4());
    if(curr_dR < TLE_min_dR_3l) TLE_min_dR_3l = curr_dR;

    short lepFlav = std::abs(leptons[i]->pdgId());

    SIP[i]             = userdatahelpers::getUserFloat(leptons[i],"SIP");
    passIsoPreFSR      = passIsoPreFSR&&(userdatahelpers::getUserFloat(leptons[i],"combRelIsoPF")<LeptonIsoHelper::isoCut(leptons[i]));

    //in the Legacy approach,  FSR-corrected iso is attached to the Z, not to the lepton!
    if (theChannel!=ZL) {
      combRelIsoPF[i]    = cand.userFloat(labels[i]+"combRelIsoPFFSRCorr"); // Note: the
      assert(SIP[i] == cand.userFloat(labels[i]+"SIP")); // Check that I don't mess up with labels[] and leptons[]
    } else {
      //FIXME: at the moment,  FSR-corrected iso is not computed for Z+L CRs
      combRelIsoPF[i]    = userdatahelpers::getUserFloat(leptons[i],"combRelIsoPF");
    }

    //Fill the info on the lepton candidates
    LepPt .push_back( leptons[i]->pt() );
    LepEta.push_back( leptons[i]->eta() );
    LepPhi.push_back( leptons[i]->phi() );
    int id =  leptons[i]->pdgId();
    if(id == 22 && (i == 1 || i == 3)) id=-22; //FIXME this assumes a standard ordering of leptons.
    LepLepId.push_back( id );
    LepSIP  .push_back( SIP[i] );
    LepTime .push_back( lepFlav==13 ? userdatahelpers::getUserFloat(leptons[i],"time") : 0. );
    LepisID .push_back( userdatahelpers::getUserFloat(leptons[i],"ID") );
    LepBDT  .push_back( lepFlav==11 ||lepFlav==22 ? userdatahelpers::getUserFloat(leptons[i],"BDT") : 0. );
    LepMissingHit.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"missingHit") : 0 );
    //LepChargedHadIso[i].push_back( userdatahelpers::getUserFloat(leptons[i],"PFChargedHadIso") );
    //LepNeutralHadIso[i].push_back( userdatahelpers::getUserFloat(leptons[i],"PFNeutralHadIso") );
    //LepPhotonIso[i].push_back( userdatahelpers::getUserFloat(leptons[i],"PFPhotonIso") );
    LepCombRelIsoPF.push_back( combRelIsoPF[i] );
    LepisLoose.push_back(userdatahelpers::hasUserFloat(leptons[i],"isLoose") == 1 ? userdatahelpers::getUserFloat(leptons[i],"isLoose") : -2);

  }

  // FSR
  for (unsigned i=0; i<fsrPhot.size(); ++i) {
    math::XYZTLorentzVector fsr = fsrPhot[i]->p4();
    fsrPt.push_back(fsr.pt());
    fsrEta.push_back(fsr.eta());
    fsrPhi.push_back(fsr.phi());
    fsrLept.push_back(fsrIndex[i]+1);
    fsrLeptID.push_back(leptons[fsrIndex[i]]->pdgId());
    fsrDR.push_back(ROOT::Math::VectorUtil::DeltaR(leptons[fsrIndex[i]]->momentum(), fsrPhot[i]->momentum()));
    int igen = MCHistoryTools::fsrMatch(fsrPhot[i], genFSR);
    double dRGenVsReco = -1.;
    double genpT = -1.;

    if (igen>=0) {
      dRGenVsReco = ROOT::Math::VectorUtil::DeltaR(genFSR[igen]->momentum(), fsrPhot[i]->momentum());
//       pTGen = genFSR[igen]->pt();
//       etaGen = genFSR[igen]->eta();
//       phiGen = genFSR[igen]->phi();
      if (dRGenVsReco<0.3) {// matching cut - FIXME
        genpT=genFSR[igen]->pt();
      }
    }
    fsrGenPt.push_back(genpT);


  }

  //convention: 0 -> 4mu   1 -> 4e   2 -> 2mu2e
  if(CRFLAG){
    ZXFakeweight = 1.;
    for(int izx=0;izx<2;izx++)
      ZXFakeweight *= getFakeWeight(Z2->daughter(izx)->pt(),Z2->daughter(izx)->eta(),Z2->daughter(izx)->pdgId(),Z1->daughter(0)->pdgId());
  }
  ZZsel = sel;

  if (theChannel!=ZL) {

    //Fill the info on categorization
    nExtraLep = cand.userFloat("nExtraLep"); // Why is this still a float at this point?
    nExtraZ=cand.userFloat("nExtraZ");

    //Fill the info on the extra leptons
    for (int iExtraLep=1; iExtraLep<=(int)nExtraLep; iExtraLep++){
      TString extraString;extraString.Form("ExtraLep%d",iExtraLep);
      if (cand.hasUserCand(extraString.Data())){
        //for(int iextra=0;iextra<4;iextra++)varExtra[iextra].Prepend(extraString.Data());
        reco::CandidatePtr candPtr=cand.userCand(extraString.Data());
        ExtraLepPt.push_back(candPtr->pt());
        ExtraLepEta.push_back(candPtr->eta());
        ExtraLepPhi.push_back(candPtr->phi());
        ExtraLepLepId.push_back(candPtr->pdgId());
      }
    }

  }

  //Compute the data/MC weight
  dataMCWeight = 1.;
  //When the trigger is not applied in the MC, apply a trigger efficiency factor instead (FIXME: here hardcoding the efficiencies computed for ICHEP2016)
  trigEffWeight = 1.;
  if(isMC) {
    for(unsigned int i=0; i<leptons.size(); ++i){
      dataMCWeight *= getAllWeight(leptons[i]);
    }
    if (applyTrigEffWeight){
      Int_t ZZFlav = abs(Z1Flav*Z2Flav);
      if(ZZFlav==121*121 || ZZFlav==121*242) //4e
	trigEffWeight = 0.992;
      if(ZZFlav==169*169) //4mu
	trigEffWeight = 0.996;
      if(ZZFlav==169*121 || ZZFlav==169*242) //2e2mu
	trigEffWeight = 0.995;
    }
  }
  //Store an overall event weight (which is normalized by gen_sumWeights)
  overallEventWeight = PUWeight * genHEPMCweight * dataMCWeight * trigEffWeight;

  /* // debug printout for weights
  cout<<"Event "<<event.id().run()<<":"<<event.luminosityBlock()<<":"<<event.id().event()<<endl;
  cout<<" pileup weight =         "<<PUWeight<<endl;
  cout<<" sign of gen. weight =   "<<genHEPMCweight/fabs(genHEPMCweight)<<endl;
  cout<<" lepton data/MC weight = "<<dataMCWeight<<endl;
  for(unsigned int i=0; i<leptons.size(); ++i)
    cout<<"   lepton ID="<<leptons[i]->pdgId()<<", pT="<<leptons[i]->pt()<<", weight="<<getAllWeight(leptons[i])<<endl;
  cout<<" trigger eff. weight =   "<<trigEffWeight<<endl;
  cout<<"product of all =         "<<overallEventWeight/fabs(genHEPMCweight)<<endl;;
  cout<<endl;
  //*/

}


void HZZ4lNtupleMaker::getCheckedUserFloat(const pat::CompositeCandidate& cand, const std::string& strval, Float_t& setval, Float_t defaultval){
  if (cand.hasUserFloat(strval)) setval = cand.userFloat(strval);
  else setval = defaultval;
}



// ------------ method called once each job just before starting event loop  ------------
void HZZ4lNtupleMaker::beginJob()
{
  edm::Service<TFileService> fs;
  myTree = new HZZ4lNtupleFactory( fs->make<TTree>(theFileName,"Event Summary"));
  const int nbins = 45;
  hCounter = fs->make<TH1F>("Counters", "Counters", nbins, 0., nbins);
  if (doreweighting) {
    hCounter_reweighted = fs->make<TH2F>("Counters_reweighted", "Counters_reweighted", nbins, 0., nbins, nReweightingSamples, 0, nReweightingSamples);
    couplingstree = fs->make<TTree>("couplings", "reweighting couplings");
    reweighting.fillcouplingstree(couplingstree);
  }
  BookAllBranches();
  buildMELABranches();
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
  hCounter->SetBinContent(8 ,gen_ZZ2emu2tau);
  hCounter->SetBinContent(9 ,gen_ZZ4tau);
  hCounter->SetBinContent(10,gen_ZZ4e_EtaAcceptance);
  hCounter->SetBinContent(11,gen_ZZ4e_LeptonAcceptance);
  hCounter->SetBinContent(14,gen_ZZ2mu2e_EtaAcceptance);
  hCounter->SetBinContent(15,gen_ZZ2mu2e_LeptonAcceptance);
  hCounter->SetBinContent(19,gen_BUGGY);
  hCounter->SetBinContent(20,gen_Unknown);

  hCounter->SetBinContent(40,gen_sumWeights); // Also stored in underflow bin; added here for convenience
  hCounter->SetBinContent(41,gen_sumGenMCWeight);
  hCounter->SetBinContent(42,gen_sumPUWeight);

  if (doreweighting) {
    for (int i = 0; i < nReweightingSamples; i++) {
      hCounter_reweighted->SetBinContent(0 , i+1, gen_sumWeights_reweighted[i]); // also stored in bin 40
      hCounter_reweighted->SetBinContent(1 , i+1, Nevt_Gen_reweighted[i]-gen_BUGGY_reweighted[i]);
      hCounter_reweighted->SetBinContent(2 , i+1, gen_ZZ4mu_reweighted[i]);
      hCounter_reweighted->SetBinContent(3 , i+1, gen_ZZ4e_reweighted[i]);
      hCounter_reweighted->SetBinContent(4 , i+1, gen_ZZ2mu2e_reweighted[i]);
      hCounter_reweighted->SetBinContent(5 , i+1, gen_ZZ2l2tau_reweighted[i]);
      hCounter_reweighted->SetBinContent(6 , i+1, gen_ZZ4mu_EtaAcceptance_reweighted[i]);
      hCounter_reweighted->SetBinContent(7 , i+1, gen_ZZ4mu_LeptonAcceptance_reweighted[i]);
      hCounter_reweighted->SetBinContent(8 , i+1, gen_ZZ2emu2tau_reweighted[i]);
      hCounter_reweighted->SetBinContent(9 , i+1, gen_ZZ4tau_reweighted[i]);
      hCounter_reweighted->SetBinContent(10, i+1, gen_ZZ4e_EtaAcceptance_reweighted[i]);
      hCounter_reweighted->SetBinContent(11, i+1, gen_ZZ4e_LeptonAcceptance_reweighted[i]);
      hCounter_reweighted->SetBinContent(14, i+1, gen_ZZ2mu2e_EtaAcceptance_reweighted[i]);
      hCounter_reweighted->SetBinContent(15, i+1, gen_ZZ2mu2e_LeptonAcceptance_reweighted[i]);
      hCounter_reweighted->SetBinContent(19, i+1, gen_BUGGY_reweighted[i]);
      hCounter_reweighted->SetBinContent(20, i+1, gen_Unknown_reweighted[i]);

      hCounter_reweighted->SetBinContent(40, i+1, gen_sumWeights_reweighted[i]); // Also stored in underflow bin; added here for convenience
      hCounter_reweighted->SetBinContent(41, i+1, gen_sumGenMCWeight_reweighted[i]);
      hCounter_reweighted->SetBinContent(42, i+1, gen_sumPUWeight);
    }
  }

  TH1 *h[2] = {hCounter, hCounter_reweighted};
  for (int i = 0; i < 2 && (doreweighting || i==0); i++) {
    h[i]->GetXaxis()->SetBinLabel(1 ,"Nevt_Gen");
    h[i]->GetXaxis()->SetBinLabel(2 ,"gen_ZZ4mu");
    h[i]->GetXaxis()->SetBinLabel(3 ,"gen_ZZ4e");
    h[i]->GetXaxis()->SetBinLabel(4 ,"gen_ZZ2mu2e");
    h[i]->GetXaxis()->SetBinLabel(5 ,"gen_ZZ2l2tau");
    h[i]->GetXaxis()->SetBinLabel(6 ,"gen_ZZ4mu_EtaAcceptance");
    h[i]->GetXaxis()->SetBinLabel(7 ,"gen_ZZ4mu_LeptonAcceptance");
    h[i]->GetXaxis()->SetBinLabel(8 ,"gen_ZZ2emu2tau");
    h[i]->GetXaxis()->SetBinLabel(9 ,"gen_ZZ4tau");
    h[i]->GetXaxis()->SetBinLabel(10,"gen_ZZ4e_EtaAcceptance");
    h[i]->GetXaxis()->SetBinLabel(11,"gen_ZZ4e_LeptonAcceptance");
    h[i]->GetXaxis()->SetBinLabel(14,"gen_ZZ2mu2e_EtaAcceptance");
    h[i]->GetXaxis()->SetBinLabel(15,"gen_ZZ2mu2e_LeptonAcceptance");
    h[i]->GetXaxis()->SetBinLabel(19,"gen_BUGGY");
    h[i]->GetXaxis()->SetBinLabel(20,"gen_Unknown");

    h[i]->GetXaxis()->SetBinLabel(40,"gen_sumWeights");
    h[i]->GetXaxis()->SetBinLabel(41,"gen_sumGenMCWeight");
    h[i]->GetXaxis()->SetBinLabel(42,"gen_sumPUWeight");
  }

  return;
}

// ------------ method called when starting to processes a run  ------------
void HZZ4lNtupleMaker::beginRun(edm::Run const& iRun, edm::EventSetup const&)
{
  // // code that helps find the indices of LHE weights
  // edm::Handle<LHERunInfoProduct> run;
  // typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
  // iRun.getByLabel( "externalLHEProducer", run );
  // LHERunInfoProduct myLHERunInfoProduct = *(run.product());
  // for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
  //   std::cout << iter->tag() << std::endl;
  //   std::vector<std::string> lines = iter->lines();
  //   for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
  //     std::cout << lines.at(iLine);
  //   }
  // }

}

// ------------ method called when ending the processing of a run  ------------
void HZZ4lNtupleMaker::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void HZZ4lNtupleMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
  Nevt_Gen_lumiBlock = 0;
  for (auto it = Nevt_Gen_lumiBlock_reweighted.begin(); it != Nevt_Gen_lumiBlock_reweighted.end(); ++it)
    *it = 0;
}

// ------------ method called when ending the processing of a luminosity block  ------------
void HZZ4lNtupleMaker::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup)
{
  Float_t Nevt_preskim = -1.;
  edm::Handle<edm::MergeableCounter> preSkimCounter;
  if (iLumi.getByToken(preSkimToken, preSkimCounter)) { // Counter before skim. Does not exist for non-skimmed samples.
    Nevt_preskim = preSkimCounter->value;
    // We do not use a filtering skim for the time being; so this is just left as a check in case we need it again in the future.
    if (Nevt_preskim>=0.) assert(Nevt_preskim == Nevt_Gen_lumiBlock);
  }

  Nevt_Gen += Nevt_Gen_lumiBlock;
  for (int i = 0; i < nReweightingSamples; i++) {
    Nevt_Gen_reweighted[i] += Nevt_Gen_lumiBlock_reweighted[i];
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


Float_t HZZ4lNtupleMaker::getAllWeight(const reco::Candidate* Lep) const
{
  Int_t   myLepID = abs(Lep->pdgId());
  if (skipMuDataMCWeight&& myLepID==13) return 1.;
  if (skipEleDataMCWeight&& myLepID==11) return 1.;
  if (myLepID==22) return 1.; // FIXME - what SFs should be used for TLEs?

  Float_t weight  = 1.;
  //Float_t errCorr = 0.;
  //Float_t errCorrSyst = 0.;

  Float_t myLepPt = Lep->pt();
  Float_t myLepEta = Lep->eta();
  Float_t mySIP = userdatahelpers::getUserFloat(Lep, "SIP"); 

  if(myLepID == 13){
    //avoid to go out of the TH boundary
    if(myLepPt > 79.) myLepPt = 79.;

    weight = hTH2D_Mu_All->GetBinContent(hTH2D_Mu_All->GetXaxis()->FindBin(myLepEta),hTH2D_Mu_All->GetYaxis()->FindBin(myLepPt));

  } else if(myLepID == 11) {

	// electron reconstruction scale factor, as a function of supercluster eta
	Float_t SCeta = userdatahelpers::getUserFloat(Lep,"SCeta");
	if(myLepPt < 20.) myLepPt = 20.;
    if(myLepPt > 199.) myLepPt = 199.;

	weight *= hTH2F_El_Reco->GetBinContent(hTH2F_El_Reco->GetXaxis()->FindBin(SCeta),hTH2F_El_Reco->GetYaxis()->FindBin(myLepPt));

    // reset pt for next lookup with different limits    
    myLepPt = Lep->pt();

    if(mySIP >= 4.0 ) {
        // No SF for RSE yet
        //return 1.;
    } else {

        if(myLepPt > 199.) myLepPt = 199.;
        Float_t myLepAbsEta = fabs(myLepEta);

        if(year >= 2016) {
            if((bool)userdatahelpers::getUserFloat(Lep,"isCrack"))
                 weight *= hTH2D_El_IdIsoSip_Cracks->GetBinContent(hTH2D_El_IdIsoSip_Cracks->FindFixBin(myLepAbsEta, myLepPt));
            else
                 weight *= hTH2D_El_IdIsoSip_notCracks->GetBinContent(hTH2D_El_IdIsoSip_notCracks->FindFixBin(myLepAbsEta, myLepPt));
        } else {
            if((bool)userdatahelpers::getUserFloat(Lep,"isCrack"))
              weight *= hTH2D_El_IdIsoSip_Cracks   ->GetBinContent(hTH2D_El_IdIsoSip_Cracks   ->GetXaxis()->FindBin(myLepPt),hTH2D_El_IdIsoSip_Cracks   ->GetYaxis()->FindBin(myLepAbsEta));
            else
              weight *= hTH2D_El_IdIsoSip_notCracks->GetBinContent(hTH2D_El_IdIsoSip_notCracks->GetXaxis()->FindBin(myLepPt),hTH2D_El_IdIsoSip_notCracks->GetYaxis()->FindBin(myLepAbsEta));
        }
    }
  } else {

    cout<<"ERROR! wrong lepton ID "<<myLepID<<endl;
    //abort();
    weight = 0.;
  }

  //FIXME
  if(myLepPt < 5. && myLepID == 13) weight = 1.;

  if(weight < 0.001 || weight > 10.){
    cout << "ERROR! LEP out of range! myLepPt = " << myLepPt << " myLepEta = " << myLepEta <<" myLepID "<<myLepID<< " weight = " << weight << endl;
    //abort();  //no correction should be zero, if you find one, stop
  }

  return weight;
}


Float_t HZZ4lNtupleMaker::getHqTWeight(double mH, double genPt) const
{
  if (skipHqTWeight) return 1.;

  //cout<<"mH = "<<mH<<", genPt = "<<genPt<<endl;
  if (mH<400 || genPt>250) return 1.;

  double weight = 1.;

  const int masses[4] = {400,600,800,1000};
  double massDiff = 1000;
  int iMass = -1;
  for (int i=0; i<4; ++i){
    double massDiffTmp = std::fabs(mH-masses[i]);
    if (massDiffTmp<massDiff){
      massDiff = massDiffTmp;
      iMass = i;
    }
  }

  if (iMass>=0) {
    weight = h_weight->GetBinContent(h_weight->FindBin(genPt));
  }
  return weight;
}


// Added by CO
Float_t HZZ4lNtupleMaker::getFakeWeight(Float_t LepPt, Float_t LepEta, Int_t LepID, Int_t LepZ1ID)
{
  if (skipFakeWeight) return 1.;

  // year 0 = 2011
  // year 1 = 2012

  Float_t weight  = 1.;

  Int_t nHisto=0;

  Float_t myLepPt   = LepPt;
  Float_t myLepEta  = fabs(LepEta);
  Int_t   myLepID   = abs(LepID);
  Int_t   myZ1LepID = abs(LepZ1ID);

  //cout << " pt = " << myLepPt << " eta = " << myLepEta << " ZID = " << myZ1LepID << " LepID = " << myLepID << endl;

  //avoid to go out of the TH boundary
  if(myLepPt > 79.) myLepPt = 79.;
  if(myLepID==13)nHisto+=2;
  if(myZ1LepID==13)nHisto+=1;


  TString Z1flavor = "Z1ee";     if(myZ1LepID==13) Z1flavor = "Z1mumu";
  TString Z2flavor = "electron"; if(myLepID==13)   Z2flavor = "muon";
  TString histo_name = "eff_"+Z1flavor+"_plus_"+Z2flavor;

  //cout << " histo = " << histo_name << endl;

  weight = h_ZXWeight[nHisto]->GetBinContent(h_ZXWeight[nHisto]->GetXaxis()->FindBin(myLepPt), h_ZXWeight[nHisto]->GetYaxis()->FindBin(myLepEta));
  /*
  h_ZXWeight[0]=FileZXWeightEle->Get("eff_Z1ee_plus_electron")->Clone();
  h_ZXWeight[2]=FileZXWeightEle->Get("eff_Z1mumu_plus_electron")->Clone();

  h_ZXWeight[5]=FileZXWeightMuo->Get("eff_Z1ee_plus_muon")->Clone();
  h_ZXWeight[7]=FileZXWeightMuo->Get("eff_Z1mumu_plus_muon")->Clone();
*/
  return weight;

} // end of getFakeWeight
void HZZ4lNtupleMaker::FillZGenInfo(Short_t Z1Id, Short_t Z2Id,
                                    const math::XYZTLorentzVector pZ1, const math::XYZTLorentzVector pZ2)
{
  GenZ1Mass= pZ1.M();
  GenZ1Pt= pZ1.Pt();
  GenZ1Eta= pZ1.Eta();
  GenZ1Phi= pZ1.Phi();
  GenZ1Flav= Z1Id;

  GenZ2Mass= pZ2.M();
  GenZ2Pt= pZ2.Pt();
  GenZ2Eta= pZ2.Eta();
  GenZ2Phi= pZ2.Phi();
  GenZ2Flav= Z2Id;

  return;
}

void HZZ4lNtupleMaker::FillLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, Short_t Lep3Id, Short_t Lep4Id,
                                      const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2,
                                      const math::XYZTLorentzVector Lep3, const math::XYZTLorentzVector Lep4)
{
  GenLep1Pt=Lep1.Pt();
  GenLep1Eta=Lep1.Eta();
  GenLep1Phi=Lep1.Phi();
  GenLep1Id=Lep1Id;

  GenLep2Pt=Lep2.Pt();
  GenLep2Eta=Lep2.Eta();
  GenLep2Phi=Lep2.Phi();
  GenLep2Id=Lep2Id;

  GenLep3Pt=Lep3.Pt();
  GenLep3Eta=Lep3.Eta();
  GenLep3Phi=Lep3.Phi();
  GenLep3Id=Lep3Id;

  GenLep4Pt=Lep4.Pt();
  GenLep4Eta=Lep4.Eta();
  GenLep4Phi=Lep4.Phi();
  GenLep4Id=Lep4Id;

  //can comment this back in if Gen angles are needed for any reason...
  //TUtil::computeAngles(zzanalysis::tlv(Lep1), Lep1Id, zzanalysis::tlv(Lep2), Lep2Id, zzanalysis::tlv(Lep3), Lep3Id, zzanalysis::tlv(Lep4), Lep4Id, Gencosthetastar, GenhelcosthetaZ1, GenhelcosthetaZ2, Genhelphi, GenphistarZ1);

  return;
}

void HZZ4lNtupleMaker::FillAssocLepGenInfo(std::vector<const reco::Candidate *>& AssocLeps)
{
  if (AssocLeps.size() >= 1) {
    GenAssocLep1Pt =AssocLeps.at(0)->p4().Pt();
    GenAssocLep1Eta=AssocLeps.at(0)->p4().Eta();
    GenAssocLep1Phi=AssocLeps.at(0)->p4().Phi();
    GenAssocLep1Id =AssocLeps.at(0)->pdgId();
  }
  if (AssocLeps.size() >= 2) {
    GenAssocLep2Pt =AssocLeps.at(1)->p4().Pt();
    GenAssocLep2Eta=AssocLeps.at(1)->p4().Eta();
    GenAssocLep2Phi=AssocLeps.at(1)->p4().Phi();
    GenAssocLep2Id =AssocLeps.at(1)->pdgId();
  }

  return;
}


void HZZ4lNtupleMaker::FillHGenInfo(const math::XYZTLorentzVector pH, float w)
{
  GenHMass=pH.M();
  GenHPt=pH.Pt();
  GenHRapidity=pH.Rapidity();

  HqTMCweight=w;

  return;
}


void HZZ4lNtupleMaker::BookAllBranches(){
   //Event variables
  myTree->Book("RunNumber",RunNumber);
  myTree->Book("EventNumber",EventNumber);
  myTree->Book("LumiNumber",LumiNumber);
  myTree->Book("NRecoMu",NRecoMu);
  myTree->Book("NRecoEle",NRecoEle);
  myTree->Book("Nvtx",Nvtx);
  myTree->Book("NObsInt",NObsInt);
  myTree->Book("NTrueInt",NTrueInt);

  myTree->Book("PFMET",PFMET);
  myTree->Book("PFMETPhi",PFMETPhi);
  myTree->Book("PFMETNoHF",PFMETNoHF);
  myTree->Book("PFMETNoHFPhi",PFMETNoHFPhi);
  myTree->Book("nCleanedJets",nCleanedJets);
  myTree->Book("nCleanedJetsPt30",nCleanedJetsPt30);
  myTree->Book("nCleanedJetsPt30_jecUp",nCleanedJetsPt30_jecUp);
  myTree->Book("nCleanedJetsPt30_jecDn",nCleanedJetsPt30_jecDn);
  myTree->Book("nCleanedJetsPt30BTagged",nCleanedJetsPt30BTagged);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF",nCleanedJetsPt30BTagged_bTagSF);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSFUp",nCleanedJetsPt30BTagged_bTagSFUp);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSFDn",nCleanedJetsPt30BTagged_bTagSFDn);
  myTree->Book("trigWord",trigWord);
  myTree->Book("ZZMass",ZZMass);
  myTree->Book("ZZMassErr",ZZMassErr);
  myTree->Book("ZZMassErrCorr",ZZMassErrCorr);
  myTree->Book("ZZMassPreFSR",ZZMassPreFSR);
  myTree->Book("ZZsel",ZZsel);
  myTree->Book("ZZPt",ZZPt);
  myTree->Book("ZZEta",ZZEta);
  myTree->Book("ZZPhi",ZZPhi);
  myTree->Book("CRflag",CRflag);
  myTree->Book("Z1Mass",Z1Mass);
  myTree->Book("Z1Pt",Z1Pt);
  myTree->Book("Z1Flav",Z1Flav);

  //Kin refitted info
  if (addKinRefit) {
    myTree->Book("ZZMassRefit",ZZMassRefit);
    myTree->Book("ZZMassRefitErr",ZZMassRefitErr);
    myTree->Book("ZZMassUnrefitErr",ZZMassUnrefitErr);
  }
  if (addVtxFit){
    myTree->Book("ZZMassCFit",ZZMassCFit);
    myTree->Book("ZZChi2CFit",ZZChi2CFit);
  }

  //Z2 variables
  myTree->Book("Z2Mass",Z2Mass);
  myTree->Book("Z2Pt",Z2Pt);
  myTree->Book("Z2Flav",Z2Flav);
  myTree->Book("costhetastar",costhetastar);
  myTree->Book("helphi",helphi);
  myTree->Book("helcosthetaZ1",helcosthetaZ1);
  myTree->Book("helcosthetaZ2",helcosthetaZ2);
  myTree->Book("phistarZ1",phistarZ1);
  myTree->Book("phistarZ2",phistarZ2);
  myTree->Book("xi",xi);
  myTree->Book("xistar",xistar);

  if (is_loose_ele_selection) {
    myTree->Book("TLE_dR_Z",TLE_dR_Z);
    myTree->Book("TLE_min_dR_3l",TLE_min_dR_3l);
  }

  myTree->Book("LepPt",LepPt);
  myTree->Book("LepEta",LepEta);
  myTree->Book("LepPhi",LepPhi);
  myTree->Book("LepLepId",LepLepId);
  myTree->Book("LepSIP",LepSIP);
  myTree->Book("LepTime",LepTime);
  myTree->Book("LepisID",LepisID);
  myTree->Book("LepisLoose",LepisLoose);
  myTree->Book("LepBDT",LepBDT);
  myTree->Book("LepMissingHit",LepMissingHit);
  //myTree->Book("LepChargedHadIso",LepChargedHadIso);
  //myTree->Book("LepNeutralHadIso",LepNeutralHadIso);
  //myTree->Book("LepPhotonIso",LepPhotonIso);
  myTree->Book("LepCombRelIsoPF",LepCombRelIsoPF);
  myTree->Book("fsrPt",fsrPt);
  myTree->Book("fsrEta",fsrEta);
  myTree->Book("fsrPhi",fsrPhi);
  myTree->Book("fsrLept",fsrLept);
  myTree->Book("passIsoPreFSR",passIsoPreFSR);
  if (addFSRDetails) {
    myTree->Book("fsrDR",fsrDR);
    myTree->Book("fsrLeptId",fsrLeptID);
    myTree->Book("fsrGenPt",fsrGenPt);
  }

  //Jet variables
  myTree->Book("JetPt",JetPt);
  myTree->Book("JetEta",JetEta);
  myTree->Book("JetPhi",JetPhi);
  myTree->Book("JetMass",JetMass);
  myTree->Book("JetBTagger",JetBTagger);
  myTree->Book("JetIsBtagged",JetIsBtagged);
  myTree->Book("JetIsBtaggedWithSF",JetIsBtaggedWithSF);
  myTree->Book("JetIsBtaggedWithSFUp",JetIsBtaggedWithSFUp);
  myTree->Book("JetIsBtaggedWithSFDn",JetIsBtaggedWithSFDn);
  myTree->Book("JetQGLikelihood",JetQGLikelihood);
  if(addQGLInputs){
    myTree->Book("JetAxis2",JetAxis2);
    myTree->Book("JetMult",JetMult);
    myTree->Book("JetPtD",JetPtD);
  }
  myTree->Book("JetSigma",JetSigma);
  myTree->Book("JetHadronFlavour",JetHadronFlavour);
  myTree->Book("DiJetMass",DiJetMass);
//   myTree->Book("DiJetMassPlus",DiJetMassPlus); // FIXME: add back once filled again
//   myTree->Book("DiJetMassMinus",DiJetMassMinus);
  myTree->Book("DiJetDEta",DiJetDEta);
  myTree->Book("DiJetFisher",DiJetFisher);
  myTree->Book("nExtraLep",nExtraLep);
  myTree->Book("nExtraZ",nExtraZ);
  myTree->Book("ExtraLepPt",ExtraLepPt);
  myTree->Book("ExtraLepEta",ExtraLepEta);
  myTree->Book("ExtraLepPhi",ExtraLepPhi);
  myTree->Book("ExtraLepLepId",ExtraLepLepId);

  myTree->Book("ZXFakeweight", ZXFakeweight);

  if (isMC){
    if (apply_K_NNLOQCD_ZZGG>0){
      myTree->Book("KFactor_QCD_ggZZ_Nominal", KFactor_QCD_ggZZ_Nominal);
      myTree->Book("KFactor_QCD_ggZZ_PDFScaleDn", KFactor_QCD_ggZZ_PDFScaleDn);
      myTree->Book("KFactor_QCD_ggZZ_PDFScaleUp", KFactor_QCD_ggZZ_PDFScaleUp);
      myTree->Book("KFactor_QCD_ggZZ_QCDScaleDn", KFactor_QCD_ggZZ_QCDScaleDn);
      myTree->Book("KFactor_QCD_ggZZ_QCDScaleUp", KFactor_QCD_ggZZ_QCDScaleUp);
      myTree->Book("KFactor_QCD_ggZZ_AsDn", KFactor_QCD_ggZZ_AsDn);
      myTree->Book("KFactor_QCD_ggZZ_AsUp", KFactor_QCD_ggZZ_AsUp);
      myTree->Book("KFactor_QCD_ggZZ_PDFReplicaDn", KFactor_QCD_ggZZ_PDFReplicaDn);
      myTree->Book("KFactor_QCD_ggZZ_PDFReplicaUp", KFactor_QCD_ggZZ_PDFReplicaUp);
    }
    if (apply_K_NLOEW_ZZQQB){
      myTree->Book("KFactor_EW_qqZZ", KFactor_EW_qqZZ);
      myTree->Book("KFactor_EW_qqZZ_unc", KFactor_EW_qqZZ_unc);
    }
    if (apply_K_NNLOQCD_ZZQQB){
      myTree->Book("KFactor_QCD_qqZZ_dPhi", KFactor_QCD_qqZZ_dPhi);
      myTree->Book("KFactor_QCD_qqZZ_M", KFactor_QCD_qqZZ_M);
      myTree->Book("KFactor_QCD_qqZZ_Pt", KFactor_QCD_qqZZ_Pt);
    }

    myTree->Book("genFinalState", genFinalState);
    myTree->Book("genProcessId", genProcessId);
    myTree->Book("genHEPMCweight", genHEPMCweight);
    myTree->Book("PUWeight", PUWeight);
    myTree->Book("dataMCWeight", dataMCWeight);
    myTree->Book("trigEffWeight", trigEffWeight);
    myTree->Book("overallEventWeight", overallEventWeight);
    myTree->Book("HqTMCweight", HqTMCweight);
    myTree->Book("xsec", xsection);
    myTree->Book("genExtInfo", genExtInfo);
    myTree->Book("GenHMass", GenHMass);
    myTree->Book("GenHPt", GenHPt);
    myTree->Book("GenHRapidity", GenHRapidity);
    myTree->Book("GenZ1Mass", GenZ1Mass);
    myTree->Book("GenZ1Pt", GenZ1Pt);
    myTree->Book("GenZ1Phi", GenZ1Phi);
    myTree->Book("GenZ1Flav", GenZ1Flav);
    myTree->Book("GenZ2Mass", GenZ2Mass);
    myTree->Book("GenZ2Pt", GenZ2Pt);
    myTree->Book("GenZ2Phi", GenZ2Phi);
    myTree->Book("GenZ2Flav", GenZ2Flav);
    myTree->Book("GenLep1Pt", GenLep1Pt);
    myTree->Book("GenLep1Eta", GenLep1Eta);
    myTree->Book("GenLep1Phi", GenLep1Phi);
    myTree->Book("GenLep1Id", GenLep1Id);
    myTree->Book("GenLep2Pt", GenLep2Pt);
    myTree->Book("GenLep2Eta", GenLep2Eta);
    myTree->Book("GenLep2Phi", GenLep2Phi);
    myTree->Book("GenLep2Id", GenLep2Id);
    myTree->Book("GenLep3Pt", GenLep3Pt);
    myTree->Book("GenLep3Eta", GenLep3Eta);
    myTree->Book("GenLep3Phi", GenLep3Phi);
    myTree->Book("GenLep3Id", GenLep3Id);
    myTree->Book("GenLep4Pt", GenLep4Pt);
    myTree->Book("GenLep4Eta", GenLep4Eta);
    myTree->Book("GenLep4Phi", GenLep4Phi);
    myTree->Book("GenLep4Id", GenLep4Id);
    myTree->Book("GenAssocLep1Pt", GenAssocLep1Pt);
    myTree->Book("GenAssocLep1Eta", GenAssocLep1Eta);
    myTree->Book("GenAssocLep1Phi", GenAssocLep1Phi);
    myTree->Book("GenAssocLep1Id", GenAssocLep1Id);
    myTree->Book("GenAssocLep2Pt", GenAssocLep2Pt);
    myTree->Book("GenAssocLep2Eta", GenAssocLep2Eta);
    myTree->Book("GenAssocLep2Phi", GenAssocLep2Phi);
    myTree->Book("GenAssocLep2Id", GenAssocLep2Id);

    if (doreweighting) myTree->Book("reweightingweights", reweightingweights);

    if (addLHEKinematics){
      myTree->Book("LHEMotherPz", LHEMotherPz);
      myTree->Book("LHEMotherE", LHEMotherE);
      myTree->Book("LHEMotherId", LHEMotherId);
      myTree->Book("LHEDaughterPt", LHEDaughterPt);
      myTree->Book("LHEDaughterEta", LHEDaughterEta);
      myTree->Book("LHEDaughterPhi", LHEDaughterPhi);
      myTree->Book("LHEDaughterMass", LHEDaughterMass);
      myTree->Book("LHEDaughterId", LHEDaughterId);
      myTree->Book("LHEAssociatedParticlePt", LHEAssociatedParticlePt);
      myTree->Book("LHEAssociatedParticleEta", LHEAssociatedParticleEta);
      myTree->Book("LHEAssociatedParticlePhi", LHEAssociatedParticlePhi);
      myTree->Book("LHEAssociatedParticleMass", LHEAssociatedParticleMass);
      myTree->Book("LHEAssociatedParticleId", LHEAssociatedParticleId);
    }

    myTree->Book("LHEPDFScale", LHEPDFScale);
    myTree->Book("LHEweight_QCDscale_muR1_muF1", LHEweight_QCDscale_muR1_muF1);
    myTree->Book("LHEweight_QCDscale_muR1_muF2", LHEweight_QCDscale_muR1_muF2);
    myTree->Book("LHEweight_QCDscale_muR1_muF0p5", LHEweight_QCDscale_muR1_muF0p5);
    myTree->Book("LHEweight_QCDscale_muR2_muF1", LHEweight_QCDscale_muR2_muF1);
    myTree->Book("LHEweight_QCDscale_muR2_muF2", LHEweight_QCDscale_muR2_muF2);
    myTree->Book("LHEweight_QCDscale_muR2_muF0p5", LHEweight_QCDscale_muR2_muF0p5);
    myTree->Book("LHEweight_QCDscale_muR0p5_muF1", LHEweight_QCDscale_muR0p5_muF1);
    myTree->Book("LHEweight_QCDscale_muR0p5_muF2", LHEweight_QCDscale_muR0p5_muF2);
    myTree->Book("LHEweight_QCDscale_muR0p5_muF0p5", LHEweight_QCDscale_muR0p5_muF0p5);
    myTree->Book("LHEweight_PDFVariation_Up", LHEweight_PDFVariation_Up);
    myTree->Book("LHEweight_PDFVariation_Dn", LHEweight_PDFVariation_Dn);
    myTree->Book("LHEweight_AsMZ_Up", LHEweight_AsMZ_Up);
    myTree->Book("LHEweight_AsMZ_Dn", LHEweight_AsMZ_Dn);
  }

  // MELA branches are booked under buildMELA
}

void HZZ4lNtupleMaker::addweight(float &weight, vector<float> &weight_reweighted, float weighttoadd) {
  //  cout << weight_reweighted.size() << endl;
  weight += weighttoadd;
  if (doreweighting)
    for (int i = 0; i < nReweightingSamples; i++)
      weight_reweighted[i] += weighttoadd * reweightingweights[i];
}


void HZZ4lNtupleMaker::buildMELABranches(){
  /***********************/
  /***********************/
  /**   Reco branches   **/
  /***********************/
  /***********************/
  for (unsigned int it=0; it<recoMElist.size(); it++){
    MELAOptionParser* me_opt = new MELAOptionParser(recoMElist.at(it));
    if (recoMElist.at(it).find("Copy")!=string::npos) recome_copyopts.push_back(me_opt);
    else recome_originalopts.push_back(me_opt);
  }
  // Resolve original options
  for (unsigned int it=0; it<recome_originalopts.size(); it++){
    MELAOptionParser* me_opt = recome_originalopts.at(it);
    myTree->BookMELABranches(me_opt, false, 0);
  }
  // Resolve copy options
  for (unsigned int it=0; it<recome_copyopts.size(); it++){
    MELAOptionParser* me_opt = recome_copyopts.at(it);
    MELAOptionParser* original_opt=0;
    // Find the original options
    for (unsigned int ih=0; ih<recome_originalopts.size(); ih++){
      if (me_opt->testCopyAlias(recome_originalopts.at(ih)->getAlias())){
        original_opt = recome_originalopts.at(ih);
        break;
      }
    }
    if (original_opt==0) continue;
    else me_opt->pickOriginalOptions(original_opt);
    myTree->BookMELABranches(me_opt, false, 0);
  }

  /**********************/
  /**********************/
  /**   LHE branches   **/
  /**********************/
  /**********************/
  for (unsigned int it=0; it<lheMElist.size(); it++){
    MELAOptionParser* lheme_opt;
    // First find out if the option has a copy specification
    // These copy options will be evaulated in a separate loop
    if (lheMElist.at(it).find("Copy")!=string::npos){
      lheme_opt = new MELAOptionParser(lheMElist.at(it));
      lheme_copyopts.push_back(lheme_opt);
      continue;
    }

    // Create a hypothesis for each option
    MELAHypothesis* lheme_hypo = new MELAHypothesis(&mela, lheMElist.at(it));
    lheme_units.push_back(lheme_hypo);

    lheme_opt = lheme_hypo->getOption();
    if (lheme_opt->isAliased()) lheme_aliased_units.push_back(lheme_hypo);

    // Create a computation for each hypothesis
    MELAComputation* lheme_computer = new MELAComputation(lheme_hypo);
    lheme_computers.push_back(lheme_computer);

    // Add the computation to a named cluster to keep track of JECUp/JECDn, or for best-pWH_SM Lep_WH computations
    addToMELACluster(lheme_computer, lheme_clusters);

    // Create the necessary branches for each computation
    myTree->BookMELABranches(lheme_opt, true, lheme_computer);
  }
  // Resolve copy options
  for (unsigned int it=0; it<lheme_copyopts.size(); it++){
    MELAOptionParser* lheme_opt = lheme_copyopts.at(it);
    MELAHypothesis* original_hypo=0;
    MELAOptionParser* original_opt=0;
    // Find the original options
    for (unsigned int ih=0; ih<lheme_aliased_units.size(); ih++){
      if (lheme_opt->testCopyAlias(lheme_aliased_units.at(ih)->getOption()->getAlias())){
        original_hypo = lheme_aliased_units.at(ih);
        original_opt = original_hypo->getOption();
        break;
      }
    }
    if (original_opt==0) continue;
    else lheme_opt->pickOriginalOptions(original_opt);
    // Create a new computation for the copy options
    MELAComputation* lheme_computer = new MELAComputation(original_hypo);
    lheme_computer->setOption(lheme_opt);
    lheme_computers.push_back(lheme_computer);

    // The rest is the same story...
    // Add the computation to a named cluster to keep track of JECUp/JECDn, or for best-pWH_SM Lep_WH computations
    addToMELACluster(lheme_computer, lheme_clusters);

    // Create the necessary branches for each computation
    myTree->BookMELABranches(lheme_opt, true, lheme_computer);
  }
  // Loop over the computations to add any contingencies to aliased hypotheses
  for (unsigned int it=0; it<lheme_computers.size(); it++) lheme_computers.at(it)->addContingencies(lheme_aliased_units);

  if (DEBUG_MB){
    std::vector<MELABranch*>* lheme_branches = myTree->getLHEMELABranches();
    for (unsigned int ib=0; ib<lheme_branches->size(); ib++) lheme_branches->at(ib)->Print();
  }
}

void HZZ4lNtupleMaker::addToMELACluster(MELAComputation* me_computer, std::vector<MELACluster*>& me_clusters){
  bool isAdded=false;
  for (unsigned int it=0; it<me_clusters.size(); it++){ if (me_clusters.at(it)->getName()==me_computer->getName()){ me_clusters.at(it)->addComputation(me_computer); isAdded=true; } }
  if (!isAdded){
    MELACluster* tmpcluster = new MELACluster(me_computer->getCluster());
    tmpcluster->addComputation(me_computer);
    me_clusters.push_back(tmpcluster);
  }
}

void HZZ4lNtupleMaker::computeMELABranches(MELACandidate* cand){
  mela.setCurrentCandidate(cand);
  updateMELAClusters_Common(); // "Common"
  mela.resetInputEvent();
}
// Common ME computations with index=0
void HZZ4lNtupleMaker::updateMELAClusters_Common(){
  MELACandidate* melaCand = mela.getCurrentCandidate();
  if (melaCand==0) return;

  for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
    MELACluster* theCluster = lheme_clusters.at(ic);
    if (theCluster->getName()=="Common"){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }
}
void HZZ4lNtupleMaker::pushRecoMELABranches(const pat::CompositeCandidate& cand){
  std::vector<MELABranch*>* recome_branches = myTree->getRecoMELABranches();
  // Pull + push...
  for (unsigned int ib=0; ib<recome_branches->size(); ib++){
    std::string branchname = recome_branches->at(ib)->bname.Data();
    if (cand.hasUserFloat(branchname)) recome_branches->at(ib)->setValue((Float_t)cand.userFloat(branchname));
    else cerr << "HZZ4lNtupleMaker::pushRecoMELABranches: Candidate does not contain the reco ME " << branchname << " it should have calculated!" << endl;
  }
}
void HZZ4lNtupleMaker::pushLHEMELABranches(){
  std::vector<MELABranch*>* lheme_branches = myTree->getLHEMELABranches();
  // Pull + push...
  for (unsigned int ib=0; ib<lheme_branches->size(); ib++) lheme_branches->at(ib)->setVal();
  // ...then reset
  for (unsigned int ic=0; ic<lheme_clusters.size(); ic++) lheme_clusters.at(ic)->reset();
}
void HZZ4lNtupleMaker::clearMELABranches(){
  for (unsigned int it=0; it<lheme_clusters.size(); it++) delete lheme_clusters.at(it);
  for (unsigned int it=0; it<lheme_computers.size(); it++) delete lheme_computers.at(it);
  for (unsigned int it=0; it<lheme_copyopts.size(); it++) delete lheme_copyopts.at(it);
  //for (unsigned int it=0; it<lheme_aliased_units.size(); it++) delete lheme_aliased_units.at(it); // DO NOT DELETE THIS, WILL BE DELETED WITH lheme_units!
  for (unsigned int it=0; it<lheme_units.size(); it++) delete lheme_units.at(it);

  for (unsigned int it=0; it<recome_copyopts.size(); it++) delete recome_copyopts.at(it);
  for (unsigned int it=0; it<recome_originalopts.size(); it++) delete recome_originalopts.at(it);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HZZ4lNtupleMaker);
