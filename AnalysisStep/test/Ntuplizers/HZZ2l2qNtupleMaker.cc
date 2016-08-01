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
#include <ZZAnalysis/AnalysisStep/interface/bitops.h>
#include <ZZAnalysis/AnalysisStep/interface/Fisher.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>
#include <ZZAnalysis/AnalysisStep/interface/JetCleaner.h>

#include "ZZ2l2qConfigHelper.h"
#include "HZZ4lNtupleFactory.h"

#include <TRandom3.h>
#include <TH2D.h>
#include "TLorentzVector.h"
#include "TSpline.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"


#include <string>

namespace {
  bool writeJets = true;     // Write jets in the tree. FIXME: make this configurable
  bool addKinRefit = true;
  bool addVtxFit = false;
  bool addFSRDetails = false;
  bool skipDataMCWeight = false; // skip computation of data/MC weight
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
  Float_t KFactorggZZ = 0;
  Float_t KFactorEWKqqZZ = 0;
  Float_t KFactorQCDqqZZ_dPhi = 0;
  Float_t KFactorQCDqqZZ_M = 0;
  Float_t KFactorQCDqqZZ_Pt = 0;
  Float_t PFMET  =  -99;
  Float_t PFMETPhi  =  -99;
  Float_t PFMETNoHF  =  -99;
  Float_t PFMETNoHFPhi  =  -99;
  Short_t nCleanedJets  =  0;
  Short_t nCleanedJetsPt30  = 0;
  Short_t nCleanedJetsPt30_jecUp  = 0;
  Short_t nCleanedJetsPt30_jecDn  = 0;
  Short_t nCleanedJetsPt30BTagged  = 0;
  Short_t trigWord  = 0;
  std::vector<float> ZZMass;
  std::vector<float> ZZMass_JecUp;
  std::vector<float> ZZMass_JecDown;
  std::vector<float> ZZMassErr;
  std::vector<float> ZZMassErrCorr;
  std::vector<float> ZZMassPreFSR;
  std::vector<short> ZZsel;
  std::vector<float> ZZPt;
  std::vector<float> ZZEta;
  std::vector<float> ZZPhi;
  std::vector<short> CRflag;
  std::vector<short> ZZCandType;         // 1 = Merged, SR
                                         // -1 = Merged, SBR
                                         // 2 = Resloved, SR
                                          // -2 = Resloved, SBR
  std::vector<float> Z1Mass;
  std::vector<float> Z1Mass_JecUp;
  std::vector<float> Z1Mass_JecDown;
  std::vector<float> Z1Pt;
  std::vector<short> Z1Flav;
  std::vector<float> Z1tau21;
  std::vector<float> ZZMassRefit;
  std::vector<float> Z1MassRefit;
  std::vector<float> ZZMassUnrefitErr;
  std::vector<float> ZZMassCFit;
  std::vector<float> ZZChi2CFit;
  std::vector<float> Z2Mass;
  std::vector<float> Z2Pt;
  std::vector<short> Z2Flav;
  std::vector<float> costhetastar;
  std::vector<float> helphi;
  std::vector<float> helcosthetaZ1;
  std::vector<float> helcosthetaZ2;
  std::vector<float> phistarZ1;
  std::vector<float> phistarZ2;
  std::vector<float> xi;
  std::vector<float> xistar;
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
  std::vector<float> fsrPt;
  std::vector<float> fsrEta;
  std::vector<float> fsrPhi;
  std::vector<float> fsrDR;
  std::vector<short> fsrLept;
  std::vector<short> fsrLeptID;
  std::vector<float> fsrGenPt;
  Bool_t passIsoPreFSR;

  // JHUGen MEs
  std::vector<float> p0plus_VAJHU;
  std::vector<float> p0minus_VAJHU;
  std::vector<float> p0hplus_VAJHU;
  std::vector<float> p2bplus_VAJHU;
  std::vector<float> p2_VAJHU;
  std::vector<float> phjj_VAJHU_highestPTJets;
  std::vector<float> pvbf_VAJHU_highestPTJets;
  std::vector<float> phjj_VAJHU_highestPTJets_up;
  std::vector<float> pvbf_VAJHU_highestPTJets_up;
  std::vector<float> phjj_VAJHU_highestPTJets_dn;
  std::vector<float> pvbf_VAJHU_highestPTJets_dn;

  // MCFM MEs
  std::vector<float> pqqZJJ_VAMCFM;
  std::vector<float> p0plus_VAMCFM;
  std::vector<float> bkg_VAMCFM;
  std::vector<float> ggzz_VAMCFM;
  std::vector<float> ggzz_p0plus_VAMCFM;
  std::vector<float> Dgg10_VAMCFM;

  /*
  std::vector<float> p0plus_VAJHU;
  std::vector<float> p0minus_VAJHU;
  std::vector<float> p0plus_VAMCFM;
  std::vector<float> p0hplus_VAJHU;
  std::vector<float> p1_VAJHU;
  std::vector<float> p1_prodIndep_VAJHU;
  std::vector<float> p1plus_VAJHU;
  std::vector<float> p1plus_prodIndep_VAJHU;
  std::vector<float> p2_VAJHU;
  std::vector<float> p2_prodIndep_VAJHU;
  std::vector<float> p2qqb_VAJHU;
  std::vector<float> p2hplus_VAJHU;
  std::vector<float> p2hminus_VAJHU;
  std::vector<float> p2bplus_VAJHU;
  std::vector<float> p2hplus_qqb_VAJHU;
  std::vector<float> p2hplus_prodIndep_VAJHU;
  std::vector<float> p2hminus_qqb_VAJHU;
  std::vector<float> p2hminus_prodIndep_VAJHU;
  std::vector<float> p2bplus_qqb_VAJHU;
  std::vector<float> p2bplus_prodIndep_VAJHU;
  std::vector<float> p2h2plus_gg_VAJHU;
  std::vector<float> p2h2plus_qqbar_VAJHU;
  std::vector<float> p2h2plus_prodIndep_VAJHU;
  std::vector<float> p2h3plus_gg_VAJHU;
  std::vector<float> p2h3plus_qqbar_VAJHU;
  std::vector<float> p2h3plus_prodIndep_VAJHU;
  std::vector<float> p2h6plus_gg_VAJHU;
  std::vector<float> p2h6plus_qqbar_VAJHU;
  std::vector<float> p2h6plus_prodIndep_VAJHU;
  std::vector<float> p2h7plus_gg_VAJHU;
  std::vector<float> p2h7plus_qqbar_VAJHU;
  std::vector<float> p2h7plus_prodIndep_VAJHU;
  std::vector<float> p2h9minus_gg_VAJHU;
  std::vector<float> p2h9minus_qqbar_VAJHU;
  std::vector<float> p2h9minus_prodIndep_VAJHU;
  std::vector<float> p2h10minus_gg_VAJHU;
  std::vector<float> p2h10minus_qqbar_VAJHU;
  std::vector<float> p2h10minus_prodIndep_VAJHU;
  std::vector<float> bkg_VAMCFM;
  std::vector<float> bkg_prodIndep_VAMCFM;
  std::vector<float> ggzz_VAMCFM;
  std::vector<float> ggzz_p0plus_VAMCFM;
  std::vector<float> ggzz_c1_VAMCFM;
  std::vector<float> ggzz_c5_VAMCFM;
  std::vector<float> ggzz_ci_VAMCFM;
  std::vector<float> pg1g4_mela;
  std::vector<float> pg1g4_VAJHU;
  std::vector<float> pg1g4_pi2_VAJHU;
  std::vector<float> pg1g2_pi2_VAJHU;
  std::vector<float> pg1g2_mela;
  std::vector<float> pg1g2_VAJHU;
  std::vector<float> p0_g1prime2_VAJHU;
  std::vector<float> pg1g1prime2_VAJHU;
  std::vector<float> Dgg10_VAMCFM;
  std::vector<float> pzzzg_VAJHU;
  std::vector<float> pzzgg_VAJHU;
  std::vector<float> pzzzg_PS_VAJHU;
  std::vector<float> pzzgg_PS_VAJHU;
  std::vector<float> p0Zgs_VAJHU;
  std::vector<float> p0gsgs_VAJHU;
  std::vector<float> p0Zgs_PS_VAJHU;
  std::vector<float> p0gsgs_PS_VAJHU;
  std::vector<float> p0Zgs_g1prime2_VAJHU;
  std::vector<float> pzzzg_g1prime2_VAJHU;
  std::vector<float> pzzzg_g1prime2_pi2_VAJHU;
  std::vector<float> p0plus_m2l2q;
  std::vector<float> bkg_m2l2q;
  std::vector<float> p0plus_m2l2q_ScaleUp;
  std::vector<float> bkg_m2l2q_ScaleUp;
  std::vector<float> p0plus_m2l2q_ScaleDown;
  std::vector<float> bkg_m2l2q_ScaleDown;
  std::vector<float> p0plus_m2l2q_ResUp;
  std::vector<float> bkg_m2l2q_ResUp;
  std::vector<float> p0plus_m2l2q_ResDown;
  std::vector<float> bkg_m2l2q_ResDown;
  std::vector<float> phjj_VAJHU_highestPTJets;
  std::vector<float> pvbf_VAJHU_highestPTJets;
  std::vector<float> phjj_VAJHU_highestPTJets_up;
  std::vector<float> pvbf_VAJHU_highestPTJets_up;
  std::vector<float> phjj_VAJHU_highestPTJets_dn;
  std::vector<float> pvbf_VAJHU_highestPTJets_dn;
  std::vector<float> phjj_VAJHU_bestDjet;
  std::vector<float> pvbf_VAJHU_bestDjet;
  std::vector<float> phjj_VAJHU_bestDjet_up;
  std::vector<float> pvbf_VAJHU_bestDjet_up;
  std::vector<float> phjj_VAJHU_bestDjet_dn;
  std::vector<float> pvbf_VAJHU_bestDjet_dn;
  std::vector<float> pAux_vbf_VAJHU;
  std::vector<float> pAux_vbf_VAJHU_up;
  std::vector<float> pAux_vbf_VAJHU_dn;
  std::vector<float> phj_VAJHU;
  std::vector<float> phj_VAJHU_up;
  std::vector<float> phj_VAJHU_dn;
  std::vector<float> pwh_hadronic_VAJHU;
  std::vector<float> pwh_hadronic_VAJHU_up;
  std::vector<float> pwh_hadronic_VAJHU_dn;
  std::vector<float> pzh_hadronic_VAJHU;
  std::vector<float> pzh_hadronic_VAJHU_up;
  std::vector<float> pzh_hadronic_VAJHU_dn;
  std::vector<float> ptth_VAJHU;
  std::vector<float> ptth_VAJHU_up;
  std::vector<float> ptth_VAJHU_dn;
  std::vector<float> pbbh_VAJHU;
  std::vector<float> pbbh_VAJHU_up;
  std::vector<float> pbbh_VAJHU_dn;
  */
  std::vector<float> JetPt;
  std::vector<float> JetEta;
  std::vector<float> JetPhi;
  std::vector<float> JetMass;
  std::vector<float> JetBTagger;
  std::vector<float> JetIsBtagged;
  std::vector<float> JetQGLikelihood;
  std::vector<float> JetSigma;
  std::vector<bool>  JetIsInZZCand;
  Float_t DiJetMass  = -99;
  //Float_t DiJetMassPlus  = -99;
  //Float_t DiJetMassMinus  = -99;
  Float_t DiJetDEta  = -99;
  Float_t DiJetFisher  = -99;
  Short_t nExtraLep  = 0;
  Short_t nExtraZ  = 0;
  std::vector<float> ExtraLepPt;
  std::vector<float> ExtraLepEta;
  std::vector<float> ExtraLepPhi;
  std::vector<short> ExtraLepLepId;
  Short_t genFinalState  = 0;
  Int_t genProcessId  = 0;
  Float_t genHEPMCweight  = 0;
  Float_t LHEweight_QCDscale_muR1_muF1  = 0;
  Float_t LHEweight_QCDscale_muR1_muF2  = 0;
  Float_t LHEweight_QCDscale_muR1_muF0p5  = 0;
  Float_t LHEweight_QCDscale_muR2_muF1  = 0;
  Float_t LHEweight_QCDscale_muR2_muF2  = 0;
  Float_t LHEweight_QCDscale_muR2_muF0p5  = 0;
  Float_t LHEweight_QCDscale_muR0p5_muF1  = 0;
  Float_t LHEweight_QCDscale_muR0p5_muF2  = 0;
  Float_t LHEweight_QCDscale_muR0p5_muF0p5  = 0;
  Short_t genExtInfo  = 0;
  Float_t xsection  = 0;
  Float_t dataMCWeight  = 0;
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
  /*Float_t GenLep3Pt  = 0;
  Float_t GenLep3Eta  = 0;
  Float_t GenLep3Phi  = 0;
  Short_t GenLep3Id  = 0;
  Float_t GenLep4Pt  = 0;
  Float_t GenLep4Eta  = 0;
  Float_t GenLep4Phi  = 0;
  Short_t GenLep4Id  = 0;*/
  Float_t GenAssocLep1Pt  = 0;
  Float_t GenAssocLep1Eta  = 0;
  Float_t GenAssocLep1Phi  = 0;
  Short_t GenAssocLep1Id  = 0;
  Float_t GenAssocLep2Pt  = 0;
  Float_t GenAssocLep2Eta  = 0;
  Float_t GenAssocLep2Phi  = 0;
  Short_t GenAssocLep2Id  = 0;
}

using namespace std;
using namespace edm;
//
// class declaration
//
class HZZ2l2qNtupleMaker : public edm::EDAnalyzer {
public:
  explicit HZZ2l2qNtupleMaker(const edm::ParameterSet&);
  ~HZZ2l2qNtupleMaker();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  void BookAllBranches();  
  virtual void FillCandidate(const pat::CompositeCandidate& higgs, bool evtPass, const edm::Event&, const Int_t CRflag, bool isMerged);
  virtual void FillJet(const pat::Jet& jet);
  virtual void endJob() ;

  void FillHGenInfo(const math::XYZTLorentzVector Hp, float w);
  void FillZGenInfo(Short_t Z1Id, Short_t Z2Id, 
		    const math::XYZTLorentzVector pZ1, const math::XYZTLorentzVector pZ2);
  void FillLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, Short_t Lep3Id, Short_t Lep4Id, 
    const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2, const math::XYZTLorentzVector Lep3, const math::XYZTLorentzVector Lep4);
  void FillAssocLepGenInfo(std::vector<const reco::Candidate *>& AssocLeps);

  Float_t getAllWeight(const reco::Candidate* Lep) const;
  Float_t getAllWeightJet(const reco::Candidate* Jet) const;
  Float_t getHqTWeight(double mH, double genPt) const;
  Float_t getFakeWeight(Float_t LepPt, Float_t LepEta, Int_t LepID, Int_t LepZ1ID);

  // ----------member data ---------------------------
  ZZ2l2qConfigHelper myHelper;
  int theChannel;
  std::string theCandLabel;
  std::string theCandLabelFat;
  TString theFileName;

  HZZ4lNtupleFactory *myTree;
  TH1F *hCounter;

  Bool_t isMC;

  bool applyTrigger;    // Keep only events passing trigger
  bool applySkim;       //   "     "      "     skim
  bool skipEmptyEvents; // Skip events whith no selected candidate (otherwise, gen info is preserved for all events)
  Float_t xsec;
  int year;

  edm::EDGetTokenT<edm::View<reco::Candidate> > genParticleToken;
  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken;
  edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > candToken;
  edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > candTokenFat;
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

  string sampleName;

  std::vector<const reco::Candidate *> genFSR;

  std::vector<std::vector<float>> ewkTable;
  TSpline3* spkfactor;

  TH2D *hTH2D_Mu_All;
  TH2F *hTH2F_El_Reco;
  TH2D *hTH2D_El_IdIsoSip_notCracks;
  TH2D *hTH2D_El_IdIsoSip_Cracks;
  TH2D* h_weight; //HqT weights
  //TH2F *h_ZXWeightMuo;
  //TH2F *h_ZXWeightEle;
  TH2D* h_ZXWeight[4];

  BTagCalibration* calib;
  BTagCalibrationReader* reader;     // cntrl
  BTagCalibrationReader* reader_up;  // sys up
  BTagCalibrationReader* reader_do;  // sys down

};

//
// constructors and destructor
//
HZZ2l2qNtupleMaker::HZZ2l2qNtupleMaker(const edm::ParameterSet& pset) :
  myHelper(pset),
  reweight(),
  hTH2D_Mu_All(0),
  hTH2F_El_Reco(0),
  hTH2D_El_IdIsoSip_notCracks(0),
  hTH2D_El_IdIsoSip_Cracks(0),
  h_weight(0)
{
  //cout<< "Beginning Constructor\n\n\n" <<endl;
  theCandLabel = pset.getUntrackedParameter<string>("CandCollection"); // Name of input ZZ collection
  theCandLabelFat = pset.getUntrackedParameter<string>("CandCollectionMerged"); // Name of input ZZ collection
  theChannel = myHelper.channel(); // Valid options: ZZ, ZLL, ZL 
  theFileName = pset.getUntrackedParameter<string>("fileName"); 
  skipEmptyEvents = pset.getParameter<bool>("skipEmptyEvents"); // Do not store 
  sampleName = pset.getParameter<string>("sampleName");
  xsec = pset.getParameter<double>("xsec");
  year = pset.getParameter<int>("setup");

  consumesMany<std::vector< PileupSummaryInfo > >();
  genParticleToken = consumes<edm::View<reco::Candidate> >(edm::InputTag("prunedGenParticles"));
  genInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  consumesMany<LHEEventProduct>();
  candToken = consumes<edm::View<pat::CompositeCandidate> >(edm::InputTag(theCandLabel));
  candTokenFat = consumes<edm::View<pat::CompositeCandidate> >(edm::InputTag(theCandLabelFat));
  triggerResultToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults"));
  vtxToken = consumes<vector<reco::Vertex> >(edm::InputTag("goodPrimaryVertices"));
  jetToken = consumes<edm::View<pat::Jet> >(edm::InputTag("cleanJets"));
  metToken     = consumes<pat::METCollection>(edm::InputTag("slimmedMETs"));
  metNoHFToken = consumes<pat::METCollection>(edm::InputTag("slimmedMETsNoHF"));
  muonToken = consumes<pat::MuonCollection>(edm::InputTag("slimmedMuons"));
  electronToken = consumes<pat::ElectronCollection>(edm::InputTag("slimmedElectrons"));

  preSkimToken = consumes<edm::MergeableCounter,edm::InLumi>(edm::InputTag("preSkimCounter"));

  isMC = myHelper.isMC();

  if (skipEmptyEvents) {
    applyTrigger=true;         // Trigger does not work in 80X
    if(isMC) applyTrigger=false;
    applySkim=true;
  } else {
    applyTrigger=false;
    applySkim=false;    
  }

  Nevt_Gen = 0;
  Nevt_Gen_lumiBlock = 0;

  //For Efficiency studies
  gen_ZZ4mu = 0;
  gen_ZZ4e = 0;
  gen_ZZ2mu2e = 0;
  gen_ZZ2l2tau = 0;
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
  edm::FileInPath ggzzFIP("ZZAnalysis/AnalysisStep/data/kfactors/Kfactor_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
  fipPath=ggzzFIP.fullPath();
  TFile* ggZZKFactorFile = TFile::Open(fipPath.data());
  spkfactor = (TSpline3*)ggZZKFactorFile->Get("sp_Kfactor")->Clone();
  ggZZKFactorFile->Close();
  
  if (!skipDataMCWeight) {
    //Scale factors for data/MC efficiency

    // only LEPTON_SETUP=2015 is supported
    /* if(year!=2015){
      cout<<"Error: data/MC weights are not supported for this LEPTON_SETUP; see HZZ2l2qNtupleMaker.cc"<<endl;
      abort();
      }*/

    TString filename;

    filename.Form("ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ScaleFactors_mu_%d.root",year);
    edm::FileInPath fipMu(filename.Data());
    fipPath = fipMu.fullPath();
    TFile *fMuWeight = TFile::Open(fipPath.data(),"READ");
    hTH2D_Mu_All = (TH2D*)fMuWeight->Get("FINAL")->Clone();
    fMuWeight->Close();

    filename.Form("ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ele_scale_factors_%d_v3.root",year);
    edm::FileInPath fipEle(filename.Data());
    fipPath = fipEle.fullPath();
    TFile *fEleWeight = TFile::Open(fipPath.data(),"READ");
    hTH2D_El_IdIsoSip_notCracks =  (TH2D*)fEleWeight->Get("ele_scale_factors")->Clone();
    hTH2D_El_IdIsoSip_Cracks =  (TH2D*)fEleWeight->Get("ele_scale_factors_gap")->Clone();
    fEleWeight->Close();

    edm::FileInPath fipEleReco("ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/egammaEffi.txt_SF2D.root");
    fipPath = fipEleReco.fullPath();
    TFile *fEleReco = TFile::Open(fipPath.data(),"READ");
    hTH2F_El_Reco = (TH2F*)fEleReco->Get("EGamma_SF2D")->Clone();
    fEleReco->Close();

    edm::FileInPath btagFIP("ZZAnalysis/AnalysisStep/data/CSVv2_4invfb_systJuly15.csv");
    fipPath=btagFIP.fullPath();
    calib = new BTagCalibration("csvv2",fipPath);
    reader = new BTagCalibrationReader(calib,BTagEntry::OP_LOOSE,"comb","central");
    reader_up = new BTagCalibrationReader(calib,BTagEntry::OP_LOOSE,"comb","up");
    reader_do = new BTagCalibrationReader(calib,BTagEntry::OP_LOOSE,"comb","down");
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

HZZ2l2qNtupleMaker::~HZZ2l2qNtupleMaker()
{
}


// ------------ method called for each event  ------------
void HZZ2l2qNtupleMaker::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  myTree->InitializeVariables();

  //----------------------------------------------------------------------
  // Analyze MC truth; collect MC weights and update counters (this is done for all generated events, 
  // including those that do not pass skim, trigger etc!)

  bool gen_ZZ2l2qInEtaAcceptance = false;   // All 4 gen leptons in eta acceptance
  bool gen_ZZ2l2qInEtaPtAcceptance = false; // All 4 gen leptons in eta,pT acceptance

  const reco::Candidate * genH = 0;
  std::vector<const reco::Candidate *> genZLeps;
  std::vector<const reco::Candidate *> genAssocLeps;

  edm::Handle<edm::View<reco::Candidate> > genParticles;
  edm::Handle<GenEventInfoProduct> genInfo;

  if (isMC) {
    // get PU weights
    vector<Handle<std::vector< PileupSummaryInfo > > >  PupInfos; //FIXME support for miniAOD v1/v2 where name changed; catch does not work...
    event.getManyByType(PupInfos);
    Handle<std::vector< PileupSummaryInfo > > PupInfo = PupInfos.front(); 
//     try {
//       cout << "TRY HZZ2l2qNtupleMaker" <<endl;
//       event.getByLabel(edm::InputTag("addPileupInfo"), PupInfo); 
//     } catch (const cms::Exception& e){
//       cout << "FAIL HZZ2l2qNtupleMaker" <<endl;
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

    MCHistoryTools mch(event, sampleName, genParticles, genInfo, false);
    genFinalState = mch.genFinalState();
    genProcessId = mch.getProcessID();
    genHEPMCweight = mch.gethepMCweight();
    genExtInfo = mch.genAssociatedFS();

    // keep track of sum of weights
    gen_sumPUWeight    += PUWeight;
    gen_sumGenMCWeight += genHEPMCweight;
    gen_sumWeights     += PUWeight*genHEPMCweight;

    // LHE weights
    vector<edm::Handle<LHEEventProduct> > EvtHandles;
    event.getManyByType(EvtHandles); // using this method because the label is not always the same (e.g. "source" in the ttH sample)
    if(EvtHandles.size()>0){
      edm::Handle<LHEEventProduct> EvtHandle = EvtHandles.front();
      if(EvtHandle.isValid() && EvtHandle->weights().size()>=9){
	LHEweight_QCDscale_muR1_muF1     = genHEPMCweight * EvtHandle->weights()[0].wgt / EvtHandle->originalXWGTUP(); // just for verification (should be 1)
	LHEweight_QCDscale_muR1_muF2     = genHEPMCweight * EvtHandle->weights()[1].wgt / EvtHandle->originalXWGTUP();
	LHEweight_QCDscale_muR1_muF0p5   = genHEPMCweight * EvtHandle->weights()[2].wgt / EvtHandle->originalXWGTUP();
	LHEweight_QCDscale_muR2_muF1     = genHEPMCweight * EvtHandle->weights()[3].wgt / EvtHandle->originalXWGTUP();
	LHEweight_QCDscale_muR2_muF2     = genHEPMCweight * EvtHandle->weights()[4].wgt / EvtHandle->originalXWGTUP();
	LHEweight_QCDscale_muR2_muF0p5   = genHEPMCweight * EvtHandle->weights()[5].wgt / EvtHandle->originalXWGTUP();
	LHEweight_QCDscale_muR0p5_muF1   = genHEPMCweight * EvtHandle->weights()[6].wgt / EvtHandle->originalXWGTUP();
	LHEweight_QCDscale_muR0p5_muF2   = genHEPMCweight * EvtHandle->weights()[7].wgt / EvtHandle->originalXWGTUP();
	LHEweight_QCDscale_muR0p5_muF0p5 = genHEPMCweight * EvtHandle->weights()[8].wgt / EvtHandle->originalXWGTUP();
      }
    }
    
    mch.genAcceptance(gen_ZZ2l2qInEtaAcceptance, gen_ZZ2l2qInEtaPtAcceptance);

    ++Nevt_Gen_lumiBlock;
    if (genFinalState == EEEE) {
      ++gen_ZZ4e;
      if (gen_ZZ2l2qInEtaAcceptance) ++gen_ZZ4e_EtaAcceptance;
      if (gen_ZZ2l2qInEtaPtAcceptance) ++gen_ZZ4e_LeptonAcceptance;
    } else if (genFinalState == MMMM) {
      ++gen_ZZ4mu;
      if (gen_ZZ2l2qInEtaAcceptance) ++gen_ZZ4mu_EtaAcceptance;
      if (gen_ZZ2l2qInEtaPtAcceptance) ++gen_ZZ4mu_LeptonAcceptance;
    } else if (genFinalState == EEMM) {
      ++gen_ZZ2mu2e;
      if (gen_ZZ2l2qInEtaAcceptance) ++gen_ZZ2mu2e_EtaAcceptance;
      if (gen_ZZ2l2qInEtaPtAcceptance) ++gen_ZZ2mu2e_LeptonAcceptance;
    } else if (genFinalState == LLTT){
      ++gen_ZZ2l2tau;
    } else if (genFinalState == BUGGY){ // handle H->ddbar 2012 generator bug!!!
      ++gen_BUGGY;
      return; // BUGGY events are skipped
    } else {
      ++gen_Unknown;
    }
    

    //Information on generated candidates, will be used later
    genH = mch.genH();
    genZLeps     = mch.sortedGenZZLeps();
    genAssocLeps = mch.genAssociatedLeps();
    genFSR       = mch.genFSR();
  }
  // End of MC history analysis ------------------------------------------


  // Get candidate collection
  edm::Handle<edm::View<pat::CompositeCandidate> > candHandle;
  if (!event.getByToken(candToken, candHandle)) return ;
  const edm::View<pat::CompositeCandidate>* cands = candHandle.product();
    
  edm::Handle<edm::View<pat::CompositeCandidate> > candHandleFat;
  if (!event.getByToken(candTokenFat, candHandleFat)) return;
  const edm::View<pat::CompositeCandidate>* candsFat = candHandleFat.product();

  if (skipEmptyEvents && cands->size() == 0 && candsFat->size() == 0) return; // Skip events with no candidate, unless skipEmptyEvents = false

  // For Z+L CRs, we want only events with exactly 1 Z+l candidate. FIXME: this has to be reviewed.
  // if (theChannel==CRZJ && cands->size() != 1) return;

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
  if (isMC) {
    if(genH != 0){
      FillHGenInfo(genH->p4(),getHqTWeight(genH->p4().M(),genH->p4().Pt()));
    }
    else if(genZLeps.size()==4){ // for 4l events take the mass of the ZZ(4l) system
      FillHGenInfo((genZLeps.at(0)->p4()+genZLeps.at(1)->p4()+genZLeps.at(2)->p4()+genZLeps.at(3)->p4()),0);
    }

    GenEventInfoProduct  genInfoP = *(genInfo.product());
    // Calculate NNLO QCD K factors for ggZZ
    //KFactorggZZ = 2.;
    //KFactorggZZ = 1.7; // Jamboree
    KFactorggZZ = (float)spkfactor->Eval(GenHMass); // Moriond2016

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

        // Calculate NNLO/NLO QCD K factors for qqZZ
        // bool sameflavor=(genZLeps.at(0)->pdgId()*genZLeps.at(1)->pdgId() == genZLeps.at(2)->pdgId()*genZLeps.at(3)->pdgId());
        KFactorQCDqqZZ_dPhi = 1.;
	// kfactor_qqZZ_qcd_dPhi( fabs(GenZ1Phi-GenZ2Phi), (sameflavor)?1:2 );
        KFactorQCDqqZZ_M    = 1.;
	// kfactor_qqZZ_qcd_M   ( GenHMass, (sameflavor)?1:2 );
        KFactorQCDqqZZ_Pt   = 1.;
	//kfactor_qqZZ_qcd_Pt  ( GenHPt, (sameflavor)?1:2 );

        // Calculate NLO EWK K factors for qqZZ
        TLorentzVector GENZ1Vec,GENZ2Vec;
        GENZ1Vec.SetPtEtaPhiM(GenZ1Pt, GenZ1Eta, GenZ1Phi, GenZ1Mass);
        GENZ2Vec.SetPtEtaPhiM(GenZ2Pt, GenZ2Eta, GenZ2Phi, GenZ2Mass);
        KFactorEWKqqZZ = EwkCorrections::getEwkCorrections(genParticles, ewkTable, genInfoP, GENZ1Vec,GENZ2Vec);
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
  }


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
    } else {
      for (unsigned i=0; i<cleanedJets.size(); ++i) {	  
	if (cleanedJets[i]!=0  &&  (!jetCleaner::isGood(*cand, *(cleanedJets[i])))) {
	  cleanedJets[i]=0;
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
    }

    // count jec up/down njets pt30
    float jec_unc = cleanedJets[i]->userFloat("jec_unc");

    float pt_up = cleanedJets[i]->pt() * (1.0 + jec_unc); 
    float pt_dn = cleanedJets[i]->pt() * (1.0 - jec_unc);
 
    if (pt_up>30) ++nCleanedJetsPt30_jecUp;
    if (pt_dn>30) ++nCleanedJetsPt30_jecDn;

    
    if (writeJets && theChannel!=ZL) FillJet(*(cleanedJets.at(i))); // No additional pT cut (for JEC studies)
  }

  // Now we can write the variables for candidates (fat first)
  int nFilled=0;
  if (theChannel!=CRZJ) {
    for( edm::View<pat::CompositeCandidate>::const_iterator cand = candsFat->begin(); cand != candsFat->end(); ++cand) {
      
      if (!( theChannel==ZL || (bool)(cand->userFloat("isBestCand")) )) continue; // Skip events other than the best cand (or CR candidates in the CR)
      
      FillCandidate(*cand, evtPassTrigger&&evtPassSkim, event, -1, true);
      
      // Fill the candidate as one entry in the tree. Do not reinitialize the event variables, as in CRs
      // there could be several candidates per event.
      // myTree->FillCurrentTree();
      ++nFilled;
    } 
  }
  for( edm::View<pat::CompositeCandidate>::const_iterator cand = cands->begin(); cand != cands->end(); ++cand) {
    size_t icand= cand-cands->begin();
    
    if (!( theChannel==ZL || CRFLAG[icand] || (bool)(cand->userFloat("isBestCand")) )) continue; // Skip events other than the best cand (or CR candidates in the CR)
    
    //For the SR, also fold information about acceptance in CRflag. 
    if (isMC && (theChannel==ZZ)) {
      if (gen_ZZ2l2qInEtaAcceptance)   set_bit(CRFLAG[icand],28);
      if (gen_ZZ2l2qInEtaPtAcceptance) set_bit(CRFLAG[icand],29);
    }
    FillCandidate(*cand, evtPassTrigger&&evtPassSkim, event, CRFLAG[icand], false);
    
    // Fill the candidate as one entry in the tree. Do not reinitialize the event variables, as in CRs
    // there could be several candidates per event.
    // myTree->FillCurrentTree();
    ++nFilled;
  }

  // If no candidate was filled but we still want to keep gen-level and weights, we need to fill one entry anyhow.
  if (skipEmptyEvents==false || nFilled>0) myTree->FillCurrentTree(); 
}


void HZZ2l2qNtupleMaker::FillJet(const pat::Jet& jet)
{
   JetPt  .push_back( jet.pt());
   JetEta .push_back( jet.eta());
   JetPhi .push_back( jet.phi());
   JetMass .push_back( jet.p4().M());
   JetBTagger .push_back( jet.userFloat("bTagger"));
   JetIsBtagged .push_back( jet.userFloat("isBtagged"));
   JetQGLikelihood .push_back( jet.userFloat("qgLikelihood"));
   JetSigma .push_back(jet.userFloat("jec_unc"));
   JetIsInZZCand.push_back(false);
}


void HZZ2l2qNtupleMaker::FillCandidate(const pat::CompositeCandidate& cand, bool evtPass, const edm::Event& event, Int_t CRFLAG, bool isMerged)
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
       
  CRflag.push_back((short)CRFLAG);

  if(theChannel!=ZL && theChannel!=CRZJ){

    //Fill the info on the Higgs candidate
    ZZMass.push_back(cand.p4().mass());
    ZZMass_JecUp.push_back(cand.userFloat("ZZMass_JecUp"));
    ZZMass_JecDown.push_back(cand.userFloat("ZZMass_JecDown"));
    // ZZMassErr.push_back(cand.userFloat("massError"));
    // ZZMassErrCorr.push_back(cand.userFloat("massErrorCorr"));
    ZZMassPreFSR.push_back(cand.userFloat("m4l"));

    ZZPt.push_back(cand.p4().pt());
    ZZEta.push_back(cand.p4().eta());
    ZZPhi.push_back(cand.p4().phi());

    if(addKinRefit){
      /*  if (isMerged) {
	ZZMassRefit.push_back(-1.);
	Z1MassRefitpush_back(-1.);
	ZZMassUnrefitErr.push_back(-1.);
	} else { */
      ZZMassRefit.push_back(cand.userFloat("ZZMassRefit"));
      Z1MassRefit.push_back(cand.userFloat("Z1MassRefit"));
      ZZMassUnrefitErr.push_back(cand.userFloat("ZZMassUnrefitErr"));
      // }
    }
    if(addVtxFit && !isMerged){
      ZZMassCFit.push_back(cand.userFloat("CFitM"));
      ZZChi2CFit.push_back(cand.userFloat("CFitChi2"));
    }

    DiJetMass = cand.userFloat("DiJetMass");
    DiJetDEta = cand.userFloat("DiJetDEta");
    DiJetFisher  = cand.userFloat("DiJetFisher");
    //    DiJetMassPlus  = cand.userFloat("DiJetMassPlus");
    //    DiJetMassMinus  = cand.userFloat("DiJetMassMinus");

    //Fill the angular variables
    helcosthetaZ1.push_back(cand.userFloat("costheta1"));
    helcosthetaZ2.push_back(cand.userFloat("costheta2"));
    helphi      .push_back(cand.userFloat("phi"));
    costhetastar.push_back(cand.userFloat("costhetastar"));
    phistarZ1     .push_back(cand.userFloat("phistar1"));
    //phistarZ2     .push_back(cand.userFloat("phistar2"));
    xi           .push_back(cand.userFloat("xi"));
    xistar       .push_back(cand.userFloat("xistar"));
  
    // JHUGen MEs
    p0plus_VAJHU.push_back(cand.userFloat("p0plus_VAJHU"));
    p0minus_VAJHU.push_back(cand.userFloat("p0minus_VAJHU"));
    p0hplus_VAJHU.push_back(cand.userFloat("p0hplus_VAJHU"));
    p2bplus_VAJHU.push_back(cand.userFloat("p2bplus_VAJHU"));
    p2_VAJHU.push_back(cand.userFloat("p2_VAJHU"));
    pvbf_VAJHU_highestPTJets.push_back(cand.userFloat("pvbf_VAJHU_highestPTJets"));
    phjj_VAJHU_highestPTJets.push_back(cand.userFloat("phjj_VAJHU_highestPTJets"));
    pvbf_VAJHU_highestPTJets_up.push_back(cand.userFloat("pvbf_VAJHU_highestPTJets_up"));
    phjj_VAJHU_highestPTJets_up.push_back(cand.userFloat("phjj_VAJHU_highestPTJets_up"));
    pvbf_VAJHU_highestPTJets_dn.push_back(cand.userFloat("pvbf_VAJHU_highestPTJets_dn"));
    phjj_VAJHU_highestPTJets_dn.push_back(cand.userFloat("phjj_VAJHU_highestPTJets_dn"));

    // MCFM MEs
    pqqZJJ_VAMCFM.push_back(cand.userFloat("pqqZJJ_VAMCFM"));
    bkg_VAMCFM.push_back(cand.userFloat("bkg_VAMCFM"));
    p0plus_VAMCFM.push_back(cand.userFloat("p0plus_VAMCFM"));
    ggzz_VAMCFM.push_back(cand.userFloat("ggzz_VAMCFM"));
    ggzz_p0plus_VAMCFM.push_back(cand.userFloat("ggzz_p0plus_VAMCFM"));
    Dgg10_VAMCFM.push_back(cand.userFloat("Dgg10_VAMCFM"));

    //Float_t ZZLD = cand.userFloat("LD");
    //Float_t ZZLDPSig = cand.userFloat("PSig");
    //Float_t ZZLDPBkg = cand.userFloat("PBkg");
    //Float_t ZZpseudoLD = cand.userFloat("pseudoLD");
    //Float_t ZZgravLD = cand.userFloat("spin2PMLD");
    //Float_t ZZMEKDLD = cand.userFloat("MEKD_LD");
    //Float_t ZZMEKDpseudoLD = cand.userFloat("MEKD_PseudoLD");
    //Float_t ZZMEKDgravLD = cand.userFloat("MEKD_GravLD");

    //Float_t p0plus_melaNorm = cand.userFloat("p0plus_melaNorm");
    //Float_t p0plus_mela = cand.userFloat("p0plus_mela");
    //Float_t p0minus_mela = cand.userFloat("p0minus_mela");
    //Float_t p0hplus_mela = cand.userFloat("p0hplus_mela"); // 0h+, analytic distribution
/* 
    p0plus_VAJHU.push_back(cand.userFloat("p0plus_VAJHU"));
    p0minus_VAJHU.push_back(cand.userFloat("p0minus_VAJHU"));
    p0plus_VAMCFM.push_back(cand.userFloat("p0plus_VAMCFM"));
    p0hplus_VAJHU.push_back(cand.userFloat("p0hplus_VAJHU")); // 0h+ (high dimensional operator), vector algebra, JHUgen
*/
    //Float_t p1_mela = cand.userFloat("p1_mela");
    //Float_t p1_prodIndep_mela = cand.userFloat("p1_prodIndep_mela");
    //Float_t p1plus_mela = cand.userFloat("p1plus_mela"); // 1+, analytic distribution 
    //Float_t p1plus_prodIndep_mela = cand.userFloat("p1plus_prodIndep_mela"); // 1+, analytic distribution 
/*
    p1_VAJHU.push_back(cand.userFloat("p1_VAJHU"));
    p1_prodIndep_VAJHU.push_back(cand.userFloat("p1_prodIndep_VAJHU"));
    p1plus_VAJHU.push_back(cand.userFloat("p1plus_VAJHU")); // 1+ (axial vector), vector algebra, JHUgen,
    p1plus_prodIndep_VAJHU.push_back(cand.userFloat("p1plus_prodIndep_VAJHU")); // 1+ (axial vector), vector algebra, JHUgen,
*/
    //Float_t p2_mela  = cand.userFloat("p2_mela");
    //Float_t p2_prodIndep_mela  = cand.userFloat("p2_prodIndep_mela");
    //Float_t p2qqb_mela = cand.userFloat("p2qqb_mela"); // graviton produced by qqbar vector algebra, analytical,
    //Float_t p2hplus_mela = cand.userFloat("p2hplus_mela"); // graviton produced by qqbar vector algebra, analytical,
    //Float_t p2hminus_mela = cand.userFloat("p2hminus_mela"); // graviton produced by qqbar vector algebra, analytical,
    //Float_t p2bplus_mela = cand.userFloat("p2bplus_mela"); // graviton produced by qqbar vector algebra, analytical,
/*
    p2_VAJHU.push_back(cand.userFloat("p2_VAJHU"));
    p2_prodIndep_VAJHU.push_back(cand.userFloat("p2_prodIndep_VAJHU"));
    p2qqb_VAJHU.push_back(cand.userFloat("p2qqb_VAJHU"));
    p2hplus_VAJHU.push_back(cand.userFloat("p2hplus_VAJHU"));
    p2hminus_VAJHU.push_back(cand.userFloat("p2hminus_VAJHU"));
    p2bplus_VAJHU.push_back(cand.userFloat("p2bplus_VAJHU"));
    p2hplus_qqb_VAJHU.push_back(cand.userFloat(           "p2hplus_qqb_VAJHU"));
    p2hplus_prodIndep_VAJHU.push_back(cand.userFloat(     "p2hplus_prodIndep_VAJHU"));
    p2hminus_qqb_VAJHU.push_back(cand.userFloat(          "p2hminus_qqb_VAJHU"));
    p2hminus_prodIndep_VAJHU.push_back(cand.userFloat(    "p2hminus_prodIndep_VAJHU"));
    p2bplus_qqb_VAJHU.push_back(cand.userFloat(           "p2bplus_qqb_VAJHU"));
    p2bplus_prodIndep_VAJHU.push_back(cand.userFloat(     "p2bplus_prodIndep_VAJHU"));
    p2h2plus_gg_VAJHU.push_back(cand.userFloat(           "p2h2plus_gg_VAJHU"));
    p2h2plus_qqbar_VAJHU.push_back(cand.userFloat(        "p2h2plus_qqbar_VAJHU"));
    p2h2plus_prodIndep_VAJHU.push_back(cand.userFloat(    "p2h2plus_prodIndep_VAJHU"));
    p2h3plus_gg_VAJHU.push_back(cand.userFloat(           "p2h3plus_gg_VAJHU"       ));
    p2h3plus_qqbar_VAJHU.push_back(cand.userFloat(        "p2h3plus_qqbar_VAJHU"    ));
    p2h3plus_prodIndep_VAJHU.push_back(cand.userFloat(    "p2h3plus_prodIndep_VAJHU"));
    p2h6plus_gg_VAJHU.push_back(cand.userFloat(           "p2h6plus_gg_VAJHU"       ));
    p2h6plus_qqbar_VAJHU.push_back(cand.userFloat(        "p2h6plus_qqbar_VAJHU"    ));
    p2h6plus_prodIndep_VAJHU.push_back(cand.userFloat(    "p2h6plus_prodIndep_VAJHU"));
    p2h7plus_gg_VAJHU.push_back(cand.userFloat(           "p2h7plus_gg_VAJHU"       ));
    p2h7plus_qqbar_VAJHU.push_back(cand.userFloat(        "p2h7plus_qqbar_VAJHU"    ));
    p2h7plus_prodIndep_VAJHU.push_back(cand.userFloat(    "p2h7plus_prodIndep_VAJHU"));
    p2h9minus_gg_VAJHU.push_back(cand.userFloat(          "p2h9minus_gg_VAJHU"       ));
    p2h9minus_qqbar_VAJHU.push_back(cand.userFloat(       "p2h9minus_qqbar_VAJHU"    ));
    p2h9minus_prodIndep_VAJHU.push_back(cand.userFloat(   "p2h9minus_prodIndep_VAJHU"));
    p2h10minus_gg_VAJHU.push_back(cand.userFloat(         "p2h10minus_gg_VAJHU"      )); 
    p2h10minus_qqbar_VAJHU.push_back(cand.userFloat(      "p2h10minus_qqbar_VAJHU"   )); 
    p2h10minus_prodIndep_VAJHU.push_back(cand.userFloat(  "p2h10minus_prodIndep_VAJHU")); 
*/
    // bkg_mela.push_back(cand.userFloat("bkg_mela"));
/*
    bkg_VAMCFM.push_back(cand.userFloat("bkg_VAMCFM"));
    bkg_prodIndep_VAMCFM.push_back(cand.userFloat("bkg_prodIndep_VAMCFM"));
    ggzz_VAMCFM.push_back(cand.userFloat("ggzz_VAMCFM"));
    ggzz_p0plus_VAMCFM.push_back(cand.userFloat("ggzz_p0plus_VAMCFM"));
    ggzz_c1_VAMCFM.push_back(cand.userFloat("ggzz_c1_VAMCFM"));
    ggzz_c5_VAMCFM.push_back(cand.userFloat("ggzz_c5_VAMCFM"));
    ggzz_ci_VAMCFM.push_back(cand.userFloat("ggzz_ci_VAMCFM"));
*/
    // bkg_VAMCFMNorm = cand.userFloat("bkg_VAMCFMNorm"));
    // p0_pt = cand.userFloat("p0_pt"));
    // p0_y = cand.userFloat("p0_y"));
    // bkg_pt = cand.userFloat("bkg_pt"));
    // bkg_y = cand.userFloat("bkg_y"));

/*
    pg1g4_mela.push_back(cand.userFloat("pg1g4_mela"));
    pg1g4_VAJHU.push_back(cand.userFloat("pg1g4_VAJHU"));
    pg1g4_pi2_VAJHU.push_back(cand.userFloat("pg1g4_pi2_VAJHU"));
    pg1g2_pi2_VAJHU.push_back(cand.userFloat("pg1g2_pi2_VAJHU"));
    pg1g2_mela.push_back(cand.userFloat("pg1g2_mela"));
    pg1g2_VAJHU.push_back(cand.userFloat("pg1g2_VAJHU"));
    p0_g1prime2_VAJHU.push_back(cand.userFloat("p0_g1prime2_VAJHU"));
    pg1g1prime2_VAJHU.push_back(cand.userFloat("pg1g1prime2_VAJHU"));
    Dgg10_VAMCFM.push_back(cand.userFloat("Dgg10_VAMCFM"));
    pzzzg_VAJHU.push_back(cand.userFloat("pzzzg_VAJHU"));
    pzzgg_VAJHU.push_back(cand.userFloat("pzzgg_VAJHU"));
    pzzzg_PS_VAJHU.push_back(cand.userFloat("pzzzg_PS_VAJHU"));
    pzzgg_PS_VAJHU.push_back(cand.userFloat("pzzgg_PS_VAJHU"));
    p0Zgs_VAJHU.push_back(cand.userFloat("p0Zgs_VAJHU"));
    p0gsgs_VAJHU.push_back(cand.userFloat("p0gsgs_VAJHU"));
    p0Zgs_PS_VAJHU.push_back(cand.userFloat("p0Zgs_PS_VAJHU"));
    p0gsgs_PS_VAJHU.push_back(cand.userFloat("p0gsgs_PS_VAJHU"));
    p0Zgs_g1prime2_VAJHU.push_back(cand.userFloat("p0Zgs_g1prime2_VAJHU"));
    pzzzg_g1prime2_VAJHU.push_back(cand.userFloat("pzzzg_g1prime2_VAJHU"));
    pzzzg_g1prime2_pi2_VAJHU.push_back(cand.userFloat("pzzzg_g1prime2_pi2_VAJHU"));

    p0plus_m2l2q.push_back(cand.userFloat("p0plus_m4l"));
    bkg_m2l2q.push_back(cand.userFloat("bkg_m4l"));
    p0plus_m2l2q_ScaleUp.push_back(cand.userFloat("p0plus_m4l_ScaleUp"));// signal m2l2q probability for systematics
    bkg_m2l2q_ScaleUp.push_back(cand.userFloat("bkg_m4l_ScaleUp"));// backgroun m2l2q probability for systematics
    p0plus_m2l2q_ScaleDown.push_back(cand.userFloat("p0plus_m4l_ScaleDown"));// signal m2l2q probability for systematics
    bkg_m2l2q_ScaleDown.push_back(cand.userFloat("bkg_m4l_ScaleDown"));// backgroun m2l2q probability for systematics
    p0plus_m2l2q_ResUp.push_back(cand.userFloat("p0plus_m4l_ResUp"));// signal m2l2q probability for systematics
    bkg_m2l2q_ResUp.push_back(cand.userFloat("bkg_m4l_ResUp"));// backgroun m2l2q probability for systematics
    p0plus_m2l2q_ResDown.push_back(cand.userFloat("p0plus_m4l_ResDown"));// signal m2l2q probability for systematics
    bkg_m2l2q_ResDown.push_back(cand.userFloat("bkg_m4l_ResDown"));// backgroun m2l2q probability for systematics

    phjj_VAJHU_highestPTJets.push_back(cand.userFloat("phjj_VAJHU_highestPTJets"));
    pvbf_VAJHU_highestPTJets.push_back(cand.userFloat("pvbf_VAJHU_highestPTJets"));
    phjj_VAJHU_highestPTJets_up.push_back(cand.userFloat("phjj_VAJHU_highestPTJets_up"));
    pvbf_VAJHU_highestPTJets_up.push_back(cand.userFloat("pvbf_VAJHU_highestPTJets_up"));
    phjj_VAJHU_highestPTJets_dn.push_back(cand.userFloat("phjj_VAJHU_highestPTJets_dn"));
    pvbf_VAJHU_highestPTJets_dn.push_back(cand.userFloat("pvbf_VAJHU_highestPTJets_dn"));
    phjj_VAJHU_bestDjet.push_back(cand.userFloat("phjj_VAJHU_bestDjet"));
    pvbf_VAJHU_bestDjet.push_back(cand.userFloat("pvbf_VAJHU_bestDjet"));
    phjj_VAJHU_bestDjet_up.push_back(cand.userFloat("phjj_VAJHU_bestDjet_up"));
    pvbf_VAJHU_bestDjet_up.push_back(cand.userFloat("pvbf_VAJHU_bestDjet_up"));
    phjj_VAJHU_bestDjet_dn.push_back(cand.userFloat("phjj_VAJHU_bestDjet_dn"));
    pvbf_VAJHU_bestDjet_dn.push_back(cand.userFloat("pvbf_VAJHU_bestDjet_dn"));

    pAux_vbf_VAJHU.push_back(cand.userFloat("pAux_vbf_VAJHU"));
    pAux_vbf_VAJHU_up.push_back(cand.userFloat("pAux_vbf_VAJHU_up"));
    pAux_vbf_VAJHU_dn.push_back(cand.userFloat("pAux_vbf_VAJHU_dn"));

    phj_VAJHU.push_back(cand.userFloat("phj_VAJHU"));
    phj_VAJHU_up.push_back(cand.userFloat("phj_VAJHU_up"));
    phj_VAJHU_dn.push_back(cand.userFloat("phj_VAJHU_dn"));
    
    pwh_hadronic_VAJHU.push_back(cand.userFloat("pwh_hadronic_VAJHU"));
    pwh_hadronic_VAJHU_up.push_back(cand.userFloat("pwh_hadronic_VAJHU_up"));
    pwh_hadronic_VAJHU_dn.push_back(cand.userFloat("pwh_hadronic_VAJHU_dn"));

    pzh_hadronic_VAJHU.push_back(cand.userFloat("pzh_hadronic_VAJHU"));
    pzh_hadronic_VAJHU_up.push_back(cand.userFloat("pzh_hadronic_VAJHU_up"));
    pzh_hadronic_VAJHU_dn.push_back(cand.userFloat("pzh_hadronic_VAJHU_dn"));

    ptth_VAJHU.push_back(cand.userFloat("ptth_VAJHU"));
    ptth_VAJHU_up.push_back(cand.userFloat("ptth_VAJHU_up"));
    ptth_VAJHU_dn.push_back(cand.userFloat("ptth_VAJHU_dn"));

    pbbh_VAJHU.push_back(cand.userFloat("pbbh_VAJHU"));
    pbbh_VAJHU_up.push_back(cand.userFloat("pbbh_VAJHU_up"));
    pbbh_VAJHU_dn.push_back(cand.userFloat("pbbh_VAJHU_dn"));
*/
  }

  //Z1 and Z2 variables
  const reco::Candidate* Z1;
  const reco::Candidate* Z2;
  vector<const reco::Candidate*> leptons;
  vector<const reco::Candidate*> fsrPhot;
  vector<short> fsrIndex;
  vector<string> labels;

  if (theChannel!=ZL && theChannel!=CRZJ) { // Regular 2l2q candidates
    Z1  = cand.daughter("Z1");
    Z2  = cand.daughter("Z2");
    userdatahelpers::getSortedJetsAndLeptons(cand, leptons, labels, fsrPhot, fsrIndex, true);
  } else {              // Special handling of Z+l candidates
    Z1   = cand.daughter(1); // the Z
    Z2   = cand.daughter(0); // This is actually the additional jet!
    userdatahelpers::getSortedJetsAndLeptons(cand, leptons, labels, fsrPhot, fsrIndex, false); // note: we get just 3 leptons in this case.
  }

  if (isMerged) Z1Mass.push_back(cand.userFloat("d0.ak8PFJetsCHSCorrPrunedMass"));
  else Z1Mass.push_back(Z1->mass());
  Z1Mass_JecUp.push_back(cand.userFloat("Z1Mass_JecUp"));
  Z1Mass_JecDown.push_back(cand.userFloat("Z1Mass_JecDown"));

  float thisZ1tau21 = 0.;
  if (isMerged) thisZ1tau21 = cand.userFloat("d0.NjettinessAK8:tau2")/cand.userFloat("d0.NjettinessAK8:tau1");
  Z1tau21.push_back(thisZ1tau21);

  short int thisZZCandType = 0;
  if (Z1Mass.back() > 70. && Z1Mass.back() < 105.) thisZZCandType = 1;
  else thisZZCandType = -1;
  if (!isMerged) thisZZCandType *= 2;
  ZZCandType.push_back(thisZZCandType);

  Z1Pt.push_back( Z1->pt());
  Z1Flav.push_back(abs(Z1->daughter(0)->pdgId()) * Z1->daughter(0)->charge() * abs(Z1->daughter(1)->pdgId()) * Z1->daughter(1)->charge()); // FIXME: temporarily changed, waiting for a fix to the mismatch of charge() and pdgId() for muons with BTT=4
  
  Z2Mass.push_back(Z2->mass());
  Z2Pt.push_back(  Z2->pt());
  Z2Flav.push_back(theChannel==ZL ? 0 : abs(Z2->daughter(0)->pdgId()) * Z2->daughter(0)->charge() * abs(Z2->daughter(1)->pdgId()) * Z2->daughter(1)->charge()); // FIXME: temporarily changed, waiting for a fix to the mismatch of charge() and pdgId() for muons with BTT=4

  Int_t sel = 0;
  if(theChannel==ZZ){

    // Precomputed selections
    bool candPass70Z2Loose = cand.userFloat("Z2Mass") && 
      cand.userFloat("MAllComb") &&
      cand.userFloat("pt1")>20 && cand.userFloat("pt2")>10. &&
      ZZMass.back()>70.;
    bool candPassFullSel70 = cand.userFloat("SR");
    bool candPassFullSel   = cand.userFloat("FullSel");
    bool candIsBest = cand.userFloat("isBestCand");
    bool passMz_zz = (Z1Mass.back()>60. && Z1Mass.back()<120. && Z2Mass.back()>60. && Z2Mass.back()<120.);   //FIXME hardcoded cut
  
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
  vector<float> SIP(2);
  vector<float> combRelIsoPF(2);
  vector<float> bTagger(2);
  vector<float> qgLik(2);
  vector<float> isBTagged(2);
  vector<float> JecUnc(2);
  passIsoPreFSR = true;
  unsigned int nJets = 2;
  if (theChannel==CRZJ) nJets = 1;
  
  for (unsigned int i=0; i<nJets; ++i){

    JetPt .push_back( leptons[i]->pt() );
    JetEta.push_back( leptons[i]->eta() );
    JetPhi.push_back( leptons[i]->phi() );
    JetMass .push_back( leptons[i]->p4().M());
    JetIsInZZCand.push_back(true);
    
    // std::cout << "PDGId " << i << " : " << leptons[i]->pdgId() << endl; 
    bTagger[i] = -900.;
    isBTagged[i] = -900.;
    qgLik[i] = -900.;
    JecUnc[i] = -900.; 

    if (!isMerged) {
      bTagger[i]           = userdatahelpers::getUserFloat(leptons[i],"bTagger");
      isBTagged[i]           = userdatahelpers::getUserFloat(leptons[i],"isBtagged");
      qgLik[i]           = userdatahelpers::getUserFloat(leptons[i],"qgLikelihood");
      JecUnc[i]           = userdatahelpers::getUserFloat(leptons[i],"jec_unc");
      
      //Fill the info on the jet candidates       
    } else {
      const pat::Jet* jet = dynamic_cast<const pat::Jet*>(leptons[i]);
      bTagger[i] = jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      isBTagged[i] = (float)(jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")>0.460);
      //Fill the info on the subjet candidates   
    }

    JetBTagger .push_back( bTagger[i] );
    JetIsBtagged .push_back( isBTagged[i] );
    JetQGLikelihood .push_back( qgLik[i] );
    JetSigma .push_back( JecUnc[i] );

  }
  for (unsigned int i=nJets; i<leptons.size(); ++i){
    // std::cout << "PDGId " << i << " : " << leptons[i]->pdgId() << endl; 
    short lepFlav = std::abs(leptons[i]->pdgId());

    SIP[i-nJets]           = userdatahelpers::getUserFloat(leptons[i],"SIP");
    passIsoPreFSR      = passIsoPreFSR&&(userdatahelpers::getUserFloat(leptons[i],"combRelIsoPF")<LeptonIsoHelper::isoCut(leptons[i]));

    //in the Legacy approach,  FSR-corrected iso is attached to the Z, not to the lepton!
    if (theChannel!=ZL) {
      combRelIsoPF[i-nJets]    = cand.userFloat(labels[i]+"combRelIsoPFFSRCorr"); // Note: the
      assert(SIP[i-nJets] == cand.userFloat(labels[i]+"SIP")); // Check that I don't mess up with labels[] and leptons[]
    } else {
      //FIXME: at the moment,  FSR-corrected iso is not computed for Z+L CRs
      combRelIsoPF[i-nJets]    = userdatahelpers::getUserFloat(leptons[i],"combRelIsoPF");
    }

    //Fill the info on the lepton candidates  
    LepPt .push_back( leptons[i]->pt() );
    LepEta.push_back( leptons[i]->eta() );
    LepPhi.push_back( leptons[i]->phi() );
    LepLepId.push_back( leptons[i]->pdgId() );
    LepSIP  .push_back( SIP[i-nJets] );
    LepTime .push_back( lepFlav==13 ? userdatahelpers::getUserFloat(leptons[i],"time") : 0. );
    LepisID .push_back( userdatahelpers::getUserFloat(leptons[i],"ID") );
    LepBDT  .push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"BDT") : 0. );
    LepMissingHit.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"missingHit") : 0 );
    //LepChargedHadIso[i].push_back( userdatahelpers::getUserFloat(leptons[i],"PFChargedHadIso") );
    //LepNeutralHadIso[i].push_back( userdatahelpers::getUserFloat(leptons[i],"PFNeutralHadIso") );
    //LepPhotonIso[i].push_back( userdatahelpers::getUserFloat(leptons[i],"PFPhotonIso") );
    LepCombRelIsoPF.push_back( combRelIsoPF[i-nJets] );
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
ZZsel.push_back(sel);

  /* if (theChannel!=ZL) {
    
    //Fill the info on categorization
    nExtraLep = cand.userFloat("nExtraLep"); // Why is this still a float at this point?
    nExtraZ=cand.userFloat("nExtraZ");
    
    //Fill the info on the extra leptons
    for (int iExtraLep=1; iExtraLep<=(int)nExtraLep; iExtraLep++){
      TString extraString;extraString.Form("ExtraLep%d",iExtraLep);
      if(cand.hasUserCand(extraString.Data())){
	//for(int iextra=0;iextra<4;iextra++)varExtra[iextra].Prepend(extraString.Data());
	reco::CandidatePtr candPtr=cand.userCand(extraString.Data());
	ExtraLepPt.push_back(candPtr->pt());
	ExtraLepEta.push_back(candPtr->eta());
	ExtraLepPhi.push_back(candPtr->phi());
	ExtraLepLepId.push_back(candPtr->pdgId());  
      }
    }

  }   */

  //Compute the data/MC weight and overall event weight
  dataMCWeight = 1.;
  if (!isMerged) {
    // for(unsigned int i=0; i<nJets; ++i){
    //  dataMCWeight *= getAllWeightJet(leptons[i]);
    // }
  }
  for(unsigned int i=nJets; i<leptons.size(); ++i){
    dataMCWeight *= getAllWeight(leptons[i]);
  }
  overallEventWeight = PUWeight * genHEPMCweight * dataMCWeight;

}



// ------------ method called once each job just before starting event loop  ------------
void HZZ2l2qNtupleMaker::beginJob()
{
  edm::Service<TFileService> fs;
  myTree = new HZZ4lNtupleFactory( fs->make<TTree>(theFileName,"Event Summary"));
  hCounter = fs->make<TH1F>("Counters", "Counters", 45, 0., 45.);
  BookAllBranches();
}

// ------------ method called once each job just after ending the event loop  ------------
void HZZ2l2qNtupleMaker::endJob()
{
  hCounter->SetBinContent(0 ,gen_sumWeights); // also stored in bin 40
  hCounter->SetBinContent(1 ,Nevt_Gen-gen_BUGGY);
  hCounter->SetBinContent(2 ,gen_ZZ4mu);
  hCounter->SetBinContent(3 ,gen_ZZ4e);
  hCounter->SetBinContent(4 ,gen_ZZ2mu2e);
  hCounter->SetBinContent(5 ,gen_ZZ2l2tau);
  hCounter->SetBinContent(6 ,gen_ZZ4mu_EtaAcceptance);
  hCounter->SetBinContent(7 ,gen_ZZ4mu_LeptonAcceptance);
  hCounter->SetBinContent(10,gen_ZZ4e_EtaAcceptance);
  hCounter->SetBinContent(11,gen_ZZ4e_LeptonAcceptance);
  hCounter->SetBinContent(14,gen_ZZ2mu2e_EtaAcceptance);
  hCounter->SetBinContent(15,gen_ZZ2mu2e_LeptonAcceptance);
  hCounter->SetBinContent(19,gen_BUGGY);
  hCounter->SetBinContent(20,gen_Unknown);

  hCounter->SetBinContent(40,gen_sumWeights); // Also stored in underflow bin; added here for convenience
  hCounter->SetBinContent(41,gen_sumGenMCWeight);
  hCounter->SetBinContent(42,gen_sumPUWeight);


  hCounter->GetXaxis()->SetBinLabel(1 ,"Nevt_Gen");
  hCounter->GetXaxis()->SetBinLabel(2 ,"gen_ZZ4mu");
  hCounter->GetXaxis()->SetBinLabel(3 ,"gen_ZZ4e");
  hCounter->GetXaxis()->SetBinLabel(4 ,"gen_ZZ2mu2e");
  hCounter->GetXaxis()->SetBinLabel(5 ,"gen_ZZ2l2tau");
  hCounter->GetXaxis()->SetBinLabel(6 ,"gen_ZZ4mu_EtaAcceptance");
  hCounter->GetXaxis()->SetBinLabel(7 ,"gen_ZZ4mu_LeptonAcceptance");
  hCounter->GetXaxis()->SetBinLabel(10,"gen_ZZ4e_EtaAcceptance");
  hCounter->GetXaxis()->SetBinLabel(11,"gen_ZZ4e_LeptonAcceptance");
  hCounter->GetXaxis()->SetBinLabel(14,"gen_ZZ2mu2e_EtaAcceptance");
  hCounter->GetXaxis()->SetBinLabel(15,"gen_ZZ2mu2e_LeptonAcceptance");
  hCounter->GetXaxis()->SetBinLabel(19,"gen_BUGGY");
  hCounter->GetXaxis()->SetBinLabel(20,"gen_Unknown");

  hCounter->GetXaxis()->SetBinLabel(40,"gen_sumWeights");
  hCounter->GetXaxis()->SetBinLabel(41,"gen_sumGenMCWeight");
  hCounter->GetXaxis()->SetBinLabel(42,"gen_sumPUWeight");

  return;
}

// ------------ method called when starting to processes a run  ------------
void HZZ2l2qNtupleMaker::beginRun(edm::Run const& iRun, edm::EventSetup const&)
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
void HZZ2l2qNtupleMaker::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void HZZ2l2qNtupleMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
  Nevt_Gen_lumiBlock = 0;
}

// ------------ method called when ending the processing of a luminosity block  ------------
void HZZ2l2qNtupleMaker::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup)
{
  Float_t Nevt_preskim = -1.;
  edm::Handle<edm::MergeableCounter> preSkimCounter;
  if (iLumi.getByToken(preSkimToken, preSkimCounter)) { // Counter before skim. Does not exist for non-skimmed samples.
    Nevt_preskim = preSkimCounter->value;
  }  

  // Nevt_gen: this is the number before any skim
  if (Nevt_preskim>=0.) {
    Nevt_Gen += Nevt_preskim; 
  } else {
    Nevt_Gen += Nevt_Gen_lumiBlock ;
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HZZ2l2qNtupleMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


Float_t HZZ2l2qNtupleMaker::getAllWeight(const reco::Candidate* Lep) const
{
  if (skipDataMCWeight) return 1.;

  Float_t weight  = 1.; 
  //Float_t errCorr = 0.;
  //Float_t errCorrSyst = 0.;

  Float_t myLepPt = Lep->pt();
  Float_t myLepEta = Lep->eta();
  Int_t   myLepID = abs(Lep->pdgId());
  
  //avoid to go out of the TH boundary
  if(myLepID == 13 && myLepPt > 79.) myLepPt = 79.;
  if(myLepID == 11 && myLepPt > 199.) myLepPt = 199.;
  if(myLepID == 11) myLepEta = fabs(myLepEta);

  if(myLepID == 13){  
   
    weight = hTH2D_Mu_All->GetBinContent(hTH2D_Mu_All->GetXaxis()->FindBin(myLepEta),hTH2D_Mu_All->GetYaxis()->FindBin(myLepPt));

  }else if(myLepID == 11){

    if((bool)userdatahelpers::getUserFloat(Lep,"isCrack"))
      weight = hTH2D_El_IdIsoSip_Cracks->GetBinContent(hTH2D_El_IdIsoSip_Cracks->FindFixBin(myLepEta, myLepPt));
    else
      weight = hTH2D_El_IdIsoSip_notCracks->GetBinContent(hTH2D_El_IdIsoSip_notCracks->FindFixBin(myLepEta, myLepPt));

    if(myLepPt < 20.) myLepPt = 20.;
    weight *= hTH2F_El_Reco->GetBinContent(hTH2F_El_Reco->GetXaxis()->FindBin(myLepEta),hTH2F_El_Reco->GetYaxis()->FindBin(myLepPt));

  }else{

    //cout<<"ERROR! wrong lepton ID "<<myLepID<<endl;
    //abort();
    weight = 0.;

  }   

  //FIXME
  if(myLepPt < 5. && myLepID == 13) weight = 1.;

  if(weight < 0.001 || weight > 10.){
    //cout << "ERROR! LEP out of range! myLepPt = " << myLepPt << " myLepEta = " << myLepEta <<" myLepID "<<myLepID<< " weight = " << weight << endl;
    //abort();  //no correction should be zero, if you find one, stop
  }

  return weight;
}

Float_t HZZ2l2qNtupleMaker::getAllWeightJet(const reco::Candidate* Jet) const
{
  if (skipDataMCWeight) return 1.;
  
  float MaxBJetPt = 669.9, MaxLJetPt = 999.9;  // value must be below the boundary
  // bool DoubleUncertainty = false;   
  float JetPt =Jet->pt();
  float weight = 1.;

  cout << "jet " << userdatahelpers::getUserFloat(Jet,"partonFlavor") << " " << JetPt << " " << Jet->eta() << endl;
  
  if (userdatahelpers::getUserFloat(Jet,"partonFlavor") == 5) {
    if (JetPt>MaxBJetPt)  { // use MaxLJetPt for light jets
      JetPt = MaxBJetPt; 
      // DoubleUncertainty = true;
    }  
    weight *= reader->eval(BTagEntry::FLAV_B, Jet->eta(), JetPt); 
  } else if (userdatahelpers::getUserFloat(Jet,"partonFlavor") == 4) {
    if (JetPt>MaxBJetPt)  { // use MaxLJetPt for light jets
      JetPt = MaxBJetPt; 
      // DoubleUncertainty = true;
    }  
    weight *= reader->eval(BTagEntry::FLAV_C, Jet->eta(), JetPt); 
  } else {
    if (JetPt>MaxLJetPt)  { // use MaxLJetPt for light jets
      JetPt = MaxLJetPt; 
      // DoubleUncertainty = true;
    }  
    weight *= reader->eval(BTagEntry::FLAV_UDSG, Jet->eta(), JetPt); 
  } 
  return weight;
}

Float_t HZZ2l2qNtupleMaker::getHqTWeight(double mH, double genPt) const
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
Float_t HZZ2l2qNtupleMaker::getFakeWeight(Float_t LepPt, Float_t LepEta, Int_t LepID, Int_t LepZ1ID)
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
void HZZ2l2qNtupleMaker::FillZGenInfo(Short_t Z1Id, Short_t Z2Id, 
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

void HZZ2l2qNtupleMaker::FillLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, Short_t Lep3Id, Short_t Lep4Id, 
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

  /* GenLep3Pt=Lep3.Pt();
  GenLep3Eta=Lep3.Eta();
  GenLep3Phi=Lep3.Phi();
  GenLep3Id=Lep3Id;

  GenLep4Pt=Lep4.Pt();
  GenLep4Eta=Lep4.Eta();
  GenLep4Phi=Lep4.Phi();
  GenLep4Id=Lep4Id; */
  
  return;
}

void HZZ2l2qNtupleMaker::FillAssocLepGenInfo(std::vector<const reco::Candidate *>& AssocLeps)
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


void HZZ2l2qNtupleMaker::FillHGenInfo(const math::XYZTLorentzVector pH, float w)
{
  GenHMass=pH.M();
  GenHPt=pH.Pt();
  GenHRapidity=pH.Rapidity();

  HqTMCweight=w;

  return;
}


void HZZ2l2qNtupleMaker::BookAllBranches(){
   //Event variables
  myTree->Book("RunNumber",RunNumber);
  myTree->Book("EventNumber",EventNumber);
  myTree->Book("LumiNumber",LumiNumber);
  myTree->Book("NRecoMu",NRecoMu);
  myTree->Book("NRecoEle",NRecoEle);
  myTree->Book("Nvtx",Nvtx);
  myTree->Book("NObsInt",NObsInt);
  myTree->Book("NTrueInt",NTrueInt);
  myTree->Book("PUWeight",PUWeight);
  myTree->Book("KFactorggZZ",KFactorggZZ);
  myTree->Book("KFactorEWKqqZZ",KFactorEWKqqZZ);
  myTree->Book("KFactorQCDqqZZ_dPhi",KFactorQCDqqZZ_dPhi);
  myTree->Book("KFactorQCDqqZZ_M",KFactorQCDqqZZ_M);
  myTree->Book("KFactorQCDqqZZ_Pt",KFactorQCDqqZZ_Pt);
  myTree->Book("PFMET",PFMET);
  myTree->Book("PFMETPhi",PFMETPhi);
  myTree->Book("PFMETNoHF",PFMETNoHF);
  myTree->Book("PFMETNoHFPhi",PFMETNoHFPhi);
  myTree->Book("nCleanedJets",nCleanedJets);
  myTree->Book("nCleanedJetsPt30",nCleanedJetsPt30);
  myTree->Book("nCleanedJetsPt30_jecUp",nCleanedJetsPt30_jecUp);
  myTree->Book("nCleanedJetsPt30_jecDn",nCleanedJetsPt30_jecDn);
  myTree->Book("nCleanedJetsPt30BTagged",nCleanedJetsPt30BTagged);
  myTree->Book("trigWord",trigWord);
  myTree->Book("ZZMass",ZZMass);
  myTree->Book("ZZMass_JecUp",ZZMass_JecUp);
  myTree->Book("ZZMass_JecDown",ZZMass_JecDown);
  myTree->Book("ZZMassErr",ZZMassErr);
  myTree->Book("ZZMassErrCorr",ZZMassErrCorr);
  myTree->Book("ZZMassPreFSR",ZZMassPreFSR);
  myTree->Book("ZZsel",ZZsel);
  myTree->Book("ZZCandType",ZZCandType);
  myTree->Book("ZZPt",ZZPt);
  myTree->Book("ZZEta",ZZEta);
  myTree->Book("ZZPhi",ZZPhi);
  myTree->Book("CRflag",CRflag);
  myTree->Book("Z1Mass",Z1Mass);
  myTree->Book("Z1Mass_JecUp",Z1Mass_JecUp);
  myTree->Book("Z1Mass_JecDown",Z1Mass_JecDown);
  myTree->Book("Z1Pt",Z1Pt);
  myTree->Book("Z1Flav",Z1Flav);
  myTree->Book("Z1tau21",Z1tau21);

  //Kin refitted info
  if (addKinRefit) {
    myTree->Book("ZZMassRefit",ZZMassRefit);
    myTree->Book("Z1MassRefit",Z1MassRefit);
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
  myTree->Book("LepPt",LepPt);
  myTree->Book("LepEta",LepEta);
  myTree->Book("LepPhi",LepPhi);
  myTree->Book("LepLepId",LepLepId);
  myTree->Book("LepSIP",LepSIP);
  myTree->Book("LepTime",LepTime);
  myTree->Book("LepisID",LepisID);
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
  
  //Discriminants

  myTree->Book("p0plus_VAJHU",p0plus_VAJHU);
  myTree->Book("p0minus_VAJHU",p0minus_VAJHU);
  myTree->Book("p0hplus_VAJHU",p0hplus_VAJHU);
  myTree->Book("p2bplus_VAJHU",p2bplus_VAJHU);
  myTree->Book("p2_VAJHU",p2_VAJHU);
  myTree->Book("phjj_VAJHU_highestPTJets", phjj_VAJHU_highestPTJets);
  myTree->Book("pvbf_VAJHU_highestPTJets", pvbf_VAJHU_highestPTJets);
  myTree->Book("phjj_VAJHU_highestPTJets_up", phjj_VAJHU_highestPTJets_up);
  myTree->Book("pvbf_VAJHU_highestPTJets_up", pvbf_VAJHU_highestPTJets_up);
  myTree->Book("phjj_VAJHU_highestPTJets_dn", phjj_VAJHU_highestPTJets_dn);
  myTree->Book("pvbf_VAJHU_highestPTJets_dn", pvbf_VAJHU_highestPTJets_dn);

  myTree->Book("pqqZJJ_VAMCFM", pqqZJJ_VAMCFM);
  myTree->Book("bkg_VAMCFM", bkg_VAMCFM);
  myTree->Book("p0plus_VAMCFM", p0plus_VAMCFM);
  myTree->Book("ggzz_VAMCFM", ggzz_VAMCFM);
  myTree->Book("ggzz_p0plus_VAMCFM", ggzz_p0plus_VAMCFM);
  myTree->Book("Dgg10_VAMCFM", Dgg10_VAMCFM);

  /*
    myTree->Book("p0plus_VAJHU",p0plus_VAJHU);
    myTree->Book("p0minus_VAJHU",p0minus_VAJHU);
    myTree->Book("p0plus_VAMCFM",p0plus_VAMCFM);
    myTree->Book("p0hplus_VAJHU",p0hplus_VAJHU);
    myTree->Book("p1_VAJHU",p1_VAJHU);
    myTree->Book("p1_prodIndep_VAJHU",p1_prodIndep_VAJHU);
    myTree->Book("p1plus_VAJHU",p1plus_VAJHU);
    myTree->Book("p1plus_prodIndep_VAJHU",p1plus_prodIndep_VAJHU);
    myTree->Book("p2_VAJHU",p2_VAJHU);
    myTree->Book("p2_prodIndep_VAJHU",p2_prodIndep_VAJHU);
    myTree->Book("p2qqb_VAJHU",p2qqb_VAJHU);
    myTree->Book("p2hplus_VAJHU",p2hplus_VAJHU);
    myTree->Book("p2hminus_VAJHU",p2hminus_VAJHU);
    myTree->Book("p2bplus_VAJHU",p2bplus_VAJHU);
    myTree->Book("p2hplus_qqb_VAJHU",p2hplus_qqb_VAJHU);
    myTree->Book("p2hplus_prodIndep_VAJHU",p2hplus_prodIndep_VAJHU);
    myTree->Book("p2hminus_qqb_VAJHU",p2hminus_qqb_VAJHU);
    myTree->Book("p2hminus_prodIndep_VAJHU",p2hminus_prodIndep_VAJHU);
    myTree->Book("p2bplus_qqb_VAJHU",p2bplus_qqb_VAJHU);
    myTree->Book("p2bplus_prodIndep_VAJHU",p2bplus_prodIndep_VAJHU);
    myTree->Book("p2h2plus_gg_VAJHU",p2h2plus_gg_VAJHU);
    myTree->Book("p2h2plus_qqbar_VAJHU"               ,p2h2plus_qqbar_VAJHU               );
    myTree->Book("p2h2plus_prodIndep_VAJHU"   ,p2h2plus_prodIndep_VAJHU   );
    myTree->Book("p2h3plus_gg_VAJHU"          ,p2h3plus_gg_VAJHU          );
    myTree->Book("p2h3plus_qqbar_VAJHU"       ,p2h3plus_qqbar_VAJHU       );
    myTree->Book("p2h3plus_prodIndep_VAJHU"   ,p2h3plus_prodIndep_VAJHU   );
    myTree->Book("p2h6plus_gg_VAJHU"          ,p2h6plus_gg_VAJHU          );
    myTree->Book("p2h6plus_qqbar_VAJHU"       ,p2h6plus_qqbar_VAJHU       );
    myTree->Book("p2h6plus_prodIndep_VAJHU"   ,p2h6plus_prodIndep_VAJHU   );
    myTree->Book("p2h7plus_gg_VAJHU"          ,p2h7plus_gg_VAJHU          );
    myTree->Book("p2h7plus_qqbar_VAJHU"       ,p2h7plus_qqbar_VAJHU       );
    myTree->Book("p2h7plus_prodIndep_VAJHU"   ,p2h7plus_prodIndep_VAJHU   );
    myTree->Book("p2h9minus_gg_VAJHU"         ,p2h9minus_gg_VAJHU         );
    myTree->Book("p2h9minus_qqbar_VAJHU"      ,p2h9minus_qqbar_VAJHU      );
    myTree->Book("p2h9minus_prodIndep_VAJHU"  ,p2h9minus_prodIndep_VAJHU  );
    myTree->Book("p2h10minus_gg_VAJHU"       ,p2h10minus_gg_VAJHU       );
    myTree->Book("p2h10minus_qqbar_VAJHU"    ,p2h10minus_qqbar_VAJHU    );
    myTree->Book("p2h10minus_prodIndep_VAJHU",p2h10minus_prodIndep_VAJHU);
    myTree->Book("bkg_VAMCFM",bkg_VAMCFM);
    myTree->Book("bkg_prodIndep_VAMCFM",bkg_prodIndep_VAMCFM);
    myTree->Book("ggzz_VAMCFM",ggzz_VAMCFM);
    myTree->Book("ggzz_p0plus_VAMCFM",ggzz_p0plus_VAMCFM);
    myTree->Book("ggzz_c1_VAMCFM",ggzz_c1_VAMCFM);
    myTree->Book("ggzz_c5_VAMCFM",ggzz_c5_VAMCFM);
    myTree->Book("ggzz_ci_VAMCFM",ggzz_ci_VAMCFM);
    myTree->Book("pg1g4_mela",pg1g4_mela);
    myTree->Book("pg1g4_VAJHU",pg1g4_VAJHU);
    myTree->Book("pg1g4_pi2_VAJHU",pg1g4_pi2_VAJHU);
    myTree->Book("pg1g2_pi2_VAJHU",pg1g2_pi2_VAJHU);
    myTree->Book("pg1g2_mela",pg1g2_mela);
    myTree->Book("pg1g2_VAJHU",pg1g2_VAJHU);
    myTree->Book("p0_g1prime2_VAJHU",p0_g1prime2_VAJHU);
    myTree->Book("pg1g1prime2_VAJHU",pg1g1prime2_VAJHU);
    myTree->Book("Dgg10_VAMCFM",Dgg10_VAMCFM);
    myTree->Book("pzzzg_VAJHU",pzzzg_VAJHU);
    myTree->Book("pzzgg_VAJHU",pzzgg_VAJHU);
    myTree->Book("pzzzg_PS_VAJHU",pzzzg_PS_VAJHU);
    myTree->Book("pzzgg_PS_VAJHU",pzzgg_PS_VAJHU);
    myTree->Book("p0Zgs_VAJHU",p0Zgs_VAJHU);
    myTree->Book("p0gsgs_VAJHU",p0gsgs_VAJHU);
    myTree->Book("p0Zgs_PS_VAJHU",p0Zgs_PS_VAJHU);
    myTree->Book("p0gsgs_PS_VAJHU",p0gsgs_PS_VAJHU);
    myTree->Book("p0Zgs_g1prime2_VAJHU",p0Zgs_g1prime2_VAJHU);
    myTree->Book("pzzzg_g1prime2_VAJHU",pzzzg_g1prime2_VAJHU);
    myTree->Book("pzzzg_g1prime2_pi2_VAJHU",pzzzg_g1prime2_pi2_VAJHU);
    myTree->Book("p0plus_m2l2q",p0plus_m2l2q);
    myTree->Book("bkg_m2l2q",bkg_m2l2q);
    myTree->Book("p0plus_m2l2q_ScaleUp",p0plus_m2l2q_ScaleUp);
    myTree->Book("bkg_m2l2q_ScaleUp",bkg_m2l2q_ScaleUp);
    myTree->Book("p0plus_m2l2q_ScaleDown",p0plus_m2l2q_ScaleDown);
    myTree->Book("bkg_m2l2q_ScaleDown",bkg_m2l2q_ScaleDown);
    myTree->Book("p0plus_m2l2q_ResUp",p0plus_m2l2q_ResUp);
    myTree->Book("bkg_m2l2q_ResUp",bkg_m2l2q_ResUp);
    myTree->Book("p0plus_m2l2q_ResDown",p0plus_m2l2q_ResDown);
    myTree->Book("bkg_m2l2q_ResDown",bkg_m2l2q_ResDown);

  //Production MELA
  myTree->Book("phjj_VAJHU_highestPTJets",phjj_VAJHU_highestPTJets);
  myTree->Book("pvbf_VAJHU_highestPTJets",pvbf_VAJHU_highestPTJets);
  myTree->Book("phjj_VAJHU_highestPTJets_up",phjj_VAJHU_highestPTJets_up);
  myTree->Book("pvbf_VAJHU_highestPTJets_up",pvbf_VAJHU_highestPTJets_up);
  myTree->Book("phjj_VAJHU_highestPTJets_dn",phjj_VAJHU_highestPTJets_dn);
  myTree->Book("pvbf_VAJHU_highestPTJets_dn",pvbf_VAJHU_highestPTJets_dn);
  myTree->Book("phjj_VAJHU_bestDjet",phjj_VAJHU_bestDjet);
  myTree->Book("pvbf_VAJHU_bestDjet",pvbf_VAJHU_bestDjet);
  myTree->Book("phjj_VAJHU_bestDjet_up",phjj_VAJHU_bestDjet_up);
  myTree->Book("pvbf_VAJHU_bestDjet_up",pvbf_VAJHU_bestDjet_up);
  myTree->Book("phjj_VAJHU_bestDjet_dn",phjj_VAJHU_bestDjet_dn);
  myTree->Book("pvbf_VAJHU_bestDjet_dn",pvbf_VAJHU_bestDjet_dn);
  myTree->Book("pAux_vbf_VAJHU",pAux_vbf_VAJHU);
  myTree->Book("pAux_vbf_VAJHU_up",pAux_vbf_VAJHU_up);
  myTree->Book("pAux_vbf_VAJHU_dn",pAux_vbf_VAJHU_dn);
  myTree->Book("phj_VAJHU",phj_VAJHU);
  myTree->Book("phj_VAJHU_up",phj_VAJHU_up);
  myTree->Book("phj_VAJHU_dn",phj_VAJHU_dn);
  myTree->Book("pwh_hadronic_VAJHU",pwh_hadronic_VAJHU);
  myTree->Book("pwh_hadronic_VAJHU_up",pwh_hadronic_VAJHU_up);
  myTree->Book("pwh_hadronic_VAJHU_dn",pwh_hadronic_VAJHU_dn);
  myTree->Book("pzh_hadronic_VAJHU",pzh_hadronic_VAJHU);
  myTree->Book("pzh_hadronic_VAJHU_up",pzh_hadronic_VAJHU_up);
  myTree->Book("pzh_hadronic_VAJHU_dn",pzh_hadronic_VAJHU_dn);
  myTree->Book("ptth_VAJHU",ptth_VAJHU);
  myTree->Book("ptth_VAJHU_up",ptth_VAJHU_up);
  myTree->Book("ptth_VAJHU_dn",ptth_VAJHU_dn);
  myTree->Book("pbbh_VAJHU",pbbh_VAJHU);
  myTree->Book("pbbh_VAJHU_up",pbbh_VAJHU_up);
  myTree->Book("pbbh_VAJHU_dn",pbbh_VAJHU_dn);
 */ 
  //Jet variables
  myTree->Book("JetPt",JetPt);
  myTree->Book("JetEta",JetEta);
  myTree->Book("JetPhi",JetPhi);
  myTree->Book("JetMass",JetMass);
  myTree->Book("JetBTagger",JetBTagger);
  myTree->Book("JetIsBtagged",JetIsBtagged);
  myTree->Book("JetQGLikelihood",JetQGLikelihood);
  myTree->Book("JetSigma",JetSigma);
  myTree->Book("JetIsInZZCand",JetIsInZZCand);
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
  myTree->Book("genFinalState",genFinalState);
  myTree->Book("genProcessId",genProcessId);
  myTree->Book("genHEPMCweight",genHEPMCweight);
  myTree->Book("LHEweight_QCDscale_muR1_muF1",LHEweight_QCDscale_muR1_muF1);
  myTree->Book("LHEweight_QCDscale_muR1_muF2",LHEweight_QCDscale_muR1_muF2);
  myTree->Book("LHEweight_QCDscale_muR1_muF0p5",LHEweight_QCDscale_muR1_muF0p5);
  myTree->Book("LHEweight_QCDscale_muR2_muF1",LHEweight_QCDscale_muR2_muF1);
  myTree->Book("LHEweight_QCDscale_muR2_muF2",LHEweight_QCDscale_muR2_muF2);
  myTree->Book("LHEweight_QCDscale_muR2_muF0p5",LHEweight_QCDscale_muR2_muF0p5);
  myTree->Book("LHEweight_QCDscale_muR0p5_muF1",LHEweight_QCDscale_muR0p5_muF1);
  myTree->Book("LHEweight_QCDscale_muR0p5_muF2",LHEweight_QCDscale_muR0p5_muF2);
  myTree->Book("LHEweight_QCDscale_muR0p5_muF0p5",LHEweight_QCDscale_muR0p5_muF0p5);
  myTree->Book("genExtInfo",genExtInfo);
  myTree->Book("xsec",xsection);
  myTree->Book("dataMCWeight",dataMCWeight);
  myTree->Book("HqTMCweight",HqTMCweight);
  myTree->Book("ZXFakeweight",ZXFakeweight);
  myTree->Book("overallEventWeight",overallEventWeight);
  myTree->Book("GenHMass",GenHMass);
  myTree->Book("GenHPt",GenHPt);
  myTree->Book("GenHRapidity",GenHRapidity);
  myTree->Book("GenZ1Mass",GenZ1Mass);
  myTree->Book("GenZ1Pt",GenZ1Pt);
  myTree->Book("GenZ1Phi",GenZ1Phi);
  myTree->Book("GenZ1Flav",GenZ1Flav);
  myTree->Book("GenZ2Mass",GenZ2Mass);
  myTree->Book("GenZ2Pt",GenZ2Pt);
  myTree->Book("GenZ2Phi",GenZ2Phi);
  myTree->Book("GenZ2Flav",GenZ2Flav);
  myTree->Book("GenLep1Pt",GenLep1Pt);
  myTree->Book("GenLep1Eta",GenLep1Eta);
  myTree->Book("GenLep1Phi",GenLep1Phi);
  myTree->Book("GenLep1Id",GenLep1Id);
  myTree->Book("GenLep2Pt",GenLep2Pt);
  myTree->Book("GenLep2Eta",GenLep2Eta);
  myTree->Book("GenLep2Phi",GenLep2Phi);
  myTree->Book("GenLep2Id",GenLep2Id);
  /* myTree->Book("GenLep3Pt",GenLep3Pt);
  myTree->Book("GenLep3Eta",GenLep3Eta);
  myTree->Book("GenLep3Phi",GenLep3Phi);
  myTree->Book("GenLep3Id",GenLep3Id);
  myTree->Book("GenLep4Pt",GenLep4Pt);
  myTree->Book("GenLep4Eta",GenLep4Eta);
  myTree->Book("GenLep4Phi",GenLep4Phi);
  myTree->Book("GenLep4Id",GenLep4Id); */
  myTree->Book("GenAssocLep1Pt",GenAssocLep1Pt);
  myTree->Book("GenAssocLep1Eta",GenAssocLep1Eta);
  myTree->Book("GenAssocLep1Phi",GenAssocLep1Phi);
  myTree->Book("GenAssocLep1Id",GenAssocLep1Id);
  myTree->Book("GenAssocLep2Pt",GenAssocLep2Pt);
  myTree->Book("GenAssocLep2Eta",GenAssocLep2Eta);
  myTree->Book("GenAssocLep2Phi",GenAssocLep2Phi);
  myTree->Book("GenAssocLep2Id",GenAssocLep2Id);

}

//define this as a plug-in
DEFINE_FWK_MODULE(HZZ2l2qNtupleMaker);
