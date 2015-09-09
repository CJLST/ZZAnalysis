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
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <DataFormats/JetReco/interface/PFJet.h>
#include <DataFormats/JetReco/interface/PFJetCollection.h>
#include <DataFormats/Math/interface/LorentzVector.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <DataFormats/Common/interface/MergeableCounter.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>
#include <ZZAnalysis/AnalysisStep/interface/MCHistoryTools.h>
#include <ZZAnalysis/AnalysisStep/interface/PUReweight.h>
#include <ZZAnalysis/AnalysisStep/interface/bitops.h>
#include <ZZAnalysis/AnalysisStep/interface/Fisher.h>
#include "ZZ4lConfigHelper.h"
#include "HZZ4lNtupleFactory.h"

#include <TRandom3.h>
#include <TH2D.h>
#include "TLorentzVector.h"

#include "ZZAnalysis/AnalysisStep/interface/PUReweight.h"

#include <string>

namespace {
  bool writePhotons = false;  // Write photons in the tree. Note: must be set also in HZZ4lNtupleFactory.cc
  bool writeJets = true;     // Write jets in the tree. FIXME: make this configurable
  bool addKinRefit = false;
  bool addVtxFit = false;
  bool skipDataMCWeight = true; // skip computation of data/MC weight 
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
  Float_t PFMET  =  -99;
  Float_t PFMETPhi  =  -99;
  Short_t nCleanedJets  =  0;
  Short_t nCleanedJetsPt30  = 0;
  Short_t nCleanedJetsPt30BTagged  = 0;
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
  Float_t ZZChi2KinFit  = 0;
  Float_t Z1MassRefit  = 0;
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
  Float_t Lep1Pt  = 0;
  Float_t Lep1Eta  = 0;
  Float_t Lep1Phi  = 0;
  Short_t Lep1LepId  = 0;
  Float_t Lep1SIP  = 0;
  Bool_t Lep1isID  = 0;
  Float_t Lep1BDT  = 0;
  char Lep1missingHit  = 0;
  Float_t Lep1combRelIsoPF  = 0;
  Float_t Lep2Pt  = 0;
  Float_t Lep2Eta  = 0;
  Float_t Lep2Phi  = 0;
  Short_t Lep2LepId  = 0;
  Float_t Lep2SIP  = 0;
  Bool_t Lep2isID  = 0;
  Float_t Lep2BDT  = 0;
  char Lep2missingHit  = 0;
  Float_t Lep2combRelIsoPF  = 0;
  Float_t Lep3Pt  = 0;
  Float_t Lep3Eta  = 0;
  Float_t Lep3Phi  = 0;
  Short_t Lep3LepId  = 0;
  Float_t Lep3SIP  = 0;
  Bool_t Lep3isID  = 0;
  Float_t Lep3BDT  = 0;
  char Lep3missingHit  = 0;
  Float_t Lep3combRelIsoPF  = 0;
  Float_t Lep4Pt  = 0;
  Float_t Lep4Eta  = 0;
  Float_t Lep4Phi  = 0;
  Short_t Lep4LepId  = 0;
  Float_t Lep4SIP  = 0;
  Bool_t Lep4isID  = 0;
  Float_t Lep4BDT  = 0;
  std::vector<float> fsrPt; 
  std::vector<float> fsrEta; 
  std::vector<float> fsrPhi;
  std::vector<short> fsrLept;
  std::vector<short> fsrLeptID;
  std::vector<float> fsrGenPt;
  char Lep4missingHit  = 0;
  Float_t Lep4combRelIsoPF  = 0;
/*Float_t Lep1chargedHadIso  = 0;
  Float_t Lep1neutralHadIso  = 0;
  Float_t Lep1photonIso  = 0;
  Float_t Lep2chargedHadIso  = 0;
  Float_t Lep2neutralHadIso  = 0;
  Float_t Lep2photonIso  = 0;
  Float_t Lep3chargedHadIso  = 0;
  Float_t Lep3neutralHadIso  = 0;
  Float_t Lep3photonIso  = 0;
  Float_t Lep4chargedHadIso  = 0;
  Float_t Lep4neutralHadIso  = 0;
  Float_t Lep4photonIso  = 0;*/
  Float_t PhotPt  = 0;
  Float_t PhotEta  = 0;
  Float_t PhotPhi  = 0;
  Float_t p0plus_VAJHU  = 0;
  Float_t p0minus_VAJHU  = 0;
  Float_t p0plus_VAMCFM  = 0;
  Float_t p0hplus_VAJHU  = 0;
  Float_t p1_VAJHU  = 0;
  Float_t p1_prodIndep_VAJHU  = 0;
  Float_t p1plus_VAJHU  = 0;
  Float_t p1plus_prodIndep_VAJHU  = 0;
  Float_t p2_VAJHU  = 0;
  Float_t p2_prodIndep_VAJHU  = 0;
  Float_t p2qqb_VAJHU  = 0;
  Float_t p2hplus_VAJHU  = 0;
  Float_t p2hminus_VAJHU  = 0;
  Float_t p2bplus_VAJHU  = 0;
  Float_t p2hplus_qqb_VAJHU  = 0;
  Float_t p2hplus_prodIndep_VAJHU  = 0;
  Float_t p2hminus_qqb_VAJHU  = 0;
  Float_t p2hminus_prodIndep_VAJHU  = 0;
  Float_t p2bplus_qqb_VAJHU  = 0;
  Float_t p2bplus_prodIndep_VAJHU  = 0;
  Float_t p2h2plus_gg_VAJHU  = 0;
  Float_t p2h2plus_qqbar_VAJHU                 =                0;
  Float_t p2h2plus_prodIndep_VAJHU     =        0;
  Float_t p2h3plus_gg_VAJHU            =        0;
  Float_t p2h3plus_qqbar_VAJHU         =        0;
  Float_t p2h3plus_prodIndep_VAJHU     =        0;
  Float_t p2h6plus_gg_VAJHU            =        0;
  Float_t p2h6plus_qqbar_VAJHU         =        0;
  Float_t p2h6plus_prodIndep_VAJHU     =        0;
  Float_t p2h7plus_gg_VAJHU            =        0;
  Float_t p2h7plus_qqbar_VAJHU         =        0;
  Float_t p2h7plus_prodIndep_VAJHU     =        0;
  Float_t p2h9minus_gg_VAJHU           =                0;
  Float_t p2h9minus_qqbar_VAJHU        =                0;
  Float_t p2h9minus_prodIndep_VAJHU    =                0;
  Float_t p2h10minus_gg_VAJHU         =                 0;
  Float_t p2h10minus_qqbar_VAJHU      =                 0;
  Float_t p2h10minus_prodIndep_VAJHU  =                 0;
  Float_t bkg_VAMCFM  = 0;
  Float_t bkg_prodIndep_VAMCFM  = 0;
  Float_t ggzz_VAMCFM  = 0;
  Float_t ggzz_p0plus_VAMCFM  = 0;
  Float_t ggzz_c1_VAMCFM  = 0;
  Float_t ggzz_c5_VAMCFM  = 0;
  Float_t ggzz_ci_VAMCFM  = 0;
  Float_t pg1g4_mela  = 0;
  Float_t pg1g4_VAJHU  = 0;
  Float_t pg1g4_pi2_VAJHU  = 0;
  Float_t pg1g2_pi2_VAJHU  = 0;
  Float_t pg1g2_mela  = 0;
  Float_t pg1g2_VAJHU  = 0;
  Float_t p0_g1prime2_VAJHU  =  0;
  Float_t pg1g1prime2_VAJHU  =  0;
  Float_t Dgg10_VAMCFM  =  0;
  Float_t pzzzg_VAJHU  =  0;
  Float_t pzzgg_VAJHU  =  0;
  Float_t pzzzg_PS_VAJHU  =  0;
  Float_t pzzgg_PS_VAJHU  =  0;
  Float_t p0Zgs_VAJHU  =  0;
  Float_t p0gsgs_VAJHU  =  0;
  Float_t p0Zgs_PS_VAJHU  =  0;
  Float_t p0gsgs_PS_VAJHU  =  0;
  Float_t p0Zgs_g1prime2_VAJHU  =  0;
  Float_t pzzzg_g1prime2_VAJHU  =  0;
  Float_t pzzzg_g1prime2_pi2_VAJHU  =  0;
  Float_t p0plus_m4l  =  0;
  Float_t bkg_m4l  =  0;
  Float_t p0plus_m4l_ScaleUp  =  0;
  Float_t bkg_m4l_ScaleUp  = 0;
  Float_t p0plus_m4l_ScaleDown  = 0;
  Float_t bkg_m4l_ScaleDown  = 0;
  Float_t p0plus_m4l_ResUp  = 0;
  Float_t bkg_m4l_ResUp  = 0;
  Float_t p0plus_m4l_ResDown  = 0;
  Float_t bkg_m4l_ResDown  = 0;
  Float_t phjj_VAJHU_old  = 0;
  Float_t pvbf_VAJHU_old  = 0;
  Float_t phjj_VAJHU_old_up  = 0;
  Float_t pvbf_VAJHU_old_up  = 0;
  Float_t phjj_VAJHU_old_dn  = 0;
  Float_t pvbf_VAJHU_old_dn  = 0;
  Float_t phjj_VAJHU_new  = 0;
  Float_t pvbf_VAJHU_new  = 0;
  Float_t phjj_VAJHU_new_up  = 0;
  Float_t pvbf_VAJHU_new_up  = 0;
  Float_t phjj_VAJHU_new_dn  = 0;
  Float_t pvbf_VAJHU_new_dn  = 0;
  Float_t pAux_vbf_VAJHU  = 0;
  Float_t pAux_vbf_VAJHU_up  = 0;
  Float_t pAux_vbf_VAJHU_dn  = 0;
  Float_t phj_VAJHU  = 0;
  Float_t phj_VAJHU_up  = 0;
  Float_t phj_VAJHU_dn  = 0;
  Float_t pwh_hadronic_VAJHU  = 0;
  Float_t pwh_hadronic_VAJHU_up  = 0;
  Float_t pwh_hadronic_VAJHU_dn  = 0;
  Float_t pzh_hadronic_VAJHU  = 0;
  Float_t pzh_hadronic_VAJHU_up  = 0;
  Float_t pzh_hadronic_VAJHU_dn  = 0;
  Float_t ptth_VAJHU  = 0;
  Float_t ptth_VAJHU_up  = 0;
  Float_t ptth_VAJHU_dn  = 0;
  Float_t pbbh_VAJHU  = 0;
  Float_t pbbh_VAJHU_up  = 0;
  Float_t pbbh_VAJHU_dn  = 0;
  std::vector<float> JetPt ;
  std::vector<float> JetEta ;
  std::vector<float> JetPhi ;
  std::vector<float> JetMass ;
  std::vector<float> JetBTagger ;
  std::vector<float> JetIsBtagged; 
  std::vector<float> JetQGLikelihood; 
  std::vector<float> JetSigma ;
  Float_t DiJetMass  = -99;
  Float_t DiJetMassPlus  = -99;
  Float_t DiJetMassMinus  = -99;
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
  Short_t genExtInfo  = 0;
  Float_t xsection  = 0;
  Float_t dataMCWeight  = 0;
  Float_t HqTMCweight  = 0;
  Float_t ZXFakeweight  = 0;
  Float_t GenHMass  = 0;
  Float_t GenHPt  = 0;
  Float_t GenHRapidity  = 0;
  Float_t GenZ1Mass  = 0;
  Float_t GenZ1Pt  = 0;
  Float_t GenZ2Mass  = 0;
  Float_t GenZ2Pt  = 0;
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
  virtual void FillCandidate(const pat::CompositeCandidate& higgs, bool evtPass, const edm::Event&, const Int_t CRflag);
  virtual void FillPhoton(const pat::Photon& photon);
  virtual void FillJet(const pat::Jet& jet);
  virtual void endJob() ;

  void FillHGenInfo(const math::XYZTLorentzVector Hp, float w);
  void FillZGenInfo(const math::XYZTLorentzVector pZ1, const math::XYZTLorentzVector pZ2);
  void FillLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, Short_t Lep3Id, Short_t Lep4Id, 
    const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2, const math::XYZTLorentzVector Lep3, 
    const math::XYZTLorentzVector Lep4, float weight);
  void FillAssocLepGenInfo(std::vector<const reco::Candidate *>& AssocLeps);

  Float_t getAllWeight(Float_t LepPt, Float_t LepEta, Int_t LepID) const;
  Float_t getHqTWeight(double mH, double genPt) const;
  Float_t getFakeWeight(Float_t LepPt, Float_t LepEta, Int_t LepID, Int_t LepZ1ID);

  // ----------member data ---------------------------
  ZZ4lConfigHelper myHelper;
  int theChannel;
  std::string theCandLabel;
  TString theFileName;

  HZZ4lNtupleFactory *myTree;
  TH1F *hCounter;

  Bool_t isMC;

  bool applyTrigger;    // Keep only events passing trigger
  bool applySkim;       //   "     "      "     skim
  bool skipEmptyEvents; // Skip events whith no selected candidate (otherwise, gen info is preserved for all events)
  Float_t xsec;
  int year;
  
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

  TH2D *hTH2D_Mu_All;// = (TH2D*)fMuWeight.Get("TH2D_ALL_2011A"); 
  TH2D *hTH2D_El_All;//  = (TH2D*)fElWeight12.Get(eleSFname.Data());
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
  reweight(),
  hTH2D_Mu_All(0),
  hTH2D_El_All(0),
  h_weight(0)
{
  //cout<< "Beginning Constructor\n\n\n" <<endl;
  theCandLabel = pset.getUntrackedParameter<string>("CandCollection"); // Name of input ZZ collection
  theChannel = myHelper.channel(); // Valid options: ZZ, ZLL, ZL 
  theFileName = pset.getUntrackedParameter<string>("fileName"); 
  skipEmptyEvents = pset.getParameter<bool>("skipEmptyEvents"); // Do not store 
  sampleName = pset.getParameter<string>("sampleName");
  xsec = pset.getParameter<double>("xsec");
  year = pset.getParameter<int>("setup");

  if (skipEmptyEvents) {
    applyTrigger=true;
    applySkim=true;
  } else {
    applyTrigger=false;
    applySkim=false;    
  }

  isMC = myHelper.isMC();

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
  
  if (!skipDataMCWeight) {
    //Scale factors for data/MC efficiency
    //FIXME: to adjust for 13 TeV
    float Run2011AFraction=0.465;
    TString yearString;yearString.Form("TH2D_ALL_%d",year);
    if(year==2011){
      TRandom3 randomGenerator(0);
      Float_t whatPeriod = randomGenerator.Uniform();
      if(whatPeriod < Run2011AFraction) yearString.Append("A");
      else yearString.Append("B");
    }
    TString filename;filename.Form("ZZAnalysis/AnalysisStep/test/Macros/scale_factors_muons%d.root",year); 
    if(year==2015){
      filename.Form("ZZAnalysis/AnalysisStep/test/Macros/scale_factors_muons%d.root",2012); //FIXME
      yearString.Form("TH2D_ALL_%d",2012);
    }
    edm::FileInPath fip(filename.Data());
    std::string fipPath=fip.fullPath();
    TFile *fMuWeight = TFile::Open(fipPath.data(),"READ");
    hTH2D_Mu_All = (TH2D*)fMuWeight->Get(yearString.Data())->Clone(); 

    filename.Form("ZZAnalysis/AnalysisStep/test/Macros/scale_factors_ele%d.root",year); 
    if(year==2015)filename.Form("ZZAnalysis/AnalysisStep/test/Macros/scale_factors_ele%d.root",2012);//FIXME  
    edm::FileInPath fipEle(filename.Data());
    fipPath=fipEle.fullPath();
    TFile *fEleWeight = TFile::Open(fipPath.data(),"READ");
    hTH2D_El_All = (TH2D*)fEleWeight->Get("h_electronScaleFactor_RecoIdIsoSip")->Clone();
    fMuWeight->Close();
    fEleWeight->Close();
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
}


// ------------ method called for each event  ------------
void HZZ4lNtupleMaker::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  // Primary vertices
  Handle<vector<reco::Vertex> >  vertexs;
  event.getByLabel("goodPrimaryVertices",vertexs);
  myTree->InitializeVariables();

  //----------------------------------------------------------------------
  // Analyze MC truth; collect MC weights and update counters (this is done for all generated events, 
  // including those that do not pass skim, trigger etc!)

  bool gen_ZZ4lInEtaAcceptance = false;   // All 4 gen leptons in eta acceptance
  bool gen_ZZ4lInEtaPtAcceptance = false; // All 4 gen leptons in eta,pT acceptance

  const reco::Candidate * genH = 0;
  std::vector<const reco::Candidate *> genZLeps;
  std::vector<const reco::Candidate *> genAssocLeps;
  if (isMC) {
    // get PU weights
    Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    event.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);

    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      if(PVI->getBunchCrossing() == 0) { 
        NObsInt  = PVI->getPU_NumInteractions();
        NTrueInt = PVI->getTrueNumInteractions();
        break;
      } 
    }

    int source = myHelper.sampleType();
    int target = myHelper.setup();
    float weight2 = reweight.weight(source,target,NTrueInt);
    PUWeight=weight2;
    MCHistoryTools mch(event, sampleName);
    genFinalState = mch.genFinalState();
    genProcessId = mch.getProcessID();
    genHEPMCweight = mch.gethepMCweight();
    genExtInfo = mch.genAssociatedFS();

    // keep track of sum of weights
    gen_sumPUWeight    += weight2;
    gen_sumGenMCWeight += genHEPMCweight;
    //weight = weight2*genHEPMCweight;
    gen_sumWeights     += weight2*genHEPMCweight;

    mch.genAcceptance(gen_ZZ4lInEtaAcceptance, gen_ZZ4lInEtaPtAcceptance);

    ++Nevt_Gen_lumiBlock;
    if (genFinalState == EEEE) {
      ++gen_ZZ4e;
      if (gen_ZZ4lInEtaAcceptance) ++gen_ZZ4e_EtaAcceptance;
      if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ4e_LeptonAcceptance;
    } else if (genFinalState == MMMM) {
      ++gen_ZZ4mu;
      if (gen_ZZ4lInEtaAcceptance) ++gen_ZZ4mu_EtaAcceptance;
      if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ4mu_LeptonAcceptance;
    } else if (genFinalState == EEMM) {
      ++gen_ZZ2mu2e;
      if (gen_ZZ4lInEtaAcceptance) ++gen_ZZ2mu2e_EtaAcceptance;
      if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ2mu2e_LeptonAcceptance;
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
  event.getByLabel(theCandLabel, candHandle);
  const edm::View<pat::CompositeCandidate>* cands = candHandle.product();

  if (skipEmptyEvents && cands->size() == 0) return; // Skip events with no candidate, unless skipEmptyEvents = false

  // For Z+L CRs, we want only events with exactly 1 Z+l candidate. FIXME: this has to be reviewed.
  if (theChannel==ZL && cands->size() != 1) return;


  // Apply MC filter (skip event)
  if (isMC && !(myHelper.passMCFilter(event))) return;

  // Apply skim
  bool evtPassSkim = myHelper.passSkim(event, trigWord);
  if (applySkim && !evtPassSkim) return;

  // Apply trigger request (skip event)
  bool evtPassTrigger = myHelper.passTrigger(event, trigWord);
  if (applyTrigger && !evtPassTrigger) return;

  //Fill MC truth information
  if (isMC) {
    if(genH != 0){
      FillHGenInfo(genH->p4(),getHqTWeight(genH->p4().M(),genH->p4().Pt()));
    }
    else if(genZLeps.size()==4){ // for 4l events take the mass of the ZZ(4l) system
      FillHGenInfo((genZLeps.at(0)->p4()+genZLeps.at(1)->p4()+genZLeps.at(2)->p4()+genZLeps.at(3)->p4()),0);
    }

    if (genFinalState!=BUGGY) {
    
      float eff_weight =1.;
      for(int nLep=0;nLep<(int)genZLeps.size();nLep++){
        int lepid = genZLeps.at(nLep)->pdgId();
        if(fabs(lepid) ==13 || fabs(lepid) ==11)eff_weight *= getAllWeight(genZLeps.at(nLep)->pt(), genZLeps.at(nLep)->eta(),genZLeps.at(nLep)->pdgId());
      }
      if (genZLeps.size()==4) {
        
        // "generated Zs" defined with standard pairing applied on gen leptons (genZLeps is sorted by MCHistoryTools)
        FillZGenInfo(genZLeps.at(0)->p4()+genZLeps.at(1)->p4(),
                             genZLeps.at(2)->p4()+genZLeps.at(3)->p4());

        // Gen leptons
        FillLepGenInfo(genZLeps.at(0)->pdgId(), genZLeps.at(1)->pdgId(), genZLeps.at(2)->pdgId(), genZLeps.at(3)->pdgId(),
                               genZLeps.at(0)->p4(), genZLeps.at(1)->p4(), genZLeps.at(2)->p4(), genZLeps.at(3)->p4(),eff_weight);      
      }

      if (genZLeps.size()==3) {
        FillLepGenInfo(genZLeps.at(0)->pdgId(), genZLeps.at(1)->pdgId(), genZLeps.at(2)->pdgId(), 0,
		       genZLeps.at(0)->p4(), genZLeps.at(1)->p4(), genZLeps.at(2)->p4(), *(new math::XYZTLorentzVector),eff_weight);
      }
      if (genZLeps.size()==2) {
        FillLepGenInfo(genZLeps.at(0)->pdgId(), genZLeps.at(1)->pdgId(), 0, 0,
		       genZLeps.at(0)->p4(), genZLeps.at(1)->p4(), *(new math::XYZTLorentzVector), *(new math::XYZTLorentzVector),eff_weight);
      }

      if (genAssocLeps.size()==1 || genAssocLeps.size()==2) {
        FillAssocLepGenInfo(genAssocLeps);
      }

    }
  }

  // Photons (store them only for events with at least 1 candidate)
  // FIXME: should rather write used FSR photons, with info on matching lepton.
  if (writePhotons && cands->size()!=0) {
    edm::Handle<vector<pat::Photon> > photons;
    if (event.getByLabel("cmgPhotonSel",photons)) {
      for( vector<pat::Photon>::const_iterator ph = photons->begin(); ph != photons->end(); ++ph) {
        if(ph->pt() > 2. && fabs(ph->eta()) < 2.8) FillPhoton(*ph);
      }
    }
  }


  // Jets
  Handle<edm::View<pat::Jet> > CleanedJets;
  event.getByLabel("cleanJets", CleanedJets);
  vector<const pat::Jet*> cleanedJets;
  vector<const pat::Jet*> cleanedJetsPt30;
  for(edm::View<pat::Jet>::const_iterator jet = CleanedJets->begin(); jet != CleanedJets->end(); ++jet){
    cleanedJets.push_back(&*jet);
    if(jet->pt()>30){
      cleanedJetsPt30.push_back(&*jet);
      if(jet->userFloat("isBtagged")) nCleanedJetsPt30BTagged++;
    }
  }

  if (cleanedJetsPt30.size()>=2) {
    DiJetDEta = cleanedJetsPt30[0]->eta()-cleanedJetsPt30[1]->eta();
    DiJetMass = (cleanedJetsPt30[0]->p4()+cleanedJetsPt30[1]->p4()).M();
    DiJetFisher = fisher(DiJetMass,DiJetDEta);      
  }

  if (writeJets){

    if(theChannel!=ZL){
      for (unsigned int i=0; i<cleanedJets.size(); i++) {
	FillJet(*(cleanedJets.at(i)));
      }

      // Note that jets variables are filled for jets above 20 GeV, to allow JES studies.
      // detajj, Mjj and ZZFisher are filled only for true dijet events (jets above 30 GeV)
      if(cleanedJets.size()>1){ 
	const pat::Jet& myjet1 = *(cleanedJets.at(0));
	const pat::Jet& myjet2 = *(cleanedJets.at(1));
        math::XYZTLorentzVector jet1 = myjet1.p4();
        math::XYZTLorentzVector jet2 = myjet2.p4();
        Float_t jesUnc1 = 0.;//myjet1.uncOnFourVectorScale();
        Float_t jesUnc2 = 0.;//myjet2.uncOnFourVectorScale();
        math::XYZTLorentzVector jetScalePlus1 = jet1*(1+jesUnc1);
        math::XYZTLorentzVector jetScaleMinus1 = jet1*(1-jesUnc1);
        math::XYZTLorentzVector jetScalePlus2 = jet2*(1+jesUnc2);
        math::XYZTLorentzVector jetScaleMinus2 = jet2*(1-jesUnc2);
        DiJetMassPlus = (jetScalePlus1+jetScalePlus2).M();
        DiJetMassMinus = (jetScaleMinus1+jetScaleMinus2).M();
        //myTree->FillDiJetInfo(Mjj,MjjPlus,MjjMinus,detajj,Fisher);
      }
    }
  }
  
  Handle<vector<reco::MET> > pfmetcoll;
  event.getByLabel("slimmedMETs", pfmetcoll);
  if (pfmetcoll.isValid()){
    PFMET = pfmetcoll->front().pt();
    PFMETPhi = pfmetcoll->front().phi();
  }
  xsection=xsec;
  //Save general event info in the tree
  RunNumber=event.id().run();
  LumiNumber=event.luminosityBlock();
  Nvtx=vertexs->size();
  nCleanedJets=cleanedJets.size();
  nCleanedJetsPt30=cleanedJetsPt30.size();
  xsection=xsec;
  EventNumber=event.id().event();


  // number of reconstructed leptons
  edm::Handle<pat::MuonCollection> muonHandle;
  event.getByLabel("slimmedMuons", muonHandle);
  for(unsigned int i = 0; i< muonHandle->size(); ++i){
    const pat::Muon* m = &((*muonHandle)[i]);
    if(m->pt()>5 && m->isPFMuon()) // these cuts are implicit in miniAOD
      NRecoMu++;
  }
  edm::Handle<pat::ElectronCollection> electronHandle;
  event.getByLabel("slimmedElectrons", electronHandle);
  for(unsigned int i = 0; i< electronHandle->size(); ++i){
    const pat::Electron* e = &((*electronHandle)[i]);
    if(e->pt()>5) // this cut is implicit in miniAOD
      NRecoEle++;
  }
  

  //Loop on the candidates
  int nFilled=0;
  for( edm::View<pat::CompositeCandidate>::const_iterator cand = cands->begin(); cand != cands->end(); ++cand) {
    Int_t CRFLAG=0;
    bool candIsBest = cand->userFloat("isBestCand");

    //    int candChannel = cand->userFloat("candChannel"); // This is currently the product of pdgId of leptons (eg 14641, 28561, 20449)
    
    if (theChannel==ZLL) {
      // AA CRs
      if(cand->userFloat("isBestCRZLLss")&&cand->userFloat("CRZLLss"))set_bit(CRFLAG,CRZLLss);      

      // A CRs
      if(cand->userFloat("isBestCRZLLos_2P2F")&&cand->userFloat("CRZLLos_2P2F"))set_bit(CRFLAG,CRZLLos_2P2F);      
      if(cand->userFloat("isBestCRZLLos_3P1F")&&cand->userFloat("CRZLLos_3P1F"))set_bit(CRFLAG,CRZLLos_3P1F);
    }
    
    //    if(theChannel==ZL){} Nothing special in this case
 
    if (!(candIsBest||CRFLAG)) continue; // Skip events other than the best cand (or CR candidates in the CR)
    
    //For the SR, also fold information about acceptance in CRflag 
    if (isMC && (theChannel==EEEE||theChannel==MMMM||theChannel==EEMM)) {
      if (gen_ZZ4lInEtaAcceptance)   set_bit(CRFLAG,28);
      if (gen_ZZ4lInEtaPtAcceptance) set_bit(CRFLAG,29);
    }
    FillCandidate(*cand, evtPassTrigger&&evtPassSkim, event, CRFLAG);


    // Fill the candidate as one entry in the tre. Do not reinitialize the event variables, as in CRs
    // there could be several candidates per event.
    myTree->FillCurrentTree();
    ++nFilled;
  }
  
  // If no candidate was filled but we still want to keep gen-level and weights, we need to fill one entry anyhow.
  if (skipEmptyEvents==false && nFilled==0)myTree->FillCurrentTree(); 
}

void HZZ4lNtupleMaker::FillPhoton(const pat::Photon& photon)
{
   PhotPt  = photon.pt();
   PhotEta = photon.eta();
   PhotPhi = photon.phi();
}


void HZZ4lNtupleMaker::FillJet(const pat::Jet& jet)
{
   JetPt  .push_back( jet.pt());
   JetEta .push_back( jet.eta());
   JetPhi .push_back( jet.phi());
   JetMass .push_back( jet.p4().M());
   JetBTagger .push_back( jet.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"));
   JetIsBtagged .push_back( jet.userFloat("isBtagged"));
   JetQGLikelihood .push_back( jet.userFloat("qgLikelihood"));
   JetSigma .push_back( 0.);//jet.uncOnFourVectorScale());
}


void HZZ4lNtupleMaker::FillCandidate(const pat::CompositeCandidate& cand, bool evtPass, const edm::Event& event, Int_t CRFLAG)
{
  //Initialize a new candidate into the tree
  myTree->createNewCandidate();

  CRflag = CRFLAG;

  //FIXME: Dobbiamo salvare info su qualita' del fit a H?
  //Chi2 e chi2constr

  //Fill the info on the Higgs candidate
   ZZMass = cand.p4().mass();
   ZZMassErr = cand.userFloat("massError");
   ZZMassErrCorr = cand.userFloat("massErrorCorr");
   ZZMassPreFSR = cand.userFloat("m4l");  
  // FIXME: Other variables to be added if addKinRefit or addVtxFit == true:
//   Float_t Z1MassRefit = cand.userFloat("mZ1Ref");
//   Float_t ZZMassRefit = cand.userFloat("m4lRef");
//   Float_t Chi2KinFit = cand.userFloat("chi2Fit");
//   Float_t ZZMassCFit = cand.userFloat("CFitM");
//   Float_t Chi2CFit = cand.userFloat("CFitChi2");
  
   ZZPt  = cand.p4().pt();
   ZZEta = cand.p4().eta();
   ZZPhi = cand.p4().phi();
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
   p0plus_VAJHU = cand.userFloat("p0plus_VAJHU");
   p0minus_VAJHU = cand.userFloat("p0minus_VAJHU");
   p0plus_VAMCFM = cand.userFloat("p0plus_VAMCFM");
   p0hplus_VAJHU = cand.userFloat("p0hplus_VAJHU"); // 0h+ (high dimensional operator), vector algebra, JHUgen
  //Float_t p1_mela = cand.userFloat("p1_mela");
  //Float_t p1_prodIndep_mela = cand.userFloat("p1_prodIndep_mela");
  //Float_t p1plus_mela = cand.userFloat("p1plus_mela"); // 1+, analytic distribution 
  //Float_t p1plus_prodIndep_mela = cand.userFloat("p1plus_prodIndep_mela"); // 1+, analytic distribution 
   p1_VAJHU = cand.userFloat("p1_VAJHU");
   p1_prodIndep_VAJHU = cand.userFloat("p1_prodIndep_VAJHU");
   p1plus_VAJHU = cand.userFloat("p1plus_VAJHU"); // 1+ (axial vector), vector algebra, JHUgen,
   p1plus_prodIndep_VAJHU = cand.userFloat("p1plus_prodIndep_VAJHU"); // 1+ (axial vector), vector algebra, JHUgen,
  //Float_t p2_mela  = cand.userFloat("p2_mela");
  //Float_t p2_prodIndep_mela  = cand.userFloat("p2_prodIndep_mela");
  //Float_t p2qqb_mela = cand.userFloat("p2qqb_mela"); // graviton produced by qqbar vector algebra, analytical,
  //Float_t p2hplus_mela = cand.userFloat("p2hplus_mela"); // graviton produced by qqbar vector algebra, analytical,
  //Float_t p2hminus_mela = cand.userFloat("p2hminus_mela"); // graviton produced by qqbar vector algebra, analytical,
  //Float_t p2bplus_mela = cand.userFloat("p2bplus_mela"); // graviton produced by qqbar vector algebra, analytical,
   p2_VAJHU = cand.userFloat("p2_VAJHU");
   p2_prodIndep_VAJHU = cand.userFloat("p2_prodIndep_VAJHU");
   p2qqb_VAJHU = cand.userFloat("p2qqb_VAJHU");
   p2hplus_VAJHU = cand.userFloat("p2hplus_VAJHU");
   p2hminus_VAJHU = cand.userFloat("p2hminus_VAJHU");
   p2bplus_VAJHU = cand.userFloat("p2bplus_VAJHU");
   p2hplus_qqb_VAJHU= cand.userFloat(           "p2hplus_qqb_VAJHU");
   p2hplus_prodIndep_VAJHU= cand.userFloat(     "p2hplus_prodIndep_VAJHU");
   p2hminus_qqb_VAJHU= cand.userFloat(          "p2hminus_qqb_VAJHU");
   p2hminus_prodIndep_VAJHU= cand.userFloat(    "p2hminus_prodIndep_VAJHU");
   p2bplus_qqb_VAJHU= cand.userFloat(           "p2bplus_qqb_VAJHU");
   p2bplus_prodIndep_VAJHU= cand.userFloat(     "p2bplus_prodIndep_VAJHU");
   p2h2plus_gg_VAJHU= cand.userFloat(           "p2h2plus_gg_VAJHU");
   p2h2plus_qqbar_VAJHU= cand.userFloat(        "p2h2plus_qqbar_VAJHU");
   p2h2plus_prodIndep_VAJHU= cand.userFloat(    "p2h2plus_prodIndep_VAJHU");
   p2h3plus_gg_VAJHU= cand.userFloat(           "p2h3plus_gg_VAJHU"       );
   p2h3plus_qqbar_VAJHU= cand.userFloat(        "p2h3plus_qqbar_VAJHU"    );
   p2h3plus_prodIndep_VAJHU= cand.userFloat(    "p2h3plus_prodIndep_VAJHU");
   p2h6plus_gg_VAJHU= cand.userFloat(           "p2h6plus_gg_VAJHU"       );
   p2h6plus_qqbar_VAJHU= cand.userFloat(        "p2h6plus_qqbar_VAJHU"    );
   p2h6plus_prodIndep_VAJHU= cand.userFloat(    "p2h6plus_prodIndep_VAJHU");
   p2h7plus_gg_VAJHU= cand.userFloat(           "p2h7plus_gg_VAJHU"       );
   p2h7plus_qqbar_VAJHU= cand.userFloat(        "p2h7plus_qqbar_VAJHU"    );
   p2h7plus_prodIndep_VAJHU= cand.userFloat(    "p2h7plus_prodIndep_VAJHU");
   p2h9minus_gg_VAJHU= cand.userFloat(          "p2h9minus_gg_VAJHU"       );
   p2h9minus_qqbar_VAJHU= cand.userFloat(       "p2h9minus_qqbar_VAJHU"    );
   p2h9minus_prodIndep_VAJHU= cand.userFloat(   "p2h9minus_prodIndep_VAJHU");
   p2h10minus_gg_VAJHU= cand.userFloat(         "p2h10minus_gg_VAJHU"      ); 
   p2h10minus_qqbar_VAJHU= cand.userFloat(      "p2h10minus_qqbar_VAJHU"   ); 
   p2h10minus_prodIndep_VAJHU= cand.userFloat(  "p2h10minus_prodIndep_VAJHU"); 
// bkg_mela = cand.userFloat("bkg_mela");
   bkg_VAMCFM = cand.userFloat("bkg_VAMCFM");
   bkg_prodIndep_VAMCFM = cand.userFloat("bkg_prodIndep_VAMCFM");
   ggzz_VAMCFM = cand.userFloat("ggzz_VAMCFM");
   ggzz_p0plus_VAMCFM = cand.userFloat("ggzz_p0plus_VAMCFM");
   ggzz_c1_VAMCFM = cand.userFloat("ggzz_c1_VAMCFM");
   ggzz_c5_VAMCFM = cand.userFloat("ggzz_c5_VAMCFM");
   ggzz_ci_VAMCFM = cand.userFloat("ggzz_ci_VAMCFM");
  // bkg_VAMCFMNorm = cand.userFloat("bkg_VAMCFMNorm");
  // p0_pt = cand.userFloat("p0_pt");
  // p0_y = cand.userFloat("p0_y");
  // bkg_pt = cand.userFloat("bkg_pt");
  // bkg_y = cand.userFloat("bkg_y");

   pg1g4_mela = cand.userFloat("pg1g4_mela");
   pg1g4_VAJHU = cand.userFloat("pg1g4_VAJHU");
   pg1g4_pi2_VAJHU = cand.userFloat("pg1g4_pi2_VAJHU");
   pg1g2_pi2_VAJHU = cand.userFloat("pg1g2_pi2_VAJHU");
   pg1g2_mela = cand.userFloat("pg1g2_mela");
   pg1g2_VAJHU = cand.userFloat("pg1g2_VAJHU");
   p0_g1prime2_VAJHU = cand.userFloat("p0_g1prime2_VAJHU");
   pg1g1prime2_VAJHU = cand.userFloat("pg1g1prime2_VAJHU");
   Dgg10_VAMCFM = cand.userFloat("Dgg10_VAMCFM");
   pzzzg_VAJHU = cand.userFloat("pzzzg_VAJHU");
   pzzgg_VAJHU = cand.userFloat("pzzgg_VAJHU");
   pzzzg_PS_VAJHU = cand.userFloat("pzzzg_PS_VAJHU");
   pzzgg_PS_VAJHU = cand.userFloat("pzzgg_PS_VAJHU");
   p0Zgs_VAJHU = cand.userFloat("p0Zgs_VAJHU");
   p0gsgs_VAJHU = cand.userFloat("p0gsgs_VAJHU");
   p0Zgs_PS_VAJHU = cand.userFloat("p0Zgs_PS_VAJHU");
   p0gsgs_PS_VAJHU = cand.userFloat("p0gsgs_PS_VAJHU");
   p0Zgs_g1prime2_VAJHU = cand.userFloat("p0Zgs_g1prime2_VAJHU");
   pzzzg_g1prime2_VAJHU = cand.userFloat("pzzzg_g1prime2_VAJHU");
   pzzzg_g1prime2_pi2_VAJHU = cand.userFloat("pzzzg_g1prime2_pi2_VAJHU");

   p0plus_m4l = cand.userFloat("p0plus_m4l");
   bkg_m4l = cand.userFloat("bkg_m4l");
   p0plus_m4l_ScaleUp = cand.userFloat("p0plus_m4l_ScaleUp");// signal m4l probability for systematics
   bkg_m4l_ScaleUp = cand.userFloat("bkg_m4l_ScaleUp");// backgroun m4l probability for systematics
   p0plus_m4l_ScaleDown = cand.userFloat("p0plus_m4l_ScaleDown");// signal m4l probability for systematics
   bkg_m4l_ScaleDown = cand.userFloat("bkg_m4l_ScaleDown");// backgroun m4l probability for systematics
   p0plus_m4l_ResUp = cand.userFloat("p0plus_m4l_ResUp");// signal m4l probability for systematics
   bkg_m4l_ResUp = cand.userFloat("bkg_m4l_ResUp");// backgroun m4l probability for systematics
   p0plus_m4l_ResDown = cand.userFloat("p0plus_m4l_ResDown");// signal m4l probability for systematics
   bkg_m4l_ResDown = cand.userFloat("bkg_m4l_ResDown");// backgroun m4l probability for systematics

   phjj_VAJHU_old = cand.userFloat("phjj_VAJHU_old");
   pvbf_VAJHU_old = cand.userFloat("pvbf_VAJHU_old");
   phjj_VAJHU_old_up = cand.userFloat("phjj_VAJHU_old_up");
   pvbf_VAJHU_old_up = cand.userFloat("pvbf_VAJHU_old_up");
   phjj_VAJHU_old_dn = cand.userFloat("phjj_VAJHU_old_dn");
   pvbf_VAJHU_old_dn = cand.userFloat("pvbf_VAJHU_old_dn");
   phjj_VAJHU_new = cand.userFloat("phjj_VAJHU_new");
   pvbf_VAJHU_new = cand.userFloat("pvbf_VAJHU_new");
   phjj_VAJHU_new_up = cand.userFloat("phjj_VAJHU_new_up");
   pvbf_VAJHU_new_up = cand.userFloat("pvbf_VAJHU_new_up");
   phjj_VAJHU_new_dn = cand.userFloat("phjj_VAJHU_new_dn");
   pvbf_VAJHU_new_dn = cand.userFloat("pvbf_VAJHU_new_dn");

   pAux_vbf_VAJHU = cand.userFloat("pAux_vbf_VAJHU");
   pAux_vbf_VAJHU_up = cand.userFloat("pAux_vbf_VAJHU_up");
   pAux_vbf_VAJHU_dn = cand.userFloat("pAux_vbf_VAJHU_dn");

   phj_VAJHU = cand.userFloat("phj_VAJHU");
   phj_VAJHU_up = cand.userFloat("phj_VAJHU_up");
   phj_VAJHU_dn = cand.userFloat("phj_VAJHU_dn");

   pwh_hadronic_VAJHU = cand.userFloat("pwh_hadronic_VAJHU");
   pwh_hadronic_VAJHU_up = cand.userFloat("pwh_hadronic_VAJHU_up");
   pwh_hadronic_VAJHU_dn = cand.userFloat("pwh_hadronic_VAJHU_dn");

   pzh_hadronic_VAJHU = cand.userFloat("pzh_hadronic_VAJHU");
   pzh_hadronic_VAJHU_up = cand.userFloat("pzh_hadronic_VAJHU_up");
   pzh_hadronic_VAJHU_dn = cand.userFloat("pzh_hadronic_VAJHU_dn");

   ptth_VAJHU = cand.userFloat("ptth_VAJHU");
   ptth_VAJHU_up = cand.userFloat("ptth_VAJHU_up");
   ptth_VAJHU_dn = cand.userFloat("ptth_VAJHU_dn");

   pbbh_VAJHU = cand.userFloat("pbbh_VAJHU");
   pbbh_VAJHU_up = cand.userFloat("pbbh_VAJHU_up");
   pbbh_VAJHU_dn = cand.userFloat("pbbh_VAJHU_dn");


  //Int_t isSignal = -1;
  //Int_t isRightPair = -1;
  //if(isMC){
  //  isSignal = cand.userFloat("MC_isRight");
  //  isRightPair = cand.userFloat("MC_isRightPair");
  //}

  //Float_t mZa = cand.userFloat("mZa");
  //Float_t mZb = cand.userFloat("mZb");
  //Float_t mLL4 = cand.userFloat("mLL4");
  //Float_t mLL6 = cand.userFloat("mLL6");
  //Float_t SIP4 = cand.userFloat("SIP4");
  //Float_t iso34 = cand.userFloat("iso34");

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
   Z1Flav = Z1->daughter(0)->pdgId()*Z1->daughter(1)->pdgId();
  
   Z2Mass = Z2->mass();
   Z2Pt =   Z2->pt();
   Z2Flav = Z2->daughter(0)->pdgId()*Z2->daughter(1)->pdgId();

  // Precomputed selections
  bool candPass70Z2Loose   = cand.userFloat("Z2Mass") && 
                             cand.userFloat("MAllComb") &&
                             cand.userFloat("pt1")>20 && cand.userFloat("pt2")>10. &&
                             ZZMass>70.;
  bool candPassFullSel70   = cand.userFloat("SR");
  bool candPassFullSel   = cand.userFloat("FullSel");
  bool candIsBest = cand.userFloat("isBestCand");
  bool passMz_zz = (Z1Mass>60. && Z1Mass<120. && Z2Mass>60. && Z2Mass<120.);   //FIXME hardcoded cut

  Int_t sel = 0;

//   if(candPassFullSel){ // FIXME: this step (fullSel, without candIsBest) is not part of the normal selection flow.
//     sel = 50;
//   }
  
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

        
  //FIXME? Not filled in the miniAOD branch, is it right?
  /*
    float Hadditional[6]={mZa, mZb, mLL4, mLL6, SIP4, iso34};                   
    TString addName[6]={};
    myTree->SetVariables(addName,Hadditional,6);
  */

  //Fill the angular variables
   helcosthetaZ1 = cand.userFloat("costheta1");
   helcosthetaZ2 = cand.userFloat("costheta2");
   helphi       = cand.userFloat("phi");
   costhetastar = cand.userFloat("costhetastar");
   phistarZ1      = cand.userFloat("phistar1");
   phistarZ2      = cand.userFloat("phistar2");
   xi            = cand.userFloat("xi");
   xistar        = cand.userFloat("xistar");

  // Retrieve the userFloat of the leptons in vectors ordered in the same way.
  vector<float> SIP(4);
  vector<float> PFChargedHadIso(4);
  vector<float> PFNeutralHadIso(4);
  vector<float> PFPhotonIso(4);
  vector<float> combRelIsoPF(4);
  vector<bool>  isID(4);
  
  for (unsigned int i=0; i<leptons.size(); ++i){
    SIP[i]             = userdatahelpers::getUserFloat(leptons[i],"SIP");
    PFChargedHadIso[i] = userdatahelpers::getUserFloat(leptons[i],"PFChargedHadIso");
    PFNeutralHadIso[i] = userdatahelpers::getUserFloat(leptons[i],"PFNeutralHadIso");
    PFPhotonIso[i]     = userdatahelpers::getUserFloat(leptons[i],"PFPhotonIso");
    isID[i]            = userdatahelpers::getUserFloat(leptons[i],"ID");

    if (theChannel==ZL) {
      combRelIsoPF[i]    = userdatahelpers::getUserFloat(leptons[i],"combRelIsoPF");
      //FIXME cannot take labels[i]+"SIP", that info only attached to the Z!!
    } else {
      combRelIsoPF[i]    = cand.userFloat(labels[i]+"combRelIsoPFFSRCorr"); // Note: the FSR-corrected iso is attached to the Z, not to the lepton!
      // Check that I don't mess up with labels[] and leptons[]
      //cout<<"SIP: "<<SIP[i]<<"  labels "<<labels[i]<<cand.userFloat(labels[i]+"SIP")<<endl;
      SIP[i] = cand.userFloat(labels[i]+"SIP");
    }

    //Fill the info on the lepton candidates  
    if(i==0){
       Lep1Pt  = leptons[i]->pt();
       Lep1Eta  = leptons[i]->eta();
       Lep1Phi  = leptons[i]->phi();
       Lep1LepId  = leptons[i]->pdgId();
       Lep1SIP  = SIP[i];
       Lep1isID  = isID[i];
       Lep1BDT  = userdatahelpers::getUserFloat(leptons[i],"BDT");
       Lep1missingHit  = userdatahelpers::getUserFloat(leptons[i],"missingHit");
       Lep1combRelIsoPF=combRelIsoPF[i];
    }else if(i==1){
       Lep2Pt  = leptons[i]->pt();
       Lep2Eta  = leptons[i]->eta();
       Lep2Phi  = leptons[i]->phi();
       Lep2LepId  = leptons[i]->pdgId();
       Lep2SIP  = SIP[i];
       Lep2isID  = isID[i];
       Lep2BDT  = userdatahelpers::getUserFloat(leptons[i],"BDT");
       Lep2missingHit  = userdatahelpers::getUserFloat(leptons[i],"missingHit");
       Lep2combRelIsoPF=combRelIsoPF[i];
    }else if(i==2){
       Lep3Pt  = leptons[i]->pt();
       Lep3Eta  = leptons[i]->eta();
       Lep3Phi  = leptons[i]->phi();
       Lep3LepId  = leptons[i]->pdgId();
       Lep3SIP  = SIP[i];
       Lep3isID  = isID[i];
       Lep3BDT  = userdatahelpers::getUserFloat(leptons[i],"BDT");
       Lep3missingHit  = userdatahelpers::getUserFloat(leptons[i],"missingHit");
       Lep3combRelIsoPF=combRelIsoPF[i];
    }else if(i==3){
       Lep4Pt  = leptons[i]->pt();
       Lep4Eta  = leptons[i]->eta();
       Lep4Phi  = leptons[i]->phi();
       Lep4LepId  = leptons[i]->pdgId();
       Lep4SIP  = SIP[i];
       Lep4isID  = isID[i];
       Lep4BDT  = userdatahelpers::getUserFloat(leptons[i],"BDT");
       Lep4missingHit  = userdatahelpers::getUserFloat(leptons[i],"missingHit");
       Lep4combRelIsoPF=combRelIsoPF[i] ;   
    }
  }

  // FSR 
  for (unsigned i=0; i<fsrPhot.size(); ++i) { 
    math::XYZTLorentzVector fsr = fsrPhot[i]->p4();
    fsrPt.push_back(fsr.pt());
    fsrEta.push_back(fsr.eta());
    fsrPhi.push_back(fsr.phi());
    fsrLept.push_back(fsrIndex[i]);
    //    fsrDR.push_back();
    fsrLeptID.push_back(leptons[fsrIndex[i]-1]->pt());
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
    for(int izx=0;izx<2;izx++)
      ZXFakeweight *= getFakeWeight(Z2->daughter(izx)->pt(),Z2->daughter(izx)->eta(),Z2->daughter(izx)->pdgId(),Z1->daughter(0)->pdgId());
  }
  ZZsel = sel;

  // FIXME: Other variables to be added if addKinRefit or addVtxFit == true:
  // Z1MassRefit, ZZMassRefit, ZZChi2KinFit, ZZMassCFit, ZZChi2CFit

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
}



// ------------ method called once each job just before starting event loop  ------------
void HZZ4lNtupleMaker::beginJob()
{
  edm::Service<TFileService> fs;
  myTree = new HZZ4lNtupleFactory( fs->make<TTree>(theFileName,"Event Summary"));
  hCounter = fs->make<TH1F>("Counters", "Counters", 45, 0., 45.);
  BookAllBranches();
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
void HZZ4lNtupleMaker::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void HZZ4lNtupleMaker::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void HZZ4lNtupleMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
  Nevt_Gen_lumiBlock = 0;
}

// ------------ method called when ending the processing of a luminosity block  ------------
void HZZ4lNtupleMaker::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup)
{
  Float_t Nevt_preskim = -1.;
  edm::Handle<edm::MergeableCounter> preSkimCounter;
  if (iLumi.getByLabel("preSkimCounter", preSkimCounter)) { // Counter before skim. Does not exist for non-skimmed samples.
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
void HZZ4lNtupleMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


Float_t HZZ4lNtupleMaker::getAllWeight(Float_t LepPt, Float_t LepEta, Int_t LepID) const
{
  if (skipDataMCWeight) return 1.;

  Float_t weight  = 1.; 
  //Float_t errCorr = 0.;
  //Float_t errCorrSyst = 0.;

  Float_t myLepPt = LepPt;
  Float_t myLepEta = LepEta;
  Int_t   myLepID = abs(LepID);
  
  //avoid to go out of the TH boundary
  if(myLepID == 13 && myLepPt > 99.) myLepPt = 99.;
  if(myLepID == 11 && myLepPt > 199.) myLepPt = 199.;
  if(myLepID == 11) myLepEta = fabs(myLepEta);

  if(myLepID == 13){                                               
      weight  = hTH2D_Mu_All->GetBinContent(hTH2D_Mu_All->GetXaxis()->FindBin(myLepPt),hTH2D_Mu_All->GetYaxis()->FindBin(LepEta));
    //  errCorr = hTH2D_Mu_All_2011A->GetBinError(hTH2D_Mu_All_2011A->GetXaxis()->FindBin(myLepPt),hTH2D_Mu_All_2011A->GetYaxis()->FindBin(LepEta));
      
  }
  else if(myLepID == 11){   

    weight  = hTH2D_El_All->GetBinContent(hTH2D_El_All->GetXaxis()->FindBin(myLepPt),hTH2D_El_All->GetYaxis()->FindBin(myLepEta));
      //errCorr = hTH2D_El_All_2011A->GetBinError(hTH2D_El_All_2011A->GetXaxis()->FindBin(myLepPt),hTH2D_El_All_2011A->GetYaxis()->FindBin(myLepEta));   
  }else {
    cout<<"ERROR! wrong lepton ID "<<myLepID<<endl;
    //abort();
    weight=0;
  }   

  
  //add the systematics on T&P corrections (for muons only, electrons have them already included)
  //if(myLepID == 13){
  //  if(myLepPt >= 15.) errCorrSyst = 0.005;
  //  else errCorrSyst = 0.015;
  //}

  //FIXME
  if(myLepPt < 5. && myLepID == 13) weight = 1.;

  if(weight < 0.001 || weight > 10.){
    cout << "ERROR! LEP out of range! myLepPt = " << myLepPt << " myLepEta = " << myLepEta <<" myLepID "<<myLepID<< " weight = " << weight << endl;
    //abort();  //no correction should be zero, if you find one, stop
  }
/*
//FIXME: in HZZ4l was working because seeder was the event number (1=first event, 10=10th event)
//Here I don't have this info. Random seed?
  static TRandom3 randomToss;
  if( ( myLepID == 13) || (myLepID == 11) ){

    //apply correlation matrix by assigning the proper seed
    Int_t CorrSeeder = Seeder;
    if(myLepID == 13){
      if(myLepPt < 20. && fabs(myLepEta) < 1.2) CorrSeeder += 100001;
      else if(myLepPt < 20. && fabs(myLepEta) >= 1.2) CorrSeeder += 100002;
    }

    randomToss.SetSeed(CorrSeeder);
    weight = randomToss.Gaus(weight,errCorr);

    //apply systematic (totally correlated in eta) for muons
    if(myLepID == 13){
      randomToss.SetSeed(Seeder);
      weight = randomToss.Gaus(weight,errCorrSyst);
    }
  }
*/
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
void HZZ4lNtupleMaker::FillZGenInfo(const math::XYZTLorentzVector pZ1, const math::XYZTLorentzVector pZ2)
{
  GenZ1Mass= pZ1.M();
  GenZ1Pt= pZ1.Pt();

  GenZ2Mass= pZ2.M();
  GenZ2Pt= pZ2.Pt();

  return;
}

void HZZ4lNtupleMaker::FillLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, Short_t Lep3Id, Short_t Lep4Id, 
				      const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2, 
				      const math::XYZTLorentzVector Lep3, const math::XYZTLorentzVector Lep4,float weight)
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

  dataMCWeight=weight;
  
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
  myTree->Book("PUWeight",PUWeight);
  myTree->Book("PFMET",PFMET);
  myTree->Book("PFMETPhi",PFMETPhi);
  myTree->Book("nCleanedJets",nCleanedJets);
  myTree->Book("nCleanedJetsPt30",nCleanedJetsPt30);
  myTree->Book("nCleanedJetsPt30BTagged",nCleanedJetsPt30BTagged);
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
    myTree->Book("ZZChi2KinFit",ZZChi2KinFit);
    myTree->Book("Z1MassRefit",Z1MassRefit);
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
  myTree->Book("Lep1Pt",Lep1Pt);
  myTree->Book("Lep1Eta",Lep1Eta);
  myTree->Book("Lep1Phi",Lep1Phi);
  myTree->Book("Lep1LepId",Lep1LepId);
  myTree->Book("Lep1SIP",Lep1SIP);
  myTree->Book("Lep1isID",Lep1isID);
  myTree->Book("Lep1BDT",Lep1BDT);
  myTree->Book("Lep1missingHit",Lep1missingHit);
  myTree->Book("Lep2Pt",Lep2Pt);
  myTree->Book("Lep2Eta",Lep2Eta);
  myTree->Book("Lep2Phi",Lep2Phi);
  myTree->Book("Lep2LepId",Lep2LepId);
  myTree->Book("Lep2SIP",Lep2SIP);
  myTree->Book("Lep2isID",Lep2isID);
  myTree->Book("Lep2BDT",Lep2BDT);
  myTree->Book("Lep2missingHit",Lep2missingHit);
  myTree->Book("Lep3Pt",Lep3Pt);
  myTree->Book("Lep3Eta",Lep3Eta);
  myTree->Book("Lep3Phi",Lep3Phi);
  myTree->Book("Lep3LepId",Lep3LepId);
  myTree->Book("Lep3SIP",Lep3SIP);
  myTree->Book("Lep3isID",Lep3isID);
  myTree->Book("Lep3BDT",Lep3BDT);
  myTree->Book("Lep3missingHit",Lep3missingHit);
  myTree->Book("Lep4Pt",Lep4Pt);
  myTree->Book("Lep4Eta",Lep4Eta);
  myTree->Book("Lep4Phi",Lep4Phi);
  myTree->Book("Lep4LepId",Lep4LepId);
  myTree->Book("Lep4SIP",Lep4SIP);
  myTree->Book("Lep4isID",Lep4isID);
  myTree->Book("Lep4BDT",Lep4BDT);
  myTree->Book("Lep4missingHit",Lep4missingHit);
  //myTree->Book("Lep1chargedHadIso",Lep1chargedHadIso);
  //myTree->Book("Lep1neutralHadIso",Lep1neutralHadIso);
  //myTree->Book("Lep1photonIso",Lep1photonIso);
  myTree->Book("Lep1combRelIsoPF",Lep1combRelIsoPF);
  //myTree->Book("Lep2chargedHadIso",Lep2chargedHadIso);
  //myTree->Book("Lep2neutralHadIso",Lep2neutralHadIso);
  //myTree->Book("Lep2photonIso",Lep2photonIso);
  myTree->Book("Lep2combRelIsoPF",Lep2combRelIsoPF);
  //myTree->Book("Lep3chargedHadIso",Lep3chargedHadIso);
  //myTree->Book("Lep3neutralHadIso",Lep3neutralHadIso);
  //myTree->Book("Lep3photonIso",Lep3photonIso);
  myTree->Book("Lep3combRelIsoPF",Lep3combRelIsoPF);
  //myTree->Book("Lep4chargedHadIso",Lep4chargedHadIso);
  //myTree->Book("Lep4neutralHadIso",Lep4neutralHadIso);
  //myTree->Book("Lep4photonIso",Lep4photonIso);
  myTree->Book("Lep4combRelIsoPF",Lep4combRelIsoPF);
  myTree->Book("fsrPt",fsrPt);
  myTree->Book("fsrEta",fsrEta);
  myTree->Book("fsrPhi",fsrPhi);
  myTree->Book("fsrLept",fsrLept);
  myTree->Book("fsrLeptId",fsrLeptID); // FIXME next ones can be skipped for mass production
  myTree->Book("fsrGenPt",fsrGenPt);

  //Photon variables
  if (writePhotons) {
    myTree->Book("PhotPt",PhotPt);
    myTree->Book("PhotEta",PhotEta);
    myTree->Book("PhotPhi",PhotPhi);
  }

  //Discriminants
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
    myTree->Book("p0plus_m4l",p0plus_m4l);
    myTree->Book("bkg_m4l",bkg_m4l);
    myTree->Book("p0plus_m4l_ScaleUp",p0plus_m4l_ScaleUp);
    myTree->Book("bkg_m4l_ScaleUp",bkg_m4l_ScaleUp);
    myTree->Book("p0plus_m4l_ScaleDown",p0plus_m4l_ScaleDown);
    myTree->Book("bkg_m4l_ScaleDown",bkg_m4l_ScaleDown);
    myTree->Book("p0plus_m4l_ResUp",p0plus_m4l_ResUp);
    myTree->Book("bkg_m4l_ResUp",bkg_m4l_ResUp);
    myTree->Book("p0plus_m4l_ResDown",p0plus_m4l_ResDown);
    myTree->Book("bkg_m4l_ResDown",bkg_m4l_ResDown);

  //Production MELA
  myTree->Book("phjj_VAJHU_old",phjj_VAJHU_old);
  myTree->Book("pvbf_VAJHU_old",pvbf_VAJHU_old);
  myTree->Book("phjj_VAJHU_old_up",phjj_VAJHU_old_up);
  myTree->Book("pvbf_VAJHU_old_up",pvbf_VAJHU_old_up);
  myTree->Book("phjj_VAJHU_old_dn",phjj_VAJHU_old_dn);
  myTree->Book("pvbf_VAJHU_old_dn",pvbf_VAJHU_old_dn);
  myTree->Book("phjj_VAJHU_new",phjj_VAJHU_new);
  myTree->Book("pvbf_VAJHU_new",pvbf_VAJHU_new);
  myTree->Book("phjj_VAJHU_new_up",phjj_VAJHU_new_up);
  myTree->Book("pvbf_VAJHU_new_up",pvbf_VAJHU_new_up);
  myTree->Book("phjj_VAJHU_new_dn",phjj_VAJHU_new_dn);
  myTree->Book("pvbf_VAJHU_new_dn",pvbf_VAJHU_new_dn);
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
  
  //Jet variables
  myTree->Book("JetPt",JetPt);
  myTree->Book("JetEta",JetEta);
  myTree->Book("JetPhi",JetPhi);
  myTree->Book("JetMass",JetMass);
  myTree->Book("JetBTagger",JetBTagger);
  myTree->Book("JetIsBtagged",JetIsBtagged);
  myTree->Book("JetQGLikelihood",JetQGLikelihood);
  myTree->Book("JetSigma",JetSigma);
  myTree->Book("DiJetMass",DiJetMass);
  myTree->Book("DiJetMassPlus",DiJetMassPlus);
  myTree->Book("DiJetMassMinus",DiJetMassMinus);
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
  myTree->Book("genExtInfo",genExtInfo);
  myTree->Book("xsec",xsection);
  myTree->Book("dataMCWeight",dataMCWeight);
  myTree->Book("HqTMCweight",HqTMCweight);
  myTree->Book("ZXFakeweight",ZXFakeweight);
  myTree->Book("GenHMass",GenHMass);
  myTree->Book("GenHPt",GenHPt);
  myTree->Book("GenHRapidity",GenHRapidity);
  myTree->Book("GenZ1Mass",GenZ1Mass);
  myTree->Book("GenZ1Pt",GenZ1Pt);
  myTree->Book("GenZ2Mass",GenZ2Mass);
  myTree->Book("GenZ2Pt",GenZ2Pt);
  myTree->Book("GenLep1Pt",GenLep1Pt);
  myTree->Book("GenLep1Eta",GenLep1Eta);
  myTree->Book("GenLep1Phi",GenLep1Phi);
  myTree->Book("GenLep1Id",GenLep1Id);
  myTree->Book("GenLep2Pt",GenLep2Pt);
  myTree->Book("GenLep2Eta",GenLep2Eta);
  myTree->Book("GenLep2Phi",GenLep2Phi);
  myTree->Book("GenLep2Id",GenLep2Id);
  myTree->Book("GenLep3Pt",GenLep3Pt);
  myTree->Book("GenLep3Eta",GenLep3Eta);
  myTree->Book("GenLep3Phi",GenLep3Phi);
  myTree->Book("GenLep3Id",GenLep3Id);
  myTree->Book("GenLep4Pt",GenLep4Pt);
  myTree->Book("GenLep4Eta",GenLep4Eta);
  myTree->Book("GenLep4Phi",GenLep4Phi);
  myTree->Book("GenLep4Id",GenLep4Id);
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
DEFINE_FWK_MODULE(HZZ4lNtupleMaker);
