#ifndef HZZ4lNtupleFactory_h
#define HZZ4lNtupleFactory_h

#include <vector>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>

class HZZ4lNtupleFactory{
  
 protected:
  
 public:
  HZZ4lNtupleFactory(TTree* outTree_input);
  ~HZZ4lNtupleFactory();
  
  //Fill the tree and initialize the variables
  void FillEvent();
  
  //Fill the tree without resetting the variables
  void FillCurrentTree(){_outTree->Fill();}
  void DumpBranches(TString filename) const;

  void createNewCandidate();
  void FillEventInfo(Int_t RunNumber, 
                     Long64_t EventNumber, 
                     Int_t LumiNumber, 
                     Int_t IndexBestCand, 
                     Int_t Nvtx, 
                     Int_t NObsInt, 
                     Float_t NTrueInt, 
                     Float_t PUweight, 
                     Float_t PFMET, 
                     Int_t nJets, 
                     Int_t nCleanedJets, 
                     Int_t nCleanedJetsPt30,
                     Int_t nCleanedJetsPt30BTagged,
                     Int_t genFinalState, 
                     Int_t genProcessId, 
                     Float_t genHEPMCweight, 
                     Short_t trigWord, 
                     Short_t genExtInfo,
                     Float_t xsec, 
                     Float_t weight
                     );
  void FillHInfo(Float_t ZZMass, 
                 Float_t ZZMassErr, 
                 Float_t ZZMassErrCorr, 
                 Float_t ZZMassPreFSR, 
                 Float_t ZZMassRefit, 
                 Float_t Chi2KinFit, 
                 Float_t ZZMassCFit, 
                 Float_t Chi2CFit, 
                 Int_t ZZsel, 
                 Float_t ZZPt, 
                 Float_t ZZEta,
                 Float_t ZZPhi, 
                 Int_t isSignal, 
                 Int_t isRightPair, 
                 Int_t CRflag=0
                 );
  void FillProbability(Float_t p0plus_VAJHU,
                       Float_t p0minus_VAJHU,
                       Float_t p0plus_VAMCFM,
                       Float_t p0hplus_VAJHU,
                       Float_t p1_VAJHU,
                       Float_t p1_prodIndep_VAJHU,
                       Float_t p1plus_VAJHU,
                       Float_t p1plus_prodIndep_VAJHU,
                       Float_t p2_VAJHU,
                       Float_t p2_prodIndep_VAJHU,
                       Float_t p2qqb_VAJHU,
                       Float_t p2hplus_VAJHU,
                       Float_t p2hminus_VAJHU,
                       Float_t p2bplus_VAJHU,
                       Float_t p2hplus_qqb_VAJHU,                                                                         
                       Float_t p2hplus_prodIndep_VAJHU,         
                       Float_t p2hminus_qqb_VAJHU,                              
                       Float_t p2hminus_prodIndep_VAJHU,        
                       Float_t p2bplus_qqb_VAJHU,                                       
                       Float_t p2bplus_prodIndep_VAJHU,         
                       Float_t p2h2plus_gg_VAJHU,               
                       Float_t p2h2plus_qqbar_VAJHU,            
                       Float_t p2h2plus_prodIndep_VAJHU,        
                       Float_t p2h3plus_gg_VAJHU,               
                       Float_t p2h3plus_qqbar_VAJHU,            
                       Float_t p2h3plus_prodIndep_VAJHU,        
                       Float_t p2h6plus_gg_VAJHU,               
                       Float_t p2h6plus_qqbar_VAJHU,            
                       Float_t p2h6plus_prodIndep_VAJHU,        
                       Float_t p2h7plus_gg_VAJHU,               
                       Float_t p2h7plus_qqbar_VAJHU,            
                       Float_t p2h7plus_prodIndep_VAJHU,        
                       Float_t p2h9minus_gg_VAJHU,              
                       Float_t p2h9minus_qqbar_VAJHU,           
                       Float_t p2h9minus_prodIndep_VAJHU,       
                       Float_t p2h10minus_gg_VAJHU,       
                       Float_t p2h10minus_qqbar_VAJHU,    
                       Float_t p2h10minus_prodIndep_VAJHU,
                       Float_t bkg_VAMCFM,
                       Float_t bkg_prodIndep_VAMCFM,
                       Float_t ggzz_VAMCFM,
                       Float_t ggzz_p0plus_VAMCFM,
                       Float_t ggzz_c1_VAMCFM,
                       Float_t ggzz_c5_VAMCFM,
                       Float_t ggzz_ci_VAMCFM,
                       Float_t phjj_VAJHU_old,
                       Float_t pvbf_VAJHU_old,
                       Float_t phjj_VAJHU_old_up,
                       Float_t pvbf_VAJHU_old_up,
                       Float_t phjj_VAJHU_old_dn,
                       Float_t pvbf_VAJHU_old_dn,
                       Float_t phjj_VAJHU_new,
                       Float_t pvbf_VAJHU_new,
                       Float_t phjj_VAJHU_new_up,
                       Float_t pvbf_VAJHU_new_up,
                       Float_t phjj_VAJHU_new_dn,
                       Float_t pvbf_VAJHU_new_dn,
                       Float_t p0_g1prime2_VAJHU,
                       Float_t pg1g1prime2_VAJHU,
                       Float_t Dgg10_VAMCFM,
                       Float_t pg1g4_mela,
                       Float_t pg1g4_VAJHU,
                       Float_t pg1g4_pi2_VAJHU,
                       Float_t pg1g2_pi2_VAJHU,
                       Float_t pg1g2_mela,
                       Float_t pg1g2_VAJHU,
                       Float_t pzzzg_VAJHU,
                       Float_t pzzgg_VAJHU,
                       Float_t pzzzg_PS_VAJHU,
                       Float_t pzzgg_PS_VAJHU,
                       Float_t p0Zgs_VAJHU,
                       Float_t p0gsgs_VAJHU,
                       Float_t p0Zgs_PS_VAJHU,
                       Float_t p0gsgs_PS_VAJHU
                       );
  void FillSuperMela(Float_t p0plus_m4l,
                     Float_t bkg_m4l,
                     Float_t p0plus_m4l_ScaleUp,
                     Float_t bkg_m4l_ScaleUp,
                     Float_t p0plus_m4l_ScaleDown,
                     Float_t bkg_m4l_ScaleDown,
                     Float_t p0plus_m4l_ResUp,
                     Float_t bkg_m4l_ResUp,
                     Float_t p0plus_m4l_ResDown,
                     Float_t bkg_m4l_ResDown
                     );
  void FillHAdditionalInfo(Float_t mZa, Float_t mZb, Float_t mLL4, Float_t mLL6, Float_t SIP4, Float_t iso34);
  void FillZInfo(Float_t ZMass, Float_t ZPt, short ZFlav, Float_t Z1MassRefit);
  void FillAngularInfo(Float_t costhetastar, Float_t phi, Float_t costheta1, Float_t costheta2, Float_t phistar1, Float_t phistar2,Float_t xi, Float_t xistar);
  void FillLepInfo(Float_t LepPt, Float_t LepEta, Float_t LepPhi, Int_t LepId, Float_t SIP, bool isID, float BDT, short parentId, int missingHit);
  void FillLepIsolInfo(Float_t LepchargedHadIso, Float_t LepneutralHadIso, Float_t LepphotonIso, Float_t LepcombRelIsoPF);
  void FillPhotonInfo(Float_t PhotPt, Float_t PhotEta, Float_t PhotPhi);
  void FillJetInfo(Float_t JetPt, Float_t JetEta, Float_t JetPhi, Float_t JetMass, Float_t JetBTagger, Float_t JetIsBtagged, Float_t JetQGLikelihood, Float_t JetSigma );
  void FillDiJetInfo(Float_t DiJetMass, Float_t DiJetMassPlus, Float_t DiJetMassMinus, Float_t DiJetDEta, Float_t DiJetFisher);
  void FillCategorizationInfo(Int_t nExtraLep, Int_t nExtraZ);
  void FillExtraLepInfo(int extraLeptonIndex, bool extraLeptonExists, reco::CandidatePtr ExtraLep);

  void FillHGenInfo(math::XYZTLorentzVector Hp);
  void FillZGenInfo(math::XYZTLorentzVector Z1p, math::XYZTLorentzVector Z2p);
  void FillLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, Short_t Lep3Id, Short_t Lep4Id,
                      math::XYZTLorentzVector Lep1, math::XYZTLorentzVector Lep2, math::XYZTLorentzVector Lep3, math::XYZTLorentzVector Lep4);
  void FillAssocLepGenInfo(std::vector<const reco::Candidate *>& AssocLeps);

  void InitializeVariables();

 private:

  TTree* _outTree;

  void InitializeBranches();

  Bool_t _firstZStored;
  Int_t _LeptonIndex;
  Int_t _LeptonIsoIndex;

  //Event variables
  Int_t _RunNumber;
  Long64_t _EventNumber;
  Int_t _LumiNumber;
  Int_t _IndexBestCand;

  //  Int_t _Nmu;
  //  Int_t _Nele;

  Int_t _Nvtx;
  Int_t _NObsInt;
  Float_t _NTrueInt;
  Float_t _PUWeight;
  Float_t _PFMET;
  Int_t _nJets;
  Int_t _nCleanedJets;
  Int_t _nCleanedJetsPt30;
  Int_t _nCleanedJetsPt30BTagged;
  Int_t _genFinalState;
  Int_t _genProcessId;
  Float_t _genHEPMCweight;
  Short_t _trigWord;
  Short_t _genExtInfo;
  Float_t _xsec;
  Float_t _gen_weight;
  
  //H variables
  Float_t _ZZMass;
  Float_t _ZZMassPreFSR;
  Float_t _ZZMassRefit;
  Float_t _ZZMassErr;
  Float_t _ZZMassErrCorr;
  Float_t _Chi2KinFit;
  Float_t _ZZMassCFit;
  Float_t _Chi2CFit;
  Int_t _ZZsel;
  Float_t _ZZPt;
  Float_t _ZZEta;
  Float_t _ZZPhi;
  Int_t _ZZgenIsSignal;
  Int_t _ZZgenIsRightPair;
  Int_t _CRflag;

  //probabilities
  Float_t _p0plus_VAJHU; // higgs, vector algebra, JHUgen
  Float_t _p0minus_VAJHU;// pseudoscalar, vector algebra, JHUgen
  Float_t _p0plus_VAMCFM;// higgs, vector algebra, MCFM
  Float_t _p0hplus_VAJHU;// 0h+ (high dimensional operator), vector algebra, JHUgen
  Float_t _p1_VAJHU;   // zprime, vector algebra, JHUgen,
  Float_t _p1_prodIndep_VAJHU;   // zprime, vector algebra, JHUgen,
  Float_t _p1plus_VAJHU;// 1+ (axial vector), vector algebra, JHUgen,
  Float_t _p1plus_prodIndep_VAJHU;// 1+ (axial vector), vector algebra, JHUgen,
  Float_t _p2_VAJHU;   // graviton, vector algebra, JHUgen,
  Float_t _p2_prodIndep_VAJHU;   // graviton, vector algebra, JHUgen,
  Float_t _p2qqb_VAJHU;   // graviton, vector algebra, JHUgen,
  Float_t _p2hplus_VAJHU;   // graviton, vector algebra, JHUgen,
  Float_t _p2hminus_VAJHU;   // graviton, vector algebra, JHUgen,
  Float_t _p2bplus_VAJHU;   // graviton, vector algebra, JHUgen,
  Float_t _p2hplus_qqb_VAJHU;        
  Float_t _p2hplus_prodIndep_VAJHU;  
  Float_t _p2hminus_qqb_VAJHU;       
  Float_t _p2hminus_prodIndep_VAJHU; 
  Float_t _p2bplus_qqb_VAJHU;        
  Float_t _p2bplus_prodIndep_VAJHU;  
  Float_t _p2h2plus_gg_VAJHU;        
  Float_t _p2h2plus_qqbar_VAJHU;     
  Float_t _p2h2plus_prodIndep_VAJHU; 
  Float_t _p2h3plus_gg_VAJHU;        
  Float_t _p2h3plus_qqbar_VAJHU;     
  Float_t _p2h3plus_prodIndep_VAJHU; 
  Float_t _p2h6plus_gg_VAJHU;        
  Float_t _p2h6plus_qqbar_VAJHU;     
  Float_t _p2h6plus_prodIndep_VAJHU; 
  Float_t _p2h7plus_gg_VAJHU;        
  Float_t _p2h7plus_qqbar_VAJHU;     
  Float_t _p2h7plus_prodIndep_VAJHU; 
  Float_t _p2h9minus_gg_VAJHU;       
  Float_t _p2h9minus_qqbar_VAJHU;    
  Float_t _p2h9minus_prodIndep_VAJHU;
  Float_t _p2h10minus_gg_VAJHU;       
  Float_t _p2h10minus_qqbar_VAJHU;    
  Float_t _p2h10minus_prodIndep_VAJHU;
  //backgrounds
  Float_t _bkg_VAMCFM; // background, vector algebra, MCFM
  Float_t _bkg_prodIndep_VAMCFM; // background, vector algebra, MCFM
  Float_t _ggzz_VAMCFM; // background, vector algebra, MCFM for ggzz
  Float_t _ggzz_p0plus_VAMCFM; // background, vector algebra, MCFM for ggzz
  Float_t _ggzz_c1_VAMCFM; // signal + background + interference w/ SM couplings, vector algebra, MCFM
  Float_t _ggzz_c5_VAMCFM; // signal + background + interference w/ 5xSM couplings, vector algebra, MCFM for ggzz
  Float_t _ggzz_ci_VAMCFM; // signal + background + interference w/ imaginary SM couplings, vector algebra, MCFM for ggzz  
  // supermela
  Float_t _p0plus_m4l;  // signal m4l probability as in datacards
  Float_t _bkg_m4l;     // backgroun m4l probability as in datacards
  Float_t _pg1g4_mela;
  Float_t _pg1g4_VAJHU;
  Float_t _pg1g4_pi2_VAJHU;
  Float_t _pg1g2_pi2_VAJHU;
  Float_t _pg1g2_mela;
  Float_t _pg1g2_VAJHU;
  Float_t _p0plus_m4l_ScaleUp;// signal m4l probability for systematics
  Float_t _bkg_m4l_ScaleUp;// backgroun m4l probability for systematics
  Float_t _p0plus_m4l_ScaleDown;// signal m4l probability for systematics
  Float_t _bkg_m4l_ScaleDown;// backgroun m4l probability for systematics
  Float_t _p0plus_m4l_ResUp;// signal m4l probability for systematics
  Float_t _bkg_m4l_ResUp;// backgroun m4l probability for systematics
  Float_t _p0plus_m4l_ResDown;// signal m4l probability for systematics
  Float_t _bkg_m4l_ResDown;// backgroun m4l probability for systematics
  // Production MELA
  Float_t _phjj_VAJHU_old; //H+jj, vector algebra, JHUGen, Legacy Jet Selection
  Float_t _pvbf_VAJHU_old; //VBF, vector algebra, JHUGen, Legacy Jet Selection
  Float_t _phjj_VAJHU_old_up; //H+jj, vector algebra, JHUGen, Legacy Jet Selection, JEC +
  Float_t _pvbf_VAJHU_old_up; //VBF, vector algebra, JHUGen, Legacy Jet Selection, JEC +
  Float_t _phjj_VAJHU_old_dn; //H+jj, vector algebra, JHUGen, Legacy Jet Selection, JEC -
  Float_t _pvbf_VAJHU_old_dn; //VBF, vector algebra, JHUGen, Legacy Jet Selection, JEC -
  Float_t _phjj_VAJHU_new; //H+jj, vector algebra, JHUGen, New Jet Selection
  Float_t _pvbf_VAJHU_new; //VBF, vector algebra, JHUGen, New Jet Selection
  Float_t _phjj_VAJHU_new_up; //H+jj, vector algebra, JHUGen, New Jet Selection, JEC +
  Float_t _pvbf_VAJHU_new_up; //VBF, vector algebra, JHUGen, New Jet Selection, JEC +
  Float_t _phjj_VAJHU_new_dn; //H+jj, vector algebra, JHUGen, New Jet Selection, JEC -
  Float_t _pvbf_VAJHU_new_dn; //VBF, vector algebra, JHUGen, New Jet Selection, JEC -
  Float_t _p0_g1prime2_VAJHU; //VBF, vector algebra, JHUGen, New Jet Selection, JEC -
  Float_t _pg1g1prime2_VAJHU; //VBF, vector algebra, JHUGen, New Jet Selection, JEC -
  Float_t _Dgg10_VAMCFM; //VBF, vector algebra, JHUGen, New Jet Selection, JEC -

  Float_t _pzzzg_VAJHU;
  Float_t _pzzgg_VAJHU;
  Float_t _pzzzg_PS_VAJHU;
  Float_t _pzzgg_PS_VAJHU;
  Float_t _p0Zgs_VAJHU;
  Float_t _p0gsgs_VAJHU;
  Float_t _p0Zgs_PS_VAJHU;
  Float_t _p0gsgs_PS_VAJHU;

  Float_t _ZZmZa;
  Float_t _ZZmZb;
  Float_t _ZZmLL4;
  Float_t _ZZmLL6;
  Float_t _ZZSIP4;

  //Z1 variables
  Float_t _Z1Mass;
  Float_t _Z1MassRefit;
  Float_t _Z1Pt;
  Short_t _Z1Flav;

  //Z2 variables
  Float_t _Z2Mass;
  Float_t _Z2Pt;
  Short_t _Z2Flav;

  //Angular variables
  Float_t _costhetastar;
  Float_t _phi;
  Float_t _costheta1;
  Float_t _costheta2;
  Float_t _phistar1;
  Float_t _phistar2;
  Float_t _xi;
  Float_t _xistar;

  //Lepton variables
  Float_t _Lep1Pt;
  Float_t _Lep1Eta;
  Float_t _Lep1Phi;
  Int_t   _Lep1LepId;
  Float_t _Lep1SIP;
  Bool_t  _Lep1isID;
  Float_t _Lep1BDT;
  Char_t  _Lep1missingHit;
  Short_t _Lep1ParentId;

  Float_t _Lep2Pt;
  Float_t _Lep2Eta;
  Float_t _Lep2Phi;
  Int_t   _Lep2LepId;
  Float_t _Lep2SIP;
  Bool_t  _Lep2isID;
  Float_t _Lep2BDT;
  Char_t  _Lep2missingHit;
  Short_t _Lep2ParentId;

  Float_t _Lep3Pt;
  Float_t _Lep3Eta;
  Float_t _Lep3Phi;
  Int_t   _Lep3LepId;
  Float_t _Lep3SIP;
  Bool_t  _Lep3isID;
  Float_t _Lep3BDT;
  Char_t  _Lep3missingHit;
  Short_t _Lep3ParentId;

  Float_t _Lep4Pt;
  Float_t _Lep4Eta;
  Float_t _Lep4Phi;
  Int_t   _Lep4LepId;
  Float_t _Lep4SIP;
  Bool_t  _Lep4isID;
  Float_t _Lep4BDT;
  Char_t  _Lep4missingHit;
  Short_t _Lep4ParentId;


  //Lepton isolation variables
  Float_t _Lep1chargedHadIso;
  Float_t _Lep1neutralHadIso;
  Float_t _Lep1photonIso;
  Float_t _Lep1combRelIsoPF;

  Float_t _Lep2chargedHadIso;
  Float_t _Lep2neutralHadIso;
  Float_t _Lep2photonIso;
  Float_t _Lep2combRelIsoPF;

  Float_t _Lep3chargedHadIso;
  Float_t _Lep3neutralHadIso;
  Float_t _Lep3photonIso;
  Float_t _Lep3combRelIsoPF;

  Float_t _Lep4chargedHadIso;
  Float_t _Lep4neutralHadIso;
  Float_t _Lep4photonIso;
  Float_t _Lep4combRelIsoPF;

  //Photon variables
  Float_t _PhotPt;
  Float_t _PhotEta;
  Float_t _PhotPhi;

  //Jet variables
  std::vector<Float_t> _JetPt;
  std::vector<Float_t> _JetEta;
  std::vector<Float_t> _JetPhi;
  std::vector<Float_t> _JetMass;
  std::vector<Float_t> _JetBTagger;
  std::vector<Float_t> _JetIsBtagged;
  std::vector<Float_t> _JetQGLikelihood;
  std::vector<Float_t> _JetSigma;
  Float_t _DiJetMass;
  Float_t _DiJetMassPlus;
  Float_t _DiJetMassMinus;
  Float_t _DiJetDEta;
  Float_t _DiJetFisher;

  //Categorization-related variables
  Int_t _nExtraLep;
  Int_t _nExtraZ;

  //Variables of extra leptons
  Float_t _ExtraLep1Pt;
  Float_t _ExtraLep1Eta;
  Float_t _ExtraLep1Phi;
  Int_t   _ExtraLep1LepId;
  Float_t _ExtraLep1SIP;
  Bool_t  _ExtraLep1isID;
  Float_t _ExtraLep1BDT;
  Char_t  _ExtraLep1missingHit;
  Float_t _ExtraLep1chargedHadIso;
  Float_t _ExtraLep1neutralHadIso;
  Float_t _ExtraLep1photonIso;
  Float_t _ExtraLep1combRelIsoPF;

  Float_t _ExtraLep2Pt;
  Float_t _ExtraLep2Eta;
  Float_t _ExtraLep2Phi;
  Int_t   _ExtraLep2LepId;
  Float_t _ExtraLep2SIP;
  Bool_t  _ExtraLep2isID;
  Float_t _ExtraLep2BDT;
  Char_t  _ExtraLep2missingHit;
  Float_t _ExtraLep2chargedHadIso;
  Float_t _ExtraLep2neutralHadIso;
  Float_t _ExtraLep2photonIso;
  Float_t _ExtraLep2combRelIsoPF;

  Float_t _ExtraLep3Pt;
  Float_t _ExtraLep3Eta;
  Float_t _ExtraLep3Phi;
  Int_t   _ExtraLep3LepId;
  Float_t _ExtraLep3SIP;
  Bool_t  _ExtraLep3isID;
  Float_t _ExtraLep3BDT;
  Char_t  _ExtraLep3missingHit;
  Float_t _ExtraLep3chargedHadIso;
  Float_t _ExtraLep3neutralHadIso;
  Float_t _ExtraLep3photonIso;
  Float_t _ExtraLep3combRelIsoPF;

  //Generation variables
  Float_t _genHMass;
  Float_t _genHPt;

  Float_t _genZ1Mass;
  Float_t _genZ1Pt;
  
  Float_t _genZ2Mass;
  Float_t _genZ2Pt;

  Float_t _genLep1Pt;
  Float_t _genLep1Eta;
  Float_t _genLep1Phi;
  Short_t _genLep1Id;

  Float_t _genLep2Pt;
  Float_t _genLep2Eta;
  Float_t _genLep2Phi;
  Short_t _genLep2Id;

  Float_t _genLep3Pt;
  Float_t _genLep3Eta;
  Float_t _genLep3Phi;
  Short_t _genLep3Id;

  Float_t _genLep4Pt;
  Float_t _genLep4Eta;
  Float_t _genLep4Phi;
  Short_t _genLep4Id;

  Float_t _genAssocLep1Pt;
  Float_t _genAssocLep1Eta;
  Float_t _genAssocLep1Phi;
  Short_t _genAssocLep1Id;

  Float_t _genAssocLep2Pt;
  Float_t _genAssocLep2Eta;
  Float_t _genAssocLep2Phi;
  Short_t _genAssocLep2Id;

};

#endif
