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

  void FillEvent();
  void DumpBranches(TString filename) const;

  void createNewCandidate();
  void FillEventInfo(const Int_t RunNumber, const Long64_t EventNumber, const Int_t LumiNumber, const Int_t IndexBestCand, Int_t Nvtx, Int_t NObsInt, Float_t NTrueInt, Float_t PUweight, Float_t PFMET, Int_t genFinalState, Int_t genProcessId, Float_t genHEPMCweight, Short_t trigWord, Short_t genExtInfo);
  void FillHInfo(const Float_t ZZMass, 
		 const Float_t ZZMassErr, 
		 const Float_t ZZMassErrCorr, 
		 const Float_t ZZMassPreFSR, 
		 const Float_t ZZMassRefit, 
		 const Float_t Chi2KinFit, 
		 const Float_t ZZMassCFit, 
		 const Float_t Chi2CFit, 
		 const Int_t ZZsel, 
		 const Float_t ZZPt, 
		 const Float_t ZZEta,
		 const Float_t ZZPhi, 
		 const Int_t isSignal, 
		 const Int_t isRightPair, 
		 const Float_t ZZFisher,
		 const Int_t CRflag=0
		 );
  void FillProbability(  //const Float_t p0plus_melaNorm,
			 //const Float_t p0plus_mela,
			 //const Float_t p0minus_mela,
			 //const Float_t p0hplus_mela,
			 const Float_t p0plus_VAJHU,
			 const Float_t p0minus_VAJHU,
			 const Float_t p0plus_VAMCFM,
			 const Float_t p0hplus_VAJHU,
			 //const Float_t p1_mela,
			 //const Float_t p1_prodIndep_mela,
			 //const Float_t p1plus_mela,
			 //const Float_t p1plus_prodIndep_mela,
			 const Float_t p1_VAJHU,
			 const Float_t p1_prodIndep_VAJHU,
			 const Float_t p1plus_VAJHU,
			 const Float_t p1plus_prodIndep_VAJHU,
			 //const Float_t p2_mela ,
			 //const Float_t p2_prodIndep_mela ,
			 //const Float_t p2qqb_mela ,
			 //const Float_t p2hplus_mela ,
			 //const Float_t p2hminus_mela ,
			 //const Float_t p2bplus_mela ,
			 const Float_t p2_VAJHU,
			 const Float_t p2_prodIndep_VAJHU,
			 const Float_t p2qqb_VAJHU,
			 const Float_t p2hplus_VAJHU,
			 const Float_t p2hminus_VAJHU,
			 const Float_t p2bplus_VAJHU,
		 	 const Float_t p2hplus_qqb_VAJHU,									  
		   const Float_t p2hplus_prodIndep_VAJHU,		
		   const Float_t p2hminus_qqb_VAJHU,				
		   const Float_t p2hminus_prodIndep_VAJHU,	
		   const Float_t p2bplus_qqb_VAJHU,					
		   const Float_t p2bplus_prodIndep_VAJHU,		
		   const Float_t p2h2plus_gg_VAJHU,      		
		   const Float_t p2h2plus_qqbar_VAJHU,   		
		   const Float_t p2h2plus_prodIndep_VAJHU,	
		   const Float_t p2h3plus_gg_VAJHU,       	
		   const Float_t p2h3plus_qqbar_VAJHU,    	
		   const Float_t p2h3plus_prodIndep_VAJHU,	
		   const Float_t p2h6plus_gg_VAJHU,       	
		   const Float_t p2h6plus_qqbar_VAJHU,    	
		   const Float_t p2h6plus_prodIndep_VAJHU,	
		   const Float_t p2h7plus_gg_VAJHU,       	
		   const Float_t p2h7plus_qqbar_VAJHU,    	
		   const Float_t p2h7plus_prodIndep_VAJHU,	
		   const Float_t p2h9minus_gg_VAJHU,       	
		   const Float_t p2h9minus_qqbar_VAJHU,    	
		   const Float_t p2h9minus_prodIndep_VAJHU,	
		   const Float_t p2h10minus_gg_VAJHU,       
		   const Float_t p2h10minus_qqbar_VAJHU,    
		   const Float_t p2h10minus_prodIndep_VAJHU,
			 //const Float_t bkg_mela,
			 const Float_t bkg_VAMCFM,
			 const Float_t bkg_prodIndep_VAMCFM,
			 const Float_t ggzz_VAMCFM,
			 const Float_t ggzz_p0plus_VAMCFM,
			 const Float_t ggzz_c1_VAMCFM,
			 const Float_t ggzz_c5_VAMCFM,
			 const Float_t ggzz_ci_VAMCFM,
			 const Float_t phjj_VAJHU_old,
			 const Float_t pvbf_VAJHU_old,
			 const Float_t phjj_VAJHU_old_up,
			 const Float_t pvbf_VAJHU_old_up,
			 const Float_t phjj_VAJHU_old_dn,
			 const Float_t pvbf_VAJHU_old_dn,
			 const Float_t phjj_VAJHU_new,
			 const Float_t pvbf_VAJHU_new,
			 const Float_t phjj_VAJHU_new_up,
			 const Float_t pvbf_VAJHU_new_up,
			 const Float_t phjj_VAJHU_new_dn,
			 const Float_t pvbf_VAJHU_new_dn,
       const Float_t p0_g1prime2_VAJHU,
       const Float_t pg1g1prime2_VAJHU,
       const Float_t Dgg10_VAMCFM,
			 const Float_t pg1g4_mela,
			 const Float_t pg1g4_VAJHU,
			 const Float_t pg1g4_pi2_VAJHU,
			 const Float_t pg1g2_pi2_VAJHU,
			 const Float_t pg1g2_mela,
			 const Float_t pg1g2_VAJHU,
			 const Float_t pzzzg_VAJHU,
			 const Float_t pzzgg_VAJHU,
			 const Float_t pzzzg_PS_VAJHU,
			 const Float_t pzzgg_PS_VAJHU,
			 const Float_t p0Zgs_VAJHU,
			 const Float_t p0gsgs_VAJHU,
			 const Float_t p0Zgs_PS_VAJHU,
			 const Float_t p0gsgs_PS_VAJHU
			 //const Float_t bkg_VAMCFMNorm,
			 //const Float_t p0_pt,
			 //const Float_t p0_y,
			 //const Float_t bkg_pt,
			 //const Float_t bkg_y
			 );
  void FillSuperMela(const Float_t p0plus_m4l,
		     const Float_t bkg_m4l,
		     const Float_t p0plus_m4l_ScaleUp,
		     const Float_t bkg_m4l_ScaleUp,
		     const Float_t p0plus_m4l_ScaleDown,
		     const Float_t bkg_m4l_ScaleDown,
		     const Float_t p0plus_m4l_ResUp,
		     const Float_t bkg_m4l_ResUp,
		     const Float_t p0plus_m4l_ResDown,
		     const Float_t bkg_m4l_ResDown
		     );
  void FillHAdditionalInfo(const Float_t mZa, const Float_t mZb, Float_t mLL4, const Float_t mLL6, const Float_t SIP4, const Float_t iso34);
  void FillZInfo(const Float_t ZMass, const Float_t ZPt, const Float_t Z1MassRefit);
  void FillAngularInfo(const Float_t costhetastar, const Float_t phi, const Float_t costheta1, const Float_t costheta2, const Float_t phistar1, const Float_t phistar2,const Float_t xi, const Float_t xistar);
  void FillLepInfo(const Float_t LepPt, const Float_t LepEta, const Float_t LepPhi, const Int_t LepId, const Float_t SIP, bool isID, float BDT, short parentId, int missingHit);
  void FillLepIsolInfo(const Float_t LepchargedHadIso, const Float_t LepneutralHadIso, const Float_t LepphotonIso, const Float_t LepcombRelIsoPF);
  void FillPhotonInfo(const Float_t PhotPt, const Float_t PhotEta, const Float_t PhotPhi);
  void FillJetInfo(const Float_t JetPt, const Float_t JetEta, const Float_t JetPhi, const Float_t JetMass, const Float_t JetBTag, const Float_t JetSigma);
  void FillDiJetInfo(const Float_t DiJetMass, const Float_t DiJetMassPlus, const Float_t DiJetMassMinus, const Float_t DiJetDEta);
  void FillCategorizationInfo(const Int_t nExtraLep, const Int_t nExtraZ, const Int_t nJets, const Int_t nCleanedJets, const Int_t nCleanedJetsPt30, const Int_t nCleanedJetsPt30BTagged);
  void FillExtraLepInfo(int extraLeptonIndex, bool extraLeptonExists, const reco::CandidatePtr ExtraLep);

  void FillHGenInfo(const math::XYZTLorentzVector Hp);
  void FillZGenInfo(const math::XYZTLorentzVector Z1p, const math::XYZTLorentzVector Z2p);
  void FillLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, Short_t Lep3Id, Short_t Lep4Id,
		      const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2, const math::XYZTLorentzVector Lep3, const math::XYZTLorentzVector Lep4);
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
  Int_t _genFinalState;
  Int_t _genProcessId;
  Float_t _genHEPMCweight;
  Short_t _trigWord;
  Short_t _genExtInfo;
  //H variables
  std::vector<Float_t> _ZZMass;
  std::vector<Float_t> _ZZMassPreFSR;
  std::vector<Float_t> _ZZMassRefit;
  std::vector<Float_t> _ZZMassErr;
  std::vector<Float_t> _ZZMassErrCorr;
  std::vector<Float_t> _Chi2KinFit;
  std::vector<Float_t> _ZZMassCFit;
  std::vector<Float_t> _Chi2CFit;
  std::vector<Int_t> _ZZsel;
  std::vector<Float_t> _ZZPt;
  std::vector<Float_t> _ZZEta;
  std::vector<Float_t> _ZZPhi;
  std::vector<Int_t> _ZZgenIsSignal;
  std::vector<Int_t> _ZZgenIsRightPair;
  std::vector<Float_t> _ZZFisher;
  std::vector<Int_t> _CRflag;

  //probabilities
  //std::vector<Float_t> _p0plus_melaNorm;// higgs, analytic distribution, normalized
  //std::vector<Float_t> _p0plus_mela;  // higgs, analytic distribution
  //std::vector<Float_t> _p0minus_mela; // pseudoscalar, analytic distribution
  //std::vector<Float_t> _p0hplus_mela;// 0h+, analytic distribution
  std::vector<Float_t> _p0plus_VAJHU; // higgs, vector algebra, JHUgen
  std::vector<Float_t> _p0minus_VAJHU;// pseudoscalar, vector algebra, JHUgen
  std::vector<Float_t> _p0plus_VAMCFM;// higgs, vector algebra, MCFM
  std::vector<Float_t> _p0hplus_VAJHU;// 0h+ (high dimensional operator), vector algebra, JHUgen
  //std::vector<Float_t> _p1_mela;    // zprime, analytic distribution
  //std::vector<Float_t> _p1_prodIndep_mela;    // zprime, analytic distribution
  //std::vector<Float_t> _p1plus_mela;// 1+, analytic distribution 
  //std::vector<Float_t> _p1plus_prodIndep_mela;// 1+, analytic distribution 
  std::vector<Float_t> _p1_VAJHU;   // zprime, vector algebra, JHUgen,
  std::vector<Float_t> _p1_prodIndep_VAJHU;   // zprime, vector algebra, JHUgen,
  std::vector<Float_t> _p1plus_VAJHU;// 1+ (axial vector), vector algebra, JHUgen,
  std::vector<Float_t> _p1plus_prodIndep_VAJHU;// 1+ (axial vector), vector algebra, JHUgen,
  //std::vector<Float_t> _p2_mela ;   // graviton, analytic distribution
  //std::vector<Float_t> _p2_prodIndep_mela ;   // graviton, analytic distribution
  //std::vector<Float_t> _p2qqb_mela;// graviton produced by qqbar vector algebra, analytical,
  //std::vector<Float_t> _p2hplus_mela;// graviton produced by qqbar vector algebra, analytical,
  //std::vector<Float_t> _p2hminus_mela;// graviton produced by qqbar vector algebra, analytical,
  //std::vector<Float_t> _p2bplus_mela;// graviton produced by qqbar vector algebra, analytical,
  std::vector<Float_t> _p2_VAJHU;   // graviton, vector algebra, JHUgen,
  std::vector<Float_t> _p2_prodIndep_VAJHU;   // graviton, vector algebra, JHUgen,
  std::vector<Float_t> _p2qqb_VAJHU;   // graviton, vector algebra, JHUgen,
  std::vector<Float_t> _p2hplus_VAJHU;   // graviton, vector algebra, JHUgen,
  std::vector<Float_t> _p2hminus_VAJHU;   // graviton, vector algebra, JHUgen,
  std::vector<Float_t> _p2bplus_VAJHU;   // graviton, vector algebra, JHUgen,
  std::vector<Float_t> _p2hplus_qqb_VAJHU;									  
  std::vector<Float_t> _p2hplus_prodIndep_VAJHU;		
  std::vector<Float_t> _p2hminus_qqb_VAJHU;				
  std::vector<Float_t> _p2hminus_prodIndep_VAJHU;	
  std::vector<Float_t> _p2bplus_qqb_VAJHU;					
  std::vector<Float_t> _p2bplus_prodIndep_VAJHU;		
  std::vector<Float_t> _p2h2plus_gg_VAJHU;      		
  std::vector<Float_t> _p2h2plus_qqbar_VAJHU;   		
  std::vector<Float_t> _p2h2plus_prodIndep_VAJHU;	
  std::vector<Float_t> _p2h3plus_gg_VAJHU;       	
  std::vector<Float_t> _p2h3plus_qqbar_VAJHU;    	
  std::vector<Float_t> _p2h3plus_prodIndep_VAJHU;	
  std::vector<Float_t> _p2h6plus_gg_VAJHU;       	
  std::vector<Float_t> _p2h6plus_qqbar_VAJHU;    	
  std::vector<Float_t> _p2h6plus_prodIndep_VAJHU;	
  std::vector<Float_t> _p2h7plus_gg_VAJHU;       	
  std::vector<Float_t> _p2h7plus_qqbar_VAJHU;    	
  std::vector<Float_t> _p2h7plus_prodIndep_VAJHU;	
  std::vector<Float_t> _p2h9minus_gg_VAJHU;       	
  std::vector<Float_t> _p2h9minus_qqbar_VAJHU;    	
  std::vector<Float_t> _p2h9minus_prodIndep_VAJHU;	
  std::vector<Float_t> _p2h10minus_gg_VAJHU;       
  std::vector<Float_t> _p2h10minus_qqbar_VAJHU;    
  std::vector<Float_t> _p2h10minus_prodIndep_VAJHU;
  //backgrounds
  //std::vector<Float_t> _bkg_mela;   // background,  analytic distribution
  std::vector<Float_t> _bkg_VAMCFM; // background, vector algebra, MCFM
  std::vector<Float_t> _bkg_prodIndep_VAMCFM; // background, vector algebra, MCFM
  std::vector<Float_t> _ggzz_VAMCFM; // background, vector algebra, MCFM for ggzz
  std::vector<Float_t> _ggzz_p0plus_VAMCFM; // background, vector algebra, MCFM for ggzz
  std::vector<Float_t> _ggzz_c1_VAMCFM; // signal + background + interference w/ SM couplings, vector algebra, MCFM
  std::vector<Float_t> _ggzz_c5_VAMCFM; // signal + background + interference w/ 5xSM couplings, vector algebra, MCFM for ggzz
  std::vector<Float_t> _ggzz_ci_VAMCFM; // signal + background + interference w/ imaginary SM couplings, vector algebra, MCFM for ggzz  
  //std::vector<Float_t> _bkg_VAMCFMNorm; // background, vector algebra, MCFM
  //pt/rapidity
  //std::vector<Float_t> _p0_pt;  // multiplicative probability for signal pt
  //std::vector<Float_t> _p0_y;   // multiplicative probability for signal y
  //std::vector<Float_t> _bkg_pt; // multiplicative probability for bkg pt
  //std::vector<Float_t> _bkg_y;  // multiplicative probability for bkg y
  // supermela
  std::vector<Float_t> _p0plus_m4l;  // signal m4l probability as in datacards
  std::vector<Float_t> _bkg_m4l;     // backgroun m4l probability as in datacards
  std::vector<Float_t> _pg1g4_mela;
  std::vector<Float_t> _pg1g4_VAJHU;
  std::vector<Float_t> _pg1g4_pi2_VAJHU;
  std::vector<Float_t> _pg1g2_pi2_VAJHU;
  std::vector<Float_t> _pg1g2_mela;
  std::vector<Float_t> _pg1g2_VAJHU;
  std::vector<Float_t> _p0plus_m4l_ScaleUp;// signal m4l probability for systematics
  std::vector<Float_t> _bkg_m4l_ScaleUp;// backgroun m4l probability for systematics
  std::vector<Float_t> _p0plus_m4l_ScaleDown;// signal m4l probability for systematics
  std::vector<Float_t> _bkg_m4l_ScaleDown;// backgroun m4l probability for systematics
  std::vector<Float_t> _p0plus_m4l_ResUp;// signal m4l probability for systematics
  std::vector<Float_t> _bkg_m4l_ResUp;// backgroun m4l probability for systematics
  std::vector<Float_t> _p0plus_m4l_ResDown;// signal m4l probability for systematics
  std::vector<Float_t> _bkg_m4l_ResDown;// backgroun m4l probability for systematics
  // Production MELA
  std::vector<Float_t> _phjj_VAJHU_old; //H+jj, vector algebra, JHUGen, Legacy Jet Selection
  std::vector<Float_t> _pvbf_VAJHU_old; //VBF, vector algebra, JHUGen, Legacy Jet Selection
  std::vector<Float_t> _phjj_VAJHU_old_up; //H+jj, vector algebra, JHUGen, Legacy Jet Selection, JEC +
  std::vector<Float_t> _pvbf_VAJHU_old_up; //VBF, vector algebra, JHUGen, Legacy Jet Selection, JEC +
  std::vector<Float_t> _phjj_VAJHU_old_dn; //H+jj, vector algebra, JHUGen, Legacy Jet Selection, JEC -
  std::vector<Float_t> _pvbf_VAJHU_old_dn; //VBF, vector algebra, JHUGen, Legacy Jet Selection, JEC -
  std::vector<Float_t> _phjj_VAJHU_new; //H+jj, vector algebra, JHUGen, New Jet Selection
  std::vector<Float_t> _pvbf_VAJHU_new; //VBF, vector algebra, JHUGen, New Jet Selection
  std::vector<Float_t> _phjj_VAJHU_new_up; //H+jj, vector algebra, JHUGen, New Jet Selection, JEC +
  std::vector<Float_t> _pvbf_VAJHU_new_up; //VBF, vector algebra, JHUGen, New Jet Selection, JEC +
  std::vector<Float_t> _phjj_VAJHU_new_dn; //H+jj, vector algebra, JHUGen, New Jet Selection, JEC -
  std::vector<Float_t> _pvbf_VAJHU_new_dn; //VBF, vector algebra, JHUGen, New Jet Selection, JEC -
  std::vector<Float_t> _p0_g1prime2_VAJHU; //VBF, vector algebra, JHUGen, New Jet Selection, JEC -
  std::vector<Float_t> _pg1g1prime2_VAJHU; //VBF, vector algebra, JHUGen, New Jet Selection, JEC -
  std::vector<Float_t> _Dgg10_VAMCFM; //VBF, vector algebra, JHUGen, New Jet Selection, JEC -

  std::vector<Float_t> _pzzzg_VAJHU;
  std::vector<Float_t> _pzzgg_VAJHU;
  std::vector<Float_t> _pzzzg_PS_VAJHU;
  std::vector<Float_t> _pzzgg_PS_VAJHU;
  std::vector<Float_t> _p0Zgs_VAJHU;
  std::vector<Float_t> _p0gsgs_VAJHU;
  std::vector<Float_t> _p0Zgs_PS_VAJHU;
  std::vector<Float_t> _p0gsgs_PS_VAJHU;

  std::vector<Float_t> _ZZmZa;
  std::vector<Float_t> _ZZmZb;
  std::vector<Float_t> _ZZmLL4;
  std::vector<Float_t> _ZZmLL6;
  std::vector<Float_t> _ZZSIP4;
  //  std::vector<Float_t> _ZZiso34;

  //Z1 variables
  std::vector<Float_t> _Z1Mass;
  std::vector<Float_t> _Z1MassRefit;
  std::vector<Float_t> _Z1Pt;

  //Z2 variables
  std::vector<Float_t> _Z2Mass;
  std::vector<Float_t> _Z2Pt;

  //Angular variables
  std::vector<Float_t> _costhetastar;
  std::vector<Float_t> _phi;
//   std::vector<Float_t> _helphiZ1;
//   std::vector<Float_t> _helphiZ2;
  std::vector<Float_t> _costheta1;
  std::vector<Float_t> _costheta2;
  std::vector<Float_t> _phistar1;
  std::vector<Float_t> _phistar2;
  std::vector<Float_t> _xi;
  std::vector<Float_t> _xistar;

  //Lepton variables
  std::vector<Float_t> _Lep1Pt;
  std::vector<Float_t> _Lep1Eta;
  std::vector<Float_t> _Lep1Phi;
  std::vector<Int_t>   _Lep1LepId;
  std::vector<Float_t> _Lep1SIP;
  std::vector<Bool_t>  _Lep1isID;
  std::vector<Float_t> _Lep1BDT;
  std::vector<Char_t>  _Lep1missingHit;
  std::vector<Short_t> _Lep1ParentId;

  std::vector<Float_t> _Lep2Pt;
  std::vector<Float_t> _Lep2Eta;
  std::vector<Float_t> _Lep2Phi;
  std::vector<Int_t>   _Lep2LepId;
  std::vector<Float_t> _Lep2SIP;
  std::vector<Bool_t>  _Lep2isID;
  std::vector<Float_t> _Lep2BDT;
  std::vector<Char_t>  _Lep2missingHit;
  std::vector<Short_t> _Lep2ParentId;

  std::vector<Float_t> _Lep3Pt;
  std::vector<Float_t> _Lep3Eta;
  std::vector<Float_t> _Lep3Phi;
  std::vector<Int_t>   _Lep3LepId;
  std::vector<Float_t> _Lep3SIP;
  std::vector<Bool_t>  _Lep3isID;
  std::vector<Float_t> _Lep3BDT;
  std::vector<Char_t>  _Lep3missingHit;
  std::vector<Short_t> _Lep3ParentId;

  std::vector<Float_t> _Lep4Pt;
  std::vector<Float_t> _Lep4Eta;
  std::vector<Float_t> _Lep4Phi;
  std::vector<Int_t>   _Lep4LepId;
  std::vector<Float_t> _Lep4SIP;
  std::vector<Bool_t>  _Lep4isID;
  std::vector<Float_t> _Lep4BDT;
  std::vector<Char_t>  _Lep4missingHit;
  std::vector<Short_t> _Lep4ParentId;


  //Lepton isolation variables
  std::vector<Float_t> _Lep1chargedHadIso;
  std::vector<Float_t> _Lep1neutralHadIso;
  std::vector<Float_t> _Lep1photonIso;
  std::vector<Float_t> _Lep1combRelIsoPF;

  std::vector<Float_t> _Lep2chargedHadIso;
  std::vector<Float_t> _Lep2neutralHadIso;
  std::vector<Float_t> _Lep2photonIso;
  std::vector<Float_t> _Lep2combRelIsoPF;

  std::vector<Float_t> _Lep3chargedHadIso;
  std::vector<Float_t> _Lep3neutralHadIso;
  std::vector<Float_t> _Lep3photonIso;
  std::vector<Float_t> _Lep3combRelIsoPF;

  std::vector<Float_t> _Lep4chargedHadIso;
  std::vector<Float_t> _Lep4neutralHadIso;
  std::vector<Float_t> _Lep4photonIso;
  std::vector<Float_t> _Lep4combRelIsoPF;

  //Photon variables
  std::vector<Float_t> _PhotPt;
  std::vector<Float_t> _PhotEta;
  std::vector<Float_t> _PhotPhi;

  //Jet variables
  std::vector<Float_t> _JetPt;
  std::vector<Float_t> _JetEta;
  std::vector<Float_t> _JetPhi;
  std::vector<Float_t> _JetMass;
  std::vector<Float_t> _JetBTag;
  std::vector<Float_t> _JetSigma;
  Float_t _DiJetMass;
  Float_t _DiJetMassPlus;
  Float_t _DiJetMassMinus;
  Float_t _DiJetDEta;

  //Categorization-related variables
  std::vector<Int_t> _nExtraLep;
  std::vector<Int_t> _nExtraZ;
  std::vector<Int_t> _nJets;
  std::vector<Int_t> _nCleanedJets;
  std::vector<Int_t> _nCleanedJetsPt30;
  std::vector<Int_t> _nCleanedJetsPt30BTagged;

  //Variables of extra leptons
  std::vector<Float_t> _ExtraLep1Pt;
  std::vector<Float_t> _ExtraLep1Eta;
  std::vector<Float_t> _ExtraLep1Phi;
  std::vector<Int_t>   _ExtraLep1LepId;
  std::vector<Float_t> _ExtraLep1SIP;
  std::vector<Bool_t>  _ExtraLep1isID;
  std::vector<Float_t> _ExtraLep1BDT;
  std::vector<Char_t>  _ExtraLep1missingHit;
  std::vector<Float_t> _ExtraLep1chargedHadIso;
  std::vector<Float_t> _ExtraLep1neutralHadIso;
  std::vector<Float_t> _ExtraLep1photonIso;
  std::vector<Float_t> _ExtraLep1combRelIsoPF;

  std::vector<Float_t> _ExtraLep2Pt;
  std::vector<Float_t> _ExtraLep2Eta;
  std::vector<Float_t> _ExtraLep2Phi;
  std::vector<Int_t>   _ExtraLep2LepId;
  std::vector<Float_t> _ExtraLep2SIP;
  std::vector<Bool_t>  _ExtraLep2isID;
  std::vector<Float_t> _ExtraLep2BDT;
  std::vector<Char_t>  _ExtraLep2missingHit;
  std::vector<Float_t> _ExtraLep2chargedHadIso;
  std::vector<Float_t> _ExtraLep2neutralHadIso;
  std::vector<Float_t> _ExtraLep2photonIso;
  std::vector<Float_t> _ExtraLep2combRelIsoPF;

  std::vector<Float_t> _ExtraLep3Pt;
  std::vector<Float_t> _ExtraLep3Eta;
  std::vector<Float_t> _ExtraLep3Phi;
  std::vector<Int_t>   _ExtraLep3LepId;
  std::vector<Float_t> _ExtraLep3SIP;
  std::vector<Bool_t>  _ExtraLep3isID;
  std::vector<Float_t> _ExtraLep3BDT;
  std::vector<Char_t>  _ExtraLep3missingHit;
  std::vector<Float_t> _ExtraLep3chargedHadIso;
  std::vector<Float_t> _ExtraLep3neutralHadIso;
  std::vector<Float_t> _ExtraLep3photonIso;
  std::vector<Float_t> _ExtraLep3combRelIsoPF;

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
