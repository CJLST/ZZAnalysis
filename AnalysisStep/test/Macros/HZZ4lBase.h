//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun 17 10:06:26 2013 by ROOT version 5.32/00
// from TTree candTree/Event Summary
// found on file: rootuples/130613/PRODFSR_8TeV/ZZ4lAnalysis_H125.root
//////////////////////////////////////////////////////////

#ifndef HZZ4lBase_h
#define HZZ4lBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <iostream>

// Fixed size dimensions of array or collections stored in the TTree if any.

class HZZ4lBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           RunNumber;
   Long64_t        EventNumber;
   Int_t           LumiNumber;
   Int_t           iBC;
   Int_t           Nvtx;
   Int_t           NObsInt;
   Float_t         NTrueInt;
   Float_t         PUWeight12;
   Float_t         PFMET;
   Int_t           genFinalState;
   Int_t           genProcessId;
   Float_t         genHEPMCweight;
   Short_t         trigWord;
   Short_t         genExtInfo;
   std::vector<float>   *ZZMass;
   std::vector<float>   *ZZMassErr;
   std::vector<float>   *ZZMassErrCorr;
   std::vector<float>   *ZZMassPreFSR;
   std::vector<int>     *ZZsel;
   std::vector<float>   *ZZPt;
   std::vector<float>   *ZZEta;
   std::vector<float>   *ZZPhi;
   std::vector<float>   *ZZFisher;
   std::vector<float>   *pg1g4_mela;
   std::vector<float>   *pg1g4_VAJHU;
   std::vector<float>   *pg1g4_pi2_VAJHU;
   std::vector<float>   *pg1g2_pi2_VAJHU;
   std::vector<float>   *pg1g2_mela;
   std::vector<float>   *pg1g2_VAJHU;
   std::vector<int>     *CRflag;
   std::vector<float>   *p0plus_VAJHU;
   std::vector<float>   *p0minus_VAJHU;
   std::vector<float>   *p0hplus_VAJHU;
   std::vector<float>   *p1_VAJHU;
   std::vector<float>   *p1_prodIndep_VAJHU;
   std::vector<float>   *p1plus_VAJHU;
   std::vector<float>   *p1plus_prodIndep_VAJHU;
   std::vector<float>   *p2_VAJHU;
   std::vector<float>   *p2_prodIndep_VAJHU;
   std::vector<float>   *p2qqb_VAJHU;
   std::vector<float>   *p2hplus_VAJHU;
   std::vector<float>   *p2hminus_VAJHU;
   std::vector<float>   *p2bplus_VAJHU;
   std::vector<float>   *p2hplus_qqb_VAJHU;
   std::vector<float>   *p2hplus_prodIndep_VAJHU;
   std::vector<float>   *p2hminus_qqb_VAJHU;
   std::vector<float>   *p2hminus_prodIndep_VAJHU;
   std::vector<float>   *p2bplus_qqb_VAJHU;
   std::vector<float>   *p2bplus_prodIndep_VAJHU;
   std::vector<float>   *p2h2plus_gg_VAJHU;
   std::vector<float>   *p2h2plus_qqbar_VAJHU;
   std::vector<float>   *p2h2plus_prodIndep_VAJHU;
   std::vector<float>   *p2h3plus_gg_VAJHU;
   std::vector<float>   *p2h3plus_qqbar_VAJHU;
   std::vector<float>   *p2h3plus_prodIndep_VAJHU;
   std::vector<float>   *p2h6plus_gg_VAJHU;
   std::vector<float>   *p2h6plus_qqbar_VAJHU;
   std::vector<float>   *p2h6plus_prodIndep_VAJHU;
   std::vector<float>   *p2h7plus_gg_VAJHU;
   std::vector<float>   *p2h7plus_qqbar_VAJHU;
   std::vector<float>   *p2h7plus_prodIndep_VAJHU;
   std::vector<float>   *p2h9minus_gg_VAJHU;
   std::vector<float>   *p2h9minus_qqbar_VAJHU;
   std::vector<float>   *p2h9minus_prodIndep_VAJHU;
   std::vector<float>   *p2h10minus_gg_VAJHU;
   std::vector<float>   *p2h10minus_qqbar_VAJHU;
   std::vector<float>   *p2h10minus_prodIndep_VAJHU;
   std::vector<float>   *bkg_VAMCFM;
   std::vector<float>   *bkg_prodIndep_VAMCFM;
   std::vector<float>   *ggzz_VAMCFM;
   std::vector<float>   *ggzz_p0plus_VAMCFM;
   std::vector<float>   *p0plus_VAMCFM;
   std::vector<float>   *ggzz_c1_VAMCFM;
   std::vector<float>   *ggzz_c5_VAMCFM;
   std::vector<float>   *ggzz_ci_VAMCFM;
   std::vector<float>   *phjj_VAJHU_old;
   std::vector<float>   *pvbf_VAJHU_old;
   std::vector<float>   *phjj_VAJHU_old_up;
   std::vector<float>   *pvbf_VAJHU_old_up;
   std::vector<float>   *phjj_VAJHU_old_dn;
   std::vector<float>   *pvbf_VAJHU_old_dn;
   std::vector<float>   *phjj_VAJHU_new;
   std::vector<float>   *pvbf_VAJHU_new;
   std::vector<float>   *phjj_VAJHU_new_up;
   std::vector<float>   *pvbf_VAJHU_new_up;
   std::vector<float>   *phjj_VAJHU_new_dn;
   std::vector<float>   *pvbf_VAJHU_new_dn;
   std::vector<float>   *p0_g1prime2_VAJHU;
   std::vector<float>   *pg1g1prime2_VAJHU;
   std::vector<float>   *Dgg10_VAMCFM;
   std::vector<float>   *pzzzg_VAJHU;
   std::vector<float>   *pzzgg_VAJHU;
   std::vector<float>   *pzzzg_PS_VAJHU;
   std::vector<float>   *pzzgg_PS_VAJHU;
   std::vector<float>   *p0Zgs_VAJHU;
   std::vector<float>   *p0gsgs_VAJHU;
   std::vector<float>   *p0Zgs_PS_VAJHU;
   std::vector<float>   *p0gsgs_PS_VAJHU;
   std::vector<float>   *p0plus_m4l;
   std::vector<float>   *bkg_m4l;
   std::vector<float>   *p0plus_m4l_ScaleUp;
   std::vector<float>   *bkg_m4l_ScaleUp;
   std::vector<float>   *p0plus_m4l_ScaleDown;
   std::vector<float>   *bkg_m4l_ScaleDown;
   std::vector<float>   *p0plus_m4l_ResUp;
   std::vector<float>   *bkg_m4l_ResUp;
   std::vector<float>   *p0plus_m4l_ResDown;
   std::vector<float>   *bkg_m4l_ResDown;
   std::vector<float>   *Z1Mass;
   std::vector<float>   *Z1Pt;
   std::vector<float>   *Z2Mass;
   std::vector<float>   *Z2Pt;
   std::vector<float>   *costhetastar;
   std::vector<float>   *helphi;
   std::vector<float>   *helcosthetaZ1;
   std::vector<float>   *helcosthetaZ2;
   std::vector<float>   *phistarZ1;
   std::vector<float>   *phistarZ2;
   std::vector<float>   *xi;
   std::vector<float>   *xistar;
   std::vector<float>   *Lep1Pt;
   std::vector<float>   *Lep1Eta;
   std::vector<float>   *Lep1Phi;
   std::vector<int>     *Lep1LepId;
   std::vector<float>   *Lep1SIP;
   std::vector<bool>    *Lep1isID;
   std::vector<float>   *Lep1BDT;
   std::vector<char>    *Lep1missingHit;
   std::vector<float>   *Lep2Pt;
   std::vector<float>   *Lep2Eta;
   std::vector<float>   *Lep2Phi;
   std::vector<int>     *Lep2LepId;
   std::vector<float>   *Lep2SIP;
   std::vector<bool>    *Lep2isID;
   std::vector<float>   *Lep2BDT;
   std::vector<char>    *Lep2missingHit;
   std::vector<float>   *Lep3Pt;
   std::vector<float>   *Lep3Eta;
   std::vector<float>   *Lep3Phi;
   std::vector<int>     *Lep3LepId;
   std::vector<float>   *Lep3SIP;
   std::vector<bool>    *Lep3isID;
   std::vector<float>   *Lep3BDT;
   std::vector<char>    *Lep3missingHit;
   std::vector<float>   *Lep4Pt;
   std::vector<float>   *Lep4Eta;
   std::vector<float>   *Lep4Phi;
   std::vector<int>     *Lep4LepId;
   std::vector<float>   *Lep4SIP;
   std::vector<bool>    *Lep4isID;
   std::vector<float>   *Lep4BDT;
   std::vector<char>    *Lep4missingHit;
   std::vector<float>   *Lep1chargedHadIso;
   std::vector<float>   *Lep1neutralHadIso;
   std::vector<float>   *Lep1photonIso;
   std::vector<float>   *Lep1combRelIsoPF;
   std::vector<float>   *Lep2chargedHadIso;
   std::vector<float>   *Lep2neutralHadIso;
   std::vector<float>   *Lep2photonIso;
   std::vector<float>   *Lep2combRelIsoPF;
   std::vector<float>   *Lep3chargedHadIso;
   std::vector<float>   *Lep3neutralHadIso;
   std::vector<float>   *Lep3photonIso;
   std::vector<float>   *Lep3combRelIsoPF;
   std::vector<float>   *Lep4chargedHadIso;
   std::vector<float>   *Lep4neutralHadIso;
   std::vector<float>   *Lep4photonIso;
   std::vector<float>   *Lep4combRelIsoPF;
   std::vector<float>   *JetPt;
   std::vector<float>   *JetEta;
   std::vector<float>   *JetPhi;
   std::vector<float>   *JetMass;
   std::vector<float>   *JetBTag;
   std::vector<float>   *JetSigma;
   std::vector<float>   *JetQG;
   std::vector<float>   *JetQGsmear;
   std::vector<float>   *JetQGL;
   Float_t         DiJetMass;
   Float_t         DiJetMassPlus;
   Float_t         DiJetMassMinus;
   Float_t         DiJetDEta;
   Float_t         GenHMass;
   Float_t         GenHPt;
   Float_t         GenZ1Mass;
   Float_t         GenZ1Pt;
   Float_t         GenZ2Mass;
   Float_t         GenZ2Pt;
   Float_t         GenLep1Pt;
   Float_t         GenLep1Eta;
   Float_t         GenLep1Phi;
   Short_t         GenLep1Id;
   Float_t         GenLep2Pt;
   Float_t         GenLep2Eta;
   Float_t         GenLep2Phi;
   Short_t         GenLep2Id;
   Float_t         GenLep3Pt;
   Float_t         GenLep3Eta;
   Float_t         GenLep3Phi;
   Short_t         GenLep3Id;
   Float_t         GenLep4Pt;
   Float_t         GenLep4Eta;
   Float_t         GenLep4Phi;
   Short_t         GenLep4Id;

   // List of branches
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_LumiNumber;   //!
   TBranch        *b_iBC;   //!
   TBranch        *b_Nvtx;   //!
   TBranch        *b_NObsInt;   //!
   TBranch        *b_NTrueInt;   //!
   TBranch        *b_PUWeight12;   //!
   TBranch        *b_PFMET;   //!
   TBranch        *b_genFinalState;   //!
   TBranch        *b_genProcessId;   //!
   TBranch        *b_genHEPMCweight; //!
   TBranch        *b_trigWord;   //!
   TBranch        *b_genExtInfo;   //!
   TBranch        *b_ZZMass;   //!
   TBranch        *b_ZZMassErr;   //!
   TBranch        *b_ZZMassErrCorr;   //!
   TBranch        *b_ZZMassPreFSR;   //!
   TBranch        *b_ZZsel;   //!
   TBranch        *b_ZZPt;   //!
   TBranch        *b_ZZEta;   //!
   TBranch        *b_ZZPhi;   //!
   TBranch        *b_ZZFisher;   //!
   TBranch        *b_pg1g4_mela; //!
   TBranch        *b_pg1g4_VAJHU; //!
   TBranch        *b_pg1g4_pi2_VAJHU; //!
   TBranch        *b_pg1g2_pi2_VAJHU; //!
   TBranch        *b_pg1g2_mela; //!
   TBranch        *b_pg1g2_VAJHU; //!
   TBranch        *b_CRflag;   //!
   TBranch        *b_p0plus_VAJHU;   //!
   TBranch        *b_p0minus_VAJHU;   //!
   TBranch        *b_p0hplus_VAJHU;   //!
   TBranch        *b_p1_VAJHU;   //!
   TBranch        *b_p1_prodIndep_VAJHU;   //!
   TBranch        *b_p1plus_VAJHU;   //!
   TBranch        *b_p1plus_prodIndep_VAJHU;   //!
   TBranch        *b_p2_VAJHU;   //!
   TBranch        *b_p2_prodIndep_VAJHU;   //!
   TBranch        *b_p2qqb_VAJHU;   //!
   TBranch        *b_p2hplus_VAJHU;   //!
   TBranch        *b_p2hminus_VAJHU;   //!
   TBranch        *b_p2bplus_VAJHU;   //!
   TBranch        *b_p2hplus_qqb_VAJHU;   //!
   TBranch        *b_p2hplus_prodIndep_VAJHU;   //!
   TBranch        *b_p2hminus_qqb_VAJHU;   //!
   TBranch        *b_p2hminus_prodIndep_VAJHU;   //!
   TBranch        *b_p2bplus_qqb_VAJHU;   //!
   TBranch        *b_p2bplus_prodIndep_VAJHU;   //!
   TBranch        *b_p2h2plus_gg_VAJHU;   //!
   TBranch        *b_p2h2plus_qqbar_VAJHU;   //!
   TBranch        *b_p2h2plus_prodIndep_VAJHU;   //!
   TBranch        *b_p2h3plus_gg_VAJHU;   //!
   TBranch        *b_p2h3plus_qqbar_VAJHU;   //!
   TBranch        *b_p2h3plus_prodIndep_VAJHU;   //!
   TBranch        *b_p2h6plus_gg_VAJHU;   //!
   TBranch        *b_p2h6plus_qqbar_VAJHU;   //!
   TBranch        *b_p2h6plus_prodIndep_VAJHU;   //!
   TBranch        *b_p2h7plus_gg_VAJHU;   //!
   TBranch        *b_p2h7plus_qqbar_VAJHU;   //!
   TBranch        *b_p2h7plus_prodIndep_VAJHU;   //!
   TBranch        *b_p2h9minus_gg_VAJHU;   //!
   TBranch        *b_p2h9minus_qqbar_VAJHU;   //!
   TBranch        *b_p2h9minus_prodIndep_VAJHU;   //!
   TBranch        *b_p2h10minus_gg_VAJHU;   //!
   TBranch        *b_p2h10minus_qqbar_VAJHU;   //!
   TBranch        *b_p2h10minus_prodIndep_VAJHU;   //!
   TBranch        *b_bkg_VAMCFM;   //!
   TBranch        *b_bkg_prodIndep_VAMCFM;   //!
   TBranch        *b_ggzz_VAMCFM;   //!
   TBranch        *b_ggzz_p0plus_VAMCFM;   //!
   TBranch        *b_p0plus_VAMCFM;   //!
   TBranch        *b_ggzz_c5_VAMCFM;   //!
   TBranch        *b_ggzz_c1_VAMCFM;   //!
   TBranch        *b_ggzz_ci_VAMCFM;   //!
   TBranch        *b_phjj_VAJHU_old;   //!
   TBranch        *b_pvbf_VAJHU_old;   //!
   TBranch        *b_phjj_VAJHU_old_up;   //!
   TBranch        *b_pvbf_VAJHU_old_up;   //!
   TBranch        *b_phjj_VAJHU_old_dn;   //!
   TBranch        *b_pvbf_VAJHU_old_dn;   //!
   TBranch        *b_phjj_VAJHU_new;   //!
   TBranch        *b_pvbf_VAJHU_new;   //!
   TBranch        *b_phjj_VAJHU_new_up;   //!
   TBranch        *b_pvbf_VAJHU_new_up;   //!
   TBranch        *b_phjj_VAJHU_new_dn;   //!
   TBranch        *b_pvbf_VAJHU_new_dn;   //!
   TBranch        *b_p0_g1prime2_VAJHU;   //!
   TBranch        *b_pg1g1prime2_VAJHU;   //!
   TBranch        *b_Dgg10_VAMCFM;   //!
   TBranch        *b_pzzzg_VAJHU;
   TBranch        *b_pzzgg_VAJHU;
   TBranch        *b_pzzzg_PS_VAJHU;
   TBranch        *b_pzzgg_PS_VAJHU;
   TBranch        *b_p0Zgs_VAJHU;
   TBranch        *b_p0gsgs_VAJHU;
   TBranch        *b_p0Zgs_PS_VAJHU;
   TBranch        *b_p0gsgs_PS_VAJHU;
   TBranch        *b_p0plus_m4l;   //!
   TBranch        *b_bkg_m4l;   //!
   TBranch        *b_p0plus_m4l_ScaleUp;   //!
   TBranch        *b_bkg_m4l_ScaleUp;   //!
   TBranch        *b_p0plus_m4l_ScaleDown;   //!
   TBranch        *b_bkg_m4l_ScaleDown;   //!
   TBranch        *b_p0plus_m4l_ResUp;   //!
   TBranch        *b_bkg_m4l_ResUp;   //!
   TBranch        *b_p0plus_m4l_ResDown;   //!
   TBranch        *b_bkg_m4l_ResDown;   //!
   TBranch        *b_Z1Mass;   //!
   TBranch        *b_Z1Pt;   //!
   TBranch        *b_Z2Mass;   //!
   TBranch        *b_Z2Pt;   //!
   TBranch        *b_costhetastar;   //!
   TBranch        *b_helphi;   //!
   TBranch        *b_helcosthetaZ1;   //!
   TBranch        *b_helcosthetaZ2;   //!
   TBranch        *b_phistarZ1;   //!
   TBranch        *b_phistarZ2;   //!
   TBranch        *b_xi;   //!
   TBranch        *b_xistar;   //!
   TBranch        *b_Lep1Pt;   //!
   TBranch        *b_Lep1Eta;   //!
   TBranch        *b_Lep1Phi;   //!
   TBranch        *b_Lep1LepId;   //!
   TBranch        *b_Lep1SIP;   //!
   TBranch        *b_Lep1isID;   //!
   TBranch        *b_Lep1BDT;   //!
   TBranch        *b_Lep1missingHit;   //!
   TBranch        *b_Lep2Pt;   //!
   TBranch        *b_Lep2Eta;   //!
   TBranch        *b_Lep2Phi;   //!
   TBranch        *b_Lep2LepId;   //!
   TBranch        *b_Lep2SIP;   //!
   TBranch        *b_Lep2isID;   //!
   TBranch        *b_Lep2BDT;   //!
   TBranch        *b_Lep2missingHit;   //!
   TBranch        *b_Lep3Pt;   //!
   TBranch        *b_Lep3Eta;   //!
   TBranch        *b_Lep3Phi;   //!
   TBranch        *b_Lep3LepId;   //!
   TBranch        *b_Lep3SIP;   //!
   TBranch        *b_Lep3isID;   //!
   TBranch        *b_Lep3BDT;   //!
   TBranch        *b_Lep3missingHit;   //!
   TBranch        *b_Lep4Pt;   //!
   TBranch        *b_Lep4Eta;   //!
   TBranch        *b_Lep4Phi;   //!
   TBranch        *b_Lep4LepId;   //!
   TBranch        *b_Lep4SIP;   //!
   TBranch        *b_Lep4isID;   //!
   TBranch        *b_Lep4BDT;   //!
   TBranch        *b_Lep4missingHit;   //!
   TBranch        *b_Lep1chargedHadIso;   //!
   TBranch        *b_Lep1neutralHadIso;   //!
   TBranch        *b_Lep1photonIso;   //!
   TBranch        *b_Lep1combRelIsoPF;   //!
   TBranch        *b_Lep2chargedHadIso;   //!
   TBranch        *b_Lep2neutralHadIso;   //!
   TBranch        *b_Lep2photonIso;   //!
   TBranch        *b_Lep2combRelIsoPF;   //!
   TBranch        *b_Lep3chargedHadIso;   //!
   TBranch        *b_Lep3neutralHadIso;   //!
   TBranch        *b_Lep3photonIso;   //!
   TBranch        *b_Lep3combRelIsoPF;   //!
   TBranch        *b_Lep4chargedHadIso;   //!
   TBranch        *b_Lep4neutralHadIso;   //!
   TBranch        *b_Lep4photonIso;   //!
   TBranch        *b_Lep4combRelIsoPF;   //!
   TBranch        *b_JetPt;   //!
   TBranch        *b_JetEta;   //!
   TBranch        *b_JetPhi;   //!
   TBranch        *b_JetMass;   //!
   TBranch        *b_JetBTag;   //!
   TBranch        *b_JetSigma;   //!
   TBranch        *b_JetQG;   //!
   TBranch        *b_JetQGsmear;   //!
   TBranch        *b_JetQGL;   //!
   TBranch        *b_DiJetMass;   //!
   TBranch        *b_DiJetMassPlus;   //!
   TBranch        *b_DiJetMassMinus;   //!
   TBranch        *b_DiJetDEta;   //!
   TBranch        *b_GenHMass;   //!
   TBranch        *b_GenHPt;   //!
   TBranch        *b_GenZ1Mass;   //!
   TBranch        *b_GenZ1Pt;   //!
   TBranch        *b_GenZ2Mass;   //!
   TBranch        *b_GenZ2Pt;   //!
   TBranch        *b_GenLep1Pt;   //!
   TBranch        *b_GenLep1Eta;   //!
   TBranch        *b_GenLep1Phi;   //!
   TBranch        *b_GenLep1Id;   //!
   TBranch        *b_GenLep2Pt;   //!
   TBranch        *b_GenLep2Eta;   //!
   TBranch        *b_GenLep2Phi;   //!
   TBranch        *b_GenLep2Id;   //!
   TBranch        *b_GenLep3Pt;   //!
   TBranch        *b_GenLep3Eta;   //!
   TBranch        *b_GenLep3Phi;   //!
   TBranch        *b_GenLep3Id;   //!
   TBranch        *b_GenLep4Pt;   //!
   TBranch        *b_GenLep4Eta;   //!
   TBranch        *b_GenLep4Phi;   //!
   TBranch        *b_GenLep4Id;   //!

   HZZ4lBase(TTree *tree=0);
   virtual ~HZZ4lBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef HZZ4lBase_cxx
HZZ4lBase::HZZ4lBase(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("rootuples/130613/PRODFSR_8TeV/ZZ4lAnalysis_H125.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("rootuples/130613/PRODFSR_8TeV/ZZ4lAnalysis_H125.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("rootuples/130613/PRODFSR_8TeV/ZZ4lAnalysis_H125.root:/ZZ4muTree");
      dir->GetObject("candTree",tree);

   }
   Init(tree);
}

HZZ4lBase::~HZZ4lBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t HZZ4lBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HZZ4lBase::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HZZ4lBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ZZMass = 0;
   ZZMassErr = 0;
   ZZMassErrCorr = 0;
   ZZMassPreFSR = 0;
   ZZsel = 0;
   ZZPt = 0;
   ZZEta = 0;
   ZZPhi = 0;
   ZZFisher = 0;
   pg1g4_mela=0; 
   pg1g4_VAJHU=0; 
   pg1g4_pi2_VAJHU=0; 
   pg1g2_pi2_VAJHU=0; 
   pg1g2_mela=0; 
   pg1g2_VAJHU=0; 
   CRflag = 0;
   p0plus_VAJHU = 0;
   p0minus_VAJHU = 0;
   p0plus_VAMCFM = 0;
   p0hplus_VAJHU = 0;
   p1_VAJHU = 0;
   p1_prodIndep_VAJHU = 0;
   p1plus_VAJHU = 0;
   p1plus_prodIndep_VAJHU = 0;
   p2_VAJHU = 0;
   p2_prodIndep_VAJHU = 0;
   p2qqb_VAJHU = 0;
   p2hplus_VAJHU = 0;
   p2hminus_VAJHU = 0;
   p2bplus_VAJHU = 0;
   p2hplus_qqb_VAJHU = 0;
   p2hplus_prodIndep_VAJHU = 0;
   p2hminus_qqb_VAJHU = 0;
   p2hminus_prodIndep_VAJHU = 0;
   p2bplus_qqb_VAJHU = 0;
   p2bplus_prodIndep_VAJHU = 0;
   p2h2plus_gg_VAJHU = 0;
   p2h2plus_qqbar_VAJHU = 0;
   p2h2plus_prodIndep_VAJHU = 0;
   p2h3plus_gg_VAJHU = 0;
   p2h3plus_qqbar_VAJHU = 0;
   p2h3plus_prodIndep_VAJHU = 0;
   p2h6plus_gg_VAJHU = 0;
   p2h6plus_qqbar_VAJHU = 0;
   p2h6plus_prodIndep_VAJHU = 0;
   p2h7plus_gg_VAJHU = 0;
   p2h7plus_qqbar_VAJHU = 0;
   p2h7plus_prodIndep_VAJHU = 0;
   p2h9minus_gg_VAJHU = 0;
   p2h9minus_qqbar_VAJHU = 0;
   p2h9minus_prodIndep_VAJHU = 0;
   p2h10minus_gg_VAJHU = 0;
   p2h10minus_qqbar_VAJHU = 0;
   p2h10minus_prodIndep_VAJHU = 0;
   bkg_VAMCFM = 0;
   bkg_prodIndep_VAMCFM = 0;
   ggzz_VAMCFM = 0;
   ggzz_p0plus_VAMCFM = 0;
   ggzz_c1_VAMCFM = 0;
   ggzz_c5_VAMCFM = 0;
   ggzz_ci_VAMCFM = 0;
   phjj_VAJHU_old = 0;
   pvbf_VAJHU_old = 0;
   phjj_VAJHU_old_up = 0;
   pvbf_VAJHU_old_up = 0;
   phjj_VAJHU_old_dn = 0;
   pvbf_VAJHU_old_dn = 0;
   phjj_VAJHU_new = 0;
   pvbf_VAJHU_new = 0;
   phjj_VAJHU_new_up = 0;
   pvbf_VAJHU_new_up = 0;
   phjj_VAJHU_new_dn = 0;
   pvbf_VAJHU_new_dn = 0;
   p0_g1prime2_VAJHU = 0;
   pg1g1prime2_VAJHU = 0;
   Dgg10_VAMCFM = 0;
   pzzzg_VAJHU=0;
   pzzgg_VAJHU=0;
   pzzzg_PS_VAJHU=0;
   pzzgg_PS_VAJHU=0;
   p0Zgs_VAJHU=0;
   p0gsgs_VAJHU=0;
   p0Zgs_PS_VAJHU=0;
   p0gsgs_PS_VAJHU=0;
   p0plus_m4l = 0;
   bkg_m4l = 0;
   p0plus_m4l_ScaleUp = 0;
   bkg_m4l_ScaleUp = 0;
   p0plus_m4l_ScaleDown = 0;
   bkg_m4l_ScaleDown = 0;
   p0plus_m4l_ResUp = 0;
   bkg_m4l_ResUp = 0;
   p0plus_m4l_ResDown = 0;
   bkg_m4l_ResDown = 0;
   Z1Mass = 0;
   Z1Pt = 0;
   Z2Mass = 0;
   Z2Pt = 0;
   costhetastar = 0;
   helphi = 0;
   helcosthetaZ1 = 0;
   helcosthetaZ2 = 0;
   phistarZ1 = 0;
   phistarZ2 = 0;
   xi = 0;
   xistar = 0;
   Lep1Pt = 0;
   Lep1Eta = 0;
   Lep1Phi = 0;
   Lep1LepId = 0;
   Lep1SIP = 0;
   Lep1isID = 0;
   Lep1BDT = 0;
   Lep1missingHit = 0;
   Lep2Pt = 0;
   Lep2Eta = 0;
   Lep2Phi = 0;
   Lep2LepId = 0;
   Lep2SIP = 0;
   Lep2isID = 0;
   Lep2BDT = 0;
   Lep2missingHit = 0;
   Lep3Pt = 0;
   Lep3Eta = 0;
   Lep3Phi = 0;
   Lep3LepId = 0;
   Lep3SIP = 0;
   Lep3isID = 0;
   Lep3BDT = 0;
   Lep3missingHit = 0;
   Lep4Pt = 0;
   Lep4Eta = 0;
   Lep4Phi = 0;
   Lep4LepId = 0;
   Lep4SIP = 0;
   Lep4isID = 0;
   Lep4BDT = 0;
   Lep4missingHit = 0;
   Lep1chargedHadIso = 0;
   Lep1neutralHadIso = 0;
   Lep1photonIso = 0;
   Lep1combRelIsoPF = 0;
   Lep2chargedHadIso = 0;
   Lep2neutralHadIso = 0;
   Lep2photonIso = 0;
   Lep2combRelIsoPF = 0;
   Lep3chargedHadIso = 0;
   Lep3neutralHadIso = 0;
   Lep3photonIso = 0;
   Lep3combRelIsoPF = 0;
   Lep4chargedHadIso = 0;
   Lep4neutralHadIso = 0;
   Lep4photonIso = 0;
   Lep4combRelIsoPF = 0;
   JetPt = 0;
   JetEta = 0;
   JetPhi = 0;
   JetMass = 0;
   JetBTag = 0;
   JetSigma = 0;
   JetQG = 0;
   JetQGsmear = 0;
   JetQGL = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("LumiNumber", &LumiNumber, &b_LumiNumber);
   fChain->SetBranchAddress("iBC", &iBC, &b_iBC);
   fChain->SetBranchAddress("Nvtx", &Nvtx, &b_Nvtx);
   fChain->SetBranchAddress("NObsInt", &NObsInt, &b_NObsInt);
   fChain->SetBranchAddress("NTrueInt", &NTrueInt, &b_NTrueInt);
   fChain->SetBranchAddress("PUWeight12", &PUWeight12, &b_PUWeight12);
   fChain->SetBranchAddress("PFMET", &PFMET, &b_PFMET);
   fChain->SetBranchAddress("genFinalState", &genFinalState, &b_genFinalState);
   fChain->SetBranchAddress("genProcessId", &genProcessId, &b_genProcessId);
   fChain->SetBranchAddress("genHEPMCweight",&genHEPMCweight,&b_genHEPMCweight);
   fChain->SetBranchAddress("trigWord", &trigWord, &b_trigWord);
   fChain->SetBranchAddress("genExtInfo", &genExtInfo, &b_genExtInfo);
   fChain->SetBranchAddress("ZZMass", &ZZMass, &b_ZZMass);
   fChain->SetBranchAddress("ZZMassErr", &ZZMassErr, &b_ZZMassErr);
   fChain->SetBranchAddress("ZZMassErrCorr", &ZZMassErrCorr, &b_ZZMassErrCorr);
   fChain->SetBranchAddress("ZZMassPreFSR", &ZZMassPreFSR, &b_ZZMassPreFSR);
   fChain->SetBranchAddress("ZZsel", &ZZsel, &b_ZZsel);
   fChain->SetBranchAddress("ZZPt", &ZZPt, &b_ZZPt);
   fChain->SetBranchAddress("ZZEta", &ZZEta, &b_ZZEta);
   fChain->SetBranchAddress("ZZPhi", &ZZPhi, &b_ZZPhi);
   fChain->SetBranchAddress("ZZFisher", &ZZFisher, &b_ZZFisher);
   fChain->SetBranchAddress("pg1g4_mela",&pg1g4_mela,&b_pg1g4_mela); 
   fChain->SetBranchAddress("pg1g4_VAJHU",&pg1g4_VAJHU,&b_pg1g4_VAJHU);
   fChain->SetBranchAddress("pg1g4_pi2_VAJHU",&pg1g4_pi2_VAJHU,&b_pg1g4_pi2_VAJHU); 
   fChain->SetBranchAddress("pg1g2_pi2_VAJHU",&pg1g2_pi2_VAJHU,&b_pg1g2_pi2_VAJHU); 
   fChain->SetBranchAddress("pg1g2_mela",&pg1g2_mela,&b_pg1g2_mela); 
   fChain->SetBranchAddress("pg1g2_VAJHU",&pg1g2_VAJHU,&b_pg1g2_VAJHU); 
   fChain->SetBranchAddress("CRflag", &CRflag, &b_CRflag);
   fChain->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU, &b_p0plus_VAJHU);
   fChain->SetBranchAddress("p0minus_VAJHU", &p0minus_VAJHU, &b_p0minus_VAJHU);
   fChain->SetBranchAddress("p0plus_VAMCFM", &p0plus_VAMCFM, &b_p0plus_VAMCFM);
   fChain->SetBranchAddress("p0hplus_VAJHU", &p0hplus_VAJHU, &b_p0hplus_VAJHU);
   fChain->SetBranchAddress("p1_VAJHU", &p1_VAJHU, &b_p1_VAJHU);
   fChain->SetBranchAddress("p1_prodIndep_VAJHU", &p1_prodIndep_VAJHU, &b_p1_prodIndep_VAJHU);
   fChain->SetBranchAddress("p1plus_VAJHU", &p1plus_VAJHU, &b_p1plus_VAJHU);
   fChain->SetBranchAddress("p1plus_prodIndep_VAJHU", &p1plus_prodIndep_VAJHU, &b_p1plus_prodIndep_VAJHU);
   fChain->SetBranchAddress("p2_VAJHU", &p2_VAJHU, &b_p2_VAJHU);
   fChain->SetBranchAddress("p2_prodIndep_VAJHU", &p2_prodIndep_VAJHU, &b_p2_prodIndep_VAJHU);
   fChain->SetBranchAddress("p2qqb_VAJHU", &p2qqb_VAJHU, &b_p2qqb_VAJHU);
   fChain->SetBranchAddress("p2hplus_VAJHU", &p2hplus_VAJHU, &b_p2hplus_VAJHU);
   fChain->SetBranchAddress("p2hminus_VAJHU", &p2hminus_VAJHU, &b_p2hminus_VAJHU);
   fChain->SetBranchAddress("p2bplus_VAJHU", &p2bplus_VAJHU, &b_p2bplus_VAJHU);
   fChain->SetBranchAddress("p2hplus_qqb_VAJHU", &p2hplus_qqb_VAJHU, &b_p2hplus_qqb_VAJHU);
   fChain->SetBranchAddress("p2hplus_prodIndep_VAJHU", &p2hplus_prodIndep_VAJHU, &b_p2hplus_prodIndep_VAJHU);
   fChain->SetBranchAddress("p2hminus_qqb_VAJHU", &p2hminus_qqb_VAJHU, &b_p2hminus_qqb_VAJHU);
   fChain->SetBranchAddress("p2hminus_prodIndep_VAJHU", &p2hminus_prodIndep_VAJHU, &b_p2hminus_prodIndep_VAJHU);
   fChain->SetBranchAddress("p2bplus_qqb_VAJHU", &p2bplus_qqb_VAJHU, &b_p2bplus_qqb_VAJHU);
   fChain->SetBranchAddress("p2bplus_prodIndep_VAJHU", &p2bplus_prodIndep_VAJHU, &b_p2bplus_prodIndep_VAJHU);
   fChain->SetBranchAddress("p2h2plus_gg_VAJHU", &p2h2plus_gg_VAJHU, &b_p2h2plus_gg_VAJHU);
   fChain->SetBranchAddress("p2h2plus_qqbar_VAJHU", &p2h2plus_qqbar_VAJHU, &b_p2h2plus_qqbar_VAJHU);
   fChain->SetBranchAddress("p2h2plus_prodIndep_VAJHU", &p2h2plus_prodIndep_VAJHU, &b_p2h2plus_prodIndep_VAJHU);
   fChain->SetBranchAddress("p2h3plus_gg_VAJHU", &p2h3plus_gg_VAJHU, &b_p2h3plus_gg_VAJHU);
   fChain->SetBranchAddress("p2h3plus_qqbar_VAJHU", &p2h3plus_qqbar_VAJHU, &b_p2h3plus_qqbar_VAJHU);
   fChain->SetBranchAddress("p2h3plus_prodIndep_VAJHU", &p2h3plus_prodIndep_VAJHU, &b_p2h3plus_prodIndep_VAJHU);
   fChain->SetBranchAddress("p2h6plus_gg_VAJHU", &p2h6plus_gg_VAJHU, &b_p2h6plus_gg_VAJHU);
   fChain->SetBranchAddress("p2h6plus_qqbar_VAJHU", &p2h6plus_qqbar_VAJHU, &b_p2h6plus_qqbar_VAJHU);
   fChain->SetBranchAddress("p2h6plus_prodIndep_VAJHU", &p2h6plus_prodIndep_VAJHU, &b_p2h6plus_prodIndep_VAJHU);
   fChain->SetBranchAddress("p2h7plus_gg_VAJHU", &p2h7plus_gg_VAJHU, &b_p2h7plus_gg_VAJHU);
   fChain->SetBranchAddress("p2h7plus_qqbar_VAJHU", &p2h7plus_qqbar_VAJHU, &b_p2h7plus_qqbar_VAJHU);
   fChain->SetBranchAddress("p2h7plus_prodIndep_VAJHU", &p2h7plus_prodIndep_VAJHU, &b_p2h7plus_prodIndep_VAJHU);
   fChain->SetBranchAddress("p2h9minus_gg_VAJHU", &p2h9minus_gg_VAJHU, &b_p2h9minus_gg_VAJHU);
   fChain->SetBranchAddress("p2h9minus_qqbar_VAJHU", &p2h9minus_qqbar_VAJHU, &b_p2h9minus_qqbar_VAJHU);
   fChain->SetBranchAddress("p2h9minus_prodIndep_VAJHU", &p2h9minus_prodIndep_VAJHU, &b_p2h9minus_prodIndep_VAJHU);
   fChain->SetBranchAddress("p2h10minus_gg_VAJHU", &p2h10minus_gg_VAJHU, &b_p2h10minus_gg_VAJHU);
   fChain->SetBranchAddress("p2h10minus_qqbar_VAJHU", &p2h10minus_qqbar_VAJHU, &b_p2h10minus_qqbar_VAJHU);
   fChain->SetBranchAddress("p2h10minus_prodIndep_VAJHU", &p2h10minus_prodIndep_VAJHU, &b_p2h10minus_prodIndep_VAJHU);
   fChain->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM, &b_bkg_VAMCFM);
   fChain->SetBranchAddress("bkg_prodIndep_VAMCFM", &bkg_prodIndep_VAMCFM, &b_bkg_prodIndep_VAMCFM);
   fChain->SetBranchAddress("ggzz_VAMCFM", &ggzz_VAMCFM, &b_ggzz_VAMCFM);
   fChain->SetBranchAddress("ggzz_p0plus_VAMCFM", &ggzz_p0plus_VAMCFM, &b_ggzz_p0plus_VAMCFM);
   fChain->SetBranchAddress("ggzz_c5_VAMCFM", &ggzz_c5_VAMCFM, &b_ggzz_c5_VAMCFM);
   fChain->SetBranchAddress("ggzz_c1_VAMCFM", &ggzz_c1_VAMCFM, &b_ggzz_c1_VAMCFM);
   fChain->SetBranchAddress("ggzz_ci_VAMCFM", &ggzz_ci_VAMCFM, &b_ggzz_ci_VAMCFM);
   fChain->SetBranchAddress("phjj_VAJHU_old", &phjj_VAJHU_old, &b_phjj_VAJHU_old);
   fChain->SetBranchAddress("pvbf_VAJHU_old", &pvbf_VAJHU_old, &b_pvbf_VAJHU_old);
   fChain->SetBranchAddress("phjj_VAJHU_old_up", &phjj_VAJHU_old_up, &b_phjj_VAJHU_old_up);
   fChain->SetBranchAddress("pvbf_VAJHU_old_up", &pvbf_VAJHU_old_up, &b_pvbf_VAJHU_old_up);
   fChain->SetBranchAddress("phjj_VAJHU_old_dn", &phjj_VAJHU_old_dn, &b_phjj_VAJHU_old_dn);
   fChain->SetBranchAddress("pvbf_VAJHU_old_dn", &pvbf_VAJHU_old_dn, &b_pvbf_VAJHU_old_dn);
   fChain->SetBranchAddress("phjj_VAJHU_new", &phjj_VAJHU_new, &b_phjj_VAJHU_new);
   fChain->SetBranchAddress("pvbf_VAJHU_new", &pvbf_VAJHU_new, &b_pvbf_VAJHU_new);
   fChain->SetBranchAddress("phjj_VAJHU_new_up", &phjj_VAJHU_new_up, &b_phjj_VAJHU_new_up);
   fChain->SetBranchAddress("pvbf_VAJHU_new_up", &pvbf_VAJHU_new_up, &b_pvbf_VAJHU_new_up);
   fChain->SetBranchAddress("phjj_VAJHU_new_dn", &phjj_VAJHU_new_dn, &b_phjj_VAJHU_new_dn);
   fChain->SetBranchAddress("pvbf_VAJHU_new_dn", &pvbf_VAJHU_new_dn, &b_pvbf_VAJHU_new_dn);
   fChain->SetBranchAddress("p0_g1prime2_VAJHU", &p0_g1prime2_VAJHU, &b_p0_g1prime2_VAJHU);
   fChain->SetBranchAddress("pg1g1prime2_VAJHU", &pg1g1prime2_VAJHU, &b_pg1g1prime2_VAJHU);
   fChain->SetBranchAddress("Dgg10_VAMCFM", &Dgg10_VAMCFM, &b_Dgg10_VAMCFM);
   if (fChain->GetBranchStatus("pzzzg_VAJHU")){
	   std::cout << "Found ZG and GG MEs!" << std::endl;

	   fChain->SetBranchAddress("pzzzg_VAJHU", &pzzzg_VAJHU, &b_pzzzg_VAJHU);
	   fChain->SetBranchAddress("pzzgg_VAJHU", &pzzgg_VAJHU, &b_pzzgg_VAJHU);
	   fChain->SetBranchAddress("pzzzg_PS_VAJHU", &pzzzg_PS_VAJHU, &b_pzzzg_PS_VAJHU);
	   fChain->SetBranchAddress("pzzgg_PS_VAJHU", &pzzgg_PS_VAJHU, &b_pzzgg_PS_VAJHU);
	   fChain->SetBranchAddress("p0Zgs_VAJHU", &p0Zgs_VAJHU, &b_p0Zgs_VAJHU);
	   fChain->SetBranchAddress("p0gsgs_VAJHU", &p0gsgs_VAJHU, &b_p0gsgs_VAJHU);
	   fChain->SetBranchAddress("p0Zgs_PS_VAJHU", &p0Zgs_PS_VAJHU, &b_p0Zgs_PS_VAJHU);
	   fChain->SetBranchAddress("p0gsgs_PS_VAJHU", &p0gsgs_PS_VAJHU, &b_p0gsgs_PS_VAJHU);
   };
   fChain->SetBranchAddress("p0plus_m4l", &p0plus_m4l, &b_p0plus_m4l);
   fChain->SetBranchAddress("bkg_m4l", &bkg_m4l, &b_bkg_m4l);
   fChain->SetBranchAddress("p0plus_m4l_ScaleUp", &p0plus_m4l_ScaleUp, &b_p0plus_m4l_ScaleUp);
   fChain->SetBranchAddress("bkg_m4l_ScaleUp", &bkg_m4l_ScaleUp, &b_bkg_m4l_ScaleUp);
   fChain->SetBranchAddress("p0plus_m4l_ScaleDown", &p0plus_m4l_ScaleDown, &b_p0plus_m4l_ScaleDown);
   fChain->SetBranchAddress("bkg_m4l_ScaleDown", &bkg_m4l_ScaleDown, &b_bkg_m4l_ScaleDown);
   fChain->SetBranchAddress("p0plus_m4l_ResUp", &p0plus_m4l_ResUp, &b_p0plus_m4l_ResUp);
   fChain->SetBranchAddress("bkg_m4l_ResUp", &bkg_m4l_ResUp, &b_bkg_m4l_ResUp);
   fChain->SetBranchAddress("p0plus_m4l_ResDown", &p0plus_m4l_ResDown, &b_p0plus_m4l_ResDown);
   fChain->SetBranchAddress("bkg_m4l_ResDown", &bkg_m4l_ResDown, &b_bkg_m4l_ResDown);
   fChain->SetBranchAddress("Z1Mass", &Z1Mass, &b_Z1Mass);
   fChain->SetBranchAddress("Z1Pt", &Z1Pt, &b_Z1Pt);
   fChain->SetBranchAddress("Z2Mass", &Z2Mass, &b_Z2Mass);
   fChain->SetBranchAddress("Z2Pt", &Z2Pt, &b_Z2Pt);
   fChain->SetBranchAddress("costhetastar", &costhetastar, &b_costhetastar);
   fChain->SetBranchAddress("helphi", &helphi, &b_helphi);
   fChain->SetBranchAddress("helcosthetaZ1", &helcosthetaZ1, &b_helcosthetaZ1);
   fChain->SetBranchAddress("helcosthetaZ2", &helcosthetaZ2, &b_helcosthetaZ2);
   fChain->SetBranchAddress("phistarZ1", &phistarZ1, &b_phistarZ1);
   fChain->SetBranchAddress("phistarZ2", &phistarZ2, &b_phistarZ2);
   fChain->SetBranchAddress("xi", &xi, &b_xi);
   fChain->SetBranchAddress("xistar", &xistar, &b_xistar);
   fChain->SetBranchAddress("Lep1Pt", &Lep1Pt, &b_Lep1Pt);
   fChain->SetBranchAddress("Lep1Eta", &Lep1Eta, &b_Lep1Eta);
   fChain->SetBranchAddress("Lep1Phi", &Lep1Phi, &b_Lep1Phi);
   fChain->SetBranchAddress("Lep1LepId", &Lep1LepId, &b_Lep1LepId);
   fChain->SetBranchAddress("Lep1SIP", &Lep1SIP, &b_Lep1SIP);
   fChain->SetBranchAddress("Lep1isID", &Lep1isID, &b_Lep1isID);
   fChain->SetBranchAddress("Lep1BDT", &Lep1BDT, &b_Lep1BDT);
   fChain->SetBranchAddress("Lep1missingHit", &Lep1missingHit, &b_Lep1missingHit);
   fChain->SetBranchAddress("Lep2Pt", &Lep2Pt, &b_Lep2Pt);
   fChain->SetBranchAddress("Lep2Eta", &Lep2Eta, &b_Lep2Eta);
   fChain->SetBranchAddress("Lep2Phi", &Lep2Phi, &b_Lep2Phi);
   fChain->SetBranchAddress("Lep2LepId", &Lep2LepId, &b_Lep2LepId);
   fChain->SetBranchAddress("Lep2SIP", &Lep2SIP, &b_Lep2SIP);
   fChain->SetBranchAddress("Lep2isID", &Lep2isID, &b_Lep2isID);
   fChain->SetBranchAddress("Lep2BDT", &Lep2BDT, &b_Lep2BDT);
   fChain->SetBranchAddress("Lep2missingHit", &Lep2missingHit, &b_Lep2missingHit);
   fChain->SetBranchAddress("Lep3Pt", &Lep3Pt, &b_Lep3Pt);
   fChain->SetBranchAddress("Lep3Eta", &Lep3Eta, &b_Lep3Eta);
   fChain->SetBranchAddress("Lep3Phi", &Lep3Phi, &b_Lep3Phi);
   fChain->SetBranchAddress("Lep3LepId", &Lep3LepId, &b_Lep3LepId);
   fChain->SetBranchAddress("Lep3SIP", &Lep3SIP, &b_Lep3SIP);
   fChain->SetBranchAddress("Lep3isID", &Lep3isID, &b_Lep3isID);
   fChain->SetBranchAddress("Lep3BDT", &Lep3BDT, &b_Lep3BDT);
   fChain->SetBranchAddress("Lep3missingHit", &Lep3missingHit, &b_Lep3missingHit);
   fChain->SetBranchAddress("Lep4Pt", &Lep4Pt, &b_Lep4Pt);
   fChain->SetBranchAddress("Lep4Eta", &Lep4Eta, &b_Lep4Eta);
   fChain->SetBranchAddress("Lep4Phi", &Lep4Phi, &b_Lep4Phi);
   fChain->SetBranchAddress("Lep4LepId", &Lep4LepId, &b_Lep4LepId);
   fChain->SetBranchAddress("Lep4SIP", &Lep4SIP, &b_Lep4SIP);
   fChain->SetBranchAddress("Lep4isID", &Lep4isID, &b_Lep4isID);
   fChain->SetBranchAddress("Lep4BDT", &Lep4BDT, &b_Lep4BDT);
   fChain->SetBranchAddress("Lep4missingHit", &Lep4missingHit, &b_Lep4missingHit);
   fChain->SetBranchAddress("Lep1chargedHadIso", &Lep1chargedHadIso, &b_Lep1chargedHadIso);
   fChain->SetBranchAddress("Lep1neutralHadIso", &Lep1neutralHadIso, &b_Lep1neutralHadIso);
   fChain->SetBranchAddress("Lep1photonIso", &Lep1photonIso, &b_Lep1photonIso);
   fChain->SetBranchAddress("Lep1combRelIsoPF", &Lep1combRelIsoPF, &b_Lep1combRelIsoPF);
   fChain->SetBranchAddress("Lep2chargedHadIso", &Lep2chargedHadIso, &b_Lep2chargedHadIso);
   fChain->SetBranchAddress("Lep2neutralHadIso", &Lep2neutralHadIso, &b_Lep2neutralHadIso);
   fChain->SetBranchAddress("Lep2photonIso", &Lep2photonIso, &b_Lep2photonIso);
   fChain->SetBranchAddress("Lep2combRelIsoPF", &Lep2combRelIsoPF, &b_Lep2combRelIsoPF);
   fChain->SetBranchAddress("Lep3chargedHadIso", &Lep3chargedHadIso, &b_Lep3chargedHadIso);
   fChain->SetBranchAddress("Lep3neutralHadIso", &Lep3neutralHadIso, &b_Lep3neutralHadIso);
   fChain->SetBranchAddress("Lep3photonIso", &Lep3photonIso, &b_Lep3photonIso);
   fChain->SetBranchAddress("Lep3combRelIsoPF", &Lep3combRelIsoPF, &b_Lep3combRelIsoPF);
   fChain->SetBranchAddress("Lep4chargedHadIso", &Lep4chargedHadIso, &b_Lep4chargedHadIso);
   fChain->SetBranchAddress("Lep4neutralHadIso", &Lep4neutralHadIso, &b_Lep4neutralHadIso);
   fChain->SetBranchAddress("Lep4photonIso", &Lep4photonIso, &b_Lep4photonIso);
   fChain->SetBranchAddress("Lep4combRelIsoPF", &Lep4combRelIsoPF, &b_Lep4combRelIsoPF);
   fChain->SetBranchAddress("JetPt", &JetPt, &b_JetPt);
   fChain->SetBranchAddress("JetEta", &JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetPhi", &JetPhi, &b_JetPhi);
   fChain->SetBranchAddress("JetMass", &JetMass, &b_JetMass);
   fChain->SetBranchAddress("JetBTag", &JetBTag, &b_JetBTag);
   fChain->SetBranchAddress("JetSigma", &JetSigma, &b_JetSigma);
   if (fChain->GetBranchStatus("JetQG")){
	   std::cout << "Found JetQG!" << std::endl;
	   fChain->SetBranchAddress("JetQG", &JetQG, &b_JetQG);
	   fChain->SetBranchAddress("JetQGsmear", &JetQGsmear, &b_JetQGsmear);
	   fChain->SetBranchAddress("JetQGL", &JetQGL, &b_JetQGL);
   };
   fChain->SetBranchAddress("DiJetMass", &DiJetMass, &b_DiJetMass);
   fChain->SetBranchAddress("DiJetMassPlus", &DiJetMassPlus, &b_DiJetMassPlus);
   fChain->SetBranchAddress("DiJetMassMinus", &DiJetMassMinus, &b_DiJetMassMinus);
   fChain->SetBranchAddress("DiJetDEta", &DiJetDEta, &b_DiJetDEta);
   fChain->SetBranchAddress("GenHMass", &GenHMass, &b_GenHMass);
   fChain->SetBranchAddress("GenHPt", &GenHPt, &b_GenHPt);
   fChain->SetBranchAddress("GenZ1Mass", &GenZ1Mass, &b_GenZ1Mass);
   fChain->SetBranchAddress("GenZ1Pt", &GenZ1Pt, &b_GenZ1Pt);
   fChain->SetBranchAddress("GenZ2Mass", &GenZ2Mass, &b_GenZ2Mass);
   fChain->SetBranchAddress("GenZ2Pt", &GenZ2Pt, &b_GenZ2Pt);
   fChain->SetBranchAddress("GenLep1Pt", &GenLep1Pt, &b_GenLep1Pt);
   fChain->SetBranchAddress("GenLep1Eta", &GenLep1Eta, &b_GenLep1Eta);
   fChain->SetBranchAddress("GenLep1Phi", &GenLep1Phi, &b_GenLep1Phi);
   fChain->SetBranchAddress("GenLep1Id", &GenLep1Id, &b_GenLep1Id);
   fChain->SetBranchAddress("GenLep2Pt", &GenLep2Pt, &b_GenLep2Pt);
   fChain->SetBranchAddress("GenLep2Eta", &GenLep2Eta, &b_GenLep2Eta);
   fChain->SetBranchAddress("GenLep2Phi", &GenLep2Phi, &b_GenLep2Phi);
   fChain->SetBranchAddress("GenLep2Id", &GenLep2Id, &b_GenLep2Id);
   fChain->SetBranchAddress("GenLep3Pt", &GenLep3Pt, &b_GenLep3Pt);
   fChain->SetBranchAddress("GenLep3Eta", &GenLep3Eta, &b_GenLep3Eta);
   fChain->SetBranchAddress("GenLep3Phi", &GenLep3Phi, &b_GenLep3Phi);
   fChain->SetBranchAddress("GenLep3Id", &GenLep3Id, &b_GenLep3Id);
   fChain->SetBranchAddress("GenLep4Pt", &GenLep4Pt, &b_GenLep4Pt);
   fChain->SetBranchAddress("GenLep4Eta", &GenLep4Eta, &b_GenLep4Eta);
   fChain->SetBranchAddress("GenLep4Phi", &GenLep4Phi, &b_GenLep4Phi);
   fChain->SetBranchAddress("GenLep4Id", &GenLep4Id, &b_GenLep4Id);
   Notify();
}

Bool_t HZZ4lBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HZZ4lBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HZZ4lBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef HZZ4lBase_cxx
