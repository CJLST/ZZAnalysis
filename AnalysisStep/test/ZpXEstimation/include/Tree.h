//////////////////////////////////////////////////////////
// This class has been automatically generated on
// $(date) by ROOT version 6.30/09
// from TTree candTree/Event Summary
// found on file: 2022EE samples
// Modified to work with 2022EE sample structure
//////////////////////////////////////////////////////////

#ifndef Tree_h
#define Tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

using namespace std;

class Tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types - branches that exist in 2022EE samples
   Int_t           RunNumber;
   Long64_t        EventNumber;
   Int_t           LumiNumber;
   Float_t         PFMET;
   Float_t         Z1Mass;
   Float_t         Z2Mass;
   Int_t           Z1Flav;
   Int_t           Z2Flav;
   Float_t         ZZMass;
   Int_t           CRflag;
   
   // Additional branches that may not exist in 2022EE samples but are used by the code
   Short_t         nCleanedJetsPt30;
   Short_t         nCleanedJetsPt30BTagged_bTagSF;
   Short_t         nExtraLep;
   Short_t         nExtraZ;
   Float_t         DiJetMass;
   Float_t         ZZPt;
   Float_t         ZZjjPt;
   Float_t         dataMCWeight;
   
   // MC-only branches
   Float_t         overallEventWeight;
   Float_t         KFactor_QCD_ggZZ_Nominal;
   Float_t         KFactor_EW_qqZZ;
   Float_t         KFactor_QCD_qqZZ_M;
   Int_t           L1prefiringWeight;
   Int_t           xsec;
   
   // JHUGen branches (may not exist in 2022EE samples)
   Float_t         p_GG_SIG_ghg2_1_ghz1_1_JHUGen;
   Float_t         p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen;
   Float_t         p_GG_SIG_ghg2_1_ghz2_1_JHUGen;
   Float_t         p_GG_SIG_ghg2_1_ghz4_1_JHUGen;
   Float_t         p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen;
   Float_t         p_GG_SIG_ghg2_1_ghza2_1_JHUGen;
   Float_t         p_GG_SIG_ghg2_1_ghza4_1_JHUGen;
   Float_t         p_GG_SIG_ghg2_1_gha2_1_JHUGen;
   Float_t         p_GG_SIG_ghg2_1_gha4_1_JHUGen;
   Float_t         p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen;
   Float_t         p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen;
   Float_t         p_GG_SIG_ghg2_1_ghz1_1_ghz2_i_JHUGen;
   Float_t         p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen;
   Float_t         p_GG_SIG_ghg2_1_ghz1_1_ghz4_i_JHUGen;
   Float_t         p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen;
   Float_t         p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4i_JHUGen;
   Float_t         p_GG_SIG_ghg2_1_ghz1_1_ghza2_1_JHUGen;
   Float_t         p_GG_SIG_ghg2_1_ghz1_1_ghza4_1_JHUGen;
   Float_t         p_GG_SIG_ghg2_1_ghz1_1_gha2_1_JHUGen;
   Float_t         p_GG_SIG_ghg2_1_ghz1_1_gha4_1_JHUGen;
   Float_t         pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
   Float_t         p_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
   Float_t         p_JQCD_SIG_ghg2_1_JHUGen_JECNominal;
   Float_t         p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal;
   Float_t         p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_JECNominal;
   Float_t         p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal;
   Float_t         p_HadWH_SIG_ghw1_1_JHUGen_JECNominal;
   Float_t         p_HadZH_SIG_ghz1_1_JHUGen_JECNominal;
	Float_t         p_HadWH_mavjj_JECNominal;
	Float_t         p_HadWH_mavjj_true_JECNominal;
   Float_t         p_HadZH_mavjj_JECNominal;
   Float_t         p_HadZH_mavjj_true_JECNominal;
   // Vector branches that exist in 2022EE samples
   vector<float>   *LepPt;
   vector<float>   *LepEta;
   vector<float>   *LepPhi;
   vector<short>   *LepLepId;
   vector<float>   *LepSIP;
   vector<float>   *Lepdxy;
   vector<float>   *Lepdz;
   vector<bool>    *LepisID;
   vector<unsigned char> *LepMissingHit;
   vector<float>   *LepCombRelIsoPF;
   vector<bool>    *Muon_ZZFullSel;
   
   // Jet branches (may not exist in 2022EE samples)
   vector<float>   *JetPt;
   vector<float>   *JetEta;
   vector<float>   *JetPhi;
   vector<float>   *JetMass;
   vector<float>   *JetQGLikelihood;

   // List of branches
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_LumiNumber;   //!
   TBranch        *b_PFMET;   //!
   TBranch        *b_Z1Mass;   //!
   TBranch        *b_Z2Mass;   //!
   TBranch        *b_Z1Flav;   //!
   TBranch        *b_Z2Flav;   //!
   TBranch        *b_ZZMass;   //!
   TBranch        *b_CRflag;   //!
   TBranch        *b_LepPt;   //!
   TBranch        *b_LepEta;   //!
   TBranch        *b_LepPhi;   //!
   TBranch        *b_LepLepId;   //!
   TBranch        *b_LepSIP;   //!
   TBranch        *b_Lepdxy;   //!
   TBranch        *b_Lepdz;   //!
   TBranch        *b_LepisID;   //!
   TBranch        *b_LepMissingHit;   //!
   TBranch        *b_LepCombRelIsoPF;   //!
   TBranch        *b_Muon_ZZFullSel;   //!
   
   // Additional branch pointers
   TBranch        *b_nCleanedJetsPt30;   //!
   TBranch        *b_nCleanedJetsPt30BTagged_bTagSF;   //!
   TBranch        *b_nExtraLep;   //!
   TBranch        *b_nExtraZ;   //!
   TBranch        *b_DiJetMass;   //!
   TBranch        *b_ZZPt;   //!
   TBranch        *b_ZZjjPt;   //!
   TBranch        *b_dataMCWeight;   //!
   TBranch        *b_JetPt;   //!
   TBranch        *b_JetEta;   //!
   TBranch        *b_JetPhi;   //!
   TBranch        *b_JetMass;   //!
   TBranch        *b_JetQGLikelihood;   //!
   
   // MC-only branch pointers
   TBranch        *b_overallEventWeight;   //!
   TBranch        *b_KFactor_QCD_ggZZ_Nominal;   //!
   TBranch        *b_KFactor_EW_qqZZ;   //!
   TBranch        *b_KFactor_QCD_qqZZ_M;   //!
   TBranch        *b_L1prefiringWeight;   //!
   TBranch        *b_xsec;   //!
   
   // JHUGen branch pointers
   TBranch        *b_p_GG_SIG_ghg2_1_ghz1_1_JHUGen;   //!
   TBranch        *b_p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen;   //!
   TBranch        *b_p_GG_SIG_ghg2_1_ghz2_1_JHUGen;   //!
   TBranch        *b_p_GG_SIG_ghg2_1_ghz4_1_JHUGen;   //!
   TBranch        *b_p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen;   //!
   TBranch        *b_p_GG_SIG_ghg2_1_ghza2_1_JHUGen;   //!
   TBranch        *b_p_GG_SIG_ghg2_1_ghza4_1_JHUGen;   //!
   TBranch        *b_p_GG_SIG_ghg2_1_gha2_1_JHUGen;   //!
   TBranch        *b_p_GG_SIG_ghg2_1_gha4_1_JHUGen;   //!
   TBranch        *b_p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen;   //!
   TBranch        *b_p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen;   //!
   TBranch        *b_p_GG_SIG_ghg2_1_ghz1_1_ghz2_i_JHUGen;   //!
   TBranch        *b_p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen;   //!
   TBranch        *b_p_GG_SIG_ghg2_1_ghz1_1_ghz4_i_JHUGen;   //!
   TBranch        *b_p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen;   //!
   TBranch        *b_p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4i_JHUGen;   //!
   TBranch        *b_p_GG_SIG_ghg2_1_ghz1_1_ghza2_1_JHUGen;   //!
   TBranch        *b_p_GG_SIG_ghg2_1_ghz1_1_ghza4_1_JHUGen;   //!
   TBranch        *b_p_GG_SIG_ghg2_1_ghz1_1_gha2_1_JHUGen;   //!
   TBranch        *b_p_GG_SIG_ghg2_1_ghz1_1_gha4_1_JHUGen;   //!
   TBranch        *b_pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal;   //!
   TBranch        *b_p_JVBF_SIG_ghv1_1_JHUGen_JECNominal;   //!
   TBranch        *b_p_JQCD_SIG_ghg2_1_JHUGen_JECNominal;   //!
   TBranch        *b_p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal;   //!
   TBranch        *b_p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_JECNominal;   //!
   TBranch        *b_p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal;   //!
   TBranch        *b_p_HadWH_SIG_ghw1_1_JHUGen_JECNominal;   //!
   TBranch        *b_p_HadZH_SIG_ghz1_1_JHUGen_JECNominal;   //!
   TBranch        *b_p_HadWH_mavjj_JECNominal;   //!
   TBranch        *b_p_HadWH_mavjj_true_JECNominal;   //!
   TBranch        *b_p_HadZH_mavjj_JECNominal;   //!
   TBranch        *b_p_HadZH_mavjj_true_JECNominal;   //!

   Tree(TTree *tree=0);
   virtual ~Tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, TString input_file_name, bool notZLregion);
   void             SetBranchAddressSafe(const char* branchName, void* addr, TBranch** branch);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Tree_cxx
Tree::Tree(TTree *tree) : fChain(0) 
{
}

Tree::~Tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Tree::LoadTree(Long64_t entry)
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

// Helper function to safely set branch address
void Tree::SetBranchAddressSafe(const char* branchName, void* addr, TBranch** branch) {
   if (fChain->GetBranch(branchName)) {
      fChain->SetBranchAddress(branchName, addr, branch);
   }
}

void Tree::Init(TTree *tree, TString input_file_name, bool notZLregion)
{
   // Set object pointer
   LepPt = 0;
   LepEta = 0;
   LepPhi = 0;
   LepLepId = 0;
   LepSIP = 0;
   Lepdxy = 0;
   Lepdz = 0;
   LepisID = 0;
   LepMissingHit = 0;
   LepCombRelIsoPF = 0;
   Muon_ZZFullSel = 0;
   JetPt = 0;
   JetEta = 0;
   JetPhi = 0;
   JetMass = 0;
   JetQGLikelihood = 0;
   
   // Initialize scalar variables to safe default values
   RunNumber = 0;
   EventNumber = 0;
   LumiNumber = 0;
   PFMET = 0.0;
   Z1Mass = 0.0;
   Z2Mass = 0.0;
   Z1Flav = 0;
   Z2Flav = 0;
   ZZMass = 0.0;
   CRflag = 0;
   nCleanedJetsPt30 = 0;
   nCleanedJetsPt30BTagged_bTagSF = 0;
   nExtraLep = 0;
   nExtraZ = 0;
   DiJetMass = 0.0;
   ZZPt = 0.0;
   ZZjjPt = 0.0;
   dataMCWeight = 1.0;
   overallEventWeight = 1.0;
   KFactor_QCD_ggZZ_Nominal = 1.0;
   KFactor_EW_qqZZ = 1.0;
   KFactor_QCD_qqZZ_M = 1.0;
   L1prefiringWeight = 1.0;
   xsec = 1.0;
   
   // Initialize JHUGen variables to safe default values
   p_GG_SIG_ghg2_1_ghz1_1_JHUGen = 0.0;
   p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen = 0.0;
   p_GG_SIG_ghg2_1_ghz2_1_JHUGen = 0.0;
   p_GG_SIG_ghg2_1_ghz4_1_JHUGen = 0.0;
   p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen = 0.0;
   p_GG_SIG_ghg2_1_ghza2_1_JHUGen = 0.0;
   p_GG_SIG_ghg2_1_ghza4_1_JHUGen = 0.0;
   p_GG_SIG_ghg2_1_gha2_1_JHUGen = 0.0;
   p_GG_SIG_ghg2_1_gha4_1_JHUGen = 0.0;
   p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen = 0.0;
   p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen = 0.0;
   p_GG_SIG_ghg2_1_ghz1_1_ghz2_i_JHUGen = 0.0;
   p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen = 0.0;
   p_GG_SIG_ghg2_1_ghz1_1_ghz4_i_JHUGen = 0.0;
   p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen = 0.0;
   p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4i_JHUGen = 0.0;
   p_GG_SIG_ghg2_1_ghz1_1_ghza2_1_JHUGen = 0.0;
   p_GG_SIG_ghg2_1_ghz1_1_ghza4_1_JHUGen = 0.0;
   p_GG_SIG_ghg2_1_ghz1_1_gha2_1_JHUGen = 0.0;
   p_GG_SIG_ghg2_1_ghz1_1_gha4_1_JHUGen = 0.0;
   pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal = 0.0;
   p_JVBF_SIG_ghv1_1_JHUGen_JECNominal = 0.0;
   p_JQCD_SIG_ghg2_1_JHUGen_JECNominal = 0.0;
   p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal = 0.0;
   p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_JECNominal = 0.0;
   p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal = 0.0;
   p_HadWH_SIG_ghw1_1_JHUGen_JECNominal = 0.0;
   p_HadZH_SIG_ghz1_1_JHUGen_JECNominal = 0.0;
   p_HadWH_mavjj_JECNominal = 0.0;
   p_HadWH_mavjj_true_JECNominal = 0.0;
   p_HadZH_mavjj_JECNominal = 0.0;
   p_HadZH_mavjj_true_JECNominal = 0.0;
   
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   // Set branches that exist in 2022EE samples
   SetBranchAddressSafe("RunNumber", &RunNumber, &b_RunNumber);
   SetBranchAddressSafe("EventNumber", &EventNumber, &b_EventNumber);
   SetBranchAddressSafe("LumiNumber", &LumiNumber, &b_LumiNumber);
   SetBranchAddressSafe("PFMET", &PFMET, &b_PFMET);
   SetBranchAddressSafe("Z1Mass", &Z1Mass, &b_Z1Mass);
   SetBranchAddressSafe("Z2Mass", &Z2Mass, &b_Z2Mass);
   SetBranchAddressSafe("Z1Flav", &Z1Flav, &b_Z1Flav);
   SetBranchAddressSafe("Z2Flav", &Z2Flav, &b_Z2Flav);
   SetBranchAddressSafe("ZZMass", &ZZMass, &b_ZZMass);
   SetBranchAddressSafe("CRflag", &CRflag, &b_CRflag);
   SetBranchAddressSafe("LepPt", &LepPt, &b_LepPt);
   SetBranchAddressSafe("LepEta", &LepEta, &b_LepEta);
   SetBranchAddressSafe("LepPhi", &LepPhi, &b_LepPhi);
   SetBranchAddressSafe("LepLepId", &LepLepId, &b_LepLepId);
   SetBranchAddressSafe("LepSIP", &LepSIP, &b_LepSIP);
   SetBranchAddressSafe("Lepdxy", &Lepdxy, &b_Lepdxy);
   SetBranchAddressSafe("Lepdz", &Lepdz, &b_Lepdz);
   SetBranchAddressSafe("LepisID", &LepisID, &b_LepisID);
   SetBranchAddressSafe("LepMissingHit", &LepMissingHit, &b_LepMissingHit);
   SetBranchAddressSafe("LepCombRelIsoPF", &LepCombRelIsoPF, &b_LepCombRelIsoPF);
   SetBranchAddressSafe("Muon_ZZFullSel", &Muon_ZZFullSel, &b_Muon_ZZFullSel);
   
   // Additional branches that may not exist in 2022EE samples
   SetBranchAddressSafe("nCleanedJetsPt30", &nCleanedJetsPt30, &b_nCleanedJetsPt30);
   SetBranchAddressSafe("nCleanedJetsPt30BTagged_bTagSF", &nCleanedJetsPt30BTagged_bTagSF, &b_nCleanedJetsPt30BTagged_bTagSF);
   SetBranchAddressSafe("nExtraLep", &nExtraLep, &b_nExtraLep);
   SetBranchAddressSafe("nExtraZ", &nExtraZ, &b_nExtraZ);
   SetBranchAddressSafe("DiJetMass", &DiJetMass, &b_DiJetMass);
   SetBranchAddressSafe("ZZPt", &ZZPt, &b_ZZPt);
   SetBranchAddressSafe("ZZjjPt", &ZZjjPt, &b_ZZjjPt);
   SetBranchAddressSafe("dataMCWeight", &dataMCWeight, &b_dataMCWeight);
   SetBranchAddressSafe("JetPt", &JetPt, &b_JetPt);
   SetBranchAddressSafe("JetEta", &JetEta, &b_JetEta);
   SetBranchAddressSafe("JetPhi", &JetPhi, &b_JetPhi);
   SetBranchAddressSafe("JetMass", &JetMass, &b_JetMass);
   SetBranchAddressSafe("JetQGLikelihood", &JetQGLikelihood, &b_JetQGLikelihood);
   
   // MC-only branches - these will only be set if they exist
   SetBranchAddressSafe("overallEventWeight", &overallEventWeight, &b_overallEventWeight);
   SetBranchAddressSafe("KFactor_QCD_ggZZ_Nominal", &KFactor_QCD_ggZZ_Nominal, &b_KFactor_QCD_ggZZ_Nominal);
   SetBranchAddressSafe("KFactor_EW_qqZZ", &KFactor_EW_qqZZ, &b_KFactor_EW_qqZZ);
   SetBranchAddressSafe("KFactor_QCD_qqZZ_M", &KFactor_QCD_qqZZ_M, &b_KFactor_QCD_qqZZ_M);
   SetBranchAddressSafe("L1prefiringWeight", &L1prefiringWeight, &b_L1prefiringWeight);
   SetBranchAddressSafe("xsec", &xsec, &b_xsec);
   
   // JHUGen branches - these will only be set if they exist
   SetBranchAddressSafe("p_GG_SIG_ghg2_1_ghz1_1_JHUGen", &p_GG_SIG_ghg2_1_ghz1_1_JHUGen, &b_p_GG_SIG_ghg2_1_ghz1_1_JHUGen);
   SetBranchAddressSafe("p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen", &p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen, &b_p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen);
   SetBranchAddressSafe("p_GG_SIG_ghg2_1_ghz2_1_JHUGen", &p_GG_SIG_ghg2_1_ghz2_1_JHUGen, &b_p_GG_SIG_ghg2_1_ghz2_1_JHUGen);
   SetBranchAddressSafe("p_GG_SIG_ghg2_1_ghz4_1_JHUGen", &p_GG_SIG_ghg2_1_ghz4_1_JHUGen, &b_p_GG_SIG_ghg2_1_ghz4_1_JHUGen);
   SetBranchAddressSafe("p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen", &p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen, &b_p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen);
   SetBranchAddressSafe("p_GG_SIG_ghg2_1_ghza2_1_JHUGen", &p_GG_SIG_ghg2_1_ghza2_1_JHUGen, &b_p_GG_SIG_ghg2_1_ghza2_1_JHUGen);
   SetBranchAddressSafe("p_GG_SIG_ghg2_1_ghza4_1_JHUGen", &p_GG_SIG_ghg2_1_ghza4_1_JHUGen, &b_p_GG_SIG_ghg2_1_ghza4_1_JHUGen);
   SetBranchAddressSafe("p_GG_SIG_ghg2_1_gha2_1_JHUGen", &p_GG_SIG_ghg2_1_gha2_1_JHUGen, &b_p_GG_SIG_ghg2_1_gha2_1_JHUGen);
   SetBranchAddressSafe("p_GG_SIG_ghg2_1_gha4_1_JHUGen", &p_GG_SIG_ghg2_1_gha4_1_JHUGen, &b_p_GG_SIG_ghg2_1_gha4_1_JHUGen);
   SetBranchAddressSafe("p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen", &p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen, &b_p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen);
   SetBranchAddressSafe("p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen", &p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen, &b_p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen);
   SetBranchAddressSafe("p_GG_SIG_ghg2_1_ghz1_1_ghz2_i_JHUGen", &p_GG_SIG_ghg2_1_ghz1_1_ghz2_i_JHUGen, &b_p_GG_SIG_ghg2_1_ghz1_1_ghz2_i_JHUGen);
   SetBranchAddressSafe("p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen", &p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen, &b_p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen);
   SetBranchAddressSafe("p_GG_SIG_ghg2_1_ghz1_1_ghz4_i_JHUGen", &p_GG_SIG_ghg2_1_ghz1_1_ghz4_i_JHUGen, &b_p_GG_SIG_ghg2_1_ghz1_1_ghz4_i_JHUGen);
   SetBranchAddressSafe("p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen", &p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen, &b_p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen);
   SetBranchAddressSafe("p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4i_JHUGen", &p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4i_JHUGen, &b_p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4i_JHUGen);
   SetBranchAddressSafe("p_GG_SIG_ghg2_1_ghz1_1_ghza2_1_JHUGen", &p_GG_SIG_ghg2_1_ghz1_1_ghza2_1_JHUGen, &b_p_GG_SIG_ghg2_1_ghz1_1_ghza2_1_JHUGen);
   SetBranchAddressSafe("p_GG_SIG_ghg2_1_ghz1_1_ghza4_1_JHUGen", &p_GG_SIG_ghg2_1_ghz1_1_ghza4_1_JHUGen, &b_p_GG_SIG_ghg2_1_ghz1_1_ghza4_1_JHUGen);
   SetBranchAddressSafe("p_GG_SIG_ghg2_1_ghz1_1_gha2_1_JHUGen", &p_GG_SIG_ghg2_1_ghz1_1_gha2_1_JHUGen, &b_p_GG_SIG_ghg2_1_ghz1_1_gha2_1_JHUGen);
   SetBranchAddressSafe("p_GG_SIG_ghg2_1_ghz1_1_gha4_1_JHUGen", &p_GG_SIG_ghg2_1_ghz1_1_gha4_1_JHUGen, &b_p_GG_SIG_ghg2_1_ghz1_1_gha4_1_JHUGen);
   SetBranchAddressSafe("pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal", &pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, &b_pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
   SetBranchAddressSafe("p_JVBF_SIG_ghv1_1_JHUGen_JECNominal", &p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, &b_p_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
   SetBranchAddressSafe("p_JQCD_SIG_ghg2_1_JHUGen_JECNominal", &p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, &b_p_JQCD_SIG_ghg2_1_JHUGen_JECNominal);
   SetBranchAddressSafe("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal", &p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, &b_p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal);
   SetBranchAddressSafe("p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_JECNominal", &p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_JECNominal, &b_p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_JECNominal);
   SetBranchAddressSafe("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal", &p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, &b_p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal);
   SetBranchAddressSafe("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal", &p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, &b_p_HadWH_SIG_ghw1_1_JHUGen_JECNominal);
   SetBranchAddressSafe("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal", &p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, &b_p_HadZH_SIG_ghz1_1_JHUGen_JECNominal);
   SetBranchAddressSafe("p_HadWH_mavjj_JECNominal", &p_HadWH_mavjj_JECNominal, &b_p_HadWH_mavjj_JECNominal);
   SetBranchAddressSafe("p_HadWH_mavjj_true_JECNominal", &p_HadWH_mavjj_true_JECNominal, &b_p_HadWH_mavjj_true_JECNominal);
   SetBranchAddressSafe("p_HadZH_mavjj_JECNominal", &p_HadZH_mavjj_JECNominal, &b_p_HadZH_mavjj_JECNominal);
   SetBranchAddressSafe("p_HadZH_mavjj_true_JECNominal", &p_HadZH_mavjj_true_JECNominal, &b_p_HadZH_mavjj_true_JECNominal);
   
   Notify();
}

Bool_t Tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the user if needed.
   // The return value is currently not used.

   return kTRUE;
}

void Tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Tree_cxx
