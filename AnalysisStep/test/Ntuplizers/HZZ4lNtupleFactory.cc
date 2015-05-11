#include <cassert>
#include <iostream>

#include "HZZ4lNtupleFactory.h"


bool addKinRefit = false;
bool addVtxFit = false;
bool writePhotons = false;  // Write photons in the tree. Must be set also in HZZ4lNtupleMaker.cc


HZZ4lNtupleFactory::HZZ4lNtupleFactory(TTree* outTree_input)
{
  //---- create output tree ----
  _outTree = outTree_input;
  InitializeVariables();
  InitializeBranches();
}

///--- Destructor ---
HZZ4lNtupleFactory::~HZZ4lNtupleFactory()
{
}

///---- Write an event to TTree ----
void HZZ4lNtupleFactory::FillEvent()
{
  _outTree->Fill();
  InitializeVariables(); // Reset all values and clean vectors
}

///---- Write to a text file branches declaration ----
void HZZ4lNtupleFactory::DumpBranches(TString filename) const
{
  //----- symply use MakeClass
  _outTree->MakeClass(filename);
  return;
}

void HZZ4lNtupleFactory::InitializeVariables()
{
  //Event variables
  _RunNumber = 0;
  _EventNumber = 0;
  _LumiNumber = 0;
  _IndexBestCand = 0;

  //  _Nmu = 0;
  //  _Nele = 0;

  _Nvtx = 0;
  _NObsInt =0;
  _NTrueInt =0;
  _PUWeight =0;
  _PFMET = -99;
  _nJets = 0;
  _nCleanedJets = 0;
  _nCleanedJetsPt30 = 0;
  _nCleanedJetsPt30BTagged = 0;
  _genFinalState = 0;
  _genProcessId = 0;
  _genHEPMCweight = 0;
  _trigWord =0;
  _genExtInfo =0;
  _xsec =0.;
  _dataMCweight=0;
  _HqTMCweight=0;

  //Info on generated particles
  _genHMass = 0.;
  _genHPt = 0.;

  _genZ1Mass = 0.;
  _genZ1Pt = 0.;
  
  _genZ2Mass = 0.;
  _genZ2Pt = 0.;

  _genLep1Pt = 0.;
  _genLep1Eta = 0.;
  _genLep1Phi = 0.;
  _genLep1Id = 0;

  _genLep2Pt = 0.;
  _genLep2Eta = 0.;
  _genLep2Phi = 0.;
  _genLep2Id = 0;

  _genLep3Pt = 0.;
  _genLep3Eta = 0.;
  _genLep3Phi = 0.;
  _genLep3Id = 0;

  _genLep4Pt = 0.;
  _genLep4Eta = 0.;
  _genLep4Phi = 0.;
  _genLep4Id = 0;

  _genAssocLep1Pt = 0.;
  _genAssocLep1Eta = 0.;
  _genAssocLep1Phi = 0.;
  _genAssocLep1Id = 0;

  _genAssocLep2Pt = 0.;
  _genAssocLep2Eta = 0.;
  _genAssocLep2Phi = 0.;
  _genAssocLep2Id = 0;


  //H variables
  _ZZMass=0;
  _ZZMassErr=0;
  _ZZMassErrCorr=0;
  _ZZMassPreFSR=0;
  _ZZMassRefit=0;
  _Chi2KinFit=0;
  _ZZMassCFit=0;
  _Chi2CFit=0;
  _ZZsel=0;
  _ZZPt=0;
  _ZZEta=0;
  _ZZPhi=0;
  _ZZgenIsSignal=0;
  _ZZgenIsRightPair=0;
  _CRflag=0;

  //_p0plus_melaNorm=0;
  //_p0plus_mela=0;
  //_p0minus_mela=0;
  //_p0hplus_mela=0; // 0h+, analytic distribution
  _p0plus_VAJHU=0;
  _p0minus_VAJHU=0;
  _p0plus_VAMCFM=0;
  _p0hplus_VAJHU=0; // 0h+ (high dimensional operator), vector algebra, JHUgen
  //_p1_mela=0;
  //_p1_prodIndep_mela=0;
  //_p1plus_mela=0; // 1+, analytic distribution 
  //_p1plus_prodIndep_mela=0; // 1+, analytic distribution 
  _p1_VAJHU=0;
  _p1_prodIndep_VAJHU=0;
  _p1plus_VAJHU=0; // 1+ (axial vector), vector algebra, JHUgen,
  _p1plus_prodIndep_VAJHU=0; // 1+ (axial vector), vector algebra, JHUgen,
  //_p2_mela=0;
  //_p2_prodIndep_mela=0;
  //_p2qqb_mela=0; // graviton produced by qqbar vector algebra, analytical,
  //_p2hplus_mela=0; // graviton produced by qqbar vector algebra, analytical,
  //_p2hminus_mela=0; // graviton produced by qqbar vector algebra, analytical,
  //_p2bplus_mela=0; // graviton produced by qqbar vector algebra, analytical,
  _p2_VAJHU=0;
  _p2_prodIndep_VAJHU=0;
  _p2qqb_VAJHU=0;
  _p2hplus_VAJHU=0;
  _p2hminus_VAJHU=0;
  _p2bplus_VAJHU=0;
  _p2hplus_qqb_VAJHU=0;
  _p2hplus_prodIndep_VAJHU=0;
  _p2hminus_qqb_VAJHU=0;
  _p2hminus_prodIndep_VAJHU=0;
  _p2bplus_qqb_VAJHU=0;
  _p2bplus_prodIndep_VAJHU=0;
  _p2h2plus_gg_VAJHU=0;
  _p2h2plus_qqbar_VAJHU=0;
  _p2h2plus_prodIndep_VAJHU=0;
  _p2h3plus_gg_VAJHU=0;
  _p2h3plus_qqbar_VAJHU=0;
  _p2h3plus_prodIndep_VAJHU=0;
  _p2h6plus_gg_VAJHU=0;
  _p2h6plus_qqbar_VAJHU=0;
  _p2h6plus_prodIndep_VAJHU=0;
  _p2h7plus_gg_VAJHU=0;
  _p2h7plus_qqbar_VAJHU=0;
  _p2h7plus_prodIndep_VAJHU=0;
  _p2h9minus_gg_VAJHU=0;
  _p2h9minus_qqbar_VAJHU=0;
  _p2h9minus_prodIndep_VAJHU=0;
  _p2h10minus_gg_VAJHU=0;
  _p2h10minus_qqbar_VAJHU=0;
  _p2h10minus_prodIndep_VAJHU=0;
  //_bkg_mela=bkg_mela;
  //_bkg_mela=0;
  _bkg_VAMCFM=0;
  _bkg_prodIndep_VAMCFM=0;
  _ggzz_VAMCFM=0;
  _ggzz_p0plus_VAMCFM=0;
  _ggzz_c1_VAMCFM=0;
  _ggzz_c5_VAMCFM=0;
  _ggzz_ci_VAMCFM=0;
  //_bkg_VAMCFMNorm=0;
  //_p0_pt=0;
  //_p0_y=0;
  //_bkg_pt=0;
  //_bkg_y=0;  
  _p0plus_m4l=0; 
  _bkg_m4l=0; 
  _pg1g4_mela=0;
  _pg1g4_VAJHU=0;
  _pg1g4_pi2_VAJHU=0;
  _pg1g2_pi2_VAJHU=0;
  _pg1g2_mela=0;
  _pg1g2_VAJHU=0;
  _p0plus_m4l=0;// signal m4l probability as in datacards
  _bkg_m4l=0;// backgroun m4l probability as in datacards
  _p0plus_m4l_ScaleUp=0;// signal m4l probability for systematics
  _bkg_m4l_ScaleUp=0;// backgroun m4l probability for systematics
  _p0plus_m4l_ScaleDown=0;// signal m4l probability for systematics
  _bkg_m4l_ScaleDown=0;// backgroun m4l probability for systematics
  _p0plus_m4l_ResUp=0;// signal m4l probability for systematics
  _bkg_m4l_ResUp=0;// backgroun m4l probability for systematics
  _p0plus_m4l_ResDown=0;// signal m4l probability for systematics
  _bkg_m4l_ResDown=0;     // backgroun m4l probability for systematics
  _phjj_VAJHU_old=0;
  _pvbf_VAJHU_old=0;
  _phjj_VAJHU_old_up=0;
  _pvbf_VAJHU_old_up=0;
  _phjj_VAJHU_old_dn=0;
  _pvbf_VAJHU_old_dn=0;
  _phjj_VAJHU_new=0;
  _pvbf_VAJHU_new=0;
  _phjj_VAJHU_new_up=0;
  _pvbf_VAJHU_new_up=0;
  _phjj_VAJHU_new_dn=0;
  _pvbf_VAJHU_new_dn=0;
  _p0_g1prime2_VAJHU=0;
  _pg1g1prime2_VAJHU=0;
  _Dgg10_VAMCFM=0;

  _pzzzg_VAJHU=0;
  _pzzgg_VAJHU=0;
  _pzzzg_PS_VAJHU=0;
  _pzzgg_PS_VAJHU=0;
  _p0Zgs_VAJHU=0;
  _p0gsgs_VAJHU=0;
  _p0Zgs_PS_VAJHU=0;
  _p0gsgs_PS_VAJHU=0;

  _ZZmZa=0;
  _ZZmZb=0;
  _ZZmLL4=0;
  _ZZmLL6=0;
  _ZZSIP4=0;
  //  _ZZiso34=0;

  //Z1 variables
  _Z1Mass=0;
  _Z1Pt=0;
  _Z1MassRefit=0;
  _Z1Flav=0;

  //Z2 variables
  _Z2Mass=0;
  _Z2Pt=0;
  _Z2Flav=0;

  //Angular variables
  _costhetastar=0;
  _phi=0;
  _costheta1=0;
  _costheta2=0;
  _phistar1=0;
  _phistar2=0;
  _xi=0;
  _xistar=0;

  //Lepton variables
  _Lep1Pt=0;
  _Lep1Eta=0;
  _Lep1Phi=0;
  _Lep1LepId=0;
  _Lep1SIP=0;
  _Lep1isID=0;
  _Lep1BDT=0;
  _Lep1missingHit=0;
  _Lep1ParentId=0;

  _Lep2Pt=0;
  _Lep2Eta=0;
  _Lep2Phi=0;
  _Lep2LepId=0;
  _Lep2SIP=0;
  _Lep2isID=0;
  _Lep2BDT=0;
  _Lep2missingHit=0;
  _Lep2ParentId=0;

  _Lep3Pt=0;
  _Lep3Eta=0;
  _Lep3Phi=0;
  _Lep3LepId=0;
  _Lep3SIP=0;
  _Lep3isID=0;
  _Lep3BDT=0;
  _Lep3missingHit=0;
  _Lep3ParentId=0;

  _Lep4Pt=0;
  _Lep4Eta=0;
  _Lep4Phi=0;
  _Lep4LepId=0;
  _Lep4SIP=0;
  _Lep4isID=0;
  _Lep4BDT=0;
  _Lep4missingHit=0;
  _Lep4ParentId=0;

  //Lepton isolation variables
  _Lep1chargedHadIso=0;
  _Lep1neutralHadIso=0;
  _Lep1photonIso=0;
  _Lep1combRelIsoPF=0;

  _Lep2chargedHadIso=0;
  _Lep2neutralHadIso=0;
  _Lep2photonIso=0;
  _Lep2combRelIsoPF=0;

  _Lep3chargedHadIso=0;
  _Lep3neutralHadIso=0;
  _Lep3photonIso=0;
  _Lep3combRelIsoPF=0;

  _Lep4chargedHadIso=0;
  _Lep4neutralHadIso=0;
  _Lep4photonIso=0;
  _Lep4combRelIsoPF=0;

  //Photon variables
  _PhotPt=0;
  _PhotEta=0;
  _PhotPhi=0; 

  //Jet variables
  _JetPt.clear();
  _JetEta.clear();
  _JetPhi.clear(); 
  _JetMass.clear(); 
  _JetBTagger.clear();
  _JetIsBtagged.clear();
  _JetQGLikelihood.clear();
  _JetSigma.clear();

  _DiJetMass=-99;
  _DiJetMassPlus=-99;
  _DiJetMassMinus=-99;
  _DiJetDEta=-99;
  _DiJetFisher=-99;
  
  //Categorization-related variables
  _nExtraLep=0;
  _nExtraZ=0;

  //Variables of extra leptons
  _ExtraLep1Pt=0;
  _ExtraLep1Eta=0;
  _ExtraLep1Phi=0;
  _ExtraLep1LepId=0;
  _ExtraLep1SIP=0;
  _ExtraLep1isID=0;
  _ExtraLep1BDT=0;
  _ExtraLep1missingHit=0;
  _ExtraLep1chargedHadIso=0;
  _ExtraLep1neutralHadIso=0;
  _ExtraLep1photonIso=0;
  _ExtraLep1combRelIsoPF=0;

  _ExtraLep2Pt=0;
  _ExtraLep2Eta=0;
  _ExtraLep2Phi=0;
  _ExtraLep2LepId=0;
  _ExtraLep2SIP=0;
  _ExtraLep2isID=0;
  _ExtraLep2BDT=0;
  _ExtraLep2missingHit=0;
  _ExtraLep2chargedHadIso=0;
  _ExtraLep2neutralHadIso=0;
  _ExtraLep2photonIso=0;
  _ExtraLep2combRelIsoPF=0;

  _ExtraLep3Pt=0;
  _ExtraLep3Eta=0;
  _ExtraLep3Phi=0;
  _ExtraLep3LepId=0;
  _ExtraLep3SIP=0;
  _ExtraLep3isID=0;
  _ExtraLep3BDT=0;
  _ExtraLep3missingHit=0;
  _ExtraLep3chargedHadIso=0;
  _ExtraLep3neutralHadIso=0;
  _ExtraLep3photonIso=0;
  _ExtraLep3combRelIsoPF=0;


  return;
}

void HZZ4lNtupleFactory::InitializeBranches()
{
  //Event variables
  _outTree->Branch("RunNumber",&_RunNumber,"RunNumber/I");
  _outTree->Branch("EventNumber",&_EventNumber,"EventNumber/L");
  _outTree->Branch("LumiNumber",&_LumiNumber,"LumiNumber/I");
  _outTree->Branch("iBC",&_IndexBestCand,"iBC/I");
  _outTree->Branch("Nvtx",&_Nvtx,"Nvtx/I");
  _outTree->Branch("NObsInt",&_NObsInt,"NObsInt/I");
  _outTree->Branch("NTrueInt",&_NTrueInt,"NTrueInt/F");
  _outTree->Branch("PUWeight12",&_PUWeight,"PUWeight12/F");
  _outTree->Branch("PFMET",&_PFMET,"PFMET/F");
  _outTree->Branch("nJets",&_nJets,"nJets/I");
  _outTree->Branch("nCleanedJets",&_nCleanedJets,"nCleanedJets/I");
  _outTree->Branch("nCleanedJetsPt30",&_nCleanedJetsPt30,"nCleanedJetsPt30/I");
  _outTree->Branch("nCleanedJetsPt30BTagged",&_nCleanedJetsPt30BTagged,"nCleanedJetsPt30BTagged/I");
  _outTree->Branch("trigWord",&_trigWord,"trigWord/S");
  
  //H variables
  _outTree->Branch("ZZMass",&_ZZMass,"ZZMass/F");
  _outTree->Branch("ZZMassErr",&_ZZMassErr);
  _outTree->Branch("ZZMassErrCorr",&_ZZMassErrCorr);
  _outTree->Branch("ZZMassPreFSR",&_ZZMassPreFSR);
  _outTree->Branch("ZZsel",&_ZZsel,"ZZsel/I");
  _outTree->Branch("ZZPt",&_ZZPt,"ZZPt/F");
  _outTree->Branch("ZZEta",&_ZZEta);
  _outTree->Branch("ZZPhi",&_ZZPhi);
  _outTree->Branch("CRflag",&_CRflag);


  //Z1 variables
  _outTree->Branch("Z1Mass",&_Z1Mass);
  _outTree->Branch("Z1Pt",&_Z1Pt);
  _outTree->Branch("Z1Flav",&_Z1Flav);  

  //Kin refitted info
  if (addKinRefit) {
    _outTree->Branch("ZZMassRefit",&_ZZMassRefit);
    _outTree->Branch("ZZChi2KinFit",&_Chi2KinFit);
    _outTree->Branch("Z1MassRefit",&_Z1MassRefit);
  }

  if (addVtxFit){
    _outTree->Branch("ZZMassCFit",&_ZZMassCFit);
    _outTree->Branch("ZZChi2CFit",&_Chi2CFit);

  }

  //Z2 variables
  _outTree->Branch("Z2Mass",&_Z2Mass);
  _outTree->Branch("Z2Pt",&_Z2Pt);
  _outTree->Branch("Z2Flav",&_Z2Flav,"ZZFlav/S");  

  //Angular variables
  _outTree->Branch("costhetastar",&_costhetastar,"costhetastar/F");
  _outTree->Branch("helphi",&_phi);
//   _outTree->Branch("helphiZ1",&_helphiZ1);
//   _outTree->Branch("helphiZ2",&_helphiZ2);
  _outTree->Branch("helcosthetaZ1",&_costheta1);
  _outTree->Branch("helcosthetaZ2",&_costheta2);
  _outTree->Branch("phistarZ1",&_phistar1);
  _outTree->Branch("phistarZ2",&_phistar2);
  _outTree->Branch("xi",&_xi);
  _outTree->Branch("xistar",&_xistar);

  //Lepton variables
  _outTree->Branch("Lep1Pt",&_Lep1Pt);
  _outTree->Branch("Lep1Eta",&_Lep1Eta);
  _outTree->Branch("Lep1Phi",&_Lep1Phi);
  _outTree->Branch("Lep1LepId",&_Lep1LepId);
  _outTree->Branch("Lep1SIP",&_Lep1SIP);
  _outTree->Branch("Lep1isID",&_Lep1isID);
  _outTree->Branch("Lep1BDT",&_Lep1BDT);
  _outTree->Branch("Lep1missingHit",&_Lep1missingHit);
  //--> Commented as these parentIDs are not computed at the moment
  //  _outTree->Branch("Lep1ParentId",&_Lep1ParentId);

  _outTree->Branch("Lep2Pt",&_Lep2Pt);
  _outTree->Branch("Lep2Eta",&_Lep2Eta);
  _outTree->Branch("Lep2Phi",&_Lep2Phi);
  _outTree->Branch("Lep2LepId",&_Lep2LepId);
  _outTree->Branch("Lep2SIP",&_Lep2SIP);
  _outTree->Branch("Lep2isID",&_Lep2isID);
  _outTree->Branch("Lep2BDT",&_Lep2BDT);
  _outTree->Branch("Lep2missingHit",&_Lep2missingHit);
  //  _outTree->Branch("Lep2ParentId",&_Lep2ParentId);

  _outTree->Branch("Lep3Pt",&_Lep3Pt);
  _outTree->Branch("Lep3Eta",&_Lep3Eta);
  _outTree->Branch("Lep3Phi",&_Lep3Phi);
  _outTree->Branch("Lep3LepId",&_Lep3LepId);
  _outTree->Branch("Lep3SIP",&_Lep3SIP);
  _outTree->Branch("Lep3isID",&_Lep3isID);
  _outTree->Branch("Lep3BDT",&_Lep3BDT);
  _outTree->Branch("Lep3missingHit",&_Lep3missingHit);
  //  _outTree->Branch("Lep3ParentId",&_Lep3ParentId);

  _outTree->Branch("Lep4Pt",&_Lep4Pt);
  _outTree->Branch("Lep4Eta",&_Lep4Eta);
  _outTree->Branch("Lep4Phi",&_Lep4Phi);
  _outTree->Branch("Lep4LepId",&_Lep4LepId);
  _outTree->Branch("Lep4SIP",&_Lep4SIP);
  _outTree->Branch("Lep4isID",&_Lep4isID);
  _outTree->Branch("Lep4BDT",&_Lep4BDT);
  _outTree->Branch("Lep4missingHit",&_Lep4missingHit);
  //  _outTree->Branch("Lep4ParentId",&_Lep4ParentId);

  //Lepton isolation variables
  _outTree->Branch("Lep1chargedHadIso",&_Lep1chargedHadIso);
  _outTree->Branch("Lep1neutralHadIso",&_Lep1neutralHadIso);
  _outTree->Branch("Lep1photonIso",&_Lep1photonIso);
  _outTree->Branch("Lep1combRelIsoPF",&_Lep1combRelIsoPF);

  _outTree->Branch("Lep2chargedHadIso",&_Lep2chargedHadIso);
  _outTree->Branch("Lep2neutralHadIso",&_Lep2neutralHadIso);
  _outTree->Branch("Lep2photonIso",&_Lep2photonIso);
  _outTree->Branch("Lep2combRelIsoPF",&_Lep2combRelIsoPF);

  _outTree->Branch("Lep3chargedHadIso",&_Lep3chargedHadIso);
  _outTree->Branch("Lep3neutralHadIso",&_Lep3neutralHadIso);
  _outTree->Branch("Lep3photonIso",&_Lep3photonIso);
  _outTree->Branch("Lep3combRelIsoPF",&_Lep3combRelIsoPF);

  _outTree->Branch("Lep4chargedHadIso",&_Lep4chargedHadIso);
  _outTree->Branch("Lep4neutralHadIso",&_Lep4neutralHadIso);
  _outTree->Branch("Lep4photonIso",&_Lep4photonIso);
  _outTree->Branch("Lep4combRelIsoPF",&_Lep4combRelIsoPF);

  //Photon variables
  if (writePhotons) {
    _outTree->Branch("PhotPt",&_PhotPt);
    _outTree->Branch("PhotEta",&_PhotEta);
    _outTree->Branch("PhotPhi",&_PhotPhi);
  }

  //Discriminants
  _outTree->Branch("p0plus_VAJHU",&_p0plus_VAJHU);
  _outTree->Branch("p0minus_VAJHU",&_p0minus_VAJHU);
  _outTree->Branch("p0plus_VAMCFM",&_p0plus_VAMCFM);
  _outTree->Branch("p0hplus_VAJHU",&_p0hplus_VAJHU); // 0h+ (high dimensional operator), vector algebra, JHUgen
  _outTree->Branch("p1_VAJHU",&_p1_VAJHU);
  _outTree->Branch("p1_prodIndep_VAJHU",&_p1_prodIndep_VAJHU);
  _outTree->Branch("p1plus_VAJHU",&_p1plus_VAJHU); // 1+ (axial vector), vector algebra, JHUgen,
  _outTree->Branch("p1plus_prodIndep_VAJHU",&_p1plus_prodIndep_VAJHU); // 1+ (axial vector), vector algebra, JHUgen,
  _outTree->Branch("p2_VAJHU",&_p2_VAJHU);
  _outTree->Branch("p2_prodIndep_VAJHU",&_p2_prodIndep_VAJHU);
  _outTree->Branch("p2qqb_VAJHU",&_p2qqb_VAJHU);
  _outTree->Branch("p2hplus_VAJHU",&_p2hplus_VAJHU);
  _outTree->Branch("p2hminus_VAJHU",&_p2hminus_VAJHU);
  _outTree->Branch("p2bplus_VAJHU",&_p2bplus_VAJHU);

  _outTree->Branch("p2hplus_qqb_VAJHU",&_p2hplus_qqb_VAJHU);                                   
  _outTree->Branch("p2hplus_prodIndep_VAJHU",&_p2hplus_prodIndep_VAJHU);             
  _outTree->Branch("p2hminus_qqb_VAJHU",&_p2hminus_qqb_VAJHU);                          
  _outTree->Branch("p2hminus_prodIndep_VAJHU",&_p2hminus_prodIndep_VAJHU);    
  _outTree->Branch("p2bplus_qqb_VAJHU",&_p2bplus_qqb_VAJHU);                                   
  _outTree->Branch("p2bplus_prodIndep_VAJHU",&_p2bplus_prodIndep_VAJHU);             
  _outTree->Branch("p2h2plus_gg_VAJHU",&_p2h2plus_gg_VAJHU);                                             
  _outTree->Branch("p2h2plus_qqbar_VAJHU"               ,               &_p2h2plus_qqbar_VAJHU);                
  _outTree->Branch("p2h2plus_prodIndep_VAJHU"   ,       &_p2h2plus_prodIndep_VAJHU);    
  _outTree->Branch("p2h3plus_gg_VAJHU"          ,       &_p2h3plus_gg_VAJHU);           
  _outTree->Branch("p2h3plus_qqbar_VAJHU"       ,       &_p2h3plus_qqbar_VAJHU);        
  _outTree->Branch("p2h3plus_prodIndep_VAJHU"   ,       &_p2h3plus_prodIndep_VAJHU);    
  _outTree->Branch("p2h6plus_gg_VAJHU"          ,       &_p2h6plus_gg_VAJHU);           
  _outTree->Branch("p2h6plus_qqbar_VAJHU"       ,       &_p2h6plus_qqbar_VAJHU);        
  _outTree->Branch("p2h6plus_prodIndep_VAJHU"   ,       &_p2h6plus_prodIndep_VAJHU);    
  _outTree->Branch("p2h7plus_gg_VAJHU"          ,       &_p2h7plus_gg_VAJHU);   
  _outTree->Branch("p2h7plus_qqbar_VAJHU"       ,       &_p2h7plus_qqbar_VAJHU);        
  _outTree->Branch("p2h7plus_prodIndep_VAJHU"   ,       &_p2h7plus_prodIndep_VAJHU);    
  _outTree->Branch("p2h9minus_gg_VAJHU"         ,               &_p2h9minus_gg_VAJHU);          
  _outTree->Branch("p2h9minus_qqbar_VAJHU"      ,               &_p2h9minus_qqbar_VAJHU);       
  _outTree->Branch("p2h9minus_prodIndep_VAJHU"  ,               &_p2h9minus_prodIndep_VAJHU);   
  _outTree->Branch("p2h10minus_gg_VAJHU"       ,                &_p2h10minus_gg_VAJHU);       
  _outTree->Branch("p2h10minus_qqbar_VAJHU"    ,                &_p2h10minus_qqbar_VAJHU);  
  _outTree->Branch("p2h10minus_prodIndep_VAJHU",                &_p2h10minus_prodIndep_VAJHU);
  _outTree->Branch("bkg_VAMCFM",&_bkg_VAMCFM);
  _outTree->Branch("bkg_prodIndep_VAMCFM",&_bkg_prodIndep_VAMCFM);
  _outTree->Branch("ggzz_VAMCFM",&_ggzz_VAMCFM);
  _outTree->Branch("ggzz_p0plus_VAMCFM",&_ggzz_p0plus_VAMCFM);
  _outTree->Branch("ggzz_c1_VAMCFM",&_ggzz_c1_VAMCFM);
  _outTree->Branch("ggzz_c5_VAMCFM",&_ggzz_c5_VAMCFM);
  _outTree->Branch("ggzz_ci_VAMCFM",&_ggzz_ci_VAMCFM);
  
  _outTree->Branch("p0plus_m4l",&_p0plus_m4l);// signal m4l probability as in datacards
  _outTree->Branch("bkg_m4l",&_bkg_m4l);// backgroun m4l probability as in datacards
  _outTree->Branch("pg1g4_mela",&_pg1g4_mela);
  _outTree->Branch("pg1g4_VAJHU",&_pg1g4_VAJHU);
  _outTree->Branch("pg1g4_pi2_VAJHU",&_pg1g4_pi2_VAJHU);
  _outTree->Branch("pg1g2_pi2_VAJHU",&_pg1g2_pi2_VAJHU);
  _outTree->Branch("pg1g2_mela",&_pg1g2_mela);
  _outTree->Branch("pg1g2_VAJHU",&_pg1g2_VAJHU);

  _outTree->Branch("p0plus_m4l_ScaleUp",&_p0plus_m4l_ScaleUp);// signal m4l probability for systematics
  _outTree->Branch("bkg_m4l_ScaleUp",&_bkg_m4l_ScaleUp);// backgroun m4l probability for systematics
  _outTree->Branch("p0plus_m4l_ScaleDown",&_p0plus_m4l_ScaleDown);// signal m4l probability for systematics
  _outTree->Branch("bkg_m4l_ScaleDown",&_bkg_m4l_ScaleDown);// backgroun m4l probability for systematics
  _outTree->Branch("p0plus_m4l_ResUp",&_p0plus_m4l_ResUp);// signal m4l probability for systematics
  _outTree->Branch("bkg_m4l_ResUp",&_bkg_m4l_ResUp);// backgroun m4l probability for systematics
  _outTree->Branch("p0plus_m4l_ResDown",&_p0plus_m4l_ResDown);// signal m4l probability for systematics
  _outTree->Branch("bkg_m4l_ResDown",&_bkg_m4l_ResDown);// backgroun m4l probability for systematics

  //Production MELA
  _outTree->Branch("phjj_VAJHU_old",&_phjj_VAJHU_old);
  _outTree->Branch("pvbf_VAJHU_old",&_pvbf_VAJHU_old);
  _outTree->Branch("phjj_VAJHU_old_up",&_phjj_VAJHU_old_up);
  _outTree->Branch("pvbf_VAJHU_old_up",&_pvbf_VAJHU_old_up);
  _outTree->Branch("phjj_VAJHU_old_dn",&_phjj_VAJHU_old_dn);
  _outTree->Branch("pvbf_VAJHU_old_dn",&_pvbf_VAJHU_old_dn);
  _outTree->Branch("phjj_VAJHU_new",&_phjj_VAJHU_new);
  _outTree->Branch("pvbf_VAJHU_new",&_pvbf_VAJHU_new);
  _outTree->Branch("phjj_VAJHU_new_up",&_phjj_VAJHU_new_up);
  _outTree->Branch("pvbf_VAJHU_new_up",&_pvbf_VAJHU_new_up);
  _outTree->Branch("phjj_VAJHU_new_dn",&_phjj_VAJHU_new_dn);
  _outTree->Branch("pvbf_VAJHU_new_dn",&_pvbf_VAJHU_new_dn);
  _outTree->Branch("p0_g1prime2_VAJHU",&_p0_g1prime2_VAJHU);
  _outTree->Branch("pg1g1prime2_VAJHU",&_pg1g1prime2_VAJHU);
  _outTree->Branch("Dgg10_VAMCFM",&_Dgg10_VAMCFM);

  _outTree->Branch("pzzzg_VAJHU",      &_pzzzg_VAJHU);
  _outTree->Branch("pzzgg_VAJHU",      &_pzzgg_VAJHU);
  _outTree->Branch("pzzzg_PS_VAJHU",   &_pzzzg_PS_VAJHU);
  _outTree->Branch("pzzgg_PS_VAJHU",   &_pzzgg_PS_VAJHU);
  _outTree->Branch("p0Zgs_VAJHU",      &_p0Zgs_VAJHU);
  _outTree->Branch("p0gsgs_VAJHU",     &_p0gsgs_VAJHU);
  _outTree->Branch("p0Zgs_PS_VAJHU",   &_p0Zgs_PS_VAJHU);
  _outTree->Branch("p0gsgs_PS_VAJHU",  &_p0gsgs_PS_VAJHU);
  
  //Jet variables
  _outTree->Branch("JetPt",&_JetPt);
  _outTree->Branch("JetEta",&_JetEta);
  _outTree->Branch("JetPhi",&_JetPhi);
  _outTree->Branch("JetMass",&_JetMass);
  _outTree->Branch("JetBTagger",&_JetBTagger);
  _outTree->Branch("JetIsBtagged",&_JetIsBtagged);
  _outTree->Branch("JetQGLikelihood",&_JetQGLikelihood);
  _outTree->Branch("JetSigma",&_JetSigma);
  _outTree->Branch("DiJetMass",&_DiJetMass,"DiJetMass/F");
  _outTree->Branch("DiJetMassPlus",&_DiJetMassPlus,"DiJetMassPlus/F");
  _outTree->Branch("DiJetMassMinus",&_DiJetMassMinus,"DiJetMassMinus/F");
  _outTree->Branch("DiJetDEta",&_DiJetDEta,"DiJetDEta/F");
  _outTree->Branch("DiJetFisher",&_DiJetFisher,"DiJetFisher/F");

  //Categorization-related variables
  _outTree->Branch("nExtraLep",&_nExtraLep);
  _outTree->Branch("nExtraZ",&_nExtraZ);

  //Variables of extra leptons
  _outTree->Branch("ExtraLep1Pt",&_ExtraLep1Pt);
  _outTree->Branch("ExtraLep1Eta",&_ExtraLep1Eta);
  _outTree->Branch("ExtraLep1Phi",&_ExtraLep1Phi);
  _outTree->Branch("ExtraLep1LepId",&_ExtraLep1LepId);

  _outTree->Branch("ExtraLep2Pt",&_ExtraLep2Pt);
  _outTree->Branch("ExtraLep2Eta",&_ExtraLep2Eta);
  _outTree->Branch("ExtraLep2Phi",&_ExtraLep2Phi);
  _outTree->Branch("ExtraLep2LepId",&_ExtraLep2LepId);

  _outTree->Branch("ExtraLep3Pt",&_ExtraLep3Pt);
  _outTree->Branch("ExtraLep3Eta",&_ExtraLep3Eta);
  _outTree->Branch("ExtraLep3Phi",&_ExtraLep3Phi);
  _outTree->Branch("ExtraLep3LepId",&_ExtraLep3LepId);

// Extended information on extra leptons.
// Currently not used, skip them for the time being
//   _outTree->Branch("ExtraLep1SIP",&_ExtraLep1SIP);
//   _outTree->Branch("ExtraLep1isID",&_ExtraLep1isID);
//   _outTree->Branch("ExtraLep1BDT",&_ExtraLep1BDT);
//   _outTree->Branch("ExtraLep1missingHit",&_ExtraLep1missingHit);
//   _outTree->Branch("ExtraLep1chargedHadIso",&_ExtraLep1chargedHadIso);
//   _outTree->Branch("ExtraLep1neutralHadIso",&_ExtraLep1neutralHadIso);
//   _outTree->Branch("ExtraLep1photonIso",&_ExtraLep1photonIso);
//   _outTree->Branch("ExtraLep1combRelIsoPF",&_ExtraLep1combRelIsoPF);
//
//   _outTree->Branch("ExtraLep2SIP",&_ExtraLep2SIP);
//   _outTree->Branch("ExtraLep2isID",&_ExtraLep2isID);
//   _outTree->Branch("ExtraLep2BDT",&_ExtraLep2BDT);
//   _outTree->Branch("ExtraLep2missingHit",&_ExtraLep2missingHit);
//   _outTree->Branch("ExtraLep2chargedHadIso",&_ExtraLep2chargedHadIso);
//   _outTree->Branch("ExtraLep2neutralHadIso",&_ExtraLep2neutralHadIso);
//   _outTree->Branch("ExtraLep2photonIso",&_ExtraLep2photonIso);
//   _outTree->Branch("ExtraLep2combRelIsoPF",&_ExtraLep2combRelIsoPF);
//
//   _outTree->Branch("ExtraLep3SIP",&_ExtraLep3SIP);
//   _outTree->Branch("ExtraLep3isID",&_ExtraLep3isID);
//   _outTree->Branch("ExtraLep3BDT",&_ExtraLep3BDT);
//   _outTree->Branch("ExtraLep3missingHit",&_ExtraLep3missingHit);
//   _outTree->Branch("ExtraLep3chargedHadIso",&_ExtraLep3chargedHadIso);
//   _outTree->Branch("ExtraLep3neutralHadIso",&_ExtraLep3neutralHadIso);
//   _outTree->Branch("ExtraLep3photonIso",&_ExtraLep3photonIso);
//   _outTree->Branch("ExtraLep3combRelIsoPF",&_ExtraLep3combRelIsoPF);

  // Gen variables
  // FIXME: don't book these fo data so we save some disk space...
  _outTree->Branch("genFinalState",&_genFinalState,"genFinalState/I");
  _outTree->Branch("genProcessId",&_genProcessId,"genProcessId/I");
  _outTree->Branch("genHEPMCweight",&_genHEPMCweight,"genHEPMCweight/F");
  _outTree->Branch("genExtInfo",&_genExtInfo,"genExtInfo/S");
  _outTree->Branch("xsec",&_xsec,"xsec/F");
  _outTree->Branch("dataMCWeight",&_dataMCweight,"dataMCweight/F");
  _outTree->Branch("HqTMCweight",&_HqTMCweight,"HqTMCweight/F");
  
  _outTree->Branch("GenHMass",&_genHMass,"GenHMass/F");
  _outTree->Branch("GenHPt",&_genHPt,"GenHPt/F");

  _outTree->Branch("GenZ1Mass",&_genZ1Mass,"GenZ1Mass/F");
  _outTree->Branch("GenZ1Pt",&_genZ1Pt,"GenZ1Pt/F");

  _outTree->Branch("GenZ2Mass",&_genZ2Mass,"GenZ2Mass/F");
  _outTree->Branch("GenZ2Pt",&_genZ2Pt,"GenZ2Pt/F");

  _outTree->Branch("GenLep1Pt",&_genLep1Pt,"GenLep1Pt/F");
  _outTree->Branch("GenLep1Eta",&_genLep1Eta,"GenLep1Eta/F");
  _outTree->Branch("GenLep1Phi",&_genLep1Phi,"GenLep1Phi/F");
  _outTree->Branch("GenLep1Id",&_genLep1Id,"GenLep1Id/S");

  _outTree->Branch("GenLep2Pt",&_genLep2Pt,"GenLep2Pt/F");
  _outTree->Branch("GenLep2Eta",&_genLep2Eta,"GenLep2Eta/F");
  _outTree->Branch("GenLep2Phi",&_genLep2Phi,"GenLep2Phi/F");
  _outTree->Branch("GenLep2Id",&_genLep2Id,"GenLep2Id/S");

  _outTree->Branch("GenLep3Pt",&_genLep3Pt,"GenLep3Pt/F");
  _outTree->Branch("GenLep3Eta",&_genLep3Eta,"GenLep3Eta/F");
  _outTree->Branch("GenLep3Phi",&_genLep3Phi,"GenLep3Phi/F");
  _outTree->Branch("GenLep3Id",&_genLep3Id,"GenLep3Id/S");

  _outTree->Branch("GenLep4Pt",&_genLep4Pt,"GenLep4Pt/F");
  _outTree->Branch("GenLep4Eta",&_genLep4Eta,"GenLep4Eta/F");
  _outTree->Branch("GenLep4Phi",&_genLep4Phi,"GenLep4Phi/F");
  _outTree->Branch("GenLep4Id",&_genLep4Id,"GenLep4Id/S");

  _outTree->Branch("GenAssocLep1Pt",&_genAssocLep1Pt,"GenAssocLep1Pt/F");
  _outTree->Branch("GenAssocLep1Eta",&_genAssocLep1Eta,"GenAssocLep1Eta/F");
  _outTree->Branch("GenAssocLep1Phi",&_genAssocLep1Phi,"GenAssocLep1Phi/F");
  _outTree->Branch("GenAssocLep1Id",&_genAssocLep1Id,"GenAssocLep1Id/S");

  _outTree->Branch("GenAssocLep2Pt",&_genAssocLep2Pt,"GenAssocLep2Pt/F");
  _outTree->Branch("GenAssocLep2Eta",&_genAssocLep2Eta,"GenAssocLep2Eta/F");
  _outTree->Branch("GenAssocLep2Phi",&_genAssocLep2Phi,"GenAssocLep2Phi/F");
  _outTree->Branch("GenAssocLep2Id",&_genAssocLep2Id,"GenAssocLep2Id/S");

  return;
}

void HZZ4lNtupleFactory::createNewCandidate()
{
  _firstZStored = false;
  _LeptonIndex = 1;
  _LeptonIsoIndex = 1;

  return;
}

void HZZ4lNtupleFactory::FillHGenInfo(math::XYZTLorentzVector pH,float w)
{
  _genHMass = pH.M();
  _genHPt = pH.Pt();
  
  _HqTMCweight=w;
  return;
}

void HZZ4lNtupleFactory::FillZGenInfo(math::XYZTLorentzVector pZ1, math::XYZTLorentzVector pZ2)
{
  _genZ1Mass = pZ1.M();
  _genZ1Pt = pZ1.Pt();

  _genZ2Mass = pZ2.M();
  _genZ2Pt = pZ2.Pt();

  return;
}

void HZZ4lNtupleFactory::FillLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, Short_t Lep3Id, Short_t Lep4Id, 
                                        math::XYZTLorentzVector Lep1, math::XYZTLorentzVector Lep2, math::XYZTLorentzVector Lep3, math::XYZTLorentzVector Lep4, float weight)
{
  _genLep1Pt = Lep1.Pt();
  _genLep1Eta = Lep1.Eta();
  _genLep1Phi = Lep1.Phi();
  _genLep1Id  = Lep1Id;

  _genLep2Pt = Lep2.Pt();
  _genLep2Eta = Lep2.Eta();
  _genLep2Phi = Lep2.Phi();
  _genLep2Id  = Lep2Id;

  _genLep3Pt = Lep3.Pt();
  _genLep3Eta = Lep3.Eta();
  _genLep3Phi = Lep3.Phi();
  _genLep3Id  = Lep3Id;

  _genLep4Pt = Lep4.Pt();
  _genLep4Eta = Lep4.Eta();
  _genLep4Phi = Lep4.Phi();
  _genLep4Id  = Lep4Id;

  _dataMCweight=weight;
  
  return;
}

void HZZ4lNtupleFactory::FillAssocLepGenInfo(std::vector<const reco::Candidate *>& AssocLeps)
{
  if (AssocLeps.size() >= 1) {
    _genAssocLep1Pt  = AssocLeps.at(0)->p4().Pt();
    _genAssocLep1Eta = AssocLeps.at(0)->p4().Eta();
    _genAssocLep1Phi = AssocLeps.at(0)->p4().Phi();
    _genAssocLep1Id  = AssocLeps.at(0)->pdgId();
  }
  if (AssocLeps.size() >= 2) {
    _genAssocLep2Pt  = AssocLeps.at(1)->p4().Pt();
    _genAssocLep2Eta = AssocLeps.at(1)->p4().Eta();
    _genAssocLep2Phi = AssocLeps.at(1)->p4().Phi();
    _genAssocLep2Id  = AssocLeps.at(1)->pdgId();
  }

  return;
}

void HZZ4lNtupleFactory::FillEventInfo(Int_t RunNumber, const Long64_t EventNumber, Int_t LumiNumber, Int_t IndexBestCand, Int_t Nvtx, 
                                       Int_t NObsInt, Float_t NTrueInt, Float_t PUweight, Float_t PFMET, Int_t nJets, Int_t nCleanedJets, 
                                       Int_t nCleanedJetsPt30, Int_t nCleanedJetsPt30BTagged, Int_t genFinalState, Int_t genProcessId, 
                                       Float_t genHEPMCweight, Short_t trigWord,  Short_t genExtInfo, Float_t xsec)
{
  _RunNumber = RunNumber;
  _EventNumber = EventNumber;
  _LumiNumber = LumiNumber;
  _IndexBestCand = IndexBestCand;
  _Nvtx = Nvtx;
  _NObsInt =NObsInt;
  _NTrueInt =NTrueInt;
  _PUWeight =PUweight;
  _xsec=xsec;
  
  _PFMET = PFMET;
  _nJets = nJets;
  _nCleanedJets = nCleanedJets;
  _nCleanedJetsPt30 = nCleanedJetsPt30;
  _nCleanedJetsPt30BTagged = nCleanedJetsPt30BTagged;
  _genFinalState = genFinalState;
  _genProcessId = genProcessId;
  _genHEPMCweight = genHEPMCweight;
  _trigWord = trigWord;
  _genExtInfo = genExtInfo;
  return;
}

void HZZ4lNtupleFactory::FillHInfo(Float_t ZZMass, Float_t ZZMassErr, Float_t ZZMassErrCorr, Float_t ZZMassPreFSR, Float_t ZZMassRefit, Float_t Chi2KinFit, Float_t ZZMassCFit, Float_t Chi2CFit, Int_t ZZsel, Float_t ZZPt, Float_t ZZEta, Float_t ZZPhi, Int_t isSignal, Int_t isRightPair, Int_t CRflag)
{
  _ZZMass=ZZMass;
  _ZZMassErr=ZZMassErr;
  _ZZMassErrCorr=ZZMassErrCorr;
  _ZZMassPreFSR=ZZMassPreFSR;
  _ZZMassRefit=ZZMassRefit;
  _Chi2KinFit=Chi2KinFit;
  _ZZMassCFit=ZZMassCFit;
  _Chi2CFit=Chi2CFit;
  _ZZsel=ZZsel;
  _ZZPt=ZZPt;
  _ZZEta=ZZEta;
  _ZZPhi=ZZPhi;
  _ZZgenIsSignal=isSignal;
  _ZZgenIsRightPair=isRightPair;
  _CRflag=CRflag;

  return;
}

void HZZ4lNtupleFactory::FillProbability(Float_t p0plus_VAJHU,
                                         Float_t p0minus_VAJHU,
                                         Float_t p0plus_VAMCFM,
                                         Float_t p0hplus_VAJHU, // 0h+ (high dimensional operator), vector algebra, JHUgen
                                         Float_t p1_VAJHU,
                                         Float_t p1_prodIndep_VAJHU,
                                         Float_t p1plus_VAJHU, // 1+ (axial vector), vector algebra, JHUgen,
                                         Float_t p1plus_prodIndep_VAJHU, // 1+ (axial vector), vector algebra, JHUgen,
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
                                         ){
  _p0plus_VAJHU=p0plus_VAJHU;
  _p0minus_VAJHU=p0minus_VAJHU;
  _p0plus_VAMCFM=p0plus_VAMCFM;
  _p0hplus_VAJHU=p0hplus_VAJHU;// 0h+ (high dimensional operator, vector algebra, JHUgen
  _p1_VAJHU=p1_VAJHU;
  _p1_prodIndep_VAJHU=p1_prodIndep_VAJHU;
  _p1plus_VAJHU=p1plus_VAJHU;// 1+ (axial vector, vector algebra, JHUgen,
  _p1plus_prodIndep_VAJHU=p1plus_prodIndep_VAJHU;// 1+ (axial vector, vector algebra, JHUgen,
  _p2_VAJHU=p2_VAJHU;
  _p2_prodIndep_VAJHU=p2_prodIndep_VAJHU;
  _p2qqb_VAJHU=p2qqb_VAJHU;
  _p2hplus_VAJHU=p2hplus_VAJHU;
  _p2hminus_VAJHU=p2hminus_VAJHU;
  _p2bplus_VAJHU=p2bplus_VAJHU;
  _p2hplus_qqb_VAJHU=                                   p2hplus_qqb_VAJHU;                                     
  _p2hplus_prodIndep_VAJHU=             p2hplus_prodIndep_VAJHU;               
  _p2hminus_qqb_VAJHU=                          p2hminus_qqb_VAJHU;                            
  _p2hminus_prodIndep_VAJHU=    p2hminus_prodIndep_VAJHU;      
  _p2bplus_qqb_VAJHU=                                   p2bplus_qqb_VAJHU;                                     
  _p2bplus_prodIndep_VAJHU=             p2bplus_prodIndep_VAJHU;               
  _p2h2plus_gg_VAJHU=                   p2h2plus_gg_VAJHU;                                               
  _p2h2plus_qqbar_VAJHU=                p2h2plus_qqbar_VAJHU;                  
  _p2h2plus_prodIndep_VAJHU=    p2h2plus_prodIndep_VAJHU;      
  _p2h3plus_gg_VAJHU=           p2h3plus_gg_VAJHU;             
  _p2h3plus_qqbar_VAJHU=        p2h3plus_qqbar_VAJHU;          
  _p2h3plus_prodIndep_VAJHU=    p2h3plus_prodIndep_VAJHU;      
  _p2h6plus_gg_VAJHU=           p2h6plus_gg_VAJHU;             
  _p2h6plus_qqbar_VAJHU=        p2h6plus_qqbar_VAJHU;          
  _p2h6plus_prodIndep_VAJHU=    p2h6plus_prodIndep_VAJHU;      
  _p2h7plus_gg_VAJHU=           p2h7plus_gg_VAJHU;     
  _p2h7plus_qqbar_VAJHU=        p2h7plus_qqbar_VAJHU;          
  _p2h7plus_prodIndep_VAJHU=    p2h7plus_prodIndep_VAJHU;      
  _p2h9minus_gg_VAJHU=          p2h9minus_gg_VAJHU;            
  _p2h9minus_qqbar_VAJHU=       p2h9minus_qqbar_VAJHU;         
  _p2h9minus_prodIndep_VAJHU=   p2h9minus_prodIndep_VAJHU;     
  _p2h10minus_gg_VAJHU=       p2h10minus_gg_VAJHU;       
  _p2h10minus_qqbar_VAJHU=    p2h10minus_qqbar_VAJHU;  
  _p2h10minus_prodIndep_VAJHU=p2h10minus_prodIndep_VAJHU;
  _bkg_VAMCFM=bkg_VAMCFM;
  _bkg_prodIndep_VAMCFM=bkg_prodIndep_VAMCFM;
  _ggzz_VAMCFM=ggzz_VAMCFM;
  _ggzz_p0plus_VAMCFM=ggzz_p0plus_VAMCFM;
  _ggzz_c1_VAMCFM=ggzz_c1_VAMCFM;
  _ggzz_c5_VAMCFM=ggzz_c5_VAMCFM;
  _ggzz_ci_VAMCFM=ggzz_ci_VAMCFM;
  _phjj_VAJHU_old=phjj_VAJHU_old;
  _pvbf_VAJHU_old=pvbf_VAJHU_old;
  _phjj_VAJHU_old_up=phjj_VAJHU_old_up;
  _pvbf_VAJHU_old_up=pvbf_VAJHU_old_up;
  _phjj_VAJHU_old_dn=phjj_VAJHU_old_dn;
  _pvbf_VAJHU_old_dn=pvbf_VAJHU_old_dn;
  _phjj_VAJHU_new=phjj_VAJHU_new;
  _pvbf_VAJHU_new=pvbf_VAJHU_new;
  _phjj_VAJHU_new_up=phjj_VAJHU_new_up;
  _pvbf_VAJHU_new_up=pvbf_VAJHU_new_up;
  _phjj_VAJHU_new_dn=phjj_VAJHU_new_dn;
  _pvbf_VAJHU_new_dn=pvbf_VAJHU_new_dn;
  _p0_g1prime2_VAJHU=p0_g1prime2_VAJHU;
  _pg1g1prime2_VAJHU=pg1g1prime2_VAJHU;
  _Dgg10_VAMCFM=Dgg10_VAMCFM;
  _pg1g4_mela=pg1g4_mela;
  _pg1g4_VAJHU=pg1g4_VAJHU;
  _pg1g4_pi2_VAJHU=pg1g4_pi2_VAJHU;
  _pg1g2_pi2_VAJHU=pg1g2_pi2_VAJHU;
  _pg1g2_mela=pg1g2_mela;
  _pg1g2_VAJHU=pg1g2_VAJHU;
  _pzzzg_VAJHU=    pzzzg_VAJHU;
  _pzzgg_VAJHU=    pzzgg_VAJHU;
  _pzzzg_PS_VAJHU= pzzzg_PS_VAJHU;
  _pzzgg_PS_VAJHU= pzzgg_PS_VAJHU;
  _p0Zgs_VAJHU=    p0Zgs_VAJHU;
  _p0gsgs_VAJHU=   p0gsgs_VAJHU;
  _p0Zgs_PS_VAJHU= p0Zgs_PS_VAJHU;
  _p0gsgs_PS_VAJHU=p0gsgs_PS_VAJHU;
  //_bkg_VAMCFMNorm=bkg_VAMCFMNorm;
}

void HZZ4lNtupleFactory::FillSuperMela(Float_t p0plus_m4l,  // signal m4l probability as in datacards
                                       Float_t bkg_m4l, // backgroun m4l probability as in datacards
                                       Float_t p0plus_m4l_ScaleUp, // signal m4l probability for systematics
                                       Float_t bkg_m4l_ScaleUp, // backgroun m4l probability for systematics
                                       Float_t p0plus_m4l_ScaleDown, // signal m4l probability for systematics
                                       Float_t bkg_m4l_ScaleDown, // backgroun m4l probability for systematics
                                       Float_t p0plus_m4l_ResUp, // signal m4l probability for systematics
                                       Float_t bkg_m4l_ResUp, // backgroun m4l probability for systematics
                                       Float_t p0plus_m4l_ResDown, // signal m4l probability for systematics
                                       Float_t bkg_m4l_ResDown){ // backgroun m4l probability for systematics

  _p0plus_m4l=p0plus_m4l;
  _bkg_m4l=bkg_m4l;
  _p0plus_m4l_ScaleUp=p0plus_m4l_ScaleUp;// signal m4l probability for systematics
  _bkg_m4l_ScaleUp=bkg_m4l_ScaleUp;// backgroun m4l probability for systematics
  _p0plus_m4l_ScaleDown=p0plus_m4l_ScaleDown;// signal m4l probability for systematics
  _bkg_m4l_ScaleDown=bkg_m4l_ScaleDown;// backgroun m4l probability for systematics
  _p0plus_m4l_ResUp=p0plus_m4l_ResUp;// signal m4l probability for systematics
  _bkg_m4l_ResUp=bkg_m4l_ResUp;// backgroun m4l probability for systematics
  _p0plus_m4l_ResDown=p0plus_m4l_ResDown;// signal m4l probability for systematics
  _bkg_m4l_ResDown=bkg_m4l_ResDown;// backgroun m4l probability for systematics
} 


void HZZ4lNtupleFactory::FillHAdditionalInfo(Float_t mZa, Float_t mZb, Float_t mLL4, Float_t mLL6, Float_t SIP4, Float_t iso34)
{
  _ZZmZa=mZa;
  _ZZmZb=mZb;
  _ZZmLL4=mLL4;
  _ZZmLL6=mLL6;
  _ZZSIP4=SIP4;
  //  _ZZiso34=iso34;

  return;
}



void HZZ4lNtupleFactory::FillZInfo(Float_t ZMass, Float_t ZPt, short ZFlav, Float_t Z1MassRefit)
{
  if(!_firstZStored){
    _Z1Mass=ZMass;
    _Z1Pt=ZPt;
    _Z1MassRefit=Z1MassRefit;
    _Z1Flav=ZFlav;
    _firstZStored = true;
  }
  else{
    _Z2Mass=ZMass;
    _Z2Pt=ZPt;
    _Z2Flav=ZFlav;
  }

  return;
}

void HZZ4lNtupleFactory::FillAngularInfo(Float_t costhetastar, Float_t phi, Float_t costheta1, Float_t costheta2, Float_t phistar1, Float_t phistar2,Float_t xi, Float_t xistar)
{
  _costhetastar=costhetastar;
  _phi=phi;
  _costheta1=costheta1;
  _costheta2=costheta2;
  _phistar1=phistar1;
  _phistar2=phistar2;
  _xi=xi;
  _xistar=xistar;

  return;
}

void HZZ4lNtupleFactory::FillLepInfo(Float_t LepPt, Float_t LepEta, Float_t LepPhi, Int_t LepId, Float_t LepSIP, bool isID, float BDT, short parentId, int missingHit)
{
  switch(_LeptonIndex){

  case 1:
    _Lep1Pt=LepPt;
    _Lep1Eta=LepEta;
    _Lep1Phi=LepPhi;
    _Lep1LepId=LepId;
    _Lep1SIP=LepSIP;
    _Lep1isID=isID;
    _Lep1BDT=BDT;
    _Lep1missingHit=(char)missingHit;
    _Lep1ParentId=parentId;
    break;

  case 2:
    _Lep2Pt=LepPt;
    _Lep2Eta=LepEta;
    _Lep2Phi=LepPhi;
    _Lep2LepId=LepId;
    _Lep2SIP=LepSIP;
    _Lep2isID=isID;
    _Lep2BDT=BDT;
    _Lep2missingHit=(char)missingHit;
    _Lep2ParentId=parentId;
    break;
  
  case 3:
    _Lep3Pt=LepPt;
    _Lep3Eta=LepEta;
    _Lep3Phi=LepPhi;
    _Lep3LepId=LepId;
    _Lep3SIP=LepSIP;
    _Lep3isID=isID;
    _Lep3BDT=BDT;
    _Lep3missingHit=(char)missingHit;
    _Lep3ParentId=parentId;
    break;

  case 4:
    _Lep4Pt=LepPt;
    _Lep4Eta=LepEta;
    _Lep4Phi=LepPhi;
    _Lep4LepId=LepId;
    _Lep4SIP=LepSIP;
    _Lep4isID=isID;
    _Lep4BDT=BDT;
    _Lep4missingHit=(char)missingHit;
    _Lep4ParentId=parentId;
    break;

  default:
    std::cout << "Error in indexing the muons! Will abort..." << std::endl;
    assert(0);
  }

  _LeptonIndex++;

  return;
}

void HZZ4lNtupleFactory::FillLepIsolInfo(Float_t LepchargedHadIso, Float_t LepneutralHadIso, Float_t LepphotonIso, Float_t LepcombRelIsoPF)
{

  switch(_LeptonIsoIndex){

  case 1:
    _Lep1chargedHadIso=LepchargedHadIso;
    _Lep1neutralHadIso=LepneutralHadIso;
    _Lep1photonIso=LepphotonIso;
    _Lep1combRelIsoPF=LepcombRelIsoPF;
    break;

  case 2:
    _Lep2chargedHadIso=LepchargedHadIso;
    _Lep2neutralHadIso=LepneutralHadIso;
    _Lep2photonIso=LepphotonIso;
    _Lep2combRelIsoPF=LepcombRelIsoPF;
    break;

  case 3:
    _Lep3chargedHadIso=LepchargedHadIso;
    _Lep3neutralHadIso=LepneutralHadIso;
    _Lep3photonIso=LepphotonIso;
    _Lep3combRelIsoPF=LepcombRelIsoPF;
    break;

  case 4:
    _Lep4chargedHadIso=LepchargedHadIso;
    _Lep4neutralHadIso=LepneutralHadIso;
    _Lep4photonIso=LepphotonIso;
    _Lep4combRelIsoPF=LepcombRelIsoPF;
    break;


  default:
    std::cout << "Error in indexing the muon isolation variables! Will abort..." << std::endl;
    assert(0);
  }

  _LeptonIsoIndex++;

  return;
}

void HZZ4lNtupleFactory::FillPhotonInfo(Float_t PhotPt, Float_t PhotEta, Float_t PhotPhi)
{

  _PhotPt=PhotPt;
  _PhotEta=PhotEta;
  _PhotPhi=PhotPhi;

  return;
}

void HZZ4lNtupleFactory::FillJetInfo(Float_t JetPt, Float_t JetEta, Float_t JetPhi, Float_t JetMass, Float_t JetBTagger, Float_t JetIsBtagged, Float_t JetQGLikelihood, Float_t JetSigma )
{

  _JetPt.push_back(JetPt);
  _JetEta.push_back(JetEta);
  _JetPhi.push_back(JetPhi);
  _JetMass.push_back(JetMass);
  _JetBTagger.push_back(JetBTagger);
  _JetIsBtagged.push_back(JetIsBtagged);
  _JetQGLikelihood.push_back(JetQGLikelihood);
  _JetSigma.push_back(JetSigma);

  return;
}

void HZZ4lNtupleFactory::FillDiJetInfo(Float_t DiJetMass, Float_t DiJetMassPlus, Float_t DiJetMassMinus, Float_t DiJetDEta, Float_t DiJetFisher)
{
  _DiJetMass = DiJetMass;
  _DiJetMassPlus = DiJetMassPlus;
  _DiJetMassMinus = DiJetMassMinus;
  _DiJetDEta = DiJetDEta;
  _DiJetFisher = DiJetFisher;

  return;
}

void HZZ4lNtupleFactory::FillCategorizationInfo(Int_t nExtraLep, Int_t nExtraZ)
{
  _nExtraLep=nExtraLep;
  _nExtraZ=nExtraZ;

  return;
}

void HZZ4lNtupleFactory::FillExtraLepInfo(int extraLeptonIndex, bool extraLeptonExists, const reco::CandidatePtr ExtraLep)
{
  Float_t Pt         = extraLeptonExists ? ExtraLep->pt()    : -9999. ;
  Float_t Eta        = extraLeptonExists ? ExtraLep->eta()   : -9999. ;
  Float_t Phi        = extraLeptonExists ? ExtraLep->phi()   : -9999. ;
  Int_t   LepId      = extraLeptonExists ? ExtraLep->pdgId() :     0  ;
//   Float_t SIP        = extraLeptonExists ? userdatahelpers::getUserFloat(&*ExtraLep,"SIP")              : -9999. ;
//   Bool_t  isID       = extraLeptonExists ? userdatahelpers::getUserFloat(&*ExtraLep,"isID")             :     0  ;
//   Float_t BDT        = extraLeptonExists ? userdatahelpers::getUserFloat(&*ExtraLep,"BDT")              : -9999. ;
//   Char_t  missingHit = extraLeptonExists ? (char)userdatahelpers::getUserFloat(&*ExtraLep,"missingHit") :     0  ;
//   Float_t chargedHadIso = extraLeptonExists ? userdatahelpers::getUserFloat(&*ExtraLep,"PFChargedHadIso") : -9999. ;
//   Float_t neutralHadIso = extraLeptonExists ? userdatahelpers::getUserFloat(&*ExtraLep,"PFNeutralHadIso") : -9999. ;
//   Float_t photonIso     = extraLeptonExists ? userdatahelpers::getUserFloat(&*ExtraLep,"PFPhotonIso")     : -9999. ;
//   Float_t combRelIsoPF  = extraLeptonExists ? userdatahelpers::getUserFloat(&*ExtraLep,"combRelIsoPF")    : -9999. ;

  switch(extraLeptonIndex){

  case 1:
    _ExtraLep1Pt        =Pt;
    _ExtraLep1Eta       =Eta;
    _ExtraLep1Phi       =Phi;
    _ExtraLep1LepId     =LepId;
//     _ExtraLep1SIP       =SIP;
//     _ExtraLep1isID      =isID;
//     _ExtraLep1BDT       =BDT;
//     _ExtraLep1missingHit=missingHit;
//     _ExtraLep1chargedHadIso=chargedHadIso;
//     _ExtraLep1neutralHadIso=neutralHadIso;
//     _ExtraLep1photonIso    =photonIso;
//     _ExtraLep1combRelIsoPF =combRelIsoPF;
    break;

  case 2:
    _ExtraLep2Pt        =Pt;
    _ExtraLep2Eta       =Eta;
    _ExtraLep2Phi       =Phi;
    _ExtraLep2LepId     =LepId;
//     _ExtraLep2SIP       =SIP;
//     _ExtraLep2isID      =isID;
//     _ExtraLep2BDT       =BDT;
//     _ExtraLep2missingHit=missingHit;
//     _ExtraLep2chargedHadIso=chargedHadIso;
//     _ExtraLep2neutralHadIso=neutralHadIso;
//     _ExtraLep2photonIso    =photonIso;
//     _ExtraLep2combRelIsoPF =combRelIsoPF;
    break;

  case 3:
    _ExtraLep3Pt        =Pt;
    _ExtraLep3Eta       =Eta;
    _ExtraLep3Phi       =Phi;
    _ExtraLep3LepId     =LepId;
//     _ExtraLep3SIP       =SIP;
//     _ExtraLep3isID      =isID;
//     _ExtraLep3BDT       =BDT;
//     _ExtraLep3missingHit=missingHit;
//     _ExtraLep3chargedHadIso=chargedHadIso;
//     _ExtraLep3neutralHadIso=neutralHadIso;
//     _ExtraLep3photonIso    =photonIso;
//     _ExtraLep3combRelIsoPF =combRelIsoPF;
    break;

  default:
    std::cout << "Error in indexing the extra leptons ! Will abort..." << std::endl;
    assert(0);
  }
  //Fill the tree once for each extra lepton (no vectors anymore)
  //_outTree->Fill();
  return;
}
