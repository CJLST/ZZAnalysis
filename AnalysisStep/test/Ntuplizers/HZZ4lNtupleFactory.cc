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
  _genFinalState = 0;
  _genProcessId = 0;
  _genHEPMCweight = 0;
  _trigWord =0;
  _genExtInfo =0;

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
  _ZZMass.clear();
  _ZZMassErr.clear();
  _ZZMassErrCorr.clear();
  _ZZMassPreFSR.clear();
  _ZZMassRefit.clear();
  _Chi2KinFit.clear();
  _ZZMassCFit.clear();
  _Chi2CFit.clear();
  _ZZsel.clear();
  _ZZPt.clear();
  _ZZEta.clear();
  _ZZPhi.clear();
  _ZZgenIsSignal.clear();
  _ZZgenIsRightPair.clear();
  _ZZFisher.clear();
  _CRflag.clear();

  //_p0plus_melaNorm.clear();
  //_p0plus_mela.clear();
  //_p0minus_mela.clear();
  //_p0hplus_mela.clear(); // 0h+, analytic distribution
  _p0plus_VAJHU.clear();
  _p0minus_VAJHU.clear();
  _p0plus_VAMCFM.clear();
  _p0hplus_VAJHU.clear(); // 0h+ (high dimensional operator), vector algebra, JHUgen
  //_p1_mela.clear();
  //_p1_prodIndep_mela.clear();
  //_p1plus_mela.clear(); // 1+, analytic distribution 
  //_p1plus_prodIndep_mela.clear(); // 1+, analytic distribution 
  _p1_VAJHU.clear();
  _p1_prodIndep_VAJHU.clear();
  _p1plus_VAJHU.clear(); // 1+ (axial vector), vector algebra, JHUgen,
  _p1plus_prodIndep_VAJHU.clear(); // 1+ (axial vector), vector algebra, JHUgen,
  //_p2_mela.clear();
  //_p2_prodIndep_mela.clear();
  //_p2qqb_mela.clear(); // graviton produced by qqbar vector algebra, analytical,
  //_p2hplus_mela.clear(); // graviton produced by qqbar vector algebra, analytical,
  //_p2hminus_mela.clear(); // graviton produced by qqbar vector algebra, analytical,
  //_p2bplus_mela.clear(); // graviton produced by qqbar vector algebra, analytical,
  _p2_VAJHU.clear();
  _p2_prodIndep_VAJHU.clear();
  _p2qqb_VAJHU.clear();
  _p2hplus_VAJHU.clear();
  _p2hminus_VAJHU.clear();
  _p2bplus_VAJHU.clear();
	_p2hplus_qqb_VAJHU.clear();					
	_p2hplus_prodIndep_VAJHU.clear();		
	_p2hminus_qqb_VAJHU.clear();				
	_p2hminus_prodIndep_VAJHU.clear();	
	_p2bplus_qqb_VAJHU.clear();					
	_p2bplus_prodIndep_VAJHU.clear();		
	_p2h2plus_gg_VAJHU.clear();      		                          
	_p2h2plus_qqbar_VAJHU.clear();   		
	_p2h2plus_prodIndep_VAJHU.clear();	
	_p2h3plus_gg_VAJHU.clear();       	
	_p2h3plus_qqbar_VAJHU.clear();    	
	_p2h3plus_prodIndep_VAJHU.clear();	
	_p2h6plus_gg_VAJHU.clear();       	
	_p2h6plus_qqbar_VAJHU.clear();    	
	_p2h6plus_prodIndep_VAJHU.clear();	
	_p2h7plus_gg_VAJHU.clear();	
	_p2h7plus_qqbar_VAJHU.clear();    	
	_p2h7plus_prodIndep_VAJHU.clear();	
	_p2h9minus_gg_VAJHU.clear();       	
	_p2h9minus_qqbar_VAJHU.clear();    	
	_p2h9minus_prodIndep_VAJHU.clear();	
	_p2h10minus_gg_VAJHU.clear();       
	_p2h10minus_qqbar_VAJHU.clear();  
	_p2h10minus_prodIndep_VAJHU.clear();
  //_bkg_mela.push_back(bkg_mela);
  //_bkg_mela.clear();
  _bkg_VAMCFM.clear();
  _bkg_prodIndep_VAMCFM.clear();
  _ggzz_VAMCFM.clear();
  _ggzz_p0plus_VAMCFM.clear();
  _ggzz_c1_VAMCFM.clear();
  _ggzz_c5_VAMCFM.clear();
  _ggzz_ci_VAMCFM.clear();
  //_bkg_VAMCFMNorm.clear();
  //_p0_pt.clear();
  //_p0_y.clear();
  //_bkg_pt.clear();
  //_bkg_y.clear();  
  _p0plus_m4l.clear(); 
  _bkg_m4l.clear(); 
  _pg1g4_mela.clear();
  _pg1g4_VAJHU.clear();
  _pg1g4_pi2_VAJHU.clear();
  _pg1g2_pi2_VAJHU.clear();
  _pg1g2_mela.clear();
  _pg1g2_VAJHU.clear();
  _p0plus_m4l.clear();// signal m4l probability as in datacards
  _bkg_m4l.clear();// backgroun m4l probability as in datacards
  _p0plus_m4l_ScaleUp.clear();// signal m4l probability for systematics
  _bkg_m4l_ScaleUp.clear();// backgroun m4l probability for systematics
  _p0plus_m4l_ScaleDown.clear();// signal m4l probability for systematics
  _bkg_m4l_ScaleDown.clear();// backgroun m4l probability for systematics
  _p0plus_m4l_ResUp.clear();// signal m4l probability for systematics
  _bkg_m4l_ResUp.clear();// backgroun m4l probability for systematics
  _p0plus_m4l_ResDown.clear();// signal m4l probability for systematics
  _bkg_m4l_ResDown.clear();     // backgroun m4l probability for systematics
  _phjj_VAJHU_old.clear();
  _pvbf_VAJHU_old.clear();
  _phjj_VAJHU_old_up.clear();
  _pvbf_VAJHU_old_up.clear();
  _phjj_VAJHU_old_dn.clear();
  _pvbf_VAJHU_old_dn.clear();
  _phjj_VAJHU_new.clear();
  _pvbf_VAJHU_new.clear();
  _phjj_VAJHU_new_up.clear();
  _pvbf_VAJHU_new_up.clear();
  _phjj_VAJHU_new_dn.clear();
  _pvbf_VAJHU_new_dn.clear();
	_p0_g1prime2_VAJHU.clear();
	_pg1g1prime2_VAJHU.clear();
	_Dgg10_VAMCFM.clear();

  _pzzzg_VAJHU.clear();
  _pzzgg_VAJHU.clear();
  _pzzzg_PS_VAJHU.clear();
  _pzzgg_PS_VAJHU.clear();
  _p0Zgs_VAJHU.clear();
  _p0gsgs_VAJHU.clear();
  _p0Zgs_PS_VAJHU.clear();
  _p0gsgs_PS_VAJHU.clear();

  _ZZmZa.clear();
  _ZZmZb.clear();
  _ZZmLL4.clear();
  _ZZmLL6.clear();
  _ZZSIP4.clear();
  //  _ZZiso34.clear();

  //Z1 variables
  _Z1Mass.clear();
  _Z1Pt.clear();
  _Z1MassRefit.clear();

  //Z2 variables
  _Z2Mass.clear();
  _Z2Pt.clear();

  //Angular variables
  _costhetastar.clear();
  _phi.clear();
//   _helphiZ1.clear();
//   _helphiZ2.clear();
  _costheta1.clear();
  _costheta2.clear();
  _phistar1.clear();
  _phistar2.clear();
  _xi.clear();
  _xistar.clear();

  //Lepton variables
  _Lep1Pt.clear();
  _Lep1Eta.clear();
  _Lep1Phi.clear();
  _Lep1LepId.clear();
  _Lep1SIP.clear();
  _Lep1isID.clear();
  _Lep1BDT.clear();
  _Lep1missingHit.clear();
  _Lep1ParentId.clear();

  _Lep2Pt.clear();
  _Lep2Eta.clear();
  _Lep2Phi.clear();
  _Lep2LepId.clear();
  _Lep2SIP.clear();
  _Lep2isID.clear();
  _Lep2BDT.clear();
  _Lep2missingHit.clear();
  _Lep2ParentId.clear();

  _Lep3Pt.clear();
  _Lep3Eta.clear();
  _Lep3Phi.clear();
  _Lep3LepId.clear();
  _Lep3SIP.clear();
  _Lep3isID.clear();
  _Lep3BDT.clear();
  _Lep3missingHit.clear();
  _Lep3ParentId.clear();

  _Lep4Pt.clear();
  _Lep4Eta.clear();
  _Lep4Phi.clear();
  _Lep4LepId.clear();
  _Lep4SIP.clear();
  _Lep4isID.clear();
  _Lep4BDT.clear();
  _Lep4missingHit.clear();
  _Lep4ParentId.clear();

  //Lepton isolation variables
  _Lep1chargedHadIso.clear();
  _Lep1neutralHadIso.clear();
  _Lep1photonIso.clear();
  _Lep1combRelIsoPF.clear();

  _Lep2chargedHadIso.clear();
  _Lep2neutralHadIso.clear();
  _Lep2photonIso.clear();
  _Lep2combRelIsoPF.clear();

  _Lep3chargedHadIso.clear();
  _Lep3neutralHadIso.clear();
  _Lep3photonIso.clear();
  _Lep3combRelIsoPF.clear();

  _Lep4chargedHadIso.clear();
  _Lep4neutralHadIso.clear();
  _Lep4photonIso.clear();
  _Lep4combRelIsoPF.clear();

  //Photon variables
  _PhotPt.clear();
  _PhotEta.clear();
  _PhotPhi.clear(); 

  //Jet variables
  _JetPt.clear();
  _JetEta.clear();
  _JetPhi.clear(); 
  _JetMass.clear(); 
  _JetBTag.clear();
  _JetSigma.clear();

  _DiJetMass=-99;
  _DiJetMassPlus=-99;
  _DiJetMassMinus=-99;
  _DiJetDEta=-99;
  
  //Categorization-related variables
  _nExtraLep.clear();
  _nExtraZ.clear();
  _nJets.clear();
  _nCleanedJets.clear();
  _nCleanedJetsPt30.clear();

  //Variables of extra leptons
  _ExtraLep1Pt.clear();
  _ExtraLep1Eta.clear();
  _ExtraLep1Phi.clear();
  _ExtraLep1LepId.clear();
  _ExtraLep1SIP.clear();
  _ExtraLep1isID.clear();
  _ExtraLep1BDT.clear();
  _ExtraLep1missingHit.clear();
  _ExtraLep1chargedHadIso.clear();
  _ExtraLep1neutralHadIso.clear();
  _ExtraLep1photonIso.clear();
  _ExtraLep1combRelIsoPF.clear();

  _ExtraLep2Pt.clear();
  _ExtraLep2Eta.clear();
  _ExtraLep2Phi.clear();
  _ExtraLep2LepId.clear();
  _ExtraLep2SIP.clear();
  _ExtraLep2isID.clear();
  _ExtraLep2BDT.clear();
  _ExtraLep2missingHit.clear();
  _ExtraLep2chargedHadIso.clear();
  _ExtraLep2neutralHadIso.clear();
  _ExtraLep2photonIso.clear();
  _ExtraLep2combRelIsoPF.clear();

  _ExtraLep3Pt.clear();
  _ExtraLep3Eta.clear();
  _ExtraLep3Phi.clear();
  _ExtraLep3LepId.clear();
  _ExtraLep3SIP.clear();
  _ExtraLep3isID.clear();
  _ExtraLep3BDT.clear();
  _ExtraLep3missingHit.clear();
  _ExtraLep3chargedHadIso.clear();
  _ExtraLep3neutralHadIso.clear();
  _ExtraLep3photonIso.clear();
  _ExtraLep3combRelIsoPF.clear();


  return;
}

void HZZ4lNtupleFactory::InitializeBranches()
{
  //Event variables
  _outTree->Branch("RunNumber",&_RunNumber,"RunNumber/I");
  _outTree->Branch("EventNumber",&_EventNumber,"EventNumber/L");
  _outTree->Branch("LumiNumber",&_LumiNumber,"LumiNumber/I");
  _outTree->Branch("iBC",&_IndexBestCand,"iBC/I");
  //  _outTree->Branch("Nmu",&_Nmu,"Nmu/I");
  //  _outTree->Branch("Nele",&_Nele,"Nele/I");
  _outTree->Branch("Nvtx",&_Nvtx,"Nvtx/I");
  _outTree->Branch("NObsInt",&_NObsInt,"NObsInt/I");
  _outTree->Branch("NTrueInt",&_NTrueInt,"NTrueInt/F");
  _outTree->Branch("PUWeight12",&_PUWeight,"PUWeight12/F");
  _outTree->Branch("PFMET",&_PFMET,"PFMET/F");
  _outTree->Branch("genFinalState",&_genFinalState,"genFinalState/I");
  _outTree->Branch("genProcessId",&_genProcessId,"genProcessId/I");
  _outTree->Branch("genHEPMCweight",&_genHEPMCweight,"genHEPMCweight/F");
  _outTree->Branch("trigWord",&_trigWord,"trigWord/S");
  _outTree->Branch("genExtInfo",&_genExtInfo,"genExtInfo/S");

  //H variables
  _outTree->Branch("ZZMass",&_ZZMass);
  _outTree->Branch("ZZMassErr",&_ZZMassErr);
  _outTree->Branch("ZZMassErrCorr",&_ZZMassErrCorr);
  _outTree->Branch("ZZMassPreFSR",&_ZZMassPreFSR);
  _outTree->Branch("ZZsel",&_ZZsel);
  _outTree->Branch("ZZPt",&_ZZPt);
  _outTree->Branch("ZZEta",&_ZZEta);
  _outTree->Branch("ZZPhi",&_ZZPhi);
  //--> Commented as these variables are not computed at the moment
  //  _outTree->Branch("ZZgenIsSignal",&_ZZgenIsSignal);
  //  _outTree->Branch("ZZgenIsRightPair",&_ZZgenIsRightPair);
  _outTree->Branch("ZZFisher",&_ZZFisher);
  _outTree->Branch("CRflag",&_CRflag);

  //_outTree->Branch("p0plus_melaNorm",&_p0plus_melaNorm);
  //_outTree->Branch("p0plus_mela",&_p0plus_mela);
  //_outTree->Branch("p0minus_mela",&_p0minus_mela);
  //_outTree->Branch("p0hplus_mela",&_p0hplus_mela); // 0h+, analytic distribution
  _outTree->Branch("p0plus_VAJHU",&_p0plus_VAJHU);
  _outTree->Branch("p0minus_VAJHU",&_p0minus_VAJHU);
  _outTree->Branch("p0plus_VAMCFM",&_p0plus_VAMCFM);
  _outTree->Branch("p0hplus_VAJHU",&_p0hplus_VAJHU); // 0h+ (high dimensional operator), vector algebra, JHUgen
  //_outTree->Branch("p1_mela",&_p1_mela);
  //_outTree->Branch("p1_prodIndep_mela",&_p1_prodIndep_mela);
  //_outTree->Branch("p1plus_mela",&_p1plus_mela); // 1+, analytic distribution 
  //_outTree->Branch("p1plus_prodIndep_mela",&_p1plus_prodIndep_mela); // 1+, analytic distribution 
  _outTree->Branch("p1_VAJHU",&_p1_VAJHU);
  _outTree->Branch("p1_prodIndep_VAJHU",&_p1_prodIndep_VAJHU);
  _outTree->Branch("p1plus_VAJHU",&_p1plus_VAJHU); // 1+ (axial vector), vector algebra, JHUgen,
  _outTree->Branch("p1plus_prodIndep_VAJHU",&_p1plus_prodIndep_VAJHU); // 1+ (axial vector), vector algebra, JHUgen,
  //_outTree->Branch("p2_mela",&_p2_mela);
  //_outTree->Branch("p2_prodIndep_mela",&_p2_prodIndep_mela);
  //_outTree->Branch("p2qqb_mela",&_p2qqb_mela); // graviton produced by qqbar vector algebra, analytical,
  //_outTree->Branch("p2hplus_mela",&_p2hplus_mela); // graviton produced by qqbar vector algebra, analytical,
  //_outTree->Branch("p2hminus_mela",&_p2hminus_mela); // graviton produced by qqbar vector algebra, analytical,
  //_outTree->Branch("p2bplus_mela",&_p2bplus_mela); // graviton produced by qqbar vector algebra, analytical,
  _outTree->Branch("p2_VAJHU",&_p2_VAJHU);
  _outTree->Branch("p2_prodIndep_VAJHU",&_p2_prodIndep_VAJHU);
  _outTree->Branch("p2qqb_VAJHU",&_p2qqb_VAJHU);
  _outTree->Branch("p2hplus_VAJHU",&_p2hplus_VAJHU);
  _outTree->Branch("p2hminus_VAJHU",&_p2hminus_VAJHU);
  _outTree->Branch("p2bplus_VAJHU",&_p2bplus_VAJHU);

	_outTree->Branch("p2hplus_qqb_VAJHU"					,		&_p2hplus_qqb_VAJHU);					
	_outTree->Branch("p2hplus_prodIndep_VAJHU"		,		&_p2hplus_prodIndep_VAJHU);		
	_outTree->Branch("p2hminus_qqb_VAJHU"				,  		&_p2hminus_qqb_VAJHU);				
	_outTree->Branch("p2hminus_prodIndep_VAJHU"	,  		&_p2hminus_prodIndep_VAJHU);	
	_outTree->Branch("p2bplus_qqb_VAJHU"					,		&_p2bplus_qqb_VAJHU);					
	_outTree->Branch("p2bplus_prodIndep_VAJHU"		,		&_p2bplus_prodIndep_VAJHU);		
	_outTree->Branch("p2h2plus_gg_VAJHU"      		,		&_p2h2plus_gg_VAJHU);      		                          
	_outTree->Branch("p2h2plus_qqbar_VAJHU"   		,		&_p2h2plus_qqbar_VAJHU);   		
	_outTree->Branch("p2h2plus_prodIndep_VAJHU"	,   	&_p2h2plus_prodIndep_VAJHU);	
	_outTree->Branch("p2h3plus_gg_VAJHU"       	,   	&_p2h3plus_gg_VAJHU);       	
	_outTree->Branch("p2h3plus_qqbar_VAJHU"    	,   	&_p2h3plus_qqbar_VAJHU);    	
	_outTree->Branch("p2h3plus_prodIndep_VAJHU"	,   	&_p2h3plus_prodIndep_VAJHU);	
	_outTree->Branch("p2h6plus_gg_VAJHU"       	,   	&_p2h6plus_gg_VAJHU);       	
	_outTree->Branch("p2h6plus_qqbar_VAJHU"    	,   	&_p2h6plus_qqbar_VAJHU);    	
	_outTree->Branch("p2h6plus_prodIndep_VAJHU"	,   	&_p2h6plus_prodIndep_VAJHU);	
	_outTree->Branch("p2h7plus_gg_VAJHU"       	,   	&_p2h7plus_gg_VAJHU);	
	_outTree->Branch("p2h7plus_qqbar_VAJHU"    	,   	&_p2h7plus_qqbar_VAJHU);    	
	_outTree->Branch("p2h7plus_prodIndep_VAJHU"	,   	&_p2h7plus_prodIndep_VAJHU);	
	_outTree->Branch("p2h9minus_gg_VAJHU"       	,		&_p2h9minus_gg_VAJHU);       	
	_outTree->Branch("p2h9minus_qqbar_VAJHU"    	,		&_p2h9minus_qqbar_VAJHU);    	
	_outTree->Branch("p2h9minus_prodIndep_VAJHU"	,		&_p2h9minus_prodIndep_VAJHU);	
	_outTree->Branch("p2h10minus_gg_VAJHU"       , 		&_p2h10minus_gg_VAJHU);       
	_outTree->Branch("p2h10minus_qqbar_VAJHU"    , 		&_p2h10minus_qqbar_VAJHU);  
	_outTree->Branch("p2h10minus_prodIndep_VAJHU", 		&_p2h10minus_prodIndep_VAJHU);
  //_outTree->Branch("bkg_mela",&_bkg_mela);
  _outTree->Branch("bkg_VAMCFM",&_bkg_VAMCFM);
  _outTree->Branch("bkg_prodIndep_VAMCFM",&_bkg_prodIndep_VAMCFM);
  _outTree->Branch("ggzz_VAMCFM",&_ggzz_VAMCFM);
  _outTree->Branch("ggzz_p0plus_VAMCFM",&_ggzz_p0plus_VAMCFM);
  _outTree->Branch("ggzz_c1_VAMCFM",&_ggzz_c1_VAMCFM);
  _outTree->Branch("ggzz_c5_VAMCFM",&_ggzz_c5_VAMCFM);
  _outTree->Branch("ggzz_ci_VAMCFM",&_ggzz_ci_VAMCFM);
  //_outTree->Branch("bkg_VAMCFMNorm",&_bkg_VAMCFMNorm);
  //_outTree->Branch("p0_pt",&_p0_pt);
  //_outTree->Branch("p0_y",&_p0_y);
  //_outTree->Branch("bkg_pt",&_bkg_pt);
  //_outTree->Branch("bkg_y",&_bkg_y);  
  
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

  // Varialbes removed - they can be retrieved from the lepton info.
  //  _outTree->Branch("ZZmZa",&_ZZmZa);
  //  _outTree->Branch("ZZmZb",&_ZZmZb);
  //  _outTree->Branch("ZZmLL4",&_ZZmLL4);
  //  _outTree->Branch("ZZmLL6",&_ZZmLL6);
  //  _outTree->Branch("ZZSIP4",&_ZZSIP4);
  //  _outTree->Branch("ZZiso34",&_ZZiso34);

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
  //Z1 variables
  _outTree->Branch("Z1Mass",&_Z1Mass);
  _outTree->Branch("Z1Pt",&_Z1Pt);

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

  //Angular variables
  _outTree->Branch("costhetastar",&_costhetastar);
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
  
  //Jet variables
  _outTree->Branch("JetPt",&_JetPt);
  _outTree->Branch("JetEta",&_JetEta);
  _outTree->Branch("JetPhi",&_JetPhi);
  _outTree->Branch("JetMass",&_JetMass);
  _outTree->Branch("JetBTag",&_JetBTag);
  _outTree->Branch("JetSigma",&_JetSigma);
  _outTree->Branch("DiJetMass",&_DiJetMass,"DiJetMass/F");
  _outTree->Branch("DiJetMassPlus",&_DiJetMassPlus,"DiJetMassPlus/F");
  _outTree->Branch("DiJetMassMinus",&_DiJetMassMinus,"DiJetMassMinus/F");
  _outTree->Branch("DiJetDEta",&_DiJetDEta,"DiJetDEta/F");

  //Categorization-related variables
  _outTree->Branch("nExtraLep",&_nExtraLep);
  _outTree->Branch("nExtraZ",&_nExtraZ);
  _outTree->Branch("nJets",&_nJets);
  _outTree->Branch("nCleanedJets",&_nCleanedJets);
  _outTree->Branch("nCleanedJetsPt30",&_nCleanedJetsPt30);

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

  //Generated particles
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

void HZZ4lNtupleFactory::FillHGenInfo(const math::XYZTLorentzVector pH)
{
  _genHMass = pH.M();
  _genHPt = pH.Pt();

  return;
}

void HZZ4lNtupleFactory::FillZGenInfo(const math::XYZTLorentzVector pZ1, const math::XYZTLorentzVector pZ2)
{
  _genZ1Mass = pZ1.M();
  _genZ1Pt = pZ1.Pt();

  _genZ2Mass = pZ2.M();
  _genZ2Pt = pZ2.Pt();

  return;
}

void HZZ4lNtupleFactory::FillLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, Short_t Lep3Id, Short_t Lep4Id, 
					const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2, const math::XYZTLorentzVector Lep3, const math::XYZTLorentzVector Lep4)
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

void HZZ4lNtupleFactory::FillEventInfo(const Int_t RunNumber, const Long64_t EventNumber, const Int_t LumiNumber, const Int_t IndexBestCand, Int_t Nvtx, 
				       Int_t NObsInt, Float_t NTrueInt, Float_t PUweight, const Float_t PFMET, Int_t genFinalState, Int_t genProcessId, Float_t genHEPMCweight, Short_t trigWord,  Short_t genExtInfo)
{
  _RunNumber = RunNumber;
  _EventNumber = EventNumber;
  _LumiNumber = LumiNumber;
  _IndexBestCand = IndexBestCand;
  _Nvtx = Nvtx;
  _NObsInt =NObsInt;
  _NTrueInt =NTrueInt;
  _PUWeight =PUweight;

  _PFMET = PFMET;
  _genFinalState = genFinalState;
  _genProcessId = genProcessId;
  _genHEPMCweight = genHEPMCweight;
  _trigWord = trigWord;
  _genExtInfo = genExtInfo;
  return;
}

void HZZ4lNtupleFactory::FillHInfo(const Float_t ZZMass, const Float_t ZZMassErr, const Float_t ZZMassErrCorr, const Float_t ZZMassPreFSR, const Float_t ZZMassRefit, const Float_t Chi2KinFit, const Float_t ZZMassCFit, const Float_t Chi2CFit, const Int_t ZZsel, const Float_t ZZPt, const Float_t ZZEta, const Float_t ZZPhi, const Int_t isSignal, const Int_t isRightPair, const Float_t ZZFisher, const Int_t CRflag)
{
  _ZZMass.push_back(ZZMass);
  _ZZMassErr.push_back(ZZMassErr);
  _ZZMassErrCorr.push_back(ZZMassErrCorr);
  _ZZMassPreFSR.push_back(ZZMassPreFSR);
  _ZZMassRefit.push_back(ZZMassRefit);
  _Chi2KinFit.push_back(Chi2KinFit);
  _ZZMassCFit.push_back(ZZMassCFit);
  _Chi2CFit.push_back(Chi2CFit);
  _ZZsel.push_back(ZZsel);
  _ZZPt.push_back(ZZPt);
  _ZZEta.push_back(ZZEta);
  _ZZPhi.push_back(ZZPhi);
  _ZZgenIsSignal.push_back(isSignal);
  _ZZgenIsRightPair.push_back(isRightPair);
  _ZZFisher.push_back(ZZFisher);
  _CRflag.push_back(CRflag);

  return;
}

void HZZ4lNtupleFactory::FillProbability(  //const Float_t p0plus_melaNorm,
					   //const Float_t p0plus_mela,
					   //const Float_t p0minus_mela,
					   //const Float_t p0hplus_mela, // 0h+, analytic distribution
					   const Float_t p0plus_VAJHU,
					   const Float_t p0minus_VAJHU,
					   const Float_t p0plus_VAMCFM,
					   const Float_t p0hplus_VAJHU, // 0h+ (high dimensional operator), vector algebra, JHUgen
					   //const Float_t p1_mela,
					   //const Float_t p1_prodIndep_mela,
					   //const Float_t p1plus_mela, // 1+, analytic distribution 
					   //const Float_t p1plus_prodIndep_mela, // 1+, analytic distribution 
					   const Float_t p1_VAJHU,
					   const Float_t p1_prodIndep_VAJHU,
					   const Float_t p1plus_VAJHU, // 1+ (axial vector), vector algebra, JHUgen,
					   const Float_t p1plus_prodIndep_VAJHU, // 1+ (axial vector), vector algebra, JHUgen,
					   //const Float_t p2_mela ,
					   //const Float_t p2_prodIndep_mela ,
					   //const Float_t p2qqb_mela, // graviton produced by qqbar vector algebra, analytical,
					   //const Float_t p2hplus_mela, // graviton produced by qqbar vector algebra, analytical,
					   //const Float_t p2hminus_mela, // graviton produced by qqbar vector algebra, analytical,
					   //const Float_t p2bplus_mela, // graviton produced by qqbar vector algebra, analytical,
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
					   ///const Float_t bkg_mela,
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
					   ){
  //_p0plus_melaNorm.push_back(p0plus_melaNorm);
  //_p0plus_mela.push_back(p0plus_mela);
  //_p0minus_mela.push_back(p0minus_mela);
  //_p0hplus_mela.push_back(p0hplus_mela);// 0h+, analytic distribution
  _p0plus_VAJHU.push_back(p0plus_VAJHU);
  _p0minus_VAJHU.push_back(p0minus_VAJHU);
  _p0plus_VAMCFM.push_back(p0plus_VAMCFM);
  _p0hplus_VAJHU.push_back(p0hplus_VAJHU);// 0h+ (high dimensional operator), vector algebra, JHUgen
  //_p1_mela.push_back(p1_mela);
  //_p1_prodIndep_mela.push_back(p1_prodIndep_mela);
  //_p1plus_mela.push_back(p1plus_mela);// 1+, analytic distribution 
  //_p1plus_prodIndep_mela.push_back(p1plus_prodIndep_mela);// 1+, analytic distribution 
  _p1_VAJHU.push_back(p1_VAJHU);
  _p1_prodIndep_VAJHU.push_back(p1_prodIndep_VAJHU);
  _p1plus_VAJHU.push_back(p1plus_VAJHU);// 1+ (axial vector), vector algebra, JHUgen,
  _p1plus_prodIndep_VAJHU.push_back(p1plus_prodIndep_VAJHU);// 1+ (axial vector), vector algebra, JHUgen,
  //_p2_mela .push_back(p2_mela );
  //_p2_prodIndep_mela .push_back(p2_prodIndep_mela );
  //_p2qqb_mela.push_back(p2qqb_mela);// graviton produced by qqbar vector algebra, analytical,
  //_p2hplus_mela.push_back(p2hplus_mela);// graviton produced by qqbar vector algebra, analytical,
  //_p2hminus_mela.push_back(p2hminus_mela);// graviton produced by qqbar vector algebra, analytical,
  //_p2bplus_mela.push_back(p2bplus_mela);// graviton produced by qqbar vector algebra, analytical,
  _p2_VAJHU.push_back(p2_VAJHU);
  _p2_prodIndep_VAJHU.push_back(p2_prodIndep_VAJHU);
  _p2qqb_VAJHU.push_back(p2qqb_VAJHU);
  _p2hplus_VAJHU.push_back(p2hplus_VAJHU);
  _p2hminus_VAJHU.push_back(p2hminus_VAJHU);
  _p2bplus_VAJHU.push_back(p2bplus_VAJHU);
	_p2hplus_qqb_VAJHU.push_back(					p2hplus_qqb_VAJHU);					
	_p2hplus_prodIndep_VAJHU.push_back(		p2hplus_prodIndep_VAJHU);		
	_p2hminus_qqb_VAJHU.push_back(				p2hminus_qqb_VAJHU);				
	_p2hminus_prodIndep_VAJHU.push_back(	p2hminus_prodIndep_VAJHU);	
	_p2bplus_qqb_VAJHU.push_back(					p2bplus_qqb_VAJHU);					
	_p2bplus_prodIndep_VAJHU.push_back(		p2bplus_prodIndep_VAJHU);		
	_p2h2plus_gg_VAJHU.push_back(      		p2h2plus_gg_VAJHU);      		                          
	_p2h2plus_qqbar_VAJHU.push_back(   		p2h2plus_qqbar_VAJHU);   		
	_p2h2plus_prodIndep_VAJHU.push_back(	p2h2plus_prodIndep_VAJHU);	
	_p2h3plus_gg_VAJHU.push_back(       	p2h3plus_gg_VAJHU);       	
	_p2h3plus_qqbar_VAJHU.push_back(    	p2h3plus_qqbar_VAJHU);    	
	_p2h3plus_prodIndep_VAJHU.push_back(	p2h3plus_prodIndep_VAJHU);	
	_p2h6plus_gg_VAJHU.push_back(       	p2h6plus_gg_VAJHU);       	
	_p2h6plus_qqbar_VAJHU.push_back(    	p2h6plus_qqbar_VAJHU);    	
	_p2h6plus_prodIndep_VAJHU.push_back(	p2h6plus_prodIndep_VAJHU);	
	_p2h7plus_gg_VAJHU.push_back(       	p2h7plus_gg_VAJHU);	
	_p2h7plus_qqbar_VAJHU.push_back(    	p2h7plus_qqbar_VAJHU);    	
	_p2h7plus_prodIndep_VAJHU.push_back(	p2h7plus_prodIndep_VAJHU);	
	_p2h9minus_gg_VAJHU.push_back(       	p2h9minus_gg_VAJHU);       	
	_p2h9minus_qqbar_VAJHU.push_back(    	p2h9minus_qqbar_VAJHU);    	
	_p2h9minus_prodIndep_VAJHU.push_back(	p2h9minus_prodIndep_VAJHU);	
	_p2h10minus_gg_VAJHU.push_back(       p2h10minus_gg_VAJHU);       
	_p2h10minus_qqbar_VAJHU.push_back(    p2h10minus_qqbar_VAJHU);  
	_p2h10minus_prodIndep_VAJHU.push_back(p2h10minus_prodIndep_VAJHU);
  //_bkg_mela.push_back(bkg_mela);
  _bkg_VAMCFM.push_back(bkg_VAMCFM);
  _bkg_prodIndep_VAMCFM.push_back(bkg_prodIndep_VAMCFM);
  _ggzz_VAMCFM.push_back(ggzz_VAMCFM);
  _ggzz_p0plus_VAMCFM.push_back(ggzz_p0plus_VAMCFM);
  _ggzz_c1_VAMCFM.push_back(ggzz_c1_VAMCFM);
  _ggzz_c5_VAMCFM.push_back(ggzz_c5_VAMCFM);
  _ggzz_ci_VAMCFM.push_back(ggzz_ci_VAMCFM);
  _phjj_VAJHU_old.push_back(phjj_VAJHU_old);
  _pvbf_VAJHU_old.push_back(pvbf_VAJHU_old);
  _phjj_VAJHU_old_up.push_back(phjj_VAJHU_old_up);
  _pvbf_VAJHU_old_up.push_back(pvbf_VAJHU_old_up);
  _phjj_VAJHU_old_dn.push_back(phjj_VAJHU_old_dn);
  _pvbf_VAJHU_old_dn.push_back(pvbf_VAJHU_old_dn);
  _phjj_VAJHU_new.push_back(phjj_VAJHU_new);
  _pvbf_VAJHU_new.push_back(pvbf_VAJHU_new);
  _phjj_VAJHU_new_up.push_back(phjj_VAJHU_new_up);
  _pvbf_VAJHU_new_up.push_back(pvbf_VAJHU_new_up);
  _phjj_VAJHU_new_dn.push_back(phjj_VAJHU_new_dn);
  _pvbf_VAJHU_new_dn.push_back(pvbf_VAJHU_new_dn);
  _p0_g1prime2_VAJHU.push_back(p0_g1prime2_VAJHU);
  _pg1g1prime2_VAJHU.push_back(pg1g1prime2_VAJHU);
  _Dgg10_VAMCFM.push_back(Dgg10_VAMCFM);
  _pg1g4_mela.push_back(pg1g4_mela);
  _pg1g4_VAJHU.push_back(pg1g4_VAJHU);
  _pg1g4_pi2_VAJHU.push_back(pg1g4_pi2_VAJHU);
  _pg1g2_pi2_VAJHU.push_back(pg1g2_pi2_VAJHU);
  _pg1g2_mela.push_back(pg1g2_mela);
  _pg1g2_VAJHU.push_back(pg1g2_VAJHU);
  _pzzzg_VAJHU.push_back(    pzzzg_VAJHU);
  _pzzgg_VAJHU.push_back(    pzzgg_VAJHU);
  _pzzzg_PS_VAJHU.push_back( pzzzg_PS_VAJHU);
  _pzzgg_PS_VAJHU.push_back( pzzgg_PS_VAJHU);
  _p0Zgs_VAJHU.push_back(    p0Zgs_VAJHU);
  _p0gsgs_VAJHU.push_back(   p0gsgs_VAJHU);
  _p0Zgs_PS_VAJHU.push_back( p0Zgs_PS_VAJHU);
  _p0gsgs_PS_VAJHU.push_back(p0gsgs_PS_VAJHU);
  //_bkg_VAMCFMNorm.push_back(bkg_VAMCFMNorm);
  //_p0_pt.push_back(p0_pt);
  //_p0_y.push_back(p0_y);
  //_bkg_pt.push_back(bkg_pt);
  //_bkg_y.push_back(bkg_y);  

  //_ZZVAKD.push_back( p0plus_VAJHU/( bkg_VAMCFMNorm + p0plus_VAJHU)  );
}

void HZZ4lNtupleFactory::FillSuperMela(	const Float_t p0plus_m4l,  // signal m4l probability as in datacards
					const Float_t bkg_m4l, // backgroun m4l probability as in datacards
					const Float_t p0plus_m4l_ScaleUp, // signal m4l probability for systematics
					const Float_t bkg_m4l_ScaleUp, // backgroun m4l probability for systematics
					const Float_t p0plus_m4l_ScaleDown, // signal m4l probability for systematics
					const Float_t bkg_m4l_ScaleDown, // backgroun m4l probability for systematics
					const Float_t p0plus_m4l_ResUp, // signal m4l probability for systematics
					const Float_t bkg_m4l_ResUp, // backgroun m4l probability for systematics
					const Float_t p0plus_m4l_ResDown, // signal m4l probability for systematics
					const Float_t bkg_m4l_ResDown){ // backgroun m4l probability for systematics

  _p0plus_m4l.push_back(p0plus_m4l);
  _bkg_m4l.push_back(bkg_m4l);
  _p0plus_m4l_ScaleUp.push_back(p0plus_m4l_ScaleUp);// signal m4l probability for systematics
  _bkg_m4l_ScaleUp.push_back(bkg_m4l_ScaleUp);// backgroun m4l probability for systematics
  _p0plus_m4l_ScaleDown.push_back(p0plus_m4l_ScaleDown);// signal m4l probability for systematics
  _bkg_m4l_ScaleDown.push_back(bkg_m4l_ScaleDown);// backgroun m4l probability for systematics
  _p0plus_m4l_ResUp.push_back(p0plus_m4l_ResUp);// signal m4l probability for systematics
  _bkg_m4l_ResUp.push_back(bkg_m4l_ResUp);// backgroun m4l probability for systematics
  _p0plus_m4l_ResDown.push_back(p0plus_m4l_ResDown);// signal m4l probability for systematics
  _bkg_m4l_ResDown.push_back(bkg_m4l_ResDown);// backgroun m4l probability for systematics
		       
} 


void HZZ4lNtupleFactory::FillHAdditionalInfo(const Float_t mZa, const Float_t mZb, Float_t mLL4, const Float_t mLL6, const Float_t SIP4, const Float_t iso34)
{
  _ZZmZa.push_back(mZa);
  _ZZmZb.push_back(mZb);
  _ZZmLL4.push_back(mLL4);
  _ZZmLL6.push_back(mLL6);
  _ZZSIP4.push_back(SIP4);
  //  _ZZiso34.push_back(iso34);

  return;
}



void HZZ4lNtupleFactory::FillZInfo(const Float_t ZMass, const Float_t ZPt, const Float_t Z1MassRefit)
{
  if(!_firstZStored){
    _Z1Mass.push_back(ZMass);
    _Z1Pt.push_back(ZPt);
    _Z1MassRefit.push_back(Z1MassRefit);
    _firstZStored = true;
  }
  else{
    _Z2Mass.push_back(ZMass);
    _Z2Pt.push_back(ZPt);
  }

  return;
}

void HZZ4lNtupleFactory::FillAngularInfo(const Float_t costhetastar, const Float_t phi, const Float_t costheta1, const Float_t costheta2, const Float_t phistar1, const Float_t phistar2,const Float_t xi, const Float_t xistar)
{
  _costhetastar.push_back(costhetastar);
  _phi.push_back(phi);
//   _helphiZ1.push_back(helphiZ1);
//   _helphiZ2.push_back(helphiZ2);
  _costheta1.push_back(costheta1);
  _costheta2.push_back(costheta2);
  _phistar1.push_back(phistar1);
  _phistar2.push_back(phistar2);
  _xi.push_back(xi);
  _xistar.push_back(xistar);

  return;
}

void HZZ4lNtupleFactory::FillLepInfo(const Float_t LepPt, const Float_t LepEta, const Float_t LepPhi, const Int_t LepId, const Float_t LepSIP, bool isID, float BDT, short parentId, int missingHit)
{
  switch(_LeptonIndex){

  case 1:
    _Lep1Pt.push_back(LepPt);
    _Lep1Eta.push_back(LepEta);
    _Lep1Phi.push_back(LepPhi);
    _Lep1LepId.push_back(LepId);
    _Lep1SIP.push_back(LepSIP);
    _Lep1isID.push_back(isID);
    _Lep1BDT.push_back(BDT);
	_Lep1missingHit.push_back((char)missingHit);
    _Lep1ParentId.push_back(parentId);
    break;

  case 2:
    _Lep2Pt.push_back(LepPt);
    _Lep2Eta.push_back(LepEta);
    _Lep2Phi.push_back(LepPhi);
    _Lep2LepId.push_back(LepId);
    _Lep2SIP.push_back(LepSIP);
    _Lep2isID.push_back(isID);
    _Lep2BDT.push_back(BDT);
	_Lep2missingHit.push_back((char)missingHit);
    _Lep2ParentId.push_back(parentId);
    break;
  
  case 3:
    _Lep3Pt.push_back(LepPt);
    _Lep3Eta.push_back(LepEta);
    _Lep3Phi.push_back(LepPhi);
    _Lep3LepId.push_back(LepId);
    _Lep3SIP.push_back(LepSIP);
    _Lep3isID.push_back(isID);
    _Lep3BDT.push_back(BDT);
	_Lep3missingHit.push_back((char)missingHit);
    _Lep3ParentId.push_back(parentId);
    break;

  case 4:
    _Lep4Pt.push_back(LepPt);
    _Lep4Eta.push_back(LepEta);
    _Lep4Phi.push_back(LepPhi);
    _Lep4LepId.push_back(LepId);
    _Lep4SIP.push_back(LepSIP);
    _Lep4isID.push_back(isID);
    _Lep4BDT.push_back(BDT);
	_Lep4missingHit.push_back((char)missingHit);
    _Lep4ParentId.push_back(parentId);
    break;

  default:
    std::cout << "Error in indexing the muons! Will abort..." << std::endl;
    assert(0);
  }

  _LeptonIndex++;

  return;
}

void HZZ4lNtupleFactory::FillLepIsolInfo(const Float_t LepchargedHadIso, const Float_t LepneutralHadIso, const Float_t LepphotonIso, const Float_t LepcombRelIsoPF)
{

  switch(_LeptonIsoIndex){

  case 1:
    _Lep1chargedHadIso.push_back(LepchargedHadIso);
    _Lep1neutralHadIso.push_back(LepneutralHadIso);
    _Lep1photonIso.push_back(LepphotonIso);
    _Lep1combRelIsoPF.push_back(LepcombRelIsoPF);
    break;

  case 2:
    _Lep2chargedHadIso.push_back(LepchargedHadIso);
    _Lep2neutralHadIso.push_back(LepneutralHadIso);
    _Lep2photonIso.push_back(LepphotonIso);
    _Lep2combRelIsoPF.push_back(LepcombRelIsoPF);
    break;

  case 3:
    _Lep3chargedHadIso.push_back(LepchargedHadIso);
    _Lep3neutralHadIso.push_back(LepneutralHadIso);
    _Lep3photonIso.push_back(LepphotonIso);
    _Lep3combRelIsoPF.push_back(LepcombRelIsoPF);
    break;

  case 4:
    _Lep4chargedHadIso.push_back(LepchargedHadIso);
    _Lep4neutralHadIso.push_back(LepneutralHadIso);
    _Lep4photonIso.push_back(LepphotonIso);
    _Lep4combRelIsoPF.push_back(LepcombRelIsoPF);
    break;


  default:
    std::cout << "Error in indexing the muon isolation variables! Will abort..." << std::endl;
    assert(0);
  }

  _LeptonIsoIndex++;

  return;
}

void HZZ4lNtupleFactory::FillPhotonInfo(const Float_t PhotPt, const Float_t PhotEta, const Float_t PhotPhi)
{

  _PhotPt.push_back(PhotPt);
  _PhotEta.push_back(PhotEta);
  _PhotPhi.push_back(PhotPhi);

  return;
}

void HZZ4lNtupleFactory::FillJetInfo(const Float_t JetPt, const Float_t JetEta, const Float_t JetPhi, const Float_t JetMass, const Float_t JetBTag, const Float_t JetSigma )
{

  _JetPt.push_back(JetPt);
  _JetEta.push_back(JetEta);
  _JetPhi.push_back(JetPhi);
  _JetMass.push_back(JetMass);
  _JetBTag.push_back(JetBTag);
  _JetSigma.push_back(JetSigma);

  return;
}

void HZZ4lNtupleFactory::FillDiJetInfo(const Float_t DiJetMass, const Float_t DiJetMassPlus, const Float_t DiJetMassMinus, const Float_t DiJetDEta)
{
  _DiJetMass = DiJetMass;
  _DiJetMassPlus = DiJetMassPlus;
  _DiJetMassMinus = DiJetMassMinus;
  _DiJetDEta = DiJetDEta;

  return;
}

void HZZ4lNtupleFactory::FillCategorizationInfo(const Int_t nExtraLep, const Int_t nExtraZ, const Int_t nJets, const Int_t nCleanedJets, const Int_t nCleanedJetsPt30)
{
  _nExtraLep.push_back(nExtraLep);
  _nExtraZ.push_back(nExtraZ);
  _nJets.push_back(nJets);
  _nCleanedJets.push_back(nCleanedJets);
  _nCleanedJetsPt30.push_back(nCleanedJetsPt30);

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
    _ExtraLep1Pt        .push_back(Pt);
    _ExtraLep1Eta       .push_back(Eta);
    _ExtraLep1Phi       .push_back(Phi);
    _ExtraLep1LepId     .push_back(LepId);
//     _ExtraLep1SIP       .push_back(SIP);
//     _ExtraLep1isID      .push_back(isID);
//     _ExtraLep1BDT       .push_back(BDT);
//     _ExtraLep1missingHit.push_back(missingHit);
//     _ExtraLep1chargedHadIso.push_back(chargedHadIso);
//     _ExtraLep1neutralHadIso.push_back(neutralHadIso);
//     _ExtraLep1photonIso    .push_back(photonIso);
//     _ExtraLep1combRelIsoPF .push_back(combRelIsoPF);
    break;

  case 2:
    _ExtraLep2Pt        .push_back(Pt);
    _ExtraLep2Eta       .push_back(Eta);
    _ExtraLep2Phi       .push_back(Phi);
    _ExtraLep2LepId     .push_back(LepId);
//     _ExtraLep2SIP       .push_back(SIP);
//     _ExtraLep2isID      .push_back(isID);
//     _ExtraLep2BDT       .push_back(BDT);
//     _ExtraLep2missingHit.push_back(missingHit);
//     _ExtraLep2chargedHadIso.push_back(chargedHadIso);
//     _ExtraLep2neutralHadIso.push_back(neutralHadIso);
//     _ExtraLep2photonIso    .push_back(photonIso);
//     _ExtraLep2combRelIsoPF .push_back(combRelIsoPF);
    break;

  case 3:
    _ExtraLep3Pt        .push_back(Pt);
    _ExtraLep3Eta       .push_back(Eta);
    _ExtraLep3Phi       .push_back(Phi);
    _ExtraLep3LepId     .push_back(LepId);
//     _ExtraLep3SIP       .push_back(SIP);
//     _ExtraLep3isID      .push_back(isID);
//     _ExtraLep3BDT       .push_back(BDT);
//     _ExtraLep3missingHit.push_back(missingHit);
//     _ExtraLep3chargedHadIso.push_back(chargedHadIso);
//     _ExtraLep3neutralHadIso.push_back(neutralHadIso);
//     _ExtraLep3photonIso    .push_back(photonIso);
//     _ExtraLep3combRelIsoPF .push_back(combRelIsoPF);
    break;

  default:
    std::cout << "Error in indexing the extra leptons ! Will abort..." << std::endl;
    assert(0);
  }

  return;
}
