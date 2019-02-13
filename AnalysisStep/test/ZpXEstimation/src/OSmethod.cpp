// Include classes
#include <ZZAnalysis/AnalysisStep/test/ZpXEstimation/include/OSmethod.h>

// Constructor
//============================================================
OSmethod::OSmethod():Tree()
{
   _current_process = -999;
   _current_final_state = -999;
   _current_category = -999;
   _current_category_stxs = -999;
   
   _s_process.push_back("Data");
   _s_process.push_back("WZ");
   _s_process.push_back("qqZZ");
   _s_process.push_back("DY");
   _s_process.push_back("ttbar");
   
   _s_flavour.push_back("ele");
   _s_flavour.push_back("mu");
   
   _s_final_state.push_back("4mu");
   _s_final_state.push_back("4e");
   _s_final_state.push_back("2e2mu");
   _s_final_state.push_back("2mu2e");
   _s_final_state.push_back("4l");
   
   _s_category.push_back("UnTagged");
   _s_category.push_back("VBF1jTagged");
   _s_category.push_back("VBF2jTagged");
   _s_category.push_back("VHLeptTagged");
   _s_category.push_back("VHHadrTagged");
   _s_category.push_back("ttHLeptTagged");
   _s_category.push_back("ttHHadrTagged");
   _s_category.push_back("VHMETTagged");
   _s_category.push_back("Inclusive");
   
   _s_category_stxs.push_back("ggH_0J_PTH_0_10");
   _s_category_stxs.push_back("ggH_0J_PTH_10_200");
   _s_category_stxs.push_back("ggH_1J_PTH_0_60");
   _s_category_stxs.push_back("ggH_1J_PTH_60_120");
   _s_category_stxs.push_back("ggH_1J_PTH_120_200");
   _s_category_stxs.push_back("ggH_2J_PTH_0_60");
   _s_category_stxs.push_back("ggH_2J_PTH_60_120");
   _s_category_stxs.push_back("ggH_2J_PTH_120_200");
   _s_category_stxs.push_back("ggH_PTH_200");
   _s_category_stxs.push_back("ggH_VBF");
   _s_category_stxs.push_back("VBF_1j");
   _s_category_stxs.push_back("VBF_2j");
   _s_category_stxs.push_back("VBF_2j_mjj_350_700_2j");
   _s_category_stxs.push_back("VBF_2j_mjj_GT700_2j");
   _s_category_stxs.push_back("VBF_2j_mjj_GT350_3j");
   _s_category_stxs.push_back("VBF_GT200_2J");
   _s_category_stxs.push_back("VH_Had");
   _s_category_stxs.push_back("VBF_rest_VH");
   _s_category_stxs.push_back("VH_lep_0_150");
   _s_category_stxs.push_back("VH_Lep_GT150");
   _s_category_stxs.push_back("ttH_Lep");
   _s_category_stxs.push_back("ttH_Had");
   _s_category_stxs.push_back("Inclusive");
   
   _s_region.push_back("2P2F");
   _s_region.push_back("3P1F");
   _s_region.push_back("OS");
	
   _s_variation.push_back("nominal");
   _s_variation.push_back("Up");
   _s_variation.push_back("Dn");
   
   DeclareFRHistos();
   DeclareDataMCHistos();
   DeclareZXHistos();
}
//============================================================



// Destructor
//====================
OSmethod::~OSmethod()
{
}
//====================


//===============================================================================
void OSmethod::FillFRHistos( TString input_file_data_name )
{
   input_file_data = TFile::Open( input_file_data_name);
   
   hCounters = (TH1F*)input_file_data->Get("CRZLTree/Counters");
   gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
   
   input_tree_data = (TTree*)input_file_data->Get("CRZLTree/candTree");
   Init( input_tree_data, input_file_data_name , false);
   
   _current_process = find_current_process(input_file_data_name);
   
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
	
	// Define some counters for control print out
	Int_t _total_events = 0;
	Int_t _failZ1MassCut = 0;
	Int_t _failLepPtCut = 0;
	Int_t _failSIPCut = 0;
	Int_t _failMETCut = 0;
	Int_t _passingSelection = 0;
	Int_t _faillingSelection = 0;
	Int_t _faillingJPsiMassCut = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
			   
      _total_events++;
		
		TLorentzVector p1,p2,p3;
		p1.SetPtEtaPhiM(LepPt->at(0), LepEta->at(0), LepPhi->at(0), 0.);
		p2.SetPtEtaPhiM(LepPt->at(1), LepEta->at(1), LepPhi->at(1), 0.);
		p3.SetPtEtaPhiM(LepPt->at(2), LepEta->at(2), LepPhi->at(2), 0.);
	   
      if ( abs(Z1Mass - 91.2) > 7. ) {_failZ1MassCut++; continue;}
      if ( (LepPt->at(0) > LepPt->at(1)) && (LepPt->at(0) < 20. || LepPt->at(1) < 10.) ) {_failLepPtCut++; continue;}
      if ( (LepPt->at(1) > LepPt->at(0)) && (LepPt->at(1) < 20. || LepPt->at(0) < 10.) ) {_failLepPtCut++; continue;}
      if ( LepSIP->at(2) > 4.) {_failSIPCut++; continue;}
      if ( PFMET > 25. ) {_failMETCut++; continue;}
		if ( (LepLepId->at(2) < 0 && LepLepId->at(0) > 0 && (p1+p3).M() < 4.) || (LepLepId->at(2) < 0 && LepLepId->at(1) > 0 && (p2+p3).M() < 4.) ) {_faillingJPsiMassCut++; continue;}
		if ( (LepLepId->at(2) > 0 && LepLepId->at(0) < 0 && (p1+p3).M() < 4.) || (LepLepId->at(2) > 0 && LepLepId->at(1) < 0 && (p2+p3).M() < 4.) ) {_faillingJPsiMassCut++; continue;}
      else
      {
         // Final event weight
         _k_factor = calculate_K_factor(input_file_data_name);
         _event_weight = (_lumi * 1000 * xsec * _k_factor * overallEventWeight) / gen_sum_weights;

         if(LepisID->at(2) && ((fabs(LepLepId->at(2)) == 11) ? LepCombRelIsoPF->at(2) < 999999. : LepCombRelIsoPF->at(2) < 0.35))
         {
		    _passingSelection++;
            if(fabs(LepLepId->at(2)) == 11 ) passing[_current_process][Settings::ele]->Fill(LepPt->at(2), (abs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
            else if(fabs(LepLepId->at(2)) == 13 ) passing[_current_process][Settings::mu]->Fill(LepPt->at(2), (abs(LepEta->at(2)) < 1.2) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
         }
         else
         {
		    _faillingSelection++;
            if(fabs(LepLepId->at(2)) == 11 ) failing[_current_process][Settings::ele]->Fill(LepPt->at(2), (abs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
            else if(fabs(LepLepId->at(2)) == 13 ) failing[_current_process][Settings::mu]->Fill(LepPt->at(2), (abs(LepEta->at(2)) < 1.2) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
         }
      }
   } // END events loop
	
	// Print Z + X expected yields for inclusive category
	if( _current_process == Settings::Data)
	{
		cout << endl;
		cout << "========================================================================================" << endl;
		cout << "[INFO] Control printout for Z+L control region." << endl;
		cout << "========================================================================================" << endl;
		cout << "[INFO] Total number of events in Z+L control region = " << _total_events << endl;
		cout << "[INFO] Events lost after  abs(Z1 - Z) < 7 GeV cut  = " << _failZ1MassCut << endl;
		cout << "[INFO] Events lost after LepPt > 20,10 GeV cut  = " << _failLepPtCut << endl;
		cout << "[INFO] Events lost after SIP < 4 cut  = " << _failSIPCut << endl;
		cout << "[INFO] Events lost after MET < 25 cut  = " << _failMETCut << endl;
		cout << "[INFO] Events lost after m_ll > 4 cut  = " << _faillingJPsiMassCut << endl;
		cout << "[INFO] Total events left = " << _passingSelection + _faillingSelection << endl;
		cout << "[INFO] Passing selection = " << _passingSelection  << endl;
		cout << "[INFO] Failling selection = " << _faillingSelection << endl;
		cout << "========================================================================================" << endl;
		cout << endl;
	}
   
   cout << "[INFO] Processing of " << input_file_data_name << " done." << endl;
}
//===============================================================================



//===============================================================================
void OSmethod::FillDataMCPlots( TString input_file_data_name )
{
   input_file_data = TFile::Open( input_file_data_name);
   
   hCounters = (TH1F*)input_file_data->Get("CRZLLTree/Counters");
   gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
   
   input_tree_data = (TTree*)input_file_data->Get("CRZLLTree/candTree");
   Init( input_tree_data, input_file_data_name , true);
   
   _current_process = find_current_process(input_file_data_name);
   
   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      if (!(test_bit(CRflag, CRZLLos_2P2F)) && !(test_bit(CRflag, CRZLLos_3P1F))) continue;
      
      _current_final_state = FindFinalState();
      
      for ( int j = 0; j < nCleanedJetsPt30; j++)
      {
         jetPt[j] = JetPt->at(j);
         jetEta[j] = JetEta->at(j);
         jetPhi[j] = JetPhi->at(j);
         jetMass[j] = JetMass->at(j);
         jetQGL[j] = JetQGLikelihood->at(j);
         jetPgOverPq[j] = 1./JetQGLikelihood->at(j) - 1.;
      }
      
		_current_category = categoryMor18(  nExtraLep,
														 nExtraZ,
														 nCleanedJetsPt30,
														 nCleanedJetsPt30BTagged_bTagSF,
														 jetQGL,
														 p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
														 p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
														 p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
														 p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
														 pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
														 p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
														 p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
														 p_HadWH_mavjj_JECNominal,
														 p_HadWH_mavjj_true_JECNominal,
														 p_HadZH_mavjj_JECNominal,
														 p_HadZH_mavjj_true_JECNominal,
														 jetPhi,
														 ZZMass,
														 PFMET,
														 false,// Use VHMET category
														 false);// Use QG tagging
      
      _current_category_stxs = stage1_reco_1p1 ( nCleanedJetsPt30,
                                                 DiJetMass,
                                                 ZZPt,
                                                 _current_category,
                                                 ZZjjPt);
      
      _k_factor = calculate_K_factor(input_file_data_name);
      _event_weight = (_lumi * 1000 * xsec * _k_factor * overallEventWeight) / gen_sum_weights;
      
      if ( test_bit(CRflag, CRZLLos_2P2F) ) histos_1D[Settings::reg2P2F][_current_process][_current_final_state][_current_category_stxs]->Fill(ZZMass, (_current_process == Settings::Data) ? 1 :  _event_weight);
      if ( test_bit(CRflag, CRZLLos_3P1F) ) histos_1D[Settings::reg3P1F][_current_process][_current_final_state][_current_category_stxs]->Fill(ZZMass, (_current_process == Settings::Data) ? 1 :  _event_weight);
		if ( Z1Flav < 0 && Z2Flav < 0 )       histos_1D[Settings::regOS][_current_process][_current_final_state][_current_category_stxs]->Fill(ZZMass, (_current_process == Settings::Data) ? 1 :  _event_weight);
   
   } // END events loop
   
   cout << "[INFO] Processing of " << input_file_data_name << " done." << endl;
}
//===============================================================================



//===============================================================================
void OSmethod::MakeHistogramsZX( TString input_file_data_name, TString  input_file_FR_name )
{
   
   FakeRates *FR = new FakeRates( input_file_FR_name );
   
   input_file_data = TFile::Open( input_file_data_name);
   input_tree_data = (TTree*)input_file_data->Get("CRZLLTree/candTree");
   Init( input_tree_data, input_file_data_name , true);
   
   
   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      if (!(test_bit(CRflag, CRZLLos_2P2F)) && !(test_bit(CRflag, CRZLLos_3P1F))) continue;
		
		if ( ZZMass < 70. ) continue;
		
      _current_final_state = FindFinalState();
      
      for ( int j = 0; j < nCleanedJetsPt30; j++)
      {
         jetPt[j] = JetPt->at(j);
         jetEta[j] = JetEta->at(j);
         jetPhi[j] = JetPhi->at(j);
         jetMass[j] = JetMass->at(j);
         jetQGL[j] = JetQGLikelihood->at(j);
         jetPgOverPq[j] = 1./JetQGLikelihood->at(j) - 1.;
      }
      
		_current_category = categoryMor18(  nExtraLep,
														 nExtraZ,
														 nCleanedJetsPt30,
														 nCleanedJetsPt30BTagged_bTagSF,
														 jetQGL,
														 p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
														 p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
														 p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
														 p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
														 pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
														 p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
														 p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
														 p_HadWH_mavjj_JECNominal,
														 p_HadWH_mavjj_true_JECNominal,
														 p_HadZH_mavjj_JECNominal,
														 p_HadZH_mavjj_true_JECNominal,
														 jetPhi,
														 ZZMass,
														 PFMET,
														 false,// Use VHMET category
														 false);// Use QG tagging
      
      _current_category_stxs = stage1_reco_1p1 ( nCleanedJetsPt30,
                                                 DiJetMass,
                                                 ZZPt,
                                                 _current_category,
                                                 ZZjjPt);

      if ( test_bit(CRflag, CRZLLos_2P2F) )
      {
         _f3    = FR->GetFakeRate(LepPt->at(2),LepEta->at(2),LepLepId->at(2));
         _f3_Up = FR->GetFakeRate_Up(LepPt->at(2),LepEta->at(2),LepLepId->at(2));
         _f3_Dn = FR->GetFakeRate_Dn(LepPt->at(2),LepEta->at(2),LepLepId->at(2));
         _f4    = FR->GetFakeRate(LepPt->at(3),LepEta->at(3),LepLepId->at(3));
         _f4_Up = FR->GetFakeRate_Up(LepPt->at(3),LepEta->at(3),LepLepId->at(3));
         _f4_Dn = FR->GetFakeRate_Dn(LepPt->at(3),LepEta->at(3),LepLepId->at(3));
			
//         cout << "===============" << endl;
//         cout << "f3 = " << _f3 << endl;
//         cout << "f4 = " << _f4 << endl;
//         cout << "weight = " << (_f3/(1-_f3))*(_f4/(1-_f4)) << endl;
//         cout << "weight_up = " << (_f3_Up/(1-_f3_Up))*(_f4_Up/(1-_f4_Up)) << endl;
//         cout << "weight_dn = " << (_f3_Dn/(1-_f3_Dn))*(_f4_Dn/(1-_f4_Dn)) << endl;
			
         h_from2P2F_SR[Settings::nominal][_current_final_state][_current_category_stxs]->Fill(ZZMass, (_f3/(1-_f3))*(_f4/(1-_f4)) );
         h_from2P2F_3P1F[Settings::nominal][_current_final_state][_current_category_stxs]->Fill(ZZMass, (_f3/(1-_f3))+(_f4/(1-_f4)) );
			
         h_from2P2F_SR[Settings::Up][_current_final_state][_current_category_stxs]->Fill(ZZMass, (_f3_Up/(1-_f3_Up))*(_f4_Up/(1-_f4_Up)) );
         h_from2P2F_3P1F[Settings::Up][_current_final_state][_current_category_stxs]->Fill(ZZMass, (_f3_Up/(1-_f3_Up))+(_f4_Up/(1-_f4_Up)) );
			
         h_from2P2F_SR[Settings::Dn][_current_final_state][_current_category_stxs]->Fill(ZZMass, (_f3_Dn/(1-_f3_Dn))*(_f4_Dn/(1-_f4_Dn)) );
         h_from2P2F_3P1F[Settings::Dn][_current_final_state][_current_category_stxs]->Fill(ZZMass, (_f3_Dn/(1-_f3_Dn))+(_f4_Dn/(1-_f4_Dn)) );
      }
      if ( test_bit(CRflag, CRZLLos_3P1F) )
      {
         if(LepisID->at(3) && ((fabs(LepLepId->at(2)) == 11) ? LepCombRelIsoPF->at(3) < 999999. : LepCombRelIsoPF->at(3) < 0.35))
         {
         	_f4    = FR->GetFakeRate(LepPt->at(2),LepEta->at(2),LepLepId->at(2));
         	_f4_Up = FR->GetFakeRate_Up(LepPt->at(2),LepEta->at(2),LepLepId->at(2));
         	_f4_Dn = FR->GetFakeRate_Dn(LepPt->at(2),LepEta->at(2),LepLepId->at(2));
         }
         else
         {
				_f4    = FR->GetFakeRate(LepPt->at(3),LepEta->at(3),LepLepId->at(3));
				_f4_Up = FR->GetFakeRate_Up(LepPt->at(3),LepEta->at(3),LepLepId->at(3));
				_f4_Dn = FR->GetFakeRate_Dn(LepPt->at(3),LepEta->at(3),LepLepId->at(3));
         }
         
         h_from3P1F_SR[Settings::nominal][_current_final_state][_current_category_stxs]->Fill(ZZMass, (_f4/(1-_f4)) );
         h_from3P1F_SR[Settings::Up][_current_final_state][_current_category_stxs]->Fill(ZZMass, (_f4_Up/(1-_f4_Up)) );
         h_from3P1F_SR[Settings::Dn][_current_final_state][_current_category_stxs]->Fill(ZZMass, (_f4_Dn/(1-_f4_Dn)) );
      }
      
   }
   
   cout << "[INFO] Processing of " << input_file_data_name << " done." << endl;
}
//===============================================================================

//===============================================================================
void OSmethod::MakeZXMCContribution( TString input_file_data_name, TString  input_file_FR_name )
{
   
   FakeRates *FR = new FakeRates( input_file_FR_name );
   input_file_data = TFile::Open( input_file_data_name);
   
   hCounters = (TH1F*)input_file_data->Get("CRZLLTree/Counters");
   gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
   
   input_tree_data = (TTree*)input_file_data->Get("CRZLLTree/candTree");
   Init( input_tree_data, input_file_data_name , true);
   
   
   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      if (!(test_bit(CRflag, CRZLLos_3P1F))) continue;
      
      _current_final_state = FindFinalState();
      
      for ( int j = 0; j < nCleanedJetsPt30; j++)
      {
         jetPt[j] = JetPt->at(j);
         jetEta[j] = JetEta->at(j);
         jetPhi[j] = JetPhi->at(j);
         jetMass[j] = JetMass->at(j);
         jetQGL[j] = JetQGLikelihood->at(j);
         jetPgOverPq[j] = 1./JetQGLikelihood->at(j) - 1.;
      }
      
		_current_category = categoryMor18(  nExtraLep,
														 nExtraZ,
														 nCleanedJetsPt30,
														 nCleanedJetsPt30BTagged_bTagSF,
														 jetQGL,
														 p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
														 p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
														 p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
														 p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
														 pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
														 p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
														 p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
														 p_HadWH_mavjj_JECNominal,
														 p_HadWH_mavjj_true_JECNominal,
														 p_HadZH_mavjj_JECNominal,
														 p_HadZH_mavjj_true_JECNominal,
														 jetPhi,
														 ZZMass,
														 PFMET,
														 false,// Use VHMET category
														 false);// Use QG tagging
      
      _current_category_stxs = stage1_reco_1p1 ( nCleanedJetsPt30,
                                                 DiJetMass,
                                                 ZZPt,
                                                 _current_category,
                                                 ZZjjPt);
      
      _k_factor = calculate_K_factor(input_file_data_name);
      _event_weight = (_lumi * 1000 * xsec * _k_factor * overallEventWeight) / gen_sum_weights;
      
      if(LepisID->at(3) && ((fabs(LepLepId->at(3)) == 11) ? LepCombRelIsoPF->at(3) < 999999. : LepCombRelIsoPF->at(3) < 0.35))
      {
			_f4    = FR->GetFakeRate(LepPt->at(2),LepEta->at(2),LepLepId->at(2));
			_f4_Up = FR->GetFakeRate_Up(LepPt->at(2),LepEta->at(2),LepLepId->at(2));
			_f4_Dn = FR->GetFakeRate_Dn(LepPt->at(2),LepEta->at(2),LepLepId->at(2));
      }
      else
      {
			_f4    = FR->GetFakeRate(LepPt->at(3),LepEta->at(3),LepLepId->at(3));
			_f4_Up = FR->GetFakeRate_Up(LepPt->at(3),LepEta->at(3),LepLepId->at(3));
			_f4_Dn = FR->GetFakeRate_Dn(LepPt->at(3),LepEta->at(3),LepLepId->at(3));
		}

      h_from3P1F_SR_ZZonly[Settings::nominal][_current_final_state][_current_category_stxs]->Fill(ZZMass, _event_weight * (_f4/(1-_f4)) );
      h_from3P1F_SR_ZZonly[Settings::Up][_current_final_state][_current_category_stxs]->Fill(ZZMass, _event_weight * (_f4_Up/(1-_f4_Up)) );
      h_from3P1F_SR_ZZonly[Settings::Dn][_current_final_state][_current_category_stxs]->Fill(ZZMass, _event_weight * (_f4_Dn/(1-_f4_Dn)) );
      
   }
   
   cout << "[INFO] Processing of " << input_file_data_name << " done." << endl;
}
//===============================================================================




//===============================================================
void OSmethod::DeclareFRHistos()
{
   for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         _histo_name = "Passing_" + _s_process.at(i_proc) + "_" + _s_flavour.at(i_flav);
         passing[i_proc][i_flav] = new TH2F(_histo_name,"", 80, 0, 80, 2, 0, 2);
         
         _histo_name = "Failing_" + _s_process.at(i_proc) + "_" + _s_flavour.at(i_flav);
         failing[i_proc][i_flav] = new TH2F(_histo_name,"", 80, 0, 80, 2, 0, 2);

      }
      
      _histo_name = "Passing_Total_" + _s_flavour.at(i_flav);
      passing[Settings::Total][i_flav] = new TH2F(_histo_name,"", 80, 0, 80, 2, 0, 2);
      _histo_name = "Failing_Total_" + _s_flavour.at(i_flav);
      failing[Settings::Total][i_flav] = new TH2F(_histo_name,"", 80, 0, 80, 2, 0, 2);

   }

}
//===============================================================

//===============================================================
void OSmethod::DeclareDataMCHistos()
{
   for (int i_reg = 0; i_reg < num_of_regions_os; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
         {
            for (int i_cat = 0; i_cat < num_of_categories_stxs; i_cat++)
            {
               _histo_name = "M4l_" + _s_region.at(i_reg) + "_" + _s_process.at(i_proc) + "_" + _s_final_state.at(i_fs) + "_" + _s_category_stxs.at(i_cat);
               _histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
               histos_1D[i_reg][i_proc][i_fs][i_cat] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
            }
         }
      }
   }
   
}
//===============================================================

//===============================================================
void OSmethod::DeclareZXHistos()
{
   for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
   {
      for (int i_cat = 0; i_cat < num_of_categories_stxs; i_cat++)
      {
      	for(int i_var = 0; i_var < num_of_fr_variations; i_var++)
      	{
				_histo_name = "h_from2P2F_SR_" + _s_variation.at(i_var) + "_" + _s_final_state.at(i_fs) + "_" + _s_category_stxs.at(i_cat);
				_histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
				h_from2P2F_SR[i_var][i_fs][i_cat] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
				
				_histo_name = "h_from2P2F_3P1F_" + _s_variation.at(i_var) + "_" + _s_final_state.at(i_fs) + "_" + _s_category_stxs.at(i_cat);
				_histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
				h_from2P2F_3P1F[i_var][i_fs][i_cat] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
				
				_histo_name = "h_from3P1F_SR_final_" + _s_variation.at(i_var) + "_" + _s_final_state.at(i_fs) + "_" + _s_category_stxs.at(i_cat);
				_histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
				h_from3P1F_SR_final[i_var][i_fs][i_cat] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
				
				_histo_name = "h_from3P1F_SR_" + _s_variation.at(i_var) + "_" + _s_final_state.at(i_fs) + "_" + _s_category_stxs.at(i_cat);
				_histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
				h_from3P1F_SR[i_var][i_fs][i_cat] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
				
				_histo_name = "h_from3P1F_SR_ZZonly_" + _s_variation.at(i_var) + "_" + _s_final_state.at(i_fs) + "_" + _s_category_stxs.at(i_cat);
				_histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
				h_from3P1F_SR_ZZonly[i_var][i_fs][i_cat] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
				
				_histo_name = "ZX_" + _s_variation.at(i_var) + "_" + _s_final_state.at(i_fs) + "_" + _s_category_stxs.at(i_cat);
				_histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
				histos_ZX[i_var][i_fs][i_cat] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
				
			}
			
      }
   }
}
//===============================================================

//===============================================================
void OSmethod::SaveFRHistos( TString file_name,  bool subtractWZ, bool remove_negative_bins)
{
   TFile* fOutHistos = TFile::Open(file_name, "recreate");
   fOutHistos->cd();
   
   // Copy data histos to total histos, if there is no WZ subtraction this is the final histo for fake rate calculation
   for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
   {
      passing[Settings::Total][i_flav]->Add(passing[Settings::Data][i_flav], 1.);
      failing[Settings::Total][i_flav]->Add(failing[Settings::Data][i_flav], 1.);
   }
   
   if (subtractWZ ) SubtractWZ(); // Subtract WZ contribution from MC estimate
   
   if ( remove_negative_bins ) // Set negative bins to zero
   {
      for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
      {
         RemoveNegativeBins2D( passing[Settings::Total][i_flav] );
         RemoveNegativeBins2D( failing[Settings::Total][i_flav] );
      }
      cout << "[INFO] Negative bins removed." << endl;
   }

   for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         passing[i_proc][i_flav]->Write();
         failing[i_proc][i_flav]->Write();
      }
      
   passing[Settings::Total][i_flav]->Write();
   failing[Settings::Total][i_flav]->Write();
      
   }
   
   fOutHistos->Close();
   delete fOutHistos;
   
   cout << "[INFO] All FakeRate histograms saved." << endl;
}
//===============================================================

//===============================================================
void OSmethod::SaveDataMCHistos( TString file_name )
{
   FillDataMCInclusive();
   
   TFile* fOutHistos = TFile::Open(file_name, "recreate");
   fOutHistos->cd();
   
   for (int i_reg = 0; i_reg < num_of_regions_os; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
         {
            for (int i_cat = 0; i_cat < num_of_categories_stxs; i_cat++)
            {
               histos_1D[i_reg][i_proc][i_fs][i_cat]->Write();
            }
         }
      }
   }
   
   fOutHistos->Close();
   delete fOutHistos;
   
   cout << "[INFO] All Data/MC histograms saved." << endl;
}
//===============================================================

//===============================================================
void OSmethod::FillDataMCInclusive( )
{
   for (int i_reg = 0; i_reg < num_of_regions_os; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < Settings::fs4l; i_fs++)
         {
            for (int i_cat = 0; i_cat < Settings::inclusive_stxs; i_cat++)
            {
               histos_1D[i_reg][i_proc][i_fs][Settings::inclusive_stxs]->Add(histos_1D[i_reg][i_proc][i_fs][i_cat]);
               histos_1D[i_reg][i_proc][Settings::fs4l][i_cat]    ->Add(histos_1D[i_reg][i_proc][i_fs][i_cat]);
            }
         }
      }
   }
   
   for (int i_reg = 0; i_reg < num_of_regions_os; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < Settings::fs4l; i_fs++)
         {
            histos_1D[i_reg][i_proc][Settings::fs4l][Settings::inclusive_stxs]->Add(histos_1D[i_reg][i_proc][i_fs][Settings::inclusive_stxs]);
         }
      }
   }
   
   cout << "[INFO] All Data/MC histograms summed." << endl;
}
//===============================================================

//===============================================================
void OSmethod::SaveZXHistos( TString file_name , bool remove_negative_bins)
{
   FillZXInclusive(remove_negative_bins);
   
   TFile* fOutHistos = TFile::Open(file_name, "recreate");
   fOutHistos->cd();
   
   for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
   {
      for (int i_cat = 0; i_cat < num_of_categories_stxs; i_cat++)
      {
			for(int i_var = 0; i_var < num_of_fr_variations; i_var++)
      	{
				h_from2P2F_SR[i_var][i_fs][i_cat]->Write();
         	h_from2P2F_3P1F[i_var][i_fs][i_cat]->Write();
         	h_from3P1F_SR_final[i_var][i_fs][i_cat]->Write();
         	h_from3P1F_SR[i_var][i_fs][i_cat]->Write();
				h_from3P1F_SR_ZZonly[i_var][i_fs][i_cat]->Write();
         	histos_ZX[i_var][i_fs][i_cat]->Write();
      	}

      }
   }
   
   fOutHistos->Close();
   delete fOutHistos;
   
   cout << "[INFO] All Z+X histograms saved." << endl;
}
//===============================================================

//===============================================================
void OSmethod::FillZXInclusive( bool remove_negative_bins )
{
	if ( remove_negative_bins )
	{
		for (int i_fs = 0; i_fs <= Settings::fs4l; i_fs++)
		{
			for (int i_cat = 0; i_cat <= Settings::inclusive_stxs; i_cat++)
			{
				for(int i_var = 0; i_var < num_of_fr_variations; i_var++)
      		{
					RemoveNegativeBins1D(h_from2P2F_SR[i_var][i_fs][i_cat]);
					RemoveNegativeBins1D(h_from2P2F_3P1F[i_var][i_fs][i_cat]);
					RemoveNegativeBins1D(h_from3P1F_SR_final[i_var][i_fs][i_cat]);
					RemoveNegativeBins1D(h_from3P1F_SR[i_var][i_fs][i_cat]);
					RemoveNegativeBins1D(h_from3P1F_SR_ZZonly[i_var][i_fs][i_cat]);
					
				}
			}
		}
	}
	
   for (int i_fs = 0; i_fs < Settings::fs4l; i_fs++)
   {
      for (int i_cat = 0; i_cat < Settings::inclusive_stxs; i_cat++)
      {
			for(int i_var = 0; i_var < num_of_fr_variations; i_var++)
      	{
				h_from2P2F_SR[i_var][i_fs][Settings::inclusive_stxs]->Add(h_from2P2F_SR[i_var][i_fs][i_cat]);
				h_from2P2F_SR[i_var][Settings::fs4l][i_cat]    ->Add(h_from2P2F_SR[i_var][i_fs][i_cat]);
				
				h_from2P2F_3P1F[i_var][i_fs][Settings::inclusive_stxs]->Add(h_from2P2F_3P1F[i_var][i_fs][i_cat]);
				h_from2P2F_3P1F[i_var][Settings::fs4l][i_cat]    ->Add(h_from2P2F_3P1F[i_var][i_fs][i_cat]);
				
				h_from3P1F_SR_final[i_var][i_fs][Settings::inclusive_stxs]->Add(h_from3P1F_SR_final[i_var][i_fs][i_cat]);
				h_from3P1F_SR_final[i_var][Settings::fs4l][i_cat]    ->Add(h_from3P1F_SR_final[i_var][i_fs][i_cat]);
				
				h_from3P1F_SR[i_var][i_fs][Settings::inclusive_stxs]->Add(h_from3P1F_SR[i_var][i_fs][i_cat]);
				h_from3P1F_SR[i_var][Settings::fs4l][i_cat]    ->Add(h_from3P1F_SR[i_var][i_fs][i_cat]);
				
				h_from3P1F_SR_ZZonly[i_var][i_fs][Settings::inclusive_stxs]->Add(h_from3P1F_SR_ZZonly[i_var][i_fs][i_cat]);
				h_from3P1F_SR_ZZonly[i_var][Settings::fs4l][i_cat]    ->Add(h_from3P1F_SR_ZZonly[i_var][i_fs][i_cat]);
         }
      }
   }
	
   for (int i_fs = 0; i_fs <= Settings::fs4l; i_fs++)
   {
      for (int i_cat = 0; i_cat <= Settings::inclusive_stxs; i_cat++)
      {
			for(int i_var = 0; i_var < num_of_fr_variations; i_var++)
      	{
				h_from3P1F_SR_final[i_var][i_fs][i_cat]->Add(h_from3P1F_SR[i_var][i_fs][i_cat], 1.);
				h_from3P1F_SR_final[i_var][i_fs][i_cat]->Add(h_from3P1F_SR_ZZonly[i_var][i_fs][i_cat], -1.);
				h_from3P1F_SR_final[i_var][i_fs][i_cat]->Add(h_from2P2F_SR[i_var][i_fs][i_cat], -2.);
         }
      }
   }
	
	if ( remove_negative_bins )
	{
		for (int i_fs = 0; i_fs <= Settings::fs4l; i_fs++)
		{
			for (int i_cat = 0; i_cat <= Settings::inclusive_stxs; i_cat++)
			{
				for(int i_var = 0; i_var < num_of_fr_variations; i_var++)
      		{
					RemoveNegativeBins1D(h_from3P1F_SR_final[i_var][i_fs][i_cat]);
				}
			}
		}
	}
	
	for (int i_fs = 0; i_fs < Settings::fs4l; i_fs++)
	{
		for(int i_var = 0; i_var < num_of_fr_variations; i_var++)
		{
			h_from2P2F_SR[i_var][Settings::fs4l][Settings::inclusive_stxs]->Add(h_from2P2F_SR[i_var][i_fs][Settings::inclusive_stxs]);
			h_from2P2F_3P1F[i_var][Settings::fs4l][Settings::inclusive_stxs]->Add(h_from2P2F_3P1F[i_var][i_fs][Settings::inclusive_stxs]);
			h_from3P1F_SR_final[i_var][Settings::fs4l][Settings::inclusive_stxs]->Add(h_from3P1F_SR_final[i_var][i_fs][Settings::inclusive_stxs]);
			h_from3P1F_SR[i_var][Settings::fs4l][Settings::inclusive_stxs]->Add(h_from3P1F_SR[i_var][i_fs][Settings::inclusive_stxs]);
			h_from3P1F_SR_ZZonly[i_var][Settings::fs4l][Settings::inclusive_stxs]->Add(h_from3P1F_SR_ZZonly[i_var][i_fs][Settings::inclusive_stxs]);
		}
		
	}
	
   for (int i_fs = 0; i_fs <= Settings::fs4l; i_fs++)
   {
      for (int i_cat = 0; i_cat <= Settings::inclusive_stxs; i_cat++)
      {
			for(int i_var = 0; i_var < num_of_fr_variations; i_var++)
      	{
				histos_ZX[i_var][i_fs][i_cat]->Add(h_from3P1F_SR_final[i_var][i_fs][i_cat], 1.);
				histos_ZX[i_var][i_fs][i_cat]->Add(h_from2P2F_SR[i_var][i_fs][i_cat], 1.);
         }
      }
   }

   
   cout << "[INFO] All Z+X histograms summed." << endl;
}
//===============================================================

//===============================================================
void OSmethod::GetFRHistos( TString file_name)
{
   TFile* histo_file = TFile::Open(file_name);
   
   for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         _histo_name = "Passing_" + _s_process.at(i_proc) + "_" + _s_flavour.at(i_flav);
         passing[i_proc][i_flav] = (TH2F*)histo_file->Get(_histo_name);
         
         _histo_name = "Failing_" + _s_process.at(i_proc) + "_" + _s_flavour.at(i_flav);
         failing[i_proc][i_flav] = (TH2F*)histo_file->Get(_histo_name);
         
      }
      
      _histo_name = "Passing_Total_" + _s_flavour.at(i_flav);
      passing[Settings::Total][i_flav] = (TH2F*)histo_file->Get(_histo_name);
      _histo_name = "Failing_Total_" + _s_flavour.at(i_flav);
      failing[Settings::Total][i_flav] = (TH2F*)histo_file->Get(_histo_name);
      
   }
   
   cout << "[INFO] All FakeRate histograms retrieved from file." << endl;
}
//===============================================================

//===============================================================
void OSmethod::GetDataMCHistos( TString file_name)
{
   TFile* histo_file = TFile::Open(file_name);
   
   for (int i_reg = 0; i_reg < num_of_regions_os; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
         {
            for (int i_cat = 0; i_cat < num_of_categories_stxs; i_cat++)
            {
               _histo_name = "M4l_" + _s_region.at(i_reg) + "_" + _s_process.at(i_proc) + "_" + _s_final_state.at(i_fs) + "_" + _s_category_stxs.at(i_cat);
               histos_1D[i_reg][i_proc][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name);
            }
         }
      }
   }
   
   cout << "[INFO] All Data/MC histograms retrieved from file." << endl;
}

//===============================================================

//===============================================================
void OSmethod::GetZXHistos( TString file_name)
{
   TFile* histo_file = TFile::Open(file_name);
   
   for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
   {
      for (int i_cat = 0; i_cat < num_of_categories_stxs; i_cat++)
      {
			for(int i_var = 0; i_var < num_of_fr_variations; i_var++)
      	{
				_histo_name = "h_from2P2F_SR_" + _s_variation.at(i_var )+ "_" + _s_final_state.at(i_fs) + "_" + _s_category_stxs.at(i_cat);
				h_from2P2F_SR[i_var][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name);
				
				_histo_name = "h_from2P2F_3P1F_" + _s_variation.at(i_var )+ "_" + _s_final_state.at(i_fs) + "_" + _s_category_stxs.at(i_cat);
				h_from2P2F_3P1F[i_var][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name);
				
				_histo_name = "h_from3P1F_SR_final_" + _s_variation.at(i_var )+ "_" + _s_final_state.at(i_fs) + "_" + _s_category_stxs.at(i_cat);
				h_from3P1F_SR_final[i_var][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name);
				
				_histo_name = "h_from3P1F_SR_" + _s_variation.at(i_var )+ "_" + _s_final_state.at(i_fs) + "_" + _s_category_stxs.at(i_cat);
				h_from3P1F_SR[i_var][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name);
				
				_histo_name = "h_from3P1F_SR_ZZonly_" + _s_variation.at(i_var )+ "_" + _s_final_state.at(i_fs) + "_" + _s_category_stxs.at(i_cat);
				h_from3P1F_SR_ZZonly[i_var][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name);
				
				_histo_name = "ZX_" + _s_variation.at(i_var )+ "_" + _s_final_state.at(i_fs) + "_" + _s_category_stxs.at(i_cat);
				histos_ZX[i_var][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name);
         }
      }
   }
   
   cout << "[INFO] All Z+X histograms retrieved from file." << endl;
}

//===============================================================



void OSmethod::PrintZXYields()
{
	double stat;
	double syst;
	double comb;
	double yield, yield_up;
	
	cout << endl;
	cout << "==============================================================" << endl;
	cout << "[INFO] Control printout for OS Z+X yields in final states "<< endl;
	cout << "==============================================================" << endl;
	for( int i_fs = 0; i_fs < Settings::fs4l ; i_fs++ )
	{
		for ( int i_cat = 0; i_cat <= Settings::inclusive_stxs; i_cat++)
		{

		   yield = histos_ZX[Settings::nominal][i_fs][i_cat]->IntegralAndError(0,histos_ZX[Settings::nominal][i_fs][i_cat]->GetSize() - 2,stat);
		   yield_up = histos_ZX[Settings::Up][i_fs][i_cat]->Integral();
		   syst = ((yield_up/yield) - 1.) * yield;
		   comb = sqrt(stat*stat + syst*syst);
			
		cout << "Category: " << _s_category_stxs.at(i_cat) << "   Final state: " << _s_final_state.at(i_fs) << endl;
		cout << yield << " +/- " << comb << " (total.)   - " << stat << " (stat.)   - " << syst << " (syst.)" << endl;
		}
		
		cout << "============================================================" << endl;
	}
}



//========================================================================================================
void OSmethod::PlotDataMC_2P2F( TString variable_name, TString folder )
{
   TCanvas *c;
   c = new TCanvas("2P2F", variable_name, 600, 600);
   
   if ( GetVarLogX( variable_name) ) c->SetLogx();
   if ( GetVarLogY( variable_name) ) c->SetLogy();
   
   for( int i_fs = 0; i_fs < Settings::fs4l ; i_fs++ )
   {
      histos_1D[Settings::reg2P2F][Settings::WZ][i_fs][Settings::inclusive_stxs]   ->SetFillColor(kMagenta);
      histos_1D[Settings::reg2P2F][Settings::qqZZ][i_fs][Settings::inclusive_stxs] ->SetFillColor(kCyan+1);
      histos_1D[Settings::reg2P2F][Settings::DY][i_fs][Settings::inclusive_stxs]   ->SetFillColor(kGreen-1);
      histos_1D[Settings::reg2P2F][Settings::ttbar][i_fs][Settings::inclusive_stxs]->SetFillColor(kBlue-2);
      
      histos_1D[Settings::reg2P2F][Settings::WZ][i_fs][Settings::inclusive_stxs]   ->SetLineColor(kMagenta);
      histos_1D[Settings::reg2P2F][Settings::qqZZ][i_fs][Settings::inclusive_stxs] ->SetLineColor(kCyan+1);
      histos_1D[Settings::reg2P2F][Settings::DY][i_fs][Settings::inclusive_stxs]   ->SetLineColor(kGreen-1);
      histos_1D[Settings::reg2P2F][Settings::ttbar][i_fs][Settings::inclusive_stxs]->SetLineColor(kBlue-2);
      
      histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive_stxs]->SetMarkerSize(0.8);
      histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive_stxs]->SetMarkerStyle(20);
      histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive_stxs]->SetBinErrorOption(TH1::kPoisson);
      histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive_stxs]->SetLineColor(kBlack);
      
      THStack *stack = new THStack( "stack", "stack" );
		stack->Add(histos_1D[Settings::reg2P2F][Settings::qqZZ][i_fs][Settings::inclusive_stxs]);
      stack->Add(histos_1D[Settings::reg2P2F][Settings::WZ][i_fs][Settings::inclusive_stxs]);
		stack->Add(histos_1D[Settings::reg2P2F][Settings::ttbar][i_fs][Settings::inclusive_stxs]);
      stack->Add(histos_1D[Settings::reg2P2F][Settings::DY][i_fs][Settings::inclusive_stxs]);
   
      stack->Draw("HIST");
      
      float data_max = histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive_stxs]->GetBinContent(histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive_stxs]->GetMaximumBin());
      float data_max_error = histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive_stxs]->GetBinErrorUp(histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive_stxs]->GetMaximumBin());
      
      stack->SetMinimum(1e-5);
      stack->SetMaximum((data_max + data_max_error)*1.1);
      
      TString _fs_label;
      if ( i_fs == Settings::fs4e) _fs_label = "m_{4#font[12]{e}} (GeV)";
      if ( i_fs == Settings::fs4mu) _fs_label = "m_{4#font[12]{#mu}} (GeV)";
      if ( i_fs == Settings::fs2e2mu) _fs_label = "m_{2#font[12]{e}2#font[12]{#mu}} (GeV)";
      if ( i_fs == Settings::fs2mu2e) _fs_label = "m_{2#font[12]{#mu}2#font[12]{e}} (GeV)";
      stack->GetXaxis()->SetTitle(_fs_label);
      stack->GetXaxis()->SetTitleSize(0.04);
      stack->GetXaxis()->SetLabelSize(0.04);
      stack->GetYaxis()->SetTitle(histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive_stxs]->GetYaxis()->GetTitle());
      stack->GetYaxis()->SetTitleSize(0.04);
      stack->GetYaxis()->SetLabelSize(0.04);
      
      stack->GetXaxis()->SetTitleOffset(1.2);
      stack->GetYaxis()->SetTitleOffset(1.25);
      
      histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive_stxs]->Draw("SAME p E1 X0");
      
      TLegend *legend;
      legend  = CreateLegend_2P2F("right",histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive_stxs],histos_1D[Settings::reg2P2F][Settings::WZ][i_fs][Settings::inclusive_stxs],histos_1D[Settings::reg2P2F][Settings::qqZZ][i_fs][Settings::inclusive_stxs],histos_1D[Settings::reg2P2F][Settings::DY][i_fs][Settings::inclusive_stxs],histos_1D[Settings::reg2P2F][Settings::ttbar][i_fs][Settings::inclusive_stxs]);
      legend->Draw();

      // Draw lumi
      CMS_lumi *lumi = new CMS_lumi;
      lumi->set_lumi(c, _lumi, 0);
      
      TString _out_file_name;
      _out_file_name = folder + "/" + variable_name + "_OS_2P2F_" + _s_final_state.at(i_fs) + "_" + _s_category_stxs.at(Settings::inclusive_stxs);
      SavePlots(c, _out_file_name);

   }
}
//========================================================================================================


//========================================================================================================
void OSmethod::PlotDataMC_3P1F( TString variable_name, TString folder )
{
   TCanvas *c;
   c = new TCanvas("3P2F", variable_name, 600, 600);
   
   if ( GetVarLogX( variable_name) ) c->SetLogx();
   if ( GetVarLogY( variable_name) ) c->SetLogy();
   
   for( int i_fs = 0; i_fs < Settings::fs4l ; i_fs++ )
   {
      histos_1D[Settings::reg3P1F][Settings::WZ][i_fs][Settings::inclusive_stxs]   ->SetFillColor(kMagenta);
      histos_1D[Settings::reg3P1F][Settings::qqZZ][i_fs][Settings::inclusive_stxs] ->SetFillColor(kCyan+1);
      histos_1D[Settings::reg3P1F][Settings::DY][i_fs][Settings::inclusive_stxs]   ->SetFillColor(kGreen-1);
      histos_1D[Settings::reg3P1F][Settings::ttbar][i_fs][Settings::inclusive_stxs]->SetFillColor(kBlue-2);
      
      histos_1D[Settings::reg3P1F][Settings::WZ][i_fs][Settings::inclusive_stxs]   ->SetLineColor(kMagenta);
      histos_1D[Settings::reg3P1F][Settings::qqZZ][i_fs][Settings::inclusive_stxs] ->SetLineColor(kCyan+1);
      histos_1D[Settings::reg3P1F][Settings::DY][i_fs][Settings::inclusive_stxs]   ->SetLineColor(kGreen-1);
      histos_1D[Settings::reg3P1F][Settings::ttbar][i_fs][Settings::inclusive_stxs]->SetLineColor(kBlue-2);
      
      h_from2P2F_3P1F[Settings::nominal][i_fs][Settings::inclusive_stxs]->SetLineColor(kRed);
      
      histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive_stxs]->SetMarkerSize(0.8);
      histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive_stxs]->SetMarkerStyle(20);
      histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive_stxs]->SetBinErrorOption(TH1::kPoisson);
      histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive_stxs]->SetLineColor(kBlack);
      
      THStack *stack = new THStack( "stack", "stack" );
		stack->Add(histos_1D[Settings::reg3P1F][Settings::qqZZ][i_fs][Settings::inclusive_stxs]);
		stack->Add(histos_1D[Settings::reg3P1F][Settings::WZ][i_fs][Settings::inclusive_stxs]);
		stack->Add(histos_1D[Settings::reg3P1F][Settings::ttbar][i_fs][Settings::inclusive_stxs]);
      stack->Add(histos_1D[Settings::reg3P1F][Settings::DY][i_fs][Settings::inclusive_stxs]);
		
      stack->Draw("HIST");
      
      h_from2P2F_3P1F[Settings::nominal][i_fs][Settings::inclusive_stxs]->Draw("HIST SAME");
      
      float data_max = histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive_stxs]->GetBinContent(histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive_stxs]->GetMaximumBin());
      float data_max_error = histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive_stxs]->GetBinErrorUp(histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive_stxs]->GetMaximumBin());
      
      stack->SetMinimum(1e-5);
      stack->SetMaximum((data_max + data_max_error)*1.35);
      
      TString _fs_label;
      if ( i_fs == Settings::fs4e) _fs_label = "m_{4#font[12]{e}} (GeV)";
      if ( i_fs == Settings::fs4mu) _fs_label = "m_{4#font[12]{#mu}} (GeV)";
      if ( i_fs == Settings::fs2e2mu) _fs_label = "m_{2#font[12]{e}2#font[12]{#mu}} (GeV)";
      if ( i_fs == Settings::fs2mu2e) _fs_label = "m_{2#font[12]{#mu}2#font[12]{e}} (GeV)";
      stack->GetXaxis()->SetTitle(_fs_label);
      stack->GetXaxis()->SetTitleSize(0.04);
      stack->GetXaxis()->SetLabelSize(0.04);
      stack->GetYaxis()->SetTitle(histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive_stxs]->GetYaxis()->GetTitle());
      stack->GetYaxis()->SetTitleSize(0.04);
      stack->GetYaxis()->SetLabelSize(0.04);
      
      stack->GetXaxis()->SetTitleOffset(1.2);
      stack->GetYaxis()->SetTitleOffset(1.25);
      
      histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive_stxs]->Draw("SAME p E1 X0");
      
      TLegend *legend;
      legend  = CreateLegend_3P1F("right",histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive_stxs],h_from2P2F_3P1F[Settings::nominal][i_fs][Settings::inclusive_stxs],histos_1D[Settings::reg3P1F][Settings::WZ][i_fs][Settings::inclusive_stxs],histos_1D[Settings::reg3P1F][Settings::qqZZ][i_fs][Settings::inclusive_stxs],histos_1D[Settings::reg3P1F][Settings::DY][i_fs][Settings::inclusive_stxs],histos_1D[Settings::reg3P1F][Settings::ttbar][i_fs][Settings::inclusive_stxs]);
      legend->Draw();
      
      // Draw lumi
      CMS_lumi *lumi = new CMS_lumi;
      lumi->set_lumi(c, _lumi, 0);
      
      TString _out_file_name;
      _out_file_name = folder + "/" + variable_name + "_OS_3P1F_" + _s_final_state.at(i_fs) + "_" + _s_category_stxs.at(Settings::inclusive_stxs);
      SavePlots(c, _out_file_name);
      
   }
}
//========================================================================================================

//========================================================================================================
void OSmethod::PlotDataMC( TString variable_name, TString folder )
{
   TCanvas *c;
   c = new TCanvas("OS", variable_name, 600, 600);
	
   if ( GetVarLogX( variable_name) ) c->SetLogx();
   if ( GetVarLogY( variable_name) ) c->SetLogy();
	
   for( int i_fs = 0; i_fs < Settings::fs4l ; i_fs++ )
   {
      histos_1D[Settings::regOS][Settings::WZ][i_fs][Settings::inclusive_stxs]   ->SetFillColor(kMagenta);
      histos_1D[Settings::regOS][Settings::qqZZ][i_fs][Settings::inclusive_stxs] ->SetFillColor(kCyan+1);
      histos_1D[Settings::regOS][Settings::DY][i_fs][Settings::inclusive_stxs]   ->SetFillColor(kGreen-1);
      histos_1D[Settings::regOS][Settings::ttbar][i_fs][Settings::inclusive_stxs]->SetFillColor(kBlue-2);
		
      histos_1D[Settings::regOS][Settings::WZ][i_fs][Settings::inclusive_stxs]   ->SetLineColor(kMagenta);
      histos_1D[Settings::regOS][Settings::qqZZ][i_fs][Settings::inclusive_stxs] ->SetLineColor(kCyan+1);
      histos_1D[Settings::regOS][Settings::DY][i_fs][Settings::inclusive_stxs]   ->SetLineColor(kGreen-1);
      histos_1D[Settings::regOS][Settings::ttbar][i_fs][Settings::inclusive_stxs]->SetLineColor(kBlue-2);
		
      histos_1D[Settings::regOS][Settings::Data][i_fs][Settings::inclusive_stxs]->SetMarkerSize(0.8);
      histos_1D[Settings::regOS][Settings::Data][i_fs][Settings::inclusive_stxs]->SetMarkerStyle(20);
      histos_1D[Settings::regOS][Settings::Data][i_fs][Settings::inclusive_stxs]->SetBinErrorOption(TH1::kPoisson);
      histos_1D[Settings::regOS][Settings::Data][i_fs][Settings::inclusive_stxs]->SetLineColor(kBlack);
		
      THStack *stack = new THStack( "stack", "stack" );
		stack->Add(histos_1D[Settings::regOS][Settings::qqZZ][i_fs][Settings::inclusive_stxs]);
		stack->Add(histos_1D[Settings::regOS][Settings::WZ][i_fs][Settings::inclusive_stxs]);
		stack->Add(histos_1D[Settings::regOS][Settings::ttbar][i_fs][Settings::inclusive_stxs]);
      stack->Add(histos_1D[Settings::regOS][Settings::DY][i_fs][Settings::inclusive_stxs]);
		
      stack->Draw("HIST");
		
      float data_max = histos_1D[Settings::regOS][Settings::Data][i_fs][Settings::inclusive_stxs]->GetBinContent(histos_1D[Settings::regOS][Settings::Data][i_fs][Settings::inclusive_stxs]->GetMaximumBin());
      float data_max_error = histos_1D[Settings::regOS][Settings::Data][i_fs][Settings::inclusive_stxs]->GetBinErrorUp(histos_1D[Settings::regOS][Settings::Data][i_fs][Settings::inclusive_stxs]->GetMaximumBin());
		
      stack->SetMinimum(1e-5);
      stack->SetMaximum((data_max + data_max_error)*1.35);
		
      TString _fs_label;
      if ( i_fs == Settings::fs4e) _fs_label = "m_{4#font[12]{e}} (GeV)";
      if ( i_fs == Settings::fs4mu) _fs_label = "m_{4#font[12]{#mu}} (GeV)";
      if ( i_fs == Settings::fs2e2mu) _fs_label = "m_{2#font[12]{e}2#font[12]{#mu}} (GeV)";
      if ( i_fs == Settings::fs2mu2e) _fs_label = "m_{2#font[12]{#mu}2#font[12]{e}} (GeV)";
      stack->GetXaxis()->SetTitle(_fs_label);
      stack->GetXaxis()->SetTitleSize(0.04);
      stack->GetXaxis()->SetLabelSize(0.04);
      stack->GetYaxis()->SetTitle(histos_1D[Settings::regOS][Settings::Data][i_fs][Settings::inclusive_stxs]->GetYaxis()->GetTitle());
      stack->GetYaxis()->SetTitleSize(0.04);
      stack->GetYaxis()->SetLabelSize(0.04);
		
      stack->GetXaxis()->SetTitleOffset(1.2);
      stack->GetYaxis()->SetTitleOffset(1.25);
		
      histos_1D[Settings::regOS][Settings::Data][i_fs][Settings::inclusive_stxs]->Draw("SAME p E1 X0");
		
      TLegend *legend;
      legend  = CreateLegend_2P2F("right",histos_1D[Settings::regOS][Settings::Data][i_fs][Settings::inclusive_stxs],histos_1D[Settings::regOS][Settings::WZ][i_fs][Settings::inclusive_stxs],histos_1D[Settings::regOS][Settings::qqZZ][i_fs][Settings::inclusive_stxs],histos_1D[Settings::regOS][Settings::DY][i_fs][Settings::inclusive_stxs],histos_1D[Settings::regOS][Settings::ttbar][i_fs][Settings::inclusive_stxs]);
      legend->Draw();
		
      // Draw lumi
      CMS_lumi *lumi = new CMS_lumi;
      lumi->set_lumi(c, _lumi, 0);
		
      TString _out_file_name;
      _out_file_name = folder + "/" + variable_name + "_OS_" + _s_final_state.at(i_fs) + "_" + _s_category_stxs.at(Settings::inclusive_stxs);
      SavePlots(c, _out_file_name);
		
   }
}
//========================================================================================================


//========================================================================================================
void OSmethod::PlotZXContributions( TString folder )
{
   TCanvas *c, *c_zx;
   TString _out_file_name;
   CMS_lumi *lumi = new CMS_lumi;

   c    = new TCanvas("c", "c", 600, 600);
   c_zx = new TCanvas("c_zx", "c_zx", 600, 600);
	
	for( int i_fs = 0; i_fs <= Settings::fs4l ; i_fs++ )
   {
		for ( int i_cat = 0; i_cat <= Settings::inclusive_stxs; i_cat++ )
      {
			c->cd();
			
			h_from3P1F_SR[Settings::nominal][i_fs][i_cat]       ->SetLineColor(kBlue);
			h_from2P2F_SR[Settings::nominal][i_fs][i_cat]       ->SetLineColor(kYellow);
			h_from3P1F_SR_final[Settings::nominal][i_fs][i_cat] ->SetLineColor(kBlack);
			h_from3P1F_SR_ZZonly[Settings::nominal][i_fs][i_cat]->SetLineColor(kRed);
			histos_ZX[Settings::nominal][i_fs][i_cat]           ->SetLineColor(kGreen);
			
			h_from3P1F_SR[Settings::nominal][i_fs][i_cat]->SetMinimum(0.0);
			
			h_from3P1F_SR[Settings::nominal][i_fs][i_cat]       ->Draw("HIST");
			h_from2P2F_SR[Settings::nominal][i_fs][i_cat]       ->Draw("HIST SAME");
			h_from3P1F_SR_final[Settings::nominal][i_fs][i_cat] ->Draw("HIST SAME");
			h_from3P1F_SR_ZZonly[Settings::nominal][i_fs][i_cat]->Draw("HIST SAME");
			histos_ZX[Settings::nominal][i_fs][i_cat]           ->Draw("HIST SAME");
			
			TString _fs_label;
			if ( i_fs == Settings::fs4e)    _fs_label = "m_{4#font[12]{e}} (GeV)";
			if ( i_fs == Settings::fs4mu)   _fs_label = "m_{4#font[12]{#mu}} (GeV)";
			if ( i_fs == Settings::fs2e2mu) _fs_label = "m_{2#font[12]{e}2#font[12]{#mu}} (GeV)";
			if ( i_fs == Settings::fs2mu2e) _fs_label = "m_{2#font[12]{#mu}2#font[12]{e}} (GeV)";
			if ( i_fs == Settings::fs4l)    _fs_label = "m_{4#font[12]{l}} (GeV)";
			h_from3P1F_SR[Settings::nominal][i_fs][i_cat]->GetXaxis()->SetTitle(_fs_label);
			h_from3P1F_SR[Settings::nominal][i_fs][i_cat]->GetXaxis()->SetTitleSize(0.04);
			h_from3P1F_SR[Settings::nominal][i_fs][i_cat]->GetXaxis()->SetLabelSize(0.04);
			h_from3P1F_SR[Settings::nominal][i_fs][i_cat]->GetYaxis()->SetTitle(h_from2P2F_SR[Settings::nominal][i_fs][i_cat]->GetYaxis()->GetTitle());
			h_from3P1F_SR[Settings::nominal][i_fs][i_cat]->GetYaxis()->SetTitleSize(0.04);
			h_from3P1F_SR[Settings::nominal][i_fs][i_cat]->GetYaxis()->SetLabelSize(0.04);
			
			h_from3P1F_SR[Settings::nominal][i_fs][i_cat]->GetXaxis()->SetTitleOffset(1.2);
			h_from3P1F_SR[Settings::nominal][i_fs][i_cat]->GetYaxis()->SetTitleOffset(1.25);
			
			TLegend *legend;
			legend  = CreateLegend_ZXcontr( "right", h_from2P2F_SR[Settings::nominal][i_fs][i_cat], h_from3P1F_SR[Settings::nominal][i_fs][i_cat],h_from3P1F_SR_ZZonly[Settings::nominal][i_fs][i_cat],h_from3P1F_SR_final[Settings::nominal][i_fs][i_cat],histos_ZX[Settings::nominal][i_fs][i_cat] );
			legend->Draw();
			
			// Draw lumi
			lumi->set_lumi(c, _lumi, 0);
			
			_out_file_name = folder + "/" + "ZX_Contributions_OS_" + _s_final_state.at(i_fs) + "_" + _s_category_stxs.at(i_cat);
			SavePlots(c, _out_file_name);
			
			c_zx->cd();
			
			histos_ZX[Settings::nominal][i_fs][i_cat]->SetLineColor(kGreen-1);
			histos_ZX[Settings::nominal][i_fs][i_cat]->SetFillColor(kGreen-1);
			histos_ZX[Settings::nominal][i_fs][i_cat]->Draw("HIST");
			lumi->set_lumi(c_zx, _lumi, 0);
			
			_out_file_name = folder + "/" + "ZX_OS_" + _s_final_state.at(i_fs) + "_" + _s_category_stxs.at(i_cat);
			SavePlots(c_zx, _out_file_name);
      }
   }
	
	
}
//========================================================================================================


//========================================================================================================
void OSmethod::FitZX( TString folder )
{
   TCanvas *c_zx;
   CMS_lumi *lumi = new CMS_lumi;
   TF1  *fit_function;
   TString _out_file_name;
   c_zx = new TCanvas("c_zx", "c_zx", 600, 600);
	
	for( int i_fs = 0; i_fs <= Settings::fs4l ; i_fs++ )
   {
		for ( int i_cat = 0; i_cat <= Settings::inclusive_stxs; i_cat++ )
      {
			c_zx->cd();
			
			gStyle->SetOptFit();
			gStyle->SetStatY(0.85);
			gStyle->SetStatX(0.95);
			gStyle->SetStatW(0.2);
			gStyle->SetStatH(0.2);

			fit_function = new TF1("fit_function","[0]*TMath::Landau(x, [1], [2]) + [3]*TMath::Landau(x, [4], [5])",70,1000);
			fit_function->SetParNames("Constant_{1}","MPV_{1}","#sigma_{1}","Constant_{2}","MPV_{2}","#sigma_{2}");
			fit_function->SetParameter(0,1.);
			fit_function->SetParameter(1,100.);
			fit_function->SetParameter(2,10.);
			fit_function->SetParameter(3,1.);
			fit_function->SetParameter(4,100.);
			fit_function->SetParameter(5,10.);
			
			histos_ZX[Settings::nominal][i_fs][i_cat]->Fit("fit_function");
			histos_ZX[Settings::nominal][i_fs][i_cat]->Draw("");
			
			// Draw lumi
			lumi->set_lumi(c_zx, _lumi, 0);
			
			
			_out_file_name = folder + "/" + "ZX_OS_fit_" + _s_final_state.at(i_fs) + "_" + _s_category_stxs.at(i_cat);
			SavePlots(c_zx, _out_file_name);
      }
   }
	gStyle->SetOptFit(0);
}
//========================================================================================================




//===============================================================
void OSmethod::SubtractWZ()
{
   for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
   {
      passing[Settings::Total][i_flav]->Add(passing[Settings::WZ][i_flav], -1.);
      failing[Settings::Total][i_flav]->Add(failing[Settings::WZ][i_flav], -1.);
   }

   cout << "[INFO] WZ contribution subtracted." << endl;
   
}
//===============================================================

//===============================================================
void OSmethod::ProduceFakeRates( TString file_name )
{
   for(int i_pT_bin = 0; i_pT_bin < _n_pT_bins - 1; i_pT_bin++ )
   {
      double temp_NP = 0;
      double temp_NF = 0;

      double temp_error_NP = 0;
      double temp_error_NF = 0;
      
      for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
      {
         if ( i_flav == Settings::ele && i_pT_bin == 0) continue; // electrons do not have 5 - 7 GeV bin
         temp_NP = passing[Settings::Total][i_flav]->IntegralAndError(passing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),passing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 1, 1, temp_error_NP);
         temp_NF = failing[Settings::Total][i_flav]->IntegralAndError(failing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),failing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 1, 1, temp_error_NF);
         
//         cout << "========================================" << endl;
//         cout << "pT bin = " << _pT_bins[i_pT_bin] << endl;
//         cout << "NP = " << temp_NP << endl;
//         cout << "error NP = " << temp_error_NP << endl;
//         cout << "NF = " << temp_NF << endl;
//         cout << "error NF = " << temp_error_NF << endl;
//         cout << "X = " << (_pT_bins[i_pT_bin] + _pT_bins[i_pT_bin + 1])/2 << endl;
//         cout << "error X = " << (_pT_bins[i_pT_bin + 1] - _pT_bins[i_pT_bin])/2 << endl;
//         cout << "Y = " << temp_NP/(temp_NP+temp_NF) << endl;
//         cout << "error Y = " << sqrt(pow((temp_NF/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NP,2) + pow((temp_NP/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NF,2)) << endl;
         
         vector_X[Settings::corrected][Settings::EB][i_flav].push_back((_pT_bins[i_pT_bin] + _pT_bins[i_pT_bin + 1])/2);
         vector_Y[Settings::corrected][Settings::EB][i_flav].push_back(temp_NP/(temp_NP+temp_NF));
         
         vector_EX[Settings::corrected][Settings::EB][i_flav].push_back((_pT_bins[i_pT_bin + 1] - _pT_bins[i_pT_bin])/2);
         vector_EY[Settings::corrected][Settings::EB][i_flav].push_back(sqrt(pow((temp_NF/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NP,2) + pow((temp_NP/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NF,2)));
         
         temp_NP = passing[Settings::Total][i_flav]->IntegralAndError(passing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),passing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 2, 2, temp_error_NP);
         temp_NF = failing[Settings::Total][i_flav]->IntegralAndError(failing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),failing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 2, 2, temp_error_NF);
         
         vector_X[Settings::corrected][Settings::EE][i_flav].push_back((_pT_bins[i_pT_bin] + _pT_bins[i_pT_bin + 1])/2);
         vector_Y[Settings::corrected][Settings::EE][i_flav].push_back(temp_NP/(temp_NP+temp_NF));
         
         vector_EX[Settings::corrected][Settings::EE][i_flav].push_back((_pT_bins[i_pT_bin + 1] - _pT_bins[i_pT_bin])/2);
         vector_EY[Settings::corrected][Settings::EE][i_flav].push_back(sqrt(pow((temp_NF/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NP,2) + pow((temp_NP/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NF,2)));
         
         // Just for fake rate plots calculate the same for histograms without WZ subtraction
         temp_NP = passing[Settings::Data][i_flav]->IntegralAndError(passing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),passing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 1, 1, temp_error_NP);
         temp_NF = failing[Settings::Data][i_flav]->IntegralAndError(failing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),failing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 1, 1, temp_error_NF);
         
         //         cout << "========================================" << endl;
         //         cout << "pT bin = " << _pT_bins[i_pT_bin] << endl;
         //         cout << "NP = " << temp_NP << endl;
         //         cout << "error NP = " << temp_error_NP << endl;
         //         cout << "NF = " << temp_NF << endl;
         //         cout << "error NF = " << temp_error_NF << endl;
         //         cout << "X = " << (_pT_bins[i_pT_bin] + _pT_bins[i_pT_bin + 1])/2 << endl;
         //         cout << "error X = " << (_pT_bins[i_pT_bin + 1] - _pT_bins[i_pT_bin])/2 << endl;
         //         cout << "Y = " << temp_NP/(temp_NP+temp_NF) << endl;
         //         cout << "error Y = " << sqrt(pow((temp_NF/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NP,2) + pow((temp_NP/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NF,2)) << endl;
         
         vector_X[Settings::uncorrected][Settings::EB][i_flav].push_back((_pT_bins[i_pT_bin] + _pT_bins[i_pT_bin + 1])/2);
         vector_Y[Settings::uncorrected][Settings::EB][i_flav].push_back(temp_NP/(temp_NP+temp_NF));
         
         vector_EX[Settings::uncorrected][Settings::EB][i_flav].push_back((_pT_bins[i_pT_bin + 1] - _pT_bins[i_pT_bin])/2);
         vector_EY[Settings::uncorrected][Settings::EB][i_flav].push_back(sqrt(pow((temp_NF/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NP,2) + pow((temp_NP/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NF,2)));
         
         temp_NP = passing[Settings::Data][i_flav]->IntegralAndError(passing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),passing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 2, 2, temp_error_NP);
         temp_NF = failing[Settings::Data][i_flav]->IntegralAndError(failing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),failing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 2, 2, temp_error_NF);
         
         vector_X[Settings::uncorrected][Settings::EE][i_flav].push_back((_pT_bins[i_pT_bin] + _pT_bins[i_pT_bin + 1])/2);
         vector_Y[Settings::uncorrected][Settings::EE][i_flav].push_back(temp_NP/(temp_NP+temp_NF));
         
         vector_EX[Settings::uncorrected][Settings::EE][i_flav].push_back((_pT_bins[i_pT_bin + 1] - _pT_bins[i_pT_bin])/2);
         vector_EY[Settings::uncorrected][Settings::EE][i_flav].push_back(sqrt(pow((temp_NF/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NP,2) + pow((temp_NP/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NF,2)));

      }
   }
   
   FR_OS_electron_EB = new TGraphErrors (vector_X[Settings::corrected][Settings::EB][Settings::ele].size(),
                                         &(vector_X[Settings::corrected][Settings::EB][Settings::ele][0]),
                                         &(vector_Y[Settings::corrected][Settings::EB][Settings::ele][0]),
                                         &(vector_EX[Settings::corrected][Settings::EB][Settings::ele][0]),
                                         &(vector_EY[Settings::corrected][Settings::EB][Settings::ele][0]));
   FR_OS_electron_EB->SetName("FR_OS_electron_EB");
   
   FR_OS_electron_EE = new TGraphErrors (vector_X[Settings::corrected][Settings::EE][Settings::ele].size(),
                                         &(vector_X[Settings::corrected][Settings::EE][Settings::ele][0]),
                                         &(vector_Y[Settings::corrected][Settings::EE][Settings::ele][0]),
                                         &(vector_EX[Settings::corrected][Settings::EE][Settings::ele][0]),
                                         &(vector_EY[Settings::corrected][Settings::EE][Settings::ele][0]));
   FR_OS_electron_EE->SetName("FR_OS_electron_EE");
   
   FR_OS_muon_EB = new TGraphErrors (vector_X[Settings::corrected][Settings::EB][Settings::mu].size(),
                                     &(vector_X[Settings::corrected][Settings::EB][Settings::mu][0]),
                                     &(vector_Y[Settings::corrected][Settings::EB][Settings::mu][0]),
                                     &(vector_EX[Settings::corrected][Settings::EB][Settings::mu][0]),
                                     &(vector_EY[Settings::corrected][Settings::EB][Settings::mu][0]));
   FR_OS_muon_EB->SetName("FR_OS_muon_EB");
   
   FR_OS_muon_EE = new TGraphErrors (vector_X[Settings::corrected][Settings::EE][Settings::mu].size(),
                                     &(vector_X[Settings::corrected][Settings::EE][Settings::mu][0]),
                                     &(vector_Y[Settings::corrected][Settings::EE][Settings::mu][0]),
                                     &(vector_EX[Settings::corrected][Settings::EE][Settings::mu][0]),
                                     &(vector_EY[Settings::corrected][Settings::EE][Settings::mu][0]));
   FR_OS_muon_EE->SetName("FR_OS_muon_EE");
   
   FR_OS_electron_EB_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EB][Settings::ele].size(),
                                         &(vector_X[Settings::uncorrected][Settings::EB][Settings::ele][0]),
                                         &(vector_Y[Settings::uncorrected][Settings::EB][Settings::ele][0]),
                                         &(vector_EX[Settings::uncorrected][Settings::EB][Settings::ele][0]),
                                         &(vector_EY[Settings::uncorrected][Settings::EB][Settings::ele][0]));
   FR_OS_electron_EB_unc->SetName("FR_OS_electron_EB_unc");
   
   FR_OS_electron_EE_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EE][Settings::ele].size(),
                                         &(vector_X[Settings::uncorrected][Settings::EE][Settings::ele][0]),
                                         &(vector_Y[Settings::uncorrected][Settings::EE][Settings::ele][0]),
                                         &(vector_EX[Settings::uncorrected][Settings::EE][Settings::ele][0]),
                                         &(vector_EY[Settings::uncorrected][Settings::EE][Settings::ele][0]));
   FR_OS_electron_EE_unc->SetName("FR_OS_electron_EE_unc");
   
   FR_OS_muon_EB_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EB][Settings::mu].size(),
                                     &(vector_X[Settings::uncorrected][Settings::EB][Settings::mu][0]),
                                     &(vector_Y[Settings::uncorrected][Settings::EB][Settings::mu][0]),
                                     &(vector_EX[Settings::uncorrected][Settings::EB][Settings::mu][0]),
                                     &(vector_EY[Settings::uncorrected][Settings::EB][Settings::mu][0]));
   FR_OS_muon_EB_unc->SetName("FR_OS_muon_EB_unc");
   
   FR_OS_muon_EE_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EE][Settings::mu].size(),
                                     &(vector_X[Settings::uncorrected][Settings::EE][Settings::mu][0]),
                                     &(vector_Y[Settings::uncorrected][Settings::EE][Settings::mu][0]),
                                     &(vector_EX[Settings::uncorrected][Settings::EE][Settings::mu][0]),
                                     &(vector_EY[Settings::uncorrected][Settings::EE][Settings::mu][0]));
   FR_OS_muon_EE_unc->SetName("FR_OS_muon_EE_unc");
   
   PlotFR();
   
   TFile* fOutHistos = TFile::Open(file_name, "recreate");
   fOutHistos->cd();
   
   FR_OS_electron_EB->Write();
   FR_OS_electron_EE->Write();
   FR_OS_muon_EB->Write();
   FR_OS_muon_EE->Write();
   
   fOutHistos->Close();
   delete fOutHistos;
   
   cout << "[INFO] Fake rates produced and stored in a file." << endl;
}
//===============================================================

//===============================================================
void OSmethod::PlotFR()
{
   TCanvas *c_ele, *c_mu;
   c_ele = new TCanvas("FR_ele", "FR_ele", 600, 600);
   c_mu  = new TCanvas("FR_mu", "FR_mu", 600, 600);
   
   mg_electrons = new TMultiGraph();
   mg_muons = new TMultiGraph();
   
   mg_electrons->Add(FR_OS_electron_EB);
   FR_OS_electron_EB->SetLineColor(kBlue);
   FR_OS_electron_EB->SetLineStyle(2);
   FR_OS_electron_EB->SetMarkerSize(0);
   FR_OS_electron_EB->SetTitle("barel corrected");
   mg_electrons->Add(FR_OS_electron_EE);
   FR_OS_electron_EE->SetLineColor(kRed);
   FR_OS_electron_EE->SetLineStyle(2);
   FR_OS_electron_EE->SetMarkerSize(0);
   FR_OS_electron_EE->SetTitle("endcap corrected");
   mg_electrons->Add(FR_OS_electron_EB_unc);
   FR_OS_electron_EB_unc->SetLineColor(kBlue);
   FR_OS_electron_EB_unc->SetLineStyle(1);
   FR_OS_electron_EB_unc->SetMarkerSize(0);
   FR_OS_electron_EB_unc->SetTitle("barel uncorrected");
   mg_electrons->Add(FR_OS_electron_EE_unc);
   FR_OS_electron_EE_unc->SetLineColor(kRed);
   FR_OS_electron_EE_unc->SetLineStyle(1);
   FR_OS_electron_EE_unc->SetMarkerSize(0);
   FR_OS_electron_EE_unc->SetTitle("endcap uncorrected");
   
   mg_muons->Add(FR_OS_muon_EB);
   FR_OS_muon_EB->SetLineColor(kBlue);
   FR_OS_muon_EB->SetLineStyle(2);
   FR_OS_muon_EB->SetMarkerSize(0);
   FR_OS_muon_EB->SetTitle("barel corrected");
   mg_muons->Add(FR_OS_muon_EE);
   FR_OS_muon_EE->SetLineColor(kRed);
   FR_OS_muon_EE->SetLineStyle(2);
   FR_OS_muon_EE->SetMarkerSize(0);
   FR_OS_muon_EE->SetTitle("endcap corrected");
   mg_muons->Add(FR_OS_muon_EB_unc);
   FR_OS_muon_EB_unc->SetLineColor(kBlue);
   FR_OS_muon_EB_unc->SetLineStyle(1);
   FR_OS_muon_EB_unc->SetMarkerSize(0);
   FR_OS_muon_EB_unc->SetTitle("barel uncorrected");
   mg_muons->Add(FR_OS_muon_EE_unc);
   FR_OS_muon_EE_unc->SetLineColor(kRed);
   FR_OS_muon_EE_unc->SetLineStyle(1);
   FR_OS_muon_EE_unc->SetMarkerSize(0);
   FR_OS_muon_EE_unc->SetTitle("endcap uncorrected");
   
   
   gStyle->SetEndErrorSize(0);
   
   TLegend *leg_ele,*leg_mu;
   CMS_lumi *lumi = new CMS_lumi;

   c_ele->cd();
   lumi->set_lumi(c_ele, _lumi, 0);
   mg_electrons->Draw("AP");
   mg_electrons->GetXaxis()->SetTitle("p_{T} [GeV]");
   mg_electrons->GetYaxis()->SetTitle("Fake Rate");
   mg_electrons->SetTitle("Electron fake rate");
   mg_electrons->SetMaximum(0.35);
   leg_ele = CreateLegend_FR("left",FR_OS_electron_EB_unc,FR_OS_electron_EB,FR_OS_electron_EE_unc,FR_OS_electron_EE);
   leg_ele->Draw();
   SavePlots(c_ele, "Plots/FR_OS_electrons");
   
   c_mu->cd();
   lumi->set_lumi(c_mu, _lumi, 0);
   mg_muons->Draw("AP");
   mg_muons->GetXaxis()->SetTitle("p_{T} [GeV]");
   mg_muons->GetYaxis()->SetTitle("Fake Rate");
   mg_muons->SetTitle("Muon fake rate");
   mg_muons->SetMaximum(0.35);
   leg_mu = CreateLegend_FR("left",FR_OS_muon_EB_unc,FR_OS_muon_EB,FR_OS_muon_EE_unc,FR_OS_muon_EE);
   leg_mu->Draw();
   SavePlots(c_mu, "Plots/FR_OS_muons");
   
}
//===============================================================

//===============================================================
void OSmethod::RemoveNegativeBins1D(TH1F *h)
{
   for (int i_bin_x = 1; i_bin_x <= h->GetXaxis()->GetNbins(); i_bin_x++)
   {
      if( h->GetBinContent(i_bin_x) < 0.) h->SetBinContent(i_bin_x, 0);
   }
   
}
//===============================================================

//===============================================================
void OSmethod::RemoveNegativeBins2D(TH2F *h)
{
   for (int i_bin_x = 1; i_bin_x <= h->GetXaxis()->GetNbins(); i_bin_x++)
   {
      for (int i_bin_y = 1; i_bin_y <= h->GetYaxis()->GetNbins(); i_bin_y++)
      {
         if( h->GetBinContent(i_bin_x,i_bin_y) < 0.) h->SetBinContent(i_bin_x,i_bin_y,0);
      }
      
   }
   
}
//===============================================================

//===============================================================
void OSmethod::Set_pT_binning(int size, float *bins)
{
   _n_pT_bins = size;

   for (int i = 0; i < size; i++)
   {
      _pT_bins[i] = bins[i];
   }
}
//===============================================================

//===============================================================
void OSmethod::SetLumi(float lumi)
{
   _lumi = lumi;
}
//===============================================================


//==========================================================
int OSmethod::find_current_process( TString input_file_name )
{
   
   int current_process = -999;
   
   // Assign dataset to correct process
   if ( input_file_name.Contains("Data") )           current_process = Settings::Data;
   if ( input_file_name.Contains("WZ") )             current_process = Settings::WZ;
   if ( input_file_name.Contains("ZZTo4l") )         current_process = Settings::qqZZ;
   if ( input_file_name.Contains("DYJetsToLL") )     current_process = Settings::DY;
   if ( input_file_name.Contains("TTJets") )         current_process = Settings::ttbar;
   if ( input_file_name.Contains("TTTo2L2Nu") )      current_process = Settings::ttbar;
   
   return current_process;
}
//==========================================================


//=============================
int OSmethod::FindFinalState()
{
   int final_state = -999;

   if ( Z1Flav == -121 )
   {
      if ( Z2Flav == -121 )
         final_state = Settings::fs4e;
      else if ( Z2Flav == -169 )
         final_state = Settings::fs2e2mu;
      else
         cerr << "[ERROR] in event " << RunNumber << ":" << LumiNumber << ":" << EventNumber << ", Z2Flav = " << Z2Flav << endl;
      }
   else if ( Z1Flav == -169 )
   {
      if ( Z2Flav == -121 )
         final_state = Settings::fs2mu2e;
      else if ( Z2Flav == -169 )
         final_state = Settings::fs4mu;
      else
         cerr << "[ERROR] in event " << RunNumber << ":" << LumiNumber << ":" << EventNumber << ", Z2Flav = " << Z2Flav << endl;
   }
   else
   {
      cerr << "[ERROR] in event " << RunNumber << ":" << LumiNumber << ":" << EventNumber << ", Z1Flav = " << Z1Flav << endl;
   }
   
   return final_state;
}
//=============================


//=================================
float OSmethod::calculate_K_factor(TString input_file_name)
{
   
   float k_factor = 1;
   
   if ( input_file_name.Contains("ZZTo4l"))
   {
      k_factor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M; // As of Moriond2016
   }
   else if ( input_file_name.Contains("ggTo"))
   {
      k_factor = KFactor_QCD_ggZZ_Nominal; // as of Moriond2016
   }
   return k_factor;
}
//=================================

//===================================================
bool OSmethod::GetVarLogX ( TString variable_name )
{
   //=============
   // M4l
   //=============
   if(variable_name == "M4l")                return bool(Plots::M4l().var_log_x);

   else
   {
      cout << "[ERROR] Wrong variable name choosen!" << endl;
      abort();
      return bool(Plots::M4l().var_log_x);
   }
}
//===================================================



//===================================================
bool OSmethod::GetVarLogY ( TString variable_name )
{
   //=============
   // M4l
   //=============
   if(variable_name == "M4l")                return bool(Plots::M4l().var_log_y);

   else
   {
      cout << "[ERROR] Wrong variable name choosen!" << endl;
      abort();
      return bool(Plots::M4l().var_log_y);
   }
}
//===================================================

//=========================================================================================================
TLegend* OSmethod::CreateLegend_FR( string position, TGraphErrors *EB_unc, TGraphErrors *EB_cor,TGraphErrors *EE_unc,TGraphErrors *EE_cor )
{
   TLegend *leg;
   leg = new TLegend( .64, .65, .97, .9 );
   if(position == "right") leg = new TLegend( .64, .65, .97, .9 );
   else if(position == "left") leg = new TLegend(.18,.65,.51,.9);
   
   leg->AddEntry( EB_unc, "barrel uncorrected", "l" );
   leg->AddEntry( EB_cor, "barrel corrected","l");
   leg->AddEntry( EE_unc, "endcap uncorrected", "l" );
   leg->AddEntry( EE_cor, "endcap corrected", "l" );
   
   return leg;
}
//=========================================================================================================

//=========================================================================================================
TLegend* OSmethod::CreateLegend_ZXcontr( string position, TH1F *h_2P2F_SR, TH1F *h_3P1F_SR,TH1F *h_3P1F_ZZ,TH1F *h_3P1F_SR_final,TH1F *total )
{
   TLegend *leg;
   leg = new TLegend( .64, .65, .97, .9 );
   if(position == "right") leg = new TLegend( .64, .65, .97, .9 );
   else if(position == "left") leg = new TLegend(.18,.65,.51,.9);
   
   leg->AddEntry( h_2P2F_SR, "2P2F", "l" );
   leg->AddEntry( h_3P1F_SR, "3P1F w/o removal","l");
   leg->AddEntry( h_3P1F_ZZ, "3P1F ZZ contr.", "l" );
   leg->AddEntry( h_3P1F_SR_final, "3P1F final", "l" );
   leg->AddEntry( total, "Z+X final", "l" );
   
   return leg;
}
//=========================================================================================================

//=========================================================================================================
TLegend* OSmethod::CreateLegend_2P2F( string position, TH1F *data, TH1F *WZ,TH1F *qqZZ,TH1F *DY,TH1F *ttbar )
{
   TLegend *leg;
   leg = new TLegend( .64, .65, .97, .9 );
   if(position == "right") leg = new TLegend( .64, .65, .97, .9 );
   else if(position == "left") leg = new TLegend(.18,.65,.51,.9);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   
   leg->AddEntry( data, "Data", "p" );
   leg->AddEntry( WZ,"WZ","f");
   leg->AddEntry( qqZZ, "Z#gamma*, ZZ", "f" );
   leg->AddEntry( DY, "Z + jets", "f" );
   leg->AddEntry( ttbar, "t#bar{t} + jets", "f" );
   
   return leg;
}
//=========================================================================================================

//=========================================================================================================
TLegend* OSmethod::CreateLegend_3P1F( string position, TH1F *data, TH1F *h_2P2F, TH1F *WZ,TH1F *qqZZ,TH1F *DY,TH1F *ttbar )
{
   TLegend *leg;
   leg = new TLegend( .64, .65, .97, .9 );
   if(position == "right") leg = new TLegend( .64, .65, .97, .9 );
   else if(position == "left") leg = new TLegend(.18,.65,.51,.9);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   
   leg->AddEntry( data, "Data", "p" );
   leg->AddEntry( h_2P2F, "2P2F extr.", "l" );
   leg->AddEntry( WZ,"WZ","f");
   leg->AddEntry( qqZZ, "Z#gamma*, ZZ", "f" );
   leg->AddEntry( DY, "Z + jets", "f" );
   leg->AddEntry( ttbar, "t#bar{t} + jets", "f" );
   
   return leg;
}
//=========================================================================================================



//=======================================
void OSmethod::SavePlots( TCanvas *c, TString name)
{
   c->SaveAs(name + ".pdf");
   c->SaveAs(name + ".root");
   c->SaveAs(name + ".eps");
   gSystem->Exec("convert -density 300 -quality 100 " + name + ".eps " + name + ".png");
}
//=======================================





