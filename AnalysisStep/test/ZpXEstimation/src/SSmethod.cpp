// Include classes
#include <ZZAnalysis/AnalysisStep/test/ZpXEstimation/include/SSmethod.h>

// Constructor
//============================================================
SSmethod::SSmethod():Tree()
{
   _current_process = -999;
   _current_final_state = -999;
   _current_category = -999;
	
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
   
   _s_region.push_back("ZLL");
   
   // Z+X SS factors
   // FIXME: recompute this for Run II, OS/SS ratio taken when computing fake rates in SS method
   _fs_ROS_SS.push_back(1.22);//4mu
   _fs_ROS_SS.push_back(0.97);//4e
   _fs_ROS_SS.push_back(1.30);//2e2mu
   _fs_ROS_SS.push_back(0.98);//2mu2e
   
   vector<float> temp;
   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
      {
         temp.push_back(0);
         _N_OS_events[i_fs][i_cat] = 0.;
         _N_SS_events[i_fs][i_cat] = 0.;
      }
      _expected_yield_SR.push_back(temp);
      _expected_yield_SR_up.push_back(temp);
      _expected_yield_SR_dn.push_back(temp);
      _number_of_events_CR.push_back(temp);
   }
   
   DeclareFRHistos();
   DeclareDataMCHistos();
   DeclareZXHistos();
}
//============================================================



// Destructor
//====================
SSmethod::~SSmethod()
{
}
//====================


//================================================================================================
void SSmethod::Calculate_SSOS_Ratio( TString input_file_data_name, TString input_file_MC_name , bool subtractMC )
{
   input_file_data = new TFile( input_file_data_name);
   input_file_MC   = new TFile( input_file_MC_name);
   
   hCounters = (TH1F*)input_file_MC->Get("CRZLLTree/Counters");
   gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
   
   //Loop over data CR to get the number of events in OS and SS
   input_tree_data = (TTree*)input_file_data->Get("CRZLLTree/candTree");
   Init( input_tree_data, input_file_data_name , true);
   _current_process = Settings::Data;
   
   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
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
														 true,// Use VHMET category
														 false);// Use QG tagging
      
      if ((test_bit(CRflag, CRZLLss))) _N_SS_events[_current_final_state][_current_category]+=1.0;
      if ((test_bit(CRflag, CRZLLos_2P2F)) || (test_bit(CRflag, CRZLLos_3P1F))) _N_OS_events[_current_final_state][_current_category]+=1.0;

      _k_factor = calculate_K_factor(input_file_data_name);
      _event_weight = (_lumi * 1000 * xsec * _k_factor * overallEventWeight) / gen_sum_weights;
   }
   
   //Loop over MC to estimate ZZTo4L events in OS
   if( subtractMC )
   {
      input_tree_MC = (TTree*)input_file_MC->Get("CRZLLTree/candTree");
      Init( input_tree_MC, input_file_MC_name , true);
      _current_process = Settings::qqZZ;
   
      if (fChain == 0) return;
   
      nentries = fChain->GetEntriesFast();
   
      nbytes = 0, nb = 0;
   
      for (Long64_t jentry=0; jentry<nentries;jentry++)
      {
         Long64_t ientry = LoadTree(jentry);
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry);
         nbytes += nb;
         
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
														 true,// Use VHMET category
														 false);// Use QG tagging
         
         _k_factor = calculate_K_factor(input_file_data_name);
         _event_weight = (_lumi * 1000 * xsec * _k_factor * overallEventWeight) / gen_sum_weights;
         
         if ((test_bit(CRflag, CRZLLos_2P2F)) || (test_bit(CRflag, CRZLLos_3P1F))) _N_OS_events[_current_final_state][_current_category]-=_event_weight;
      
      }
   }
   
   //Calculate inclusive numbers
   for (  int i_cat = 0; i_cat < num_of_categories - 1; i_cat++  )
   {
      for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++  )
      {
         _N_SS_events[Settings::fs4l][i_cat]     += _N_SS_events[i_fs][i_cat];   //calculate N events for inclusive 4l final state
         _N_OS_events[Settings::fs4l][i_cat]     += _N_OS_events[i_fs][i_cat];
         _N_SS_events[i_fs][Settings::inclusive] += _N_SS_events[i_fs][i_cat];   //calculate N events for inclusive category
         _N_OS_events[i_fs][Settings::inclusive] += _N_OS_events[i_fs][i_cat];
         
         if (false)//( MERGE_2E2MU )
         {
            _N_SS_events[Settings::fs2e2mu][i_cat]     += _N_SS_events[Settings::fs2mu2e][i_cat];   //merge 2e2mu and 2mu2e final state
            _N_OS_events[Settings::fs2e2mu][i_cat]     += _N_OS_events[Settings::fs2mu2e][i_cat];
            _N_SS_events[Settings::fs2mu2e][i_cat]      = 0.;
            _N_OS_events[Settings::fs2mu2e][i_cat]      = 0.;
         }
      }
   }
   for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++  )
   {
      _N_SS_events[Settings::fs4l][Settings::inclusive] += _N_SS_events[i_fs][Settings::inclusive];
      _N_OS_events[Settings::fs4l][Settings::inclusive] += _N_OS_events[i_fs][Settings::inclusive];
   }
   
   // Print Z + X expected yields for inclusive category
   cout << endl;
   cout << "========================================================================================" << endl;
   cout << "[INFO] Control printout of OS/SS ratio calculation." << endl;
   cout << "========================================================================================" << endl;
   for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++ )
   {
      if (false) continue;//( MERGE_2E2MU && i_fs == Settings::fs2mu2e) continue;
      cout << "Category: " << _s_category.at(Settings::inclusive) << "   Final state: " << _s_final_state.at(i_fs) << endl;
      cout << _N_OS_events[i_fs][Settings::inclusive]/_N_SS_events[i_fs][Settings::inclusive] << " +/- " << sqrt(1./_N_OS_events[i_fs][Settings::inclusive]
      + 1./_N_SS_events[i_fs][Settings::inclusive]) << endl;
   }
   
   cout << "[INFO] Total = " << _N_OS_events[Settings::fs4l][Settings::inclusive]/_N_SS_events[Settings::fs4l][Settings::inclusive] << " +/- " <<
   sqrt(1./_N_OS_events[Settings::fs4l][Settings::inclusive] + 1./_N_SS_events[Settings::fs4l][Settings::inclusive]) << endl;
   cout << "========================================================================================" << endl;
   cout << endl;
   
   if(true)
	{
		_fs_ROS_SS[Settings::fs4mu]   = _N_OS_events[Settings::fs4mu][Settings::inclusive]/_N_SS_events[Settings::fs4mu][Settings::inclusive];//4mu
		_fs_ROS_SS[Settings::fs4e]    = _N_OS_events[Settings::fs4e][Settings::inclusive]/_N_SS_events[Settings::fs4e][Settings::inclusive];//4e
		_fs_ROS_SS[Settings::fs2e2mu] = _N_OS_events[Settings::fs2e2mu][Settings::inclusive]/_N_SS_events[Settings::fs2e2mu][Settings::inclusive];//2e2mu
		_fs_ROS_SS[Settings::fs2mu2e] = _N_OS_events[Settings::fs2mu2e][Settings::inclusive]/_N_SS_events[Settings::fs2mu2e][Settings::inclusive];//2mu2e
	}

	
   cout << "[INFO] OS/SS ratios calculated." << endl;
   
   
}
//================================================================================================

//===============================================================================
void SSmethod::FillFRHistos( TString input_file_data_name )
{
   input_file_data = new TFile( input_file_data_name);
   
   hCounters = (TH1F*)input_file_data->Get("CRZLTree/Counters");
   gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
   
   input_tree_data = (TTree*)input_file_data->Get("CRZLTree/candTree");
   Init( input_tree_data, input_file_data_name , false);
   
   _current_process = find_current_process(input_file_data_name);
   
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   // Define some counters for control print out
	Int_t _total_events[num_of_final_states];
	Int_t _failZ1MassCut[num_of_final_states];
	Int_t _failLepPtCut[num_of_final_states];
	Int_t _failSIPCut[num_of_final_states];
	Int_t _failMETCut[num_of_final_states];
	Int_t _passingSelection[num_of_final_states];
	Int_t _faillingSelection[num_of_final_states];
	
	for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
	{
		_total_events[i_fs] = 0.;
		_failZ1MassCut[i_fs] = 0.;
		_failLepPtCut[i_fs] = 0.;
		_failSIPCut[i_fs] = 0.;
		_failMETCut[i_fs] = 0.;
		_passingSelection[i_fs] = 0.;
		_faillingSelection[i_fs] = 0.;
	}

   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;

		(fabs(LepLepId->at(2)) == 11) ? _total_events[Settings::ele]++ : _total_events[Settings::mu]++;
	   
		if ( Z1Mass < 40. ) {(fabs(LepLepId->at(2)) == 11) ? _failZ1MassCut[Settings::ele]++ : _failZ1MassCut[Settings::mu]++; continue;}
	   if ( Z1Mass > 120. ) {(fabs(LepLepId->at(2)) == 11) ? _failZ1MassCut[Settings::ele]++ : _failZ1MassCut[Settings::mu]++; continue;}
	   if ( (LepPt->at(0) > LepPt->at(1)) && (LepPt->at(0) < 20. || LepPt->at(1) < 10.) ) {(fabs(LepLepId->at(2)) == 11) ? _failLepPtCut[Settings::ele]++ : _failLepPtCut[Settings::mu]++; continue;}
	   if ( (LepPt->at(1) > LepPt->at(0)) && (LepPt->at(1) < 20. || LepPt->at(0) < 10.) ) {(fabs(LepLepId->at(2)) == 11) ? _failLepPtCut[Settings::ele]++ : _failLepPtCut[Settings::mu]++; continue;}
	   if ( LepSIP->at(2) > 4.) {(fabs(LepLepId->at(2)) == 11) ? _failSIPCut[Settings::ele]++ : _failSIPCut[Settings::mu]++; continue;}
	   if ( PFMET > 25. ) {(fabs(LepLepId->at(2)) == 11) ? _failMETCut[Settings::ele]++ : _failMETCut[Settings::mu]++; continue;}
      else
	   {
         // Final event weight
         _k_factor = calculate_K_factor(input_file_data_name);
         _event_weight = (_lumi * 1000 * xsec * _k_factor * overallEventWeight) / gen_sum_weights;

         if(LepisID->at(2) && ((fabs(LepLepId->at(2)) == 11) ? LepCombRelIsoPF->at(2) < 999999. : LepCombRelIsoPF->at(2) < 0.35))
         {
				(fabs(LepLepId->at(2)) == 11) ? _passingSelection[Settings::ele]++ : _passingSelection[Settings::mu]++;
            if(fabs(LepLepId->at(2)) == 11 ) passing[_current_process][Settings::ele]->Fill(LepPt->at(2), (abs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
            else if(fabs(LepLepId->at(2)) == 13 ) passing[_current_process][Settings::mu]->Fill(LepPt->at(2), (abs(LepEta->at(2)) < 1.2) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
         }
         else
         {
				(fabs(LepLepId->at(2)) == 11) ? _faillingSelection[Settings::ele]++ : _faillingSelection[Settings::mu]++;
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
		cout << "[INFO] Control printout for electrons in Z+L control region." << endl;
		cout << "========================================================================================" << endl;
		cout << "[INFO] Total number of events in Z+L control region = " << _total_events[Settings::ele] << endl;
		cout << "[INFO] Events after 40 < Z1 < 120 GeV cut  = " << _total_events[Settings::ele] - _failZ1MassCut[Settings::ele] << endl;
		cout << "[INFO] Events after LepPt > 20,10 GeV cut  = " << _total_events[Settings::ele] - _failZ1MassCut[Settings::ele] - _failLepPtCut[Settings::ele] << endl;
		cout << "[INFO] Events after SIP < 4 cut  = " << _total_events[Settings::ele] - _failZ1MassCut[Settings::ele] - _failLepPtCut[Settings::ele] - _failSIPCut[Settings::ele] << endl;
		cout << "[INFO] Events after MET < 25 cut  = " << _total_events[Settings::ele] - _failZ1MassCut[Settings::ele] - _failLepPtCut[Settings::ele] - _failSIPCut[Settings::ele] - _failMETCut[Settings::ele] << endl;
		cout << "[INFO] Total events left = " << _passingSelection[Settings::ele] + _faillingSelection[Settings::ele] << endl;
		cout << "[INFO] Passing selection = " << _passingSelection[Settings::ele]  << endl;
		cout << "[INFO] Failling selection = " << _faillingSelection[Settings::ele] << endl;
		cout << "========================================================================================" << endl;
		cout << endl;
		
		cout << endl;
		cout << "========================================================================================" << endl;
		cout << "[INFO] Control printout for muons in Z+L control region." << endl;
		cout << "========================================================================================" << endl;
		cout << "[INFO] Total number of events in Z+L control region = " << _total_events[Settings::mu] << endl;
		cout << "[INFO] Events after 40 < Z1 < 120 GeV cut  = " << _total_events[Settings::mu] - _failZ1MassCut[Settings::mu] << endl;
		cout << "[INFO] Events after LepPt > 20,10 GeV cut  = " << _total_events[Settings::mu] - _failZ1MassCut[Settings::mu] - _failLepPtCut[Settings::mu] << endl;
		cout << "[INFO] Events after SIP < 4 cut  = " << _total_events[Settings::mu] - _failZ1MassCut[Settings::mu] - _failLepPtCut[Settings::mu] - _failSIPCut[Settings::mu] << endl;
		cout << "[INFO] Events after MET < 25 cut  = " << _total_events[Settings::mu] - _failZ1MassCut[Settings::mu] - _failLepPtCut[Settings::mu] - _failSIPCut[Settings::mu] - _failMETCut[Settings::mu] << endl;
		cout << "[INFO] Total events left = " << _passingSelection[Settings::mu] + _faillingSelection[Settings::mu] << endl;
		cout << "[INFO] Passing selection = " << _passingSelection[Settings::mu]  << endl;
		cout << "[INFO] Failling selection = " << _faillingSelection[Settings::mu] << endl;
		cout << "========================================================================================" << endl;
		cout << endl;
	}
	
   cout << "[INFO] Processing of " << input_file_data_name << " done." << endl;
}
//===============================================================================



//===============================================================================
void SSmethod::FillDataMCPlots( TString input_file_data_name )
{
   input_file_data = new TFile( input_file_data_name);
   
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
      
      if (!(test_bit(CRflag, CRZLLss))) continue;
      
      
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
														 true,// Use VHMET category
														 false);// Use QG tagging

      
      _k_factor = calculate_K_factor(input_file_data_name);
      _event_weight = (_lumi * 1000 * xsec * _k_factor * overallEventWeight) / gen_sum_weights;
   
      histos_1D[Settings::regZLL][_current_process][_current_final_state][_current_category]->Fill(ZZMass,(_current_process == Settings::Data) ? 1 :  _event_weight);

   } // END events loop
   
   cout << "[INFO] Processing of " << input_file_data_name << " done." << endl;
}
//===============================================================================



//===============================================================================
void SSmethod::MakeHistogramsZX( TString input_file_data_name, TString  input_file_FR_name )
{
   
   FakeRates *FR = new FakeRates( input_file_FR_name );
   
   input_file_data = new TFile( input_file_data_name);
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
      
      if ( !CRflag ) continue;
      if ( !test_bit(CRflag, CRZLLss) ) continue;
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
														 true,// Use VHMET category
														 false);// Use QG tagging
      
      
      _k_factor = calculate_K_factor(input_file_data_name);
      _event_weight = (_lumi * 1000 * xsec * _k_factor * overallEventWeight) / gen_sum_weights;
      
      // Calculate yield
      _yield_SR = _fs_ROS_SS.at(_current_final_state)*FR->GetFakeRate(LepPt->at(2),LepEta->at(2),LepLepId->at(2))*FR->GetFakeRate(LepPt->at(3),LepEta->at(3),LepLepId->at(3));
      _yield_SR_up = _fs_ROS_SS.at(_current_final_state)*FR->GetFakeRate_Up(LepPt->at(2),LepEta->at(2),LepLepId->at(2))*FR->GetFakeRate_Up(LepPt->at(3),LepEta->at(3),LepLepId->at(3));
      _yield_SR_dn = _fs_ROS_SS.at(_current_final_state)*FR->GetFakeRate_Dn(LepPt->at(2),LepEta->at(2),LepLepId->at(2))*FR->GetFakeRate_Dn(LepPt->at(3),LepEta->at(3),LepLepId->at(3));
      
      
      _expected_yield_SR[_current_final_state][_current_category] += _yield_SR;
      _expected_yield_SR_up[_current_final_state][_current_category] += _yield_SR_up;
      _expected_yield_SR_dn[_current_final_state][_current_category] += _yield_SR_dn;
      _number_of_events_CR[_current_final_state][_current_category]++;
      //cout << _current_process << " " <<  _current_final_state << " " << _current_category << endl;
		
      // Fill m4l Z+X histograms
      histos_ZX[Settings::regZLL][_current_process][_current_final_state][_current_category]->Fill(ZZMass,(_current_process == Settings::Data) ? _yield_SR :  _yield_SR*_event_weight);
      

   } // End events loop
   
   for (  int i_cat = 0; i_cat < num_of_categories - 1; i_cat++  )
   {
      for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++  )
      {
         _expected_yield_SR[Settings::fs4l][i_cat]       += _expected_yield_SR[i_fs][i_cat];   //calculate expected yield for inclusive 4l final state
         _number_of_events_CR[Settings::fs4l][i_cat]     += _number_of_events_CR[i_fs][i_cat];
         _expected_yield_SR[i_fs][Settings::inclusive]   += _expected_yield_SR[i_fs][i_cat];   //calculate expected yield for inclusive category
         _expected_yield_SR_up[i_fs][Settings::inclusive]   += _expected_yield_SR_up[i_fs][i_cat];
         _expected_yield_SR_dn[i_fs][Settings::inclusive]   += _expected_yield_SR_dn[i_fs][i_cat];
         _number_of_events_CR[i_fs][Settings::inclusive] += _number_of_events_CR[i_fs][i_cat];
         
         if (false)//( MERGE_2E2MU )
         {
            _expected_yield_SR[Settings::fs2e2mu][i_cat]       += _expected_yield_SR[Settings::fs2mu2e][i_cat];   //merge 2e2mu and 2mu2e final state
            _expected_yield_SR_up[Settings::fs2e2mu][i_cat]       += _expected_yield_SR_up[Settings::fs2mu2e][i_cat];
            _expected_yield_SR_dn[Settings::fs2e2mu][i_cat]       += _expected_yield_SR_dn[Settings::fs2mu2e][i_cat];
            _number_of_events_CR[Settings::fs2e2mu][i_cat]     += _number_of_events_CR[Settings::fs2mu2e][i_cat];
            _expected_yield_SR[Settings::fs2mu2e][i_cat]        = 0.;
            _expected_yield_SR_up[Settings::fs2mu2e][i_cat]        = 0.;
            _expected_yield_SR_dn[Settings::fs2mu2e][i_cat]        = 0.;
            _number_of_events_CR[Settings::fs2mu2e][i_cat]      = 0.;
         }
      }
   }
   for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++  )
   {
      _expected_yield_SR[Settings::fs4l][Settings::inclusive] += _expected_yield_SR[i_fs][Settings::inclusive];
      _expected_yield_SR_up[Settings::fs4l][Settings::inclusive] += _expected_yield_SR_dn[i_fs][Settings::inclusive];
      _expected_yield_SR_up[Settings::fs4l][Settings::inclusive] += _expected_yield_SR_dn[i_fs][Settings::inclusive];
   }
   
   // Print Z + X expected yields and uncertainties
   cout << endl;
   cout << "===================================================================================================================================" << endl;
   cout << "[INFO] Control printout for Z+X yields in final states derived with this fake rate file: " << input_file_FR_name << endl;
   cout << "===================================================================================================================================" << endl;
   for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++ )
   {
      if (false) continue;//( MERGE_2E2MU && i_fs == Settings::fs2mu2e) continue;
      for ( int i_cat = 0; i_cat < num_of_categories; i_cat++)
      {
        float stat    = _expected_yield_SR[i_fs][i_cat]/sqrt(_number_of_events_CR[i_fs][i_cat]);
        float syst_up = _expected_yield_SR[i_fs][i_cat]*((_expected_yield_SR_up[i_fs][i_cat]/_expected_yield_SR[i_fs][i_cat]) - 1.);
        float syst_dn = _expected_yield_SR[i_fs][i_cat]*(1. - (_expected_yield_SR_dn[i_fs][i_cat]/_expected_yield_SR[i_fs][i_cat]));
        float comb_up = sqrt(stat*stat + syst_up*syst_up);
        float comb_dn = sqrt(stat*stat + syst_dn*syst_dn);
			
			
			cout << "Category: " << _s_category.at(i_cat) << "   Final state: " << _s_final_state.at(i_fs) << endl;
         cout << _expected_yield_SR[i_fs][i_cat] << " +/- " << comb_up << "(total):" << "  - " << stat << " (stat., evt: " <<
         _number_of_events_CR[i_fs][i_cat] << ")" << "   - " << syst_up << " (syst.)" <<  "  " << (1. - comb_dn/_expected_yield_SR[i_fs][i_cat]) << "/" << (1. + comb_up/_expected_yield_SR[i_fs][i_cat]) << endl;
		}
	cout << "==================================================================================================================================" << endl;
   }
   
   cout << "[INFO] Total = " << _expected_yield_SR[Settings::fs4l][Settings::inclusive] << endl;
   cout << "==================================================================================================================================" << endl;
   cout << endl;
	
   
   
   cout << "[INFO] Z+X histograms filled." << endl;
}
//===============================================================================




//===============================================================
void SSmethod::DeclareFRHistos()
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
void SSmethod::DeclareDataMCHistos()
{
   for (int i_reg = 0; i_reg < num_of_regions_ss; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
         {
            for (int i_cat = 0; i_cat < num_of_categories; i_cat++)
            {
               _histo_name = "M4l_" + _s_region.at(i_reg) + "_" + _s_process.at(i_proc) + "_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
               _histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
               histos_1D[i_reg][i_proc][i_fs][i_cat] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
            }
         }
      }
   }
   
}
//===============================================================

//===============================================================
void SSmethod::DeclareZXHistos()
{
   for (int i_reg = 0; i_reg < num_of_regions_ss; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
         {
            for (int i_cat = 0; i_cat < num_of_categories; i_cat++)
            {
               _histo_name = "ZX_M4l_" + _s_region.at(i_reg) + "_" + _s_process.at(i_proc) + "_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
               _histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
               histos_ZX[i_reg][i_proc][i_fs][i_cat] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
            }
         }
      }
   }
}
//===============================================================

//===============================================================
void SSmethod::SaveFRHistos( TString file_name, bool subtractWZ, bool remove_negative_bins)
{
   TFile* fOutHistos = new TFile(file_name, "recreate");
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
void SSmethod::SaveDataMCHistos( TString file_name )
{
   FillDataMCInclusive();
   
   TFile* fOutHistos = new TFile(file_name, "recreate");
   fOutHistos->cd();
   
   for (int i_reg = 0; i_reg < num_of_regions_ss; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
         {
            for (int i_cat = 0; i_cat < num_of_categories; i_cat++)
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
void SSmethod::FillDataMCInclusive( )
{
   for (int i_reg = 0; i_reg < num_of_regions_ss; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < Settings::fs4l; i_fs++)
         {
            for (int i_cat = 0; i_cat < Settings::inclusive; i_cat++)
            {
               histos_1D[i_reg][i_proc][i_fs][Settings::inclusive]->Add(histos_1D[i_reg][i_proc][i_fs][i_cat]);
               histos_1D[i_reg][i_proc][Settings::fs4l][i_cat]    ->Add(histos_1D[i_reg][i_proc][i_fs][i_cat]);
            }
         }
      }
   }
   
   for (int i_reg = 0; i_reg < num_of_regions_ss; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < Settings::fs4l; i_fs++)
         {
            histos_1D[i_reg][i_proc][Settings::fs4l][Settings::inclusive]->Add(histos_1D[i_reg][i_proc][i_fs][Settings::inclusive]);
         }
      }
   }
   
   cout << "[INFO] All Data/MC histograms summed." << endl;
}
//===============================================================

//===============================================================
void SSmethod::SaveZXHistos( TString file_name )
{
   FillZXInclusive();
   
   TFile* fOutHistos = new TFile(file_name, "recreate");
   fOutHistos->cd();
   
   for (int i_reg = 0; i_reg < num_of_regions_ss; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
         {
            for (int i_cat = 0; i_cat < num_of_categories; i_cat++)
            {
               histos_ZX[i_reg][i_proc][i_fs][i_cat]->Write();
            }
         }
      }
   }
   
   fOutHistos->Close();
   delete fOutHistos;
   
   cout << "[INFO] All Z+X histograms saved." << endl;
}
//===============================================================

//===============================================================
void SSmethod::FillZXInclusive( )
{
   for (int i_reg = 0; i_reg < num_of_regions_ss; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < Settings::fs4l; i_fs++)
         {
            for (int i_cat = 0; i_cat < Settings::inclusive; i_cat++)
            {
               histos_ZX[i_reg][i_proc][i_fs][Settings::inclusive]->Add(histos_ZX[i_reg][i_proc][i_fs][i_cat]);
               histos_ZX[i_reg][i_proc][Settings::fs4l][i_cat]    ->Add(histos_ZX[i_reg][i_proc][i_fs][i_cat]);
            }
         }
      }
   }
   
   for (int i_reg = 0; i_reg < num_of_regions_ss; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < Settings::fs4l; i_fs++)
         {
            histos_ZX[i_reg][i_proc][Settings::fs4l][Settings::inclusive]->Add(histos_ZX[i_reg][i_proc][i_fs][Settings::inclusive]);
         }
      }
   }
   
   cout << "[INFO] All Z+X histograms summed." << endl;
}
//===============================================================

//===============================================================
void SSmethod::GetFRHistos( TString file_name)
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
void SSmethod::GetDataMCHistos( TString file_name)
{
   TFile* histo_file = TFile::Open(file_name);
   
   for (int i_reg = 0; i_reg < num_of_regions_ss; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
         {
            for (int i_cat = 0; i_cat < num_of_categories; i_cat++)
            {
               _histo_name = "M4l_" + _s_region.at(i_reg) + "_" + _s_process.at(i_proc) + "_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
               histos_1D[i_reg][i_proc][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name);
            }
         }
      }
   }
   
   cout << "[INFO] All Data/MC histograms retrieved from file." << endl;
}

//===============================================================

//===============================================================
void SSmethod::GetZXHistos( TString file_name)
{
   TFile* histo_file = TFile::Open(file_name);
   
   for (int i_reg = 0; i_reg < num_of_regions_ss; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
         {
            for (int i_cat = 0; i_cat < num_of_categories; i_cat++)
            {
               _histo_name = "ZX_M4l_" + _s_region.at(i_reg) + "_" + _s_process.at(i_proc) + "_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
               histos_ZX[i_reg][i_proc][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name);
            }
         }
      }
   }
   
   cout << "[INFO] All Z+X histograms retrieved from file." << endl;
}

//===============================================================


//===============================================================
void SSmethod::ProduceFakeRates( TString file_name , TString input_file_data_name /*= "DONT_CORRECT"*/)
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
	
   FR_SS_electron_EB_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EB][Settings::ele].size(),
															&(vector_X[Settings::uncorrected][Settings::EB][Settings::ele][0]),
															&(vector_Y[Settings::uncorrected][Settings::EB][Settings::ele][0]),
															&(vector_EX[Settings::uncorrected][Settings::EB][Settings::ele][0]),
															&(vector_EY[Settings::uncorrected][Settings::EB][Settings::ele][0]));
   FR_SS_electron_EB_unc->SetName("FR_SS_electron_EB_unc");
	
   FR_SS_electron_EE_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EE][Settings::ele].size(),
															&(vector_X[Settings::uncorrected][Settings::EE][Settings::ele][0]),
															&(vector_Y[Settings::uncorrected][Settings::EE][Settings::ele][0]),
															&(vector_EX[Settings::uncorrected][Settings::EE][Settings::ele][0]),
															&(vector_EY[Settings::uncorrected][Settings::EE][Settings::ele][0]));
   FR_SS_electron_EE_unc->SetName("FR_SS_electron_EE_unc");
	
   FR_SS_muon_EB_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EB][Settings::mu].size(),
													  &(vector_X[Settings::uncorrected][Settings::EB][Settings::mu][0]),
													  &(vector_Y[Settings::uncorrected][Settings::EB][Settings::mu][0]),
													  &(vector_EX[Settings::uncorrected][Settings::EB][Settings::mu][0]),
													  &(vector_EY[Settings::uncorrected][Settings::EB][Settings::mu][0]));
   FR_SS_muon_EB_unc->SetName("FR_SS_muon_EB_unc");
	
   FR_SS_muon_EE_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EE][Settings::mu].size(),
													  &(vector_X[Settings::uncorrected][Settings::EE][Settings::mu][0]),
													  &(vector_Y[Settings::uncorrected][Settings::EE][Settings::mu][0]),
													  &(vector_EX[Settings::uncorrected][Settings::EE][Settings::mu][0]),
													  &(vector_EY[Settings::uncorrected][Settings::EE][Settings::mu][0]));
   FR_SS_muon_EE_unc->SetName("FR_SS_muon_EE_unc");
	
	
   FR_SS_muon_EB = new TGraphErrors (vector_X[Settings::corrected][Settings::EB][Settings::mu].size(),
                                     &(vector_X[Settings::corrected][Settings::EB][Settings::mu][0]),
                                     &(vector_Y[Settings::corrected][Settings::EB][Settings::mu][0]),
                                     &(vector_EX[Settings::corrected][Settings::EB][Settings::mu][0]),
                                     &(vector_EY[Settings::corrected][Settings::EB][Settings::mu][0]));
   FR_SS_muon_EB->SetName("FR_SS_muon_EB");
   
   FR_SS_muon_EE = new TGraphErrors (vector_X[Settings::corrected][Settings::EE][Settings::mu].size(),
                                     &(vector_X[Settings::corrected][Settings::EE][Settings::mu][0]),
                                     &(vector_Y[Settings::corrected][Settings::EE][Settings::mu][0]),
                                     &(vector_EX[Settings::corrected][Settings::EE][Settings::mu][0]),
                                     &(vector_EY[Settings::corrected][Settings::EE][Settings::mu][0]));
   FR_SS_muon_EE->SetName("FR_SS_muon_EE");
	
   // Electron fake rates must be corrected using average number of missing hits
	if ( input_file_data_name != "DONT_CORRECT" ) CorrectElectronFakeRate(input_file_data_name);
	
   FR_SS_electron_EB = new TGraphErrors (vector_X[Settings::corrected][Settings::EB][Settings::ele].size(),
													  &(vector_X[Settings::corrected][Settings::EB][Settings::ele][0]),
													  &(vector_Y[Settings::corrected][Settings::EB][Settings::ele][0]),
													  &(vector_EX[Settings::corrected][Settings::EB][Settings::ele][0]),
													  &(vector_EY[Settings::corrected][Settings::EB][Settings::ele][0]));
   FR_SS_electron_EB->SetName("FR_SS_electron_EB");
	
   FR_SS_electron_EE = new TGraphErrors (vector_X[Settings::corrected][Settings::EE][Settings::ele].size(),
													  &(vector_X[Settings::corrected][Settings::EE][Settings::ele][0]),
													  &(vector_Y[Settings::corrected][Settings::EE][Settings::ele][0]),
													  &(vector_EX[Settings::corrected][Settings::EE][Settings::ele][0]),
													  &(vector_EY[Settings::corrected][Settings::EE][Settings::ele][0]));
   FR_SS_electron_EE->SetName("FR_SS_electron_EE");

   

   
   PlotFR();
   
   TFile* fOutHistos = new TFile(file_name, "recreate");
   fOutHistos->cd();
   
   FR_SS_electron_EB->Write();
   FR_SS_electron_EE->Write();
   FR_SS_muon_EB->Write();
   FR_SS_muon_EE->Write();
   
   FR_SS_electron_EB_unc->Write();
   FR_SS_electron_EE_unc->Write();
   FR_SS_muon_EB_unc->Write();
   FR_SS_muon_EE_unc->Write();
   
   fOutHistos->Close();
   delete fOutHistos;
   
   cout << "[INFO] Fake rates produced and stored in a file." << endl;
}
//===============================================================


//========================================================================
void SSmethod::CorrectElectronFakeRate( TString input_file_data_name )
{
	TGraphErrors *FR_MissingHits_graph[num_of_eta_bins][99];
	
	Calculate_FR_nMissingHits(input_file_data_name, FR_MissingHits_graph);
	Fit_FRnMH_graphs(FR_MissingHits_graph);
	cout << "[INFO] All graphs fitted." << endl;
	Correct_Final_FR( input_file_data_name );
	cout << "[INFO] Electron fake rates have been corrected." << endl;
}
//========================================================================


//========================================================================
void SSmethod::Calculate_FR_nMissingHits( TString input_file_data_name, TGraphErrors *FR_MissingHits_graph[99][99] )
{
	input_file_data = new TFile( input_file_data_name);
	
	hCounters = (TH1F*)input_file_data->Get("CRZLTree/Counters");
	gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
	
	input_tree_data = (TTree*)input_file_data->Get("CRZLTree/candTree");
	Init( input_tree_data, input_file_data_name , false);
	
	_current_process = find_current_process(input_file_data_name);
	
	for ( int i_pt = 0; i_pt < _n_pT_bins-2; i_pt++)
	{
		for ( int i_eta = 0; i_eta < num_of_eta_bins; i_eta++)
		{
			for ( int i_ZMass = 0; i_ZMass < num_of_z_mass_windows; i_ZMass++ )
			{
				_N_MissingHits[i_ZMass][i_eta][i_pt] = 0.;
				_N_Passing[i_ZMass][i_eta][i_pt] = 0.;
				_N_Failling[i_ZMass][i_eta][i_pt] = 0.;
			}
		}
	}
	
	if (fChain == 0) return;
	
	Long64_t nentries = fChain->GetEntriesFast();
	
	Long64_t nbytes = 0, nb = 0;
	
	for (Long64_t jentry=0; jentry<nentries;jentry++)
	{
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);
		nbytes += nb;
		
		if ( abs(LepLepId->at(2)) != 11 ) continue; // only electrons
		if ( (LepPt->at(0) > LepPt->at(1)) && (LepPt->at(0) < 20. || LepPt->at(1) < 10.) ) continue;
		if ( (LepPt->at(1) > LepPt->at(0)) && (LepPt->at(1) < 20. || LepPt->at(0) < 10.) ) continue;
		if ( LepSIP->at(2) > 4.) continue;
		if ( PFMET > 25. ) continue;
		else
		{
			_current_pT_bin = Find_Ele_pT_bin ( LepPt->at(2) );
			_current_eta_bin = Find_Ele_eta_bin ( LepEta->at(2));
			
			if ( (Z1Mass > 40.) && (Z1Mass < 120.) )
			{
				_N_MissingHits[Settings::_40_MZ1_120][_current_eta_bin][_current_pT_bin] += LepMissingHit->at(2);
				if(LepisID->at(2) && ((fabs(LepLepId->at(2)) == 11) ? LepCombRelIsoPF->at(2) < 999999. : LepCombRelIsoPF->at(2) < 0.35)) _N_Passing[Settings::_40_MZ1_120][_current_eta_bin][_current_pT_bin] += 1.;
				else _N_Failling[Settings::_40_MZ1_120][_current_eta_bin][_current_pT_bin] += 1.;
			}
			
			if ( abs( Z1Mass - 91.2 ) < 7. )
			{
				_N_MissingHits[Settings::_MZ1mMZtrue_7][_current_eta_bin][_current_pT_bin] += LepMissingHit->at(2);
				if(LepisID->at(2) && ((fabs(LepLepId->at(2)) == 11) ? LepCombRelIsoPF->at(2) < 999999. : LepCombRelIsoPF->at(2) < 0.35)) _N_Passing[Settings::_MZ1mMZtrue_7][_current_eta_bin][_current_pT_bin] += 1.;
				else _N_Failling[Settings::_MZ1mMZtrue_7][_current_eta_bin][_current_pT_bin] += 1.;
			}
			
			if ( (Z1Mass > 60.) && (Z1Mass < 120.) )
			{
				_N_MissingHits[Settings::_60_MZ1_120][_current_eta_bin][_current_pT_bin] += LepMissingHit->at(2);
				if(LepisID->at(2) && ((fabs(LepLepId->at(2)) == 11) ? LepCombRelIsoPF->at(2) < 999999. : LepCombRelIsoPF->at(2) < 0.35)) _N_Passing[Settings::_60_MZ1_120][_current_eta_bin][_current_pT_bin] += 1.;
				else _N_Failling[Settings::_60_MZ1_120][_current_eta_bin][_current_pT_bin] += 1.;
			}
			
			TLorentzVector p1,p2,p3;
			p1.SetPtEtaPhiM(LepPt->at(0), LepEta->at(0), LepPhi->at(0), 0.);
			p2.SetPtEtaPhiM(LepPt->at(1), LepEta->at(1), LepPhi->at(1), 0.);
			p3.SetPtEtaPhiM(LepPt->at(2), LepEta->at(2), LepPhi->at(2), 0.);
			
			if ( abs( ((p1+p2)+p3).M() - 91.2 ) < 5. )//3 lepton mass
			{
				_N_MissingHits[Settings::_MZ1EmMZtrue_5][_current_eta_bin][_current_pT_bin] += LepMissingHit->at(2);
				if(LepisID->at(2) && ((fabs(LepLepId->at(2)) == 11) ? LepCombRelIsoPF->at(2) < 999999. : LepCombRelIsoPF->at(2) < 0.35)) _N_Passing[Settings::_MZ1EmMZtrue_5][_current_eta_bin][_current_pT_bin] += 1.;
				else _N_Failling[Settings::_MZ1EmMZtrue_5][_current_eta_bin][_current_pT_bin] += 1.;
			}
		}
	} // END events loop
	
	//Fill vectors to produce TGraphs
	
	vector<Float_t> vector_x[num_of_eta_bins][_n_pT_bins-2];
	vector<Float_t> vector_y[num_of_eta_bins][_n_pT_bins-2];
	vector<Float_t> vector_ex[num_of_eta_bins][_n_pT_bins-2];
	vector<Float_t> vector_ey[num_of_eta_bins][_n_pT_bins-2];
	
	for ( int i_pt = 0; i_pt < _n_pT_bins-2; i_pt++)
	{
		for ( int i_eta = 0; i_eta < num_of_eta_bins; i_eta++)
		{
			for ( int i_ZMass = 0; i_ZMass < num_of_z_mass_windows; i_ZMass++ )
			{
				vector_x[i_eta][i_pt].push_back(_N_MissingHits[i_ZMass][i_eta][i_pt]/(_N_Passing[i_ZMass][i_eta][i_pt] + _N_Failling[i_ZMass][i_eta][i_pt]));
				vector_ex[i_eta][i_pt].push_back(sqrt(pow((1./pow(_N_Failling[i_ZMass][i_eta][i_pt]+_N_Passing[i_ZMass][i_eta][i_pt],1)),2)*_N_MissingHits[i_ZMass][i_eta][i_pt] + pow((_N_MissingHits[i_ZMass][i_eta][i_pt]/pow(_N_Failling[i_ZMass][i_eta][i_pt]+_N_Passing[i_ZMass][i_eta][i_pt],2)),2)*(_N_Failling[i_ZMass][i_eta][i_pt]+_N_Passing[i_ZMass][i_eta][i_pt])));
				
				vector_y[i_eta][i_pt].push_back(_N_Passing[i_ZMass][i_eta][i_pt]/(_N_Passing[i_ZMass][i_eta][i_pt] + _N_Failling[i_ZMass][i_eta][i_pt]));
				vector_ey[i_eta][i_pt].push_back(sqrt(pow((_N_Failling[i_ZMass][i_eta][i_pt]/pow(_N_Failling[i_ZMass][i_eta][i_pt]+_N_Passing[i_ZMass][i_eta][i_pt],2)),2)*_N_Passing[i_ZMass][i_eta][i_pt] + pow((_N_Passing[i_ZMass][i_eta][i_pt]/pow(_N_Failling[i_ZMass][i_eta][i_pt]+_N_Passing[i_ZMass][i_eta][i_pt],2)),2)*_N_Failling[i_ZMass][i_eta][i_pt]));
				
//				cout << "========================================" << endl;
//				cout << "[INFO] Z+L printout." << endl;
//				cout << "========================================" << endl;
//				cout << "Control region number: " << i_ZMass << endl;
//				cout << "pT bin = " << _pT_bins[i_pt+1] << " - " <<  _pT_bins[i_pt + 2] << endl;
//				cout << "eta bin = " << i_eta << endl;
//				cout << "NP = " << _N_Passing[i_ZMass][i_eta][i_pt] << endl;
//				cout << "NF = " << _N_Failling[i_ZMass][i_eta][i_pt] << endl;
//				cout << "MH = " << _N_MissingHits[i_ZMass][i_eta][i_pt] << endl;
//				cout << "avg_MH = " << vector_x[i_eta][i_pt][i_ZMass] << endl;
//				cout << "avg_MH error = " << vector_ex[i_eta][i_pt][i_ZMass] << endl;
//				cout << "FR = " << vector_y[i_eta][i_pt][i_ZMass] << endl;
//				cout << "FR error = " << vector_ey[i_eta][i_pt][i_ZMass] << endl;
			}
		}
		
	}
	
	for ( int i_pt = 0; i_pt < _n_pT_bins-2; i_pt++)
	{
		for ( int i_eta = 0; i_eta < num_of_eta_bins; i_eta++)
		{
			FR_MissingHits_graph[i_eta][i_pt] = new TGraphErrors(vector_x[i_eta][i_pt].size(),
																					&(vector_x[i_eta][i_pt][0]),
																					&(vector_y[i_eta][i_pt][0]),
																					&(vector_ex[i_eta][i_pt][0]),
																					&(vector_ey[i_eta][i_pt][0]));
		}
		
	}
	
}
//========================================================================


//============================================================================
void SSmethod::Fit_FRnMH_graphs(TGraphErrors *FR_MissingHits_graph[99][99])
{
	for ( int i_pt = 0; i_pt < _n_pT_bins-2; i_pt++)
	{
		for ( int i_eta = 0; i_eta < num_of_eta_bins; i_eta++)
		{
			TString func_name;
			func_name.Form("FR_MissingHits_func_eta_%d_pT_%d",i_eta,i_pt);
			Ele_FR_correction_function[i_eta][i_pt] = new TF1(func_name,"[0]*x+[1]",0,3);
			Ele_FR_correction_function[i_eta][i_pt]->SetParameter(0,1.);
			Ele_FR_correction_function[i_eta][i_pt]->SetParameter(1,0.);
			
			FR_MissingHits_graph[i_eta][i_pt]->Fit(Ele_FR_correction_function[i_eta][i_pt], "Q");
			
			TString graph_name;
			graph_name.Form("FR_MissingHits_graph_eta_%d_pT_%d",i_eta,i_pt);
			FR_MissingHits_graph[i_eta][i_pt]->SetName(graph_name);
			FR_MissingHits_graph[i_eta][i_pt]->GetXaxis()->SetTitle("<# Missing Hits>");
			FR_MissingHits_graph[i_eta][i_pt]->GetYaxis()->SetTitle("Fake Rate");
			TCanvas *c1 = new TCanvas(graph_name,graph_name,900,900);
			c1->cd();
			FR_MissingHits_graph[i_eta][i_pt]->Draw("AP");
			system("mkdir -p Fits");
			SavePlots(c1, "Fits/" + graph_name);
		}
	}
}
//============================================================================


//=============================================================
void SSmethod::Correct_Final_FR( TString input_file_data_name)
{
	input_file_data = new TFile( input_file_data_name);
	
	hCounters = (TH1F*)input_file_data->Get("CRZLLTree/Counters");
	gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
	
	input_tree_data = (TTree*)input_file_data->Get("CRZLLTree/candTree");
	Init( input_tree_data, input_file_data_name , true);
	
	_current_process = find_current_process(input_file_data_name);
	
	if (fChain == 0) return;
	
	Long64_t nentries = fChain->GetEntriesFast();
	
	Long64_t nbytes = 0, nb = 0;
	
	float _N_MissingHits_ZLL[num_of_eta_bins][_n_pT_bins];
	float _N_Passing_ZLL[num_of_eta_bins][_n_pT_bins];
	float _N_Failling_ZLL[num_of_eta_bins][_n_pT_bins];
	
	float _avg_MissingHits_ZLL[num_of_eta_bins][_n_pT_bins];
	
	for ( int i_pt = 0; i_pt <= _n_pT_bins-2; i_pt++)
	{
		for ( int i_eta = 0; i_eta < num_of_eta_bins; i_eta++)
		{
			_N_MissingHits_ZLL[i_eta][i_pt] = 0.;
			_N_Passing_ZLL[i_eta][i_pt] = 0.;
			_N_Failling_ZLL[i_eta][i_pt] = 0.;
		}
	}
	
	for (Long64_t jentry=0; jentry<nentries;jentry++)
	{
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);
		nbytes += nb;
		
		if (!(test_bit(CRflag, CRZLLss))) continue;
		
		if ( abs(Z2Flav) != 121) continue; // only electrons
		if ( abs(Z1Flav) != 121) continue; // only 4e
//		if ( (LepPt->at(0) > LepPt->at(1)) && (LepPt->at(0) < 20. || LepPt->at(1) < 10.) ) continue;
//		if ( (LepPt->at(1) > LepPt->at(0)) && (LepPt->at(1) < 20. || LepPt->at(0) < 10.) ) continue;
		else
		{
			_current_pT_bin = Find_Ele_pT_bin ( LepPt->at(2) );
			_current_eta_bin = Find_Ele_eta_bin ( LepEta->at(2));

			_N_MissingHits_ZLL[_current_eta_bin][_current_pT_bin] += LepMissingHit->at(2);
			if(LepisID->at(2) && ((fabs(LepLepId->at(2)) == 11) ? LepCombRelIsoPF->at(2) < 999999. : LepCombRelIsoPF->at(2) < 0.35)) _N_Passing_ZLL[_current_eta_bin][_current_pT_bin] += 1.;
			else _N_Failling_ZLL[_current_eta_bin][_current_pT_bin] += 1.;
			
			_current_pT_bin = Find_Ele_pT_bin ( LepPt->at(3) );
			_current_eta_bin = Find_Ele_eta_bin ( LepEta->at(3));
			
			_N_MissingHits_ZLL[_current_eta_bin][_current_pT_bin] += LepMissingHit->at(3);
			if(LepisID->at(3) && ((fabs(LepLepId->at(3)) == 11) ? LepCombRelIsoPF->at(3) < 999999. : LepCombRelIsoPF->at(3) < 0.35)) _N_Passing_ZLL[_current_eta_bin][_current_pT_bin] += 1.;
			else _N_Failling_ZLL[_current_eta_bin][_current_pT_bin] += 1.;
		}
		
	} // END events loop
	
	for ( int i_pt = 0; i_pt < _n_pT_bins-2; i_pt++)
	{
		Float_t sigma_avgMH = 0;
		_avg_MissingHits_ZLL[Settings::EB][i_pt] = _N_MissingHits_ZLL[Settings::EB][i_pt]/(_N_Passing_ZLL[Settings::EB][i_pt] + _N_Failling_ZLL[Settings::EB][i_pt]);
		sigma_avgMH = sqrt(pow((1./pow(_N_Failling_ZLL[Settings::EB][i_pt]+_N_Passing_ZLL[Settings::EB][i_pt],1)),2)*_N_MissingHits_ZLL[Settings::EB][i_pt] + pow((_N_MissingHits_ZLL[Settings::EB][i_pt]/pow(_N_Failling_ZLL[Settings::EB][i_pt]+_N_Passing_ZLL[Settings::EB][i_pt],2)),2)*(_N_Failling_ZLL[Settings::EB][i_pt]+_N_Passing_ZLL[Settings::EB][i_pt]));
		
		vector_X[Settings::corrected][Settings::EB][Settings::ele][i_pt] = ((_pT_bins[i_pt + 1] + _pT_bins[i_pt + 2])/2);
		vector_Y[Settings::corrected][Settings::EB][Settings::ele][i_pt] = (Ele_FR_correction_function[Settings::EB][i_pt]->Eval(_avg_MissingHits_ZLL[Settings::EB][i_pt]));
		
		vector_EX[Settings::corrected][Settings::EB][Settings::ele][i_pt] = ((_pT_bins[i_pt + 2] - _pT_bins[i_pt + 1])/2);
		vector_EY[Settings::corrected][Settings::EB][Settings::ele][i_pt] = (Ele_FR_correction_function[Settings::EB][i_pt]->Eval(_avg_MissingHits_ZLL[Settings::EB][i_pt]) - Ele_FR_correction_function[Settings::EB][i_pt]->Eval(_avg_MissingHits_ZLL[Settings::EB][i_pt] - sigma_avgMH));
		
		
		_avg_MissingHits_ZLL[Settings::EE][i_pt] = _N_MissingHits_ZLL[Settings::EE][i_pt]/(_N_Passing_ZLL[Settings::EE][i_pt] + _N_Failling_ZLL[Settings::EE][i_pt]);
		sigma_avgMH = sqrt(pow((1./pow(_N_Failling_ZLL[Settings::EE][i_pt]+_N_Passing_ZLL[Settings::EE][i_pt],1)),2)*_N_MissingHits_ZLL[Settings::EE][i_pt] + pow((_N_MissingHits_ZLL[Settings::EE][i_pt]/pow(_N_Failling_ZLL[Settings::EE][i_pt]+_N_Passing_ZLL[Settings::EE][i_pt],2)),2)*(_N_Failling_ZLL[Settings::EE][i_pt]+_N_Passing_ZLL[Settings::EE][i_pt]));
		
		vector_X[Settings::corrected][Settings::EE][Settings::ele][i_pt] = ((_pT_bins[i_pt + 1] + _pT_bins[i_pt + 2])/2);
		vector_Y[Settings::corrected][Settings::EE][Settings::ele][i_pt] = (Ele_FR_correction_function[Settings::EE][i_pt]->Eval(_avg_MissingHits_ZLL[Settings::EE][i_pt]));
		
		vector_EX[Settings::corrected][Settings::EE][Settings::ele][i_pt] = ((_pT_bins[i_pt + 2] - _pT_bins[i_pt + 1])/2);
		vector_EY[Settings::corrected][Settings::EE][Settings::ele][i_pt] = (Ele_FR_correction_function[Settings::EE][i_pt]->Eval(_avg_MissingHits_ZLL[Settings::EE][i_pt]) - Ele_FR_correction_function[Settings::EE][i_pt]->Eval(_avg_MissingHits_ZLL[Settings::EE][i_pt] - sigma_avgMH));
		
//		cout << "========================================" << endl;
//		cout << "[INFO] Z+LL printout." << endl;
//		cout << "========================================" << endl;
//		cout << "pT bin = " << _pT_bins[i_pt + 1] << " - " <<  _pT_bins[i_pt + 2] << endl;
//		cout << "eta bin = " << Settings::EB << endl;
//		cout << "NP = " << _N_Passing_ZLL[Settings::EB][i_pt] << endl;
//		cout << "NF = " << _N_Failling_ZLL[Settings::EB][i_pt] << endl;
//		cout << "avg_MH = " << _avg_MissingHits_ZLL[Settings::EB][i_pt] << endl;
//		cout << "FR = " << _N_Passing_ZLL[Settings::EB][i_pt]/(_N_Passing_ZLL[Settings::EB][i_pt]+_N_Failling_ZLL[Settings::EB][i_pt]) << endl;
//		cout << "corr FR = " << vector_Y[Settings::corrected][Settings::EB][Settings::ele][i_pt] << endl;
//		cout << "========================================" << endl;
//		cout << "pT bin = " << _pT_bins[i_pt + 1] << " - " <<  _pT_bins[i_pt + 2] << endl;
//		cout << "eta bin = " << Settings::EE << endl;
//		cout << "NP = " << _N_Passing_ZLL[Settings::EE][i_pt] << endl;
//		cout << "NF = " << _N_Failling_ZLL[Settings::EE][i_pt] << endl;
//		cout << "avg_MH = " << _avg_MissingHits_ZLL[Settings::EE][i_pt] << endl;
//		cout << "FR = " << _N_Passing_ZLL[Settings::EE][i_pt]/(_N_Passing_ZLL[Settings::EE][i_pt]+_N_Failling_ZLL[Settings::EE][i_pt]) << endl;
//		cout << "corr FR = " << vector_Y[Settings::corrected][Settings::EE][Settings::ele][i_pt] << endl;
	}
	
}
//=============================================================


//===============================================================
void SSmethod::SubtractWZ()
{
   passing[Settings::Total][Settings::mu]->Add(passing[Settings::WZ][Settings::mu], -1.);
   failing[Settings::Total][Settings::mu]->Add(failing[Settings::WZ][Settings::mu], -1.);
   
   cout << "[INFO] WZ contribution subtracted." << endl;
   
}
//===============================================================

//===============================================================
void SSmethod::PlotFR()
{
   TCanvas *c_ele, *c_mu;
   c_ele = new TCanvas("FR_ele", "FR_ele", 600, 600);
   c_mu  = new TCanvas("FR_mu", "FR_mu", 600, 600);
   
   mg_electrons = new TMultiGraph();
   mg_muons = new TMultiGraph();
   
   mg_electrons->Add(FR_SS_electron_EB);
   FR_SS_electron_EB->SetLineColor(kBlue);
   FR_SS_electron_EB->SetLineStyle(2);
   FR_SS_electron_EB->SetMarkerSize(0);
   FR_SS_electron_EB->SetTitle("barel corrected");
   mg_electrons->Add(FR_SS_electron_EE);
   FR_SS_electron_EE->SetLineColor(kRed);
   FR_SS_electron_EE->SetLineStyle(2);
   FR_SS_electron_EE->SetMarkerSize(0);
   FR_SS_electron_EE->SetTitle("endcap corrected");
   mg_electrons->Add(FR_SS_electron_EB_unc);
   FR_SS_electron_EB_unc->SetLineColor(kBlue);
   FR_SS_electron_EB_unc->SetLineStyle(1);
   FR_SS_electron_EB_unc->SetMarkerSize(0);
   FR_SS_electron_EB_unc->SetTitle("barel uncorrected");
   mg_electrons->Add(FR_SS_electron_EE_unc);
   FR_SS_electron_EE_unc->SetLineColor(kRed);
   FR_SS_electron_EE_unc->SetLineStyle(1);
   FR_SS_electron_EE_unc->SetMarkerSize(0);
   FR_SS_electron_EE_unc->SetTitle("endcap uncorrected");
   
   mg_muons->Add(FR_SS_muon_EB);
   FR_SS_muon_EB->SetLineColor(kBlue);
   FR_SS_muon_EB->SetLineStyle(2);
   FR_SS_muon_EB->SetMarkerSize(0);
   FR_SS_muon_EB->SetTitle("barel corrected");
   mg_muons->Add(FR_SS_muon_EE);
   FR_SS_muon_EE->SetLineColor(kRed);
   FR_SS_muon_EE->SetLineStyle(2);
   FR_SS_muon_EE->SetMarkerSize(0);
   FR_SS_muon_EE->SetTitle("endcap corrected");
   mg_muons->Add(FR_SS_muon_EB_unc);
   FR_SS_muon_EB_unc->SetLineColor(kBlue);
   FR_SS_muon_EB_unc->SetLineStyle(1);
   FR_SS_muon_EB_unc->SetMarkerSize(0);
   FR_SS_muon_EB_unc->SetTitle("barel uncorrected");
   mg_muons->Add(FR_SS_muon_EE_unc);
   FR_SS_muon_EE_unc->SetLineColor(kRed);
   FR_SS_muon_EE_unc->SetLineStyle(1);
   FR_SS_muon_EE_unc->SetMarkerSize(0);
   FR_SS_muon_EE_unc->SetTitle("endcap uncorrected");

   gStyle->SetEndErrorSize(0);
   
   TLegend *leg_ele,*leg_mu;
   
   c_ele->cd();
   mg_electrons->Draw("AP");
	mg_electrons->GetXaxis()->SetTitle("p_{T} [GeV]");
	mg_electrons->GetYaxis()->SetTitle("Fake Rate");
	mg_electrons->SetTitle("Electron fake rate");
   mg_electrons->SetMaximum(0.35);
   leg_ele = CreateLegend_FR("left",FR_SS_electron_EB_unc,FR_SS_electron_EB,FR_SS_electron_EE_unc,FR_SS_electron_EE);
   leg_ele->Draw();
   SavePlots(c_ele, "Plots/FR_SS_electrons");
   
   c_mu->cd();
   mg_muons->Draw("AP");
	mg_muons->GetXaxis()->SetTitle("p_{T} [GeV]");
	mg_muons->GetYaxis()->SetTitle("Fake Rate");
	mg_muons->SetTitle("Muon fake rate");
   mg_muons->SetMaximum(0.35);
   leg_mu = CreateLegend_FR("left",FR_SS_muon_EB_unc,FR_SS_muon_EB,FR_SS_muon_EE_unc,FR_SS_muon_EE);
   leg_mu->Draw();
   SavePlots(c_mu, "Plots/FR_SS_muons");
   
}
//===============================================================

//========================================================================================================
void SSmethod::PlotDataMC( TString variable_name, TString folder )
{
   TCanvas *c;
   c = new TCanvas("ZLLss", variable_name, 600, 600);
   
   if ( GetVarLogX( variable_name) ) c->SetLogx();
   if ( GetVarLogY( variable_name) ) c->SetLogy();
   
   for( int i_fs = 0; i_fs <= Settings::fs4l ; i_fs++ )
   {
      for ( int i_cat = 0; i_cat <= Settings::inclusive; i_cat++ )
      {
         histos_1D[Settings::regZLL][Settings::WZ][i_fs][i_cat]   ->SetFillColor(kMagenta);
         histos_1D[Settings::regZLL][Settings::qqZZ][i_fs][i_cat] ->SetFillColor(kCyan+1);
         histos_1D[Settings::regZLL][Settings::DY][i_fs][i_cat]   ->SetFillColor(kGreen-1);
         histos_1D[Settings::regZLL][Settings::ttbar][i_fs][i_cat]->SetFillColor(kBlue-2);
         
         histos_1D[Settings::regZLL][Settings::WZ][i_fs][i_cat]   ->SetLineColor(kMagenta);
         histos_1D[Settings::regZLL][Settings::qqZZ][i_fs][i_cat] ->SetLineColor(kCyan+1);
         histos_1D[Settings::regZLL][Settings::DY][i_fs][i_cat]   ->SetLineColor(kGreen-1);
         histos_1D[Settings::regZLL][Settings::ttbar][i_fs][i_cat]->SetLineColor(kBlue-2);
         
         histos_1D[Settings::regZLL][Settings::Data][i_fs][i_cat]->SetMarkerSize(0.8);
         histos_1D[Settings::regZLL][Settings::Data][i_fs][i_cat]->SetMarkerStyle(20);
         histos_1D[Settings::regZLL][Settings::Data][i_fs][i_cat]->SetBinErrorOption(TH1::kPoisson);
         histos_1D[Settings::regZLL][Settings::Data][i_fs][i_cat]->SetLineColor(kBlack);
         
         THStack *stack = new THStack( "stack", "stack" );
			stack->Add(histos_1D[Settings::regZLL][Settings::qqZZ][i_fs][i_cat]);
			stack->Add(histos_1D[Settings::regZLL][Settings::WZ][i_fs][i_cat]);
			stack->Add(histos_1D[Settings::regZLL][Settings::ttbar][i_fs][i_cat]);
         stack->Add(histos_1D[Settings::regZLL][Settings::DY][i_fs][i_cat]);
			
         
         stack->Draw("HIST");
         
         float data_max = histos_1D[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetBinContent(histos_1D[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetMaximumBin());
         float data_max_error = histos_1D[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetBinErrorUp(histos_1D[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetMaximumBin());
         
         stack->SetMinimum(1e-5);
         stack->SetMaximum((data_max + data_max_error)*1.1);
         
         TString _fs_label;
			if ( i_fs == Settings::fs4l) _fs_label = "m_{4#font[12]{l}} (GeV)";
         if ( i_fs == Settings::fs4e) _fs_label = "m_{4#font[12]{e}} (GeV)";
         if ( i_fs == Settings::fs4mu) _fs_label = "m_{4#font[12]{#mu}} (GeV)";
         if ( i_fs == Settings::fs2e2mu) _fs_label = "m_{2#font[12]{e}2#font[12]{#mu}} (GeV)";
         if ( i_fs == Settings::fs2mu2e) _fs_label = "m_{2#font[12]{#mu}2#font[12]{e}} (GeV)";
         stack->GetXaxis()->SetTitle(_fs_label);
         stack->GetXaxis()->SetTitleSize(0.04);
         stack->GetXaxis()->SetLabelSize(0.04);
         stack->GetYaxis()->SetTitle(histos_1D[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetYaxis()->GetTitle());
         stack->GetYaxis()->SetTitleSize(0.04);
         stack->GetYaxis()->SetLabelSize(0.04);
         
         stack->GetXaxis()->SetTitleOffset(1.2);
         stack->GetYaxis()->SetTitleOffset(1.25);
         
         histos_1D[Settings::regZLL][Settings::Data][i_fs][i_cat]->Draw("SAME p E1 X0");
         
         TLegend *legend;
         legend  = CreateLegend_ZLL("right",histos_1D[Settings::regZLL][Settings::Data][i_fs][i_cat],histos_1D[Settings::regZLL][Settings::WZ][i_fs][i_cat],histos_1D[Settings::regZLL][Settings::qqZZ][i_fs][i_cat],histos_1D[Settings::regZLL][Settings::DY][i_fs][i_cat],histos_1D[Settings::regZLL][Settings::ttbar][i_fs][i_cat]);
         legend->Draw();
         
         // Draw lumi
         CMS_lumi *lumi = new CMS_lumi;
         lumi->set_lumi(c, _lumi, 0);
         
         TString _out_file_name;
         _out_file_name = folder + "/" + variable_name + "_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
         SavePlots(c, _out_file_name);
         
      }
      
   }
}
//========================================================================================================


//========================================================================================================
void SSmethod::PlotZX( TString variable_name, TString folder )
{
   TCanvas *c;
   c = new TCanvas("ZLLss", variable_name, 600, 600);
   
   if ( GetVarLogX( variable_name) ) c->SetLogx();
   if ( GetVarLogY( variable_name) ) c->SetLogy();
   
   for( int i_fs = 0; i_fs <= Settings::fs4l ; i_fs++ )
   {
      for ( int i_cat = 0; i_cat <= Settings::inclusive; i_cat++ )
      {
         histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->SetFillColor(kGreen-1);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->SetLineColor(kGreen-1);
         
         histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->Draw("HIST");
         
         TString _fs_label;
         if ( i_fs == Settings::fs4e)    _fs_label = "m_{4#font[12]{e}} (GeV)";
         if ( i_fs == Settings::fs4mu)   _fs_label = "m_{4#font[12]{#mu}} (GeV)";
         if ( i_fs == Settings::fs2e2mu) _fs_label = "m_{2#font[12]{e}2#font[12]{#mu}} (GeV)";
         if ( i_fs == Settings::fs2mu2e) _fs_label = "m_{2#font[12]{#mu}2#font[12]{e}} (GeV)";
			if ( i_fs == Settings::fs4l)    _fs_label = "m_{4#font[12]{l}} (GeV)";
         histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetXaxis()->SetTitle(_fs_label);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetXaxis()->SetTitleSize(0.04);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetXaxis()->SetLabelSize(0.04);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetYaxis()->SetTitle(histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetYaxis()->GetTitle());
         histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetYaxis()->SetTitleSize(0.04);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetYaxis()->SetLabelSize(0.04);
         
         histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetXaxis()->SetTitleOffset(1.2);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetYaxis()->SetTitleOffset(1.25);
         
         // Draw lumi
         CMS_lumi *lumi = new CMS_lumi;
         lumi->set_lumi(c, _lumi, 0);
         
         TString _out_file_name;
         _out_file_name = folder + "/" + variable_name + "_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
         SavePlots(c, _out_file_name);
         
      }
      
   }
}
//========================================================================================================


//========================================================================================================
void SSmethod::FitZX( TString variable_name, TString folder )
{
   TCanvas *c;
   TF1  *fit_function;
   CMS_lumi *lumi = new CMS_lumi;
	
	c = new TCanvas("Fits_ZLLss", variable_name, 600, 600);
	
   if ( GetVarLogX( variable_name) ) c->SetLogx();
   if ( GetVarLogY( variable_name) ) c->SetLogy();
	
   for( int i_fs = 0; i_fs <= Settings::fs4l ; i_fs++ )
   {
      for ( int i_cat = 0; i_cat <= Settings::inclusive; i_cat++ )
      {
         TString _fs_label;
         if ( i_fs == Settings::fs4e)    _fs_label = "m_{4#font[12]{e}} (GeV)";
         if ( i_fs == Settings::fs4mu)   _fs_label = "m_{4#font[12]{#mu}} (GeV)";
         if ( i_fs == Settings::fs2e2mu) _fs_label = "m_{2#font[12]{e}2#font[12]{#mu}} (GeV)";
         if ( i_fs == Settings::fs2mu2e) _fs_label = "m_{2#font[12]{#mu}2#font[12]{e}} (GeV)";
         if ( i_fs == Settings::fs4l)    _fs_label = "m_{4#font[12]{l}} (GeV)";
         histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetXaxis()->SetTitle(_fs_label);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetXaxis()->SetTitleSize(0.04);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetXaxis()->SetLabelSize(0.04);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetYaxis()->SetTitle(histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetYaxis()->GetTitle());
         histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetYaxis()->SetTitleSize(0.04);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetYaxis()->SetLabelSize(0.04);
			
         histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetXaxis()->SetTitleOffset(1.2);
         histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->GetYaxis()->SetTitleOffset(1.25);
			
			gStyle->SetOptFit();
			gStyle->SetStatY(0.9);
			gStyle->SetStatX(0.95);
			gStyle->SetStatW(0.2);
			gStyle->SetStatH(0.1);
			
			fit_function = new TF1("fit_function","[0]*TMath::Landau(x, [1], [2])",70,1000);
			fit_function->SetParNames("Constant","MPV","#sigma");
			fit_function->SetParameter(0,1.);
			fit_function->SetParameter(1,100.);
			fit_function->SetParameter(2,10.);

			histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->Fit("fit_function");
			histos_ZX[Settings::regZLL][Settings::Data][i_fs][i_cat]->Draw("");
			
         // Draw lumi
         lumi->set_lumi(c, _lumi, 0);
			
         TString _out_file_name;
         _out_file_name = folder + "/" + variable_name + "_ZX_SS_fit_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
         SavePlots(c, _out_file_name);
      }
	}
	gStyle->SetOptFit(0);
}
//========================================================================================================



//===============================================================
void SSmethod::RemoveNegativeBins1D(TH1F *h)
{
   for (int i_bin_x = 1; i_bin_x <= h->GetXaxis()->GetNbins(); i_bin_x++)
   {
      if( h->GetBinContent(i_bin_x) < 0.) h->SetBinContent(i_bin_x, 0);
   }
   
}
//===============================================================

//===============================================================
void SSmethod::RemoveNegativeBins2D(TH2F *h)
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
void SSmethod::Set_pT_binning(int size, float *bins)
{
   _n_pT_bins = size;

   for (int i = 0; i < size; i++)
   {
      _pT_bins[i] = bins[i];
   }
}
//===============================================================

//===============================================================
void SSmethod::SetLumi(float lumi)
{
   _lumi = lumi;
}
//===============================================================


//==========================================================
int SSmethod::find_current_process( TString input_file_name )
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
int SSmethod::FindFinalState()
{
   int final_state = -999;

   if ( Z1Flav == -121 )
   {
      if ( abs(Z2Flav) == 121 )
         final_state = Settings::fs4e;
      else if ( abs(Z2Flav) == 169 )
         final_state = Settings::fs2e2mu;
      else
         cerr << "[ERROR] in event " << RunNumber << ":" << LumiNumber << ":" << EventNumber << ", Z2Flav = " << Z2Flav << endl;
      }
   else if ( Z1Flav == -169 )
   {
      if ( abs(Z2Flav) == 121 )
         final_state = Settings::fs2mu2e;
      else if ( abs(Z2Flav) == 169 )
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

//=========================================
int SSmethod::Find_Ele_pT_bin( Float_t pT )
{
	int bin = 0;
	
	for ( int i = 2; i <= _n_pT_bins ; i++)
	{
		if ( (pT > _pT_bins[i-1]) && (pT < _pT_bins[i]) ) bin = i - 2;
	}
	if ( pT > 80. ) bin = _n_pT_bins - 2;

	//cout << "PT = " << pT << " bin = " << bin << endl;
	return bin;
}
//=========================================

//=========================================
int SSmethod::Find_Ele_eta_bin( Float_t eta )
{
	int bin = 0;
	
	if (abs(eta) < 1.479) bin = Settings::EB;
	else bin = Settings::EE;
	
	//cout << "eta = " << eta << " bin = " << bin << endl;
	return bin;
}
//=========================================


//=================================
float SSmethod::calculate_K_factor(TString input_file_name)
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
bool SSmethod::GetVarLogX ( TString variable_name )
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
bool SSmethod::GetVarLogY ( TString variable_name )
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
TLegend* SSmethod::CreateLegend_FR( string position, TGraphErrors *EB_unc, TGraphErrors *EB_cor,TGraphErrors *EE_unc,TGraphErrors *EE_cor )
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
TLegend* SSmethod::CreateLegend_ZLL( string position, TH1F *data, TH1F *WZ,TH1F *qqZZ,TH1F *DY,TH1F *ttbar )
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

//=======================================
void SSmethod::SavePlots( TCanvas *c, TString name)
{
   c->SaveAs(name + ".pdf");
   c->SaveAs(name + ".root");
   c->SaveAs(name + ".eps");
   gSystem->Exec("convert -density 300 -quality 100 " + name + ".eps " + name + ".png");
}
//=======================================





