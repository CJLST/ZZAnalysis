// Include classes
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/include/Yields.h>

// Constructor
//=======================
Yields::Yields( double lumi ):Tree()
{
   yields_histos = new Histograms( lumi );
   histo_map["Yields"] = yields_histos;

   _lumi = lumi;
   _merge_2e2mu = true;
   _current_process = -999;
   _k_factor = 1;
   _current_final_state = -999;
   _current_category = -999;
   
   //Get colours needed for plots
   _tclr->GetColor("#000000");
   _tclr->GetColor("#ffafaf");
   _tclr->GetColor("#ff9090");
   _tclr->GetColor("#99ccff");
   _tclr->GetColor("#8bc5ff");
   _tclr->GetColor("#3366ff");
   _tclr->GetColor("#2c5Df0");
   _tclr->GetColor("#669966");
   _tclr->GetColor("#6dae6d");
   _tclr->GetColor("#9a6666");
   _tclr->GetColor("#cc0000");
   _tclr->GetColor("#000099");
   _tclr->GetColor("#003300");
   _tclr->GetColor("#5f3f3f");
   
   // Z+X SS factors
   /*****************************************
 *        4e,    4mu,       2e2mu,    2mu2e
 * 2016: 1.00245, 0.998863, 1.03338, 0.998852
 *
 * 2017: 1.01198, 1.03949, 1.013128, 1.00257
 *
 * 2018: 1.00568, 1.02926, 1.03226, 1.00432
 * *******************************************/
   _fs_ROS_SS.push_back(1.00245);//4e
   _fs_ROS_SS.push_back(0.998863);//4mu
   _fs_ROS_SS.push_back(1.03338);//2e2mu
   _fs_ROS_SS.push_back(0.998852);//2mu2e

// Place holder. At some point we might want to 
// have a more general constructor that takes the
// year as input.
/*
   if(year == 2016) {
       _fs_ROS_SS.push_back(1.00245);//4e
       _fs_ROS_SS.push_back(0.998863);//4mu
       _fs_ROS_SS.push_back(1.03338);//2e2mu
       _fs_ROS_SS.push_back(0.998852);//2mu2e
   }
   else if (year == 2017) {
       _fs_ROS_SS.push_back(1.01198);//4e
       _fs_ROS_SS.push_back(1.03949);//4mu
       _fs_ROS_SS.push_back(1.013128);//2e2mu
       _fs_ROS_SS.push_back(1.00257);//2mu2e
   }
   else {
       _fs_ROS_SS.push_back(1.00568);//4e
       _fs_ROS_SS.push_back(1.02926);//4mu
       _fs_ROS_SS.push_back(1.03226);//2e2mu
       _fs_ROS_SS.push_back(1.00432);//2mu2e
   }
*/   
   vector<float> temp;
   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
      {
         temp.push_back(0);
      }
      _expected_yield_SR.push_back(temp);
      _number_of_events_CR.push_back(temp);
   }
}
//============================================================



// Destructor
//====================
Yields::~Yields()
{
}
//====================


//===================================================================
void Yields::Split_2e2mu()
{
   _merge_2e2mu = false;
   yields_histos->Split_2e2mu();
}
//===================================================================


//=====================================================
void Yields::MakeHistograms( TString input_file_name , int year)
{
   
   input_file = new TFile(input_file_name);
   
   hCounters = (TH1F*)input_file->Get("ZZTree/Counters");
   n_gen_events = (Long64_t)hCounters->GetBinContent(1);
   gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
   
   input_tree = (TTree*)input_file->Get("ZZTree/candTree");
   Init( input_tree, input_file_name );
   
   if(year == 2018) _lumi = 59.7;
   if(year == 2017) _lumi = 41.5;
   if(year == 2016) _lumi = 35.9;
   
   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      // Check number of leptons in event
      if ( LepEta->size() != 4 )
      {
         cout << "[ERROR] in event " << RunNumber << ":" << LumiNumber << ":" << EventNumber << ", stored " << LepEta->size() << " leptons instead of 4" << endl;
         continue;
      }
      
      if ( !(ZZsel >= 90) ) continue;
      if ( LepSIP->at(0) > 8. || LepSIP->at(1) > 8. || LepSIP->at(2) > 8. || LepSIP->at(3) > 8.) continue;//SIP Sanity check
      // Find current process
      gen_assoc_lep_id_.push_back(GenAssocLep1Id);
      gen_assoc_lep_id_.push_back(GenAssocLep2Id);
      _n_gen_assoc_lep = CountAssociatedLeptons();
      _current_process = find_current_process( input_file_name, genExtInfo , _n_gen_assoc_lep);
      gen_assoc_lep_id_.clear();
      
      // Final states
      _current_final_state = FindFinalState();
      
      // Find current category
      for ( int j = 0; j < nCleanedJetsPt30; j++)
      {
         jetPt[j] = JetPt->at(j);
         jetEta[j] = JetEta->at(j);
         jetPhi[j] = JetPhi->at(j);
         jetMass[j] = JetMass->at(j);
         jetQGL[j] = JetQGLikelihood->at(j);
         jetPgOverPq[j] = 1./JetQGLikelihood->at(j) - 1.;
      }
      
      _current_category = categoryMor18(nExtraLep,
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
      
      // K factors
      _k_factor = calculate_K_factor(input_file_name);

      // Final event weight
      _event_weight = (_lumi * 1000 * xsec * _k_factor * overallEventWeight * L1prefiringWeight) / gen_sum_weights;
      if ( input_file_name.Contains("ggH") ) _event_weight *= ggH_NNLOPS_weight; // reweight POWHEG ggH to NNLOPS
   
      // Fill M4l histograms
       yields_histos->FillYields( ZZMass, _event_weight, _current_final_state, _current_category, _current_process );
   } // end for loop
   
   cout << "[INFO] Histograms for " << input_file_name << " filled." << endl;
}
//=====================================================



//===============================================================================
void Yields::Calculate_SS_ZX_Yields( TString input_file_data_name, TString  input_file_FR_name )
{
   
   FakeRates *FR = new FakeRates( input_file_FR_name );
   
   input_file_data = new TFile(input_file_data_name);
   input_tree_data = (TTree*)input_file_data->Get("CRZLLTree/candTree");
   Init( input_tree_data, input_file_data_name );
   
   
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
      
      _current_final_state = FindFinalStateZX();
      
      _current_category = categoryMor18(nExtraLep,
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
      
      
      // Calculate yield
      _yield_SR = _fs_ROS_SS.at(_current_final_state)*FR->GetFakeRate(LepPt->at(2),LepEta->at(2),LepLepId->at(2))*FR->GetFakeRate(LepPt->at(3),LepEta->at(3),LepLepId->at(3));
      
      _expected_yield_SR[_current_final_state][_current_category] += _yield_SR; // this number needs to be used when renormalizing histograms that have some cut/blinding
      _number_of_events_CR[_current_final_state][_current_category]++;
      
   } // END events loop
   
   for (  int i_cat = 0; i_cat < num_of_categories - 1; i_cat++  )
   {
      for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++  )
      {
         _expected_yield_SR[Settings::fs4l][i_cat]       += _expected_yield_SR[i_fs][i_cat];   //calculate expected yield for inclusive 4l final state
         _number_of_events_CR[Settings::fs4l][i_cat]     += _number_of_events_CR[i_fs][i_cat];
         _expected_yield_SR[i_fs][Settings::inclusive]   += _expected_yield_SR[i_fs][i_cat];   //calculate expected yield for inclusive category
         _number_of_events_CR[i_fs][Settings::inclusive] += _number_of_events_CR[i_fs][i_cat];
         
         if ( _merge_2e2mu )
         {
            _expected_yield_SR[Settings::fs2e2mu][i_cat]       += _expected_yield_SR[Settings::fs2mu2e][i_cat];   //merge 2e2mu and 2mu2e final state
            _number_of_events_CR[Settings::fs2e2mu][i_cat]     += _number_of_events_CR[Settings::fs2mu2e][i_cat];
            _expected_yield_SR[Settings::fs2mu2e][i_cat]        = 0.;
            _number_of_events_CR[Settings::fs2mu2e][i_cat]      = 0.;
         }
         
      }
   }
   for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++  )
   {
      _expected_yield_SR[Settings::fs4l][Settings::inclusive] += _expected_yield_SR[i_fs][Settings::inclusive];
   }
   
   // Print Z + X expected yields for inclusive category
   cout << endl;
   cout << "========================================================================================" << endl;
   cout << "[INFO] Control printout." << endl << "!!! Numbers shoud be identical to yields from SS method !!!" << endl;
   for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++ )
   {
      if ( _merge_2e2mu && i_fs == Settings::fs2mu2e) continue;
      cout << "Category: " << Settings::inclusive << "   Final state: " << i_fs << endl;
      cout << _expected_yield_SR[i_fs][Settings::inclusive] << " +/- " <<
      _expected_yield_SR[i_fs][Settings::inclusive]/sqrt(_number_of_events_CR[i_fs][Settings::inclusive]) << " (stat., evt: " <<
      _number_of_events_CR[i_fs][Settings::inclusive] << ")" << " +/- " << _expected_yield_SR[i_fs][Settings::inclusive]*0.50 << " (syst.)" << endl;
   }
   
   cout << "[INFO] Total = " << _expected_yield_SR[Settings::fs4l][Settings::inclusive] << endl;
   cout << "========================================================================================" << endl;
   cout << endl;
   
   cout << "[INFO] Z+X yields calculated using SS method." << endl;
}
//===============================================================================


//=====================================================
void Yields::ProduceDataROOTFiles( TString input_file_name , TString output_folder_name )
{

   input_file = new TFile(input_file_name);

   hCounters = (TH1F*)input_file->Get("ZZTree/Counters");
   
   input_tree = (TTree*)input_file->Get("ZZTree/candTree");
   Init( input_tree, input_file_name );
   
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   
   _s_final_state.push_back("4e");
   _s_final_state.push_back("4mu");
   _s_final_state.push_back("2e2mu");
   
   _s_category.push_back("UntaggedMor18");
   _s_category.push_back("VBF1JetTaggedMor18");
   _s_category.push_back("VBF2JetTaggedMor18");
   _s_category.push_back("VHLeptTaggedMor18");
   _s_category.push_back("VHHadrTaggedMor18");
   _s_category.push_back("ttHLeptTaggedMor18");
   _s_category.push_back("ttHHadrTaggedMor18");
   
   TString file_name;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      if ( LepEta->size() != 4 )
      {
         cout << "[ERROR] in event " << RunNumber << ":" << LumiNumber << ":" << EventNumber << ", stored " << LepEta->size() << " leptons instead of 4" << endl;
         continue;
      }

      if ( !(ZZsel >= 90) ) continue;
      
      // Final states
      _current_final_state = FindFinalState();
      
      // Categories
      for ( int j = 0; j < nCleanedJetsPt30; j++)
      {
         jetPt[j] = JetPt->at(j);
         jetEta[j] = JetEta->at(j);
         jetPhi[j] = JetPhi->at(j);
         jetMass[j] = JetMass->at(j);
         jetQGL[j] = JetQGLikelihood->at(j);
         jetPgOverPq[j] = 1./JetQGLikelihood->at(j) - 1.;
      }

      _current_category = categoryMor18(nExtraLep,
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
      
      // Calculate kinematic discriminants
      if(_current_category == Settings::VBF_2j_tagged)
      {
          _kd[_current_final_state][_current_category].push_back(D_bkg_VBFdec(  p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal,
                              p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal,
                              p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal,
                              p_JJVBF_BKG_MCFM_JECNominal,
                              p_HadZH_BKG_MCFM_JECNominal,
                              p_HadWH_BKG_MCFM_JECNominal,
                              p_JJQCD_BKG_MCFM_JECNominal,
                              p_HadZH_mavjj_JECNominal,
                              p_HadZH_mavjj_true_JECNominal,
                              p_HadWH_mavjj_JECNominal,
                              p_HadWH_mavjj_true_JECNominal,
                              pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal,
                              pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal,
                              pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal,
                              pConst_JJVBF_BKG_MCFM_JECNominal,
                              pConst_HadZH_BKG_MCFM_JECNominal,
                              pConst_HadWH_BKG_MCFM_JECNominal,
                              pConst_JJQCD_BKG_MCFM_JECNominal,
                              Z1Flav*Z2Flav,
                              ZZMass));
      }
      
      else if(_current_category == Settings::VH_hadron_tagged)
      {
          _kd[_current_final_state][_current_category].push_back(D_bkg_VHdec(   p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal,
                              p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal,
                              p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal,
                              p_JJVBF_BKG_MCFM_JECNominal,
                              p_HadZH_BKG_MCFM_JECNominal,
                              p_HadWH_BKG_MCFM_JECNominal,
                              p_JJQCD_BKG_MCFM_JECNominal,
                              p_HadZH_mavjj_JECNominal,
                              p_HadZH_mavjj_true_JECNominal,
                              p_HadWH_mavjj_JECNominal,
                              p_HadWH_mavjj_true_JECNominal,
                              pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal,
                              pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal,
                              pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal,
                              pConst_JJVBF_BKG_MCFM_JECNominal,
                              pConst_HadZH_BKG_MCFM_JECNominal,
                              pConst_HadWH_BKG_MCFM_JECNominal,
                              pConst_JJQCD_BKG_MCFM_JECNominal,
                              Z1Flav*Z2Flav,
                              ZZMass));
      }
      else _kd[_current_final_state][_current_category].push_back(p_GG_SIG_ghg2_1_ghz1_1_JHUGen / ( p_GG_SIG_ghg2_1_ghz1_1_JHUGen + p_QQB_BKG_MCFM*getDbkgkinConstant(Z1Flav*Z2Flav,ZZMass) ));
      
      _mass[_current_final_state][_current_category].push_back(ZZMass);
      
   } // end for loop
   
   system("mkdir -p " + output_folder_name);
   
   for ( int i_fs = 0; i_fs < num_of_final_states - 2; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories - 2; i_cat++ )
      {
         file_name = output_folder_name + "/data_obs_" + _s_category.at(i_cat) + "_" + _s_final_state.at(i_fs) + ".root";
         data_root_file = new TFile(file_name, "recreate");
         data_obs = new TTree("data_obs","data_obs");
         data_obs->Branch("mass4l",&mass4l,"mass4l/D");
         if( i_cat == Settings::VBF_2j_tagged ) data_obs->Branch("kdCoarse",&kd,"kdCoarse/D");
         else if( i_cat == Settings::VH_hadron_tagged ) data_obs->Branch("kdCoarse",&kd,"kdCoarse/D");
         else data_obs->Branch("kd",&kd,"kd/D");
         
         for (unsigned int i_event = 0; i_event < _mass[i_fs][i_cat].size(); i_event++ )
         {
            mass4l = _mass[i_fs][i_cat].at(i_event);
            kd = _kd[i_fs][i_cat].at(i_event);
            data_obs->Fill();
         }
         
         data_obs->Write();
         data_root_file->Close();
      }
   }
   
   cout << "[INFO] Histograms for " << input_file_name << " filled." << endl;
}
//=====================================================



//=========================================
void Yields::GetHistos( TString file_name )
{
   histo_map[file_name]->GetYieldsHistos( "ROOT_files/" + file_name + ".root" );
   
   cout << "[INFO] Got all histograms." << endl;
}
//=========================================



//===========================
void Yields::FillInclusive()
{
   yields_histos->FillInclusiveYields();
   
   cout << "[INFO] Summing of histograms finished." << endl;
}
//===========================



//==================
void Yields::Save()
{
   system("mkdir -p ROOT_files");
   yields_histos->SaveYieldHistos("ROOT_files/Yields.root");
   
   cout << "[INFO] All yield histograms are saved in a root file." << endl;
}
//==================



//==================
void Yields::Delete()
{
   yields_histos->DeleteYieldsHistos();
   
   cout << "[INFO] Memory clean-up done." << endl;
}
//==================

//==================
void Yields::FillGraphs( TString file_name, float M4l_down, float M4l_up , TString option )
{
   histo_map[file_name]->FillYieldGraphs( M4l_down, M4l_up , option);
   
   cout << "[INFO] Yield vs mH graphs are filled." << endl;   
}
//==================


//==================
void Yields::PrepareYamlFiles( TString file_name, TString sqrt, float M4l_down, float M4l_up )
{
   histo_map[file_name]->PrepareYamlFiles( sqrt, M4l_down, M4l_up, _expected_yield_SR );
   
   cout << "[INFO] Prepared Yaml files." << endl;   
}
//==================



//==================
void Yields::Print( TString file_name )
{
   histo_map[file_name]->PrintYields( _expected_yield_SR );
   
}
//==================



//==================
void Yields::Print( TString file_name, float M4l_down, float M4l_up  )
{
   histo_map[file_name]->PrintYields( M4l_down, M4l_up, _expected_yield_SR);
   
}
//==================

//==================
void Yields::PrintLatexTables( TString file_name, float M4l_down, float M4l_up  )
{
   histo_map[file_name]->PrintLatexTables( M4l_down, M4l_up, _expected_yield_SR);
   
}
//==================


//==========================================================
int Yields::find_current_process( TString input_file_name , int genExtInfo, int n_gen_assoc_lep)
{

   int current_process = -999;
   
   // Assign dataset to correct process
   if ( input_file_name.Contains("Data") )           current_process = Settings::yData;
   
   if ( input_file_name.Contains("ggH120") )         current_process = Settings::yH120ggH;
   if ( input_file_name.Contains("ggH124") )         current_process = Settings::yH124ggH;
   if ( input_file_name.Contains("ggH125") )         current_process = Settings::yH125ggH;
   if ( input_file_name.Contains("ggH126") )         current_process = Settings::yH126ggH;
   if ( input_file_name.Contains("ggH130") )         current_process = Settings::yH130ggH;
   
   if ( input_file_name.Contains("VBFH120") )        current_process = Settings::yH120VBF;
   if ( input_file_name.Contains("VBFH124") )        current_process = Settings::yH124VBF;
   if ( input_file_name.Contains("VBFH125") )        current_process = Settings::yH125VBF;
   if ( input_file_name.Contains("VBFH126") )        current_process = Settings::yH126VBF;
   if ( input_file_name.Contains("VBFH130") )        current_process = Settings::yH130VBF;
   
   if ( input_file_name.Contains("WplusH120") && genExtInfo > 10)       current_process = Settings::yH120WHlep;
   if ( input_file_name.Contains("WplusH124") && genExtInfo > 10)       current_process = Settings::yH124WHlep;
   if ( input_file_name.Contains("WplusH125") && genExtInfo > 10)       current_process = Settings::yH125WHlep;
   if ( input_file_name.Contains("WplusH126") && genExtInfo > 10)       current_process = Settings::yH126WHlep;
   if ( input_file_name.Contains("WplusH130") && genExtInfo > 10)       current_process = Settings::yH130WHlep;
   
   if ( input_file_name.Contains("WminusH120") && genExtInfo > 10)      current_process = Settings::yH120WHlep;
   if ( input_file_name.Contains("WminusH124") && genExtInfo > 10)      current_process = Settings::yH124WHlep;
   if ( input_file_name.Contains("WminusH125") && genExtInfo > 10)      current_process = Settings::yH125WHlep;
   if ( input_file_name.Contains("WminusH126") && genExtInfo > 10)      current_process = Settings::yH126WHlep;
   if ( input_file_name.Contains("WminusH130") && genExtInfo > 10)      current_process = Settings::yH130WHlep;
      
   if ( input_file_name.Contains("ZH120") && genExtInfo > 10)           current_process = Settings::yH120ZHlep;
   if ( input_file_name.Contains("ZH124") && genExtInfo > 10)           current_process = Settings::yH124ZHlep;
   if ( input_file_name.Contains("ZH125") && genExtInfo > 10)           current_process = Settings::yH125ZHlep;
   if ( input_file_name.Contains("ZH126") && genExtInfo > 10)           current_process = Settings::yH126ZHlep;
   if ( input_file_name.Contains("ZH130") && genExtInfo > 10)           current_process = Settings::yH130ZHlep;
   
   if ( input_file_name.Contains("WplusH120") && genExtInfo <= 10)      current_process = Settings::yH120WHhad;
   if ( input_file_name.Contains("WplusH124") && genExtInfo <= 10)      current_process = Settings::yH124WHhad;
   if ( input_file_name.Contains("WplusH125") && genExtInfo <= 10)      current_process = Settings::yH125WHhad;
   if ( input_file_name.Contains("WplusH126") && genExtInfo <= 10)      current_process = Settings::yH126WHhad;
   if ( input_file_name.Contains("WplusH130") && genExtInfo <= 10)      current_process = Settings::yH130WHhad;
   
   if ( input_file_name.Contains("WminusH120") && genExtInfo <= 10)     current_process = Settings::yH120WHhad;
   if ( input_file_name.Contains("WminusH124") && genExtInfo <= 10)     current_process = Settings::yH124WHhad;
   if ( input_file_name.Contains("WminusH125") && genExtInfo <= 10)     current_process = Settings::yH125WHhad;
   if ( input_file_name.Contains("WminusH126") && genExtInfo <= 10)     current_process = Settings::yH126WHhad;
   if ( input_file_name.Contains("WminusH130") && genExtInfo <= 10)     current_process = Settings::yH130WHhad;
   
   if ( input_file_name.Contains("ZH120") && genExtInfo <= 10)          current_process = Settings::yH120ZHhad;
   if ( input_file_name.Contains("ZH124") && genExtInfo <= 10)          current_process = Settings::yH124ZHhad;
   if ( input_file_name.Contains("ZH125") && genExtInfo <= 10)          current_process = Settings::yH125ZHhad;
   if ( input_file_name.Contains("ZH126") && genExtInfo <= 10)          current_process = Settings::yH126ZHhad;
   if ( input_file_name.Contains("ZH130") && genExtInfo <= 10)          current_process = Settings::yH130ZHhad;
   
   if ( input_file_name.Contains("ttH120") && n_gen_assoc_lep > 0)      current_process = Settings::yH120ttHlep;
   if ( input_file_name.Contains("ttH124") && n_gen_assoc_lep > 0)      current_process = Settings::yH124ttHlep;
   if ( input_file_name.Contains("ttH125") && n_gen_assoc_lep > 0)      current_process = Settings::yH125ttHlep;
   if ( input_file_name.Contains("ttH126") && n_gen_assoc_lep > 0)      current_process = Settings::yH126ttHlep;
   if ( input_file_name.Contains("ttH130") && n_gen_assoc_lep > 0)      current_process = Settings::yH130ttHlep;
   
   if ( input_file_name.Contains("ttH120") && n_gen_assoc_lep == 0)     current_process = Settings::yH120ttHhad;
   if ( input_file_name.Contains("ttH124") && n_gen_assoc_lep == 0)     current_process = Settings::yH124ttHhad;
   if ( input_file_name.Contains("ttH125") && n_gen_assoc_lep == 0)     current_process = Settings::yH125ttHhad;
   if ( input_file_name.Contains("ttH126") && n_gen_assoc_lep == 0)     current_process = Settings::yH126ttHhad;
   if ( input_file_name.Contains("ttH130") && n_gen_assoc_lep == 0)     current_process = Settings::yH130ttHhad;
   
   if ( input_file_name.Contains("bbH120") )     current_process = Settings::yH120bbH;
   if ( input_file_name.Contains("bbH124") )     current_process = Settings::yH124bbH;
   if ( input_file_name.Contains("bbH125") )     current_process = Settings::yH125bbH;
   if ( input_file_name.Contains("bbH126") )     current_process = Settings::yH126bbH;
   if ( input_file_name.Contains("bbH130") )     current_process = Settings::yH130bbH;
   
   if ( input_file_name.Contains("tqH120") )     current_process = Settings::yH120tqH;
   if ( input_file_name.Contains("tqH124") )     current_process = Settings::yH124tqH;
   if ( input_file_name.Contains("tqH125") )     current_process = Settings::yH125tqH;
   if ( input_file_name.Contains("tqH126") )     current_process = Settings::yH126tqH;
   if ( input_file_name.Contains("tqH130") )     current_process = Settings::yH130tqH;
      
   if ( input_file_name.Contains("ZZTo4l") )         current_process = Settings::yqqZZ;
   if ( input_file_name.Contains("ggTo4e") )         current_process = Settings::yggZZ;
   if ( input_file_name.Contains("ggTo4mu") )        current_process = Settings::yggZZ;
   if ( input_file_name.Contains("ggTo4tau") )       current_process = Settings::yggZZ;
   if ( input_file_name.Contains("ggTo2e2mu") )      current_process = Settings::yggZZ;
   if ( input_file_name.Contains("ggTo2e2tau") )     current_process = Settings::yggZZ;
   if ( input_file_name.Contains("ggTo2mu2tau") )    current_process = Settings::yggZZ;
   if ( input_file_name.Contains("DYJetsToLL_M50") ) current_process = Settings::yDY;
   if ( input_file_name.Contains("TTJets") )         current_process = Settings::yttbar;
   if ( input_file_name.Contains("TTTo2L2Nu") )      current_process = Settings::yttbar;
   // End assign dataset to correct process

   return current_process;
}
//==========================================================



//=================================
float Yields::calculate_K_factor(TString input_file_name)
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



//===========================
int Yields::FindFinalState()
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
   
   if ( _merge_2e2mu && final_state == Settings::fs2mu2e ) final_state = Settings::fs2e2mu;

   return final_state;
}
//===========================

//=============================
int Yields::FindFinalStateZX()
{
   int final_state = -999;
   
   if ( Z1Flav == -121 )
   {
      if ( Z2Flav == +121 )
         final_state = Settings::fs4e;
      else if ( Z2Flav == +169 )
         final_state = Settings::fs2e2mu;
      else
         cerr << "[ERROR] in event " << RunNumber << ":" << LumiNumber << ":" << EventNumber << ", Z2Flav = " << Z2Flav << endl;
   }
   else if ( Z1Flav == -169 )
   {
      if ( Z2Flav == +121 )
         final_state = Settings::fs2mu2e;
      else if ( Z2Flav == +169 )
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


//=============================
int Yields::CountAssociatedLeptons()
{
   int n_gen_assoc_leptons = 0;
   for ( int i_gen_assoc_lep = 0; i_gen_assoc_lep < 2; i_gen_assoc_lep++ )
      {
         if ( abs(gen_assoc_lep_id_.at(i_gen_assoc_lep)) == 11 || abs(gen_assoc_lep_id_.at(i_gen_assoc_lep)) == 13 )
         {
            n_gen_assoc_leptons++;
         }
      }
   return n_gen_assoc_leptons;
}
//=============================

