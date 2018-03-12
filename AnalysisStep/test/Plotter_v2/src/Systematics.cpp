// Include classes
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/include/Systematics.h>

// Constructor
//============================================================
Systematics::Systematics():Tree()
{
   _current_process = -999;
   _k_factor = 1;
   _current_final_state = -999;
   _current_category = -999;

   vector<float> temp;
   for ( int i_prod = 0; i_prod < num_of_production_modes; i_prod++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
      {
         temp.push_back(0);
      }
      _expected_yield.push_back(temp);
      _expected_yield_UP.push_back(temp);
      _expected_yield_DN.push_back(temp);
   }
	
	_s_category.push_back("UnTagged");
   _s_category.push_back("VBF1jTagged");
   _s_category.push_back("VBF2jTagged");
   _s_category.push_back("VHLeptTagged");
   _s_category.push_back("VHHadrTagged");
   _s_category.push_back("ttHLeptTagged");
   _s_category.push_back("ttHHadrTagged");
   _s_category.push_back("VHMETTagged");
   _s_category.push_back("Inclusive");
	
   _s_production_mode.push_back("ggH");
   _s_production_mode.push_back("VBF");
   _s_production_mode.push_back("WH_lep");
   _s_production_mode.push_back("WH_had");
   _s_production_mode.push_back("ZH_lep");
   _s_production_mode.push_back("ZH_had");
   _s_production_mode.push_back("ttH_lep");
   _s_production_mode.push_back("ttH_had");
   _s_production_mode.push_back("bbH");
   _s_production_mode.push_back("tqH");
}
//============================================================



// Destructor
//====================
Systematics::~Systematics()
{
}
//====================



//==========================================================
void Systematics::PrintSystematics_JEC( TString input_file_name)
{
	input_file = new TFile(input_file_name);

   hCounters = (TH1F*)input_file->Get("ZZTree/Counters");
   n_gen_events = (Long64_t)hCounters->GetBinContent(1);
   gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
	
   input_tree = (TTree*)input_file->Get("ZZTree/candTree");
   Init( input_tree, input_file_name );
	
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
	
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

      // Find current production mode
      gen_assoc_lep_id_.push_back(GenAssocLep1Id);
   	gen_assoc_lep_id_.push_back(GenAssocLep2Id);
      _n_gen_assoc_lep = CountAssociatedLeptons();
		_current_production_mode = find_current_production_mode( input_file_name, genExtInfo , _n_gen_assoc_lep);
      gen_assoc_lep_id_.clear();
		
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
													 true,// Use VHMET category
													 false);// Use QG tagging
		
		_current_category_UP = categoryMor18(nExtraLep,
											 nExtraZ,
											 nCleanedJetsPt30_jecUp,
											 nCleanedJetsPt30BTagged_bTagSF_jecUp,
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
											 PFMET_jesUp,
											 true,// Use VHMET category
											 false);// Use QG tagging
		
		_current_category_DN = categoryMor18(nExtraLep,
													 nExtraZ,
													 nCleanedJetsPt30_jecDn,
													 nCleanedJetsPt30BTagged_bTagSF_jecDn,
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
													 PFMET_jesDn,
													 true,// Use VHMET category
													 false);// Use QG tagging
      // K factors
      _k_factor = calculate_K_factor(input_file_name);

      // Final event weight
      _event_weight    = ( xsec * _k_factor * overallEventWeight) / gen_sum_weights;
      _event_weight_UP = ( xsec * _k_factor * overallEventWeight) / gen_sum_weights;
      _event_weight_DN = ( xsec * _k_factor * overallEventWeight) / gen_sum_weights;
      //if ( input_file_name.Contains("ggH") ) _event_weight *= ggH_NNLOPS_weight; // reweight POWHEG ggH to NNLOPS

      // Fill M4l histograms
       _expected_yield[_current_production_mode][_current_category]       += _event_weight;
       _expected_yield_UP[_current_production_mode][_current_category_UP] += _event_weight_UP;
       _expected_yield_DN[_current_production_mode][_current_category_DN] += _event_weight_DN;
   } // end for loop
	
	SumInclusive();
   cout << "[INFO] Systematics for " << _s_production_mode.at(_current_production_mode) << endl;
	
	
//	for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
//	{
//		cout << "Nominal = " << _expected_yield[_current_production_mode][i_cat]/_expected_yield[_current_production_mode][Settings::inclusive] << endl;
//		cout << "Up = " << _expected_yield_UP[_current_production_mode][i_cat]/_expected_yield_UP[_current_production_mode][Settings::inclusive] << endl;
//		cout << "Down = " << _expected_yield_DN[_current_production_mode][i_cat]/_expected_yield_DN[_current_production_mode][Settings::inclusive] << endl;
//		cout << endl;
//	}
	
	for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
	{
		cout << _s_category.at(i_cat) << " : " << (_expected_yield_DN[_current_production_mode][i_cat]/_expected_yield_DN[_current_production_mode][Settings::inclusive])/(_expected_yield[_current_production_mode][i_cat]/_expected_yield[_current_production_mode][Settings::inclusive]) << " / " <<
		   (_expected_yield_UP[_current_production_mode][i_cat]/_expected_yield_UP[_current_production_mode][Settings::inclusive])/(_expected_yield[_current_production_mode][i_cat]/_expected_yield[_current_production_mode][Settings::inclusive]) << endl;
	}
	
}
//==========================================================


//=================================
void Systematics::SumInclusive()
{
	for ( int i_prod = 0; i_prod < num_of_production_modes; i_prod++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
      {
		 _expected_yield[i_prod][Settings::inclusive]    += _expected_yield[i_prod][i_cat];
       _expected_yield_UP[i_prod][Settings::inclusive] += _expected_yield_UP[i_prod][i_cat];
       _expected_yield_DN[i_prod][Settings::inclusive] += _expected_yield_DN[i_prod][i_cat];
      }
   }
}
//=================================


//==========================================================
int Systematics::find_current_production_mode( TString input_file_name , int genExtInfo, int n_gen_assoc_lep)
{
	
   int current_production_mode = -999;

   // Assign dataset to correct process
   if ( input_file_name.Contains("ggH125") )         							current_production_mode = Settings::ggH;
   if ( input_file_name.Contains("VBFH125") )        							current_production_mode = Settings::VBF;
   if ( input_file_name.Contains("WplusH125")  && genExtInfo > 10)   	current_production_mode = Settings::WH_lep; //prepare for splitting VH and ttH into lep and had
   if ( input_file_name.Contains("WminusH125") && genExtInfo > 10)   	current_production_mode = Settings::WH_lep;
   if ( input_file_name.Contains("ZH125")      && genExtInfo > 10)   	current_production_mode = Settings::ZH_lep;
   if ( input_file_name.Contains("ttH125")     && n_gen_assoc_lep > 0)  current_production_mode = Settings::ttH_lep;
	if ( input_file_name.Contains("WplusH125")  && genExtInfo <= 10)  	current_production_mode = Settings::WH_had; //prepare for splitting VH and ttH into lep and had
	if ( input_file_name.Contains("WminusH125") && genExtInfo <= 10)  	current_production_mode = Settings::WH_had;
	if ( input_file_name.Contains("ZH125")      && genExtInfo <= 10)  	current_production_mode = Settings::ZH_had;
   if ( input_file_name.Contains("ttH125")     && n_gen_assoc_lep == 0) current_production_mode = Settings::ttH_had;
   if ( input_file_name.Contains("bbH125") )         							current_production_mode = Settings::bbH;
   if ( input_file_name.Contains("tqH125") )         							current_production_mode = Settings::tqH;

   // End assign dataset to correct process
	
   return current_production_mode;
}
//==========================================================

//=================================
float Systematics::calculate_K_factor(TString input_file_name)
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
int Systematics::FindFinalState()
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
//===========================

//=============================
int Systematics::CountAssociatedLeptons()
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

