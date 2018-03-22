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
      _expected_yield_PU.push_back(temp);
      _expected_yield_PU_UP.push_back(temp);
      _expected_yield_PU_DN.push_back(temp);
		
      _expected_yield_JEC.push_back(temp);
      _expected_yield_JEC_UP.push_back(temp);
      _expected_yield_JEC_DN.push_back(temp);
		
		_expected_yield_BTag.push_back(temp);
      _expected_yield_BTag_UP.push_back(temp);
      _expected_yield_BTag_DN.push_back(temp);
		
		_expected_yield_muR.push_back(temp);
      _expected_yield_muR_UP.push_back(temp);
      _expected_yield_muR_DN.push_back(temp);
		
		_expected_yield_muF.push_back(temp);
      _expected_yield_muF_UP.push_back(temp);
      _expected_yield_muF_DN.push_back(temp);
		
		_expected_yield_As.push_back(temp);
      _expected_yield_As_UP.push_back(temp);
      _expected_yield_As_DN.push_back(temp);
		
		_expected_yield_PDF.push_back(temp);
      _expected_yield_PDF_UP.push_back(temp);
      _expected_yield_PDF_DN.push_back(temp);
		
		_expected_yield_EWCorr.push_back(temp);
      _expected_yield_EWCorr_UP.push_back(temp);
      _expected_yield_EWCorr_DN.push_back(temp);
		
		_expected_yield_PythiaScale.push_back(temp);
      _expected_yield_PythiaScale_UP.push_back(temp);
      _expected_yield_PythiaScale_DN.push_back(temp);
		
		_expected_yield_PythiaTune.push_back(temp);
      _expected_yield_PythiaTune_UP.push_back(temp);
      _expected_yield_PythiaTune_DN.push_back(temp);
		
      _expected_yield_THU_ggH_Mu.push_back(temp);
      _expected_yield_THU_ggH_Mu_UP.push_back(temp);
      _expected_yield_THU_ggH_Mu_DN.push_back(temp);
		
		_expected_yield_THU_ggH_Res.push_back(temp);
      _expected_yield_THU_ggH_Res_UP.push_back(temp);
      _expected_yield_THU_ggH_Res_DN.push_back(temp);
		
      _expected_yield_THU_ggH_Mig01.push_back(temp);
      _expected_yield_THU_ggH_Mig01_UP.push_back(temp);
      _expected_yield_THU_ggH_Mig01_DN.push_back(temp);
		
      _expected_yield_THU_ggH_Mig12.push_back(temp);
      _expected_yield_THU_ggH_Mig12_UP.push_back(temp);
      _expected_yield_THU_ggH_Mig12_DN.push_back(temp);
		
      _expected_yield_THU_ggH_VBF2j.push_back(temp);
      _expected_yield_THU_ggH_VBF2j_UP.push_back(temp);
      _expected_yield_THU_ggH_VBF2j_DN.push_back(temp);
		
      _expected_yield_THU_ggH_VBF3j.push_back(temp);
      _expected_yield_THU_ggH_VBF3j_UP.push_back(temp);
      _expected_yield_THU_ggH_VBF3j_DN.push_back(temp);
		
      _expected_yield_THU_ggH_PT60.push_back(temp);
      _expected_yield_THU_ggH_PT60_UP.push_back(temp);
      _expected_yield_THU_ggH_PT60_DN.push_back(temp);
		
      _expected_yield_THU_ggH_PT120.push_back(temp);
      _expected_yield_THU_ggH_PT120_UP.push_back(temp);
      _expected_yield_THU_ggH_PT120_DN.push_back(temp);
		
      _expected_yield_THU_ggH_qmtop.push_back(temp);
      _expected_yield_THU_ggH_qmtop_UP.push_back(temp);
      _expected_yield_THU_ggH_qmtop_DN.push_back(temp);
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
   _s_production_mode.push_back("qqH");
   _s_production_mode.push_back("WH_lep");
	_s_production_mode.push_back("WH_had");
   _s_production_mode.push_back("ZH_lep");
	_s_production_mode.push_back("ZH_had");
   _s_production_mode.push_back("ttH_lep");
   _s_production_mode.push_back("ttH_had");
   _s_production_mode.push_back("tqH");
   _s_production_mode.push_back("bbH");
   _s_production_mode.push_back("qqZZ");
   _s_production_mode.push_back("ggZZ");
}
//============================================================



// Destructor
//====================
Systematics::~Systematics()
{
}
//====================



//==========================================================
void Systematics::FillSystematics( TString input_file_name)
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
      if ( ZZMass < 105. || ZZMass > 140. ) continue;

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
		
		//============================================================
		// PileUp
		//============================================================
		
		// K factors
      _k_factor 			    = calculate_K_factor(input_file_name);
		
      // Final event weight
      _event_weight    = ( xsec * _k_factor * overallEventWeight * ( 1. )                  						     ) / gen_sum_weights;
      _event_weight_UP = ( xsec * _k_factor * overallEventWeight * (PUWeight == 0. ? 0. : PUWeight_Up/PUWeight)  ) / gen_sum_weights;
      _event_weight_DN = ( xsec * _k_factor * overallEventWeight * (PUWeight == 0. ? 0. : PUWeight_Dn/PUWeight)  ) / gen_sum_weights;
      //if ( input_file_name.Contains("ggH") ) _event_weight *= ggH_NNLOPS_weight; // reweight POWHEG ggH to NNLOPS

      // Sum yields
       _expected_yield_PU[_current_production_mode][_current_category]           += _event_weight;
       _expected_yield_PU_UP[_current_production_mode][_current_category]        += _event_weight_UP;
       _expected_yield_PU_DN[_current_production_mode][_current_category]        += _event_weight_DN;
		
       _expected_yield_PU[_current_production_mode][Settings::inclusive]         += _event_weight;
       _expected_yield_PU_UP[_current_production_mode][Settings::inclusive]      += _event_weight_UP;
       _expected_yield_PU_DN[_current_production_mode][Settings::inclusive]      += _event_weight_DN;
		
		//============================================================
		// JEC
		//============================================================
		
		_current_category_JEC_UP = categoryMor18(nExtraLep,
													 nExtraZ,
													 nCleanedJetsPt30_jecUp,
													 nCleanedJetsPt30BTagged_bTagSF_jecUp,
													 jetQGL,
													 p_JJQCD_SIG_ghg2_1_JHUGen_JECUp,
													 p_JQCD_SIG_ghg2_1_JHUGen_JECUp,
													 p_JJVBF_SIG_ghv1_1_JHUGen_JECUp,
													 p_JVBF_SIG_ghv1_1_JHUGen_JECUp,
													 pAux_JVBF_SIG_ghv1_1_JHUGen_JECUp,
													 p_HadWH_SIG_ghw1_1_JHUGen_JECUp,
													 p_HadZH_SIG_ghz1_1_JHUGen_JECUp,
													 p_HadWH_mavjj_JECUp,
													 p_HadWH_mavjj_true_JECUp,
													 p_HadZH_mavjj_JECUp,
													 p_HadZH_mavjj_true_JECUp,
													 jetPhi,
													 ZZMass,
													 PFMET_jesUp,
													 true,// Use VHMET category
													 false);// Use QG tagging
		
		_current_category_JEC_DN = categoryMor18(nExtraLep,
													 nExtraZ,
													 nCleanedJetsPt30_jecDn,
													 nCleanedJetsPt30BTagged_bTagSF_jecDn,
													 jetQGL,
													 p_JJQCD_SIG_ghg2_1_JHUGen_JECDn,
													 p_JQCD_SIG_ghg2_1_JHUGen_JECDn,
													 p_JJVBF_SIG_ghv1_1_JHUGen_JECDn,
													 p_JVBF_SIG_ghv1_1_JHUGen_JECDn,
													 pAux_JVBF_SIG_ghv1_1_JHUGen_JECDn,
													 p_HadWH_SIG_ghw1_1_JHUGen_JECDn,
													 p_HadZH_SIG_ghz1_1_JHUGen_JECDn,
													 p_HadWH_mavjj_JECDn,
													 p_HadWH_mavjj_true_JECDn,
													 p_HadZH_mavjj_JECDn,
													 p_HadZH_mavjj_true_JECDn,
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

      // Sum yields
       _expected_yield_JEC[_current_production_mode][_current_category]           += _event_weight;
       _expected_yield_JEC_UP[_current_production_mode][_current_category_JEC_UP] += _event_weight_UP;
       _expected_yield_JEC_DN[_current_production_mode][_current_category_JEC_DN] += _event_weight_DN;
		
       _expected_yield_JEC[_current_production_mode][Settings::inclusive]    += _event_weight;
       _expected_yield_JEC_UP[_current_production_mode][Settings::inclusive] += _event_weight_UP;
       _expected_yield_JEC_DN[_current_production_mode][Settings::inclusive] += _event_weight_DN;
		
		//============================================================
		// BTag
		//============================================================
		
		_current_category_BTag_UP = categoryMor18(nExtraLep,
													 nExtraZ,
													 nCleanedJetsPt30,
													 nCleanedJetsPt30BTagged_bTagSFUp,
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
		
		_current_category_BTag_DN = categoryMor18(nExtraLep,
													 nExtraZ,
													 nCleanedJetsPt30,
													 nCleanedJetsPt30BTagged_bTagSFDn,
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
      // K factors
      _k_factor = calculate_K_factor(input_file_name);

      // Final event weight
      _event_weight    = ( xsec * _k_factor * overallEventWeight) / gen_sum_weights;
      _event_weight_UP = ( xsec * _k_factor * overallEventWeight) / gen_sum_weights;
      _event_weight_DN = ( xsec * _k_factor * overallEventWeight) / gen_sum_weights;
      //if ( input_file_name.Contains("ggH") ) _event_weight *= ggH_NNLOPS_weight; // reweight POWHEG ggH to NNLOPS

      // Sum yields
       _expected_yield_BTag[_current_production_mode][_current_category]           += _event_weight;
       _expected_yield_BTag_UP[_current_production_mode][_current_category_BTag_UP] += _event_weight_UP;
       _expected_yield_BTag_DN[_current_production_mode][_current_category_BTag_DN] += _event_weight_DN;
		
		 _expected_yield_BTag[_current_production_mode][Settings::inclusive]           += _event_weight;
       _expected_yield_BTag_UP[_current_production_mode][Settings::inclusive] += _event_weight_UP;
       _expected_yield_BTag_DN[_current_production_mode][Settings::inclusive] += _event_weight_DN;
		
		
		//============================================================
		// muR Scale
		//============================================================
		
		// K factors
      _k_factor 			    = calculate_K_factor(input_file_name);
		
      // Final event weight
      _event_weight    = ( xsec * _k_factor * overallEventWeight * ( 1. )                      											  ) / gen_sum_weights;
      _event_weight_UP = ( xsec * _k_factor * overallEventWeight * (LHEweight_QCDscale_muR2_muF1/LHEweight_QCDscale_muR1_muF1)  ) / gen_sum_weights;
      _event_weight_DN = ( xsec * _k_factor * overallEventWeight * (LHEweight_QCDscale_muR0p5_muF1/LHEweight_QCDscale_muR1_muF1)) / gen_sum_weights;
      //if ( input_file_name.Contains("ggH") ) _event_weight *= ggH_NNLOPS_weight; // reweight POWHEG ggH to NNLOPS

      // Sum yields
       _expected_yield_muR[_current_production_mode][_current_category]           += _event_weight;
       _expected_yield_muR_UP[_current_production_mode][_current_category]        += _event_weight_UP;
       _expected_yield_muR_DN[_current_production_mode][_current_category]        += _event_weight_DN;
		
       _expected_yield_muR[_current_production_mode][Settings::inclusive]           += _event_weight;
       _expected_yield_muR_UP[_current_production_mode][Settings::inclusive]        += _event_weight_UP;
       _expected_yield_muR_DN[_current_production_mode][Settings::inclusive]        += _event_weight_DN;
		
		//============================================================
		// muF Scale
		//============================================================
		
		// K factors
      _k_factor 			    = calculate_K_factor(input_file_name);
		
      // Final event weight
      _event_weight    = ( xsec * _k_factor * overallEventWeight * ( 1. )                      											  ) / gen_sum_weights;
      _event_weight_UP = ( xsec * _k_factor * overallEventWeight * (LHEweight_QCDscale_muR1_muF2/LHEweight_QCDscale_muR1_muF1)  ) / gen_sum_weights;
      _event_weight_DN = ( xsec * _k_factor * overallEventWeight * (LHEweight_QCDscale_muR1_muF0p5/LHEweight_QCDscale_muR1_muF1)) / gen_sum_weights;
      //if ( input_file_name.Contains("ggH") ) _event_weight *= ggH_NNLOPS_weight; // reweight POWHEG ggH to NNLOPS

      // Sum yields
       _expected_yield_muF[_current_production_mode][_current_category]           += _event_weight;
       _expected_yield_muF_UP[_current_production_mode][_current_category]        += _event_weight_UP;
       _expected_yield_muF_DN[_current_production_mode][_current_category]        += _event_weight_DN;
		
       _expected_yield_muF[_current_production_mode][Settings::inclusive]           += _event_weight;
       _expected_yield_muF_UP[_current_production_mode][Settings::inclusive]        += _event_weight_UP;
       _expected_yield_muF_DN[_current_production_mode][Settings::inclusive]        += _event_weight_DN;
		
		//============================================================
		// Alpha strong
		//============================================================
		
		// K factors
      _k_factor 			    = calculate_K_factor(input_file_name);
		
      // Final event weight
      _event_weight    = ( xsec * _k_factor * overallEventWeight * ( 1. )              ) / gen_sum_weights;
      _event_weight_UP = ( xsec * _k_factor * overallEventWeight * (LHEweight_AsMZ_Up) ) / gen_sum_weights;
      _event_weight_DN = ( xsec * _k_factor * overallEventWeight * (LHEweight_AsMZ_Dn) ) / gen_sum_weights;
      //if ( input_file_name.Contains("ggH") ) _event_weight *= ggH_NNLOPS_weight; // reweight POWHEG ggH to NNLOPS

      // Sum yields
       _expected_yield_As[_current_production_mode][_current_category]           += _event_weight;
       _expected_yield_As_UP[_current_production_mode][_current_category]        += _event_weight_UP;
       _expected_yield_As_DN[_current_production_mode][_current_category]        += _event_weight_DN;
		
       _expected_yield_As[_current_production_mode][Settings::inclusive]         += _event_weight;
       _expected_yield_As_UP[_current_production_mode][Settings::inclusive]      += _event_weight_UP;
       _expected_yield_As_DN[_current_production_mode][Settings::inclusive]      += _event_weight_DN;
		
		//============================================================
		// PDF Variation
		//============================================================
		
		// K factors
      _k_factor 			    = calculate_K_factor(input_file_name);
		
      // Final event weight
      _event_weight    = ( xsec * _k_factor * overallEventWeight * ( 1. )                      ) / gen_sum_weights;
      _event_weight_UP = ( xsec * _k_factor * overallEventWeight * (LHEweight_PDFVariation_Up) ) / gen_sum_weights;
      _event_weight_DN = ( xsec * _k_factor * overallEventWeight * (LHEweight_PDFVariation_Dn) ) / gen_sum_weights;
      //if ( input_file_name.Contains("ggH") ) _event_weight *= ggH_NNLOPS_weight; // reweight POWHEG ggH to NNLOPS

      // Sum yields
       _expected_yield_PDF[_current_production_mode][_current_category]           += _event_weight;
       _expected_yield_PDF_UP[_current_production_mode][_current_category]        += _event_weight_UP;
       _expected_yield_PDF_DN[_current_production_mode][_current_category]        += _event_weight_DN;
		
       _expected_yield_PDF[_current_production_mode][Settings::inclusive]         += _event_weight;
       _expected_yield_PDF_UP[_current_production_mode][Settings::inclusive]      += _event_weight_UP;
       _expected_yield_PDF_DN[_current_production_mode][Settings::inclusive]      += _event_weight_DN;
		
		//============================================================
		// EW corrections
		//============================================================
		
		// K factors
      _k_factor 			    = calculate_K_factor(input_file_name);
      if(input_file_name.Contains("ZZTo4l"))
      {
			_k_factor_EWCorr_UP = KFactor_EW_qqZZ * ( 1. + KFactor_EW_qqZZ_unc ) * KFactor_QCD_qqZZ_M;
			_k_factor_EWCorr_DN = KFactor_EW_qqZZ * ( 1. - KFactor_EW_qqZZ_unc ) * KFactor_QCD_qqZZ_M;
		}
		
      // Final event weight
      _event_weight    = ( xsec * _k_factor           * overallEventWeight ) / gen_sum_weights;
      _event_weight_UP = ( xsec * _k_factor_EWCorr_UP * overallEventWeight ) / gen_sum_weights;
      _event_weight_DN = ( xsec * _k_factor_EWCorr_DN * overallEventWeight ) / gen_sum_weights;
      //if ( input_file_name.Contains("ggH") ) _event_weight *= ggH_NNLOPS_weight; // reweight POWHEG ggH to NNLOPS

      // Sum yields
       _expected_yield_EWCorr[_current_production_mode][_current_category]           += _event_weight;
       _expected_yield_EWCorr_UP[_current_production_mode][_current_category]        += _event_weight_UP;
       _expected_yield_EWCorr_DN[_current_production_mode][_current_category]        += _event_weight_DN;
		
       _expected_yield_EWCorr[_current_production_mode][Settings::inclusive]         += _event_weight;
       _expected_yield_EWCorr_UP[_current_production_mode][Settings::inclusive]      += _event_weight_UP;
       _expected_yield_EWCorr_DN[_current_production_mode][Settings::inclusive]      += _event_weight_DN;
		
		
		//============================================================
		// Pythia Scale
		//============================================================
		
		// K factors
      _k_factor 			    = calculate_K_factor(input_file_name);
		
      // Final event weight
      _event_weight    = ( xsec * _k_factor * overallEventWeight * ( 1. )                                                ) / gen_sum_weights;
      _event_weight_UP = ( xsec * _k_factor * overallEventWeight * (PythiaWeight_isr_muR4    * PythiaWeight_fsr_muR4   ) ) / gen_sum_weights;
      _event_weight_DN = ( xsec * _k_factor * overallEventWeight * (PythiaWeight_isr_muR0p25 * PythiaWeight_fsr_muR0p25) ) / gen_sum_weights;
      //if ( input_file_name.Contains("ggH") ) _event_weight *= ggH_NNLOPS_weight; // reweight POWHEG ggH to NNLOPS

      // Sum yields
       _expected_yield_PythiaScale[_current_production_mode][_current_category]           += _event_weight;
       _expected_yield_PythiaScale_UP[_current_production_mode][_current_category]        += _event_weight_UP;
       _expected_yield_PythiaScale_DN[_current_production_mode][_current_category]        += _event_weight_DN;
		
       _expected_yield_PythiaScale[_current_production_mode][Settings::inclusive]         += _event_weight;
       _expected_yield_PythiaScale_UP[_current_production_mode][Settings::inclusive]      += _event_weight_UP;
       _expected_yield_PythiaScale_DN[_current_production_mode][Settings::inclusive]      += _event_weight_DN;
		
		//============================================================
		// Pythia Tune
		//============================================================
		
		// K factors
      _k_factor 			    = calculate_K_factor(input_file_name);
		
      // Final event weight
      _event_weight    = ( xsec * _k_factor * overallEventWeight * ( 1. )                                                ) / gen_sum_weights;
      //if ( input_file_name.Contains("ggH") ) _event_weight *= ggH_NNLOPS_weight; // reweight POWHEG ggH to NNLOPS

      // Sum yields
       _expected_yield_PythiaTune[_current_production_mode][_current_category]           += _event_weight;
       _expected_yield_PythiaTune[_current_production_mode][Settings::inclusive]         += _event_weight;
		//Up/Dn variations have dedicated samples and are filled in FillSystematics_tuneUpDn() function

		
		//============================================================
		// THU_ggH uncertainties
		//============================================================
		if(input_file_name.Contains("ggH125"))
      {
			// K factors
			_k_factor 			    = calculate_K_factor(input_file_name);
			
			// Final event weight
			_event_weight    = ( xsec * _k_factor           * overallEventWeight ) / gen_sum_weights;
			_event_weight_UP = ( xsec * _k_factor * overallEventWeight * (qcd_ggF_uncertSF->at(0))) / gen_sum_weights;
			_event_weight_DN = ( xsec * _k_factor * overallEventWeight * (1. - (qcd_ggF_uncertSF->at(0) - 1.))) / gen_sum_weights;
			
			//if ( input_file_name.Contains("ggH") ) _event_weight *= ggH_NNLOPS_weight; // reweight POWHEG ggH to NNLOPS

			// Sum yields
			 _expected_yield_THU_ggH_Mu[_current_production_mode][_current_category]           += _event_weight;
			 _expected_yield_THU_ggH_Mu_UP[_current_production_mode][_current_category]        += _event_weight_UP;
			 _expected_yield_THU_ggH_Mu_DN[_current_production_mode][_current_category]        += _event_weight_DN;
			
			 _expected_yield_THU_ggH_Mu[_current_production_mode][Settings::inclusive]         += _event_weight;
			 _expected_yield_THU_ggH_Mu_UP[_current_production_mode][Settings::inclusive]      += _event_weight_UP;
			 _expected_yield_THU_ggH_Mu_DN[_current_production_mode][Settings::inclusive]      += _event_weight_DN;
			
			// K factors
			_k_factor 			    = calculate_K_factor(input_file_name);
			
			// Final event weight
			_event_weight    = ( xsec * _k_factor           * overallEventWeight ) / gen_sum_weights;
			_event_weight_UP = ( xsec * _k_factor * overallEventWeight * (qcd_ggF_uncertSF->at(1))) / gen_sum_weights;
			_event_weight_DN = ( xsec * _k_factor * overallEventWeight * (1. - (qcd_ggF_uncertSF->at(1) - 1.))) / gen_sum_weights;
			
			//if ( input_file_name.Contains("ggH") ) _event_weight *= ggH_NNLOPS_weight; // reweight POWHEG ggH to NNLOPS

			// Sum yields
			 _expected_yield_THU_ggH_Res[_current_production_mode][_current_category]           += _event_weight;
			 _expected_yield_THU_ggH_Res_UP[_current_production_mode][_current_category]        += _event_weight_UP;
			 _expected_yield_THU_ggH_Res_DN[_current_production_mode][_current_category]        += _event_weight_DN;
			
			 _expected_yield_THU_ggH_Res[_current_production_mode][Settings::inclusive]         += _event_weight;
			 _expected_yield_THU_ggH_Res_UP[_current_production_mode][Settings::inclusive]      += _event_weight_UP;
			 _expected_yield_THU_ggH_Res_DN[_current_production_mode][Settings::inclusive]      += _event_weight_DN;
			
			// K factors
			_k_factor 			    = calculate_K_factor(input_file_name);
			
			// Final event weight
			_event_weight    = ( xsec * _k_factor           * overallEventWeight ) / gen_sum_weights;
			_event_weight_UP = ( xsec * _k_factor * overallEventWeight * (qcd_ggF_uncertSF->at(2))) / gen_sum_weights;
			_event_weight_DN = ( xsec * _k_factor * overallEventWeight * (1. - (qcd_ggF_uncertSF->at(2) - 1.))) / gen_sum_weights;
			
			//if ( input_file_name.Contains("ggH") ) _event_weight *= ggH_NNLOPS_weight; // reweight POWHEG ggH to NNLOPS

			// Sum yields
			 _expected_yield_THU_ggH_Mig01[_current_production_mode][_current_category]           += _event_weight;
			 _expected_yield_THU_ggH_Mig01_UP[_current_production_mode][_current_category]        += _event_weight_UP;
			 _expected_yield_THU_ggH_Mig01_DN[_current_production_mode][_current_category]        += _event_weight_DN;
			
			 _expected_yield_THU_ggH_Mig01[_current_production_mode][Settings::inclusive]         += _event_weight;
			 _expected_yield_THU_ggH_Mig01_UP[_current_production_mode][Settings::inclusive]      += _event_weight_UP;
			 _expected_yield_THU_ggH_Mig01_DN[_current_production_mode][Settings::inclusive]      += _event_weight_DN;
			
			// K factors
			_k_factor 			    = calculate_K_factor(input_file_name);
			
			// Final event weight
			_event_weight    = ( xsec * _k_factor           * overallEventWeight ) / gen_sum_weights;
			_event_weight_UP = ( xsec * _k_factor * overallEventWeight * (qcd_ggF_uncertSF->at(3))) / gen_sum_weights;
			_event_weight_DN = ( xsec * _k_factor * overallEventWeight * (1. - (qcd_ggF_uncertSF->at(3) - 1.))) / gen_sum_weights;
			
			//if ( input_file_name.Contains("ggH") ) _event_weight *= ggH_NNLOPS_weight; // reweight POWHEG ggH to NNLOPS

			// Sum yields
			 _expected_yield_THU_ggH_Mig12[_current_production_mode][_current_category]           += _event_weight;
			 _expected_yield_THU_ggH_Mig12_UP[_current_production_mode][_current_category]        += _event_weight_UP;
			 _expected_yield_THU_ggH_Mig12_DN[_current_production_mode][_current_category]        += _event_weight_DN;
			
			 _expected_yield_THU_ggH_Mig12[_current_production_mode][Settings::inclusive]         += _event_weight;
			 _expected_yield_THU_ggH_Mig12_UP[_current_production_mode][Settings::inclusive]      += _event_weight_UP;
			 _expected_yield_THU_ggH_Mig12_DN[_current_production_mode][Settings::inclusive]      += _event_weight_DN;
			
			// K factors
			_k_factor 			    = calculate_K_factor(input_file_name);
			
			// Final event weight
			_event_weight    = ( xsec * _k_factor           * overallEventWeight ) / gen_sum_weights;
			_event_weight_UP = ( xsec * _k_factor * overallEventWeight * (qcd_ggF_uncertSF->at(4))) / gen_sum_weights;
			_event_weight_DN = ( xsec * _k_factor * overallEventWeight * (1. - (qcd_ggF_uncertSF->at(4) - 1.))) / gen_sum_weights;
			
			//if ( input_file_name.Contains("ggH") ) _event_weight *= ggH_NNLOPS_weight; // reweight POWHEG ggH to NNLOPS

			// Sum yields
			 _expected_yield_THU_ggH_VBF2j[_current_production_mode][_current_category]           += _event_weight;
			 _expected_yield_THU_ggH_VBF2j_UP[_current_production_mode][_current_category]        += _event_weight_UP;
			 _expected_yield_THU_ggH_VBF2j_DN[_current_production_mode][_current_category]        += _event_weight_DN;
			
			 _expected_yield_THU_ggH_VBF2j[_current_production_mode][Settings::inclusive]         += _event_weight;
			 _expected_yield_THU_ggH_VBF2j_UP[_current_production_mode][Settings::inclusive]      += _event_weight_UP;
			 _expected_yield_THU_ggH_VBF2j_DN[_current_production_mode][Settings::inclusive]      += _event_weight_DN;
			
			// K factors
			_k_factor 			    = calculate_K_factor(input_file_name);
			
			// Final event weight
			_event_weight    = ( xsec * _k_factor           * overallEventWeight ) / gen_sum_weights;
			_event_weight_UP = ( xsec * _k_factor * overallEventWeight * (qcd_ggF_uncertSF->at(5))) / gen_sum_weights;
			_event_weight_DN = ( xsec * _k_factor * overallEventWeight * (1. - (qcd_ggF_uncertSF->at(5) - 1.))) / gen_sum_weights;
			
			//if ( input_file_name.Contains("ggH") ) _event_weight *= ggH_NNLOPS_weight; // reweight POWHEG ggH to NNLOPS

			// Sum yields
			 _expected_yield_THU_ggH_VBF3j[_current_production_mode][_current_category]           += _event_weight;
			 _expected_yield_THU_ggH_VBF3j_UP[_current_production_mode][_current_category]        += _event_weight_UP;
			 _expected_yield_THU_ggH_VBF3j_DN[_current_production_mode][_current_category]        += _event_weight_DN;
			
			 _expected_yield_THU_ggH_VBF3j[_current_production_mode][Settings::inclusive]         += _event_weight;
			 _expected_yield_THU_ggH_VBF3j_UP[_current_production_mode][Settings::inclusive]      += _event_weight_UP;
			 _expected_yield_THU_ggH_VBF3j_DN[_current_production_mode][Settings::inclusive]      += _event_weight_DN;
			
			// K factors
			_k_factor 			    = calculate_K_factor(input_file_name);
			
			// Final event weight
			_event_weight    = ( xsec * _k_factor           * overallEventWeight ) / gen_sum_weights;
			_event_weight_UP = ( xsec * _k_factor * overallEventWeight * (qcd_ggF_uncertSF->at(6))) / gen_sum_weights;
			_event_weight_DN = ( xsec * _k_factor * overallEventWeight * (1. - (qcd_ggF_uncertSF->at(6) - 1.))) / gen_sum_weights;
			
			//if ( input_file_name.Contains("ggH") ) _event_weight *= ggH_NNLOPS_weight; // reweight POWHEG ggH to NNLOPS

			// Sum yields
			 _expected_yield_THU_ggH_PT60[_current_production_mode][_current_category]           += _event_weight;
			 _expected_yield_THU_ggH_PT60_UP[_current_production_mode][_current_category]        += _event_weight_UP;
			 _expected_yield_THU_ggH_PT60_DN[_current_production_mode][_current_category]        += _event_weight_DN;
			
			 _expected_yield_THU_ggH_PT60[_current_production_mode][Settings::inclusive]         += _event_weight;
			 _expected_yield_THU_ggH_PT60_UP[_current_production_mode][Settings::inclusive]      += _event_weight_UP;
			 _expected_yield_THU_ggH_PT60_DN[_current_production_mode][Settings::inclusive]      += _event_weight_DN;
			
			// K factors
			_k_factor 			    = calculate_K_factor(input_file_name);
			
			// Final event weight
			_event_weight    = ( xsec * _k_factor           * overallEventWeight ) / gen_sum_weights;
			_event_weight_UP = ( xsec * _k_factor * overallEventWeight * (qcd_ggF_uncertSF->at(7))) / gen_sum_weights;
			_event_weight_DN = ( xsec * _k_factor * overallEventWeight * (1. - (qcd_ggF_uncertSF->at(7) - 1.))) / gen_sum_weights;
			
			//if ( input_file_name.Contains("ggH") ) _event_weight *= ggH_NNLOPS_weight; // reweight POWHEG ggH to NNLOPS

			// Sum yields
			 _expected_yield_THU_ggH_PT120[_current_production_mode][_current_category]           += _event_weight;
			 _expected_yield_THU_ggH_PT120_UP[_current_production_mode][_current_category]        += _event_weight_UP;
			 _expected_yield_THU_ggH_PT120_DN[_current_production_mode][_current_category]        += _event_weight_DN;
			
			 _expected_yield_THU_ggH_PT120[_current_production_mode][Settings::inclusive]         += _event_weight;
			 _expected_yield_THU_ggH_PT120_UP[_current_production_mode][Settings::inclusive]      += _event_weight_UP;
			 _expected_yield_THU_ggH_PT120_DN[_current_production_mode][Settings::inclusive]      += _event_weight_DN;
			
			// K factors
			_k_factor 			    = calculate_K_factor(input_file_name);
			
			// Final event weight
			_event_weight    = ( xsec * _k_factor           * overallEventWeight ) / gen_sum_weights;
			_event_weight_UP = ( xsec * _k_factor * overallEventWeight * (qcd_ggF_uncertSF->at(8))) / gen_sum_weights;
			_event_weight_DN = ( xsec * _k_factor * overallEventWeight * (1. - (qcd_ggF_uncertSF->at(8) - 1.))) / gen_sum_weights;
			
			//if ( input_file_name.Contains("ggH") ) _event_weight *= ggH_NNLOPS_weight; // reweight POWHEG ggH to NNLOPS

			// Sum yields
			 _expected_yield_THU_ggH_qmtop[_current_production_mode][_current_category]           += _event_weight;
			 _expected_yield_THU_ggH_qmtop_UP[_current_production_mode][_current_category]        += _event_weight_UP;
			 _expected_yield_THU_ggH_qmtop_DN[_current_production_mode][_current_category]        += _event_weight_DN;
			
			 _expected_yield_THU_ggH_qmtop[_current_production_mode][Settings::inclusive]         += _event_weight;
			 _expected_yield_THU_ggH_qmtop_UP[_current_production_mode][Settings::inclusive]      += _event_weight_UP;
			 _expected_yield_THU_ggH_qmtop_DN[_current_production_mode][Settings::inclusive]      += _event_weight_DN;
			
       }
   } // end for loop
	
cout << "[INFO] Systematics for " << input_file_name << " filled." << endl;
}
//==========================================================


//==========================================================
void Systematics::FillSystematics_tuneUpDn( TString input_file_name)
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
      if ( ZZMass < 105. || ZZMass > 140. ) continue;

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
		
		//============================================================
		// Pythia Tune
		//============================================================
		
		// K factors
      _k_factor 			    = calculate_K_factor(input_file_name);
		
      // Final event weight
      _event_weight    = ( xsec * _k_factor * overallEventWeight * ( 1. )                                                ) / gen_sum_weights;
      //if ( input_file_name.Contains("ggH") ) _event_weight *= ggH_NNLOPS_weight; // reweight POWHEG ggH to NNLOPS

      // Sum yields
       if(input_file_name.Contains("tuneup"))   _expected_yield_PythiaTune_UP[_current_production_mode][_current_category]        += _event_weight;
       if(input_file_name.Contains("tunedown")) _expected_yield_PythiaTune_DN[_current_production_mode][_current_category]        += _event_weight;
		
       if(input_file_name.Contains("tuneup"))   _expected_yield_PythiaTune_UP[_current_production_mode][Settings::inclusive]      += _event_weight;
       if(input_file_name.Contains("tunedown")) _expected_yield_PythiaTune_DN[_current_production_mode][Settings::inclusive]      += _event_weight;
		
   } // end for loop
	
cout << "[INFO] Systematics for " << input_file_name << " filled." << endl;
}
//==========================================================



//========================================
void Systematics::PrintSystematics_PU( )
{
	
	for ( int i_prod = 0; i_prod < num_of_production_modes; i_prod++ )
   {
   	for ( int i_cat = 0; i_cat < num_of_categories ; i_cat++ )
		{
			cout << _s_category.at(i_cat) << endl;
			cout << "Nominal = " << _expected_yield_PU[i_prod][i_cat] << endl;
			cout << "Up = " << _expected_yield_PU_UP[i_prod][i_cat] << endl;
			cout << "Down = " << _expected_yield_PU_DN[i_prod][i_cat] << endl;
			cout << endl;
		}
   }
	
	for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
   {
   	cout << "[INFO] PileUp systematics for " << _s_category.at(i_cat) << endl;
		for ( int i_prod = 0; i_prod < num_of_production_modes; i_prod++ )
		{
			cout << _s_production_mode.at(i_prod)  << " : " << (_expected_yield_PU_UP[i_prod][i_cat]/_expected_yield_PU_UP[i_prod][Settings::inclusive])/(_expected_yield_PU[i_prod][i_cat]/_expected_yield_PU[i_prod][Settings::inclusive]) << "/" <<(_expected_yield_PU_DN[i_prod][i_cat]/_expected_yield_PU_DN[i_prod][Settings::inclusive])/(_expected_yield_PU[i_prod][i_cat]/_expected_yield_PU[i_prod][Settings::inclusive]) << endl;
		}
   }

	cout << "[INFO] PileUp systematics for inclusive: " << endl;
	for ( int i_prod = 0; i_prod < num_of_production_modes; i_prod++ )
	{
		cout << _s_production_mode.at(i_prod)  << " : " << (_expected_yield_PU_UP[i_prod][Settings::inclusive]/_expected_yield_PU[i_prod][Settings::inclusive]) << "/" << (_expected_yield_PU_DN[i_prod][Settings::inclusive]/_expected_yield_PU[i_prod][Settings::inclusive]) << endl;
	}
}
//========================================



//========================================
void Systematics::PrintSystematics_JEC( )
{
	
	for ( int i_prod = 0; i_prod < num_of_production_modes; i_prod++ )
   {
   	for ( int i_cat = 0; i_cat < num_of_categories ; i_cat++ )
		{
			cout << _s_category.at(i_cat) << endl;
			cout << "Nominal = " << _expected_yield_JEC[i_prod][i_cat] << endl;
			cout << "Up = " << _expected_yield_JEC_UP[i_prod][i_cat] << endl;
			cout << "Down = " << _expected_yield_JEC_DN[i_prod][i_cat] << endl;
			cout << endl;
		}
   }
	
	for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
   {
   	cout << "[INFO] JEC systematics for " << _s_category.at(i_cat) << endl;
		for ( int i_prod = 0; i_prod < num_of_production_modes; i_prod++ )
		{
			cout << _s_production_mode.at(i_prod)  << " : " <<(_expected_yield_JEC_UP[i_prod][i_cat]/_expected_yield_JEC_UP[i_prod][Settings::inclusive])/(_expected_yield_JEC[i_prod][i_cat]/_expected_yield_JEC[i_prod][Settings::inclusive]) << "/" <<(_expected_yield_JEC_DN[i_prod][i_cat]/_expected_yield_JEC_DN[i_prod][Settings::inclusive])/(_expected_yield_JEC[i_prod][i_cat]/_expected_yield_JEC[i_prod][Settings::inclusive]) << endl;
		}
   }

	
}
//========================================

//========================================
void Systematics::PrintSystematics_BTag( )
{
	
//	for ( int i_prod = 0; i_prod < num_of_production_modes; i_prod++ )
//   {
//   	for ( int i_cat = 0; i_cat < num_of_categories ; i_cat++ )
//		{
//			cout << _s_category.at(i_cat) << endl;
//			cout << "Nominal = " << _expected_yield_BTag[i_prod][i_cat] << endl;
//			cout << "Up = " << _expected_yield_BTag_UP[i_prod][i_cat] << endl;
//			cout << "Down = " << _expected_yield_BTag_DN[i_prod][i_cat] << endl;
//			cout << endl;
//		}
//   }
	
	for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
   {
   	cout << "[INFO] BTag systematics for " << _s_category.at(i_cat) << endl;
		for ( int i_prod = 0; i_prod < num_of_production_modes; i_prod++ )
		{
			cout << _s_production_mode.at(i_prod)  << " : " << (_expected_yield_BTag_UP[i_prod][i_cat]/_expected_yield_BTag_UP[i_prod][Settings::inclusive])/(_expected_yield_BTag[i_prod][i_cat]/_expected_yield_BTag[i_prod][Settings::inclusive]) << "/" <<(_expected_yield_BTag_DN[i_prod][i_cat]/_expected_yield_BTag_DN[i_prod][Settings::inclusive])/(_expected_yield_BTag[i_prod][i_cat]/_expected_yield_BTag[i_prod][Settings::inclusive]) << endl;
		}
   }

	
}
//========================================

//========================================
void Systematics::PrintSystematics_muRmuFScale( )
{
	
//	for ( int i_prod = 0; i_prod < num_of_production_modes; i_prod++ )
//   {
//   	if ( i_prod != Settings::ggH ) continue;
//   	cout << "===========================" << endl;
//   	cout << _s_production_mode.at(i_prod) << endl;
//   	cout << "===========================" << endl;
//
//   	for ( int i_cat = 0; i_cat < num_of_categories ; i_cat++ )
//		{
//			cout << _s_category.at(i_cat) << endl;
//			cout << "muR Nominal = " << _expected_yield_muR[i_prod][i_cat] << endl;
//			cout << "muR Up = " << _expected_yield_muR_UP[i_prod][i_cat] << endl;
//			cout << "muR Down = " << _expected_yield_muR_DN[i_prod][i_cat] << endl;
//			cout << endl;
//
//			cout << _s_category.at(i_cat) << endl;
//			cout << "muF Nominal = " << _expected_yield_muF[i_prod][i_cat] << endl;
//			cout << "muF Up = " << _expected_yield_muF_UP[i_prod][i_cat] << endl;
//			cout << "muF Down = " << _expected_yield_muF_DN[i_prod][i_cat] << endl;
//			cout << endl;
//
//		}
//   }
	
	float muR_Up = 1.;
	float muF_Up = 1.;
	
	float muR_Dn = 1.;
	float muF_Dn = 1.;
	
	int combination = 0;
	
	for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
   {
   	cout << "[INFO] QCD scale systematics for " << _s_category.at(i_cat) << endl;
		for ( int i_prod = 0; i_prod < num_of_production_modes; i_prod++ )
		{
			muR_Up = (_expected_yield_muR_UP[i_prod][i_cat]/_expected_yield_muR_UP[i_prod][Settings::inclusive])/(_expected_yield_muR[i_prod][i_cat]/_expected_yield_muR[i_prod][Settings::inclusive]);
			muR_Dn = (_expected_yield_muR_DN[i_prod][i_cat]/_expected_yield_muR_DN[i_prod][Settings::inclusive])/(_expected_yield_muR[i_prod][i_cat]/_expected_yield_muR[i_prod][Settings::inclusive]);
			
			muF_Up = (_expected_yield_muF_UP[i_prod][i_cat]/_expected_yield_muF_UP[i_prod][Settings::inclusive])/(_expected_yield_muF[i_prod][i_cat]/_expected_yield_muF[i_prod][Settings::inclusive]);
			muF_Dn = (_expected_yield_muF_DN[i_prod][i_cat]/_expected_yield_muF_DN[i_prod][Settings::inclusive])/(_expected_yield_muF[i_prod][i_cat]/_expected_yield_muF[i_prod][Settings::inclusive]);
			
			if( i_cat == Settings::untagged && muR_Up > 1. && muF_Up > 1.) combination = 0;
			if( i_cat == Settings::untagged && muR_Up > 1. && muF_Up < 1.) combination = 1;
			if( i_cat == Settings::untagged && muR_Up < 1. && muF_Up > 1.) combination = 2;
			if( i_cat == Settings::untagged && muR_Up < 1. && muF_Up < 1.) combination = 3;
			
//			cout << muR_Up << " " << muR_Dn << endl;
//			cout << muF_Up << " " << muF_Dn << endl;
//			cout << combination << endl;
			
			if(combination == 0) cout << _s_production_mode.at(i_prod)  << " : " << 1. - sqrt(pow(1. - muR_Dn,2)+pow(1. - muF_Dn,2)) << "/" << 1. + sqrt(pow(muR_Up - 1.,2)+pow(muF_Up - 1.,2)) << endl;
			if(combination == 1) cout << _s_production_mode.at(i_prod)  << " : " << 1. - sqrt(pow(1. - muR_Dn,2)+pow(1. - muF_Up,2)) << "/" << 1. + sqrt(pow(muR_Up - 1.,2)+pow(muF_Dn - 1.,2)) << endl;
			if(combination == 2) cout << _s_production_mode.at(i_prod)  << " : " << 1. - sqrt(pow(1. - muR_Up,2)+pow(1. - muF_Dn,2)) << "/" << 1. + sqrt(pow(muR_Dn - 1.,2)+pow(muF_Up - 1.,2)) << endl;
			if(combination == 3) cout << _s_production_mode.at(i_prod)  << " : " << 1. - sqrt(pow(1. - muR_Up,2)+pow(1. - muF_Up,2)) << "/" << 1. + sqrt(pow(muR_Dn - 1.,2)+pow(muF_Dn - 1.,2)) << endl;
		}
   }
}
//========================================


//========================================
void Systematics::PrintSystematics_PDFScale( )
{
	
//	for ( int i_prod = 0; i_prod < num_of_production_modes; i_prod++ )
//   {
//   	if ( i_prod != Settings::ggH ) continue;
//   	cout << "===========================" << endl;
//   	cout << _s_production_mode.at(i_prod) << endl;
//   	cout << "===========================" << endl;
//
//   	for ( int i_cat = 0; i_cat < num_of_categories ; i_cat++ )
//		{
//			cout << _s_category.at(i_cat) << endl;
//			cout << "muR Nominal = " << _expected_yield_muR[i_prod][i_cat] << endl;
//			cout << "muR Up = " << _expected_yield_As_Up[i_prod][i_cat] << endl;
//			cout << "muR Down = " << _expected_yield_As_Dn[i_prod][i_cat] << endl;
//			cout << endl;
//
//			cout << _s_category.at(i_cat) << endl;
//			cout << "muF Nominal = " << _expected_yield_muF[i_prod][i_cat] << endl;
//			cout << "muF Up = " << _expected_yield_PDF_Up[i_prod][i_cat] << endl;
//			cout << "muF Down = " << _expected_yield_PDF_Dn[i_prod][i_cat] << endl;
//			cout << endl;
//
//		}
//   }
	
	float As_Up = 1.;
	float PDF_Up = 1.;
	
	float As_Dn = 1.;
	float PDF_Dn = 1.;
	
	int combination = 0;
	
	for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
   {
   	cout << "[INFO] PDF scale systematics for " << _s_category.at(i_cat) << endl;
		for ( int i_prod = 0; i_prod < num_of_production_modes; i_prod++ )
		{
			As_Up = (_expected_yield_As_UP[i_prod][i_cat]/_expected_yield_As_UP[i_prod][Settings::inclusive])/(_expected_yield_As[i_prod][i_cat]/_expected_yield_As[i_prod][Settings::inclusive]);
			As_Dn = (_expected_yield_As_DN[i_prod][i_cat]/_expected_yield_As_DN[i_prod][Settings::inclusive])/(_expected_yield_As[i_prod][i_cat]/_expected_yield_As[i_prod][Settings::inclusive]);
			
			PDF_Up = (_expected_yield_PDF_UP[i_prod][i_cat]/_expected_yield_PDF_UP[i_prod][Settings::inclusive])/(_expected_yield_PDF[i_prod][i_cat]/_expected_yield_PDF[i_prod][Settings::inclusive]);
			PDF_Dn = (_expected_yield_PDF_DN[i_prod][i_cat]/_expected_yield_PDF_DN[i_prod][Settings::inclusive])/(_expected_yield_PDF[i_prod][i_cat]/_expected_yield_PDF[i_prod][Settings::inclusive]);
			
			if( i_cat == Settings::untagged && As_Up > 1. && PDF_Up > 1.) combination = 0;
			if( i_cat == Settings::untagged && As_Up > 1. && PDF_Up < 1.) combination = 1;
			if( i_cat == Settings::untagged && As_Up < 1. && PDF_Up > 1.) combination = 2;
			if( i_cat == Settings::untagged && As_Up < 1. && PDF_Up < 1.) combination = 3;
			
//			cout << As_Up << " " << As_Dn << endl;
//			cout << PDF_Up << " " << PDF_Dn << endl;
//			cout << combination << endl;
			
			if(combination == 0) cout << _s_production_mode.at(i_prod)  << " : " << 1. - sqrt(pow(1. - As_Dn,2)+pow(1. - PDF_Dn,2)) << "/" << 1. + sqrt(pow(As_Up - 1.,2)+pow(PDF_Up - 1.,2)) << endl;
			if(combination == 1) cout << _s_production_mode.at(i_prod)  << " : " << 1. - sqrt(pow(1. - As_Dn,2)+pow(1. - PDF_Up,2)) << "/" << 1. + sqrt(pow(As_Up - 1.,2)+pow(PDF_Dn - 1.,2)) << endl;
			if(combination == 2) cout << _s_production_mode.at(i_prod)  << " : " << 1. - sqrt(pow(1. - As_Up,2)+pow(1. - PDF_Dn,2)) << "/" << 1. + sqrt(pow(As_Dn - 1.,2)+pow(PDF_Up - 1.,2)) << endl;
			if(combination == 3) cout << _s_production_mode.at(i_prod)  << " : " << 1. - sqrt(pow(1. - As_Up,2)+pow(1. - PDF_Up,2)) << "/" << 1. + sqrt(pow(As_Dn - 1.,2)+pow(PDF_Dn - 1.,2)) << endl;
		}
   }
	
}
//========================================

//========================================
void Systematics::PrintSystematics_EWCorr( )
{
	
//	for ( int i_prod = 0; i_prod < num_of_production_modes; i_prod++ )
//   {
//   	for ( int i_cat = 0; i_cat < num_of_categories ; i_cat++ )
//		{
//			cout << _s_category.at(i_cat) << endl;
//			cout << "Nominal = " << _expected_yield_BTag[i_prod][i_cat] << endl;
//			cout << "Up = " << _expected_yield_BTag_UP[i_prod][i_cat] << endl;
//			cout << "Down = " << _expected_yield_BTag_DN[i_prod][i_cat] << endl;
//			cout << endl;
//		}
//   }
	
	for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
   {
   	cout << "[INFO] EWCorr systematics for " << _s_category.at(i_cat) << endl;
		cout << _s_production_mode.at(Settings::qqToZZ)  << " : " << (_expected_yield_EWCorr_DN[Settings::qqToZZ][i_cat]/_expected_yield_EWCorr_DN[Settings::qqToZZ][Settings::inclusive])/(_expected_yield_EWCorr[Settings::qqToZZ][i_cat]/_expected_yield_EWCorr[Settings::qqToZZ][Settings::inclusive]) << "/" <<
		   (_expected_yield_EWCorr_UP[Settings::qqToZZ][i_cat]/_expected_yield_EWCorr_UP[Settings::qqToZZ][Settings::inclusive])/(_expected_yield_EWCorr[Settings::qqToZZ][i_cat]/_expected_yield_EWCorr[Settings::qqToZZ][Settings::inclusive]) << endl;
   }

	cout << "[INFO] EWCorr systematics for inclusive " << endl;
		cout << _s_production_mode.at(Settings::qqToZZ)  << " : " << (_expected_yield_EWCorr_DN[Settings::qqToZZ][Settings::inclusive]/_expected_yield_EWCorr[Settings::qqToZZ][Settings::inclusive]) << "/" <<
		   (_expected_yield_EWCorr_UP[Settings::qqToZZ][Settings::inclusive]/_expected_yield_EWCorr[Settings::qqToZZ][Settings::inclusive]) << endl;
	
}
//========================================

//========================================
void Systematics::PrintSystematics_PythiaScale( )
{
	
//	for ( int i_prod = 0; i_prod < num_of_production_modes; i_prod++ )
//   {
//   	for ( int i_cat = 0; i_cat < num_of_categories ; i_cat++ )
//		{
//			cout << _s_category.at(i_cat) << endl;
//			cout << "Nominal = " << _expected_yield_PythiaScale[i_prod][i_cat] << endl;
//			cout << "Up = " << _expected_yield_PythiaScale_UP[i_prod][i_cat] << endl;
//			cout << "Down = " << _expected_yield_PythiaScale_DN[i_prod][i_cat] << endl;
//			cout << endl;
//		}
//   }
	
	for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
   {
   	cout << "[INFO] Pythia Scale systematics for " << _s_category.at(i_cat) << endl;
		for ( int i_prod = 0; i_prod < num_of_production_modes; i_prod++ )
		{
			cout << _s_production_mode.at(i_prod)  << " : " << (_expected_yield_PythiaScale_DN[i_prod][i_cat]/_expected_yield_PythiaScale_DN[i_prod][Settings::inclusive])/(_expected_yield_PythiaScale[i_prod][i_cat]/_expected_yield_PythiaScale[i_prod][Settings::inclusive]) << "/" <<
		   (_expected_yield_PythiaScale_UP[i_prod][i_cat]/_expected_yield_PythiaScale_UP[i_prod][Settings::inclusive])/(_expected_yield_PythiaScale[i_prod][i_cat]/_expected_yield_PythiaScale[i_prod][Settings::inclusive]) << endl;
		}
   }

	cout << "[INFO] Pythia Scale systematics for inclusive: " << endl;
	for ( int i_prod = 0; i_prod < num_of_production_modes; i_prod++ )
	{
		cout << _s_production_mode.at(i_prod)  << " : " << (_expected_yield_PythiaScale_DN[i_prod][Settings::inclusive]/_expected_yield_PythiaScale[i_prod][Settings::inclusive]) << "/" <<
		(_expected_yield_PythiaScale_UP[i_prod][Settings::inclusive]/_expected_yield_PythiaScale[i_prod][Settings::inclusive]) << endl;
	}
}
//========================================

//========================================
void Systematics::PrintSystematics_PythiaTune( )
{
	
//	for ( int i_prod = 0; i_prod < num_of_production_modes; i_prod++ )
//   {
//   	for ( int i_cat = 0; i_cat < num_of_categories ; i_cat++ )
//		{
//			cout << _s_category.at(i_cat) << endl;
//			cout << "Nominal = " << _expected_yield_PythiaTune[i_prod][i_cat] << endl;
//			cout << "Up = " << _expected_yield_PythiaTune_UP[i_prod][i_cat] << endl;
//			cout << "Down = " << _expected_yield_PythiaTune_DN[i_prod][i_cat] << endl;
//			cout << endl;
//		}
//   }
	float Dn_ratio = 1.;
	float Up_ratio = 1.;
	float Nominal_ratio = 1.;
	
	for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
   {
   	cout << "[INFO] Pythia Tune systematics for " << _s_category.at(i_cat) << endl;
		for ( int i_prod = 0; i_prod < num_of_production_modes; i_prod++ )
		{
			//Deal with low statistics in tuneup/tunedown samples
			if( _expected_yield_PythiaTune_DN[i_prod][i_cat] == 0 ) Dn_ratio = 0.;
			else Dn_ratio = _expected_yield_PythiaTune_DN[i_prod][i_cat]/_expected_yield_PythiaTune_DN[i_prod][Settings::inclusive];
			if( _expected_yield_PythiaTune_UP[i_prod][i_cat] == 0 ) Up_ratio = 0.;
			else Up_ratio = _expected_yield_PythiaTune_UP[i_prod][i_cat]/_expected_yield_PythiaTune_UP[i_prod][Settings::inclusive];
			Nominal_ratio = _expected_yield_PythiaTune[i_prod][i_cat]/_expected_yield_PythiaTune[i_prod][Settings::inclusive];
			
			
			cout << _s_production_mode.at(i_prod)  << " : " << (Dn_ratio)/(Nominal_ratio) << "/" << (Up_ratio)/(Nominal_ratio) << endl;
		}
   }

}
//========================================

//========================================
void Systematics::PrintSystematics_QCDScale( )
{
	float muR_Up = 1.;
	float muF_Up = 1.;
	float muR_Dn = 1.;
	float muF_Dn = 1.;
	
	float mu_Up  = 1.;
	float mu_Dn  = 1.;
	
	float PythiaScale_Up = 1.;
	float PythiaTune_Up  = 1.;
	float PythiaScale_Dn = 1.;
	float PythiaTune_Dn  = 1.;
	
	float Pythia_Up = 1.;
	float Pythia_Dn = 1.;
	
	float comb_Up = 1.;
	float comb_Dn = 1.;
	
	int combination_mu     = 0;
	int combination_Pythia = 0;
	
	for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
   {
   	cout << "[INFO] QCD scale systematics for " << _s_category.at(i_cat) << endl;
		for ( int i_prod = 0; i_prod < num_of_production_modes; i_prod++ )
		{
			// Combine muR and muF
			muR_Up = (_expected_yield_muR_UP[i_prod][i_cat]/_expected_yield_muR_UP[i_prod][Settings::inclusive])/(_expected_yield_muR[i_prod][i_cat]/_expected_yield_muR[i_prod][Settings::inclusive]);
			muR_Dn = (_expected_yield_muR_DN[i_prod][i_cat]/_expected_yield_muR_DN[i_prod][Settings::inclusive])/(_expected_yield_muR[i_prod][i_cat]/_expected_yield_muR[i_prod][Settings::inclusive]);
			
			muF_Up = (_expected_yield_muF_UP[i_prod][i_cat]/_expected_yield_muF_UP[i_prod][Settings::inclusive])/(_expected_yield_muF[i_prod][i_cat]/_expected_yield_muF[i_prod][Settings::inclusive]);
			muF_Dn = (_expected_yield_muF_DN[i_prod][i_cat]/_expected_yield_muF_DN[i_prod][Settings::inclusive])/(_expected_yield_muF[i_prod][i_cat]/_expected_yield_muF[i_prod][Settings::inclusive]);
			
			if( i_cat == Settings::untagged && muR_Up > 1. && muF_Up > 1.) combination_mu = 0;
			if( i_cat == Settings::untagged && muR_Up > 1. && muF_Up < 1.) combination_mu = 1;
			if( i_cat == Settings::untagged && muR_Up < 1. && muF_Up > 1.) combination_mu = 2;
			if( i_cat == Settings::untagged && muR_Up < 1. && muF_Up < 1.) combination_mu = 3;
			
			if(combination_mu == 0) {mu_Dn = 1. - sqrt(pow(1. - muR_Dn,2)+pow(1. - muF_Dn,2));  mu_Up = 1. + sqrt(pow(muR_Up - 1.,2)+pow(muF_Up - 1.,2));}
			if(combination_mu == 1) {mu_Dn = 1. - sqrt(pow(1. - muR_Dn,2)+pow(1. - muF_Up,2));  mu_Up = 1. + sqrt(pow(muR_Up - 1.,2)+pow(muF_Dn - 1.,2));}
			if(combination_mu == 2) {mu_Dn = 1. - sqrt(pow(1. - muR_Up,2)+pow(1. - muF_Dn,2));  mu_Up = 1. + sqrt(pow(muR_Dn - 1.,2)+pow(muF_Up - 1.,2));}
			if(combination_mu == 3) {mu_Dn = 1. - sqrt(pow(1. - muR_Up,2)+pow(1. - muF_Up,2));  mu_Up = 1. + sqrt(pow(muR_Dn - 1.,2)+pow(muF_Dn - 1.,2));}
			
			
			// Combine Pythia tune and Pythia scale
			PythiaScale_Up = (_expected_yield_PythiaScale_UP[i_prod][i_cat]/_expected_yield_PythiaScale_UP[i_prod][Settings::inclusive])/(_expected_yield_PythiaScale[i_prod][i_cat]/_expected_yield_PythiaScale[i_prod][Settings::inclusive]);
			if(_expected_yield_PythiaScale_UP[i_prod][i_cat] == 0) PythiaScale_Up = 1.;
			PythiaScale_Dn = (_expected_yield_PythiaScale_DN[i_prod][i_cat]/_expected_yield_PythiaScale_DN[i_prod][Settings::inclusive])/(_expected_yield_PythiaScale[i_prod][i_cat]/_expected_yield_PythiaScale[i_prod][Settings::inclusive]);
			if(_expected_yield_PythiaScale_DN[i_prod][i_cat] == 0) PythiaScale_Dn = 1.;
			
			PythiaTune_Up = (_expected_yield_PythiaTune_UP[i_prod][i_cat]/_expected_yield_PythiaTune_UP[i_prod][Settings::inclusive])/(_expected_yield_PythiaTune[i_prod][i_cat]/_expected_yield_PythiaTune[i_prod][Settings::inclusive]);
			if(PythiaTune_Up == 0) PythiaTune_Up = 1.;
			if(PythiaTune_Up > 2.0) PythiaTune_Up = 1.5;
			PythiaTune_Dn = (_expected_yield_PythiaTune_DN[i_prod][i_cat]/_expected_yield_PythiaTune_DN[i_prod][Settings::inclusive])/(_expected_yield_PythiaTune[i_prod][i_cat]/_expected_yield_PythiaTune[i_prod][Settings::inclusive]);
			if(PythiaTune_Dn== 0) PythiaTune_Dn = 1.;
			if(PythiaTune_Dn > 2.0) PythiaTune_Dn = 1.5;
			
			if( i_cat == Settings::untagged && PythiaScale_Up > 1. && PythiaTune_Up > 1.) combination_Pythia = 0;
			if( i_cat == Settings::untagged && PythiaScale_Up > 1. && PythiaTune_Up < 1.) combination_Pythia = 1;
			if( i_cat == Settings::untagged && PythiaScale_Up < 1. && PythiaTune_Up > 1.) combination_Pythia = 2;
			if( i_cat == Settings::untagged && PythiaScale_Up < 1. && PythiaTune_Up < 1.) combination_Pythia = 3;
			
			if(combination_Pythia == 0) {Pythia_Dn = 1. - sqrt(pow(1. - PythiaScale_Dn,2)+pow(1. - PythiaTune_Dn,2));  Pythia_Up = 1. + sqrt(pow(PythiaScale_Up - 1.,2)+pow(PythiaTune_Up - 1.,2));}
			if(combination_Pythia == 1) {Pythia_Dn = 1. - sqrt(pow(1. - PythiaScale_Dn,2)+pow(1. - PythiaTune_Up,2));  Pythia_Up = 1. + sqrt(pow(PythiaScale_Up - 1.,2)+pow(PythiaTune_Dn - 1.,2));}
			if(combination_Pythia == 2) {Pythia_Dn = 1. - sqrt(pow(1. - PythiaScale_Up,2)+pow(1. - PythiaTune_Dn,2));  Pythia_Up = 1. + sqrt(pow(PythiaScale_Dn - 1.,2)+pow(PythiaTune_Up - 1.,2));}
			if(combination_Pythia == 3) {Pythia_Dn = 1. - sqrt(pow(1. - PythiaScale_Up,2)+pow(1. - PythiaTune_Up,2));  Pythia_Up = 1. + sqrt(pow(PythiaScale_Dn - 1.,2)+pow(PythiaTune_Dn - 1.,2));}
			
			// Combine everything
			comb_Up = 1. + sqrt(pow(Pythia_Up - 1.,2)+pow(mu_Up - 1.,2));
			comb_Dn = 1. - sqrt(pow(1. - Pythia_Dn,2)+pow(1. - mu_Dn,2));

			cout << _s_production_mode.at(i_prod)  << " : " << comb_Dn << "/" << comb_Up << endl;
			
		}
   }
}
//========================================


//========================================
void Systematics::PrintSystematics_THU_ggH( )
{
	
//	for ( int i_prod = 0; i_prod < num_of_production_modes; i_prod++ )
//   {
//   	for ( int i_cat = 0; i_cat < num_of_categories ; i_cat++ )
//		{
//			cout << _s_category.at(i_cat) << endl;
//			cout << "Nominal = " << _expected_yield_BTag[i_prod][i_cat] << endl;
//			cout << "Up = " << _expected_yield_BTag_UP[i_prod][i_cat] << endl;
//			cout << "Down = " << _expected_yield_BTag_DN[i_prod][i_cat] << endl;
//			cout << endl;
//		}
//   }
	

	cout << "[INFO] THU_ggH_Mu systematics for " ;
		for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
   {
    	cout <<_s_category.at(i_cat) << endl;
		cout << _s_production_mode.at(Settings::ggH)  << " : " << (_expected_yield_THU_ggH_Mu_UP[Settings::ggH][i_cat]/_expected_yield_THU_ggH_Mu[Settings::ggH][i_cat]) << endl;
	}
	
	cout << "[INFO] THU_ggH_Res systematics for " ;
		for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
   {
   	cout <<_s_category.at(i_cat) << endl;
		cout << _s_production_mode.at(Settings::ggH)  << " : " << (_expected_yield_THU_ggH_Res_UP[Settings::ggH][i_cat]/_expected_yield_THU_ggH_Res[Settings::ggH][i_cat]) << endl;
	}
	
	cout << "[INFO] THU_ggH_Mig01 systematics for " ;
		for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
   {
   	cout <<_s_category.at(i_cat) << endl;
		cout << _s_production_mode.at(Settings::ggH)  << " : " << (_expected_yield_THU_ggH_Mig01_UP[Settings::ggH][i_cat]/_expected_yield_THU_ggH_Mig01_UP[Settings::ggH][Settings::inclusive])/(_expected_yield_THU_ggH_Mig01[Settings::ggH][i_cat]/_expected_yield_THU_ggH_Mig01[Settings::ggH][Settings::inclusive]) << endl;
	}
	
	cout << "[INFO] THU_ggH_Mig12 systematics for " ;
		for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
   {
   	cout <<_s_category.at(i_cat) << endl;
		cout << _s_production_mode.at(Settings::ggH)  << " : " << (_expected_yield_THU_ggH_Mig12_UP[Settings::ggH][i_cat]/_expected_yield_THU_ggH_Mig12_UP[Settings::ggH][Settings::inclusive])/(_expected_yield_THU_ggH_Mig12[Settings::ggH][i_cat]/_expected_yield_THU_ggH_Mig12[Settings::ggH][Settings::inclusive]) << endl;
	}
	
	cout << "[INFO] THU_ggH_VBF2j systematics for " ;
		for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
   {
   	cout <<_s_category.at(i_cat) << endl;
		cout << _s_production_mode.at(Settings::ggH)  << " : " << (_expected_yield_THU_ggH_VBF2j_UP[Settings::ggH][i_cat]/_expected_yield_THU_ggH_VBF2j_UP[Settings::ggH][Settings::inclusive])/(_expected_yield_THU_ggH_VBF2j[Settings::ggH][i_cat]/_expected_yield_THU_ggH_VBF2j[Settings::ggH][Settings::inclusive]) << endl;
	}
	
	cout << "[INFO] THU_ggH_VBF3j systematics for " ;
		for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
   {
   	cout <<_s_category.at(i_cat) << endl;
		cout << _s_production_mode.at(Settings::ggH)  << " : " << (_expected_yield_THU_ggH_VBF3j_UP[Settings::ggH][i_cat]/_expected_yield_THU_ggH_VBF3j_UP[Settings::ggH][Settings::inclusive])/(_expected_yield_THU_ggH_VBF3j[Settings::ggH][i_cat]/_expected_yield_THU_ggH_VBF3j[Settings::ggH][Settings::inclusive]) << endl;
	}
	
	cout << "[INFO] THU_ggH_PT60 systematics for " ;
		for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
   {
   	cout <<_s_category.at(i_cat) << endl;
		cout << _s_production_mode.at(Settings::ggH)  << " : " << (_expected_yield_THU_ggH_PT60_UP[Settings::ggH][i_cat]/_expected_yield_THU_ggH_PT60_UP[Settings::ggH][Settings::inclusive])/(_expected_yield_THU_ggH_PT60[Settings::ggH][i_cat]/_expected_yield_THU_ggH_PT60[Settings::ggH][Settings::inclusive]) << endl;
	}
	
	cout << "[INFO] THU_ggH_PT120 systematics for " ;
		for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
   {
   	cout <<_s_category.at(i_cat) << endl;
		cout << _s_production_mode.at(Settings::ggH)  << " : " << (_expected_yield_THU_ggH_PT120_UP[Settings::ggH][i_cat]/_expected_yield_THU_ggH_PT120_UP[Settings::ggH][Settings::inclusive])/(_expected_yield_THU_ggH_PT120[Settings::ggH][i_cat]/_expected_yield_THU_ggH_PT120[Settings::ggH][Settings::inclusive]) << endl;
	}
	
	cout << "[INFO] THU_ggH_qmtop systematics for " ;
		for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
   {
   	cout <<_s_category.at(i_cat) << endl;
		cout << _s_production_mode.at(Settings::ggH)  << " : " << (_expected_yield_THU_ggH_qmtop_UP[Settings::ggH][i_cat]/_expected_yield_THU_ggH_qmtop_UP[Settings::ggH][Settings::inclusive])/(_expected_yield_THU_ggH_qmtop[Settings::ggH][i_cat]/_expected_yield_THU_ggH_qmtop[Settings::ggH][Settings::inclusive]) << endl;
	}
	


	
}
//========================================


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
   if ( input_file_name.Contains("ggTo") )          							current_production_mode = Settings::ggToZZ;
   if ( input_file_name.Contains("ZZTo") )         							current_production_mode = Settings::qqToZZ;

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

