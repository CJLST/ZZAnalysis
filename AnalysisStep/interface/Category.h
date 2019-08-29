#ifndef CATEGORY_H
#define CATEGORY_H



//---------- RunI categorization 

enum CategoryLegacy {
  ZeroOneJet = 0,
  Dijet      = 1
};

extern "C" int categoryLegacy( int nCleanedJetsPt30 );



//---------- Moriond 2016 categorization 

enum CategoryMor16 {
  UntaggedMor16  = 0,
  VBFTaggedMor16 = 1
};

extern "C" int categoryMor16(
			     int nCleanedJetsPt30,
			     float p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
			     float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal
			     );



//---------- ICHEP 2016 categorization

enum CategoryIchep16 {
  UntaggedIchep16     = 0,
  VBF1jTaggedIchep16  = 1,
  VBF2jTaggedIchep16  = 2, 
  VHLeptTaggedIchep16 = 3, 
  VHHadrTaggedIchep16 = 4, 
  ttHTaggedIchep16    = 5
};

extern "C" int categoryIchep16(
			       int nExtraLep,
			       int nExtraZ,
			       int nCleanedJetsPt30, 
			       int nCleanedJetsPt30BTagged_bTagSF,
			       float* jetQGLikelihood,
			       float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
			       float p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
			       float p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
			       float p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
			       float pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
			       float p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
			       float p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
					 float p_HadWH_mavjj_JECNominal,
					 float p_HadWH_mavjj_true_JECNominal,
					 float p_HadZH_mavjj_JECNominal,
					 float p_HadZH_mavjj_true_JECNominal,
			       float* jetPhi,
			       float ZZMass,
			       bool useQGTagging = false
			       );



//---------- Moriond 2017 categorization

enum CategoryMor17 {
  UntaggedMor17     = 0,
  VBF1jTaggedMor17  = 1,
  VBF2jTaggedMor17  = 2, 
  VHLeptTaggedMor17 = 3, 
  VHHadrTaggedMor17 = 4,
  ttHTaggedMor17    = 5,
  VHMETTaggedMor17  = 6
};

extern "C" int categoryMor17(
			     int nExtraLep,
			     int nExtraZ,
			     int nCleanedJetsPt30, 
			     int nCleanedJetsPt30BTagged_bTagSF,
			     float* jetQGLikelihood,
			     float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
			     float p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
			     float p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
			     float p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
			     float pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
			     float p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
			     float p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
				  float p_HadWH_mavjj_JECNominal,
				  float p_HadWH_mavjj_true_JECNominal,
				  float p_HadZH_mavjj_JECNominal,
				  float p_HadZH_mavjj_true_JECNominal,
			     float* jetPhi,
			     float ZZMass,
			     float PFMET,
			     bool useVHMETTagged = true,
			     bool useQGTagging = false
			     );

//---------- Moriond 2018 categorization

enum CategoryMor18 {
  UntaggedMor18      = 0,
  VBF1jTaggedMor18   = 1,
  VBF2jTaggedMor18   = 2,
  VHLeptTaggedMor18  = 3,
  VHHadrTaggedMor18  = 4,
  ttHLeptTaggedMor18 = 5,
  ttHHadrTaggedMor18 = 6,
  VHMETTaggedMor18   = 7
};

extern "C" int categoryMor18(
			     int nExtraLep,
			     int nExtraZ,
			     int nCleanedJetsPt30,
			     int nCleanedJetsPt30BTagged_bTagSF,
			     float* jetQGLikelihood,
			     float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
			     float p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
			     float p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
			     float p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
			     float pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
			     float p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
			     float p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
				  float p_HadWH_mavjj_JECNominal,
				  float p_HadWH_mavjj_true_JECNominal,
				  float p_HadZH_mavjj_JECNominal,
				  float p_HadZH_mavjj_true_JECNominal,
			     float* jetPhi,
			     float ZZMass,
			     float PFMET,
			     bool useVHMETTagged = true,
			     bool useQGTagging = false
			     );


//---------- Run II Legacy categorization

enum CategoryLegacyRunII {
  ggH_0J_PTH_0_10     = 0,
  ggH_0J_PTH_10_200   = 1,
  ggH_1J_PTH_0_60     = 2,
  ggH_1J_PTH_60_120   = 3,
  ggH_1J_PTH_120_200  = 4,
  ggH_2J_PTH_0_60     = 5,
  ggH_2J_PTH_60_120   = 6,
  ggH_2J_PTH_120_200  = 7,
  ggH_PTH_200         = 8,
  ggH_VBF             = 9,
  VBF_1j              = 10,
  VBF_2j              = 11,
  VBF_2j_mjj_350_700_2j = 12,
  VBF_2j_mjj_GT700_2j   = 13,
  VBF_2j_mjj_GT350_3j   = 14,
  VBF_GT200_2J          = 15,
  VH_Had                = 16,
  VBF_rest_VH           = 17,
  VH_lep_0_150          = 18,
  VH_Lep_GT150          = 19,
  ttH_Lep               = 20,
  ttH_Had               = 21
};


extern "C" int stage1_reco_1p1(
                     int Njets,
                     float mjj,
                     float H_pt,
                     int category,
                     float pt_hjj
                     );


//---------- Anomalous couplings 2019 categorization
//VBF1j, VH Lep and MET and ttH are merged into untagged or boosted,
//but keep the numbers for compatibility with Mor18

enum CategoryAC19 {
  UntaggedAC19      = UntaggedMor18,
  VBF1jTaggedAC19   = VBF1jTaggedMor18,
  VBF2jTaggedAC19   = VBF2jTaggedMor18,
  VHLeptTaggedAC19  = VHLeptTaggedMor18,
  VHHadrTaggedAC19  = VHHadrTaggedMor18,
  ttHLeptTaggedAC19 = ttHLeptTaggedMor18,
  ttHHadrTaggedAC19 = ttHHadrTaggedMor18,
  VHMETTaggedAC19   = VHMETTaggedMor18,
  BoostedAC19       = 8
};

extern "C" int categoryAC19(
			     int nExtraLep,
			     int nExtraZ,
			     int nCleanedJetsPt30,
			     int nCleanedJetsPt30BTagged_bTagSF,
			     float* jetQGLikelihood,
			     float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
			     float p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
			     float p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
			     float p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
			     float pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
			     float p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
			     float p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
				  float p_HadWH_mavjj_JECNominal,
				  float p_HadWH_mavjj_true_JECNominal,
				  float p_HadZH_mavjj_JECNominal,
				  float p_HadZH_mavjj_true_JECNominal,
			     float* jetPhi,
			     float ZZMass,
			     float ZZPt,
			     float PFMET,
			     bool useVHMETTagged = true,
			     bool useQGTagging = false
			     );


#endif
