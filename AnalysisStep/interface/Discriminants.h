#ifndef DISCRIMINANTS_H
#define DISCRIMINANTS_H



extern "C" float D_bkg_kin(
    float p_GG_SIG_ghg2_1_ghz1_1_JHUGen,
    float p_QQB_BKG_MCFM,
    int   ZZflav,
    float ZZMass);

extern "C" float D_bkg(
    float p_GG_SIG_ghg2_1_ghz1_1_JHUGen,
    float p_m4l_SIG,
    float p_QQB_BKG_MCFM,
    float p_m4l_BKG,
    int   ZZflav,
    float ZZMass);

extern "C" float D_bkg_VBFdec(
    float p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal,
    float p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal,
    float p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal,
    float p_JJVBF_BKG_MCFM_JECNominal,
    float p_HadZH_BKG_MCFM_JECNominal,
    float p_HadWH_BKG_MCFM_JECNominal,
    float p_JJQCD_BKG_MCFM_JECNominal,
    float p_HadZH_mavjj_JECNominal,
    float p_HadZH_mavjj_true_JECNominal,
    float p_HadWH_mavjj_JECNominal,
    float p_HadWH_mavjj_true_JECNominal,
    float pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal,
    float pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal,
    float pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal,
    float pConst_JJVBF_BKG_MCFM_JECNominal,
    float pConst_HadZH_BKG_MCFM_JECNominal,
    float pConst_HadWH_BKG_MCFM_JECNominal,
    float pConst_JJQCD_BKG_MCFM_JECNominal,
	 int   ZZFlav,
    float ZZMass);


extern "C" float D_bkg_VHdec(
    float p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal,
    float p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal,
    float p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal,
    float p_JJVBF_BKG_MCFM_JECNominal,
    float p_HadZH_BKG_MCFM_JECNominal,
    float p_HadWH_BKG_MCFM_JECNominal,
    float p_JJQCD_BKG_MCFM_JECNominal,
    float p_HadZH_mavjj_JECNominal,
    float p_HadZH_mavjj_true_JECNominal,
    float p_HadWH_mavjj_JECNominal,
    float p_HadWH_mavjj_true_JECNominal,
    float pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal,
    float pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal,
    float pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal,
    float pConst_JJVBF_BKG_MCFM_JECNominal,
    float pConst_HadZH_BKG_MCFM_JECNominal,
    float pConst_HadWH_BKG_MCFM_JECNominal,
    float pConst_JJQCD_BKG_MCFM_JECNominal,
	 int   ZZFlav,
    float ZZMass);

// D0-
extern "C" float D_g4(
    float p_GG_SIG_ghg2_1_ghz1_1_JHUGen,
    float p_GG_SIG_ghg2_1_ghz4_1_JHUGen);


// matrix-element only discriminants

extern "C" float DVBF2j_ME(
    float p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
    float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
    float ZZMass 
			   );
extern "C" float DVBF1j_ME(
    float p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
    float pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
    float p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
    float ZZMass
			   );
extern "C" float DWHh_ME(
    float p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
    float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
	 float p_HadWH_mavjj_JECNominal,
	 float p_HadWH_mavjj_true_JECNominal,
    float ZZMass
			 );
extern "C" float DZHh_ME(
    float p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
    float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
	 float p_HadZH_mavjj_JECNominal,
	 float p_HadZH_mavjj_true_JECNominal,
    float ZZMass
			 );


// discriminants using matrix elements and q/g tagging

float jetPgOverPq(float jetQGLikelihood, float jetPhi);

extern "C" float DVBF2j_ME_QG(
    float p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
    float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
    float ZZMass,
    float* jetQGLikelihood,
    float* jetPhi
			      );
extern "C" float DVBF1j_ME_QG(
    float p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
    float pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
    float p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
    float ZZMass,
    float* jetQGLikelihood,
    float* jetPhi
			      );
extern "C" float DWHh_ME_QG(
    float p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
    float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
	 float p_HadWH_mavjj_JECNominal,
	 float p_HadWH_mavjj_true_JECNominal,
    float ZZMass,
    float* jetQGLikelihood,
    float* jetPhi
			    );
extern "C" float DZHh_ME_QG(
    float p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
    float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
	 float p_HadZH_mavjj_JECNominal,
	 float p_HadWH_mavjj_true_JECNominal,
    float ZZMass,
    float* jetQGLikelihood,
    float* jetPhi
			    );


#endif
