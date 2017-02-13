#ifndef DISCRIMINANTS_H
#define DISCRIMINANTS_H


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
    float ZZMass
			 );
extern "C" float DZHh_ME(
    float p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
    float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
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
    float ZZMass,
    float* jetQGLikelihood,
    float* jetPhi
			    );
extern "C" float DZHh_ME_QG(
    float p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
    float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
    float ZZMass,
    float* jetQGLikelihood,
    float* jetPhi
			    );


#endif
