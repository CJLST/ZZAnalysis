#include <ZZAnalysis/AnalysisStep/interface/Discriminants.h>
#include <ZZAnalysis/AnalysisStep/interface/cConstants.h>

#include <cmath>

#include "TMath.h"
#include "TRandom3.h"


extern "C" float D_bkg_kin(
    float p_GG_SIG_ghg2_1_ghz1_1_JHUGen,
    float p_QQB_BKG_MCFM,
    int   ZZFlav,
    float ZZMass) 
{
  return p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(p_GG_SIG_ghg2_1_ghz1_1_JHUGen + p_QQB_BKG_MCFM*getDbkgkinConstant(ZZFlav,ZZMass)); 
}


extern "C" float D_bkg(
    float p_GG_SIG_ghg2_1_ghz1_1_JHUGen,
    float p_m4l_SIG,
    float p_QQB_BKG_MCFM,
    float p_m4l_BKG,
    int   ZZFlav,
    float ZZMass)
{
  return p_GG_SIG_ghg2_1_ghz1_1_JHUGen*p_m4l_SIG/(p_GG_SIG_ghg2_1_ghz1_1_JHUGen*p_m4l_SIG + p_QQB_BKG_MCFM*p_m4l_BKG*getDbkgConstant(ZZFlav,ZZMass));
}


extern "C" float D_g4( 
    float p_GG_SIG_ghg2_1_ghz1_1_JHUGen,
    float p_GG_SIG_ghg2_1_ghz4_1_JHUGen) 
{
  return p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(p_GG_SIG_ghg2_1_ghz1_1_JHUGen + pow(2.521, 2)*p_GG_SIG_ghg2_1_ghz4_1_JHUGen); //Note the hardcoded c-constant!
}




extern "C" float DVBF2j_ME(
    float p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
    float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
    float ZZMass 
			   )
{
  float c_Mela2j = getDVBF2jetsConstant(ZZMass);
  return 1./(1.+ c_Mela2j*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal);
}

extern "C" float DVBF1j_ME(
    float p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
    float pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
    float p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
    float ZZMass
			   )
{
  float c_Mela1j = getDVBF1jetConstant(ZZMass);
  return 1./(1.+ c_Mela1j*p_JQCD_SIG_ghg2_1_JHUGen_JECNominal/(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal*pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal));
}

extern "C" float DWHh_ME(
    float p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
    float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
	 float p_HadWH_mavjj_JECNominal,
    float ZZMass
			 )
{
  float c_MelaWH = getDWHhConstant(ZZMass);
  return 1./(1.+ c_MelaWH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/(p_HadWH_mavjj_JECNominal*p_HadWH_SIG_ghw1_1_JHUGen_JECNominal));
}

extern "C" float DZHh_ME(
    float p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
    float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
	 float p_HadZH_mavjj_JECNominal,
    float ZZMass
			 )
{
  float c_MelaZH = getDZHhConstant(ZZMass);
  return 1./(1.+ c_MelaZH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/(p_HadZH_mavjj_JECNominal*p_HadZH_SIG_ghz1_1_JHUGen_JECNominal));
}

float jetPgOverPq(float jetQGLikelihood, float jetPhi)
{
  if(jetQGLikelihood<0.){
    TRandom3 rand;
    rand.SetSeed(abs(static_cast<int>(sin(jetPhi)*100000)));
    return 1./rand.Uniform() - 1.;
  }else{    
    return 1./jetQGLikelihood - 1.;
  }
}

extern "C" float DVBF2j_ME_QG(
    float p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
    float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
    float ZZMass,
    float* jetQGLikelihood,
    float* jetPhi
			      )
{
  float DVBF2jME = DVBF2j_ME(p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
  float GOverQ = TMath::Power( jetPgOverPq(jetQGLikelihood[0],jetPhi[0]) * jetPgOverPq(jetQGLikelihood[1],jetPhi[1]) , 1./3. );
  return 1./(1.+ (1./DVBF2jME - 1.) * GOverQ);
}

extern "C" float DVBF1j_ME_QG(
    float p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
    float pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
    float p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
    float ZZMass,
    float* jetQGLikelihood,
    float* jetPhi
			      )
{
  float DVBF1jME = DVBF1j_ME(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
  float GOverQ = TMath::Power( jetPgOverPq(jetQGLikelihood[0],jetPhi[0]) , 1./3. );
  return 1./(1.+ (1./DVBF1jME - 1.) * GOverQ);
}

extern "C" float DWHh_ME_QG(
    float p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
    float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
	 float p_HadWH_mavjj_JECNominal,
    float ZZMass,
    float* jetQGLikelihood,
    float* jetPhi
			    )
{
  float DWHhME = DWHh_ME(p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_HadWH_mavjj_JECNominal, ZZMass);
  float GOverQ = TMath::Power( jetPgOverPq(jetQGLikelihood[0],jetPhi[0]) * jetPgOverPq(jetQGLikelihood[1],jetPhi[1]) , 1./3. );
  return 1./(1.+ (1./DWHhME - 1.) * GOverQ);
}

extern "C" float DZHh_ME_QG(
    float p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
    float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
	 float p_HadZH_mavjj_JECNominal,
    float ZZMass,
    float* jetQGLikelihood,
    float* jetPhi
			    )
{
  float DZHhME = DZHh_ME(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_HadZH_mavjj_JECNominal, ZZMass);
  float GOverQ = TMath::Power( jetPgOverPq(jetQGLikelihood[0],jetPhi[0]) * jetPgOverPq(jetQGLikelihood[1],jetPhi[1]) , 1./3. );
  return 1./(1.+ (1./DZHhME - 1.) * GOverQ);
}

