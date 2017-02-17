#include <ZZAnalysis/AnalysisStep/interface/Category.h>
#include <ZZAnalysis/AnalysisStep/interface/Discriminants.h>
#include <ZZAnalysis/AnalysisStep/interface/cConstants.h>

#include <cmath>
#include <iostream>

#include "TMath.h"
#include "TRandom3.h"

using namespace std;

#define VERBOSE 1



extern "C" int categoryLegacy( int nCleanedJetsPt30 )
{
  if(VERBOSE) cout << "WARNING: using deprecated categorization function 'categoryLegacy'" << endl;

  if(nCleanedJetsPt30>=2)
    return Dijet;
  else
    return ZeroOneJet;
}


extern "C" int categoryMor16(
			     int nCleanedJetsPt30,
			     float p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
			     float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal
			     )
{
  if(VERBOSE) cout << "WARNING: using deprecated categorization function 'categoryMor16'" << endl;

  float vbfMela = p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal / ( p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal + p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal );

  if(nCleanedJetsPt30>=2 && vbfMela>0.5)
    return VBFTaggedMor16;
  else
    return UntaggedMor16;
}


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
			       float* jetPhi,
			       float ZZMass,
			       bool useQGTagging 
			       )
{
  if(VERBOSE) cout << "WARNING: using deprecated categorization function 'categoryIchep16'" << endl;

  float D_VBF2j = -2;
  float D_VBF1j = -2;
  float D_WHh   = -2;
  float D_ZHh   = -2;
  if(useQGTagging){
    if(nCleanedJetsPt30==1){
      D_VBF1j = DVBF1j_ME_QG(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
    }else if(nCleanedJetsPt30>=2){
      D_VBF2j = DVBF2j_ME_QG(p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
      D_WHh   = DWHh_ME_QG(p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
      D_ZHh   = DZHh_ME_QG(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
    }
  }else{
    if(nCleanedJetsPt30==1){
      D_VBF1j = DVBF1j_ME(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
    }else if(nCleanedJetsPt30>=2){
      D_VBF2j = DVBF2j_ME(p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
      D_WHh   = DWHh_ME(p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
      D_ZHh   = DZHh_ME(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
    }
  }

  float WP_VBF2j = getDVBF2jetsWP(ZZMass, useQGTagging);
  float WP_VBF1j = getDVBF1jetWP(ZZMass, useQGTagging);
  float WP_WHh = getDWHhWP(ZZMass, useQGTagging);
  float WP_ZHh = getDZHhWP(ZZMass, useQGTagging);

  if( nExtraLep==0 && nCleanedJetsPt30==1 && D_VBF1j>WP_VBF1j ){

    return VBF1jTaggedIchep16;

  }else if( nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && D_VBF2j>WP_VBF2j ){

    return VBF2jTaggedIchep16;

  }else if(   ( nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && (D_WHh>WP_WHh||D_ZHh>WP_ZHh) )
	   || ( nExtraLep==0 && (nCleanedJetsPt30==2||nCleanedJetsPt30==3) && nCleanedJetsPt30BTagged_bTagSF>=2 ) ){

    return VHHadrTaggedIchep16;

  }else if(   ( nCleanedJetsPt30<=3 && nCleanedJetsPt30BTagged_bTagSF==0 && (nExtraLep==1||nExtraZ>=1) )
	   || ( nCleanedJetsPt30==0 && nExtraLep>=1 ) ){

    return VHLeptTaggedIchep16;

  }else if(   ( nCleanedJetsPt30>=4 && nCleanedJetsPt30BTagged_bTagSF>=1 )
	   || nExtraLep>=1 ){

    return ttHTaggedIchep16;

  }else{

    return UntaggedIchep16;

  }

}


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
			     float* jetPhi,
			     float ZZMass,
			     float PFMET,
			     bool useVHMETTagged,
			     bool useQGTagging
			     )
{

  float D_VBF2j = -2;
  float D_VBF1j = -2;
  float D_WHh   = -2;
  float D_ZHh   = -2;
  if(useQGTagging){
    if(nCleanedJetsPt30==1){
      D_VBF1j = DVBF1j_ME_QG(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
    }else if(nCleanedJetsPt30>=2){
      D_VBF2j = DVBF2j_ME_QG(p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
      D_WHh   = DWHh_ME_QG(p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
      D_ZHh   = DZHh_ME_QG(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
    }
  }else{
    if(nCleanedJetsPt30==1){
      D_VBF1j = DVBF1j_ME(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
    }else if(nCleanedJetsPt30>=2){
      D_VBF2j = DVBF2j_ME(p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
      D_WHh   = DWHh_ME(p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
      D_ZHh   = DZHh_ME(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
    }
  }

  float WP_VBF2j = getDVBF2jetsWP(ZZMass, useQGTagging);
  float WP_VBF1j = getDVBF1jetWP(ZZMass, useQGTagging);
  float WP_WHh = getDWHhWP(ZZMass, useQGTagging);
  float WP_ZHh = getDZHhWP(ZZMass, useQGTagging);

  if( nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && D_VBF2j>WP_VBF2j ){

    return VBF2jTaggedMor17;

  }else if( nExtraLep==0 && (nCleanedJetsPt30==2||nCleanedJetsPt30==3||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && (D_WHh>WP_WHh||D_ZHh>WP_ZHh) ){

    return VHHadrTaggedMor17;

  }else if(   ( nCleanedJetsPt30<=3 && nCleanedJetsPt30BTagged_bTagSF==0 && (nExtraLep==1||nExtraZ>=1) )
	   || ( nCleanedJetsPt30==0 && nExtraLep>=1 ) ){

    return VHLeptTaggedMor17;

  }else if(   ( nCleanedJetsPt30>=4 && nCleanedJetsPt30BTagged_bTagSF>=1 )
	   || nExtraLep>=1 ){

    return ttHTaggedMor17;
  
  }else if( useVHMETTagged && nExtraLep==0 && (nCleanedJetsPt30==0||nCleanedJetsPt30==1) && PFMET>100 ){

    return VHMETTaggedMor17;

  }else if( nExtraLep==0 && nCleanedJetsPt30==1 && D_VBF1j>WP_VBF1j ){

    return VBF1jTaggedMor17;

  }else{

    return UntaggedMor17;

  }

}
