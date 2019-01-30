#include <ZZAnalysis/AnalysisStep/interface/Category.h>
#include <ZZAnalysis/AnalysisStep/interface/Discriminants.h>
#include <ZZAnalysis/AnalysisStep/interface/cConstants.h>

#include <cmath>
#include <iostream>

#include "TMath.h"
#include "TRandom3.h"
#include "TH1F.h"

using namespace std;

#define VERBOSE 1


float bins_hpt4[]={0,60,120,200};
TH1F *hpt_bin=new TH1F("hpt_bin","",3, bins_hpt4);

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
					 float p_HadWH_mavjj_JECNominal,
					 float p_HadWH_mavjj_true_JECNominal,
					 float p_HadZH_mavjj_JECNominal,
					 float p_HadZH_mavjj_true_JECNominal,
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
      D_WHh   = DWHh_ME_QG(p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_HadWH_mavjj_JECNominal, p_HadWH_mavjj_true_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
      D_ZHh   = DZHh_ME_QG(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_HadZH_mavjj_JECNominal, p_HadZH_mavjj_true_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
    }
  }else{
    if(nCleanedJetsPt30==1){
      D_VBF1j = DVBF1j_ME(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
    }else if(nCleanedJetsPt30>=2){
      D_VBF2j = DVBF2j_ME(p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
      D_WHh   = DWHh_ME(p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_HadWH_mavjj_JECNominal, p_HadWH_mavjj_true_JECNominal, ZZMass);
      D_ZHh   = DZHh_ME(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_HadZH_mavjj_JECNominal, p_HadZH_mavjj_true_JECNominal, ZZMass);
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
				  float p_HadWH_mavjj_JECNominal,
				  float p_HadWH_mavjj_true_JECNominal,
				  float p_HadZH_mavjj_JECNominal,
				  float p_HadZH_mavjj_true_JECNominal,
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
      D_WHh   = DWHh_ME_QG(p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_HadWH_mavjj_JECNominal, p_HadWH_mavjj_true_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
      D_ZHh   = DZHh_ME_QG(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_HadZH_mavjj_JECNominal, p_HadZH_mavjj_true_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
    }
  }else{
    if(nCleanedJetsPt30==1){
      D_VBF1j = DVBF1j_ME(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
    }else if(nCleanedJetsPt30>=2){
      D_VBF2j = DVBF2j_ME(p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
      D_WHh   = DWHh_ME(p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_HadWH_mavjj_JECNominal, p_HadWH_mavjj_true_JECNominal, ZZMass);
      D_ZHh   = DZHh_ME(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_HadZH_mavjj_JECNominal, p_HadZH_mavjj_true_JECNominal, ZZMass);
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
      D_WHh   = DWHh_ME_QG(p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_HadWH_mavjj_JECNominal, p_HadWH_mavjj_true_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
      D_ZHh   = DZHh_ME_QG(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_HadZH_mavjj_JECNominal, p_HadZH_mavjj_true_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
    }
  }else{
    if(nCleanedJetsPt30==1){
      D_VBF1j = DVBF1j_ME(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
    }else if(nCleanedJetsPt30>=2){
      D_VBF2j = DVBF2j_ME(p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
      D_WHh   = DWHh_ME(p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_HadWH_mavjj_JECNominal, p_HadWH_mavjj_true_JECNominal, ZZMass);
      D_ZHh   = DZHh_ME(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_HadZH_mavjj_JECNominal, p_HadZH_mavjj_true_JECNominal, ZZMass);
    }
  }

  float WP_VBF2j = getDVBF2jetsWP(ZZMass, useQGTagging);
  float WP_VBF1j = getDVBF1jetWP(ZZMass, useQGTagging);
  float WP_WHh = getDWHhWP(ZZMass, useQGTagging);
  float WP_ZHh = getDZHhWP(ZZMass, useQGTagging);

  if( nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && D_VBF2j>WP_VBF2j ){

    return VBF2jTaggedMor18;

  }else if( nExtraLep==0 && (nCleanedJetsPt30==2||nCleanedJetsPt30==3||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && (D_WHh>WP_WHh||D_ZHh>WP_ZHh) ){

    return VHHadrTaggedMor18;

  }else if(   ( nCleanedJetsPt30<=3 && nCleanedJetsPt30BTagged_bTagSF==0 && (nExtraLep==1||nExtraZ>=1) )
	   || ( nCleanedJetsPt30==0 && nExtraLep>=1 ) ){

    return VHLeptTaggedMor18;

  }else if( nCleanedJetsPt30>=4 && nCleanedJetsPt30BTagged_bTagSF>=1 && nExtraLep ==0){

    return ttHHadrTaggedMor18;
  
  }else if( nExtraLep>=1 ){
  
  	return ttHLeptTaggedMor18;
	
  }else if( useVHMETTagged && nExtraLep==0 && (nCleanedJetsPt30==0||nCleanedJetsPt30==1) && PFMET>100 ){

    return VHMETTaggedMor18;

  }else if( nExtraLep==0 && nCleanedJetsPt30==1 && D_VBF1j>WP_VBF1j ){

    return VBF1jTaggedMor18;

  }else{

    return UntaggedMor18;

  }

}
// inputs: Njets, mjj, H_pt, category number returned from categoryMor18 without MET, pt_Hjj, pass a string argument in the last so that it returns you the name of the category 
extern "C" int stage1_reco_1p1_sync(int Njets, float mjj, float H_pt, int category, float pt_hjj,string &reco_catName){
	int vbfTopo=0;
	if (Njets<2) vbfTopo=0; 
	vbfTopo = mjj > 350.0; 
	if(category>4 ){
		switch (category){
			case 6: reco_catName = "TTH_Had"; break; 
			case 5: reco_catName = "TTH_Lep"; break; 
		}
		return  15+category; 
	}
	else if(category==3){
		if(H_pt<150 ){reco_catName = "VH_Lep_0_150"; return 18;}
		else {reco_catName = "VH_Lep_GT150"; return 19;}
	}
	else if (category ==1){
		reco_catName = "VBF_1j"; return 10;
	}
	else if(category==2){
		if(vbfTopo){
			if (H_pt>200 )   {
				reco_catName = "VBF_GT200_2J"; return 15;
			}
			else{
				if (pt_hjj>25)
				{reco_catName = "VBF_3j"; return 14;}
				else{
					if (mjj > 350 && mjj < 700 ){
						reco_catName = "VBF_2j_mjj_350_700_2j"; return 12;
					}
					else if (mjj > 700 ){
						reco_catName = "VBF_2j_mjj_GT700_2j"; return 13;
					}
				}
			}
		}
		else {
			reco_catName = "VBF_2j"; return 11;
		}
	}
	else if (category == 4){
		if ( 60 < mjj && mjj < 120)
		{reco_catName = "VH_Had"; return 16;}
		else
		{reco_catName = "VBF_rest_VH"; return 17;}
	}
	else {
		if (H_pt>200 )   {
			reco_catName = "ggH_PTH_200"; return 8;
		}
		else{
			if (Njets==0)        {
				if(H_pt<10)
				{reco_catName = "ggH_0j_0_10"; return 0;}
				else
				{reco_catName = "ggH_0j_10_200"; return 1;}
			}
			else if (Njets==1)   {
				int binpt = hpt_bin->FindBin(H_pt);
				switch(binpt){
					case 1: reco_catName = "ggH_1j_0_60"; break;
					case 2: reco_catName = "ggH_1j_60_120"; break;
					case 3: reco_catName = "ggH_1j_120_200"; break;
				}
				return binpt+1;
			} 
			else if ( Njets>=2) {
				if(vbfTopo) {
					reco_catName = "ggH_VBF";return 9;
				}
				int binpt = hpt_bin->FindBin(H_pt);
				switch(binpt){
					case 1: reco_catName = "ggH_2j_0_60"; break;
					case 2: reco_catName = "ggH_2j_60_120"; break;
					case 3: reco_catName = "ggH_2j_120_200"; break;
				}
				return binpt+4;
			}
		}
	}
	return -1;
}
// stage1p1_reco:
// ggH_0J_PTH_0_10 =0
// ggH_0J_PTH_10_200 =1
// ggH_1J_PTH_0_60 =2
// ggH_1J_PTH_60_120 =3
// ggH_1J_PTH_120_200=4
// ggH_2J_PTH_0_60 =5
// ggH_2J_PTH_60_120 =6
// ggH_2J_PTH_120_200=7
// ggH_PTH_200=8
// ggH_VBF=9

// VBF_1j = 10;
// VBF_2j= 11;
// VBF_2j_mjj_350_700_2j = 12;
// VBF_2j_mjj_GT700_2j = 13;
// VBF_2j_mjj_GT350_3j = 14;
// VBF_GT200_2J = 15;

// VH_Had = 16;
// VBF_rest_VH =17;

// VH_lep_0_150 =18;
// VH_Lep_GT150=19;

// ttH_Lep = 20 ;
// ttH_Hap = 21;
