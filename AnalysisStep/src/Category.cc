#include <ZZAnalysis/AnalysisStep/interface/Category.h>
#include <ZZAnalysis/AnalysisStep/interface/cConstants.h>

#include <cmath>

#include "TMath.h"
#include "TRandom3.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >  LV;


int flagDijetVH(
		int nCleanedJetsPt30, 
		float* jetpt,
		float* jeteta,
		float* jetphi,
		float* jetmass
		)
{

  bool found = false;

  if(nCleanedJetsPt30>=2){

    for(int j1=0; j1<nCleanedJetsPt30; j1++){
      if( std::abs(jeteta[j1])<2.4 && jetpt[j1]>40. ){

	for(int j2=j1+1; j2<nCleanedJetsPt30; j2++){
	  if( std::abs(jeteta[j2])<2.4 && jetpt[j2]>40. ){

	    LV jet1 (jetpt[j1],jeteta[j1],jetphi[j1],jetmass[j1]);
	    LV jet2 (jetpt[j2],jeteta[j2],jetphi[j2],jetmass[j2]);
	    float mjj = (jet1+jet2).mass();

	    if( 60.<mjj && mjj<120. ){
	      found = true;
	      break;
	    }

	  }
	}

	if(found) break;
      }
    }

  }
  
  return found;

}


//int category(
extern "C" int category(
	     int nExtraLep,
	     float ZZPt,
	     float ZZMass,
	     int nCleanedJetsPt30, 
	     int nCleanedJetsPt30BTagged,
	     float* jetpt,
	     float* jeteta,
	     float* jetphi,
	     float* jetmass,
	     float DiJetFisher
	     )
{

  int category = -1;

  if( nExtraLep==0 && nCleanedJetsPt30>=2 && nCleanedJetsPt30BTagged<=1 && DiJetFisher>0.5 ){

    category = VBFTagged; 

  }else if( ( nExtraLep==0 && nCleanedJetsPt30>=2 && ZZPt>ZZMass && flagDijetVH(nCleanedJetsPt30,jetpt,jeteta,jetphi,jetmass) )
            || ( nExtraLep==0 && nCleanedJetsPt30==2 && nCleanedJetsPt30BTagged==2 ) ){

    category = VHHadrTagged;

  }else if( nExtraLep>=1 && nCleanedJetsPt30<=2 && nCleanedJetsPt30BTagged==0 ){

    category = VHLeptTagged;

  }else if( nExtraLep>=1 || (nCleanedJetsPt30>=3 && nCleanedJetsPt30BTagged>=1) ){

    category = ttHTagged;

  }else if(nCleanedJetsPt30>=1){

    category = OneJetTagged;

  }else{

    category = Untagged;

  }

  return category;

}


extern "C" int categoryMor16(
			     int nCleanedJetsPt30,
			     float p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
			     float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal
			     )
{
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
	     int nCleanedJetsPt30BTagged,
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

  float WP_VBF2j = getDVBF2jetsWP(ZZMass, useQGTagging);
  float WP_VBF1j = getDVBF1jetWP(ZZMass, useQGTagging);
  float WP_WHh = getDWHhWP(ZZMass, useQGTagging);
  float WP_ZHh = getDZHhWP(ZZMass, useQGTagging);

  float c_Mela2j = getDVBF2jetsConstant(ZZMass);
  float c_Mela1j = getDVBF1jetConstant(ZZMass);
  float c_MelaWH = getDWHhConstant(ZZMass);
  float c_MelaZH = getDZHhConstant(ZZMass);

  float jetPgOverPq[nCleanedJetsPt30];
  for(int j=0; j<nCleanedJetsPt30; j++){
    if(jetQGLikelihood[j]<0. && j<2){
      TRandom3 rand;
      rand.SetSeed(abs(static_cast<int>(sin(jetPhi[j])*100000)));
      jetPgOverPq[j] = 1./rand.Uniform() - 1.;
    }else{    
      jetPgOverPq[j] = 1./jetQGLikelihood[j] - 1.;
    }
  }

  float D_VBF2j = (nCleanedJetsPt30>=2) ? 1/(1+ c_Mela2j*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal * ( useQGTagging ? TMath::Power(jetPgOverPq[0]*jetPgOverPq[1],1/3.) : 1. ) ) : -2 ;
  float D_VBF1j = (nCleanedJetsPt30>=1) ? 1/(1+ c_Mela1j*p_JQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JVBF_SIG_ghv1_1_JHUGen_JECNominal*pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal * (useQGTagging ? TMath::Power(jetPgOverPq[0],1/3.) : 1. ) ) : -2 ;
  float D_WHh = (nCleanedJetsPt30>=2) ? 1/(1+ c_MelaWH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadWH_SIG_ghw1_1_JHUGen_JECNominal * (useQGTagging ? TMath::Power(jetPgOverPq[0]*jetPgOverPq[1],1/3.) : 1. ) ) : -2 ;
  float D_ZHh = (nCleanedJetsPt30>=2) ? 1/(1+ c_MelaZH*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_HadZH_SIG_ghz1_1_JHUGen_JECNominal * (useQGTagging ? TMath::Power(jetPgOverPq[0]*jetPgOverPq[1],1/3.) : 1. ) ) : -2 ;

  if( nExtraLep==0 && nCleanedJetsPt30==1 && D_VBF1j>WP_VBF1j ){

    return VBF1jTaggedIchep16;

  }else if( nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged==0)) && D_VBF2j>WP_VBF2j ){

    return VBF2jTaggedIchep16;

  }else if(   ( nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged==0)) && (D_WHh>WP_WHh||D_ZHh>WP_ZHh) )
	   || ( nExtraLep==0 && (nCleanedJetsPt30==2||nCleanedJetsPt30==3) && nCleanedJetsPt30BTagged>=2 ) ){

    return VHHadrTaggedIchep16;

  }else if(   ( nCleanedJetsPt30<=3 && nCleanedJetsPt30BTagged==0 && (nExtraLep==1||nExtraZ>=1) )
	   || ( nCleanedJetsPt30==0 && nExtraLep>=1 ) ){

    return VHLeptTaggedIchep16;

  }else if(   ( nCleanedJetsPt30>=4 && nCleanedJetsPt30BTagged>=1 )
	   || nExtraLep>=1 ){

    return ttHTaggedIchep16;

  }else{

    return UntaggedIchep16;

  }

}


extern "C" int categoryLegacy( int nCleanedJetsPt30 )
{
  if(nCleanedJetsPt30>=2)
    return Dijet;
  else
    return ZeroOneJet;
}
