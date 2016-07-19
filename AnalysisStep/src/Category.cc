#include <ZZAnalysis/AnalysisStep/interface/Category.h>

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
			     float pvbf_VAJHU_highestPTJets,
			     float phjj_VAJHU_highestPTJets
			     )
{
  float vbfMela = pvbf_VAJHU_highestPTJets / ( phjj_VAJHU_highestPTJets + pvbf_VAJHU_highestPTJets );

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
	     float phjj_VAJHU_highestPTJets,
	     float phj_VAJHU,
	     float pvbf_VAJHU_highestPTJets,
	     float pAux_vbf_VAJHU,
	     float pwh_hadronic_VAJHU,
	     float pzh_hadronic_VAJHU,
             float* jetPhi,
	     float ZZMass,
	     bool useQGTagging 
	     )
{

  float WP_VBF2j, WP_VBF1j, WP_WHh, WP_ZHh;
  if(useQGTagging){
    WP_VBF2j = 0.926;
    WP_VBF1j = 0.829;
    WP_WHh = 0.973;
    WP_ZHh = 0.996;
  }else{
    //WP_VBF2j = 0.938; // non mass-dependent cut
    WP_VBF2j = 0.536+665./(ZZMass+1530.); // mass-dependent cut
    WP_VBF1j = 0.815;
    WP_WHh = 0.959;
    WP_ZHh = 0.9946;
  }

  float c_Mela2j = 0.06;
  float c_Mela1j = 0.3;
  float c_MelaWH = 100000.;
  float c_MelaZH = 10000.;

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

  float D_VBF2j = (nCleanedJetsPt30>=2) ? 1/(1+ c_Mela2j*phjj_VAJHU_highestPTJets/pvbf_VAJHU_highestPTJets * ( useQGTagging ? TMath::Power(jetPgOverPq[0]*jetPgOverPq[1],1/3.) : 1. ) ) : -2 ;
  float D_VBF1j = (nCleanedJetsPt30>=1) ? 1/(1+ (c_Mela1j*phj_VAJHU)/(pvbf_VAJHU_highestPTJets*pAux_vbf_VAJHU) * (useQGTagging ? TMath::Power(jetPgOverPq[0],1/3.) : 1. ) ) : -2 ;
  float D_WHh = (nCleanedJetsPt30>=2) ? 1/(1+ c_MelaWH*phjj_VAJHU_highestPTJets/pwh_hadronic_VAJHU * (useQGTagging ? TMath::Power(jetPgOverPq[0]*jetPgOverPq[1],1/3.) : 1. ) ) : -2 ;
  float D_ZHh = (nCleanedJetsPt30>=2) ? 1/(1+ c_MelaZH*phjj_VAJHU_highestPTJets/pzh_hadronic_VAJHU * (useQGTagging ? TMath::Power(jetPgOverPq[0]*jetPgOverPq[1],1/3.) : 1. ) ) : -2 ;

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
