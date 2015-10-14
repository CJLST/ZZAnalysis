#include <ZZAnalysis/AnalysisStep/interface/Category.h>

#include <cmath>

#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >  LV;


int flagDijetVH(
		int nJets, 
		float* jetpt,
		float* jeteta,
		float* jetphi,
		float* jetmass
		)
{

  bool found = false;

  if(nJets>=2){

    for(int j1=0; j1<nJets; j1++){
      if( std::abs(jeteta[j1])<2.4 && jetpt[j1]>40. ){

	for(int j2=j1+1; j2<nJets; j2++){
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
	     int nExtraLeptons,
	     float ZZPt,
	     float ZZMass,
	     int nJets, 
	     int nBTaggedJets,
	     float* jetpt,
	     float* jeteta,
	     float* jetphi,
	     float* jetmass,
	     float Fisher
	     )
{

  int category = -1;

  if( nExtraLeptons==0 && nJets>=2 && nBTaggedJets<=1 && Fisher>0.5 ){

    category = VBFTagged; 

  }else if( ( nExtraLeptons==0 && nJets>=2 && ZZPt>ZZMass && flagDijetVH(nJets,jetpt,jeteta,jetphi,jetmass) )
            || ( nExtraLeptons==0 && nJets==2 && nBTaggedJets==2 ) ){

    category = VHHadrTagged;

  }else if( nExtraLeptons>=1 && nJets<=2 && nBTaggedJets==0 ){

    category = VHLeptTagged;

  }else if( nExtraLeptons>=1 || (nJets>=3 && nBTaggedJets>=1) ){

    category = ttHTagged;

  }else if(nJets>=1){

    category = OneJetTagged;

  }else{

    category = Untagged;

  }

  return category;

}
