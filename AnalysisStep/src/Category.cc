#include <ZZAnalysis/AnalysisStep/interface/Category.h>
#include <cmath>

int category(
	     int nExtraLeptons,
	     float ZZPt,
	     float ZZMass,
	     int nJets, 
	     int nBTaggedJets,
	     float ptj1,
	     float ptj2,
	     float etaj1,
	     float etaj2,
	     float mjj,
	     float Fisher
	     )
{

  int category = -1;
  // 0 = Untagged  
  // 1 = 1-jet tagged  
  // 2 = VBF tagged  
  // 3 = VH-leptonic tagged  
  // 4 = VH-hadronic tagged  
  // 5 = ttH tagged  

  if( nExtraLeptons==0 && nJets>=2 && nBTaggedJets<=1 && Fisher>0.5 ){

    category = 2; // VBF tagged

  }else if( ( nExtraLeptons==0 && nJets>=2 && 
	         ptj1>40 && ptj2>40 && fabs(etaj1)<2.4 && fabs(etaj2)<2.4 && 
              60<mjj && mjj<120 && ZZPt>ZZMass )
            || ( nExtraLeptons==0 && nJets==2 && nBTaggedJets==2 ) ){

    category = 4; // VH-hadronic tagged

  }else if( nExtraLeptons>=1 && nJets<=2 && nBTaggedJets==0 ){

    category = 3; // VH-leptonic tagged

  }else if( nExtraLeptons>=1 || (nJets>=3 && nBTaggedJets>=1) ){

    category = 5; // ttH tagged

  }else if(nJets>=1){

    category = 1; // 1-jet tagged

  }else{

    category = 0; // Untagged

  }

  return category;

}
