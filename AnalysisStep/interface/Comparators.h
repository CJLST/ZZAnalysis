#ifndef COMPARATORS_H
#define COMPARATORS_H

/** \class Comparators
 *
 *  No description available.
 *
 *  $Date: 2012/06/06 00:27:43 $
 *  $Revision: 1.3 $
 *  \author S. Casasso - TORINO
 *  \author N. Amapane - CERN

 */

#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>

#include <algorithm>

namespace Comparators {

  const float ZmassValue = 91.1876;

  // Needed to sort best candidates after pre-selection
  struct candWithIdx{
    int idxCand;
    double mZ1;
    double ptSumZ2;

    // candWithIdx(int i, pat::CompositeCandidate* theCand):
    candWithIdx(int i, double mZ1, double ptSumZ2):
      idxCand(i),
      mZ1(mZ1),
      ptSumZ2(ptSumZ2)
    {}    

  };

  struct bestZ1bestZ2 {
    bool operator() (Comparators::candWithIdx* candWithIdx1, Comparators::candWithIdx* candWithIdx2){
      
      double mZ1_1 = candWithIdx1->mZ1;
      double mZ1_2 = candWithIdx2->mZ1;
      
      double ptSumZ2_1 = candWithIdx1->ptSumZ2;
      double ptSumZ2_2 = candWithIdx2->ptSumZ2;
      
      if ( fabs(mZ1_1-mZ1_2) < 1e-04 ){
	if (ptSumZ2_1 > ptSumZ2_2) {
	  return true;
	}
	else {
	  return false;
	}
      }
      else {
	if ( mZ1_1-mZ1_2 > 0. ) {
	  return true;
	}
	else {
	  return false;
	}
      }      

   };
  }; //bestZ1bestZ1

  
}// end namespace Comparators
  

#endif
