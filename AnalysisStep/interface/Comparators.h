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

  struct bestZ1bestZ2 {
    
    pat::CompositeCandidateCollection candCollection;
    bestZ1bestZ2(pat::CompositeCandidateCollection candCollection):
      candCollection(candCollection)
    {}

    bool operator() (int i, int j){
      
      pat::CompositeCandidate& cand_i = candCollection[i];
      pat::CompositeCandidate& cand_j = candCollection[j];

      int idxZ1_i = 0; 
      int idxZ2_i = 1;
      if(fabs(cand_i.daughter(0)->masterClone().get()->mass()-Comparators::ZmassValue)>=fabs(cand_i.daughter(1)->masterClone().get()->mass()-Comparators::ZmassValue)) swap(idxZ1_i,idxZ2_i);
      const reco::Candidate* Z1_i = cand_i.daughter(idxZ1_i)->masterClone().get();
      
      int idxZ1_j = 0; 
      int idxZ2_j = 1;
      if(fabs(cand_j.daughter(0)->masterClone().get()->mass()-Comparators::ZmassValue)>=fabs(cand_j.daughter(1)->masterClone().get()->mass()-Comparators::ZmassValue)) swap(idxZ1_j,idxZ2_j);
      const reco::Candidate* Z1_j = cand_j.daughter(idxZ1_j)->masterClone().get();
      
      if ( fabs(Z1_i->mass()-Z1_j->mass()) < 1e-04 ){
	const reco::Candidate* Z2_i = cand_i.daughter(idxZ2_i)->masterClone().get();
	const reco::Candidate* Z2_j = cand_j.daughter(idxZ2_j)->masterClone().get();
	double ptSumZ2_i = Z2_i->daughter(0)->masterClone().get()->pt() + Z2_i->daughter(1)->masterClone().get()->pt();
	double ptSumZ2_j = Z2_j->daughter(0)->masterClone().get()->pt() + Z2_j->daughter(1)->masterClone().get()->pt();
	if (ptSumZ2_i > ptSumZ2_j) {
	  return true;
	}
	else {
	  return false;
	}
      }
      else {
	if ( Z1_i->mass()-Z1_j->mass() > 0. ) {
	  return true;
	}
	else {
	  return false;
	}
      }      

    }; //operator()    
    

  }; //bestZ1bestZ1

  
}// end namespace Comparators
  

#endif
