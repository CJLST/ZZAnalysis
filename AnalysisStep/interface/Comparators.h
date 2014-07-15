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
    
    const pat::CompositeCandidateCollection& candCollection;
    bestZ1bestZ2(const pat::CompositeCandidateCollection& candCollection):
      candCollection(candCollection)
    {}

    bool operator() (int i, int j){
      
      const pat::CompositeCandidate& cand_i = candCollection[i];
      const pat::CompositeCandidate& cand_j = candCollection[j];

      const reco::Candidate* Z1_i = cand_i.daughter("Z1")->masterClone().get();      
      const reco::Candidate* Z1_j = cand_j.daughter("Z1")->masterClone().get();
      
      if ( fabs(Z1_i->mass()-Z1_j->mass()) < 1e-04 ){
	const reco::Candidate* Z2_i = cand_i.daughter("Z2")->masterClone().get();
	const reco::Candidate* Z2_j = cand_j.daughter("Z2")->masterClone().get();
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



  struct bestKD {

    const pat::CompositeCandidateCollection& candCollection;
    bestKD(const pat::CompositeCandidateCollection& candCollection):
      candCollection(candCollection)
    {}

    bool operator() (int i, int j){
      
      const pat::CompositeCandidate& cand_i = candCollection[i];
      const pat::CompositeCandidate& cand_j = candCollection[j];

      double KD_i = cand_i.userFloat("p0plus_VAJHU")/( cand_i.userFloat("p0plus_VAJHU") + cand_i.userFloat("bkg_VAMCFM") );
      double KD_j = cand_j.userFloat("p0plus_VAJHU")/( cand_j.userFloat("p0plus_VAJHU") + cand_j.userFloat("bkg_VAMCFM") );
      
      return (KD_i>KD_j);
      
    }// bestKD

    
  };

  
}// end namespace Comparators
  

#endif
