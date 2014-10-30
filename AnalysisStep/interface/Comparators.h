#ifndef COMPARATORS_H
#define COMPARATORS_H

/** \class Comparators
 *
 *  No description available.
 *
 *  $Date: 2012/06/06 00:27:43 $
 *  $Revision: 1.3 $
 *  \author S. Casasso - TORINO
 *  \author N. Amapane - TORINO

 */

#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>

#include <algorithm>
#include <limits>

namespace Comparators {

  const float ZmassValue = 91.1876;

  enum ComparatorTypes {byBestZ1bestZ2 = 0, byBestKD=1, byBestDbkgVH};

  // Abstract base class
  struct BestCandComparator {
    const pat::CompositeCandidateCollection& candCollection;
    
    BestCandComparator(const pat::CompositeCandidateCollection& candCollection, ComparatorTypes type):
      candCollection(candCollection),
      theType(type)
    {}

    bool bestZ1bestZ2(const pat::CompositeCandidate& cand_i, const pat::CompositeCandidate& cand_j){

      float mZ1_i = cand_i.daughter("Z1")->mass();
      float mZ1_j = cand_j.daughter("Z1")->mass();
      if ( fabs(mZ1_i-mZ1_j) < 1e-04 ){ // same Z1: choose the candidate with highest-pT Z2 leptons
	const reco::Candidate* Z2_i = cand_i.daughter("Z2");
	const reco::Candidate* Z2_j = cand_j.daughter("Z2");
	double ptSumZ2_i = Z2_i->daughter(0)->pt() + Z2_i->daughter(1)->pt();
	double ptSumZ2_j = Z2_j->daughter(0)->pt() + Z2_j->daughter(1)->pt();

	// cout <<  "Comparator: compare Z2 pTs: " << ptSumZ2_i << " " << ptSumZ2_j << endl;
	return ( ptSumZ2_i > ptSumZ2_j );
      }
      else { // choose the candidate with Z1 closest to nominal mass 
	// cout << "Comparator: compare Z1 masses: " << mZ1_i << " " << mZ1_j << endl;
	return ( fabs(mZ1_i-ZmassValue)<fabs(mZ1_j-ZmassValue) );
      }      
      
    }


    virtual bool operator() (int i, int j){

      const pat::CompositeCandidate& cand_i = candCollection[i];
      const pat::CompositeCandidate& cand_j = candCollection[j];

      if (theType==byBestZ1bestZ2) { // "Legacy" best cand logic	
	return bestZ1bestZ2(cand_i,cand_j);
      } //end of byBestZ1bestZ2
  
      else if (theType==byBestKD) {	
	double KD_i = cand_i.userFloat("p0plus_VAJHU")/( cand_i.userFloat("p0plus_VAJHU") + cand_i.userFloat("bkg_VAMCFM") );
	double KD_j = cand_j.userFloat("p0plus_VAJHU")/( cand_j.userFloat("p0plus_VAJHU") + cand_j.userFloat("bkg_VAMCFM") );

	bool kdEqual = KD_i == KD_j;

	if (kdEqual) return bestZ1bestZ2(cand_i,cand_j); //same 4 leptons, different pairing
	else return (KD_i>KD_j);

      } // end of byBestKD

      else if (theType==byBestDbkgVH) {
	double ps_i = cand_i.userFloat("pzh_VAJHU");
	if (ps_i<0) ps_i = cand_i.userFloat("p0plus_VAJHU");
	double ps_j = cand_j.userFloat("pzh_VAJHU");
	if (ps_j<0) ps_j = cand_j.userFloat("p0plus_VAJHU");
	
	double KD_i = ps_i/( ps_i + cand_i.userFloat("bkg_VAMCFM") );
	double KD_j = ps_j/( ps_j + cand_j.userFloat("bkg_VAMCFM") );
      
	return (KD_i>KD_j);
      } // end of byBestKD
      
      else {
	abort();
      }

    } // end of operator()    
    
    ComparatorTypes theType;
    
  }; // BestCandComparator


  
}// end namespace Comparators
  

#endif
