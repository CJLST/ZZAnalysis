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

  enum ComparatorTypes {byBestZ1bestZ2 = 0, byBestKD=1, byBestKD_VH=2, byBestPsig=3, byMHWindow=4};

  struct BestCandComparator {
    const pat::CompositeCandidateCollection& candCollection;
    
    BestCandComparator(const pat::CompositeCandidateCollection& candCollection, ComparatorTypes type):
      candCollection(candCollection),
      theType(type)
    {}

    virtual ~BestCandComparator() = default;

    // Helper function, implements the "legacy" Z1/Z2 selection logic
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

    // Sorting based on discriminants cannot distinguish between SF candidates with the same leptons and FSR,
    // but different pairing. We need to check for these cases.
    // Note that for SS CRs, this function will always return false.
    bool isEquivalent(const pat::CompositeCandidate& cand_i, const pat::CompositeCandidate& cand_j){
      if (std::abs(cand_i.mass()-cand_j.mass())<1e-4) { // tolerance: 100 keV
	// Check same FS, and that candidate is SF, OS. 
	int c_i = cand_i.userFloat("candChannel");
	int c_j = cand_i.userFloat("candChannel");
	if (c_i==c_j && (c_i==28561||c_j==14641)) return true;
      }
      return false;      
    }
    

    bool inMHWindow(const pat::CompositeCandidate& cand){
      const double m_min=120.;
      const double m_max=130.;
      double m =cand.mass();
      return (m>m_min&&m<m_max);
    }

    virtual bool operator() (int i, int j){

      const pat::CompositeCandidate& cand_i = candCollection[i];
      const pat::CompositeCandidate& cand_j = candCollection[j];


      // "Legacy" best cand logic	
      if (theType==byBestZ1bestZ2) { 
	return bestZ1bestZ2(cand_i,cand_j);
      } //end of byBestZ1bestZ2
  

      // Choose by best KD; for equivalent candidates (same leptons and FSR) that differ only on pairing,
      // choose based on legacy logic
      else if (theType==byBestKD) {	
	if (isEquivalent(cand_i,cand_j)) return bestZ1bestZ2(cand_i,cand_j); //same 4 leptons, different pairing
	double KD_i = cand_i.userFloat("p_GG_SIG_ghg2_1_ghz1_1_JHUGen")/( cand_i.userFloat("p_GG_SIG_ghg2_1_ghz1_1_JHUGen") + cand_i.userFloat("p_QQB_BKG_MCFM") );
	double KD_j = cand_j.userFloat("p_GG_SIG_ghg2_1_ghz1_1_JHUGen")/( cand_j.userFloat("p_GG_SIG_ghg2_1_ghz1_1_JHUGen") + cand_j.userFloat("p_QQB_BKG_MCFM") );
	return (KD_i>KD_j);
      } // end of byBestKD


      else if (theType==byBestKD_VH) {
	// FIXME: This is just a temporary implementation for tests!!
	// FIXME: Note that pzh_VAJHU can in some case be defined only for one of the two candidates (because of FSR/isolation). 
	// The comparison of KD does not make much sense in that case...
	double ps_i = cand_i.userFloat("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
	double pszh_i = 1.;
	if (pszh_i>0) ps_i *= pszh_i;

	double ps_j = cand_j.userFloat("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
	double pszh_j = 1.;
	if (pszh_j>0) ps_j *= pszh_j;
	
	double KD_i = ps_i/( ps_i + cand_i.userFloat("p_QQB_BKG_MCFM") );
	double KD_j = ps_j/( ps_j + cand_j.userFloat("p_QQB_BKG_MCFM") );

	if (isEquivalent(cand_i,cand_j)) return bestZ1bestZ2(cand_i,cand_j);
	else return (KD_i>KD_j);
      } // end of byBestKD_VH
      

      else if (theType==byBestPsig) {
	if (isEquivalent(cand_i,cand_j)) return bestZ1bestZ2(cand_i,cand_j); //same 4 leptons, different pairing
	return (cand_i.userFloat("p_GG_SIG_ghg2_1_ghz1_1_JHUGen")>cand_j.userFloat("p_GG_SIG_ghg2_1_ghz1_1_JHUGen"));
      }

      else if (theType==byMHWindow) {
	// Precedence for a candidate in the H mass window
	bool i_inMHWin=inMHWindow(cand_i);
	bool j_inMHWin=inMHWindow(cand_j);
	if (i_inMHWin&&!j_inMHWin) return true;
	if (j_inMHWin&&!i_inMHWin) return false;
	// If neither, or both candidate are in the mass window, use the old logic. Could choose by best psig as well.
	return bestZ1bestZ2(cand_i,cand_j); //same 4 leptons, different pairing
      }

      else {
	abort();
      }

    } // end of operator()    
    
    ComparatorTypes theType;
    
  }; // BestCandComparator


  
}// end namespace Comparators
  

#endif
