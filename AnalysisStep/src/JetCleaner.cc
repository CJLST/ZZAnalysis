#include <ZZAnalysis/AnalysisStep/interface/JetCleaner.h>
#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>

using namespace std;

bool jetCleaner::isGood(const pat::CompositeCandidate& cand, const pat::Jet& jet) {
  for(unsigned i=0; i<cand.numberOfDaughters(); ++i){
    if (cand.daughter(i)->numberOfDaughters()==0) { //FSR or single leptons. //FIXME: for Z+L, FSR is not currently attached as daughter!
      if (reco::deltaR(*(cand.daughter(i)), jet) < 0.4) {
	return false;
      }
    } else {
      for (unsigned j=0; j<cand.daughter(i)->numberOfDaughters(); ++j){
	if (reco::deltaR(*(cand.daughter(i)->daughter(j)), jet) < 0.4) {
	  return false;
	}
      }
    }
  }
  return true;
}  

