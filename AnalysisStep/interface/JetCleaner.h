#ifndef JetCleaner_h
#define JetCleaner_h

/**
 *
 *  Clean jets from leptons of a CompositeCandidate, and their FSR (assuming the CC has either Z or leptons abd FSR 
 *  as daughters, as in our ZZ and Z+l events)
 *
 */

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"


namespace jetCleaner {
  bool isGood(const pat::CompositeCandidate& cand, const pat::Jet& jet);
};
#endif

