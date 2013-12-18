#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/Merger.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>

// A merger which produces a collection of pat::CompositeCandidate
// (instead of a collection of reco::CompositeCandidate)
typedef Merger<pat::CompositeCandidateCollection, pat::CompositeCandidateCollection> PATCompositeCandidateMerger;


DEFINE_FWK_MODULE( PATCompositeCandidateMerger );


#include <CommonTools/UtilAlgos/interface/StringCutObjectSelector.h>
#include <CommonTools/CandAlgos/interface/CandCombiner.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>


// A CandCombiner which produces a collection of pat::CompositeCandidate
// (instead of a collection of reco::CompositeCandidate)
namespace reco {
  namespace modules {
    typedef CandCombiner<
      StringCutObjectSelector<reco::Candidate, true>,
      AnyPairSelector,
      combiner::helpers::ShallowClone,
      pat::CompositeCandidateCollection
      > PATCandViewShallowCloneCombiner;

DEFINE_FWK_MODULE( PATCandViewShallowCloneCombiner );

  }
}

