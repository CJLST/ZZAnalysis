#ifndef VBFCandidateJetSelector_h
#define VBFCandidateJetSelector_h

/** \class VBFCandidateJetSelector
 *
 *  Sete
 *
 *  $Date: 2013/05/16 14:29:13 $
 *  $Revision: 1.3 $
 *
 */

#include <vector>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include "DataFormats/Math/interface/deltaR.h"

class VBFCandidateJetSelector {
 public:
	
  
  /// Constructor
	VBFCandidateJetSelector() {};

  /// Destructor
	virtual ~VBFCandidateJetSelector() {};

	typedef std::vector<const pat::Jet*> container ;
	
	std::vector<const pat::Jet*> cleanJets(const pat::CompositeCandidate& cand,
						 edm::Handle<edm::View<pat::Jet> > jets, int year);
	
private:       
	container selected_ ;

};
#endif




