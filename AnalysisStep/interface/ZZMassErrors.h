#ifndef ZZMassErrors_h
#define ZZMassErrors_h

/** \class ZZMassErrors
 *
 *  No description available.
 *
 *  $Date: 2012/05/05 09:08:33 $
 *  $Revision: 1.1 $
 *  \author N. Amapane - CERN
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/PatCandidates/interface/Electron.h"


class ZZMassErrors {
public:
  /// Constructor
  ZZMassErrors(const edm::EventSetup& eventSetup);

  /// Destructor
  virtual ~ZZMassErrors(){};
  
  // Operations
  float massError(const reco::Candidate* Z1Lp, 
		  const reco::Candidate* Z1Lm, 
		  const reco::Candidate* Z2Lp, 
		  const reco::Candidate* Z2Lm);

  double pError(const pat::Electron* ele);

 private:
  edm::ESHandle<MagneticField> magfield;
  
};
#endif

