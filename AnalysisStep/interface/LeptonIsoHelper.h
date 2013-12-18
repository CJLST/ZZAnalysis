#ifndef LeptonIsoHelper_h
#define LeptonIsoHelper_h

/** \class LeptonIsoHelper
 *
 *  Helper for computing lepon isolation
 *
 *  $Date: 2012/06/10 17:40:50 $
 *  $Revision: 1.3 $
 *  \author N. Amapane
 */

#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>

namespace LeptonIsoHelper {

  /// Set the rho tag and the EA targed based on setup, for muons and electrons
  edm::InputTag getMuRhoTag(int sampleType, int setup);

  edm::InputTag getEleRhoTag(int sampleType, int setup);
  
  /// Compute combRelIsoPF for a mu
  float combRelIsoPF(int sampleType, int setup, double rho, const pat::Muon& mu, float fsr=0);
  
  /// Compute combRelIsoPF for an ele
  float combRelIsoPF(int sampleType, int setup, double rho, const pat::Electron& ele, float fsr=0);
  
  /// Generic version, assuming that Lep is a PATObject; calls one of the above
  float combRelIsoPF(int sampleType, int setup, double rho, const reco::Candidate* lep, float fsr=0);
}
#endif


