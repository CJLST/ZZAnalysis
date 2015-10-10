#ifndef LeptonIsoHelper_h
#define LeptonIsoHelper_h

/** \class LeptonIsoHelper
 *
 *  Helper for computing lepton isolation
 *
 *  \author N. Amapane
 */

#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidate.h>

namespace LeptonIsoHelper {

  extern int defaultCorrTypeMu;
  extern int defaultCorrTypeEle;

  /// Set the rho tag and the EA targed based on setup, for muons and electrons
  edm::InputTag getMuRhoTag(int sampleType, int setup);

  edm::InputTag getEleRhoTag(int sampleType, int setup);
  
  /// Compute combRelIsoPF for a mu
  float combRelIsoPF(int sampleType, int setup, double rho, const pat::Muon& mu, float fsr=0, int correctionType=defaultCorrTypeMu);
  
  /// Compute combRelIsoPF for an ele
  float combRelIsoPF(int sampleType, int setup, double rho, const pat::Electron& ele, float fsr=0, int correctionType=defaultCorrTypeEle);
  
  /// Generic version, assuming that Lep is a PATObject; calls one of the above
  float combRelIsoPF(int sampleType, int setup, double rho, const reco::Candidate* lep, float fsr=0, int correctionType=-1);

  // Compute FSR isolation
  void fsrIso(const reco::PFCandidate* photon, edm::Handle<edm::View<pat::PackedCandidate> > pfcands, double& ptSumNe, double& ptSumCh, double& ptSumChByWorstPV);

  // Value of the iso cut 
  float isoCut(const reco::Candidate* d);

}
#endif


