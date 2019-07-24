#ifndef PhotonIDHelper_h
#define PhotonIDHelper_h

/** \class PhotonIDHelper
 *
 *  Helper for computing photon identification
 *
 *  \author T.Sculac
 */

#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidate.h>

namespace PhotonIDHelper {

   bool isCutBasedID_Loose(int year, const pat::Photon& photon);

}
#endif


