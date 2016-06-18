#ifndef DaughterDataEmbedder_h
#define DaughterDataEmbedder_h

/** \class DaughterDataEmbedder
 *
 *  No description available.
 *
 *  $Date: 2012/06/06 00:27:43 $
 *  $Revision: 1.3 $
 *  \author N. Amapane - CERN
 */


#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <ZZAnalysis/AnalysisStep/interface/PhotonFwd.h>

#include <string>
#include <vector>
#include <sstream>


namespace userdatahelpers {

  /// Retrieve the userFloat "name" from a reco::Candidate c
  float getUserFloat(const reco::Candidate* c, const char* name);

  /// Test if the userFloat "name" from a reco::Candidate c exists
  int hasUserFloat(const reco::Candidate* c, const char* name);


  /// Retrieve matched photons stored as userData
  const PhotonPtrVector* getUserPhotons(const reco::Candidate* c);

  /// Add the userFloats of daughters into cand, with proper prefix 
  void embedDaughterData(pat::CompositeCandidate& cand);

  /// Add the userFloats of daughter d into cand, with proper prefix
  template <typename T>
  void embedDaughterData(pat::CompositeCandidate& cand, unsigned i, const T* d) {
    using namespace std;

    string base;
    stringstream str;
    str << "d" << i << ".";
    str >> base;

    /*
    if(d->isElectron()) edm::LogPrint("") << "Filling Electron daughter";
    if(d->isPhoton()) edm::LogPrint("") << "Filling Photon daughter";
    if(d->isElectron() || d->isPhoton()) { 
        edm::LogPrint("") << "eta " << d->eta() << " pt " << d->pt(); 
        edm::LogPrint("") << "isBDT " << d->userFloat("isBDT");
    }*/
    const vector<string> & userLabels = d->userFloatNames();
    for (vector<string>::const_iterator name = userLabels.begin(); name!= userLabels.end(); ++name){      
      string newname = base + *name;
      //if(d->isElectron() || d->isPhoton())
      //  edm::LogPrint("") << "In DaughterHelper adding "<< newname << " value " << d->userFloat(*name);
      cand.addUserFloat(newname, d->userFloat(*name));
    }
    // edm::LogPrint("") << "--------------------------------------------";
  
  }


  /// Get the daughters daughters sorted by Z1/Z2 and charge so that the order is Z1Lp, Z1Ln, Z2Lp, Z2Ln
  /// Keeps the original sorting for the same-sign collections used for CRs
  /// Use with is4l=false for a consistent behaviour with the Z+l collection
  /// Results: 
  /// labels are the prefixes for cand.userFloat() string of each lepton.
  /// fsrIndex are indices of the leptons corresponding to each photon after sorting
  void getSortedLeptons(const pat::CompositeCandidate& cand, 
			std::vector<const reco::Candidate*>& leptons, 
			std::vector<std::string>& labels, 
			std::vector<const reco::Candidate*>& fsr,
			std::vector<short>& fsrIndex,
			bool is4l=true);

  void getSortedZLeptons(const pat::CompositeCandidate& cand, 
			std::vector<const reco::Candidate*>& leptons, 
			std::vector<std::string>& labels, 
			std::vector<const reco::Candidate*>& fsr,
			std::vector<short>& fsrIndex);


}
#endif
