#ifndef EwkCorrections_h
#define EwkCorrections_h


#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

//Load here all the dataformat that we will need
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/HepMCCandidate/interface/PdfInfo.h"

//need for the good lumi filter
#include "DataFormats/Provenance/interface/LuminosityBlockID.h"
#include "DataFormats/Provenance/interface/LuminosityBlockRange.h"
#include "FWCore/Utilities/interface/Algorithms.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

//#include "UserCode/llvv_fwk/interface/MacroUtils.h"
//#include "UserCode/llvv_fwk/interface/LumiUtils.h"
//#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"

// Electron ID
#include "RecoEgamma/ElectronIdentification/interface/VersionedPatElectronSelector.h"

#include <vector>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TGraph.h"
#include <Math/VectorUtil.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

inline bool sort_CandidatesByPt(const reco::Candidate *a, const reco::Candidate *b) { return a->pt()>b->pt(); }

namespace EwkCorrections
{
  std::vector<std::vector<float>> readFile_and_loadEwkTable(TString dtag);
  std::vector<float> findCorrection(const std::vector<std::vector<float>> & Table_EWK, float sqrt_s_hat, float t_hat);
  double getEwkCorrections(const edm::Handle<edm::View<reco::Candidate> > & particles, 
                           const std::vector<std::vector<float>> & Table, 
                           const GenEventInfoProduct & eventInfo,
                           TLorentzVector Z1, TLorentzVector Z2);
}

#endif
