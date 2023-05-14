#ifndef EwkCorrections_h
#define EwkCorrections_h


#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

//Load here all the dataformat that we will need
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/View.h" 

//need for the good lumi filter
#include "FWCore/Utilities/interface/Algorithms.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

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
