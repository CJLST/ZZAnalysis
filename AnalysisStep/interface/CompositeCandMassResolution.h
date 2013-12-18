#ifndef CompositeCandMassResolution_h
#define CompositeCandMassResolution_h

namespace reco { class Candidate; class Muon; class GsfElectron; class Track; class PFCandidate; class LeafCandidate; }
namespace edm { class EventSetup; }

#include <vector>
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include <TMatrixDSym.h>

/* Mutuated from 
   /UserCode/Mangano/WWAnalysis/AnalysisStep/src/CompositeCandMassResolution.cc

 */

class CompositeCandMassResolution  {
    public:
        CompositeCandMassResolution() {}
        ~CompositeCandMassResolution() {}
        void init(const edm::EventSetup &iSetup);
        double getMassResolution(const reco::Candidate &c) const ;
        double getMassResolutionWithComponents(const reco::Candidate &c, std::vector<double> &errI) const;
    private:
        void   fillP3Covariance(const reco::Candidate &c, TMatrixDSym &bigCov, int offset) const ;
        void   fillP3Covariance(const reco::GsfElectron &c, TMatrixDSym &bigCov, int offset) const ;
        void   fillP3Covariance(const reco::Muon &c, TMatrixDSym &bigCov, int offset) const ;
        void   fillP3Covariance(const reco::PFCandidate &c, TMatrixDSym &bigCov, int offset) const ;
        void   fillP3Covariance(const reco::Candidate &c, const reco::Track &t, TMatrixDSym &bigCov, int offset) const ;
	void   fillP3Covariance(const reco::LeafCandidate &c, TMatrixDSym &bigCov, int offset) const ;

        edm::ESHandle<MagneticField> magfield_;

        // 1 if this is a lead, recursive number of leafs if composite
        void getLeaves(const reco::Candidate &c, std::vector<const reco::Candidate *> &out) const ;

        double getMassResolution_(const reco::Candidate &c, std::vector<double> &errI, bool doComponents) const;
};

#endif

