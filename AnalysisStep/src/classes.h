#include <DataFormats/PatCandidates/interface/PFParticle.h>
#include <DataFormats/PatCandidates/interface/UserData.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonSFHelper.h>
#include <ZZAnalysis/AnalysisStep/interface/kFactors.h>
#include <ZZAnalysis/AnalysisStep/interface/WeightCalculatorFromHistogram.h>
#include <ZZAnalysis/AnalysisStep/interface/MuonScaRe.h>
#include <vector>


edm::Ptr<pat::PFParticle> dummy1;
pat::UserHolder<std::vector<edm::Ptr<pat::PFParticle> > > dummy2;

#include <JHUGenMELA/MELA/interface/Mela.h>
