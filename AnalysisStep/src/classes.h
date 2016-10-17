#include <DataFormats/PatCandidates/interface/PFParticle.h>
#include <DataFormats/PatCandidates/interface/UserData.h>
#include <ZZAnalysis/AnalysisStep/interface/ExtendedBranch.h>

edm::Ptr<pat::PFParticle> dummy1;
pat::UserHolder<std::vector<edm::Ptr<pat::PFParticle> > > dummy2;

BranchHelpers::ExtendedBranch<Bool_t> bdummy_b(0,"",(Bool_t)false);
BranchHelpers::ExtendedBranch<Char_t> bdummy_c(0, "", (Char_t)' ');
BranchHelpers::ExtendedBranch<Short_t> bdummy_s(0, "", (Short_t)0);
BranchHelpers::ExtendedBranch<Int_t> bdummy_i(0, "", (Int_t)0);
BranchHelpers::ExtendedBranch<Long64_t> bdummy_l(0, "", (Long64_t)0);
BranchHelpers::ExtendedBranch<Float_t> bdummy_f(0, "", (Float_t)0);
BranchHelpers::ExtendedBranch<Double_t> bdummy_d(0, "", (Double_t)0);
