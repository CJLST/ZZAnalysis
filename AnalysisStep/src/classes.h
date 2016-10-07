#include <DataFormats/PatCandidates/interface/PFParticle.h>
#include <DataFormats/PatCandidates/interface/UserData.h>
#include <ZZAnalysis/AnalysisStep/interface/ExtendedBranch.h>

edm::Ptr<pat::PFParticle> dummy1;
pat::UserHolder<std::vector<edm::Ptr<pat::PFParticle> > > dummy2;

BranchHelpers::ExtendedBranch<Bool_t> bdummy_b;
BranchHelpers::ExtendedBranch<Char_t> bdummy_c;
BranchHelpers::ExtendedBranch<Short_t> bdummy_s;
BranchHelpers::ExtendedBranch<Int_t> bdummy_i;
BranchHelpers::ExtendedBranch<Long64_t> bdummy_l;
BranchHelpers::ExtendedBranch<Float_t> bdummy_f;
BranchHelpers::ExtendedBranch<Double_t> bdummy_d;
BranchHelpers::ExtendedBranch<vectorBool_t> bdummy_B;
BranchHelpers::ExtendedBranch<vectorChar_t> bdummy_C;
BranchHelpers::ExtendedBranch<vectorShort_t> bdummy_S;
BranchHelpers::ExtendedBranch<vectorInt_t> bdummy_I;
BranchHelpers::ExtendedBranch<vectorFloat_t> bdummy_F;
BranchHelpers::ExtendedBranch<vectorDouble_t> bdummy_D;
