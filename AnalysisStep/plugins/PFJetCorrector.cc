#if CMSSW_VERSION>500

#include "AnalysisDataFormats/CMGTools/interface/PFJet.h"
#include "CMGTools/Common/plugins/JetEnergyCorrector.h"

typedef cmg::JetEnergyCorrector<cmg::PFJet> PFJetCorrector;

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(PFJetCorrector);

#endif
