#include "../interface/PileUpWeight.h"

float PileUpWeight::weight(float input, PileUpWeight::PUvar var) {

  if (h_nominal != nullptr && var == PUvar::NOMINAL) {
    return h_nominal->GetBinContent(h_nominal->FindBin(input));
  } else if (h_down != nullptr && var == PUvar::VARDOWN) {
    return h_down->GetBinContent(h_down->FindBin(input));
  } else if (h_up != nullptr && var == PUvar::VARUP) {
    return h_up->GetBinContent(h_up->FindBin(input));
  } else {
    return -1.;
  }
}


PileUpWeight::PileUpWeight(int MC, int target) { 

 if (MC==2016&&target==2016) {
    edm::FileInPath fip("ZZAnalysis/AnalysisStep/data/PileUpWeights/puWeightsMoriond17_v2.root");

    TFile *fPUWeight = TFile::Open(fip.fullPath().data(),"READ");

    h_nominal.reset((TH1*)fPUWeight->Get("weights")->Clone());
    h_up.reset((TH1*)fPUWeight->Get("weights_varUp")->Clone());
    h_down.reset((TH1*)fPUWeight->Get("weights_varDn")->Clone());

    fPUWeight->Close();
  }
  if(h_nominal == nullptr) {
     edm::LogError("PU reweight") << "Did not find reweighting histogram to reweight MC=" << MC << " to data=" << target;
     edm::LogError("PU reweight") << "Legacy configurations are in PUReweight.cc";
  }
}
