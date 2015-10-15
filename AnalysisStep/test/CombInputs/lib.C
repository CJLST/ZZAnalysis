#include "TROOT.h"
#include "TSystem.h"

#include <iostream>


void lib(){
gSystem->Load("libRooFit.so");
gSystem->Load("libHiggsAnalysisCombinedLimit.so");
gSystem->Load("libZZAnalysisAnalysisStep.so");
gSystem->AddIncludePath("-I$ROOFITSYS");
gSystem->AddIncludePath("-I$CMSSW_BASE");
}
void lib();

