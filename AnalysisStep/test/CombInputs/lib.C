#include "TROOT.h"
#include "TSystem.h"

#include <iostream>


void lib(){
gSystem->Load("libRooFit.so");
gSystem->Load("libHiggsAnalysisCombinedLimit.so");
gSystem->Load("libZZAnalysisAnalysisStep.so");
gSystem->AddIncludePath("-I/cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/roofit/5.34.18-cms12/include/");
gSystem->AddIncludePath("-I/afs/cern.ch/user/m/mkiani/work/Analysis/CMSSW_7_2_4/src");
}
void lib();

