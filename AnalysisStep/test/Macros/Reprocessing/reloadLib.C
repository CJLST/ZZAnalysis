{

  gSystem->AddIncludePath("-I$ROOFITSYS/include/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/");
  gSystem->Load("libZZMatrixElementMELA.so");
  gSystem->Load("$CMSSW_BASE/src/ZZMatrixElement/MELA/data/$SCRAM_ARCH/libmcfm_6p8.so");
  gROOT->LoadMacro("$CMSSW_BASE/src/ZZMatrixElement/MELA/interface/Mela.h++");
//  gROOT->LoadMacro("$CMSSW_BASE/src/ZZMatrixElement/MELA/interface/TVar.hh+");
//  gROOT->LoadMacro("/afs/cern.ch/user/u/usarica/work/snowmass2013/src/RooSpinZero_7DComplex_withAccep_withFepspr.cc+");
//  gROOT->LoadMacro("/afs/cern.ch/user/u/usarica/work/snowmass2013/src/ScalarPdfFactory_withFepspr.cc+");
/*  gROOT->LoadMacro("/afs/cern.ch/user/u/usarica/work/snowmass2013/src/RooSpinZero_KDInt_ZH_fepspr.cc++");
  gROOT->LoadMacro("/afs/cern.ch/user/u/usarica/work/snowmass2013/src/RooSpinZero_KD_ZH_fepspr.cc++");
  gROOT->LoadMacro("./Pdfs/RooSpinZero_KDInt_fepspr_withBkg.cc++");
  gROOT->LoadMacro("./Pdfs/RooSpinZero_KD_fepspr_withBkg.cc++");
  gROOT->LoadMacro("./Pdfs/RooSpinZero_KDInt_withBkg.cc++");
  gROOT->LoadMacro("./Pdfs/RooSpinZero_KD_withBkg.cc++");
*/

// loading tdr style for plots
//  gROOT->ProcessLine(".L  tdrstyle.C");
//  setTDRStyle();  

}
