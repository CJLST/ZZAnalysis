void analyzer() {

  if (! TString(gSystem->GetLibraries()).Contains("ZZTreeReader")) {
    gROOT->LoadMacro("$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/test/Plotter/root_lib/TFileServiceWrapper.cc+");
    gROOT->LoadMacro("$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/test/Macros/HZZ4lBase.C+");
    gROOT->LoadMacro("$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/test/Plotter/ZZTreeReader.cc+");
  }


  //HCP TREES
 
  TString inputDir="/data_CMS/cms/salerno/CJLSTReduceTree/130702/"; //FIXME folder on polui
  TString outDir = "HZZ_plots/130702/";

  vector<TString> files;
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_DoubleEle_1963.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_DoubleMu_1963.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_DoubleOr_1963.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_H120.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_H125.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_H126.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_VBFH126.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_WH126.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_ZH126.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_ttH126.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_H350.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_ZZTo2e2mu.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_ZZTo2e2tau.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_ZZTo2mu2tau.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_ZZTo4e.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_ZZTo4mu.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_ZZTo4tau.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_ggZZ2l2l.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_ggZZ4l.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_DYJetsToLLTuneZ2M10B.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_DYJetsToLLTuneZ2M10NoB.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_DYJetsToLLTuneZ2M50B.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_DYJetsToLLTuneZ2M50NoB.root");
  files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_TTTo2L2Nu2B.root");
  //files.push_back("PRODFSR_8TeV/ZZ4lAnalysis_WZ.root");

  files.push_back("PRODFSR/ZZ4lAnalysis_DoubleEle.root");
  files.push_back("PRODFSR/ZZ4lAnalysis_DoubleMu.root");
  files.push_back("PRODFSR/ZZ4lAnalysis_DoubleOr.root");
  files.push_back("PRODFSR/ZZ4lAnalysis_H120.root");
  files.push_back("PRODFSR/ZZ4lAnalysis_H125.root");
  files.push_back("PRODFSR/ZZ4lAnalysis_H126.root");
  //files.push_back("PRODFSR/ZZ4lAnalysis_VBFH126.root");  //FIXME the 7 TeV does not use in the plot.r use 8TeV samples
  files.push_back("PRODFSR/ZZ4lAnalysis_WH126.root");
  files.push_back("PRODFSR/ZZ4lAnalysis_ZH126.root");
  files.push_back("PRODFSR/ZZ4lAnalysis_ttH126.root");
  files.push_back("PRODFSR/ZZ4lAnalysis_H350.root");
  files.push_back("PRODFSR/ZZ4lAnalysis_ZZTo2e2mu.root");
  files.push_back("PRODFSR/ZZ4lAnalysis_ZZTo2e2tau.root");
  files.push_back("PRODFSR/ZZ4lAnalysis_ZZTo2mu2tau.root");
  files.push_back("PRODFSR/ZZ4lAnalysis_ZZTo4e.root");
  files.push_back("PRODFSR/ZZ4lAnalysis_ZZTo4mu.root");
  files.push_back("PRODFSR/ZZ4lAnalysis_ZZTo4tau.root");
  files.push_back("PRODFSR/ZZ4lAnalysis_ggZZ2l2l.root");
  files.push_back("PRODFSR/ZZ4lAnalysis_ggZZ4l.root");
  files.push_back("PRODFSR/ZZ4lAnalysis_DYJetsToLLTuneZ2M10B.root");
  files.push_back("PRODFSR/ZZ4lAnalysis_DYJetsToLLTuneZ2M10NoB.root");
  files.push_back("PRODFSR/ZZ4lAnalysis_DYJetsToLLTuneZ2M50B.root");
  files.push_back("PRODFSR/ZZ4lAnalysis_DYJetsToLLTuneZ2M50NoB.root");
  files.push_back("PRODFSR/ZZ4lAnalysis_TTTo2L2Nu2B.root");
  //files.push_back("PRODFSR/ZZ4lAnalysis_WZ.root");


  TString finalStates[] = {"4mu", "4e", "2e2mu"};

  for (vector<TString>::const_iterator file = files.begin(); file!=files.end(); ++file) {
    TString filename = outDir + *file;
    cout << "Processing " << filename << endl;
    TString mkdircmd = "mkdir -p `dirname ";
    gSystem->Exec(mkdircmd+filename+"`");    

    TFile* outputFile = TFileServiceWrapper::init(filename);
    

    for (int i=0; i<3; ++i) {      
      TFile* input = TFile::Open(inputDir+(*file));
      TString plotPath = "Plots" + finalStates[i];
      cout << " " << plotPath <<endl;
      TFileServiceWrapper::mkdir(plotPath);
      TString chainName = "ZZ" + finalStates[i]+"Tree/candTree";
      TChain *tree = new TChain(chainName);
      tree->Add(inputDir+(*file));
      ZZTreeReader reader(tree);
      reader.Loop();
      TH1F* nEventComplete = (TH1F*) (input->Get(plotPath+"/nEventComplete"))->Clone();
      input->Close();
      nEventComplete->Write();
    }

    outputFile->Close();
  }
}
