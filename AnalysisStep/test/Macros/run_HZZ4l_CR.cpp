#include "HZZ4l.h"

#include <TH1F.h>

#include <iostream>
#include <cstdlib>
#include <string>

int main (int argc, char ** argv) 
{
// Should match ZZAnalysis...root names

	TChain *tree_CRZLLTree = new TChain("CRZLLTree/crTree");

  //Get the TTree with the candidates
  std::string samplename(argv[2]);
  std::string inputfilename(argv[3]);
  std::string outputfilename(argv[4]);

  HZZ4l analyzer(tree_CRZLLTree, samplename);

  //Set sqrt(s)
  bool is8TeV(atoi(argv[1]));
  analyzer.set8TeV(is8TeV);

  //check if we want to do non-SM shapes
  // 0 means SM
  // 1 means EWK-singlet
  // 2 means 2HDM
  int BSM_flag = 0;
  if(outputfilename.find("kappa") != std::string::npos) BSM_flag = 1;
  if(outputfilename.find("2HDM") != std::string::npos)  BSM_flag = 2;
  analyzer.setBSM(BSM_flag);

  analyzer.identifySample(outputfilename);

  //Let the analyzer know this is a Control Region
  analyzer.setCR(true);

  tree_CRZLLTree->Add(inputfilename.c_str());

  //Get the normalization. It's the first bin of the histo
  TFile fIn(inputfilename.c_str());
  std::string histoName = "CRZLLTree/Counters";
  TH1F *nEventComplete = (TH1F*)fIn.Get(histoName.c_str());
  Int_t Nevt_Gen = nEventComplete->GetBinContent(1);

  std::string fileOut_CRZLLTree;     
 
  std::cout << "################################################" << std::endl;
  analyzer.Loop(0, outputfilename);

  return 0;
}
