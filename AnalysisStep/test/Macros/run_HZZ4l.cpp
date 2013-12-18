#include "HZZ4l.h"

#include <TH1F.h>

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <string>


int main (int argc, char ** argv) 
{
  int theChannel(atoi(argv[1]));

  std::string channelname;

  switch(theChannel){

  case 0:
    channelname = "ZZ4muTree/";
    break;

  case 1:
    channelname = "ZZ4eTree/";
    break;

  case 2:
    channelname = "ZZ2e2muTree/";
    break;

  default:
    std::cout << "Need to specify the correct final state! " << std::endl;
    abort();
  }


  //Get the TTree with the candidates
  std::string chainName = channelname + "candTree";
  TChain *tree = new TChain(chainName.c_str());

  std::string samplename(argv[3]);
  std::string inputfilename(argv[4]);
  tree->Add(inputfilename.c_str());

  std::string outputfilename(argv[5]);

  HZZ4l analyzer(tree, samplename);

  //Set sqrt(s)
  bool is8TeV(atoi(argv[2]));
  analyzer.set8TeV(is8TeV);

  //check if we want to do non-SM shapes
  // 0 means SM
  // 1 means EWK-singlet
  // 2 means 2HDM
  int BSM_flag = 0;
  if(outputfilename.find("kappa") != std::string::npos) BSM_flag = 1;
  if(outputfilename.find("2HDM") != std::string::npos)  BSM_flag = 2;
  analyzer.setBSM(BSM_flag);

  //Get the normalization. It's the first bin of the histo
  TFile fIn(inputfilename.c_str());
  std::string histoName = channelname + "Counters";
  TH1F *nEventComplete = (TH1F*)fIn.Get(histoName.c_str());
  Int_t Nevt_Gen = nEventComplete->GetBinContent(1);
  Float_t Nevt_Gen_weighted = nEventComplete->GetBinContent(0); // Weighted #of events, can be used only when no filter was applied

  //SM case
  if(BSM_flag == 0){
    std::cout << "################################################" << std::endl;
    analyzer.Loop(theChannel, outputfilename.c_str());
  }

  //EWK-singlet
  else if(BSM_flag == 1){
    std::cout << "Doing EWK singlet root files!" << std::endl;
    std::string::size_type pos = outputfilename.find("kappa");

    for(float kappa = 0.2; kappa < 1.1 ; kappa = kappa + 0.2){
	std::cout << "################################################" << std::endl;
      std::cout << "Now doing kappa' = " << kappa << std::endl;
      std::string outputfilename_tmp = outputfilename;
      std::stringstream ss (std::stringstream::in | std::stringstream::out);
      ss << kappa;
      std::string replace_string = "kappa_" + ss.str();
      if(pos != std::string::npos) outputfilename_tmp.replace(pos, 5, replace_string); 
      std::cout << outputfilename << " -> " << outputfilename_tmp << std::endl;
      analyzer.setKappa(kappa);
      analyzer.Loop(theChannel, outputfilename_tmp.c_str());
    }

  }


  return 0;
}
