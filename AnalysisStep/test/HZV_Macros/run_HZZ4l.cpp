#include "HZZ4l.h"

#include <TH1F.h>

#include <iostream>
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

  bool is8TeV(atoi(argv[2]));

  std::string filename(argv[3]);
  tree->Add(filename.c_str());

  std::string outputfilename(argv[4]);

  bool isVH(false);
  std::vector<std::string> VHfinalStates;
  VHfinalStates.push_back("WH");
  VHfinalStates.push_back("ZH");
  VHfinalStates.push_back("ttH");
  if (filename.find("VH")!=std::string::npos) isVH=true;  

  //Get the normalization. It's the first bin of the histo
  TFile fIn(filename.c_str());
  std::string histoName = channelname + "Counters";
  TH1F *nEventComplete = (TH1F*)fIn.Get(histoName.c_str());
  Int_t Nevt_Gen = nEventComplete->GetBinContent(1);

  HZZ4l analyzer(tree);
  if (!isVH) {
    std::cout << "################################################" << std::endl;
    std::cout << "Name of the output file : " << outputfilename << std::endl;
    std::cout << "Number of total entries : " << tree->GetEntries() << std::endl ;
    std::cout << "Number of generated events : " << Nevt_Gen << std::endl;
    analyzer.Loop(Nevt_Gen, theChannel, outputfilename.c_str(), is8TeV, 0);
  }
  else{
    std::string::size_type pos = outputfilename.find("VH");
    for (unsigned int iVH=0; iVH<VHfinalStates.size(); ++iVH){
      std::cout<<filename<<" -> "<<outputfilename<<VHfinalStates[iVH]<<"  "<<pos<<"  "<<std::endl;
      if ( pos != std::string::npos ) outputfilename.replace( pos, 2, VHfinalStates[iVH] ); 
      std::cout << "################################################" << std::endl;
      std::cout << "Name of the output file : " << outputfilename << std::endl;
      std::cout << "Number of total entries : " << tree->GetEntries() << std::endl ;
      std::cout << "Number of generated events : " << Nevt_Gen << std::endl;
      analyzer.Loop(Nevt_Gen, theChannel, outputfilename.c_str(), is8TeV, iVH+1);
    }
  }

  return 0;
}
