#include "HZZ4l.h"

#include <TH1F.h>
#include <TH2F.h>

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
  // Should match ZZAnalysis...root names
	char* myPrimarySample_SpinZero[kNumFiles]={
		"jhuGenV3H126",
		"0PH126",
		"0MH126",
		"0Mf05ph0H126",
		"0Mf05ph90H126",
		"0Mf05ph180H126",
		"0Mf05ph270H126",
		"0Mf01ph0H126", 
		"0Mf01ph90H126",
		"0Mf01ph180H126", 
		"0Mf01ph270H126",

		"0PMH125.6",
		"0PH125.6",
		"0MH125.6",
		"0L1H125.6",

		"0PHf05ph0H125.6",
		"0Mf05ph0H125.6",
		"0PHf05ph0Mf05ph0H125.6",
		"0L1f05ph180H125.6",

		"0PHf033ph0Mf033ph0H125.6",
		"0PHf01ph0H125.6",
		"0Mf01ph0H125.6",
		"0PHf01ph0Mf01ph0H125.6",

		"0L1f05ph0H125.6",
		"0L1f01ph0H125.6",

		"0PHf05ph90H125.6",
		"0Mf05ph90H125.6",
		"0PHf05ph0Mf05ph90H125.6",

		"0PHf033ph0Mf033ph90H125.6",
		"0PHf01ph90H125.6",
		"0Mf01ph90H125.6",
		"0PHf01ph0Mf01ph90H125.6",

		"0PHf05ph180H125.6",
		"0Mf05ph180H125.6",
		"0PHf05ph180Mf05ph0H125.6"
	};



  //Get the TTree with the candidates
  std::string chainName = channelname + "candTree";
  TChain *tree = new TChain(chainName.c_str());

  std::string samplename(argv[3]);
  std::string inputfilename(argv[4]);
  tree->Add(inputfilename.c_str());
  // tree->Show(15); 

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

  bool HZZ4l_flag=false;
  int HZZ4l_code=0;
  float HZZ4L_HMass=125.6;
  for(int f=0;f<kNumFiles;f++){
	  if( outputfilename.find( myPrimarySample_SpinZero[f] ) != std::string::npos){
		  HZZ4l_code=f;
		  HZZ4l_flag=true;
		  if(f<11) HZZ4L_HMass=126;
	  };
  };
  if(HZZ4l_flag) cout << "HZZ4l spin code is " << HZZ4l_code << endl;
  analyzer.setHZZ4l(HZZ4l_flag,HZZ4l_code,HZZ4L_HMass);

  bool qqZZ_flag=false;
  bool ggZZ_flag=false;
  if( 
	  (
		  outputfilename.find( "ggZZ4l" ) != std::string::npos 
		  || outputfilename.find( "ggZZ2l2l" ) != std::string::npos
		  || outputfilename.find( "ggTo4l_Continuum" ) != std::string::npos
		  || outputfilename.find( "ggTo2l2l_Continuum" ) != std::string::npos
	  ) && !(
		  outputfilename.find( "ggTo4l_H125.6" ) != std::string::npos 
		  || outputfilename.find( "ggTo2l2l_H125.6" ) != std::string::npos
		  || outputfilename.find( "ggTo4l_ContinuumInterfH125.6" ) != std::string::npos
		  || outputfilename.find( "ggTo2l2l_ContinuumInterfH125.6" ) != std::string::npos
	  )
	) ggZZ_flag=true;
  if( 
	  (
		  outputfilename.find( "_ZZTo" ) != std::string::npos 
		  || outputfilename.find( "_ZZ95-160To" ) != std::string::npos
	  )
	) qqZZ_flag=true;
  analyzer.setGGQQB(qqZZ_flag,ggZZ_flag);

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
