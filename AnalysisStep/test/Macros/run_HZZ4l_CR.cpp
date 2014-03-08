#include "HZZ4l.h"

#include <TH1F.h>

#include <iostream>
#include <cstdlib>
#include <string>

int main (int argc, char ** argv) 
{
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

  //Let the analyzer know this is a Control Region
  analyzer.setCR(true);

  tree_CRZLLTree->Add(inputfilename.c_str());

  bool HZZ4l_flag=false;
  bool HZZ4l_NoLepInt=false;
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
  if(outputfilename.find( "MCFM67" ) != std::string::npos) HZZ4l_NoLepInt=true;
  if(HZZ4l_NoLepInt) cout << "Sample with need for ggZZ lepton interference is found" << endl;
  analyzer.setHZZ4l_NoLepInt(HZZ4l_NoLepInt);

  bool qqZZ_flag=false;
  bool ggZZ_flag=false;
  if( 
	  (
		  outputfilename.find( "ggZZ4l" ) != std::string::npos 
		  || outputfilename.find( "ggZZ2l2l" ) != std::string::npos
		  || outputfilename.find( "ggTo4l_Continuum" ) != std::string::npos
		  || outputfilename.find( "ggTo2l2l_Continuum" ) != std::string::npos
		  || outputfilename.find( "ggTo4e_Contin" ) != std::string::npos
		  || outputfilename.find( "ggTo4mu_Contin" ) != std::string::npos
		  || outputfilename.find( "ggTo2e2mu_Contin" ) != std::string::npos
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
  std::string histoName = "CRZLLTree/Counters";
  TH1F *nEventComplete = (TH1F*)fIn.Get(histoName.c_str());
  Int_t Nevt_Gen = nEventComplete->GetBinContent(1);

  std::string fileOut_CRZLLTree;     
 
  std::cout << "################################################" << std::endl;
  analyzer.Loop(0, outputfilename);

  return 0;
}
