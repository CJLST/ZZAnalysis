#include "HZZ4l.h"

#include <TH1F.h>

#include <iostream>
#include <cstdlib>
#include <string>

int main (int argc, char ** argv) 
{
  TChain *tree_CRZLLTree = new TChain("CRZLLTree/crTree");
  TChain *tree_CRZLLHiSIPTree = new TChain("CRZLLHiSIPTree/crTree");

  TChain *tree_CRMMMMssTree = new TChain("CRMMMMssTree/crTree");
  TChain *tree_CREEEEssTree = new TChain("CREEEEssTree/crTree");
  TChain *tree_CREEMMssTree = new TChain("CREEMMssTree/crTree");
  TChain *tree_CRMMEEssTree = new TChain("CRMMEEssTree/crTree");

  TChain *tree_CRMMMMosTree = new TChain("CRMMMMosTree/crTree");
  TChain *tree_CREEEEosTree = new TChain("CREEEEosTree/crTree");
  TChain *tree_CREEMMosTree = new TChain("CREEMMosTree/crTree");
  TChain *tree_CRMMEEosTree = new TChain("CRMMEEosTree/crTree");

  bool is8TeV(atoi(argv[1]));

  //Get the TTree with the candidates
  std::string filename(argv[2]);

  std::string outputfilename(argv[3]);

  bool isVH(false);
  std::vector<std::string> VHfinalStates;
  VHfinalStates.push_back("WH");
  VHfinalStates.push_back("ZH");
  VHfinalStates.push_back("ttH");
  if (filename.find("VH")!=std::string::npos) isVH=true;  


  tree_CRZLLTree->Add(filename.c_str());
  tree_CRZLLHiSIPTree->Add(filename.c_str());

  tree_CRMMMMssTree->Add(filename.c_str());
  tree_CREEEEssTree->Add(filename.c_str());
  tree_CREEMMssTree->Add(filename.c_str());
  tree_CRMMEEssTree->Add(filename.c_str());

  tree_CRMMMMosTree->Add(filename.c_str());
  tree_CREEEEosTree->Add(filename.c_str());
  tree_CREEMMosTree->Add(filename.c_str());
  tree_CRMMEEosTree->Add(filename.c_str());

  //Get the normalization. It's the first bin of the histo
  TFile fIn(filename.c_str());
  std::string histoName = "CRZLLTree/Counters";
  TH1F *nEventComplete = (TH1F*)fIn.Get(histoName.c_str());
  Int_t Nevt_Gen = nEventComplete->GetBinContent(1);

  HZZ4l analyzer_CRZLLTree(tree_CRZLLTree);
  HZZ4l analyzer_CRZLLHiSIPTree(tree_CRZLLHiSIPTree);

  HZZ4l analyzer_CRMMMMssTree(tree_CRMMMMssTree);
  HZZ4l analyzer_CREEEEssTree(tree_CREEEEssTree);
  HZZ4l analyzer_CREEMMssTree(tree_CREEMMssTree);
  HZZ4l analyzer_CRMMEEssTree(tree_CRMMEEssTree);

  HZZ4l analyzer_CRMMMMosTree(tree_CRMMMMosTree);
  HZZ4l analyzer_CREEEEosTree(tree_CREEEEosTree);
  HZZ4l analyzer_CREEMMosTree(tree_CREEMMosTree);
  HZZ4l analyzer_CRMMEEosTree(tree_CRMMEEosTree);


  std::string fileOut_CRZLLTree;     
  std::string fileOut_CRZLLHiSIPTree;

  std::string fileOut_CRMMMMssTree;  
  std::string fileOut_CREEEEssTree;  
  std::string fileOut_CREEMMssTree;  
  std::string fileOut_CRMMEEssTree;  

  std::string fileOut_CRMMMMosTree;  
  std::string fileOut_CREEEEosTree;  
  std::string fileOut_CREEMMosTree;  
  std::string fileOut_CRMMEEosTree;  

  if (!isVH) {
    fileOut_CRZLLTree      = outputfilename + "_CRZLLTree.root"; 
    /*    
    fileOut_CRZLLHiSIPTree = outputfilename + "_CRZLLHiSIPTree.root"; 
    
    fileOut_CRMMMMssTree   = outputfilename + "_CRMMMMssTree.root";   
    fileOut_CREEEEssTree   = outputfilename + "_CREEEEssTree.root";   
    fileOut_CREEMMssTree   = outputfilename + "_CREEMMssTree.root";   
    fileOut_CRMMEEssTree   = outputfilename + "_CRMMEEssTree.root";   
    
    fileOut_CRMMMMosTree   = outputfilename + "_CRMMMMosTree.root";   
    fileOut_CREEEEosTree   = outputfilename + "_CREEEEosTree.root";   
    fileOut_CREEMMosTree   = outputfilename + "_CREEMMosTree.root";   
    fileOut_CRMMEEosTree   = outputfilename + "_CRMMEEosTree.root";       
    */
    std::cout << "################################################" << std::endl;
    std::cout << "Name of the output file : " << outputfilename << std::endl;
    std::cout << "Number of generated events : " << Nevt_Gen << std::endl;

    analyzer_CRZLLTree.Loop(Nevt_Gen, 0, fileOut_CRZLLTree.c_str(), is8TeV, 0);
    /*    analyzer_CRZLLHiSIPTree.Loop(Nevt_Gen, 0, fileOut_CRZLLHiSIPTree.c_str(), is8TeV, 0);
    
    analyzer_CRMMMMssTree.Loop(Nevt_Gen, 0, fileOut_CRMMMMssTree.c_str(), is8TeV, 0);
    analyzer_CREEEEssTree.Loop(Nevt_Gen, 1, fileOut_CREEEEssTree.c_str(), is8TeV, 0);
    analyzer_CREEMMssTree.Loop(Nevt_Gen, 2, fileOut_CREEMMssTree.c_str(), is8TeV, 0);
    analyzer_CRMMEEssTree.Loop(Nevt_Gen, 2, fileOut_CRMMEEssTree.c_str(), is8TeV, 0);
    
    analyzer_CRMMMMosTree.Loop(Nevt_Gen, 0, fileOut_CRMMMMosTree.c_str(), is8TeV, 0);
    analyzer_CREEEEosTree.Loop(Nevt_Gen, 1, fileOut_CREEEEosTree.c_str(), is8TeV, 0);
    analyzer_CREEMMosTree.Loop(Nevt_Gen, 2, fileOut_CREEMMosTree.c_str(), is8TeV, 0);
    analyzer_CRMMEEosTree.Loop(Nevt_Gen, 2, fileOut_CRMMEEosTree.c_str(), is8TeV, 0);*/
  }
  else{
    std::string::size_type pos = outputfilename.find("VH");
    for (unsigned int iVH=0; iVH<VHfinalStates.size(); ++iVH){
      std::cout<<filename<<" -> "<<outputfilename<<VHfinalStates[iVH]<<"  "<<pos<<"  "<<std::endl;
      if ( pos != std::string::npos ) outputfilename.replace( pos, 2, VHfinalStates[iVH] ); 

      fileOut_CRZLLTree      = outputfilename + "_CRZLLTree.root";     
      /*      fileOut_CRZLLHiSIPTree = outputfilename + "_CRZLLHiSIPTree.root"; 
      
      fileOut_CRMMMMssTree   = outputfilename + "_CRMMMMssTree.root";   
      fileOut_CREEEEssTree   = outputfilename + "_CREEEEssTree.root";   
      fileOut_CREEMMssTree   = outputfilename + "_CREEMMssTree.root";   
      fileOut_CRMMEEssTree   = outputfilename + "_CRMMEEssTree.root";   
      
      fileOut_CRMMMMosTree   = outputfilename + "_CRMMMMosTree.root";   
      fileOut_CREEEEosTree   = outputfilename + "_CREEEEosTree.root";   
      fileOut_CREEMMosTree   = outputfilename + "_CREEMMosTree.root";   
      fileOut_CRMMEEosTree   = outputfilename + "_CRMMEEosTree.root";     */  

      std::cout << "################################################" << std::endl;
      std::cout << "Name of the output file : " << outputfilename << std::endl;
      std::cout << "Number of generated events : " << Nevt_Gen << std::endl;
      
      analyzer_CRZLLTree.Loop(Nevt_Gen, 0, fileOut_CRZLLTree.c_str(), is8TeV, iVH+1);
      /*     analyzer_CRZLLHiSIPTree.Loop(Nevt_Gen, 0, fileOut_CRZLLHiSIPTree.c_str(), is8TeV, iVH+1);
      
      analyzer_CRMMMMssTree.Loop(Nevt_Gen, 0, fileOut_CRMMMMssTree.c_str(), is8TeV, iVH+1);
      analyzer_CREEEEssTree.Loop(Nevt_Gen, 1, fileOut_CREEEEssTree.c_str(), is8TeV, iVH+1);
      analyzer_CREEMMssTree.Loop(Nevt_Gen, 2, fileOut_CREEMMssTree.c_str(), is8TeV, iVH+1);
      analyzer_CRMMEEssTree.Loop(Nevt_Gen, 2, fileOut_CRMMEEssTree.c_str(), is8TeV, iVH+1);
      
      analyzer_CRMMMMosTree.Loop(Nevt_Gen, 0, fileOut_CRMMMMosTree.c_str(), is8TeV, iVH+1);
      analyzer_CREEEEosTree.Loop(Nevt_Gen, 1, fileOut_CREEEEosTree.c_str(), is8TeV, iVH+1);
      analyzer_CREEMMosTree.Loop(Nevt_Gen, 2, fileOut_CREEMMosTree.c_str(), is8TeV, iVH+1);
      analyzer_CRMMEEosTree.Loop(Nevt_Gen, 2, fileOut_CRMMEEosTree.c_str(), is8TeV, iVH+1);*/
    }
  }


  return 0;
}
