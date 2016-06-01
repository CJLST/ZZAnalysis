#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>

namespace {
  const unsigned nFS = 29;
}


std::string finalState(int iFS) {
  if (iFS<0 || iFS>=nFS) return "None";
  const std::string finalStates[nFS] = {"MMMM",      // 0
					"EEEE",      // 1
					"EEMM",      // 2
					"ZZ",        // 3
					"LLTT",      // 4
					"CRMMMMss",  // 5 
					"CRMMMMos",  // 6 
					"CREEEEss",  // 7 
					"CREEEEos",  //	8
					"CREEMMss",  // 9
					"CREEMMos",  // 10
					"CRMMEEss",  // 11
					"CRMMEEos",  // 12
					"ZL",        // 13
					"CRZLLHiSIP",   //14
					"CRZLLHiSIPMM", //15
					"CRZLLHiSIPKin",//16
					"CRZLL",     //17
					"ZLL",       //18
					"CRZ2mLL",   //19
					"CRZ2eLL",   //20
					"CRZLLss",   //21
					"CRZLLos_2P2F", //22
					"CRZLLos_3P1F", //23
					"CRZLLos_2P2F_ZZOnShell", //24
					"CRZLLos_3P1F_ZZOnShell", //25
					"ZZOnShell", //26
					"llTT",      //27
					"TTTT"       //28
  };
  return finalStates[iFS];			     	
}


std::string finalStateNiceName(int iFS) {
  if (iFS<0||iFS>2) return finalState(iFS);
  const std::string finalStates[3] = {"4mu", "4e", "2e2mu"};
  
  return finalStates[iFS];
}



Channel finalState(std::string sFS) {
  for (unsigned i=0; i<nFS; ++i) {
    if (sFS==finalState(i)) return (Channel) i;
  }
  return NONE;
}
