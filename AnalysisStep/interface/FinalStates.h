#ifndef FinalStates_h
#define FinalStates_h
// Some standard labels for final states.

#include <string>

enum Channel {MMMM=0, 
	      EEEE=1, 
	      EEMM=2, 
	      ZZ=3,
	      LLTT=4, // tau final states
	      CRMMMMss=5,
	      CRMMMMos=6,
	      CREEEEss=7,
	      CREEEEos=8,
	      CREEMMss=9,
	      CREEMMos=10,
	      CRMMEEss=11,
	      CRMMEEos=12,
	      ZL=13,
	      CRZLLHiSIP=14,
	      CRZLLHiSIPMM=15,
	      CRZLLHiSIPKin=16,
	      CRZLL=17,
	      ZLL=18,
	      NONE = 99, BUGGY=666};

//Return string corresponding to integer code 
std::string finalState(int iFS);

//Return "4e", "4mu, "2e2mu
std::string finalStateNiceName(int iFS);

//Return integer code corresponding to name
Channel finalState(std::string sFS);

#endif
