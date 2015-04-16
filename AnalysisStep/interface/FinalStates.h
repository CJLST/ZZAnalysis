#ifndef FinalStates_h
#define FinalStates_h
// Some standard labels for final states.

#include <string>

enum Channel {MMMM=0, 
	      EEEE=1, 
	      EEMM=2, 
	      ZZ=3,
	      LLTT=4, // tau final states
	      ZLL=18, // Generic label for a CR (used in job configuration)
	      CRMMMMss=5, // AA CRs (ss loose letptons; os ones are for x-check)
	      CRMMMMos=6,
	      CREEEEss=7,
	      CREEEEos=8,
	      CREEMMss=9,
	      CREEMMos=10,
	      CRMMEEss=11,
	      CRMMEEos=12,
	      CRZLLss=21,  //merged AA ss CRs
	      ZL=13, // Fake rate CR (Z+loose lepton)
	      CRZLLHiSIP=14,   // Old inverted SIP CRs
	      CRZLLHiSIPMM=15,
	      CRZLLHiSIPKin=16,
	      CRZLL=17,   // Old CR for Z2 with no SIP
	      CRZLLos_2P2F=22,  //CR for Z+2l opposite sign, 2P2F
	      CRZLLos_3P1F=23,
	      NONE = 99, BUGGY=666};

//Return string corresponding to integer code 
std::string finalState(int iFS);

//Return "4e", "4mu, "2e2mu
std::string finalStateNiceName(int iFS);

//Return integer code corresponding to name
Channel finalState(std::string sFS);

#endif
