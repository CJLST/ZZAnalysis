#ifndef FinalStates_h
#define FinalStates_h
// Some standard labels for final states.

#include <string>

enum Channel {MMMM=0, 
	      EEEE=1, 
	      EEMM=2, 
	      ZZ=3,
	      LLTT=4, // tau final states (L=e,mu,tau); obsolete.
	      llTT=27,// tau final states (L=e,mu)
	      TTTT=28,// 4tau
	      ZLL=18, // Generic label for a CR (used in job configuration)
	      ZZOnShell = 26,
	      CRMMMMss=5, // AA CRs (ss loose letptons; os ones are for x-check)
	      CRMMMMos=6,
	      CREEEEss=7,
	      CREEEEos=8,
	      CREEMMss=9,
	      CREEMMos=10,
	      CRMMEEss=11,
	      CRMMEEos=12,
	      CRZLLss=21,  //merged AA ss CRs
	      CRZLLos_2P2F=22,  //CR for Z+2l opposite sign, 2P2F
	      CRZLLos_3P1F=23,  //CR for Z+2l opposite sign, 3P1F
	      CRZLLos_2P2F_ZZOnShell=24,  //CR for Z+2l opposite sign, 2P2F, 2 on-shell Zs
	      CRZLLos_3P1F_ZZOnShell=25,  //CR for Z+2l opposite sign, 3P1F, 2 on-shell Zs
	      ZL=13, // Fake rate CR (Z+loose lepton)
	      CRZLLHiSIP=14,   // Old inverted SIP CRs
	      CRZLLHiSIPMM=15,
	      CRZLLHiSIPKin=16,
	      CRZLL=17,   // Old CR for Z2 with no SIP
	      CRZ2mLL=19, // A CR: Z->mumu + l+l- (with at least 1F)
	      CRZ2eLL=20, // A CR: Z->ee   + l+l- (with at least 1F)
	      NONE = 99, BUGGY=666};

//Return string corresponding to integer code 
std::string finalState(int iFS);

//Return "4e", "4mu, "2e2mu
std::string finalStateNiceName(int iFS);

//Return integer code corresponding to name
Channel finalState(std::string sFS);

#endif
