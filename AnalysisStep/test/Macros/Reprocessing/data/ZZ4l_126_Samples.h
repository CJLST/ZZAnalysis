#ifndef ZZ4L_126_SAMPLES_H
#define ZZ4L_126_SAMPLES_H

#include <iostream>
#include <string>
#include "TMath.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TCanvas.h"

const float PI_VAL = TMath::Pi();
const int kNumFiles=11;
const int nFinalStates=5;

enum sample {
	kfg2_0_fg4_0,
	kfg2_1_fg4_0,
	kfg2_0_fg4_1,
	kfLambda1_1,

	kfg2_05_fg4_0,
	kfg2_0_fg4_05,
	kfg2_05_fg4_05,
	kfLambda1_m05,

	kfg2_33_fg4_33,
	kfg2_01_fg4_0,
	kfg2_0_fg4_01,
	kfg2_01_fg4_01,

	kfLambda1_05,
	kfLambda1_03, // No sample here
	kfLambda1_01,

	kfg2_05_fg4_0_p290,
	kfg2_0_fg4_05_p390,
	kfg2_05_fg4_05_p390,

	kfg2_33_fg4_33_p390,
	kfg2_01_fg4_0_p290,
	kfg2_0_fg4_01_p390,
	kfg2_01_fg4_01_p390,

	kfg2_05_fg4_0_p2Pi,
	kfg2_0_fg4_05_p3Pi,
	kfg2_05_fg4_05_p2Pi,

	kNumSamples
};
string sampleName[kNumSamples] = {
	"fg2_0_fg4_0",
	"fg2_1_fg4_0",
	"fg2_0_fg4_1",
	"fLambda1_1",

	"fg2_05_fg4_0",
	"fg2_0_fg4_05",
	"fg2_05_fg4_05",
	"fLambda1_m05",

	"fg2_33_fg4_33",
	"fg2_01_fg4_0",
	"fg2_0_fg4_01",
	"fg2_01_fg4_01",

	"fLambda1_05",
	"fLambda1_03", // No sample here
	"fLambda1_01",

	"fg2_05_fg4_0_p290",
	"fg2_0_fg4_05_p390",
	"fg2_05_fg4_05_p390",

	"fg2_33_fg4_33_p390",
	"fg2_01_fg4_0_p290",
	"fg2_0_fg4_01_p390",
	"fg2_01_fg4_01_p390",

	"fg2_05_fg4_0_p2Pi",
	"fg2_0_fg4_05_p3Pi",
	"fg2_05_fg4_05_p2Pi",
};
char* sampleName_label[kNumSamples]={
	"0^{+}_{m} 126 GeV",
	"0^{+}_{h} 126 GeV",
	"0^{-} 126 GeV",
	"0^{+}_{m} Pure Q^{2} 126 GeV",

	"f_{a2}=0.5 #phi_{a2}=0 126 GeV",
	"f_{a3}=0.5 #phi_{a3}=0 126 GeV",
	"f_{a2}=f_{a3}=0.5 #phi_{a2}=#phi_{a3}=0 126 GeV",
	"f_{#Lambda1}=-0.5 126 GeV",

	"f_{a2}=f_{a3}=0.33 #phi_{a2}=#phi_{a3}=0 126 GeV",
	"f_{a2}=0.1 #phi_{a2}=0 126 GeV",
	"f_{a3}=0.1 #phi_{a3}=0 126 GeV",
	"f_{a2}=f_{a3}=0.1 #phi_{a2}=#phi_{a3}=0 126 GeV",

	"f_{#Lambda1}=0.5 126 GeV",
	"f_{#Lambda1}=0.3 126 GeV",
	"f_{#Lambda1}=0.1 126 GeV",

	"f_{a2}=0.5 #phi_{a2}=#pi/2 126 GeV",
	"f_{a3}=0.5 #phi_{a3}=#pi/2 126 GeV",
	"f_{a2}=f_{a3}=0.5 #phi_{a3}=#pi/2 126 GeV",

	"f_{a2}=0.1 #phi_{a2}=#pi/2 126 GeV",
	"f_{a3}=0.1 #phi_{a3}=#pi/2 126 GeV",
	"f_{a2}=f_{a3}=0.1 #phi_{a3}=#pi/2 126 GeV",

	"f_{a2}=0.5 #phi_{a2}=#pi 126 GeV",
	"f_{a3}=0.5 #phi_{a3}=#pi 126 GeV",
	"f_{a2}=f_{a3}=0.5 #phi_{a2}=#pi 126 GeV"
	"f_{a2}=f_{a3}=0.33 #phi_{a3}=#pi/2 126 GeV",
};
char* sample_file[kNumFiles]={
	"HZZ4lTree_powheg15jhuGenV3H126",
	"HZZ4lTree_powheg15jhuGenV3ScaHH126",
	"HZZ4lTree_powheg15jhuGenV3PseHH126",
	"HZZ4lTree_powheg15jhuGenV3-0Mf05ph0H126",
	"HZZ4lTree_powheg15jhuGenV3-0Mf05ph90H126",
	"HZZ4lTree_powheg15jhuGenV3-0Mf05ph180H126",
	"HZZ4lTree_powheg15jhuGenV3-0Mf05ph270H126",
	"HZZ4lTree_powheg15jhuGenV3-0Mf01ph0H126", 
	"HZZ4lTree_powheg15jhuGenV3-0Mf01ph90H126",
	"HZZ4lTree_powheg15jhuGenV3-0Mf01ph180H126", 
	"HZZ4lTree_powheg15jhuGenV3-0Mf01ph270H126"
};
char* sample_file_label[kNumFiles]={
	"0^{+}_{m} 126 GeV",
	"0^{+}_{h} 126 GeV",
	"0^{-} 126 GeV",
	"f_{a3}=0.5 #phi_{a3}=0 126 GeV",
	"f_{a3}=0.5 #phi_{a3}=#pi/2 126 GeV",
	"f_{a3}=0.5 #phi_{a3}=#pi 126 GeV",
	"f_{a3}=0.5 #phi_{a3}=3#pi/2 126 GeV",
	"f_{a3}=0.1 #phi_{a3}=0 126 GeV",
	"f_{a3}=0.1 #phi_{a3}=#pi/2 126 GeV",
	"f_{a3}=0.1 #phi_{a3}=#pi 126 GeV",
	"f_{a3}=0.1 #phi_{a3}=3#pi/2 126 GeV"
};
char* sample_filePrimary[kNumFiles]={
	"ZZ4lAnalysis_powheg15jhuGenV3H126",
	"ZZ4lAnalysis_powheg15jhuGenV3ScaHH126",
	"ZZ4lAnalysis_powheg15jhuGenV3PseHH126",
	"ZZ4lAnalysis_powheg15jhuGenV3-0Mf05ph0H126",
	"ZZ4lAnalysis_powheg15jhuGenV3-0Mf05ph90H126",
	"ZZ4lAnalysis_powheg15jhuGenV3-0Mf05ph180H126",
	"ZZ4lAnalysis_powheg15jhuGenV3-0Mf05ph270H126",
	"ZZ4lAnalysis_powheg15jhuGenV3-0Mf01ph0H126", 
	"ZZ4lAnalysis_powheg15jhuGenV3-0Mf01ph90H126",
	"ZZ4lAnalysis_powheg15jhuGenV3-0Mf01ph180H126", 
	"ZZ4lAnalysis_powheg15jhuGenV3-0Mf01ph270H126"
};
const float gi_phi2_phi4_files[kNumFiles][9]={
	{	1.0,	0,	0,	0,		0,	0,	0,	0,	0},
	{	0,	1.0,	0,	0,		0,	0,	1.0,	0,	0},
	{	0,	0,	0,	1.0,		0,	0,	0,	1.0,	0},

	{	1.0,	0,	0,	2.498,		0,	0,	0,	0.5,	0},
	{	1.0,	0,	0,	2.498,		0,	PI_VAL/2.0,	0,	0.5,	0},
	{	1.0,	0,	0,	2.498,		0,	PI_VAL,	0,	0.5,	0},
	{	1.0,	0,	0,	2.498,		0,	PI_VAL*1.5,	0,	0.5,	0},

	{	1.0,	0,	0,	0.8327,		0,	0,	0,	0.1,	0},
	{	1.0,	0,	0,	0.8327,		0,	PI_VAL/2.0,	0,	0.1,	0},
	{	1.0,	0,	0,	0.8327,		0,	PI_VAL,	0,	0.1,	0},
	{	1.0,	0,	0,	0.8327,		0,	PI_VAL*1.5,	0,	0.1,	0}
};

double sample_xsec_ratio[11]={ // Assuming g_i below
	1.0,
	0.379,
	0.1608,
	6.94081E-09, // g1=1, g1'=-1, lambda=10e4
	3.698,	2.0,	0.3244,	// phase 0 //
	3.935147, // g1=12004.14, g1'=-12003.14, lambda=10e4
	3.698,	2.0,	0.3244 // phase pi/2 //
};
double epspr_mixedSample_lambda1pr=10000.0;

const float gi_phi2_phi4[kNumSamples][9]={ // g1-4; phia2,3; fa2, 3; fepspr
	{	1.0,	0,	0,	0,		0,	0,	0,	0,	0}, // Pure SM
	{	0,	1.0,	0,	0,		0,	0,	1.0,	0,	0}, // fa2=1
	{	0,	0,	0,	1.0,		0,	0,	0,	1.0,	0}, // fa3=1
	{	0,	0,	0,	0,		0,	0,	0,	0,	1.0}, // fL1=1

	{	1.0,	1.624,	0,	0,		0,	0,	0.5,	0,	0}, // fa2=0.5
	{	1.0,	0,	0,	2.484,		0,	0,	0,	0.5,	0}, // fa3=0.5
	{	0,	0.654,	0,	1.0,		0,	0,	0.5,	0.5,	0}, // fa2=fa3=0.5
	{	1.0,	0,	0,	0,		0,	0,	0,	0,	12003.14}, // fLambda1=-0.5, for T3 templates

	{	1.0,	1.624,	0,	2.484,		0,	0,	1.0/3.0,	1.0/3.0,	0}, // fa2=fa3=1/3
	{	1.0,	0.541,	0,	0,		0,	0,	0.1,	0,	0}, // fa2=0.1
	{	1.0,	0,	0,	0.828,		0,	0,	0,	0.1,	0}, // fa3=0.1
	{	1.0,	0.574,	0,	0.878,		0,	0,	0.1,	0.1,	0}, // fa2=fa3=0.1

	{	1.0,	0,	0,	0,		0,	0,	0,	0,	-12003.14}, // fLambda1=0.5
	{	1.0,	0,	0,	0,		0,	0,	0,	0,	-7857.90}, // g1pp=-7857.90, g1=1: flambda1=0.3
	{	1.0,	0,	0,	0,		0,	0,	0,	0,	-4001.05}, // g1pp=-4001.05, g1=1: flambda1=0.1 

	{	1.0,	1.624,	0,	0,		PI_VAL/2.0,	0,	0.5,	0,	0}, // fa2=0.5, phia2=90
	{	1.0,	0,	0,	2.484,		0,	PI_VAL/2.0,	0,	0.5,	0}, // fa3=0.5, phia3=90
	{	0,	0.654,	0,	1.0,		0,	PI_VAL/2.0,	0.5,	0.5,	0}, // fa2=fa3=0.5, phia3=90

	{	1.0,	1.624,	0,	2.484,		0,	PI_VAL/2.0,	1.0/3.0,	1.0/3.0,	0}, // fa2=fa3=1/3, phia3=90
	{	1.0,	0.541,	0,	0,		PI_VAL/2.0,	0,	0.1,	0,	0}, // fa2=0.1, phia2=90
	{	1.0,	0,	0,	0.828,		0,	PI_VAL/2.0,	0,	0.1,	0}, // fa3=0.1, phia3=90
	{	1.0,	0.574,	0,	0.878,		0,	PI_VAL/2.0,	0.1,	0.1,	0}, // fa2=0.1, fa3=0.1, phia3=90

	{	1.0,	1.624,	0,	0,		PI_VAL,	0,	0.5,	0,	0}, // fa2=-0.5
	{	1.0,	0,	0,	2.484,		0,	PI_VAL,	0,	0.5,	0}, // fa3=-0.5
	{	0,	0.654,	0,	1.0,		PI_VAL,	0,	0.5,	0.5,	0} // fa2=-0.5, fa3=0.5
};

string user_dir="/afs/cern.ch/work/u/usarica/HZZ4l-126_0-FullSimAnalysis/";
char* user_folder[5]={
	"4mu",
	"4e",
	"2mu2e",
	"CR",
	"data"
};

const int kNumBkg=8;
char hzz4lprefix[]="HZZ4lTree_";
char* sample_BackgroundFile[kNumBkg]={
	"ggZZ4l",
	"ggZZ2l2l",
	"ZZTo4mu",
	"ZZTo4e",
	"ZZTo2e2mu",
	"ZZTo2mu2tau",
	"ZZTo2e2tau",
	"ZZTo4tau"
};

#endif
