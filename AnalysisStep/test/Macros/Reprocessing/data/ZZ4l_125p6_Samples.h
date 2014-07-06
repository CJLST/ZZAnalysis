#ifndef ZZ4L_125.6_SAMPLES_H
#define ZZ4L_125.6_SAMPLES_H

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
const int kNumFiles=24;
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

	kfg2_05_fLambda1_m05,
	kfg2_05_fLambda1_m05_p2270,
	kfg2_05_fLambda1_05,
	kfg2_33_fLambda1_33_p2Pi,

	kfg4_05_fLambda1_m05,
	kfg4_05_fLambda1_m05_p3270,
	kfg4_05_fLambda1_05,
	kfg4_33_fLambda1_33_p3Pi,

	kfZG_1_fGG_0,
	kfZG_0_fGG_1,
	kfZG_05_fGG_0,
	kfZG_0_fGG_05,

	kfMZG_1_fMGG_0,
	kfMZG_0_fMGG_1,
	kfMZG_05_fMGG_0,
	kfMZG_0_fMGG_05,

	kfZG_05_fMZG_05,
	kfGG_05_fMGG_05,

	kfZG_SM_fGG_SM,

	kfZG_Lambda1_1,
	kfZG_Lambda1_05,

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

	"fg2_05_fLambda1_m05",
	"fg2_05_fLambda1_m05_p2270",
	"fg2_05_fLambda1_05",
	"fg2_33_fLambda1_33_p2Pi",

	"fg4_05_fLambda1_m05",
	"fg4_05_fLambda1_m05_p3270",
	"fg4_05_fLambda1_05",
	"fg4_33_fLambda1_33_p3Pi",

	"fZG_1_fGG_0",
	"fZG_0_fGG_1",
	"fZG_05_fGG_0",
	"fZG_0_fGG_05",

	"fMZG_1_fMGG_0",
	"fMZG_0_fMGG_1",
	"fMZG_05_fMGG_0",
	"fMZG_0_fMGG_05",

	"fZG_05_fMZG_05",
	"fGG_05_fMGG_05",

	"fZG_SM_fGG_SM",

	"fZG_Lambda1_1",
	"fZG_Lambda1_05"
};
char* sampleName_label[kNumSamples]={
	"0^{+}_{m} 125.6 GeV",
	"0^{+}_{h} 125.6 GeV",
	"0^{-} 125.6 GeV",
	"0^{+}_{m} Pure Q^{2} 125.6 GeV",

	"f_{a2}=0.5 #phi_{a2}=0 125.6 GeV",
	"f_{a3}=0.5 #phi_{a3}=0 125.6 GeV",
	"f_{a2}=f_{a3}=0.5 #phi_{a2}=#phi_{a3}=0 125.6 GeV",
	"f_{#Lambda1}=-0.5 125.6 GeV",

	"f_{a2}=f_{a3}=0.33 #phi_{a2}=#phi_{a3}=0 125.6 GeV",
	"f_{a2}=0.1 #phi_{a2}=0 125.6 GeV",
	"f_{a3}=0.1 #phi_{a3}=0 125.6 GeV",
	"f_{a2}=f_{a3}=0.1 #phi_{a2}=#phi_{a3}=0 125.6 GeV",

	"f_{#Lambda1}=0.5 125.6 GeV",
	"f_{#Lambda1}=0.3 125.6 GeV",
	"f_{#Lambda1}=0.1 125.6 GeV",

	"f_{a2}=0.5 #phi_{a2}=#pi/2 125.6 GeV",
	"f_{a3}=0.5 #phi_{a3}=#pi/2 125.6 GeV",
	"f_{a2}=f_{a3}=0.5 #phi_{a3}=#pi/2 125.6 GeV",

	"f_{a2}=f_{a3}=0.33 #phi_{a3}=#pi/2 125.6 GeV",
	"f_{a2}=0.1 #phi_{a2}=#pi/2 125.6 GeV",
	"f_{a3}=0.1 #phi_{a3}=#pi/2 125.6 GeV",
	"f_{a2}=f_{a3}=0.1 #phi_{a3}=#pi/2 125.6 GeV",

	"f_{a2}=0.5 #phi_{a2}=#pi 125.6 GeV",
	"f_{a3}=0.5 #phi_{a3}=#pi 125.6 GeV",
	"f_{a2}=f_{a3}=0.5 #phi_{a2}=#pi 125.6 GeV",

	"f_{a2}=f_{#Lambda1}=0.5 #phi_{#Lambda1}=#pi 125.6 GeV",
	"f_{a2}=f_{#Lambda1}=0.5 #phi_{a2}=3#pi/2 #phi_{#Lambda1}=#pi 125.6 GeV",
	"f_{a2}=f_{#Lambda1}=0.5 #phi_{#Lambda1}=0 125.6 GeV",
	"f_{a2}=f_{#Lambda1}=0.33 #phi_{a2}=#pi 125.6 GeV",

	"f_{a3}=f_{#Lambda1}=0.5 #phi_{#Lambda1}=#pi 125.6 GeV",
	"f_{a3}=f_{#Lambda1}=0.5 #phi_{a3}=3#pi/2 #phi_{#Lambda1}=#pi 125.6 GeV",
	"f_{a3}=f_{#Lambda1}=0.5 #phi_{#Lambda1}=0 125.6 GeV",
	"f_{a3}=f_{#Lambda1}=0.33 #phi_{a3}=#pi 125.6 GeV",

	"f_{Z#gamma}=1, f_{#gamma#gamma}=0 125.6 GeV",
	"f_{Z#gamma}=0, f_{#gamma#gamma}=1 125.6 GeV",
	"f_{Z#gamma}=0.5, f_{#gamma#gamma}=0 125.6 GeV",
	"f_{Z#gamma}=0, f_{#gamma#gamma}=0.5 125.6 GeV",

	"f_{0M-Z#gamma}=1, f_{0M-#gamma#gamma}=0 125.6 GeV",
	"f_{0M-Z#gamma}=0, f_{0M-#gamma#gamma}=1 125.6 GeV",
	"f_{0M-Z#gamma}=0.5, f_{0M-#gamma#gamma}=0 125.6 GeV",
	"f_{0M-Z#gamma}=0, f_{0M-#gamma#gamma}=0.5 125.6 GeV",

	"f_{Z#gamma}=0.5, f_{0M-Z#gamma}=0.5 125.6 GeV",
	"f_{#gamma#gamma}=0.5, f_{0M-#gamma#gamma}=0.5 125.6 GeV",

	"f_{Z#gamma}, f_{#gamma#gamma} Full SM 125.6 GeV",

	"f_{Z#gamma, #Lambda1}=1 125.6 GeV",
	"f_{Z#gamma, #Lambda1}=0.5 125.6 GeV"
};
char* sample_file[kNumFiles]={
	"HZZ4lTree_powheg15jhuGenV3-0PMH125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHH125.6",
	"HZZ4lTree_powheg15jhuGenV3-0MH125.6",
	"HZZ4lTree_powheg15jhuGenV3-0L1H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf05ph0H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0Mf05ph0H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf05ph0Mf05ph0H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0L1f05ph180H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf033ph0Mf033ph0H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf01ph0H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0Mf01ph0H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf01ph0Mf01ph0H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0L1f05ph0H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0L1f01ph0H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf05ph90H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0Mf05ph90H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf05ph0Mf05ph90H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf033ph0Mf033ph90H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf01ph90H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0Mf01ph90H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf01ph0Mf01ph90H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf05ph180H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0Mf05ph180H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf05ph180Mf05ph0H125.6"
};
char* sample_file_label[kNumFiles]={
	"0^{+}_{m} 125.6 GeV",
	"0^{+}_{h} 125.6 GeV",
	"0^{-} 125.6 GeV",
	"0^{+}_{m} Pure Q^{2} 125.6 GeV",

	"f_{a2}=0.5 #phi_{a2}=0 125.6 GeV",
	"f_{a3}=0.5 #phi_{a3}=0 125.6 GeV",
	"f_{a2}=f_{a3}=0.5 #phi_{a2}=#phi_{a3}=0 125.6 GeV",
	"f_{#Lambda1}=-0.5 125.6 GeV",

	"f_{a2}=f_{a3}=0.33 #phi_{a2}=#phi_{a3}=0 125.6 GeV",
	"f_{a2}=0.1 #phi_{a2}=0 125.6 GeV",
	"f_{a3}=0.1 #phi_{a3}=0 125.6 GeV",
	"f_{a2}=f_{a3}=0.1 #phi_{a2}=#phi_{a3}=0 125.6 GeV",

	"f_{#Lambda1}=0.5 125.6 GeV",
	"f_{#Lambda1}=0.1 125.6 GeV",

	"f_{a2}=0.5 #phi_{a2}=#pi/2 125.6 GeV",
	"f_{a3}=0.5 #phi_{a3}=#pi/2 125.6 GeV",
	"f_{a2}=f_{a3}=0.5 #phi_{a3}=#pi/2 125.6 GeV",

	"f_{a2}=f_{a3}=0.33 #phi_{a3}=#pi/2 125.6 GeV",
	"f_{a2}=0.1 #phi_{a2}=#pi/2 125.6 GeV",
	"f_{a3}=0.1 #phi_{a3}=#pi/2 125.6 GeV",
	"f_{a2}=f_{a3}=0.1 #phi_{a3}=#pi/2 125.6 GeV",

	"f_{a2}=0.5 #phi_{a2}=#pi 125.6 GeV",
	"f_{a3}=0.5 #phi_{a3}=#pi 125.6 GeV",
	"f_{a2}=f_{a3}=0.5 #phi_{a2}=#pi 125.6 GeV"
};
char* sample_filePrimary[kNumFiles]={
	"ZZ4lAnalysis_powheg15jhuGenV3-0PMH125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0PHH125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0MH125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0L1H125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0PHf05ph0H125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0Mf05ph0H125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0PHf05ph0Mf05ph0H125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0L1f05ph180H125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0PHf033ph0Mf033ph0H125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0PHf01ph0H125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0Mf01ph0H125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0PHf01ph0Mf01ph0H125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0L1f05ph0H125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0L1f01ph0H125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0PHf05ph90H125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0Mf05ph90H125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0PHf05ph0Mf05ph90H125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0PHf033ph0Mf033ph90H125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0PHf01ph90H125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0Mf01ph90H125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0PHf01ph0Mf01ph90H125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0PHf05ph180H125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0Mf05ph180H125.6",
	"ZZ4lAnalysis_powheg15jhuGenV3-0PHf05ph180Mf05ph0H125.6"
};
char* sample_fileLHE[kNumFiles]={
	"Higgs0PMToZZTo4L_M-125p6_",
	"Higgs0PHToZZTo4L_M-125p6_",
	"Higgs0MToZZTo4L_M-125p6_",
	"Higgs0L1ToZZTo4L_M-125p6_",
	"Higgs0PHf05ph0ToZZTo4L_M-125p6_",
	"Higgs0Mf05ph0ToZZTo4L_M-125p6_",
	"Higgs0PHf05ph0Mf05ph0ToZZTo4L_M-125p6_",
	"Higgs0L1f05ph180ToZZTo4L_M-125p6_",
	"Higgs0PHf033ph0Mf033ph0ToZZTo4L_M-125p6_",
	"Higgs0PHf01ph0ToZZTo4L_M-125p6_",
	"Higgs0Mf01ph0ToZZTo4L_M-125p6_",
	"Higgs0PHf01ph0Mf01ph0ToZZTo4L_M-125p6_",
	"Higgs0L1f05ph0ToZZTo4L_M-125p6_",
	"Higgs0L1f01ph0ToZZTo4L_M-125p6_",
	"Higgs0PHf05ph90ToZZTo4L_M-125p6_",
	"Higgs0Mf05ph90ToZZTo4L_M-125p6_",
	"Higgs0PHf05ph0Mf05ph90ToZZTo4L_M-125p6_",
	"Higgs0PHf033ph0Mf033ph90ToZZTo4L_M-125p6_",
	"Higgs0PHf01ph90ToZZTo4L_M-125p6_",
	"Higgs0Mf01ph90ToZZTo4L_M-125p6_",
	"Higgs0PHf01ph0Mf01ph90ToZZTo4L_M-125p6_",
	"Higgs0PHf05ph180ToZZTo4L_M-125p6_",
	"Higgs0Mf05ph180ToZZTo4L_M-125p6_",
	"Higgs0PHf05ph180Mf05ph0ToZZTo4L_M-125p6_"
};
const double gi_phi2_phi4_files[kNumFiles][14]={ // g1-4; phia2,3; MASS, WIDTH; g1pp; gZG, gGG; gZG-PS, gGG-PS
	{	1.0,	0,	0,	0,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // Pure SM
	{	0,	1.0,	0,	0,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=1
	{	0,	0,	0,	1.0,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa3=1
	{	0,	0,	0,	0,		0,	0,	125.6,	0.00415,	1.0, 0, 0, 0, 0, 0}, // fL1=1
											
	{	1.0,	1.638,	0,	0,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=0.5
	{	1.0,	0,	0,	2.521,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa3=0.5
	{	0,	0.650,	0,	1.0,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=fa3=0.5
	{	1.0,	0,	0,	0,		0,	0,	125.6,	0.00415,	12046.01, 0, 0, 0, 0, 0}, // fLambda1=-0.5, for T3 templates
											
	{	1.0,	1.638,	0,	2.521,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=fa3=1/3
	{	1.0,	0.546,	0,	0,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=0.1
	{	1.0,	0,	0,	0.840,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa3=0.1
	{	1.0,	0.579,	0,	0.891,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=fa3=0.1
											
	{	1.0,	0,	0,	0,		0,	0,	125.6,	0.00415,	-12046.01, 0, 0, 0, 0, 0}, // fLambda1=0.5
	// NOTE: No FLambda1=0.3 sample at 125.6 GeV is available										
	{	1.0,	0,	0,	0,		0,	0,	125.6,	0.00415,	-4015.337, 0, 0, 0, 0, 0}, // flambda1=0.1 
											
	{	1.0,	1.638,	0,	0,		PI_VAL/2.0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=0.5, phia2=90
	{	1.0,	0,	0,	2.521,		0,	PI_VAL/2.0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa3=0.5, phia3=90
	{	0,	0.650,	0,	1.0,		0,	PI_VAL/2.0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=fa3=0.5, phia3=90
											
	{	1.0,	1.638,	0,	2.521,		0,	PI_VAL/2.0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=fa3=1/3, phia3=90
	{	1.0,	0.546,	0,	0,		PI_VAL/2.0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=0.1, phia2=90
	{	1.0,	0,	0,	0.840,		0,	PI_VAL/2.0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa3=0.1, phia3=90
	{	1.0,	0.579,	0,	0.891,		0,	PI_VAL/2.0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=0.1, fa3=0.1, phia3=90
											
	{	1.0,	1.638,	0,	0,		PI_VAL,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=-0.5
	{	1.0,	0,	0,	2.521,		0,	PI_VAL,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa3=-0.5
	{	0,	0.650,	0,	1.0,		PI_VAL,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0} // fa2=-0.5, fa3=0.5
};

double sample_mixedxsec_ratio[4]={ // Assuming g_i below
	3.660,	2.0,	0.3146, // 0+m0+h, 0+m0-, 0+h0-
	3.966863 // fL1=-0.5
};
double epspr_mixedSample_lambda1pr=10000.0;

const double gi_phi2_phi4[kNumSamples][14]={ // g1-4; phia2,3; MASS, WIDTH; g1pp; gZG, gGG; gZG-PS, gGG-PS
	{	1.0,	0,	0,	0,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // Pure SM 0
	{	0,	1.0,	0,	0,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=1 1
	{	0,	0,	0,	1.0,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa3=1 2
	{	0,	0,	0,	0,		0,	0,	125.6,	0.00415,	1.0, 0, 0, 0, 0, 0}, // fL1=1 3
											
	{	1.0,	1.638,	0,	0,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=0.5 4
	{	1.0,	0,	0,	2.521,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa3=0.5 5
	{	0,	0.650,	0,	1.0,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=fa3=0.5 6
	{	1.0,	0,	0,	0,		0,	0,	125.6,	0.00415,	12046.01, 0, 0, 0, 0, 0}, // fLambda1=-0.5, for T3 templates 7
											
	{	1.0,	1.638,	0,	2.521,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=fa3=1/3 8
	{	1.0,	0.546,	0,	0,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=0.1 9
	{	1.0,	0,	0,	0.840,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa3=0.1 10
	{	1.0,	0.579,	0,	0.891,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=fa3=0.1 11
											
	{	1.0,	0,	0,	0,		0,	0,	125.6,	0.00415,	-12046.01, 0, 0, 0, 0, 0}, // fLambda1=0.5 12
	{	1.0,	0,	0,	0,		0,	0,	125.6,	0.00415,	-7885.965, 0, 0, 0, 0, 0}, // flambda1=0.3 13
	{	1.0,	0,	0,	0,		0,	0,	125.6,	0.00415,	-4015.337, 0, 0, 0, 0, 0}, // flambda1=0.1 14
											
	{	1.0,	1.638,	0,	0,		PI_VAL/2.0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=0.5, phia2=90 15
	{	1.0,	0,	0,	2.521,		0,	PI_VAL/2.0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa3=0.5, phia3=90 16
	{	0,	0.650,	0,	1.0,		0,	PI_VAL/2.0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=fa3=0.5, phia3=90 17
											
	{	1.0,	1.638,	0,	2.521,		0,	PI_VAL/2.0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=fa3=1/3, phia3=90 18
	{	1.0,	0.546,	0,	0,		PI_VAL/2.0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=0.1, phia2=90 19
	{	1.0,	0,	0,	0.840,		0,	PI_VAL/2.0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa3=0.1, phia3=90 20
	{	1.0,	0.579,	0,	0.891,		0,	PI_VAL/2.0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=0.1, fa3=0.1, phia3=90 21
											
	{	1.0,	1.638,	0,	0,		PI_VAL,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=-0.5 22
	{	1.0,	0,	0,	2.521,		0,	PI_VAL,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa3=-0.5 23
	{	0,	0.650,	0,	1.0,		PI_VAL,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 0}, // fa2=-0.5, fa3=0.5 24
											
	{	0,	1.638,	0,	0,		0,	0,	125.6,	0.00415,	12046.01, 0, 0, 0, 0, 0}, // fa2=0.5, fLambda1=-0.5, for templates 25
	{	0,	1.638,	0,	0,		PI_VAL*3.0/2.0,	0,	125.6,	0.00415,	12046.01, 0, 0, 0, 0, 0}, // fa2=-0.5i, fLambda1=-0.5, for templates 26
	{	0,	1.638,	0,	0,		0,	0,	125.6,	0.00415,	-12046.01, 0, 0, 0, 0, 0}, // fa2=0.5, fLambda1=0.5 27
	{	1.0,	1.638,	0,	0,		PI_VAL,	0,	125.6,	0.00415,	-12046.01, 0, 0, 0, 0, 0}, // fa2=-0.33, fLambda1=0.33 28
											
	{	0,	0,	0,	2.521,		0,	0,	125.6,	0.00415,	12046.01, 0, 0, 0, 0, 0}, // fa3=0.5, fLambda1=-0.5, for templates 29
	{	0,	0,	0,	2.521,		PI_VAL*3.0/2.0,	0,	125.6,	0.00415,	12046.01, 0, 0, 0, 0, 0}, // fa3=-0.5i, fLambda1=-0.5, for templates 30
	{	0,	0,	0,	2.521,		0,	0,	125.6,	0.00415,	-12046.01, 0, 0, 0, 0, 0}, // fa3=0.5, fLambda1=0.5 31
	{	1.0,	0,	0,	2.521,		PI_VAL,	0,	125.6,	0.00415,	-12046.01, 0, 0, 0, 0, 0}, // fa3=-0.33, fLambda1=0.33 32
											
	{	0,	0,	0,	0,		0,	0,	125.6,	0.00415,	0, 1, 0, 0, 0, 0}, // Pure PM-ZG 33
	{	0,	0,	0,	0,		0,	0,	125.6,	0.00415,	0, 0, 1, 0, 0, 0}, // Pure PM-GG 34
	{	1.0,	0,	0,	0,		0,	0,	125.6,	0.00415,	0, 0.0473, 0, 0, 0, 0}, // fPM-ZG=0.5 35
	{	1.0,	0,	0,	0,		0,	0,	125.6,	0.00415,	0, 0, -0.0531, 0, 0, 0}, // fPM-GG=0.5 36
										
	{	0,	0,	0,	0,		0,	0,	125.6,	0.00415,	0, 0, 0, 1, 0, 0}, // Pure M-ZG 37
	{	0,	0,	0,	0,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 1, 0}, // Pure M-GG 38
	{	1.0,	0,	0,	0,		0,	0,	125.6,	0.00415,	0, 0, 0, 0.052161, 0, 0}, // fM-ZG=0.5 39
	{	1.0,	0,	0,	0,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, -0.053666, 0}, // fM-GG=0.5 40
										
	{	0,	0,	0,	0,		0,	0,	125.6,	0.00415,	0, 0.0473, 0, 0.052161, 0, 0}, // fZG=fM-ZG=0.5 41
	{	0,	0,	0,	0,		0,	0,	125.6,	0.00415,	0, 0, -0.0531, 0, -0.053666, 0}, // fGG=fM-GG=0.5 42
										
	{	1.0,	0,	0,	0,		0,	0,	125.6,	0.00415,	0, 0.00175, -0.002, 0, 0, 0}, // SM ZZ, ZG. GG 43

	{	0,	0,	0,	0,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, 1.0}, // Pure ZG-L1 44
	{	1.0,	0,	0,	0,		0,	0,	125.6,	0.00415,	0, 0, 0, 0, 0, -7591.914} // fZG-L1=0.5 45

};

string user_dir_raw = "/afs/cern.ch/work/u/usarica/ggtoH-PWGSamples-125_6/";
string user_dir="/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/";
string user_dir_gritsan="/afs/cern.ch/work/g/gritsan/public/backup/newProduction/";
string user_dir_newProduction="/afs/cern.ch/work/l/lfeng/public_write/usarica/Production/";
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
