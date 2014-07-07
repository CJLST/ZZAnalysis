#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <ctime>
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TChain.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "./data/ZZ4l_125p6_Samples.h"

using namespace std;

const int kSignalSamples=24;
//const int kggZZSamples=6;
const int kggZZSamples=7;
const int kqqZZSamples=6;
const int kqqZZSamples_Dedicated=5;
const int kZXSamples=1;
char* sampleName_Sig[kSignalSamples]={
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
char* sampleName_ggZZ[kggZZSamples]={
	"HZZ4lTree_ggZZ4l",
	"HZZ4lTree_ggZZ2l2l",
	"HZZ4lTree_ggTo4mu_Contin-MCFM67",
	"HZZ4lTree_ggTo4e_Contin-MCFM67",
	"HZZ4lTree_ggTo2e2mu_Contin-MCFM67",
	"HZZ4lTree_ggTo2l2l_Continuum",
	"HZZ4lTree_ggTo4l_Continuum"
};
char* sampleName_qqZZ[kqqZZSamples]={
	"HZZ4lTree_ZZTo4mu",
	"HZZ4lTree_ZZTo4e",
	"HZZ4lTree_ZZTo2e2mu",
	"HZZ4lTree_ZZTo2mu2tau",
	"HZZ4lTree_ZZTo2e2tau",
	"HZZ4lTree_ZZTo4tau"
};
char* sampleName_qqZZDedicated[kqqZZSamples_Dedicated]={
	"HZZ4lTree_ZZ95-160To2e2mu",
	"HZZ4lTree_ZZ95-160To2mu2tau",
	"HZZ4lTree_ZZ95-160To4e",
	"HZZ4lTree_ZZ95-160To4mu",
	"HZZ4lTree_ZZ95-160To4tau"
};
char* sampleName_ZX[kZXSamples]={
	"HZZ4lTree_DoubleOr_CRZLLTree"
};

float evaluateAlphaSShift(float mzz, int analytic=1, int systematic=0, float mPOLE = 125.6){
	float scale_alphaS = 1;
	if(analytic==1){
		float parA = 1.28618e-2;
		float parErrA = 1.41313e-3;
		float parB = 7.27744e-1;
		float parErrB = 3.07639e-2;
		float myA = parA;
		float myB = parB;
		float myDiv = 2.0;

		scale_alphaS = pow( ( (1.0 + myA*mPOLE/2.0) / (1.0 + myA*mzz/myDiv) ),myB);

		if(systematic==1){
			parA = 4.56257e-2;
			parErrA = 2.20198e-3;
			parB = 4.74865e-1;
			parErrB = 5.83891e-3;
			myA = parA;
			myB = parB;

			myDiv=4.0;
			myA -= parErrA;
			myB -= parErrB;

			scale_alphaS = pow( ( (1.0 + myA*mzz/2.0) / (1.0 + myA*mzz/myDiv) ),myB);
		}
		else if(systematic==-1){
			parA = 3.04304e-2;
			parErrA = 1.89309e-3;
			parB = 4.18857e-1;
			parErrB = 5.32661e-3;
			myA = parA;
			myB = parB;

			myDiv=1.0;
			myA += parErrA;
			myB += parErrB;

			scale_alphaS = pow( ( (1.0 + myA*mzz/2.0) / (1.0 + myA*mzz/myDiv) ),myB);
		};
	};
	return scale_alphaS;
};

void shuffleSignalBkg_125p6(int flavor, int erg_tev){
	srand( time(0) );

	string cinput_KDFactor="./data/HZZ4l-KDFactorGraph";
	if(erg_tev==7) cinput_KDFactor = cinput_KDFactor + "_7TeV";
	cinput_KDFactor = cinput_KDFactor + ".root";
	TFile* finput_KDFactor = new TFile(cinput_KDFactor.c_str(),"read");
	string tgkfname = "KDScale_";
	tgkfname = tgkfname + "AllFlavors_UnNormalized";
	TGraphAsymmErrors* tgkf = (TGraphAsymmErrors*) finput_KDFactor->Get(tgkfname.c_str());
	string tgkfname_up = tgkfname + "_Down"; // K+ is down
	TGraphAsymmErrors* tgkf_up = (TGraphAsymmErrors*) finput_KDFactor->Get(tgkfname_up.c_str());
	string tgkfname_down = tgkfname + "_Up"; // K- is up
	TGraphAsymmErrors* tgkf_down = (TGraphAsymmErrors*) finput_KDFactor->Get(tgkfname_down.c_str());

	string OUTPUT_NAME = hzz4lprefix;
	OUTPUT_NAME = OUTPUT_NAME + "H125p6_ShuffledSignalBkg";
	TString TREE_NAME = "SelectedTree";

	char erg_dir[1000];
	sprintf(erg_dir,"LHC_%iTeV",erg_tev);
	char cerg[1000];
	sprintf(cerg,"%iTeV",erg_tev);
	string cinput_common = user_dir + erg_dir;
//	string cinput_common = user_dir_newProduction + cerg;
	string coutput_common = cinput_common;

	const int numSamples=46;
	int myNumSamples=numSamples;
	float MC_weight_spin0[numSamples]={0};
	float MC_weight_4GeVcut_spin0[numSamples]={0};

	float GenHMass=0;
	float GenZ1Mass=0;
	float GenZ2Mass=0;
	float ZZMass=0;
	float Z1Mass=0;
	float Z2Mass=0;
	float mycosthetastar,myhelphi,myhelcosthetaZ1,myhelcosthetaZ2,myphistarZ1,myphistarZ2;
	short Z1ids,Z2ids;	

	float ZXfake_weightProper=0;
	Long64_t EventNumber=-1;
	int RunNumber=-1;
	int EventSample=-1;

	int genFinalState=0;
	float MC_weight=1;
	float MC_weight_noxsec=1;
	float MC_weight_xsec=0;
	float MC_weight_4GeVcut=1;
	float MC_weight_xsec_4GeVcut=0;
	float MC_weight_norm=1;
	float MC_weight_PUWeight=1;
	float MC_weight_powhegWeight=1;
	float MC_weight_dataMC=1;
	float MC_weight_HqT=1;

	int kNumQQBGGWeights=2;
	float MC_weight_QQBGGProper[2]={0};

	float MC_weight_down=1;
	float MC_weight_up=1;
	float MC_weight_Kfactor=1;
	float MC_weight_Kfactor_down=1;
	float MC_weight_Kfactor_up=1;
	float MC_weight_PDF_down=0;
	float MC_weight_PDF_up=0;
	float MC_weight_Kfactor_Norm_down=1;
	float MC_weight_Kfactor_Norm_up=1;
	float MC_weight_PDF_Norm_down=0;
	float MC_weight_PDF_Norm_up=0;
	float MC_weight_alphaS_down=1;
	float MC_weight_alphaS_up=1;
	float MC_weight_alphaS=1;

	float D_ggZZ=-99;
	float D_bkg_kin=-99;
	float D_bkg=-99;
	float D_bkg_ScaleUp=-99;
	float D_bkg_ScaleDown=-99;
	float D_bkg_ResUp=-99;
	float D_bkg_ResDown=-99;
	float D_bkg_prodIndep=-99;
	float D_bkg_prodIndep_ScaleUp=-99;
	float D_bkg_prodIndep_ScaleDown=-99;
	float D_bkg_prodIndep_ResUp=-99;
	float D_bkg_prodIndep_ResDown=-99;
	float D_Gamma_r1=-99;
	float D_Gamma_r5=-99;
	float D_Gamma_r10=-99;
	float D_Gamma_r15=-99;
	float D_Gamma_r20=-99;
	float D_Gamma_r25=-99;
	float D_Gamma_gg_r1=-99;
	float D_Gamma_gg_r5=-99;
	float D_Gamma_gg_r10=-99;
	float D_Gamma_gg_r15=-99;
	float D_Gamma_gg_r20=-99;
	float D_Gamma_gg_r25=-99;
	float D_Gamma=-99;
	float D_Gamma_int=-99;
	float D_g1q2=-99,D_g1q2int=-99;
	float D_g2=-99,D_g2int=-99,D_g2int_perp=-99;
	float D_g4=-99,D_g4int=-99,D_g4int_perp=-99;
	float D_ZG=-99, D_GG=-99;
	float D_ZG_PS=-99, D_GG_PS=-99;
	float D_ZGint=-99, D_GGint=-99;
	float D_ZG_PSint=-99, D_GG_PSint=-99;
	float D_ZG_L1=-99,D_ZG_L1int=-99,D_ZG_L1int_perp=-99;

	// Spin 1 and 2 discriminants
	float p1plusProdIndepKD = -99;
	float p1minusProdIndepKD = -99;
	float p1plusKD =-99;
	float p1minusKD =-99;
	float graviKD = -99;
	float qqgraviKD = -99;
	float p2h2plusKD = -99;
	float p2h2plus_qqb_KD = -99;
	float p2h3plusKD = -99;
	float p2h3plus_qqb_KD = -99;
	float p2hplusKD = -99;
	float p2hplus_qqb_KD = -99;
	float p2bplusKD = -99;
	float p2bplus_qqb_KD = -99;
	float p2h6plusKD = -99;
	float p2h6plus_qqb_KD = -99;
	float p2h7plusKD = -99;
	float p2h7plus_qqb_KD = -99;
	float p2hminusKD = -99;
	float p2hminus_qqb_KD = -99;
	float p2h9minusKD = -99;
	float p2h9minus_qqb_KD = -99;
	float p2h10minusKD = -99;
	float p2h10minus_qqb_KD = -99;
	float p2mProdIndepKD = -99;
	float p2h2plusProdIndepKD = -99;
	float p2h3plusProdIndepKD = -99;
	float p2hplusProdIndepKD = -99;
	float p2bplusProdIndepKD = -99;
	float p2h6plusProdIndepKD = -99;
	float p2h7plusProdIndepKD = -99;
	float p2hminusProdIndepKD = -99;
	float p2h9minusProdIndepKD = -99;
	float p2h10minusProdIndepKD = -99;

	float MC_CV_weight[numSamples];
	float MC_CV_4GeVcut_weight[numSamples];

	TFile* finput[kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples+kqqZZSamples_Dedicated];
	TTree* myTree[kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples+kqqZZSamples_Dedicated];

	int nevents[kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples+kqqZZSamples_Dedicated]={0};
	int naccevents[kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples+kqqZZSamples_Dedicated]={0};

	int nTotal_Sig=0;
	int nTotal_ggZZ=0;
	int nTotal_qqZZ=0;
	int nTotal_qqZZ_Dedicated=0;
	int nTotal_ZX=0;

	double NGenTotal_weighted[kSignalSamples][numSamples];
	double NGenTotal_unweighted[kSignalSamples]={0};
	double NGenTotal_weighted_perFlavor[kSignalSamples][numSamples][nFinalStates];
	double fGen_BRweight_perFlavor[numSamples][nFinalStates]={{0}};
	double fGen_BRweight[numSamples]={0};
	double NGenTotal_unweighted_perFlavor[kSignalSamples][nFinalStates];
	double sum_NGenTotal_unweighted=0;
	double sum_NGenTotal_weighted[numSamples]={0};

	double NGenTotal_4GeVcut_weighted[kSignalSamples][numSamples];
	double NGenTotal_4GeVcut_unweighted[kSignalSamples]={0};
	double NGenTotal_4GeVcut_weighted_perFlavor[kSignalSamples][numSamples][nFinalStates];
	double fGen_4GeVcut_BRweight_perFlavor[numSamples][nFinalStates]={{0}};
	double fGen_4GeVcut_BRweight[numSamples]={0};
	double NGenTotal_4GeVcut_unweighted_perFlavor[kSignalSamples][nFinalStates];
	double sum_NGenTotal_4GeVcut_unweighted=0;
	double sum_NGenTotal_4GeVcut_weighted[numSamples]={0};

	for(int smp=0;smp<kSignalSamples;smp++){
		for(int hypo=0;hypo<numSamples;hypo++){
			NGenTotal_weighted[smp][hypo]=0;
			NGenTotal_4GeVcut_weighted[smp][hypo]=0;
			for (int fl = 0; fl < nFinalStates; fl++){
				NGenTotal_weighted_perFlavor[smp][hypo][fl] = 0;
				NGenTotal_4GeVcut_weighted_perFlavor[smp][hypo][fl] = 0;
			};
		};
		for (int fl = 0; fl < nFinalStates; fl++){
			NGenTotal_unweighted_perFlavor[smp][fl] = 0;
			NGenTotal_4GeVcut_unweighted_perFlavor[smp][fl] = 0;
		};
	};

	for(int smp=0;smp<(kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples+kqqZZSamples_Dedicated);smp++){
		string cinput = cinput_common;
		cinput = cinput + "/" + user_folder[flavor] + "/";
		if(smp<kSignalSamples) cinput = cinput + sampleName_Sig[smp] + "_Reprocessed.root";
		else if(smp>=kSignalSamples && smp<kggZZSamples+kSignalSamples) cinput = cinput + sampleName_ggZZ[smp-kSignalSamples] + "_Reprocessed.root";
		else if(smp>=kggZZSamples+kSignalSamples && smp<kqqZZSamples+kggZZSamples+kSignalSamples) cinput = cinput + sampleName_qqZZ[smp-kSignalSamples-kggZZSamples] + "_Reprocessed.root";
		else if(smp>=kqqZZSamples+kggZZSamples+kSignalSamples && smp<kZXSamples+kqqZZSamples+kggZZSamples+kSignalSamples) cinput = cinput + sampleName_ZX[smp-kSignalSamples-kggZZSamples-kqqZZSamples] + "_Reprocessed.root";
		else if(smp>=kZXSamples+kqqZZSamples+kggZZSamples+kSignalSamples) cinput = cinput + sampleName_qqZZDedicated[smp-kSignalSamples-kggZZSamples-kqqZZSamples-kZXSamples] + "_Reprocessed.root";

		cout << cinput << '\t' << smp << endl;

		finput[smp] = new TFile(cinput.c_str(),"read");
		myTree[smp] = (TTree*) finput[smp]->Get(TREE_NAME);
		nevents[smp] = myTree[smp]->GetEntries();
		naccevents[smp] = nevents[smp];

		if(smp<kSignalSamples){
			string cinputGen = cinput_common + "/GenSignal/" + sampleName_Sig[smp] + "_GenLevel_ReWeighed.root";
			cout << cinputGen << endl;
			TFile* ftemp = new TFile(cinputGen.c_str(),"read");
			TTree* ttemp = (TTree*) ftemp->Get("GenTree");
			float MC_gen_spin0[numSamples]={0};
			int genFlavor=0;
			float genMZ1=0;
			float genMZ2=0;
			int nGenEntries = ttemp->GetEntries();
			ttemp->SetBranchAddress("genFinalState",&genFlavor);
			ttemp->SetBranchAddress("GenZ1Mass",&genMZ1);
			ttemp->SetBranchAddress("GenZ2Mass",&genMZ2);
			ttemp->SetBranchAddress("MC_weight_spin0",MC_gen_spin0);
			TH2F* htrue = (TH2F*) ftemp->Get("hCounters_spin0_RW");
			TH2F* htrue_4GeVcut = (TH2F*) ftemp->Get("hCounters_spin0_4GeVcut_RW");
			cout << htrue->GetName() << '\t' << htrue_4GeVcut->GetName() << endl;
			TH1F* hreco = (TH1F*) finput[smp]->Get("Counters");
			double nTrueProcessed = hreco->GetBinContent(1);
			for(int myf=0;myf<nFinalStates;myf++){
				NGenTotal_unweighted_perFlavor[smp][myf] = htrue->GetBinContent(myf+1,1);
				NGenTotal_4GeVcut_unweighted_perFlavor[smp][myf] = htrue_4GeVcut->GetBinContent(myf+1,1);
			};
			for(int myf=0;myf<nFinalStates;myf++){
				NGenTotal_unweighted[smp] += NGenTotal_unweighted_perFlavor[smp][myf];
				NGenTotal_4GeVcut_unweighted[smp] += NGenTotal_4GeVcut_unweighted_perFlavor[smp][myf];
			};
			for(int myf=0;myf<nFinalStates;myf++){
				for(int p=0;p<numSamples;p++){
					NGenTotal_weighted_perFlavor[smp][p][myf] = htrue->GetBinContent(myf+1,p+2);
					NGenTotal_4GeVcut_weighted_perFlavor[smp][p][myf] = htrue_4GeVcut->GetBinContent(myf+1,p+2);
					if(NGenTotal_weighted_perFlavor[smp][p][myf]!=NGenTotal_weighted_perFlavor[smp][p][myf]){
						NGenTotal_weighted_perFlavor[smp][p][myf]=0;
						NGenTotal_4GeVcut_weighted_perFlavor[smp][p][myf]=0;
						for(int ev=0;ev<nGenEntries;ev++){
							ttemp->GetEntry(ev);
							if (genFlavor == myf && MC_gen_spin0[p] == MC_gen_spin0[p]){
								NGenTotal_weighted_perFlavor[smp][p][myf] += MC_gen_spin0[p];
								if (genMZ1>4 && genMZ2 > 4) NGenTotal_4GeVcut_weighted_perFlavor[smp][p][myf] += MC_gen_spin0[p];
							};
						};
					};
				};
			};
			for(int myf=0;myf<nFinalStates;myf++){
				for(int p=0;p<numSamples;p++){
					NGenTotal_weighted[smp][p] += NGenTotal_weighted_perFlavor[smp][p][myf];
					NGenTotal_4GeVcut_weighted[smp][p] += NGenTotal_4GeVcut_weighted_perFlavor[smp][p][myf];
//					NGenTotal_weighted[smp][p] += NGenTotal_weighted_perFlavor[smp][p][myf]/NGenTotal_weighted_perFlavor[smp][p][2]*NGenTotal_unweighted_perFlavor[smp][2];
				};
			};

			for (int p = 0; p < myNumSamples; p++){
				NGenTotal_weighted[smp][p] *= nTrueProcessed / NGenTotal_unweighted[smp];
				NGenTotal_4GeVcut_weighted[smp][p] *= nTrueProcessed / NGenTotal_unweighted[smp];
			};
			NGenTotal_unweighted[smp] *= nTrueProcessed/NGenTotal_unweighted[smp];
			NGenTotal_4GeVcut_unweighted[smp] *= nTrueProcessed/NGenTotal_unweighted[smp];
			for(int myf=0;myf<nFinalStates;myf++){
				for(int p=0;p<myNumSamples;p++){
					NGenTotal_weighted_perFlavor[smp][p][myf] *= nTrueProcessed/NGenTotal_unweighted[smp];
					NGenTotal_4GeVcut_weighted_perFlavor[smp][p][myf] *= nTrueProcessed/NGenTotal_unweighted[smp];
				};
				NGenTotal_unweighted_perFlavor[smp][myf] *= nTrueProcessed/NGenTotal_unweighted[smp];
				NGenTotal_4GeVcut_unweighted_perFlavor[smp][myf] *= nTrueProcessed/NGenTotal_unweighted[smp];
			};

			sum_NGenTotal_unweighted += NGenTotal_unweighted[smp];
			sum_NGenTotal_4GeVcut_unweighted += NGenTotal_4GeVcut_unweighted[smp];
			for (int p = 0; p < myNumSamples; p++){
				sum_NGenTotal_weighted[p] += NGenTotal_weighted[smp][p];
				sum_NGenTotal_4GeVcut_weighted[p] += NGenTotal_4GeVcut_weighted[smp][p];
			};
			
			delete hreco;
			delete htrue_4GeVcut;
			delete htrue;
			delete ttemp;

			ftemp->Close();
		};

		if(smp<kSignalSamples){
			nTotal_Sig += nevents[smp];
		}
		else if(smp>=kSignalSamples && smp<kggZZSamples+kSignalSamples){
			nTotal_ggZZ += nevents[smp];
		}
		else if(smp>=kggZZSamples+kSignalSamples && smp<kqqZZSamples+kggZZSamples+kSignalSamples){
			nTotal_qqZZ += nevents[smp];
		}
		else if(smp>=kqqZZSamples+kggZZSamples+kSignalSamples && smp<kZXSamples+kqqZZSamples+kggZZSamples+kSignalSamples){
			nTotal_ZX += nevents[smp];
		}
		else if(smp>=kZXSamples+kqqZZSamples+kggZZSamples+kSignalSamples){
			nTotal_qqZZ_Dedicated += nevents[smp];
		};

		if(smp!=0
			&& smp!=kSignalSamples
			&& smp!=kggZZSamples+kSignalSamples
			&& smp!=kqqZZSamples+kggZZSamples+kSignalSamples
			&& smp!=kZXSamples+kqqZZSamples+kggZZSamples+kSignalSamples) naccevents[smp] += naccevents[smp-1];
		cout << nevents[smp] << endl;
		cout << naccevents[smp] << endl;

		if(smp>=kqqZZSamples+kggZZSamples+kSignalSamples && smp<kZXSamples+kqqZZSamples+kggZZSamples+kSignalSamples) myTree[smp]->SetBranchAddress("ZXfake_weightProper",&ZXfake_weightProper);
		else{
			myTree[smp]->SetBranchAddress("MC_weight_spin0",MC_weight_spin0);
			myTree[smp]->SetBranchAddress("genFinalState", &genFinalState);
			myTree[smp]->SetBranchAddress("GenHMass", &GenHMass);
			myTree[smp]->SetBranchAddress("GenZ1Mass", &GenZ1Mass);
			myTree[smp]->SetBranchAddress("GenZ2Mass", &GenZ2Mass);
			myTree[smp]->SetBranchAddress("MC_weight",&MC_weight);
			myTree[smp]->SetBranchAddress("MC_weight_noxsec",&MC_weight_noxsec);
			myTree[smp]->SetBranchAddress("MC_weight_xsec",&MC_weight_xsec);
			myTree[smp]->SetBranchAddress("MC_weight_PUWeight",&MC_weight_PUWeight);
			myTree[smp]->SetBranchAddress("MC_weight_powhegWeight",&MC_weight_powhegWeight);
			myTree[smp]->SetBranchAddress("MC_weight_dataMC",&MC_weight_dataMC);
			myTree[smp]->SetBranchAddress("MC_weight_HqT",&MC_weight_HqT);

			myTree[smp]->SetBranchAddress("MC_weight_Kfactor",&MC_weight_Kfactor);
			myTree[smp]->SetBranchAddress("MC_weight_Kfactor_down",&MC_weight_Kfactor_down);
			myTree[smp]->SetBranchAddress("MC_weight_Kfactor_up",&MC_weight_Kfactor_up);
//			myTree[smp]->SetBranchAddress("MC_weight_down",&MC_weight_down);
//			myTree[smp]->SetBranchAddress("MC_weight_up",&MC_weight_up);

			if( ( smp>=kSignalSamples && smp<kqqZZSamples+kggZZSamples+kSignalSamples ) || ( smp>=kZXSamples+kqqZZSamples+kggZZSamples+kSignalSamples && smp<kqqZZSamples_Dedicated+kZXSamples+kqqZZSamples+kggZZSamples+kSignalSamples ) ){
				myTree[smp]->SetBranchAddress("kNumQQBGGWeights",&kNumQQBGGWeights);
				myTree[smp]->SetBranchAddress("MC_weight_QQBGGProper",MC_weight_QQBGGProper);
			};

			myTree[smp]->SetBranchAddress("GenHMass", &GenHMass);
		};

		myTree[smp]->SetBranchAddress("ZZMass", &ZZMass);
		myTree[smp]->SetBranchAddress("Z1Mass", &Z1Mass);
		myTree[smp]->SetBranchAddress("Z2Mass", &Z2Mass);
		myTree[smp]->SetBranchAddress("helcosthetaZ1", &myhelcosthetaZ1);
		myTree[smp]->SetBranchAddress("helcosthetaZ2", &myhelcosthetaZ2);
		myTree[smp]->SetBranchAddress("helphi", &myhelphi);
		myTree[smp]->SetBranchAddress("costhetastar", &mycosthetastar);
		myTree[smp]->SetBranchAddress("phistarZ1", &myphistarZ1);
		myTree[smp]->SetBranchAddress("Z1ids", &Z1ids);
		myTree[smp]->SetBranchAddress("Z2ids", &Z2ids);

		myTree[smp]->SetBranchAddress("D_g1_vs_g2_phi0",&D_g2);
		myTree[smp]->SetBranchAddress("D_g1_vs_g4_phi0",&D_g4);
		myTree[smp]->SetBranchAddress("D_g2int_phi0",&D_g2int);
		myTree[smp]->SetBranchAddress("D_g4int_phi0",&D_g4int);
		myTree[smp]->SetBranchAddress("D_g2int_phi90",&D_g2int_perp);
		myTree[smp]->SetBranchAddress("D_g4int_phi90",&D_g4int_perp);
		myTree[smp]->SetBranchAddress("D_g1Q2_phi0",&D_g1q2);
		myTree[smp]->SetBranchAddress("D_g1Q2int_phi0",&D_g1q2int);
		myTree[smp]->SetBranchAddress("D_ZG",&D_ZG);
		myTree[smp]->SetBranchAddress("D_GG",&D_GG);
		myTree[smp]->SetBranchAddress("D_ZG_PS",&D_ZG_PS);
		myTree[smp]->SetBranchAddress("D_GG_PS",&D_GG_PS);
		myTree[smp]->SetBranchAddress("D_ZGint",&D_ZGint);
		myTree[smp]->SetBranchAddress("D_GGint",&D_GGint);
		myTree[smp]->SetBranchAddress("D_ZG_PSint",&D_ZG_PSint);
		myTree[smp]->SetBranchAddress("D_GG_PSint",&D_GG_PSint);
		myTree[smp]->SetBranchAddress("D_ZG_L1",&D_ZG_L1);
		myTree[smp]->SetBranchAddress("D_ZG_L1int_phi0",&D_ZG_L1int);
		myTree[smp]->SetBranchAddress("D_ZG_L1int_phi90",&D_ZG_L1int_perp);
		myTree[smp]->SetBranchAddress("D_bkg_kin",&D_bkg_kin);
		myTree[smp]->SetBranchAddress("D_bkg",&D_bkg);
		myTree[smp]->SetBranchAddress("D_bkg_ScaleUp",&D_bkg_ScaleUp);
		myTree[smp]->SetBranchAddress("D_bkg_ScaleDown",&D_bkg_ScaleDown);
		myTree[smp]->SetBranchAddress("D_bkg_ResUp",&D_bkg_ResUp);
		myTree[smp]->SetBranchAddress("D_bkg_ResDown",&D_bkg_ResDown);

		if (smp >= kSignalSamples){
			myTree[smp]->SetBranchAddress("D_bkg_prodIndep", &D_bkg_prodIndep);
			myTree[smp]->SetBranchAddress("D_bkg_prodIndep_ScaleUp", &D_bkg_prodIndep_ScaleUp);
			myTree[smp]->SetBranchAddress("D_bkg_prodIndep_ScaleDown", &D_bkg_prodIndep_ScaleDown);
			myTree[smp]->SetBranchAddress("D_bkg_prodIndep_ResUp", &D_bkg_prodIndep_ResUp);
			myTree[smp]->SetBranchAddress("D_bkg_prodIndep_ResDown", &D_bkg_prodIndep_ResDown);

			myTree[smp]->SetBranchAddress("p1plusProdIndepKD", &p1plusProdIndepKD);
			myTree[smp]->SetBranchAddress("p1minusProdIndepKD", &p1minusProdIndepKD);
			myTree[smp]->SetBranchAddress("p1plusKD", &p1plusKD);
			myTree[smp]->SetBranchAddress("p1minusKD", &p1minusKD);

			myTree[smp]->SetBranchAddress("graviKD", &graviKD);
			myTree[smp]->SetBranchAddress("qqgraviKD", &qqgraviKD);
			myTree[smp]->SetBranchAddress("p2h2plusKD", &p2h2plusKD);
			myTree[smp]->SetBranchAddress("p2h2plus_qqb_KD", &p2h2plus_qqb_KD);
			myTree[smp]->SetBranchAddress("p2h3plusKD", &p2h3plusKD);
			myTree[smp]->SetBranchAddress("p2h3plus_qqb_KD", &p2h3plus_qqb_KD);
			myTree[smp]->SetBranchAddress("p2hplusKD", &p2hplusKD);
			myTree[smp]->SetBranchAddress("p2hplus_qqb_KD", &p2hplus_qqb_KD);
			myTree[smp]->SetBranchAddress("p2bplusKD", &p2bplusKD);
			myTree[smp]->SetBranchAddress("p2bplus_qqb_KD", &p2bplus_qqb_KD);
			myTree[smp]->SetBranchAddress("p2h6plusKD", &p2h6plusKD);
			myTree[smp]->SetBranchAddress("p2h6plus_qqb_KD", &p2h6plus_qqb_KD);
			myTree[smp]->SetBranchAddress("p2h7plusKD", &p2h7plusKD);
			myTree[smp]->SetBranchAddress("p2h7plus_qqb_KD", &p2h7plus_qqb_KD);
			myTree[smp]->SetBranchAddress("p2hminusKD", &p2hminusKD);
			myTree[smp]->SetBranchAddress("p2hminus_qqb_KD", &p2hminus_qqb_KD);
			myTree[smp]->SetBranchAddress("p2h9minusKD", &p2h9minusKD);
			myTree[smp]->SetBranchAddress("p2h9minus_qqb_KD", &p2h9minus_qqb_KD);
			myTree[smp]->SetBranchAddress("p2h10minusKD", &p2h10minusKD);
			myTree[smp]->SetBranchAddress("p2h10minus_qqb_KD", &p2h10minus_qqb_KD);
			myTree[smp]->SetBranchAddress("p2mProdIndepKD", &p2mProdIndepKD);
			myTree[smp]->SetBranchAddress("p2h2plusProdIndepKD", &p2h2plusProdIndepKD);
			myTree[smp]->SetBranchAddress("p2h3plusProdIndepKD", &p2h3plusProdIndepKD);
			myTree[smp]->SetBranchAddress("p2hplusProdIndepKD", &p2hplusProdIndepKD);
			myTree[smp]->SetBranchAddress("p2bplusProdIndepKD", &p2bplusProdIndepKD);
			myTree[smp]->SetBranchAddress("p2h6plusProdIndepKD", &p2h6plusProdIndepKD);
			myTree[smp]->SetBranchAddress("p2h7plusProdIndepKD", &p2h7plusProdIndepKD);
			myTree[smp]->SetBranchAddress("p2hminusProdIndepKD", &p2hminusProdIndepKD);
			myTree[smp]->SetBranchAddress("p2h9minusProdIndepKD", &p2h9minusProdIndepKD);
			myTree[smp]->SetBranchAddress("p2h10minusProdIndepKD", &p2h10minusProdIndepKD);
		};


/*		myTree[smp]->SetBranchAddress("D_ggZZ",&D_ggZZ);
		myTree[smp]->SetBranchAddress("D_Gamma_r1",&D_Gamma_r1);
		myTree[smp]->SetBranchAddress("D_Gamma_r5",&D_Gamma_r5);
		myTree[smp]->SetBranchAddress("D_Gamma_r10",&D_Gamma_r10);
		myTree[smp]->SetBranchAddress("D_Gamma_r15",&D_Gamma_r15);
		myTree[smp]->SetBranchAddress("D_Gamma_r20",&D_Gamma_r20);
		myTree[smp]->SetBranchAddress("D_Gamma_r25",&D_Gamma_r25);
		myTree[smp]->SetBranchAddress("D_Gamma_gg_r1",&D_Gamma_gg_r1);
		myTree[smp]->SetBranchAddress("D_Gamma_gg_r5",&D_Gamma_gg_r5);
		myTree[smp]->SetBranchAddress("D_Gamma_gg_r10",&D_Gamma_gg_r10);
		myTree[smp]->SetBranchAddress("D_Gamma_gg_r15",&D_Gamma_gg_r15);
		myTree[smp]->SetBranchAddress("D_Gamma_gg_r20",&D_Gamma_gg_r20);
		myTree[smp]->SetBranchAddress("D_Gamma_gg_r25",&D_Gamma_gg_r25);
		myTree[smp]->SetBranchAddress("D_Gamma",&D_Gamma);
		myTree[smp]->SetBranchAddress("D_Gamma_int",&D_Gamma_int);
*/	};

	for(int p=0;p<numSamples;p++){
/*
	float NGenTotal_weighted[kSignalSamples][numSamples];
	float NGenTotal_unweighted[kSignalSamples]={0};
	float NGenTotal_weighted_perFlavor[kSignalSamples][numSamples][nFinalStates];
	float fGen_BRweight_perFlavor[numSamples][nFinalStates]={{0}};
	float fGen_BRweight[numSamples]={0};
	float sum_NGenTotal_unweighted=0;
	float sum_NGenTotal_weighted[numSamples]={0};
*/
		double flavor_contribution[nFinalStates]={0};
		double flavor_contribution_4GeVcut[nFinalStates]={0};
		double sum_nevents=0;
		double sum_nevents_4GeVcut=0;
		for(int smp=0;smp<kSignalSamples;smp++){
			double total_weighted = NGenTotal_weighted[smp][p];
			double total_unweighted = NGenTotal_unweighted[smp];
			double total_4GeVcut_weighted = NGenTotal_4GeVcut_weighted[smp][p];
			double total_4GeVcut_unweighted = NGenTotal_4GeVcut_unweighted[smp];
			for(int f=0;f<nFinalStates;f++){
				flavor_contribution[f] += NGenTotal_weighted_perFlavor[smp][p][f] * total_unweighted/total_weighted;
				sum_nevents += NGenTotal_weighted_perFlavor[smp][p][f] * total_unweighted/total_weighted;
				flavor_contribution_4GeVcut[f] += NGenTotal_4GeVcut_weighted_perFlavor[smp][p][f] * total_4GeVcut_unweighted/total_4GeVcut_weighted;
				sum_nevents_4GeVcut += NGenTotal_4GeVcut_weighted_perFlavor[smp][p][f] * total_4GeVcut_unweighted/total_4GeVcut_weighted;
			};
		};
		for (int f = 0; f < nFinalStates; f++){
			fGen_BRweight_perFlavor[p][f] = flavor_contribution[f] / sum_nevents;
			fGen_4GeVcut_BRweight_perFlavor[p][f] = flavor_contribution_4GeVcut[f] / sum_nevents_4GeVcut;
		};
	};
	for(int p=0;p<numSamples;p++){
		cout << "BR 2e2mu in hypothesis " << p << ": " << fGen_BRweight_perFlavor[p][2] << endl;
		fGen_BRweight[p] = fGen_BRweight_perFlavor[0][2]/fGen_BRweight_perFlavor[p][2];
		cout << "BR 2e2mu (mZ>4 GeV) in hypothesis " << p << ": " << fGen_4GeVcut_BRweight_perFlavor[p][2] << endl;
		fGen_4GeVcut_BRweight[p] = fGen_4GeVcut_BRweight_perFlavor[0][2]/fGen_4GeVcut_BRweight_perFlavor[p][2];
	};

	string coutput = coutput_common + "/" + user_folder[flavor] + "/" + OUTPUT_NAME + ".root";
	TFile* foutput = new TFile(coutput.c_str(),"recreate");
	cout << "XSEC Ratios for Gen. Level No Interf" << endl;
	cout << "Hypothesis\tXSEC NO CUT\tXSEC WITH MZ>4 GEV" << endl;
	TH1F* hGenXsecRatios = new TH1F("hGen2e2muXsecRatios","2e2mu Xsec Ratios wrt SM for Each Hypothesis",0,myNumSamples,myNumSamples);
	TH1F* hGenXsecRatios_4GeVcut = new TH1F("hGen2e2muXsecRatios_4GeVcut","2e2mu Xsec Ratios wrt SM for Each Hypothesis After mZ>4 cut",0,myNumSamples,myNumSamples);
	double xsecSMScale_4GeVcut=1;
	for (int p = 0; p < myNumSamples; p++){
		double nWeighted=0,nUnweighted=0;
		double nWeighted_4GeVcut=0,nUnweighted_4GeVcut=0;

		double nWeightedSMRatio_4GeVcut=0,nUnweighted_Total=0;

		double xsecRatio=1;
		double xsecRatio_4GeVcut=1;

		for (int smp = 0; smp < (kSignalSamples); smp++){
			nWeighted += NGenTotal_unweighted_perFlavor[smp][2]*NGenTotal_weighted_perFlavor[smp][p][2]/NGenTotal_weighted_perFlavor[smp][0][2];
			nUnweighted += NGenTotal_unweighted_perFlavor[smp][2];
			nWeighted_4GeVcut += NGenTotal_4GeVcut_unweighted_perFlavor[smp][2]*NGenTotal_4GeVcut_weighted_perFlavor[smp][p][2]/NGenTotal_4GeVcut_weighted_perFlavor[smp][0][2];
			nUnweighted_4GeVcut += NGenTotal_4GeVcut_unweighted_perFlavor[smp][2];

/*
			for(int p=0;p<13;p++) cout << p << "\t" << nHypo[p] << "\t" <<  (nHypoPredicted[p]*fGen_BRweight[p]) << "\t" << nHypoPredicted[p] << "\t" << fGen_BRweight[p] << endl;
			for(int p=13;p<14;p++) cout << p << "\t" << nHypo[p] << "\tn/a\tn/a\tn/a" << endl;
			for(int p=14;p<myNumSamples;p++) cout << p << "\t" << nHypo[p] << "\t" << (nHypoPredicted[p-1]*fGen_BRweight[p-1]) << "\t" << nHypoPredicted[p-1] << "\t" << fGen_BRweight[p-1] << endl;
*/
			nWeightedSMRatio_4GeVcut += NGenTotal_unweighted[smp] * NGenTotal_4GeVcut_weighted[smp][0] / NGenTotal_weighted[smp][p];
			nUnweighted_Total += NGenTotal_unweighted[smp];
		};
		xsecRatio = nWeighted/nUnweighted;
		xsecRatio_4GeVcut = nWeighted_4GeVcut/nUnweighted_4GeVcut;

		if(p==0) xsecSMScale_4GeVcut = nWeightedSMRatio_4GeVcut/nUnweighted_Total;

		hGenXsecRatios->SetBinContent(p+1,xsecRatio);
		hGenXsecRatios_4GeVcut->SetBinContent(p+1,xsecRatio_4GeVcut);

		cout << p << '\t' << xsecRatio << '\t' << xsecRatio_4GeVcut << endl;
	};
	foutput->WriteTObject(hGenXsecRatios);
	foutput->WriteTObject(hGenXsecRatios_4GeVcut);
	delete hGenXsecRatios_4GeVcut;
	delete hGenXsecRatios;

	TTree* shuffledTree = new TTree(TREE_NAME,TREE_NAME);
	shuffledTree->SetAutoSave(3000000000);
	shuffledTree->Branch("EventSample",&EventSample);
	shuffledTree->Branch("kNumSamples",&myNumSamples);
	shuffledTree->Branch("MC_CV_weight",MC_CV_weight,"MC_CV_weight[kNumSamples]/F");
//	shuffledTree->Branch("MC_CVwrt2mu2e_weight",MC_CVwrt2mu2e_weight,"MC_CVwrt2mu2e_weight[kNumSamples]/F");
	shuffledTree->Branch("MC_weight_spin0",MC_weight_spin0,"MC_weight_spin0[kNumSamples]/F");
	shuffledTree->Branch("MC_weight",&MC_weight);
	shuffledTree->Branch("MC_weight_xsec",&MC_weight_xsec);
	shuffledTree->Branch("MC_CV_4GeVcut_weight",MC_CV_4GeVcut_weight,"MC_CV_4GeVcut_weight[kNumSamples]/F");
	shuffledTree->Branch("MC_weight_4GeVcut_spin0",MC_weight_4GeVcut_spin0,"MC_weight_4GeVcut_spin0[kNumSamples]/F");
	shuffledTree->Branch("MC_weight_4GeVcut",&MC_weight_4GeVcut);
	shuffledTree->Branch("MC_weight_xsec_4GeVcut",&MC_weight_xsec_4GeVcut);
	shuffledTree->Branch("MC_weight_Kfactor",&MC_weight_Kfactor);
	shuffledTree->Branch("MC_weight_Kfactor_down",&MC_weight_Kfactor_down);
	shuffledTree->Branch("MC_weight_Kfactor_up",&MC_weight_Kfactor_up);
	shuffledTree->Branch("GenHMass", &GenHMass);
//	shuffledTree->Branch("GenZ1Mass", &GenZ1Mass);
//	shuffledTree->Branch("GenZ2Mass", &GenZ2Mass);
	shuffledTree->Branch("ZZMass", &ZZMass);
	shuffledTree->Branch("Z1Mass", &Z1Mass);
	shuffledTree->Branch("Z2Mass", &Z2Mass);
	shuffledTree->Branch("helcosthetaZ1", &myhelcosthetaZ1);
	shuffledTree->Branch("helcosthetaZ2", &myhelcosthetaZ2);
	shuffledTree->Branch("helphi", &myhelphi);
	shuffledTree->Branch("costhetastar", &mycosthetastar);
	shuffledTree->Branch("phistarZ1", &myphistarZ1);
	shuffledTree->Branch("Z1ids", &Z1ids);
	shuffledTree->Branch("Z2ids", &Z2ids);
	shuffledTree->Branch("D_g1_vs_g2_phi0",&D_g2);
	shuffledTree->Branch("D_g1_vs_g4_phi0",&D_g4);
	shuffledTree->Branch("D_g2int_phi0",&D_g2int);
	shuffledTree->Branch("D_g4int_phi0",&D_g4int);
	shuffledTree->Branch("D_g2int_phi90",&D_g2int_perp);
	shuffledTree->Branch("D_g4int_phi90",&D_g4int_perp);
	shuffledTree->Branch("D_g1Q2_phi0",&D_g1q2);
	shuffledTree->Branch("D_g1Q2int_phi0",&D_g1q2int);
	shuffledTree->Branch("D_ZG",&D_ZG);
	shuffledTree->Branch("D_GG",&D_GG);
	shuffledTree->Branch("D_ZG_PS",&D_ZG_PS);
	shuffledTree->Branch("D_GG_PS",&D_GG_PS);
	shuffledTree->Branch("D_ZGint",&D_ZGint);
	shuffledTree->Branch("D_GGint",&D_GGint);
	shuffledTree->Branch("D_ZG_PSint",&D_ZG_PSint);
	shuffledTree->Branch("D_GG_PSint",&D_GG_PSint);
	shuffledTree->Branch("D_ZG_L1",&D_ZG_L1);
	shuffledTree->Branch("D_ZG_L1int_phi0",&D_ZG_L1int);
	shuffledTree->Branch("D_ZG_L1int_phi90",&D_ZG_L1int_perp);
	shuffledTree->Branch("D_bkg_kin",&D_bkg_kin);
	shuffledTree->Branch("D_bkg",&D_bkg);
	shuffledTree->Branch("D_bkg_ScaleUp",&D_bkg_ScaleUp);
	shuffledTree->Branch("D_bkg_ScaleDown",&D_bkg_ScaleDown);
	shuffledTree->Branch("D_bkg_ResUp",&D_bkg_ResUp);
	shuffledTree->Branch("D_bkg_ResDown",&D_bkg_ResDown);
/*	shuffledTree->Branch("D_ggZZ",&D_ggZZ);
	shuffledTree->Branch("D_Gamma_r1",&D_Gamma_r1);
	shuffledTree->Branch("D_Gamma_r5",&D_Gamma_r5);
	shuffledTree->Branch("D_Gamma_r10",&D_Gamma_r10);
	shuffledTree->Branch("D_Gamma_r15",&D_Gamma_r15);
	shuffledTree->Branch("D_Gamma_r20",&D_Gamma_r20);
	shuffledTree->Branch("D_Gamma_r25",&D_Gamma_r25);
	shuffledTree->Branch("D_Gamma_gg_r1",&D_Gamma_gg_r1);
	shuffledTree->Branch("D_Gamma_gg_r5",&D_Gamma_gg_r5);
	shuffledTree->Branch("D_Gamma_gg_r10",&D_Gamma_gg_r10);
	shuffledTree->Branch("D_Gamma_gg_r15",&D_Gamma_gg_r15);
	shuffledTree->Branch("D_Gamma_gg_r20",&D_Gamma_gg_r20);
	shuffledTree->Branch("D_Gamma_gg_r25",&D_Gamma_gg_r25);
	shuffledTree->Branch("D_Gamma",&D_Gamma);
	shuffledTree->Branch("D_Gamma_int",&D_Gamma_int);
*/
	TString TREE_NAME_GGZZ = TREE_NAME;
	TREE_NAME_GGZZ = TREE_NAME_GGZZ + "_ggZZ";
	TTree* shuffledggZZBkg = new TTree(TREE_NAME_GGZZ,TREE_NAME_GGZZ);
	shuffledggZZBkg->SetAutoSave(3000000000);
	shuffledggZZBkg->Branch("MC_weight",&MC_weight);
	shuffledggZZBkg->Branch("MC_weight_noxsec",&MC_weight_noxsec);
	shuffledggZZBkg->Branch("MC_weight_down",&MC_weight_down);
	shuffledggZZBkg->Branch("MC_weight_up",&MC_weight_up);
	shuffledggZZBkg->Branch("MC_weight_xsec",&MC_weight_xsec);
	shuffledggZZBkg->Branch("MC_weight_Kfactor",&MC_weight_Kfactor);
	shuffledggZZBkg->Branch("MC_weight_Kfactor_down",&MC_weight_Kfactor_down);
	shuffledggZZBkg->Branch("MC_weight_Kfactor_up",&MC_weight_Kfactor_up);
	shuffledggZZBkg->Branch("kNumQQBGGWeights",&kNumQQBGGWeights);
	shuffledggZZBkg->Branch("MC_weight_QQBGGProper",MC_weight_QQBGGProper,"MC_weight_QQBGGProper[kNumQQBGGWeights]/F");
	shuffledggZZBkg->Branch("ZZMass", &ZZMass);
	shuffledggZZBkg->Branch("Z1Mass", &Z1Mass);
	shuffledggZZBkg->Branch("Z2Mass", &Z2Mass);
	shuffledggZZBkg->Branch("helcosthetaZ1", &myhelcosthetaZ1);
	shuffledggZZBkg->Branch("helcosthetaZ2", &myhelcosthetaZ2);
	shuffledggZZBkg->Branch("helphi", &myhelphi);
	shuffledggZZBkg->Branch("costhetastar", &mycosthetastar);
	shuffledggZZBkg->Branch("phistarZ1", &myphistarZ1);
	shuffledggZZBkg->Branch("Z1ids", &Z1ids);
	shuffledggZZBkg->Branch("Z2ids", &Z2ids);
	shuffledggZZBkg->Branch("D_g1_vs_g2_phi0",&D_g2);
	shuffledggZZBkg->Branch("D_g1_vs_g4_phi0",&D_g4);
	shuffledggZZBkg->Branch("D_g2int_phi0",&D_g2int);
	shuffledggZZBkg->Branch("D_g4int_phi0",&D_g4int);
	shuffledggZZBkg->Branch("D_g2int_phi90",&D_g2int_perp);
	shuffledggZZBkg->Branch("D_g4int_phi90",&D_g4int_perp);
	shuffledggZZBkg->Branch("D_g1Q2_phi0",&D_g1q2);
	shuffledggZZBkg->Branch("D_g1Q2int_phi0",&D_g1q2int);
	shuffledggZZBkg->Branch("D_ZG",&D_ZG);
	shuffledggZZBkg->Branch("D_GG",&D_GG);
	shuffledggZZBkg->Branch("D_ZG_PS",&D_ZG_PS);
	shuffledggZZBkg->Branch("D_GG_PS",&D_GG_PS);
	shuffledggZZBkg->Branch("D_ZGint",&D_ZGint);
	shuffledggZZBkg->Branch("D_GGint",&D_GGint);
	shuffledggZZBkg->Branch("D_ZG_PSint",&D_ZG_PSint);
	shuffledggZZBkg->Branch("D_GG_PSint",&D_GG_PSint);
	shuffledggZZBkg->Branch("D_ZG_L1",&D_ZG_L1);
	shuffledggZZBkg->Branch("D_ZG_L1int_phi0",&D_ZG_L1int);
	shuffledggZZBkg->Branch("D_ZG_L1int_phi90",&D_ZG_L1int_perp);
	shuffledggZZBkg->Branch("D_bkg_kin",&D_bkg_kin);
	shuffledggZZBkg->Branch("D_bkg",&D_bkg);
	shuffledggZZBkg->Branch("D_bkg_ScaleUp",&D_bkg_ScaleUp);
	shuffledggZZBkg->Branch("D_bkg_ScaleDown",&D_bkg_ScaleDown);
	shuffledggZZBkg->Branch("D_bkg_ResUp",&D_bkg_ResUp);
	shuffledggZZBkg->Branch("D_bkg_ResDown",&D_bkg_ResDown);
/*	shuffledggZZBkg->Branch("D_ggZZ",&D_ggZZ);
	shuffledggZZBkg->Branch("D_Gamma_r1",&D_Gamma_r1);
	shuffledggZZBkg->Branch("D_Gamma_r5",&D_Gamma_r5);
	shuffledggZZBkg->Branch("D_Gamma_r10",&D_Gamma_r10);
	shuffledggZZBkg->Branch("D_Gamma_r15",&D_Gamma_r15);
	shuffledggZZBkg->Branch("D_Gamma_r20",&D_Gamma_r20);
	shuffledggZZBkg->Branch("D_Gamma_r25",&D_Gamma_r25);
	shuffledggZZBkg->Branch("D_Gamma_gg_r1",&D_Gamma_gg_r1);
	shuffledggZZBkg->Branch("D_Gamma_gg_r5",&D_Gamma_gg_r5);
	shuffledggZZBkg->Branch("D_Gamma_gg_r10",&D_Gamma_gg_r10);
	shuffledggZZBkg->Branch("D_Gamma_gg_r15",&D_Gamma_gg_r15);
	shuffledggZZBkg->Branch("D_Gamma_gg_r20",&D_Gamma_gg_r20);
	shuffledggZZBkg->Branch("D_Gamma_gg_r25",&D_Gamma_gg_r25);
	shuffledggZZBkg->Branch("D_Gamma",&D_Gamma);
	shuffledggZZBkg->Branch("D_Gamma_int",&D_Gamma_int);
*/
	shuffledggZZBkg->Branch("D_bkg_prodIndep", &D_bkg_prodIndep);
	shuffledggZZBkg->Branch("D_bkg_prodIndep_ScaleUp", &D_bkg_prodIndep_ScaleUp);
	shuffledggZZBkg->Branch("D_bkg_prodIndep_ScaleDown", &D_bkg_prodIndep_ScaleDown);
	shuffledggZZBkg->Branch("D_bkg_prodIndep_ResUp", &D_bkg_prodIndep_ResUp);
	shuffledggZZBkg->Branch("D_bkg_prodIndep_ResDown", &D_bkg_prodIndep_ResDown);
	shuffledggZZBkg->Branch("p1plusProdIndepKD", &p1plusProdIndepKD);
	shuffledggZZBkg->Branch("p1minusProdIndepKD", &p1minusProdIndepKD);
	shuffledggZZBkg->Branch("p1plusKD", &p1plusKD);
	shuffledggZZBkg->Branch("p1minusKD", &p1minusKD);
	shuffledggZZBkg->Branch("graviKD", &graviKD);
	shuffledggZZBkg->Branch("qqgraviKD", &qqgraviKD);
	shuffledggZZBkg->Branch("p2h2plusKD", &p2h2plusKD);
	shuffledggZZBkg->Branch("p2h2plus_qqb_KD", &p2h2plus_qqb_KD);
	shuffledggZZBkg->Branch("p2h3plusKD", &p2h3plusKD);
	shuffledggZZBkg->Branch("p2h3plus_qqb_KD", &p2h3plus_qqb_KD);
	shuffledggZZBkg->Branch("p2hplusKD", &p2hplusKD);
	shuffledggZZBkg->Branch("p2hplus_qqb_KD", &p2hplus_qqb_KD);
	shuffledggZZBkg->Branch("p2bplusKD", &p2bplusKD);
	shuffledggZZBkg->Branch("p2bplus_qqb_KD", &p2bplus_qqb_KD);
	shuffledggZZBkg->Branch("p2h6plusKD", &p2h6plusKD);
	shuffledggZZBkg->Branch("p2h6plus_qqb_KD", &p2h6plus_qqb_KD);
	shuffledggZZBkg->Branch("p2h7plusKD", &p2h7plusKD);
	shuffledggZZBkg->Branch("p2h7plus_qqb_KD", &p2h7plus_qqb_KD);
	shuffledggZZBkg->Branch("p2hminusKD", &p2hminusKD);
	shuffledggZZBkg->Branch("p2hminus_qqb_KD", &p2hminus_qqb_KD);
	shuffledggZZBkg->Branch("p2h9minusKD", &p2h9minusKD);
	shuffledggZZBkg->Branch("p2h9minus_qqb_KD", &p2h9minus_qqb_KD);
	shuffledggZZBkg->Branch("p2h10minusKD", &p2h10minusKD);
	shuffledggZZBkg->Branch("p2h10minus_qqb_KD", &p2h10minus_qqb_KD);
	shuffledggZZBkg->Branch("p2mProdIndepKD", &p2mProdIndepKD);
	shuffledggZZBkg->Branch("p2h2plusProdIndepKD", &p2h2plusProdIndepKD);
	shuffledggZZBkg->Branch("p2h3plusProdIndepKD", &p2h3plusProdIndepKD);
	shuffledggZZBkg->Branch("p2hplusProdIndepKD", &p2hplusProdIndepKD);
	shuffledggZZBkg->Branch("p2bplusProdIndepKD", &p2bplusProdIndepKD);
	shuffledggZZBkg->Branch("p2h6plusProdIndepKD", &p2h6plusProdIndepKD);
	shuffledggZZBkg->Branch("p2h7plusProdIndepKD", &p2h7plusProdIndepKD);
	shuffledggZZBkg->Branch("p2hminusProdIndepKD", &p2hminusProdIndepKD);
	shuffledggZZBkg->Branch("p2h9minusProdIndepKD", &p2h9minusProdIndepKD);
	shuffledggZZBkg->Branch("p2h10minusProdIndepKD", &p2h10minusProdIndepKD);

	TString TREE_NAME_QQZZ = TREE_NAME;
	TREE_NAME_QQZZ = TREE_NAME_QQZZ + "_qqZZ";
	TTree* shuffledqqZZBkg = new TTree(TREE_NAME_QQZZ,TREE_NAME_QQZZ);
	shuffledqqZZBkg->SetAutoSave(3000000000);
	shuffledqqZZBkg->Branch("MC_weight",&MC_weight);
	shuffledqqZZBkg->Branch("MC_weight_noxsec",&MC_weight_noxsec);
	shuffledqqZZBkg->Branch("MC_weight_down",&MC_weight_down);
	shuffledqqZZBkg->Branch("MC_weight_up",&MC_weight_up);
	shuffledqqZZBkg->Branch("MC_weight_xsec",&MC_weight_xsec);
	shuffledqqZZBkg->Branch("MC_weight_Kfactor",&MC_weight_Kfactor);
	shuffledqqZZBkg->Branch("MC_weight_Kfactor_down",&MC_weight_Kfactor_down);
	shuffledqqZZBkg->Branch("MC_weight_Kfactor_up",&MC_weight_Kfactor_up);
	shuffledqqZZBkg->Branch("kNumQQBGGWeights",&kNumQQBGGWeights);
	shuffledqqZZBkg->Branch("MC_weight_QQBGGProper",MC_weight_QQBGGProper,"MC_weight_QQBGGProper[kNumQQBGGWeights]/F");
	shuffledqqZZBkg->Branch("ZZMass", &ZZMass);
	shuffledqqZZBkg->Branch("Z1Mass", &Z1Mass);
	shuffledqqZZBkg->Branch("Z2Mass", &Z2Mass);
	shuffledqqZZBkg->Branch("helcosthetaZ1", &myhelcosthetaZ1);
	shuffledqqZZBkg->Branch("helcosthetaZ2", &myhelcosthetaZ2);
	shuffledqqZZBkg->Branch("helphi", &myhelphi);
	shuffledqqZZBkg->Branch("costhetastar", &mycosthetastar);
	shuffledqqZZBkg->Branch("phistarZ1", &myphistarZ1);
	shuffledqqZZBkg->Branch("Z1ids", &Z1ids);
	shuffledqqZZBkg->Branch("Z2ids", &Z2ids);
	shuffledqqZZBkg->Branch("D_g1_vs_g2_phi0",&D_g2);
	shuffledqqZZBkg->Branch("D_g1_vs_g4_phi0",&D_g4);
	shuffledqqZZBkg->Branch("D_g2int_phi0",&D_g2int);
	shuffledqqZZBkg->Branch("D_g4int_phi0",&D_g4int);
	shuffledqqZZBkg->Branch("D_g2int_phi90",&D_g2int_perp);
	shuffledqqZZBkg->Branch("D_g4int_phi90",&D_g4int_perp);
	shuffledqqZZBkg->Branch("D_g1Q2_phi0",&D_g1q2);
	shuffledqqZZBkg->Branch("D_g1Q2int_phi0",&D_g1q2int);
	shuffledqqZZBkg->Branch("D_ZG",&D_ZG);
	shuffledqqZZBkg->Branch("D_GG",&D_GG);
	shuffledqqZZBkg->Branch("D_ZG_PS",&D_ZG_PS);
	shuffledqqZZBkg->Branch("D_GG_PS",&D_GG_PS);
	shuffledqqZZBkg->Branch("D_ZGint",&D_ZGint);
	shuffledqqZZBkg->Branch("D_GGint",&D_GGint);
	shuffledqqZZBkg->Branch("D_ZG_PSint",&D_ZG_PSint);
	shuffledqqZZBkg->Branch("D_GG_PSint",&D_GG_PSint);
	shuffledqqZZBkg->Branch("D_ZG_L1",&D_ZG_L1);
	shuffledqqZZBkg->Branch("D_ZG_L1int_phi0",&D_ZG_L1int);
	shuffledqqZZBkg->Branch("D_ZG_L1int_phi90",&D_ZG_L1int_perp);
	shuffledqqZZBkg->Branch("D_bkg_kin",&D_bkg_kin);
	shuffledqqZZBkg->Branch("D_bkg",&D_bkg);
	shuffledqqZZBkg->Branch("D_bkg_ScaleUp",&D_bkg_ScaleUp);
	shuffledqqZZBkg->Branch("D_bkg_ScaleDown",&D_bkg_ScaleDown);
	shuffledqqZZBkg->Branch("D_bkg_ResUp",&D_bkg_ResUp);
	shuffledqqZZBkg->Branch("D_bkg_ResDown",&D_bkg_ResDown);
/*	shuffledqqZZBkg->Branch("D_ggZZ",&D_ggZZ);
	shuffledqqZZBkg->Branch("D_Gamma_r1",&D_Gamma_r1);
	shuffledqqZZBkg->Branch("D_Gamma_r5",&D_Gamma_r5);
	shuffledqqZZBkg->Branch("D_Gamma_r10",&D_Gamma_r10);
	shuffledqqZZBkg->Branch("D_Gamma_r15",&D_Gamma_r15);
	shuffledqqZZBkg->Branch("D_Gamma_r20",&D_Gamma_r20);
	shuffledqqZZBkg->Branch("D_Gamma_r25",&D_Gamma_r25);
	shuffledqqZZBkg->Branch("D_Gamma_gg_r1",&D_Gamma_gg_r1);
	shuffledqqZZBkg->Branch("D_Gamma_gg_r5",&D_Gamma_gg_r5);
	shuffledqqZZBkg->Branch("D_Gamma_gg_r10",&D_Gamma_gg_r10);
	shuffledqqZZBkg->Branch("D_Gamma_gg_r15",&D_Gamma_gg_r15);
	shuffledqqZZBkg->Branch("D_Gamma_gg_r20",&D_Gamma_gg_r20);
	shuffledqqZZBkg->Branch("D_Gamma_gg_r25",&D_Gamma_gg_r25);
	shuffledqqZZBkg->Branch("D_Gamma",&D_Gamma);
	shuffledqqZZBkg->Branch("D_Gamma_int",&D_Gamma_int);
*/
	shuffledqqZZBkg->Branch("D_bkg_prodIndep", &D_bkg_prodIndep);
	shuffledqqZZBkg->Branch("D_bkg_prodIndep_ScaleUp", &D_bkg_prodIndep_ScaleUp);
	shuffledqqZZBkg->Branch("D_bkg_prodIndep_ScaleDown", &D_bkg_prodIndep_ScaleDown);
	shuffledqqZZBkg->Branch("D_bkg_prodIndep_ResUp", &D_bkg_prodIndep_ResUp);
	shuffledqqZZBkg->Branch("D_bkg_prodIndep_ResDown", &D_bkg_prodIndep_ResDown);
	shuffledqqZZBkg->Branch("p1plusProdIndepKD", &p1plusProdIndepKD);
	shuffledqqZZBkg->Branch("p1minusProdIndepKD", &p1minusProdIndepKD);
	shuffledqqZZBkg->Branch("p1plusKD", &p1plusKD);
	shuffledqqZZBkg->Branch("p1minusKD", &p1minusKD);
	shuffledqqZZBkg->Branch("graviKD", &graviKD);
	shuffledqqZZBkg->Branch("qqgraviKD", &qqgraviKD);
	shuffledqqZZBkg->Branch("p2h2plusKD", &p2h2plusKD);
	shuffledqqZZBkg->Branch("p2h2plus_qqb_KD", &p2h2plus_qqb_KD);
	shuffledqqZZBkg->Branch("p2h3plusKD", &p2h3plusKD);
	shuffledqqZZBkg->Branch("p2h3plus_qqb_KD", &p2h3plus_qqb_KD);
	shuffledqqZZBkg->Branch("p2hplusKD", &p2hplusKD);
	shuffledqqZZBkg->Branch("p2hplus_qqb_KD", &p2hplus_qqb_KD);
	shuffledqqZZBkg->Branch("p2bplusKD", &p2bplusKD);
	shuffledqqZZBkg->Branch("p2bplus_qqb_KD", &p2bplus_qqb_KD);
	shuffledqqZZBkg->Branch("p2h6plusKD", &p2h6plusKD);
	shuffledqqZZBkg->Branch("p2h6plus_qqb_KD", &p2h6plus_qqb_KD);
	shuffledqqZZBkg->Branch("p2h7plusKD", &p2h7plusKD);
	shuffledqqZZBkg->Branch("p2h7plus_qqb_KD", &p2h7plus_qqb_KD);
	shuffledqqZZBkg->Branch("p2hminusKD", &p2hminusKD);
	shuffledqqZZBkg->Branch("p2hminus_qqb_KD", &p2hminus_qqb_KD);
	shuffledqqZZBkg->Branch("p2h9minusKD", &p2h9minusKD);
	shuffledqqZZBkg->Branch("p2h9minus_qqb_KD", &p2h9minus_qqb_KD);
	shuffledqqZZBkg->Branch("p2h10minusKD", &p2h10minusKD);
	shuffledqqZZBkg->Branch("p2h10minus_qqb_KD", &p2h10minus_qqb_KD);
	shuffledqqZZBkg->Branch("p2mProdIndepKD", &p2mProdIndepKD);
	shuffledqqZZBkg->Branch("p2h2plusProdIndepKD", &p2h2plusProdIndepKD);
	shuffledqqZZBkg->Branch("p2h3plusProdIndepKD", &p2h3plusProdIndepKD);
	shuffledqqZZBkg->Branch("p2hplusProdIndepKD", &p2hplusProdIndepKD);
	shuffledqqZZBkg->Branch("p2bplusProdIndepKD", &p2bplusProdIndepKD);
	shuffledqqZZBkg->Branch("p2h6plusProdIndepKD", &p2h6plusProdIndepKD);
	shuffledqqZZBkg->Branch("p2h7plusProdIndepKD", &p2h7plusProdIndepKD);
	shuffledqqZZBkg->Branch("p2hminusProdIndepKD", &p2hminusProdIndepKD);
	shuffledqqZZBkg->Branch("p2h9minusProdIndepKD", &p2h9minusProdIndepKD);
	shuffledqqZZBkg->Branch("p2h10minusProdIndepKD", &p2h10minusProdIndepKD);

	TString TREE_NAME_ZX = TREE_NAME;
	TREE_NAME_ZX = TREE_NAME_ZX + "_ZX";
	TTree* shuffledZXBkg = new TTree(TREE_NAME_ZX,TREE_NAME_ZX);
	shuffledZXBkg->SetAutoSave(3000000000);
	shuffledZXBkg->Branch("ZXfake_weightProper",&ZXfake_weightProper);
	shuffledZXBkg->Branch("ZZMass", &ZZMass);
	shuffledZXBkg->Branch("Z1Mass", &Z1Mass);
	shuffledZXBkg->Branch("Z2Mass", &Z2Mass);
	shuffledZXBkg->Branch("helcosthetaZ1", &myhelcosthetaZ1);
	shuffledZXBkg->Branch("helcosthetaZ2", &myhelcosthetaZ2);
	shuffledZXBkg->Branch("helphi", &myhelphi);
	shuffledZXBkg->Branch("costhetastar", &mycosthetastar);
	shuffledZXBkg->Branch("phistarZ1", &myphistarZ1);
	shuffledZXBkg->Branch("Z1ids", &Z1ids);
	shuffledZXBkg->Branch("Z2ids", &Z2ids);
	shuffledZXBkg->Branch("D_g1_vs_g2_phi0",&D_g2);
	shuffledZXBkg->Branch("D_g1_vs_g4_phi0",&D_g4);
	shuffledZXBkg->Branch("D_g2int_phi0",&D_g2int);
	shuffledZXBkg->Branch("D_g4int_phi0",&D_g4int);
	shuffledZXBkg->Branch("D_g2int_phi90",&D_g2int_perp);
	shuffledZXBkg->Branch("D_g4int_phi90",&D_g4int_perp);
	shuffledZXBkg->Branch("D_g1Q2_phi0",&D_g1q2);
	shuffledZXBkg->Branch("D_g1Q2int_phi0",&D_g1q2int);
	shuffledZXBkg->Branch("D_ZG",&D_ZG);
	shuffledZXBkg->Branch("D_GG",&D_GG);
	shuffledZXBkg->Branch("D_ZG_PS",&D_ZG_PS);
	shuffledZXBkg->Branch("D_GG_PS",&D_GG_PS);
	shuffledZXBkg->Branch("D_ZGint",&D_ZGint);
	shuffledZXBkg->Branch("D_GGint",&D_GGint);
	shuffledZXBkg->Branch("D_ZG_PSint",&D_ZG_PSint);
	shuffledZXBkg->Branch("D_GG_PSint",&D_GG_PSint);
	shuffledZXBkg->Branch("D_ZG_L1",&D_ZG_L1);
	shuffledZXBkg->Branch("D_ZG_L1int_phi0",&D_ZG_L1int);
	shuffledZXBkg->Branch("D_ZG_L1int_phi90",&D_ZG_L1int_perp);
	shuffledZXBkg->Branch("D_bkg_kin",&D_bkg_kin);
	shuffledZXBkg->Branch("D_bkg",&D_bkg);
	shuffledZXBkg->Branch("D_bkg_ScaleUp",&D_bkg_ScaleUp);
	shuffledZXBkg->Branch("D_bkg_ScaleDown",&D_bkg_ScaleDown);
	shuffledZXBkg->Branch("D_bkg_ResUp",&D_bkg_ResUp);
	shuffledZXBkg->Branch("D_bkg_ResDown",&D_bkg_ResDown);
/*	shuffledZXBkg->Branch("D_ggZZ",&D_ggZZ);
	shuffledZXBkg->Branch("D_Gamma_r1",&D_Gamma_r1);
	shuffledZXBkg->Branch("D_Gamma_r5",&D_Gamma_r5);
	shuffledZXBkg->Branch("D_Gamma_r10",&D_Gamma_r10);
	shuffledZXBkg->Branch("D_Gamma_r15",&D_Gamma_r15);
	shuffledZXBkg->Branch("D_Gamma_r20",&D_Gamma_r20);
	shuffledZXBkg->Branch("D_Gamma_r25",&D_Gamma_r25);
	shuffledZXBkg->Branch("D_Gamma_gg_r1",&D_Gamma_gg_r1);
	shuffledZXBkg->Branch("D_Gamma_gg_r5",&D_Gamma_gg_r5);
	shuffledZXBkg->Branch("D_Gamma_gg_r10",&D_Gamma_gg_r10);
	shuffledZXBkg->Branch("D_Gamma_gg_r15",&D_Gamma_gg_r15);
	shuffledZXBkg->Branch("D_Gamma_gg_r20",&D_Gamma_gg_r20);
	shuffledZXBkg->Branch("D_Gamma_gg_r25",&D_Gamma_gg_r25);
	shuffledZXBkg->Branch("D_Gamma",&D_Gamma);
	shuffledZXBkg->Branch("D_Gamma_int",&D_Gamma_int);
*/
	shuffledZXBkg->Branch("D_bkg_prodIndep", &D_bkg_prodIndep);
	shuffledZXBkg->Branch("D_bkg_prodIndep_ScaleUp", &D_bkg_prodIndep_ScaleUp);
	shuffledZXBkg->Branch("D_bkg_prodIndep_ScaleDown", &D_bkg_prodIndep_ScaleDown);
	shuffledZXBkg->Branch("D_bkg_prodIndep_ResUp", &D_bkg_prodIndep_ResUp);
	shuffledZXBkg->Branch("D_bkg_prodIndep_ResDown", &D_bkg_prodIndep_ResDown);
	shuffledZXBkg->Branch("p1plusProdIndepKD", &p1plusProdIndepKD);
	shuffledZXBkg->Branch("p1minusProdIndepKD", &p1minusProdIndepKD);
	shuffledZXBkg->Branch("p1plusKD", &p1plusKD);
	shuffledZXBkg->Branch("p1minusKD", &p1minusKD);
	shuffledZXBkg->Branch("graviKD", &graviKD);
	shuffledZXBkg->Branch("qqgraviKD", &qqgraviKD);
	shuffledZXBkg->Branch("p2h2plusKD", &p2h2plusKD);
	shuffledZXBkg->Branch("p2h2plus_qqb_KD", &p2h2plus_qqb_KD);
	shuffledZXBkg->Branch("p2h3plusKD", &p2h3plusKD);
	shuffledZXBkg->Branch("p2h3plus_qqb_KD", &p2h3plus_qqb_KD);
	shuffledZXBkg->Branch("p2hplusKD", &p2hplusKD);
	shuffledZXBkg->Branch("p2hplus_qqb_KD", &p2hplus_qqb_KD);
	shuffledZXBkg->Branch("p2bplusKD", &p2bplusKD);
	shuffledZXBkg->Branch("p2bplus_qqb_KD", &p2bplus_qqb_KD);
	shuffledZXBkg->Branch("p2h6plusKD", &p2h6plusKD);
	shuffledZXBkg->Branch("p2h6plus_qqb_KD", &p2h6plus_qqb_KD);
	shuffledZXBkg->Branch("p2h7plusKD", &p2h7plusKD);
	shuffledZXBkg->Branch("p2h7plus_qqb_KD", &p2h7plus_qqb_KD);
	shuffledZXBkg->Branch("p2hminusKD", &p2hminusKD);
	shuffledZXBkg->Branch("p2hminus_qqb_KD", &p2hminus_qqb_KD);
	shuffledZXBkg->Branch("p2h9minusKD", &p2h9minusKD);
	shuffledZXBkg->Branch("p2h9minus_qqb_KD", &p2h9minus_qqb_KD);
	shuffledZXBkg->Branch("p2h10minusKD", &p2h10minusKD);
	shuffledZXBkg->Branch("p2h10minus_qqb_KD", &p2h10minus_qqb_KD);
	shuffledZXBkg->Branch("p2mProdIndepKD", &p2mProdIndepKD);
	shuffledZXBkg->Branch("p2h2plusProdIndepKD", &p2h2plusProdIndepKD);
	shuffledZXBkg->Branch("p2h3plusProdIndepKD", &p2h3plusProdIndepKD);
	shuffledZXBkg->Branch("p2hplusProdIndepKD", &p2hplusProdIndepKD);
	shuffledZXBkg->Branch("p2bplusProdIndepKD", &p2bplusProdIndepKD);
	shuffledZXBkg->Branch("p2h6plusProdIndepKD", &p2h6plusProdIndepKD);
	shuffledZXBkg->Branch("p2h7plusProdIndepKD", &p2h7plusProdIndepKD);
	shuffledZXBkg->Branch("p2hminusProdIndepKD", &p2hminusProdIndepKD);
	shuffledZXBkg->Branch("p2h9minusProdIndepKD", &p2h9minusProdIndepKD);
	shuffledZXBkg->Branch("p2h10minusProdIndepKD", &p2h10minusProdIndepKD);


	TString TREE_NAME_QQZZ_DEDICATED = TREE_NAME;
	TREE_NAME_QQZZ_DEDICATED = TREE_NAME_QQZZ_DEDICATED + "_qqZZ_Dedicated";
	TTree* shuffledqqZZ_DedicatedBkg = new TTree(TREE_NAME_QQZZ_DEDICATED,TREE_NAME_QQZZ_DEDICATED);
	shuffledqqZZ_DedicatedBkg->SetAutoSave(3000000000);
	shuffledqqZZ_DedicatedBkg->Branch("MC_weight",&MC_weight);
	shuffledqqZZ_DedicatedBkg->Branch("MC_weight_noxsec",&MC_weight_noxsec);
	shuffledqqZZ_DedicatedBkg->Branch("MC_weight_down",&MC_weight_down);
	shuffledqqZZ_DedicatedBkg->Branch("MC_weight_up",&MC_weight_up);
	shuffledqqZZ_DedicatedBkg->Branch("MC_weight_xsec",&MC_weight_xsec);
	shuffledqqZZ_DedicatedBkg->Branch("MC_weight_Kfactor",&MC_weight_Kfactor);
	shuffledqqZZ_DedicatedBkg->Branch("MC_weight_Kfactor_down",&MC_weight_Kfactor_down);
	shuffledqqZZ_DedicatedBkg->Branch("MC_weight_Kfactor_up",&MC_weight_Kfactor_up);
	shuffledqqZZ_DedicatedBkg->Branch("kNumQQBGGWeights",&kNumQQBGGWeights);
	shuffledqqZZ_DedicatedBkg->Branch("MC_weight_QQBGGProper",MC_weight_QQBGGProper,"MC_weight_QQBGGProper[kNumQQBGGWeights]/F");
	shuffledqqZZ_DedicatedBkg->Branch("ZZMass", &ZZMass);
	shuffledqqZZ_DedicatedBkg->Branch("Z1Mass", &Z1Mass);
	shuffledqqZZ_DedicatedBkg->Branch("Z2Mass", &Z2Mass);
	shuffledqqZZ_DedicatedBkg->Branch("helcosthetaZ1", &myhelcosthetaZ1);
	shuffledqqZZ_DedicatedBkg->Branch("helcosthetaZ2", &myhelcosthetaZ2);
	shuffledqqZZ_DedicatedBkg->Branch("helphi", &myhelphi);
	shuffledqqZZ_DedicatedBkg->Branch("costhetastar", &mycosthetastar);
	shuffledqqZZ_DedicatedBkg->Branch("phistarZ1", &myphistarZ1);
	shuffledqqZZ_DedicatedBkg->Branch("Z1ids", &Z1ids);
	shuffledqqZZ_DedicatedBkg->Branch("Z2ids", &Z2ids);
	shuffledqqZZ_DedicatedBkg->Branch("D_g1_vs_g2_phi0",&D_g2);
	shuffledqqZZ_DedicatedBkg->Branch("D_g1_vs_g4_phi0",&D_g4);
	shuffledqqZZ_DedicatedBkg->Branch("D_g2int_phi0",&D_g2int);
	shuffledqqZZ_DedicatedBkg->Branch("D_g4int_phi0",&D_g4int);
	shuffledqqZZ_DedicatedBkg->Branch("D_g2int_phi90",&D_g2int_perp);
	shuffledqqZZ_DedicatedBkg->Branch("D_g4int_phi90",&D_g4int_perp);
	shuffledqqZZ_DedicatedBkg->Branch("D_g1Q2_phi0",&D_g1q2);
	shuffledqqZZ_DedicatedBkg->Branch("D_g1Q2int_phi0",&D_g1q2int);
	shuffledqqZZ_DedicatedBkg->Branch("D_ZG",&D_ZG);
	shuffledqqZZ_DedicatedBkg->Branch("D_GG",&D_GG);
	shuffledqqZZ_DedicatedBkg->Branch("D_ZG_PS",&D_ZG_PS);
	shuffledqqZZ_DedicatedBkg->Branch("D_GG_PS",&D_GG_PS);
	shuffledqqZZ_DedicatedBkg->Branch("D_ZGint",&D_ZGint);
	shuffledqqZZ_DedicatedBkg->Branch("D_GGint",&D_GGint);
	shuffledqqZZ_DedicatedBkg->Branch("D_ZG_PSint",&D_ZG_PSint);
	shuffledqqZZ_DedicatedBkg->Branch("D_GG_PSint",&D_GG_PSint);
	shuffledqqZZ_DedicatedBkg->Branch("D_ZG_L1",&D_ZG_L1);
	shuffledqqZZ_DedicatedBkg->Branch("D_ZG_L1int_phi0",&D_ZG_L1int);
	shuffledqqZZ_DedicatedBkg->Branch("D_ZG_L1int_phi90",&D_ZG_L1int_perp);
	shuffledqqZZ_DedicatedBkg->Branch("D_bkg_kin",&D_bkg_kin);
	shuffledqqZZ_DedicatedBkg->Branch("D_bkg",&D_bkg);
	shuffledqqZZ_DedicatedBkg->Branch("D_bkg_ScaleUp",&D_bkg_ScaleUp);
	shuffledqqZZ_DedicatedBkg->Branch("D_bkg_ScaleDown",&D_bkg_ScaleDown);
	shuffledqqZZ_DedicatedBkg->Branch("D_bkg_ResUp",&D_bkg_ResUp);
	shuffledqqZZ_DedicatedBkg->Branch("D_bkg_ResDown",&D_bkg_ResDown);
/*	shuffledqqZZ_DedicatedBkg->Branch("D_ggZZ",&D_ggZZ);
	shuffledqqZZ_DedicatedBkg->Branch("D_Gamma_r1",&D_Gamma_r1);
	shuffledqqZZ_DedicatedBkg->Branch("D_Gamma_r5",&D_Gamma_r5);
	shuffledqqZZ_DedicatedBkg->Branch("D_Gamma_r10",&D_Gamma_r10);
	shuffledqqZZ_DedicatedBkg->Branch("D_Gamma_r15",&D_Gamma_r15);
	shuffledqqZZ_DedicatedBkg->Branch("D_Gamma_r20",&D_Gamma_r20);
	shuffledqqZZ_DedicatedBkg->Branch("D_Gamma_r25",&D_Gamma_r25);
	shuffledqqZZ_DedicatedBkg->Branch("D_Gamma_gg_r1",&D_Gamma_gg_r1);
	shuffledqqZZ_DedicatedBkg->Branch("D_Gamma_gg_r5",&D_Gamma_gg_r5);
	shuffledqqZZ_DedicatedBkg->Branch("D_Gamma_gg_r10",&D_Gamma_gg_r10);
	shuffledqqZZ_DedicatedBkg->Branch("D_Gamma_gg_r15",&D_Gamma_gg_r15);
	shuffledqqZZ_DedicatedBkg->Branch("D_Gamma_gg_r20",&D_Gamma_gg_r20);
	shuffledqqZZ_DedicatedBkg->Branch("D_Gamma_gg_r25",&D_Gamma_gg_r25);
	shuffledqqZZ_DedicatedBkg->Branch("D_Gamma",&D_Gamma);
	shuffledqqZZ_DedicatedBkg->Branch("D_Gamma_int",&D_Gamma_int);
*/
	shuffledqqZZ_DedicatedBkg->Branch("D_bkg_prodIndep", &D_bkg_prodIndep);
	shuffledqqZZ_DedicatedBkg->Branch("D_bkg_prodIndep_ScaleUp", &D_bkg_prodIndep_ScaleUp);
	shuffledqqZZ_DedicatedBkg->Branch("D_bkg_prodIndep_ScaleDown", &D_bkg_prodIndep_ScaleDown);
	shuffledqqZZ_DedicatedBkg->Branch("D_bkg_prodIndep_ResUp", &D_bkg_prodIndep_ResUp);
	shuffledqqZZ_DedicatedBkg->Branch("D_bkg_prodIndep_ResDown", &D_bkg_prodIndep_ResDown);
	shuffledqqZZ_DedicatedBkg->Branch("p1plusProdIndepKD", &p1plusProdIndepKD);
	shuffledqqZZ_DedicatedBkg->Branch("p1minusProdIndepKD", &p1minusProdIndepKD);
	shuffledqqZZ_DedicatedBkg->Branch("p1plusKD", &p1plusKD);
	shuffledqqZZ_DedicatedBkg->Branch("p1minusKD", &p1minusKD);
	shuffledqqZZ_DedicatedBkg->Branch("graviKD", &graviKD);
	shuffledqqZZ_DedicatedBkg->Branch("qqgraviKD", &qqgraviKD);
	shuffledqqZZ_DedicatedBkg->Branch("p2h2plusKD", &p2h2plusKD);
	shuffledqqZZ_DedicatedBkg->Branch("p2h2plus_qqb_KD", &p2h2plus_qqb_KD);
	shuffledqqZZ_DedicatedBkg->Branch("p2h3plusKD", &p2h3plusKD);
	shuffledqqZZ_DedicatedBkg->Branch("p2h3plus_qqb_KD", &p2h3plus_qqb_KD);
	shuffledqqZZ_DedicatedBkg->Branch("p2hplusKD", &p2hplusKD);
	shuffledqqZZ_DedicatedBkg->Branch("p2hplus_qqb_KD", &p2hplus_qqb_KD);
	shuffledqqZZ_DedicatedBkg->Branch("p2bplusKD", &p2bplusKD);
	shuffledqqZZ_DedicatedBkg->Branch("p2bplus_qqb_KD", &p2bplus_qqb_KD);
	shuffledqqZZ_DedicatedBkg->Branch("p2h6plusKD", &p2h6plusKD);
	shuffledqqZZ_DedicatedBkg->Branch("p2h6plus_qqb_KD", &p2h6plus_qqb_KD);
	shuffledqqZZ_DedicatedBkg->Branch("p2h7plusKD", &p2h7plusKD);
	shuffledqqZZ_DedicatedBkg->Branch("p2h7plus_qqb_KD", &p2h7plus_qqb_KD);
	shuffledqqZZ_DedicatedBkg->Branch("p2hminusKD", &p2hminusKD);
	shuffledqqZZ_DedicatedBkg->Branch("p2hminus_qqb_KD", &p2hminus_qqb_KD);
	shuffledqqZZ_DedicatedBkg->Branch("p2h9minusKD", &p2h9minusKD);
	shuffledqqZZ_DedicatedBkg->Branch("p2h9minus_qqb_KD", &p2h9minus_qqb_KD);
	shuffledqqZZ_DedicatedBkg->Branch("p2h10minusKD", &p2h10minusKD);
	shuffledqqZZ_DedicatedBkg->Branch("p2h10minus_qqb_KD", &p2h10minus_qqb_KD);
	shuffledqqZZ_DedicatedBkg->Branch("p2mProdIndepKD", &p2mProdIndepKD);
	shuffledqqZZ_DedicatedBkg->Branch("p2h2plusProdIndepKD", &p2h2plusProdIndepKD);
	shuffledqqZZ_DedicatedBkg->Branch("p2h3plusProdIndepKD", &p2h3plusProdIndepKD);
	shuffledqqZZ_DedicatedBkg->Branch("p2hplusProdIndepKD", &p2hplusProdIndepKD);
	shuffledqqZZ_DedicatedBkg->Branch("p2bplusProdIndepKD", &p2bplusProdIndepKD);
	shuffledqqZZ_DedicatedBkg->Branch("p2h6plusProdIndepKD", &p2h6plusProdIndepKD);
	shuffledqqZZ_DedicatedBkg->Branch("p2h7plusProdIndepKD", &p2h7plusProdIndepKD);
	shuffledqqZZ_DedicatedBkg->Branch("p2hminusProdIndepKD", &p2hminusProdIndepKD);
	shuffledqqZZ_DedicatedBkg->Branch("p2h9minusProdIndepKD", &p2h9minusProdIndepKD);
	shuffledqqZZ_DedicatedBkg->Branch("p2h10minusProdIndepKD", &p2h10minusProdIndepKD);

	cout << "Starting to shuffle signal trees... Total N_events is " << sum_NGenTotal_unweighted << "." << endl;
	cout << "External weights for hypo=SM:" << endl;
	cout << "N_total Unweighted\nFile\tWeight" << endl;
	for(int l=0;l<kSignalSamples;l++) cout << l << '\t' << NGenTotal_unweighted[l] << endl;
	cout << "N_total Weighted\nFile\tWeight" << endl;
	for(int l=0;l<kSignalSamples;l++) cout << l << '\t' << NGenTotal_weighted[l][0] << endl;

	double nSM[kSignalSamples]={0};
	double nSM_4GeVcut[kSignalSamples]={0};

	double nHypoPredicted[kSignalSamples]={0};
	double nHypo[numSamples+2]={0};
	double nHypoPredicted_4GeVcut[kSignalSamples]={0};
	double nHypo_4GeVcut[numSamples+2]={0};
	double nHypo_125p6MassWin[numSamples]={0};
	double nHypo_125p6RestrictedMassWin[numSamples]={0};
	double nHypo_4GeVcut_125p6MassWin[numSamples]={0};
	double nHypo_4GeVcut_125p6RestrictedMassWin[numSamples]={0};

	double nGGZZ=0;
	double nQQZZ=0;
	double nQQZZ_Dedicated=0;
	double nGGZZ160=0;
	double nQQZZ160=0;
	double nQQZZ_Dedicated160=0;
	double nZX=0;
	double nGGZZReWgt=0;
	double nGGZZReWgt160=0;
	double nGGZZReWgt_in160=0;
	double nQQZZReWgt=0;
	double nQQZZReWgt_in160=0;
	double nQQZZ_DedicatedReWgt=0;
	
	int nUsed[(kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples+kqqZZSamples_Dedicated)]={0};

	cout << "4 GeV cut xsec scale is " << xsecSMScale_4GeVcut << endl;

	int ctr=0;
	while(ctr<nTotal_Sig){
		int coin = rand() % naccevents[(kSignalSamples-1)];
		int lucky=-1;
		for(int smp=0;smp<(kSignalSamples);smp++){
			if(coin<naccevents[smp] && (nevents[smp]-nUsed[smp])!=0){lucky=smp;break;};
		};
		if(lucky<0) continue;

		myTree[lucky]->GetEntry(nUsed[lucky]);
		nUsed[lucky] += 1;
		EventSample = lucky;
		MC_weight_xsec = MC_weight/MC_weight_noxsec;
		MC_weight_xsec_4GeVcut = MC_weight_xsec*xsecSMScale_4GeVcut;
		MC_weight_4GeVcut = MC_weight_xsec_4GeVcut*MC_weight_noxsec;

		if(lucky<kSignalSamples){
			for(int p=0;p<myNumSamples;p++){
				MC_weight_4GeVcut_spin0[p] = MC_weight_spin0[p]*( NGenTotal_4GeVcut_unweighted[lucky]/NGenTotal_4GeVcut_weighted[lucky][p]*fGen_4GeVcut_BRweight[p] );
				MC_weight_spin0[p] *= ( NGenTotal_unweighted[lucky]/NGenTotal_weighted[lucky][p]*fGen_BRweight[p] );
				double myxsec = MC_weight_xsec * MC_weight_spin0[p]  / sum_NGenTotal_unweighted ;
				double myxsec_4GeVcut = MC_weight_xsec_4GeVcut * MC_weight_4GeVcut_spin0[p]  / sum_NGenTotal_4GeVcut_unweighted ;
				MC_CV_weight[p] = myxsec*MC_weight_PUWeight*MC_weight_powhegWeight*MC_weight_dataMC*MC_weight_HqT;
				MC_CV_4GeVcut_weight[p] = myxsec_4GeVcut*MC_weight_PUWeight*MC_weight_powhegWeight*MC_weight_dataMC*MC_weight_HqT;
				if(!(GenZ1Mass>4 && GenZ2Mass>4)) MC_CV_4GeVcut_weight[p]=0;

				if(MC_CV_weight[p]>=1.48e-5) MC_CV_weight[p] = pow(1.48e-5,2)/MC_CV_weight[p];
				if(MC_CV_4GeVcut_weight[p]>=1.48e-5) MC_CV_4GeVcut_weight[p] = pow(1.48e-5,2)/MC_CV_4GeVcut_weight[p];

				double ratioAsIf2mu2e = 0.5;
				if(flavor==2) ratioAsIf2mu2e = 1;
				ratioAsIf2mu2e *= NGenTotal_weighted_perFlavor[lucky][p][2];
				ratioAsIf2mu2e = NGenTotal_weighted_perFlavor[lucky][p][flavor] / ratioAsIf2mu2e;

				nHypo[p] += MC_CV_weight[p];
				nHypo_4GeVcut[p] += MC_CV_4GeVcut_weight[p];
				if (ZZMass < 140.6 && ZZMass >= 105.6){
					nHypo_125p6MassWin[p] += MC_CV_weight[p];
					nHypo_4GeVcut_125p6MassWin[p] += MC_CV_4GeVcut_weight[p];
					if (ZZMass < 135 && ZZMass >= 115){
						nHypo_125p6RestrictedMassWin[p] += MC_CV_weight[p];
						nHypo_4GeVcut_125p6RestrictedMassWin[p] += MC_CV_4GeVcut_weight[p];
					};
				};
			};
			nHypo[myNumSamples+1] += MC_weight_xsec*MC_weight_PUWeight*MC_weight_powhegWeight*MC_weight_dataMC*MC_weight_HqT/sum_NGenTotal_unweighted;
			nHypoPredicted[lucky] += MC_weight;
			nSM[lucky] += MC_CV_weight[0];
//			if (GenZ1Mass > 4 && GenZ2Mass > 4){
//				nHypo_4GeVcut[myNumSamples + 1] += MC_weight_xsec_4GeVcut*MC_weight_PUWeight*MC_weight_powhegWeight*MC_weight_dataMC*MC_weight_HqT / sum_NGenTotal_4GeVcut_unweighted;
//				nHypoPredicted_4GeVcut[lucky] += MC_weight_4GeVcut;
//				nSM_4GeVcut[lucky] += MC_CV_4GeVcut_weight[0];
//			};
			nHypo_4GeVcut[myNumSamples + 1] += MC_weight_xsec_4GeVcut*MC_weight_PUWeight*MC_weight_powhegWeight*MC_weight_dataMC*MC_weight_HqT / sum_NGenTotal_4GeVcut_unweighted;
			nHypoPredicted_4GeVcut[lucky] += MC_weight_4GeVcut;
			nSM_4GeVcut[lucky] += MC_CV_4GeVcut_weight[0];
		};
		shuffledTree->Fill();
		for(int smp=lucky;smp<(kSignalSamples);smp++) naccevents[smp] -= 1;
		ctr++;
	};
	foutput->WriteTObject(shuffledTree);

	ctr=0;
	while(ctr<nTotal_ggZZ){
		int coin = rand() % naccevents[(kSignalSamples+kggZZSamples-1)];
		int lucky=-1;
		for(int smp=kSignalSamples;smp<(kSignalSamples+kggZZSamples);smp++){
			if(coin<naccevents[smp] && (nevents[smp]-nUsed[smp])!=0){lucky=smp;break;};
		};
		if(lucky<0) continue;

		myTree[lucky]->GetEntry(nUsed[lucky]);
		nUsed[lucky] += 1;

		MC_weight_xsec = MC_weight/MC_weight_noxsec;
		for(int p=0;p<myNumSamples;p++){
			MC_CV_weight[p] = MC_weight;
		};
		nHypo[myNumSamples] += MC_CV_weight[0];

		MC_weight_Kfactor = tgkf->Eval(GenHMass);
		MC_weight_Kfactor_up = tgkf_up->Eval(GenHMass);
		MC_weight_Kfactor_down = tgkf_down->Eval(GenHMass);
		MC_weight_alphaS_down = evaluateAlphaSShift(GenHMass, 1, -1, 125.6);
		MC_weight_alphaS_up = evaluateAlphaSShift(GenHMass, 1, 1, 125.6);
		MC_weight_Kfactor_up *= MC_weight_alphaS_up;
		MC_weight_Kfactor_down *= MC_weight_alphaS_down;

		if(ZZMass>=220) nGGZZ += MC_weight;
		if(ZZMass<140.6 && ZZMass>=105.6) nGGZZ160 += MC_weight;

		if(ZZMass>=220) nGGZZReWgt += MC_weight_QQBGGProper[0]*MC_weight_noxsec;
		if(ZZMass<140.6 && ZZMass>=105.6) nGGZZReWgt_in160 += MC_weight_QQBGGProper[0]*MC_weight_noxsec;
		if(ZZMass<140.6 && ZZMass>=105.6) nGGZZReWgt160 += MC_weight_QQBGGProper[0]*MC_weight_noxsec;
		if(ZZMass>=220) nQQZZReWgt += MC_weight_QQBGGProper[0]*MC_weight_noxsec;
		if(ZZMass<140.6 && ZZMass>=105.6) nQQZZReWgt_in160 += MC_weight_QQBGGProper[1]*MC_weight_noxsec;
		if(ZZMass<140.6 && ZZMass>=105.6) nQQZZ_DedicatedReWgt += MC_weight_QQBGGProper[1]*MC_weight_noxsec;

		shuffledggZZBkg->Fill();
		for(int smp=lucky;smp<(kSignalSamples+kggZZSamples);smp++) naccevents[smp] -= 1;

		ctr++;
	};
	foutput->WriteTObject(shuffledggZZBkg);

	ctr=0;
	while(ctr<nTotal_qqZZ){
		int coin = rand() % naccevents[(kSignalSamples+kggZZSamples+kqqZZSamples-1)];
		int lucky=-1;
		for(int smp=kSignalSamples+kggZZSamples;smp<(kSignalSamples+kggZZSamples+kqqZZSamples);smp++){
			if(coin<naccevents[smp] && (nevents[smp]-nUsed[smp])!=0){lucky=smp;break;};
		};
		if(lucky<0) continue;

		myTree[lucky]->GetEntry(nUsed[lucky]);
		nUsed[lucky] += 1;

		MC_weight_xsec = MC_weight/MC_weight_noxsec;
		for(int p=0;p<myNumSamples;p++){
			MC_CV_weight[p] = MC_weight;
		};
		nHypo[myNumSamples] += MC_CV_weight[0];

		MC_weight_Kfactor = tgkf->Eval(GenHMass);
		MC_weight_Kfactor_up = tgkf_up->Eval(GenHMass);
		MC_weight_Kfactor_down = tgkf_down->Eval(GenHMass);
		MC_weight_alphaS_down = evaluateAlphaSShift(GenHMass, 1, -1, 125.6);
		MC_weight_alphaS_up = evaluateAlphaSShift(GenHMass, 1, 1, 125.6);
		MC_weight_Kfactor_up *= MC_weight_alphaS_up;
		MC_weight_Kfactor_down *= MC_weight_alphaS_down;

		if(ZZMass>=220) nQQZZ += MC_weight;
		if(ZZMass<140.6 && ZZMass>=105.6) nQQZZ160 += MC_weight;

		if(ZZMass>=220) nGGZZReWgt += MC_weight_QQBGGProper[0]*MC_weight_noxsec;
		if(ZZMass<140.6 && ZZMass>=105.6) nGGZZReWgt_in160 += MC_weight_QQBGGProper[0]*MC_weight_noxsec;
		if(ZZMass<140.6 && ZZMass>=105.6) nGGZZReWgt160 += MC_weight_QQBGGProper[0]*MC_weight_noxsec;
		if(ZZMass>=220) nQQZZReWgt += MC_weight_QQBGGProper[0]*MC_weight_noxsec;
		if(ZZMass<140.6 && ZZMass>=105.6) nQQZZReWgt_in160 += MC_weight_QQBGGProper[1]*MC_weight_noxsec;
		if(ZZMass<140.6 && ZZMass>=105.6) nQQZZ_DedicatedReWgt += MC_weight_QQBGGProper[1]*MC_weight_noxsec;

		shuffledqqZZBkg->Fill();
		for(int smp=lucky;smp<(kSignalSamples+kggZZSamples+kqqZZSamples);smp++) naccevents[smp] -= 1;

		ctr++;
	};
	foutput->WriteTObject(shuffledqqZZBkg);

	ctr=0;
	while(ctr<nTotal_ZX){
		int coin = rand() % naccevents[(kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples-1)];
		int lucky=-1;
		for(int smp=kSignalSamples+kggZZSamples+kqqZZSamples;smp<(kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples);smp++){
			if(coin<naccevents[smp] && (nevents[smp]-nUsed[smp])!=0){lucky=smp;break;};
		};
		if(lucky<0) continue;

		myTree[lucky]->GetEntry(nUsed[lucky]);
		nUsed[lucky] += 1;

		nHypo[myNumSamples] += ZXfake_weightProper;
		if(ZZMass>=220) nZX += ZXfake_weightProper;
		shuffledZXBkg->Fill();
		for(int smp=lucky;smp<(kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples);smp++) naccevents[smp] -= 1;

		ctr++;
	};
	foutput->WriteTObject(shuffledZXBkg);

	ctr=0;
	while(ctr<nTotal_qqZZ_Dedicated){
		int coin = rand() % naccevents[(kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples+kqqZZSamples_Dedicated-1)];
		int lucky=-1;
		for(int smp=kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples;smp<(kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples+kqqZZSamples_Dedicated);smp++){
			if(coin<naccevents[smp] && (nevents[smp]-nUsed[smp])!=0){lucky=smp;break;};
		};
		if(lucky<0) continue;


		myTree[lucky]->GetEntry(nUsed[lucky]);
		nUsed[lucky] += 1;

		MC_weight_xsec = MC_weight/MC_weight_noxsec;
		for(int p=0;p<myNumSamples;p++){
			MC_CV_weight[p] = MC_weight;
		};
		nHypo[myNumSamples] += MC_CV_weight[0];

		MC_weight_Kfactor = tgkf->Eval(GenHMass);
		MC_weight_Kfactor_up = tgkf_up->Eval(GenHMass);
		MC_weight_Kfactor_down = tgkf_down->Eval(GenHMass);
		MC_weight_alphaS_down = evaluateAlphaSShift(GenHMass, 1, -1, 125.6);
		MC_weight_alphaS_up = evaluateAlphaSShift(GenHMass, 1, 1, 125.6);
		MC_weight_Kfactor_up *= MC_weight_alphaS_up;
		MC_weight_Kfactor_down *= MC_weight_alphaS_down;

		if(ZZMass>=220) nQQZZ_Dedicated += MC_weight;
		if(ZZMass<140.6 && ZZMass>=105.6) nQQZZ_Dedicated160 += MC_weight;

		if(ZZMass>=220) nGGZZReWgt += MC_weight_QQBGGProper[0]*MC_weight_noxsec;
		if(ZZMass<140.6 && ZZMass>=105.6) nGGZZReWgt_in160 += MC_weight_QQBGGProper[0]*MC_weight_noxsec;
		if(ZZMass<140.6 && ZZMass>=105.6) nGGZZReWgt160 += MC_weight_QQBGGProper[0]*MC_weight_noxsec;
		if(ZZMass>=220) nQQZZReWgt += MC_weight_QQBGGProper[0]*MC_weight_noxsec;
		if(ZZMass<140.6 && ZZMass>=105.6) nQQZZReWgt_in160 += MC_weight_QQBGGProper[1]*MC_weight_noxsec;
		if(ZZMass<140.6 && ZZMass>=105.6) nQQZZ_DedicatedReWgt += MC_weight_QQBGGProper[1]*MC_weight_noxsec;

		shuffledqqZZ_DedicatedBkg->Fill();
		for(int smp=lucky;smp<(kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples+kqqZZSamples_Dedicated);smp++) naccevents[smp] -= 1;

		ctr++;
	};
	foutput->WriteTObject(shuffledqqZZ_DedicatedBkg);


	delete shuffledqqZZ_DedicatedBkg;
	delete shuffledZXBkg;
	delete shuffledqqZZBkg;
	delete shuffledggZZBkg;
	delete shuffledTree;

	cout << "Sample SM integrals are " << endl;
	for(int f=0;f<kSignalSamples;f++) cout << f << ": " << nSM[f] << endl;
	cout << "Sample SM 4GeV cut integrals are " << endl;
	for(int f=0;f<kSignalSamples;f++) cout << f << ": " << nSM_4GeVcut[f] << endl;
	float sumnSM=0;
	float sumnSM_4GeVcut=0;
	for(int f=0;f<kSignalSamples;f++){
		sumnSM += nSM[f];
		sumnSM_4GeVcut += nSM_4GeVcut[f];
	};
	cout << "Sum SM (no cut, 4 GeV cut): " << sumnSM << '\t' << sumnSM_4GeVcut << endl;
	cout << "Average hypothesis integral is " << nHypo[myNumSamples+1] << endl;
	cout << "Average hypothesis integral with 4 GeV cut is " << nHypo_4GeVcut[myNumSamples+1] << endl;

	cout << "Hypothesis integrals: " << endl;
	cout << "Index\tCombined\tSample Prediction\tRaw Sample\t2e2mu/e,mu of Hypothesis" << endl;
	for(int p=0;p<13;p++) cout << p << "\t" << nHypo[p] << "\t" <<  (nHypoPredicted[p]*fGen_BRweight[p]) << "\t" << nHypoPredicted[p] << "\t" << fGen_BRweight[p] << endl;
	for(int p=13;p<14;p++) cout << p << "\t" << nHypo[p] << "\tn/a\tn/a\tn/a" << endl;
	for(int p=14;p<kSignalSamples+1;p++) cout << p << "\t" << nHypo[p] << "\t" << (nHypoPredicted[p-1]*fGen_BRweight[p-1]) << "\t" << nHypoPredicted[p-1] << "\t" << fGen_BRweight[p-1] << endl;
	for(int p=kSignalSamples+1;p<myNumSamples;p++) cout << p << "\t" << nHypo[p] << "\tn/a\tn/a\t" << fGen_BRweight[p-1] << endl;
	cout << "Hypothesis integrals with 4 GeV cut: " << endl;
	cout << "Index\tCombined\tSample Prediction\tRaw Sample\t2e2mu/e,mu of Hypothesis" << endl;
	for(int p=0;p<13;p++) cout << p << "\t" << nHypo_4GeVcut[p] << "\t" <<  (nHypoPredicted_4GeVcut[p]*fGen_4GeVcut_BRweight[p]) << "\t" << nHypoPredicted_4GeVcut[p] << "\t" << fGen_4GeVcut_BRweight[p] << endl;
	for(int p=13;p<14;p++) cout << p << "\t" << nHypo_4GeVcut[p] << "\tn/a\tn/a\tn/a" << endl;
	for(int p=14;p<kSignalSamples+1;p++) cout << p << "\t" << nHypo_4GeVcut[p] << "\t" << (nHypoPredicted_4GeVcut[p-1]*fGen_4GeVcut_BRweight[p-1]) << "\t" << nHypoPredicted_4GeVcut[p-1] << "\t" << fGen_4GeVcut_BRweight[p-1] << endl;
	for(int p=kSignalSamples+1;p<myNumSamples;p++) cout << p << "\t" << nHypo_4GeVcut[p] << "\tn/a\tn/a\t" << fGen_4GeVcut_BRweight[p-1] << endl;
	cout << "Total Un-weighted Bkg: " << nHypo[myNumSamples] << endl;

	cout << "Predicted Mass Window Signal Integrals: " << endl;
	cout << "Index\tNo cut\tWith cut\tNo cut restricted\tWith cut restricted" << endl;
	for(int p=0;p<myNumSamples;p++) cout << p << "\t" << nHypo_125p6MassWin[p] << "\t" << nHypo_4GeVcut_125p6MassWin[p] << "\t" << nHypo_125p6RestrictedMassWin[p] << "\t" << nHypo_4GeVcut_125p6RestrictedMassWin[p] << endl;

	cout << "nGGZZ Raw (mZZ>=220): " << nGGZZ/3.0 << endl;
	cout << "nGGZZ + qqZZZ (mZZ>=220): " << nGGZZReWgt << endl;
	cout << "nGGZZ Raw (105.6<=mZZ<140.6): " << nGGZZ160/3.0 << endl;
	cout << "nGGZZ + qqZZZ (105.6<=mZZ<140.6): " << nGGZZReWgt_in160 << endl;
	cout << "nGGZZ + qqZZ Dedicated (105.6<=mZZ<140.6): " << nGGZZReWgt160 << endl;
	cout << "nQQZZ Raw (mZZ>=220): " << nQQZZ << endl;
	cout << "nQQZZ + ggZZ (mZZ>=220): " << nQQZZReWgt << endl;
	cout << "nQQZZ Raw (105.6<=mZZ<140.6): " << nQQZZ160 << endl;
	cout << "nQQZZ Dedicated Raw (105.6<=mZZ<140.6): " << nQQZZ_Dedicated160 << endl;
	cout << "nQQZZ Dedicated + ggZZ (105.6<=mZZ<140.6): " << nQQZZ_DedicatedReWgt << endl;
	cout << "nQQZZ + ggZZ (105.6<=mZZ<140.6): " << nQQZZReWgt_in160 << endl;
	cout << "nZX (mZZ>=220): " << nZX << endl;

	foutput->Close();
	for(int smp=0;smp<(kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples+kqqZZSamples_Dedicated);smp++){
		delete myTree[smp];
		finput[smp]->Close();
	};

	finput_KDFactor->Close();
};
