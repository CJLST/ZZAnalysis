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


void countShuffleSignal_125p6(int flavor, int erg_tev){
	srand( time(0) );

	string INPUT_NAME = hzz4lprefix;
	INPUT_NAME = INPUT_NAME + "H125p6_ShuffledSignalBkg";
	TString TREE_NAME = "SelectedTree";

	char erg_dir[1000];
	sprintf(erg_dir,"LHC_%iTeV",erg_tev);
	char cerg[1000];
	sprintf(cerg,"%iTeV",erg_tev);
	string cinput_common = user_dir + erg_dir;
//	string cinput_common = user_dir_newProduction + cerg;
	string coutput_common = cinput_common;

	const int numSamples=44;
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
//	float MC_CVwrt2mu2e_weight[numSamples];

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

	TFile* finput[kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples+kqqZZSamples_Dedicated];
	TTree* myTree[kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples+kqqZZSamples_Dedicated];

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

	for (int smp = 0; smp < kSignalSamples; smp++){
		string cinput = cinput_common;
		cinput = cinput + "/" + user_folder[flavor] + "/";
		if (smp < kSignalSamples) cinput = cinput + sampleName_Sig[smp] + "_Reprocessed.root";

//		cout << cinput << '\t' << smp << endl;

		finput[smp] = new TFile(cinput.c_str(), "read");

		if (smp < kSignalSamples){
			string cinputGen = cinput_common + "/GenSignal/" + sampleName_Sig[smp] + "_GenLevel_ReWeighed.root";
//			cout << cinputGen << endl;
			TFile* ftemp = new TFile(cinputGen.c_str(), "read");
			TTree* ttemp = (TTree*)ftemp->Get("GenTree");
			float MC_gen_spin0[numSamples] = { 0 };
			int genFlavor = 0;
			float genMZ1 = 0;
			float genMZ2 = 0;
			int nGenEntries = ttemp->GetEntries();
			ttemp->SetBranchAddress("genFinalState", &genFlavor);
			ttemp->SetBranchAddress("GenZ1Mass", &genMZ1);
			ttemp->SetBranchAddress("GenZ2Mass", &genMZ2);
			ttemp->SetBranchAddress("MC_weight_spin0", MC_gen_spin0);
			TH2F* htrue = (TH2F*)ftemp->Get("hCounters_spin0_RW");
			TH2F* htrue_4GeVcut = (TH2F*)ftemp->Get("hCounters_spin0_4GeVcut_RW");
//			cout << htrue->GetName() << '\t' << htrue_4GeVcut->GetName() << endl;
			TH1F* hreco = (TH1F*)finput[smp]->Get("Counters");
			double nTrueProcessed = hreco->GetBinContent(1);
			for (int myf = 0; myf < nFinalStates; myf++){
				NGenTotal_unweighted_perFlavor[smp][myf] = htrue->GetBinContent(myf + 1, 1);
				NGenTotal_4GeVcut_unweighted_perFlavor[smp][myf] = htrue_4GeVcut->GetBinContent(myf + 1, 1);
			};
			for (int myf = 0; myf < nFinalStates; myf++){
				NGenTotal_unweighted[smp] += NGenTotal_unweighted_perFlavor[smp][myf];
				NGenTotal_4GeVcut_unweighted[smp] += NGenTotal_4GeVcut_unweighted_perFlavor[smp][myf];
			};
			for (int myf = 0; myf < nFinalStates; myf++){
				for (int p = 0; p < numSamples; p++){
					NGenTotal_weighted_perFlavor[smp][p][myf] = htrue->GetBinContent(myf + 1, p + 2);
					NGenTotal_4GeVcut_weighted_perFlavor[smp][p][myf] = htrue_4GeVcut->GetBinContent(myf + 1, p + 2);
					if (NGenTotal_weighted_perFlavor[smp][p][myf] != NGenTotal_weighted_perFlavor[smp][p][myf]){
						NGenTotal_weighted_perFlavor[smp][p][myf] = 0;
						NGenTotal_4GeVcut_weighted_perFlavor[smp][p][myf] = 0;
						for (int ev = 0; ev<nGenEntries; ev++){
							ttemp->GetEntry(ev);
							if (genFlavor == myf && MC_gen_spin0[p] == MC_gen_spin0[p]){
								NGenTotal_weighted_perFlavor[smp][p][myf] += MC_gen_spin0[p];
								if (genMZ1>4 && genMZ2 > 4) NGenTotal_4GeVcut_weighted_perFlavor[smp][p][myf] += MC_gen_spin0[p];
							};
						};
					};
				};
			};
			for (int myf = 0; myf < nFinalStates; myf++){
				for (int p = 0; p < numSamples; p++){
					NGenTotal_weighted[smp][p] += NGenTotal_weighted_perFlavor[smp][p][myf];
					NGenTotal_4GeVcut_weighted[smp][p] += NGenTotal_4GeVcut_weighted_perFlavor[smp][p][myf];
					//					NGenTotal_weighted[smp][p] += NGenTotal_weighted_perFlavor[smp][p][myf]/NGenTotal_weighted_perFlavor[smp][p][2]*NGenTotal_unweighted_perFlavor[smp][2];
				};
			};

			for (int p = 0; p < myNumSamples; p++){
				NGenTotal_weighted[smp][p] *= nTrueProcessed / NGenTotal_unweighted[smp];
				NGenTotal_4GeVcut_weighted[smp][p] *= nTrueProcessed / NGenTotal_unweighted[smp];
			};
			NGenTotal_unweighted[smp] *= nTrueProcessed / NGenTotal_unweighted[smp];
			NGenTotal_4GeVcut_unweighted[smp] *= nTrueProcessed / NGenTotal_unweighted[smp];
			for (int myf = 0; myf < nFinalStates; myf++){
				for (int p = 0; p < myNumSamples; p++){
					NGenTotal_weighted_perFlavor[smp][p][myf] *= nTrueProcessed / NGenTotal_unweighted[smp];
					NGenTotal_4GeVcut_weighted_perFlavor[smp][p][myf] *= nTrueProcessed / NGenTotal_unweighted[smp];
				};
				NGenTotal_unweighted_perFlavor[smp][myf] *= nTrueProcessed / NGenTotal_unweighted[smp];
				NGenTotal_4GeVcut_unweighted_perFlavor[smp][myf] *= nTrueProcessed / NGenTotal_unweighted[smp];
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
		finput[smp]->Close();
	};
	for(int p=0;p<numSamples;p++){
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
//		cout << "BR 2e2mu in hypothesis " << p << ": " << fGen_BRweight_perFlavor[p][2] << endl;
		fGen_BRweight[p] = fGen_BRweight_perFlavor[0][2]/fGen_BRweight_perFlavor[p][2];
//		cout << "BR 2e2mu (mZ>4 GeV) in hypothesis " << p << ": " << fGen_4GeVcut_BRweight_perFlavor[p][2] << endl;
		fGen_4GeVcut_BRweight[p] = fGen_4GeVcut_BRweight_perFlavor[0][2]/fGen_4GeVcut_BRweight_perFlavor[p][2];
	};


	string cinput_tree = cinput_common + "/" + user_folder[flavor] + "/" + INPUT_NAME + ".root";
	TFile* finput_tree = new TFile(cinput_tree.c_str(),"read");

	TTree* shuffledTree = (TTree*) finput_tree->Get(TREE_NAME);
	shuffledTree->SetBranchAddress("EventSample",&EventSample);
	shuffledTree->SetBranchAddress("kNumSamples",&myNumSamples);
	shuffledTree->SetBranchAddress("MC_CV_weight",MC_CV_weight);
/*//	shuffledTree->SetBranchAddress("MC_CVwrt2mu2e_weight",MC_CVwrt2mu2e_weight,"MC_CVwrt2mu2e_weight[kNumSamples]/F");
	shuffledTree->SetBranchAddress("MC_weight_spin0",MC_weight_spin0,"MC_weight_spin0[kNumSamples]/F");
*/	shuffledTree->SetBranchAddress("MC_weight",&MC_weight);
//	shuffledTree->SetBranchAddress("MC_weight_xsec",&MC_weight_xsec);
	shuffledTree->SetBranchAddress("MC_CV_4GeVcut_weight",MC_CV_4GeVcut_weight);
//	shuffledTree->SetBranchAddress("MC_weight_4GeVcut_spin0",MC_weight_4GeVcut_spin0,"MC_weight_4GeVcut_spin0[kNumSamples]/F");
	shuffledTree->SetBranchAddress("MC_weight_4GeVcut",&MC_weight_4GeVcut);
/*	shuffledTree->SetBranchAddress("MC_weight_xsec_4GeVcut",&MC_weight_xsec_4GeVcut);
	shuffledTree->SetBranchAddress("MC_weight_Kfactor",&MC_weight_Kfactor);
	shuffledTree->SetBranchAddress("MC_weight_Kfactor_down",&MC_weight_Kfactor_down);
	shuffledTree->SetBranchAddress("MC_weight_Kfactor_up",&MC_weight_Kfactor_up);
*///	shuffledTree->SetBranchAddress("GenHMass", &GenHMass);
//	shuffledTree->SetBranchAddress("GenZ1Mass", &GenZ1Mass);
//	shuffledTree->SetBranchAddress("GenZ2Mass", &GenZ2Mass);
	shuffledTree->SetBranchAddress("ZZMass", &ZZMass);
/*	shuffledTree->SetBranchAddress("Z1Mass", &Z1Mass);
	shuffledTree->SetBranchAddress("Z2Mass", &Z2Mass);
	shuffledTree->SetBranchAddress("helcosthetaZ1", &myhelcosthetaZ1);
	shuffledTree->SetBranchAddress("helcosthetaZ2", &myhelcosthetaZ2);
	shuffledTree->SetBranchAddress("helphi", &myhelphi);
	shuffledTree->SetBranchAddress("costhetastar", &mycosthetastar);
	shuffledTree->SetBranchAddress("phistarZ1", &myphistarZ1);
*/	shuffledTree->SetBranchAddress("Z1ids", &Z1ids);
	shuffledTree->SetBranchAddress("Z2ids", &Z2ids);
/*	shuffledTree->SetBranchAddress("D_g1_vs_g2_phi0",&D_g2);
	shuffledTree->SetBranchAddress("D_g1_vs_g4_phi0",&D_g4);
	shuffledTree->SetBranchAddress("D_g2int_phi0",&D_g2int);
	shuffledTree->SetBranchAddress("D_g4int_phi0",&D_g4int);
	shuffledTree->SetBranchAddress("D_g2int_phi90",&D_g2int_perp);
	shuffledTree->SetBranchAddress("D_g4int_phi90",&D_g4int_perp);
	shuffledTree->SetBranchAddress("D_g1Q2_phi0",&D_g1q2);
	shuffledTree->SetBranchAddress("D_g1Q2int_phi0",&D_g1q2int);
	shuffledTree->SetBranchAddress("D_ZG",&D_ZG);
	shuffledTree->SetBranchAddress("D_GG",&D_GG);
	shuffledTree->SetBranchAddress("D_ZG_PS",&D_ZG_PS);
	shuffledTree->SetBranchAddress("D_GG_PS",&D_GG_PS);
	shuffledTree->SetBranchAddress("D_ZGint",&D_ZGint);
	shuffledTree->SetBranchAddress("D_GGint",&D_GGint);
	shuffledTree->SetBranchAddress("D_ZG_PSint",&D_ZG_PSint);
	shuffledTree->SetBranchAddress("D_GG_PSint",&D_GG_PSint);
	shuffledTree->SetBranchAddress("D_bkg_kin",&D_bkg_kin);
	shuffledTree->SetBranchAddress("D_bkg",&D_bkg);
	shuffledTree->SetBranchAddress("D_bkg_ScaleUp",&D_bkg_ScaleUp);
	shuffledTree->SetBranchAddress("D_bkg_ScaleDown",&D_bkg_ScaleDown);
	shuffledTree->SetBranchAddress("D_bkg_ResUp",&D_bkg_ResUp);
	shuffledTree->SetBranchAddress("D_bkg_ResDown",&D_bkg_ResDown);
	shuffledTree->SetBranchAddress("D_ggZZ",&D_ggZZ);
	shuffledTree->SetBranchAddress("D_Gamma_r1",&D_Gamma_r1);
	shuffledTree->SetBranchAddress("D_Gamma_r5",&D_Gamma_r5);
	shuffledTree->SetBranchAddress("D_Gamma_r10",&D_Gamma_r10);
	shuffledTree->SetBranchAddress("D_Gamma_r15",&D_Gamma_r15);
	shuffledTree->SetBranchAddress("D_Gamma_r20",&D_Gamma_r20);
	shuffledTree->SetBranchAddress("D_Gamma_r25",&D_Gamma_r25);
	shuffledTree->SetBranchAddress("D_Gamma_gg_r1",&D_Gamma_gg_r1);
	shuffledTree->SetBranchAddress("D_Gamma_gg_r5",&D_Gamma_gg_r5);
	shuffledTree->SetBranchAddress("D_Gamma_gg_r10",&D_Gamma_gg_r10);
	shuffledTree->SetBranchAddress("D_Gamma_gg_r15",&D_Gamma_gg_r15);
	shuffledTree->SetBranchAddress("D_Gamma_gg_r20",&D_Gamma_gg_r20);
	shuffledTree->SetBranchAddress("D_Gamma_gg_r25",&D_Gamma_gg_r25);
	shuffledTree->SetBranchAddress("D_Gamma",&D_Gamma);
	shuffledTree->SetBranchAddress("D_Gamma_int",&D_Gamma_int);
*/
	cout << "Starting to count shuffle signal trees... Total N_events is " << sum_NGenTotal_unweighted << "." << endl;

	double nHypoPredicted[kSignalSamples][3] = { { 0 } };
	double nHypoPredicted_125p6MassWin[kSignalSamples][3] = { { 0 } };
	double nHypoPredicted_125p6RestrictedMassWin[kSignalSamples][3] = { { 0 } };

	double nHypo[numSamples][3] = { { 0 } };
	double nHypo_4GeVcut[numSamples][3] = { { 0 } };

	double nHypo_125p6MassWin[numSamples][3] = { { 0 } };
	double nHypo_4GeVcut_125p6MassWin[numSamples][3] = { { 0 } };

	double nHypo_125p6RestrictedMassWin[numSamples][3] = { { 0 } };
	double nHypo_4GeVcut_125p6RestrictedMassWin[numSamples][3] = { { 0 } };

	int ctr=0;
	for (int ev = 0; ev < shuffledTree->GetEntries(); ev++){
		shuffledTree->GetEntry(ev);

		for(int p=0;p<myNumSamples;p++){
			nHypo[p][0] += MC_CV_weight[p];
			nHypo_4GeVcut[p][0] += MC_CV_4GeVcut_weight[p];
			nHypo[p][1] += pow(MC_CV_weight[p],2);
			nHypo_4GeVcut[p][1] += pow(MC_CV_4GeVcut_weight[p],2);
			nHypo[p][2] += 1;
			if (MC_CV_4GeVcut_weight[p]>0) nHypo_4GeVcut[p][2] += 1;
			if (ZZMass < 140.6 && ZZMass >= 105.6){
				nHypo_125p6MassWin[p][0] += MC_CV_weight[p];
				nHypo_4GeVcut_125p6MassWin[p][0] += MC_CV_4GeVcut_weight[p];
				nHypo_125p6MassWin[p][1] += pow(MC_CV_weight[p],2);
				nHypo_4GeVcut_125p6MassWin[p][1] += pow(MC_CV_4GeVcut_weight[p],2);
				nHypo_125p6MassWin[p][2] += 1;
				if (MC_CV_4GeVcut_weight[p]>0) nHypo_4GeVcut_125p6MassWin[p][2] += 1;
				if (ZZMass < 135 && ZZMass >= 115){
					nHypo_125p6RestrictedMassWin[p][0] += MC_CV_weight[p];
					nHypo_4GeVcut_125p6RestrictedMassWin[p][0] += MC_CV_4GeVcut_weight[p];
					nHypo_125p6RestrictedMassWin[p][1] += pow(MC_CV_weight[p],2);
					nHypo_4GeVcut_125p6RestrictedMassWin[p][1] += pow(MC_CV_4GeVcut_weight[p],2);
					nHypo_125p6RestrictedMassWin[p][2] += 1;
					if (MC_CV_4GeVcut_weight[p]>0) nHypo_4GeVcut_125p6RestrictedMassWin[p][2] += 1;
				};
			};
		};
		nHypoPredicted[EventSample][0] += MC_weight;
		nHypoPredicted[EventSample][1] += pow(MC_weight,2);
		nHypoPredicted[EventSample][2] += 1;

		if (ZZMass < 140.6 && ZZMass >= 105.6){
			nHypoPredicted_125p6MassWin[EventSample][0] += MC_weight;
			nHypoPredicted_125p6MassWin[EventSample][1] += pow(MC_weight,2);
			if (MC_CV_4GeVcut_weight[0]>0) nHypoPredicted_125p6MassWin[EventSample][2] += 1;
			if (ZZMass < 135 && ZZMass >= 115){
				nHypoPredicted_125p6RestrictedMassWin[EventSample][0] += MC_weight;
				nHypoPredicted_125p6RestrictedMassWin[EventSample][1] += pow(MC_weight,2);
				if (MC_CV_4GeVcut_weight[0]>0) nHypoPredicted_125p6RestrictedMassWin[EventSample][2] += 1;
			};
		};
	};

	delete shuffledTree;
	finput_tree->Close();

	cout << "Hypothesis integrals: " << endl;
	cout << "Index\tCombined\tDedicated Sample" << endl;
	for(int p=0;p<13;p++) cout << p << "\t" << nHypo[p][0] << "\t" <<  nHypoPredicted[p][0]*fGen_BRweight[p] << endl;
	for(int p=13;p<14;p++) cout << p << "\t" << nHypo[p][0] << "\tn/a" << endl;
	for(int p=14;p<kSignalSamples+1;p++) cout << p << "\t" << nHypo[p][0] << "\t" << nHypoPredicted[p-1][0]*fGen_BRweight[p-1] << endl;
	for(int p=kSignalSamples+1;p<myNumSamples;p++) cout << p << "\t" << nHypo[p][0] << "\tn/a" << endl;
	cout << "Hypothesis integrals with 4 GeV cut: " << endl;
	cout << "Index\tCombined" << endl;
	for(int p=0;p<myNumSamples;p++) cout << p << "\t" << nHypo_4GeVcut[p][0] << endl;
	cout << "Predicted Mass Window Signal Integrals (Combined): " << endl;
	cout << "Index\tNo cut\tWith cut\tNo cut restricted\tWith cut restricted" << endl;
	for(int p=0;p<myNumSamples;p++) cout << p << "\t" << nHypo_125p6MassWin[p][0] << "\t" << nHypo_4GeVcut_125p6MassWin[p][0] << "\t" << nHypo_125p6RestrictedMassWin[p][0] << "\t" << nHypo_4GeVcut_125p6RestrictedMassWin[p][0] << endl;

	cout << "Hypothesis integral fractional uncertainties: " << endl;
	cout << "Index\tCombined\tDedicated Sample" << endl;
	for(int p=0;p<13;p++) cout << p << "\t" << sqrt(nHypo[p][1]*sum_NGenTotal_unweighted/nHypo[p][2])/nHypo[p][0] << "\t" <<  sqrt(nHypoPredicted[p][1]*fGen_BRweight[p]*NGenTotal_unweighted[p]/nHypoPredicted[p][2])/(nHypoPredicted[p][0]*fGen_BRweight[p]) << endl;
	for(int p=13;p<14;p++) cout << p << "\t" << sqrt(nHypo[p][1]*sum_NGenTotal_unweighted/nHypo[p][2])/nHypo[p][0] << "\tn/a" << endl;
	for(int p=14;p<kSignalSamples+1;p++) cout << p << "\t" << sqrt(nHypo[p][1]*sum_NGenTotal_unweighted/nHypo[p][2])/nHypo[p][0] << "\t" << sqrt(nHypoPredicted[p-1][1]*fGen_BRweight[p-1]*NGenTotal_unweighted[p-1]/nHypoPredicted[p-1][2])/(nHypoPredicted[p-1][0]*fGen_BRweight[p-1]) << endl;
	for(int p=kSignalSamples+1;p<myNumSamples;p++) cout << p << "\t" << sqrt(nHypo[p][1]*sum_NGenTotal_unweighted/nHypo[p][2])/nHypo[p][0] << "\tn/a" << endl;
	cout << "Hypothesis integrals with 4 GeV cut: " << endl;
	cout << "Index\tCombined" << endl;
	for(int p=0;p<myNumSamples;p++) cout << p << "\t" << sqrt(nHypo_4GeVcut[p][1]*sum_NGenTotal_4GeVcut_unweighted/nHypo_4GeVcut[p][2])/nHypo_4GeVcut[p][0] << endl;
	cout << "Predicted Mass Window Signal Integral (Combined) Fractional Uncertainties: " << endl;
	cout << "Index\tNo cut\tWith cut\tNo cut restricted\tWith cut restricted" << endl;
	for(int p=0;p<myNumSamples;p++) cout << p << "\t" << sqrt(nHypo_125p6MassWin[p][1]*sum_NGenTotal_unweighted/nHypo_125p6MassWin[p][2])/nHypo_125p6MassWin[p][0] << "\t" << sqrt(nHypo_4GeVcut_125p6MassWin[p][1]*sum_NGenTotal_4GeVcut_unweighted/nHypo_4GeVcut_125p6MassWin[p][2])/nHypo_4GeVcut_125p6MassWin[p][0] << "\t" << sqrt(nHypo_125p6RestrictedMassWin[p][1]*sum_NGenTotal_unweighted/nHypo_125p6RestrictedMassWin[p][2])/nHypo_125p6RestrictedMassWin[p][0] << "\t" << sqrt(nHypo_4GeVcut_125p6RestrictedMassWin[p][1]*sum_NGenTotal_4GeVcut_unweighted/nHypo_4GeVcut_125p6RestrictedMassWin[p][2])/nHypo_4GeVcut_125p6RestrictedMassWin[p][0] << endl;

	cout << "Hypothesis N_events: " << endl;
	cout << "Index\tCombined\tDedicated Sample" << endl;
	for(int p=0;p<13;p++) cout << p << "\t" << nHypo[p][2] << "\t" <<  nHypoPredicted[p][2] << endl;
	for(int p=13;p<14;p++) cout << p << "\t" << nHypo[p][2] << "\tn/a" << endl;
	for(int p=14;p<kSignalSamples+1;p++) cout << p << "\t" << nHypo[p][2] << "\t" << nHypoPredicted[p-1][2] << endl;
	for(int p=kSignalSamples+1;p<myNumSamples;p++) cout << p << "\t" << nHypo[p][2] << "\tn/a" << endl;
	cout << "Hypothesis N_events with 4 GeV cut: " << endl;
	cout << "Index\tCombined" << endl;
	for(int p=0;p<myNumSamples;p++) cout << p << "\t" << nHypo_4GeVcut[p][2] << endl;
	cout << "Predicted Mass Window Signal N_events (Combined): " << endl;
	cout << "Index\tNo cut\tWith cut\tNo cut restricted\tWith cut restricted" << endl;
	for(int p=0;p<myNumSamples;p++) cout << p << "\t" << nHypo_125p6MassWin[p][2] << "\t" << nHypo_4GeVcut_125p6MassWin[p][2] << "\t" << nHypo_125p6RestrictedMassWin[p][2] << "\t" << nHypo_4GeVcut_125p6RestrictedMassWin[p][2] << endl;
};
