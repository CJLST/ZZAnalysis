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

const int kSignalSamples=21;
const int kggZZSamples=2;
const int kqqZZSamples=6;
const int kZXSamples=2;
char* sampleName_Sig[kSignalSamples]={
	"HZZ4lTree_powheg15jhuGenV3-0PMH125.6",
	"HZZ4lTree_powheg15jhuGenV3-0L1H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0L1f01ph0H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0L1f05ph0H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0L1f05ph180H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0MH125.6",
	"HZZ4lTree_powheg15jhuGenV3-0Mf01ph90H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0Mf05ph180H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0Mf05ph90H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf01ph0H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf01ph0Mf01ph0H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf01ph0Mf01ph90H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf01ph90H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf033ph0Mf033ph0H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf033ph0Mf033ph90H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf05ph0H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf05ph0Mf05ph0H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf05ph0Mf05ph90H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf05ph180H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf05ph180Mf05ph0H125.6",
	"HZZ4lTree_powheg15jhuGenV3-0PHf05ph90H125.6"
};
char* sampleName_ggZZ[kggZZSamples]={
	"HZZ4lTree_ggZZ4l",
	"HZZ4lTree_ggZZ2l2l"
};
char* sampleName_qqZZ[kqqZZSamples]={
	"HZZ4lTree_ZZTo4mu",
	"HZZ4lTree_ZZTo4e",
	"HZZ4lTree_ZZTo2e2mu",
	"HZZ4lTree_ZZTo2mu2tau",
	"HZZ4lTree_ZZTo2e2tau",
	"HZZ4lTree_ZZTo4tau"
};
char* sampleName_ZX[kZXSamples]={
	"HZZ4lTree_DYJetsToLLTuneZ2M50-B",
	"HZZ4lTree_DYJetsToLLTuneZ2M50-NoB"
};

void shuffleSignalBkg_125p6(int flavor, int erg_tev){
	srand( time(0) );

	string OUTPUT_NAME = hzz4lprefix;
	OUTPUT_NAME = OUTPUT_NAME + "H125p6_ShuffledSignalBkg";
	TString TREE_NAME = "SelectedTree";

	char erg_dir[1000];
	sprintf(erg_dir,"LHC_%iTeV",erg_tev);
	string cinput_common = user_dir + erg_dir;
	string coutput_common = cinput_common;

	const int numSamples=25;
	int myNumSamples=numSamples;
	float MC_weight_spin0[numSamples]={0};

	int genFinalState=0;
	float ZZMass=0;
	float Z1Mass=0;
	float Z2Mass=0;

	float MC_weight=1;
	float MC_weight_noxsec=1;
	float MC_weight_xsec=0;
	float MC_weight_norm=1;
	float MC_weight_PUWeight=1;
	float MC_weight_powhegWeight=1;
	float MC_weight_dataMC=1;
	float MC_weight_HqT=1;

	float MC_weight_Kfactor=1;
	float MC_weight_down=1;
	float MC_weight_up=1;

	float D_ggZZ=-99;
	float D_bkg=-99;
	float D_bkg_kin=-99;
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

	float MC_CV_weight[numSamples];
	float MC_CVwrt2mu2e_weight[numSamples];

	TFile* finput[kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples];
	TTree* myTree[kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples];

	int nevents[kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples]={0};
	int naccevents[kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples]={0};

	int nTotal_Sig=0;
	int nTotal_ggZZ=0;
	int nTotal_qqZZ=0;
	int nTotal_ZX=0;

	float NGenTotal_weighted[kSignalSamples][numSamples];
	float NGenTotal_unweighted[kSignalSamples]={0};
	float NGenTotal_weighted_perFlavor[kSignalSamples][numSamples][nFinalStates];
	float NGenTotal_unweighted_perFlavor[kSignalSamples][nFinalStates];
	float sum_NGenTotal_unweighted=0;

	for(int smp=0;smp<kSignalSamples;smp++){
		for(int hypo=0;hypo<numSamples;hypo++){
			NGenTotal_weighted[smp][hypo]=0;
			for(int fl=0;fl<nFinalStates;fl++) NGenTotal_weighted_perFlavor[smp][hypo][fl]=0;
		};
		for(int fl=0;fl<nFinalStates;fl++) NGenTotal_unweighted_perFlavor[smp][fl]=0;
	};

	for(int smp=0;smp<(kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples);smp++){
		string cinput = cinput_common;
		cinput = cinput + "/" + user_folder[flavor] + "/";
		if(smp<kSignalSamples) cinput = cinput + sampleName_Sig[smp] + "_Reprocessed.root";
		else if(smp>=kSignalSamples && smp<kggZZSamples+kSignalSamples) cinput = cinput + sampleName_ggZZ[smp-kSignalSamples] + "_Reprocessed.root";
		else if(smp>=kggZZSamples+kSignalSamples && smp<kqqZZSamples+kggZZSamples+kSignalSamples) cinput = cinput + sampleName_qqZZ[smp-kSignalSamples-kggZZSamples] + "_Reprocessed.root";
		else if(smp>=kqqZZSamples+kggZZSamples+kSignalSamples) cinput = cinput + sampleName_ZX[smp-kSignalSamples-kggZZSamples-kqqZZSamples] + "_Reprocessed.root";

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
			int nGenEntries = ttemp->GetEntries();
			ttemp->SetBranchAddress("MC_weight_spin0",MC_gen_spin0);
			ttemp->SetBranchAddress("genFinalState",&genFlavor);
			TH2F* htrue = (TH2F*) ftemp->Get("hCounters_spin0_RW");
			TH1F* hreco = (TH1F*) finput[smp]->Get("Counters");
			float nTrueProcessed = hreco->GetBinContent(1);
			for(int myf=0;myf<nFinalStates;myf++){
				for(int p=0;p<numSamples;p++){
					NGenTotal_weighted_perFlavor[smp][p][myf] = htrue->GetBinContent(myf+1,p+2);
					if(NGenTotal_weighted_perFlavor[smp][p][myf]!=NGenTotal_weighted_perFlavor[smp][p][myf]){
						NGenTotal_weighted_perFlavor[smp][p][myf]=0;
						for(int ev=0;ev<nGenEntries;ev++){
							ttemp->GetEntry(ev);
							if(genFlavor==myf && MC_gen_spin0[p]==MC_gen_spin0[p]) NGenTotal_weighted_perFlavor[smp][p][myf] += MC_gen_spin0[p];
						};
					};
					NGenTotal_weighted[smp][p] += NGenTotal_weighted_perFlavor[smp][p][myf];
				};
				NGenTotal_unweighted_perFlavor[smp][myf] = htrue->GetBinContent(myf+1,1);
				NGenTotal_unweighted[smp] += NGenTotal_unweighted_perFlavor[smp][myf];
			};

			for(int p=0;p<myNumSamples;p++) NGenTotal_weighted[smp][p] *= nTrueProcessed/NGenTotal_unweighted[smp];
			NGenTotal_unweighted[smp] *= nTrueProcessed/NGenTotal_unweighted[smp];
			for(int myf=0;myf<nFinalStates;myf++){
				for(int p=0;p<myNumSamples;p++){
					NGenTotal_weighted_perFlavor[smp][p][myf] *= nTrueProcessed/NGenTotal_unweighted[smp];
				};
				NGenTotal_unweighted_perFlavor[smp][myf] *= nTrueProcessed/NGenTotal_unweighted[smp];
			};

			sum_NGenTotal_unweighted += NGenTotal_unweighted[smp];
			
			delete hreco;
			delete htrue;
			delete ttemp;

			ftemp->Close();

			nTotal_Sig += nevents[smp];
		}
		else if(smp>=kSignalSamples && smp<kggZZSamples+kSignalSamples){
			nTotal_ggZZ += nevents[smp];
		}
		else if(smp>=kggZZSamples+kSignalSamples && smp<kqqZZSamples+kggZZSamples+kSignalSamples){
			nTotal_qqZZ += nevents[smp];
		}
		else if(smp>=kqqZZSamples+kggZZSamples+kSignalSamples){
			nTotal_ZX += nevents[smp];
		};

		if(smp!=0 && smp!=kSignalSamples && smp!=kggZZSamples+kSignalSamples && smp!=kqqZZSamples+kggZZSamples+kSignalSamples) naccevents[smp] += naccevents[smp-1];
		cout << nevents[smp] << endl;
		cout << naccevents[smp] << endl;

		myTree[smp]->SetBranchAddress("MC_weight_spin0",MC_weight_spin0);
		myTree[smp]->SetBranchAddress("genFinalState", &genFinalState);
		myTree[smp]->SetBranchAddress("MC_weight",&MC_weight);
		myTree[smp]->SetBranchAddress("MC_weight_noxsec",&MC_weight_noxsec);
		myTree[smp]->SetBranchAddress("MC_weight_xsec",&MC_weight_xsec);
		myTree[smp]->SetBranchAddress("MC_weight_PUWeight",&MC_weight_PUWeight);
		myTree[smp]->SetBranchAddress("MC_weight_powhegWeight",&MC_weight_powhegWeight);
		myTree[smp]->SetBranchAddress("MC_weight_dataMC",&MC_weight_dataMC);
		myTree[smp]->SetBranchAddress("MC_weight_HqT",&MC_weight_HqT);

		myTree[smp]->SetBranchAddress("MC_weight_Kfactor",&MC_weight_Kfactor);
		myTree[smp]->SetBranchAddress("MC_weight_down",&MC_weight_down);
		myTree[smp]->SetBranchAddress("MC_weight_up",&MC_weight_up);

		myTree[smp]->SetBranchAddress("ZZMass", &ZZMass);
		myTree[smp]->SetBranchAddress("Z1Mass", &Z1Mass);
		myTree[smp]->SetBranchAddress("Z2Mass", &Z2Mass);

		myTree[smp]->SetBranchAddress("D_g1_vs_g2_phi0",&D_g2);
		myTree[smp]->SetBranchAddress("D_g1_vs_g4_phi0",&D_g4);
		myTree[smp]->SetBranchAddress("D_g2int_phi0",&D_g2int);
		myTree[smp]->SetBranchAddress("D_g4int_phi0",&D_g4int);
		myTree[smp]->SetBranchAddress("D_g2int_phi90",&D_g2int_perp);
		myTree[smp]->SetBranchAddress("D_g4int_phi90",&D_g4int_perp);
		myTree[smp]->SetBranchAddress("D_g1Q2_phi0",&D_g1q2);
		myTree[smp]->SetBranchAddress("D_g1Q2int_phi0",&D_g1q2int);
		myTree[smp]->SetBranchAddress("D_bkg",&D_bkg);
		myTree[smp]->SetBranchAddress("D_bkg_kin",&D_bkg_kin);

		myTree[smp]->SetBranchAddress("D_ggZZ",&D_ggZZ);
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
	};

	string coutput = coutput_common + "/" + user_folder[flavor] + "/" + OUTPUT_NAME + ".root";
	TFile* foutput = new TFile(coutput.c_str(),"recreate");

	TTree* shuffledTree = new TTree(TREE_NAME,TREE_NAME);
	shuffledTree->SetAutoSave(3000000000);
	shuffledTree->Branch("kNumSamples",&myNumSamples);
	shuffledTree->Branch("MC_CV_weight",MC_CV_weight,"MC_CV_weight[kNumSamples]/F");
	shuffledTree->Branch("MC_CVwrt2mu2e_weight",MC_CVwrt2mu2e_weight,"MC_CVwrt2mu2e_weight[kNumSamples]/F");
	shuffledTree->Branch("MC_weight_spin0",MC_weight_spin0,"MC_weight_spin0[kNumSamples]/F");
	shuffledTree->Branch("MC_weight",&MC_weight);
	shuffledTree->Branch("MC_weight_down",&MC_weight_down);
	shuffledTree->Branch("MC_weight_up",&MC_weight_up);
//	shuffledTree->Branch("MC_weight_noxsec",&MC_weight_noxsec);
	shuffledTree->Branch("MC_weight_xsec",&MC_weight_xsec);
	shuffledTree->Branch("MC_weight_Kfactor",&MC_weight_Kfactor);
	shuffledTree->Branch("ZZMass", &ZZMass);
	shuffledTree->Branch("Z1Mass", &Z1Mass);
	shuffledTree->Branch("Z2Mass", &Z2Mass);
	shuffledTree->Branch("D_g1_vs_g2_phi0",&D_g2);
	shuffledTree->Branch("D_g1_vs_g4_phi0",&D_g4);
	shuffledTree->Branch("D_g2int_phi0",&D_g2int);
	shuffledTree->Branch("D_g4int_phi0",&D_g4int);
	shuffledTree->Branch("D_g2int_phi90",&D_g2int_perp);
	shuffledTree->Branch("D_g4int_phi90",&D_g4int_perp);
	shuffledTree->Branch("D_g1Q2_phi0",&D_g1q2);
	shuffledTree->Branch("D_g1Q2int_phi0",&D_g1q2int);
	shuffledTree->Branch("D_bkg",&D_bkg);
	shuffledTree->Branch("D_bkg_kin",&D_bkg_kin);
	shuffledTree->Branch("D_ggZZ",&D_ggZZ);
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

	TString TREE_NAME_GGZZ = TREE_NAME;
	TREE_NAME_GGZZ = TREE_NAME_GGZZ + "_ggZZ";
	TTree* shuffledggZZBkg = new TTree(TREE_NAME_GGZZ,TREE_NAME_GGZZ);
	shuffledggZZBkg->SetAutoSave(3000000000);
	shuffledggZZBkg->Branch("MC_weight",&MC_weight);
	shuffledggZZBkg->Branch("MC_weight_down",&MC_weight_down);
	shuffledggZZBkg->Branch("MC_weight_up",&MC_weight_up);
//	shuffledggZZBkg->Branch("MC_weight_noxsec",&MC_weight_noxsec);
	shuffledggZZBkg->Branch("MC_weight_xsec",&MC_weight_xsec);
	shuffledggZZBkg->Branch("MC_weight_Kfactor",&MC_weight_Kfactor);
	shuffledggZZBkg->Branch("ZZMass", &ZZMass);
	shuffledggZZBkg->Branch("Z1Mass", &Z1Mass);
	shuffledggZZBkg->Branch("Z2Mass", &Z2Mass);
	shuffledggZZBkg->Branch("D_g1_vs_g2_phi0",&D_g2);
	shuffledggZZBkg->Branch("D_g1_vs_g4_phi0",&D_g4);
	shuffledggZZBkg->Branch("D_g2int_phi0",&D_g2int);
	shuffledggZZBkg->Branch("D_g4int_phi0",&D_g4int);
	shuffledggZZBkg->Branch("D_g2int_phi90",&D_g2int_perp);
	shuffledggZZBkg->Branch("D_g4int_phi90",&D_g4int_perp);
	shuffledggZZBkg->Branch("D_g1Q2_phi0",&D_g1q2);
	shuffledggZZBkg->Branch("D_g1Q2int_phi0",&D_g1q2int);
	shuffledggZZBkg->Branch("D_bkg",&D_bkg);
	shuffledggZZBkg->Branch("D_bkg_kin",&D_bkg_kin);
	shuffledggZZBkg->Branch("D_ggZZ",&D_ggZZ);
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

	TString TREE_NAME_QQZZ = TREE_NAME;
	TREE_NAME_QQZZ = TREE_NAME_QQZZ + "_qqZZ";
	TTree* shuffledqqZZBkg = new TTree(TREE_NAME_QQZZ,TREE_NAME_QQZZ);
	shuffledqqZZBkg->SetAutoSave(3000000000);
	shuffledqqZZBkg->Branch("MC_weight",&MC_weight);
	shuffledqqZZBkg->Branch("MC_weight_down",&MC_weight_down);
	shuffledqqZZBkg->Branch("MC_weight_up",&MC_weight_up);
//	shuffledqqZZBkg->Branch("MC_weight_noxsec",&MC_weight_noxsec);
	shuffledqqZZBkg->Branch("MC_weight_xsec",&MC_weight_xsec);
	shuffledqqZZBkg->Branch("MC_weight_Kfactor",&MC_weight_Kfactor);
	shuffledqqZZBkg->Branch("ZZMass", &ZZMass);
	shuffledqqZZBkg->Branch("Z1Mass", &Z1Mass);
	shuffledqqZZBkg->Branch("Z2Mass", &Z2Mass);
	shuffledqqZZBkg->Branch("D_g1_vs_g2_phi0",&D_g2);
	shuffledqqZZBkg->Branch("D_g1_vs_g4_phi0",&D_g4);
	shuffledqqZZBkg->Branch("D_g2int_phi0",&D_g2int);
	shuffledqqZZBkg->Branch("D_g4int_phi0",&D_g4int);
	shuffledqqZZBkg->Branch("D_g2int_phi90",&D_g2int_perp);
	shuffledqqZZBkg->Branch("D_g4int_phi90",&D_g4int_perp);
	shuffledqqZZBkg->Branch("D_g1Q2_phi0",&D_g1q2);
	shuffledqqZZBkg->Branch("D_g1Q2int_phi0",&D_g1q2int);
	shuffledqqZZBkg->Branch("D_bkg",&D_bkg);
	shuffledqqZZBkg->Branch("D_bkg_kin",&D_bkg_kin);
	shuffledqqZZBkg->Branch("D_ggZZ",&D_ggZZ);
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

	TString TREE_NAME_ZX = TREE_NAME;
	TREE_NAME_ZX = TREE_NAME_ZX + "_ZX";
	TTree* shuffledZXBkg = new TTree(TREE_NAME_ZX,TREE_NAME_ZX);
	shuffledZXBkg->SetAutoSave(3000000000);
	shuffledZXBkg->Branch("MC_weight",&MC_weight);
	shuffledZXBkg->Branch("MC_weight_down",&MC_weight_down);
	shuffledZXBkg->Branch("MC_weight_up",&MC_weight_up);
//	shuffledZXBkg->Branch("MC_weight_noxsec",&MC_weight_noxsec);
	shuffledZXBkg->Branch("MC_weight_xsec",&MC_weight_xsec);
	shuffledZXBkg->Branch("MC_weight_Kfactor",&MC_weight_Kfactor);
	shuffledZXBkg->Branch("ZZMass", &ZZMass);
	shuffledZXBkg->Branch("Z1Mass", &Z1Mass);
	shuffledZXBkg->Branch("Z2Mass", &Z2Mass);
	shuffledZXBkg->Branch("D_g1_vs_g2_phi0",&D_g2);
	shuffledZXBkg->Branch("D_g1_vs_g4_phi0",&D_g4);
	shuffledZXBkg->Branch("D_g2int_phi0",&D_g2int);
	shuffledZXBkg->Branch("D_g4int_phi0",&D_g4int);
	shuffledZXBkg->Branch("D_g2int_phi90",&D_g2int_perp);
	shuffledZXBkg->Branch("D_g4int_phi90",&D_g4int_perp);
	shuffledZXBkg->Branch("D_g1Q2_phi0",&D_g1q2);
	shuffledZXBkg->Branch("D_g1Q2int_phi0",&D_g1q2int);
	shuffledZXBkg->Branch("D_bkg",&D_bkg);
	shuffledZXBkg->Branch("D_bkg_kin",&D_bkg_kin);
	shuffledZXBkg->Branch("D_ggZZ",&D_ggZZ);
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

	cout << "Starting to shuffle signal trees... Total N_events is " << sum_NGenTotal_unweighted << "." << endl;
	cout << "External weights for hypo=SM:" << endl;
	cout << "N_total Unweighted\nFile\tWeight" << endl;
	for(int l=0;l<kSignalSamples;l++) cout << l << '\t' << NGenTotal_unweighted[l] << endl;
	cout << "N_total Weighted\nFile\tWeight" << endl;
	for(int l=0;l<kSignalSamples;l++) cout << l << '\t' << NGenTotal_weighted[l][0] << endl;

	float nSM[kSignalSamples]={0};
	float nSM_2mu2e[kSignalSamples]={0};
	float nHypo[numSamples+2]={0};
	float nHypo_AsIf2e2mu[numSamples+2]={0};
	float nGGZZ=0;
	float nQQZZ=0;
	float nZX=0;
	
	int nUsed[(kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples)]={0};

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
//		EventSample = lucky;
		MC_weight_xsec = MC_weight/MC_weight_noxsec;
		if(lucky<kSignalSamples){
			for(int p=0;p<myNumSamples;p++){
				MC_weight_spin0[p] *= ( NGenTotal_unweighted[lucky]/NGenTotal_weighted[lucky][p] );
				float myxsec = MC_weight_xsec * MC_weight_spin0[p]  / sum_NGenTotal_unweighted ;
				MC_CV_weight[p] = myxsec*MC_weight_PUWeight*MC_weight_powhegWeight*MC_weight_dataMC*MC_weight_HqT;

				float ratioAsIf2mu2e = 0.5;
				if(flavor==2) ratioAsIf2mu2e = 1;
				ratioAsIf2mu2e *= NGenTotal_weighted_perFlavor[lucky][p][2];
				ratioAsIf2mu2e = NGenTotal_weighted_perFlavor[lucky][p][flavor] / ratioAsIf2mu2e;

				MC_CVwrt2mu2e_weight[p] = MC_CV_weight[p] * ratioAsIf2mu2e;
				nHypo[p] += MC_CV_weight[p];
				nHypo_AsIf2e2mu[p] += MC_CVwrt2mu2e_weight[p];
			};
			nHypo[myNumSamples+1] += MC_weight_xsec*MC_weight_PUWeight*MC_weight_powhegWeight*MC_weight_dataMC*MC_weight_HqT/sum_NGenTotal_unweighted;
			nHypo_AsIf2e2mu[myNumSamples+1] += MC_weight_xsec*MC_weight_PUWeight*MC_weight_powhegWeight*MC_weight_dataMC*MC_weight_HqT/sum_NGenTotal_unweighted;
			nSM[lucky] += MC_CV_weight[0];
			nSM_2mu2e[lucky] += MC_CVwrt2mu2e_weight[0];
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
		if(ZZMass>=240) nGGZZ += MC_weight;
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
		if(ZZMass>=240) nQQZZ += MC_weight;
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

		MC_weight_xsec = MC_weight/MC_weight_noxsec;
		for(int p=0;p<myNumSamples;p++){
			MC_CV_weight[p] = MC_weight;
		};
		nHypo[myNumSamples] += MC_CV_weight[0];
		if(ZZMass>=240) nZX += MC_weight;
		shuffledZXBkg->Fill();
		for(int smp=lucky;smp<(kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples);smp++) naccevents[smp] -= 1;

		ctr++;
	};
	foutput->WriteTObject(shuffledZXBkg);

	delete shuffledZXBkg;
	delete shuffledqqZZBkg;
	delete shuffledggZZBkg;
	delete shuffledTree;

	cout << "Sample SM integrals are " << endl;
	for(int f=0;f<kSignalSamples;f++) cout << f << ": " << nSM[f] << endl;
	cout << "Sample SM 2e2mu integrals are " << endl;
	for(int f=0;f<kSignalSamples;f++) cout << f << ": " << nSM_2mu2e[f] << endl;
	float sumnSM=0;
	float sumnSM_2mu2e=0;
	for(int f=0;f<kSignalSamples;f++){
		sumnSM += nSM[f];
		sumnSM_2mu2e += nSM_2mu2e[f];
	};
	cout << "Sum SM (with interf., without interf): " << sumnSM << '\t' << sumnSM_2mu2e << endl;
	cout << "Average hypothesis integral is " << nHypo[myNumSamples+1] << endl;
	cout << "Average hypothesis 2e2mu integral is " << nHypo_AsIf2e2mu[myNumSamples+1] << endl;
	cout << "Hypothesis integrals: " << endl;
	for(int p=0;p<myNumSamples+1;p++) cout << p << ": " << nHypo[p] << endl;
	cout << "nGGZZ (mZZ>=240): " << nGGZZ << endl;
	cout << "nQQZZ (mZZ>=240): " << nQQZZ << endl;
	cout << "nZX (mZZ>=240): " << nZX << endl;

	foutput->Close();
	for(int smp=0;smp<(kSignalSamples+kggZZSamples+kqqZZSamples+kZXSamples);smp++){
		delete myTree[smp];
		finput[smp]->Close();
	};
};
