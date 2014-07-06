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
#include "./data/ZZ4l_126_Samples.h"
#include "./data/FitDimensionsList.h"

using namespace std;


void shuffleSignalBkg(int flavor, int erg_tev){
	srand( time(0) );

	string OUTPUT_NAME = hzz4lprefix;
	OUTPUT_NAME = OUTPUT_NAME + "H126_ShuffledSignalBkg";
	TString TREE_NAME = "SelectedTree";
	TString TREE_GEN_NAME = "GenTree";

	char erg_dir[1000];
	sprintf(erg_dir,"LHC_%iTeV",erg_tev);
	TString cinput_common = user_dir + erg_dir;
	TString coutput_common = cinput_common;

	char cinput_KDFactor[1000];
	sprintf(cinput_KDFactor,"./data/HZZ4l-KDFactorGraph.root",user_dir.c_str());
	TFile* finput_KDFactor = new TFile(cinput_KDFactor,"read");
	string tgkfname = "KDScale_";
	tgkfname = tgkfname + "AllFlavors";
	TGraphAsymmErrors* tgkf = (TGraphAsymmErrors*) finput_KDFactor->Get(tgkfname.c_str());

	char* cGenArray[8]={"GenHMass","GenZ1Mass","GenZ2Mass","Gencosthetastar","GenhelcosthetaZ1","GenhelcosthetaZ2","Genhelphi","GenphistarZ1"};
	float kGenArray[8]={0};

	int genFinalState=0;
	float MC_weight,MC_weight_noxsec,MC_weight_xsec,PUWeight,powheg_weight,dataMCWeight,HqTWeight;
	float GenZZMass=0;
	float ZZMass=0;
	float MC_CV_weight[kNumSamples];
	float MC_spin0_wt[kNumSamples];
	float fitarray[kFitTotalD];
	float fitextraarray[kFitTotalExtraD];
	float myDbkg=0;
	int numSamples=kNumSamples;

	TFile* finput[kNumFiles+kNumBkg];
	TTree* myTree[kNumFiles+kNumBkg];
	int nevents[kNumFiles+kNumBkg]={0};
	int naccevents[kNumFiles+kNumBkg]={0};
	int nTotal=0;
	int nTotal_bkg=0;
	float NGenTotal_weighted[kNumFiles][kNumSamples];
	float NGenTotal_unweighted[kNumFiles]={0};
	float sum_NGenTotal_unweighted=0;
	for(int smp=0;smp<kNumFiles;smp++){
		for(int hypo=0;hypo<kNumSamples;hypo++) NGenTotal_weighted[smp][hypo]=0;
	};

	TFile* fGeninput[kNumFiles];
	TTree* myGenTree[kNumFiles];
	int nGenevents[kNumFiles]={0};
	int nGenaccevents[kNumFiles]={0};
	int nGenTotal=0;

	for(int smp=0;smp<(kNumFiles+kNumBkg);smp++){
		TString cinput = cinput_common + "/" + user_folder[flavor] + "/" + sample_file[smp] + ".root";
		if(smp>=kNumFiles) cinput = cinput_common + "/" + user_folder[flavor] + "/" + hzz4lprefix + sample_BackgroundFile[smp-kNumFiles] + ".root";

		cout << cinput << '\t' << smp << endl;

		finput[smp] = new TFile(cinput,"read");
		myTree[smp] = (TTree*) finput[smp]->Get(TREE_NAME);
		nevents[smp] = myTree[smp]->GetEntries();
		naccevents[smp] = nevents[smp];

		if(smp<kNumFiles){
			TString cinputGen = cinput_common + "/" + user_folder[flavor] + "/" + sample_file[smp] + "_GenLevel_ReWeighed.root";
			TFile* ftemp = new TFile(cinputGen,"read");
			TH2F* htrue = (TH2F*) ftemp->Get("hCounters_spin0_RW");
			for(int myf=0;myf<5;myf++){
				for(int p=0;p<kNumSamples;p++) NGenTotal_weighted[smp][p] += htrue->GetBinContent(myf+1,p+2);
				NGenTotal_unweighted[smp] += htrue->GetBinContent(myf+1,1);
			};
			sum_NGenTotal_unweighted += NGenTotal_unweighted[smp];
			delete htrue;
			ftemp->Close();

			nTotal += nevents[smp];

			fGeninput[smp] = new TFile(cinputGen,"read");
			myGenTree[smp] = (TTree*) fGeninput[smp]->Get(TREE_GEN_NAME);
			nGenevents[smp] = myGenTree[smp]->GetEntries();
			nGenaccevents[smp] = nGenevents[smp];
			nGenTotal += nGenevents[smp];
		}
		else{
			nTotal_bkg += nevents[smp];
		};
		if(smp!=0 && smp!=kNumFiles) naccevents[smp] += naccevents[smp-1];
		if(smp!=0 && smp<kNumFiles) nGenaccevents[smp] += nGenaccevents[smp-1];
		cout << nevents[smp] << endl;
		cout << naccevents[smp] << endl;
		cout << nGenaccevents[smp] << endl;

		myTree[smp]->SetBranchAddress("genFinalState",&genFinalState);
		myTree[smp]->SetBranchAddress("MC_weight",&MC_weight);
		myTree[smp]->SetBranchAddress("MC_weight_noxsec",&MC_weight_noxsec);
		myTree[smp]->SetBranchAddress("MC_weight_PUWeight",&PUWeight);
		myTree[smp]->SetBranchAddress("MC_weight_powhegWeight",&powheg_weight);
		myTree[smp]->SetBranchAddress("MC_weight_dataMC",&dataMCWeight);
		myTree[smp]->SetBranchAddress("MC_weight_HqT",&HqTWeight);
		myTree[smp]->SetBranchAddress("GenHMass",&GenZZMass);
		myTree[smp]->SetBranchAddress("ZZMass",&ZZMass);
		myTree[smp]->SetBranchAddress("kNumSamples",&numSamples);
		myTree[smp]->SetBranchAddress("MC_weight_spin0",MC_spin0_wt);

		myTree[smp]->SetBranchAddress(strDBkg,&myDbkg);
		for(int git=0;git<kFitTotalD;git++){
			if(git==kFitTotalD/2-1 || git==kFitTotalD-1) continue;
			fitarray[git]=0;
			int branch_set = myTree[smp]->SetBranchAddress(strFitDim[git],(fitarray+git));
			if(branch_set!=0){
				cerr << "Could NOT find branch named " << strFitDim[git] << "!!! Setting fitarray[" << git << "] = 0." << endl;
				fitarray[git]=0;
			};
		};
		for(int git = (kFitTotalExtraD-1);git>=(kFitTotalExtraD/2-1);git--){
			fitextraarray[git]=0;
			int branch_set = myTree[smp]->SetBranchAddress(strFitExtraDim[git],(fitextraarray+git));
			if(branch_set!=0){
				cerr << "Could NOT find branch named " << strFitExtraDim[git] << "!!! Setting fitextraarray[" << git << "] = 0." << endl;
				fitextraarray[git]=0;
			};
		};
		if(smp<kNumFiles){
			myGenTree[smp]->SetBranchAddress("genFinalState",&genFinalState);
			myGenTree[smp]->SetBranchAddress("GenZZMass",&ZZMass);
			myGenTree[smp]->SetBranchAddress("kNumSamples",&numSamples);
			myGenTree[smp]->SetBranchAddress("MC_weight_spin0",MC_spin0_wt);

			for(int git=0;git<8;git++){
				kGenArray[git]=0;
				int branch_set = myGenTree[smp]->SetBranchAddress(cGenArray[git],(kGenArray+git));
				if(branch_set!=0){
					cerr << "Could NOT find branch named " << cGenArray[git] << "!!! Setting GenArray[" << git << "] = 0." << endl;
					kGenArray[git]=0;
				};
			};
			for(int git=0;git<kFitTotalD;git++){
				if(git==kFitTotalD/2-1 || git==kFitTotalD-1) continue;
				fitarray[git]=0;
				int branch_set = myGenTree[smp]->SetBranchAddress(strFitDim[git],(fitarray+git));
				if(branch_set!=0){
					cerr << "Could NOT find branch named " << strFitDim[git] << "!!! Setting fitarray[" << git << "] = 0." << endl;
					fitarray[git]=0;
				};
			};
			for(int git = (kFitTotalExtraD-1);git>=(kFitTotalExtraD/2-1);git--){
				fitextraarray[git]=0;
				int branch_set = myGenTree[smp]->SetBranchAddress(strFitExtraDim[git],(fitextraarray+git));
				if(branch_set!=0){
					cerr << "Could NOT find branch named " << strFitExtraDim[git] << "!!! Setting fitextraarray[" << git << "] = 0." << endl;
					fitextraarray[git]=0;
				};
			};
		};
	};

	TString coutput = coutput_common + "/" + user_folder[flavor] + "/" + OUTPUT_NAME + ".root";
	TFile* foutput = new TFile(coutput,"recreate");

	int EventSample=0;

	TTree* shuffledTree = new TTree(TREE_NAME,TREE_NAME);
	shuffledTree->SetAutoSave(3000000000);
	shuffledTree->Branch("EventSample",&EventSample);
	shuffledTree->Branch("genFinalState",&genFinalState);
	shuffledTree->Branch("MC_weight_xsec",&MC_weight_xsec);
	shuffledTree->Branch("ZZMass",&ZZMass);
	shuffledTree->Branch("kNumSamples",&numSamples);
	shuffledTree->Branch("MC_CV_weight",MC_CV_weight,"MC_CV_weight[kNumSamples]/F");
	shuffledTree->Branch("MC_weight_spin0",MC_spin0_wt,"MC_weight_spin0[kNumSamples]/F");
	shuffledTree->Branch(strDBkg,&myDbkg);
	for(int git=0;git<kFitTotalD;git++){
		if(git==kFitTotalD/2-1 || git==kFitTotalD-1) continue;
		shuffledTree->Branch(strFitDim[git],(fitarray+git));
	};
	for(int git = (kFitTotalExtraD-1);git>=(kFitTotalExtraD/2-1);git--){
		shuffledTree->Branch(strFitExtraDim[git],(fitextraarray+git));
	};

	TString TREE_NAME_BKG = TREE_NAME;
	TREE_NAME_BKG = TREE_NAME_BKG + "_Bkg";
	TTree* shuffledBkg = new TTree(TREE_NAME_BKG,TREE_NAME_BKG);
	shuffledBkg->SetAutoSave(3000000000);
	shuffledBkg->Branch("EventSample",&EventSample);
	shuffledBkg->Branch("genFinalState",&genFinalState);
	shuffledBkg->Branch("MC_weight_xsec",&MC_weight_xsec);
	shuffledBkg->Branch("ZZMass",&ZZMass);
	shuffledBkg->Branch("kNumSamples",&numSamples);
	shuffledBkg->Branch("MC_CV_weight",MC_CV_weight,"MC_CV_weight[kNumSamples]/F");
	shuffledBkg->Branch("MC_weight_spin0",MC_spin0_wt,"MC_weight_spin0[kNumSamples]/F");
	shuffledBkg->Branch(strDBkg,&myDbkg);
	for(int git=0;git<kFitTotalD;git++){
		if(git==kFitTotalD/2-1 || git==kFitTotalD-1) continue;
		shuffledBkg->Branch(strFitDim[git],(fitarray+git));
	};
	for(int git = (kFitTotalExtraD-1);git>=(kFitTotalExtraD/2-1);git--){
		shuffledBkg->Branch(strFitExtraDim[git],(fitextraarray+git));
	};

	TTree* shuffledGenTree = new TTree(TREE_GEN_NAME,TREE_GEN_NAME);
	shuffledGenTree->SetAutoSave(3000000000);
	shuffledGenTree->Branch("EventSample",&EventSample);
	shuffledGenTree->Branch("genFinalState",&genFinalState);
	shuffledGenTree->Branch("GenHMass",&ZZMass);
	shuffledGenTree->Branch("kNumSamples",&numSamples);
	shuffledGenTree->Branch("MC_weight_spin0",MC_spin0_wt,"MC_weight_spin0[kNumSamples]/F");
	for(int git=0;git<kFitTotalD;git++){
		if(git==kFitTotalD/2-1 || git==kFitTotalD-1) continue;
		shuffledGenTree->Branch(strFitDim[git],(fitarray+git));
	};
	for(int git = (kFitTotalExtraD-1);git>=(kFitTotalExtraD/2-1);git--){
		shuffledGenTree->Branch(strFitExtraDim[git],(fitextraarray+git));
	};
	for(int git=0;git<8;git++){
		shuffledGenTree->Branch(cGenArray[git],(kGenArray+git));
	};


	cout << "Starting to shuffle signal trees... Total N_events is " << sum_NGenTotal_unweighted << "." << endl;
	cout << "External weights for hypo=SM:" << endl;
	cout << "N_total Unweighted\nFile\tWeight" << endl;
	for(int l=0;l<kNumFiles;l++) cout << l << '\t' << NGenTotal_unweighted[l] << endl;
	cout << "N_total Weighted\nFile\tWeight" << endl;
	for(int l=0;l<kNumFiles;l++) cout << l << '\t' << NGenTotal_weighted[l][0] << endl;

	float nSM[kNumFiles]={0};
	float nSM_2mu2e[kNumFiles]={0};
	float nHypo[kNumSamples+2]={0};
	int nUsed[(kNumFiles+kNumBkg)]={0};
	float nGenSM[kNumFiles]={0};
	float nGenHypo[kNumSamples+2]={0};
	int nGenUsed[kNumFiles]={0};
	int ctr=0;
	while(ctr<nTotal){
		int coin = rand() % naccevents[(kNumFiles-1)];
		int lucky=-1;
		for(int smp=0;smp<(kNumFiles);smp++){
			if(coin<naccevents[smp] && (nevents[smp]-nUsed[smp])!=0){lucky=smp;break;};
		};
		if(lucky<0) continue;

		myTree[lucky]->GetEntry(nUsed[lucky]);
		nUsed[lucky] += 1;
		EventSample = lucky;
		MC_weight_xsec = MC_weight/MC_weight_noxsec;
		float myKDwgt = tgkf->Eval(GenZZMass);
		MC_weight_xsec *= myKDwgt; // Weight for signal lineshape; comes from K factor
		if(lucky<kNumFiles){
			for(int p=0;p<kNumSamples;p++){
				MC_spin0_wt[p] *= ( NGenTotal_unweighted[lucky]/NGenTotal_weighted[lucky][p] );
				float myxsec = MC_weight_xsec * MC_spin0_wt[p]  / sum_NGenTotal_unweighted ;
				MC_CV_weight[p] = myxsec*PUWeight*powheg_weight*dataMCWeight*HqTWeight;
				nHypo[p] += MC_CV_weight[p];
			};
			nHypo[kNumSamples+1] += MC_weight_xsec*PUWeight*powheg_weight*dataMCWeight*HqTWeight/sum_NGenTotal_unweighted;
			nSM[lucky] += MC_CV_weight[0];
			if(genFinalState==2) nSM_2mu2e[lucky] += MC_CV_weight[0];
		};
		shuffledTree->Fill();
		for(int smp=lucky;smp<(kNumFiles);smp++) naccevents[smp] -= 1;
		ctr++;
	};
	foutput->WriteTObject(shuffledTree);

	ctr=0;
	while(ctr<nTotal_bkg){
		int coin = rand() % naccevents[(kNumFiles+kNumBkg-1)];
		int lucky=-1;
		for(int smp=kNumFiles;smp<(kNumFiles+kNumBkg);smp++){
			if(coin<naccevents[smp] && (nevents[smp]-nUsed[smp])!=0){lucky=smp;break;};
		};
		if(lucky<0) continue;

		myTree[lucky]->GetEntry(nUsed[lucky]);
		nUsed[lucky] += 1;
		EventSample = lucky;
		float myKDwgt = tgkf->Eval(GenZZMass);
		if(lucky<kNumFiles+2){
			MC_weight *= myKDwgt;
		};
		MC_weight_xsec = MC_weight/MC_weight_noxsec;
		for(int p=0;p<kNumSamples;p++){
			MC_CV_weight[p] = MC_weight;
		};
		nHypo[kNumSamples] += MC_CV_weight[0];
		shuffledBkg->Fill();
		for(int smp=lucky;smp<(kNumFiles+kNumBkg);smp++) naccevents[smp] -= 1;

		ctr++;
	};
	foutput->WriteTObject(shuffledBkg);

	ctr=0;
	while(ctr<nGenTotal){
		int coin = rand() % nGenaccevents[kNumFiles-1];
		int lucky=-1;
		for(int smp=0;smp<kNumFiles;smp++){
			if(coin<nGenaccevents[smp] && (nGenevents[smp]-nGenUsed[smp])!=0){lucky=smp;break;};
		};
		if(lucky<0) continue;

		myGenTree[lucky]->GetEntry(nGenUsed[lucky]);
		nGenUsed[lucky] += 1;
		EventSample = lucky;
		float myKDwgt = tgkf->Eval(GenZZMass);
		for(int p=0;p<kNumSamples;p++){
			MC_spin0_wt[p] *= ( ( NGenTotal_unweighted[lucky]/NGenTotal_weighted[lucky][p] ) / sum_NGenTotal_unweighted )*myKDwgt;
		};
		shuffledGenTree->Fill();
		for(int smp=lucky;smp<kNumFiles;smp++) nGenaccevents[smp] -= 1;
		ctr++;
	};
	foutput->WriteTObject(shuffledGenTree);

	delete shuffledTree;
	delete shuffledGenTree;
	delete shuffledBkg;

	cout << "Sample SM integrals are " << endl;
	for(int f=0;f<kNumFiles;f++) cout << f << ": " << nSM[f] << endl;
	cout << "Sample SM 2e2mu integrals are " << endl;
	for(int f=0;f<kNumFiles;f++) cout << f << ": " << nSM_2mu2e[f] << endl;
	cout << "Average hypothesis integral is " << nHypo[kNumSamples+1] << endl;
	cout << "Hypothesis integrals: " << endl;
	for(int p=0;p<kNumSamples+1;p++) cout << p << ": " << nHypo[p] << endl;

	foutput->Close();
	for(int smp=0;smp<kNumFiles;smp++){
		delete myGenTree[smp];
		fGeninput[smp]->Close();
	};
	for(int smp=0;smp<(kNumFiles+kNumBkg);smp++){
		delete myTree[smp];
		finput[smp]->Close();
	};
	finput_KDFactor->Close();
};
