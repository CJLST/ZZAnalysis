#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include <ZZMatrixElement/MELA/interface/Mela.h>
//#include <ZZMatrixElement/MELA/src/computeAngles.h>
#include "./data/ZZ4l_125p6_Samples.h"
#include "./data/FitDimensionsList.h"

using namespace RooFit;

using namespace std;

namespace {
  const bool doMuEffCorr    = true;
  const bool doEleEffCorr   = true;
  const bool doHighmassCorr = true; 
  const bool doHqTCorr      = false;
  const bool saveJets       = true;
  const bool applySystMu    = false;
  const bool applySystEle   = false;

  const int pwhg_flag = 0;    // 0 means standard high mass weights
                                // 1 means CPS+ + Interference
                                // 2 means CPS- + Interference
                                // 3 means CPS  + Interference+
                                // 4 means CPS  + Interference-

  const float Run2011AFraction = 0.465;
  const float Zmass = 91.1876;
  const float gamZ = 2.5;
  const float R1Val = 0.15;
  const float R2Val = 0.15;
  const int para = 1;
  const bool use_acc=false;
  const float M_muon = 0.105658389;
  const float M_electron = 0.00051099907;
  const float M_tau = 1.777;
  const int PDG_electron=11,PDG_muon=13,PDG_tau=15;
}


float getJHUGenMELAWeight(Mela& myMela, int lepId[4], float angularOrdered[8], double selfDHvvcoupl[31][2]){
	float myprob=1.0;
	int myflavor=-1;
	if(abs(lepId[0])==abs(lepId[1]) &&
		abs(lepId[0])==abs(lepId[2]) &&
		abs(lepId[0])==abs(lepId[3])){
			if(abs(lepId[0])==11) myflavor=1;
			else myflavor=2;
	}
	else myflavor=3;

	if(myflavor>=0) myMela.computeP(angularOrdered[0],angularOrdered[1],angularOrdered[2],angularOrdered[3],
		angularOrdered[4],angularOrdered[5],angularOrdered[6],angularOrdered[7],
	    myflavor,
	    selfDHvvcoupl,
	    myprob
		);
	return myprob;
}
/*
void testSpin0MEDivergence(int iSample, int iHypo, float& MEVal){
	if (iSample == kfLambda1_1){
		if(iHypo == kfZG_0_fGG_1 && MEVal>2.0e14) MEVal = pow(2.0e14,2)/MEVal;
		if(iHypo == kfMZG_0_fMGG_1 && MEVal>2.0e14) MEVal = pow(2.0e14,2)/MEVal;
		if(iHypo == kfZG_0_fGG_05 && MEVal>8.0e11) MEVal = pow(8.0e11,2)/MEVal;
		if(iHypo == kfMZG_0_fMGG_05 && MEVal>6.5e11) MEVal = pow(6.5e11,2)/MEVal;
	};
	if (iSample == kfLambda1_m05){
		if(iHypo == kfZG_0_fGG_1 && MEVal>4.0e5) MEVal = pow(4.0e5,2)/MEVal;
		if(iHypo == kfMZG_0_fMGG_1 && MEVal>4.0e5) MEVal = pow(4.0e5,2)/MEVal;
		if(iHypo == kfZG_0_fGG_05 && MEVal>1000) MEVal = pow(1000.0,2)/MEVal;
//		if(iHypo == kfMZG_0_fMGG_05 && MEVal>150) MEVal = pow(150.0,2)/MEVal;
	};
	if (iSample == kfLambda1_05){
		if(iHypo == kfZG_1_fGG_0 && MEVal>5.0e7) MEVal = pow(5.0e7,2)/MEVal;
		if(iHypo == kfMZG_1_fMGG_0 && MEVal>5.0e7) MEVal = pow(5.0e7,2)/MEVal;
		if(iHypo == kfZG_0_fGG_1 && MEVal>2.0e7) MEVal = pow(2.0e7,2)/MEVal;
		if(iHypo == kfMZG_0_fMGG_1 && MEVal>2.0e7) MEVal = pow(2.0e7,2)/MEVal;
		if(iHypo == kfZG_0_fGG_05 && MEVal>6000) MEVal = pow(6000.0,2)/MEVal;
	};
	if (iSample == kfg2_0_fg4_1){
		if(iHypo == kfZG_0_fGG_1 && MEVal>2.0e7) MEVal = pow(2.0e7,2)/MEVal;
		if(iHypo == kfMZG_0_fMGG_1 && MEVal>2.0e7) MEVal = pow(2.0e7,2)/MEVal;
		if(iHypo == kfZG_0_fGG_05 && MEVal>30000) MEVal = pow(30000.0,2)/MEVal;
//		if(iHypo == kfMZG_0_fMGG_05 && MEVal>8000) MEVal = pow(8000.0,2)/MEVal;
	};
	if (iSample == kfg2_0_fg4_05){
		if(iHypo == kfZG_0_fGG_1 && MEVal>1.5e6) MEVal = pow(1.5e6,2)/MEVal;
		if(iHypo == kfMZG_0_fMGG_1 && MEVal>1.5e6) MEVal = pow(1.5e6,2)/MEVal;
		if(iHypo == kfZG_0_fGG_05 && MEVal>4000) MEVal = pow(4000.0,2)/MEVal;
//		if(iHypo == kfMZG_0_fMGG_05 && MEVal>900) MEVal = pow(900.0,2)/MEVal;
	};
	if (iSample == kfg2_0_fg4_01){
		if(iHypo == kfZG_0_fGG_1 && MEVal>2.0e6) MEVal = pow(2.0e6,2)/MEVal;
		if(iHypo == kfMZG_0_fMGG_1 && MEVal>2.0e6) MEVal = pow(2.0e6,2)/MEVal;
		if(iHypo == kfZG_0_fGG_05 && MEVal>6000) MEVal = pow(6000.0,2)/MEVal;
//		if(iHypo == kfMZG_0_fMGG_05 && MEVal>1200) MEVal = pow(1200.0,2)/MEVal;
	};
	if (iSample == kfg2_0_fg4_01_p390-1){
		if(iHypo == kfZG_0_fGG_1 && MEVal>1.0e6) MEVal = pow(1.0e6,2)/MEVal;
		if(iHypo == kfMZG_0_fMGG_1 && MEVal>1.0e6) MEVal = pow(1.0e6,2)/MEVal;
		if(iHypo == kfZG_0_fGG_05 && MEVal>3000) MEVal = pow(3000.0,2)/MEVal;
//		if(iHypo == kfMZG_0_fMGG_05 && MEVal>600) MEVal = pow(600.0,2)/MEVal;
	};
	if (iSample == kfg2_0_fg4_05_p3Pi-1){
		if(iHypo == kfZG_0_fGG_1 && MEVal>2.0e6) MEVal = pow(2.0e6,2)/MEVal;
		if(iHypo == kfMZG_0_fMGG_1 && MEVal>2.0e6) MEVal = pow(2.0e6,2)/MEVal;
		if(iHypo == kfZG_0_fGG_05 && MEVal>6000) MEVal = pow(6000.0,2)/MEVal;
//		if(iHypo == kfMZG_0_fMGG_05 && MEVal>1200) MEVal = pow(1200.0,2)/MEVal;
	};
	if (iSample == kfg2_05_fg4_0){
		if(iHypo == kfZG_0_fGG_1 && MEVal>3.0e6) MEVal = pow(3.0e6,2)/MEVal;
		if(iHypo == kfMZG_0_fMGG_1 && MEVal>3.0e6) MEVal = pow(3.0e6,2)/MEVal;
		if(iHypo == kfZG_0_fGG_05 && MEVal>1000) MEVal = pow(1000.0,2)/MEVal;
//		if(iHypo == kfMZG_0_fMGG_05 && MEVal>150) MEVal = pow(150.0,2)/MEVal;
	};
	if (iSample == kfg2_05_fg4_0_p290-1){
		if(iHypo == kfZG_0_fGG_1 && MEVal>1.0e6) MEVal = pow(1.0e6,2)/MEVal;
		if(iHypo == kfMZG_0_fMGG_1 && MEVal>1.0e6) MEVal = pow(1.0e6,2)/MEVal;
		if(iHypo == kfZG_0_fGG_05 && MEVal>3000) MEVal = pow(3000.0,2)/MEVal;
//		if(iHypo == kfMZG_0_fMGG_05 && MEVal>600) MEVal = pow(600.0,2)/MEVal;
	};
	if (iSample == kfg2_05_fg4_0_p2Pi-1){
		if(iHypo == kfZG_0_fGG_1 && MEVal>1.0e6) MEVal = pow(1.0e6,2)/MEVal;
		if(iHypo == kfMZG_0_fMGG_1 && MEVal>1.0e6) MEVal = pow(1.0e6,2)/MEVal;
		if(iHypo == kfZG_0_fGG_05 && MEVal>20000) MEVal = pow(20000.0,2)/MEVal;
//		if(iHypo == kfMZG_0_fMGG_05 && MEVal>3000) MEVal = pow(3000.0,2)/MEVal;
	};
	if (iSample == kfg2_05_fg4_05){
		if(iHypo == kfZG_0_fGG_1 && MEVal>5.0e6) MEVal = pow(5.0e6,2)/MEVal;
		if(iHypo == kfMZG_0_fMGG_1 && MEVal>5.0e6) MEVal = pow(5.0e6,2)/MEVal;
		if(iHypo == kfZG_0_fGG_05 && MEVal>20000) MEVal = pow(20000.0,2)/MEVal;
//		if(iHypo == kfMZG_0_fMGG_05 && MEVal>6000) MEVal = pow(6000.0,2)/MEVal;
	};
	if (iSample == kfg2_05_fg4_05_p2Pi-1){
		if(iHypo == kfZG_0_fGG_1 && MEVal>3.0e6) MEVal = pow(3.0e6,2)/MEVal;
		if(iHypo == kfMZG_0_fMGG_1 && MEVal>3.0e6) MEVal = pow(3.0e6,2)/MEVal;
		if(iHypo == kfZG_0_fGG_05 && MEVal>20000) MEVal = pow(20000.0,2)/MEVal;
//		if(iHypo == kfMZG_0_fMGG_05 && MEVal>2000) MEVal = pow(2000.0,2)/MEVal;
	};
	if (iSample == kfg2_33_fg4_33_p390-1){
		if(iHypo == kfZG_0_fGG_1 && MEVal>3.0e5) MEVal = pow(3.0e5,2)/MEVal;
		if(iHypo == kfMZG_0_fMGG_1 && MEVal>3.0e5) MEVal = pow(3.0e5,2)/MEVal;
		if(iHypo == kfZG_0_fGG_05 && MEVal>800) MEVal = pow(800.0,2)/MEVal;
//		if(iHypo == kfMZG_0_fMGG_05 && MEVal>150) MEVal = pow(150.0,2)/MEVal;
	};
	if (iSample == kfg2_01_fg4_01_p390-1){
		if(iHypo == kfZG_0_fGG_1 && MEVal>1.6e6) MEVal = pow(1.6e6,2)/MEVal;
		if(iHypo == kfMZG_0_fMGG_1 && MEVal>1.6e6) MEVal = pow(1.6e6,2)/MEVal;
		if(iHypo == kfZG_0_fGG_05 && MEVal>5000) MEVal = pow(5000.0,2)/MEVal;
//		if(iHypo == kfMZG_0_fMGG_05 && MEVal>500) MEVal = pow(500.0,2)/MEVal;
	};
};
*/
void reweightIdeal_125p6(int erg_tev, int smp_min=0, int smp_max=kNumFiles, float mPOLE=125.6){
	char GenLevel_Location[]="/GenSignal";
	char TREE_NAME[] = "GenTree";
	char erg_dir[1000];
	sprintf(erg_dir,"LHC_%iTeV",erg_tev);
	char cerg[1000];
	sprintf(cerg,"%iTeV",erg_tev);
	TString cinput_common = user_dir + erg_dir;
//	TString coutput_common = cinput_common;
	TString coutput_common = user_dir_newProduction + cerg;

	Mela mela(erg_tev,mPOLE);
	mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);

	const int firstNonZZHypothesis = kfZG_1_fGG_0;
	cout << "Will apply a 4 GeV cut to hypotheses " << firstNonZZHypothesis << "+..." << endl;
	if(smp_max>kNumFiles) smp_max=kNumFiles; 
	for(int smp=smp_min; smp<smp_max; smp++){
		TString cinput = cinput_common + GenLevel_Location + "/" + sample_file[smp] + "_Generated.root";
		TString coutput = coutput_common + "/" + GenLevel_Location + "/" + sample_file[smp] + "_GenLevel_ReWeighed.root";

		TFile* finput = new TFile(cinput,"read");
		TTree* tree = (TTree*) finput->Get(TREE_NAME);
		tree->SetAutoSave(3000000000);

		int genFinalState;
		float GenHMass,GenZ1Mass,GenZ2Mass,GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id,myGencosthetastar,myGenhelphi,myGenhelcosthetaZ1,myGenhelcosthetaZ2,myGenphistarZ1,myGenphistarZ2,sample_probPdf;
		tree->SetBranchAddress("genFinalState", &genFinalState);
		tree->SetBranchAddress("GenZZMass", &GenHMass);
		tree->SetBranchAddress("GenZ1Mass", &GenZ1Mass);
		tree->SetBranchAddress("GenZ2Mass", &GenZ2Mass);
		tree->SetBranchAddress("GenLep1Id", &GenLep1Id);
		tree->SetBranchAddress("GenLep2Id", &GenLep2Id);
		tree->SetBranchAddress("GenLep3Id", &GenLep3Id);
		tree->SetBranchAddress("GenLep4Id", &GenLep4Id);
		tree->SetBranchAddress("Gencosthetastar",&myGencosthetastar);
		tree->SetBranchAddress("Genhelphi",&myGenhelphi);
		tree->SetBranchAddress("GenhelcosthetaZ1",&myGenhelcosthetaZ1);
		tree->SetBranchAddress("GenhelcosthetaZ2",&myGenhelcosthetaZ2);
		tree->SetBranchAddress("GenphistarZ1",&myGenphistarZ1);
		tree->SetBranchAddress("GenphistarZ2",&myGenphistarZ2);

		TFile* foutput = new TFile(coutput,"recreate");

		mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);

		float NGen_weighted=0;
		float MC_weight_samples[kNumSamples];
		int numSamples = kNumSamples;
		for(int s=0;s<kNumSamples;s++){
			MC_weight_samples[s]=1;
		};
		foutput->cd();
		TH2F* hCount_new = new TH2F("hCounters_spin0_RW","Counters with Spin0 Re-weights",nFinalStates,0,nFinalStates,numSamples+1,0,numSamples+1);
		TH2F* hCount_4GeVcut_new = new TH2F("hCounters_spin0_4GeVcut_RW","Counters with Spin0 Re-weights with 4GeV m1,2 cut",nFinalStates,0,nFinalStates,numSamples+1,0,numSamples+1);

		TTree* mytree = (TTree*) tree->CloneTree(0);
		mytree->SetAutoSave(3000000000);
		mytree->Branch("kNumSamples",&numSamples);
		mytree->Branch("MC_weight_spin0",MC_weight_samples,"MC_weight_spin0[kNumSamples]/F");
		mytree->Branch("sampleprob",&sample_probPdf);

		const float g1q2scale = 3.28;
		const float g1g2scale_old = 2.1;
		const float g1g4scale_old = 6.0;
		const float g1g2scale = 2.134;
		const float g1g4scale = 5.548;
		const float g2g4scale = 1.539;
		const float g1q2intscale = 1.0;
		const float g1g2intscale = 1.0;
		const float g1g4intscale = 0.5/0.2;
		const float g1gZGintscale = 0.0688;
		const float g1gGGintscale = -0.0898;
		const float g4gZGscale = 0.0351;
		const float g4gGGscale = -0.0372;
		float Gen_D_g1q2=-99,Gen_D_g1q2int=-99;
		float Gen_D_g2=-99,Gen_D_g2int=-99,Gen_D_g2int_perp=-99;
		float Gen_D_g4=-99,Gen_D_g4int=-99,Gen_D_g4int_perp=-99;
		float Gen_D_ZG, Gen_D_GG;
		float Gen_D_ZGint, Gen_D_GGint;
		float Gen_D_ZG_PS, Gen_D_GG_PS;
		float Gen_D_ZG_PSint, Gen_D_GG_PSint;
/*
		mytree->Branch("D_g1Q2_phi0",&Gen_D_g1q2);
		mytree->Branch("D_g1_vs_g2_phi0",&Gen_D_g2);
		mytree->Branch("D_g1_vs_g4_phi0",&Gen_D_g4);
//		mytree->Branch("D_g1Q2intPdf_phi0",&Gen_D_g1q2int);
		mytree->Branch("D_g2intPdf_phi0",&Gen_D_g2int);
		mytree->Branch("D_g4intPdf_phi0",&Gen_D_g4int);
		mytree->Branch("D_g2intPdf_phi90",&Gen_D_g2int_perp);
		mytree->Branch("D_g4intPdf_phi90",&Gen_D_g4int_perp);
		mytree->Branch("Gen_D_ZG",&Gen_D_ZG);
		mytree->Branch("Gen_D_GG",&Gen_D_GG);
		mytree->Branch("Gen_D_ZGint",&Gen_D_ZGint);
		mytree->Branch("Gen_D_GGint",&Gen_D_GGint);
		mytree->Branch("Gen_D_ZG_PS",&Gen_D_ZG_PS);
		mytree->Branch("Gen_D_GG_PS",&Gen_D_GG_PS);
		mytree->Branch("Gen_D_ZG_PSint",&Gen_D_ZG_PSint);
		mytree->Branch("Gen_D_GG_PSint",&Gen_D_GG_PSint);
*/
		float N_generated[nFinalStates][kNumSamples+1];
		for(int xb=0;xb<kNumSamples+1;xb++){ for(int yb=0;yb<nFinalStates;yb++) N_generated[yb][xb]=0;};
		float N_generated_4GeVcut[nFinalStates][kNumSamples+1];
		for(int xb=0;xb<kNumSamples+1;xb++){ for(int yb=0;yb<nFinalStates;yb++) N_generated_4GeVcut[yb][xb]=0;};

		const int size_zzcoupl=SIZE_HVV; // Directly from Mela package: TVar.hh
		double selfDHvvcoupl[size_zzcoupl][2];
		for(int gx=0;gx<size_zzcoupl;gx++){
			selfDHvvcoupl[gx][0]=0;
			selfDHvvcoupl[gx][1]=0;
		};

		cout << "MELA Gen. ME width: " << gi_phi2_phi4_files[smp][p_GH] << endl;
		int nEntries = tree->GetEntries();
		cout << "Entries " << nEntries << endl;
		for(int ev=0;ev<nEntries;ev++){
//		for(int ev=0;ev<100;ev++){
			tree->GetEntry(ev);
			sample_probPdf=1.0;

			if((ev%50000)==0) cout << "Event: " <<  ev << endl;

			int lepIdOrdered[4]={ GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id };
			float angularOrdered[8]={GenHMass,GenZ1Mass,GenZ2Mass,myGencosthetastar,myGenhelcosthetaZ1,myGenhelcosthetaZ2,myGenhelphi,myGenphistarZ1};

			if(abs(lepIdOrdered[0])==15 ||
				abs(lepIdOrdered[1])==15 ||
				abs(lepIdOrdered[2])==15 ||
				abs(lepIdOrdered[3])==15
				) genFinalState=4;
			else if(abs(lepIdOrdered[0])==abs(lepIdOrdered[1]) &&
					abs(lepIdOrdered[0])==abs(lepIdOrdered[2]) &&
					abs(lepIdOrdered[0])==abs(lepIdOrdered[3])
					){
					if(abs(lepIdOrdered[0])==11) genFinalState=1;
					else if(abs(lepIdOrdered[0])==13) genFinalState=0;
			}
			else genFinalState=2;

			for(int gx=0;gx<size_zzcoupl;gx++){
				selfDHvvcoupl[gx][0]=0;
				selfDHvvcoupl[gx][1]=0;
			};

			selfDHvvcoupl[0][0] = gi_phi2_phi4_files[smp][p_Re_g1]; // Re-g1
			selfDHvvcoupl[0][1] = gi_phi2_phi4_files[smp][p_Im_g1]; // Im-g1
			selfDHvvcoupl[1][0] = gi_phi2_phi4_files[smp][p_Re_g2]; // Re-g2
			selfDHvvcoupl[1][1] = gi_phi2_phi4_files[smp][p_Im_g2]; // Im-g2
			selfDHvvcoupl[3][0] = gi_phi2_phi4_files[smp][p_Re_g4]; // Re-g4
			selfDHvvcoupl[3][1] = gi_phi2_phi4_files[smp][p_Im_g4]; // Im-g4
			selfDHvvcoupl[11][0] = gi_phi2_phi4_files[smp][p_Re_L1];  // Re-g1prime2
			selfDHvvcoupl[11][1] = gi_phi2_phi4_files[smp][p_Im_L1];  // Im-g1prime2
			selfDHvvcoupl[4][0] = gi_phi2_phi4_files[smp][p_Re_a2ZG]; // Re-g_zgsg_2
			selfDHvvcoupl[4][1] = gi_phi2_phi4_files[smp][p_Im_a2ZG]; // Im-g_zgsg_2
			selfDHvvcoupl[7][0] = gi_phi2_phi4_files[smp][p_Re_a2GG]; // Re-g_gsgs_2
			selfDHvvcoupl[7][1] = gi_phi2_phi4_files[smp][p_Im_a2GG]; // Im-g_gsgs_2
			selfDHvvcoupl[6][0] = gi_phi2_phi4_files[smp][p_Re_a3ZG]; // Re-g_zgs_4
			selfDHvvcoupl[6][1] = gi_phi2_phi4_files[smp][p_Im_a3ZG]; // Im-g_zgs_4
			selfDHvvcoupl[9][0] = gi_phi2_phi4_files[smp][p_Re_a3GG]; // Re-g_gsgs_4
			selfDHvvcoupl[9][1] = gi_phi2_phi4_files[smp][p_Im_a3GG]; // Im-g_gsgs_4
			selfDHvvcoupl[30][0] = gi_phi2_phi4_files[smp][p_Re_ZGL1]; // Re-g_zgs_1prime2
			selfDHvvcoupl[30][1] = gi_phi2_phi4_files[smp][p_Im_ZGL1]; // Im-g_zgs_1prime2
			selfDHvvcoupl[13][0] = gi_phi2_phi4_files[smp][p_Re_LQ]; // Re-g1prime4 (g1-LambdaQ)
			selfDHvvcoupl[13][1] = gi_phi2_phi4_files[smp][p_Im_LQ]; // Im-g1prime4 (g1-LambdaQ)

			mela.setMelaHiggsWidth(gi_phi2_phi4_files[smp][p_GH]);
			sample_probPdf = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl); // MELA automatically resets width to default after calculation as part of protection for discriminants
			bool isCountable=true;
			if(sample_probPdf==0 || sample_probPdf!=sample_probPdf){
				sample_probPdf=1.0;
				isCountable=false;
			};

			float weight_probPdf=1.0;

			Gen_D_ZG = -99;
			Gen_D_GG = -99;
			Gen_D_ZG_PS = -99;
			Gen_D_GG_PS = -99;
			Gen_D_g2 = -99;
			Gen_D_g4 = -99;
			Gen_D_g1q2 = -99;
			Gen_D_g2int = -99;
			Gen_D_g4int = -99;
			Gen_D_g1q2int = -99;
			Gen_D_g2int_perp = -99;
			Gen_D_g4int_perp = -99;
			Gen_D_ZGint=-99;
			Gen_D_GGint=-99;
			Gen_D_ZG_PSint=-99;
			Gen_D_GG_PSint=-99;

			float g1probPdf_true=1.0;
			float g2probPdf_true=1.0;
			float g4probPdf_true=1.0;
			float g1q2probPdf_true=1.0;
			float gZGprobPdf_true=1.0;
			float gGGprobPdf_true=1.0;
			float gPSZGprobPdf_true=1.0;
			float gPSGGprobPdf_true=1.0;

			float g1g2probPdf_true=1.0;
			float g1g4probPdf_true=1.0;
			float g1g1q2probPdf_true=1.0;
			float g1gZGprobPdf_true=1.0;
			float g1gGGprobPdf_true=1.0;
			float g4gZGprobPdf_true=1.0;
			float g4gGGprobPdf_true=1.0;

			float g1g2perpprobPdf_true=1.0;
			float g1g4perpprobPdf_true=1.0;

			for(int hypo=0;hypo<kNumSamples;hypo++){
				weight_probPdf=1.0;

				for(int gx=0;gx<size_zzcoupl;gx++){
					selfDHvvcoupl[gx][0]=0;
					selfDHvvcoupl[gx][1]=0;
				};

				selfDHvvcoupl[0][0] = gi_phi2_phi4[hypo][p_Re_g1]; // Re-g1
				selfDHvvcoupl[0][1] = gi_phi2_phi4[hypo][p_Im_g1]; // Im-g1
				selfDHvvcoupl[1][0] = gi_phi2_phi4[hypo][p_Re_g2]; // Re-g2
				selfDHvvcoupl[1][1] = gi_phi2_phi4[hypo][p_Im_g2]; // Im-g2
				selfDHvvcoupl[3][0] = gi_phi2_phi4[hypo][p_Re_g4]; // Re-g4
				selfDHvvcoupl[3][1] = gi_phi2_phi4[hypo][p_Im_g4]; // Im-g4
				selfDHvvcoupl[11][0] = gi_phi2_phi4[hypo][p_Re_L1];  // Re-g1prime2
				selfDHvvcoupl[11][1] = gi_phi2_phi4[hypo][p_Im_L1];  // Im-g1prime2
				selfDHvvcoupl[4][0] = gi_phi2_phi4[hypo][p_Re_a2ZG]; // Re-g_zgsg_2
				selfDHvvcoupl[4][1] = gi_phi2_phi4[hypo][p_Im_a2ZG]; // Im-g_zgsg_2
				selfDHvvcoupl[7][0] = gi_phi2_phi4[hypo][p_Re_a2GG]; // Re-g_gsgs_2
				selfDHvvcoupl[7][1] = gi_phi2_phi4[hypo][p_Im_a2GG]; // Im-g_gsgs_2
				selfDHvvcoupl[6][0] = gi_phi2_phi4[hypo][p_Re_a3ZG]; // Re-g_zgs_4
				selfDHvvcoupl[6][1] = gi_phi2_phi4[hypo][p_Im_a3ZG]; // Im-g_zgs_4
				selfDHvvcoupl[9][0] = gi_phi2_phi4[hypo][p_Re_a3GG]; // Re-g_gsgs_4
				selfDHvvcoupl[9][1] = gi_phi2_phi4[hypo][p_Im_a3GG]; // Im-g_gsgs_4
				selfDHvvcoupl[30][0] = gi_phi2_phi4[hypo][p_Re_ZGL1]; // Re-g_zgs_1prime2
				selfDHvvcoupl[30][1] = gi_phi2_phi4[hypo][p_Im_ZGL1]; // Im-g_zgs_1prime2
				selfDHvvcoupl[13][0] = gi_phi2_phi4[hypo][p_Re_LQ]; // Re-g1prime4 (g1-LambdaQ)
				selfDHvvcoupl[13][1] = gi_phi2_phi4[hypo][p_Im_LQ]; // Im-g1prime4 (g1-LambdaQ)

				if( genFinalState<=4 ){
					if (hypo<firstNonZZHypothesis || (GenZ1Mass>4 && GenZ2Mass > 4)){
						mela.setMelaHiggsWidth(gi_phi2_phi4[hypo][p_GH]); // MELA automatically resets width to default after calculation as part of protection for discriminants
						weight_probPdf = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
					}
					else weight_probPdf=0;
				}
				else{
					weight_probPdf = 1.0;
				};
				if(weight_probPdf==weight_probPdf){
					MC_weight_samples[hypo] = weight_probPdf/sample_probPdf;
//					testSpin0MEDivergence(smp,hypo,MC_weight_samples[hypo]);				
				}
				else{
					MC_weight_samples[hypo] = 0;
					isCountable=false;
				};
			};
			for(int hypo=0;hypo<kNumSamples;hypo++){
				if(genFinalState<=4 && isCountable){
					N_generated[genFinalState][hypo+1] += MC_weight_samples[hypo];
					if(GenZ1Mass>4 && GenZ2Mass>4) N_generated_4GeVcut[genFinalState][hypo+1] += MC_weight_samples[hypo];
				};
			};
			if(genFinalState<=4 && isCountable){
				N_generated[genFinalState][0] += 1.0;
				if(GenZ1Mass>4 && GenZ2Mass>4) N_generated_4GeVcut[genFinalState][0] += 1.0;
			};

			mytree->Fill();
		};
		for(int binx=0;binx<nFinalStates;binx++){
			for(int biny=0;biny<numSamples+1;biny++){
				hCount_new->SetBinContent(binx+1,biny+1,N_generated[binx][biny]);
				hCount_4GeVcut_new->SetBinContent(binx+1,biny+1,N_generated_4GeVcut[binx][biny]);
			};
		};
		foutput->WriteTObject(hCount_new);
		foutput->WriteTObject(hCount_4GeVcut_new);
		foutput->WriteTObject(mytree);

		foutput->Close();
		finput->Close();
	};
};

