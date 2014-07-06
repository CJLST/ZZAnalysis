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
//#include "/afs/cern.ch/user/u/usarica/work/snowmass2013/src/ScalarPdfFactory_withFepspr.cc"
#include <ZZMatrixElement/MELA/interface/Mela.h>
#include <ZZMatrixElement/MELA/src/computeAngles.h>
#include "./data/ZZ4l_126_Samples.h"
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


float getJHUGenMELAWeight(Mela& myMela, int lepId[4], float angularOrdered[8], double selfDHvvcoupl[20][2]){
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

void reweightIdeal(int flavor, int erg_tev, float mPOLE=126.0){
	char GenLevel_Location[]="/GenSignal";
	char TREE_NAME[] = "GenTree";
	char erg_dir[1000];
	sprintf(erg_dir,"LHC_%iTeV",erg_tev);
	TString cinput_common = user_dir + erg_dir;
	TString coutput_common = cinput_common;

	for(int smp=0; smp<kNumFiles; smp++){
		if(flavor>2) continue;

		TString cinput = cinput_common + GenLevel_Location + "/" + sample_file[smp] + "_Generated.root";
		TString coutput = coutput_common + "/" + user_folder[flavor] + "/" + sample_file[smp] + "_GenLevel_ReWeighed.root";

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

		Mela mela(erg_tev,mPOLE);
		mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);

		float NGen_weighted=0;
		float MC_weight_samples[kNumSamples];
		int numSamples = kNumSamples;
		for(int s=0;s<kNumSamples;s++) MC_weight_samples[s]=1;
		TH2F* hCount_new = new TH2F("hCounters_spin0_RW","Counters with Spin0 Re-weights",nFinalStates,0,nFinalStates,numSamples+1,0,numSamples+1);
		TH2F* hkDs_true[8] = {
			new TH2F("D_g1Q2_phi0_True","",30,0.0,1.0,numSamples+1,0,numSamples+1),
			new TH2F("D_g1_vs_g2_phi0_True","",30,0.0,1.0,numSamples+1,0,numSamples+1),
			new TH2F("D_g1_vs_g4_phi0_True","",30,0.0,1.0,numSamples+1,0,numSamples+1),
			new TH2F("D_g1Q2intPdf_phi0_True","",30,-1.0,1.0,numSamples+1,0,numSamples+1),
			new TH2F("D_g2intPdf_phi0_True","",30,-1.0,1.0,numSamples+1,0,numSamples+1),
			new TH2F("D_g4intPdf_phi0_True","",30,-1.0,1.0,numSamples+1,0,numSamples+1),
			new TH2F("D_g2intPdf_phi90_True","",30,-1.0,1.0,numSamples+1,0,numSamples+1),
			new TH2F("D_g4intPdf_phi90_True","",30,-1.0,1.0,numSamples+1,0,numSamples+1)
		};
		for(int ht=0;ht<8;ht++){
			for(int binx=1;binx<=30;binx++){
				for(int biny=1;biny<=numSamples+1;biny++) hkDs_true[ht]->SetBinContent(binx,biny,0);
			};
		};

		foutput->cd();
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
		float Gen_D_g1q2=-99,Gen_D_g1q2int=-99;
		float Gen_D_g2=-99,Gen_D_g2int=-99,Gen_D_g2int_perp=-99;
		float Gen_D_g4=-99,Gen_D_g4int=-99,Gen_D_g4int_perp=-99;
		mytree->Branch("D_g1Q2_phi0",&Gen_D_g1q2);
		mytree->Branch("D_g1_vs_g2_phi0",&Gen_D_g2);
		mytree->Branch("D_g1_vs_g4_phi0",&Gen_D_g4);
		mytree->Branch("D_g1Q2intPdf_phi0",&Gen_D_g1q2int);
		mytree->Branch("D_g2intPdf_phi0",&Gen_D_g2int);
		mytree->Branch("D_g4intPdf_phi0",&Gen_D_g4int);
		mytree->Branch("D_g2intPdf_phi90",&Gen_D_g2int_perp);
		mytree->Branch("D_g4intPdf_phi90",&Gen_D_g4int_perp);

		float N_generated[nFinalStates][kNumSamples+1];
		for(int xb=0;xb<kNumSamples+1;xb++){ for(int yb=0;yb<nFinalStates;yb++) N_generated[yb][xb]=0;};

		double selfDHvvcoupl[20][2];
		for(int gx=0;gx<20;gx++){
			selfDHvvcoupl[gx][0]=0;
			selfDHvvcoupl[gx][1]=0;
		};

		int nEntries = tree->GetEntries();
		cout << "Entries " << nEntries << endl;
		for(int ev=0;ev<nEntries;ev++){
			tree->GetEntry(ev);
			sample_probPdf=1.0;

			if((ev%10000)==0) cout << "Event: " <<  ev << endl;

			int lepIdOrdered[4]={ GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id };
			float angularOrdered[8]={GenHMass,GenZ1Mass,GenZ2Mass,myGencosthetastar,myGenhelcosthetaZ1,myGenhelcosthetaZ2,myGenhelphi,myGenphistarZ1};

			for(int gx=0;gx<20;gx++){
				selfDHvvcoupl[gx][0]=0;
				selfDHvvcoupl[gx][1]=0;
			};
			selfDHvvcoupl[0][0] = (gi_phi2_phi4_files[smp][0]) * (1.0-abs( gi_phi2_phi4_files[smp][8] ) );
			selfDHvvcoupl[1][0] = (gi_phi2_phi4_files[smp][1]) * cos( gi_phi2_phi4_files[smp][4] );
			selfDHvvcoupl[1][1] = (gi_phi2_phi4_files[smp][1]) * sin( gi_phi2_phi4_files[smp][4] );
			selfDHvvcoupl[3][0] = (gi_phi2_phi4_files[smp][3]) * cos( gi_phi2_phi4_files[smp][5] );
			selfDHvvcoupl[3][1] = (gi_phi2_phi4_files[smp][3]) * sin( gi_phi2_phi4_files[smp][5] );
			selfDHvvcoupl[5][0] = ( gi_phi2_phi4_files[smp][8]);
			sample_probPdf = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
			if(sample_probPdf==0) sample_probPdf=1.0;

			float weight_probPdf=1.0;

			Gen_D_g2 = -99;
			Gen_D_g4 = -99;
			Gen_D_g1q2 = -99;
			Gen_D_g2int = -99;
			Gen_D_g4int = -99;
			Gen_D_g1q2int = -99;
			Gen_D_g2int_perp = -99;
			Gen_D_g4int_perp = -99;

			float g1probPdf_true=1.0;
			float g2probPdf_true=1.0;
			float g4probPdf_true=1.0;
			float g1q2probPdf_true=1.0;

			float g1g2probPdf_true=1.0;
			float g1g4probPdf_true=1.0;
			float g1g1q2probPdf_true=1.0;

			float g1g2perpprobPdf_true=1.0;
			float g1g4perpprobPdf_true=1.0;

			for(int hypo=0;hypo<kNumSamples;hypo++){
				weight_probPdf=1.0;

				for(int gx=0;gx<20;gx++){
					selfDHvvcoupl[gx][0]=0;
					selfDHvvcoupl[gx][1]=0;
				};
				selfDHvvcoupl[0][0] = (gi_phi2_phi4[hypo][0]);
				selfDHvvcoupl[1][0] = (gi_phi2_phi4[hypo][1]) * cos( gi_phi2_phi4[hypo][4] );
				selfDHvvcoupl[1][1] = (gi_phi2_phi4[hypo][1]) * sin( gi_phi2_phi4[hypo][4] );
				selfDHvvcoupl[3][0] = (gi_phi2_phi4[hypo][3]) * cos( gi_phi2_phi4[hypo][5] );
				selfDHvvcoupl[3][1] = (gi_phi2_phi4[hypo][3]) * sin( gi_phi2_phi4[hypo][5] );
				selfDHvvcoupl[5][0] = (gi_phi2_phi4[hypo][8]);

				if( (genFinalState<=4 )
					|| hypo==0
					|| hypo==1
					|| hypo==2
					|| hypo==3
					|| hypo==4
					|| hypo==5
					|| hypo==7
					|| hypo==13
					|| hypo==14
					){
						weight_probPdf = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

						if(hypo==0) g1probPdf_true = weight_probPdf;
						if(hypo==1) g2probPdf_true = weight_probPdf;
						if(hypo==2) g4probPdf_true = weight_probPdf;
						if(hypo==3) g1q2probPdf_true = weight_probPdf;

						if(hypo==4) g1g2probPdf_true = weight_probPdf;
						if(hypo==5) g1g4probPdf_true = weight_probPdf;
						if(hypo==7) g1g1q2probPdf_true = weight_probPdf;

						if(hypo==13) g1g2perpprobPdf_true = weight_probPdf;
						if(hypo==14) g1g4perpprobPdf_true = weight_probPdf;
				};

				if( genFinalState>4 ) weight_probPdf = 1.0;
				MC_weight_samples[hypo] = weight_probPdf/sample_probPdf;
				if(genFinalState<=4) N_generated[genFinalState][hypo+1] += MC_weight_samples[hypo];
			};
			if(genFinalState<=4) N_generated[genFinalState][0] += 1.0;

			Gen_D_g2 = g1probPdf_true/(g1probPdf_true + g2probPdf_true*g1g2scale_old);
			Gen_D_g4 = g1probPdf_true/(g1probPdf_true + g4probPdf_true*g1g4scale_old);
			Gen_D_g1q2 = g1probPdf_true/(g1probPdf_true + g1q2probPdf_true*g1q2scale );

			Gen_D_g2int = ( g1g2probPdf_true - g1probPdf_true - g2probPdf_true*pow(gi_phi2_phi4[4][1],2.0) ) * ( g1g2intscale/gi_phi2_phi4[4][1] ) / (g1probPdf_true + g2probPdf_true*g1g2scale_old);
			Gen_D_g4int = ( g1g4probPdf_true - g1probPdf_true - g4probPdf_true*pow(gi_phi2_phi4[5][3],2.0) ) * ( g1g4intscale/gi_phi2_phi4[5][3] ) / (g1probPdf_true + g4probPdf_true*g1g4scale_old);
			Gen_D_g1q2int = ( g1g1q2probPdf_true - g1probPdf_true - g1q2probPdf_true ) / (g1probPdf_true + g1q2probPdf_true);

			Gen_D_g2int_perp = ( g1g2perpprobPdf_true - g1probPdf_true - g2probPdf_true*pow(gi_phi2_phi4[13][1],2.0) ) * ( g1g2intscale/gi_phi2_phi4[13][1] ) / (g1probPdf_true + g2probPdf_true*g1g2scale_old);
			Gen_D_g4int_perp = ( g1g4perpprobPdf_true - g1probPdf_true - g4probPdf_true*pow(gi_phi2_phi4[14][3],2.0) ) * ( g1g4intscale/gi_phi2_phi4[14][3] ) / (g1probPdf_true + g4probPdf_true*g1g4scale_old);

			if((ev%100000)==0) cout << Gen_D_g1q2 << '\t' << g1g2probPdf_true << '\t' << g1probPdf_true << '\t' << g2probPdf_true << endl;

			if(flavor!=genFinalState) continue;

			for(int s=0;s<kNumSamples+1;s++){
				for(int ht=0;ht<8;ht++){
					float someKD;
					if(ht==0) someKD=Gen_D_g1q2;
					if(ht==1) someKD=Gen_D_g2;
					if(ht==2) someKD=Gen_D_g4;
					if(ht==3) someKD=Gen_D_g1q2int;
					if(ht==4) someKD=Gen_D_g2int;
					if(ht==5) someKD=Gen_D_g4int;
					if(ht==6) someKD=Gen_D_g2int_perp;
					if(ht==7) someKD=Gen_D_g4int_perp;

					int KDbin=hkDs_true[ht]->GetXaxis()->FindBin(someKD);
					float bincontent = hkDs_true[ht]->GetBinContent(KDbin,s+1);
					if(s==0) bincontent += 1.0;
					else bincontent += MC_weight_samples[s-1];
					hkDs_true[ht]->SetBinContent(KDbin,s+1,bincontent);
				};
			};

			mytree->Fill();
		};
		for(int binx=0;binx<nFinalStates;binx++){
			for(int biny=0;biny<kNumSamples+1;biny++) hCount_new->SetBinContent(binx+1,biny+1,N_generated[binx][biny]);
		};
		for(int ht=0;ht<8;ht++){
//			hkDs_true[ht]->SetOption("colz");
			foutput->WriteTObject(hkDs_true[ht]);
			delete hkDs_true[ht];
		};
		foutput->WriteTObject(hCount_new);
		foutput->WriteTObject(mytree);

		foutput->Close();
		finput->Close();
	};
};

