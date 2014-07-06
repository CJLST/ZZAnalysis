#include <iostream>
#include <cmath>
#include <string>
#include "TChain.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "/afs/cern.ch/user/u/usarica/work/snowmass2013/src/ScalarPdfFactory_withFepspr.cc"
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

		RooRealVar* mzz_rrv = new RooRealVar("mzz","mzz",GenHMass,1.0e-20,1000.0);
		RooRealVar* mz1_rrv = new RooRealVar("mz1","mz1",1.0e-20,1000.0);
		RooRealVar* mz2_rrv = new RooRealVar("mz2","mz2",1.0e-20,1000.0);
		RooRealVar* hs_rrv = new RooRealVar("costhetastar","costhetastar",-1.,1.);
		RooRealVar* h1_rrv = new RooRealVar("costhetaZ1","costhetaZ1",-1.,1.);
		RooRealVar* h2_rrv = new RooRealVar("costhetaZ2","costhetaZ2",-1.,1.);
		RooRealVar* phi_rrv = new RooRealVar("phi","phi",-PI_VAL,PI_VAL);
		RooRealVar* phi1_rrv = new RooRealVar("phistarZ1","phistarZ1",-PI_VAL,PI_VAL);

		ScalarPdfFactory_withFepspr* sample_scalar = new ScalarPdfFactory_withFepspr(mz1_rrv,mz2_rrv,hs_rrv,h1_rrv,h2_rrv,phi_rrv,phi1_rrv,mzz_rrv,para,use_acc,true,mPOLE);
		ScalarPdfFactory_withFepspr* weight_scalar = new ScalarPdfFactory_withFepspr(mz1_rrv,mz2_rrv,hs_rrv,h1_rrv,h2_rrv,phi_rrv,phi1_rrv,mzz_rrv,para,use_acc,true,mPOLE);
		sample_scalar->_modelParams.fepspr->setVal(gi_phi2_phi4_files[smp][8]);
		sample_scalar->_modelParams.g1Val->setVal(gi_phi2_phi4_files[smp][0]);
		sample_scalar->_modelParams.g1ValIm->setVal(0);
		sample_scalar->_modelParams.g2Val->setVal(gi_phi2_phi4_files[smp][1]*cos(gi_phi2_phi4_files[smp][4]));
		sample_scalar->_modelParams.g2ValIm->setVal(gi_phi2_phi4_files[smp][1]*sin(gi_phi2_phi4_files[smp][4]));
		sample_scalar->_modelParams.g3Val->setVal(gi_phi2_phi4_files[smp][2]);
		sample_scalar->_modelParams.g3ValIm->setVal(0.0);
		sample_scalar->_modelParams.g4Val->setVal(gi_phi2_phi4_files[smp][3]*cos(gi_phi2_phi4_files[smp][5]));
		sample_scalar->_modelParams.g4ValIm->setVal(gi_phi2_phi4_files[smp][3]*sin(gi_phi2_phi4_files[smp][5]));

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
/*		mytree->Branch("Gen_D_g1Q2_phi0",&Gen_D_g1q2);
		mytree->Branch("Gen_D_g1_vs_g2_phi0",&Gen_D_g2);
		mytree->Branch("Gen_D_g1_vs_g4_phi0",&Gen_D_g4);
		mytree->Branch("Gen_D_g1Q2intPdf_phi0",&Gen_D_g1q2int);
		mytree->Branch("Gen_D_g2intPdf_phi0",&Gen_D_g2int);
		mytree->Branch("Gen_D_g4intPdf_phi0",&Gen_D_g4int);
		mytree->Branch("Gen_D_g2intPdf_phi90",&Gen_D_g2int_perp);
		mytree->Branch("Gen_D_g4intPdf_phi90",&Gen_D_g4int_perp);
*/		mytree->Branch("D_g1Q2_phi0",&Gen_D_g1q2);
		mytree->Branch("D_g1_vs_g2_phi0",&Gen_D_g2);
		mytree->Branch("D_g1_vs_g4_phi0",&Gen_D_g4);
		mytree->Branch("D_g1Q2intPdf_phi0",&Gen_D_g1q2int);
		mytree->Branch("D_g2intPdf_phi0",&Gen_D_g2int);
		mytree->Branch("D_g4intPdf_phi0",&Gen_D_g4int);
		mytree->Branch("D_g2intPdf_phi90",&Gen_D_g2int_perp);
		mytree->Branch("D_g4intPdf_phi90",&Gen_D_g4int_perp);

		float N_generated[nFinalStates][kNumSamples+1];
		for(int xb=0;xb<kNumSamples+1;xb++){ for(int yb=0;yb<nFinalStates;yb++) N_generated[yb][xb]=0;};

		int nEntries = tree->GetEntries();
		cout << "Entries " << nEntries << endl;
		for(int ev=0;ev<nEntries;ev++){
			tree->GetEntry(ev);

			mzz_rrv->setVal(mPOLE);
			mz1_rrv->setVal(GenZ1Mass);
			mz2_rrv->setVal(GenZ2Mass);
			hs_rrv->setVal(myGencosthetastar);
			h1_rrv->setVal(myGenhelcosthetaZ1);
			h2_rrv->setVal(myGenhelcosthetaZ2);
			phi_rrv->setVal(myGenhelphi);
			phi1_rrv->setVal(myGenphistarZ1);

			if(genFinalState<=4) sample_probPdf = sample_scalar->PDF->getVal();
			else sample_probPdf=1.0;
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

				weight_scalar->_modelParams.fepspr->setVal(gi_phi2_phi4[hypo][8]);
				weight_scalar->_modelParams.g1Val->setVal(gi_phi2_phi4[hypo][0]);
				weight_scalar->_modelParams.g1ValIm->setVal(0);
				weight_scalar->_modelParams.g2Val->setVal(gi_phi2_phi4[hypo][1]*cos(gi_phi2_phi4[hypo][4]));
				weight_scalar->_modelParams.g2ValIm->setVal(gi_phi2_phi4[hypo][1]*sin(gi_phi2_phi4[hypo][4]));
				weight_scalar->_modelParams.g3Val->setVal(gi_phi2_phi4[hypo][2]);
				weight_scalar->_modelParams.g3ValIm->setVal(0.0);
				weight_scalar->_modelParams.g4Val->setVal(gi_phi2_phi4[hypo][3]*cos(gi_phi2_phi4[hypo][5]));
				weight_scalar->_modelParams.g4ValIm->setVal(gi_phi2_phi4[hypo][3]*sin(gi_phi2_phi4[hypo][5]));

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
						weight_probPdf = weight_scalar->PDF->getVal();

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
			Gen_D_g1q2int = ( g1g1q2probPdf_true - g1probPdf_true*pow( ( 1.0-abs(gi_phi2_phi4[7][8]) ),2.0) - g1q2probPdf_true*pow( abs(gi_phi2_phi4[7][8]) ,2.0) )*(g1q2intscale/( ( 1.0-abs(gi_phi2_phi4[7][8]) )*abs(gi_phi2_phi4[7][8]) ) ) / (g1probPdf_true + g1q2probPdf_true*g1q2scale);

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
			hkDs_true[ht]->SetOption("colz");
			foutput->WriteTObject(hkDs_true[ht]);
			delete hkDs_true[ht];
		};
		foutput->WriteTObject(hCount_new);
		foutput->WriteTObject(mytree);

		foutput->Close();
		finput->Close();
	};
};
