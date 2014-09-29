#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include "TMath.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include <ZZMatrixElement/MELA/interface/Mela.h>
#include "./data/HZZ4l_LeptonInterference.h"

string user_dir="/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/";
string user_dir_newProduction="/afs/cern.ch/work/l/lfeng/public_write/usarica/Production/";
const int kAllSamples=79;
const int kGGHSamples=35;
const int kGGMCFMSamples=49;
const int kGGSamples=55;
const int kQQBZZSamples=61;
const int kQQBZZSamples_Dedicated=66;
const int kVBF_Phantom=78;
char* sampleName[kAllSamples]={
	"HZZ4lTree_powheg15jhuGenV3H126",
	"HZZ4lTree_powheg15jhuGenV3-0MH126",
	"HZZ4lTree_powheg15jhuGenV3-0PHH126",
	"HZZ4lTree_powheg15jhuGenV3-0Mf01ph0H126",
	"HZZ4lTree_powheg15jhuGenV3-0Mf01ph180H126",
	"HZZ4lTree_powheg15jhuGenV3-0Mf01ph270H126",
	"HZZ4lTree_powheg15jhuGenV3-0Mf01ph90H126",
	"HZZ4lTree_powheg15jhuGenV3-0Mf05ph0H126",
	"HZZ4lTree_powheg15jhuGenV3-0Mf05ph180H126",
	"HZZ4lTree_powheg15jhuGenV3-0Mf05ph270H126",
	"HZZ4lTree_powheg15jhuGenV3-0Mf05ph90H126",
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
	"HZZ4lTree_powheg15jhuGenV3-0PHf05ph180Mf05ph0H125.6",

	"HZZ4lTree_ggZZ4l",
	"HZZ4lTree_ggZZ2l2l",

	"HZZ4lTree_ggTo4mu_Contin-MCFM67",
	"HZZ4lTree_ggTo4e_Contin-MCFM67",
	"HZZ4lTree_ggTo2e2mu_Contin-MCFM67",

	"HZZ4lTree_ggTo4mu_SMH-MCFM67",
	"HZZ4lTree_ggTo4e_SMH-MCFM67",
	"HZZ4lTree_ggTo2e2mu_SMH-MCFM67",

	"HZZ4lTree_ggTo4mu_SMHContinInterf-MCFM67",
	"HZZ4lTree_ggTo4e_SMHContinInterf-MCFM67",
	"HZZ4lTree_ggTo2e2mu_SMHContinInterf-MCFM67",

	"HZZ4lTree_ggTo4mu_BSMHContinInterf-MCFM67",
	"HZZ4lTree_ggTo4e_BSMHContinInterf-MCFM67",
	"HZZ4lTree_ggTo2e2mu_BSMHContinInterf-MCFM67",

	"HZZ4lTree_ggTo2l2l_H125.6",
	"HZZ4lTree_ggTo2l2l_Continuum",
	"HZZ4lTree_ggTo2l2l_ContinuumInterfH125.6",
	"HZZ4lTree_ggTo4l_H125.6",
	"HZZ4lTree_ggTo4l_Continuum",
	"HZZ4lTree_ggTo4l_ContinuumInterfH125.6",
		"HZZ4lTree_ZZTo4mu",
		"HZZ4lTree_ZZTo4e",
		"HZZ4lTree_ZZTo2e2mu",
		"HZZ4lTree_ZZTo2mu2tau",
		"HZZ4lTree_ZZTo2e2tau",
		"HZZ4lTree_ZZTo4tau",
"HZZ4lTree_ZZ95-160To2e2mu",
"HZZ4lTree_ZZ95-160To2mu2tau",
"HZZ4lTree_ZZ95-160To4e",
"HZZ4lTree_ZZ95-160To4mu",
"HZZ4lTree_ZZ95-160To4tau",

	"HZZ4lTree_ZZTo4muJJ_Contin",
	"HZZ4lTree_ZZTo4eJJ_Contin",
	"HZZ4lTree_ZZTo2e2muJJ_Contin",

	"HZZ4lTree_ZZTo4muJJ_SMHContinInterf_H125.6",
	"HZZ4lTree_ZZTo4eJJ_SMHContinInterf_H125.6",
	"HZZ4lTree_ZZTo2e2muJJ_SMHContinInterf_H125.6",

	"HZZ4lTree_ZZTo4muJJ_BSM10HContinInterf_H125.6",
	"HZZ4lTree_ZZTo4eJJ_BSM10HContinInterf_H125.6",
	"HZZ4lTree_ZZTo2e2muJJ_BSM10HContinInterf_H125.6",

	"HZZ4lTree_ZZTo4muJJ_BSM25HContinInterf_H125.6",
	"HZZ4lTree_ZZTo4eJJ_BSM25HContinInterf_H125.6",
	"HZZ4lTree_ZZTo2e2muJJ_BSM25HContinInterf_H125.6",

	"HZZ4lTree_DoubleOr_CRZLLTree"
};
char* data_files[3]={
	"HZZ4lTree_DoubleMu",
	"HZZ4lTree_DoubleEle",
	"HZZ4lTree_DoubleOr"
};
char* user_folder[5]={
	"4mu",
	"4e",
	"2mu2e",
	"CR",
	"data"
};




void setCTotalBkgGraphs(int folder, TFile* fcontainer, const int tgCSize, TGraph* tgC[]){
	string tgname = "C_TotalBkgM4l_";

	float rValues[6]={1,5,10,15,20,25};

	for(int mR=0;mR<tgCSize;mR++){

		float myWidth = 1;

		char crValue[20];

		int rCode = mR % 6;
		if(mR<6){
			myWidth = rValues[rCode];
			sprintf(crValue,"D_Gamma_r%.0f",myWidth);
		}
		else if(mR<12){
			myWidth = rValues[rCode];
			sprintf(crValue,"D_Gamma_gg_r%.0f",myWidth);
		}
		else if(mR==12) sprintf(crValue,"D_Gamma");
		else if(mR==13) sprintf(crValue,"D_Gamma_int");

		string ctgM4L = tgname;
		ctgM4L = ctgM4L + user_folder[folder] + "_";
		ctgM4L = ctgM4L + crValue;
		tgC[mR] = (TGraph*) fcontainer->Get(ctgM4L.c_str());
	};
};
void getCTotalBkg (float MZZ, const int tgCSize, float myCTotalBkg[], TGraph* tgC[]){
	for(int mR=0;mR<tgCSize;mR++) myCTotalBkg[mR] = tgC[mR]->Eval(MZZ);
};
void constructVariousDGamma(float c_ggzz, float bkg_VAMCFM_noscale, float ggzz_VAMCFM_noscale, float ggHZZ_prob_pure_noscale, float ggHZZ_prob_int_noscale, const int tgCSize, float myCTotalBkg[], float myDGamma[]){
	float rValues[6]={1,5,10,15,20,25};
	for(int mR=0;mR<tgCSize;mR++){
		float total_sig_ME;
		float total_bkg_ME;
		float ctotal_bkg = myCTotalBkg[mR];
		float myWidth = 1;
		int rCode = mR % 6;
		if(mR<12) myWidth = rValues[rCode];

		if(mR<6){
			total_sig_ME = (myWidth * ggHZZ_prob_pure_noscale + sqrt(myWidth) * ggHZZ_prob_int_noscale);
			total_bkg_ME = (c_ggzz*ggzz_VAMCFM_noscale + bkg_VAMCFM_noscale)*ctotal_bkg;
		}
		else if(mR<12){
			total_sig_ME = (myWidth * ggHZZ_prob_pure_noscale + sqrt(myWidth) * ggHZZ_prob_int_noscale + ggzz_VAMCFM_noscale);
			total_bkg_ME = bkg_VAMCFM_noscale*ctotal_bkg;
		}
		else if(mR==12){
			total_sig_ME = ggHZZ_prob_pure_noscale;
			total_bkg_ME = (c_ggzz*ggzz_VAMCFM_noscale + bkg_VAMCFM_noscale)*ctotal_bkg;
		}
		else if(mR==13){
			total_sig_ME = ggHZZ_prob_int_noscale;
			total_bkg_ME = (c_ggzz*ggzz_VAMCFM_noscale + bkg_VAMCFM_noscale)*ctotal_bkg;
		};
		float kd_denominator = (total_sig_ME+total_bkg_ME);
		if(mR==13) kd_denominator = (ggHZZ_prob_pure_noscale + total_bkg_ME);

		float kd = total_sig_ME/kd_denominator;
		if(mR<6){
			if(kd<=1.0) kd = exp(kd-1.0);
			else kd = (2.0-exp(1.0-kd));
			kd -= 1.0;
		};
		myDGamma[mR] = kd;
	};
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
float evaluatePDFUncertainty(float mzz, TGraph* tg[], int systematic=0, float mPOLE = 125.6){
	float dPDF=0;
/*
	int tgIndex=-1;
	if(systematic==1) tgIndex=1;
	else if(systematic==-1) tgIndex=0;
	if(tgIndex>=0 && tgIndex<2) dPDF = tg[tgIndex]->Eval(mzz);
*/
//	if(systematic==1) dPDF = 1.91007e-5*(mzz-mPOLE);
//	if(systematic==-1) dPDF = -1.91007e-5*(mzz-mPOLE);
	if(systematic==1) dPDF = 0.01751 - 7.514738e-5*mzz;
	if(systematic==-1) dPDF = -0.01751 + 7.514738e-5*mzz;

	return dPDF;
};
float getMCFMMELAWeight(Mela& myMela, int lepId[4], float angularOrdered[8], bool isWeighted=false){
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
	    myprob,
		isWeighted
		);
	return myprob;
};
TGraph* make_HZZ_LeptonInterferenceGraph(){
//void make_HZZ_LeptonInterferenceGraph={
	string cinput="./data/HZZ_LeptonInterferenceGraph.root";
	float x[leptonInterf_YR3_Size];
	float y[leptonInterf_YR3_Size];
	for(int a=0;a<leptonInterf_YR3_Size;a++){
		x[a] = leptonInterf_YR3[a][0];
		y[a] = leptonInterf_YR3[a][1];
	};
	TGraph* tg = new TGraph(leptonInterf_YR3_Size,x,y);
	tg->SetName("tgHZZ_LeptonInterference");
	tg->SetTitle("H#rightarrowZZ#rightarrow4l Lepton Interference Weight on 4e, 4#mu wrt. 2e2#mu");

//	TFile* foutput = new TFile(cinput.c_str(),"recreate");
//	foutput->WriteTObject(tg);
//	foutput->Close();
//	delete tg;
	return tg;
};

void make_QQBGG_Counters(int erg_tev, int folder){
	char OUTPUT_NAME[]="HZZ4l_QQBGGSampleRecoCounter";
	char TREE_NAME[] = "SelectedTree";
	char erg_dir[1000];
	sprintf(erg_dir,"LHC_%iTeV",erg_tev);
	char cerg[1000];
	sprintf(cerg,"%iTeV",erg_tev);
	string cinput_common = user_dir;
	cinput_common = cinput_common + erg_dir;
//	string cinput_common = user_dir_newProduction;
//	cinput_common = cinput_common + cerg;
	string coutput_common = cinput_common;

	string coutput_main = coutput_common;
	coutput_main = coutput_main + "/" + user_folder[folder] + "/";
	char coutput[1000];
	sprintf(coutput,"%s%s%s",coutput_main.c_str(),OUTPUT_NAME,".root");
	TFile* foutput = new TFile(coutput,"recreate");
	const int kNumGGTypes=3;
	const int kNumQQTypes=2;
	const int nbinsx = (kNumGGTypes+kNumQQTypes);
	TH1F* heff = new TH1F("hQQBGG_EfficiencyMatcher","Match gg-QQ Weights: Bins are ggZZ, qqZZ and dedicated qqZZ x ggZZ, qqZZ",2*nbinsx+2,0,2*nbinsx+2);
	TChain* tree[nbinsx];
	for(int t=0;t<nbinsx;t++) tree[t] = new TChain(TREE_NAME);
	double nReco[2*nbinsx] = {0};
	double sumMCwgt[2*nbinsx] = {0};
	double sumMCwgt_noxsec[2*nbinsx] = {0};
	double sumMCwgtMEratio[2*nbinsx] = {0};
	double sumMCwgt_noxsecMEratio[2*nbinsx] = {0};
	double avgXSEC[2*nbinsx]={0};
	double avgXSEC_weighted[2*nbinsx]={0};
	double chosenXSEC[2]={0};

//	float DedicatedQQBSample_GenZZcuts[2]={105,155};
	float DedicatedQQBSample_GenZZcuts[2]={105,141};
	float GenHMass=1;
	float ZZMass=1;
	float MC_weight=1;
	float MC_weight_noxsec=1;
	float MC_weight_QQBGG_VAMCFM=1;

	for(int smp=kGGHSamples; smp<kQQBZZSamples_Dedicated; smp++){
		if(smp>=kGGHSamples+5 && smp<kGGSamples && !(smp==kGGHSamples+15 || smp==kGGHSamples+18)) continue;
//		if(folder!=2 && smp==kGGHSamples+18) continue; // ggTo4l fix
		cout << "Sample: " << smp << endl;

		string cinput_main = cinput_common;
		cinput_main = cinput_main + "/" + user_folder[folder] + "/";

		char cinput[1000];
		sprintf(cinput,"%s%s%s",cinput_main.c_str(),sampleName[smp],".root");
		cout << cinput << endl;

		if(smp<kGGHSamples+2) tree[0]->Add(cinput);
		else if(smp<kGGHSamples+5) tree[1]->Add(cinput);
		else if(smp<kGGSamples) tree[2]->Add(cinput);
		else if(smp<kQQBZZSamples) tree[3]->Add(cinput);
		else if(smp<kQQBZZSamples_Dedicated) tree[4]->Add(cinput);
	};
	for(int t=0;t<nbinsx;t++){
		tree[t]->SetBranchAddress("GenHMass",&GenHMass);
		tree[t]->SetBranchAddress("ZZMass",&ZZMass);
		tree[t]->SetBranchAddress("MC_weight",&MC_weight);
		tree[t]->SetBranchAddress("MC_weight_noxsec",&MC_weight_noxsec);
		tree[t]->SetBranchAddress("MC_weight_QQBGG_VAMCFM",&MC_weight_QQBGG_VAMCFM);

		int nEntries = tree[t]->GetEntries();
		for(int ev=0;ev<nEntries;ev++){
			tree[t]->GetEntry(ev);

			double MC_weight_alphaS = evaluateAlphaSShift(GenHMass, 1, 0,125.6);

			nReco[2*t] += 1;
			if(t==0 || t==2){
				sumMCwgt[2*t] += MC_weight*MC_weight_alphaS;
				sumMCwgtMEratio[2*t] += MC_weight*MC_weight_QQBGG_VAMCFM*MC_weight_alphaS;
				sumMCwgt_noxsec[2*t] += MC_weight_noxsec*MC_weight_alphaS;
				sumMCwgt_noxsecMEratio[2*t] += MC_weight_noxsec*MC_weight_QQBGG_VAMCFM*MC_weight_alphaS;
			}
			else{
				sumMCwgt[2*t] += MC_weight;
				sumMCwgtMEratio[2*t] += MC_weight*MC_weight_QQBGG_VAMCFM;
				sumMCwgt_noxsec[2*t] += MC_weight_noxsec;
				sumMCwgt_noxsecMEratio[2*t] += MC_weight_noxsec*MC_weight_QQBGG_VAMCFM;
			};
//			if( GenHMass<DedicatedQQBSample_GenZZcuts[1] && GenHMass>=DedicatedQQBSample_GenZZcuts[0] && ZZMass<DedicatedQQBSample_GenZZcuts[1] && ZZMass>=DedicatedQQBSample_GenZZcuts[0] ){
			if( ZZMass<DedicatedQQBSample_GenZZcuts[1] && ZZMass>=DedicatedQQBSample_GenZZcuts[0] ){
				nReco[2*t+1] += 1;
				if(t==0 || t==2){
					sumMCwgt[2*t+1] += MC_weight*MC_weight_alphaS;
					sumMCwgt_noxsec[2*t+1] += MC_weight_noxsec*MC_weight_alphaS;
					sumMCwgtMEratio[2*t+1] += MC_weight*MC_weight_QQBGG_VAMCFM*MC_weight_alphaS;
					sumMCwgt_noxsecMEratio[2*t+1] += MC_weight_noxsec*MC_weight_QQBGG_VAMCFM*MC_weight_alphaS;
				}
				else{
					sumMCwgt[2*t+1] += MC_weight;
					sumMCwgt_noxsec[2*t+1] += MC_weight_noxsec;
					sumMCwgtMEratio[2*t+1] += MC_weight*MC_weight_QQBGG_VAMCFM;
					sumMCwgt_noxsecMEratio[2*t+1] += MC_weight_noxsec*MC_weight_QQBGG_VAMCFM;
				};
			};
		};
		for(int m=0;m<2;m++){
			avgXSEC[2*t+m] = sumMCwgt[2*t+m] / sumMCwgt_noxsec[2*t+m];
			avgXSEC_weighted[2*t+m] = sumMCwgtMEratio[2*t+m] / sumMCwgt_noxsec[2*t+m];
		};
	};
	for(int t=0;t<nbinsx;t++) delete tree[t];
	int chosenGGIndex = 0;
	int chosenQQIndex = 0;
	for(int t=0;t<kNumGGTypes;t++){
		if(nReco[2*t+1]>chosenXSEC[0]){
			chosenXSEC[0]=nReco[2*t+1];
			chosenGGIndex=t;
		};
	};
	for(int t=kNumGGTypes;t<nbinsx;t++){
		if(nReco[2*t+1]>chosenXSEC[1]){
			chosenXSEC[1]=nReco[2*t+1];
			chosenQQIndex=t;
		};
	};
	chosenXSEC[0] = avgXSEC[2*chosenGGIndex+1];
	chosenXSEC[1] = avgXSEC[2*chosenQQIndex+1];
	heff->SetBinContent(2*nbinsx+1,chosenXSEC[0]);
	heff->SetBinContent(2*nbinsx+2,chosenXSEC[1]);

	cout << "ggZZ chosen: " << chosenGGIndex << " at " << chosenXSEC[0] << endl;
	cout << "qqZZ chosen: " << chosenQQIndex << " at " << chosenXSEC[1] << endl;

	double totalGenerated_Dedicated = nReco[1] + nReco[3] + nReco[5] + nReco[7] + nReco[9];

	heff->SetBinContent(1, (sumMCwgt_noxsec[2*chosenGGIndex+1]/sumMCwgt_noxsec[1])*nReco[1]/( totalGenerated_Dedicated ) );
	heff->SetBinContent(2, (sumMCwgt_noxsec[2*chosenGGIndex+1]/sumMCwgt_noxsec[3])*nReco[3]/( totalGenerated_Dedicated ) );
	heff->SetBinContent(3, (sumMCwgt_noxsec[2*chosenGGIndex+1]/sumMCwgt_noxsec[5])*nReco[5]/( totalGenerated_Dedicated ) );
	heff->SetBinContent(4, (sumMCwgt_noxsec[2*chosenGGIndex+1]/sumMCwgt_noxsecMEratio[7])*nReco[7]/( totalGenerated_Dedicated ) );
	heff->SetBinContent(5, (sumMCwgt_noxsec[2*chosenGGIndex+1]/sumMCwgt_noxsecMEratio[9])*nReco[9]/( totalGenerated_Dedicated ) );

	heff->SetBinContent(6, (sumMCwgt_noxsec[2*chosenQQIndex+1]/sumMCwgt_noxsecMEratio[1])*nReco[1]/( totalGenerated_Dedicated ) );
	heff->SetBinContent(7, (sumMCwgt_noxsec[2*chosenQQIndex+1]/sumMCwgt_noxsecMEratio[3])*nReco[3]/( totalGenerated_Dedicated ) );
	heff->SetBinContent(8, (sumMCwgt_noxsec[2*chosenQQIndex+1]/sumMCwgt_noxsecMEratio[5])*nReco[5]/( totalGenerated_Dedicated ) );
	heff->SetBinContent(9, (sumMCwgt_noxsec[2*chosenQQIndex+1]/sumMCwgt_noxsec[7])*nReco[7]/( totalGenerated_Dedicated ) );
	heff->SetBinContent(10, (sumMCwgt_noxsec[2*chosenQQIndex+1]/sumMCwgt_noxsec[9])*nReco[9]/( totalGenerated_Dedicated ) );

	for(int b=1;b<=2*nbinsx+2;b++) cout << heff->GetBinContent(b) << endl;
	double sum=0;
	sum += heff->GetBinContent(1)*sumMCwgt_noxsec[1]*chosenXSEC[0];
	sum += heff->GetBinContent(2)*sumMCwgt_noxsec[3]*chosenXSEC[0];
	sum += heff->GetBinContent(3)*sumMCwgt_noxsec[5]*chosenXSEC[0];
	sum += heff->GetBinContent(4)*sumMCwgt_noxsecMEratio[7]*chosenXSEC[0];
	sum += heff->GetBinContent(5)*sumMCwgt_noxsecMEratio[9]*chosenXSEC[0];
	cout << sum << endl;
	foutput->WriteTObject(heff);
	foutput->Close();
};
double getQQZZEWKCorrection(float GenHMass, int systematics = 0){
	double A = -0.0782143;
	double MA = 125.369;
	double EWKcorr = A*log(GenHMass / MA);
	double m4lThreshold = 2.0*(91.1876 - 2.4952*2.0);

	double B = 1.89547;
	double MB = 68.1253;
	double OB = 1.53646;
	double QCDcorr = OB - B*exp(-GenHMass / MB);

	if (EWKcorr < -1) EWKcorr = -1;
	if (QCDcorr < 0) QCDcorr = 0;
	if (EWKcorr>1.0){
		EWKcorr = 1.0;
		QCDcorr = 1.0;
	};

	double result=1;
	if (systematics == 0) result = (1.0 + EWKcorr);
	else if (systematics == 1) result = (1.0 + EWKcorr*(2.0 - QCDcorr));
	else if (systematics == -1) result = (1.0 + EWKcorr*QCDcorr);
	else result = 1;
	if(GenHMass<m4lThreshold) result=1;
	return result;
};

void reprocess_AllSamples(int erg_tev, int folder=2, int sample_min=0, int sample_max=kAllSamples, float mPOLE=125.6){
	char TREE_NAME[] = "SelectedTree";
	char erg_dir[1000];
	sprintf(erg_dir,"LHC_%iTeV",erg_tev);
	char cerg[1000];
	sprintf(cerg,"%iTeV",erg_tev);
	string cinput_common = user_dir;
	cinput_common = cinput_common + erg_dir;
	string coutput_common = cinput_common;
/*	string cinput_common = user_dir_newProduction;
	cinput_common = cinput_common + cerg;
	string coutput_common = user_dir_newProduction;
	coutput_common = coutput_common + cerg;
*/
	bool isData=false;
	if(folder==4) isData=true;
	if(isData && sample_min<0) sample_min=0;
	if(isData && sample_max<0) sample_max=3;
	if(isData && sample_min>2) sample_min=0;
	if(isData && sample_max>3) sample_max=3;

	TGraph* tg_interf = make_HZZ_LeptonInterferenceGraph();

	float c_ggzz = 1;
	char cinput_CggZZ[]="./data/ZZ4l-C_ggZZGraph.root";
	TFile* finput_CggZZ = new TFile(cinput_CggZZ,"read");
	string tgcggzzname = "C_ZZGG_AllFlavors";
	TGraphAsymmErrors* tgcggzz = (TGraphAsymmErrors*) finput_CggZZ->Get(tgcggzzname.c_str());
	if(folder<3) c_ggzz = tgcggzz->Eval(folder);
	else c_ggzz = tgcggzz->Eval(2); // 2mu2e
	finput_CggZZ->Close();

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

	char* cinput_PDFError[2]={
		"./data/CTEQ6LPDFErrorDown.root",
		"./data/CTEQ6LPDFErrorUp.root",
	};
	string tgPDFErrorname = "Graph";
	TFile* finput_PDFError[2];
	TGraph* tgPDFError[2];
	for(int f=0;f<2;f++){
		finput_PDFError[f] = new TFile( cinput_PDFError[f] ,"read");
		tgPDFError[f] = (TGraph*) finput_PDFError[f]->Get(tgPDFErrorname.c_str());
	};

	const int tgCTotalBkgSize=14;
	char cinput_ctotbkg[]="./data/ZZ4l-C_TotalBkgM4lGraph.root";
	TFile* finput_ctotbkg = TFile::Open(cinput_ctotbkg,"read");
	TGraph* tgtotalbkg[tgCTotalBkgSize];
	if(!isData) setCTotalBkgGraphs(folder,finput_ctotbkg, tgCTotalBkgSize, tgtotalbkg);

	if(sample_min<kQQBZZSamples_Dedicated && sample_max>=kGGSamples) make_QQBGG_Counters(erg_tev,folder);
	char INPUT_QQBGG_PROPER_NAME[]="HZZ4l_QQBGGSampleRecoCounter.root";
	string cqqbggProper = coutput_common;
	cqqbggProper = cqqbggProper + "/" + user_folder[folder] + "/";
	char cQQBGGfile[1000];
	sprintf(cQQBGGfile,"%s%s",cqqbggProper.c_str(),INPUT_QQBGG_PROPER_NAME);
	TFile* fQQBGGProper = new TFile(cQQBGGfile,"read");
	TH1F* hQQBGGProper = (TH1F*) fQQBGGProper->Get("hQQBGG_EfficiencyMatcher");

	Mela mela(erg_tev,mPOLE);
	mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
	TGraph* vaScale_4e = mela.vaScale_4e;
	TGraph* vaScale_4mu = mela.vaScale_4mu;
	TGraph* vaScale_2e2mu = mela.vaScale_2e2mu;
	TGraph* DggZZ_scalefactor = mela.DggZZ_scalefactor;

	for(int smp=sample_min; smp<sample_max; smp++){
		bool resurrection=false;
//		if(smp==kGGHSamples+18 && !isData) resurrection=true;
//		if(resurrection && folder==2) continue;
		if(isData){
			setCTotalBkgGraphs(smp,finput_ctotbkg, tgCTotalBkgSize, tgtotalbkg);
			cout << tgtotalbkg[0]->GetName() << endl;
		};
//		Some K factor-related note-taking:
		double sum_MC_weight_noK = 0;
		double sum_MC_weight_wK = 0;
		double sum_MC_weight_wKu = 0;
		double sum_MC_weight_wKd = 0;
		double sum_MC_weight_wPu = 0;
		double sum_MC_weight_wPd = 0;
		double Ku_MH = 1, PDFu_MH = 1;
		double Kd_MH = 1, PDFd_MH = 1;
		double K_MH = 1;

		string cinput_main = cinput_common;
		if(smp<(kAllSamples-1) || isData) cinput_main = cinput_main + "/" + user_folder[folder] + "/";
		else cinput_main = cinput_main + "/" + user_folder[3] + "/";
		string coutput_main = cinput_common;
		coutput_main = coutput_main + "/" + user_folder[folder] + "/";
//		if(resurrection) cinput_main = cinput_common + "/" + user_folder[2] + "/";

		char coutput[1000];
		if(!isData) sprintf(coutput,"%s%s%s",coutput_main.c_str(),sampleName[smp],"_Reprocessed.root");
		else sprintf(coutput,"%s%s%s",coutput_main.c_str(),data_files[smp],"_Reprocessed.root");
		char cinput[1000];
		if(!isData) sprintf(cinput,"%s%s%s",cinput_main.c_str(),sampleName[smp],".root");
		else sprintf(cinput,"%s%s%s",cinput_main.c_str(),data_files[smp],".root");
//		if(resurrection) sprintf(cinput,"%s%s%s",cinput_main.c_str(),sampleName[kGGHSamples+15],".root");
		cout << "Sample: " << smp << '\t' << cinput << endl;

		TChain* tree = new TChain(TREE_NAME);
		tree->Add(cinput);

		int RunNumber=-1;
		Long64_t EventNumber=-1;
		int RunNumber_CRreference=-1;
		Long64_t EventNumber_CRreference=-1;
		double FakeRate=-1;
		double FakeRate_w=-1;
		float ZXfake_weight=0;
		float ZXfake_weightProper=0;
		float ZXfake_weightReference=0;
		int CRflavor=-1;
		int CRreferenceflag=-1;
		int CRflag=-1;
		float CRZZMass=-1;
		TChain* CRtree;
		int nCRreferenceEntries=0;
		if(smp==(kAllSamples-1) && !isData){
			CRtree = new TChain(TREE_NAME);
			char cCRreference[1000];
			sprintf(cCRreference,"%s%s%s",cinput_main.c_str(),sampleName[smp],"_Reference.root");
			CRtree->Add(cCRreference);

			nCRreferenceEntries = CRtree->GetEntries();
			CRtree->SetBranchAddress("RunNumber",&RunNumber_CRreference);
			CRtree->SetBranchAddress("EventNumber",&EventNumber_CRreference);
			CRtree->SetBranchAddress("FakeRate",&FakeRate);
			CRtree->SetBranchAddress("FakeRate_w",&FakeRate_w);
			CRtree->SetBranchAddress("flavor",&CRflavor);
			CRtree->SetBranchAddress("flag",&CRreferenceflag);
			CRtree->SetBranchAddress("ZZMass",&CRZZMass);
		};

		int genFinalState;
		float mycosthetastar,myhelphi,myhelcosthetaZ1,myhelcosthetaZ2,myphistarZ1,myphistarZ2;
		float ZZMass,Z1Mass,Z2Mass;
		float GenHMass=0,GenZ1Mass=0,GenZ2Mass=0;
		int GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id;
		float m_l1minus_pT, m_l1plus_pT, m_l2minus_pT, m_l2plus_pT;
		float m_l1minus_eta, m_l1plus_eta, m_l2minus_eta, m_l2plus_eta;
		float m_l1minus_phi, m_l1plus_phi, m_l2minus_phi, m_l2plus_phi;
		float m_l1minus_E, m_l1plus_E, m_l2minus_E, m_l2plus_E;

		int Lep1Id,Lep2Id,Lep3Id,Lep4Id;
		short Z1ids,Z2ids;

		float p0plus_VAJHU;
		float p0hplus_VAJHU;
		float p0minus_VAJHU;
		float p0_g1prime2_VAJHU;
		float pg1g1prime2_VAJHU;
		float pg1g2_VAJHU;
		float pg1g4_VAJHU;
		float pg1g2_pi2_VAJHU;
		float pg1g4_pi2_VAJHU;

		float pzzzg_VAJHU;
		float pzzgg_VAJHU;
		float pzzzg_PS_VAJHU;
		float pzzgg_PS_VAJHU;
		float p0Zgs_VAJHU;
		float p0gsgs_VAJHU;
		float p0Zgs_PS_VAJHU;
		float p0gsgs_PS_VAJHU;
		float p0Zgs_g1prime2_VAJHU;
		float pzzzgs_g1prime2_VAJHU;
		float pzzzgs_g1prime2_pi2_VAJHU;
		
		float p1plusProdIndepVA, p1minusProdIndepVA;
		float p1plus_VAJHU,p1_VAJHU;
		float p2minimalVA, p2minimalVA_qqb, p2h2plusVA, p2h2plusVA_qqb;
		float p2h3plusVA, p2h3plusVA_qqb, p2hplusVA, p2hplusVA_qqb;
		float p2bplusVA, p2bplusVA_qqb, p2h6plusVA, p2h6plusVA_qqb, p2h7plusVA, p2h7plusVA_qqb, p2hminusVA, p2hminusVA_qqb, p2h9minusVA, p2h9minusVA_qqb, p2h10minusVA, p2h10minusVA_qqb;
		float p2mProdIndepVA, p2h2plusProdIndepVA, p2h3plusProdIndepVA, p2hplusProdIndepVA, p2bplusProdIndepVA, p2h6plusProdIndepVA, p2h7plusProdIndepVA, p2hminusProdIndepVA, p2h9minusProdIndepVA, p2h10minusProdIndepVA;

		float bkg_VAMCFM, bkg_prodIndep_VAMCFM;
		float bkg_VAMCFM_OLD,ggzz_p0plus_VAMCFM_OLD,ggzz_VAMCFM_OLD,p0plus_VAMCFM_OLD;
		float bkg_VAMCFM_NEW,ggzz_p0plus_VAMCFM_NEW,ggzz_VAMCFM_NEW,p0plus_VAMCFM_NEW;
		float ggzz_VAMCFM,ggzz_c1_VAMCFM,ggzz_c5_VAMCFM;
		float bkg_VAMCFM_noscale , ggzz_VAMCFM_noscale , ggHZZ_prob_pure , ggHZZ_prob_int , ggHZZ_prob_pure_noscale , ggHZZ_prob_int_noscale;
		float qqScale , ggScale;
		float p0plus_m4l, p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown, p0plus_m4l_ResUp, p0plus_m4l_ResDown;
		float bkg_m4l, bkg_m4l_ScaleUp, bkg_m4l_ScaleDown, bkg_m4l_ResUp, bkg_m4l_ResDown;
		float is_selected1=0;
		float is_selected2=0;
		float is_selected3=0;
		float is_selected=0;

		float MC_weight=1;
		float MC_weight_noxsec=1;
		float MC_weight_xsec=0;
		float MC_weight_ggZZLepInt=1;

		float MC_weight_norm=1;
		float MC_weight_PUWeight=1;
		float MC_weight_powhegWeight=1;
		float MC_weight_dataMC=1;
		float MC_weight_HqT=1;

		const int numSamples=52;
		const int firstNonZZhypo=39;
		int kNumSamples=numSamples;
		float MC_weight_spin0[numSamples]={0};
//		float MC_ME_as2mu2e_spin0[numSamples]={0};

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
		float MC_weight_QQZZEWK = 1;
		float MC_weight_QQZZEWK_down = 1;
		float MC_weight_QQZZEWK_up = 1;

		float Dgg10_VAMCFM=-99;
		float D_ggZZ=-99;
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
		float D_ZG=-99, D_GG=-99;
		float D_ZG_PS=-99, D_GG_PS=-99;
		float D_ZGint=-99, D_GGint=-99;
		float D_ZG_PSint=-99, D_GG_PSint=-99;
		float D_ZG_L1=-99,D_ZG_L1int=-99,D_ZG_L1int_perp=-99;

		float D_bkg=-99;
		float D_bkg_prodIndep=-99;
		float D_bkg_ScaleUp=-99;
		float D_bkg_prodIndep_ScaleUp=-99;
		float D_bkg_ScaleDown=-99;
		float D_bkg_prodIndep_ScaleDown=-99;
		float D_bkg_ResUp=-99;
		float D_bkg_prodIndep_ResUp=-99;
		float D_bkg_ResDown=-99;
		float D_bkg_prodIndep_ResDown=-99;

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

// Jet Variables
		float DiJetMass,DiJetMassPlus,DiJetMassMinus;
		short NJets30;
		float phjj_VAJHU_old,phjj_VAJHU_old_up,phjj_VAJHU_old_dn;
		float pvbf_VAJHU_old,pvbf_VAJHU_old_up,pvbf_VAJHU_old_dn;
		std::vector<double> myJetPt;
		std::vector<double> myJetSigma;
		std::vector<double> * JetPt=0;
		std::vector<double> * JetSigma=0;
		TBranch* bJetPt=0;
		TBranch* bJetSigma=0;

		float Djet_VAJHU=-99;
		float Djet_VAJHU_up=-99;
		float Djet_VAJHU_dn=-99;

// Spin 1 and 2 MEs
		if ((smp >= kGGSamples || smp == 0 || smp == 11 || (smp >= kGGHSamples && smp < kGGHSamples + 5) || (smp == kGGMCFMSamples + 1 || smp == kGGMCFMSamples + 4)) || isData){
			tree->SetBranchAddress("p1plus_prodIndep_VAJHU", &p1plusProdIndepVA);
			tree->SetBranchAddress("p1_prodIndep_VAJHU", &p1minusProdIndepVA);
			tree->SetBranchAddress("p1plus_VAJHU", &p1plus_VAJHU);
			tree->SetBranchAddress("p1_VAJHU", &p1_VAJHU);

			tree->SetBranchAddress("p2_VAJHU", &p2minimalVA);
			tree->SetBranchAddress("p2qqb_VAJHU", &p2minimalVA_qqb);
			tree->SetBranchAddress("p2h2plus_gg_VAJHU", &p2h2plusVA);
			tree->SetBranchAddress("p2h2plus_qqbar_VAJHU", &p2h2plusVA_qqb);
			tree->SetBranchAddress("p2h3plus_gg_VAJHU", &p2h3plusVA);
			tree->SetBranchAddress("p2h3plus_qqbar_VAJHU", &p2h3plusVA_qqb);
			tree->SetBranchAddress("p2hplus_VAJHU", &p2hplusVA);
			tree->SetBranchAddress("p2hplus_qqb_VAJHU", &p2hplusVA_qqb);
			tree->SetBranchAddress("p2bplus_VAJHU", &p2bplusVA);
			tree->SetBranchAddress("p2bplus_qqb_VAJHU", &p2bplusVA_qqb);
			tree->SetBranchAddress("p2h6plus_gg_VAJHU", &p2h6plusVA);
			tree->SetBranchAddress("p2h6plus_qqbar_VAJHU", &p2h6plusVA_qqb);
			tree->SetBranchAddress("p2h7plus_gg_VAJHU", &p2h7plusVA);
			tree->SetBranchAddress("p2h7plus_qqbar_VAJHU", &p2h7plusVA_qqb);
			tree->SetBranchAddress("p2hminus_VAJHU", &p2hminusVA);
			tree->SetBranchAddress("p2hminus_qqb_VAJHU", &p2hminusVA_qqb);
			tree->SetBranchAddress("p2h9minus_gg_VAJHU", &p2h9minusVA);
			tree->SetBranchAddress("p2h9minus_qqbar_VAJHU", &p2h9minusVA_qqb);
			tree->SetBranchAddress("p2h10minus_gg_VAJHU", &p2h10minusVA);
			tree->SetBranchAddress("p2h10minus_qqbar_VAJHU", &p2h10minusVA_qqb);
			tree->SetBranchAddress("p2_prodIndep_VAJHU", &p2mProdIndepVA);
			tree->SetBranchAddress("p2h2plus_prodIndep_VAJHU", &p2h2plusProdIndepVA);
			tree->SetBranchAddress("p2h3plus_prodIndep_VAJHU", &p2h3plusProdIndepVA);
			tree->SetBranchAddress("p2hplus_prodIndep_VAJHU", &p2hplusProdIndepVA);
			tree->SetBranchAddress("p2bplus_prodIndep_VAJHU", &p2bplusProdIndepVA);
			tree->SetBranchAddress("p2h6plus_prodIndep_VAJHU", &p2h6plusProdIndepVA);
			tree->SetBranchAddress("p2h7plus_prodIndep_VAJHU", &p2h7plusProdIndepVA);
			tree->SetBranchAddress("p2hminus_prodIndep_VAJHU", &p2hminusProdIndepVA);
			tree->SetBranchAddress("p2h9minus_prodIndep_VAJHU", &p2h9minusProdIndepVA);
			tree->SetBranchAddress("p2h10minus_prodIndep_VAJHU", &p2h10minusProdIndepVA);
			tree->SetBranchAddress("bkg_prodIndep_VAMCFM",&bkg_prodIndep_VAMCFM);
		};

		if(smp<(kAllSamples-1) && !isData){
			tree->SetBranchAddress("genFinalState", &genFinalState);
			tree->SetBranchAddress("MC_weight",&MC_weight);
			tree->SetBranchAddress("MC_weight_noxsec",&MC_weight_noxsec);
			tree->SetBranchAddress("MC_weight_norm",&MC_weight_norm);
			tree->SetBranchAddress("MC_weight_PUWeight",&MC_weight_PUWeight);
			tree->SetBranchAddress("MC_weight_powhegWeight",&MC_weight_powhegWeight);
			tree->SetBranchAddress("MC_weight_dataMC",&MC_weight_dataMC);
			tree->SetBranchAddress("MC_weight_HqT",&MC_weight_HqT);
			if (smp < kGGHSamples){
				tree->SetBranchAddress("kNumSamples", &kNumSamples);
				tree->SetBranchAddress("MC_weight_spin0", MC_weight_spin0);
			};
			tree->SetBranchAddress("GenHMass", &GenHMass);
			tree->SetBranchAddress("GenZ1Mass", &GenZ1Mass);
			tree->SetBranchAddress("GenZ2Mass", &GenZ2Mass);
		}
		else{ // Data and CR
			tree->SetBranchAddress("RunNumber",&RunNumber);
			tree->SetBranchAddress("EventNumber",&EventNumber);
//			tree->SetBranchAddress("Dgg10_VAMCFM",&Dgg10_VAMCFM);
			if(!isData){
				tree->SetBranchAddress("ZXfake_weight",&ZXfake_weight);
				tree->SetBranchAddress("CRflag",&CRflag);
			};
		};

		float MC_weight_QQBGG_VAMCFM=1;
		float MC_weight_QQBGGProper[2]={0};
		int kNumQQBGGWeights=4;
		if( (
			smp>=kGGHSamples && smp<kQQBZZSamples_Dedicated
			&& (smp<kGGHSamples+5 || smp>=kGGSamples || smp==kGGHSamples+15 || smp==kGGHSamples+18) 
			)
			&& !isData ){
			tree->SetBranchAddress("MC_weight_QQBGG_VAMCFM",&MC_weight_QQBGG_VAMCFM);
		};

		tree->SetBranchAddress("ZZMass", &ZZMass);
		tree->SetBranchAddress("Z1Mass", &Z1Mass);
		tree->SetBranchAddress("Z2Mass", &Z2Mass);
		tree->SetBranchAddress("helcosthetaZ1", &myhelcosthetaZ1);
		tree->SetBranchAddress("helcosthetaZ2", &myhelcosthetaZ2);
		tree->SetBranchAddress("helphi", &myhelphi);
		tree->SetBranchAddress("costhetastar", &mycosthetastar);
		tree->SetBranchAddress("phistarZ1", &myphistarZ1);
		tree->SetBranchAddress("Z1ids", &Z1ids);
		tree->SetBranchAddress("Z2ids", &Z2ids);
		tree->SetBranchAddress("Lep1ID", &Lep1Id);
		tree->SetBranchAddress("Lep2ID", &Lep2Id);
		tree->SetBranchAddress("Lep3ID", &Lep3Id);
		tree->SetBranchAddress("Lep4ID", &Lep4Id);

		tree->SetBranchAddress("p0plus_VAJHU",&p0plus_VAJHU);
		tree->SetBranchAddress("p0hplus_VAJHU",&p0hplus_VAJHU);
		tree->SetBranchAddress("p0minus_VAJHU",&p0minus_VAJHU);
		tree->SetBranchAddress("p0_g1prime2_VAJHU",&p0_g1prime2_VAJHU);
		tree->SetBranchAddress("pg1g1prime2_VAJHU",&pg1g1prime2_VAJHU);
		tree->SetBranchAddress("pg1g2_VAJHU",&pg1g2_VAJHU);
		tree->SetBranchAddress("pg1g4_VAJHU",&pg1g4_VAJHU);
		tree->SetBranchAddress("pg1g2_pi2_VAJHU",&pg1g2_pi2_VAJHU);
		tree->SetBranchAddress("pg1g4_pi2_VAJHU",&pg1g4_pi2_VAJHU);

		tree->SetBranchAddress("pzzzg_VAJHU",&pzzzg_VAJHU);
		tree->SetBranchAddress("pzzgg_VAJHU",&pzzgg_VAJHU);
		tree->SetBranchAddress("pzzzg_PS_VAJHU",&pzzzg_PS_VAJHU);
		tree->SetBranchAddress("pzzgg_PS_VAJHU",&pzzgg_PS_VAJHU);
		tree->SetBranchAddress("p0Zgs_VAJHU",&p0Zgs_VAJHU);
		tree->SetBranchAddress("p0gsgs_VAJHU",&p0gsgs_VAJHU);
		tree->SetBranchAddress("p0Zgs_PS_VAJHU",&p0Zgs_PS_VAJHU);
		tree->SetBranchAddress("p0gsgs_PS_VAJHU",&p0gsgs_PS_VAJHU);
		tree->SetBranchAddress("p0Zgs_g1prime2_VAJHU",&p0Zgs_g1prime2_VAJHU);
		tree->SetBranchAddress("pzzzgs_g1prime2_VAJHU",&pzzzgs_g1prime2_VAJHU);
		tree->SetBranchAddress("pzzzgs_g1prime2_pi2_VAJHU",&pzzzgs_g1prime2_pi2_VAJHU);

		tree->SetBranchAddress("p0plus_VAMCFM",&p0plus_VAMCFM_OLD);
		tree->SetBranchAddress("bkg_VAMCFM",&bkg_VAMCFM_OLD);
		tree->SetBranchAddress("ggzz_p0plus_VAMCFM",&ggzz_p0plus_VAMCFM_OLD);
		tree->SetBranchAddress("ggzz_VAMCFM",&ggzz_VAMCFM_OLD);

		tree->SetBranchAddress("p0plus_m4l",&p0plus_m4l);
		tree->SetBranchAddress("bkg_m4l",&bkg_m4l);
		tree->SetBranchAddress("p0plus_m4l_ScaleUp",&p0plus_m4l_ScaleUp);
		tree->SetBranchAddress("bkg_m4l_ScaleUp",&bkg_m4l_ScaleUp);
		tree->SetBranchAddress("p0plus_m4l_ScaleDown",&p0plus_m4l_ScaleDown);
		tree->SetBranchAddress("bkg_m4l_ScaleDown",&bkg_m4l_ScaleDown);
		tree->SetBranchAddress("p0plus_m4l_ResUp",&p0plus_m4l_ResUp);
		tree->SetBranchAddress("bkg_m4l_ResUp",&bkg_m4l_ResUp);
		tree->SetBranchAddress("p0plus_m4l_ResDown",&p0plus_m4l_ResDown);
		tree->SetBranchAddress("bkg_m4l_ResDown",&bkg_m4l_ResDown);

		tree->SetBranchAddress("NJets30",&NJets30);
		tree->SetBranchAddress("DiJetMass",&DiJetMass);
		tree->SetBranchAddress("DiJetMassPlus",&DiJetMassPlus);
		tree->SetBranchAddress("DiJetMassMinus",&DiJetMassMinus);
		tree->SetBranchAddress("JetPt",&JetPt,&bJetPt);
		tree->SetBranchAddress("JetSigma",&JetSigma,&bJetSigma);
		tree->SetBranchAddress("phjj_VAJHU_old",&phjj_VAJHU_old);
		tree->SetBranchAddress("phjj_VAJHU_old_up",&phjj_VAJHU_old_up);
		tree->SetBranchAddress("phjj_VAJHU_old_dn",&phjj_VAJHU_old_dn);
		tree->SetBranchAddress("pvbf_VAJHU_old",&pvbf_VAJHU_old);
		tree->SetBranchAddress("pvbf_VAJHU_old_up",&pvbf_VAJHU_old_up);
		tree->SetBranchAddress("pvbf_VAJHU_old_dn",&pvbf_VAJHU_old_dn);
		tree->SetBranchAddress("pvbf_VAJHU_old_dn",&pvbf_VAJHU_old_dn);


		TFile* foutput = new TFile(coutput,"recreate");
		TFile* finput = new TFile(cinput,"read");
		foutput->cd();
		TH1F* hCounters = (TH1F*) finput->Get("hCounters");
		foutput->WriteTObject(hCounters);
		TH2F* hSpin0Counters;
		TH2F* hSpin0Counters_4GeVcut;
		if(smp<kGGHSamples && !isData){
			foutput->cd();
			hSpin0Counters = (TH2F*) finput->Get("hCounters_spin0_RW");
			hSpin0Counters_4GeVcut = (TH2F*) finput->Get("hCounters_spin0_4GeVcutME_RW");
			foutput->WriteTObject(hSpin0Counters);
			foutput->WriteTObject(hSpin0Counters_4GeVcut);
		};
		finput->Close();
		foutput->cd();
		TTree* mytree = new TTree("SelectedTree","SelectedTree");
		mytree->SetAutoSave(6000000000);
		if(smp<(kAllSamples-1) && !isData){
			if (smp < kGGHSamples){
				mytree->Branch("kNumSamples", &kNumSamples);
				mytree->Branch("MC_weight_spin0", MC_weight_spin0, "MC_weight_spin0[kNumSamples]/F");
			};
//			if(smp<kGGHSamples) mytree->Branch("MC_ME_as2mu2e_spin0",MC_ME_as2mu2e_spin0,"MC_ME_as2mu2e_spin0[kNumSamples]/F");
			mytree->Branch("genFinalState", &genFinalState);
			mytree->Branch("MC_weight",&MC_weight);
			mytree->Branch("MC_weight_down",&MC_weight_down);
			mytree->Branch("MC_weight_up",&MC_weight_up);
			mytree->Branch("MC_weight_norm",&MC_weight_norm);
			mytree->Branch("MC_weight_PUWeight",&MC_weight_PUWeight);
			mytree->Branch("MC_weight_powhegWeight",&MC_weight_powhegWeight);
			mytree->Branch("MC_weight_dataMC",&MC_weight_dataMC);
			mytree->Branch("MC_weight_HqT",&MC_weight_HqT);
			mytree->Branch("MC_weight_noxsec",&MC_weight_noxsec);
			mytree->Branch("MC_weight_xsec",&MC_weight_xsec);
			mytree->Branch("MC_weight_Kfactor",&MC_weight_Kfactor);
			mytree->Branch("MC_weight_Kfactor_down",&MC_weight_Kfactor_down);
			mytree->Branch("MC_weight_Kfactor_up",&MC_weight_Kfactor_up);
			mytree->Branch("MC_weight_PDF_down",&MC_weight_PDF_down);
			mytree->Branch("MC_weight_PDF_up",&MC_weight_PDF_up);
			mytree->Branch("MC_weight_Kfactor_Norm_down",&MC_weight_Kfactor_Norm_down);
			mytree->Branch("MC_weight_Kfactor_Norm_up",&MC_weight_Kfactor_Norm_up);
			mytree->Branch("MC_weight_PDF_Norm_down",&MC_weight_PDF_Norm_down);
			mytree->Branch("MC_weight_PDF_Norm_up",&MC_weight_PDF_Norm_up);

			if (smp >= kGGHSamples){
				mytree->Branch("MC_weight_QQZZEWK", &MC_weight_QQZZEWK);
				mytree->Branch("MC_weight_QQZZEWK_down", &MC_weight_QQZZEWK_down);
				mytree->Branch("MC_weight_QQZZEWK_up", &MC_weight_QQZZEWK_up);
			};

			if (smp >= kGGHSamples + 2 && smp<kGGMCFMSamples) mytree->Branch("MC_weight_ggZZLepInt", &MC_weight_ggZZLepInt); // MCFM ggZZ mZZ-dependent lepton interference

			if( (
				smp>=kGGHSamples && smp<kQQBZZSamples_Dedicated
				&& (smp<kGGHSamples+5 || smp>=kGGSamples || smp==kGGHSamples+15 || smp==kGGHSamples+18) 
				)
				&& !isData ){
				mytree->Branch("kNumQQBGGWeights",&kNumQQBGGWeights);
				mytree->Branch("MC_weight_QQBGGProper",MC_weight_QQBGGProper,"MC_weight_QQBGGProper[kNumQQBGGWeights]/F");
			};

			mytree->Branch("GenHMass", &GenHMass);
			mytree->Branch("GenZ1Mass", &GenZ1Mass);
			mytree->Branch("GenZ2Mass", &GenZ2Mass);
		}
		else{
			mytree->Branch("RunNumber",&RunNumber);
			mytree->Branch("EventNumber",&EventNumber);
//			mytree->Branch("Dgg10_VAMCFM", &Dgg10_VAMCFM);
			mytree->Branch("helcosthetaZ1", &myhelcosthetaZ1);
			mytree->Branch("helcosthetaZ2", &myhelcosthetaZ2);
			mytree->Branch("helphi", &myhelphi);
			mytree->Branch("costhetastar", &mycosthetastar);
			mytree->Branch("phistarZ1", &myphistarZ1);
			if(!isData){
				mytree->Branch("ZXfake_weight",&ZXfake_weight);
				mytree->Branch("ZXfake_weightProper",&ZXfake_weightProper);
//				mytree->Branch("ZXfake_weightReference",&ZXfake_weightReference);
				mytree->Branch("flag",&CRreferenceflag);
			};
		};

		mytree->Branch("ZZMass", &ZZMass);
		mytree->Branch("Z1Mass", &Z1Mass);
		mytree->Branch("Z2Mass", &Z2Mass);
		mytree->Branch("helcosthetaZ1", &myhelcosthetaZ1);
		mytree->Branch("helcosthetaZ2", &myhelcosthetaZ2);
		mytree->Branch("helphi", &myhelphi);
		mytree->Branch("costhetastar", &mycosthetastar);
		mytree->Branch("phistarZ1", &myphistarZ1);

		mytree->Branch("Z1ids", &Z1ids);
		mytree->Branch("Z2ids", &Z2ids);
		mytree->Branch("Lep1ID", &Lep1Id);
		mytree->Branch("Lep2ID", &Lep2Id);
		mytree->Branch("Lep3ID", &Lep3Id);
		mytree->Branch("Lep4ID", &Lep4Id);

		mytree->Branch("D_g1_vs_g2_phi0",&D_g2);
		mytree->Branch("D_g1_vs_g4_phi0",&D_g4);
		mytree->Branch("D_g2int_phi0",&D_g2int);
		mytree->Branch("D_g4int_phi0",&D_g4int);
		mytree->Branch("D_g2int_phi90",&D_g2int_perp);
		mytree->Branch("D_g4int_phi90",&D_g4int_perp);
		mytree->Branch("D_g1Q2_phi0",&D_g1q2);
		mytree->Branch("D_g1Q2int_phi0",&D_g1q2int);

		mytree->Branch("D_ZG",&D_ZG);
		mytree->Branch("D_GG",&D_GG);
		mytree->Branch("D_ZG_PS",&D_ZG_PS);
		mytree->Branch("D_GG_PS",&D_GG_PS);
		mytree->Branch("D_ZGint",&D_ZGint);
		mytree->Branch("D_GGint",&D_GGint);
		mytree->Branch("D_ZG_PSint",&D_ZG_PSint);
		mytree->Branch("D_GG_PSint",&D_GG_PSint);
		mytree->Branch("D_ZG_L1",&D_ZG_L1);
		mytree->Branch("D_ZG_L1int_phi0",&D_ZG_L1int);
		mytree->Branch("D_ZG_L1int_phi90",&D_ZG_L1int_perp);

		mytree->Branch("D_bkg",&D_bkg);
		mytree->Branch("D_bkg_kin",&D_bkg_kin);
		mytree->Branch("D_bkg_ScaleUp",&D_bkg_ScaleUp);
		mytree->Branch("D_bkg_ScaleDown",&D_bkg_ScaleDown);
		mytree->Branch("D_bkg_ResUp",&D_bkg_ResUp);
		mytree->Branch("D_bkg_ResDown",&D_bkg_ResDown);

		mytree->Branch("D_ggZZ",&D_ggZZ);
/*		mytree->Branch("D_Gamma_r1",&D_Gamma_r1);
		mytree->Branch("D_Gamma_r5",&D_Gamma_r5);
		mytree->Branch("D_Gamma_r10",&D_Gamma_r10);
		mytree->Branch("D_Gamma_r15",&D_Gamma_r15);
		mytree->Branch("D_Gamma_r20",&D_Gamma_r20);
		mytree->Branch("D_Gamma_r25",&D_Gamma_r25);
*/
//		mytree->Branch("D_Gamma_gg_r1",&D_Gamma_gg_r1);
//		mytree->Branch("D_Gamma_gg_r5",&D_Gamma_gg_r5);
		mytree->Branch("D_Gamma_gg_r10",&D_Gamma_gg_r10);
/*		mytree->Branch("D_Gamma_gg_r15",&D_Gamma_gg_r15);
		mytree->Branch("D_Gamma_gg_r20",&D_Gamma_gg_r20);
		mytree->Branch("D_Gamma_gg_r25",&D_Gamma_gg_r25);
		mytree->Branch("D_Gamma",&D_Gamma);
		mytree->Branch("D_Gamma_int",&D_Gamma_int);
*/

// Jet Variables
		mytree->Branch("NJets30",&NJets30);
		mytree->Branch("DiJetMass",&DiJetMass);
		mytree->Branch("DiJetMassPlus",&DiJetMassPlus);
		mytree->Branch("DiJetMassMinus",&DiJetMassMinus);
		mytree->Branch("Djet_VAJHU",&Djet_VAJHU);
		mytree->Branch("Djet_VAJHU_up",&Djet_VAJHU_up);
		mytree->Branch("Djet_VAJHU_dn",&Djet_VAJHU_dn);
		mytree->Branch("JetPt",&myJetPt);
		mytree->Branch("JetSigma",&myJetSigma);
		
// Spin 1 and 2 discriminants
		if ((smp >= kGGSamples || smp == 0 || smp == 11 || (smp >= kGGHSamples && smp<kGGHSamples + 5) || (smp == kGGMCFMSamples + 1 || smp == kGGMCFMSamples + 4)) || isData){
			mytree->Branch("D_bkg_prodIndep",&D_bkg_prodIndep);
			mytree->Branch("D_bkg_prodIndep_ScaleUp",&D_bkg_prodIndep_ScaleUp);
			mytree->Branch("D_bkg_prodIndep_ScaleDown",&D_bkg_prodIndep_ScaleDown);
			mytree->Branch("D_bkg_prodIndep_ResUp",&D_bkg_prodIndep_ResUp);
			mytree->Branch("D_bkg_prodIndep_ResDown",&D_bkg_prodIndep_ResDown);

			mytree->Branch("p1plusProdIndepKD", &p1plusProdIndepKD);
			mytree->Branch("p1minusProdIndepKD", &p1minusProdIndepKD);
			mytree->Branch("p1plusKD", &p1plusKD);
			mytree->Branch("p1minusKD", &p1minusKD);

			mytree->Branch("graviKD",&graviKD);
			mytree->Branch("qqgraviKD",&qqgraviKD);
			mytree->Branch("p2h2plusKD",&p2h2plusKD);
			mytree->Branch("p2h2plus_qqb_KD",&p2h2plus_qqb_KD);
			mytree->Branch("p2h3plusKD",&p2h3plusKD);
			mytree->Branch("p2h3plus_qqb_KD",&p2h3plus_qqb_KD);
			mytree->Branch("p2hplusKD",&p2hplusKD);
			mytree->Branch("p2hplus_qqb_KD",&p2hplus_qqb_KD);
			mytree->Branch("p2bplusKD",&p2bplusKD);
			mytree->Branch("p2bplus_qqb_KD",&p2bplus_qqb_KD);
			mytree->Branch("p2h6plusKD",&p2h6plusKD);
			mytree->Branch("p2h6plus_qqb_KD",&p2h6plus_qqb_KD);
			mytree->Branch("p2h7plusKD",&p2h7plusKD);
			mytree->Branch("p2h7plus_qqb_KD",&p2h7plus_qqb_KD);
			mytree->Branch("p2hminusKD",&p2hminusKD);
			mytree->Branch("p2hminus_qqb_KD",&p2hminus_qqb_KD);
			mytree->Branch("p2h9minusKD",&p2h9minusKD);
			mytree->Branch("p2h9minus_qqb_KD",&p2h9minus_qqb_KD);
			mytree->Branch("p2h10minusKD",&p2h10minusKD);
			mytree->Branch("p2h10minus_qqb_KD",&p2h10minus_qqb_KD);
			mytree->Branch("p2mProdIndepKD",&p2mProdIndepKD);
			mytree->Branch("p2h2plusProdIndepKD",&p2h2plusProdIndepKD);
			mytree->Branch("p2h3plusProdIndepKD",&p2h3plusProdIndepKD);
			mytree->Branch("p2hplusProdIndepKD",&p2hplusProdIndepKD);
			mytree->Branch("p2bplusProdIndepKD",&p2bplusProdIndepKD);
			mytree->Branch("p2h6plusProdIndepKD",&p2h6plusProdIndepKD);
			mytree->Branch("p2h7plusProdIndepKD",&p2h7plusProdIndepKD);
			mytree->Branch("p2hminusProdIndepKD",&p2hminusProdIndepKD);
			mytree->Branch("p2h9minusProdIndepKD",&p2h9minusProdIndepKD);
			mytree->Branch("p2h10minusProdIndepKD",&p2h10minusProdIndepKD);
		};
/*		if(isData){
			mytree->Branch("p0plus_VAMCFM_OLD",&p0plus_VAMCFM_OLD);
			mytree->Branch("bkg_VAMCFM_OLD",&bkg_VAMCFM_OLD);
			mytree->Branch("ggzz_p0plus_VAMCFM_OLD",&ggzz_p0plus_VAMCFM_OLD);
			mytree->Branch("ggzz_VAMCFM_OLD",&ggzz_VAMCFM_OLD);
			mytree->Branch("p0plus_VAMCFM_NEW",&p0plus_VAMCFM_NEW);
			mytree->Branch("bkg_VAMCFM_NEW",&bkg_VAMCFM_NEW);
			mytree->Branch("ggzz_p0plus_VAMCFM_NEW",&ggzz_p0plus_VAMCFM_NEW);
			mytree->Branch("ggzz_VAMCFM_NEW",&ggzz_VAMCFM_NEW);
		};
*/
/*
		mytree->Branch("p0plus_m4l",&p0plus_m4l);
		mytree->Branch("bkg_m4l",&bkg_m4l);
		mytree->Branch("p0plus_m4l_ScaleUp",&p0plus_m4l_ScaleUp);
		mytree->Branch("bkg_m4l_ScaleUp",&bkg_m4l_ScaleUp);
		mytree->Branch("p0plus_m4l_ScaleDown",&p0plus_m4l_ScaleDown);
		mytree->Branch("bkg_m4l_ScaleDown",&bkg_m4l_ScaleDown);
		mytree->Branch("p0plus_m4l_ResUp",&p0plus_m4l_ResUp);
		mytree->Branch("bkg_m4l_ResUp",&bkg_m4l_ResUp);
		mytree->Branch("p0plus_m4l_ResDown",&p0plus_m4l_ResDown);
		mytree->Branch("bkg_m4l_ResDown",&bkg_m4l_ResDown);
*/
		int nEntries = tree->GetEntries();
		cout << "Entries " << nEntries << endl;
		for(int ev=0;ev<nEntries;ev++){
			tree->GetEntry(ev);
			if(Z2Mass<=12) continue;
			if (smp < kGGHSamples && !isData){
				int iZGhypo = firstNonZZhypo;
				if (!(GenZ1Mass>4 && GenZ2Mass>4)){
					for (int p = iZGhypo; p < kNumSamples; p++) MC_weight_spin0[p] = 0;
				};
			};

			bool CRmatched = false;
			bool debug = false;
			if((ev%100000)==0){
				cout << "Event: " <<  ev << "... XSEC = ";
			};
			if(smp==(kAllSamples-1) && !isData){
				if(EventNumber!=EventNumber) continue;
				for(int CRrefev=0;CRrefev<nCRreferenceEntries;CRrefev++){
					CRtree->GetEntry(CRrefev);

//					if(EventNumber!=EventNumber_CRreference) continue;
					if(ZZMass!=CRZZMass) continue;
					if(RunNumber!=RunNumber_CRreference) continue;
					if(folder == 0 && CRflavor!=1) continue;
					if(folder == 1 && CRflavor!=0) continue;
					if(folder == 2 && !(CRflavor==2 || CRflavor==3) ) continue;
					if (CRflavor>3) continue;

					ZXfake_weightProper = FakeRate_w;
					ZXfake_weightReference = FakeRate;

					if((CRrefev%100)==0) cout << "CR reference event: " << CRrefev << endl;

					CRmatched=true;
					break;
				};
				if(!CRmatched) continue;
			};

/*
			if(!isData){
//				if(folder==2 && !resurrection){
				if(folder==2){
					Lep1Id=-11;
					Lep2Id=11;
					Lep3Id=-13;
					Lep4Id=13;
				}
				else if(folder==1){
					Lep1Id=-11;
					Lep2Id=11;
					Lep3Id=-11;
					Lep4Id=11;
				}
				else{
					Lep1Id=-13;
					Lep2Id=13;
					Lep3Id=-13;
					Lep4Id=13;
				};
			}
			else{
				if(smp==2){
					Lep1Id=-11;
					Lep2Id=11;
					Lep3Id=-13;
					Lep4Id=13;
				}
				else if(smp==1){
					Lep1Id=-11;
					Lep2Id=11;
					Lep3Id=-11;
					Lep4Id=11;
				}
				else{
					Lep1Id=-13;
					Lep2Id=13;
					Lep3Id=-13;
					Lep4Id=13;
				};

			};
*/
			if (isData || smp==(kAllSamples-1) || (smp>=kQQBZZSamples_Dedicated && smp<kVBF_Phantom) ){
				if(Z1ids==-13*13){
					Lep1Id=-13;
					Lep2Id=13;
				}
				else if(Z1ids==-11*11){
					Lep1Id=-11;
					Lep2Id=11;
				};
				if(Z2ids==-13*13){
					Lep3Id=-13;
					Lep4Id=13;
				}
				else if(Z2ids==-11*11){
					Lep3Id=-11;
					Lep4Id=11;
				};
			};

			if(MC_weight<0){
				if (ev == 0){
					cout << "WARNING: XSEC IS LESS THAN 0" << endl;
					cout << "ev " << ev << " MC_weight: " << MC_weight << endl;
				};
				if(smp>=11 && smp<kGGHSamples){
					if (erg_tev == 7) MC_weight *= 1.9786800;//-2.0055388;
					if(erg_tev == 8) MC_weight *= -2.5198800;//-2.562639;
					if(MC_weight<0) cout << "Not fixed!" << endl;
				}
				else cout << "Sample: " << smp << " weight < 0 !!!" << endl; // Keep bothering the user for all events...
			};
//			if(resurrection){
//				if(folder==0) MC_weight *= 0.58185465;
//				if(folder==1) MC_weight *= 0.401331435;
//			};
//const int kQQBZZSamples_Dedicated=66;
//const int kVBF_Phantom=78;
//			if( smp>=(kGGHSamples+2) && smp<kGGMCFMSamples && ((smp-(kGGHSamples+2))%3)!=2 && erg_tev==8) MC_weight = MC_weight/2.0; // TEMPORARY FIX TO 4e, 4mu MCFM XSEC
			MC_weight_xsec = MC_weight/MC_weight_noxsec;
/*			if (smp >= kQQBZZSamples_Dedicated && smp < kVBF_Phantom){ // Workaround for Phantom
				MC_weight_xsec = 1.0/MC_weight_xsec*1.0e6;
				MC_weight = MC_weight_noxsec * MC_weight_xsec;
			};
*/			if((ev%100000)==0) cout << MC_weight_xsec << endl;
			if(MC_weight_noxsec==0 && !isData) continue;

// Jet clearing for filling...
			myJetPt.clear();
			myJetSigma.clear();
			for (int i = 0; i < JetPt->size(); i++){
				myJetPt.push_back(JetPt->at(i));
				myJetSigma.push_back(JetSigma->at(i));
			};


			int lepIdOrdered[4]={ Lep1Id,Lep2Id,Lep3Id,Lep4Id };
			float angularOrdered[8]={ZZMass,Z1Mass,Z2Mass,mycosthetastar,myhelcosthetaZ1,myhelcosthetaZ2,myhelphi,myphistarZ1};

			std::vector<int> partId; 
			for(int it=0;it<4;it++) partId.push_back(lepIdOrdered[it]);

			MC_weight_Kfactor = 1;
			MC_weight_Kfactor_up = 1;
			MC_weight_Kfactor_down = 1;
			MC_weight_PDF_down=0;
			MC_weight_PDF_up=0;
			MC_weight_ggZZLepInt=1;

			if(smp<kGGSamples && !isData){
				MC_weight_Kfactor = tgkf->Eval(GenHMass);
				MC_weight_Kfactor_up = tgkf_up->Eval(GenHMass);
				MC_weight_Kfactor_down = tgkf_down->Eval(GenHMass);

				MC_weight_alphaS = evaluateAlphaSShift(GenHMass, 1, 0,mPOLE);
				MC_weight_alphaS_down = evaluateAlphaSShift(GenHMass, 1, -1,mPOLE);
				MC_weight_alphaS_up = evaluateAlphaSShift(GenHMass, 1, 1,mPOLE);
				MC_weight_alphaS_down *= MC_weight_alphaS;
				MC_weight_alphaS_up *= MC_weight_alphaS;

				if(smp<kGGMCFMSamples && smp>=kGGHSamples+2){
					MC_weight_alphaS_down /= MC_weight_alphaS;
					MC_weight_alphaS_up /= MC_weight_alphaS;
					MC_weight_alphaS = 1;
				};

				float MC_weight_Kfactor_MH = tgkf->Eval(mPOLE);
				float MC_weight_Kfactor_up_MH = tgkf_up->Eval(mPOLE);
				float MC_weight_Kfactor_down_MH = tgkf_down->Eval(mPOLE);

				float MC_weight_alphaS_MH = evaluateAlphaSShift(mPOLE, 1, 0,mPOLE);
				float MC_weight_alphaS_down_MH = evaluateAlphaSShift(mPOLE, 1, -1,mPOLE);
				float MC_weight_alphaS_up_MH = evaluateAlphaSShift(mPOLE, 1, 1,mPOLE);
				MC_weight_alphaS_down_MH *= MC_weight_alphaS_MH;
				MC_weight_alphaS_up_MH *= MC_weight_alphaS_MH;
				if(smp<kGGMCFMSamples && smp>=kGGHSamples+2){
					MC_weight_alphaS_down_MH /= MC_weight_alphaS_MH;
					MC_weight_alphaS_up_MH /= MC_weight_alphaS_MH;
					MC_weight_alphaS_MH = 1;

					if(genFinalState<2) MC_weight_ggZZLepInt = tg_interf->Eval(GenHMass);
				};

				float MC_weight_PDF_down_MH = 0;
				float MC_weight_PDF_up_MH = 0;

				if(smp>=kGGHSamples){
					MC_weight_PDF_down = evaluatePDFUncertainty(GenHMass,tgPDFError,-1);
					MC_weight_PDF_up = evaluatePDFUncertainty(GenHMass,tgPDFError,1);
					MC_weight_PDF_down_MH = evaluatePDFUncertainty(mPOLE,tgPDFError,-1);
					MC_weight_PDF_up_MH = evaluatePDFUncertainty(mPOLE,tgPDFError,1);
/*
					MC_weight_PDF_down *= 0.01;
					MC_weight_PDF_up *= 0.01;
					MC_weight_PDF_down_MH *= 0.01;
					MC_weight_PDF_up_MH *= 0.01;

					MC_weight_PDF_down = MC_weight_PDF_down/(1.0 + MC_weight_PDF_down);
					MC_weight_PDF_up = MC_weight_PDF_up/(1.0+MC_weight_PDF_up);
					MC_weight_PDF_down_MH = MC_weight_PDF_down_MH/(1.0 + MC_weight_PDF_down_MH);
					MC_weight_PDF_up_MH = MC_weight_PDF_up_MH/(1.0+MC_weight_PDF_up_MH);
*/
				};

				MC_weight_Kfactor_up *= MC_weight_alphaS_up;
				MC_weight_Kfactor_down *= MC_weight_alphaS_down;
				MC_weight_Kfactor *= MC_weight_alphaS;
				MC_weight_Kfactor_up_MH *= MC_weight_alphaS_up_MH;
				MC_weight_Kfactor_down_MH *= MC_weight_alphaS_down_MH;
				MC_weight_Kfactor_MH *= MC_weight_alphaS_MH;

				if(smp<kGGHSamples){
					MC_weight_Kfactor_up /= MC_weight_Kfactor;
					MC_weight_Kfactor_down /= MC_weight_Kfactor;
					MC_weight_Kfactor = 1;

					MC_weight_Kfactor_up_MH /= MC_weight_Kfactor_MH;
					MC_weight_Kfactor_down_MH /= MC_weight_Kfactor_MH;
					MC_weight_Kfactor_MH = 1;
				};
				MC_weight_up = 1.0 + sqrt( pow(MC_weight_Kfactor_up/MC_weight_Kfactor - 1 , 2) + pow(MC_weight_PDF_up,2) );
				MC_weight_down = 1.0 - sqrt( pow(MC_weight_Kfactor_down/MC_weight_Kfactor - 1 , 2) + pow(MC_weight_PDF_down,2) );
				float MC_weight_up_MH = 1.0 + sqrt( pow(MC_weight_Kfactor_up_MH/MC_weight_Kfactor_MH - 1 , 2) + pow(MC_weight_PDF_up_MH,2) );
				float MC_weight_down_MH = 1.0 - sqrt( pow(MC_weight_Kfactor_down_MH/MC_weight_Kfactor_MH - 1 , 2) + pow(MC_weight_PDF_down_MH,2) );

				MC_weight_up /= MC_weight_up_MH;
				MC_weight_down /= MC_weight_down_MH;

				MC_weight_up *= MC_weight;
				MC_weight_down *= MC_weight;

				MC_weight_Kfactor_Norm_up = MC_weight_Kfactor_up/MC_weight_Kfactor_up_MH;
				MC_weight_Kfactor_Norm_down = MC_weight_Kfactor_down/MC_weight_Kfactor_down_MH;
				MC_weight_PDF_Norm_up = (1.0 + (MC_weight_PDF_up) )/(1.0 + (MC_weight_PDF_up_MH) );
				MC_weight_PDF_Norm_down = (1.0 + (MC_weight_PDF_down) )/(1.0 + (MC_weight_PDF_down_MH) );

				if (((smp < kGGHSamples + 11) || (smp<kGGMCFMSamples+6 && smp>=kGGMCFMSamples)) && ZZMass >= (mPOLE - 20.0) && ZZMass<(mPOLE + 15.0)){
					sum_MC_weight_noK += MC_weight;
					sum_MC_weight_wK += MC_weight*MC_weight_Kfactor;
					sum_MC_weight_wKu += MC_weight*MC_weight_Kfactor_up;
					sum_MC_weight_wKd += MC_weight*MC_weight_Kfactor_down;
					sum_MC_weight_wPu += MC_weight*(1.0+MC_weight_PDF_up);
					sum_MC_weight_wPd += MC_weight*(1.0+MC_weight_PDF_down);
					K_MH = MC_weight_Kfactor_MH;
					Ku_MH = MC_weight_Kfactor_up_MH;
					Kd_MH = MC_weight_Kfactor_down_MH;
					PDFu_MH = (1.0 + MC_weight_PDF_up_MH);
					PDFd_MH = (1.0 + MC_weight_PDF_down_MH);
				};
			}
			else if(!isData){
				for(int si=0;si<kNumSamples;si++) MC_weight_spin0[si]=1;
				MC_weight_up = MC_weight;
				MC_weight_down = MC_weight;
			};

			if (smp >= kGGHSamples && smp < kQQBZZSamples_Dedicated && !isData){
				MC_weight_QQZZEWK = getQQZZEWKCorrection(GenHMass, 0);
				MC_weight_QQZZEWK_down = getQQZZEWKCorrection(GenHMass, -1);
				MC_weight_QQZZEWK_up = getQQZZEWKCorrection(GenHMass, 1);
			};

			if( (
				smp>=kGGHSamples && smp<kQQBZZSamples_Dedicated
				&& (smp<kGGHSamples+5 || smp>=kGGSamples || smp==kGGHSamples+15 || smp==kGGHSamples+18) 
				)
				&& !isData ){
				double chosenXSEC[2] = {0};
				if(smp<kGGHSamples+2){
					chosenXSEC[0] = hQQBGGProper->GetBinContent(11);
					chosenXSEC[1] = hQQBGGProper->GetBinContent(12);
					MC_weight_QQBGGProper[0] = hQQBGGProper->GetBinContent(1)*chosenXSEC[0]*MC_weight_alphaS;
					MC_weight_QQBGGProper[1] = MC_weight_QQBGG_VAMCFM*hQQBGGProper->GetBinContent(6)*chosenXSEC[1]*MC_weight_alphaS;
				}
				else if(smp<kGGHSamples+5){
					chosenXSEC[0] = hQQBGGProper->GetBinContent(11);
					chosenXSEC[1] = hQQBGGProper->GetBinContent(12);
					MC_weight_QQBGGProper[0] = hQQBGGProper->GetBinContent(2)*chosenXSEC[0];
					MC_weight_QQBGGProper[1] = MC_weight_QQBGG_VAMCFM*hQQBGGProper->GetBinContent(7)*chosenXSEC[1];
				}
				else if(smp<kGGSamples){
					chosenXSEC[0] = hQQBGGProper->GetBinContent(11);
					chosenXSEC[1] = hQQBGGProper->GetBinContent(12);
					MC_weight_QQBGGProper[0] = hQQBGGProper->GetBinContent(3)*chosenXSEC[0]*MC_weight_alphaS;
					MC_weight_QQBGGProper[1] = MC_weight_QQBGG_VAMCFM*hQQBGGProper->GetBinContent(8)*chosenXSEC[1]*MC_weight_alphaS;
				}
				else if(smp<kQQBZZSamples){
					chosenXSEC[0] = hQQBGGProper->GetBinContent(11);
					chosenXSEC[1] = hQQBGGProper->GetBinContent(12);
					MC_weight_QQBGGProper[0] = MC_weight_QQBGG_VAMCFM*hQQBGGProper->GetBinContent(4)*chosenXSEC[0];
					MC_weight_QQBGGProper[1] = hQQBGGProper->GetBinContent(9)*chosenXSEC[1];
				}
				else if(smp<kQQBZZSamples_Dedicated){
					chosenXSEC[0] = hQQBGGProper->GetBinContent(11);
					chosenXSEC[1] = hQQBGGProper->GetBinContent(12);
					MC_weight_QQBGGProper[0] = MC_weight_QQBGG_VAMCFM*hQQBGGProper->GetBinContent(5)*chosenXSEC[0];
					MC_weight_QQBGGProper[1] = hQQBGGProper->GetBinContent(10)*chosenXSEC[1];
				};
			};

/*			if(resurrection){
				mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
				bkg_VAMCFM_OLD = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, true);
				mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
				ggzz_VAMCFM_OLD = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, false);
				mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
				ggzz_p0plus_VAMCFM_OLD = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, false);
				mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
				p0plus_VAMCFM_OLD = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, false);
			};
*/
/*			if(isData){
				mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
				bkg_VAMCFM_NEW = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, true);
				mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
				ggzz_VAMCFM_NEW = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, false);
				mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
				ggzz_p0plus_VAMCFM_NEW = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, false);
				mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
				p0plus_VAMCFM_NEW = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, false);
			};
*/
			bkg_VAMCFM=bkg_VAMCFM_OLD;
			ggzz_VAMCFM_noscale=ggzz_VAMCFM_OLD;
			ggHZZ_prob_pure_noscale = p0plus_VAMCFM_OLD;
			ggHZZ_prob_int_noscale = ggzz_p0plus_VAMCFM_OLD;

			ggHZZ_prob_int_noscale = ggHZZ_prob_int_noscale - ggHZZ_prob_pure_noscale - ggzz_VAMCFM_noscale;

			ggScale=1.0;
			qqScale=1.0;
			if(abs(lepIdOrdered[0])==abs(lepIdOrdered[1]) &&
				abs(lepIdOrdered[0])==abs(lepIdOrdered[2]) &&
				abs(lepIdOrdered[0])==abs(lepIdOrdered[3])){
					if(abs(lepIdOrdered[0])==11){
						if(ZZMass > 900) qqScale = vaScale_4e->Eval(900.);
						else if (ZZMass <  100 ) qqScale = vaScale_4e->Eval(100.);
						else qqScale = vaScale_4e->Eval(ZZMass);

						if(ZZMass > 900) ggScale = vaScale_4e->Eval(900.);
						else if (ZZMass <  110 ) ggScale = vaScale_4e->Eval(110.);
						else ggScale = vaScale_4e->Eval(ZZMass);
					}
					else{
						if(ZZMass > 900) qqScale = vaScale_4mu->Eval(900.);
						else if (ZZMass <  100 ) qqScale = vaScale_4mu->Eval(100.);
						else qqScale = vaScale_4mu->Eval(ZZMass);

						if(ZZMass > 900) ggScale = vaScale_4mu->Eval(900.);
						else if (ZZMass <  110 ) ggScale = vaScale_4mu->Eval(110.);
						else ggScale = vaScale_4mu->Eval(ZZMass);
					};
			}
			else{
				if(ZZMass > 900) qqScale = vaScale_2e2mu->Eval(900.);
				else if (ZZMass <  100 ) qqScale = vaScale_2e2mu->Eval(100.);
				else qqScale = vaScale_2e2mu->Eval(ZZMass);

				if(ZZMass > 900) ggScale = vaScale_2e2mu->Eval(900.);
				else if (ZZMass <  110 ) ggScale = vaScale_2e2mu->Eval(110.);
				else ggScale = vaScale_2e2mu->Eval(ZZMass);
			};
			if(ZZMass > 900) ggScale /= DggZZ_scalefactor->Eval(900.);
			else if (ZZMass <  110 ) ggScale /= DggZZ_scalefactor->Eval(110.);
			else ggScale /= DggZZ_scalefactor->Eval(ZZMass);

			bkg_VAMCFM_noscale = bkg_VAMCFM/qqScale;

			ggzz_VAMCFM = ggzz_VAMCFM_noscale*ggScale;
			ggHZZ_prob_pure = ggHZZ_prob_pure_noscale*ggScale;
			ggHZZ_prob_int = ggHZZ_prob_int_noscale*ggScale;

// Construct the discriminants
			float CTotalBkg[tgCTotalBkgSize] = {0};
			float myDGamma[tgCTotalBkgSize] = {0};
			getCTotalBkg (ZZMass, tgCTotalBkgSize, CTotalBkg, tgtotalbkg);
			constructVariousDGamma(c_ggzz, bkg_VAMCFM_noscale, ggzz_VAMCFM_noscale, ggHZZ_prob_pure_noscale, ggHZZ_prob_int_noscale, tgCTotalBkgSize, CTotalBkg, myDGamma);
//			if(isData) cout << bkg_VAMCFM_noscale << '\t' << ggzz_VAMCFM_noscale << '\t' << ggHZZ_prob_pure_noscale << '\t' << ggHZZ_prob_int_noscale << endl;

			D_ggZZ = ggzz_VAMCFM / (ggzz_VAMCFM + bkg_VAMCFM);

			D_Gamma_r1 = myDGamma[0];
			D_Gamma_r5 = myDGamma[2];
			D_Gamma_r10 = myDGamma[2];
			D_Gamma_r15 = myDGamma[3];
			D_Gamma_r20 = myDGamma[4];
			D_Gamma_r25 = myDGamma[5];
			D_Gamma_gg_r1 = myDGamma[6];
			D_Gamma_gg_r5 = myDGamma[7];
			D_Gamma_gg_r10 = myDGamma[8];
			D_Gamma_gg_r15 = myDGamma[9];
			D_Gamma_gg_r20 = myDGamma[10];
			D_Gamma_gg_r25 = myDGamma[11];
			D_Gamma = myDGamma[12];
			D_Gamma_int = myDGamma[13];

			D_bkg_kin = p0plus_VAJHU / ( p0plus_VAJHU + bkg_VAMCFM );
			D_bkg = p0plus_VAJHU*p0plus_m4l / ( p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l );
			D_bkg_ScaleUp = (p0plus_VAJHU*p0plus_m4l_ScaleUp) / (p0plus_VAJHU*p0plus_m4l_ScaleUp + bkg_VAMCFM_OLD*bkg_m4l_ScaleUp);
			D_bkg_ScaleDown = (p0plus_VAJHU*p0plus_m4l_ScaleDown) / (p0plus_VAJHU*p0plus_m4l_ScaleDown + bkg_VAMCFM_OLD*bkg_m4l_ScaleDown);
			D_bkg_ResUp = (p0plus_VAJHU*p0plus_m4l_ResUp) / (p0plus_VAJHU*p0plus_m4l_ResUp + bkg_VAMCFM_OLD*bkg_m4l_ResUp);
			D_bkg_ResDown = (p0plus_VAJHU*p0plus_m4l_ResDown) / (p0plus_VAJHU*p0plus_m4l_ResDown + bkg_VAMCFM_OLD*bkg_m4l_ResDown);

			if(phjj_VAJHU_old>=0 && pvbf_VAJHU_old>=0) Djet_VAJHU = pvbf_VAJHU_old / ( pvbf_VAJHU_old + phjj_VAJHU_old ) ;
			else{
				if((phjj_VAJHU_old<0 || pvbf_VAJHU_old<0) && !(phjj_VAJHU_old<0 && pvbf_VAJHU_old<0)) cout << "WARNING: Inconsistent jet presence in nominal Djet ME!" << endl;
				Djet_VAJHU=-1;
			};
			if(phjj_VAJHU_old_up>=0 && pvbf_VAJHU_old_up>=0) Djet_VAJHU_up = pvbf_VAJHU_old_up / ( pvbf_VAJHU_old_up + phjj_VAJHU_old_up ) ;
			else{
				if((phjj_VAJHU_old_up<0 || pvbf_VAJHU_old_up<0) && !(phjj_VAJHU_old_up<0 && pvbf_VAJHU_old_up<0)) cout << "WARNING: Inconsistent jet presence in Djet_up ME!" << endl;
				Djet_VAJHU_up=-1;
			};
			if(phjj_VAJHU_old_dn>=0 && pvbf_VAJHU_old_dn>=0) Djet_VAJHU_dn = pvbf_VAJHU_old_dn / ( pvbf_VAJHU_old_dn + phjj_VAJHU_old_dn ) ;
			else{
				if((phjj_VAJHU_old_dn<0 || pvbf_VAJHU_old_dn<0) && !(phjj_VAJHU_old_dn<0 && pvbf_VAJHU_old_dn<0)) cout << "WARNING: Inconsistent jet presence in Djet_dn ME!" << endl;
				Djet_VAJHU_dn=-1;
			};

// Spin 1 and 2 discriminants
			if ((smp >= kGGSamples || smp == 0 || smp == 11 || (smp >= kGGHSamples && smp<kGGHSamples + 5) || (smp == kGGMCFMSamples + 1 || smp == kGGMCFMSamples + 4)) || isData){
				D_bkg_prodIndep = (p0plus_VAJHU*p0plus_m4l)/(p0plus_VAJHU*p0plus_m4l + bkg_prodIndep_VAMCFM*bkg_m4l);
				D_bkg_prodIndep_ScaleUp   = (p0plus_VAJHU*p0plus_m4l_ScaleUp) / (p0plus_VAJHU*p0plus_m4l_ScaleUp + bkg_prodIndep_VAMCFM*bkg_m4l_ScaleUp);
				D_bkg_prodIndep_ScaleDown = (p0plus_VAJHU*p0plus_m4l_ScaleDown) / (p0plus_VAJHU*p0plus_m4l_ScaleDown + bkg_prodIndep_VAMCFM*bkg_m4l_ScaleDown);
				D_bkg_prodIndep_ResUp     = (p0plus_VAJHU*p0plus_m4l_ResUp) / (p0plus_VAJHU*p0plus_m4l_ResUp + bkg_prodIndep_VAMCFM*bkg_m4l_ResUp);
				D_bkg_prodIndep_ResDown     = (p0plus_VAJHU*p0plus_m4l_ResDown) / (p0plus_VAJHU*p0plus_m4l_ResDown + bkg_prodIndep_VAMCFM*bkg_m4l_ResDown);

				p1plusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p1plusProdIndepVA);
				p1minusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p1minusProdIndepVA);
				p1plusKD = p0plus_VAJHU / (p0plus_VAJHU + p1plus_VAJHU);
				p1minusKD = p0plus_VAJHU / (p0plus_VAJHU + p1_VAJHU);

				graviKD = p0plus_VAJHU / (p0plus_VAJHU + p2minimalVA);
				qqgraviKD = p0plus_VAJHU / (p0plus_VAJHU + p2minimalVA_qqb);
				p2h2plusKD = p0plus_VAJHU / (p0plus_VAJHU + p2h2plusVA);
				p2h2plus_qqb_KD = p0plus_VAJHU / (p0plus_VAJHU + p2h2plusVA_qqb);
				p2h3plusKD = p0plus_VAJHU / (p0plus_VAJHU + p2h3plusVA);
				p2h3plus_qqb_KD = p0plus_VAJHU / (p0plus_VAJHU + p2h3plusVA_qqb);
				p2hplusKD = p0plus_VAJHU / (p0plus_VAJHU + p2hplusVA);
				p2hplus_qqb_KD = p0plus_VAJHU / (p0plus_VAJHU + p2hplusVA_qqb);
				p2bplusKD = p0plus_VAJHU / (p0plus_VAJHU + p2bplusVA);
				p2bplus_qqb_KD = p0plus_VAJHU / (p0plus_VAJHU + p2bplusVA_qqb);
				p2h6plusKD = p0plus_VAJHU / (p0plus_VAJHU + p2h6plusVA);
				p2h6plus_qqb_KD = p0plus_VAJHU / (p0plus_VAJHU + p2h6plusVA_qqb);
				p2h7plusKD = p0plus_VAJHU / (p0plus_VAJHU + p2h7plusVA);
				p2h7plus_qqb_KD = p0plus_VAJHU / (p0plus_VAJHU + p2h7plusVA_qqb);
				p2hminusKD = p0plus_VAJHU / (p0plus_VAJHU + p2hminusVA);
				p2hminus_qqb_KD = p0plus_VAJHU / (p0plus_VAJHU + p2hminusVA_qqb);
				p2h9minusKD = p0plus_VAJHU / (p0plus_VAJHU + p2h9minusVA);
				p2h9minus_qqb_KD = p0plus_VAJHU / (p0plus_VAJHU + p2h9minusVA_qqb);
				p2h10minusKD = p0plus_VAJHU / (p0plus_VAJHU + p2h10minusVA);
				p2h10minus_qqb_KD = p0plus_VAJHU / (p0plus_VAJHU + p2h10minusVA_qqb);
				p2mProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p2mProdIndepVA);
				p2h2plusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p2h2plusProdIndepVA);
				p2h3plusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p2h3plusProdIndepVA);
				p2hplusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p2hplusProdIndepVA);
				p2bplusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p2bplusProdIndepVA);
				p2h6plusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p2h6plusProdIndepVA);
				p2h7plusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p2h7plusProdIndepVA);
				p2hminusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p2hminusProdIndepVA);
				p2h9minusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p2h9minusProdIndepVA);
				p2h10minusProdIndepKD = p0plus_VAJHU / (p0plus_VAJHU + p2h10minusProdIndepVA);
			};

			double cadd_g1g2 = 1.0;
			double cadd_g1g4 = 1.0;
			double cadd_g4ZGGG = 1.0;
			double cadd_ZG_L1 = 2.6;
			if(abs(lepIdOrdered[0])==abs(lepIdOrdered[1]) &&
				abs(lepIdOrdered[0])==abs(lepIdOrdered[2]) &&
				abs(lepIdOrdered[0])==abs(lepIdOrdered[3])){
//			if (folder < 2){
				cadd_g1g2 = pow(1.638,2)/2.3;
				cadd_g1g4 = pow(2.521,2)/7.0;
				cadd_g4ZGGG = 1.0/7.0;
			}
			else{
				cadd_g1g2 = pow(1.638, 2) / 2.1;
				cadd_g1g4 = pow(2.521, 2) / 6.0;
				cadd_g4ZGGG = 1.0/6.0;
			};
			D_g1q2 = p0plus_VAJHU / (p0plus_VAJHU + p0_g1prime2_VAJHU);
			D_g1q2int = pg1g1prime2_VAJHU/(p0plus_VAJHU + p0_g1prime2_VAJHU);
			D_g2 = p0plus_VAJHU/(p0plus_VAJHU + p0hplus_VAJHU);
			D_g4 = p0plus_VAJHU/(p0plus_VAJHU + p0minus_VAJHU);
			D_g2int = pg1g2_VAJHU / (p0plus_VAJHU + p0hplus_VAJHU*cadd_g1g2);
			D_g4int = pg1g4_VAJHU / (p0plus_VAJHU + p0minus_VAJHU*cadd_g1g4);
			D_g2int_perp = pg1g2_pi2_VAJHU / (p0plus_VAJHU + p0hplus_VAJHU*cadd_g1g2);
			D_g4int_perp = pg1g4_pi2_VAJHU / (p0plus_VAJHU + p0minus_VAJHU*cadd_g1g4);

			D_ZG = p0plus_VAJHU / (p0plus_VAJHU + p0Zgs_VAJHU);
			D_GG = p0plus_VAJHU / (p0plus_VAJHU + p0gsgs_VAJHU);
			D_ZGint = pzzzg_VAJHU / (p0plus_VAJHU + p0Zgs_VAJHU);
			D_GGint = pzzgg_VAJHU / (p0plus_VAJHU + p0gsgs_VAJHU);
/*
			D_ZG_PS = p0minus_VAJHU*cadd_g4ZGGG / (p0minus_VAJHU*cadd_g4ZGGG + p0Zgs_PS_VAJHU);
			D_GG_PS = p0minus_VAJHU*cadd_g4ZGGG / (p0minus_VAJHU*cadd_g4ZGGG + p0gsgs_PS_VAJHU);
			D_ZG_PSint = pzzzg_PS_VAJHU / (p0minus_VAJHU*cadd_g4ZGGG + p0Zgs_PS_VAJHU);
			D_GG_PSint = pzzgg_PS_VAJHU / (p0minus_VAJHU*cadd_g4ZGGG + p0gsgs_PS_VAJHU);
*/
			D_ZG_PS = p0plus_VAJHU / (p0plus_VAJHU + p0Zgs_PS_VAJHU);
			D_GG_PS = p0plus_VAJHU / (p0plus_VAJHU + p0gsgs_PS_VAJHU);
			D_ZG_PSint = pzzzg_PS_VAJHU / (p0plus_VAJHU + p0Zgs_PS_VAJHU);
			D_GG_PSint = pzzgg_PS_VAJHU / (p0plus_VAJHU + p0gsgs_PS_VAJHU);

			D_ZG_L1 = p0plus_VAJHU / (p0plus_VAJHU + p0Zgs_g1prime2_VAJHU);
			D_ZG_L1int = pzzzgs_g1prime2_VAJHU*sqrt(cadd_ZG_L1) / (p0plus_VAJHU + p0Zgs_g1prime2_VAJHU*cadd_ZG_L1);
			D_ZG_L1int_perp = pzzzgs_g1prime2_pi2_VAJHU*sqrt(cadd_ZG_L1) / (p0plus_VAJHU + p0Zgs_g1prime2_VAJHU*cadd_ZG_L1);

			mytree->Fill();
		};
		if ((smp < kGGHSamples + 11) || (smp<kGGMCFMSamples + 6 && smp >= kGGMCFMSamples)){
			cout << coutput << endl;
			cout << sum_MC_weight_noK << '\t' << sum_MC_weight_wK << '\t' << sum_MC_weight_wKd << '\t' << sum_MC_weight_wKu << '\t' << sum_MC_weight_wPd << '\t' << sum_MC_weight_wPu << endl;
			cout << K_MH << '\t' << Kd_MH << '\t' << Ku_MH << '\t' << PDFd_MH << '\t' << PDFu_MH << endl;
		};

		foutput->WriteTObject(mytree);
		foutput->Close();
		if(isData){
			for(int g=0;g<tgCTotalBkgSize;g++) delete tgtotalbkg[g];
		};
		if(smp==(kAllSamples-1) && !isData) delete CRtree;
		delete tree;
	};
	fQQBGGProper->Close();
	finput_ctotbkg->Close();
	for(int f=0;f<2;f++) finput_PDFError[f]->Close();
	finput_KDFactor->Close();
	delete tg_interf;
};
