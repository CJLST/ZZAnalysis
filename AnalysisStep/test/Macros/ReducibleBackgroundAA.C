#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include <ZZAnalysis/AnalysisStep/interface/bitops.h>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>


//7TeV
double fakeRateEle7TeV(double pt, double eta) 
{
	double factor = 0;
	if (fabs(eta) < 1.449) {
		if  ( pt > 3  && pt < 5)
			factor = 0. ;
		if	( pt > 5  && pt < 7)
			factor = 0. ;
		if	( pt > 7  && pt < 10)
			factor = 0.0196523 ;
		if	( pt > 10  && pt < 15)
			factor = 0.0402844 ;
		if	( pt > 15  && pt < 20)
			factor = 0.0597707 ;
		if	( pt > 20  && pt < 25)
			factor = 0.0596864 ;
		if	( pt > 25  && pt < 30)
			factor = 0.0613398 ;
		if	( pt > 30 && pt < 40)
			factor = 0.0602504 ;
		if	( pt > 40 && pt < 50)
			factor = 0.0549928 ;
		if	( pt > 50 )
			factor = 0.0896104 ;
	}
	
	if (fabs(eta) > 1.449) {
		if  ( pt > 3  && pt < 5)
			factor = 0 ;
		if	( pt > 5  && pt < 7)
			factor = 0.0 ;
		if	( pt > 7  && pt < 10)
			factor = 0.0437974 ;
		if	( pt > 10  && pt < 15)
			factor = 0.0852655 ;
		if	( pt > 15  && pt < 20)
			factor = 0.136328 ;
		if	( pt > 20  && pt < 35)
			factor = 0.139205 ;
		if	( pt > 25  && pt < 30)
			factor = 0.177725 ;
		if	( pt > 30 && pt < 40)
			factor = 0.155844 ;
		if	( pt > 40 && pt < 50)
			factor = 0.140426 ;
		if	( pt > 50 )
			factor = 0.210526 ;
	}
	return factor;
}

//8TeV
double fakeRateEle8TeV(double pt, double eta) 
{
	double factor = 0;
	if (fabs(eta) < 1.449) {
	if  ( pt > 3  && pt < 5)
		factor = 0. ;
	if	( pt > 5  && pt < 7)
		factor = 0. ;
	if	( pt > 7  && pt < 10)
		factor = 0.0235961 ;
	if	( pt > 10  && pt < 15)
		factor = 0.0321762 ;
	if	( pt > 15  && pt < 20)
		factor = 0.0379071 ;
	if	( pt > 20  && pt < 25)
		factor = 0.0495186  ;
	if	( pt > 25  && pt < 30)
		factor = 0.0551827 ;
	if	( pt > 30 && pt < 40)
		factor = 0.048103;
	if	( pt > 40 && pt < 50)
		factor = 0.0503876  ;
	if	( pt > 50 )
		factor = 0.0774487  ;
	}
	
	if (fabs(eta) > 1.449) {
		if  ( pt > 3  && pt < 5)
			factor = 0 ;
		if	( pt > 5  && pt < 7)
			factor = 0 ;
		if	( pt > 7  && pt < 10)
			factor = 0.0355005  ;
		if	( pt > 10  && pt < 15)
			factor = 0.0617876 ;
		if	( pt > 15  && pt < 20)
			factor = 0.0929597 ;
		if	( pt > 20  && pt < 35)
			factor = 0.124096 ;
		if	( pt > 25  && pt < 30)
			factor = 0.104603 ;
		if	( pt > 30 && pt < 40)
			factor = 0.14652 ;
		if	( pt > 40 && pt < 50)
			factor = 0.182836 ;
		if	( pt > 50 )
			factor = 0.190184 ;
	}
	return factor;
}

//7TeV
double fakeRateMu7TeV(double pt, double eta) 
{
	double factor = 0;
	if (fabs(eta) < 1.2) {
		if  ( pt > 3  && pt < 5)
			factor = 0  ;
		if	( pt > 5  && pt < 7)
			factor = 0.132577 ;
		if	( pt > 7  && pt < 10)
			factor = 0.116645 ;
		if	( pt > 10  && pt < 15)
			factor = 0.0992579 ;
		if	( pt > 15  && pt < 20)
			factor = 0.0680412;
		if	( pt > 20  && pt < 35)
			factor = 0.0909091 ;
		if	( pt > 25  && pt < 30)
			factor = 0.103448 ;
		if	( pt > 30 && pt < 40)
			factor = 0.112 ;
		if	( pt > 40 && pt < 50)
			factor = 0.08 ;
		if	( pt > 50 )
			factor = 0.163934 ;
	}
	if (fabs(eta) > 1.2) {
		if  ( pt > 3  && pt < 5)
			factor = 0 ;
		if	( pt > 5  && pt < 7)
			factor = 0.139007 ;
		if	( pt > 7  && pt < 10)
			factor = 0.119403 ;
		if	( pt > 10  && pt < 15)
			factor = 0.106017 ;
		if	( pt > 15  && pt < 20)
			factor = 0.0716981 ;
		if	( pt > 20  && pt < 35)
			factor = 0.118056 ;
		if	( pt > 25  && pt < 30)
			factor = 0.135135 ;
		if	( pt > 30 && pt < 40)
			factor =0.157895;
		if	( pt > 40 && pt < 50)
			factor = 0.272727 ;
		if	( pt > 50 )
			factor = 0.324324;
	}
	return factor;
}

//8TeV
double fakeRateMu8TeV(double pt, double eta) 
{
	double factor = 0;
	if (fabs(eta) < 1.2) {
		if  ( pt > 3  && pt < 5)
			factor = 0  ;
		if	( pt > 5  && pt < 7)
			factor = 0.117108 ;
		if	( pt > 7  && pt < 10)
			factor = 0.107685 ;
		if	( pt > 10  && pt < 15)
			factor = 0.0766689 ;
		if	( pt > 15  && pt < 20)
			factor = 0.0789845 ;
		if	( pt > 20  && pt < 35)
			factor = 0.0773639 ;
		if	( pt > 25  && pt < 30)
			factor = 0.115 ;
		if	( pt > 30 && pt < 40)
			factor = 0.10396 ;
		if	( pt > 40 && pt < 50)
			factor = 0.114583 ;
		if	( pt > 50 )
			factor = 0.121951 ;
	}
	if (fabs(eta) > 1.2) {
		if  ( pt > 3  && pt < 5)
			factor = 0 ;
		if	( pt > 5  && pt < 7)
			factor = 0.133566 ;
		if	( pt > 7  && pt < 10)
			factor = 0.121115 ;
		if	( pt > 10  && pt < 15)
			factor = 0.107735 ;
		if	( pt > 15  && pt < 20)
			factor = 0.0962963 ;
		if	( pt > 20  && pt < 35)
			factor = 0.0663717 ;
		if	( pt > 25  && pt < 30)
			factor = 0.11811 ;
		if	( pt > 30 && pt < 40)
			factor = 0.072 ;
		if	( pt > 40 && pt < 50)
			factor = 0.16129 ;
		if	( pt > 50 )
			factor = 0.157895 ;
	}
	return factor;
	
}

void ReducibleBackgroundAA()
{

	gROOT->SetStyle("Plain");
	
	//decide which energy
	bool is7TeV = true;
	bool is8TeV = false;
	//decide which final state
	bool is2e2e = true;
	bool is2mu2mu = false;
	bool is2mu2e = false;
	bool is2e2mu = false;
	
	//load the tree + need branches
	if (is7TeV) {
		if (is2e2mu || is2e2e) 
			TFile * myfile = new TFile("/* 7TeV DoubleEle root*/") ;
		if (is2mu2e || is2mu2mu) 
			TFile * myfile = new TFile("/* 7TeV DoubleMu root*/") ;
	}
	
	if (is8TeV) {
		if (is2e2mu || is2e2e) 
			TFile * myfile = new TFile("/* 8TeV DoubleEle root*/") ;
		if (is2mu2e || is2mu2mu) 
			TFile * myfile = new TFile("/* 8TeV DoubleMu root*/") ;
	}	

	TTree * mytree = myfile->Get("ZLLCand/crTree") ;
	int finalState =0;
	if (is2e2e) finalState=CREEEEss;
	if (is2mu2mu) finalState=CRMMMMss;
	if (is2mu2e) finalState=CRMMEEss;
	if (is2e2mu) finalState=CREEMMss;
	std::vector<float>  *ZZMass;
	std::vector<float>  *Z2Mass;
	std::vector<float>  *Z1Mass;
	std::vector<float>  *ZZLD;
	std::vector<float>  *Lep3Pt;
	std::vector<float>  *Lep3Eta;
	std::vector<float>  *Lep3SIP;
	std::vector<int>  *Lep3LepId;
	std::vector<float>  *Lep3combRelIsoPF;
	std::vector<float>  *Lep4Pt;
	std::vector<float>  *Lep4Eta;
	std::vector<float>  *Lep4SIP;
	std::vector<int>  *Lep4LepId;
	std::vector<float>  *Lep4combRelIsoPF;
	int CRflag=0;
	ZZMass = 0;
	Z2Mass = 0;
	Z1Mass = 0;
	ZZLD = 0;
	Lep3Pt = 0;
	Lep3Eta = 0;
	Lep3SIP = 0;
	Lep3LepId = 0;
	Lep3combRelIsoPF = 0;
	Lep4Pt = 0;
	Lep4Eta = 0;
	Lep4SIP = 0;
	Lep4LepId = 0;
	Lep4combRelIsoPF = 0;
	TBranch        *b_ZZMass;   //!
	TBranch        *b_Z2Mass;
	TBranch        *b_Z1Mass;
	TBranch        *b_ZZLD;
	TBranch        *b_Lep3Pt;   //!
	TBranch        *b_Lep3Eta;   //!
	TBranch        *b_Lep3SIP;   //!
	TBranch        *b_Lep3LepId;   //!
	TBranch        *b_Lep3combRelIsoPF;   //!
	TBranch        *b_Lep4Pt;   //!
	TBranch        *b_Lep4Eta;   //!
	TBranch        *b_Lep4SIP;   //!
	TBranch        *b_Lep4LepId;   //!
	TBranch        *b_Lep4combRelIsoPF;   //!
	TBranch        *b_CRflag;   //!
	mytree->SetBranchAddress("ZZMass", &ZZMass,&b_ZZMass);
	mytree->SetBranchAddress("Z2Mass", &Z2Mass,&b_Z2Mass);
	mytree->SetBranchAddress("Z1Mass", &Z1Mass,&b_Z1Mass);
	mytree->SetBranchAddress("ZZLD", &ZZLD,&b_ZZLD);
	mytree->SetBranchAddress("Lep3Pt", &Lep3Pt,&b_Lep3Pt);
	mytree->SetBranchAddress("Lep3Eta", &Lep3Eta,&b_Lep3Eta);
	mytree->SetBranchAddress("Lep3SIP",&Lep3SIP,&b_Lep3SIP);
	mytree->SetBranchAddress("Lep3LepId",&Lep3LepId,&b_Lep3LepId);
	mytree->SetBranchAddress("Lep3combRelIsoPF",&Lep3combRelIsoPF,&b_Lep3combRelIsoPF);
	mytree->SetBranchAddress("Lep4Pt", &Lep4Pt,&b_Lep4Pt);
	mytree->SetBranchAddress("Lep4Eta", &Lep4Eta,&b_Lep4Eta);
	mytree->SetBranchAddress("Lep4SIP",&Lep4SIP,&b_Lep4SIP);
	mytree->SetBranchAddress("Lep4LepId",&Lep4LepId,&b_Lep4LepId);
	mytree->SetBranchAddress("Lep4combRelIsoPF",&Lep4combRelIsoPF,&b_Lep4combRelIsoPF);
	mytree->SetBranchAddress("CRflag", &CRflag,&b_CRflag);

	//loop on the data 
	double expectedNumberOfEvents = 0 ;
	double NumberOfEvents = 0 ;
	for (Int_t iEvt = 0; iEvt < mytree->GetEntries() ; ++iEvt) {
		mytree->GetEntry(iEvt);
		if(!CRflag)continue;
		if(!test_bit(CRflag,finalState))continue;
		if (is7TeV) {
			if (is2e2e) {
				expectedNumberOfEvents += 0.9*(fakeRateEle7TeV(Lep3Pt->at(index),Lep3Eta->at(index))*fakeRateEle7TeV(Lep4Pt->at(index),Lep4Eta->at(index))) ;
				++NumberOfEvents;
			}
			if (is2mu2e) {
				expectedNumberOfEvents += 0.9*(fakeRateEle7TeV(Lep3Pt->at(index),Lep3Eta->at(index))*fakeRateEle7TeV(Lep4Pt->at(index),Lep4Eta->at(index))) ;
				++NumberOfEvents;
			}
			if (is2e2mu) {
				expectedNumberOfEvents += 0.9*(fakeRateMu7TeV(Lep3Pt->at(index),Lep3Eta->at(index))*fakeRateMu7TeV(Lep4Pt->at(index),Lep4Eta->at(index))) ;
				++NumberOfEvents;
			}
			if (is2mu2mu) {
				expectedNumberOfEvents += 1.2*(fakeRateMu7TeV(Lep3Pt->at(index),Lep3Eta->at(index))*fakeRateMu7TeV(Lep4Pt->at(index),Lep4Eta->at(index))) ;
				++NumberOfEvents;
			}
		}
		
		if (is8TeV) {
			if (is2e2e) {
				expectedNumberOfEvents += 0.9*(fakeRateEle8TeV(Lep3Pt->at(index),Lep3Eta->at(index))*fakeRateEle8TeV(Lep4Pt->at(index),Lep4Eta->at(index))) ;
				++NumberOfEvents;
			}
			if (is2mu2e) {
				expectedNumberOfEvents += 0.9*(fakeRateEle8TeV(Lep3Pt->at(index),Lep3Eta->at(index))*fakeRateEle8TeV(Lep4Pt->at(index),Lep4Eta->at(index))) ;
				++NumberOfEvents;
			}
			if (is2e2mu) {
				expectedNumberOfEvents += 0.9*(fakeRateMu8TeV(Lep3Pt->at(index),Lep3Eta->at(index))*fakeRateMu8TeV(Lep4Pt->at(index),Lep4Eta->at(index))) ;
				++NumberOfEvents;
			}
			if (is2mu2mu) {
				expectedNumberOfEvents += 1.2*(fakeRateMu8TeV(Lep3Pt->at(index),Lep3Eta->at(index))*fakeRateMu8TeV(Lep4Pt->at(index),Lep4Eta->at(index))) ;
				++NumberOfEvents;
			}
		}
		
		
	}

	//output 
	if (is7TeV) 
		cout << "--- 7TeV --- " <<  endl; 
	if (is8TeV) 
		cout << "--- 8TeV --- " <<  endl; 
    if (is2e2e) {
		cout 
		<< "Z1->ee + ee : " <<  expectedNumberOfEvents 
		<< " +/- " <<  expectedNumberOfEvents/sqrt(NumberOfEvents)<< " (stat., evt: " << NumberOfEvents << ")" 
		<< " +/- " <<  expectedNumberOfEvents*0.50<< " (syt.)" 
		<< endl ;
	}
    if (is2e2mu) {
		cout 
		<< "Z1->ee + mumu : " <<  expectedNumberOfEvents 
		<< " +/- " <<  expectedNumberOfEvents/sqrt(NumberOfEvents)<< " (stat., evt: " << NumberOfEvents << ")" 
		<< " +/- " <<  expectedNumberOfEvents*0.50<< " (syt.)" 
		<< endl ;
	}
    if (is2mu2e) {
		cout 
		<< "Z1->mumu + ee : " <<  expectedNumberOfEvents 
		<< " +/- " <<  expectedNumberOfEvents/sqrt(NumberOfEvents)<< " (stat., evt: " << NumberOfEvents << ")" 
		<< " +/- " <<  expectedNumberOfEvents*0.50<< " (syt.)" 
		<< endl ;
	}
    if (is2mu2mu) {
		cout 
		<< "Z1->mumu + mumu : " <<  expectedNumberOfEvents 
		<< " +/- " <<  expectedNumberOfEvents/sqrt(NumberOfEvents)<< " (stat., evt: " << NumberOfEvents << ")" 
		<< " +/- " <<  expectedNumberOfEvents*0.50<< " (syt.)" 
		<< endl ;
	}
	
}







