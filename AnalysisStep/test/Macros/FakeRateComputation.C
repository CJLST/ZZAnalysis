#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"

void FakeRateComputation()
{
	
	gROOT->SetStyle("Plain");
	
    //choose which fake rate computed (only one)
	bool isMuon = false ;
	bool isEle  = true ;
	
	TFile * myfile = new TFile("8TeV/ZZ4lAnalysis_DoubleMu_5213.root");
	TTree * mytree = myfile->Get("CRZLTree/crTree") ;
	
	//4e tree
	std::vector<float>  *ZZMass;
	std::vector<float>  *Z2Mass;
	std::vector<float>  *Lep3Pt;
	std::vector<float>  *Lep3Eta;
	std::vector<float>  *Lep3SIP;
	std::vector<short>  *Lep3isID;
	std::vector<int>  *Lep3LepId;
	std::vector<float>  *Lep3combRelIsoPF;
	float  PFMET =  0;
	ZZMass = 0;
	Z2Mass = 0;
	Lep3Pt = 0;
	Lep3Eta = 0;
	Lep3SIP = 0;
	Lep3LepId = 0;
	Lep3isID = 0;
	Lep3combRelIsoPF = 0;
	PFMET = 0;
	TBranch        *b_ZZMass;   //!
	TBranch        *b_Z2Mass;
	TBranch        *b_Lep3Pt;   //!
	TBranch        *b_Lep3Eta;   //!
	TBranch        *b_Lep3SIP;   //!
	TBranch        *b_Lep3isID;   //!
	TBranch        *b_Lep3LepId;   //!
	TBranch        *b_Lep3combRelIsoPF;   //!
	TBranch        *b_PFMET;   //!
	mytree->SetBranchAddress("ZZMass", &ZZMass,&b_ZZMass);
	mytree->SetBranchAddress("Z2Mass", &Z2Mass,&b_Z2Mass);
	mytree->SetBranchAddress("Lep3Pt", &Lep3Pt,&b_Lep3Pt);
	mytree->SetBranchAddress("Lep3Eta", &Lep3Eta,&b_Lep3Eta);
	mytree->SetBranchAddress("Lep3SIP",&Lep3SIP,&b_Lep3SIP);
	mytree->SetBranchAddress("Lep3LepId",&Lep3LepId,&b_Lep3LepId);
	mytree->SetBranchAddress("Lep3isID",&Lep3isID,&b_Lep3isID);
	mytree->SetBranchAddress("Lep3combRelIsoPF",&Lep3combRelIsoPF,&b_Lep3combRelIsoPF);
	mytree->SetBranchAddress("PFMET",&PFMET,&b_PFMET);
	
	char hname[20];
	char htitle[80];

	TH1F *hpT[2][2];
	double pT_bins[11] = { 3, 5, 7, 10, 15, 20, 25, 30, 40, 50, 80};
	for (int i=0;i<2;i++) {
		for (int j=0;j<2;j++) {
			sprintf(hname,"h%d_%d",i,j);
			sprintf(htitle,"hist for pt: %d-%d",i,j);
			hpT[i][j] = new TH1F(hname,htitle,10,pT_bins);
		    hpT[i][j]->GetXaxis()->SetTitle("p_{T} [GeV]");
			hpT[i][j]->GetYaxis()->SetTitle("Fake Rate");
			hpT[i][j]->Sumw2();
		}
	}
	TH1F *hRatio[2];
	for (int i=0;i<2;i++) {
		sprintf(hname,"h%d",i);
		sprintf(htitle,"hist for Ratio: %d",i);
		hRatio[i] = new TH1F(hname,htitle,10,pT_bins);
		hRatio[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
		hRatio[i]->GetYaxis()->SetTitle("Fake Rate");
		hRatio[i]->GetYaxis()->SetLabelSize(0.10);
		hRatio[i]->GetXaxis()->SetLabelSize(0.10);
		hRatio[i]->Sumw2();
	}
	
	
	
	for (Int_t iEvt = 0; iEvt < mytree->GetEntries() ; ++iEvt) {
		
		mytree->GetEntry(iEvt);
		
		int index = 0 ; 
        if (PFMET < 25) 
		{
			if ( 
				(isMuon && fabs(Lep3Eta->at(index)) < 1.2   && abs(Lep3LepId->at(index)) == 13 ) ||
				(isEle  && fabs(Lep3Eta->at(index)) < 1.499 && abs(Lep3LepId->at(index)) == 11 ) 
			   )
			{
				hpT[0][0]->Fill(Lep3Pt->at(index));
				if (Lep3isID->at(index) > 0.5 && Lep3combRelIsoPF->at(index) < 0.4) 
					hpT[0][1]->Fill(Lep3Pt->at(index));
			}
			if ( 
				(isMuon && fabs(Lep3Eta->at(index)) > 1.2   && abs(Lep3LepId->at(index)) == 13 ) ||
				(isEle  && fabs(Lep3Eta->at(index)) > 1.499 && abs(Lep3LepId->at(index)) == 11 ) 
			   )
				{
				hpT[1][0]->Fill(Lep3Pt->at(index));
				if (Lep3isID->at(index) > 0.5 && Lep3combRelIsoPF->at(index) < 0.4) 
					hpT[1][1]->Fill(Lep3Pt->at(index));
			}
		}
	}
	
	//fake rate histo (num/den)
	for (int j=0;j<2;j++) {
		hRatio[j] = hpT[j][1];
		hRatio[j]->Divide(hpT[j][0]);
		hRatio[j]->SetMaximum(1);
		hRatio[j]->SetLineColor(4);
		hRatio[j]->SetLineWidth(1);
		hRatio[j]->SetMarkerColor(4);
	}

    //plotting
	TLegend *legFake = new TLegend(0.5,0.7,0.75,0.85,NULL,"brNDC");
	legFake->SetTextSize(0.04);
	legFake->SetTextFont(42);
	legFake->SetLineWidth(2);
	legFake->SetTextSize(0.04);
	legFake->SetLineColor(1);
	legFake->SetLineStyle(1);
	legFake->SetLineWidth(1);
	legFake->SetFillColor(0);
	legFake->SetFillStyle(0);
	legFake->Draw();
	entry=legFake->AddEntry("","Barrel","pl");
	entry->SetMarkerColor(4);
	entry->SetLineColor(4);
	entry->SetLineWidth(3);
	entry->SetMarkerSize(1.2);
	entry->SetMarkerStyle(20);
	TLegendEntry *entry=legFake->AddEntry("","Endcap","pl");
	entry->SetMarkerColor(2);
	entry->SetLineColor(2);
	entry->SetLineWidth(3);
	entry->SetMarkerSize(1.2);
	entry->SetMarkerStyle(20);	
	
	TCanvas * canFakepT = new TCanvas("canFakepT","canFakepT");
	hRatio[0]->SetMarkerColor(4);
	hRatio[0]->SetLineColor(4);
	hRatio[0]->Draw();
	hRatio[1]->SetMarkerColor(2);
	hRatio[1]->SetLineColor(2);
	hRatio[1]->Draw("same");
	legFake->Draw();

	for (int j=1;j<11;j++) {
		if (isMuon) cout << "Muon barrel: " << "pT " << pT_bins[j-1] << " , " << hRatio[0]->GetBinContent(j) << " " << endl;
		if (isEle) cout << "Ele   barrel: " << "pT " << pT_bins[j-1] << " , " << hRatio[0]->GetBinContent(j) << " " << endl;
	}	
	
	for (int j=1;j<11;j++) {
		if (isMuon) cout << "Muon endcap: "  << "pT " << pT_bins[j-1] << " , " << hRatio[1]->GetBinContent(j) << " " << endl;
		if (isEle) cout  << "Ele  endcap: " << "pT " << pT_bins[j-1] << " , " << hRatio[1]->GetBinContent(j) << " " << endl;
	}	
	
	
}







