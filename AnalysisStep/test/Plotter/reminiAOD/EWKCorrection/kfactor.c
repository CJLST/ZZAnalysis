void kfactor(){
	TChain *t = new TChain("ZZTree/candTree"); 
	t->Add("root://lxcms03://data3/Higgs/170222/ZZTo4l/ZZ4lAnalysis.root");
	TH2F *h = new TH2F("h","",60,200,2500,140,0,1.4);
	t->Draw("KFactor_EW_qqZZ:GenHMass>>h");
	TH1F *ewk=(TH1F*)h->ProfileX("kfactor_ewk");
	ewk->Draw();
	TFile *fnew = new TFile("kfactor.root","recreate");
	ewk->Write();
	fnew->Close();
	
}
