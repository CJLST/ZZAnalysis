
void plotPU(){
  TFile* fMoriond17 = TFile::Open("puWeightsMoriond17_v2.root");  
  //  TFile* fData2017 = TFile::Open("puWeight_Spring2016MC_to_2017Data_294927-301141.root");
  //  TFile* fData2017 = TFile::Open("puWeight_Spring2016MC_to_2017Data_294927-301997.root");
  TFile* fData2017 = TFile::Open("puWeight_Spring2016MC_to_2017Data_294927-305636.root");
  
  auto hMC = (TH1F*) fMoriond17->Get("MC out-of-the-box");
  auto hData2016 = (TH1F*) fMoriond17->Get("Data");
  auto hData2017 = (TH1F*) fData2017->Get("Data");

  
  hMC->SetMaximum(0.065);
  hMC->SetStats(false);
  hMC->SetFillColor(kBlue-10);
  hData2016->SetMarkerStyle(21);
  hData2016->SetMarkerSize(0.8);

  hData2017->SetMarkerStyle(21);
  hData2017->SetMarkerSize(0.8);
  hData2017->SetMarkerColor(kRed);
  

  hMC->Draw("histo");
  hData2016->Draw("sameP");
  hData2017->Draw("sameP");
  leg = new TLegend(0.56,0.71,0.88,0.88);
  leg->AddEntry(hMC, "Spring2016 MC", "f");
  leg->AddEntry(hData2016,"Data 2016", "p");
  leg->AddEntry(hData2017,"Data 2017, 35.88/fb", "p");
  leg->SetLineColor(kWhite);
  leg->Draw();

}
