{

  gStyle->SetOptStat(0);

  TFile fSig10("trees/DarkBoson130502/PRODFSR_8TeV/2mu2e/HZZ4lTree_HZV10.root");
  fSig10.cd();

  TFile fHZZ("trees/DarkBoson130502/PRODFSR_8TeV/2mu2e/HZZ4lTree_H126.root");
  fHZZ.cd();

  TFile fZZ("trees/DarkBoson130502/PRODFSR_8TeV/2mu2e/HZZ4lTree_ZZTo2e2mu.root");
  fZZ.cd();

  TFile fggZZ("trees/DarkBoson130502/PRODFSR_8TeV/2mu2e/HZZ4lTree_ggZZ2l2l.root");
  fggZZ.cd();

  TTree *tree_Sig10 = fSig10.Get("SelectedTree");
  TTree *tree_HZZ = fHZZ.Get("SelectedTree");
  TTree *tree_ZZ = fZZ.Get("SelectedTree");
  TTree *tree_ggZZ = fggZZ.Get("SelectedTree");

  Float_t Z2Mass = 0.;
  Float_t MC_weight = 0.;

  tree_Sig10->SetBranchAddress("Z2Mass",&Z2Mass);
  tree_Sig10->SetBranchAddress("MC_weight",&MC_weight);

  tree_HZZ->SetBranchAddress("Z2Mass",&Z2Mass);
  tree_HZZ->SetBranchAddress("MC_weight",&MC_weight);

  tree_ZZ->SetBranchAddress("Z2Mass",&Z2Mass);
  tree_ZZ->SetBranchAddress("MC_weight",&MC_weight);

  tree_ggZZ->SetBranchAddress("Z2Mass",&Z2Mass);
  tree_ggZZ->SetBranchAddress("MC_weight",&MC_weight);

  TH1F hSig10("hSig10","hSig10",100,0.,105.);
  TH1F hHZZ("hHZZ","hHZZ",100,0.,105.);
  TH1F hZZ("hZZ","hZZ",100,0.,105.);
  TH1F hggZZ("hggZZ","hggZZ",100,0.,105.);

  Long64_t nentries = tree_Sig10->GetEntriesFast();
  for(Long64_t jentry=0; jentry<nentries;jentry++) {
    tree_Sig10->GetEntry(jentry);

    hSig10.Fill(Z2Mass,MC_weight);
  }

  nentries = tree_HZZ->GetEntriesFast();
  for(Long64_t jentry=0; jentry<nentries;jentry++) {
    tree_HZZ->GetEntry(jentry);

    hHZZ.Fill(Z2Mass,MC_weight);
  }

  nentries = tree_ZZ->GetEntriesFast();
  for(Long64_t jentry=0; jentry<nentries;jentry++) {
    tree_ZZ->GetEntry(jentry);

    hZZ.Fill(Z2Mass,MC_weight);
  }

  nentries = tree_ggZZ->GetEntriesFast();
  for(Long64_t jentry=0; jentry<nentries;jentry++) {
    tree_ggZZ->GetEntry(jentry);

    hggZZ.Fill(Z2Mass,MC_weight);
  }

  hSig10.SetTitle("");
  hSig10.GetXaxis()->SetTitle("m_{Z_{2}} [GeV/c^{2}]");
  hSig10.SetFillColor(3);

  hHZZ.SetTitle("");
  hHZZ.GetXaxis()->SetTitle("m_{Z_{2}} [GeV/c^{2}]");
  hHZZ.SetFillColor(4);

  hZZ.SetTitle("");
  hZZ.GetXaxis()->SetTitle("m_{Z_{2}} [GeV/c^{2}]");
  hZZ.SetFillColor(2);

  hggZZ.SetTitle("");
  hggZZ.GetXaxis()->SetTitle("m_{Z_{2}} [GeV/c^{2}]");
  hggZZ.SetFillColor(6);

  hSig10.Add(&hHZZ);
  hSig10.Add(&hZZ);
  hSig10.Add(&hggZZ);

  hHZZ.Add(&hZZ);
  hHZZ.Add(&hggZZ);

  hggZZ.Add(&hZZ);

  TLegend legend(0.5,0.6,0.7,0.8,"");
  legend.AddEntry(&hSig10,"V_{d} 10 GeV","f");
  legend.AddEntry(&hHZZ,"Higgs 126 GeV","f");
  legend.AddEntry(&hZZ,"qqZZ","f");
  legend.AddEntry(&hggZZ,"ggZZ","f");

  TCanvas c1;
  c1.cd();
  hSig10.Draw();
  hHZZ.Draw("SAME");
  hggZZ.Draw("SAME");
  hZZ.Draw("SAME");
  legend.Draw();
  c1.SaveAs("mZ2_2mu2e.gif");


}
