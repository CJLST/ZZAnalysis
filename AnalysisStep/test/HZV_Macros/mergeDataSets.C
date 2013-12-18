{

  RooRealVar ZZMass("ZZMass","ZZMass",100.,200.);
  RooRealVar Z1Mass("Z1Mass","Z1Mass",20.,120.);
  RooRealVar Z2Mass("Z2Mass","Z2Mass",0.,120.);
  RooRealVar MC_weight("MC_weight","MC_weight",0.,10.);

  RooArgSet vars;
  vars.add(RooArgSet(ZZMass,Z1Mass,Z2Mass,MC_weight));

  TFile f1("trees/DarkBoson130414/PRODFSR_8TeV/4mu/HZZ4lTree_H126.root");
  f1.cd();
  TTree *selected_1 = f1.Get("SelectedTree");
  RooDataSet data_1("data_H4mu","data",selected_1,vars,"","MC_weight");

  TFile f2("trees/DarkBoson130414/PRODFSR_8TeV/4mu/HZZ4lTree_ZZTo4mu.root");
  f2.cd();
  TTree *selected_2 = f2.Get("SelectedTree");
  RooDataSet data_2("data_ZZ4mu","data",selected_2,vars,"","MC_weight");

  TFile f3("trees/DarkBoson130414/PRODFSR_8TeV/4mu/HZZ4lTree_ggZZ4l.root");
  f3.cd();
  TTree *selected_3 = f3.Get("SelectedTree");
  RooDataSet data_3("data_ggZZ4mu","data",selected_3,vars,"","MC_weight");

  RooDataSet dataMu("dataMu","dataMu",vars,"MC_weight");
  dataMu.append(data_1);
  dataMu.append(data_2);
  dataMu.append(data_3);

  TFile f4("trees/DarkBoson130414/PRODFSR_8TeV/4e/HZZ4lTree_H126.root");
  f4.cd();
  TTree *selected_4 = f4.Get("SelectedTree");
  RooDataSet data_4("data_H4e","data",selected_4,vars,"","MC_weight");

  TFile f5("trees/DarkBoson130414/PRODFSR_8TeV/4e/HZZ4lTree_ZZTo4e.root");
  f5.cd();
  TTree *selected_5 = f5.Get("SelectedTree");
  RooDataSet data_5("data_ZZ4e","data",selected_5,vars,"","MC_weight");

  TFile f6("trees/DarkBoson130414/PRODFSR_8TeV/4e/HZZ4lTree_ggZZ4l.root");
  f6.cd();
  TTree *selected_6 = f6.Get("SelectedTree");
  RooDataSet data_6("data_ggZZ4e","data",selected_6,vars,"","MC_weight");

  RooDataSet dataE("dataE","dataE",vars,"MC_weight");
  dataE.append(data_4);
  dataE.append(data_5);
  dataE.append(data_6);

  TFile f7("trees/DarkBoson130414/PRODFSR_8TeV/2mu2e/HZZ4lTree_H126.root");
  f7.cd();
  TTree *selected_7 = f7.Get("SelectedTree");
  RooDataSet data_7("data_H2mu2e","data",selected_7,vars,"","MC_weight");

  TFile f8("trees/DarkBoson130414/PRODFSR_8TeV/2mu2e/HZZ4lTree_ZZTo2e2mu.root");
  f8.cd();
  TTree *selected_8 = f8.Get("SelectedTree");
  RooDataSet data_8("data_ZZ2mu2e","data",selected_8,vars,"","MC_weight");

  TFile f9("trees/DarkBoson130414/PRODFSR_8TeV/2mu2e/HZZ4lTree_ggZZ2l2l.root");
  f9.cd();
  TTree *selected_9 = f9.Get("SelectedTree");
  RooDataSet data_9("data_ggZZ2mu2e","data",selected_9,vars,"","MC_weight");

  RooDataSet dataMuE("dataMuE","dataMuE",vars,"MC_weight");
  dataMuE.append(data_7);
  dataMuE.append(data_8);
  dataMuE.append(data_9);

  TFile fout("MC_Z2_all.root","RECREATE");
  fout.cd();
  dataMu.Write();
  dataE.Write();
  dataMuE.Write();
  data_1.Write();
  data_2.Write();
  data_3.Write();
  data_4.Write();
  data_5.Write();
  data_6.Write();
  data_7.Write();
  data_8.Write();
  data_9.Write();
  fout.Close();

}
