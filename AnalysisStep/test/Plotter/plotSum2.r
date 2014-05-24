void plotSum2() {

  if (! TString(gSystem->GetLibraries()).Contains("Number")) {
    gROOT->LoadMacro("$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/test/Plotter/root_lib/TFileServiceWrapper.cc+");
    gROOT->LoadMacro("$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/interface/Histograms.h");
    gROOT->LoadMacro("$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/test/Plotter/root_lib/Number.cc+");
    gROOT->LoadMacro("$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/test/Plotter/root_lib/XSecReader.cc+");
    gROOT->LoadMacro("$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/test/Plotter/root_lib/LumiNormalization.cc+");
    gROOT->LoadMacro("$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/test/Plotter/root_lib/HistoStack.cc+");
    gROOT->LoadMacro("macros.C");
    gROOT->LoadMacro("stackplot.C");
  }

  gStyle->Reset();
  TStyle * style = getStyle("ZZ");
  //style->SetOptTitle(1);
  //style->SetPadGridX(false);
  //style->SetPadGridY(false);
  //style->SetOptStat(0);

  style->SetMarkerSize(0.8);
  //style->SetTitleOffset(0.85,"X");
  //style->SetPadBottomMargin(0.16);
  //style->SetPadRightMargin(0.05);
  style->cd();

  TString xlabel_M4l = "m_{4l} [GeV]";
  TString xlabel_Mll = "m_{ll} [GeV]";
  TString xlabel_MZ1 = "m_{Z1} [GeV]";
  TString xlabel_MZ2 = "m_{Z2} [GeV]";
  TString xlabel_Mllll = "m_{4l} [GeV]"; //+#gamma
  TString ylabel_M   = "GeV";
  TString ylabel_P   = "GeV";
  TString xlabel_SIP = "SIP3D";
  TString xlabel_ISO = "R_{iso}";
  TString xlabel_CombRelIso2 = "R_{iso,i}+R_{iso,j}";

  f1=TFile::Open("sum_7TeV.root");
  f2=TFile::Open("sum_8TeV.root");

  THStack* s1 = (THStack*) f1->Get("ZZMass/MC");
  TH1F*    d1 = (TH1F*)    f1->Get("ZZMass/data");
  THStack* s2 = (THStack*) f2->Get("ZZMass/MC");
  TH1F*    d2 = (TH1F*)    f2->Get("ZZMass/data");
  TH1F*    dsum = (TH1F*) d1->Clone();
  dsum->Add(d2);
  THStack* ssum = add(s1,s2);
  cout << (long) ssum << " " << (long) dsum <<endl;



  // 10 GeV binning
  TCanvas* c = newCanvas("ZZMass_7Plus8TeV_70-600_10GeV");
  drawStack(ssum, dsum, 20, "drawOverflow", "", xlabel_M4l,ylabel_M,70,600,0,70); //25 for ICHEP

  TCanvas* c = newCanvas("ZZMass_7Plus8TeV_70-800_10GeV");
  drawStack(ssum, dsum, 20, "", "", xlabel_M4l,ylabel_M,70,800,0,70); //25 for ICHEP

  TCanvas* c = newCanvas("ZZMass_7Plus8TeV_70-1000_10GeV");
  drawStack(ssum, dsum, 20, "", "", xlabel_M4l,ylabel_M,70,1000,0,70); //25 for ICHEP


  // 5 GeV binning
  THStack* ssum4 = skipBins(ssum,2.5);
  TH1F* dsum4 = skipBins(dsum,2.5);
  TCanvas* c = newCanvas("ZZMass_7Plus8TeV_70-600_5GeV");
  drawStack(ssum4, dsum4, 10, "drawOverflow" , "", xlabel_M4l,ylabel_M,70,600,0,60); //20 for ICHEP
  //  return;

  //--- zoom version
  THStack* s1 = (THStack*) f1->Get("ZZMass_zoom/MC");
  TH1F*    d1 = (TH1F*)    f1->Get("ZZMass_zoom/data");
  THStack* s2 = (THStack*) f2->Get("ZZMass_zoom/MC");
  TH1F*    d2 = (TH1F*)    f2->Get("ZZMass_zoom/data");
  TH1F*    dsum = (TH1F*) d1->Clone();
  dsum->Add(d2);
  THStack* ssum = add(s1,s2);

  // 2GeV binning starting g from 100.
  THStack* ssum2 = skipBins(ssum,1.);
  TH1F* dsum2 = skipBins(dsum,1.);
  TCanvas* c = newCanvas("ZZMass_7Plus8TeV_70-170_2GeV");
  drawStack(ssum2, dsum2, 4, "", "", xlabel_M4l,ylabel_M,70,170,0,35);

  // 3 GeV binning starting from 100.5...
  THStack* ssum3 = skipBins(ssum,1.5);
  TH1F* dsum3 = skipBins(dsum,1.5);
  TCanvas* c = newCanvas("ZZMass_7Plus8TeV_70-180_3GeV");
  drawStack(ssum3, dsum3, 6, "" , "", xlabel_M4l,ylabel_M,70.5,180,0,35); //13 for ICHEP
  TCanvas* c = newCanvas("ZZMass_7Plus8TeV_100-180_3GeV");
  drawStack(ssum3, dsum3, 6, "" , "", xlabel_M4l,ylabel_M,100.5,180,0,20); //8 for ICHEP


//   drawopt="";
//    TCanvas* c = newCanvas("ZZMass_7Plus8TeV_70-170");
//    drawStack(ssum2, dsum2, 4, drawopt, "", xlabel_M4l,ylabel_M,70,170,0,12);  

//    TCanvas* c = newCanvas("ZZMass_7Plus8TeV_100-170");
//    drawStack(ssum2, dsum2, 4, drawopt, "", xlabel_M4l,ylabel_M,100,170,0,8);

//   //With zero bins
//   drawopt="draw0Bins";
//   TCanvas* c = newCanvas("ZZMass_7Plus8TeV_70-170_zb");
//   drawStack(ssum2, dsum2, 4, drawopt, "", xlabel_M4l,ylabel_M,70,170,0,12);

//   TCanvas* c = newCanvas("ZZMass_7Plus8TeV_100-170_zb");
//   drawStack(ssum2, dsum2, 4, drawopt, "", xlabel_M4l,ylabel_M,100,170,0,8);


  // mela>0.5 plot
  if (false) {
   TCanvas* c = newCanvas("ZZMass_7Plus8TeV_100-180_Mela05");
   drawStack(ssum3, dsum3, 6, "", "", xlabel_M4l,ylabel_M,100.5,180,0,8);
   TPaveText *ll = new TPaveText(0.19, 0.83, 0.39, 0.93, "NDC");
   ll->SetTextSize(0.035);
   ll->SetTextFont(42);
   ll->SetFillColor(0);
   ll->SetBorderSize(0);
   ll->SetMargin(0.01);
   ll->SetTextAlign(12); // align left
   TString text = "MELA>0.5";
   ll->AddText(0.01,0.5,text);
   ll->Draw();
   return;
  }



//   //With zero bins
//   TCanvas* c = newCanvas("ZZMass_7Plus8TeV_70-170_3GeV_zb");
//   drawStack(ssum3, dsum3, 6, "draw0Bins" , "", xlabel_M4l,ylabel_M,70.5,169.5,0,13);

//   TCanvas* c = newCanvas("ZZMass_7Plus8TeV_100-170_3GeV_zb");
//   drawStack(ssum3, dsum3, 6, "draw0Bins" , "", xlabel_M4l,ylabel_M,100.5,169.5,0,8);

// 4 GeV binning - nobody likes this
//   TCanvas* c = newCanvas("ZZMass_7Plus8TeV_100-200_4Gev");
//   drawStack(ssum, dsum, 8, "", "", xlabel_M4l,ylabel_M,100,198,0,12);

//   TCanvas* c = newCanvas("ZZMass_7Plus8TeV_70-200_4Gev");
//   drawStack(ssum, dsum, 8, "", "", xlabel_M4l,ylabel_M,70,198,0,15);


  //--- Z1
  THStack* s1 = (THStack*) f1->Get("Z1/MC");
  TH1F*    d1 = (TH1F*)    f1->Get("Z1/data");
  THStack* s2 = (THStack*) f2->Get("Z1/MC");
  TH1F*    d2 = (TH1F*)    f2->Get("Z1/data");
  TH1F*    dsum = (TH1F*) d1->Clone();
  dsum->Add(d2);
  THStack* ssum = add(s1,s2);

//   THStack* ssum2 = skipBins(ssum);
//   TH1F* dsum2 = skipBins(dsum);

  TCanvas* c = newCanvas("Z1Mass_7Plus8TeV");
  drawStack(ssum, dsum, 4, "", "", xlabel_MZ1,ylabel_M,40,119,0,250); //85 for ICHEP
  
  // Set of plots tuned for 121-131 mass range
  if (true) {
    TCanvas* c = newCanvas("Z1Mass_7Plus8TeV");
    drawStack(ssum, dsum, 4, "", "", xlabel_MZ1,ylabel_M,40,119,0,10);

    TCanvas* c = newCanvas("Z1Mass_7Plus8TeV_3GeV");
    drawStack(ssum, dsum, 6, "", "", xlabel_MZ1,ylabel_M,40,119,0,6);

    TCanvas* c = newCanvas("Z1Mass_7Plus8TeV_4GeV");
    drawStack(ssum, dsum, 8, "", "", xlabel_MZ1,ylabel_M,40,119,0,7);
  }
  

  //--- Z2
  THStack* s1 = (THStack*) f1->Get("Z2/MC");
  TH1F*    d1 = (TH1F*)    f1->Get("Z2/data");
  THStack* s2 = (THStack*) f2->Get("Z2/MC");
  TH1F*    d2 = (TH1F*)    f2->Get("Z2/data");
  TH1F*    dsum = (TH1F*) d1->Clone();
  dsum->Add(d2);
  THStack* ssum = add(s1,s2);

//   THStack* ssum2 = skipBins(ssum);
//   TH1F* dsum2 = skipBins(dsum);

  TCanvas* c = newCanvas("Z2Mass_7Plus8TeV");
  drawStack(ssum, dsum, 4, "", "", xlabel_MZ2,ylabel_M,12.5,119,0,80); //30 for ICHEP

  // Set of plots tuned for 121-131 mass range
  if (true) {
    drawStack(ssum, dsum, 4, "", "", xlabel_MZ2,ylabel_M,12.5,119,0,8);

    TCanvas* c = newCanvas("Z2Mass_7Plus8TeV_1GeV");
    drawStack(ssum, dsum, 2, "", "", xlabel_MZ2,ylabel_M,12,60,0,7);

    TCanvas* c = newCanvas("Z2Mass_7Plus8TeV_3GeV");
    drawStack(ssum, dsum, 6, "", "", xlabel_MZ2,ylabel_M,12.5,119,0,9);

    TCanvas* c = newCanvas("Z2Mass_7Plus8TeV_4GeV");
    drawStack(ssum, dsum, 8, "", "", xlabel_MZ2,ylabel_M,12.5,119,0,10);
  }
  

  //--- LD
  THStack* s1 = (THStack*) f1->Get("LD_lowmass/MC");
  TH1F*    d1 = (TH1F*)    f1->Get("LD_lowmass/data");
  THStack* s2 = (THStack*) f2->Get("LD_lowmass/MC");
  TH1F*    d2 = (TH1F*)    f2->Get("LD_lowmass/data");
  TH1F*    dsum = (TH1F*) d1->Clone();
  dsum->Add(d2);
  THStack* ssum = add(s1,s2);
  TCanvas* c = newCanvas("LD_lowmass_7Plus8TeV");
  drawStack(ssum, dsum, 3, "", "", "D_{bkg}^{kin}","",0,1,0,9); //5 for ICHEP
  
  TCanvas* c = newCanvas("LD_lowmass_7Plus8TeV_ns");
  drawStack(ssum, dsum, 3, "noStack", "", "D_{bkg}^{kin}","",0,1,0,9); //for ICHEP



  if (false) {
    THStack* s1 = (THStack*) f1->Get("LD_himass/MC");
    TH1F*    d1 = (TH1F*)    f1->Get("LD_himass/data");
    THStack* s2 = (THStack*) f2->Get("LD_himass/MC");
    TH1F*    d2 = (TH1F*)    f2->Get("LD_himass/data");
    TH1F*    dsum = (TH1F*) d1->Clone();
    dsum->Add(d2);
    THStack* ssum = add(s1,s2);
    TCanvas* c = newCanvas("LD_himass_7Plus8TeV");
    drawStack(ssum, dsum, 2, "", "", "MELA","",0,1,0,20);
  }

  if (false) {
    THStack* s1 = (THStack*) f1->Get("pseudoLD/MC");
    TH1F*    d1 = (TH1F*)    f1->Get("pseudoLD/data");
    THStack* s2 = (THStack*) f2->Get("pseudoLD/MC");
    TH1F*    d2 = (TH1F*)    f2->Get("pseudoLD/data");
    TH1F*    dsum = (TH1F*) d1->Clone();
    dsum->Add(d2);
    THStack* ssum = add(s1,s2);
    TCanvas* c = newCanvas("pseudoLD_7Plus8TeV");
    drawStack(ssum, dsum, 2, "", "", "pseudo-MELA","",0,1,0,5);    
  }
  
 }
 
