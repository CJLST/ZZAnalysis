void drawPulls() {

  //   if (! TString(gSystem->GetLibraries()).Contains("DTDetId_cc")) {
  gROOT->LoadMacro("macros.C");
  //   }

  gStyle->Reset();
  TStyle * style = getStyle("ZZ");
  style->cd();
  gStyle->SetOptFit(1);
  gStyle->SetFitFormat("5.4g");
  gStyle->SetFuncColor(2);
  gStyle->SetFuncStyle(1);
  gStyle->SetFuncWidth(1);


  int sqrts = 1; //0: 7 TeV. 1: 8TeV  
  TString channel = "4mu";
  //TString channel = "4e";
  //TString channel = "2e2mu";

  bool doH  = true;
  bool doZZ = false;
 
  
  TChain ch;

  //  TFile* input = TFile::Open(inputDir+(*file));  

  TString inputDir = "root://lxcms00//data3/2012/HZZ_root/061012/PRODFSR";

  vector<TString> files;

  TString sqrtsName;
  
  if (sqrts == 0){
    sqrtsName="_7TeV";
    files.push_back("/ZZ4lAnalysis_H120.root");
    files.push_back("/ZZ4lAnalysis_H130.root");
    files.push_back("/ZZ4lAnalysis_H140.root");
    files.push_back("/ZZ4lAnalysis_H150.root");
    files.push_back("/ZZ4lAnalysis_H160.root");
    files.push_back("/ZZ4lAnalysis_H170.root");
    files.push_back("/ZZ4lAnalysis_H180.root");
    files.push_back("/ZZ4lAnalysis_H190.root");
    files.push_back("/ZZ4lAnalysis_H200.root");
    files.push_back("/ZZ4lAnalysis_H210.root");
    files.push_back("/ZZ4lAnalysis_H220.root");
    files.push_back("/ZZ4lAnalysis_H250.root");
    files.push_back("/ZZ4lAnalysis_H275.root");
    files.push_back("/ZZ4lAnalysis_H300.root");
    files.push_back("/ZZ4lAnalysis_H325.root");
    files.push_back("/ZZ4lAnalysis_H350.root");
    files.push_back("/ZZ4lAnalysis_H400.root");
    files.push_back("/ZZ4lAnalysis_H425.root");
    files.push_back("/ZZ4lAnalysis_H450.root");
    files.push_back("/ZZ4lAnalysis_H475.root");
    files.push_back("/ZZ4lAnalysis_H525.root");
    files.push_back("/ZZ4lAnalysis_H550.root");
    files.push_back("/ZZ4lAnalysis_H575.root");
    files.push_back("/ZZ4lAnalysis_H600.root");
    files.push_back("/ZZ4lAnalysis_H650.root");
    files.push_back("/ZZ4lAnalysis_H700.root");
    files.push_back("/ZZ4lAnalysis_H750.root");
    files.push_back("/ZZ4lAnalysis_H800.root");
    //    files.push_back("/ZZ4lAnalysis_H850.root");
    files.push_back("/ZZ4lAnalysis_H900.root");
    files.push_back("/ZZ4lAnalysis_H950.root");
    //    files.push_back("/ZZ4lAnalysis_H1000.root");
  }

  if (sqrts == 1){
    inputDir+="_8TeV";
    sqrtsName="_8TeV";
    files.push_back("/ZZ4lAnalysis_H115.root");
    files.push_back("/ZZ4lAnalysis_H116.root");
    files.push_back("/ZZ4lAnalysis_H117.root");
    files.push_back("/ZZ4lAnalysis_H118.root");
    files.push_back("/ZZ4lAnalysis_H119.root");
    files.push_back("/ZZ4lAnalysis_H120.root");
    files.push_back("/ZZ4lAnalysis_H121.root");
    files.push_back("/ZZ4lAnalysis_H122.root");
    files.push_back("/ZZ4lAnalysis_H123.root");
    files.push_back("/ZZ4lAnalysis_H124.root");
    files.push_back("/ZZ4lAnalysis_H125.root");
    files.push_back("/ZZ4lAnalysis_H126.root");
    files.push_back("/ZZ4lAnalysis_H127.root");
    files.push_back("/ZZ4lAnalysis_H128.root");
    files.push_back("/ZZ4lAnalysis_H129.root");
    files.push_back("/ZZ4lAnalysis_H130.root");
    files.push_back("/ZZ4lAnalysis_H145.root");
    files.push_back("/ZZ4lAnalysis_H150.root");
    files.push_back("/ZZ4lAnalysis_H180.root");
    files.push_back("/ZZ4lAnalysis_H200.root");
    files.push_back("/ZZ4lAnalysis_H250.root");
    files.push_back("/ZZ4lAnalysis_H300.root");
    files.push_back("/ZZ4lAnalysis_H325.root");
    files.push_back("/ZZ4lAnalysis_H350.root");
    //    files.push_back("/ZZ4lAnalysis_H400.root");
    files.push_back("/ZZ4lAnalysis_H450.root");
    files.push_back("/ZZ4lAnalysis_H500.root");
    //    files.push_back("/ZZ4lAnalysis_H550.root");
    files.push_back("/ZZ4lAnalysis_H600.root");
    files.push_back("/ZZ4lAnalysis_H650.root");
    files.push_back("/ZZ4lAnalysis_H700.root");
    files.push_back("/ZZ4lAnalysis_H750.root");
    files.push_back("/ZZ4lAnalysis_H800.root");
    files.push_back("/ZZ4lAnalysis_H850.root");
    files.push_back("/ZZ4lAnalysis_H900.root");
    files.push_back("/ZZ4lAnalysis_H950.root");
    files.push_back("/ZZ4lAnalysis_H1000.root");
  }
 
  TString chainName = "ZZ";
  ch = new TChain(chainName+channel+"Tree/candTree");
  for (vector<TString>::const_iterator file = files.begin(); file!=files.end(); ++file) {
    //    TFile* input = TFile::Open(inputDir+(*file));
    ch.Add(inputDir+(*file));
    //    input->Close();
  }
  
  //  std::cout<<"Number of entries in the chain = "<<ch.GetEntries()<<std::endl;

//   c = newCanvas();
//   ch.Draw("(ZZMass-genM4l)/ZZMassErr>>h(300,-30,30", "sel>=100");
//   TH1F* h1 = h->Clone();
//   drawGFit(h1, 1, 3, -20, 20);

  TString cname = channel+sqrtsName;


  TH2F *h2d = new TH2F("h2d","Pulls vs generated H mass",1000,1,1001,300,-30,30);
  newCanvas(cname+"_PullvsM4l");
  //Note: x Binning is such that ibin = low edge!
  //  ch.Draw("(ZZMass-GenHMass)/ZZMassErr:ZZMass>>h2d(1000,1,1001,300,-30,30", "sel>=100");
  ch.Draw("(ZZMass[iBC]-GenHMass)/ZZMassErr[iBC]:GenHMass>>h2d","iBC>=0 && ZZsel[iBC]>=100");
  //  gStyle->SetOptTitle(0);
  h2d->SetTitle(cname);
  h2d->GetXaxis()->SetTitle("m_{H} (GeV)");
  h2d->GetYaxis()->SetTitle("Pull");
  pull2D = h2d->Clone();
  pull2D->Draw();
//   cout << (long) pull2D <<endl;

//   newCanvas(channel+"_Pull");
//   ch.Draw("(ZZMass-GenHMass)/ZZMassErr>>h(300,-30,30)", "sel>=100");
//   TH1F* h1 = h->Clone();
//   TF1* f = drawGFit(h1, 1.5, 4, -10, 10);

//   return;
  

//    const int ibin =21;
//    float bins[ibin+1] = {110, 125, 135, 145, 155, 165, 175, 185, 195, 205, 215, 225, 240, 260, 290, 325, 375, 425, 475, 525, 575, 625};

   const int ibin =13;
   float bins[ibin+1] = {110, 135, 145, 165, 185, 240, 260, 310, 375, 425, 475, 525, 575, 625};


  
  float fX[ibin+1];
  float fs[ibin+1];
  float fm[ibin+1];
  float fEx[ibin+1];
  float fEm[ibin+1];
  float fEs[ibin+1];

  for (int i=1;i<ibin+1;++i){
    gStyle->SetOptTitle(1);
    //    cout <<  bins[i]-bins[i-1] << endl;
    TString name = "_pull_";
    TString cname2 = cname+"_pull_"+long(bins[i-1]) + "-" + long(bins[i]);
    newCanvas(cname2);
    hpr = h2d->ProjectionY(cname, bins[i-1],bins[i]);
    hpr->SetTitle(cname2);
    TF1* f = drawGFit(hpr, 1.3, 4, -10, 10);
    float m = f->GetParameter("Mean");
    float s = f->GetParameter("Sigma");  
    cout << long(bins[i-1]) << "-" <<long(bins[i]) << " : "  << m << " \t" << s << endl;
    fX[i-1]  = (bins[i]+bins[i-1])/2.;
    fEx[i-1] = (bins[i]-bins[i-1])/2.;
    fm[i-1]  = m;
    fs[i-1]  = s;
    fEm[i-1] =  f->GetParError(f->GetParNumber("Mean"));
    fEs[i-1] =  f->GetParError(f->GetParNumber("Sigma"));
  }

  newCanvas(cname+"_width");
  gPad->SetGrid(1,1);
  gStyle->SetGridColor(15);
  TGraphErrors* gs = new TGraphErrors(ibin, fX, fs, fEx, fEs);
  gs->SetMaximum(1.4);
  gs->SetMinimum(0.9);
  gs->SetLineColor(kRed);
  gs->SetMarkerColor(kRed);
  gs->SetTitle(cname);
  gs->GetXaxis()->SetTitle("m_{H} (GeV)");
  gs->GetYaxis()->SetTitle("Width");
  gs->Draw("AP");


  newCanvas(cname+"_mean");
  gPad->SetGrid(1,1);
  gStyle->SetGridColor(15);
  TGraphErrors* gm = new TGraphErrors(ibin, fX, fm, fEx, fEm);
  gm->SetMaximum(0.5);
  gm->SetMinimum(-0.5);
  gm->SetTitle(cname);
  gm->GetXaxis()->SetTitle("m_{H} (GeV)");
  gm->GetYaxis()->SetTitle("Mean");
//   gm->SetLineColor();
//   gm->SetMarkerColor();
  gm->Draw("AP");


}
