{

  TFile fIn("MC_Z2_all.root");
  fIn.cd();

  const Int_t nBins = 60;
  const Double_t mZ2min = 3.;
  const Double_t mZ2max = 70.;
  const Double_t binWidth = (mZ2max - mZ2min)/nBins;

  //output files
  ofstream ofsN("card_Nevt_4mu.txt",fstream::out | fstream::app);

  RooRealVar ZZMass("ZZMass","ZZMass",100.,150.);
  RooRealVar Z2Mass("Z2Mass","Z2Mass",mZ2min,mZ2max);

  RooDataSet *dataMu = (RooDataSet*)fIn.Get("dataMu");
  RooDataSet *data_1 = (RooDataSet*)fIn.Get("data_H4mu");
  RooDataSet *data_2 = (RooDataSet*)fIn.Get("data_ZZ4mu");
  RooDataSet *data_3 = (RooDataSet*)fIn.Get("data_ggZZ4mu");

  RooDataSet *data         = (RooDataSet*)dataMu->reduce(RooFit::SelectVars(RooArgSet(Z2Mass)),RooFit::Cut("ZZMass > 110. && ZZMass < 140."));
  RooDataSet *data_H4mu    = (RooDataSet*)data_1->reduce(RooFit::SelectVars(RooArgSet(Z2Mass)),RooFit::Cut("ZZMass > 110. && ZZMass < 140."));
  RooDataSet *data_ZZ4mu   = (RooDataSet*)data_2->reduce(RooFit::SelectVars(RooArgSet(Z2Mass)),RooFit::Cut("ZZMass > 110. && ZZMass < 140."));
  RooDataSet *data_ggZZ4mu = (RooDataSet*)data_3->reduce(RooFit::SelectVars(RooArgSet(Z2Mass)),RooFit::Cut("ZZMass > 110. && ZZMass < 140."));

  Z2Mass.setBins(nBins);
  ZZMass.setBins(25);

  //Dark boson signal
  RooRealVar mean_S("mean_S","mean_S",mZ2min,mZ2max);
  RooRealVar sigma_S("sigma_S","sigma_S",0.,10.);

  RooGaussian signalPDF("signalPDF","signalPDF",Z2Mass,mean_S,sigma_S);

  //Higgs Signal
  RooDataHist *h_H4mu = (RooDataHist*)data_H4mu->binnedClone();
  RooHistPdf H4muPDF("H4muPDF","H4muPDF",RooArgSet(Z2Mass),*h_H4mu,4);

  //qqZZ signal
  RooDataHist *h_ZZ4mu = (RooDataHist*)data_ZZ4mu->binnedClone();
  RooHistPdf ZZ4muPDF("ZZ4muPDF","ZZ4muPDF",RooArgSet(Z2Mass),*h_ZZ4mu,4);

  //ggZZ signal
  RooDataHist *h_ggZZ4mu = (RooDataHist*)data_ggZZ4mu->binnedClone();
  RooHistPdf ggZZ4muPDF("ggZZ4muPDF","ggZZ4muPDF",RooArgSet(Z2Mass),*h_ggZZ4mu,4);

  RooRealVar N_Dark("N_Dark","N_Dark",0.,100.);
  RooRealVar N_H("N_H","N_H",0.,100.);
  RooRealVar N_bkg("N_bkg","N_bkg",0.,100.);

  RooRealVar frac_bkg("frac_bkg","frac_bkg",0.,1.);
  RooAddPdf PDFbkg("PDFbkg","PDFbkg",RooArgList(ZZ4muPDF,ggZZ4muPDF),RooArgList(frac_bkg));

  RooAddPdf PDFTot("PDFTot","PDFTot",RooArgList(signalPDF,H4muPDF,PDFbkg),RooArgList(N_Dark,N_H,N_bkg));

  for(int i=1;i<nBins;i++){

    RooWorkspace w("w");

    Double_t m_mean = mZ2min + 0.5*binWidth + i*binWidth;
    Double_t m_sigma = 1.;

    mean_S.setVal(m_mean);
    sigma_S.setVal(m_sigma);

    N_Dark.setVal(2.);
    N_H.setVal(1.);
    N_bkg.setVal(5.);
    frac_bkg.setVal(0.8);

    mean_S.setConstant(kTRUE);
    sigma_S.setConstant(kTRUE);

    //PDFTot.fitTo(*data,RooFit::SumW2Error(kFALSE),RooFit::Extended(kTRUE));

    RooPlot *Z2plot = Z2Mass.frame();
    data->plotOn(Z2plot);
    PDFTot.plotOn(Z2plot);

    TCanvas c1;
    c1.cd();
    Z2plot->Draw();
    if(i == 1) c1.SaveAs("fit.ps(");
    else if(i == nBins - 1) c1.SaveAs("fit.ps)");
    else c1.SaveAs("fit.ps");

    ofsN << mean_S.getVal() << " " << N_Dark.getVal() << " " << N_H.getVal() << " " << N_bkg.getVal() << " " << frac_bkg.getVal() << endl;

    w.import(*data);
    w.import(PDFTot);

    char outname[100];
    sprintf(outname,"workspaces/ws_fitMZ2_%i.root",i);
    TFile fOut(outname,"RECREATE");
    fOut.cd();
    w.Write();
    fOut.Close();

  }

    RooPlot *ZZplot = ZZMass.frame();
    dataMu->plotOn(ZZplot);

    TCanvas c2;
    c2.cd();
    ZZplot->Draw();
    c2.SaveAs("fit_ZZMass.gif");

}
