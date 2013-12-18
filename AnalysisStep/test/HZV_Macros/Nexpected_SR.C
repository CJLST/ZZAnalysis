{

  RooRealVar ZZMass("ZZMass","ZZMass",100.,600.);

  RooRealVar meanL("meanL","meanL",140.3,100.,150.);
  RooRealVar sigmaL("sigmaL","sigmaL",21.7,0.,100.);
  RooLandau ZjetsPDF("ZjetsPDF","ZjetsPDF",ZZMass,meanL,sigmaL);

  ZZMass.setRange("SR",110.,140.);
  ZZMass.setRange("full",100.,600.);

  RooAbsReal *Integral_Zbb_SR = ZjetsPDF.createIntegral(RooArgList(ZZMass),"SR");
  RooAbsReal *Integral_Zbb_full = ZjetsPDF.createIntegral(RooArgList(ZZMass),"full");

  //this is 4mu
  double r4mu = Integral_Zbb_SR->getVal()/Integral_Zbb_full->getVal();
  cout << "Ratio of Zbb in signal region for 4mu = " << r4mu << endl;

  //this is 4e
  meanL.setVal(148.9);
  sigmaL.setVal(20.2);
  double r4e = Integral_Zbb_SR->getVal()/Integral_Zbb_full->getVal();
  cout << "Ratio of Zbb in signal region for 4e = " << r4e << endl;

  //this is 2mu2e
  meanL.setVal(146.4);
  sigmaL.setVal(19.6);
  double r2mu2e = Integral_Zbb_SR->getVal()/Integral_Zbb_full->getVal();
  cout << "Ratio of Zbb in signal region for 2mu2e = " << r2mu2e << endl;

  TFile fIn("MC_Z2_all.root");
  fIn.cd();

  //This is qqZZ
  RooDataSet *data_ZZ4mu_tmp = (RooDataSet*)fIn.Get("data_ZZ4mu");
  RooDataSet *data_ZZ4mu = (RooDataSet*)data_ZZ4mu_tmp->reduce(RooFit::SelectVars(RooArgSet(ZZMass)));
  RooDataHist *h_ZZ4mu = (RooDataHist*)data_ZZ4mu->binnedClone();
  RooHistPdf ZZ4muPDF("ZZ4muPDF","ZZ4muPDF",RooArgSet(ZZMass),*h_ZZ4mu,4);

  RooDataSet *data_ZZ4e_tmp = (RooDataSet*)fIn.Get("data_ZZ4e");
  RooDataSet *data_ZZ4e = (RooDataSet*)data_ZZ4e_tmp->reduce(RooFit::SelectVars(RooArgSet(ZZMass)));
  RooDataHist *h_ZZ4e = (RooDataHist*)data_ZZ4e->binnedClone();
  RooHistPdf ZZ4ePDF("ZZ4ePDF","ZZ4ePDF",RooArgSet(ZZMass),*h_ZZ4e,4);

  RooDataSet *data_ZZ2mu2e_tmp = (RooDataSet*)fIn.Get("data_ZZ2mu2e");
  RooDataSet *data_ZZ2mu2e = (RooDataSet*)data_ZZ2mu2e_tmp->reduce(RooFit::SelectVars(RooArgSet(ZZMass)));
  RooDataHist *h_ZZ2mu2e = (RooDataHist*)data_ZZ2mu2e->binnedClone();
  RooHistPdf ZZ2mu2ePDF("ZZ2mu2ePDF","ZZ2mu2ePDF",RooArgSet(ZZMass),*h_ZZ2mu2e,4);

  RooAbsReal *Integral_qqZZ_SR_4mu = ZZ4muPDF.createIntegral(RooArgList(ZZMass),"SR");
  RooAbsReal *Integral_qqZZ_full_4mu = ZZ4muPDF.createIntegral(RooArgList(ZZMass),"full");

  RooAbsReal *Integral_qqZZ_SR_4e = ZZ4ePDF.createIntegral(RooArgList(ZZMass),"SR");
  RooAbsReal *Integral_qqZZ_full_4e = ZZ4ePDF.createIntegral(RooArgList(ZZMass),"full");

  RooAbsReal *Integral_qqZZ_SR_2mu2e = ZZ2mu2ePDF.createIntegral(RooArgList(ZZMass),"SR");
  RooAbsReal *Integral_qqZZ_full_2mu2e = ZZ2mu2ePDF.createIntegral(RooArgList(ZZMass),"full");

  double r4mu_qqZZ = Integral_qqZZ_SR_4mu->getVal()/Integral_qqZZ_full_4mu->getVal();
  double r4e_qqZZ = Integral_qqZZ_SR_4e->getVal()/Integral_qqZZ_full_4e->getVal();
  double r2mu2e_qqZZ = Integral_qqZZ_SR_2mu2e->getVal()/Integral_qqZZ_full_2mu2e->getVal();

  cout << "Ratio of qqZZ in signal region for 4mu = " << r4mu_qqZZ << endl;
  cout << "Ratio of qqZZ in signal region for 4e = " << r4e_qqZZ << endl;
  cout << "Ratio of qqZZ in signal region for 2mu2e = " << r2mu2e_qqZZ << endl;

  //This is ggZZ
  RooDataSet *data_ggZZ4mu_tmp = (RooDataSet*)fIn.Get("data_ggZZ4mu");
  RooDataSet *data_ggZZ4mu = (RooDataSet*)data_ggZZ4mu_tmp->reduce(RooFit::SelectVars(RooArgSet(ZZMass)));
  RooDataHist *h_ggZZ4mu = (RooDataHist*)data_ggZZ4mu->binnedClone();
  RooHistPdf ggZZ4muPDF("ggZZ4muPDF","ggZZ4muPDF",RooArgSet(ZZMass),*h_ggZZ4mu,4);

  RooDataSet *data_ggZZ4e_tmp = (RooDataSet*)fIn.Get("data_ggZZ4e");
  RooDataSet *data_ggZZ4e = (RooDataSet*)data_ggZZ4e_tmp->reduce(RooFit::SelectVars(RooArgSet(ZZMass)));
  RooDataHist *h_ggZZ4e = (RooDataHist*)data_ggZZ4e->binnedClone();
  RooHistPdf ggZZ4ePDF("ggZZ4ePDF","ggZZ4ePDF",RooArgSet(ZZMass),*h_ggZZ4e,4);

  RooDataSet *data_ggZZ2mu2e_tmp = (RooDataSet*)fIn.Get("data_ggZZ2mu2e");
  RooDataSet *data_ggZZ2mu2e = (RooDataSet*)data_ggZZ2mu2e_tmp->reduce(RooFit::SelectVars(RooArgSet(ZZMass)));
  RooDataHist *h_ggZZ2mu2e = (RooDataHist*)data_ggZZ2mu2e->binnedClone();
  RooHistPdf ggZZ2mu2ePDF("ggZZ2mu2ePDF","ggZZ2mu2ePDF",RooArgSet(ZZMass),*h_ggZZ2mu2e,4);

  RooAbsReal *Integral_ggZZ_SR_4mu = ggZZ4muPDF.createIntegral(RooArgList(ZZMass),"SR");
  RooAbsReal *Integral_ggZZ_full_4mu = ggZZ4muPDF.createIntegral(RooArgList(ZZMass),"full");

  RooAbsReal *Integral_ggZZ_SR_4e = ggZZ4ePDF.createIntegral(RooArgList(ZZMass),"SR");
  RooAbsReal *Integral_ggZZ_full_4e = ggZZ4ePDF.createIntegral(RooArgList(ZZMass),"full");

  RooAbsReal *Integral_ggZZ_SR_2mu2e = ggZZ2mu2ePDF.createIntegral(RooArgList(ZZMass),"SR");
  RooAbsReal *Integral_ggZZ_full_2mu2e = ggZZ2mu2ePDF.createIntegral(RooArgList(ZZMass),"full");

  double r4mu_ggZZ = Integral_ggZZ_SR_4mu->getVal()/Integral_ggZZ_full_4mu->getVal();
  double r4e_ggZZ = Integral_ggZZ_SR_4e->getVal()/Integral_ggZZ_full_4e->getVal();
  double r2mu2e_ggZZ = Integral_ggZZ_SR_2mu2e->getVal()/Integral_ggZZ_full_2mu2e->getVal();

  cout << "Ratio of ggZZ in signal region for 4mu = " << r4mu_ggZZ << endl;
  cout << "Ratio of ggZZ in signal region for 4e = " << r4e_ggZZ << endl;
  cout << "Ratio of ggZZ in signal region for 2mu2e = " << r2mu2e_ggZZ << endl;





}
