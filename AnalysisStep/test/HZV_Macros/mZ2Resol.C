#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

/*
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TSystem.h"

#include "RooGlobalFunc.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooFFTConvPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooCBShape.h"
#include "RooPlot.h"
*/


//----------> SET INPUT VARIABLES in Config.h
//#include "Config.h"
//<----------

using namespace RooFit ;
using namespace std;

//Declaration
void signalFits(int channel, int sqrts);
//float WidthValue(float mHStarWidth);

void mZ2Resol()
{
  //gSystem->Load("../CreateDatacards/CMSSW_5_2_5/lib/slc5_amd64_gcc462/libHiggsAnalysisCombinedLimit.so");
  signalFits_HSM(1,8);
  //signalFits_HSM(2,8);
  //signalFits_HSM(3,8);
  return;
}


void signalFits(int channel, int sqrts)
{

  gSystem->Exec("mkdir -p Z2Figs8TeV_HZV");

  string schannel;
  if (channel == 1) schannel = "4mu";
  if (channel == 2) schannel = "4e";
  if (channel == 3) schannel = "2e2mu";
  cout << "Final state = " << schannel << " and sqrt(s) = " << sqrts << endl;

  const int nPoints=3;
  double masses[nPoints]={10.,20.,30.};

  //Prepare to store all the shape parameters for the mass points
  const int arraySize=200;
  assert(arraySize >= nPoints);

  //Double_t masses[arraySize];
  Double_t a_meanLandau[arraySize];  Double_t a_meanLandau_err[arraySize]; 
  Double_t a_sigmaLandau[arraySize]; Double_t a_sigmaLandau_err[arraySize];
  Double_t a_meanGauss[arraySize]; Double_t a_meanGauss_err[arraySize];
  Double_t a_sigmaGauss[arraySize];     Double_t a_sigmaGauss_err[arraySize];
  Double_t a_gf[arraySize];   Double_t a_gf_err[arraySize];
  Double_t a_gl[arraySize];   Double_t a_gl_err[arraySize];

  Double_t a_fitCovQual[arraySize];
  Double_t a_fitEDM[arraySize];
  Double_t a_fitStatus[arraySize];
  
  char outfile[192];
  sprintf(outfile,"Z2Figs%iTeV_HZV",sqrts);
  
  //Loop over the mass points
  for (int iPoint=0; iPoint<nPoints; iPoint++){
    //Pick the correct mass points and paths
    TString filePath="root://lxcms02//data/Higgs/rootuplesOut/DarkBoson130502/PRODFSR_8TeV/";
    filePath.Append(schannel=="2e2mu"?"2mu2e":schannel);
    TString massname;massname.Form("/HZZ4lTree_HZV%d.root",(int)masses[iPoint]);
    filePath.Append(massname.Data());
    TFile *f = TFile::Open(filePath.Data()) ;
    TTree *tree= (TTree*) f->Get("SelectedTree");
    if(tree==NULL){
      cout << "Impossible to retrieve the tree for mass point " << masses[iPoint]  <<" GeV " << endl;
      abort();
    }

    double start=masses[iPoint]-5*iPoint;
    double end=masses[iPoint]+10;
    //Get variables	
    RooRealVar Z2Mass("Z2Mass","Z2Mass",start,end);
    RooRealVar ZZMass("Z2MassErr","Z2MassErr",0,3.5);
    //RooRealVar MC_weight("MC_weight","MC_weight",-10.,10.);
    //RooRealVar invZZMass("-1*Z2MassErr","invZ2MassErr",-6,0);

    //RooDataSet *set2 = new RooDataSet("data","data", tree, RooArgSet(ZZMass,Z2Mass,MC_weight), "", "MC_weight");
    RooDataSet *set2 = new RooDataSet("data","data", tree, RooArgSet(ZZMass,Z2Mass), "", "");
    RooDataHist *set = (RooDataHist*)set2->binnedClone("datahist","datahist");

    //PDFs
    //Landau
    double endsl=1.;
    //if(start<30)endsl=0.05;
    RooRealVar ml("ml","mean landau",0.3,0.,2.) ;//0.3
    RooRealVar sl("sl","sigma landau",0.01,endsl) ;//0.5
    RooLandau Landau("Landau","Landau",ZZMass,ml,sl) ;

    //CB
    RooRealVar meanCB("meanCB","meanCB",0.,-20.,20.);
    RooRealVar sigmaCB("sigmaCB","sigmaCB",1.,0.01,100.);
    RooRealVar alphaCB("alphaCB","alphaCB",-3.,-10.,10.);
    RooRealVar nCB("nCB","nCB",1.,-10.,10.);
    RooCBShape CBShape("CBShape","crystal ball",ZZMass,meanCB,sigmaCB,alphaCB,nCB);

    //Gauss
    double gsigma=0.087;//+0.01*start;
    double gsigmamin=0.01;//+0.009*start;
    double gsigmamax=0.1;//+0.011*start;
    if(channel==3)gsigmamin=0.045+masses[iPoint]/10000.;
    RooRealVar mg("mg","mg",0.1+0.004*masses[iPoint],0.7) ;//0.15-2
    RooRealVar sg("sg","sg",gsigma,gsigmamin,gsigmamax);//0.01-2
    RooGaussian gauss("gauss","gauss",ZZMass,mg,sg) ;

    //Expo
    RooRealVar ec("ec","ec",-0.1,-1,0.1);
    RooExponential expo("expo","expo",ZZMass,ec);

    //LogNormal
    RooRealVar lns("lns","lns",1,0.5,2);
    RooRealVar lnm("lnm","lnm",0.2,0.1,1);
    RooLognormal lognormal("lognormal","lognormal",ZZMass,lns,lnm);


    //RooRealVar gf("gf","gf",0.2,0,0.5);
    //RooRealVar gl("gl","gl",0.8,0.5,1);

    RooRealVar gf("gf","gf",0.5);
    RooRealVar gl("gl","gl",1);//0.8-1

    //RooFFTConvPdf *sigPDF;
    RooAddPdf *sigPDF;
    //sigPDF = new RooFFTConvPdf("sigPDF","sigPDF",ZZMass, SignalTheor,massRes);
    //sigPDF = new RooFFTConvPdf("sigPDF","sigPDF",ZZMass,gauss,SignalTheor);
    //sigPDF = new RooAddPdf("sigPDF","gauss+SignalTheor",RooArgList(gauss,SignalTheor),RooArgList(gf,gl));
    //sigPDF = new RooAddPdf("sigPDF","gauss+CBShape",RooArgList(gauss,CBShape),RooArgList(gf,gl));
    sigPDF = new RooAddPdf("sigPDF","gauss+expo",RooArgList(gauss,Landau),RooArgList(gf,gl));
    //sigPDF->setBufferFraction(0.2);

    //Fit the shape
    //RooFitResult *fitRes = sigPDF->fitTo(*set,Save(1), SumW2Error(kTRUE));
    RooFitResult *fitRes = Landau->fitTo(*set,Save(1), SumW2Error(kTRUE));
    //RooFitResult *fitRes = gauss->fitTo(*set,Save(1), SumW2Error(kTRUE));
 
    int i=iPoint;
    //Fill parameters
    a_fitEDM[i] = fitRes->edm();
    a_fitCovQual[i] = fitRes->covQual();
    a_fitStatus[i] = fitRes->status();

    a_meanLandau[i]  = ml.getVal();
    a_sigmaLandau[i] = sl.getVal();
    //a_meanGauss[i]  = mg.getVal();
    //a_sigmaGauss[i]  = sg.getVal();
    //a_gf[i]     = gf.getVal();
    //a_gl[i]     = gl.getVal();

    a_meanLandau_err[i]  = ml.getError();
    a_sigmaLandau_err[i] = sl.getError();
    //a_meanGauss_err[i]  = mg.getError();
    //a_sigmaGauss_err[i]  = sg.getError();
    //a_gf_err[i]     = gf.getError();
    //a_gl_err[i]     = gl.getError();
    
    //Plot in the figures directory
    RooPlot *xplot = ZZMass.frame();
    set->plotOn(xplot);
    //sigPDF->plotOn(xplot);
    Landau->plotOn(xplot);

    TCanvas canv;
    canv.cd();
    xplot->Draw();
    
    string tmp_plotFileTitle;
    tmp_plotFileTitle.insert(0,outfile);
    tmp_plotFileTitle += "/fitMass_H";
    char tmp2_plotFileTitle[200];
    sprintf(tmp2_plotFileTitle,"%i_%iTeV_",masses[iPoint],sqrts);
    string plotFileTitle = tmp_plotFileTitle + tmp2_plotFileTitle + schannel + ".gif";

    canv.SaveAs(plotFileTitle.c_str());
  }

  Double_t x_err[arraySize];
  
  TGraph* gr_meanLandau  = new TGraph(nPoints, masses, a_meanLandau);
  TGraph* gr_sigmaLandau = new TGraph(nPoints, masses, a_sigmaLandau);
  //TGraph* gr_meanGauss = new TGraph(nPoints, masses, a_meanGauss);
  //TGraph* gr_sigmaGauss  = new TGraph(nPoints, masses, a_sigmaGauss);
  //TGraph* gr_gf   = new TGraph(nPoints, masses, a_gf);
  //TGraph* gr_gl   = new TGraph(nPoints, masses, a_gl);

  gr_meanLandau->SetMarkerStyle(20);
  gr_sigmaLandau->SetMarkerStyle(20);
  //gr_meanGauss->SetMarkerStyle(20);
  //gr_sigmaGauss->SetMarkerStyle(20);
  //gr_gl->SetMarkerStyle(20);
  //gr_gf->SetMarkerStyle(20);

  gr_meanLandau->Fit("pol1");  
  gr_sigmaLandau->Fit("pol1"); 
  //gr_meanGauss->Fit("pol1");   
  //gr_sigmaGauss->Fit("pol1");  
  //gr_gl->Fit("pol1");	       
  //gr_gf->Fit("pol1");          


  TF1 *fit_meanLandau  = gr_meanLandau->GetListOfFunctions()->First();
  TF1 *fit_sigmaLandau = gr_sigmaLandau->GetListOfFunctions()->First();
  //TF1 *fit_meanGauss = gr_meanGauss->GetListOfFunctions()->First();
  //TF1 *fit_sigmaGauss     = gr_sigmaGauss->GetListOfFunctions()->First();
  //TF1 *fit_gf   = gr_gf->GetListOfFunctions()->First();
  //TF1 *fit_gl   = gr_gl->GetListOfFunctions()->First();

  gr_meanLandau->GetXaxis()->SetTitle("Mean value of the Landau function");
  gr_sigmaLandau->GetXaxis()->SetTitle("Sigma of the Landau function");
  //gr_meanGauss->GetXaxis()->SetTitle("Mean Value of the Gauss function");
  //gr_sigmaGauss->GetXaxis()->SetTitle("Sigma of the Gauss function");
  //gr_gf->GetXaxis()->SetTitle("Gauss fraction");
  //gr_gl->GetXaxis()->SetTitle("Landau fraction");

  gr_meanLandau->SetTitle("");
  gr_sigmaLandau->SetTitle("");
  //gr_meanGauss->SetTitle("");
  //gr_sigmaGauss->SetTitle("");
  //gr_gf->SetTitle("");
  //gr_gl->SetTitle("");


  char tmp_outCardName[200];
  sprintf(tmp_outCardName,"CardFragments/MZ2Err_%s_%iTeV_",schannel.c_str(),sqrts);
  string outCardName = tmp_outCardName + schannel + ".txt";
  ofstream ofsCard(outCardName.c_str(),fstream::out);
  ofsCard << "## MZ2 functions --- no spaces! ##" << endl;
  //  ofsCard << "MZ2Shape Gauss m " << fit_sigmaGauss->GetParameter(0) << "+(" << fit_sigmaGauss->GetParameter(1) << "*@0)+(" << fit_sigmaGauss->GetParameter(2) << "*@0*@0)+(" << fit_sigmaGauss->GetParameter(3) << "*@0*@0*@0)" << endl;
  //  ofsCard << "MZ2Shape Gauss s " << fit_meanGauss->GetParameter(0) << "+(" << fit_meanGauss->GetParameter(1) << "*@0)+(" << fit_meanGauss->GetParameter(2) << "*@0*@0)+(" << fit_meanGauss->GetParameter(3) << "*@0*@0*@0)" << endl;
  //  ofsCard << "MZ2Shape Landau m " << fit_meanLandau->GetParameter(0) << "+(" << fit_meanLandau->GetParameter(1) << "*@0)+(" << fit_meanLandau->GetParameter(2) << "*@0*@0)+(" << fit_meanLandau->GetParameter(3) << "*@0*@0*@0)" << endl;
  //  ofsCard << "MZ2Shape Landau m " << fit_sigmaLandau->GetParameter(0) << "+(" << fit_sigmaLandau->GetParameter(1) << "*@0)+(" << fit_sigmaLandau->GetParameter(2) << "*@0*@0)+(" << fit_sigmaLandau->GetParameter(3) << "*@0*@0*@0)" << endl;
  //  ofsCard << "MZ2Shape gf " << fit_gf->GetParameter(0) << "+(" << fit_gf->GetParameter(1) << "*@0)+(" << fit_gf->GetParameter(2) << "*@0*@0)+(" << fit_gf->GetParameter(3) << "*@0*@0*@0)" << endl;
  //  ofsCard << "MZ2Shape gl " << fit_gl->GetParameter(0) << "+(" << fit_gl->GetParameter(1) << "*@0)+(" << fit_gl->GetParameter(2) << "*@0*@0)+(" << fit_gl->GetParameter(3) << "*@0*@0*@0)" << endl;
  //  ofsCard << endl;

  //  ofsCard << "MZ2Shape Gauss m " << fit_sigmaGauss->GetParameter(0) << "+(" << fit_sigmaGauss->GetParameter(1) << "*@0)"  << endl;
  //  ofsCard << "MZ2Shape Gauss s " << fit_meanGauss->GetParameter(0) << "+(" << fit_meanGauss->GetParameter(1) << "*@0)"  << endl;
  ofsCard << "MZ2Shape Landau m " << fit_meanLandau->GetParameter(0) << "+(" << fit_meanLandau->GetParameter(1) << "*@0)" << endl;
  ofsCard << "MZ2Shape Landau s " << fit_sigmaLandau->GetParameter(0) << "+(" << fit_sigmaLandau->GetParameter(1) << "*@0)" << endl;
  //  ofsCard << "MZ2Shape gf " << fit_gf->GetParameter(0) << endl;
  //  ofsCard << "MZ2Shape gl " << fit_gl->GetParameter(0) << endl;
  ofsCard << endl;


  TCanvas *cg = new TCanvas();
  cg->Divide(1,2);
  cg->cd(1);
  gr_meanLandau->Draw("ALP");
  fit_meanLandau->Draw("SAME");
  cg->cd(2);
  gr_sigmaLandau->Draw("ALP");
  fit_sigmaLandau->Draw("SAME");
  //cg->cd(3);
  //gr_meanGauss->Draw("ALP");
  //fit_meanGauss->Draw("SAME");
  //cg->cd(4);
  //gr_sigmaGauss->Draw("ALP");
  //fit_sigmaGauss->Draw("SAME");
  //cg->cd(5);
  //gr_gf->Draw("ALP");
  //fit_gf->Draw("SAME");
  //cg->cd(6);
  //gr_gl->Draw("ALP");
  //fit_gl->Draw("SAME");
  
}

//The actual job
void signalFits_HSM(int channel, int sqrts)
{

  gSystem->Exec("mkdir -p Z2Figs8TeV_qqLandau");

  string schannel;
  if (channel == 1) schannel = "4mu";
  if (channel == 2) schannel = "4e";
  if (channel == 3) schannel = "2e2mu";
  cout << "Final state = " << schannel << " and sqrt(s) = " << sqrts << endl;

  //Pick the correct mass points and paths
  TString filePath="root://lxcms02//data/Higgs/rootuplesOut/DarkBoson130502/PRODFSR_8TeV/";

  filePath.Append(schannel=="2e2mu"?"2mu2e":schannel);
  //filePath.Append("/HZZ4lTree_H.root");//merge 124-127
  filePath.Append("/HZZ4lTree_ZZTo");
  filePath.Append(schannel);
  filePath.Append(".root");
  //Prepare to store all the shape parameters for the mass points
  const int arraySize=200;
  
  Double_t masses[arraySize];  Double_t massesErr[arraySize];
  Double_t a_meanLandau[arraySize];  Double_t a_meanLandau_err[arraySize]; 
  Double_t a_sigmaLandau[arraySize]; Double_t a_sigmaLandau_err[arraySize];
  Double_t a_meanGauss[arraySize]; Double_t a_meanGauss_err[arraySize];
  Double_t a_sigmaGauss[arraySize];     Double_t a_sigmaGauss_err[arraySize];
  Double_t a_gf[arraySize];   Double_t a_gf_err[arraySize];
  Double_t a_gl[arraySize];   Double_t a_gl_err[arraySize];

  Double_t a_fitCovQual[arraySize];
  Double_t a_fitEDM[arraySize];
  Double_t a_fitStatus[arraySize];
  
  char outfile[192];
  sprintf(outfile,"Z2Figs%iTeV_qqLandau",sqrts);
  
  int step=4;              
  int i=0;

  //Loop over the mass points
  for (int start=4; start<60; start=start+step){
    //    if(start>=20){
    //      step=2;
    if(start==24)step=6;
    if(start>=30)step=15;
    //    }else step=4;
    
    masses[i]=(Double_t)start+(step)/2.;
    massesErr[i]=(step)/2.;
    TFile *f = TFile::Open(filePath.Data()) ;
    TTree *tree= (TTree*) f->Get("SelectedTree");
    if(tree==NULL){
      cout << "Impossible to retrieve the tree for mass point " << start  <<" GeV " << endl;
      abort();
    }

    //Get variables	
    RooRealVar Z2Mass("Z2Mass","Z2Mass",(double)start,(double)start+(double)step);
    RooRealVar HMass("ZZMass","ZZMass",120,130);
    RooRealVar ZZMass("Z2MassErr","Z2MassErr",0,3.5);
    RooRealVar MC_weight("MC_weight","MC_weight",0.,10.);
    RooRealVar invZZMass("-1*Z2MassErr","invZ2MassErr",-6,0);

    RooDataSet *set2 = new RooDataSet("data","data", tree, RooArgSet(ZZMass,Z2Mass,HMass,MC_weight), "", "MC_weight");
    RooDataHist *set = (RooDataHist*)set2->binnedClone("datahist","datahist");

    //PDFs
    //Landau
    double endsl=1.;
    //if(start<30&&channel!=2)endsl=0.05;
    //if(start<20&&channel==3)endsl=0.04;
    //if(start<16&&channel==3)endsl=0.03;
    RooRealVar ml("ml","mean landau",0.3,0,1) ;//0.3
    RooRealVar sl("sl","sigma landau",0.001,endsl) ;//0.5
    RooLandau Landau("Landau","Landau",ZZMass,ml,sl) ;

    //CB
    RooRealVar meanCB("meanCB","meanCB",0.,-20.,20.);
    RooRealVar sigmaCB("sigmaCB","sigmaCB",1.,0.01,100.);
    RooRealVar alphaCB("alphaCB","alphaCB",-3.,-10.,10.);
    RooRealVar nCB("nCB","nCB",1.,-10.,10.);
    RooCBShape CBShape("CBShape","crystal ball",ZZMass,meanCB,sigmaCB,alphaCB,nCB);

    //Gauss
    double gsigma=0.087+0.01*start;
    double gsigmamin=0.05+0.009*start;
    double gsigmamax=0.1+0.011*start;
    double gmeanmin=0;
    RooRealVar mg("mg","mg",0.01,2) ;//0-2
    RooRealVar sg("sg","sg",gsigma,gsigmamin,gsigmamax);//0.01-2
    RooGaussian gauss("gauss","gauss",ZZMass,mg,sg) ;

    //Expo
    RooRealVar ec("ec","ec",-0.1,-1,0.1);
    RooExponential expo("expo","expo",ZZMass,ec);

    //LogNormal
    RooRealVar lns("lns","lns",1,0.5,2);
    RooRealVar lnm("lnm","lnm",0.2,0.1,1);
    RooLognormal lognormal("lognormal","lognormal",ZZMass,lns,lnm);


    //RooRealVar gf("gf","gf",0.2,0,0.5);
    //RooRealVar gl("gl","gl",0.8,0.5,1);

    RooRealVar gf("gf","gf",0.5);
    RooRealVar gl("gl","gl",1);//0.8-1

    //RooFFTConvPdf *sigPDF;
    RooAddPdf *sigPDF;
    //sigPDF = new RooFFTConvPdf("sigPDF","sigPDF",ZZMass, SignalTheor,massRes);
    //sigPDF = new RooFFTConvPdf("sigPDF","sigPDF",ZZMass,gauss,SignalTheor);
    //sigPDF = new RooAddPdf("sigPDF","gauss+SignalTheor",RooArgList(gauss,SignalTheor),RooArgList(gf,gl));
    //sigPDF = new RooAddPdf("sigPDF","gauss+CBShape",RooArgList(gauss,CBShape),RooArgList(gf,gl));
    sigPDF = new RooAddPdf("sigPDF","gauss+expo",RooArgList(gauss,Landau),RooArgList(gf,gl));
    //sigPDF->setBufferFraction(0.2);

    //Fit the shape
    //RooFitResult *fitRes = sigPDF->fitTo(*set,Save(1), SumW2Error(kTRUE));
    RooFitResult *fitRes = Landau->fitTo(*set,Save(1), SumW2Error(kTRUE));


    //Fill parameters
    a_fitEDM[i] = fitRes->edm();
    a_fitCovQual[i] = fitRes->covQual();
    a_fitStatus[i] = fitRes->status();

    a_meanLandau[i]  = ml.getVal();
    a_sigmaLandau[i] = sl.getVal();
    //a_meanGauss[i]  = mg.getVal();
    //a_sigmaGauss[i]  = sg.getVal();
    //a_gf[i]     = gf.getVal();
    //a_gl[i]     = gl.getVal();

    a_meanLandau_err[i]  = ml.getError();
    a_sigmaLandau_err[i] = sl.getError();
    //a_meanGauss_err[i]  = mg.getError();
    //a_sigmaGauss_err[i]  = sg.getError();
    //a_gf_err[i]     = gf.getError();
    //a_gl_err[i]     = gl.getError();

    //Plot in the figures directory
    RooPlot *xplot = ZZMass.frame();
    set->plotOn(xplot);
    //sigPDF->plotOn(xplot);
    Landau->plotOn(xplot);

    TCanvas canv;
    canv.cd();
    xplot->Draw();
    
    string tmp_plotFileTitle;
    tmp_plotFileTitle.insert(0,outfile);
    tmp_plotFileTitle += "/fitMass_H";
    char tmp2_plotFileTitle[200];
    sprintf(tmp2_plotFileTitle,"%i_%iTeV_",start,sqrts);
    string plotFileTitle = tmp_plotFileTitle + tmp2_plotFileTitle + schannel + ".gif";

    canv.SaveAs(plotFileTitle.c_str());
    i++;
  }

  Int_t nPoints=i;
  Double_t x_err[arraySize];

  TGraph* gr_meanLandau  = new TGraph(nPoints, masses, a_meanLandau);
  TGraph* gr_sigmaLandau = new TGraph(nPoints, masses, a_sigmaLandau);
  //TGraph* gr_meanGauss = new TGraph(nPoints, masses, a_meanGauss);
  //TGraph* gr_sigmaGauss  = new TGraph(nPoints, masses, a_sigmaGauss);
  //TGraph* gr_gf   = new TGraph(nPoints, masses, a_gf);
  //TGraph* gr_gl   = new TGraph(nPoints, masses, a_gl);

  gr_meanLandau->SetMarkerStyle(20);
  gr_sigmaLandau->SetMarkerStyle(20);
  //gr_meanGauss->SetMarkerStyle(20);
  //gr_sigmaGauss->SetMarkerStyle(20);
  //gr_gl->SetMarkerStyle(20);
  //gr_gf->SetMarkerStyle(20);

  gr_meanLandau->Fit("pol3");  
  gr_sigmaLandau->Fit("pol3"); 
  //gr_meanGauss->Fit("pol3");   
  //gr_sigmaGauss->Fit("pol3");  
  //gr_gl->Fit("pol3");	       
  //gr_gf->Fit("pol3");          


  TF1 *fit_meanLandau  = gr_meanLandau->GetListOfFunctions()->First();
  TF1 *fit_sigmaLandau = gr_sigmaLandau->GetListOfFunctions()->First();
  //TF1 *fit_meanGauss = gr_meanGauss->GetListOfFunctions()->First();
  //TF1 *fit_sigmaGauss     = gr_sigmaGauss->GetListOfFunctions()->First();
  //TF1 *fit_gf   = gr_gf->GetListOfFunctions()->First();
  //TF1 *fit_gl   = gr_gl->GetListOfFunctions()->First();

  gr_meanLandau->GetXaxis()->SetTitle("Mean value of the Landau function");
  gr_sigmaLandau->GetXaxis()->SetTitle("Sigma of the Landau function");
  // gr_meanGauss->GetXaxis()->SetTitle("Mean Value of the Gauss function");
  //gr_sigmaGauss->GetXaxis()->SetTitle("Sigma of the Gauss function");
  //gr_gf->GetXaxis()->SetTitle("Gauss fraction");
  //gr_gl->GetXaxis()->SetTitle("Landau fraction");

  gr_meanLandau->SetTitle("");
  gr_sigmaLandau->SetTitle("");
  //gr_meanGauss->SetTitle("");
  //gr_sigmaGauss->SetTitle("");
  //gr_gf->SetTitle("");
  //gr_gf->SetTitle("");


  char tmp_outCardName[200];
  sprintf(tmp_outCardName,"CardFragments/MZ2Err_%iTeV_",sqrts);
  string outCardName = tmp_outCardName + schannel + ".txt";
  ofstream ofsCard(outCardName.c_str(),fstream::out);
  ofsCard << "## MZ2 functions --- no spaces! ##" << endl;
  //ofsCard << "MZ2Shape Gauss m " << fit_sigmaGauss->GetParameter(0) << "+(" << fit_sigmaGauss->GetParameter(1) << "*@0)+(" << fit_sigmaGauss->GetParameter(2) << "*@0*@0)+(" << fit_sigmaGauss->GetParameter(3) << "*@0*@0*@0)" << endl;
  //ofsCard << "MZ2Shape Gauss s " << fit_meanGauss->GetParameter(0) << "+(" << fit_meanGauss->GetParameter(1) << "*@0)+(" << fit_meanGauss->GetParameter(2) << "*@0*@0)+(" << fit_meanGauss->GetParameter(3) << "*@0*@0*@0)" << endl;
  ofsCard << "MZ2Shape Landau m " << fit_meanLandau->GetParameter(0) << "+(" << fit_meanLandau->GetParameter(1) << "*@0)+(" << fit_meanLandau->GetParameter(2) << "*@0*@0)+(" << fit_meanLandau->GetParameter(3) << "*@0*@0*@0)" << endl;
  ofsCard << "MZ2Shape Landau s " << fit_sigmaLandau->GetParameter(0) << "+(" << fit_sigmaLandau->GetParameter(1) << "*@0)+(" << fit_sigmaLandau->GetParameter(2) << "*@0*@0)+(" << fit_sigmaLandau->GetParameter(3) << "*@0*@0*@0)" << endl;
  //ofsCard << "MZ2Shape gf " << fit_gf->GetParameter(0) << "+(" << fit_gf->GetParameter(1) << "*@0)+(" << fit_gf->GetParameter(2) << "*@0*@0)+(" << fit_gf->GetParameter(3) << "*@0*@0*@0)" << endl;
  //ofsCard << "MZ2Shape gl " << fit_gl->GetParameter(0) << "+(" << fit_gl->GetParameter(1) << "*@0)+(" << fit_gl->GetParameter(2) << "*@0*@0)+(" << fit_gl->GetParameter(3) << "*@0*@0*@0)" << endl;
  ofsCard << endl;

  TCanvas *cg = new TCanvas();
  cg->Divide(1,2);
  cg->cd(1);
  gr_meanLandau->Draw("ALP");
  fit_meanLandau->Draw("SAME");
  cg->cd(2);
  gr_sigmaLandau->Draw("ALP");
  fit_sigmaLandau->Draw("SAME");
  //cg->cd(3);
  //gr_meanGauss->Draw("ALP");
  //fit_meanGauss->Draw("SAME");
  //cg->cd(4);
  //gr_sigmaGauss->Draw("ALP");
  //fit_sigmaGauss->Draw("SAME");
  //cg->cd(5);
  //gr_gf->Draw("ALP");
  //fit_gf->Draw("SAME");
  //cg->cd(6);
  //gr_gl->Draw("ALP");
  //fit_gl->Draw("SAME");

}
