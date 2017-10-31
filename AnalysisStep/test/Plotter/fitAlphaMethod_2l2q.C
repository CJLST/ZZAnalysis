/* 
 * usage: 
 * -specify parameters at the end of this file
 * -run with:
 *   root -l plotDataVsMC.C++
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Math/DistFunc.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TROOT.h"
#include "TLine.h"
#include "TSpline.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"
#include "TMatrix.h"
#include "TFitResult.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include "tdrstyle.C"
// #include "CMS_lumi.C"
#include "plotUtils.C"
#include "ZZAnalysis/AnalysisStep/src/kFactors.C"

#include <ZZAnalysis/AnalysisStep/src/Category.cc>
#include <ZZAnalysis/AnalysisStep/src/bitops.cc>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>
#include <ZZAnalysis/AnalysisStep/test/Plotter/fit_functions3.C>

using namespace std;

float corrFactor[4] = {1.,1.,1.,1.};

const int nChannels = 8;
const int maxBinsForAlpha = 11;

string channelSPart1[nChannels] = {"resolved","merged","resolved","merged","resolved","merged","resolved","merged"};
string channelSPart2[nChannels] = {"","","btag","btag","vbf","vbf","nobtag","nobtag"};
float minX[nChannels] = {400.,750.,400.,650.,400.,800.,400.,750.};
float maxX[nChannels] = {2400.,2400.,2400.,2400.,2400.,2400.,2400.,2400.};
int binsForAlpha[nChannels][maxBinsForAlpha] = {
                                               {1,2,4,6,8,10,12,14,17,20,23},
                                               {2,4,6,8,10,12,14,17,20,23,-1},
                                               {2,5,8,12,-1,-1,-1,-1,-1,-1,-1},
                                               {2,5,8,10,12,-1,-1,-1,-1,-1,-1},
                                               {2,5,8,12,-1,-1,-1,-1,-1,-1,-1},
                                               {2,5,8,12,-1,-1,-1,-1,-1,-1,-1},
                                               {1,2,4,6,8,10,12,14,17,20,23},
                                               {2,4,6,8,10,12,14,17,20,23,-1}  } ;

int binsForAlpha2[nChannels][maxBinsForAlpha] = {
                                               {1,3,5,7,9,11,13,15,18,21,23},
                                               {1,3,5,7,9,11,13,18,21,23,-1},
                                               {3,6,9,12,-1,-1,-1,-1,-1,-1,-1},
                                               {3,6,9,11,12,-1,-1,-1,-1,-1,-1},
                                               {3,6,9,12,-1,-1,-1,-1,-1,-1,-1},
                                               {3,6,9,12,-1,-1,-1,-1,-1,-1,-1},
                                               {1,3,5,7,9,11,13,15,18,21,23},
                                               {1,3,5,7,9,11,13,15,18,23,-1}  } ;


TH1F* rebinHisto(TH1F* origHisto, int theChannel) {
  TH1F* result = (TH1F*)origHisto->Clone();
  std::vector<float> binEnt;
  std::vector<float> binErr;
  float totBin = 0.;
  float totErr = 0.;
  
  int countBinForAlpha = 0;
  
  for (int i=1; i<=origHisto->GetNbinsX(); i++) {
    totBin += origHisto->GetBinContent(i);
    totErr += (origHisto->GetBinError(i))*(origHisto->GetBinError(i));
    if (i == binsForAlpha[theChannel][countBinForAlpha]) {
      binEnt.push_back(totBin);            totBin = 0;
      binErr.push_back(sqrt(totErr));      totErr = 0;
      // cout << countBinForAlpha << " " << i << " " << binsForAlpha[theChannel][countBinForAlpha] << endl;
      countBinForAlpha++;
    }
  }
  
  countBinForAlpha = 0;
  
  for (int i=1; i<=origHisto->GetNbinsX(); i++) {
    result->SetBinContent(i,binEnt.at(countBinForAlpha));
    result->SetBinError(i,binErr.at(countBinForAlpha));
    if (i == binsForAlpha[theChannel][countBinForAlpha]) countBinForAlpha++;
  }
  
  binEnt.clear();
  binErr.clear();
  totBin = 0.;
  totErr = 0.; 

  return result;
}

void fitAlphaMethod_2l2q(string dirout = "fitAlphaExtended", string theNtuple = "TMVAAndRoofitInputs.root", bool extended = false, bool altBinning = false)
{
  
  float minLim = 450.;
  if (extended) minLim = 400.;
  for (int io = 0; io < 4; io++) {minX[2*io] = minLim;}

  setTDRStyle();
  // gStyle->SetOptStat(1111111);
  TFile ffunc("results_alternat.root","READ");
  TFile ffunc2("results_alternatbin.root","READ");

  string infileName( theNtuple );
  TFile* inFile = TFile::Open( infileName.c_str() );
  TCanvas c1("c1","c1",10,10,650,650);
  TF1* ffit[nChannels];
  TF1* ffitup[nChannels];
  TF1* ffitdown[nChannels];
  TFile outf("results.root","RECREATE");

  for (int tc=0; tc<nChannels; tc++) {
    char histoName[1000];
    sprintf(histoName,"hmass_%sSB%s_Data",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    TH1F* fitHist = (TH1F*)inFile->Get(histoName);
    sprintf(histoName,"hmass_%sSB%s_TT",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    TH1F* ttWZZZHist = (TH1F*)inFile->Get(histoName); 
    sprintf(histoName,"hmass_%sSB%s_DY",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    TH1F* alphaDen = (TH1F*)inFile->Get(histoName);
    sprintf(histoName,"hmass_%sSR%s_DY",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    TH1F* alphaNum = (TH1F*)inFile->Get(histoName);
    sprintf(histoName,"hmass_up_%sSB%s_DY",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    TH1F* alphaDen_up = (TH1F*)inFile->Get(histoName);
    sprintf(histoName,"hmass_up_%sSR%s_DY",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    TH1F* alphaNum_up = (TH1F*)inFile->Get(histoName);
    sprintf(histoName,"hmass_down_%sSB%s_DY",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    TH1F* alphaDen_down = (TH1F*)inFile->Get(histoName);
    sprintf(histoName,"hmass_down_%sSR%s_DY",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    TH1F* alphaNum_down = (TH1F*)inFile->Get(histoName);
    
    sprintf(histoName,"results_%s%s.txt",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    ofstream out(histoName);  

    if (altBinning) {
      for (int tb=0; tb<maxBinsForAlpha; tb++) {
	binsForAlpha[tc][tb] = binsForAlpha2[tc][tb];
      }
    }

    // "rebin" alpha histograms if not enough events
    TH1F* alphaNumRebin = rebinHisto(alphaNum,tc);
    TH1F* alphaDenRebin = rebinHisto(alphaDen,tc);
 
    alphaNumRebin->Divide(alphaNumRebin,alphaDenRebin);
    alphaNumRebin->SetLineColor(4);  
    alphaNumRebin->SetFillColor(4);
    alphaNumRebin->SetFillStyle(0);
    alphaNumRebin->SetMarkerStyle(0);

    // Syst variation: JEC
    TH1F* alphaNumRebin_jecup = rebinHisto(alphaNum_up,tc);
    TH1F* alphaNumRebin_jecdown = rebinHisto(alphaNum_down,tc);  
    
    TH1F* alphaDenRebin_jecup = rebinHisto(alphaDen_up,tc);
    TH1F* alphaDenRebin_jecdown = rebinHisto(alphaDen_down,tc);

    alphaNumRebin_jecup->Divide(alphaNumRebin_jecup,alphaDenRebin_jecup);
    alphaNumRebin_jecup->SetLineColor(2);  
    alphaNumRebin_jecdown->Divide(alphaNumRebin_jecdown,alphaDenRebin_jecdown);
    alphaNumRebin_jecdown->SetLineColor(2); 

    float reset;
    for (int i=1; i<=alphaNumRebin->GetNbinsX(); i++) {
      reset = alphaNumRebin->GetBinContent(i) + sqrt(pow(alphaNumRebin_jecup->GetBinContent(i)-alphaNumRebin->GetBinContent(i),2) + pow(alphaNumRebin->GetBinError(i),2));
      alphaNumRebin_jecup->SetBinContent(i,reset);
      reset = alphaNumRebin->GetBinContent(i) - sqrt(pow(alphaNumRebin_jecdown->GetBinContent(i)-alphaNumRebin->GetBinContent(i),2) + pow(alphaNumRebin->GetBinError(i),2));
      alphaNumRebin_jecdown->SetBinContent(i,reset);
    }
    
    // Syst. variation: binning
    for (int tb=0; tb<maxBinsForAlpha; tb++) {
      binsForAlpha[tc][tb] = binsForAlpha2[tc][tb];
    }

    TH1F* alphaNumRebin_alt = rebinHisto(alphaNum,tc);
    TH1F* alphaDenRebin_alt = rebinHisto(alphaDen,tc);
 
    alphaNumRebin_alt->Divide(alphaNumRebin_alt,alphaDenRebin_alt);
    alphaNumRebin_alt->SetLineColor(3);
    alphaDenRebin_alt->SetLineWidth(3);

    /* for (int i=1; i<=alphaNumRebin->GetNbinsX(); i++) {
      reset = alphaNumRebin->GetBinContent(i) + sqrt(pow(alphaNumRebin_alt->GetBinContent(i)-alphaNumRebin->GetBinContent(i),2) + pow(alphaNumRebin_jecup->GetBinContent(i)-alphaNumRebin->GetBinContent(i),2) + pow(alphaNumRebin->GetBinError(i),2));
      alphaNumRebin_alt->SetBinContent(i,reset);
      reset = alphaNumRebin->GetBinContent(i) - sqrt(pow(alphaNumRebin_alt->GetBinContent(i)-alphaNumRebin->GetBinContent(i),2) + pow(alphaNumRebin_jecdown->GetBinContent(i)-alphaNumRebin->GetBinContent(i),2) + pow(alphaNumRebin->GetBinError(i),2));
      alphaDenRebin_alt->SetBinContent(i,reset);
      }*/
    
    c1.cd();
    
    // TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.35,0.95,0.97);
    // pad1->Draw();
    // TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.02,0.95,0.35);
    // pad2->Draw();

    // pad1->cd();
    alphaNumRebin->SetMinimum(-0.7);
    alphaNumRebin->SetMaximum(1.5);
    
    alphaNumRebin->Draw("e2hist");
    alphaNumRebin_jecup->Draw("histsame");
    alphaNumRebin_jecdown->Draw("histsame");
    alphaNumRebin_alt->Draw("histsame");
    // alphaDenRebin_alt->Draw("histsame");
    alphaNum->Divide(alphaNum,alphaDen);
    alphaNum->Draw("same");

    TLegend *legend = new TLegend(0.35,0.20,0.85,0.35,NULL,"brNDC");
    legend->SetBorderSize(     0);
    legend->SetFillColor (     0);
    legend->SetTextAlign (    12);
    legend->SetTextFont  (    42);
    legend->SetTextSize  (0.03);
	      
    legend->AddEntry(alphaNumRebin, "stat. unc." , "l");
    legend->AddEntry(alphaNumRebin_jecup, "stat. + JEC unc." , "l");
    legend->AddEntry(alphaNumRebin_alt, "alternative binning" , "l");
    legend->Draw("same");

    /* pad2->cd();
    TH1F* ratio = (TH1F*)alphaNum->Clone();
    ratio->Add(alphaNum,alphaNumRebin,1.,-1.);
    ratio->Divide(ratio,alphaNum);
    ratio->SetMinimum(-0.8);
    ratio->SetMaximum(0.8); 
    ratio->GetYaxis()->SetTitle("(alpha-average)/alpha");
    ratio->Draw("e");
    TLine line(ratio->GetXaxis()->GetBinLowEdge(1),0.,
	       ratio->GetXaxis()->GetBinUpEdge(ratio->GetNbinsX()),0.);
    line.SetLineColor(kRed);
    line.SetLineStyle(kDashed);
    line.Draw("same"); */

    sprintf(histoName,"~/www/graviton/%s/alphaAlphaMethod_%s%s.png",dirout.c_str(),channelSPart1[tc].c_str(),channelSPart2[tc].c_str()); 
    c1.SaveAs(histoName);
    
    fitHist->Add(fitHist,ttWZZZHist,1.,-1.);
    fitHist->Multiply(fitHist,alphaNumRebin);
    for (int ib=1; ib<=fitHist->GetNbinsX(); ib++) {
      if (fitHist->GetBinContent(ib) <= 0.1) {
	fitHist->SetBinContent(ib,0.1);
        fitHist->SetBinError(ib,alphaNumRebin->GetBinError(ib)/alphaNumRebin->GetBinContent(ib));
      }
    } 

    sprintf(histoName,"ffit_%sSR%s",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    // if (tc == 1) ffit[tc] = new TF1(histoName,myfunction2,minX[tc],maxX[tc],3); else
    ffit[tc] = new TF1(histoName,myfunction,minX[tc],maxX[tc],6);
    ffit[tc]->SetParNames("constant","slope","constant2","slope2","constant3","slope3");
    ffit[tc]->SetLineColor(2);
    ffit[tc]->SetLineStyle(kDashed);
    ffit[tc]->SetLineWidth(3);
    int bin1 = fitHist->FindBin(minX[tc]); int bin2 = fitHist->FindBin(maxX[tc]);
    ffit[tc]->SetParameter(0,fitHist->Integral(bin1,bin2));
    ffit[tc]->SetParameter(1,0.005);
    ffit[tc]->SetParameter(2,fitHist->Integral(bin1,bin2)/10.);
    ffit[tc]->SetParameter(3,0.01);
    ffit[tc]->SetParameter(4,fitHist->Integral(bin1,bin2)/100.);
    ffit[tc]->SetParameter(5,0.001);
    if (tc%2 != 0 || (!extended && tc==2) || (!extended && tc==4)) {
      ffit[tc]->FixParameter(2,0.);
      ffit[tc]->FixParameter(3,0.01);
      ffit[tc]->FixParameter(4,0.);
      ffit[tc]->FixParameter(5,0.1);
    }      

    TCanvas c2("c2","c2",10,10,600,850);
    c2.cd();

    TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.35,0.95,0.97);
    pad1->Draw();
    TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.02,0.95,0.35);
    pad2->Draw();

    pad1->cd();

    TFitResultPtr frp = fitHist->Fit(histoName,"S","",minX[tc],maxX[tc]);
    TMatrixD cov = frp->GetCovarianceMatrix();
    double* covElem = cov.GetMatrixArray(); 

    if (tc%2 == 0 ) {
      float intLim = (extended ? minLim : minLim+50.);
      out << "The integral of the function in [" << int(intLim) << "-2500] is:" << endl;
      out << "   " << ffit[tc]->Integral(intLim,2500)/(50.*corrFactor[0]) << endl;
    }

    out << "The integral of the function in [700-2500] is:" << endl;
    out << "   " << ffit[tc]->Integral(700,2500)/(50.*corrFactor[0]) << endl;
    out << "The parameters are: " << endl;
    for (int j=0; j<6; j++) 
      out << ffit[tc]->GetParName(j) << "\t" << ffit[tc]->GetParameter(j) << endl;
    out << "The covariance matrix is: " << endl;
    for (int i=0; i<36; i++) {
      out << covElem[i] << "\t";
      if ((i+1)%6 == 0) out << endl;
    }      
    out << endl;  
    out.close();

    // Build upper and lower 1sigma function
    sprintf(histoName,"ffitup_%sSR%s",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    ffitup[tc] = new TF1(histoName,myfunctionErrUp,minX[tc],maxX[tc],42);
    for (int j=0; j<6; j++) {
      ffitup[tc]->FixParameter(j,ffit[tc]->GetParameter(j));
      for (int i=0; i<6; i++) {
	ffitup[tc]->FixParameter(6+j+6*i,covElem[6*j+i]);   // array is filled per column
                                                            // not per row
      }
    }
    ffitup[tc]->SetLineColor(4);
    ffitup[tc]->SetLineStyle(kDashed);
    ffitup[tc]->SetLineWidth(3);  
    ffitup[tc]->Draw("same");

    sprintf(histoName,"ffitdown_%sSR%s",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    ffitdown[tc] = new TF1(histoName,myfunctionErrDown,minX[tc],maxX[tc],42);
    for (int j=0; j<6; j++) {
      ffitdown[tc]->FixParameter(j,ffit[tc]->GetParameter(j));
      for (int i=0; i<6; i++) {
	ffitdown[tc]->FixParameter(6+j+6*i,covElem[4*j+i]);   // array is filled per column
                                                              // not per row
      }
    }
    ffitdown[tc]->SetLineColor(4);
    ffitdown[tc]->SetLineStyle(kDashed);
    ffitdown[tc]->SetLineWidth(3);

    ffitdown[tc]->Draw("same");

    sprintf(histoName,"ffit_%sSR%s",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    TF1* ffit_temp = (TF1*)ffunc.Get(histoName);
    if (ffit_temp) {
      ffit_temp->SetLineColor(kGreen+2);
      ffit_temp->Draw("same");
    }

    TF1* ffit_temp2 = (TF1*)ffunc2.Get(histoName);
    if (ffit_temp2) {
      ffit_temp2->SetLineColor(kOrange-2);
      ffit_temp2->Draw("same");
    }

    TLegend *legend2 = new TLegend(0.55,0.55,0.9,0.9,NULL,"brNDC");
    legend2->SetBorderSize(     0);
    legend2->SetFillColor (     0);
    legend2->SetTextAlign (    12);
    legend2->SetTextFont  (    42);
    legend2->SetTextSize  (0.03);
	      
    legend2->AddEntry(ffit[tc], "nominal" , "l");
    legend2->AddEntry(ffitdown[tc], "#pm 1#sigma parameters" , "l");
    legend2->AddEntry(ffit_temp2, "alternative binning" , "l"); 
    legend2->AddEntry(ffit_temp, "alternative function" , "l");
    legend2->Draw("same");
    
    pad2->cd();
    gPad->SetLogy(0);
    TH1F* funcHist = (TH1F*)ttWZZZHist->Clone();
    TH1F* errHist  = (TH1F*)ttWZZZHist->Clone();
    for (int ib=1; ib<=fitHist->GetNbinsX(); ib++) {
      funcHist->SetBinContent(ib,ffit[tc]->Eval(fitHist->GetBinCenter(ib)));
      errHist->SetBinContent(ib,fitHist->GetBinError(ib));
    }
    TH1F* ratio2 = (TH1F*)funcHist->Clone();
    ratio2->Add(fitHist,funcHist,1.,-1.);
    ratio2->Divide(ratio2,errHist);
    for (int ib=1; ib<=ratio2->GetNbinsX(); ib++) {
      ratio2->SetBinError(ib,1.);
      if (ratio2->GetBinCenter(ib) < minX[tc] || ratio2->GetBinCenter(ib) > maxX[tc]) ratio2->SetBinContent(ib,-999.);
    }
    ratio2->SetMinimum(-5.);
    ratio2->SetMaximum(5.); 
    ratio2->GetYaxis()->SetTitle("(Data-fit)/#sigma_{data}");
    ratio2->Draw("e");
    TLine line2(ratio2->GetXaxis()->GetBinLowEdge(1),0.,
		ratio2->GetXaxis()->GetBinUpEdge(ratio2->GetNbinsX()),0.);
    line2.SetLineColor(kRed);
    line2.SetLineStyle(kDashed);
    line2.Draw("same");  

    sprintf(histoName,"~/www/graviton/%s/bkgshapeAlphaMethod_%s%s.png",dirout.c_str(),channelSPart1[tc].c_str(),channelSPart2[tc].c_str()); 
    c2.SaveAs(histoName);
    if (tc == 0) {
      pad1->cd();
      gPad->SetLogy(1);
      sprintf(histoName,"~/www/graviton/%s/bkgshapeAlphaMethod_%s%s_log.png",dirout.c_str(),channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
      c2.SaveAs(histoName);
      gPad->SetLogy(0);
    }
  }

  outf.cd();
  for (int tc=0; tc<nChannels; tc++) {
    ffit[tc]->Write();
    ffitup[tc]->Write();
    ffitdown[tc]->Write();
  }
  outf.Close();
}
