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
#include <ZZAnalysis/AnalysisStep/test/Plotter/fit_functions_DijetFunc.C>

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


void fitAlphaMethodDijetFunc_2l2q(string dirout = "fitAlphaExtended", string theNtuple = "TMVAAndRoofitInputs.root", bool extended = true, bool altBinning = false)
{
  
  float minLim = 450.;
  if (extended) minLim = 400.;
  for (int io = 0; io < 4; io++) {minX[2*io] = minLim;}

  if (altBinning) {
    for (int tc=0; tc<nChannels; tc++) {
      for (int tb=0; tb<maxBinsForAlpha; tb++) {
	binsForAlpha[tc][tb] = binsForAlpha2[tc][tb];
      }
    }
  }

  setTDRStyle();
  // gStyle->SetOptStat(1111111);

  string infileName( theNtuple );
  TFile* inFile = TFile::Open( infileName.c_str() );
  TCanvas c1("c1","c1",10,10,650,900);
  TF1* ffit[nChannels];
  TF1* ffitup[nChannels];
  TF1* ffitdown[nChannels];
  TFile outf("results_alternat.root","RECREATE");

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
    sprintf(histoName,"results_alternat_%s%s.txt",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    ofstream out(histoName);  

    // "rebin" alpha histograms if not enough events
    TH1F* alphaNumRebin = (TH1F*)alphaNum->Clone();
    TH1F* alphaDenRebin = (TH1F*)alphaDen->Clone();

    std::vector<float> binEnt;
    std::vector<float> binErr;
    float totBin = 0.;
    float totErr = 0.;
    
    int countBinForAlpha = 0;

    for (int i=1; i<=alphaNum->GetNbinsX(); i++) {
      totBin += alphaNum->GetBinContent(i);
      totErr += (alphaNum->GetBinError(i))*(alphaNum->GetBinError(i));
      if (i == binsForAlpha[tc][countBinForAlpha]) {
        binEnt.push_back(totBin);            totBin = 0;
        binErr.push_back(sqrt(totErr));      totErr = 0;
        // cout << countBinForAlpha << " " << i << " " << binsForAlpha[tc][countBinForAlpha] << endl;
	countBinForAlpha++;
      }
    }

    countBinForAlpha = 0;
   
    for (int i=1; i<=alphaNum->GetNbinsX(); i++) {
      alphaNumRebin->SetBinContent(i,binEnt.at(countBinForAlpha));
      alphaNumRebin->SetBinError(i,binErr.at(countBinForAlpha));
      if (i == binsForAlpha[tc][countBinForAlpha]) countBinForAlpha++;
    }
   
    binEnt.clear();
    binErr.clear();
    totBin = 0.;
    totErr = 0.;    
    countBinForAlpha = 0;

    for (int i=1; i<=alphaDen->GetNbinsX(); i++) {
      totBin += alphaDen->GetBinContent(i);
      totErr += (alphaDen->GetBinError(i))*(alphaDen->GetBinError(i));
      if (i == binsForAlpha[tc][countBinForAlpha]) {
        binEnt.push_back(totBin);            totBin = 0;
        binErr.push_back(sqrt(totErr));      totErr = 0;
	countBinForAlpha++;
      }
    }

    countBinForAlpha = 0;
   
    for (int i=1; i<=alphaDen->GetNbinsX(); i++) {
      alphaDenRebin->SetBinContent(i,binEnt.at(countBinForAlpha));
      alphaDenRebin->SetBinError(i,binErr.at(countBinForAlpha));
      if (i == binsForAlpha[tc][countBinForAlpha]) countBinForAlpha++;
    }
    
    // end rebin  
 
    alphaNum->Divide(alphaNum,alphaDen);
    alphaNumRebin->Divide(alphaNumRebin,alphaDenRebin);
    alphaNumRebin->SetLineColor(4);  

    c1.cd();
    c1.cd();
    
    TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.35,0.95,0.97);
    pad1->Draw();
    TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.02,0.95,0.35);
    pad2->Draw();

    pad1->cd();
    alphaNum->SetMinimum(0.);
    alphaNum->SetMaximum(1.);
    
    alphaNum->Draw();
    alphaNumRebin->Draw("histsame");    

    pad2->cd();
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
    line.Draw("same");

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
    ffit[tc] = new TF1(histoName,myfunction,minX[tc],maxX[tc],3);
    ffit[tc]->SetParNames("constant","a","b");
    ffit[tc]->SetLineColor(2);
    ffit[tc]->SetLineStyle(kDashed);
    ffit[tc]->SetLineWidth(3);
    int bin1 = fitHist->FindBin(minX[tc]); int bin2 = fitHist->FindBin(maxX[tc]);
    ffit[tc]->SetParameter(0,fitHist->Integral(bin1,bin2));
    ffit[tc]->SetParameter(1,300.);
    ffit[tc]->SetParameter(2,0.03);
    if (tc%2 != 0 || (!extended && tc==2) || (!extended && tc==4)) {
      ffit[tc]->FixParameter(2,ffit[tc-1]->GetParameter(2));
      if (tc==3 || tc==5) ffit[tc]->FixParameter(1,ffit[tc-1]->GetParameter(1));
    }      

    c1.cd();
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
    for (int j=0; j<3; j++) 
      out << ffit[tc]->GetParName(j) << "\t" << ffit[tc]->GetParameter(j) << endl;
    out << "The covariance matrix is: " << endl;
    for (int i=0; i<9; i++) {
      out << covElem[i] << "\t";
      if ((i+1)%3 == 0) out << endl;
    }      
    out << endl;  
    out.close();

    // Build upper and lower 1sigma function
    sprintf(histoName,"ffitup_%sSR%s",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    ffitup[tc] = new TF1(histoName,myfunctionErrUp,minX[tc],maxX[tc],12);
    for (int j=0; j<3; j++) {
      ffitup[tc]->FixParameter(j,ffit[tc]->GetParameter(j));
      for (int i=0; i<3; i++) {
	ffitup[tc]->FixParameter(3+j+3*i,covElem[3*j+i]);   // array is filled per column
                                                            // not per row
      }
    }
    ffitup[tc]->SetLineColor(4);
    ffitup[tc]->SetLineStyle(kDashed);
    ffitup[tc]->SetLineWidth(3);  
    ffitup[tc]->Draw("same");

    sprintf(histoName,"ffitdown_%sSR%s",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    ffitdown[tc] = new TF1(histoName,myfunctionErrDown,minX[tc],maxX[tc],12);
    for (int j=0; j<3; j++) {
      ffitdown[tc]->FixParameter(j,ffit[tc]->GetParameter(j));
      for (int i=0; i<3; i++) {
	ffitdown[tc]->FixParameter(3+j+3*i,covElem[3*j+i]);   // array is filled per column
                                                              // not per row
      }
    }
    ffitdown[tc]->SetLineColor(4);
    ffitdown[tc]->SetLineStyle(kDashed);
    ffitdown[tc]->SetLineWidth(3);
    
    ffitdown[tc]->Draw("same");
    
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
    c1.SaveAs(histoName);
    if (tc == 0) {
      pad1->cd();
      gPad->SetLogy(1);
      sprintf(histoName,"~/www/graviton/%s/bkgshapeAlphaMethod_%s%s_log.png",dirout.c_str(),channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
      c1.SaveAs(histoName);
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
