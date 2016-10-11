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
#include <ZZAnalysis/AnalysisStep/test/Plotter/fit_functions.C>

using namespace std;

// float corrFactor[4] = {1.,1.,1.,1.};

const int nChannels = 8;
const int maxBinsForAlpha = 11;

string channelSPart1[nChannels] = {"resolved","merged","resolved","merged","resolved","merged","resolved","merged"};
string channelSPart2[nChannels] = {"","","btag","btag","vbf","vbf","nobtag","nobtag"};
float minX[nChannels] = {400.,750.,400.,650.,400.,800.,400.,750.};
float maxX[nChannels] = {2000.,2000.,2000.,1600.,1600.,1600.,2000.,2000};
int binsForAlpha[nChannels][maxBinsForAlpha] = {
                                               {1,2,4,6,8,10,12,14,17,20,23},
                                               {2,4,6,8,10,12,14,17,20,23,-1},
                                               {2,5,8,12,-1,-1,-1,-1,-1,-1,-1},
                                               {2,5,8,10,12,-1,-1,-1,-1,-1,-1},
                                               {2,5,8,12,-1,-1,-1,-1,-1,-1,-1},
                                               {2,5,8,12,-1,-1,-1,-1,-1,-1,-1},
                                               {1,2,4,6,8,10,12,14,17,20,23},
                                               {2,4,6,8,10,12,14,17,20,23,-1}  } ;

void fitAlphaMethod_2l2q(string dirout = "fitAlphaExtended", string theNtuple = "TMVAAndRoofitInputs.root")
{
  
  setTDRStyle();
  // gStyle->SetOptStat(1111111);

  string infileName( theNtuple );
  TFile* inFile = TFile::Open( infileName.c_str() );
  TCanvas c1;
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
    sprintf(histoName,"results_%s%s.txt",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
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
        cout << countBinForAlpha << " " << i << " " << binsForAlpha[tc][countBinForAlpha] << endl;
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
    alphaNum->SetMinimum(0.);
    alphaNum->SetMaximum(1.);
    
    alphaNum->Draw();
    alphaNumRebin->Draw("histsame");    

    sprintf(histoName,"~/www/graviton/%s/alphaAlphaMethod_%s%s.png",dirout.c_str(),channelSPart1[tc].c_str(),channelSPart2[tc].c_str()); 
    c1.SaveAs(histoName);
    
    fitHist->Add(fitHist,ttWZZZHist,1.,-1.);
    fitHist->Multiply(fitHist,alphaNumRebin);

    sprintf(histoName,"ffit_%sSR%s",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    // if (tc == 1) ffit[tc] = new TF1(histoName,myfunction2,minX[tc],maxX[tc],3); else
    ffit[tc] = new TF1(histoName,myfunction,minX[tc],maxX[tc],4);
    ffit[tc]->SetParNames("constant","slope","constant2","slope2");
    ffit[tc]->SetLineColor(2);
    ffit[tc]->SetLineStyle(kDashed);
    ffit[tc]->SetLineWidth(3);
    int bin1 = fitHist->FindBin(minX[tc]); int bin2 = fitHist->FindBin(maxX[tc]);
    ffit[tc]->SetParameter(0,fitHist->Integral(bin1,bin2));
    ffit[tc]->SetParameter(1,0.005);
    ffit[tc]->SetParameter(2,fitHist->Integral(bin1,bin2)/10.);
    ffit[tc]->SetParameter(3,0.01);
    if (tc%2 != 0) {
      ffit[tc]->FixParameter(2,0.);
      ffit[tc]->FixParameter(3,0.01);
    }      

    c1.cd();
    TFitResultPtr frp = fitHist->Fit(histoName,"S","",minX[tc],maxX[tc]);
    TMatrixD cov = frp->GetCovarianceMatrix();
    double* covElem = cov.GetMatrixArray(); 

    if (tc%2 == 0 ) {
      out << "The integral of the function in [400-2500] is:" << endl;
      out << "   " << ffit[tc]->Integral(400,2500)/(50./* *corrFactor[0]*/) << endl;
    }

    out << "The integral of the function in [600-2500] is:" << endl;
    out << "   " << ffit[tc]->Integral(600,2500)/(50./* *corrFactor[0]*/) << endl;
    out << "The parameters are: " << endl;
    for (int j=0; j<4; j++) 
      out << ffit[tc]->GetParName(j) << "\t" << ffit[tc]->GetParameter(j)/* /corrFactor[j]*/ << endl;
    out << "The covariance matrix is: " << endl;
    for (int i=0; i<16; i++) {
      out << covElem[i] << "\t";
      if ((i+1)%4 == 0) out << endl;
    }      
    out << endl;  
    out.close();

    // Build upper and lower 1sigma function
    sprintf(histoName,"ffitup_%sSR%s",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    ffitup[tc] = new TF1(histoName,myfunctionErrUp,minX[tc],maxX[tc],20);
    for (int j=0; j<4; j++) {
      ffitup[tc]->FixParameter(j,ffit[tc]->GetParameter(j));
      for (int i=0; i<4; i++) {
	ffitup[tc]->FixParameter(4+j+4*i,covElem[4*j+i]);   // array is filled per column
                                                            // not per row
      }
    }
    ffitup[tc]->SetLineColor(4);
    ffitup[tc]->SetLineStyle(kDashed);
    ffitup[tc]->SetLineWidth(3);  
    ffitup[tc]->Draw("same");

    sprintf(histoName,"ffitdown_%sSR%s",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    ffitdown[tc] = new TF1(histoName,myfunctionErrDown,minX[tc],maxX[tc],20);
    for (int j=0; j<4; j++) {
      ffitdown[tc]->FixParameter(j,ffit[tc]->GetParameter(j));
      for (int i=0; i<4; i++) {
	ffitdown[tc]->FixParameter(4+j+4*i,covElem[4*j+i]);   // array is filled per column
                                                              // not per row
      }
    }
    ffitdown[tc]->SetLineColor(4);
    ffitdown[tc]->SetLineStyle(kDashed);
    ffitdown[tc]->SetLineWidth(3);
    
    ffitdown[tc]->Draw("same");
    
    sprintf(histoName,"~/www/graviton/%s/bkgshapeAlphaMethod_%s%s.png",dirout.c_str(),channelSPart1[tc].c_str(),channelSPart2[tc].c_str()); 
    c1.SaveAs(histoName);
    if (tc == 0) {
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
