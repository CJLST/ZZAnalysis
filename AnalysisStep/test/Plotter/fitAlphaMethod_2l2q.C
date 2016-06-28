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
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include "tdrstyle.C"
// #include "CMS_lumi.C"
#include "plotUtils.C"
#include "ZZAnalysis/AnalysisStep/src/kFactors.C"

#include <ZZAnalysis/AnalysisStep/src/Category.cc>
#include <ZZAnalysis/AnalysisStep/src/bitops.cc>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>

using namespace std;

const int nChannels = 6;
string channelSPart1[nChannels] = {"resolved","merged","resolved","merged","resolved","merged"};
string channelSPart2[nChannels] = {"","","btag","btag","vbf","vbf"};
float minX[nChannels] = {450.,500.,450.,600.,450.,600.};
float maxX[nChannels] = {2000.,1600.,1200.,1200.,900.,900.};
int binsForAlpha[nChannels] = {1,3,4,10,4,4};

Double_t myfunction(Double_t *x, Double_t *par)
{
   Float_t xx =x[0];
   Double_t f = par[0]*exp(-par[1]*xx);
   return f;
}

Double_t myfunction2(Double_t *x, Double_t *par)
{
   Float_t xx =x[0];
   Double_t f = (par[0]+par[2]*xx)*exp(-par[1]*xx);
   return f;
}

void fitAlphaMethod_2l2q(string dirout = "fitAlpha", string theNtuple = "TMVAAndRoofitInputs.root")
{
  
  setTDRStyle();
  // gStyle->SetOptStat(1111111);

  string infileName( theNtuple );
  TFile* inFile = TFile::Open( infileName.c_str() );
  TCanvas c1;
  TF1* ffit[nChannels];
  TF1* ffitup[nChannels];
  TF1* ffitdown[nChannels];

  for (int tc=0; tc<nChannels; tc++) {
    char histoName[100];
    sprintf(histoName,"hmass_%sSB%s_Data",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    TH1F* fitHist = (TH1F*)inFile->Get(histoName);
    sprintf(histoName,"hmass_%sSB%s_TT",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    TH1F* ttWZZZHist = (TH1F*)inFile->Get(histoName); 
    sprintf(histoName,"hmass_%sSB%s_DY",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    TH1F* alphaDen = (TH1F*)inFile->Get(histoName);
    sprintf(histoName,"hmass_%sSR%s_Data",channelSPart1[tc].c_str(),channelSPart2[tc].c_str());
    TH1F* alphaNum = (TH1F*)inFile->Get(histoName);

    // "rebin" alpha histograms if not enough events
    TH1F* alphaNumRebin = (TH1F*)alphaNum->Clone();
    TH1F* alphaDenRebin = (TH1F*)alphaDen->Clone();

    std::vector<float> binEnt;
    std::vector<float> binErr;
    float totBin = 0.;
    float totErr = 0.;
    
    for (int i=1; i<=alphaNum->GetNbinsX(); i++) {
      totBin += alphaNum->GetBinContent(i);
      totErr += (alphaNum->GetBinError(i))*(alphaNum->GetBinError(i));
      if (i%binsForAlpha[tc] == 0) {
        binEnt.push_back(totBin);            totBin = 0;
        binErr.push_back(sqrt(totErr));      totErr = 0;
      }
    }
    
    for (int j=0; j<(int)binEnt.size(); j++) {
      for (int i=j*binsForAlpha[tc]+1; i<=(j+1)*binsForAlpha[tc]; i++) {
        alphaNumRebin->SetBinContent(i,binEnt.at(j));
        alphaNumRebin->SetBinError(i,binErr.at(j));
      }
    }

    binEnt.clear();
    binErr.clear();
    totBin = 0.;
    totErr = 0.;
    
    for (int i=1; i<=alphaDen->GetNbinsX(); i++) {
      totBin += alphaDen->GetBinContent(i);
      totErr += (alphaDen->GetBinError(i))*(alphaDen->GetBinError(i));
      if (i%binsForAlpha[tc] == 0) {
        binEnt.push_back(totBin);            totBin = 0;
        binErr.push_back(sqrt(totErr));      totErr = 0;
      }
    }
    
    for (int j=0; j<(int)binEnt.size(); j++) {
      for (int i=j*binsForAlpha[tc]+1; i<=(j+1)*binsForAlpha[tc]; i++) {
        alphaDenRebin->SetBinContent(i,binEnt.at(j));
        alphaDenRebin->SetBinError(i,binErr.at(j));
      }
    }
    // end rebin  
 
    alphaNum->Divide(alphaNum,alphaDen);
    alphaNumRebin->Divide(alphaNumRebin,alphaDenRebin);
    alphaNumRebin->SetLineColor(4);  

    c1.cd();
    alphaNum->Draw();
    alphaNumRebin->Draw("histsame");    

    sprintf(histoName,"~/www/graviton/%s/alphaAlphaMethod_%s%s.png",dirout.c_str(),channelSPart1[tc].c_str(),channelSPart2[tc].c_str()); 
    c1.SaveAs(histoName);
    
    fitHist->Add(fitHist,ttWZZZHist,1.,-1.);
    fitHist->Multiply(fitHist,alphaDenRebin);

    sprintf(histoName,"ffit%d",tc);
    if (tc == 1) ffit[tc] = new TF1(histoName,myfunction2,minX[tc],maxX[tc],3);
    else ffit[tc] = new TF1(histoName,myfunction,minX[tc],maxX[tc],2);
    ffit[tc]->SetParNames("constant","slope");
    ffit[tc]->SetLineColor(2);
    ffit[tc]->SetLineStyle(kDashed);
    ffit[tc]->SetLineWidth(3);
    ffit[tc]->SetParameter(0,fitHist->Integral());
    ffit[tc]->SetParameter(1,0.005);

    c1.cd();
    fitHist->Fit(histoName,"","",minX[tc],maxX[tc]);

    sprintf(histoName,"ffitup%d",tc);
    if (tc == 1) ffitup[tc] = new TF1(histoName,myfunction2,minX[tc],maxX[tc],3);
    else ffitup[tc] = new TF1(histoName,myfunction,minX[tc],maxX[tc],2);
    
    ffitup[tc]->SetParameter(0,TMath::Max(0.,ffit[tc]->GetParameter(0)+ffit[tc]->GetParError(0)));
    // ffitup[tc]->SetParameter(1,ffit[tc]->GetParameter(1)+ffit[tc]->GetParError(1));
    ffitup[tc]->SetParameter(1,ffit[tc]->GetParameter(1));
    if (tc == 1) ffitup[tc]->SetParameter(2,ffit[tc]->GetParameter(2)); 
    ffitup[tc]->SetParNames("constant","slope");
    ffitup[tc]->SetLineColor(4);
    ffitup[tc]->SetLineStyle(kDashed);
    ffitup[tc]->SetLineWidth(3);  

    ffitup[tc]->Draw("same");

    sprintf(histoName,"ffitdown%d",tc);
    if (tc == 1) ffitdown[tc] = new TF1(histoName,myfunction2,minX[tc],maxX[tc],3);
    else ffitdown[tc] = new TF1(histoName,myfunction,minX[tc],maxX[tc],2);
    
    ffitdown[tc]->SetParameter(0,TMath::Max(0.,ffit[tc]->GetParameter(0)-ffit[tc]->GetParError(0)));
    // ffitdown[tc]->SetParameter(1,ffit[tc]->GetParameter(1)-ffit[tc]->GetParError(1));
    ffitdown[tc]->SetParameter(1,ffit[tc]->GetParameter(1));
    if (tc == 1) ffitdown[tc]->SetParameter(2,ffit[tc]->GetParameter(2)); 
    ffitdown[tc]->SetParNames("constant","slope");
    ffitdown[tc]->SetLineColor(4);
    ffitdown[tc]->SetLineStyle(kDashed);
    ffitdown[tc]->SetLineWidth(3);
    
    ffitdown[tc]->Draw("same");
    
    sprintf(histoName,"~/www/graviton/%s/bkgshapeAlphaMethod_%s%s.png",dirout.c_str(),channelSPart1[tc].c_str(),channelSPart2[tc].c_str()); 
    c1.SaveAs(histoName);
  }
}
