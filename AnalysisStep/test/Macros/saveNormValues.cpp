#include "../Plotter/root_lib/XSecReader.h"
#include <TString.h>
#include <TChain.h>

#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TH1F.h>

using namespace std;

int main (int argc, char ** argv) 
{
  bool is8TeV = false;

  XSecReader *xsecRead;
  if(is8TeV) xsecRead = new XSecReader("../Plotter/Xsection8TeV_YR3.txt","../Plotter/Luminosity.txt");
  else xsecRead = new XSecReader("../Plotter/Xsection7TeV_YR3.txt","../Plotter/Luminosity.txt");

  ofstream outputFile;
  outputFile.open("saveNormValues.txt");
  if(!outputFile){
    std::cout << "Error: no output file!" << std::endl;
    return 0;
  }

  //We start with ggZZ
  //Proper weights:
  //   4mu   : w_4mu_ggzz4l
  //   4e    : w_4e_ggzz4l
  //   2e2mu : w_2e2mu_ggzz2l2l
  //Improper weights (must be normalized to proper weights):
  //   4mu   : w_4mu_ggzz2l2l
  //   4e    : w_4e_ggzz2l2l
  //   2e2mu : w_2e2mu_ggzz4l

  string histoName4e = "ZZ4eTree/Counters";
  string histoName4mu = "ZZ4muTree/Counters";
  string histoName2e2mu = "ZZ2e2muTree/Counters";

  //************************************************
  //******************   ggZZ   ********************
  //************************************************

  TFile fIn_ggZZ4l("rootuples/130702/PRODFSR/ZZ4lAnalysis_ggZZ4l.root");
  TH1F *histo_ggZZ4l = (TH1F*)fIn_ggZZ4l.Get(histoName4e.c_str());
  Double_t Nevt_Gen_ggZZ4l  = histo_ggZZ4l->GetBinContent(1);
  Double_t gen_ZZ4mu_ggZZ4l = histo_ggZZ4l->GetBinContent(2); //gen_ZZ4mu
  Double_t gen_ZZ4e_ggZZ4l  = histo_ggZZ4l->GetBinContent(3); //gen_ZZ4e

  if(is8TeV){
    gen_ZZ4mu_ggZZ4l = Nevt_Gen_ggZZ4l/3.;
    gen_ZZ4e_ggZZ4l = Nevt_Gen_ggZZ4l/3.;
  }

  TFile fIn_ggZZ2l2l("rootuples/130702/PRODFSR/ZZ4lAnalysis_ggZZ2l2l.root");
  TH1F *histo_ggZZ2l2l = (TH1F*)fIn_ggZZ2l2l.Get(histoName2e2mu.c_str());
  Double_t Nevt_Gen_ggZZ2l2l  = histo_ggZZ2l2l->GetBinContent(1);
  Double_t gen_ZZ2e2mu_ggZZ2l2l = histo_ggZZ2l2l->GetBinContent(4); //gen_ZZ2e2mu

  if(is8TeV) gen_ZZ2e2mu_ggZZ2l2l = Nevt_Gen_ggZZ2l2l/3.;
 
 //proper weights
  float w_4mu_ggzz4l = xsecRead->getWeight("ggZZ4l", "1fb-1","all", true)/Nevt_Gen_ggZZ4l;
  float w_4e_ggzz4l = xsecRead->getWeight("ggZZ4l", "1fb-1","all", true)/Nevt_Gen_ggZZ4l;
  float w_2e2mu_ggzz2l2l = xsecRead->getWeight("ggZZ2l2l", "1fb-1","all", true)/Nevt_Gen_ggZZ2l2l;

  float w1_4mu_ggZZ4l = 1/gen_ZZ4mu_ggZZ4l;
  float w1_4e_ggZZ4l = 1/gen_ZZ4e_ggZZ4l;
  float w1_2e2mu_ggZZ2l2l = 1/gen_ZZ2e2mu_ggZZ2l2l;

  outputFile << "4mu_ggZZ4l "   << w1_4mu_ggZZ4l << endl;
  outputFile << "4e_ggZZ4l "    << w1_4e_ggZZ4l << endl;
  outputFile << "2e2mu_ggZZ2l2l " << w1_2e2mu_ggZZ2l2l << endl;

  //improper weights
  float w_4mu_ggzz2l2l = xsecRead->getWeight("ggZZ2l2l", "1fb-1","all", true)/Nevt_Gen_ggZZ2l2l;
  float w_4e_ggzz2l2l = xsecRead->getWeight("ggZZ2l2l", "1fb-1","all", true)/Nevt_Gen_ggZZ2l2l;
  float w_2e2mu_ggzz4l = xsecRead->getWeight("ggZZ4l", "1fb-1","all", true)/Nevt_Gen_ggZZ4l;

  float w1_4mu_ggZZ2l2l = w_4mu_ggzz2l2l/(w_4mu_ggzz4l * gen_ZZ4mu_ggZZ4l);
  float w1_4e_ggZZ2l2l = w_4e_ggzz2l2l/(w_4e_ggzz4l * gen_ZZ4e_ggZZ4l);
  float w1_2e2mu_ggZZ4l = w_2e2mu_ggzz4l/(w_2e2mu_ggzz2l2l * gen_ZZ2e2mu_ggZZ2l2l);

  std::cout << "w1_4mu_ggZZ2l2l = " << w1_4mu_ggZZ2l2l << std::endl; 

  outputFile << "4mu_ggZZ2l2l "   << w1_4mu_ggZZ2l2l << endl;
  outputFile << "4e_ggZZ2l2l "    << w1_4e_ggZZ2l2l << endl;
  outputFile << "2e2mu_ggZZ4l "   << w1_2e2mu_ggZZ4l << endl;


  //************************************************
  //******************   qqZZ   ********************
  //************************************************

  //ZZ->4e

  TFile fIn_qqZZ4e("rootuples/130702/PRODFSR/ZZ4lAnalysis_ZZTo4e.root");
  TH1F *histo_qqZZ4e = (TH1F*)fIn_qqZZ4e.Get(histoName4e.c_str());
  Double_t Nevt_Gen_qqZZ4e  = histo_qqZZ4e->GetBinContent(1);
  Double_t gen_ZZ4e_qqZZ4e  = histo_qqZZ4e->GetBinContent(3); //gen_ZZ4e

  //ZZ->4mu

  TFile fIn_qqZZ4mu("rootuples/130702/PRODFSR/ZZ4lAnalysis_ZZTo4mu.root");
  TH1F *histo_qqZZ4mu = (TH1F*)fIn_qqZZ4mu.Get(histoName4mu.c_str());
  Double_t Nevt_Gen_qqZZ4mu  = histo_qqZZ4mu->GetBinContent(1);
  Double_t gen_ZZ4mu_qqZZ4mu = histo_qqZZ4mu->GetBinContent(2); //gen_ZZ4mu

  //ZZ->2e2mu

  TFile fIn_qqZZ2e2mu("rootuples/130702/PRODFSR/ZZ4lAnalysis_ZZTo2e2mu.root");
  TH1F *histo_qqZZ2e2mu = (TH1F*)fIn_qqZZ2e2mu.Get(histoName2e2mu.c_str());
  Double_t Nevt_Gen_qqZZ2e2mu  = histo_qqZZ2e2mu->GetBinContent(1);
  Double_t gen_ZZ2e2mu_qqZZ2e2mu = histo_qqZZ2e2mu->GetBinContent(4); //gen_ZZ2e2mu

  //ZZ->2mu2tau

  TFile fIn_qqZZ2mu2tau("rootuples/130702/PRODFSR/ZZ4lAnalysis_ZZTo2mu2tau.root");
  TH1F *histo_qqZZ2mu2tau = (TH1F*)fIn_qqZZ2mu2tau.Get(histoName4mu.c_str());
  Double_t Nevt_Gen_qqZZ2mu2tau    = histo_qqZZ2mu2tau->GetBinContent(1);
  Double_t gen_ZZ4mu_qqZZ2mu2tau   = histo_qqZZ2mu2tau->GetBinContent(2); //gen_ZZ4mu
  Double_t gen_ZZ4e_qqZZ2mu2tau    = histo_qqZZ2mu2tau->GetBinContent(3); //gen_ZZ4e
  Double_t gen_ZZ2e2mu_qqZZ2mu2tau = histo_qqZZ2mu2tau->GetBinContent(4); //gen_ZZ2e2mu

  //ZZ->2e2tau

  TFile fIn_qqZZ2e2tau("rootuples/130702/PRODFSR/ZZ4lAnalysis_ZZTo2e2tau.root");
  TH1F *histo_qqZZ2e2tau = (TH1F*)fIn_qqZZ2e2tau.Get(histoName4e.c_str());
  Double_t Nevt_Gen_qqZZ2e2tau    = histo_qqZZ2e2tau->GetBinContent(1);
  Double_t gen_ZZ4mu_qqZZ2e2tau   = histo_qqZZ2e2tau->GetBinContent(2); //gen_ZZ4mu
  Double_t gen_ZZ4e_qqZZ2e2tau    = histo_qqZZ2e2tau->GetBinContent(3); //gen_ZZ4e
  Double_t gen_ZZ2e2e_qqZZ2e2tau = histo_qqZZ2e2tau->GetBinContent(4); //gen_ZZ2e2mu

  float w_4mu_qqzz4mu    = xsecRead->getWeight("ZZTo4mu",  "1fb-1","all",  true)/Nevt_Gen_qqZZ4mu;
  float w_4e_qqzz4e      = xsecRead->getWeight("ZZTo4e",   "1fb-1","all",   true)/Nevt_Gen_qqZZ4e;
  float w_2e2mu_qqzz2e2mu     = xsecRead->getWeight("ZZTo2e2mu",  "1fb-1","all",  true)/Nevt_Gen_qqZZ2e2mu;

  float w_2mu2tau_qqzz2mu2tau = xsecRead->getWeight("ZZTo2mu2tau","1fb-1","all",true)/Nevt_Gen_qqZZ2mu2tau;
  float w_2e2tau_qqzz2e2tau   = xsecRead->getWeight("ZZTo2e2tau", "1fb-1","all", true)/Nevt_Gen_qqZZ2e2tau;

  //proper weights
  float w1_4mu_qqZZ4mu     = 1/gen_ZZ4mu_qqZZ4mu;
  float w1_4e_qqZZ4e       = 1/gen_ZZ4e_qqZZ4e;
  float w1_2e2mu_qqZZ2e2mu = 1/gen_ZZ2e2mu_qqZZ2e2mu;

  outputFile << "4mu_ZZTo4mu "     << w1_4mu_qqZZ4mu << endl;
  outputFile << "4e_ZZTo4e "       << w1_4e_qqZZ4e << endl;
  outputFile << "2e2mu_ZZTo2e2mu " << w1_2e2mu_qqZZ2e2mu << endl;

  //improper weights - w1_CHANNEL_SAMPLE

  //4mu

  float w1_4mu_qqZZ4e      = w_4e_qqzz4e          /(w_4mu_qqzz4mu * gen_ZZ4mu_qqZZ4mu);
  float w1_4mu_qqZZ2e2mu   = w_2e2mu_qqzz2e2mu    /(w_4mu_qqzz4mu * gen_ZZ4mu_qqZZ4mu);
  float w1_4mu_qqZZ2mu2tau = w_2mu2tau_qqzz2mu2tau/(w_4mu_qqzz4mu * gen_ZZ4mu_qqZZ4mu);
  float w1_4mu_qqZZ2e2tau  = w_2e2tau_qqzz2e2tau  /(w_4mu_qqzz4mu * gen_ZZ4mu_qqZZ4mu);

  outputFile << "4mu_ZZTo4e "      << w1_4mu_qqZZ4e << endl;
  outputFile << "4mu_ZZTo2e2mu "   << w1_4mu_qqZZ2e2mu << endl;
  outputFile << "4mu_ZZTo2mu2tau " << w1_4mu_qqZZ2mu2tau << endl;
  outputFile << "4mu_ZZTo2e2tau "  << w1_4mu_qqZZ2e2tau << endl;

  //4e

  float w1_4e_qqZZ4mu     = w_4mu_qqzz4mu        /(w_4e_qqzz4e * gen_ZZ4e_qqZZ4e);
  float w1_4e_qqZZ2e2mu   = w_2e2mu_qqzz2e2mu    /(w_4e_qqzz4e * gen_ZZ4e_qqZZ4e);
  float w1_4e_qqZZ2mu2tau = w_2mu2tau_qqzz2mu2tau/(w_4e_qqzz4e * gen_ZZ4e_qqZZ4e);
  float w1_4e_qqZZ2e2tau  = w_2e2tau_qqzz2e2tau  /(w_4e_qqzz4e * gen_ZZ4e_qqZZ4e);

  outputFile << "4e_ZZTo4mu "     << w1_4e_qqZZ4mu << endl;
  outputFile << "4e_ZZTo2e2mu "   << w1_4e_qqZZ2e2mu << endl;
  outputFile << "4e_ZZTo2mu2tau " << w1_4e_qqZZ2mu2tau << endl;
  outputFile << "4e_ZZTo2e2tau "  << w1_4e_qqZZ2e2tau << endl;

  //2e2mu

  float w1_2e2mu_qqZZ4mu     = w_4mu_qqzz4mu        /(w_2e2mu_qqzz2e2mu * gen_ZZ2e2mu_qqZZ2e2mu);
  float w1_2e2mu_qqZZ4e      = w_4e_qqzz4e          /(w_2e2mu_qqzz2e2mu * gen_ZZ2e2mu_qqZZ2e2mu);
  float w1_2e2mu_qqZZ2mu2tau = w_2mu2tau_qqzz2mu2tau/(w_2e2mu_qqzz2e2mu * gen_ZZ2e2mu_qqZZ2e2mu);
  float w1_2e2mu_qqZZ2e2tau  = w_2e2tau_qqzz2e2tau  /(w_2e2mu_qqzz2e2mu * gen_ZZ2e2mu_qqZZ2e2mu);

  outputFile << "2e2mu_ZZTo4mu "     << w1_2e2mu_qqZZ4mu << endl;
  outputFile << "2e2mu_ZZTo4e "      << w1_2e2mu_qqZZ4e << endl;
  outputFile << "2e2mu_ZZTo2mu2tau " << w1_2e2mu_qqZZ2mu2tau << endl;
  outputFile << "2e2mu_ZZTo2e2tau "  << w1_2e2mu_qqZZ2e2tau << endl;

  outputFile.close();

  return 0;
}
