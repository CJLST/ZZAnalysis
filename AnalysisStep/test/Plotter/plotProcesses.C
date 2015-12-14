/* 
 * usage: 
 * -specify 'inputPath' at the end of this file
 * -run with:
 *   root -l plotProcesses.C++
 */

#include <iostream>
#include <fstream>
#include <utility>
#include <map>
#include <vector>
#include <cmath>

#include "TCanvas.h"
#include "TColor.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "plotUtils.C"

#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >  LV;

using namespace std;

#define MASSZ 91.1876

#define DEBUG 0
#define FINALSTATE 3 // 0:4mu, 1:4e, 2:2e2mu, 3:4mu+4e+2e2mu, 4: none('emptyEvents')

#define PRINTPLOTINFO 1
string INFO = "13TeV samples";
string INFO2 = "#sqrt{s} = 13 TeV";

#define VARIABLELIST 0 // put 0 for plots shown on September 8th

#define BASKETLIST   5 // choice of category definitions:
                       //  - categories shown on Sep 23rd 2014 : nb. 1
                       //  - categories shown on Nov 7th 2014  : nb. 0 and 3 (subcat: 2)
                       //  - categories shown on Dec 19th 2014 : nb. 3 (subcat: 2) and 5 (subcat: 4)

#define requireExactly4GoodLeps 0 
#define requireAtLeast5GoodLeps 0 
#define requireExactly5GoodLeps 0 
#define requireExactly6GoodLeps 0 
#define requireHLepsAreInEtaPtAcc    0 // impose that all 4 gen leptons from the Higgs are in the acceptance
#define requireHLepsAreGood     0 // impose that all 4 gen leptons from the Higgs are matched to good leptons (from the candidate or extra leptons)

#define excludeH2l2X            1 // completely exclude ZH,H->2l2X and ttH,H->2l2X events from the study
#define treatH2l2XAsBkgd        1 // treat these events as a background in categorization-related plots 
#define printYields             0 // print yield numbers on basket efficiency plots

#define doProdComp       0 // plots shown on September 8th
#define doProdCompMatch4 0
#define doMatch4OrNot    0
#define doMatchHLeps     0
#define doMatchAllLeps   0
#define doMatchWHZHttH   0 // plots shown on September 8th
#define do2DPlots        0 // plots shown on September 8th
#define doBaskets        1 // categorization plots
#define doYieldStudy     0
#define doTTHCheck       0
#define doPlotsByGenChan 0

#define nSamples 8
#define nHiggsSamples 7
#define nGenChannels 11
#define nRecoChannels 7
#define nVariables 16
#define nBaskets (BASKETLIST==0?7:BASKETLIST==1?9:BASKETLIST==2?24:BASKETLIST==3?6:BASKETLIST==4?38:BASKETLIST==5?7:BASKETLIST==6?10:1)

#define nMatchHLepsStatuses 6
#define nMatchAllLepsStatuses 5
#define nMatchWHStatuses 5
#define nMatchZHStatuses 7
#define nMatchttHStatuses 9

#define groupWminusWplus 1
#define nAssocWDecays 2
#define nAssocZDecays 3
#define nAssocttDecays 4
#define nAssocDecays (nAssocWDecays + nAssocZDecays + nAssocttDecays)





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////////////////////// MISCELLANEOUS FUNCTIONS ////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


bool test_bit_16( short mask, unsigned int iBit ) { return (mask >> iBit) & 1; }

Double_t deltaR(Double_t e1, Double_t p1, Double_t e2, Double_t p2) {
  Double_t deltaPhi = acos(cos(p1-p2));
  return TMath::Sqrt((e1-e2)*(e1-e2) + deltaPhi*deltaPhi);
}

string eventID(Int_t nRun, Int_t nEvent, Int_t nLumi) {
  return string(Form("run%i_lumi%i_event%i",nRun,nEvent,nLumi));
}

void printStatus(Int_t n, Int_t interval, Int_t total, string label) {
  if(n==0)              cout<<"  "<<1  <<" / "<<total<<" "<<label<<endl;
  if((n+1)%interval==0) cout<<"  "<<n+1<<" / "<<total<<" "<<label<<endl;
  if(n==total-1)        cout<<"  "<<n+1<<" / "<<total<<" "<<label<<endl;
}

string fixWidth(string str, unsigned n, bool atBeginning) {
  if(str.length()>n){
    cout<<"Error in function fixWidth('"<<str<<"',"<<n<<") : string '"<<str<<"' is too long."<<endl;
    return ""; 
  }else{
    string res = "";
    if(atBeginning) res += str;
    for(int i=0; i<(int)n-(int)(str.length()); i++) res += " ";
    if(!atBeginning) res += str;
    return res;
  }
}

string repeat(string pattern, int n) {
  string res = "";
  for(int i=0; i<n; i++) res += pattern;
  return res;
}

string rounding2(float x) {
  return string(Form("%.2f",x));
}
string rounding3(float x) {
  return string(Form("%.3f",x));
}

string percentage(float frac) {
  return string(Form("%.1f",100*frac));
}

void prepareBasketlabels(TH1F* h, vector<string> basketLabel){
  for(int i=1; i<=h->GetNbinsX(); i++)  h->GetXaxis()->SetBinLabel(i,basketLabel[i-1].c_str());
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetXaxis()->LabelsOption("h");
  h->GetXaxis()->SetNdivisions(-h->GetNbinsX());
}

void prepareBasketlabelsHorizontal(TH2F* h2, vector<string> basketLabel){
  int nbins = h2->GetNbinsY();
  for(int i=1; i<=nbins; i++)  h2->GetYaxis()->SetBinLabel(nbins+1-i,basketLabel[i-1].c_str());
  h2->GetYaxis()->SetLabelSize(0.05);
  h2->GetYaxis()->LabelsOption("h");
  h2->GetYaxis()->SetNdivisions(-h2->GetNbinsY());
}

void TriBulle(float* array, int length, bool ascendingOrder = true){
  int i = length;
  bool swap = true;
  while(i>0 && swap){
    swap = false;
    for(int j=0; j<i-1; j++){
      if( (array[j]>array[j+1] && ascendingOrder) || (array[j]<array[j+1] && !ascendingOrder) ){
	float tmp = array[j];
	array[j] = array[j+1];
	array[j+1] = tmp;
	swap = true;
      }
    }
    i--;
  }
}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////// PLOTTING FUNCTIONS ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void DrawByProdmodes(TCanvas* c, TH1F* hSource[nSamples][nRecoChannels], string* prodName, Bool_t* isPresent, Color_t* colors, Bool_t logY = false) {

  gStyle->SetOptTitle(0);

  c->cd();
  if(logY) c->SetLogy();
  c->SetTicks(0,0);

  TLegend* lgd = new TLegend(0.69,0.67,0.89,0.89);
  lgd->SetFillStyle(0);
  lgd->SetBorderSize(0);

  TH1F* h[nSamples];
  for(int p=0; p<nSamples; p++){
    if(!isPresent[p]) continue;
    h[p] = (TH1F*)hSource[p][FINALSTATE]->Clone(); // includes the 3 decay channels !!
  }

  Float_t max = 0.;
  for(int p=0; p<nSamples; p++){
    if(!isPresent[p]) continue;

    h[p]->Scale(1/h[p]->Integral());
    Float_t maxtemp = h[p]->GetMaximum();
    if(maxtemp>max) max = maxtemp;

  }

  for(int p=0; p<nSamples; p++){
    if(!isPresent[p]) continue;
    string pn = prodName[p];

    h[p]->SetLineColor(colors[p]);
    h[p]->SetLineWidth(2);
    h[p]->SetFillStyle(0);
    h[p]->SetStats(0);
    //if(!logY) h[p]->SetMinimum(0);

    if(p==0){
      h[p]->SetMaximum(1.1*max);
      h[p]->GetYaxis()->SetTitle("normalized to 1");
      h[p]->Draw("hist"); 
    }else{
      h[p]->Draw("histsame");
    }

    lgd->AddEntry( h[p], pn.c_str(), "l" );
  }

  lgd->Draw();

  gPad->RedrawAxis();

}

void DrawMatch4OrNot(TCanvas* c, TH1F* h, TH1F* hMatch, string title, Color_t color, Color_t colorNoMatch, Bool_t logY = false) {

  gStyle->SetOptTitle(1);

  c->cd();
  if(logY) c->SetLogy();
  c->SetTicks(0,0);

//   Int_t scale = 1/h->Integral();
//   h->Scale(scale);
//   hMatch->Scale(scale);

  h->SetLineColor(colorNoMatch);
  h->SetFillColor(colorNoMatch);
  hMatch->SetLineColor(color);
  hMatch->SetFillColor(color);
  h->SetTitle(title.c_str()); 
  hMatch->SetTitle(title.c_str());
  h->SetStats(0);
  hMatch->SetStats(0);
  //if(!logY) h->SetMinimum(0);
  
  h->Draw(); 
  hMatch->Draw("sames");
  
  TLegend* lgd = new TLegend(0.78,0.70,0.98,0.82);
  lgd->SetFillColor(0);
  lgd->AddEntry(hMatch,"match 4 gen","f");
  lgd->AddEntry(h,"other cases","f");
  lgd->Draw();

  gPad->RedrawAxis();

}

void DrawMatchHLeps(TCanvas* c, TH1F* h, TH1F** hMatch, string title, Color_t* color, string* matchHLepsKeys, Bool_t logY = false) {

  gStyle->SetOptTitle(1);

  c->cd();
  if(logY) c->SetLogy();
  c->SetTicks(0,0);

  h->SetLineColor(1);
  h->SetFillColor(1);
  h->SetTitle(title.c_str()); 
  h->SetStats(0);
  //if(!logY) h->SetMinimum(0);
  h->Draw(); 

  TH1F* hStacks[nMatchHLepsStatuses];
  for(int m=0; m<nMatchHLepsStatuses; m++){
    hStacks[m] = (TH1F*)hMatch[m]->Clone();
    if(m>0) hStacks[m]->Add(hStacks[m-1]);
    hStacks[m]->SetLineColor(color[m]);
    hStacks[m]->SetFillColor(color[m]);
    hStacks[m]->SetTitle(title.c_str()); 
    hStacks[m]->SetStats(0);
  }
  for(int m=nMatchHLepsStatuses-1; m>=0; m--){
    hStacks[m]->Draw("same");
  }
  
  TLegend* lgd = new TLegend(0.78,0.45,0.98,0.82);
  lgd->SetFillColor(0);
  for(int m=0; m<nMatchHLepsStatuses; m++){
    lgd->AddEntry(hStacks[m],matchHLepsKeys[m].c_str(),"f");
  }
  lgd->Draw();

  gPad->RedrawAxis();

}

void DrawMatchAllLeps(TCanvas* c, TH1F* h, TH1F** hMatch, string title, Color_t* color, string* matchAllLepsKeys, Bool_t logY = false) {

  gStyle->SetOptTitle(1);

  c->cd();
  if(logY) c->SetLogy();
  c->SetTicks(0,0);

  h->SetLineColor(1);
  h->SetFillColor(1);
  h->SetTitle(title.c_str()); 
  h->SetStats(0);
  //if(!logY) h->SetMinimum(0);
  h->Draw(); 

  TH1F* hStacks[nMatchAllLepsStatuses];
  for(int m=0; m<nMatchAllLepsStatuses; m++){
    hStacks[m] = (TH1F*)hMatch[m]->Clone();
    if(m>0) hStacks[m]->Add(hStacks[m-1]);
    hStacks[m]->SetLineColor(color[m]);
    hStacks[m]->SetFillColor(color[m]);
    hStacks[m]->SetTitle(title.c_str()); 
    hStacks[m]->SetStats(0);
  }
  for(int m=nMatchAllLepsStatuses-1; m>=0; m--){
    hStacks[m]->Draw("same");
  }
  
  TLegend* lgd = new TLegend(0.73,0.50,0.98,0.82);
  lgd->SetFillStyle(0);
  for(int m=0; m<nMatchAllLepsStatuses; m++){
    lgd->AddEntry(hStacks[m],matchAllLepsKeys[m].c_str(),"f");
  }
  lgd->Draw();

  gPad->RedrawAxis();

}

void DrawMatchCustom(TCanvas* c, TH1F* h, TH1F** hMatch, string title, Color_t* color, string* matchKeys, Int_t nStatuses, Float_t legLeft, Float_t legBottom, Bool_t logY = false) {

  gStyle->SetOptTitle(1);

  c->cd();
  if(logY) c->SetLogy();
  c->SetTicks(0,0);

  h->SetLineColor(1);
  h->SetFillColor(1);
  h->SetTitle(title.c_str()); 
  h->SetStats(0);
  //if(!logY) h->SetMinimum(0);
  h->Draw(); 

  Float_t denom = h->Integral();
  string percentages[20];
  for(int m=0; m<nStatuses; m++){
    if(excludeH2l2X && matchKeys[m].find("2l2X")!=string::npos) continue;
    percentages[m] = percentage(hMatch[m]->Integral()/denom);
  }

  TH1F* hStacks[20];
  for(int m=0; m<nStatuses; m++){
    hStacks[m] = (TH1F*)hMatch[m]->Clone();
    if(m>0) hStacks[m]->Add(hStacks[m-1]);
    hStacks[m]->SetLineColor(color[m]);
    hStacks[m]->SetFillColor(color[m]);
    hStacks[m]->SetTitle(title.c_str()); 
    hStacks[m]->SetStats(0);
  }
  for(int m=nStatuses-1; m>=0; m--){
    hStacks[m]->Draw("same");
  }
  
  TLegend* lgd = new TLegend(legLeft,legBottom,legLeft+0.55,0.89);
  lgd->SetFillStyle(0);
  lgd->SetBorderSize(0);
  for(int m=0; m<nStatuses; m++){
    if(excludeH2l2X && matchKeys[m].find("2l2X")!=string::npos) continue;
    lgd->AddEntry(hStacks[m],matchKeys[m].c_str(),"f");
  }
  lgd->Draw();
  
  TPaveText* pav = new TPaveText(legLeft+0.03,legBottom,legLeft+0.11,0.89,"brNDC");
  pav->SetFillStyle(0);
  pav->SetBorderSize(0);
  pav->SetTextColor(0);
  for(int m=0; m<nStatuses; m++){
    if(excludeH2l2X && matchKeys[m].find("2l2X")!=string::npos) continue;
    pav->AddText((percentages[m]+" %").c_str());
  }
  pav->Draw();

  gPad->RedrawAxis();

}

void Draw2D(TCanvas* c, TH2F* h, string title) {
  gStyle->SetOptTitle(1);
  setColZGradient_1();
  c->cd();
  h->SetTitle(title.c_str());
  h->SetStats(0);
  h->Draw("COLZ"); 
  if(PRINTPLOTINFO) printInfo(INFO2,0.1,0.9,0.75,0.94);
}

void DrawHorizontal(TH1 *h, vector<string> basketLabel, Bool_t same = false, Int_t nDiv = 0, Bool_t without1stBin = false) {
  Double_t ymin = 0.; //h->GetMinimum();
  Double_t ymax = h->GetMaximum(); // 1.1*h->GetMaximum(); 
  TAxis *axis   = h->GetXaxis();
  Double_t xmin = axis->GetXmin();
  Double_t xmax = axis->GetXmax();
  Int_t nbins   = axis->GetNbins();
  TH2F *h2 = new TH2F(Form("h2_%i",rand()),Form("%s;%s;%s",h->GetTitle(),h->GetYaxis()->GetTitle(),""),10,ymin,ymax,nbins,xmin,xmax);
  h2->SetBit(kCanDelete);
  //h2->SetDirectory(0);
  h2->SetTickLength(0.01);
  h2->GetXaxis()->SetLabelSize(0.03);
  if(without1stBin) h2->GetYaxis()->SetRangeUser(0,nBaskets-2);
  h2->SetStats(0);
  prepareBasketlabelsHorizontal(h2,basketLabel);
  h2->Draw(same?"same":"");
  if(nDiv!=0) h2->GetXaxis()->SetNdivisions(nDiv);
  TBox box;
  Int_t color = h->GetFillColor();
  if (color == 0) color = 1;
  box.SetFillColor(color);
  Double_t dy,x1,y1,x2,y2;
  for (Int_t i=1;i<=nbins;i++) {
    if(without1stBin && i==1) continue;
    dy = axis->GetBinWidth(nbins+1-i);
    x1 = 0.;
    y1 = axis->GetBinCenter(nbins+1-i)-0.4*dy;
    x2 = h->GetBinContent(i);
    y2 = axis->GetBinCenter(nbins+1-i)+0.4*dy;
    box.DrawBox(x1,y1,x2,y2);
  }
}

void DrawBasketEfficiencies(TCanvas* c, TH1F* h1, string pn, TH1F** hAssocDecay, string* assocDecayName, Bool_t H2l2XAsBkgd, vector<string> basketLabel, Bool_t logY = false) {

  if(pn=="WminusH") { cout<<"ERROR in method DrawBasketEfficiencies"<<endl; return; }

  bool without1stBin = false;

  gStyle->SetOptTitle(1);
  gStyle->SetTitleX(0.25);

  c->cd();
  if(logY) c->SetLogy();
  //c->SetTicks(0,0);
  c->SetLeftMargin(0.25);

  bool doAssoc = (pn=="WH"||pn=="WplusH"||pn=="WminusH"||pn=="ZH"||pn=="ttH");
  Int_t nDecays = (pn=="WH"||pn=="WplusH"||pn=="WminusH") ? nAssocWDecays : pn=="ZH" ? nAssocZDecays : pn=="ttH" ? nAssocttDecays : 0;
  Int_t startAt = (pn=="WH"||pn=="WplusH"||pn=="WminusH") ? 0 : pn=="ZH" ? nAssocWDecays : pn=="ttH" ? nAssocWDecays+nAssocZDecays : 0;

  TH1F* h = (TH1F*)h1->Clone();
  if(H2l2XAsBkgd){
    if(pn=="ZH") h->Add(hAssocDecay[nAssocWDecays+nAssocZDecays-1],-1);
    if(pn=="ttH") h->Add(hAssocDecay[nAssocWDecays+nAssocZDecays+nAssocttDecays-1],-1);
  }
  if(doAssoc) h->SetLineColor(1);
  if(doAssoc) h->SetFillColor(1);
  h->SetTitle(pn.c_str()); 
  h->SetStats(0);
  //if(!logY) h->SetMinimum(0);
  //h->GetXaxis()->SetRangeUser(1,nBaskets-1);
  //h->Draw(); 
  DrawHorizontal(h,basketLabel,false,508,without1stBin);

  TH1F* hStacks[nDecays];
  for(int a=0; a<nDecays; a++){
    if(H2l2XAsBkgd && (pn=="ZH"||pn=="ttH") && a==nDecays-1) continue;
    hStacks[a] = (TH1F*)hAssocDecay[startAt+a]->Clone();
    if(a>0) hStacks[a]->Add(hStacks[a-1]);
    hStacks[a]->SetStats(0);
  }
  bool first = true;
  for(int a=nDecays-1; a>=0; a--){
    if(H2l2XAsBkgd && (pn=="ZH"||pn=="ttH") && a==nDecays-1) continue;
    //hStacks[a]->Draw("same");
    DrawHorizontal(hStacks[a],basketLabel,!first,0,without1stBin);
    first = false;
  }
  
  if(doAssoc){
    TLegend* lgd = new TLegend((H2l2XAsBkgd?0.7:0.55),0.13,0.86,(((pn=="WH"||pn=="WplusH"||pn=="WminusH")||(pn=="ZH"&&H2l2XAsBkgd))? 0.24 : (pn=="ZH"||(pn=="ttH"&&H2l2XAsBkgd))? 0.3 : pn=="ttH"? 0.36 : 0.));
    //lgd->SetFillColor(0);
    lgd->SetFillStyle(0);
    lgd->SetBorderSize(0);
    for(int a=0; a<nDecays; a++){
      if(H2l2XAsBkgd && (pn=="ZH"||pn=="ttH") && a==nDecays-1) continue;
      lgd->AddEntry(hStacks[a],assocDecayName[startAt+a].c_str(),"f");
    }
    lgd->Draw();
  }

  if(printYields){
    TPaveText* pavNExp = new TPaveText(0.905,0.1,0.995,0.9,"brNDC");
    pavNExp->SetFillStyle(0);
    pavNExp->SetBorderSize(0);
    for(int b=2; b<=h->GetNbinsX(); b++){
      pavNExp->AddText(Form("%.4f",h->GetBinContent(b)));
    }
    pavNExp->Draw();
  }

  if(PRINTPLOTINFO) printInfo(INFO,0.23,0.9,0.75,0.94);

  gPad->RedrawAxis();

}

vector<TH1F*> StackForPurity(vector<TH1F*> histos) {

  vector<TH1F*> result;

  int n = histos.size();
  for(int i=0; i<n; i++){
    TH1F* hTemp = (TH1F*)histos[i]->Clone();
    if(i>0) hTemp->Add(result[i-1]);
    result.push_back(hTemp);
  }

  int nBins = histos[0]->GetNbinsX();
  for(int b=1; b<=nBins; b++){
    Float_t norm = result[n-1]->GetBinContent(b);
    for(int i=0; i<n; i++){
      result[i]->SetBinContent(b,result[i]->GetBinContent(b)/norm);
    }
  }

  return result;
}

void DrawBasketPurities(TCanvas* c, TH1F** h, string* prodName, Bool_t* isPresent, TH1F** hAssocDecay, string* assocDecayName, Bool_t withAssocDecays, Bool_t H2l2XAsBkgd, vector<string> basketLabel, float lumi) {

  float offset = 0.1;

  gStyle->SetOptTitle(0);
  //gStyle->SetFrameLineWidth(2);

  c->cd();
  c->SetLeftMargin(0.15+offset);
  c->SetRightMargin(0.35-offset);

  vector<TH1F*> histos;
  
  TLegend* lgd = new TLegend(0.67+offset,withAssocDecays?(H2l2XAsBkgd?0.45:0.35):0.65,((0.95+offset<1.)?0.95+offset:1.),0.9);
  lgd->SetFillColor(0);
  lgd->SetBorderSize(0);
  
  TPaveText* pavNExp = new TPaveText(0.15+offset,0.1,0.4+offset,0.9,"brNDC");
  pavNExp->SetFillStyle(0);
  pavNExp->SetBorderSize(0);
  pavNExp->SetTextColor(0);
  TH1F* hNExpectedEvts = (TH1F*)h[0]->Clone();

  for(int p=0; p<nHiggsSamples; p++){
    if(!isPresent[p]) continue;
    string pn = prodName[p];
    if(groupWminusWplus && pn=="WminusH") continue;
    if(groupWminusWplus && pn=="WplusH") pn = "WH";
    bool doAssoc = ((pn=="WH"||pn=="WplusH"||pn=="WminusH")||pn=="ZH"||pn=="ttH");
    Int_t nDecays = (pn=="WH"||pn=="WplusH"||pn=="WminusH") ? nAssocWDecays : pn=="ZH" ? nAssocZDecays : pn=="ttH" ? nAssocttDecays : 0;
    Int_t startAt = (pn=="WH"||pn=="WplusH"||pn=="WminusH") ? 0 : pn=="ZH" ? nAssocWDecays : pn=="ttH" ? nAssocWDecays+nAssocZDecays : 0;

    TH1F* hTemp = (TH1F*)h[p]->Clone();
    if(H2l2XAsBkgd){
      if(pn=="ZH") hTemp->Add(hAssocDecay[nAssocWDecays+nAssocZDecays-1],-1);
      if(pn=="ttH") hTemp->Add(hAssocDecay[nAssocWDecays+nAssocZDecays+nAssocttDecays-1],-1);
    }

    if(withAssocDecays && doAssoc){
      for(int a=0; a<nDecays; a++){
	bool isH2l2X = ((pn=="ZH"||pn=="ttH") && a==nDecays-1); 
	if(!(H2l2XAsBkgd && isH2l2X)){
	  if(!(pn=="WminusH")){
	    histos.push_back(hAssocDecay[startAt+a]);
	    lgd->AddEntry(hAssocDecay[startAt+a],assocDecayName[startAt+a].c_str(),"f");
	  }
	}
      }
    }else{
      histos.push_back(hTemp);
      lgd->AddEntry(hTemp,pn.c_str(),"f");
    }

    if(p!=0) hNExpectedEvts->Add(hTemp);
  }

  vector<TH1F*> stackedHistos = StackForPurity(histos); 
  int n = stackedHistos.size();
  for(int i=0; i<n; i++){
    stackedHistos[n-1-i]->SetStats(0);
    if(i==0){
      stackedHistos[n-1-i]->GetYaxis()->SetRangeUser(0.,1.);
      stackedHistos[n-1-i]->GetYaxis()->SetTitle("signal fraction");
      //stackedHistos[n-1-i]->Draw(); 
      DrawHorizontal(stackedHistos[n-1-i],basketLabel,false,-210);
    }else{
      //stackedHistos[n-1-i]->Draw("same");
      DrawHorizontal(stackedHistos[n-1-i],basketLabel,true,0);
    }
  }

  lgd->Draw();

  for(int b=1; b<=hNExpectedEvts->GetNbinsX(); b++){
    pavNExp->AddText(Form("%.3f exp. events in %.0f fb^{-1}",hNExpectedEvts->GetBinContent(b),lumi));
  }
  pavNExp->Draw();

  if(PRINTPLOTINFO) printInfo(INFO,0.23,0.9,0.75,0.94);
  
  gPad->RedrawAxis();

}

void DrawBasketSOSPB(TCanvas* c, TH1F* h, vector<string> basketLabel) {

  gStyle->SetOptTitle(0);

  c->cd();
  c->SetLeftMargin(0.25);

  h->SetFillColor(kBlue+2);
  h->GetYaxis()->SetRangeUser(0.,1.);
  h->GetYaxis()->SetTitle("S/(S+B)");
  DrawHorizontal(h,basketLabel,false,-210);

  if(PRINTPLOTINFO) printInfo(INFO,0.23,0.9,0.75,0.94);

}

void DrawSorted(TCanvas* c, TH1F** hGen, TH1F** hReco, string pn, string channel, Float_t maxY = -1.) {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  c->cd();

  Color_t colorGen [4] = { kBlue, kRed+1, kGreen+2, kGray+1 };
  Color_t colorReco[4] = { kBlue, kRed+1, kGreen+2, kGray+1 };

  Float_t SF = 1/hGen[0]->Integral();

  for(int sl=0; sl<4; sl++){

    hGen[sl]->Scale(SF);
    hGen[sl]->SetLineColor(colorReco[sl]);
    hGen[sl]->SetLineWidth(2);
    hGen[sl]->SetLineStyle(2);
    hGen[sl]->SetFillColor(0);
    hGen[sl]->SetFillStyle(0);
    hGen[sl]->GetYaxis()->SetTitle("a.u.");
    //hGen[sl]->GetYaxis()->SetTitleOffset(0.5);
    hGen[sl]->GetYaxis()->SetLabelSize(0.03);
    if(maxY>0.) hGen[sl]->SetMaximum(maxY);

    hReco[sl]->Scale(SF);
    hReco[sl]->SetLineColor(colorGen[sl]);
    hReco[sl]->SetLineWidth(2);
    hReco[sl]->SetLineStyle(2);
    hReco[sl]->SetFillColor(colorGen[sl]);
    hReco[sl]->SetFillStyle(3002);
    hReco[sl]->GetYaxis()->SetTitle("a.u.");
    //hReco[sl]->GetYaxis()->SetTitleOffset(0.5);
    hReco[sl]->GetYaxis()->SetLabelSize(0.03);

  }

  hGen[3]->Draw();
  for(int sl=2;sl>=0;sl--) hGen [sl]->Draw("SAME");
  for(int sl=3;sl>=0;sl--) hReco[sl]->Draw("SAME");

  TLegend *leg = new TLegend(0.45,0.6,0.88,0.88);
  leg->SetLineColor(0);
  leg->SetFillColor(0);  
  TString legtitle; legtitle.Form("H#rightarrow ZZ* #rightarrow %s",channel.c_str());
  //TString massstring; massstring.Form("m_{H} = %d GeV",mass);
  leg->AddEntry((TObject*)0, legtitle.Data(), "");
  //leg->AddEntry((TObject*)0, massstring.Data(), "");
  TLegendEntry *before = leg->AddEntry("Before Analysis Selection","Before Analysis Selection","f");
  TLegendEntry *after  = leg->AddEntry("After Analysis Selection","After Analysis Selection","f");
  after->SetLineStyle(2);
  after->SetLineWidth(2);
  after->SetFillColor(kBlack);
  after->SetLineColor(kBlack);
  after->SetFillStyle(3002);
  before->SetLineStyle(2);
  before->SetLineWidth(2);
  before->SetFillColor(0);
  before->SetLineColor(kBlack);

  leg->Draw("SAME");

  gPad->RedrawAxis();

  if(PRINTPLOTINFO) printInfo(INFO2,0.1,0.9,0.75,0.94);

}

void Draw1D(TCanvas* c, TH1F* h) {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  c->cd();

  //Float_t SF = 1/h->Integral();
  //h->Scale(SF);
  h->SetLineColor(kBlue);
  h->SetLineWidth(2);
  h->SetLineStyle(1);
  h->SetFillColor(0);
  h->SetFillStyle(0);
  //h->GetYaxis()->SetTitle("a.u.");
  //h->GetYaxis()->SetTitleOffset(0.5);
  //h->GetYaxis()->SetLabelSize(0.03);

  h->Draw();

  if(PRINTPLOTINFO) printInfo(INFO2,0.1,0.9,0.75,0.94);

}





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// MAIN MACRO //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void run(	
	 string inDir,	   
	 string outDir,
	 float lumi,
	 
	 // path of the input primary tree (just let "" for trees you don't have)
	 string file_ggH     = "",
	 string file_VBF     = "",
	 string file_WplusH  = "",
	 string file_WminusH = "",
	 string file_ZH      = "",
	 string file_ttH     = "",
	 string file_bbH     = "",
	 string file_ZZ4l    = "",
	 
	 string tagIn = "",
	 string tagOut = ""
	 
	 )
{

  srand(time(0));

  ofstream txtOut;
  TString txtOutName = (TString)outDir.c_str()+"/matchingInfo.txt";
  txtOut.open(txtOutName);


  // ------------------------------------------------------------
  // ---------------------- Definitions -------------------------
  // ------------------------------------------------------------

  string prodName[nSamples] = {
    "ggH",
    "VBF",
    "WplusH",
    "WminusH",
    "ZH",
    "ttH",
    "bbH",
    "ZZTo4l",
  };
  string fileName[nSamples] = {
    file_ggH    ,
    file_VBF    ,
    file_WplusH ,
    file_WminusH,
    file_ZH     ,
    file_ttH    ,
    file_bbH    ,
    file_ZZ4l   , 
  };
  Bool_t isPresent[nSamples]; for(int p=0; p<nSamples; p++) isPresent[p] = (fileName[p]!="");
  //Bool_t allSignalsArePresent = (file_ggH!="" && file_VBF!="" && file_WplusH!="" && file_WminusH!="" && file_ZH!="" && file_ttH!="" && file_bbH!="");
  Bool_t allSignalsArePresent = (file_ggH!="" && file_VBF!="" && file_WplusH!="" && file_WminusH!="" && file_ZH!="" && file_ttH!="");
  Bool_t ZZBkgdIsPresent = file_ZZ4l!= "";
  Color_t colors[nSamples] = { kBlue, kGreen+2, kRed, kRed, kOrange+1, kMagenta-7, kYellow+1, kBlack };

  Color_t allColorsMP[nMatchHLepsStatuses][nHiggsSamples] = {
    { kBlue+2 , kGreen+3, kRed+2 , kRed+2 , kOrange+4, kMagenta+2 , kYellow+3 },
    { kBlue-4 , kGreen+2, kRed-4 , kRed-4 , kOrange+9, kMagenta   , kYellow+2 },
    { kBlue-7 , kGreen+1, kRed-7 , kRed-7 , kOrange+7, kMagenta-7 , kYellow+1 },
    { kBlue-9 , kGreen  , kRed-9 , kRed-9 , kOrange+1, kMagenta-9 , kYellow-6 },
    { kBlue-10, kGreen-9, kRed-10, kRed-10, kOrange-9, kMagenta-10, kYellow-8 },
    { kGray+1 , kGray+1 , kGray+1, kGray+1, kGray+1  , kGray+1    , kGray+1   },
  };
  Color_t allColorsPM1[nHiggsSamples][nMatchHLepsStatuses];
  for(int p=0; p<nHiggsSamples; p++){
    for(int m=0; m<nMatchHLepsStatuses; m++){
      allColorsPM1[p][m] = allColorsMP[m][p];
    }
  }
  Color_t allColorsPM2[nHiggsSamples][nMatchHLepsStatuses];
  for(int p=0; p<nHiggsSamples; p++){
    for(int m=0; m<nMatchHLepsStatuses; m++){
      if(m<3)
	allColorsPM2[p][m] = allColorsMP[m][p];
      else
	allColorsPM2[p][m] = allColorsMP[m+1][p];
    }
  }

  string genChannels[nGenChannels] = {
    "4mu",
    "4e",
    "2e2mu",
    "4tau",
    "2e2tau",
    "2mu2tau",
    "2L2X1",
    "4l",
    "2L2tau",
    "2L2X2",
    "all",
  };
  string recoChannels[nRecoChannels] = {
    "4mu",
    "4e",
    "2e2mu",
    "4l",
    "not4l1",
    "not4l2",
  };

  string nbZDaughtersFromH[4] = { "2", "1", "0", "ambig." };

  string WHdecays[nAssocWDecays] = {
    "H->ZZ->4l, W->X     ",
    "H->ZZ->4l, W->lnu   ",
  };
  string ZHdecays[nAssocZDecays] = {
    "H->ZZ->4l, Z->X     ",
    "H->ZZ->4l, Z->2l    ",
    "H->ZZ->2l2X, Z->2l  ",
  };
  string ttHdecays[nAssocttDecays] = {
    "H->ZZ->4l, tt->X    ",
    "H->ZZ->4l, tt->lX   ",
    "H->ZZ->4l, tt->2lX  ",
    "H->ZZ->2l2X, tt->2lX",
  };
  string assocDecayName1[nAssocDecays] = {
    "H#rightarrowZZ#rightarrow4l, W#rightarrowX",
    "H#rightarrowZZ#rightarrow4l, W#rightarrowl#nu",
    "H#rightarrowZZ#rightarrow4l, Z#rightarrowX",
    "H#rightarrowZZ#rightarrow4l, Z#rightarrow2l",
    "H#rightarrowZZ#rightarrow2l2X, Z#rightarrow2l",
    "H#rightarrowZZ#rightarrow4l, t#bar{t}#rightarrow0l+X",
    "H#rightarrowZZ#rightarrow4l, t#bar{t}#rightarrow1l+X",
    "H#rightarrowZZ#rightarrow4l, t#bar{t}#rightarrow2l+X",
    "H#rightarrowZZ#rightarrow2l2X, t#bar{t}#rightarrow2l+X",
  };
  string assocDecayName2[nAssocDecays] = {
    " W#rightarrowX",
    " W#rightarrowl#nu",
    " Z#rightarrowX",
    " Z#rightarrow2l",
    " ERROR",
    " t#bar{t}#rightarrow0l+X",
    " t#bar{t}#rightarrow1l+X",
    " t#bar{t}#rightarrow2l+X",
    " ERROR",
  };
  string assocDecayName3[nAssocDecays] = {
    "WH, H#rightarrow4l, W#rightarrowX",
    "WH, H#rightarrow4l, W#rightarrowl#nu",
    "ZH, H#rightarrow4l, Z#rightarrowX",
    "ZH, H#rightarrow4l, Z#rightarrow2l",
    "ZH, H#rightarrow2l2X, Z#rightarrow2l",
    "ttH, H#rightarrow4l, t#bar{t}#rightarrow0l+X",
    "ttH, H#rightarrow4l, t#bar{t}#rightarrow1l+X",
    "ttH, H#rightarrow4l, t#bar{t}#rightarrow2l+X",
    "ttH, H#rightarrow2l2X, t#bar{t}#rightarrow2l+X",
  };
  string assocDecayName4[nAssocDecays] = {
    "WH, W#rightarrowX",
    "WH, W#rightarrowl#nu",
    "ZH, Z#rightarrowX",
    "ZH, Z#rightarrow2l",
    "ERROR",
    "ttH, t#bar{t}#rightarrow0l+X",
    "ttH, t#bar{t}#rightarrow1l+X",
    "ttH, t#bar{t}#rightarrow2l+X",
    "ERROR",
  };
  Color_t assocDecayColor[nAssocDecays] = {
    kRed,
    kRed+1,
    kOrange+1, 
    kOrange+7,
    kOrange-8,
//     kViolet-4,
//     kViolet-5,
//     kViolet-6,
//     kViolet-8, 
    kMagenta-7, 
    kMagenta-3, 
    kMagenta+2, 
    kMagenta-8,
  };

  string varName[nVariables] = {
    "M4l",
    "MZ1",
    "MZ2",
    "KD",
    "Djet",
    "Pt4l",
    "NGenLep",
    "NGenLepInEtaPtAcc",
    "NGenLepNotInEtaPtAcc",
    "NGenHLepNotInEtaPtAcc",
    "NGenAssocLepNotInEtaPtAcc",
    "NGenLepMinusNGoodLep",
    "NGenLepInEtaPtAccMinusNGoodLep",
    "NExtraLep",
    "NExtraZ",
    "NJets",
  };
  Bool_t plotThisVar[2][nVariables] = {
    {1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,},
    {0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,},
  };
  string varLabel[nVariables] = {
    "m_{4l} (GeV)",
    "m_{Z_{1}} (GeV)",
    "m_{Z_{2}} (GeV)",
    "D_{bkg}^{kin}",
    "D_{jet}",
    "p_{T}^{4l} (GeV)",
    "# gen leptons",
    "# gen leptons in acceptance",
    "# gen leptons not in acceptance",
    "# gen leptons from H not in acceptance",
    "# gen associated leptons not in acceptance",
    "# gen leptons - # good leptons",
    "# gen leptons in acceptance - # good leptons",
    "# extra leptons",
    "# extra Z candidates",
    "# jets",
  };
  Int_t varNbin[nVariables] = { 
    100,  
    75,  
    75,  
    50,  
    50,  
    50,
    7,
    7,
    7,
    5,
    3,
    6,
    6,
    6,
    6,
    15,
  };
  Float_t varMin[nVariables] = {  
    50,  
    0,   
    0,   
    0,   
    0,   
    0,  
    0,   
    0, 
    0,
    0,   
    0,
    -3,
    -3, 
    0,
    0,
    0,
  };
  Float_t varMax[nVariables] = { 
    850, 
    150, 
    150,   
    1,   
    2, 
    500,
    7,
    7,
    7,
    5,
    3,
    3,
    3,
    6,
    6,
    15,
  };

  const Int_t n2DHist = 2;
  Int_t varXindex[n2DHist] = { 0, 2 };
  Int_t varYindex[n2DHist] = { 3, 3 };
  const Int_t nDecays = 5;
  string decayLabel[nDecays] = { ", all", ", H#rightarrow4l", ", H#rightarrow2l2X", ", H#rightarrow4l, 4 from H", ", H#rightarrow4l, <4 from H" };
  string decayInfix[nDecays] = { "HtoAny", "Hto4l", "Hto2l2X", "Hto4l4fromH", "Hto4lnot4fromH" };

  string matchHLepsKeys[nMatchHLepsStatuses] = {
    "4 matches",
    "3 matches",
    "2 matches",
    "1 match",
    "0 match",
    "ambiguous",
  };
  string matchHLepsInfix[nMatchHLepsStatuses] = {
    "4",
    "3",
    "2",
    "1",
    "0",
    "Ambig",
  };

  string matchAllLepsKeys[nMatchAllLepsStatuses] = {
    "4 from H",
    "3 from H, 1 assoc.",
    "2 from H, 2 assoc.",
    "< 4 matches",
    "ambiguous",
  };
  string matchAllLepsInfix[nMatchAllLepsStatuses] = {
    "4h0a",
    "3h1a",
    "2h2a",
    "Other",
    "Ambig",
  };

  string matchWH[nMatchWHStatuses] = {
    "H#rightarrowZZ#rightarrow4l, W#rightarrowX ; 4 from H",
    "H#rightarrowZZ#rightarrow4l, W#rightarrowl#nu ; 4 from H",
    "H#rightarrowZZ#rightarrow4l, W#rightarrowl#nu ; 3 from H, 1 from W",
    "< 4 matches",
    "ambiguous",
  };
  string matchZH[nMatchZHStatuses] = {
    "H#rightarrowZZ#rightarrow4l, Z#rightarrowX ; 4 from H",
    "H#rightarrowZZ#rightarrow4l, Z#rightarrow2l ; 4 from H",
    "H#rightarrowZZ#rightarrow4l, Z#rightarrow2l ; 3 from H, 1 from Z",
    "H#rightarrowZZ#rightarrow4l, Z#rightarrow2l ; 2 from H, 2 from Z",
    "H#rightarrowZZ#rightarrow2l2X, Z#rightarrow2l ; 2 from H, 2 from Z",
    "< 4 matches",
    "ambiguous",
  };
  string matchttH[nMatchttHStatuses] = {
    "H#rightarrowZZ#rightarrow4l, t#bar{t}#rightarrowX ; 4 from H",
    "H#rightarrowZZ#rightarrow4l, t#bar{t}#rightarrow1l+X ; 4 from H",
    "H#rightarrowZZ#rightarrow4l, t#bar{t}#rightarrow1l+X ; 3 from H, 1 from t#bar{t}",
    "H#rightarrowZZ#rightarrow4l, t#bar{t}#rightarrow2l+X ; 4 from H",
    "H#rightarrowZZ#rightarrow4l, t#bar{t}#rightarrow2l+X ; 3 from H, 1 from t#bar{t}",
    "H#rightarrowZZ#rightarrow4l, t#bar{t}#rightarrow2l+X ; 2 from H, 2 from t#bar{t}",
    "H#rightarrowZZ#rightarrow2l2X, t#bar{t}#rightarrow2l+X ; 2 from H, 2 from t#bar{t}",
    "< 4 matches",
    "ambiguous",
  };

  Color_t colorsMatchWH[nMatchWHStatuses] = {
    kGreen+2,
    kAzure,
    kAzure-4,
    kGray,
    kGray+1,
  };
  Color_t colorsMatchZH[nMatchZHStatuses] = {
    kGreen+2,
    kViolet+2,
    kViolet+1,
    kViolet-9,
    kOrange+7,
    kGray,
    kGray+1,
  };
  Color_t colorsMatchttH[nMatchttHStatuses] = {
    kGreen+2,
    kAzure,
    kAzure-4,
    kViolet+2,
    kViolet+1,
    kViolet-9,
    kOrange+7,
    kGray,
    kGray+1,
  };

  vector<string> basketLabel[7];
  basketLabel[0].push_back("all");
  basketLabel[0].push_back("#splitline{   0/1 jet}{4 leptons}");
  basketLabel[0].push_back("#splitline{   0/1 jet}{5 leptons}");
  basketLabel[0].push_back("#splitline{    0/1 jet}{#geq6 leptons}");
  basketLabel[0].push_back("#splitline{  #geq2 jets}{4 leptons}");
  basketLabel[0].push_back("#splitline{  #geq2 jets}{5 leptons}");
  basketLabel[0].push_back("#splitline{   #geq2 jets}{#geq6 leptons}");
  basketLabel[1].push_back("all");
  basketLabel[1].push_back("#splitline{0/1 jet}{4 leptons}");
  //basketLabel[1].push_back("#splitline{0/1 jet, 5 leptons}{E_{T}^{miss}>45, p_{T}^{l5}>20}");
  basketLabel[1].push_back("#splitline{0/1 jet, 5 leptons}{E_{T}^{miss}>45}");
  basketLabel[1].push_back("#splitline{0/1 jet, 5 leptons}{(else)}");
  basketLabel[1].push_back("#splitline{0/1 jet}{#geq6 leptons}");
  basketLabel[1].push_back("#splitline{#geq2 jets, 4 leptons}{D_{jet}>0.5}");
  basketLabel[1].push_back("#splitline{#geq2 jets, 4 leptons}{(else)}");
  basketLabel[1].push_back("#splitline{#geq2 jets}{5 leptons}");
  basketLabel[1].push_back("#splitline{#geq2 jets}{#geq6 leptons}");
  basketLabel[2].push_back("all");
  basketLabel[2].push_back("0j 4l");
  basketLabel[2].push_back("0j 5l");
  basketLabel[2].push_back("0j #geq6l");
  basketLabel[2].push_back("1j 0b 4l");
  basketLabel[2].push_back("1j 0b 5l");
  basketLabel[2].push_back("1j 0b #geq6l");
  basketLabel[2].push_back("1j 1b 4l");
  basketLabel[2].push_back("1j 1b 5l");
  basketLabel[2].push_back("1j 1b #geq6l");
  basketLabel[2].push_back("#geq2j 0b 4l D_{jet}");
  basketLabel[2].push_back("#geq2j 0b 4l VH-h");
  basketLabel[2].push_back("#geq2j 0b 4l else");
  basketLabel[2].push_back("#geq2j 0b 5l");
  basketLabel[2].push_back("#geq2j 0b #geq6l");
  basketLabel[2].push_back("#geq2j 1b 4l D_{jet}");
  basketLabel[2].push_back("#geq2j 1b 4l VH-h");
  basketLabel[2].push_back("#geq2j 1b 4l else");
  basketLabel[2].push_back("#geq2j 1b 5l");
  basketLabel[2].push_back("#geq2j 1b #geq6l");
  basketLabel[2].push_back("#geq2j #geq2b 4l VH-h");
  basketLabel[2].push_back("#geq2j #geq2b 4l else");
  basketLabel[2].push_back("#geq2j #geq2b 5l");
  basketLabel[2].push_back("#geq2j #geq2b #geq6l");
  basketLabel[3].push_back("all");
  basketLabel[3].push_back("Untagged");
  basketLabel[3].push_back("VBF tagged");
  basketLabel[3].push_back("#splitline{VH-leptonic}{    tagged}");
  basketLabel[3].push_back("#splitline{VH-hadronic}{     tagged}");
  basketLabel[3].push_back("ttH tagged");
  basketLabel[4].push_back("all");
  basketLabel[4].push_back("0j 4l");
  basketLabel[4].push_back("0j 5l");
  basketLabel[4].push_back("0j #geq6l");
  basketLabel[4].push_back("1j 0b 4l");
  basketLabel[4].push_back("1j 0b 5l");
  basketLabel[4].push_back("1j 0b #geq6l");
  basketLabel[4].push_back("1j 1b 4l");
  basketLabel[4].push_back("1j 1b 5l");
  basketLabel[4].push_back("1j 1b #geq6l");
  basketLabel[4].push_back("2j 0b 4l D_{jet}");
  basketLabel[4].push_back("2j 0b 4l VH-h");
  basketLabel[4].push_back("2j 0b 4l else");
  basketLabel[4].push_back("2j 0b 5l");
  basketLabel[4].push_back("2j 0b #geq6l");
  basketLabel[4].push_back("2j 1b 4l D_{jet}");
  basketLabel[4].push_back("2j 1b 4l VH-h");
  basketLabel[4].push_back("2j 1b 4l else");
  basketLabel[4].push_back("2j 1b 5l");
  basketLabel[4].push_back("2j 1b #geq6l");
  basketLabel[4].push_back("2j 2b 4l VH-h");
  basketLabel[4].push_back("2j 2b 4l else");
  basketLabel[4].push_back("2j 2b 5l");
  basketLabel[4].push_back("2j 2b #geq6l");
  basketLabel[4].push_back("#geq3j 0b 4l D_{jet}");
  basketLabel[4].push_back("#geq3j 0b 4l VH-h");
  basketLabel[4].push_back("#geq3j 0b 4l else");
  basketLabel[4].push_back("#geq3j 0b 5l");
  basketLabel[4].push_back("#geq3j 0b #geq6l");
  basketLabel[4].push_back("#geq3j 1b 4l D_{jet}");
  basketLabel[4].push_back("#geq3j 1b 4l VH-h");
  basketLabel[4].push_back("#geq3j 1b 4l else");
  basketLabel[4].push_back("#geq3j 1b 5l");
  basketLabel[4].push_back("#geq3j 1b #geq6l");
  basketLabel[4].push_back("#geq3j #geq2b 4l VH-h");
  basketLabel[4].push_back("#geq3j #geq2b 4l else");
  basketLabel[4].push_back("#geq3j #geq2b 5l");
  basketLabel[4].push_back("#geq3j #geq2b #geq6l");
  basketLabel[5].push_back("all");
  basketLabel[5].push_back("Untagged");
  basketLabel[5].push_back("1-jet tagged");
  basketLabel[5].push_back("VBF tagged");
  basketLabel[5].push_back("#splitline{VH-leptonic}{    tagged}");
  basketLabel[5].push_back("#splitline{VH-hadronic}{     tagged}");
  basketLabel[5].push_back("ttH tagged");
  basketLabel[6].push_back("all");
  basketLabel[6].push_back("4l 0j");
  basketLabel[6].push_back("4l 1j");
  basketLabel[6].push_back("4l #geq2j");
  basketLabel[6].push_back("5l 0j");
  basketLabel[6].push_back("5l 1j");
  basketLabel[6].push_back("5l #geq2j");
  basketLabel[6].push_back("#geq6l 0j");
  basketLabel[6].push_back("#geq6l 1j");
  basketLabel[6].push_back("#geq6l #geq2j");

  TFile* inputFile[nSamples];
  TTree* inputTree[nSamples];
  TH1F* hCounters[nSamples];
  Long64_t NGenEvt[nSamples];
  Float_t preselExp[nSamples];
  Double_t gen_sumWeights[nSamples];
  Float_t partialSampleWeight[nSamples];

  Int_t nRun;
  Long64_t nEvent;
  Int_t nLumi;
  Float_t PFMET;
  Short_t NRecoMu;
  Short_t NRecoEle;
  Short_t Nvtx;
  Short_t NObsInt;
  Float_t NTrueInt;
  Short_t trigWord;
  Float_t overallEventWeight;
  Float_t xsec;
  Short_t ZZsel;
  Float_t ZZMass;
  Float_t ZZPt;
  Float_t DiJetFisher;
  Float_t p0plus_VAJHU;
  Float_t bkg_VAMCFM;
  Float_t Z1Mass;
  Float_t Z2Mass;
  vector<Float_t> *CandLepPt = 0;
  vector<Float_t> *CandLepEta = 0;
  vector<Float_t> *CandLepPhi = 0;
  vector<Short_t> *CandLepId = 0;
  Short_t nExtraLep;
  vector<Float_t> *ExtraLepPt = 0;
  vector<Float_t> *ExtraLepEta = 0;
  vector<Float_t> *ExtraLepPhi = 0;
  vector<Float_t> *ExtraLepId = 0;
  Short_t nExtraZ;
  Short_t nJets;
  Short_t nJetsBTagged;
  vector<Float_t> *JetPt = 0;
  vector<Float_t> *JetEta = 0;
  vector<Float_t> *JetPhi = 0;
  vector<Float_t> *JetMass = 0;
  vector<Float_t> *JetIsBtagged = 0;
  vector<Float_t> *JetQGLikelihood = 0;
  Float_t DiJetMass;
  Float_t GenHPt;
  Float_t GenHRapidity;
  Float_t GenLep1Pt;
  Float_t GenLep1Eta;
  Float_t GenLep1Phi;
  Short_t GenLep1Id;
  Float_t GenLep2Pt;
  Float_t GenLep2Eta;
  Float_t GenLep2Phi;
  Short_t GenLep2Id;
  Float_t GenLep3Pt;
  Float_t GenLep3Eta;
  Float_t GenLep3Phi;
  Short_t GenLep3Id;
  Float_t GenLep4Pt;
  Float_t GenLep4Eta;
  Float_t GenLep4Phi;
  Short_t GenLep4Id;
  Float_t GenAssocLep1Pt;
  Float_t GenAssocLep1Eta;
  Float_t GenAssocLep1Phi;
  Short_t GenAssocLep1Id;
  Float_t GenAssocLep2Pt;
  Float_t GenAssocLep2Eta;
  Float_t GenAssocLep2Phi;
  Short_t GenAssocLep2Id;

  Int_t nbStored[nSamples][nGenChannels];
  Int_t nbHLepsAreInEtaAcc[nSamples][nGenChannels];
  Int_t nbHLepsAreInEtaPtAcc[nSamples][nGenChannels];
  Int_t nb4GenLeps[nSamples][nGenChannels];
  Int_t nb4RecoLeps[nSamples][nGenChannels];
  Int_t nbWithBCFullSel100G[nSamples][nGenChannels];
  Int_t nbPassTriggerG[nSamples][nGenChannels];
  Int_t nbPassTriggerGNo1E[nSamples][nGenChannels];
  Int_t nbPassTriggerGWithBCFullSel100G[nSamples][nGenChannels];
  Int_t nbPassTriggerGNo1EWithBCFullSel100G[nSamples][nGenChannels];
  Int_t nbHLepsAreInEtaPtAcc4RecoLeps[nSamples][nGenChannels];
  Int_t nbHLepsAreInEtaPtAccWithBCFullSel100G[nSamples][nGenChannels];
  Int_t nbHLepsAreInEtaPtAccPassTriggerG[nSamples][nGenChannels];
  Int_t nbHLepsAreInEtaPtAccPassTriggerGNo1E[nSamples][nGenChannels];
  Int_t nbHLepsAreInEtaPtAccPassTriggerGWithBCFullSel100G[nSamples][nGenChannels];
  Int_t nbHLepsAreInEtaPtAccPassTriggerGNo1EWithBCFullSel100G[nSamples][nGenChannels];
  Float_t yieldStored[nSamples][nGenChannels];
  Float_t yieldHLepsAreInEtaAcc[nSamples][nGenChannels];
  Float_t yieldHLepsAreInEtaPtAcc[nSamples][nGenChannels];
  Float_t yield4GenLeps[nSamples][nGenChannels];
  Float_t yield4RecoLeps[nSamples][nGenChannels];
  Float_t yieldWithBCFullSel100G[nSamples][nGenChannels];
  Float_t yieldPassTriggerG[nSamples][nGenChannels];
  Float_t yieldPassTriggerGNo1E[nSamples][nGenChannels];
  Float_t yieldPassTriggerGWithBCFullSel100G[nSamples][nGenChannels];
  Float_t yieldPassTriggerGNo1EWithBCFullSel100G[nSamples][nGenChannels];
  Float_t yieldHLepsAreInEtaPtAcc4RecoLeps[nSamples][nGenChannels];
  Float_t yieldHLepsAreInEtaPtAccWithBCFullSel100G[nSamples][nGenChannels];
  Float_t yieldHLepsAreInEtaPtAccPassTriggerG[nSamples][nGenChannels];
  Float_t yieldHLepsAreInEtaPtAccPassTriggerGNo1E[nSamples][nGenChannels];
  Float_t yieldHLepsAreInEtaPtAccPassTriggerGWithBCFullSel100G[nSamples][nGenChannels];
  Float_t yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCFullSel100G[nSamples][nGenChannels];
  TH1F* h1GenptHLepsAreInEtaAcc[nSamples][nGenChannels][4];
  TH1F* h1RecoptWithBCFullSel100GPassTriggerG[nSamples][nGenChannels][4];
  TH1F* h1GenetaHLepsAreInPtAcc[nSamples][nGenChannels][4];
  TH1F* h1RecoetaWithBCFullSel100GPassTriggerG[nSamples][nGenChannels][4];
  TH1F* h1GenHPt[nSamples];
  TH1F* h1GenHRapidity[nSamples];
  TH2F* h2GenHRapidityVsPt[nSamples];
  TH1F* h1NRecoLepHLepsAreInEtaPtAcc[nSamples][nGenChannels];
  TH1F* h1NRecoMuoHLepsAreInEtaPtAcc[nSamples][nGenChannels];
  TH1F* h1NRecoEleHLepsAreInEtaPtAcc[nSamples][nGenChannels];

  const Int_t nTTHGenStates = 100;
  const Int_t nTTHGenStatesUsed = 58;
  Int_t nbCheckTTH[nTTHGenStates];
  Float_t yieldCheckTTH[nTTHGenStates];
  for(int nt=0; nt<nTTHGenStates; nt++){
    nbCheckTTH[nt] = 0;
    yieldCheckTTH[nt] = 0.;
  }
  string labelCheckTTH[nTTHGenStates] = {
    "  All events       ",
    "  * H->4l          ",
    "    H->2l2t        ",
    "    H->4t          ",
    "    H->2l2X        ",
    "    H->2t2X        ",
    "  * TTH->4lX       ",
    "    TTH->5lX       ",
    "    TTH->6lX       ",
    "    trash          ",
    "  * TTH->4lX, 0e   ",
    "    TTH->4lX, 1e   ",
    "    TTH->4lX, 2e   ",
    "    TTH->4lX, 3e   ",
    "    TTH->4lX, 4e   ",
    "    TTH->4lX, trash",
    "  * TTH->4lX, 0m   ",
    "    TTH->4lX, 1m   ",
    "    TTH->4lX, 2m   ",
    "    TTH->4lX, 3m   ",
    "    TTH->4lX, 4m   ",
    "    TTH->4lX, trash",
    "  * TTH->4lX, 0t   ",
    "    TTH->4lX, 1t   ",
    "    TTH->4lX, 2t   ",
    "    TTH->4lX, 3t   ",
    "    TTH->4lX, 4t   ",
    "    TTH->4lX, trash",
    "  * TTH->4lX, 2+2- ",
    "    TTH->4lX, trash",
    "  * TTH->4lX, 4e   ",
    "    TTH->4lX, 4m   ",
    "    TTH->4lX, 4t   ",
    "    TTH->4lX, 2e2m ",
    "    TTH->4lX, 2e2t ",
    "    TTH->4lX, 2m2t ",
    "    TTH->4lX, else ",
    "  * H->4e          ",
    "    H->4m          ",
    "    H->4t          ",
    "    H->2e2m        ",
    "    H->2e2t        ",
    "    H->2m2t        ",
    "    H->2e2X        ",
    "    H->2m2X        ",
    "    H->2t2X        ",
    "    other          ",
    "  * TT->2eX        ",
    "    TT->2mX        ",
    "    TT->2tX        ",
    "    TT->emX        ",
    "    TT->etX        ",
    "    TT->mtX        ",
    "    TT->eX         ",
    "    TT->mX         ",
    "    TT->tX         ",
    "    TT->X          ",
    "    other          ",
    "","","","","","","","","","",
    "","","","","","","","","","",
    "","","","","","","","","","",
    "","","","","","","","","","",
    "","",
  };

  Int_t nbWithBC[nSamples][nRecoChannels];
  Int_t nbWithBCFullSel70[nSamples][nRecoChannels];
  Int_t nbWithBCFullSel100[nSamples][nRecoChannels];
  Int_t nbWithBCFullSel100Exactly4GoodLeps[nSamples][nRecoChannels];
  Int_t nbWithBCFullSel100AtLeast5GoodLeps[nSamples][nRecoChannels];
  Int_t nbWithBCFullSel100Exactly5GoodLeps[nSamples][nRecoChannels];
  Int_t nbWithBCFullSel100Exactly6GoodLeps[nSamples][nRecoChannels];
  Int_t nbWithBCFullSel100HLepsAreInEtaPtAcc[nSamples][nRecoChannels];
  Int_t nbWithBCFullSel100HLepsAreGood[nSamples][nRecoChannels];
  Int_t nbWithBCFullSel100All4LepRight[nSamples][nRecoChannels];
  Int_t nbWithBCFullSel100MatchHLeps[nMatchHLepsStatuses][nSamples][nRecoChannels];
  Int_t nbWithBCFullSel100MatchAllLeps[nMatchAllLepsStatuses][nSamples][nRecoChannels];
  Int_t nbWithBCFullSel100MatchWH[nMatchWHStatuses][nRecoChannels];
  Int_t nbWithBCFullSel100MatchZH[nMatchZHStatuses][nRecoChannels];
  Int_t nbWithBCFullSel100MatchttH[nMatchttHStatuses][nRecoChannels];
  Float_t yieldWithBC[nSamples][nRecoChannels];
  Float_t yieldWithBCFullSel70[nSamples][nRecoChannels];
  Float_t yieldWithBCFullSel100[nSamples][nRecoChannels];

  TH1F* hBCFullSel100[nVariables][nSamples][nRecoChannels];
  TH1F* hBCFullSel100MatchHLeps[nVariables][nMatchHLepsStatuses][nSamples][nRecoChannels];
  TH1F* hBCFullSel100MatchAllLeps[nVariables][nMatchAllLepsStatuses][nSamples][nRecoChannels];
  TH1F* hBCFullSel100MatchWH[nVariables][nMatchWHStatuses][nRecoChannels];
  TH1F* hBCFullSel100MatchZH[nVariables][nMatchZHStatuses][nRecoChannels];
  TH1F* hBCFullSel100MatchttH[nVariables][nMatchttHStatuses][nRecoChannels];

  TH2F* h2DBCFullSel100[n2DHist][nSamples][nRecoChannels];
  TH2F* h2DBCFullSel100Decays[n2DHist][nSamples][nDecays][nRecoChannels];

  TH1F* hBCFullSel100Baskets[nSamples][nRecoChannels];
  TH1F* hBCFullSel100BasketsAssocDecays[nAssocDecays][nRecoChannels];
  TH1F* hBCFullSel100MasswindowBaskets[nSamples][nRecoChannels];
  TH1F* hBCFullSel100MasswindowBasketsAssocDecays[nAssocDecays][nRecoChannels];

  for(int a=0; a<nAssocDecays; a++){
    for(int rc=0; rc<nRecoChannels; rc++){
      hBCFullSel100BasketsAssocDecays[a][rc] = new TH1F(Form("hBCFullSel100_baskets_assocDecays_%i_%s",a,recoChannels[rc].c_str()),Form(";;# exp. events in %.0f fb^{-1}",lumi),nBaskets,0,nBaskets);
      prepareBasketlabels(hBCFullSel100BasketsAssocDecays[a][rc], basketLabel[BASKETLIST]);
      hBCFullSel100BasketsAssocDecays[a][rc]->SetLineColor(assocDecayColor[a]);
      hBCFullSel100BasketsAssocDecays[a][rc]->SetFillColor(assocDecayColor[a]);
      hBCFullSel100MasswindowBasketsAssocDecays[a][rc] = new TH1F(Form("hBCFullSel100Masswindow_baskets_assocDecays_%i_%s",a,recoChannels[rc].c_str()),Form(";;# exp. events in %.0f fb^{-1}",lumi),nBaskets,0,nBaskets);
      prepareBasketlabels(hBCFullSel100MasswindowBasketsAssocDecays[a][rc], basketLabel[BASKETLIST]);
      hBCFullSel100MasswindowBasketsAssocDecays[a][rc]->SetLineColor(assocDecayColor[a]);
      hBCFullSel100MasswindowBasketsAssocDecays[a][rc]->SetFillColor(assocDecayColor[a]);
    }
  }

  Int_t nbTotalWH[nAssocWDecays] = {0,0};
  Int_t nbHLepsAreInEtaPtAccWH[nAssocWDecays] = {0,0};
  Int_t nbHLepsAreGoodWH[nAssocWDecays] = {0,0};
  Int_t nbAll4LepRightWH[nAssocWDecays] = {0,0};
  Int_t nbZ1DaughtersFromHWH[nAssocWDecays][4]; for(int i=0; i<nAssocWDecays; i++) for(int j=0; j<4; j++) nbZ1DaughtersFromHWH[i][j] = 0;
  Int_t nbZ2DaughtersFromHWH[nAssocWDecays][4]; for(int i=0; i<nAssocWDecays; i++) for(int j=0; j<4; j++) nbZ2DaughtersFromHWH[i][j] = 0;

  Int_t nbTotalZH[nAssocZDecays] = {0,0,0};
  Int_t nbHLepsAreInEtaPtAccZH[nAssocZDecays] = {0,0,0};
  Int_t nbHLepsAreGoodZH[nAssocZDecays] = {0,0,0};
  Int_t nbAll4LepRightZH[nAssocZDecays] = {0,0,0};
  Int_t nbZ1DaughtersFromHZH[nAssocZDecays][4]; for(int i=0; i<nAssocZDecays; i++) for(int j=0; j<4; j++) nbZ1DaughtersFromHZH[i][j] = 0;
  Int_t nbZ2DaughtersFromHZH[nAssocZDecays][4]; for(int i=0; i<nAssocZDecays; i++) for(int j=0; j<4; j++) nbZ2DaughtersFromHZH[i][j] = 0;

  Int_t nbTotalttH[nAssocttDecays] = {0,0,0,0};
  Int_t nbHLepsAreInEtaPtAccttH[nAssocttDecays] = {0,0,0,0};
  Int_t nbHLepsAreGoodttH[nAssocttDecays] = {0,0,0,0};
  Int_t nbAll4LepRightttH[nAssocttDecays] = {0,0,0,0};
  Int_t nbZ1DaughtersFromHttH[nAssocttDecays][4]; for(int i=0; i<nAssocttDecays; i++) for(int j=0; j<4; j++) nbZ1DaughtersFromHttH[i][j] = 0;
  Int_t nbZ2DaughtersFromHttH[nAssocttDecays][4]; for(int i=0; i<nAssocttDecays; i++) for(int j=0; j<4; j++) nbZ2DaughtersFromHttH[i][j] = 0;



  // ------------------------------------------------------------
  // ---------------------- Processing --------------------------
  // ------------------------------------------------------------

  for(int p=0; p<nSamples; p++){
    if(!isPresent[p]) continue;

    cout<<prodName[p]<<endl;
    txtOut<<prodName[p]<<endl;

    for(int gc=0; gc<nGenChannels; gc++){

      string suffix = prodName[p]+"_"+genChannels[gc];
      string title = prodName[p]+", gen. channel "+genChannels[gc];

      nbStored[p][gc] = 0;
      nbHLepsAreInEtaAcc[p][gc] = 0;
      nbHLepsAreInEtaPtAcc[p][gc] = 0;
      nb4GenLeps[p][gc] = 0;
      nb4RecoLeps[p][gc] = 0;
      nbWithBCFullSel100G[p][gc] = 0;
      nbPassTriggerG[p][gc] = 0;
      nbPassTriggerGNo1E[p][gc] = 0;
      nbPassTriggerGWithBCFullSel100G[p][gc] = 0;
      nbPassTriggerGNo1EWithBCFullSel100G[p][gc] = 0;
      nbHLepsAreInEtaPtAcc4RecoLeps[p][gc] = 0;
      nbHLepsAreInEtaPtAccWithBCFullSel100G[p][gc] = 0;
      nbHLepsAreInEtaPtAccPassTriggerG[p][gc] = 0;
      nbHLepsAreInEtaPtAccPassTriggerGNo1E[p][gc] = 0;
      nbHLepsAreInEtaPtAccPassTriggerGWithBCFullSel100G[p][gc] = 0;
      nbHLepsAreInEtaPtAccPassTriggerGNo1EWithBCFullSel100G[p][gc] = 0;
      yieldStored[p][gc] = 0.;
      yieldHLepsAreInEtaAcc[p][gc] = 0.;
      yieldHLepsAreInEtaPtAcc[p][gc] = 0.;
      yield4GenLeps[p][gc] = 0.;
      yield4RecoLeps[p][gc] = 0.;
      yieldWithBCFullSel100G[p][gc] = 0.;
      yieldPassTriggerG[p][gc] = 0.;
      yieldPassTriggerGNo1E[p][gc] = 0.;
      yieldPassTriggerGWithBCFullSel100G[p][gc] = 0.;
      yieldPassTriggerGNo1EWithBCFullSel100G[p][gc] = 0.;
      yieldHLepsAreInEtaPtAcc4RecoLeps[p][gc] = 0.;
      yieldHLepsAreInEtaPtAccWithBCFullSel100G[p][gc] = 0.;
      yieldHLepsAreInEtaPtAccPassTriggerG[p][gc] = 0.;
      yieldHLepsAreInEtaPtAccPassTriggerGNo1E[p][gc] = 0.;
      yieldHLepsAreInEtaPtAccPassTriggerGWithBCFullSel100G[p][gc] = 0.;
      yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCFullSel100G[p][gc] = 0.;

      for(int sl=0; sl<4; sl++){
	h1GenptHLepsAreInEtaAcc[p][gc][sl] = new TH1F(Form("h1GenptHLepsAreInEtaAcc_%s_%i",suffix.c_str(),sl),(title+";p_{T};").c_str(),80,0,80);
	h1RecoptWithBCFullSel100GPassTriggerG[p][gc][sl] = new TH1F(Form("h1RecoptWithBCFullSel100GPassTriggerG_%s_%i",suffix.c_str(),sl),(title+";p_{T};").c_str(),100,0,100);
	h1GenetaHLepsAreInPtAcc[p][gc][sl] = new TH1F(Form("h1GenetaHLepsAreInPtAcc_%s_%i",suffix.c_str(),sl),(title+";|#eta|;").c_str(),100,0,5);
	h1RecoetaWithBCFullSel100GPassTriggerG[p][gc][sl] = new TH1F(Form("h1RecoetaWithBCFullSel100GPassTriggerG_%s_%i",suffix.c_str(),sl),(title+";|#eta|;").c_str(),100,0,5);
      }

      h1NRecoLepHLepsAreInEtaPtAcc[p][gc] = new TH1F(Form("h1NRecoLepHLepsAreInEtaPtAcc_%s",suffix.c_str()),(title+";nb. of reco. leptons;").c_str(),10,0,10);
      h1NRecoMuoHLepsAreInEtaPtAcc[p][gc] = new TH1F(Form("h1NRecoMuoHLepsAreInEtaPtAcc_%s",suffix.c_str()),(title+";nb. of reco. muons;").c_str(),10,0,10);
      h1NRecoEleHLepsAreInEtaPtAcc[p][gc] = new TH1F(Form("h1NRecoEleHLepsAreInEtaPtAcc_%s",suffix.c_str()),(title+";nb. of reco. electrons;").c_str(),10,0,10);

    }

    h1GenHPt[p] = new TH1F(Form("h1GenHPt_%s",prodName[p].c_str()),";p_{T}^{H};",100,0,100);
    h1GenHRapidity[p] = new TH1F(Form("h1GenHRapidity_%s",prodName[p].c_str()),";y^{H};",100,0,5);
    h2GenHRapidityVsPt[p] = new TH2F(Form("h2GenHRapidityVsPt_%s",prodName[p].c_str()),";p_{T}^{H};y^{H}",20,0,100,20,0,5);

    for(int rc=0; rc<nRecoChannels; rc++){
      
      string suffix = "_"+prodName[p]+"_"+recoChannels[rc];
      string title = prodName[p]+", channel "+recoChannels[rc];

      nbWithBC[p][rc] = 0;
      nbWithBCFullSel70[p][rc] = 0;
      nbWithBCFullSel100[p][rc] = 0;
      nbWithBCFullSel100Exactly4GoodLeps[p][rc] = 0;
      nbWithBCFullSel100AtLeast5GoodLeps[p][rc] = 0;
      nbWithBCFullSel100Exactly5GoodLeps[p][rc] = 0;
      nbWithBCFullSel100Exactly6GoodLeps[p][rc] = 0;
      nbWithBCFullSel100HLepsAreInEtaPtAcc[p][rc] = 0;
      nbWithBCFullSel100HLepsAreGood[p][rc] = 0;
      nbWithBCFullSel100All4LepRight[p][rc] = 0;
      for(int m=0; m<nMatchHLepsStatuses; m++) nbWithBCFullSel100MatchHLeps[m][p][rc] = 0;
      for(int m=0; m<nMatchAllLepsStatuses; m++) nbWithBCFullSel100MatchAllLeps[m][p][rc] = 0;
      if(prodName[p]=="WplusH") for(int m=0; m<nMatchWHStatuses; m++) nbWithBCFullSel100MatchWH[m][rc] = 0;
      if(prodName[p]=="ZH") for(int m=0; m<nMatchZHStatuses; m++) nbWithBCFullSel100MatchZH[m][rc] = 0;
      if(prodName[p]=="ttH") for(int m=0; m<nMatchttHStatuses; m++) nbWithBCFullSel100MatchttH[m][rc] = 0;

      yieldWithBC[p][rc] = 0.;
      yieldWithBCFullSel70[p][rc] = 0.;
      yieldWithBCFullSel100[p][rc] = 0.;

      for(int v=0; v<nVariables; v++){
	if(!plotThisVar[VARIABLELIST][v]) continue;
	hBCFullSel100[v][p][rc] = new TH1F(("hBCFullSel100_"+varName[v]+suffix).c_str(),Form("%s;%s;# exp. events in %.0f fb^{-1}",title.c_str(),varLabel[v].c_str(),lumi),varNbin[v],varMin[v],varMax[v]);
	for(int m=0; m<nMatchHLepsStatuses; m++)
	  hBCFullSel100MatchHLeps[v][m][p][rc] = new TH1F(("hBCFullSel100MatchHLeps_"+varName[v]+"_"+matchHLepsInfix[m]+suffix).c_str(),Form("%s;%s;# exp. events in %.0f fb^{-1}",title.c_str(),varLabel[v].c_str(),lumi),varNbin[v],varMin[v],varMax[v]);
	for(int m=0; m<nMatchAllLepsStatuses; m++)
	  hBCFullSel100MatchAllLeps[v][m][p][rc] = new TH1F(("hBCFullSel100MatchAllLeps_"+varName[v]+"_"+matchAllLepsInfix[m]+suffix).c_str(),Form("%s;%s;# exp. events in %.0f fb^{-1}",title.c_str(),varLabel[v].c_str(),lumi),varNbin[v],varMin[v],varMax[v]);
	if(prodName[p]=="WplusH")
	  for(int m=0; m<nMatchWHStatuses; m++)
	    hBCFullSel100MatchWH[v][m][rc] = new TH1F(Form("hBCFullSel100MatchWH_%s_%i_%s",varName[v].c_str(),m,recoChannels[rc].c_str()),Form("%s;%s;# exp. events in %.0f fb^{-1}",title.c_str(),varLabel[v].c_str(),lumi),varNbin[v],varMin[v],varMax[v]);
	if(prodName[p]=="ZH")
	  for(int m=0; m<nMatchZHStatuses; m++)
	    hBCFullSel100MatchZH[v][m][rc] = new TH1F(Form("hBCFullSel100MatchZH_%s_%i_%s",varName[v].c_str(),m,recoChannels[rc].c_str()),Form("%s;%s;# exp. events in %.0f fb^{-1}",title.c_str(),varLabel[v].c_str(),lumi),varNbin[v],varMin[v],varMax[v]);
	if(prodName[p]=="ttH")
	  for(int m=0; m<nMatchttHStatuses; m++)
	    hBCFullSel100MatchttH[v][m][rc] = new TH1F(Form("hBCFullSel100MatchttH_%s_%i_%s",varName[v].c_str(),m,recoChannels[rc].c_str()),Form("%s;%s;# exp. events in %.0f fb^{-1}",title.c_str(),varLabel[v].c_str(),lumi),varNbin[v],varMin[v],varMax[v]);
      }

      for(int v2=0; v2<n2DHist; v2++){
	h2DBCFullSel100[v2][p][rc] = new TH2F(("h2DBCFullSel100_"+varName[varXindex[v2]]+"_"+varName[varYindex[v2]]+suffix).c_str(),(title+";"+varLabel[varXindex[v2]]+";"+varLabel[varYindex[v2]]).c_str(),varNbin[varXindex[v2]],varMin[varXindex[v2]],varMax[varXindex[v2]],varNbin[varYindex[v2]],varMin[varYindex[v2]],varMax[varYindex[v2]]);
	for(int d=0; d<nDecays; d++){
	  h2DBCFullSel100Decays[v2][p][d][rc] = new TH2F(("h2DBCFullSel100_"+decayInfix[d]+"_"+varName[varXindex[v2]]+"_"+varName[varYindex[v2]]+suffix).c_str(),(title+";"+varLabel[varXindex[v2]]+";"+varLabel[varYindex[v2]]).c_str(),varNbin[varXindex[v2]],varMin[varXindex[v2]],varMax[varXindex[v2]],varNbin[varYindex[v2]],varMin[varYindex[v2]],varMax[varYindex[v2]]);
	}
      }
      
      hBCFullSel100Baskets[p][rc] = new TH1F(("hBCFullSel100_baskets_"+suffix).c_str(),Form("%s;;# exp. events in %.0f fb^{-1}",title.c_str(),lumi),nBaskets,0,nBaskets);
      prepareBasketlabels(hBCFullSel100Baskets[p][rc], basketLabel[BASKETLIST]);
      hBCFullSel100Baskets[p][rc]->SetLineColor(colors[p]);
      hBCFullSel100Baskets[p][rc]->SetFillColor(colors[p]);
      hBCFullSel100MasswindowBaskets[p][rc] = new TH1F(("hBCFullSel100Masswindow_baskets_"+suffix).c_str(),Form("%s;;# exp. events in %.0f fb^{-1}",title.c_str(),lumi),nBaskets,0,nBaskets);
      prepareBasketlabels(hBCFullSel100MasswindowBaskets[p][rc], basketLabel[BASKETLIST]);
      hBCFullSel100MasswindowBaskets[p][rc]->SetLineColor(colors[p]);
      hBCFullSel100MasswindowBaskets[p][rc]->SetFillColor(colors[p]);
      
    }

    inputFile[p] = TFile::Open((inDir+tagIn+"/"+fileName[p]).c_str());
    hCounters[p] = (TH1F*)inputFile[p]->Get("ZZTree/Counters");
    NGenEvt[p] = (Long64_t)hCounters[p]->GetBinContent(1);
    cout<<"  number of generated events: "<<NGenEvt[p]<<endl;
    gen_sumWeights[p] = (Long64_t)hCounters[p]->GetBinContent(40);
    partialSampleWeight[p] = lumi * 1000 / gen_sumWeights[p] ;

    inputTree[p] = (TTree*)inputFile[p]->Get("ZZTree/candTree");
    inputTree[p]->SetBranchAddress("RunNumber", &nRun);
    inputTree[p]->SetBranchAddress("EventNumber", &nEvent);
    inputTree[p]->SetBranchAddress("LumiNumber", &nLumi);
    inputTree[p]->SetBranchAddress("PFMET", &PFMET);
    inputTree[p]->SetBranchAddress("NRecoMu", &NRecoMu);
    inputTree[p]->SetBranchAddress("NRecoEle", &NRecoEle);
    inputTree[p]->SetBranchAddress("Nvtx", &Nvtx);
    inputTree[p]->SetBranchAddress("NObsInt", &NObsInt);
    inputTree[p]->SetBranchAddress("NTrueInt", &NTrueInt);
    inputTree[p]->SetBranchAddress("trigWord", &trigWord);
    inputTree[p]->SetBranchAddress("overallEventWeight", &overallEventWeight);
    inputTree[p]->SetBranchAddress("xsec", &xsec);
    inputTree[p]->SetBranchAddress("ZZsel", &ZZsel);
    inputTree[p]->SetBranchAddress("ZZMass", &ZZMass);
    inputTree[p]->SetBranchAddress("ZZPt", &ZZPt);
    inputTree[p]->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU);
    inputTree[p]->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM);
    inputTree[p]->SetBranchAddress("Z1Mass", &Z1Mass);
    inputTree[p]->SetBranchAddress("Z2Mass", &Z2Mass);
    inputTree[p]->SetBranchAddress("LepPt", &CandLepPt);
    inputTree[p]->SetBranchAddress("LepEta", &CandLepEta);
    inputTree[p]->SetBranchAddress("LepPhi", &CandLepPhi);
    inputTree[p]->SetBranchAddress("LepLepId", &CandLepId);
    inputTree[p]->SetBranchAddress("nExtraLep", &nExtraLep);
    inputTree[p]->SetBranchAddress("ExtraLepPt", &ExtraLepPt);
    inputTree[p]->SetBranchAddress("ExtraLepEta", &ExtraLepEta);
    inputTree[p]->SetBranchAddress("ExtraLepPhi", &ExtraLepPhi);
    inputTree[p]->SetBranchAddress("ExtraLepLepId", &ExtraLepId);
    inputTree[p]->SetBranchAddress("nExtraZ", &nExtraZ);
    inputTree[p]->SetBranchAddress("nCleanedJetsPt30", &nJets);
    inputTree[p]->SetBranchAddress("nCleanedJetsPt30BTagged", &nJetsBTagged);
    inputTree[p]->SetBranchAddress("JetPt", &JetPt);
    inputTree[p]->SetBranchAddress("JetEta", &JetEta);
    inputTree[p]->SetBranchAddress("JetPhi", &JetPhi);
    inputTree[p]->SetBranchAddress("JetMass", &JetMass);
    inputTree[p]->SetBranchAddress("JetIsBtagged", &JetIsBtagged);
    inputTree[p]->SetBranchAddress("JetQGLikelihood", &JetQGLikelihood);
    inputTree[p]->SetBranchAddress("DiJetMass", &DiJetMass);
    inputTree[p]->SetBranchAddress("DiJetFisher", &DiJetFisher);
    inputTree[p]->SetBranchAddress("GenHPt", &GenHPt);
    inputTree[p]->SetBranchAddress("GenHRapidity", &GenHRapidity);
    inputTree[p]->SetBranchAddress("GenLep1Pt", &GenLep1Pt);
    inputTree[p]->SetBranchAddress("GenLep1Eta", &GenLep1Eta);
    inputTree[p]->SetBranchAddress("GenLep1Phi", &GenLep1Phi);
    inputTree[p]->SetBranchAddress("GenLep1Id", &GenLep1Id);
    inputTree[p]->SetBranchAddress("GenLep2Pt", &GenLep2Pt);
    inputTree[p]->SetBranchAddress("GenLep2Eta", &GenLep2Eta);
    inputTree[p]->SetBranchAddress("GenLep2Phi", &GenLep2Phi);
    inputTree[p]->SetBranchAddress("GenLep2Id", &GenLep2Id);
    inputTree[p]->SetBranchAddress("GenLep3Pt", &GenLep3Pt);
    inputTree[p]->SetBranchAddress("GenLep3Eta", &GenLep3Eta);
    inputTree[p]->SetBranchAddress("GenLep3Phi", &GenLep3Phi);
    inputTree[p]->SetBranchAddress("GenLep3Id", &GenLep3Id);
    inputTree[p]->SetBranchAddress("GenLep4Pt", &GenLep4Pt);
    inputTree[p]->SetBranchAddress("GenLep4Eta", &GenLep4Eta);
    inputTree[p]->SetBranchAddress("GenLep4Phi", &GenLep4Phi);
    inputTree[p]->SetBranchAddress("GenLep4Id", &GenLep4Id);
    inputTree[p]->SetBranchAddress("GenAssocLep1Pt", &GenAssocLep1Pt);
    inputTree[p]->SetBranchAddress("GenAssocLep1Eta", &GenAssocLep1Eta);
    inputTree[p]->SetBranchAddress("GenAssocLep1Phi", &GenAssocLep1Phi);
    inputTree[p]->SetBranchAddress("GenAssocLep1Id", &GenAssocLep1Id);
    inputTree[p]->SetBranchAddress("GenAssocLep2Pt", &GenAssocLep2Pt);
    inputTree[p]->SetBranchAddress("GenAssocLep2Eta", &GenAssocLep2Eta);
    inputTree[p]->SetBranchAddress("GenAssocLep2Phi", &GenAssocLep2Phi);
    inputTree[p]->SetBranchAddress("GenAssocLep2Id", &GenAssocLep2Id);

    preselExp[p] = lumi * 1000 * xsec ;

    Long64_t entries = inputTree[p]->GetEntries();

    for (Long64_t z=0; z<entries; ++z){

      if(DEBUG && z>1000) break;

      printStatus(z,20000,entries,"entries");

      inputTree[p]->GetEntry(z);

      Double_t eventWeight = partialSampleWeight[p] * xsec * overallEventWeight ;

      Bool_t FullSel70 = ZZsel>=90;
      Bool_t FullSel100 = ZZsel>=100;
      Bool_t passTrigger = test_bit_16(trigWord,0);
      Bool_t passTriggerNo1E = test_bit_16(trigWord,8);

      Short_t GenHLepId[4] = {GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id};
      Float_t GenHLepPt[4] = {GenLep1Pt,GenLep2Pt,GenLep3Pt,GenLep4Pt};
      Float_t GenHLepEta[4] = {GenLep1Eta,GenLep2Eta,GenLep3Eta,GenLep4Eta};
      Float_t GenHLepPhi[4] = {GenLep1Phi,GenLep2Phi,GenLep3Phi,GenLep4Phi};
      Short_t GenAssocLepId[2] = {GenAssocLep1Id,GenAssocLep2Id};
      Float_t GenAssocLepPt[2] = {GenAssocLep1Pt,GenAssocLep2Pt};
      Float_t GenAssocLepEta[2] = {GenAssocLep1Eta,GenAssocLep2Eta};
      Float_t GenAssocLepPhi[2] = {GenAssocLep1Phi,GenAssocLep2Phi};


      // ---------------------- reco decay channel --------------------------

      Int_t nCandEle = 0;
      Int_t nCandMu = 0;
      for(Int_t iCandLep=0; iCandLep<4; iCandLep++){
	if(abs(CandLepId->at(iCandLep))==11) nCandEle++;
	if(abs(CandLepId->at(iCandLep))==13) nCandMu++;
      }
      Int_t rc = 4;
      if     (nCandEle==0 && nCandMu==4) rc = 0;
      else if(nCandEle==4 && nCandMu==0) rc = 1;
      else if(nCandEle==2 && nCandMu==2) rc = 2;
      //if(rc==4) cout<<"error in channel definition"<<endl;
      Int_t rc2 = 5;
      if(rc==0||rc==1||rc==2) rc2 = 3;
      Int_t rc3 = 6;


      // ---------------------- Gen lepton counting --------------------------
 
      Int_t nGenHLep = 0;
      Int_t nGenHLepInEtaAcc   = 0;
      Int_t nGenHLepInPtAcc    = 0;
      Int_t nGenHLepInEtaPtAcc = 0;
      Int_t nGenHEle = 0;
      Int_t nGenHMu  = 0;
      Int_t nGenHTau = 0;
      Int_t nGenHLEP = 0;
      Int_t nGenAssocLep = 0;
      Int_t nGenAssocLepInEtaAcc   = 0;
      Int_t nGenAssocLepInPtAcc    = 0;
      Int_t nGenAssocLepInEtaPtAcc = 0;
      Int_t nGenAssocEle = 0;
      Int_t nGenAssocMu  = 0;
      Int_t nGenAssocTau = 0;
      Int_t nGenAssocLEP = 0;
      Int_t nGenLEPPlus  = 0;
      Int_t nGenLEPMinus = 0;
      Bool_t  GenHLepIsInEtaAcc  [4];
      Bool_t  GenHLepIsInPtAcc   [4];
      Bool_t  GenHLepIsInEtaPtAcc[4];
      Bool_t  GenAssocLepIsInEtaAcc  [4];
      Bool_t  GenAssocLepIsInPtAcc   [4];
      Bool_t  GenAssocLepIsInEtaPtAcc[4];
      for(Int_t iGenHLep=0; iGenHLep<4; iGenHLep++){
	if(abs(GenHLepId[iGenHLep])==11){
	  nGenHEle++; nGenHLep++; nGenHLEP++;
	  if(GenHLepId[iGenHLep]>0) nGenLEPPlus++; else nGenLEPMinus++;
	  GenHLepIsInEtaAcc[iGenHLep] = fabs(GenHLepEta[iGenHLep])<2.5 ;
	  GenHLepIsInPtAcc [iGenHLep] = GenHLepPt[iGenHLep]>7. ;
	  GenHLepIsInEtaPtAcc[iGenHLep] = GenHLepIsInEtaAcc[iGenHLep] && GenHLepIsInPtAcc[iGenHLep] ;
	  if(GenHLepIsInEtaAcc  [iGenHLep]) nGenHLepInEtaAcc  ++;
	  if(GenHLepIsInPtAcc   [iGenHLep]) nGenHLepInPtAcc   ++;
	  if(GenHLepIsInEtaPtAcc[iGenHLep]) nGenHLepInEtaPtAcc++;
	}else if(abs(GenHLepId[iGenHLep])==13){
	  nGenHMu ++; nGenHLep++; nGenHLEP++;
	  if(GenHLepId[iGenHLep]>0) nGenLEPPlus++; else nGenLEPMinus++;
	  GenHLepIsInEtaAcc[iGenHLep] = fabs(GenHLepEta[iGenHLep])<2.4 ;
	  GenHLepIsInPtAcc [iGenHLep] = GenHLepPt[iGenHLep]>5. ;
	  GenHLepIsInEtaPtAcc[iGenHLep] = GenHLepIsInEtaAcc[iGenHLep] && GenHLepIsInPtAcc[iGenHLep] ;
	  if(GenHLepIsInEtaAcc  [iGenHLep]) nGenHLepInEtaAcc  ++;
	  if(GenHLepIsInPtAcc   [iGenHLep]) nGenHLepInPtAcc   ++;
	  if(GenHLepIsInEtaPtAcc[iGenHLep]) nGenHLepInEtaPtAcc++;
	}else if(abs(GenHLepId[iGenHLep])==15){
	  nGenHTau++; nGenHLEP++;
	  if(GenHLepId[iGenHLep]>0) nGenLEPPlus++; else nGenLEPMinus++;
	}
      }
      for(Int_t iGenAssocLep=0; iGenAssocLep<2; iGenAssocLep++){
	if(abs(GenAssocLepId[iGenAssocLep])==11){
	  nGenAssocEle++; nGenAssocLep++; nGenAssocLEP++;
	  if(GenAssocLepId[iGenAssocLep]>0) nGenLEPPlus++; else nGenLEPMinus++;
	  GenAssocLepIsInEtaAcc[iGenAssocLep] = fabs(GenAssocLepEta[iGenAssocLep])<2.5 ;
	  GenAssocLepIsInPtAcc [iGenAssocLep] = GenAssocLepPt[iGenAssocLep]>7. ;
	  GenAssocLepIsInEtaPtAcc[iGenAssocLep] = GenAssocLepIsInEtaAcc[iGenAssocLep] && GenAssocLepIsInPtAcc[iGenAssocLep] ;
	  if(GenAssocLepIsInEtaAcc  [iGenAssocLep]) nGenAssocLepInEtaAcc  ++;
	  if(GenAssocLepIsInPtAcc   [iGenAssocLep]) nGenAssocLepInPtAcc   ++;
	  if(GenAssocLepIsInEtaPtAcc[iGenAssocLep]) nGenAssocLepInEtaPtAcc++;
	}else if(abs(GenAssocLepId[iGenAssocLep])==13){
	  nGenAssocMu ++; nGenAssocLep++; nGenAssocLEP++;
	  if(GenAssocLepId[iGenAssocLep]>0) nGenLEPPlus++; else nGenLEPMinus++;
	  GenAssocLepIsInEtaAcc[iGenAssocLep] = fabs(GenAssocLepEta[iGenAssocLep])<2.4 ;
	  GenAssocLepIsInPtAcc [iGenAssocLep] = GenAssocLepPt[iGenAssocLep]>5. ;
	  GenAssocLepIsInEtaPtAcc[iGenAssocLep] = GenAssocLepIsInEtaAcc[iGenAssocLep] && GenAssocLepIsInPtAcc[iGenAssocLep] ;
	  if(GenAssocLepIsInEtaAcc  [iGenAssocLep]) nGenAssocLepInEtaAcc  ++;
	  if(GenAssocLepIsInPtAcc   [iGenAssocLep]) nGenAssocLepInPtAcc   ++;
	  if(GenAssocLepIsInEtaPtAcc[iGenAssocLep]) nGenAssocLepInEtaPtAcc++;
	}else if(abs(GenAssocLepId[iGenAssocLep])==15){
	  nGenAssocTau++; nGenAssocLEP++;
	  if(GenAssocLepId[iGenAssocLep]>0) nGenLEPPlus++; else nGenLEPMinus++;
	}
      }
      Int_t nGenLep = nGenHLep + nGenAssocLep;
      Int_t nGenLepInEtaAcc   = nGenHLepInEtaAcc   + nGenAssocLepInEtaAcc  ;
      Int_t nGenLepInPtAcc    = nGenHLepInPtAcc    + nGenAssocLepInPtAcc   ;
      Int_t nGenLepInEtaPtAcc = nGenHLepInEtaPtAcc + nGenAssocLepInEtaPtAcc;
      Int_t nGenEle = nGenHEle + nGenAssocEle;
      Int_t nGenMu  = nGenHMu  + nGenAssocMu ;
      Int_t nGenTau = nGenHTau + nGenAssocTau;
      Int_t nGenLEP = nGenHLEP + nGenAssocLEP;
 
      Int_t currentGenHDecay = -1;
      if     (nGenHMu ==4               ) currentGenHDecay = 0;
      else if(nGenHEle==4               ) currentGenHDecay = 1;
      else if(nGenHEle==2 && nGenHMu ==2) currentGenHDecay = 2;
      else if(nGenHTau==4               ) currentGenHDecay = 3;
      else if(nGenHEle==2 && nGenHTau==2) currentGenHDecay = 4;
      else if(nGenHMu ==2 && nGenHTau==2) currentGenHDecay = 5;
      else                                currentGenHDecay = 6;
      Int_t gc1 = currentGenHDecay;
      Int_t gc2 = 9;
      if(0<=currentGenHDecay && currentGenHDecay<=2) gc2 = 7;
      if(3<=currentGenHDecay && currentGenHDecay<=5) gc2 = 8;
      Int_t gc3 = 10;


      // ---------------------- Control counters --------------------------

      nbStored[p][gc1]++;
      nbStored[p][gc2]++;
      nbStored[p][gc3]++;
      yieldStored[p][gc1] += eventWeight;
      yieldStored[p][gc2] += eventWeight;
      yieldStored[p][gc3] += eventWeight;
      if(nGenHLepInEtaAcc==4){
	nbHLepsAreInEtaAcc[p][gc1]++;
	nbHLepsAreInEtaAcc[p][gc2]++;
	nbHLepsAreInEtaAcc[p][gc3]++;
	yieldHLepsAreInEtaAcc[p][gc1] += eventWeight;
	yieldHLepsAreInEtaAcc[p][gc2] += eventWeight;
	yieldHLepsAreInEtaAcc[p][gc3] += eventWeight;
	Float_t SortedGenHLepPt[4];
	for(int sl=0; sl<4; sl++) SortedGenHLepPt[sl] = GenHLepPt[sl];
	TriBulle(SortedGenHLepPt,4,false);
	for(int sl=0; sl<4; sl++){
	  h1GenptHLepsAreInEtaAcc[p][gc1][sl]->Fill(SortedGenHLepPt[sl],eventWeight);
	  h1GenptHLepsAreInEtaAcc[p][gc2][sl]->Fill(SortedGenHLepPt[sl],eventWeight);
	  h1GenptHLepsAreInEtaAcc[p][gc3][sl]->Fill(SortedGenHLepPt[sl],eventWeight);
	}
      }
      if(nGenHLepInPtAcc==4){
	Float_t SortedGenHLepAbseta[4];
	for(int sl=0; sl<4; sl++) SortedGenHLepAbseta[sl] = fabs(GenHLepEta[sl]);
	TriBulle(SortedGenHLepAbseta,4,false);
	for(int sl=0; sl<4; sl++){
	  h1GenetaHLepsAreInPtAcc[p][gc1][sl]->Fill(SortedGenHLepAbseta[sl],eventWeight);
	  h1GenetaHLepsAreInPtAcc[p][gc2][sl]->Fill(SortedGenHLepAbseta[sl],eventWeight);
	  h1GenetaHLepsAreInPtAcc[p][gc3][sl]->Fill(SortedGenHLepAbseta[sl],eventWeight);
	}
      }
      if(nGenHLepInEtaPtAcc==4){
	nbHLepsAreInEtaPtAcc[p][gc1]++;
	nbHLepsAreInEtaPtAcc[p][gc2]++;
	nbHLepsAreInEtaPtAcc[p][gc3]++;
	yieldHLepsAreInEtaPtAcc[p][gc1] += eventWeight;
	yieldHLepsAreInEtaPtAcc[p][gc2] += eventWeight;
	yieldHLepsAreInEtaPtAcc[p][gc3] += eventWeight;
	h1NRecoLepHLepsAreInEtaPtAcc[p][gc1]->Fill(NRecoMu+NRecoEle,eventWeight);
	h1NRecoLepHLepsAreInEtaPtAcc[p][gc2]->Fill(NRecoMu+NRecoEle,eventWeight);
	h1NRecoLepHLepsAreInEtaPtAcc[p][gc3]->Fill(NRecoMu+NRecoEle,eventWeight);
	h1NRecoMuoHLepsAreInEtaPtAcc[p][gc1]->Fill(NRecoMu,eventWeight);
	h1NRecoMuoHLepsAreInEtaPtAcc[p][gc2]->Fill(NRecoMu,eventWeight);
	h1NRecoMuoHLepsAreInEtaPtAcc[p][gc3]->Fill(NRecoMu,eventWeight);
	h1NRecoEleHLepsAreInEtaPtAcc[p][gc1]->Fill(NRecoEle,eventWeight);
	h1NRecoEleHLepsAreInEtaPtAcc[p][gc2]->Fill(NRecoEle,eventWeight);
	h1NRecoEleHLepsAreInEtaPtAcc[p][gc3]->Fill(NRecoEle,eventWeight);
      }
      if(nGenLep>=4){
	nb4GenLeps[p][gc1]++;
	nb4GenLeps[p][gc2]++;
	nb4GenLeps[p][gc3]++;
	yield4GenLeps[p][gc1] += eventWeight;
	yield4GenLeps[p][gc2] += eventWeight;
	yield4GenLeps[p][gc3] += eventWeight;
      }
      if(NRecoMu+NRecoEle>=4){
	nb4RecoLeps[p][gc1]++;
	nb4RecoLeps[p][gc2]++;
	nb4RecoLeps[p][gc3]++;
	yield4RecoLeps[p][gc1] += eventWeight;
	yield4RecoLeps[p][gc2] += eventWeight;
	yield4RecoLeps[p][gc3] += eventWeight;
      }
      if(FullSel100){
	nbWithBCFullSel100G[p][gc1]++;
	nbWithBCFullSel100G[p][gc2]++;
	nbWithBCFullSel100G[p][gc3]++;
	yieldWithBCFullSel100G[p][gc1] += eventWeight;
	yieldWithBCFullSel100G[p][gc2] += eventWeight;
	yieldWithBCFullSel100G[p][gc3] += eventWeight;
      }
      if(passTrigger){
	nbPassTriggerG[p][gc1]++;
	nbPassTriggerG[p][gc2]++;
	nbPassTriggerG[p][gc3]++;
	yieldPassTriggerG[p][gc1] += eventWeight;
	yieldPassTriggerG[p][gc2] += eventWeight;
	yieldPassTriggerG[p][gc3] += eventWeight;
      }
      if(passTriggerNo1E){
	nbPassTriggerGNo1E[p][gc1]++;
	nbPassTriggerGNo1E[p][gc2]++;
	nbPassTriggerGNo1E[p][gc3]++;
	yieldPassTriggerGNo1E[p][gc1] += eventWeight;
	yieldPassTriggerGNo1E[p][gc2] += eventWeight;
	yieldPassTriggerGNo1E[p][gc3] += eventWeight;
      }
      if(FullSel100 && passTrigger){
	nbPassTriggerGWithBCFullSel100G[p][gc1]++;
	nbPassTriggerGWithBCFullSel100G[p][gc2]++;
	nbPassTriggerGWithBCFullSel100G[p][gc3]++;
	yieldPassTriggerGWithBCFullSel100G[p][gc1] += eventWeight;
	yieldPassTriggerGWithBCFullSel100G[p][gc2] += eventWeight;
	yieldPassTriggerGWithBCFullSel100G[p][gc3] += eventWeight;
	Float_t SortedCandLepPt[4];
	for(int sl=0; sl<4; sl++) SortedCandLepPt[sl] = CandLepPt->at(sl);
	TriBulle(SortedCandLepPt,4,false);
	for(int sl=0; sl<4; sl++){
	  h1RecoptWithBCFullSel100GPassTriggerG[p][gc1][sl]->Fill(SortedCandLepPt[sl],eventWeight);
	  h1RecoptWithBCFullSel100GPassTriggerG[p][gc2][sl]->Fill(SortedCandLepPt[sl],eventWeight);
	  h1RecoptWithBCFullSel100GPassTriggerG[p][gc3][sl]->Fill(SortedCandLepPt[sl],eventWeight);
	}
	Float_t SortedCandLepAbseta[4];
	for(int sl=0; sl<4; sl++) SortedCandLepAbseta[sl] = fabs(CandLepEta->at(sl));
	TriBulle(SortedCandLepAbseta,4,false);
	for(int sl=0; sl<4; sl++){
	  h1RecoetaWithBCFullSel100GPassTriggerG[p][gc1][sl]->Fill(SortedCandLepAbseta[sl],eventWeight);
	  h1RecoetaWithBCFullSel100GPassTriggerG[p][gc2][sl]->Fill(SortedCandLepAbseta[sl],eventWeight);
	  h1RecoetaWithBCFullSel100GPassTriggerG[p][gc3][sl]->Fill(SortedCandLepAbseta[sl],eventWeight);
	}
      }
      if(FullSel100 && passTriggerNo1E){
	nbPassTriggerGNo1EWithBCFullSel100G[p][gc1]++;
	nbPassTriggerGNo1EWithBCFullSel100G[p][gc2]++;
	nbPassTriggerGNo1EWithBCFullSel100G[p][gc3]++;
	yieldPassTriggerGNo1EWithBCFullSel100G[p][gc1] += eventWeight;
	yieldPassTriggerGNo1EWithBCFullSel100G[p][gc2] += eventWeight;
	yieldPassTriggerGNo1EWithBCFullSel100G[p][gc3] += eventWeight;
      }
      if(nGenHLepInEtaPtAcc==4 && NRecoMu+NRecoEle>=4){
	nbHLepsAreInEtaPtAcc4RecoLeps[p][gc1]++;
	nbHLepsAreInEtaPtAcc4RecoLeps[p][gc2]++;
	nbHLepsAreInEtaPtAcc4RecoLeps[p][gc3]++;
	yieldHLepsAreInEtaPtAcc4RecoLeps[p][gc1] += eventWeight;
	yieldHLepsAreInEtaPtAcc4RecoLeps[p][gc2] += eventWeight;
	yieldHLepsAreInEtaPtAcc4RecoLeps[p][gc3] += eventWeight;
      }
      if(nGenHLepInEtaPtAcc==4 && FullSel100){
	nbHLepsAreInEtaPtAccWithBCFullSel100G[p][gc1]++;
	nbHLepsAreInEtaPtAccWithBCFullSel100G[p][gc2]++;
	nbHLepsAreInEtaPtAccWithBCFullSel100G[p][gc3]++;
	yieldHLepsAreInEtaPtAccWithBCFullSel100G[p][gc1] += eventWeight;
	yieldHLepsAreInEtaPtAccWithBCFullSel100G[p][gc2] += eventWeight;
	yieldHLepsAreInEtaPtAccWithBCFullSel100G[p][gc3] += eventWeight;
      }
      if(nGenHLepInEtaPtAcc==4 && passTrigger){
	nbHLepsAreInEtaPtAccPassTriggerG[p][gc1]++;
	nbHLepsAreInEtaPtAccPassTriggerG[p][gc2]++;
	nbHLepsAreInEtaPtAccPassTriggerG[p][gc3]++;
	yieldHLepsAreInEtaPtAccPassTriggerG[p][gc1] += eventWeight;
	yieldHLepsAreInEtaPtAccPassTriggerG[p][gc2] += eventWeight;
	yieldHLepsAreInEtaPtAccPassTriggerG[p][gc3] += eventWeight;
      }
      if(nGenHLepInEtaPtAcc==4 && passTriggerNo1E){
	nbHLepsAreInEtaPtAccPassTriggerGNo1E[p][gc1]++;
	nbHLepsAreInEtaPtAccPassTriggerGNo1E[p][gc2]++;
	nbHLepsAreInEtaPtAccPassTriggerGNo1E[p][gc3]++;
	yieldHLepsAreInEtaPtAccPassTriggerGNo1E[p][gc1] += eventWeight;
	yieldHLepsAreInEtaPtAccPassTriggerGNo1E[p][gc2] += eventWeight;
	yieldHLepsAreInEtaPtAccPassTriggerGNo1E[p][gc3] += eventWeight;
      }
      if(nGenHLepInEtaPtAcc==4 && FullSel100 && passTrigger){
	nbHLepsAreInEtaPtAccPassTriggerGWithBCFullSel100G[p][gc1]++;
	nbHLepsAreInEtaPtAccPassTriggerGWithBCFullSel100G[p][gc2]++;
	nbHLepsAreInEtaPtAccPassTriggerGWithBCFullSel100G[p][gc3]++;
	yieldHLepsAreInEtaPtAccPassTriggerGWithBCFullSel100G[p][gc1] += eventWeight;
	yieldHLepsAreInEtaPtAccPassTriggerGWithBCFullSel100G[p][gc2] += eventWeight;
	yieldHLepsAreInEtaPtAccPassTriggerGWithBCFullSel100G[p][gc3] += eventWeight;
      }
      if(nGenHLepInEtaPtAcc==4 && FullSel100 && passTriggerNo1E){
	nbHLepsAreInEtaPtAccPassTriggerGNo1EWithBCFullSel100G[p][gc1]++;
	nbHLepsAreInEtaPtAccPassTriggerGNo1EWithBCFullSel100G[p][gc2]++;
	nbHLepsAreInEtaPtAccPassTriggerGNo1EWithBCFullSel100G[p][gc3]++;
	yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCFullSel100G[p][gc1] += eventWeight;
	yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCFullSel100G[p][gc2] += eventWeight;
	yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCFullSel100G[p][gc3] += eventWeight;
      }
      h1GenHPt[p]->Fill(GenHPt,eventWeight);
      h1GenHRapidity[p]->Fill(GenHRapidity,eventWeight);
      h2GenHRapidityVsPt[p]->Fill(GenHPt,GenHRapidity,eventWeight);
      //       nbStored[p][rc]++;
      //       nbStored[p][rc2]++;
      //       nbStored[p][rc3]++;
      //       yieldStored[p][rc] += eventWeight;
      //       yieldStored[p][rc2] += eventWeight;
      //       yieldStored[p][rc3] += eventWeight;

      //if(prodName[p]=="ttH" && FullSel100){
      if(prodName[p]=="ttH"){
	nbCheckTTH[0]++; yieldCheckTTH[0]+=eventWeight;
	if      (nGenHLep==4){ //H->4l
	  nbCheckTTH[1]++; yieldCheckTTH[1]+=eventWeight;
	}else if(nGenHLep==2 && nGenHTau==2){ //H->2l2t
	  nbCheckTTH[2]++; yieldCheckTTH[2]+=eventWeight;
	}else if(nGenHLep==0 && nGenHTau==4){ //H->4t
	  nbCheckTTH[3]++; yieldCheckTTH[3]+=eventWeight;
	}else if(nGenHLep==2 && nGenHTau==0){ //H->2l2X 
	  nbCheckTTH[4]++; yieldCheckTTH[4]+=eventWeight;
	}else if(nGenHLep==0 && nGenHTau==2){ //H->2t2X
	  nbCheckTTH[5]++; yieldCheckTTH[5]+=eventWeight;
	}else{ //trash
	  cout<<"Error in ttH check"<<endl;
	}
	if      (nGenLEP==4){ //ttH->4lX
	  nbCheckTTH[6]++; yieldCheckTTH[6]+=eventWeight;
	}else if(nGenLEP==5){ //ttH->5lX
	  nbCheckTTH[7]++; yieldCheckTTH[7]+=eventWeight;
	}else if(nGenLEP==6){ //ttH->6lX
	  nbCheckTTH[8]++; yieldCheckTTH[8]+=eventWeight;
	}else{ //trash
	  nbCheckTTH[9]++; yieldCheckTTH[9]+=eventWeight;
	}
	if(nGenLEP==4){
	  if      (nGenEle==0){
	    nbCheckTTH[10]++; yieldCheckTTH[10]+=eventWeight;
	  }else if(nGenEle==1){
	    nbCheckTTH[11]++; yieldCheckTTH[11]+=eventWeight;
	  }else if(nGenEle==2){
	    nbCheckTTH[12]++; yieldCheckTTH[12]+=eventWeight;
	  }else if(nGenEle==3){
	    nbCheckTTH[13]++; yieldCheckTTH[13]+=eventWeight;
	  }else if(nGenEle==4){
	    nbCheckTTH[14]++; yieldCheckTTH[14]+=eventWeight;
	  }else{ //trash
	    nbCheckTTH[15]++; yieldCheckTTH[15]+=eventWeight;
	  }
	  if      (nGenMu ==0){
	    nbCheckTTH[16]++; yieldCheckTTH[16]+=eventWeight;
	  }else if(nGenMu ==1){
	    nbCheckTTH[17]++; yieldCheckTTH[17]+=eventWeight;
	  }else if(nGenMu ==2){
	    nbCheckTTH[18]++; yieldCheckTTH[18]+=eventWeight;
	  }else if(nGenMu ==3){
	    nbCheckTTH[19]++; yieldCheckTTH[19]+=eventWeight;
	  }else if(nGenMu ==4){
	    nbCheckTTH[20]++; yieldCheckTTH[20]+=eventWeight;
	  }else{ //trash
	    nbCheckTTH[21]++; yieldCheckTTH[21]+=eventWeight;
	  }
	  if      (nGenTau==0){
	    nbCheckTTH[22]++; yieldCheckTTH[22]+=eventWeight;
	  }else if(nGenTau==1){
	    nbCheckTTH[23]++; yieldCheckTTH[23]+=eventWeight;
	  }else if(nGenTau==2){
	    nbCheckTTH[24]++; yieldCheckTTH[24]+=eventWeight;
	  }else if(nGenTau==3){
	    nbCheckTTH[25]++; yieldCheckTTH[25]+=eventWeight;
	  }else if(nGenTau==4){
	    nbCheckTTH[26]++; yieldCheckTTH[26]+=eventWeight;
	  }else{ //trash
	    nbCheckTTH[27]++; yieldCheckTTH[27]+=eventWeight;
	  }
	  if(nGenLEPPlus==2 && nGenLEPMinus==2){
	    nbCheckTTH[28]++; yieldCheckTTH[28]+=eventWeight;
	  }else{ //trash
	    nbCheckTTH[29]++; yieldCheckTTH[29]+=eventWeight;
	  }
	  if      (nGenEle==4 && nGenMu ==0 && nGenTau==0){
	    nbCheckTTH[30]++; yieldCheckTTH[30]+=eventWeight;
	  }else if(nGenEle==0 && nGenMu ==4 && nGenTau==0){
	    nbCheckTTH[31]++; yieldCheckTTH[31]+=eventWeight;
	  }else if(nGenEle==0 && nGenMu ==0 && nGenTau==4){
	    nbCheckTTH[32]++; yieldCheckTTH[32]+=eventWeight;
	  }else if(nGenEle==2 && nGenMu ==2 && nGenTau==0){
	    nbCheckTTH[33]++; yieldCheckTTH[33]+=eventWeight;
	  }else if(nGenEle==2 && nGenMu ==0 && nGenTau==2){
	    nbCheckTTH[34]++; yieldCheckTTH[34]+=eventWeight;
	  }else if(nGenEle==0 && nGenMu ==2 && nGenTau==2){
	    nbCheckTTH[35]++; yieldCheckTTH[35]+=eventWeight;
	  }else{
	    nbCheckTTH[36]++; yieldCheckTTH[36]+=eventWeight;
	  }
	}
	//if(1){
	//if(nGenHLEP<4){
	if(nGenHLep<4){
	  if      (nGenHEle==4 && nGenHMu ==0 && nGenHTau==0){
	    nbCheckTTH[37]++; yieldCheckTTH[37]+=eventWeight;
	  }else if(nGenHEle==0 && nGenHMu ==4 && nGenHTau==0){
	    nbCheckTTH[38]++; yieldCheckTTH[38]+=eventWeight;
	  }else if(nGenHEle==0 && nGenHMu ==0 && nGenHTau==4){
	    nbCheckTTH[39]++; yieldCheckTTH[39]+=eventWeight;
	  }else if(nGenHEle==2 && nGenHMu ==2 && nGenHTau==0){
	    nbCheckTTH[40]++; yieldCheckTTH[40]+=eventWeight;
	  }else if(nGenHEle==2 && nGenHMu ==0 && nGenHTau==2){
	    nbCheckTTH[41]++; yieldCheckTTH[41]+=eventWeight;
	  }else if(nGenHEle==0 && nGenHMu ==2 && nGenHTau==2){
	    nbCheckTTH[42]++; yieldCheckTTH[42]+=eventWeight;
	  }else if(nGenHEle==2 && nGenHMu ==0 && nGenHTau==0){
	    nbCheckTTH[43]++; yieldCheckTTH[43]+=eventWeight;
	  }else if(nGenHEle==0 && nGenHMu ==2 && nGenHTau==0){
	    nbCheckTTH[44]++; yieldCheckTTH[44]+=eventWeight;
	  }else if(nGenHEle==0 && nGenHMu ==0 && nGenHTau==2){
	    nbCheckTTH[45]++; yieldCheckTTH[45]+=eventWeight;
	  }else{
	    nbCheckTTH[46]++; yieldCheckTTH[46]+=eventWeight;
	  }
	  if      (nGenAssocEle==2 && nGenAssocMu ==0 && nGenAssocTau==0){
	    nbCheckTTH[47]++; yieldCheckTTH[47]+=eventWeight;
	  }else if(nGenAssocEle==0 && nGenAssocMu ==2 && nGenAssocTau==0){
	    nbCheckTTH[48]++; yieldCheckTTH[48]+=eventWeight;
	  }else if(nGenAssocEle==0 && nGenAssocMu ==0 && nGenAssocTau==2){
	    nbCheckTTH[49]++; yieldCheckTTH[49]+=eventWeight;
	  }else if(nGenAssocEle==1 && nGenAssocMu ==1 && nGenAssocTau==0){
	    nbCheckTTH[50]++; yieldCheckTTH[50]+=eventWeight;
	  }else if(nGenAssocEle==1 && nGenAssocMu ==0 && nGenAssocTau==1){
	    nbCheckTTH[51]++; yieldCheckTTH[51]+=eventWeight;
	  }else if(nGenAssocEle==0 && nGenAssocMu ==1 && nGenAssocTau==1){
	    nbCheckTTH[52]++; yieldCheckTTH[52]+=eventWeight;
	  }else if(nGenAssocEle==1 && nGenAssocMu ==0 && nGenAssocTau==0){
	    nbCheckTTH[53]++; yieldCheckTTH[53]+=eventWeight;
	  }else if(nGenAssocEle==0 && nGenAssocMu ==1 && nGenAssocTau==0){
	    nbCheckTTH[54]++; yieldCheckTTH[54]+=eventWeight;
	  }else if(nGenAssocEle==0 && nGenAssocMu ==0 && nGenAssocTau==1){
	    nbCheckTTH[55]++; yieldCheckTTH[55]+=eventWeight;
	  }else if(nGenAssocEle==0 && nGenAssocMu ==0 && nGenAssocTau==0){
	    nbCheckTTH[56]++; yieldCheckTTH[56]+=eventWeight;
	  }else{
	    nbCheckTTH[57]++; yieldCheckTTH[57]+=eventWeight;
	  }
	}
      }


      // ---------------------- Cuts related to gen leptons --------------------------

      // Don't want to consider such events.
      // (They are sometimes selected because some e/mu from tau decays are reconstructed as good leptons.)
      if(nGenLep<4 && p<nHiggsSamples) continue; 

      if(
	 (p==7 && !(nGenHLep==4&&nGenAssocLep==0)) ||
	 (p<nHiggsSamples && !( nGenHLep==4 || (nGenHLep==2&&nGenAssocLep==2&&(prodName[p]=="ZH"||prodName[p]=="ttH")) ) )
	 ){
	cout<<"ERROR with nb. of gen leptons (sample "<<prodName[p]<<", event "<<z<<", nGenHLep="<<nGenHLep<<", nGenAssocLep="<<nGenAssocLep<<")"<<endl;
	continue;
      }

      ///////////////////////////
      if(excludeH2l2X)
	if(nGenHLep!=4 && p<nHiggsSamples) continue;
      ///////////////////////////


      // ---------------------- Gen associated decay --------------------------

      Int_t currentAssocDecay = -1;
      if(prodName[p]=="WH"||prodName[p]=="WplusH"||prodName[p]=="WminusH"){
	if(nGenHLep==4 && nGenAssocLep==0) currentAssocDecay = 0;
	else if(nGenHLep==4 && nGenAssocLep==1) currentAssocDecay = 1;
      }else if(prodName[p]=="ZH"){
	if(nGenHLep==4 && nGenAssocLep==0) currentAssocDecay = 2;
	else if(nGenHLep==4 && nGenAssocLep==2) currentAssocDecay = 3;
	else if(nGenHLep==2 && nGenAssocLep==2) currentAssocDecay = 4;
      }else if(prodName[p]=="ttH"){
	if(nGenHLep==4 && nGenAssocLep==0) currentAssocDecay = 5;
	else if(nGenHLep==4 && nGenAssocLep==1) currentAssocDecay = 6;
	else if(nGenHLep==4 && nGenAssocLep==2) currentAssocDecay = 7;
	else if(nGenHLep==2 && nGenAssocLep==2) currentAssocDecay = 8;
      }
      if((prodName[p]=="WH"||prodName[p]=="WplusH"||prodName[p]=="WminusH" || prodName[p]=="ZH" || prodName[p]=="ttH") && currentAssocDecay == -1)
	cout<<"ERROR with assoc decays: process="<<prodName[p]<<", nGenHLep="<<nGenHLep<<", nGenAssocLep="<<nGenAssocLep<<endl;


      // ---------------------- Successive selection steps --------------------------

      nbWithBC[p][rc]++;
      nbWithBC[p][rc2]++;
      yieldWithBC[p][rc] += eventWeight;
      yieldWithBC[p][rc2] += eventWeight;
	  
      if(FullSel70){
	nbWithBCFullSel70[p][rc]++;
	nbWithBCFullSel70[p][rc2]++;
	yieldWithBCFullSel70[p][rc] += eventWeight;
	yieldWithBCFullSel70[p][rc2] += eventWeight;
      }
	  
      if(FullSel100){

	nbWithBCFullSel100[p][rc]++;
	nbWithBCFullSel100[p][rc2]++;
	yieldWithBCFullSel100[p][rc] += eventWeight;
	yieldWithBCFullSel100[p][rc2] += eventWeight;


	// ---------------------- Gen to Good lepton matching --------------------------

	Int_t nRecoLepMatchedToGenHLep[4] = {0,0,0,0};
	Int_t nCandLepMatchedToGenHLep[4] = {0,0,0,0};
	Int_t nRecoLepMatchedToGenAssocLep[2] = {0,0};
	Int_t nCandLepMatchedToGenAssocLep[2] = {0,0};
	Int_t nGenLepMatchedToCandLep[4] = {0,0,0,0};
	Int_t nGenLepMatchedToRecoLep[7] = {0,0,0,0,0,0,0};
	Int_t nGenHLepMatchedToZ1Lep[4] = {0,0};
	Int_t nGenHLepMatchedToZ2Lep[4] = {0,0};
	for(Int_t iGenHLep=0; iGenHLep<4; iGenHLep++){
	  if(abs(GenHLepId[iGenHLep])==11 || abs(GenHLepId[iGenHLep])==13){
	    for(Int_t iCandLep=0; iCandLep<4; iCandLep++){
	      if(deltaR(GenHLepEta[iGenHLep],GenHLepPhi[iGenHLep],CandLepEta->at(iCandLep),CandLepPhi->at(iCandLep)) < 0.1){
		nRecoLepMatchedToGenHLep[iGenHLep]++;
		nCandLepMatchedToGenHLep[iGenHLep]++;
		nGenLepMatchedToRecoLep[iCandLep]++;
		nGenLepMatchedToCandLep[iCandLep]++;
		if(iCandLep<2){
		  nGenHLepMatchedToZ1Lep[iCandLep]++;
		}else{
		  nGenHLepMatchedToZ2Lep[iCandLep-2]++;
		}
	      }
	    }
	    for(Int_t iExtraLep=0; iExtraLep<nExtraLep; iExtraLep++){
	      if(deltaR(GenHLepEta[iGenHLep],GenHLepPhi[iGenHLep],ExtraLepEta->at(iExtraLep),ExtraLepPhi->at(iExtraLep)) < 0.1){
		nRecoLepMatchedToGenHLep[iGenHLep]++;
		nGenLepMatchedToRecoLep[4+iExtraLep]++;
	      }
	    }
	  }
	}
	for(Int_t iGenAssocLep=0; iGenAssocLep<2; iGenAssocLep++){
	  if(abs(GenAssocLepId[iGenAssocLep])==11 || abs(GenAssocLepId[iGenAssocLep])==13){
	    for(Int_t iCandLep=0; iCandLep<4; iCandLep++){
	      if(deltaR(GenAssocLepEta[iGenAssocLep],GenAssocLepPhi[iGenAssocLep],CandLepEta->at(iCandLep),CandLepPhi->at(iCandLep)) < 0.1){
		nRecoLepMatchedToGenAssocLep[iGenAssocLep]++;
		nCandLepMatchedToGenAssocLep[iGenAssocLep]++;
		nGenLepMatchedToRecoLep[iCandLep]++;
		nGenLepMatchedToCandLep[iCandLep]++;
	      }
	    }
	    for(Int_t iExtraLep=0; iExtraLep<nExtraLep; iExtraLep++){
	      if(deltaR(GenAssocLepEta[iGenAssocLep],GenAssocLepPhi[iGenAssocLep],ExtraLepEta->at(iExtraLep),ExtraLepPhi->at(iExtraLep)) < 0.1){
		nRecoLepMatchedToGenAssocLep[iGenAssocLep]++;
		nGenLepMatchedToRecoLep[4+iExtraLep]++;
	      }
	    }
	  }
	}


	// ---------------------- Exploit matching information --------------------------
	    
	Bool_t foundMatchingAmbiguity = false;
	for(Int_t iGenHLep=0; iGenHLep<4; iGenHLep++) if(nRecoLepMatchedToGenHLep[iGenHLep]>1){ foundMatchingAmbiguity = true; break; }
	for(Int_t iGenAssocLep=0; iGenAssocLep<2; iGenAssocLep++) if(nRecoLepMatchedToGenAssocLep[iGenAssocLep]>1){ foundMatchingAmbiguity = true; break; }
	for(Int_t iRecoLep=0; iRecoLep<4+nExtraLep; iRecoLep++) if(nGenLepMatchedToRecoLep[iRecoLep]>1){ foundMatchingAmbiguity = true; break; }

	Int_t nOnes = 0;
	for(Int_t iGenHLep=0; iGenHLep<4; iGenHLep++) if(nRecoLepMatchedToGenHLep[iGenHLep]==1) nOnes++;

	Int_t nOnesHLeps = 0;
	for(Int_t iGenHLep=0; iGenHLep<4; iGenHLep++) if(nCandLepMatchedToGenHLep[iGenHLep]==1) nOnesHLeps++;
	Int_t nOnesAssocLeps = 0;
	for(Int_t iGenAssocLep=0; iGenAssocLep<2; iGenAssocLep++) if(nCandLepMatchedToGenAssocLep[iGenAssocLep]==1) nOnesAssocLeps++;
	Int_t currentMatchHLepsStatus = -1;
	if(foundMatchingAmbiguity){
	  currentMatchHLepsStatus = 5;
	}else{
	  if(nOnesHLeps==4) currentMatchHLepsStatus = 0;
	  if(nOnesHLeps==3) currentMatchHLepsStatus = 1;
	  if(nOnesHLeps==2) currentMatchHLepsStatus = 2;
	  if(nOnesHLeps==1) currentMatchHLepsStatus = 3;
	  if(nOnesHLeps==0) currentMatchHLepsStatus = 4;
	}
	Int_t currentMatchAllLepsStatus = -1;
	if(foundMatchingAmbiguity){
	  currentMatchAllLepsStatus = 4;
	}else{
	  if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchAllLepsStatus = 0;
	  if(nOnesHLeps==3 && nOnesAssocLeps==1) currentMatchAllLepsStatus = 1;
	  if(nOnesHLeps==2 && nOnesAssocLeps==2) currentMatchAllLepsStatus = 2;
	  if(nOnesHLeps+nOnesAssocLeps<4) currentMatchAllLepsStatus = 3;
	}
	Int_t currentMatchWHStatus = -1;
	Int_t currentMatchZHStatus = -1;
	Int_t currentMatchttHStatus = -1;
	if(prodName[p]=="WH"||prodName[p]=="WplusH"||prodName[p]=="WminusH"){
	  if(foundMatchingAmbiguity){
	    currentMatchWHStatus = 4;
	  }else{
	    if(nOnesHLeps+nOnesAssocLeps<4){
	      currentMatchWHStatus = 3;
	    }else{
	      if(nGenHLep==4 && nGenAssocLep==0){
		if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchWHStatus = 0;
		else cout<<"error nOnes"<<endl;
	      }else if(nGenHLep==4 && nGenAssocLep==1){
		if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchWHStatus = 1;
		else if(nOnesHLeps==3 && nOnesAssocLeps==1) currentMatchWHStatus = 2;
		else cout<<"error nOnes"<<endl;
	      }else{
		cout<<"error nGen"<<endl;
	      }
	    }
	  } 
	}
	if(prodName[p]=="ZH"){
	  if(foundMatchingAmbiguity){
	    currentMatchZHStatus = 6;
	  }else{
	    if(nOnesHLeps+nOnesAssocLeps<4){
	      currentMatchZHStatus = 5;
	    }else{
	      if(nGenHLep==4 && nGenAssocLep==0){
		if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchZHStatus = 0;
		else cout<<"error nOnes"<<endl;
	      }else if(nGenHLep==4 && nGenAssocLep==2){
		if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchZHStatus = 1;
		else if(nOnesHLeps==3 && nOnesAssocLeps==1) currentMatchZHStatus = 2;
		else if(nOnesHLeps==2 && nOnesAssocLeps==2) currentMatchZHStatus = 3;
		else cout<<"error nOnes"<<endl;
	      }else if(nGenHLep==2 && nGenAssocLep==2){
		if(nOnesHLeps==2 && nOnesAssocLeps==2) currentMatchZHStatus = 4;
		else cout<<"error nOnes"<<endl;
	      }else{
		cout<<"error nGen"<<endl;
	      }
	    }
	  } 
	}
	if(prodName[p]=="ttH"){
	  if(foundMatchingAmbiguity){
	    currentMatchttHStatus = 8;
	  }else{
	    if(nOnesHLeps+nOnesAssocLeps<4){
	      currentMatchttHStatus = 7;
	    }else{
	      if(nGenHLep==4 && nGenAssocLep==0){
		if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchttHStatus = 0;
		else cout<<"error nOnes"<<endl;
	      }else if(nGenHLep==4 && nGenAssocLep==1){
		if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchttHStatus = 1;
		else if(nOnesHLeps==3 && nOnesAssocLeps==1) currentMatchttHStatus = 2;
		else cout<<"error nOnes"<<endl;
	      }else if(nGenHLep==4 && nGenAssocLep==2){
		if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchttHStatus = 3;
		else if(nOnesHLeps==3 && nOnesAssocLeps==1) currentMatchttHStatus = 4;
		else if(nOnesHLeps==2 && nOnesAssocLeps==2) currentMatchttHStatus = 5;
		else cout<<"error nOnes"<<endl;
	      }else if(nGenHLep==2 && nGenAssocLep==2){
		if(nOnesHLeps==2 && nOnesAssocLeps==2) currentMatchttHStatus = 6;
		else cout<<"error nOnes"<<endl;
	      }else{
		cout<<"error nGen"<<endl;
	      }
	    }
	  } 
	} 

	Int_t currentZ1MatchStatus = -1;
	if(nGenHLepMatchedToZ1Lep[0]>1 || nGenHLepMatchedToZ1Lep[1]>1) currentZ1MatchStatus = 3;
	else currentZ1MatchStatus = 2 - (nGenHLepMatchedToZ1Lep[0] + nGenHLepMatchedToZ1Lep[1]);
	Int_t currentZ2MatchStatus = -1;
	if(nGenHLepMatchedToZ2Lep[0]>1 || nGenHLepMatchedToZ2Lep[1]>1) currentZ2MatchStatus = 3;
	else currentZ2MatchStatus = 2 - (nGenHLepMatchedToZ2Lep[0] + nGenHLepMatchedToZ2Lep[1]);


	// ---------------------- Cuts impacting counters/histograms --------------------------

	Bool_t Exactly4GoodLeps = nExtraLep==0 ;
	Bool_t AtLeast5GoodLeps = nExtraLep>=1 ;
	Bool_t Exactly5GoodLeps = nExtraLep==1 ;
	Bool_t Exactly6GoodLeps = nExtraLep==2 ;
	Bool_t HLepsAreInEtaPtAcc    = nGenHLepInEtaPtAcc==4 ;
	Bool_t HLepsAreGood     = !foundMatchingAmbiguity && nOnes==4 ;

	if(requireExactly4GoodLeps && !Exactly4GoodLeps) continue;
	if(requireAtLeast5GoodLeps && !AtLeast5GoodLeps) continue;
	if(requireExactly5GoodLeps && !Exactly5GoodLeps) continue;
	if(requireExactly6GoodLeps && !Exactly6GoodLeps) continue;
	if(requireHLepsAreInEtaPtAcc    && !HLepsAreInEtaPtAcc   ) continue;
	if(requireHLepsAreGood     && !HLepsAreGood    ) continue;


	// ---------------------- Fill histograms and increment counters --------------------------

	if(Exactly4GoodLeps){ nbWithBCFullSel100Exactly4GoodLeps[p][rc]++; nbWithBCFullSel100Exactly4GoodLeps[p][rc2]++; }
	if(AtLeast5GoodLeps){ nbWithBCFullSel100AtLeast5GoodLeps[p][rc]++; nbWithBCFullSel100AtLeast5GoodLeps[p][rc2]++; }
	if(Exactly5GoodLeps){ nbWithBCFullSel100Exactly5GoodLeps[p][rc]++; nbWithBCFullSel100Exactly5GoodLeps[p][rc2]++; }
	if(Exactly6GoodLeps){ nbWithBCFullSel100Exactly6GoodLeps[p][rc]++; nbWithBCFullSel100Exactly6GoodLeps[p][rc2]++; }
	if(HLepsAreInEtaPtAcc   ){ nbWithBCFullSel100HLepsAreInEtaPtAcc   [p][rc]++; nbWithBCFullSel100HLepsAreInEtaPtAcc   [p][rc2]++; }
	if(HLepsAreGood    ){ nbWithBCFullSel100HLepsAreGood    [p][rc]++; nbWithBCFullSel100HLepsAreGood    [p][rc2]++; }

	Float_t varVal[nVariables] = {
	  ZZMass,
	  Z1Mass,
	  Z2Mass,
	  p0plus_VAJHU / ( p0plus_VAJHU + bkg_VAMCFM ),
	  DiJetFisher,
	  ZZPt,
	  (Float_t)nGenLep,
	  (Float_t)nGenLepInEtaPtAcc,
	  (Float_t)(nGenLep-nGenLepInEtaPtAcc),
	  (Float_t)(nGenHLep-nGenHLepInEtaPtAcc),
	  (Float_t)(nGenAssocLep-nGenAssocLepInEtaPtAcc),
	  (Float_t)(nGenLep-(4+nExtraLep)),
	  (Float_t)(nGenLepInEtaPtAcc-(4+nExtraLep)),
	  (Float_t)nExtraLep,
	  (Float_t)nExtraZ,
	  (Float_t)nJets,
	};

	for(int v=0; v<nVariables; v++){
	  if(!plotThisVar[VARIABLELIST][v]) continue;
	  hBCFullSel100[v][p][rc]->Fill(varVal[v],eventWeight);
	  hBCFullSel100[v][p][rc2]->Fill(varVal[v],eventWeight);
	}

	for(int v2=0; v2<n2DHist; v2++){
	  h2DBCFullSel100[v2][p][rc]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],eventWeight);
	  h2DBCFullSel100[v2][p][rc2]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],eventWeight);
	  h2DBCFullSel100Decays[v2][p][0][rc]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],eventWeight);
	  h2DBCFullSel100Decays[v2][p][0][rc2]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],eventWeight);
	  if(nGenHLep==4){
	    h2DBCFullSel100Decays[v2][p][1][rc]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],eventWeight);
	    h2DBCFullSel100Decays[v2][p][1][rc2]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],eventWeight);
	    if(currentMatchHLepsStatus==0){
	      h2DBCFullSel100Decays[v2][p][3][rc]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],eventWeight);
	      h2DBCFullSel100Decays[v2][p][3][rc2]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],eventWeight);
	    }else if(currentMatchHLepsStatus<5){
	      h2DBCFullSel100Decays[v2][p][4][rc]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],eventWeight);
	      h2DBCFullSel100Decays[v2][p][4][rc2]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],eventWeight);
	    }
	  }else{
	    h2DBCFullSel100Decays[v2][p][2][rc]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],eventWeight);
	    h2DBCFullSel100Decays[v2][p][2][rc2]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],eventWeight);
	  }
	}

	nbWithBCFullSel100MatchHLeps[currentMatchHLepsStatus][p][rc]++;
	nbWithBCFullSel100MatchHLeps[currentMatchHLepsStatus][p][rc2]++;
	for(int v=0; v<nVariables; v++){
	  if(!plotThisVar[VARIABLELIST][v]) continue;
	  hBCFullSel100MatchHLeps[v][currentMatchHLepsStatus][p][rc]->Fill(varVal[v],eventWeight);
	  hBCFullSel100MatchHLeps[v][currentMatchHLepsStatus][p][rc2]->Fill(varVal[v],eventWeight);
	}

	nbWithBCFullSel100MatchAllLeps[currentMatchAllLepsStatus][p][rc]++;
	nbWithBCFullSel100MatchAllLeps[currentMatchAllLepsStatus][p][rc2]++;
	for(int v=0; v<nVariables; v++){
	  if(!plotThisVar[VARIABLELIST][v]) continue;
	  hBCFullSel100MatchAllLeps[v][currentMatchAllLepsStatus][p][rc]->Fill(varVal[v],eventWeight);
	  hBCFullSel100MatchAllLeps[v][currentMatchAllLepsStatus][p][rc2]->Fill(varVal[v],eventWeight);
	}
	    
	if(prodName[p]=="WH"||prodName[p]=="WplusH"||prodName[p]=="WminusH"){
	  nbWithBCFullSel100MatchWH[currentMatchWHStatus][rc]++;
	  nbWithBCFullSel100MatchWH[currentMatchWHStatus][rc2]++;
	  for(int v=0; v<nVariables; v++){
	    if(!plotThisVar[VARIABLELIST][v]) continue;
	    hBCFullSel100MatchWH[v][currentMatchWHStatus][rc]->Fill(varVal[v],eventWeight);
	    hBCFullSel100MatchWH[v][currentMatchWHStatus][rc2]->Fill(varVal[v],eventWeight);
	  }
	}
	if(prodName[p]=="ZH"){
	  nbWithBCFullSel100MatchZH[currentMatchZHStatus][rc]++;
	  nbWithBCFullSel100MatchZH[currentMatchZHStatus][rc2]++;
	  for(int v=0; v<nVariables; v++){
	    if(!plotThisVar[VARIABLELIST][v]) continue;
	    hBCFullSel100MatchZH[v][currentMatchZHStatus][rc]->Fill(varVal[v],eventWeight);
	    hBCFullSel100MatchZH[v][currentMatchZHStatus][rc2]->Fill(varVal[v],eventWeight);
	  }
	}
	if(prodName[p]=="ttH"){
	  nbWithBCFullSel100MatchttH[currentMatchttHStatus][rc]++;
	  nbWithBCFullSel100MatchttH[currentMatchttHStatus][rc2]++;
	  for(int v=0; v<nVariables; v++){
	    if(!plotThisVar[VARIABLELIST][v]) continue;
	    hBCFullSel100MatchttH[v][currentMatchttHStatus][rc]->Fill(varVal[v],eventWeight);
	    hBCFullSel100MatchttH[v][currentMatchttHStatus][rc2]->Fill(varVal[v],eventWeight);
	  }
	}

	if(currentMatchAllLepsStatus==0){
	  nbWithBCFullSel100All4LepRight[p][rc]++;
	  nbWithBCFullSel100All4LepRight[p][rc2]++;
	}

	if(prodName[p]=="WH"||prodName[p]=="WplusH"||prodName[p]=="WminusH"){
	  Int_t currentWDecay = -1;
	  if(nGenHLep==4 && nGenAssocLep==0) currentWDecay = 0;
	  else if(nGenHLep==4 && nGenAssocLep==1) currentWDecay = 1;
	  else cout<<"error"<<endl;
	  nbTotalWH[currentWDecay]++;
	  if(nGenHLepInEtaPtAcc==4) nbHLepsAreInEtaPtAccWH[currentWDecay]++;
	  if(HLepsAreGood) nbHLepsAreGoodWH[currentWDecay]++;
	  if(currentMatchAllLepsStatus==0) nbAll4LepRightWH[currentWDecay]++;
	  nbZ1DaughtersFromHWH[currentWDecay][currentZ1MatchStatus]++;
	  nbZ2DaughtersFromHWH[currentWDecay][currentZ2MatchStatus]++;
	}
	if(prodName[p]=="ZH"){
	  Int_t currentZDecay = -1;
	  if(nGenHLep==4 && nGenAssocLep==0) currentZDecay = 0;
	  else if(nGenHLep==4 && nGenAssocLep==2) currentZDecay = 1;
	  else if(nGenHLep==2 && nGenAssocLep==2) currentZDecay = 2;
	  else cout<<"error"<<endl;
	  nbTotalZH[currentZDecay]++;
	  if(nGenHLepInEtaPtAcc==4) nbHLepsAreInEtaPtAccZH[currentZDecay]++;
	  if(HLepsAreGood) nbHLepsAreGoodZH[currentZDecay]++;
	  if(currentMatchAllLepsStatus==0) nbAll4LepRightZH[currentZDecay]++;
	  nbZ1DaughtersFromHZH[currentZDecay][currentZ1MatchStatus]++;
	  nbZ2DaughtersFromHZH[currentZDecay][currentZ2MatchStatus]++;
	}
	if(prodName[p]=="ttH"){
	  Int_t currentttDecay = -1;
	  if(nGenHLep==4 && nGenAssocLep==0) currentttDecay = 0;
	  else if(nGenHLep==4 && nGenAssocLep==1) currentttDecay = 1;
	  else if(nGenHLep==4 && nGenAssocLep==2) currentttDecay = 2;
	  else if(nGenHLep==2 && nGenAssocLep==2) currentttDecay = 3;
	  else cout<<"error"<<endl;
	  nbTotalttH[currentttDecay]++;
	  if(nGenHLepInEtaPtAcc==4) nbHLepsAreInEtaPtAccttH[currentttDecay]++;
	  if(HLepsAreGood) nbHLepsAreGoodttH[currentttDecay]++;
	  if(currentMatchAllLepsStatus==0) nbAll4LepRightttH[currentttDecay]++;
	  nbZ1DaughtersFromHttH[currentttDecay][currentZ1MatchStatus]++;
	  nbZ2DaughtersFromHttH[currentttDecay][currentZ2MatchStatus]++;
	}


	// ---------------------- Categorization stuff --------------------------

	Bool_t HadVHTag = 0;
	//* // method that is used in all plots shown before June 8th 2015
	if(nJets>=2){
	  if( 60.<DiJetMass && DiJetMass<120. && ZZPt>ZZMass ){
	    for(int j1=0; j1<nJets; j1++){
	      for(int j2=j1; j2<nJets; j2++){
		if( abs(JetEta->at(j1))<2.4 && abs(JetEta->at(j2))<2.4 
		    && JetPt->at(j1)>40. && JetPt->at(j2)>40. ){
		  HadVHTag = true;
		  break;
		}
	      }
	      if(HadVHTag) break;
	    }
	  }
	}
	//*/
	/* method used in first synchronization (but worse purity !)
	if(nJets>=2){
	  for(int j1=0; j1<nJets; j1++){
	    if( abs(JetEta->at(j1))<2.4 && JetPt->at(j1)>40. ){
	      for(int j2=j1+1; j2<nJets; j2++){
		if( abs(JetEta->at(j2))<2.4 && JetPt->at(j2)>40. ){
		  LV jet1 (JetPt->at(j1),JetEta->at(j1),JetPhi->at(j1),JetMass->at(j1));
		  LV jet2 (JetPt->at(j2),JetEta->at(j2),JetPhi->at(j2),JetMass->at(j2));
		  float mjj = (jet1+jet2).mass();
		  if( 60.<mjj && mjj<120. ){
		    HadVHTag = true;
		    break;
		  }
		}
	      }
	      if(HadVHTag) break;
	    }
	  }
	}
	//*/

	Int_t currentBasket = -1;

	if(BASKETLIST==0){

	  if(nJets==0 || nJets==1){
	    if(nExtraLep==0) currentBasket = 1;
	    else if(nExtraLep==1) currentBasket = 2;
	    else if(nExtraLep>=2) currentBasket = 3;
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else if(nJets>=2){
	    if(nExtraLep==0) currentBasket = 4;
	    else if(nExtraLep==1) currentBasket = 5;
	    else if(nExtraLep>=2) currentBasket = 6;
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else{
	    cout<<"WARNING : inconsistent nJets"<<endl;
	  }

	}else if(BASKETLIST==1){

	  if(nJets==0 || nJets==1){
	    if(nExtraLep==0) currentBasket = 1;
	    else if(nExtraLep==1){ 
	      //if(ExtraLep1Pt>20. && PFMET>45.) currentBasket = 2;
	      if(PFMET>45.) currentBasket = 2;
	      else currentBasket = 3;
	    }else if(nExtraLep>=2) currentBasket = 4;
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else if(nJets>=2){
	    if(nExtraLep==0){
	      if(DiJetFisher>0.5) currentBasket = 5;
	      else currentBasket = 6;	      
	    }else if(nExtraLep==1) currentBasket = 7;
	    else if(nExtraLep>=2) currentBasket = 8;
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else{
	    cout<<"WARNING : inconsistent nJets"<<endl;
	  }

	}else if(BASKETLIST==2){

	  if(nJets==0){
	    if(nExtraLep==0) currentBasket = 1;
	    else if(nExtraLep==1) currentBasket = 2;
	    else if(nExtraLep>=2) currentBasket = 3;
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else if(nJets==1){
	    if(nJetsBTagged==0){
	      //if(nExtraLep==0 && abs(JetEta->at(0))>2.4 && ZZPt>50.) currentBasket = 24; else
	      if(nExtraLep==0) currentBasket = 4;
	      else if(nExtraLep==1) currentBasket = 5;
	      else if(nExtraLep>=2) currentBasket = 6;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else if(nJetsBTagged==1){
	      if(nExtraLep==0) currentBasket = 7;
	      else if(nExtraLep==1) currentBasket = 8;
	      else if(nExtraLep>=2) currentBasket = 9;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else cout<<"WARNING : inconsistent nJetsBTagged"<<endl;
	  }else if(nJets>=2){
	    if(nJetsBTagged==0){
	      if(nExtraLep==0 && DiJetFisher>0.5) currentBasket = 10;
	      else if(nExtraLep==0 && HadVHTag) currentBasket = 11;
	      else if(nExtraLep==0) currentBasket = 12;
	      else if(nExtraLep==1) currentBasket = 13;
	      else if(nExtraLep>=2) currentBasket = 14;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else if(nJetsBTagged==1){
	      if(nExtraLep==0 && DiJetFisher>0.5) currentBasket = 15;
	      else if(nExtraLep==0 && HadVHTag) currentBasket = 16;
	      else if(nExtraLep==0) currentBasket = 17;
	      else if(nExtraLep==1) currentBasket = 18;
	      else if(nExtraLep>=2) currentBasket = 19;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else if(nJetsBTagged>=2){
	      if(nExtraLep==0 && HadVHTag) currentBasket = 20;
	      else if(nExtraLep==0) currentBasket = 21;
	      else if(nExtraLep==1) currentBasket = 22;
	      else if(nExtraLep>=2) currentBasket = 23;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else cout<<"WARNING : inconsistent nJetsBTagged"<<endl;
	  }else{
	    cout<<"WARNING : inconsistent nJets"<<endl;
	  }

	}else if(BASKETLIST==3){

	  int un  = 1;
	  int vbf = 2;
	  int vhl = 3;
	  int vhh = 4;
	  int tth = 5;
	  if(nJets==0){
	    if(nExtraLep==0) currentBasket = un;
	    else if(nExtraLep==1) currentBasket = vhl;
	    else if(nExtraLep>=2) currentBasket = vhl;
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else if(nJets==1){
	    if(nJetsBTagged==0){
	      if(nExtraLep==0) currentBasket = un;
	      else if(nExtraLep==1) currentBasket = vhl;
	      else if(nExtraLep>=2) currentBasket = vhl;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else if(nJetsBTagged==1){
	      if(nExtraLep==0) currentBasket = un;
	      else if(nExtraLep==1) currentBasket = tth; // can do better here
	      else if(nExtraLep>=2) currentBasket = tth; // can do better here
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else cout<<"WARNING : inconsistent nJetsBTagged"<<endl;
	  }else if(nJets>=2){
	    if(nJetsBTagged==0){
	      if(nExtraLep==0 && DiJetFisher>0.5) currentBasket = vbf; // can do better here
	      else if(nExtraLep==0 && HadVHTag) currentBasket = vhh; // can do better here
	      else if(nExtraLep==0) currentBasket = un; // can do better here
	      else if(nExtraLep==1) currentBasket = vhl; // can also choose vhl ?
	      else if(nExtraLep>=2) currentBasket = vhl; // can also choose vhl ?
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else if(nJetsBTagged==1){
	      if(nExtraLep==0 && DiJetFisher>0.5) currentBasket = vbf; // can do better here
	      else if(nExtraLep==0 && HadVHTag) currentBasket = vhh; // can do better here
	      else if(nExtraLep==0) currentBasket = un; // can do better here
	      else if(nExtraLep==1) currentBasket = tth;
	      else if(nExtraLep>=2) currentBasket = tth;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else if(nJetsBTagged>=2){
	      if(nExtraLep==0 && HadVHTag) currentBasket = vhh; // can do better here
	      else if(nExtraLep==0) currentBasket = tth; // can also choose un ?
	      else if(nExtraLep==1) currentBasket = tth;
	      else if(nExtraLep>=2) currentBasket = tth;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else cout<<"WARNING : inconsistent nJetsBTagged"<<endl;
	  }else{
	    cout<<"WARNING : inconsistent nJets"<<endl;
	  }

	}else if(BASKETLIST==4){

	  if(nJets==0){
	    if(nExtraLep==0) currentBasket = 1;
	    else if(nExtraLep==1) currentBasket = 2;
	    else if(nExtraLep>=2) currentBasket = 3;
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else if(nJets==1){
	    if(nJetsBTagged==0){
	      //if(nExtraLep==0 && abs(JetEta->at(0))>2.4 && ZZPt>50.) currentBasket = 24; else
	      if(nExtraLep==0) currentBasket = 4;
	      else if(nExtraLep==1) currentBasket = 5;
	      else if(nExtraLep>=2) currentBasket = 6;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else if(nJetsBTagged==1){
	      if(nExtraLep==0) currentBasket = 7;
	      else if(nExtraLep==1) currentBasket = 8;
	      else if(nExtraLep>=2) currentBasket = 9;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else cout<<"WARNING : inconsistent nJetsBTagged"<<endl;
	  }else if(nJets==2){
	    if(nJetsBTagged==0){
	      if(nExtraLep==0 && DiJetFisher>0.5) currentBasket = 10;
	      else if(nExtraLep==0 && HadVHTag) currentBasket = 11;
	      else if(nExtraLep==0) currentBasket = 12;
	      else if(nExtraLep==1) currentBasket = 13;
	      else if(nExtraLep>=2) currentBasket = 14;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else if(nJetsBTagged==1){
	      if(nExtraLep==0 && DiJetFisher>0.5) currentBasket = 15;
	      else if(nExtraLep==0 && HadVHTag) currentBasket = 16;
	      else if(nExtraLep==0) currentBasket = 17;
	      else if(nExtraLep==1) currentBasket = 18;
	      else if(nExtraLep>=2) currentBasket = 19;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else if(nJetsBTagged==2){
	      if(nExtraLep==0 && HadVHTag) currentBasket = 20;
	      else if(nExtraLep==0) currentBasket = 21;
	      else if(nExtraLep==1) currentBasket = 22;
	      else if(nExtraLep>=2) currentBasket = 23;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else cout<<"WARNING : inconsistent nJetsBTagged"<<endl;
	  }else if(nJets>=3){
	    if(nJetsBTagged==0){
	      if(nExtraLep==0 && DiJetFisher>0.5) currentBasket = 24;
	      else if(nExtraLep==0 && HadVHTag) currentBasket = 25;
	      else if(nExtraLep==0) currentBasket = 26;
	      else if(nExtraLep==1) currentBasket = 27;
	      else if(nExtraLep>=2) currentBasket = 28;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else if(nJetsBTagged==1){
	      if(nExtraLep==0 && DiJetFisher>0.5) currentBasket = 29;
	      else if(nExtraLep==0 && HadVHTag) currentBasket = 30;
	      else if(nExtraLep==0) currentBasket = 31;
	      else if(nExtraLep==1) currentBasket = 32;
	      else if(nExtraLep>=2) currentBasket = 33;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else if(nJetsBTagged>=2){
	      if(nExtraLep==0 && HadVHTag) currentBasket = 34;
	      else if(nExtraLep==0) currentBasket = 35;
	      else if(nExtraLep==1) currentBasket = 36;
	      else if(nExtraLep>=2) currentBasket = 37;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else cout<<"WARNING : inconsistent nJetsBTagged"<<endl;
	  }else{
	    cout<<"WARNING : inconsistent nJets"<<endl;
	  }

	}else if(BASKETLIST==5){

	  int un = 1;
	  int vbf1j = 2;
	  int vbf2j = 3;
	  int vhl = 4;
	  int vhh = 5;
	  int tth = 6;
	  if(nJets==0){
	    if(nExtraLep==0) currentBasket = un;
	    else if(nExtraLep==1) currentBasket = vhl;
	    else if(nExtraLep>=2) currentBasket = vhl;
	    else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }else if(nJets==1){
	    if(nJetsBTagged==0){
	      if(nExtraLep==0) currentBasket = vbf1j;
	      else if(nExtraLep==1) currentBasket = vhl;
	      else if(nExtraLep>=2) currentBasket = vhl;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else if(nJetsBTagged==1){
	      if(nExtraLep==0) currentBasket = vbf1j;
	      else if(nExtraLep==1) currentBasket = tth; // can do better here
	      else if(nExtraLep>=2) currentBasket = tth; // can do better here
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else cout<<"WARNING : inconsistent nJetsBTagged"<<endl;
	  }else if(nJets==2){
	    if(nJetsBTagged==0){
	      if(nExtraLep==0 && DiJetFisher>0.5) currentBasket = vbf2j;
	      else if(nExtraLep==0 && HadVHTag) currentBasket = vhh;
	      else if(nExtraLep==0) currentBasket = vbf1j;
	      else if(nExtraLep==1) currentBasket = vhl;
	      else if(nExtraLep>=2) currentBasket = vhl;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else if(nJetsBTagged==1){
	      if(nExtraLep==0 && DiJetFisher>0.5) currentBasket = vbf2j;
	      else if(nExtraLep==0 && HadVHTag) currentBasket = vhh;
	      else if(nExtraLep==0) currentBasket = vbf1j;
	      else if(nExtraLep==1) currentBasket = tth;
	      else if(nExtraLep>=2) currentBasket = tth;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else if(nJetsBTagged==2){
	      if(nExtraLep==0 && HadVHTag) currentBasket = vhh;
	      else if(nExtraLep==0) currentBasket = vhh;
	      else if(nExtraLep==1) currentBasket = tth;
	      else if(nExtraLep>=2) currentBasket = tth;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else cout<<"WARNING : inconsistent nJetsBTagged"<<endl;
	  }else if(nJets>=3){
	    if(nJetsBTagged==0){
	      if(nExtraLep==0 && DiJetFisher>0.5) currentBasket = vbf2j;
	      else if(nExtraLep==0 && HadVHTag) currentBasket = vhh;
	      else if(nExtraLep==0) currentBasket = vbf1j;
	      else if(nExtraLep==1) currentBasket = tth;
	      else if(nExtraLep>=2) currentBasket = tth;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else if(nJetsBTagged==1){
	      if(nExtraLep==0 && DiJetFisher>0.5) currentBasket = vbf2j;
	      else if(nExtraLep==0 && HadVHTag) currentBasket = vhh;
	      else if(nExtraLep==0) currentBasket = tth;
	      else if(nExtraLep==1) currentBasket = tth;
	      else if(nExtraLep>=2) currentBasket = tth;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else if(nJetsBTagged>=2){
	      if(nExtraLep==0 && HadVHTag) currentBasket = vhh;
	      else if(nExtraLep==0) currentBasket = tth;
	      else if(nExtraLep==1) currentBasket = tth;
	      else if(nExtraLep>=2) currentBasket = tth;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else cout<<"WARNING : inconsistent nJetsBTagged"<<endl;
	  }else{
	    cout<<"WARNING : inconsistent nJets"<<endl;
	  }

	}else if(BASKETLIST==6){

	  if(nExtraLep==0){ 
	    if(nJets==0) currentBasket = 1;
	    else if(nJets==1) currentBasket = 2;
	    else if(nJets>=2) currentBasket = 3;
	    else cout<<"WARNING : inconsistent nJets"<<endl;
	  }else if(nExtraLep==1){
	    if(nJets==0) currentBasket = 4;
	    else if(nJets==1) currentBasket = 5;
	    else if(nJets>=2) currentBasket = 6;
	    else cout<<"WARNING : inconsistent nJets"<<endl;
	  }else if(nExtraLep>=2){
	    if(nJets==0) currentBasket = 7;
	    else if(nJets==1) currentBasket = 8;
	    else if(nJets>=2) currentBasket = 9;
	    else cout<<"WARNING : inconsistent nJets"<<endl;
	  }else{
	    cout<<"WARNING : inconsistent nExtraLep"<<endl;
	  }

	}

	hBCFullSel100Baskets[p][rc]->Fill(currentBasket,eventWeight);
	hBCFullSel100Baskets[p][rc]->Fill(0.,eventWeight);
	hBCFullSel100Baskets[p][rc2]->Fill(currentBasket,eventWeight);
	hBCFullSel100Baskets[p][rc2]->Fill(0.,eventWeight);
	if(currentAssocDecay>=0){
	  hBCFullSel100BasketsAssocDecays[currentAssocDecay][rc]->Fill(currentBasket,eventWeight);
	  hBCFullSel100BasketsAssocDecays[currentAssocDecay][rc]->Fill(0.,eventWeight);
	  hBCFullSel100BasketsAssocDecays[currentAssocDecay][rc2]->Fill(currentBasket,eventWeight);
	  hBCFullSel100BasketsAssocDecays[currentAssocDecay][rc2]->Fill(0.,eventWeight);	      
	}
	if(106<ZZMass && ZZMass<141){
	  hBCFullSel100MasswindowBaskets[p][rc]->Fill(currentBasket,eventWeight);
	  hBCFullSel100MasswindowBaskets[p][rc]->Fill(0.,eventWeight);
	  hBCFullSel100MasswindowBaskets[p][rc2]->Fill(currentBasket,eventWeight);
	  hBCFullSel100MasswindowBaskets[p][rc2]->Fill(0.,eventWeight);
	  if(currentAssocDecay>=0){
	    hBCFullSel100MasswindowBasketsAssocDecays[currentAssocDecay][rc]->Fill(currentBasket,eventWeight);
	    hBCFullSel100MasswindowBasketsAssocDecays[currentAssocDecay][rc]->Fill(0.,eventWeight);
	    hBCFullSel100MasswindowBasketsAssocDecays[currentAssocDecay][rc2]->Fill(currentBasket,eventWeight);
	    hBCFullSel100MasswindowBasketsAssocDecays[currentAssocDecay][rc2]->Fill(0.,eventWeight);	      
	  }
	}

	
      } // end if(FullSel100)

    } // end for entries


    // ---------------------- Printing --------------------------

    Int_t widthColumn2 = 7;

    for(int rc=0; rc<nRecoChannels; rc++){
      txtOut<<" "<<recoChannels[rc]<<endl;
      txtOut<<"  has BC :     "<<fixWidth(Form("%i",nbWithBC[p][rc]),widthColumn2,false)<<endl;
      txtOut<<"  FullSel70 :  "<<fixWidth(Form("%i",nbWithBCFullSel70[p][rc]),widthColumn2,false)<<endl;
      txtOut<<"  FullSel100 : "<<fixWidth(Form("%i",nbWithBCFullSel100[p][rc]),widthColumn2,false)<<endl;
      if(requireExactly4GoodLeps) txtOut<<"  [ From this point, require that there is exactly 4 good leptons ]"<<endl;
      if(requireAtLeast5GoodLeps) txtOut<<"  [ From this point, require that there is at least 5 good leptons ]"<<endl;
      if(requireExactly5GoodLeps) txtOut<<"  [ From this point, require that there is exactly 5 good leptons ]"<<endl;
      if(requireExactly6GoodLeps) txtOut<<"  [ From this point, require that there is exactly 6 good leptons ]"<<endl;
      if(requireHLepsAreInEtaPtAcc   ) txtOut<<"  [ From this point, require that the 4 gen-leptons from the H are in the acceptance ]"<<endl;
      if(requireHLepsAreGood    ) txtOut<<"  [ From this point, require that the 4 gen-leptons from the H are reconstructed as good leptons ]"<<endl;
      txtOut<<"  FullSel100, there is ==4 / >=5 / ==5 / ==6 good leptons : "<<nbWithBCFullSel100Exactly4GoodLeps[p][rc]<<" / "<<nbWithBCFullSel100AtLeast5GoodLeps[p][rc]<<" / "<<nbWithBCFullSel100Exactly5GoodLeps[p][rc]<<" / "<<nbWithBCFullSel100Exactly6GoodLeps[p][rc]<<endl;
      txtOut<<"  FullSel100, the 4 gen-leptons from the H are in the acceptance : "<<nbWithBCFullSel100HLepsAreInEtaPtAcc[p][rc]<<endl;
      txtOut<<"  FullSel100, the 4 gen-leptons from the H are reconstructed as good leptons : "<<nbWithBCFullSel100HLepsAreGood[p][rc]<<endl;
      txtOut<<"  FullSel100, the 4 gen-leptons from the H are the 4 good leptons of the best candidate : "<<nbWithBCFullSel100All4LepRight[p][rc]<<endl;
    }

    int widthColumn1 = 22;
    int widthOtherColumns = 7;
    string separator = repeat("-",widthColumn1+4*(2+widthOtherColumns));
    if(prodName[p]=="WH"||prodName[p]=="WplusH"||prodName[p]=="WminusH"){
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<fixWidth("# of Z1 leptons from H",widthColumn1,true);
      for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(nbZDaughtersFromH[j],widthOtherColumns,false);
      txtOut<<endl;
      for(int i=0; i<nAssocWDecays; i++){
	txtOut<<" "<<fixWidth(WHdecays[i],widthColumn1,true);
	for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(percentage((float)nbZ1DaughtersFromHWH[i][j]/nbTotalWH[i])+" %",widthOtherColumns,false);
	txtOut<<endl;
      }
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<fixWidth("# of Z2 leptons from H",widthColumn1,true);
      for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(nbZDaughtersFromH[j],widthOtherColumns,false);
      txtOut<<endl;
      for(int i=0; i<nAssocWDecays; i++){
	txtOut<<" "<<fixWidth(WHdecays[i],widthColumn1,true);
	for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(percentage((float)nbZ2DaughtersFromHWH[i][j]/nbTotalWH[i])+" %",widthOtherColumns,false);
	txtOut<<endl;
      }
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<"the 4 gen-leptons from the H are in the acceptance :"<<endl;
      for(int i=0; i<nAssocWDecays; i++) txtOut<<" "<<fixWidth(WHdecays[i],widthColumn1,true)<<"   "<<percentage((float)nbHLepsAreInEtaPtAccWH[i]/nbTotalWH[i])+" %"<<endl;
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<"the 4 gen-leptons from the H are reconstructed as good leptons :"<<endl;
      for(int i=0; i<nAssocWDecays; i++) txtOut<<" "<<fixWidth(WHdecays[i],widthColumn1,true)<<"   "<<percentage((float)nbHLepsAreGoodWH[i]/nbTotalWH[i])+" %"<<endl;
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<"the 4 gen-leptons from the H are the 4 good leptons of the best candidate :"<<endl;
      for(int i=0; i<nAssocWDecays; i++) txtOut<<" "<<fixWidth(WHdecays[i],widthColumn1,true)<<"   "<<percentage((float)nbAll4LepRightWH[i]/nbTotalWH[i])+" %"<<endl;
    }
    if(prodName[p]=="ZH"){
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<fixWidth("# of Z1 leptons from H",widthColumn1,true);
      for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(nbZDaughtersFromH[j],widthOtherColumns,false);
      txtOut<<endl;
      for(int i=0; i<nAssocZDecays; i++){
	txtOut<<" "<<fixWidth(ZHdecays[i],widthColumn1,true);
	for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(percentage((float)nbZ1DaughtersFromHZH[i][j]/nbTotalZH[i])+" %",widthOtherColumns,false);
	txtOut<<endl;
      }
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<fixWidth("# of Z2 leptons from H",widthColumn1,true);
      for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(nbZDaughtersFromH[j],widthOtherColumns,false);
      txtOut<<endl;
      for(int i=0; i<nAssocZDecays; i++){
	txtOut<<" "<<fixWidth(ZHdecays[i],widthColumn1,true);
	for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(percentage((float)nbZ2DaughtersFromHZH[i][j]/nbTotalZH[i])+" %",widthOtherColumns,false);
	txtOut<<endl;
      }
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<"the 4 gen-leptons from the H are in the acceptance :"<<endl;
      for(int i=0; i<nAssocZDecays; i++) txtOut<<" "<<fixWidth(ZHdecays[i],widthColumn1,true)<<"   "<<percentage((float)nbHLepsAreInEtaPtAccZH[i]/nbTotalZH[i])+" %"<<endl;
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<"the 4 gen-leptons from the H are reconstructed as good leptons :"<<endl;
      for(int i=0; i<nAssocZDecays; i++) txtOut<<" "<<fixWidth(ZHdecays[i],widthColumn1,true)<<"   "<<percentage((float)nbHLepsAreGoodZH[i]/nbTotalZH[i])+" %"<<endl;
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<"the 4 gen-leptons from the H are the 4 good leptons of the best candidate :"<<endl;
      for(int i=0; i<nAssocZDecays; i++) txtOut<<" "<<fixWidth(ZHdecays[i],widthColumn1,true)<<"   "<<percentage((float)nbAll4LepRightZH[i]/nbTotalZH[i])+" %"<<endl;
    }
    if(prodName[p]=="ttH"){
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<fixWidth("# of Z1 leptons from H",widthColumn1,true);
      for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(nbZDaughtersFromH[j],widthOtherColumns,false);
      txtOut<<endl;
      for(int i=0; i<nAssocttDecays; i++){
	txtOut<<" "<<fixWidth(ttHdecays[i],widthColumn1,true);
	for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(percentage((float)nbZ1DaughtersFromHttH[i][j]/nbTotalttH[i])+" %",widthOtherColumns,false);
	txtOut<<endl;
      }
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<fixWidth("# of Z2 leptons from H",widthColumn1,true);
      for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(nbZDaughtersFromH[j],widthOtherColumns,false);
      txtOut<<endl;
      for(int i=0; i<nAssocttDecays; i++){
	txtOut<<" "<<fixWidth(ttHdecays[i],widthColumn1,true);
	for(int j=0; j<4; j++) txtOut<<"  "<<fixWidth(percentage((float)nbZ2DaughtersFromHttH[i][j]/nbTotalttH[i])+" %",widthOtherColumns,false);
	txtOut<<endl;
      }
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<"the 4 gen-leptons from the H are in the acceptance :"<<endl;
      for(int i=0; i<nAssocttDecays; i++) txtOut<<" "<<fixWidth(ttHdecays[i],widthColumn1,true)<<"   "<<percentage((float)nbHLepsAreInEtaPtAccttH[i]/nbTotalttH[i])+" %"<<endl;
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<"the 4 gen-leptons from the H are reconstructed as good leptons :"<<endl;
      for(int i=0; i<nAssocttDecays; i++) txtOut<<" "<<fixWidth(ttHdecays[i],widthColumn1,true)<<"   "<<percentage((float)nbHLepsAreGoodttH[i]/nbTotalttH[i])+" %"<<endl;
      txtOut<<" "<<separator<<endl;
      txtOut<<" "<<"the 4 gen-leptons from the H are the 4 good leptons of the best candidate :"<<endl;
      for(int i=0; i<nAssocttDecays; i++) txtOut<<" "<<fixWidth(ttHdecays[i],widthColumn1,true)<<"   "<<percentage((float)nbAll4LepRightttH[i]/nbTotalttH[i])+" %"<<endl;
    }
    
    txtOut<<endl;

  } // end for production modes

  txtOut.close();


  for(int p=0; p<nSamples; p++){
    if(!isPresent[p]) continue;
    //if(p!=0) continue;

    if(doYieldStudy){

      cout<<"For process "<<prodName[p]<<": "<<endl;
      cout<<" Nb generated / total yield (L*xsec)   :     "<<NGenEvt[p]<<"  /  "<<preselExp[p]<<endl;
      cout<<" Nb / Yield of stored events           :     "<<nbStored[p][10]<<"  /  "<<yieldStored[p][10]<<endl;
      cout<<" Nb / Yield of stored events (H->4l)   :     "<<nbStored[p][7]<<"  /  "<<yieldStored[p][7]<<endl;
      cout<<"   same, 4 gen l in eta acceptance     :     "<<nbHLepsAreInEtaAcc[p][7]<<"  /  "<<yieldHLepsAreInEtaAcc[p][7]<<endl;
      cout<<"   same, 4 gen l in eta+pt acceptance  :     "<<nbHLepsAreInEtaPtAcc[p][7]<<"  /  "<<yieldHLepsAreInEtaPtAcc[p][7]<<endl;
      cout<<"   same, >=4 reconstructed leptons     :     "<<nb4RecoLeps[p][7]<<"  /  "<<yield4RecoLeps[p][7]<<endl;
      cout<<"   same, 4 gen in eta+pt acc. + 4 reco :     "<<nbHLepsAreInEtaPtAcc4RecoLeps[p][7]<<"  /  "<<yieldHLepsAreInEtaPtAcc4RecoLeps[p][7]<<endl;
      //cout<<" Nb / Yield of stored events (H->2L2t) :     "<<nbStored[p][8]<<"  /  "<<yieldStored[p][8]<<endl;
      //cout<<" Nb / Yield of stored events (H->else) :     "<<nbStored[p][9]<<"  /  "<<yieldStored[p][9]<<endl;
      //cout<<" Nb / Yield of stored events (HX->4lX) :     "<<nb4GenLeps[p][10]<<"  /  "<<yield4GenLeps[p][10]<<endl;
      //cout<<" Nb / Yield of events with BC          :     "<<nbWithBC[p][3]<<"  /  "<<yieldWithBC[p][3]<<endl;
      cout<<" Nb / Yield of events with FullSel70   :     "<<nbWithBCFullSel70[p][3]<<"  /  "<<yieldWithBCFullSel70[p][3]<<endl;
      cout<<" Nb / Yield of events with FullSel100  :     "<<nbWithBCFullSel100[p][3]<<"  /  "<<yieldWithBCFullSel100[p][3]<<endl;

      Int_t wi = 6;
      cout<<" Generated final state             ";
      cout<<fixWidth("   4mu",2*wi+1,1)<<"   ";
      cout<<fixWidth("   4e",2*wi+1,1)<<"   ";
      cout<<fixWidth("  2e2mu",2*wi+1,1)<<"   ";
      cout<<fixWidth("   all",2*wi+1,1)<<"   ";
      cout<<endl;
      cout<<" "<<repeat("-",92)<<endl;
      cout<<"                                   ";
      cout<<fixWidth("yield",wi,1)<<"  "<<fixWidth(" eff.",wi,1)<<"   ";
      cout<<fixWidth("yield",wi,1)<<"  "<<fixWidth(" eff.",wi,1)<<"   ";
      cout<<fixWidth("yield",wi,1)<<"  "<<fixWidth(" eff.",wi,1)<<"   ";
      cout<<fixWidth("yield",wi,1)<<"  "<<fixWidth(" eff.",wi,1)<<"   ";
      cout<<endl;
      cout<<" All generated events              ";
      cout<<fixWidth(rounding2(yieldStored[p][0]),wi,0)<<"  "<<repeat(" ",wi)<<"   ";
      cout<<fixWidth(rounding2(yieldStored[p][1]),wi,0)<<"  "<<repeat(" ",wi)<<"   ";
      cout<<fixWidth(rounding2(yieldStored[p][2]),wi,0)<<"  "<<repeat(" ",wi)<<"   ";
      cout<<fixWidth(rounding2(yieldStored[p][7]),wi,0)<<"  "<<repeat(" ",wi)<<"   ";
      cout<<endl;
      cout<<" 4 gen l in eta acceptance         ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaAcc[p][0]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaAcc[p][0]/yieldStored[p][0])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaAcc[p][1]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaAcc[p][1]/yieldStored[p][1])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaAcc[p][2]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaAcc[p][2]/yieldStored[p][2])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaAcc[p][7]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaAcc[p][7]/yieldStored[p][7])+"%",wi,0)<<"   ";
      cout<<endl;
      cout<<" 4 gen l in eta+pt acceptance      ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc[p][0]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc[p][0]/yieldHLepsAreInEtaAcc[p][0])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc[p][1]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc[p][1]/yieldHLepsAreInEtaAcc[p][1])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc[p][2]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc[p][2]/yieldHLepsAreInEtaAcc[p][2])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc[p][7]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc[p][7]/yieldHLepsAreInEtaAcc[p][7])+"%",wi,0)<<"   ";
      cout<<endl;
      cout<<"  (same but eff wrt. all)          ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc[p][0]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc[p][0]/yieldStored[p][0])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc[p][1]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc[p][1]/yieldStored[p][1])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc[p][2]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc[p][2]/yieldStored[p][2])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc[p][7]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc[p][7]/yieldStored[p][7])+"%",wi,0)<<"   ";
      cout<<endl;
      cout<<" eta+pt acc + >=4 reco leptons     ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc4RecoLeps[p][0]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc4RecoLeps[p][0]/yieldHLepsAreInEtaPtAcc[p][0])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc4RecoLeps[p][1]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc4RecoLeps[p][1]/yieldHLepsAreInEtaPtAcc[p][1])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc4RecoLeps[p][2]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc4RecoLeps[p][2]/yieldHLepsAreInEtaPtAcc[p][2])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAcc4RecoLeps[p][7]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAcc4RecoLeps[p][7]/yieldHLepsAreInEtaPtAcc[p][7])+"%",wi,0)<<"   ";
      cout<<endl;
      cout<<" eta+pt acc + FullSel100           ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccWithBCFullSel100G[p][0]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccWithBCFullSel100G[p][0]/yieldHLepsAreInEtaPtAcc4RecoLeps[p][0])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccWithBCFullSel100G[p][1]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccWithBCFullSel100G[p][1]/yieldHLepsAreInEtaPtAcc4RecoLeps[p][1])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccWithBCFullSel100G[p][2]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccWithBCFullSel100G[p][2]/yieldHLepsAreInEtaPtAcc4RecoLeps[p][2])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccWithBCFullSel100G[p][7]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccWithBCFullSel100G[p][7]/yieldHLepsAreInEtaPtAcc4RecoLeps[p][7])+"%",wi,0)<<"   ";
      cout<<endl;
      cout<<" eta+pt acc + FullSel100 + trigger ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGWithBCFullSel100G[p][0]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGWithBCFullSel100G[p][0]/yieldHLepsAreInEtaPtAccWithBCFullSel100G[p][0])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGWithBCFullSel100G[p][1]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGWithBCFullSel100G[p][1]/yieldHLepsAreInEtaPtAccWithBCFullSel100G[p][1])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGWithBCFullSel100G[p][2]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGWithBCFullSel100G[p][2]/yieldHLepsAreInEtaPtAccWithBCFullSel100G[p][2])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGWithBCFullSel100G[p][7]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGWithBCFullSel100G[p][7]/yieldHLepsAreInEtaPtAccWithBCFullSel100G[p][7])+"%",wi,0)<<"   ";
      cout<<endl;
      cout<<"  (same without single ele path)   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCFullSel100G[p][0]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCFullSel100G[p][0]/yieldHLepsAreInEtaPtAccWithBCFullSel100G[p][0])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCFullSel100G[p][1]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCFullSel100G[p][1]/yieldHLepsAreInEtaPtAccWithBCFullSel100G[p][1])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCFullSel100G[p][2]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCFullSel100G[p][2]/yieldHLepsAreInEtaPtAccWithBCFullSel100G[p][2])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCFullSel100G[p][7]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGNo1EWithBCFullSel100G[p][7]/yieldHLepsAreInEtaPtAccWithBCFullSel100G[p][7])+"%",wi,0)<<"   ";
      cout<<endl;
      cout<<" eta+pt acc + trigger (% wrt. acc) ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerG[p][0]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerG[p][0]/yieldHLepsAreInEtaPtAcc[p][0])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerG[p][1]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerG[p][1]/yieldHLepsAreInEtaPtAcc[p][1])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerG[p][2]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerG[p][2]/yieldHLepsAreInEtaPtAcc[p][2])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerG[p][7]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerG[p][7]/yieldHLepsAreInEtaPtAcc[p][7])+"%",wi,0)<<"   ";
      cout<<endl;
      cout<<"  (same without single ele path)   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGNo1E[p][0]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGNo1E[p][0]/yieldHLepsAreInEtaPtAcc[p][0])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGNo1E[p][1]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGNo1E[p][1]/yieldHLepsAreInEtaPtAcc[p][1])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGNo1E[p][2]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGNo1E[p][2]/yieldHLepsAreInEtaPtAcc[p][2])+"%",wi,0)<<"   ";
      cout<<fixWidth(rounding2(yieldHLepsAreInEtaPtAccPassTriggerGNo1E[p][7]),wi,0)<<"  "<<fixWidth(percentage(yieldHLepsAreInEtaPtAccPassTriggerGNo1E[p][7]/yieldHLepsAreInEtaPtAcc[p][7])+"%",wi,0)<<"   ";
      cout<<endl;
      cout<<" "<<repeat("-",92)<<endl;

    }

    if(doTTHCheck && prodName[p]=="ttH"){
      cout<<"Checking generated final states of ttH wrong signal (H -> not 4l) : "<<endl;
      for(int nt=0; nt<nTTHGenStatesUsed; nt++){
	//cout<<labelCheckTTH[nt]<<nbCheckTTH[nt]<<"  /  "<<yieldCheckTTH[nt]<<endl;
	cout<<labelCheckTTH[nt]<<rounding3(yieldCheckTTH[nt])<<endl;
      }
    }

  }
  

  // ------------------------------------------------------------
  // -------------------------- Plots ---------------------------
  // ------------------------------------------------------------

  if(doProdComp){
    TCanvas* cBCFullSel100[nVariables];
    for(int v=0; v<nVariables; v++){
      if(!plotThisVar[VARIABLELIST][v]) continue;
      cBCFullSel100[v] = new TCanvas(Form("cBCFullSel100_%s",varName[v].c_str()),Form("cBCFullSel100_%s",varName[v].c_str()),500,500);
      DrawByProdmodes(cBCFullSel100[v],hBCFullSel100[v],prodName,isPresent,colors,v==0);
      SaveCanvas(outDir,cBCFullSel100[v]);
    }
  }

  if(doProdCompMatch4){
    TCanvas* cBCFullSel100Match4[nVariables];
    for(int v=0; v<nVariables; v++){
      if(!plotThisVar[VARIABLELIST][v]) continue;
      cBCFullSel100Match4[v] = new TCanvas(Form("cBCFullSel100Match4_%s",varName[v].c_str()),Form("cBCFullSel100Match4_%s",varName[v].c_str()),500,500);
      DrawByProdmodes(cBCFullSel100Match4[v],hBCFullSel100MatchHLeps[v][0],prodName,isPresent,colors,v==0);
      SaveCanvas(outDir,cBCFullSel100Match4[v]);
    }
  }

  if(doMatch4OrNot){
    TCanvas* cBCFullSel100Match4OrNot[nVariables][nHiggsSamples];
    for(int p=0; p<nHiggsSamples; p++){
      if(!isPresent[p]) continue;
      string pn = prodName[p];
      for(int v=0; v<nVariables; v++){
	if(!plotThisVar[VARIABLELIST][v]) continue;
	cBCFullSel100Match4OrNot[v][p] = new TCanvas(Form("cBCFullSel100Match4OrNot_%s_%s",varName[v].c_str(),pn.c_str()),Form("cBCFullSel100Match4OrNot_%s_%s",varName[v].c_str(),pn.c_str()),500,500);
	DrawMatch4OrNot(cBCFullSel100Match4OrNot[v][p],hBCFullSel100[v][p][FINALSTATE],hBCFullSel100MatchHLeps[v][0][p][FINALSTATE],pn,allColorsMP[1][p],allColorsMP[4][p],v==0);
	SaveCanvas(outDir,cBCFullSel100Match4OrNot[v][p]);
      }
    }
  }

  if(doMatchHLeps){
    TCanvas* cBCFullSel100MatchHLeps[nVariables][nHiggsSamples];
    for(int p=0; p<nHiggsSamples; p++){
      if(!isPresent[p]) continue;
      string pn = prodName[p];
      for(int v=0; v<nVariables; v++){
	if(!plotThisVar[VARIABLELIST][v]) continue;
	TH1F* h[nMatchHLepsStatuses]; for(int m=0; m<nMatchHLepsStatuses; m++) h[m] = (TH1F*)hBCFullSel100MatchHLeps[v][m][p][FINALSTATE];
	cBCFullSel100MatchHLeps[v][p] = new TCanvas(Form("cBCFullSel100MatchHLeps_%s_%s",varName[v].c_str(),pn.c_str()),Form("cBCFullSel100MatchHLeps_%s_%s",varName[v].c_str(),pn.c_str()),500,500);
	DrawMatchHLeps(cBCFullSel100MatchHLeps[v][p],hBCFullSel100[v][p][FINALSTATE],h,pn,allColorsPM1[p],matchHLepsKeys,v==0);
	SaveCanvas(outDir,cBCFullSel100MatchHLeps[v][p]);
      }
    }
  }

  if(doMatchAllLeps){
    TCanvas* cBCFullSel100MatchAllLeps[nVariables][nHiggsSamples];
    for(int p=0; p<nHiggsSamples; p++){
      if(!isPresent[p]) continue;
      if(p==0||p==1) continue;
      string pn = prodName[p];
      for(int v=0; v<nVariables; v++){
	if(!plotThisVar[VARIABLELIST][v]) continue;
	TH1F* h[nMatchAllLepsStatuses]; for(int m=0; m<nMatchAllLepsStatuses; m++) h[m] = (TH1F*)hBCFullSel100MatchAllLeps[v][m][p][FINALSTATE];
	cBCFullSel100MatchAllLeps[v][p] = new TCanvas(Form("cBCFullSel100MatchAllLeps_%s_%s",varName[v].c_str(),pn.c_str()),Form("cBCFullSel100MatchAllLeps_%s_%s",varName[v].c_str(),pn.c_str()),500,500);
	DrawMatchAllLeps(cBCFullSel100MatchAllLeps[v][p],hBCFullSel100[v][p][FINALSTATE],h,pn,allColorsPM2[p],matchAllLepsKeys,v==0);
	SaveCanvas(outDir,cBCFullSel100MatchAllLeps[v][p]);
      }
    }
  }
  
  if(doMatchWHZHttH){

    if(file_WplusH!="" && file_WminusH!=""){
      TCanvas* cBCFullSel100MatchWH[nVariables];
      for(int v=0; v<nVariables; v++){
	if(!plotThisVar[VARIABLELIST][v]) continue;
	TH1F* hMatchWH[nMatchWHStatuses]; for(int m=0; m<nMatchWHStatuses; m++) hMatchWH[m] = (TH1F*)hBCFullSel100MatchWH[v][m][FINALSTATE];
	cBCFullSel100MatchWH[v] = new TCanvas(Form("cBCFullSel100MatchWH_%s",varName[v].c_str()),Form("cBCFullSel100MatchWH_%s",varName[v].c_str()),500,500);
	if(v==1 || v==3)
	  DrawMatchCustom(cBCFullSel100MatchWH[v],hBCFullSel100[v][2][FINALSTATE],hMatchWH,prodName[2],colorsMatchWH,matchWH,nMatchWHStatuses,0.11,0.6,v==0);
	else
	  DrawMatchCustom(cBCFullSel100MatchWH[v],hBCFullSel100[v][2][FINALSTATE],hMatchWH,prodName[2],colorsMatchWH,matchWH,nMatchWHStatuses,0.34,0.6,v==0);
	SaveCanvas(outDir,cBCFullSel100MatchWH[v]);
      }
    }

    if(file_ZH!=""){
      TCanvas* cBCFullSel100MatchZH[nVariables];
      for(int v=0; v<nVariables; v++){
	if(!plotThisVar[VARIABLELIST][v]) continue;
	TH1F* hMatchZH[nMatchZHStatuses]; for(int m=0; m<nMatchZHStatuses; m++) hMatchZH[m] = (TH1F*)hBCFullSel100MatchZH[v][m][FINALSTATE];
	cBCFullSel100MatchZH[v] = new TCanvas(Form("cBCFullSel100MatchZH_%s",varName[v].c_str()),Form("cBCFullSel100MatchZH_%s",varName[v].c_str()),500,500);
	if(v==1 || v==3)
	  DrawMatchCustom(cBCFullSel100MatchZH[v],hBCFullSel100[v][3][FINALSTATE],hMatchZH,prodName[3],colorsMatchZH,matchZH,nMatchZHStatuses,0.11,0.5,v==0);
	else
	  DrawMatchCustom(cBCFullSel100MatchZH[v],hBCFullSel100[v][3][FINALSTATE],hMatchZH,prodName[3],colorsMatchZH,matchZH,nMatchZHStatuses,0.34,0.5,v==0);
	SaveCanvas(outDir,cBCFullSel100MatchZH[v]);
      }
    }

    if(file_ttH!=""){
      TCanvas* cBCFullSel100MatchttH[nVariables];
      for(int v=0; v<nVariables; v++){
	if(!plotThisVar[VARIABLELIST][v]) continue;
	TH1F* hMatchttH[nMatchttHStatuses]; for(int m=0; m<nMatchttHStatuses; m++) hMatchttH[m] = (TH1F*)hBCFullSel100MatchttH[v][m][FINALSTATE];
	cBCFullSel100MatchttH[v] = new TCanvas(Form("cBCFullSel100MatchttH_%s",varName[v].c_str()),Form("cBCFullSel100MatchttH_%s",varName[v].c_str()),500,500);
	if(v==1 || v==3)
	  DrawMatchCustom(cBCFullSel100MatchttH[v],hBCFullSel100[v][4][FINALSTATE],hMatchttH,prodName[4],colorsMatchttH,matchttH,nMatchttHStatuses,0.11,0.4,v==0);
	else
	  DrawMatchCustom(cBCFullSel100MatchttH[v],hBCFullSel100[v][4][FINALSTATE],hMatchttH,prodName[4],colorsMatchttH,matchttH,nMatchttHStatuses,0.34,0.4,v==0);
	SaveCanvas(outDir,cBCFullSel100MatchttH[v]);
      }
    }

  }
  
  if(do2DPlots){

    TCanvas* c2DBCFullSel100[n2DHist][nSamples];
    for(int v2=0; v2<n2DHist; v2++){
      for(int p=0; p<nSamples; p++){
	if(!isPresent[p]) continue;
	string name = "c2DBCFullSel100_"+varName[varXindex[v2]]+"_"+varName[varYindex[v2]]+"_"+prodName[p];
	c2DBCFullSel100[v2][p] = new TCanvas(name.c_str(),name.c_str(),600,400);
	Draw2D(c2DBCFullSel100[v2][p],h2DBCFullSel100[v2][p][FINALSTATE],prodName[p]);
	SaveCanvas(outDir,c2DBCFullSel100[v2][p]);
      }
    }

    TCanvas* c2DBCFullSel100Decays[n2DHist][nSamples][nDecays];
    for(int v2=0; v2<n2DHist; v2++){
      for(int p=0; p<nSamples; p++){
	if(!isPresent[p]) continue;
	if(p!=2&&p!=3&&p!=4) continue;
	for(int d=0; d<nDecays; d++){
	  if(p==2&&(d==0||d==2)) continue;
	  string name = "c2DBCFullSel100_"+decayInfix[d]+"_"+varName[varXindex[v2]]+"_"+varName[varYindex[v2]]+"_"+prodName[p];
	  c2DBCFullSel100Decays[v2][p][d] = new TCanvas(name.c_str(),name.c_str(),600,400);
	  Draw2D(c2DBCFullSel100Decays[v2][p][d],h2DBCFullSel100Decays[v2][p][d][FINALSTATE],prodName[p]+decayLabel[d]);
	  SaveCanvas(outDir,c2DBCFullSel100Decays[v2][p][d]);
	}
      }
    }

  }
  
  if(doBaskets){

    TCanvas* cBCFullSel100BasketsAll = new TCanvas("cBCFullSel100_baskets","cBCFullSel100_baskets",500,500);
    DrawByProdmodes(cBCFullSel100BasketsAll,hBCFullSel100Baskets,prodName,isPresent,colors,false);
    SaveCanvas(outDir,cBCFullSel100BasketsAll);
    
    TCanvas* cBCFullSel100BasketEfficiency[nSamples];
    TH1F* hProdModes[nSamples]; for(int p=0; p<nSamples; p++) hProdModes[p] = (TH1F*)hBCFullSel100Baskets[p][FINALSTATE];
    TH1F* hProdModesMasswindow[nSamples]; for(int p=0; p<nSamples; p++) hProdModesMasswindow[p] = (TH1F*)hBCFullSel100MasswindowBaskets[p][FINALSTATE];
    TH1F* hAssocDecays[nAssocDecays]; for(int a=0; a<nAssocDecays; a++) hAssocDecays[a] = (TH1F*)hBCFullSel100BasketsAssocDecays[a][FINALSTATE];
    TH1F* hAssocDecaysMasswindow[nAssocDecays]; for(int a=0; a<nAssocDecays; a++) hAssocDecaysMasswindow[a] = (TH1F*)hBCFullSel100MasswindowBasketsAssocDecays[a][FINALSTATE];
    for(int p=0; p<nSamples; p++){
      if(!isPresent[p]) continue;
      string pn = prodName[p];
      if(groupWminusWplus && pn=="WminusH"){
	continue;
      }else if(groupWminusWplus && pn=="WplusH"){
	hProdModes[p]->Add(hProdModes[p+1]);
	pn = "WH";
      }
      cBCFullSel100BasketEfficiency[p] = new TCanvas(Form("cBCFullSel100_basketEfficiency_%s",pn.c_str()),Form("cBCFullSel100_basketEfficiency_%s",pn.c_str()),500,((BASKETLIST==2||BASKETLIST==4)?1000:500));
      DrawBasketEfficiencies(cBCFullSel100BasketEfficiency[p],hProdModes[p],pn,hAssocDecays,treatH2l2XAsBkgd?assocDecayName2:assocDecayName1,treatH2l2XAsBkgd,basketLabel[BASKETLIST],false);
      //DrawBasketEfficiencies(cBCFullSel100BasketEfficiency[p],hProdModesMasswindow[p],pn,hAssocDecaysMasswindow,treatH2l2XAsBkgd?assocDecayName2:assocDecayName1,treatH2l2XAsBkgd,basketLabel[BASKETLIST],false);
      SaveCanvas(outDir,cBCFullSel100BasketEfficiency[p]);
    }

    if(allSignalsArePresent){
      TCanvas* cBCFullSel100BasketPurity = new TCanvas("cBCFullSel100_basketPurity","cBCFullSel100_basketPurity",800,((BASKETLIST==2||BASKETLIST==4)?1000:500));
      DrawBasketPurities(cBCFullSel100BasketPurity,hProdModes,prodName,isPresent,hAssocDecays,treatH2l2XAsBkgd?assocDecayName4:assocDecayName3,false,treatH2l2XAsBkgd,basketLabel[BASKETLIST],lumi);
      //DrawBasketPurities(cBCFullSel100BasketPurity,hProdModesMasswindow,prodName,isPresent,hAssocDecaysMasswindow,treatH2l2XAsBkgd?assocDecayName4:assocDecayName3,false,treatH2l2XAsBkgd,basketLabel[BASKETLIST],lumi);
      SaveCanvas(outDir,cBCFullSel100BasketPurity);
      TCanvas* cBCFullSel100BasketPuritySplit = new TCanvas("cBCFullSel100_basketPuritySplit","cBCFullSel100_basketPuritySplit",800,((BASKETLIST==2||BASKETLIST==4)?1000:500));
      DrawBasketPurities(cBCFullSel100BasketPuritySplit,hProdModes,prodName,isPresent,hAssocDecays,treatH2l2XAsBkgd?assocDecayName4:assocDecayName3,true,treatH2l2XAsBkgd,basketLabel[BASKETLIST],lumi);
      //DrawBasketPurities(cBCFullSel100BasketPuritySplit,hProdModesMasswindow,prodName,isPresent,hAssocDecaysMasswindow,treatH2l2XAsBkgd?assocDecayName4:assocDecayName3,true,treatH2l2XAsBkgd,basketLabel[BASKETLIST],lumi);
      SaveCanvas(outDir,cBCFullSel100BasketPuritySplit);
    }

    if(allSignalsArePresent && ZZBkgdIsPresent){

      TH1F* hSumSgnl = (TH1F*)hBCFullSel100Baskets[0][FINALSTATE]->Clone();
      for(int p=1; p<nHiggsSamples; p++){
	if(!isPresent[p]) continue;
	hSumSgnl->Add(hBCFullSel100Baskets[p][FINALSTATE]);
      }
      TH1F* hSumBkgd = (TH1F*)hBCFullSel100Baskets[nHiggsSamples][FINALSTATE]->Clone();
      if(treatH2l2XAsBkgd){
	hSumSgnl->Add(hBCFullSel100BasketsAssocDecays[nAssocWDecays+nAssocZDecays-1][FINALSTATE],-1);
	hSumSgnl->Add(hBCFullSel100BasketsAssocDecays[nAssocWDecays+nAssocZDecays+nAssocttDecays-1][FINALSTATE],-1);
	hSumBkgd->Add(hBCFullSel100BasketsAssocDecays[nAssocWDecays+nAssocZDecays-1][FINALSTATE]);
	hSumBkgd->Add(hBCFullSel100BasketsAssocDecays[nAssocWDecays+nAssocZDecays+nAssocttDecays-1][FINALSTATE]);
      }

      TH1F* hDenom = (TH1F*)hSumSgnl->Clone();
      hDenom->Add(hSumBkgd);
      TH1F* hSOSPB = (TH1F*)hSumSgnl->Clone();
      hSOSPB->Divide(hDenom);
      TCanvas* cBCFullSel100BasketSOSPB = new TCanvas("cBCFullSel100_basketSOSPB","cBCFullSel100_basketSOSPB",500,((BASKETLIST==2||BASKETLIST==4)?1000:500));
      DrawBasketSOSPB(cBCFullSel100BasketSOSPB,hSOSPB,basketLabel[BASKETLIST]);
      SaveCanvas(outDir,cBCFullSel100BasketSOSPB);

      TH1F* hSumSgnlMasswindow = (TH1F*)hBCFullSel100MasswindowBaskets[0][FINALSTATE]->Clone();
      for(int p=1; p<nHiggsSamples; p++){
	if(!isPresent[p]) continue;
	hSumSgnlMasswindow->Add(hBCFullSel100MasswindowBaskets[p][FINALSTATE]);
      }
      TH1F* hSumBkgdMasswindow = (TH1F*)hBCFullSel100MasswindowBaskets[nHiggsSamples][FINALSTATE]->Clone();
      if(treatH2l2XAsBkgd){
	hSumSgnlMasswindow->Add(hBCFullSel100MasswindowBasketsAssocDecays[nAssocWDecays+nAssocZDecays-1][FINALSTATE],-1);
	hSumSgnlMasswindow->Add(hBCFullSel100MasswindowBasketsAssocDecays[nAssocWDecays+nAssocZDecays+nAssocttDecays-1][FINALSTATE],-1);
	hSumBkgdMasswindow->Add(hBCFullSel100MasswindowBasketsAssocDecays[nAssocWDecays+nAssocZDecays-1][FINALSTATE]);
	hSumBkgdMasswindow->Add(hBCFullSel100MasswindowBasketsAssocDecays[nAssocWDecays+nAssocZDecays+nAssocttDecays-1][FINALSTATE]);
      }

      TH1F* hDenomMasswindow = (TH1F*)hSumSgnlMasswindow->Clone();
      hDenomMasswindow->Add(hSumBkgdMasswindow);
      TH1F* hSOSPBMasswindow = (TH1F*)hSumSgnlMasswindow->Clone();
      hSOSPBMasswindow->Divide(hDenomMasswindow);
      TCanvas* cBCFullSel100MasswindowBasketSOSPB = new TCanvas("cBCFullSel100Masswindow_basketSOSPB","cBCFullSel100Masswindow_basketSOSPB",500,((BASKETLIST==2||BASKETLIST==4)?1000:500));
      DrawBasketSOSPB(cBCFullSel100MasswindowBasketSOSPB,hSOSPBMasswindow,basketLabel[BASKETLIST]);
      SaveCanvas(outDir,cBCFullSel100MasswindowBasketSOSPB);

    }

    for(int p=0; p<nSamples; p++) 
      if(isPresent[p]) cout<<"process "<<prodName[p]<<": "<<0.5*hBCFullSel100Baskets[p][FINALSTATE]->Integral()<<endl;

    /*
    for(int p=0; p<nSamples; p++) 
      if(isPresent[p]) cout<<"process "<<prodName[p]<<", VH-leptonic: "<<hBCFullSel100Baskets[p][FINALSTATE]->GetBinContent(5)<<endl;
    for(int a=0; a<nAssocDecays; a++)
      cout<<"assoc decay "<<assocDecayName2[a]<<", VH-leptonic: "<<hBCFullSel100BasketsAssocDecays[a][FINALSTATE]->GetBinContent(5)<<endl;
    //*/

  }

  if(doPlotsByGenChan){

    TCanvas* cSortedPt [nHiggsSamples][nGenChannels];
    TCanvas* cSortedEta[nHiggsSamples][nGenChannels];
    TCanvas* cNRecoLep[nHiggsSamples][nGenChannels];
    TCanvas* cNRecoMuo[nHiggsSamples][nGenChannels];
    TCanvas* cNRecoEle[nHiggsSamples][nGenChannels];
    TCanvas* cGenHPt[nHiggsSamples];
    TCanvas* cGenHRapidity[nHiggsSamples];
    TCanvas* cGenHRapidityVsPt[nHiggsSamples];

    for(int p=0; p<nHiggsSamples; p++){
      if(!isPresent[p]) continue;
      if(!(prodName[p]=="ggH")) continue;

      for(int gc=0; gc<3; gc++){
	string suffix = prodName[p]+"_"+genChannels[gc];

	cSortedPt[p][gc] = new TCanvas(("cSortedPt_"+suffix).c_str(),("cSortedPt_"+suffix).c_str(),500,500);
	DrawSorted(cSortedPt[p][gc],h1GenptHLepsAreInEtaAcc[p][gc],h1RecoptWithBCFullSel100GPassTriggerG[p][gc],prodName[p],genChannels[gc],0.1);
	SaveCanvas(outDir,cSortedPt[p][gc],tagOut);
	cSortedEta[p][gc] = new TCanvas(("cSortedEta_"+suffix).c_str(),("cSortedEta_"+suffix).c_str(),500,500);
	DrawSorted(cSortedEta[p][gc],h1GenetaHLepsAreInPtAcc[p][gc],h1RecoetaWithBCFullSel100GPassTriggerG[p][gc],prodName[p],genChannels[gc],0.1);
	SaveCanvas(outDir,cSortedEta[p][gc],tagOut);

	cNRecoLep[p][gc] = new TCanvas(("cNRecoLep_"+suffix).c_str(),("cNRecoLep_"+suffix).c_str(),500,500);
	Draw1D(cNRecoLep[p][gc],h1NRecoLepHLepsAreInEtaPtAcc[p][gc]);
	SaveCanvas(outDir,cNRecoLep[p][gc],tagOut);
	cNRecoMuo[p][gc] = new TCanvas(("cNRecoMuo_"+suffix).c_str(),("cNRecoMuo_"+suffix).c_str(),500,500);
	Draw1D(cNRecoMuo[p][gc],h1NRecoMuoHLepsAreInEtaPtAcc[p][gc]);
	SaveCanvas(outDir,cNRecoMuo[p][gc],tagOut);
	cNRecoEle[p][gc] = new TCanvas(("cNRecoEle_"+suffix).c_str(),("cNRecoEle_"+suffix).c_str(),500,500);
	Draw1D(cNRecoEle[p][gc],h1NRecoEleHLepsAreInEtaPtAcc[p][gc]);
	SaveCanvas(outDir,cNRecoEle[p][gc],tagOut);

      }

      cGenHPt[p] = new TCanvas(("cGenHPt_"+prodName[p]).c_str(),("cGenHPt_"+prodName[p]).c_str(),500,500);
      Draw1D(cGenHPt[p],h1GenHPt[p]);
      SaveCanvas(outDir,cGenHPt[p],tagOut);
      cGenHRapidity[p] = new TCanvas(("cGenHRapidity_"+prodName[p]).c_str(),("cGenHRapidity_"+prodName[p]).c_str(),500,500);
      Draw1D(cGenHRapidity[p],h1GenHRapidity[p]);
      SaveCanvas(outDir,cGenHRapidity[p],tagOut);
      cGenHRapidityVsPt[p] = new TCanvas(("cGenHRapidityVsPt_"+prodName[p]).c_str(),("cGenHRapidityVsPt_"+prodName[p]).c_str(),500,500);
      Draw2D(cGenHRapidityVsPt[p],h2GenHRapidityVsPt[p],"");
      SaveCanvas(outDir,cGenHRapidityVsPt[p],tagOut);

    }

  }

}





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// MAIN MACRO //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void plotProcesses() {

  string inputPath = ""; // without "/" at the end
  string outputPath = "PlotsProcesses/";
  gSystem->Exec(("mkdir -p "+outputPath).c_str());

  float lumi = 15.;
		  
  run( 
      inputPath,
      outputPath,
      lumi,

      "ggH125/ZZ4lAnalysis.root",
      "VBFH125/ZZ4lAnalysis.root",
      "WplusH125/ZZ4lAnalysis.root",
      "WminusH125/ZZ4lAnalysis.root",
      "ZH125/ZZ4lAnalysis.root",
      "ttH125/ZZ4lAnalysis.root",
      "",
      "",//"ZZTo4l/ZZ4lAnalysis.root",

      "",
      ""
       );
  
}
