#include <iostream>
#include <fstream>
#include <utility>
#include <map>
#include <vector>
#include <cmath>

#include "TStyle.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TColor.h"
#include "TSystem.h"

using namespace std;

#define MASSZ 91.1876
#define LUMI 20.

#define DEBUG 0

#define VARIABLELIST 0 // put 0 for plots shown on September 8th

#define requireHLepsAreInAcc    0 // impose that all 4 gen leptons from the Higgs are in the acceptance
#define requireExactly4GoodLeps 0 
#define require5GoodLepsOrMore  0 
#define requireExactly5GoodLeps 0 
#define requireExactly6GoodLeps 0 
#define requireHLepsAreGood     0 // impose that all 4 gen leptons from the Higgs are matched to good leptons (from the candidate or extra leptons)

#define excludeH2l2X            1 // completely exclude ZH,H->2l2X and ttH,H->2l2X events from the study
#define acceptanceIncludesPt    1 // if not, acceptance is just on eta
#define treatH2l2XAsBkgd        0 // treat these events as a background for the purity and S/(S+B) plots 

#define doProdComp       1 // plots shown on September 8th
#define doProdCompMatch4 0
#define doMatch4OrNot    0
#define doMatchHLeps     0
#define doMatchAllLeps   0
#define doMatchWHZHttH   1 // plots shown on September 8th
#define do2DPlots        0 // plots shown on September 8th
#define doBaskets        0 // plots shown on September 23rd

#define nSamples 8
#define nHiggsSamples 5
#define nChannels 4
#define nVariables 16
#define nBaskets 8

#define nMatchHLepsStatuses 6
#define nMatchAllLepsStatuses 5
#define nMatchWHStatuses 5
#define nMatchZHStatuses 7
#define nMatchttHStatuses 9

#define nAssocWDecays 2
#define nAssocZDecays 3
#define nAssocttDecays 4
#define nAssocDecays (nAssocWDecays + nAssocZDecays + nAssocttDecays)




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////////////////////// MISCELLANEOUS FUNCTIONS ////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


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

string percentage(float frac) {
  return string(Form("%.1f",100*frac));
}

void SaveCanvas(string directory, TCanvas* c, const char* tag = "") {
  c->SaveAs(Form("%s%s_%s.root",directory.c_str(),c->GetName(),tag));
  c->SaveAs(Form("%s%s_%s.pdf" ,directory.c_str(),c->GetName(),tag));  
  c->SaveAs(Form("%s%s_%s.png" ,directory.c_str(),c->GetName(),tag));
}

void set_plot_style() {
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

void prepareBasketlabels(TH1F* h, string* basketLabel){
  for(int i=1; i<=h->GetNbinsX(); i++)  h->GetXaxis()->SetBinLabel(i,basketLabel[i-1].c_str());
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetXaxis()->LabelsOption("h");
  h->GetXaxis()->SetNdivisions(-h->GetNbinsX());
}

void prepareBasketlabelsHorizontal(TH2F* h2, string* basketLabel){
  int nbins = h2->GetNbinsY();
  for(int i=1; i<=nbins; i++)  h2->GetYaxis()->SetBinLabel(nbins+1-i,basketLabel[i-1].c_str());
  h2->GetYaxis()->SetLabelSize(0.04);
  h2->GetYaxis()->LabelsOption("h");
  h2->GetYaxis()->SetNdivisions(-h2->GetNbinsY());
}

void DrawHorizontal(TH1 *h, string* basketLabel, Bool_t same = false, Int_t nDiv = 0, Bool_t without1stBin = false) {
  Double_t ymin = 0.; //h->GetMinimum();
  Double_t ymax = h->GetMaximum(); // 1.1*h->GetMaximum(); 
  TAxis *axis   = h->GetXaxis();
  Double_t xmin = axis->GetXmin();
  Double_t xmax = axis->GetXmax();
  Int_t nbins   = axis->GetNbins();
  TH2F *h2 = new TH2F("h2",Form("%s;%s;%s",h->GetTitle(),h->GetYaxis()->GetTitle(),""),10,ymin,ymax,nbins,xmin,xmax);
  h2->SetBit(kCanDelete);
  //h2->SetDirectory(0);
  h2->SetTickLength(0.01);
  h2->GetXaxis()->SetLabelSize(0.03);
  if(without1stBin) h2->GetYaxis()->SetRangeUser(0,nBaskets-1);
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




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////// PLOTTING FUNCTIONS ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void DrawByProdmodes(TCanvas* c, TH1F* hSource[nSamples][nChannels], string* prodName, Bool_t* isPresent, Color_t* colors, Bool_t mergeZZBkgd, Bool_t logY = false) {

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
    h[p] = (TH1F*)hSource[p][3]->Clone(); // includes the 3 decay channels !!
    if(mergeZZBkgd&&(p==6||p==7)) h[5]->Add(h[p]);
  }

  Float_t max = 0.;
  for(int p=0; p<nSamples; p++){
    if(!isPresent[p]) continue;
    if(mergeZZBkgd&&(p==6||p==7)) continue;

    h[p]->Scale(1/h[p]->Integral());
    Float_t maxtemp = h[p]->GetMaximum();
    if(maxtemp>max) max = maxtemp;

  }

  for(int p=0; p<nSamples; p++){
    if(!isPresent[p]) continue;
    if(mergeZZBkgd&&(p==6||p==7)) continue;
    string pn = prodName[p];

    h[p]->SetLineColor(colors[p]);
    h[p]->SetLineWidth(2);
    h[p]->SetFillStyle(0);
    h[p]->SetStats(0);
    //if(!logY) h[p]->SetMinimum(0);

    if(p==0){
      h[p]->SetMaximum(1.1*max);
      h[p]->GetYaxis()->SetTitle("normalized to 1");
      h[p]->Draw(); 
    }else{
      h[p]->Draw("sames");
    }

    lgd->AddEntry( h[p], (mergeZZBkgd&&p==5)?"qqtoZZ":pn.c_str(), "l" );
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
  string percentages[nStatuses];
  for(int m=0; m<nStatuses; m++){
    if(excludeH2l2X && matchKeys[m].find("2l2X")!=string::npos) continue;
    percentages[m] = percentage(hMatch[m]->Integral()/denom);
  }

  TH1F* hStacks[nStatuses];
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
  set_plot_style();
  c->cd();
  h->SetTitle(title.c_str());
  h->SetStats(0);
  h->Draw("COLZ"); 
}

void DrawBasketEfficiencies(TCanvas* c, TH1F* h1, string pn, TH1F** hAssocDecay, string* assocDecayName, string* basketLabel, Bool_t logY = false) {

  bool without1stBin = true;

  gStyle->SetOptTitle(1);
  gStyle->SetTitleX(0.25);

  c->cd();
  if(logY) c->SetLogy();
  //c->SetTicks(0,0);
  c->SetLeftMargin(0.25);

  bool doAssoc = (pn=="WH"||pn=="ZH"||pn=="ttH");
  Int_t nDecays = pn=="WH" ? nAssocWDecays : pn=="ZH" ? nAssocZDecays : pn=="ttH" ? nAssocttDecays : 0;
  Int_t startAt = pn=="WH" ? 0 : pn=="ZH" ? nAssocWDecays : pn=="ttH" ? nAssocWDecays+nAssocZDecays : 0;

  TH1F* h = (TH1F*)h1->Clone();
  if(doAssoc) h->SetLineColor(1);
  if(doAssoc) h->SetFillColor(1);
  h->SetTitle(pn.c_str()); 
  h->SetStats(0);
  //if(!logY) h->SetMinimum(0);
  //h->GetXaxis()->SetRangeUser(1,nBaskets);
  //h->Draw(); 
  DrawHorizontal(h,basketLabel,false,508,without1stBin);

  TH1F* hStacks[nDecays];
  for(int a=0; a<nDecays; a++){
    hStacks[a] = (TH1F*)hAssocDecay[startAt+a]->Clone();
    if(a>0) hStacks[a]->Add(hStacks[a-1]);
    hStacks[a]->SetStats(0);
  }
  for(int a=nDecays-1; a>=0; a--){
    //hStacks[a]->Draw("same");
    DrawHorizontal(hStacks[a],basketLabel,true,0,without1stBin);
  }
  
  if(doAssoc){
    TLegend* lgd = new TLegend(0.6,0.12,0.88,(pn=="WH"? 0.24 : pn=="ZH"? 0.3 : pn=="ttH"? 0.36 : 0.));
    //lgd->SetFillColor(0);
    lgd->SetFillStyle(0);
    lgd->SetBorderSize(0);
    for(int a=0; a<nDecays; a++){
      lgd->AddEntry(hStacks[a],assocDecayName[startAt+a].c_str(),"f");
    }
    lgd->Draw();
  }

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

void DrawBasketPurities(TCanvas* c, TH1F** h, string* prodName, TH1F** hAssocDecay, string* assocDecayName, Bool_t withAssocDecays, Bool_t H2l2XAsBkgd, string* basketLabel) {

  gStyle->SetOptTitle(0);
  //gStyle->SetFrameLineWidth(2);

  c->cd();
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.35);

  vector<TH1F*> histos;
  
  TLegend* lgd = new TLegend(0.67,withAssocDecays?(H2l2XAsBkgd?0.45:0.35):0.65,0.95,0.9);
  lgd->SetFillColor(0);
  lgd->SetBorderSize(0);
  
  TPaveText* pavNExp = new TPaveText(0.15,0.1,0.4,0.9,"brNDC");
  pavNExp->SetFillStyle(0);
  pavNExp->SetBorderSize(0);
  pavNExp->SetTextColor(0);
  TH1F* hNExpectedEvts = (TH1F*)h[0]->Clone();

  for(int p=0; p<nHiggsSamples; p++){
    string pn = prodName[p];
    bool doAssoc = (pn=="WH"||pn=="ZH"||pn=="ttH");
    Int_t nDecays = pn=="WH" ? nAssocWDecays : pn=="ZH" ? nAssocZDecays : pn=="ttH" ? nAssocttDecays : 0;
    Int_t startAt = pn=="WH" ? 0 : pn=="ZH" ? nAssocWDecays : pn=="ttH" ? nAssocWDecays+nAssocZDecays : 0;

    TH1F* hTemp = (TH1F*)h[p]->Clone();
    if(H2l2XAsBkgd){
      if(pn=="ZH") hTemp->Add(hAssocDecay[nAssocWDecays+nAssocZDecays-1],-1);
      if(pn=="ttH") hTemp->Add(hAssocDecay[nAssocWDecays+nAssocZDecays+nAssocttDecays-1],-1);
    }

    if(withAssocDecays && doAssoc){
      for(int a=0; a<nDecays; a++){
	bool isH2l2X = ((pn=="ZH"||pn=="ttH") && a==nDecays-1); 
	if(!(H2l2XAsBkgd && isH2l2X)){
	  histos.push_back(hAssocDecay[startAt+a]);
	  lgd->AddEntry(hAssocDecay[startAt+a],assocDecayName[startAt+a].c_str(),"f");
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
    pavNExp->AddText(Form("%.2f exp. events in %.0f fb^{-1}",hNExpectedEvts->GetBinContent(b),LUMI));
  }
  pavNExp->Draw();
  
  gPad->RedrawAxis();

}

void DrawBasketSOSPB(TCanvas* c, TH1F* h, string* basketLabel) {

  gStyle->SetOptTitle(0);

  c->cd();
  c->SetLeftMargin(0.25);

  h->SetFillColor(kBlue+2);
  h->GetYaxis()->SetRangeUser(0.,1.);
  h->GetYaxis()->SetTitle("S/(S+B)");
  DrawHorizontal(h,basketLabel,false,-210);

}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// MAIN MACRO //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void plotProductionModeComparison(

				  // directory to save the plots
				  string outDir,

				  // path of the input primary tree (just let "" for trees you don't have)
				  string file_ggH,
				  string file_VBF     = "",
				  string file_WH      = "",
				  string file_ZH      = "",
				  string file_ttH     = "",
				  string file_ZZ4mu   = "",
				  string file_ZZ4e    = "",
				  string file_ZZ2e2mu = ""

				  )
{

  gSystem->Exec("mkdir -p "+(TString)outDir.c_str());

  ofstream txtOut;
  TString txtOutName = (TString)outDir.c_str()+"/matchingInfo.txt";
  txtOut.open(txtOutName);


  // ------------------------------------------------------------
  // ---------------------- Definitions -------------------------
  // ------------------------------------------------------------

  string prodName[nSamples] = {
    "ggH",
    "VBF",
    "WH",
    "ZH",
    "ttH",
    "ZZ4mu",
    "ZZ4e",
    "ZZ2e2mu",
  };
  Bool_t isPresent[nSamples] = {
    file_ggH != "",
    file_VBF != "",
    file_WH  != "",
    file_ZH  != "",
    file_ttH != "",
    file_ZZ4mu   != "", 
    file_ZZ4e    != "",
    file_ZZ2e2mu != "",
  };
  Bool_t allSignalsArePresent = (file_ggH!="" && file_VBF!="" && file_WH!="" && file_ZH!="" && file_ttH!="");
  Bool_t allZZBkgdsArePresent = (file_ZZ4mu!="" && file_ZZ4e!="" && file_ZZ2e2mu!="");
  Color_t colors[nSamples] = { kBlue, kGreen+2, kRed, kOrange+1, kMagenta, kRed+4, kRed+3, kRed-2 };

  Color_t allColorsMP[nMatchHLepsStatuses][nHiggsSamples] = {
    { kBlue+2 , kGreen+3, kRed+2 , kOrange+4, kMagenta+2  },
    { kBlue-4 , kGreen+2, kRed-4 , kOrange+9, kMagenta    },
    { kBlue-7 , kGreen+1, kRed-7 , kOrange+7, kMagenta-7  },
    { kBlue-9 , kGreen  , kRed-9 , kOrange+1, kMagenta-9  },
    { kBlue-10, kGreen-9, kRed-10, kOrange-9, kMagenta-10 },
    { kGray+1 , kGray+1 , kGray+1, kGray+1  , kGray+1     },
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

  string channels[nChannels] = {
    "4mu",
    "4e",
    "2e2mu",
    "4mu + 4e + 2e2mu",
  };
  string channelsShort[nChannels] = {
    "4mu",
    "4e",
    "2e2mu",
    "Sumchans",
  };
  string subdir[nChannels] = {
    "ZZ4muTree",
    "ZZ4eTree",
    "ZZ2e2muTree",
    "",
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
  string assocDecayName[nAssocDecays] = {
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
  string assocDecayName3[nAssocDecays] = {
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
    kMagenta, 
    kMagenta+1, 
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
    "NGenLepInAcc",
    "NGenLepNotInAcc",
    "NGenHLepNotInAcc",
    "NGenAssocLepNotInAcc",
    "NGenLepMinusNGoodLep",
    "NGenLepInAccMinusNGoodLep",
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

  string basketLabel[nBaskets+1] = {
    "all",
    "#splitline{0/1 jet}{4 leptons}",
    //"#splitline{0/1 jet, 5 leptons}{E_{T}^{miss}>45, p_{T}^{l5}>20}",
    "#splitline{0/1 jet, 5 leptons}{E_{T}^{miss}>45}",
    "#splitline{0/1 jet, 5 leptons}{(else)}",
    "#splitline{0/1 jet}{#geq6 leptons}",
    "#splitline{#geq2 jets, 4 leptons}{D_{jet}>0.5}",
    "#splitline{#geq2 jets, 4 leptons}{(else)}",
    "#splitline{#geq2 jets}{5 leptons}",
    "#splitline{#geq2 jets}{#geq6 leptons}",
  };

//   string basketLabel[nBaskets+1] = {
//     "all",
//     "0/1 jet, 0 extra lepton",
//     "0/1 jet, 1 extra lepton",
//     "0/1 jet, #geq2 extra leptons",
//     "#geq2 jets, 0 extra lepton",
//     "#geq2 jets, 1 extra lepton",
//     "#geq2 jets, #geq2 extra leptons",
//   };

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

  TChain* chain[nSamples][nChannels];
  for(int p=0; p<nSamples; p++){
    if(!isPresent[p]) continue;
    for(int c=0; c<nChannels-1; c++){
      chain[p][c] = new TChain(Form("%s/candTree",subdir[c].c_str()));
      if(prodName[p]=="ggH") chain[p][c]->Add(file_ggH.c_str());
      if(prodName[p]=="VBF") chain[p][c]->Add(file_VBF.c_str());
      if(prodName[p]=="WH" ) chain[p][c]->Add(file_WH .c_str());
      if(prodName[p]=="ZH" ) chain[p][c]->Add(file_ZH .c_str());
      if(prodName[p]=="ttH") chain[p][c]->Add(file_ttH.c_str());
      if(prodName[p]=="ZZ4mu"  ) chain[p][c]->Add(file_ZZ4mu.c_str());
      if(prodName[p]=="ZZ4e"   ) chain[p][c]->Add(file_ZZ4e.c_str());
      if(prodName[p]=="ZZ2e2mu") chain[p][c]->Add(file_ZZ2e2mu.c_str());
    }
  }
  Long64_t NGenEvt[nSamples];
  for(int p=0; p<nSamples; p++){
    if(!isPresent[p]) continue;
    TFile* fTemp = 0;
    if(prodName[p]=="ggH") fTemp = TFile::Open(file_ggH.c_str());
    if(prodName[p]=="VBF") fTemp = TFile::Open(file_VBF.c_str());
    if(prodName[p]=="WH" ) fTemp = TFile::Open(file_WH .c_str());
    if(prodName[p]=="ZH" ) fTemp = TFile::Open(file_ZH .c_str());
    if(prodName[p]=="ttH") fTemp = TFile::Open(file_ttH.c_str());
    if(prodName[p]=="ZZ4mu"  ) fTemp = TFile::Open(file_ZZ4mu  .c_str());
    if(prodName[p]=="ZZ4e"   ) fTemp = TFile::Open(file_ZZ4e   .c_str());
    if(prodName[p]=="ZZ2e2mu") fTemp = TFile::Open(file_ZZ2e2mu.c_str());
    TH1F* hCounters = (TH1F*)fTemp->Get("ZZ4muTree/Counters");
    NGenEvt[p] = (Long64_t)hCounters->GetBinContent(1);
    cout<<"Number of generated events for "+prodName[p]+" is "<<NGenEvt[p]<<endl;
  }
  cout<<endl;

  Float_t xSec[nSamples] = { // X-sec (pb) * BR
    0.0053185200 ,
    0.0004355280 ,
    0.0186014400 * 0.010384 ,
    0.0109639200 * 0.028544 ,
    0.0034135200 * 0.029658 ,
    0.07691 ,
    0.07691 ,
    0.1767 ,
  };
  Float_t weights[nSamples];
  for(int p=0; p<nSamples; p++){
    if(!isPresent[p]) continue;
    weights[p] = LUMI * 1000 * xSec[p] / NGenEvt[p] ;
  }

  Int_t nRun;
  Long64_t nEvent;
  Int_t nLumi;
  Float_t PFMET;
  Int_t iBC;
  Int_t Nvtx;
  Int_t NObsInt;
  Float_t NTrueInt;
  vector<Int_t> *ZZsel = 0;
  vector<Float_t> *ZZMass = 0;
  vector<Float_t> *ZZPt = 0;
  vector<Float_t> *ZZFisher = 0;
  vector<Float_t> *p0plus_VAJHU = 0;
  vector<Float_t> *bkg_VAMCFM = 0;
  vector<Float_t> *Z1Mass = 0;
  vector<Float_t> *Z2Mass = 0;
  vector<Float_t> *CandLep1Pt = 0;
  vector<Float_t> *CandLep1Eta = 0;
  vector<Float_t> *CandLep1Phi = 0;
  vector<Int_t> *CandLep1LepId = 0;
  vector<Float_t> *CandLep2Pt = 0;
  vector<Float_t> *CandLep2Eta = 0;
  vector<Float_t> *CandLep2Phi = 0;
  vector<Int_t> *CandLep2LepId = 0;
  vector<Float_t> *CandLep3Pt = 0;
  vector<Float_t> *CandLep3Eta = 0;
  vector<Float_t> *CandLep3Phi = 0;
  vector<Int_t> *CandLep3LepId = 0;
  vector<Float_t> *CandLep4Pt = 0;
  vector<Float_t> *CandLep4Eta = 0;
  vector<Float_t> *CandLep4Phi = 0;
  vector<Int_t> *CandLep4LepId = 0;
  vector<Float_t> *ExtraLep1Pt = 0;
  vector<Float_t> *ExtraLep1Eta = 0;
  vector<Float_t> *ExtraLep1Phi = 0;
  vector<Int_t> *ExtraLep1LepId = 0;
  vector<Float_t> *ExtraLep2Pt = 0;
  vector<Float_t> *ExtraLep2Eta = 0;
  vector<Float_t> *ExtraLep2Phi = 0;
  vector<Int_t> *ExtraLep2LepId = 0;
  vector<Float_t> *ExtraLep3Pt = 0;
  vector<Float_t> *ExtraLep3Eta = 0;
  vector<Float_t> *ExtraLep3Phi = 0;
  vector<Int_t> *ExtraLep3LepId = 0;
  vector<Int_t> *nExtraLep = 0;
  vector<Int_t> *nExtraZ = 0;
  vector<Int_t> *nJets = 0;
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

  map<string,vector<string> > overlapMapStored[nSamples];
  map<string,vector<string> > overlapMapWithBC[nSamples];
  map<string,vector<string> > overlapMapWithBCFullSel70[nSamples];
  map<string,vector<string> > overlapMapWithBCFullSel100[nSamples];

  Int_t nbStored[nSamples][nChannels];
  Int_t nbWithBC[nSamples][nChannels];
  Int_t nbWithBCFullSel70[nSamples][nChannels];
  Int_t nbWithBCFullSel100[nSamples][nChannels];
  Int_t nbWithBCFullSel100HLepsAreInAcc[nSamples][nChannels];
  Int_t nbWithBCFullSel100HLepsAreGood[nSamples][nChannels];
  Int_t nbWithBCFullSel100All4LepRight[nSamples][nChannels];
  Int_t nbWithBCFullSel100MatchHLeps[nMatchHLepsStatuses][nSamples][nChannels];
  Int_t nbWithBCFullSel100MatchAllLeps[nMatchAllLepsStatuses][nSamples][nChannels];
  Int_t nbWithBCFullSel100MatchWH[nMatchWHStatuses][nChannels];
  Int_t nbWithBCFullSel100MatchZH[nMatchZHStatuses][nChannels];
  Int_t nbWithBCFullSel100MatchttH[nMatchttHStatuses][nChannels];

  TH1F* hBCFullSel100[nVariables][nSamples][nChannels];
  TH1F* hBCFullSel100MatchHLeps[nVariables][nMatchHLepsStatuses][nSamples][nChannels];
  TH1F* hBCFullSel100MatchAllLeps[nVariables][nMatchAllLepsStatuses][nSamples][nChannels];
  TH1F* hBCFullSel100MatchWH[nVariables][nMatchWHStatuses][nChannels];
  TH1F* hBCFullSel100MatchZH[nVariables][nMatchZHStatuses][nChannels];
  TH1F* hBCFullSel100MatchttH[nVariables][nMatchttHStatuses][nChannels];

  TH2F* h2DBCFullSel100[n2DHist][nSamples][nChannels];
  TH2F* h2DBCFullSel100Decays[n2DHist][nSamples][nDecays][nChannels];

  TH1F* hBCFullSel100Baskets[nSamples][nChannels];
  TH1F* hBCFullSel100BasketsAssocDecays[nAssocDecays][nChannels];
  TH1F* hBCFullSel100MasswindowBaskets[nSamples][nChannels];
  TH1F* hBCFullSel100MasswindowBasketsAssocDecays[nAssocDecays][nChannels];

  for(int a=0; a<nAssocDecays; a++){
    for(int c=0; c<nChannels; c++){
      hBCFullSel100BasketsAssocDecays[a][c] = new TH1F(Form("hBCFullSel100_baskets_assocDecays_%i_%s",a,channelsShort[c].c_str()),Form(";;# exp. events in %.0f fb^{-1}",LUMI),nBaskets+1,0,nBaskets+1);
      prepareBasketlabels(hBCFullSel100BasketsAssocDecays[a][c], basketLabel);
      hBCFullSel100BasketsAssocDecays[a][c]->SetLineColor(assocDecayColor[a]);
      hBCFullSel100BasketsAssocDecays[a][c]->SetFillColor(assocDecayColor[a]);
      hBCFullSel100MasswindowBasketsAssocDecays[a][c] = new TH1F(Form("hBCFullSel100Masswindow_baskets_assocDecays_%i_%s",a,channelsShort[c].c_str()),Form(";;# exp. events in %.0f fb^{-1}",LUMI),nBaskets+1,0,nBaskets+1);
      prepareBasketlabels(hBCFullSel100MasswindowBasketsAssocDecays[a][c], basketLabel);
      hBCFullSel100MasswindowBasketsAssocDecays[a][c]->SetLineColor(assocDecayColor[a]);
      hBCFullSel100MasswindowBasketsAssocDecays[a][c]->SetFillColor(assocDecayColor[a]);
    }
  }

  Int_t nbTotalWH[nAssocWDecays] = {0,0};
  Int_t nbHLepsAreInAccWH[nAssocWDecays] = {0,0};
  Int_t nbHLepsAreGoodWH[nAssocWDecays] = {0,0};
  Int_t nbAll4LepRightWH[nAssocWDecays] = {0,0};
  Int_t nbZ1DaughtersFromHWH[nAssocWDecays][4]; for(int i=0; i<nAssocWDecays; i++) for(int j=0; j<4; j++) nbZ1DaughtersFromHWH[i][j] = 0;
  Int_t nbZ2DaughtersFromHWH[nAssocWDecays][4]; for(int i=0; i<nAssocWDecays; i++) for(int j=0; j<4; j++) nbZ2DaughtersFromHWH[i][j] = 0;

  Int_t nbTotalZH[nAssocZDecays] = {0,0,0};
  Int_t nbHLepsAreInAccZH[nAssocZDecays] = {0,0,0};
  Int_t nbHLepsAreGoodZH[nAssocZDecays] = {0,0,0};
  Int_t nbAll4LepRightZH[nAssocZDecays] = {0,0,0};
  Int_t nbZ1DaughtersFromHZH[nAssocZDecays][4]; for(int i=0; i<nAssocZDecays; i++) for(int j=0; j<4; j++) nbZ1DaughtersFromHZH[i][j] = 0;
  Int_t nbZ2DaughtersFromHZH[nAssocZDecays][4]; for(int i=0; i<nAssocZDecays; i++) for(int j=0; j<4; j++) nbZ2DaughtersFromHZH[i][j] = 0;

  Int_t nbTotalttH[nAssocttDecays] = {0,0,0,0};
  Int_t nbHLepsAreInAccttH[nAssocttDecays] = {0,0,0,0};
  Int_t nbHLepsAreGoodttH[nAssocttDecays] = {0,0,0,0};
  Int_t nbAll4LepRightttH[nAssocttDecays] = {0,0,0,0};
  Int_t nbZ1DaughtersFromHttH[nAssocttDecays][4]; for(int i=0; i<nAssocttDecays; i++) for(int j=0; j<4; j++) nbZ1DaughtersFromHttH[i][j] = 0;
  Int_t nbZ2DaughtersFromHttH[nAssocttDecays][4]; for(int i=0; i<nAssocttDecays; i++) for(int j=0; j<4; j++) nbZ2DaughtersFromHttH[i][j] = 0;



  // ------------------------------------------------------------
  // ---------------------- Processing --------------------------
  // ------------------------------------------------------------

  for(int p=0; p<nSamples; p++){
    if(!isPresent[p]) continue;

    txtOut<<prodName[p]<<endl;

    for(int c=0; c<nChannels; c++){
      
      string suffix = "_"+prodName[p]+"_"+channelsShort[c];
      string title = prodName[p]+", channel "+channels[c];

      nbStored[p][c] = 0;
      nbWithBC[p][c] = 0;
      nbWithBCFullSel70[p][c] = 0;
      nbWithBCFullSel100[p][c] = 0;
      nbWithBCFullSel100HLepsAreInAcc[p][c] = 0;
      nbWithBCFullSel100HLepsAreGood[p][c] = 0;
      nbWithBCFullSel100All4LepRight[p][c] = 0;
      for(int m=0; m<nMatchHLepsStatuses; m++) nbWithBCFullSel100MatchHLeps[m][p][c] = 0;
      for(int m=0; m<nMatchAllLepsStatuses; m++) nbWithBCFullSel100MatchAllLeps[m][p][c] = 0;
      if(prodName[p]=="WH") for(int m=0; m<nMatchWHStatuses; m++) nbWithBCFullSel100MatchWH[m][c] = 0;
      if(prodName[p]=="ZH") for(int m=0; m<nMatchZHStatuses; m++) nbWithBCFullSel100MatchZH[m][c] = 0;
      if(prodName[p]=="ttH") for(int m=0; m<nMatchttHStatuses; m++) nbWithBCFullSel100MatchttH[m][c] = 0;

      for(int v=0; v<nVariables; v++){
	if(!plotThisVar[VARIABLELIST][v]) continue;
	hBCFullSel100[v][p][c] = new TH1F(("hBCFullSel100_"+varName[v]+suffix).c_str(),Form("%s;%s;# exp. events in %.0f fb^{-1}",title.c_str(),varLabel[v].c_str(),LUMI),varNbin[v],varMin[v],varMax[v]);
	for(int m=0; m<nMatchHLepsStatuses; m++)
	  hBCFullSel100MatchHLeps[v][m][p][c] = new TH1F(("hBCFullSel100MatchHLeps_"+varName[v]+"_"+matchHLepsInfix[m]+suffix).c_str(),Form("%s;%s;# exp. events in %.0f fb^{-1}",title.c_str(),varLabel[v].c_str(),LUMI),varNbin[v],varMin[v],varMax[v]);
	for(int m=0; m<nMatchAllLepsStatuses; m++)
	  hBCFullSel100MatchAllLeps[v][m][p][c] = new TH1F(("hBCFullSel100MatchAllLeps_"+varName[v]+"_"+matchAllLepsInfix[m]+suffix).c_str(),Form("%s;%s;# exp. events in %.0f fb^{-1}",title.c_str(),varLabel[v].c_str(),LUMI),varNbin[v],varMin[v],varMax[v]);
	if(prodName[p]=="WH")
	  for(int m=0; m<nMatchWHStatuses; m++)
	    hBCFullSel100MatchWH[v][m][c] = new TH1F(Form("hBCFullSel100MatchWH_%s_%i_%s",varName[v].c_str(),m,channelsShort[c].c_str()),Form("%s;%s;# exp. events in %.0f fb^{-1}",title.c_str(),varLabel[v].c_str(),LUMI),varNbin[v],varMin[v],varMax[v]);
	if(prodName[p]=="ZH")
	  for(int m=0; m<nMatchZHStatuses; m++)
	    hBCFullSel100MatchZH[v][m][c] = new TH1F(Form("hBCFullSel100MatchZH_%s_%i_%s",varName[v].c_str(),m,channelsShort[c].c_str()),Form("%s;%s;# exp. events in %.0f fb^{-1}",title.c_str(),varLabel[v].c_str(),LUMI),varNbin[v],varMin[v],varMax[v]);
	if(prodName[p]=="ttH")
	  for(int m=0; m<nMatchttHStatuses; m++)
	    hBCFullSel100MatchttH[v][m][c] = new TH1F(Form("hBCFullSel100MatchttH_%s_%i_%s",varName[v].c_str(),m,channelsShort[c].c_str()),Form("%s;%s;# exp. events in %.0f fb^{-1}",title.c_str(),varLabel[v].c_str(),LUMI),varNbin[v],varMin[v],varMax[v]);
      }

      for(int v2=0; v2<n2DHist; v2++){
	h2DBCFullSel100[v2][p][c] = new TH2F(("h2DBCFullSel100_"+varName[varXindex[v2]]+"_"+varName[varYindex[v2]]+suffix).c_str(),(title+";"+varLabel[varXindex[v2]]+";"+varLabel[varYindex[v2]]).c_str(),varNbin[varXindex[v2]],varMin[varXindex[v2]],varMax[varXindex[v2]],varNbin[varYindex[v2]],varMin[varYindex[v2]],varMax[varYindex[v2]]);
	for(int d=0; d<nDecays; d++){
	  h2DBCFullSel100Decays[v2][p][d][c] = new TH2F(("h2DBCFullSel100_"+decayInfix[d]+"_"+varName[varXindex[v2]]+"_"+varName[varYindex[v2]]+suffix).c_str(),(title+";"+varLabel[varXindex[v2]]+";"+varLabel[varYindex[v2]]).c_str(),varNbin[varXindex[v2]],varMin[varXindex[v2]],varMax[varXindex[v2]],varNbin[varYindex[v2]],varMin[varYindex[v2]],varMax[varYindex[v2]]);
	}
      }
      
      hBCFullSel100Baskets[p][c] = new TH1F(("hBCFullSel100_baskets_"+suffix).c_str(),Form("%s;;# exp. events in %.0f fb^{-1}",title.c_str(),LUMI),nBaskets+1,0,nBaskets+1);
      prepareBasketlabels(hBCFullSel100Baskets[p][c], basketLabel);
      hBCFullSel100Baskets[p][c]->SetLineColor(colors[p]);
      hBCFullSel100Baskets[p][c]->SetFillColor(colors[p]);
      hBCFullSel100MasswindowBaskets[p][c] = new TH1F(("hBCFullSel100Masswindow_baskets_"+suffix).c_str(),Form("%s;;# exp. events in %.0f fb^{-1}",title.c_str(),LUMI),nBaskets+1,0,nBaskets+1);
      prepareBasketlabels(hBCFullSel100MasswindowBaskets[p][c], basketLabel);
      hBCFullSel100MasswindowBaskets[p][c]->SetLineColor(colors[p]);
      hBCFullSel100MasswindowBaskets[p][c]->SetFillColor(colors[p]);
      
    }

    for(int c=0; c<nChannels-1; c++){

      chain[p][c]->SetBranchAddress("RunNumber", &nRun);
      chain[p][c]->SetBranchAddress("EventNumber", &nEvent);
      chain[p][c]->SetBranchAddress("LumiNumber", &nLumi);
      chain[p][c]->SetBranchAddress("PFMET", &PFMET);
      chain[p][c]->SetBranchAddress("iBC", &iBC);
      chain[p][c]->SetBranchAddress("Nvtx", &Nvtx);
      chain[p][c]->SetBranchAddress("NObsInt", &NObsInt);
      chain[p][c]->SetBranchAddress("NTrueInt", &NTrueInt);
      chain[p][c]->SetBranchAddress("ZZsel", &ZZsel);
      chain[p][c]->SetBranchAddress("ZZMass", &ZZMass);
      chain[p][c]->SetBranchAddress("ZZPt", &ZZPt);
      chain[p][c]->SetBranchAddress("ZZFisher", &ZZFisher);
      chain[p][c]->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU);
      chain[p][c]->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM);
      chain[p][c]->SetBranchAddress("Z1Mass", &Z1Mass);
      chain[p][c]->SetBranchAddress("Z2Mass", &Z2Mass);
      chain[p][c]->SetBranchAddress("Lep1Pt", &CandLep1Pt);
      chain[p][c]->SetBranchAddress("Lep1Eta", &CandLep1Eta);
      chain[p][c]->SetBranchAddress("Lep1Phi", &CandLep1Phi);
      chain[p][c]->SetBranchAddress("Lep1LepId", &CandLep1LepId);
      chain[p][c]->SetBranchAddress("Lep2Pt", &CandLep2Pt);
      chain[p][c]->SetBranchAddress("Lep2Eta", &CandLep2Eta);
      chain[p][c]->SetBranchAddress("Lep2Phi", &CandLep2Phi);
      chain[p][c]->SetBranchAddress("Lep2LepId", &CandLep2LepId);
      chain[p][c]->SetBranchAddress("Lep3Pt", &CandLep3Pt);
      chain[p][c]->SetBranchAddress("Lep3Eta", &CandLep3Eta);
      chain[p][c]->SetBranchAddress("Lep3Phi", &CandLep3Phi);
      chain[p][c]->SetBranchAddress("Lep3LepId", &CandLep3LepId);
      chain[p][c]->SetBranchAddress("Lep4Pt", &CandLep4Pt);
      chain[p][c]->SetBranchAddress("Lep4Eta", &CandLep4Eta);
      chain[p][c]->SetBranchAddress("Lep4Phi", &CandLep4Phi);
      chain[p][c]->SetBranchAddress("Lep4LepId", &CandLep4LepId);
      chain[p][c]->SetBranchAddress("ExtraLep1Pt", &ExtraLep1Pt);
      chain[p][c]->SetBranchAddress("ExtraLep1Eta", &ExtraLep1Eta);
      chain[p][c]->SetBranchAddress("ExtraLep1Phi", &ExtraLep1Phi);
      chain[p][c]->SetBranchAddress("ExtraLep1LepId", &ExtraLep1LepId);
      chain[p][c]->SetBranchAddress("ExtraLep2Pt", &ExtraLep2Pt);
      chain[p][c]->SetBranchAddress("ExtraLep2Eta", &ExtraLep2Eta);
      chain[p][c]->SetBranchAddress("ExtraLep2Phi", &ExtraLep2Phi);
      chain[p][c]->SetBranchAddress("ExtraLep2LepId", &ExtraLep2LepId);
      chain[p][c]->SetBranchAddress("ExtraLep3Pt", &ExtraLep3Pt);
      chain[p][c]->SetBranchAddress("ExtraLep3Eta", &ExtraLep3Eta);
      chain[p][c]->SetBranchAddress("ExtraLep3Phi", &ExtraLep3Phi);
      chain[p][c]->SetBranchAddress("ExtraLep3LepId", &ExtraLep3LepId);
      chain[p][c]->SetBranchAddress("nExtraLep", &nExtraLep);
      chain[p][c]->SetBranchAddress("nExtraZ", &nExtraZ);
      //chain[p][c]->SetBranchAddress("nJets", &nJets);
      //chain[p][c]->SetBranchAddress("nCleanedJets", &nJets);
      chain[p][c]->SetBranchAddress("nCleanedJetsPt30", &nJets);
      chain[p][c]->SetBranchAddress("GenLep1Pt", &GenLep1Pt);
      chain[p][c]->SetBranchAddress("GenLep1Eta", &GenLep1Eta);
      chain[p][c]->SetBranchAddress("GenLep1Phi", &GenLep1Phi);
      chain[p][c]->SetBranchAddress("GenLep1Id", &GenLep1Id);
      chain[p][c]->SetBranchAddress("GenLep2Pt", &GenLep2Pt);
      chain[p][c]->SetBranchAddress("GenLep2Eta", &GenLep2Eta);
      chain[p][c]->SetBranchAddress("GenLep2Phi", &GenLep2Phi);
      chain[p][c]->SetBranchAddress("GenLep2Id", &GenLep2Id);
      chain[p][c]->SetBranchAddress("GenLep3Pt", &GenLep3Pt);
      chain[p][c]->SetBranchAddress("GenLep3Eta", &GenLep3Eta);
      chain[p][c]->SetBranchAddress("GenLep3Phi", &GenLep3Phi);
      chain[p][c]->SetBranchAddress("GenLep3Id", &GenLep3Id);
      chain[p][c]->SetBranchAddress("GenLep4Pt", &GenLep4Pt);
      chain[p][c]->SetBranchAddress("GenLep4Eta", &GenLep4Eta);
      chain[p][c]->SetBranchAddress("GenLep4Phi", &GenLep4Phi);
      chain[p][c]->SetBranchAddress("GenLep4Id", &GenLep4Id);
      chain[p][c]->SetBranchAddress("GenAssocLep1Pt", &GenAssocLep1Pt);
      chain[p][c]->SetBranchAddress("GenAssocLep1Eta", &GenAssocLep1Eta);
      chain[p][c]->SetBranchAddress("GenAssocLep1Phi", &GenAssocLep1Phi);
      chain[p][c]->SetBranchAddress("GenAssocLep1Id", &GenAssocLep1Id);
      chain[p][c]->SetBranchAddress("GenAssocLep2Pt", &GenAssocLep2Pt);
      chain[p][c]->SetBranchAddress("GenAssocLep2Eta", &GenAssocLep2Eta);
      chain[p][c]->SetBranchAddress("GenAssocLep2Phi", &GenAssocLep2Phi);
      chain[p][c]->SetBranchAddress("GenAssocLep2Id", &GenAssocLep2Id);

      Long64_t entries = chain[p][c]->GetEntries();

      for (Long64_t z=0; z<entries; ++z){

	if(DEBUG && z>1000) break;

	printStatus(z,20000,entries,"entries");

	chain[p][c]->GetEntry(z);


	// ---------------------- Gen lepton counting --------------------------
 
	Int_t nGenHLep = 0;
	Int_t nGenHLepInAcc = 0;
	Int_t nGenAssocLep = 0;
	Int_t nGenAssocLepInAcc = 0;
	Int_t   GenHLepId[4] = {GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id};
	Float_t GenHLepPt[4] = {GenLep1Pt,GenLep2Pt,GenLep3Pt,GenLep4Pt};
	Float_t GenHLepEta[4] = {GenLep1Eta,GenLep2Eta,GenLep3Eta,GenLep4Eta};
	Float_t GenHLepPhi[4] = {GenLep1Phi,GenLep2Phi,GenLep3Phi,GenLep4Phi};
	Bool_t  GenHLepIsInAcc[4];
	Int_t   GenAssocLepId[2] = {GenAssocLep1Id,GenAssocLep2Id};
	Float_t GenAssocLepPt[2] = {GenAssocLep1Pt,GenAssocLep2Pt};
	Float_t GenAssocLepEta[2] = {GenAssocLep1Eta,GenAssocLep2Eta};
	Float_t GenAssocLepPhi[2] = {GenAssocLep1Phi,GenAssocLep2Phi};
	Bool_t  GenAssocLepIsInAcc[4];
	for(Int_t iGenHLep=0; iGenHLep<4; iGenHLep++){
	  if(abs(GenHLepId[iGenHLep])==11 || abs(GenHLepId[iGenHLep])==13){
	    nGenHLep++;
	    GenHLepIsInAcc[iGenHLep] = (abs(GenHLepId[iGenHLep])==11 && (acceptanceIncludesPt?GenHLepPt[iGenHLep]>7.:true) && fabs(GenHLepEta[iGenHLep])<2.5) || 
	                               (abs(GenHLepId[iGenHLep])==13 && (acceptanceIncludesPt?GenHLepPt[iGenHLep]>5.:true) && fabs(GenHLepEta[iGenHLep])<2.4) ;
	    if(GenHLepIsInAcc[iGenHLep]) nGenHLepInAcc++;
	  }
	}
	for(Int_t iGenAssocLep=0; iGenAssocLep<2; iGenAssocLep++){
	  if(abs(GenAssocLepId[iGenAssocLep])==11 || abs(GenAssocLepId[iGenAssocLep])==13){
	    nGenAssocLep++;
	    GenAssocLepIsInAcc[iGenAssocLep] = (abs(GenAssocLepId[iGenAssocLep])==11 && (acceptanceIncludesPt?GenAssocLepPt[iGenAssocLep]>7.:true) && fabs(GenAssocLepEta[iGenAssocLep])<2.5) || 
	                                       (abs(GenAssocLepId[iGenAssocLep])==13 && (acceptanceIncludesPt?GenAssocLepPt[iGenAssocLep]>5.:true) && fabs(GenAssocLepEta[iGenAssocLep])<2.4) ;
	    if(GenAssocLepIsInAcc[iGenAssocLep]) nGenAssocLepInAcc++;
	  }
	}
	Int_t nGenLep = nGenHLep + nGenAssocLep;
	Int_t nGenLepInAcc = nGenHLepInAcc + nGenAssocLepInAcc;

	// Don't want to consider such events.
	// (They are sometimes selected because some e/mu from tau decays are reconstructed as good leptons.)
	if(nGenLep<4) continue;

	///////////////////////////
	if(excludeH2l2X)
	  if(nGenHLep!=4) continue;
	///////////////////////////

	Int_t currentAssocDecay = -1;
	if(prodName[p]=="WH"){
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


	// ---------------------- Successive selection steps --------------------------

	string evtID = eventID(nRun,nLumi,nEvent);
	overlapMapStored[p][evtID].push_back(channels[c]);

	nbStored[p][c]++;
	nbStored[p][3]++;

	if(iBC>=0){

	  overlapMapWithBC[p][evtID].push_back(channels[c]);

	  nbWithBC[p][c]++;
	  nbWithBC[p][3]++;

	  Bool_t FullSel70 = ZZsel->at(iBC)>=90;
	  Bool_t FullSel100 = ZZsel->at(iBC)>=100;
	  
	  if(FullSel70){
	    overlapMapWithBCFullSel70[p][evtID].push_back(channels[c]);
	    nbWithBCFullSel70[p][c]++;
	    nbWithBCFullSel70[p][3]++;
	  }
	  
	  if(FullSel100){

	    overlapMapWithBCFullSel100[p][evtID].push_back(channels[c]);

	    nbWithBCFullSel100[p][c]++;
	    nbWithBCFullSel100[p][3]++;


	    // ---------------------- Gen to Good lepton matching --------------------------

	    Float_t CandLepEta[4] = {CandLep1Eta->at(iBC),CandLep2Eta->at(iBC),CandLep3Eta->at(iBC),CandLep4Eta->at(iBC)};
	    Float_t CandLepPhi[4] = {CandLep1Phi->at(iBC),CandLep2Phi->at(iBC),CandLep3Phi->at(iBC),CandLep4Phi->at(iBC)};
	    Float_t ExtraLepEta[3] = {ExtraLep1Eta->at(iBC),ExtraLep2Eta->at(iBC),ExtraLep3Eta->at(iBC)};
	    Float_t ExtraLepPhi[3] = {ExtraLep1Phi->at(iBC),ExtraLep2Phi->at(iBC),ExtraLep3Phi->at(iBC)};
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
		  if(deltaR(GenHLepEta[iGenHLep],GenHLepPhi[iGenHLep],CandLepEta[iCandLep],CandLepPhi[iCandLep]) < 0.1){
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
		for(Int_t iExtraLep=0; iExtraLep<nExtraLep->at(iBC); iExtraLep++){
		  if(deltaR(GenHLepEta[iGenHLep],GenHLepPhi[iGenHLep],ExtraLepEta[iExtraLep],ExtraLepPhi[iExtraLep]) < 0.1){
		    nRecoLepMatchedToGenHLep[iGenHLep]++;
		    nGenLepMatchedToRecoLep[4+iExtraLep]++;
		  }
		}
	      }
	    }
	    for(Int_t iGenAssocLep=0; iGenAssocLep<2; iGenAssocLep++){
	      if(abs(GenAssocLepId[iGenAssocLep])==11 || abs(GenAssocLepId[iGenAssocLep])==13){
		for(Int_t iCandLep=0; iCandLep<4; iCandLep++){
		  if(deltaR(GenAssocLepEta[iGenAssocLep],GenAssocLepPhi[iGenAssocLep],CandLepEta[iCandLep],CandLepPhi[iCandLep]) < 0.1){
		    nRecoLepMatchedToGenAssocLep[iGenAssocLep]++;
		    nCandLepMatchedToGenAssocLep[iGenAssocLep]++;
		    nGenLepMatchedToRecoLep[iCandLep]++;
		    nGenLepMatchedToCandLep[iCandLep]++;
		  }
		}
		for(Int_t iExtraLep=0; iExtraLep<nExtraLep->at(iBC); iExtraLep++){
		  if(deltaR(GenAssocLepEta[iGenAssocLep],GenAssocLepPhi[iGenAssocLep],ExtraLepEta[iExtraLep],ExtraLepPhi[iExtraLep]) < 0.1){
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
	    for(Int_t iRecoLep=0; iRecoLep<4+nExtraLep->at(iBC); iRecoLep++) if(nGenLepMatchedToRecoLep[iRecoLep]>1){ foundMatchingAmbiguity = true; break; }

	    Int_t nOnes = 0;
	    for(Int_t iGenHLep=0; iGenHLep<4; iGenHLep++) if(nRecoLepMatchedToGenHLep[iGenHLep]==1) nOnes++;
	    Bool_t HLepsAreGood = !foundMatchingAmbiguity && nOnes==4;

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
	    if(prodName[p]=="WH"){
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

	    if(requireHLepsAreInAcc)
	      if(nGenHLepInAcc!=4) continue;
	    if(requireExactly4GoodLeps)
	      if(!(nExtraLep->at(iBC)==0)) continue;
	    if(require5GoodLepsOrMore)
	      if(!(nExtraLep->at(iBC)>=1)) continue;
	    if(requireExactly5GoodLeps)
	      if(!(nExtraLep->at(iBC)==1)) continue;
	    if(requireExactly6GoodLeps)
	      if(!(nExtraLep->at(iBC)==2)) continue;
	    if(requireHLepsAreGood)
	      if(!HLepsAreGood) continue;


	    // ---------------------- Fill histograms and increment counters --------------------------

	    if(nGenHLepInAcc==4){
	      nbWithBCFullSel100HLepsAreInAcc[p][c]++;
	      nbWithBCFullSel100HLepsAreInAcc[p][3]++;
	    }
	    if(HLepsAreGood){
	      nbWithBCFullSel100HLepsAreGood[p][c]++;
	      nbWithBCFullSel100HLepsAreGood[p][3]++;
	    }

	    Float_t varVal[nVariables] = {
	      ZZMass->at(iBC),
	      Z1Mass->at(iBC),
	      Z2Mass->at(iBC),
	      p0plus_VAJHU->at(iBC) / ( p0plus_VAJHU->at(iBC) + bkg_VAMCFM->at(iBC) ),
	      ZZFisher->at(iBC),
	      ZZPt->at(iBC),
	      nGenLep,
	      nGenLepInAcc,
	      nGenLep-nGenLepInAcc,
	      nGenHLep-nGenHLepInAcc,
	      nGenAssocLep-nGenAssocLepInAcc,
	      nGenLep-(4+nExtraLep->at(iBC)),
	      nGenLepInAcc-(4+nExtraLep->at(iBC)),
	      nExtraLep->at(iBC),
	      nExtraZ->at(iBC),
	      nJets->at(iBC),
	    };

	    for(int v=0; v<nVariables; v++){
	      if(!plotThisVar[VARIABLELIST][v]) continue;
	      hBCFullSel100[v][p][c]->Fill(varVal[v],weights[p]);
	      hBCFullSel100[v][p][3]->Fill(varVal[v],weights[p]);
	    }

	    for(int v2=0; v2<n2DHist; v2++){
	      h2DBCFullSel100[v2][p][c]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],weights[p]);
	      h2DBCFullSel100[v2][p][3]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],weights[p]);
	      h2DBCFullSel100Decays[v2][p][0][c]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],weights[p]);
	      h2DBCFullSel100Decays[v2][p][0][3]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],weights[p]);
	      if(nGenHLep==4){
		h2DBCFullSel100Decays[v2][p][1][c]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],weights[p]);
		h2DBCFullSel100Decays[v2][p][1][3]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],weights[p]);
		if(currentMatchHLepsStatus==0){
		  h2DBCFullSel100Decays[v2][p][3][c]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],weights[p]);
		  h2DBCFullSel100Decays[v2][p][3][3]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],weights[p]);
		}else if(currentMatchHLepsStatus<5){
		  h2DBCFullSel100Decays[v2][p][4][c]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],weights[p]);
		  h2DBCFullSel100Decays[v2][p][4][3]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],weights[p]);
		}
	      }else{
		h2DBCFullSel100Decays[v2][p][2][c]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],weights[p]);
		h2DBCFullSel100Decays[v2][p][2][3]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]],weights[p]);
	      }
	    }

	    nbWithBCFullSel100MatchHLeps[currentMatchHLepsStatus][p][c]++;
	    nbWithBCFullSel100MatchHLeps[currentMatchHLepsStatus][p][3]++;
	    for(int v=0; v<nVariables; v++){
	      if(!plotThisVar[VARIABLELIST][v]) continue;
	      hBCFullSel100MatchHLeps[v][currentMatchHLepsStatus][p][c]->Fill(varVal[v],weights[p]);
	      hBCFullSel100MatchHLeps[v][currentMatchHLepsStatus][p][3]->Fill(varVal[v],weights[p]);
	    }

	    nbWithBCFullSel100MatchAllLeps[currentMatchAllLepsStatus][p][c]++;
	    nbWithBCFullSel100MatchAllLeps[currentMatchAllLepsStatus][p][3]++;
	    for(int v=0; v<nVariables; v++){
	      if(!plotThisVar[VARIABLELIST][v]) continue;
	      hBCFullSel100MatchAllLeps[v][currentMatchAllLepsStatus][p][c]->Fill(varVal[v],weights[p]);
	      hBCFullSel100MatchAllLeps[v][currentMatchAllLepsStatus][p][3]->Fill(varVal[v],weights[p]);
	    }
	    
	    if(prodName[p]=="WH"){
	      nbWithBCFullSel100MatchWH[currentMatchWHStatus][c]++;
	      nbWithBCFullSel100MatchWH[currentMatchWHStatus][3]++;
	      for(int v=0; v<nVariables; v++){
		if(!plotThisVar[VARIABLELIST][v]) continue;
		hBCFullSel100MatchWH[v][currentMatchWHStatus][c]->Fill(varVal[v],weights[p]);
		hBCFullSel100MatchWH[v][currentMatchWHStatus][3]->Fill(varVal[v],weights[p]);
	      }
	    }
	    if(prodName[p]=="ZH"){
	      nbWithBCFullSel100MatchZH[currentMatchZHStatus][c]++;
	      nbWithBCFullSel100MatchZH[currentMatchZHStatus][3]++;
	      for(int v=0; v<nVariables; v++){
		if(!plotThisVar[VARIABLELIST][v]) continue;
		hBCFullSel100MatchZH[v][currentMatchZHStatus][c]->Fill(varVal[v],weights[p]);
		hBCFullSel100MatchZH[v][currentMatchZHStatus][3]->Fill(varVal[v],weights[p]);
	      }
	    }
	    if(prodName[p]=="ttH"){
	      nbWithBCFullSel100MatchttH[currentMatchttHStatus][c]++;
	      nbWithBCFullSel100MatchttH[currentMatchttHStatus][3]++;
	      for(int v=0; v<nVariables; v++){
		if(!plotThisVar[VARIABLELIST][v]) continue;
		hBCFullSel100MatchttH[v][currentMatchttHStatus][c]->Fill(varVal[v],weights[p]);
		hBCFullSel100MatchttH[v][currentMatchttHStatus][3]->Fill(varVal[v],weights[p]);
	      }
	    }

	    if(currentMatchAllLepsStatus==0){
	      nbWithBCFullSel100All4LepRight[p][c]++;
	      nbWithBCFullSel100All4LepRight[p][3]++;
	    }

	    if(prodName[p]=="WH"){
	      Int_t currentWDecay = -1;
	      if(nGenHLep==4 && nGenAssocLep==0) currentWDecay = 0;
	      else if(nGenHLep==4 && nGenAssocLep==1) currentWDecay = 1;
	      else cout<<"error"<<endl;
	      nbTotalWH[currentWDecay]++;
	      if(nGenHLepInAcc==4) nbHLepsAreInAccWH[currentWDecay]++;
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
	      if(nGenHLepInAcc==4) nbHLepsAreInAccZH[currentZDecay]++;
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
	      if(nGenHLepInAcc==4) nbHLepsAreInAccttH[currentttDecay]++;
	      if(HLepsAreGood) nbHLepsAreGoodttH[currentttDecay]++;
	      if(currentMatchAllLepsStatus==0) nbAll4LepRightttH[currentttDecay]++;
	      nbZ1DaughtersFromHttH[currentttDecay][currentZ1MatchStatus]++;
	      nbZ2DaughtersFromHttH[currentttDecay][currentZ2MatchStatus]++;
	    }


	    // ---------------------- Categorization stuff --------------------------

	    Int_t currentBasket = -1;
	    if(nJets->at(iBC)==0 || nJets->at(iBC)==1){
	      if(nExtraLep->at(iBC)==0) currentBasket = 1;
	      else if(nExtraLep->at(iBC)==1){ 
		//if(ExtraLep1Pt->at(iBC)>20. && PFMET>45.) currentBasket = 2;
		if(PFMET>45.) currentBasket = 2;
		else currentBasket = 3;
	      }else if(nExtraLep->at(iBC)>=2) currentBasket = 4;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else if(nJets->at(iBC)>=2){
	      if(nExtraLep->at(iBC)==0){
		if(ZZFisher->at(iBC)>0.5) currentBasket = 5;
		else currentBasket = 6;	      
	      }else if(nExtraLep->at(iBC)==1) currentBasket = 7;
	      else if(nExtraLep->at(iBC)>=2) currentBasket = 8;
	      else cout<<"WARNING : inconsistent nExtraLep"<<endl;
	    }else{
	      cout<<"WARNING : inconsistent nJets"<<endl;
	    }

	    hBCFullSel100Baskets[p][c]->Fill(currentBasket,weights[p]);
	    hBCFullSel100Baskets[p][c]->Fill(0.,weights[p]);
	    hBCFullSel100Baskets[p][3]->Fill(currentBasket,weights[p]);
	    hBCFullSel100Baskets[p][3]->Fill(0.,weights[p]);
	    if(currentAssocDecay>=0){
	      hBCFullSel100BasketsAssocDecays[currentAssocDecay][c]->Fill(currentBasket,weights[p]);
	      hBCFullSel100BasketsAssocDecays[currentAssocDecay][c]->Fill(0.,weights[p]);
	      hBCFullSel100BasketsAssocDecays[currentAssocDecay][3]->Fill(currentBasket,weights[p]);
	      hBCFullSel100BasketsAssocDecays[currentAssocDecay][3]->Fill(0.,weights[p]);	      
	    }
	    if(106<ZZMass->at(iBC) && ZZMass->at(iBC)<141){
	      hBCFullSel100MasswindowBaskets[p][c]->Fill(currentBasket,weights[p]);
	      hBCFullSel100MasswindowBaskets[p][c]->Fill(0.,weights[p]);
	      hBCFullSel100MasswindowBaskets[p][3]->Fill(currentBasket,weights[p]);
	      hBCFullSel100MasswindowBaskets[p][3]->Fill(0.,weights[p]);
	      if(currentAssocDecay>=0){
		hBCFullSel100MasswindowBasketsAssocDecays[currentAssocDecay][c]->Fill(currentBasket,weights[p]);
		hBCFullSel100MasswindowBasketsAssocDecays[currentAssocDecay][c]->Fill(0.,weights[p]);
		hBCFullSel100MasswindowBasketsAssocDecays[currentAssocDecay][3]->Fill(currentBasket,weights[p]);
		hBCFullSel100MasswindowBasketsAssocDecays[currentAssocDecay][3]->Fill(0.,weights[p]);	      
	      }
	    }

	
	  } // end if(FullSel100)

	} // end if(iBC>=0)

      } // end for entries

    } // end for channels


    // ---------------------- Printing --------------------------

    Int_t widthColumn2 = 6;

    for(int c=0; c<nChannels; c++){
      txtOut<<" "<<channels[c]<<endl;
      txtOut<<"  stored :            "<<fixWidth(Form("%i",nbStored[p][c]),widthColumn2,false)<<endl;
      txtOut<<"  iBC>0 :             "<<fixWidth(Form("%i",nbWithBC[p][c]),widthColumn2,false)<<endl;
      txtOut<<"  iBC>0, FullSel70 :  "<<fixWidth(Form("%i",nbWithBCFullSel70[p][c]),widthColumn2,false)<<endl;
      txtOut<<"  iBC>0, FullSel100 : "<<fixWidth(Form("%i",nbWithBCFullSel100[p][c]),widthColumn2,false)<<endl;
      if(requireHLepsAreInAcc)    txtOut<<"  [ From this point, require that the 4 gen-leptons from the H are in the acceptance ]"<<endl;
      if(requireExactly4GoodLeps) txtOut<<"  [ From this point, require that there is exactly 4 good leptons ]"<<endl;
      if(require5GoodLepsOrMore)  txtOut<<"  [ From this point, require that there is at least 5 good leptons ]"<<endl;
      if(requireExactly5GoodLeps) txtOut<<"  [ From this point, require that there is exactly 5 good leptons ]"<<endl;
      if(requireExactly6GoodLeps) txtOut<<"  [ From this point, require that there is exactly 6 good leptons ]"<<endl;
      if(requireHLepsAreGood)     txtOut<<"  [ From this point, require that the 4 gen-leptons from the H are reconstructed as good leptons ]"<<endl;
      txtOut<<"  iBC>0, FullSel100, the 4 gen-leptons from the H are in the acceptance : "<<fixWidth(Form("%i",nbWithBCFullSel100HLepsAreInAcc[p][c]),widthColumn2,false)<<endl;
      txtOut<<"  iBC>0, FullSel100, the 4 gen-leptons from the H are reconstructed as good leptons : "<<fixWidth(Form("%i",nbWithBCFullSel100HLepsAreGood[p][c]),widthColumn2,false)<<endl;
      txtOut<<"  iBC>0, FullSel100, the 4 gen-leptons from the H are the 4 good leptons of the best candidate : "<<fixWidth(Form("%i",nbWithBCFullSel100All4LepRight[p][c]),widthColumn2,false)<<endl;
    }

    txtOut<<" "<<"# events :"<<endl;
    txtOut<<"  stored :            "<<fixWidth(Form("%i",(int)overlapMapStored[p].size()),widthColumn2,false)<<endl;
    txtOut<<"  iBC>0 :             "<<fixWidth(Form("%i",(int)overlapMapWithBC[p].size()),widthColumn2,false)<<endl;
    txtOut<<"  iBC>0, FullSel70 :  "<<fixWidth(Form("%i",(int)overlapMapWithBCFullSel70[p].size()),widthColumn2,false)<<endl;
    txtOut<<"  iBC>0, FullSel100 : "<<fixWidth(Form("%i",(int)overlapMapWithBCFullSel100[p].size()),widthColumn2,false)<<endl;

    txtOut<<" "<<"# events stored in 2:3 decay channels :"<<endl;
    Int_t nbOverlapStored2 = 0;
    Int_t nbOverlapStored3 = 0;
    for(map<string,vector<string> >::iterator it=overlapMapStored[p].begin(); it!=overlapMapStored[p].end(); it++){
      if(it->second.size()==2) nbOverlapStored2++;
      if(it->second.size()==3) nbOverlapStored3++;
    }
    txtOut<<"  stored :            "<<nbOverlapStored2<<":"<<nbOverlapStored3<<" ("<<percentage((float)(nbOverlapStored2+nbOverlapStored3)/(int)overlapMapStored[p].size())<<" %)"<<endl;
    Int_t nbOverlapWithBC2 = 0;
    Int_t nbOverlapWithBC3 = 0;
    for(map<string,vector<string> >::iterator it=overlapMapWithBC[p].begin(); it!=overlapMapWithBC[p].end(); it++){
      if(it->second.size()==2) nbOverlapWithBC2++;
      if(it->second.size()==3) nbOverlapWithBC3++;
    }
    txtOut<<"  iBC>0 :             "<<nbOverlapWithBC2<<":"<<nbOverlapWithBC3<<" ("<<percentage((float)(nbOverlapWithBC2+nbOverlapWithBC3)/(int)overlapMapWithBC[p].size())<<" %)"<<endl;
    Int_t nbOverlapWithBCFullSel702 = 0;
    Int_t nbOverlapWithBCFullSel703 = 0;
    for(map<string,vector<string> >::iterator it=overlapMapWithBCFullSel70[p].begin(); it!=overlapMapWithBCFullSel70[p].end(); it++){
      if(it->second.size()==2) nbOverlapWithBCFullSel702++;
      if(it->second.size()==3) nbOverlapWithBCFullSel703++;
    }
    txtOut<<"  iBC>0, FullSel70 :  "<<nbOverlapWithBCFullSel702<<":"<<nbOverlapWithBCFullSel703<<" ("<<percentage((float)(nbOverlapWithBCFullSel702+nbOverlapWithBCFullSel703)/(int)overlapMapWithBCFullSel70[p].size())<<" %)"<<endl;
    Int_t nbOverlapWithBCFullSel1002 = 0;
    Int_t nbOverlapWithBCFullSel1003 = 0;
    for(map<string,vector<string> >::iterator it=overlapMapWithBCFullSel100[p].begin(); it!=overlapMapWithBCFullSel100[p].end(); it++){
      if(it->second.size()==2) nbOverlapWithBCFullSel1002++;
      if(it->second.size()==3) nbOverlapWithBCFullSel1003++;
    }
    txtOut<<"  iBC>0, FullSel100 : "<<nbOverlapWithBCFullSel1002<<":"<<nbOverlapWithBCFullSel1003<<" ("<<percentage((float)(nbOverlapWithBCFullSel1002+nbOverlapWithBCFullSel1003)/(int)overlapMapWithBCFullSel100[p].size())<<" %)"<<endl;

    int widthColumn1 = 22;
    int widthOtherColumns = 7;
    string separator = repeat("-",widthColumn1+4*(2+widthOtherColumns));
    if(prodName[p]=="WH"){
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
      for(int i=0; i<nAssocWDecays; i++) txtOut<<" "<<fixWidth(WHdecays[i],widthColumn1,true)<<"   "<<percentage((float)nbHLepsAreInAccWH[i]/nbTotalWH[i])+" %"<<endl;
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
      for(int i=0; i<nAssocZDecays; i++) txtOut<<" "<<fixWidth(ZHdecays[i],widthColumn1,true)<<"   "<<percentage((float)nbHLepsAreInAccZH[i]/nbTotalZH[i])+" %"<<endl;
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
      for(int i=0; i<nAssocttDecays; i++) txtOut<<" "<<fixWidth(ttHdecays[i],widthColumn1,true)<<"   "<<percentage((float)nbHLepsAreInAccttH[i]/nbTotalttH[i])+" %"<<endl;
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



  // ------------------------------------------------------------
  // -------------------------- Plots ---------------------------
  // ------------------------------------------------------------

  if(doProdComp){
    TCanvas* cBCFullSel100[nVariables];
    for(int v=0; v<nVariables; v++){
      if(!plotThisVar[VARIABLELIST][v]) continue;
      cBCFullSel100[v] = new TCanvas(Form("cBCFullSel100_%s",varName[v].c_str()),Form("cBCFullSel100_%s",varName[v].c_str()),500,500);
      DrawByProdmodes(cBCFullSel100[v],hBCFullSel100[v],prodName,isPresent,colors,allZZBkgdsArePresent,v==0);
      SaveCanvas(outDir,cBCFullSel100[v]);
    }
  }

  if(doProdCompMatch4){
    TCanvas* cBCFullSel100Match4[nVariables];
    for(int v=0; v<nVariables; v++){
      if(!plotThisVar[VARIABLELIST][v]) continue;
      cBCFullSel100Match4[v] = new TCanvas(Form("cBCFullSel100Match4_%s",varName[v].c_str()),Form("cBCFullSel100Match4_%s",varName[v].c_str()),500,500);
      DrawByProdmodes(cBCFullSel100Match4[v],hBCFullSel100MatchHLeps[v][0],prodName,isPresent,colors,allZZBkgdsArePresent,v==0);
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
	DrawMatch4OrNot(cBCFullSel100Match4OrNot[v][p],hBCFullSel100[v][p][3],hBCFullSel100MatchHLeps[v][0][p][3],pn,allColorsMP[1][p],allColorsMP[4][p],v==0);
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
	TH1F* h[nMatchHLepsStatuses]; for(int m=0; m<nMatchHLepsStatuses; m++) h[m] = (TH1F*)hBCFullSel100MatchHLeps[v][m][p][3];
	cBCFullSel100MatchHLeps[v][p] = new TCanvas(Form("cBCFullSel100MatchHLeps_%s_%s",varName[v].c_str(),pn.c_str()),Form("cBCFullSel100MatchHLeps_%s_%s",varName[v].c_str(),pn.c_str()),500,500);
	DrawMatchHLeps(cBCFullSel100MatchHLeps[v][p],hBCFullSel100[v][p][3],h,pn,allColorsPM1[p],matchHLepsKeys,v==0);
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
	TH1F* h[nMatchAllLepsStatuses]; for(int m=0; m<nMatchAllLepsStatuses; m++) h[m] = (TH1F*)hBCFullSel100MatchAllLeps[v][m][p][3];
	cBCFullSel100MatchAllLeps[v][p] = new TCanvas(Form("cBCFullSel100MatchAllLeps_%s_%s",varName[v].c_str(),pn.c_str()),Form("cBCFullSel100MatchAllLeps_%s_%s",varName[v].c_str(),pn.c_str()),500,500);
	DrawMatchAllLeps(cBCFullSel100MatchAllLeps[v][p],hBCFullSel100[v][p][3],h,pn,allColorsPM2[p],matchAllLepsKeys,v==0);
	SaveCanvas(outDir,cBCFullSel100MatchAllLeps[v][p]);
      }
    }
  }

  if(doMatchWHZHttH){

    if(file_WH!=""){
      TCanvas* cBCFullSel100MatchWH[nVariables];
      for(int v=0; v<nVariables; v++){
	if(!plotThisVar[VARIABLELIST][v]) continue;
	TH1F* hMatchWH[nMatchWHStatuses]; for(int m=0; m<nMatchWHStatuses; m++) hMatchWH[m] = (TH1F*)hBCFullSel100MatchWH[v][m][3];
	cBCFullSel100MatchWH[v] = new TCanvas(Form("cBCFullSel100MatchWH_%s",varName[v].c_str()),Form("cBCFullSel100MatchWH_%s",varName[v].c_str()),500,500);
	if(v==1 || v==3)
	  DrawMatchCustom(cBCFullSel100MatchWH[v],hBCFullSel100[v][2][3],hMatchWH,prodName[2],colorsMatchWH,matchWH,nMatchWHStatuses,0.11,0.6,v==0);
	else
	  DrawMatchCustom(cBCFullSel100MatchWH[v],hBCFullSel100[v][2][3],hMatchWH,prodName[2],colorsMatchWH,matchWH,nMatchWHStatuses,0.34,0.6,v==0);
	SaveCanvas(outDir,cBCFullSel100MatchWH[v]);
      }
    }

    if(file_ZH!=""){
      TCanvas* cBCFullSel100MatchZH[nVariables];
      for(int v=0; v<nVariables; v++){
	if(!plotThisVar[VARIABLELIST][v]) continue;
	TH1F* hMatchZH[nMatchZHStatuses]; for(int m=0; m<nMatchZHStatuses; m++) hMatchZH[m] = (TH1F*)hBCFullSel100MatchZH[v][m][3];
	cBCFullSel100MatchZH[v] = new TCanvas(Form("cBCFullSel100MatchZH_%s",varName[v].c_str()),Form("cBCFullSel100MatchZH_%s",varName[v].c_str()),500,500);
	if(v==1 || v==3)
	  DrawMatchCustom(cBCFullSel100MatchZH[v],hBCFullSel100[v][3][3],hMatchZH,prodName[3],colorsMatchZH,matchZH,nMatchZHStatuses,0.11,0.5,v==0);
	else
	  DrawMatchCustom(cBCFullSel100MatchZH[v],hBCFullSel100[v][3][3],hMatchZH,prodName[3],colorsMatchZH,matchZH,nMatchZHStatuses,0.34,0.5,v==0);
	SaveCanvas(outDir,cBCFullSel100MatchZH[v]);
      }
    }

    if(file_ttH!=""){
      TCanvas* cBCFullSel100MatchttH[nVariables];
      for(int v=0; v<nVariables; v++){
	if(!plotThisVar[VARIABLELIST][v]) continue;
	TH1F* hMatchttH[nMatchttHStatuses]; for(int m=0; m<nMatchttHStatuses; m++) hMatchttH[m] = (TH1F*)hBCFullSel100MatchttH[v][m][3];
	cBCFullSel100MatchttH[v] = new TCanvas(Form("cBCFullSel100MatchttH_%s",varName[v].c_str()),Form("cBCFullSel100MatchttH_%s",varName[v].c_str()),500,500);
	if(v==1 || v==3)
	  DrawMatchCustom(cBCFullSel100MatchttH[v],hBCFullSel100[v][4][3],hMatchttH,prodName[4],colorsMatchttH,matchttH,nMatchttHStatuses,0.11,0.4,v==0);
	else
	  DrawMatchCustom(cBCFullSel100MatchttH[v],hBCFullSel100[v][4][3],hMatchttH,prodName[4],colorsMatchttH,matchttH,nMatchttHStatuses,0.34,0.4,v==0);
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
	Draw2D(c2DBCFullSel100[v2][p],h2DBCFullSel100[v2][p][3],prodName[p]);
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
	  Draw2D(c2DBCFullSel100Decays[v2][p][d],h2DBCFullSel100Decays[v2][p][d][3],prodName[p]+decayLabel[d]);
	  SaveCanvas(outDir,c2DBCFullSel100Decays[v2][p][d]);
	}
      }
    }

  }
  
  if(doBaskets){

    TCanvas* cBCFullSel100BasketsAll = new TCanvas("cBCFullSel100_baskets","cBCFullSel100_baskets",500,500);
    DrawByProdmodes(cBCFullSel100BasketsAll,hBCFullSel100Baskets,prodName,isPresent,colors,allZZBkgdsArePresent,false);
    SaveCanvas(outDir,cBCFullSel100BasketsAll);
    
    TCanvas* cBCFullSel100BasketEfficiency[nHiggsSamples+1];
    TH1F* hProdModes[nHiggsSamples]; for(int p=0; p<nHiggsSamples; p++) hProdModes[p] = (TH1F*)hBCFullSel100Baskets[p][3];
    TH1F* hAssocDecays[nAssocDecays]; for(int a=0; a<nAssocDecays; a++) hAssocDecays[a] = (TH1F*)hBCFullSel100BasketsAssocDecays[a][3];
    for(int p=0; p<nHiggsSamples; p++){
      string pn = prodName[p];
      cBCFullSel100BasketEfficiency[p] = new TCanvas(Form("cBCFullSel100_basketEfficiency_%s",pn.c_str()),Form("cBCFullSel100_basketEfficiency_%s",pn.c_str()),500,500);
      DrawBasketEfficiencies(cBCFullSel100BasketEfficiency[p],hBCFullSel100Baskets[p][3],pn,hAssocDecays,assocDecayName,basketLabel,false);
      SaveCanvas(outDir,cBCFullSel100BasketEfficiency[p]);
    }

    if(allSignalsArePresent){
      TCanvas* cBCFullSel100BasketPurity = new TCanvas("cBCFullSel100_basketPurity","cBCFullSel100_basketPurity",800,500);
      DrawBasketPurities(cBCFullSel100BasketPurity,hProdModes,prodName,hAssocDecays,treatH2l2XAsBkgd?assocDecayName3:assocDecayName2,false,treatH2l2XAsBkgd,basketLabel);
      SaveCanvas(outDir,cBCFullSel100BasketPurity);
      TCanvas* cBCFullSel100BasketPuritySplit = new TCanvas("cBCFullSel100_basketPuritySplit","cBCFullSel100_basketPuritySplit",800,500);
      DrawBasketPurities(cBCFullSel100BasketPuritySplit,hProdModes,prodName,hAssocDecays,treatH2l2XAsBkgd?assocDecayName3:assocDecayName2,true,treatH2l2XAsBkgd,basketLabel);
      SaveCanvas(outDir,cBCFullSel100BasketPuritySplit);
    }
    if(allSignalsArePresent && allZZBkgdsArePresent){

      TH1F* hSumSgnl = (TH1F*)hBCFullSel100Baskets[0][3]->Clone();
      for(int p=1; p<nHiggsSamples; p++) hSumSgnl->Add(hBCFullSel100Baskets[p][3]);
      TH1F* hSumBkgd = (TH1F*)hBCFullSel100Baskets[5][3]->Clone();
      hSumBkgd->Add(hBCFullSel100Baskets[6][3]);
      hSumBkgd->Add(hBCFullSel100Baskets[7][3]);
      if(treatH2l2XAsBkgd){
	hSumSgnl->Add(hBCFullSel100BasketsAssocDecays[nAssocWDecays+nAssocZDecays-1][3],-1);
	hSumSgnl->Add(hBCFullSel100BasketsAssocDecays[nAssocWDecays+nAssocZDecays+nAssocttDecays-1][3],-1);
	hSumBkgd->Add(hBCFullSel100BasketsAssocDecays[nAssocWDecays+nAssocZDecays-1][3]);
	hSumBkgd->Add(hBCFullSel100BasketsAssocDecays[nAssocWDecays+nAssocZDecays+nAssocttDecays-1][3]);
      }

      TH1F* hDenom = (TH1F*)hSumSgnl->Clone();
      hDenom->Add(hSumBkgd);
      TH1F* hSOSPB = (TH1F*)hSumSgnl->Clone();
      hSOSPB->Divide(hDenom);
      TCanvas* cBCFullSel100BasketSOSPB = new TCanvas("cBCFullSel100_basketSOSPB","cBCFullSel100_basketSOSPB",500,500);
      DrawBasketSOSPB(cBCFullSel100BasketSOSPB,hSOSPB,basketLabel);
      SaveCanvas(outDir,cBCFullSel100BasketSOSPB);

      string pn = "qqtoZZ";
      cBCFullSel100BasketEfficiency[nHiggsSamples] = new TCanvas(Form("cBCFullSel100_basketEfficiency_%s",pn.c_str()),Form("cBCFullSel100_basketEfficiency_%s",pn.c_str()),500,500);
      DrawBasketEfficiencies(cBCFullSel100BasketEfficiency[nHiggsSamples],hSumBkgd,pn,hAssocDecays,assocDecayName,basketLabel,false);
      SaveCanvas(outDir,cBCFullSel100BasketEfficiency[nHiggsSamples]);


      TH1F* hSumSgnlMasswindow = (TH1F*)hBCFullSel100MasswindowBaskets[0][3]->Clone();
      for(int p=1; p<nHiggsSamples; p++) hSumSgnlMasswindow->Add(hBCFullSel100MasswindowBaskets[p][3]);
      TH1F* hSumBkgdMasswindow = (TH1F*)hBCFullSel100MasswindowBaskets[5][3]->Clone();
      hSumBkgdMasswindow->Add(hBCFullSel100MasswindowBaskets[6][3]);
      hSumBkgdMasswindow->Add(hBCFullSel100MasswindowBaskets[7][3]);
      if(treatH2l2XAsBkgd){
	hSumSgnlMasswindow->Add(hBCFullSel100MasswindowBasketsAssocDecays[nAssocWDecays+nAssocZDecays-1][3],-1);
	hSumSgnlMasswindow->Add(hBCFullSel100MasswindowBasketsAssocDecays[nAssocWDecays+nAssocZDecays+nAssocttDecays-1][3],-1);
	hSumBkgdMasswindow->Add(hBCFullSel100MasswindowBasketsAssocDecays[nAssocWDecays+nAssocZDecays-1][3]);
	hSumBkgdMasswindow->Add(hBCFullSel100MasswindowBasketsAssocDecays[nAssocWDecays+nAssocZDecays+nAssocttDecays-1][3]);
      }

      TH1F* hDenomMasswindow = (TH1F*)hSumSgnlMasswindow->Clone();
      hDenomMasswindow->Add(hSumBkgdMasswindow);
      TH1F* hSOSPBMasswindow = (TH1F*)hSumSgnlMasswindow->Clone();
      hSOSPBMasswindow->Divide(hDenomMasswindow);
      TCanvas* cBCFullSel100MasswindowBasketSOSPB = new TCanvas("cBCFullSel100Masswindow_basketSOSPB","cBCFullSel100Masswindow_basketSOSPB",500,500);
      DrawBasketSOSPB(cBCFullSel100MasswindowBasketSOSPB,hSOSPBMasswindow,basketLabel);
      SaveCanvas(outDir,cBCFullSel100MasswindowBasketSOSPB);

    }

  }

}


