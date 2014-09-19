#include <iostream>
#include <utility>
#include <map>
#include <vector>
#include <cmath>

#include "TStyle.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TColor.h"

using namespace std;

#define MASSZ 91.1876

#define DEBUG 0

#define doProdComp       1
#define doProdCompMatch4 0
#define doMatch4OrNot    0
#define doMatchHLeps     0
#define doMatchAllLeps   0
#define doMatchWHZHttH   1
#define do2DPlots        1

#define nSamples 6
#define nHiggsSamples 5
#define nCategories 4
#define nVariables 6
#define nMatchHLepsStatuses 6
#define nMatchAllLepsStatuses 5
#define nMatchWHStatuses 5
#define nMatchZHStatuses 7
#define nMatchttHStatuses 9




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




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////// PLOTTING FUNCTIONS ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void DrawByProdmodes(TCanvas* c, TH1F* hSource[nSamples][nCategories], string* prodName, Bool_t* isPresent, Color_t* colors, Bool_t logY = false) {

  gStyle->SetOptTitle(0);

  c->cd();
  if(logY) c->SetLogy();
  c->SetTicks(0,0);

  TLegend* lgd = new TLegend(0.69,0.67,0.89,0.89);
  lgd->SetFillStyle(0);
  lgd->SetBorderSize(0);

  TH1F* h[nSamples];

  Float_t max = 0.;
  for(int p=0; p<nSamples; p++){
    if(!isPresent[p]) continue;
    h[p] = (TH1F*)hSource[p][3]->Clone(); // includes the 3 decay channels !!
    h[p]->Scale(1/h[p]->Integral());
    Float_t maxtemp = h[p]->GetMaximum();
    if(maxtemp>max) max = maxtemp;
  }

  for(int p=0; p<nSamples; p++){
    if(!isPresent[p]) continue;
    string pn = prodName[p];

    h[p]->SetLineColor(colors[p]);
    h[p]->SetLineWidth(2);
    h[p]->SetStats(0);
    //if(!logY) h[p]->SetMinimum(0);
    if(p==0) h[p]->SetMaximum(1.1*max);

    if(p==0) h[p]->Draw(); else h[p]->Draw("sames");
    
    lgd->AddEntry(h[p],pn.c_str(),"l");
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
  lgd->SetFillColor(0);
  for(int m=0; m<nMatchAllLepsStatuses; m++){
    lgd->AddEntry(hStacks[m],matchAllLepsKeys[m].c_str(),"f");
  }
  lgd->Draw();

  gPad->RedrawAxis();

}

void DrawMatchCustom(TCanvas* c, TH1F* h, TH1F** hMatch, string title, Color_t* color, string* matchWHKeys, Int_t nStatuses, Float_t legLeft, Float_t legBottom, Bool_t logY = false) {

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
    lgd->AddEntry(hStacks[m],matchWHKeys[m].c_str(),"f");
  }
  lgd->Draw();
  
  TPaveText* pav = new TPaveText(legLeft+0.03,legBottom,legLeft+0.11,0.89,"brNDC");
  pav->SetFillStyle(0);
  pav->SetBorderSize(0);
  pav->SetTextColor(0);
  for(int m=0; m<nStatuses; m++){
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


  // ------------------------------------------------------------
  // ---------------------- Definitions -------------------------
  // ------------------------------------------------------------

  string prodName[nSamples] = {
    "ggH",
    "VBF",
    "WH",
    "ZH",
    "ttH",
    "qqtoZZ",
  };
  Bool_t isPresent[nSamples] = {
    file_ggH != "",
    file_VBF != "",
    file_WH  != "",
    file_ZH  != "",
    file_ttH != "",
    file_ZZ4mu != "" && file_ZZ4e != "" && file_ZZ2e2mu != "",
  };

  Color_t colors[nSamples] = {kBlue, kGreen+2, kRed, kOrange+1, kMagenta+1, kRed-1};
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

  string categories[nCategories] = {
    "4mu",
    "4e",
    "2e2mu",
    "4mu + 4e + 2e2mu",
  };
  string categoriesShort[nCategories] = {
    "4mu",
    "4e",
    "2e2mu",
    "Sumcat",
  };
  string subdir[nCategories] = {
    "ZZ4muTree",
    "ZZ4eTree",
    "ZZ2e2muTree",
    "",
  };

  string varName[nVariables] = {
    "M4l",
    "MZ1",
    "MZ2",
    "KD",
    "Djet",
    "Pt4l",
  };
  string varLabel[nVariables] = {
    "m_{4l} (GeV)",
    "m_{Z_{1}} (GeV)",
    "m_{Z_{2}} (GeV)",
    "D_{bkg}^{kin}",
    "D_{jet}",
    "p_{T}^{4l} (GeV)",
  };
  Int_t varNbin[nVariables]  = { 100,  75,  75,  50,  50,  50 };
  Float_t varMin[nVariables] = {  50,   0,   0,   0,   0,   0 };
  Float_t varMax[nVariables] = { 850, 150, 150,   1,   2, 500 };

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

  TChain* chain[nSamples][nCategories];
  for(int p=0; p<nSamples; p++){
    if(!isPresent[p]) continue;
    for(int c=0; c<nCategories-1; c++){
      chain[p][c] = new TChain(Form("%s/candTree",subdir[c].c_str()));
      if(prodName[p]=="ggH") chain[p][c]->Add(file_ggH.c_str());
      if(prodName[p]=="VBF") chain[p][c]->Add(file_VBF.c_str());
      if(prodName[p]=="WH" ) chain[p][c]->Add(file_WH .c_str());
      if(prodName[p]=="ZH" ) chain[p][c]->Add(file_ZH .c_str());
      if(prodName[p]=="ttH") chain[p][c]->Add(file_ttH.c_str());
      if(prodName[p]=="qqtoZZ"){
	chain[p][c]->Add(file_ZZ4mu.c_str());
	chain[p][c]->Add(file_ZZ4e.c_str());
	chain[p][c]->Add(file_ZZ2e2mu.c_str());
      }
    }
  }

  Int_t nRun;
  Long64_t nEvent;
  Int_t nLumi;
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
  vector<Float_t> *Lep1Pt = 0;
  vector<Float_t> *Lep1Eta = 0;
  vector<Float_t> *Lep1Phi = 0;
  vector<Int_t> *Lep1LepId = 0;
  vector<Float_t> *Lep2Pt = 0;
  vector<Float_t> *Lep2Eta = 0;
  vector<Float_t> *Lep2Phi = 0;
  vector<Int_t> *Lep2LepId = 0;
  vector<Float_t> *Lep3Pt = 0;
  vector<Float_t> *Lep3Eta = 0;
  vector<Float_t> *Lep3Phi = 0;
  vector<Int_t> *Lep3LepId = 0;
  vector<Float_t> *Lep4Pt = 0;
  vector<Float_t> *Lep4Eta = 0;
  vector<Float_t> *Lep4Phi = 0;
  vector<Int_t> *Lep4LepId = 0;
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

  Int_t nbStored[nSamples][nCategories];
  Int_t nbWithBC[nSamples][nCategories];
  Int_t nbWithBCFullSel70[nSamples][nCategories];
  Int_t nbWithBCFullSel100[nSamples][nCategories];
  Int_t nbWithBCFullSel100MatchHLeps[nMatchHLepsStatuses][nSamples][nCategories];
  Int_t nbWithBCFullSel100MatchAllLeps[nMatchAllLepsStatuses][nSamples][nCategories];
  Int_t nbWithBCFullSel100MatchWH[nMatchWHStatuses][nCategories];
  Int_t nbWithBCFullSel100MatchZH[nMatchZHStatuses][nCategories];
  Int_t nbWithBCFullSel100MatchttH[nMatchttHStatuses][nCategories];

  TH1F* hBCFullSel100[nVariables][nSamples][nCategories];
  TH1F* hBCFullSel100MatchHLeps[nVariables][nMatchHLepsStatuses][nSamples][nCategories];
  TH1F* hBCFullSel100MatchAllLeps[nVariables][nMatchAllLepsStatuses][nSamples][nCategories];
  TH1F* hBCFullSel100MatchWH[nVariables][nMatchWHStatuses][nCategories];
  TH1F* hBCFullSel100MatchZH[nVariables][nMatchZHStatuses][nCategories];
  TH1F* hBCFullSel100MatchttH[nVariables][nMatchttHStatuses][nCategories];

  TH2F* h2DBCFullSel100[n2DHist][nSamples][nCategories];
  TH2F* h2DBCFullSel100Decays[n2DHist][nSamples][nDecays][nCategories];

  string nbZDaughtersFromH[4] = { "2", "1", "0", "ambig." };
  string WHdecays[2] = {
    "H->ZZ->4l, W->X     ",
    "H->ZZ->4l, W->lnu   ",
  };
  Int_t nbZ1DaughtersFromHWH[2][4]; for(int i=0; i<2; i++) for(int j=0; j<4; j++) nbZ1DaughtersFromHWH[i][j] = 0;
  Int_t nbZ2DaughtersFromHWH[2][4]; for(int i=0; i<2; i++) for(int j=0; j<4; j++) nbZ2DaughtersFromHWH[i][j] = 0;
  string ZHdecays[3] = {
    "H->ZZ->4l, Z->X     ",
    "H->ZZ->4l, Z->2l    ",
    "H->ZZ->2l2X, Z->2l  ",
  };
  Int_t nbZ1DaughtersFromHZH[3][4]; for(int i=0; i<3; i++) for(int j=0; j<4; j++) nbZ1DaughtersFromHZH[i][j] = 0;
  Int_t nbZ2DaughtersFromHZH[3][4]; for(int i=0; i<3; i++) for(int j=0; j<4; j++) nbZ2DaughtersFromHZH[i][j] = 0;
  string ttHdecays[4] = {
    "H->ZZ->4l, tt->X    ",
    "H->ZZ->4l, tt->lX   ",
    "H->ZZ->4l, tt->2lX  ",
    "H->ZZ->2l2X, tt->2lX",
  };
  Int_t nbZ1DaughtersFromHttH[4][4]; for(int i=0; i<4; i++) for(int j=0; j<4; j++) nbZ1DaughtersFromHttH[i][j] = 0;
  Int_t nbZ2DaughtersFromHttH[4][4]; for(int i=0; i<4; i++) for(int j=0; j<4; j++) nbZ2DaughtersFromHttH[i][j] = 0;



  // ------------------------------------------------------------
  // ---------------------- Processing --------------------------
  // ------------------------------------------------------------

  for(int p=0; p<nSamples; p++){
    if(!isPresent[p]) continue;

    cout<<prodName[p]<<endl;

    for(int c=0; c<nCategories; c++){
      
      string suffix = "_"+prodName[p]+"_"+categoriesShort[c];
      string title = prodName[p]+", category "+categories[c];

      nbStored[p][c] = 0;
      nbWithBC[p][c] = 0;
      nbWithBCFullSel70[p][c] = 0;
      nbWithBCFullSel100[p][c] = 0;
      for(int m=0; m<nMatchHLepsStatuses; m++) nbWithBCFullSel100MatchHLeps[m][p][c] = 0;
      for(int m=0; m<nMatchAllLepsStatuses; m++) nbWithBCFullSel100MatchAllLeps[m][p][c] = 0;
      if(prodName[p]=="WH") for(int m=0; m<nMatchWHStatuses; m++) nbWithBCFullSel100MatchWH[m][c] = 0;
      if(prodName[p]=="ZH") for(int m=0; m<nMatchZHStatuses; m++) nbWithBCFullSel100MatchZH[m][c] = 0;
      if(prodName[p]=="ttH") for(int m=0; m<nMatchttHStatuses; m++) nbWithBCFullSel100MatchttH[m][c] = 0;

      for(int v=0; v<nVariables; v++){
	hBCFullSel100[v][p][c] = new TH1F(("hBCFullSel100_"+varName[v]+suffix).c_str(),(title+";"+varLabel[v]+";entries").c_str(),varNbin[v],varMin[v],varMax[v]);
	for(int m=0; m<nMatchHLepsStatuses; m++)
	  hBCFullSel100MatchHLeps[v][m][p][c] = new TH1F(("hBCFullSel100MatchHLeps_"+varName[v]+"_"+matchHLepsInfix[m]+suffix).c_str(),(title+";"+varLabel[v]+";entries").c_str(),varNbin[v],varMin[v],varMax[v]);
	for(int m=0; m<nMatchAllLepsStatuses; m++)
	  hBCFullSel100MatchAllLeps[v][m][p][c] = new TH1F(("hBCFullSel100MatchAllLeps_"+varName[v]+"_"+matchAllLepsInfix[m]+suffix).c_str(),(title+";"+varLabel[v]+";entries").c_str(),varNbin[v],varMin[v],varMax[v]);
	if(prodName[p]=="WH")
	  for(int m=0; m<nMatchWHStatuses; m++)
	    hBCFullSel100MatchWH[v][m][c] = new TH1F(Form("hBCFullSel100MatchWH_%s_%i_%s",varName[v].c_str(),m,categoriesShort[c].c_str()),(title+";"+varLabel[v]+";entries").c_str(),varNbin[v],varMin[v],varMax[v]);
	if(prodName[p]=="ZH")
	  for(int m=0; m<nMatchZHStatuses; m++)
	    hBCFullSel100MatchZH[v][m][c] = new TH1F(Form("hBCFullSel100MatchZH_%s_%i_%s",varName[v].c_str(),m,categoriesShort[c].c_str()),(title+";"+varLabel[v]+";entries").c_str(),varNbin[v],varMin[v],varMax[v]);
	if(prodName[p]=="ttH")
	  for(int m=0; m<nMatchttHStatuses; m++)
	    hBCFullSel100MatchttH[v][m][c] = new TH1F(Form("hBCFullSel100MatchttH_%s_%i_%s",varName[v].c_str(),m,categoriesShort[c].c_str()),(title+";"+varLabel[v]+";entries").c_str(),varNbin[v],varMin[v],varMax[v]);
      }

      for(int v2=0; v2<n2DHist; v2++){
	h2DBCFullSel100[v2][p][c] = new TH2F(("h2DBCFullSel100_"+varName[varXindex[v2]]+"_"+varName[varYindex[v2]]+suffix).c_str(),(title+";"+varLabel[varXindex[v2]]+";"+varLabel[varYindex[v2]]).c_str(),varNbin[varXindex[v2]],varMin[varXindex[v2]],varMax[varXindex[v2]],varNbin[varYindex[v2]],varMin[varYindex[v2]],varMax[varYindex[v2]]);
	for(int d=0; d<nDecays; d++){
	  h2DBCFullSel100Decays[v2][p][d][c] = new TH2F(("h2DBCFullSel100_"+decayInfix[d]+"_"+varName[varXindex[v2]]+"_"+varName[varYindex[v2]]+suffix).c_str(),(title+";"+varLabel[varXindex[v2]]+";"+varLabel[varYindex[v2]]).c_str(),varNbin[varXindex[v2]],varMin[varXindex[v2]],varMax[varXindex[v2]],varNbin[varYindex[v2]],varMin[varYindex[v2]],varMax[varYindex[v2]]);
	}
      }

    }

    for(int c=0; c<nCategories-1; c++){

      chain[p][c]->SetBranchAddress("RunNumber", &nRun);
      chain[p][c]->SetBranchAddress("EventNumber", &nEvent);
      chain[p][c]->SetBranchAddress("LumiNumber", &nLumi);
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
      chain[p][c]->SetBranchAddress("Lep1Pt", &Lep1Pt);
      chain[p][c]->SetBranchAddress("Lep1Eta", &Lep1Eta);
      chain[p][c]->SetBranchAddress("Lep1Phi", &Lep1Phi);
      chain[p][c]->SetBranchAddress("Lep1LepId", &Lep1LepId);
      chain[p][c]->SetBranchAddress("Lep2Pt", &Lep2Pt);
      chain[p][c]->SetBranchAddress("Lep2Eta", &Lep2Eta);
      chain[p][c]->SetBranchAddress("Lep2Phi", &Lep2Phi);
      chain[p][c]->SetBranchAddress("Lep2LepId", &Lep2LepId);
      chain[p][c]->SetBranchAddress("Lep3Pt", &Lep3Pt);
      chain[p][c]->SetBranchAddress("Lep3Eta", &Lep3Eta);
      chain[p][c]->SetBranchAddress("Lep3Phi", &Lep3Phi);
      chain[p][c]->SetBranchAddress("Lep3LepId", &Lep3LepId);
      chain[p][c]->SetBranchAddress("Lep4Pt", &Lep4Pt);
      chain[p][c]->SetBranchAddress("Lep4Eta", &Lep4Eta);
      chain[p][c]->SetBranchAddress("Lep4Phi", &Lep4Phi);
      chain[p][c]->SetBranchAddress("Lep4LepId", &Lep4LepId);
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

	printStatus(z,10000,entries,"entries");

	chain[p][c]->GetEntry(z);
 
	Int_t nGenHLep = 0;
	Int_t nGenAssocLep = 0;
	Float_t GenHLepId[4] = {GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id};
	Float_t GenAssocLepId[2] = {GenAssocLep1Id,GenAssocLep2Id};
	for(Int_t iGen=0; iGen<4; iGen++){
	  if(GenHLepId[iGen]!=0){
	    nGenHLep++;
	  }
	}
	for(Int_t iGenAssoc=0; iGenAssoc<2; iGenAssoc++){
	  if(GenAssocLepId[iGenAssoc]!=0){
	    nGenAssocLep++;
	  }
	}

	/*
	if(nGenHLep!=4) continue;
	//*/

	string evtID = eventID(nRun,nLumi,nEvent);
	overlapMapStored[p][evtID].push_back(categories[c]);

	nbStored[p][c]++;
	nbStored[p][3]++;

	if(iBC>=0){

	  overlapMapWithBC[p][evtID].push_back(categories[c]);

	  nbWithBC[p][c]++;
	  nbWithBC[p][3]++;

	  Bool_t FullSel70 = ZZsel->at(iBC)>=90;
	  Bool_t FullSel100 = ZZsel->at(iBC)>=100;
	  
	  if(FullSel70){
	    overlapMapWithBCFullSel70[p][evtID].push_back(categories[c]);
	    nbWithBCFullSel70[p][c]++;
	    nbWithBCFullSel70[p][3]++;
	  }
	  
	  if(FullSel100){

	    overlapMapWithBCFullSel100[p][evtID].push_back(categories[c]);

	    nbWithBCFullSel100[p][c]++;
	    nbWithBCFullSel100[p][3]++;

	    Float_t varVal[nVariables] = {
	      ZZMass->at(iBC),
	      Z1Mass->at(iBC),
	      Z2Mass->at(iBC),
	      p0plus_VAJHU->at(iBC) / ( p0plus_VAJHU->at(iBC) + bkg_VAMCFM->at(iBC) ),
	      ZZFisher->at(iBC),
	      ZZPt->at(iBC),
	    };

	    for(int v=0; v<nVariables; v++){
	      hBCFullSel100[v][p][c]->Fill(varVal[v]);
	      hBCFullSel100[v][p][3]->Fill(varVal[v]);
	    }

	    Int_t nMatchedToGenHLep[4] = {0,0,0,0};
	    Int_t nMatchedToGenAssocLep[2] = {0,0};
	    Int_t nMatchedToRecoLep[4] = {0,0,0,0};
	    Int_t nGenHLepMatchedToZ1Lep[4] = {0,0};
	    Int_t nGenHLepMatchedToZ2Lep[4] = {0,0};
	    Float_t GenHLepEta[4] = {GenLep1Eta,GenLep2Eta,GenLep3Eta,GenLep4Eta};
	    Float_t GenHLepPhi[4] = {GenLep1Phi,GenLep2Phi,GenLep3Phi,GenLep4Phi};
	    Float_t GenAssocLepEta[2] = {GenAssocLep1Eta,GenAssocLep2Eta};
	    Float_t GenAssocLepPhi[2] = {GenAssocLep1Phi,GenAssocLep2Phi};
	    Float_t RecoLepEta[4] = {Lep1Eta->at(iBC),Lep2Eta->at(iBC),Lep3Eta->at(iBC),Lep4Eta->at(iBC)};
	    Float_t RecoLepPhi[4] = {Lep1Phi->at(iBC),Lep2Phi->at(iBC),Lep3Phi->at(iBC),Lep4Phi->at(iBC)};
	    for(Int_t iGen=0; iGen<4; iGen++){
	      if(GenHLepId[iGen]!=0){
		for(Int_t iReco=0; iReco<4; iReco++){
		  if(deltaR(GenHLepEta[iGen],GenHLepPhi[iGen],RecoLepEta[iReco],RecoLepPhi[iReco]) < 0.1){
		    nMatchedToGenHLep[iGen]++;
		    nMatchedToRecoLep[iReco]++;
		    if(iReco<2){
		      nGenHLepMatchedToZ1Lep[iReco]++;
		    }else{
		      nGenHLepMatchedToZ2Lep[iReco-2]++;
		    }
		  }
		}
	      }
	    }
	    for(Int_t iGenAssoc=0; iGenAssoc<2; iGenAssoc++){
	      if(GenAssocLepId[iGenAssoc]!=0){
		for(Int_t iReco=0; iReco<4; iReco++){
		  if(deltaR(GenAssocLepEta[iGenAssoc],GenAssocLepPhi[iGenAssoc],RecoLepEta[iReco],RecoLepPhi[iReco]) < 0.1){
		    nMatchedToGenAssocLep[iGenAssoc]++;
		    nMatchedToRecoLep[iReco]++;
		  }
		}
	      }
	    }

	    Int_t currentMatchHLepsStatus = -1;
	    Bool_t foundAmbiguityHLeps = false;
	    for(Int_t iGen=0; iGen<4; iGen++) if(nMatchedToGenHLep[iGen]>1){ foundAmbiguityHLeps = true; break; }
	    for(Int_t iReco=0; iReco<4; iReco++) if(nMatchedToRecoLep[iReco]>1){ foundAmbiguityHLeps = true; break; }
	    if(foundAmbiguityHLeps){
	      currentMatchHLepsStatus = 5;
	    }else{
	      Int_t nOnesHLeps = 0;
	      for(Int_t iGen=0; iGen<4; iGen++) if(nMatchedToGenHLep[iGen]==1) nOnesHLeps++;
	      if(nOnesHLeps==4) currentMatchHLepsStatus = 0;
	      if(nOnesHLeps==3) currentMatchHLepsStatus = 1;
	      if(nOnesHLeps==2) currentMatchHLepsStatus = 2;
	      if(nOnesHLeps==1) currentMatchHLepsStatus = 3;
	      if(nOnesHLeps==0) currentMatchHLepsStatus = 4;
	    }
	    Int_t currentMatchAllLepsStatus = -1;
	    Bool_t foundAmbiguityAllLeps = false;
	    for(Int_t iGen=0; iGen<4; iGen++) if(nMatchedToGenHLep[iGen]>1){ foundAmbiguityAllLeps = true; break; }
	    for(Int_t iGenAssoc=0; iGenAssoc<2; iGenAssoc++) if(nMatchedToGenAssocLep[iGenAssoc]>1){ foundAmbiguityAllLeps = true; break; }
	    for(Int_t iReco=0; iReco<4; iReco++) if(nMatchedToRecoLep[iReco]>1){ foundAmbiguityAllLeps = true; break; }
	    Int_t nOnesHLeps = 0;
	    Int_t nOnesAssocLeps = 0;
	    for(Int_t iGen=0; iGen<4; iGen++) if(nMatchedToGenHLep[iGen]==1) nOnesHLeps++;
	    for(Int_t iGenAssoc=0; iGenAssoc<4; iGenAssoc++) if(nMatchedToGenAssocLep[iGenAssoc]==1) nOnesAssocLeps++;
	    if(foundAmbiguityAllLeps){
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
	      if(foundAmbiguityAllLeps){
		currentMatchWHStatus = 4;
	      }else{
		if(nOnesHLeps+nOnesAssocLeps<4){
		  currentMatchWHStatus = 3;
		}else{
		  if(nGenHLep==4 && nGenAssocLep==0){
		    if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchWHStatus = 0;
		    else cout<<"error"<<endl;
		  }else if(nGenHLep==4 && nGenAssocLep==1){
		    if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchWHStatus = 1;
		    else if(nOnesHLeps==3 && nOnesAssocLeps==1) currentMatchWHStatus = 2;
		    else cout<<"error"<<endl;
		  }else{
		    cout<<"error"<<endl;
		  }
		}
	      } 
	    }
	    if(prodName[p]=="ZH"){
	      if(foundAmbiguityAllLeps){
		currentMatchZHStatus = 6;
	      }else{
		if(nOnesHLeps+nOnesAssocLeps<4){
		  currentMatchZHStatus = 5;
		}else{
		  if(nGenHLep==4 && nGenAssocLep==0){
		    if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchZHStatus = 0;
		    else cout<<"error"<<endl;
		  }else if(nGenHLep==4 && nGenAssocLep==2){
		    if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchZHStatus = 1;
		    else if(nOnesHLeps==3 && nOnesAssocLeps==1) currentMatchZHStatus = 2;
		    else if(nOnesHLeps==2 && nOnesAssocLeps==2) currentMatchZHStatus = 3;
		    else cout<<"error"<<endl;
		  }else if(nGenHLep==2 && nGenAssocLep==2){
		    if(nOnesHLeps==2 && nOnesAssocLeps==2) currentMatchZHStatus = 4;
		    else cout<<"error"<<endl;
		  }else{
		    cout<<"error"<<endl;
		  }
		}
	      } 
	    }
	    if(prodName[p]=="ttH"){
	      if(foundAmbiguityAllLeps){
		currentMatchttHStatus = 8;
	      }else{
		if(nOnesHLeps+nOnesAssocLeps<4){
		  currentMatchttHStatus = 7;
		}else{
		  if(nGenHLep==4 && nGenAssocLep==0){
		    if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchttHStatus = 0;
		    else cout<<"error"<<endl;
		  }else if(nGenHLep==4 && nGenAssocLep==1){
		    if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchttHStatus = 1;
		    else if(nOnesHLeps==3 && nOnesAssocLeps==1) currentMatchttHStatus = 2;
		    else cout<<"error"<<endl;
		  }else if(nGenHLep==4 && nGenAssocLep==2){
		    if(nOnesHLeps==4 && nOnesAssocLeps==0) currentMatchttHStatus = 3;
		    else if(nOnesHLeps==3 && nOnesAssocLeps==1) currentMatchttHStatus = 4;
		    else if(nOnesHLeps==2 && nOnesAssocLeps==2) currentMatchttHStatus = 5;
		    else cout<<"error"<<endl;
		  }else if(nGenHLep==2 && nGenAssocLep==2){
		    if(nOnesHLeps==2 && nOnesAssocLeps==2) currentMatchttHStatus = 6;
		    else cout<<"error"<<endl;
		  }else{
		    cout<<"error"<<endl;
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
	    if(prodName[p]=="WH"){
	      if(nGenHLep==4 && nGenAssocLep==0){
		nbZ1DaughtersFromHWH[0][currentZ1MatchStatus]++;
		nbZ2DaughtersFromHWH[0][currentZ2MatchStatus]++;
	      }else if(nGenHLep==4 && nGenAssocLep==1){
		nbZ1DaughtersFromHWH[1][currentZ1MatchStatus]++;
		nbZ2DaughtersFromHWH[1][currentZ2MatchStatus]++;
	      }else{
		cout<<"error"<<endl;
	      }
	    }
	    if(prodName[p]=="ZH"){
	      if(nGenHLep==4 && nGenAssocLep==0){
		nbZ1DaughtersFromHZH[0][currentZ1MatchStatus]++;
		nbZ2DaughtersFromHZH[0][currentZ2MatchStatus]++;
	      }else if(nGenHLep==4 && nGenAssocLep==2){
		nbZ1DaughtersFromHZH[1][currentZ1MatchStatus]++;
		nbZ2DaughtersFromHZH[1][currentZ2MatchStatus]++;
	      }else if(nGenHLep==2 && nGenAssocLep==2){
		nbZ1DaughtersFromHZH[2][currentZ1MatchStatus]++;
		nbZ2DaughtersFromHZH[2][currentZ2MatchStatus]++;
	      }else{
		cout<<"error"<<endl;
	      }
	    }
	    if(prodName[p]=="ttH"){
	      if(nGenHLep==4 && nGenAssocLep==0){
		nbZ1DaughtersFromHttH[0][currentZ1MatchStatus]++;
		nbZ2DaughtersFromHttH[0][currentZ2MatchStatus]++;
	      }else if(nGenHLep==4 && nGenAssocLep==1){
		nbZ1DaughtersFromHttH[1][currentZ1MatchStatus]++;
		nbZ2DaughtersFromHttH[1][currentZ2MatchStatus]++;
	      }else if(nGenHLep==4 && nGenAssocLep==2){
		nbZ1DaughtersFromHttH[2][currentZ1MatchStatus]++;
		nbZ2DaughtersFromHttH[2][currentZ2MatchStatus]++;
	      }else if(nGenHLep==2 && nGenAssocLep==2){
		nbZ1DaughtersFromHttH[3][currentZ1MatchStatus]++;
		nbZ2DaughtersFromHttH[3][currentZ2MatchStatus]++;
	      }else{
		cout<<"error"<<endl;
	      }
	    }
	    // 	cout<<nGenHLep;
	    // 	cout<<" ";
	    // 	cout<<nGenAssocLep;
	    // 	cout<<" ";
	    // 	for(Int_t iGen =0; iGen <4; iGen ++) cout<<nMatchedToGenHLep [iGen ];
	    // 	cout<<" ";
	    // 	for(Int_t iGenAssoc =0; iGenAssoc <2; iGenAssoc ++) cout<<nMatchedToGenAssocLep[iGenAssoc];
	    // 	cout<<" ";
	    // 	for(Int_t iReco=0; iReco<4; iReco++) cout<<nMatchedToRecoLep[iReco];
	    // 	cout<<" ";
	    // 	cout<<currentMatchHLepsStatus;
	    // 	cout<<" ";
	    // 	cout<<currentMatchAllLepsStatus;
	    // 	cout<<endl;


	    for(int v2=0; v2<n2DHist; v2++){
	      h2DBCFullSel100[v2][p][c]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]]);
	      h2DBCFullSel100[v2][p][3]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]]);
	      h2DBCFullSel100Decays[v2][p][0][c]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]]);
	      h2DBCFullSel100Decays[v2][p][0][3]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]]);
	      if(nGenHLep==4){
		h2DBCFullSel100Decays[v2][p][1][c]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]]);
		h2DBCFullSel100Decays[v2][p][1][3]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]]);
		if(currentMatchHLepsStatus==0){
		  h2DBCFullSel100Decays[v2][p][3][c]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]]);
		  h2DBCFullSel100Decays[v2][p][3][3]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]]);
		}else if(currentMatchHLepsStatus<5){
		  h2DBCFullSel100Decays[v2][p][4][c]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]]);
		  h2DBCFullSel100Decays[v2][p][4][3]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]]);
		}
	      }else{
		h2DBCFullSel100Decays[v2][p][2][c]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]]);
		h2DBCFullSel100Decays[v2][p][2][3]->Fill(varVal[varXindex[v2]],varVal[varYindex[v2]]);
	      }
	    }

	    nbWithBCFullSel100MatchHLeps[currentMatchHLepsStatus][p][c]++;
	    nbWithBCFullSel100MatchHLeps[currentMatchHLepsStatus][p][3]++;
	    for(int v=0; v<nVariables; v++){
	      hBCFullSel100MatchHLeps[v][currentMatchHLepsStatus][p][c]->Fill(varVal[v]);
	      hBCFullSel100MatchHLeps[v][currentMatchHLepsStatus][p][3]->Fill(varVal[v]);
	    }

	    nbWithBCFullSel100MatchAllLeps[currentMatchAllLepsStatus][p][c]++;
	    nbWithBCFullSel100MatchAllLeps[currentMatchAllLepsStatus][p][3]++;
	    for(int v=0; v<nVariables; v++){
	      hBCFullSel100MatchAllLeps[v][currentMatchAllLepsStatus][p][c]->Fill(varVal[v]);
	      hBCFullSel100MatchAllLeps[v][currentMatchAllLepsStatus][p][3]->Fill(varVal[v]);
	    }
	    
	    if(prodName[p]=="WH"){
	      nbWithBCFullSel100MatchWH[currentMatchWHStatus][c]++;
	      nbWithBCFullSel100MatchWH[currentMatchWHStatus][3]++;
	      for(int v=0; v<nVariables; v++){
		hBCFullSel100MatchWH[v][currentMatchWHStatus][c]->Fill(varVal[v]);
		hBCFullSel100MatchWH[v][currentMatchWHStatus][3]->Fill(varVal[v]);
	      }
	    }
	    if(prodName[p]=="ZH"){
	      nbWithBCFullSel100MatchZH[currentMatchZHStatus][c]++;
	      nbWithBCFullSel100MatchZH[currentMatchZHStatus][3]++;
	      for(int v=0; v<nVariables; v++){
		hBCFullSel100MatchZH[v][currentMatchZHStatus][c]->Fill(varVal[v]);
		hBCFullSel100MatchZH[v][currentMatchZHStatus][3]->Fill(varVal[v]);
	      }
	    }
	    if(prodName[p]=="ttH"){
	      nbWithBCFullSel100MatchttH[currentMatchttHStatus][c]++;
	      nbWithBCFullSel100MatchttH[currentMatchttHStatus][3]++;
	      for(int v=0; v<nVariables; v++){
		hBCFullSel100MatchttH[v][currentMatchttHStatus][c]->Fill(varVal[v]);
		hBCFullSel100MatchttH[v][currentMatchttHStatus][3]->Fill(varVal[v]);
	      }
	    }
	
	  }

	}

      } // end for entries

    } // end for categories

    Int_t widthColumn2 = 6;

    for(int c=0; c<nCategories; c++){
      cout<<" "<<categories[c]<<endl;
      cout<<"  stored :            "<<fixWidth(Form("%i",nbStored[p][c]),widthColumn2,false)<<endl;
      cout<<"  iBC>0 :             "<<fixWidth(Form("%i",nbWithBC[p][c]),widthColumn2,false)<<endl;
      cout<<"  iBC>0, FullSel70 :  "<<fixWidth(Form("%i",nbWithBCFullSel70[p][c]),widthColumn2,false)<<endl;
      cout<<"  iBC>0, FullSel100 : "<<fixWidth(Form("%i",nbWithBCFullSel100[p][c]),widthColumn2,false)<<endl;
      cout<<"  iBC>0, FullSel100, match to 4 gen leptons : "<<fixWidth(Form("%i",nbWithBCFullSel100MatchHLeps[0][p][c]),widthColumn2,false)<<endl;
    }

    cout<<" "<<"# events :"<<endl;
    cout<<"  stored :            "<<fixWidth(Form("%i",(int)overlapMapStored[p].size()),widthColumn2,false)<<endl;
    cout<<"  iBC>0 :             "<<fixWidth(Form("%i",(int)overlapMapWithBC[p].size()),widthColumn2,false)<<endl;
    cout<<"  iBC>0, FullSel70 :  "<<fixWidth(Form("%i",(int)overlapMapWithBCFullSel70[p].size()),widthColumn2,false)<<endl;
    cout<<"  iBC>0, FullSel100 : "<<fixWidth(Form("%i",(int)overlapMapWithBCFullSel100[p].size()),widthColumn2,false)<<endl;

    cout<<" "<<"# events stored in 2:3 decay channels :"<<endl;
    Int_t nbOverlapStored2 = 0;
    Int_t nbOverlapStored3 = 0;
    for(map<string,vector<string> >::iterator it=overlapMapStored[p].begin(); it!=overlapMapStored[p].end(); it++){
      if(it->second.size()==2) nbOverlapStored2++;
      if(it->second.size()==3) nbOverlapStored3++;
    }
    cout<<"  stored :            "<<nbOverlapStored2<<":"<<nbOverlapStored3<<" ("<<percentage((float)(nbOverlapStored2+nbOverlapStored3)/(int)overlapMapStored[p].size())<<" %)"<<endl;
    Int_t nbOverlapWithBC2 = 0;
    Int_t nbOverlapWithBC3 = 0;
    for(map<string,vector<string> >::iterator it=overlapMapWithBC[p].begin(); it!=overlapMapWithBC[p].end(); it++){
      if(it->second.size()==2) nbOverlapWithBC2++;
      if(it->second.size()==3) nbOverlapWithBC3++;
    }
    cout<<"  iBC>0 :             "<<nbOverlapWithBC2<<":"<<nbOverlapWithBC3<<" ("<<percentage((float)(nbOverlapWithBC2+nbOverlapWithBC3)/(int)overlapMapWithBC[p].size())<<" %)"<<endl;
    Int_t nbOverlapWithBCFullSel702 = 0;
    Int_t nbOverlapWithBCFullSel703 = 0;
    for(map<string,vector<string> >::iterator it=overlapMapWithBCFullSel70[p].begin(); it!=overlapMapWithBCFullSel70[p].end(); it++){
      if(it->second.size()==2) nbOverlapWithBCFullSel702++;
      if(it->second.size()==3) nbOverlapWithBCFullSel703++;
    }
    cout<<"  iBC>0, FullSel70 :  "<<nbOverlapWithBCFullSel702<<":"<<nbOverlapWithBCFullSel703<<" ("<<percentage((float)(nbOverlapWithBCFullSel702+nbOverlapWithBCFullSel703)/(int)overlapMapWithBCFullSel70[p].size())<<" %)"<<endl;
    Int_t nbOverlapWithBCFullSel1002 = 0;
    Int_t nbOverlapWithBCFullSel1003 = 0;
    for(map<string,vector<string> >::iterator it=overlapMapWithBCFullSel100[p].begin(); it!=overlapMapWithBCFullSel100[p].end(); it++){
      if(it->second.size()==2) nbOverlapWithBCFullSel1002++;
      if(it->second.size()==3) nbOverlapWithBCFullSel1003++;
    }
    cout<<"  iBC>0, FullSel100 : "<<nbOverlapWithBCFullSel1002<<":"<<nbOverlapWithBCFullSel1003<<" ("<<percentage((float)(nbOverlapWithBCFullSel1002+nbOverlapWithBCFullSel1003)/(int)overlapMapWithBCFullSel100[p].size())<<" %)"<<endl;

    int widthColumn1 = 22;
    int widthOtherColumns = 7;
    if(prodName[p]=="WH"){
      cout<<" "<<repeat("-",widthColumn1+4*(2+widthOtherColumns))<<endl;
      cout<<" "<<fixWidth("# of Z1 leptons from H",widthColumn1,true);
      for(int j=0; j<4; j++) cout<<"  "<<fixWidth(nbZDaughtersFromH[j],widthOtherColumns,false);
      cout<<endl;
      for(int i=0; i<2; i++){
	cout<<" ";
	cout<<fixWidth(WHdecays[i],widthColumn1,true);
	int total = 0; for(int j=0; j<4; j++) total += nbZ1DaughtersFromHWH[i][j];
	for(int j=0; j<4; j++) cout<<"  "<<fixWidth(percentage((float)nbZ1DaughtersFromHWH[i][j]/total)+" %",widthOtherColumns,false);
	cout<<endl;
      }
      cout<<" "<<repeat("-",widthColumn1+4*(2+widthOtherColumns))<<endl;
      cout<<" "<<fixWidth("# of Z2 leptons from H",widthColumn1,true);
      for(int j=0; j<4; j++) cout<<"  "<<fixWidth(nbZDaughtersFromH[j],widthOtherColumns,false);
      cout<<endl;
      for(int i=0; i<2; i++){
	cout<<" ";
	cout<<fixWidth(WHdecays[i],widthColumn1,true);
	int total = 0; for(int j=0; j<4; j++) total += nbZ2DaughtersFromHWH[i][j];
	for(int j=0; j<4; j++) cout<<"  "<<fixWidth(percentage((float)nbZ2DaughtersFromHWH[i][j]/total)+" %",widthOtherColumns,false);
	cout<<endl;
      }
    }
    if(prodName[p]=="ZH"){
      cout<<" "<<repeat("-",widthColumn1+4*(2+widthOtherColumns))<<endl;
      cout<<" "<<fixWidth("# of Z1 leptons from H",widthColumn1,true);
      for(int j=0; j<4; j++) cout<<"  "<<fixWidth(nbZDaughtersFromH[j],widthOtherColumns,false);
      cout<<endl;
      for(int i=0; i<3; i++){
	cout<<" ";
	cout<<fixWidth(ZHdecays[i],widthColumn1,true);
	int total = 0; for(int j=0; j<4; j++) total += nbZ1DaughtersFromHZH[i][j];
	for(int j=0; j<4; j++) cout<<"  "<<fixWidth(percentage((float)nbZ1DaughtersFromHZH[i][j]/total)+" %",widthOtherColumns,false);
	cout<<endl;
      }
      cout<<" "<<repeat("-",widthColumn1+4*(2+widthOtherColumns))<<endl;
      cout<<" "<<fixWidth("# of Z2 leptons from H",widthColumn1,true);
      for(int j=0; j<4; j++) cout<<"  "<<fixWidth(nbZDaughtersFromH[j],widthOtherColumns,false);
      cout<<endl;
      for(int i=0; i<3; i++){
	cout<<" ";
	cout<<fixWidth(ZHdecays[i],widthColumn1,true);
	int total = 0; for(int j=0; j<4; j++) total += nbZ2DaughtersFromHZH[i][j];
	for(int j=0; j<4; j++) cout<<"  "<<fixWidth(percentage((float)nbZ2DaughtersFromHZH[i][j]/total)+" %",widthOtherColumns,false);
	cout<<endl;
      }
    }
    if(prodName[p]=="ttH"){
      cout<<" "<<repeat("-",widthColumn1+4*(2+widthOtherColumns))<<endl;
      cout<<" "<<fixWidth("# of Z1 leptons from H",widthColumn1,true);
      for(int j=0; j<4; j++) cout<<"  "<<fixWidth(nbZDaughtersFromH[j],widthOtherColumns,false);
      cout<<endl;
      for(int i=0; i<4; i++){
	cout<<" ";
	cout<<fixWidth(ttHdecays[i],widthColumn1,true);
	int total = 0; for(int j=0; j<4; j++) total += nbZ1DaughtersFromHttH[i][j];
	for(int j=0; j<4; j++) cout<<"  "<<fixWidth(percentage((float)nbZ1DaughtersFromHttH[i][j]/total)+" %",widthOtherColumns,false);
	cout<<endl;
      }
      cout<<" "<<repeat("-",widthColumn1+4*(2+widthOtherColumns))<<endl;
      cout<<" "<<fixWidth("# of Z2 leptons from H",widthColumn1,true);
      for(int j=0; j<4; j++) cout<<"  "<<fixWidth(nbZDaughtersFromH[j],widthOtherColumns,false);
      cout<<endl;
      for(int i=0; i<4; i++){
	cout<<" ";
	cout<<fixWidth(ttHdecays[i],widthColumn1,true);
	int total = 0; for(int j=0; j<4; j++) total += nbZ2DaughtersFromHttH[i][j];
	for(int j=0; j<4; j++) cout<<"  "<<fixWidth(percentage((float)nbZ2DaughtersFromHttH[i][j]/total)+" %",widthOtherColumns,false);
	cout<<endl;
      }
    }
    
    cout<<endl;

  } // end for production modes



  // ------------------------------------------------------------
  // -------------------------- Plots ---------------------------
  // ------------------------------------------------------------

  if(doProdComp){
    TCanvas* cBCFullSel100[nVariables];
    for(int v=0; v<nVariables; v++){
      cBCFullSel100[v] = new TCanvas(Form("cBCFullSel100_%s",varName[v].c_str()),Form("cBCFullSel100_%s",varName[v].c_str()),500,500);
      DrawByProdmodes(cBCFullSel100[v],hBCFullSel100[v],prodName,isPresent,colors,v==0);
      SaveCanvas(outDir,cBCFullSel100[v]);
    }
  }

  if(doProdCompMatch4){
    TCanvas* cBCFullSel100Match4[nVariables];
    for(int v=0; v<nVariables; v++){
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
  
}


