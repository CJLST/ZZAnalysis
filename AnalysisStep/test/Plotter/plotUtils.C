#include "Math/DistFunc.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TStyle.h"


void setColZGradient_1() {
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}
void setColZGradient_2() {
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.20, 0.55, 0.88, 1.00 };
    Double_t red[NRGBs]   = { 0.35, 0.00, 0.87, 1.00, 0.70 };
    Double_t green[NRGBs] = { 0.30, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 1.00, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}
void setColZGradient_2_Compact() {
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.03, 0.12, 0.50, 1.00 };
    Double_t red[NRGBs]   = { 0.35, 0.00, 0.87, 1.00, 0.70 };
    Double_t green[NRGBs] = { 0.30, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 1.00, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}
void setColZGradient_Gray() {
    const Int_t NRGBs = 2;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 1.00 };
    Double_t red[NRGBs]   = { 0.90, 0.20 };
    Double_t green[NRGBs] = { 0.90, 0.20 };
    Double_t blue[NRGBs]  = { 0.90, 0.20 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}
void setColZGradient_Gray_Compact() {
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.03, 0.12, 0.50, 1.00 };
    Double_t red[NRGBs]   = { 0.90, 0.75, 0.60, 0.40, 0.20 };
    Double_t green[NRGBs] = { 0.90, 0.75, 0.60, 0.40, 0.20 };
    Double_t blue[NRGBs]  = { 0.90, 0.75, 0.60, 0.40, 0.20 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}
void setColZGradient_OneColor(int col) {
    const Int_t NRGBs = 2;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 1.00 };
    Double_t red[NRGBs]   = { 1., col==1?0.50:col==2?223./256.:col==3?255./256.:col==4?200./256.:col==5? 25./256.:col==6? 52./256.:col==7?116./256.:0.40 };
    Double_t green[NRGBs] = { 1., col==1?0.50:col==2? 48./256.:col==3?128./256.:col==4?167./256.:col==5?121./256.:col==6?182./256.:col==7?186./256.:0.40 };
    Double_t blue[NRGBs]  = { 1., col==1?0.50:col==2?164./256.:col==3?  0./256.:col==4?  0./256.:col==5?218./256.:col==6?152./256.:col==7?255./256.:0.40 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}


void printInfo(string info, Double_t x1, Double_t y1, Double_t x2, Double_t y2){
  TPaveText* pav = new TPaveText(x1,y1,x2,y2,"brNDC");
  pav->SetFillStyle(0);
  pav->SetBorderSize(0);
  pav->SetTextAlign(12);
  pav->AddText(info.c_str());
  pav->Draw();
}


void SaveCanvas(string directory, TCanvas* c, string tag = "") {
  c->SaveAs(Form("%s%s_%s.root",directory.c_str(),c->GetName(),tag.c_str()));
  //c->SaveAs(Form("%s%s_%s.C"   ,directory.c_str(),c->GetName(),tag.c_str())); // triggers a segfault !?
  c->SaveAs(Form("%s%s_%s.pdf" ,directory.c_str(),c->GetName(),tag.c_str()));  
  c->SaveAs(Form("%s%s_%s.png" ,directory.c_str(),c->GetName(),tag.c_str()));
}


TGraphAsymmErrors* getDataGraph(TH1F* h, bool drawZeroBins=false) {

  float fX[1000];
  float fY[1000];
  float fEXlow[1000];
  float fEXhigh[1000];
  float fEYlow[1000];
  float fEYhigh[1000];  

  TAxis *xaxis = ((TH1*)h)->GetXaxis();
  double q=(1-0.6827)/2.;

  int ibin=0;
  for (Int_t i=0; i<h->GetNbinsX(); ++i) {
    float yy = h->GetBinContent(i+1);
    if (drawZeroBins || yy > 0){
      float xx = xaxis->GetBinCenter(i+1);
      fX[ibin] = xx;
      fY[ibin] = yy;
      fEXlow [ibin] = 0.;
      fEXhigh[ibin] = 0.;
      double N = yy;
      fEYlow [ibin] = (N==0)?0:(N-ROOT::Math::chisquared_quantile_c(1-q,2*N)/2.);
      fEYhigh[ibin] = ROOT::Math::chisquared_quantile_c(q,2*(N+1))/2.-N;
      ++ibin;
    }
  }

  TGraphAsymmErrors* g = new TGraphAsymmErrors(ibin, fX, fY, fEXlow, fEXhigh, fEYlow, fEYhigh);
  h->TAttLine::Copy(*g);
  h->TAttFill::Copy(*g);
  h->TAttMarker::Copy(*g);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.9);
//   g->SetName(h->GetName());
//   g->SetTitle(h->GetTitle());

  return g;
}


TGraphAsymmErrors* getDataOverMCGraph(TGraphAsymmErrors* gData, TH1F* hMC) {
 
  int nPoints = gData->GetN();
  float fX[nPoints];
  float fY[nPoints];
  float fEXlow[nPoints];
  float fEXhigh[nPoints];
  float fEYlow[nPoints];
  float fEYhigh[nPoints];  

  int ipt=0;
  for(Int_t i=0; i<nPoints; ++i) {
    Float_t denom = hMC->GetBinContent(hMC->FindBin(gData->GetX()[i]));
    if(denom!=0.){
      fX[ipt] = gData->GetX()[i];
      fY[ipt] = gData->GetY()[i]/denom;
      fEXlow [ipt] = gData->GetEXlow ()[i];
      fEXhigh[ipt] = gData->GetEXhigh()[i];
      fEYlow [ipt] = gData->GetEYlow ()[i]/denom;
      fEYhigh[ipt] = gData->GetEYhigh()[i]/denom;
      ++ipt;
    }
  }

  TGraphAsymmErrors* gRatio = new TGraphAsymmErrors(ipt, fX, fY, fEXlow, fEXhigh, fEYlow, fEYhigh);
  gData->TAttLine::Copy(*gRatio);
  gData->TAttFill::Copy(*gRatio);
  gData->TAttMarker::Copy(*gRatio);

  return gRatio;
}


void restrictXAxis(TH1F* h, Int_t countmax){
  for(Int_t i=1; i<=h->GetNbinsX(); i++)
    h->GetXaxis()->SetBinLabel(i,Form("%i",(Int_t)(h->GetXaxis()->GetBinCenter(i))));
  h->GetXaxis()->SetBinLabel(countmax+1,Form("#geq%i",(Int_t)(h->GetXaxis()->GetBinCenter(countmax+1))));
  for(Int_t i=countmax+2; i<=h->GetNbinsX()+1; i++)
    h->SetBinContent(countmax+1,h->GetBinContent(countmax+1)+h->GetBinContent(i));
  h->GetXaxis()->SetLabelSize(0.06);
  h->GetXaxis()->SetNdivisions(-h->GetNbinsX());
  h->GetXaxis()->SetRangeUser(0., (float)(countmax+1));
}


TH1F* Smooth(TH1F* hIn, Int_t factor){
  TH1F* hOut = (TH1F*)hIn->Clone();
  TH1F* hInRebinned = (TH1F*)hIn->Rebin(factor,"hInRebinned");
  hInRebinned->Scale(1./(Float_t)factor);
  for(int ibin=1; ibin<=hIn->GetNbinsX(); ibin++)
    hOut->SetBinContent(ibin, hInRebinned->GetBinContent(1+(ibin-1)/factor));
  return hOut;
}


string ReplaceString(string subject, 
		     const string& search,
		     const string& replace) {
  if(search.empty()) return subject;
  size_t pos = 0;
  while((pos = subject.find(search, pos)) != std::string::npos) {
    subject.replace(pos, search.length(), replace);
    pos += replace.length();
  }
  return subject;
}
