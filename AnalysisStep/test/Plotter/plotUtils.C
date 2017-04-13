#include "Math/DistFunc.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TStyle.h"


void setColZGradient_Rainbow1() {
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}
void setColZGradient_Rainbow2() {
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.20, 0.55, 0.88, 1.00 };
    Double_t red[NRGBs]   = { 0.35, 0.00, 0.87, 1.00, 0.70 };
    Double_t green[NRGBs] = { 0.30, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 1.00, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}
void setColZGradient_Rainbow2_Compact() {
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
void setColZGradient_OneColor(int col, bool shift = false) {
    const Int_t NRGBs = 2;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 1.00 };
    Double_t red[NRGBs]   = { 1., col==1?0.50:col==2?223./256.:col==3?255./256.:col==4?200./256.:col==5? 25./256.:col==6? 52./256.:col==7?116./256.:0.40 };
    Double_t green[NRGBs] = { 1., col==1?0.50:col==2? 48./256.:col==3?128./256.:col==4?167./256.:col==5?121./256.:col==6?182./256.:col==7?186./256.:0.40 };
    Double_t blue[NRGBs]  = { 1., col==1?0.50:col==2?164./256.:col==3?  0./256.:col==4?  0./256.:col==5?218./256.:col==6?152./256.:col==7?255./256.:0.40 };
    if(shift){
      for(int i=0; i<NRGBs; i++){
	red  [i] -= 0.10;
	green[i] -= 0.10;
	blue [i] -= 0.10;
      }
    }
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}
void setColZGradient_TwoColors() {
    const Int_t NRGBs = 2;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 1.00 };
    Double_t red[NRGBs]   = { 60./256., 1.00 };
    Double_t green[NRGBs] = { 140./256., 1.00 };
    Double_t blue[NRGBs]  = { 1.00, 0.50 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}


TPaveText* printInfo(string info, Double_t x1, Double_t y1, Double_t x2, Double_t y2){
  TPaveText* pav = new TPaveText(x1,y1,x2,y2,"brNDC");
  pav->SetFillStyle(0);
  pav->SetBorderSize(0);
  pav->SetTextAlign(12);
  pav->AddText(info.c_str());
  pav->Draw();
  return pav;
}


void SaveCanvas(string directory, TCanvas* c, string tag = "") {
  string prefix = directory + string(c->GetName()) + "_" + tag;
  c->SaveAs((prefix+".root").c_str());
  c->SaveAs((prefix+".C").c_str()); //sometimes causes a segfault, apprently due to empty TGraphs
  c->SaveAs((prefix+".pdf").c_str());
  //c->SaveAs((prefix+".png").c_str()); //gives bad quality, it's much better to produce a .png from the .eps by doing: convert -density 150 -quality 100 FILE.eps FILE.png
  c->SaveAs((prefix+".eps").c_str());
  gSystem->Exec(("convert -density 150 -quality 100 "+prefix+".eps "+prefix+".png").c_str());
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


void NormalizePerSlice(TH2F* h, Bool_t sliceInX = true){
  Int_t nbinsx = h->GetNbinsX();
  Int_t nbinsy = h->GetNbinsY();
  if(sliceInX)
    for(Int_t bx=0; bx<=nbinsx+1; bx++){
      Float_t sliceNorm = h->Integral(bx,bx,0,nbinsy+1);
      for(Int_t by=0; by<=nbinsy+1; by++)
	h->SetBinContent(bx,by,h->GetBinContent(bx,by)/sliceNorm);
    }
  else
    for(Int_t by=0; by<=nbinsy+1; by++){
      Float_t sliceNorm = h->Integral(0,nbinsx+1,by,by);
      for(Int_t bx=0; bx<=nbinsx+1; bx++)
	h->SetBinContent(bx,by,h->GetBinContent(bx,by)/sliceNorm);
    }
  TAxis* axisToChange = sliceInX ? h->GetXaxis() : h->GetYaxis() ;
  axisToChange->SetTitle(Form("%s (normalized per slice)",axisToChange->GetTitle()));
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


TGraph* doROC (
	       TH1F* hSgnl,
	       TH1F* hBkgd,
	       Bool_t inverted = false,
	       Bool_t un = false, // include underflow bin
	       Bool_t ov = false, // include overflow bin
	       string title = "",
	       string titleSgnl = "signal efficiency",
	       string titleBkgd = "background efficiency"
	       )
{
  if (hSgnl!=0 && hBkgd!=0) {

    Int_t nbBinsSgnl = hSgnl->GetNbinsX();
    Int_t nbBinsBkgd = hBkgd->GetNbinsX(); 
  
    if (nbBinsSgnl==nbBinsBkgd) {

      Int_t nbBins = nbBinsSgnl;
      if(hSgnl->GetBinContent(0)!=0)
	cout<<"Warning in function doROC : underflow bin of histogram "<<hSgnl->GetName()<<" is not empty ("<<(un?"":"not ")<<"including it)"<<endl;
      if(hBkgd->GetBinContent(0)!=0)
	cout<<"Warning in function doROC : underflow bin of histogram "<<hBkgd->GetName()<<" is not empty ("<<(un?"":"not ")<<"including it)"<<endl;
      if(hSgnl->GetBinContent(nbBins+1)!=0)
	cout<<"Warning in function doROC : overflow bin of histogram "<<hSgnl->GetName()<<" is not empty ("<<(ov?"":"not ")<<"including it)"<<endl;
      if(hBkgd->GetBinContent(nbBins+1)!=0)
	cout<<"Warning in function doROC : overflow bin of histogram "<<hBkgd->GetName()<<" is not empty ("<<(ov?"":"not ")<<"including it)"<<endl;

      Double_t vSgnl[nbBins+1+un+ov];
      Double_t vBkgd[nbBins+1+un+ov];
      Double_t intSgnl = hSgnl->Integral(1-un,nbBins+ov);
      Double_t intBkgd = hBkgd->Integral(1-un,nbBins+ov);

      vSgnl[0] = 0;
      vBkgd[0] = 0;
      for(Int_t i=1-un; i<=nbBins+ov; i++){
	//vSgnl[i+un] = hSgnl->Integral(1-un,i)/intSgnl; //if(i==1-un)cout<<vSgnl[i+un]<<endl;
	//vBkgd[i+un] = hBkgd->Integral(1-un,i)/intBkgd; //if(i==1-un)cout<<vBkgd[i+un]<<endl;
	vSgnl[i+un] = vSgnl[i-1+un] + hSgnl->GetBinContent(i)/intSgnl;
	vBkgd[i+un] = vBkgd[i-1+un] + hBkgd->GetBinContent(i)/intBkgd;
      }

      if(inverted){
	for(Int_t i=0; i<=nbBins+1+un+ov; i++){
	  vSgnl[i] = 1-vSgnl[i];
	  vBkgd[i] = 1-vBkgd[i];
	}
      }

      TGraph* gRoc = new TGraph(nbBins+1+un+ov,vBkgd,vSgnl);  
      gRoc->SetTitle(title.c_str());
      gRoc->GetXaxis()->SetTitle(titleBkgd.c_str());
      gRoc->GetYaxis()->SetTitle(titleSgnl.c_str());
      gRoc->GetXaxis()->SetLimits(0.,1.);
      gRoc->SetMinimum(0);
      gRoc->SetMaximum(1);
      
      return gRoc; 
      
    } else {
      cout << Form("Function doROC has raised an error : the histograms %s and %s don't have the same number of bins",hSgnl->GetName(),hBkgd->GetName()) << endl;
      return new TGraph();
    }

  } else {
    cout << Form("Function doROC has raised an error : histogram is not allocated (%s or %s)",hSgnl->GetName(),hBkgd->GetName()) << endl;
    return new TGraph();
  }
}


TGraph* doEFF (
	       TH1F* h,
	       Bool_t inverted = false,
	       string titleYaxis = "signal efficiency",
	       Bool_t un = false, // include underflow bin
	       Bool_t ov = false, // include overflow bin
	       Float_t minXaxis = -999.,
	       Float_t maxXaxis = -999.,
	       string title = "",
	       string titleXaxis = "working point"
	       )
{
  if (h!=0) {

    Int_t nbBins = h->GetNbinsX();
  
    if(h->GetBinContent(0)!=0)
      cout<<"Warning in function doEFF : underflow bin of histogram "<<h->GetName()<<" is not empty ("<<(un?"":"not ")<<"including it)"<<endl;
    if(h->GetBinContent(nbBins+1)!=0)
      cout<<"Warning in function doEFF : overflow bin of histogram "<<h->GetName()<<" is not empty ("<<(ov?"":"not ")<<"including it)"<<endl;

    Double_t thr[nbBins+1];
    Double_t eff[nbBins+1];
    Double_t integral = h->Integral(1-un,nbBins+ov);

    eff[0] = un * h->GetBinContent(0)/integral;
    thr[0] = h->GetBinLowEdge(1);
    for(Int_t i=1; i<=nbBins; i++){
      eff[i] = eff[i-1] + h->GetBinContent(i)/integral;
      thr[i] = h->GetBinLowEdge(i+1);
    }  

    if(inverted){
      for(Int_t i=0; i<=nbBins+1; i++){
	eff[i] = 1-eff[i];
      }
    }

    TGraph* gRoc = new TGraph(nbBins+1,thr,eff);  
    gRoc->SetTitle(title.c_str());
    gRoc->GetXaxis()->SetTitle(titleXaxis.c_str());
    gRoc->GetYaxis()->SetTitle(titleYaxis.c_str());
    gRoc->GetXaxis()->SetLimits((minXaxis!=-999.?minXaxis:h->GetBinLowEdge(2-un)),(maxXaxis!=-999.?maxXaxis:h->GetBinLowEdge(nbBins+ov+1)));
    gRoc->SetMinimum(0);
    gRoc->SetMaximum(1);
      
    return gRoc; 
      
  } else {
    cout << Form("Function doEFF has raised an error : histogram is not allocated (%s)",h->GetName()) << endl;
    return new TGraph();
  }
}


TGraph* doWP (
	      TH1F* hSgnl,
	      TH1F* hBkgd, 
	      Double_t WP,
	      Bool_t inverted = false,
	      Bool_t un = false, // include underflow bin
	      Bool_t ov = false // include overflow bin
	      )
{  
  if (hSgnl!=0 && hBkgd!=0) {

    Int_t nbBinsSgnl = hSgnl->GetNbinsX();
    Int_t nbBinsBkgd = hBkgd->GetNbinsX(); 
  
    if (nbBinsSgnl==nbBinsBkgd) {

      Int_t nbBins = nbBinsSgnl;
      if(hSgnl->GetBinContent(0)!=0)
	cout<<"Warning in function doWP : underflow bin of histogram "<<hSgnl->GetName()<<" is not empty ("<<(un?"":"not ")<<"including it)"<<endl;
      if(hBkgd->GetBinContent(0)!=0)
	cout<<"Warning in function doWP : underflow bin of histogram "<<hBkgd->GetName()<<" is not empty ("<<(un?"":"not ")<<"including it)"<<endl;
      if(hSgnl->GetBinContent(nbBins+1)!=0)
	cout<<"Warning in function doWP : overflow bin of histogram "<<hSgnl->GetName()<<" is not empty ("<<(ov?"":"not ")<<"including it)"<<endl;
      if(hBkgd->GetBinContent(nbBins+1)!=0)
	cout<<"Warning in function doWP : overflow bin of histogram "<<hBkgd->GetName()<<" is not empty ("<<(ov?"":"not ")<<"including it)"<<endl;

      Double_t effSgnl[1];
      Double_t effBkgd[1];
      Double_t intSgnl = hSgnl->Integral(1-un,nbBins+ov);
      Double_t intBkgd = hBkgd->Integral(1-un,nbBins+ov);

      for(Int_t i=1-un; i<=nbBins+ov; i++){
	if((Float_t)hSgnl->GetBinLowEdge(i+1)<=(Float_t)WP && (Float_t)WP<(Float_t)hSgnl->GetBinLowEdge(i+2)){ // casting is essential here !
	  effSgnl[0] = hSgnl->Integral(1-un,i)/intSgnl;
	  effBkgd[0] = hBkgd->Integral(1-un,i)/intBkgd;
	  if(inverted){
	    effSgnl[0] = 1.-effSgnl[0];
	    effBkgd[0] = 1.-effBkgd[0];
	  }
	  break;
	}
      }    

      TGraph* gRoc = new TGraph(1,effBkgd,effSgnl);  

      return gRoc; 
      
    } else {
      cout << Form("Function doWP has raised an error : the histograms %s and %s don't have the same number of bins",hSgnl->GetName(),hBkgd->GetName()) << endl;
      return new TGraph();
    }

  } else {
    cout << Form("Function doWP has raised an error : histogram is not allocated (%s or %s)",hSgnl->GetName(),hBkgd->GetName()) << endl;
    return new TGraph();
  }
}
