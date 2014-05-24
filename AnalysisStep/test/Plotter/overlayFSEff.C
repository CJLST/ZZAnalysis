// Overlay different final states for efficiency plots, TO BE RUN in AnalysisInputs after having 
//   run the efficiency fits:
// root -b -q -l "overlayFSEff.C(\"ggH\",true,true)"
// parameters should be self explaining
// output is stored in sigFigs directory


#include "TFile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"
#include "TStyle.h"

void overlayFSEff(TString sprocess = "ggH", Bool_t is8TeV = true, Bool_t isTot = true)
{

  if (!(sprocess=="ggH" || sprocess=="qqH" || sprocess=="ZH" || sprocess=="WH" || sprocess=="ttH")) {
    std::cout<<"Invalid process string, choose ggH, qqH, WH, ZH, ttH"<<std::endl;
    exit(1);
  }

  TString baseDir = "sigFigs"+(is8TeV?TString("8TeV"):TString("7TeV"));
  TString ratioOrNot = isTot?"":"_ratio";
  TString sratio = isTot?"TotalEfficiency":"Ratio";
  TString sprocessExt;

  if (sprocess=="ggH") sprocessExt = "gluon-gluon fusion";
  else if (sprocess=="qqH") sprocessExt = "vector boson fusion"; 
  else if (sprocess=="ZH") sprocessExt = "associated production ZH";
  else if (sprocess=="WH") sprocessExt = "associated production WH";
  else sprocessExt = "associated production ttH";
  
  TFile* f_2e2mu = TFile::Open(baseDir+"/eff_"+sprocess+"_2e2mu"+ratioOrNot+".root","READ");
  TFile* f_4e = TFile::Open(baseDir+"/eff_"+sprocess+"_4e"+ratioOrNot+".root","READ");
  TFile* f_4mu = TFile::Open(baseDir+"/eff_"+sprocess+"_4mu"+ratioOrNot+".root","READ");

  TGraphErrors* tge_2e2mu = (TGraphErrors*)f_2e2mu->Get(sratio);
  TGraphErrors* tge_4mu = (TGraphErrors*)f_4mu->Get(sratio);
  TGraphErrors* tge_4e = (TGraphErrors*)f_4e->Get(sratio);

  tge_2e2mu->SetMarkerStyle(21); 
  tge_4mu->SetMarkerStyle(22);
  tge_4e->SetMarkerStyle(20);
  tge_2e2mu->SetMarkerColor(4); 
  tge_4mu->SetMarkerColor(2);
  tge_4e->SetMarkerColor(3);
  tge_2e2mu->SetLineColor(4); 
  tge_4mu->SetLineColor(2);
  tge_4e->SetLineColor(3);

  TF1* func_2e2mu, *func_4e, *func_4mu;

  if (isTot) {
    func_2e2mu = tge_2e2mu->GetFunction("polyFunctot");
    func_4mu   = tge_4mu->GetFunction("polyFunctot");
    func_4e    = tge_4e->GetFunction("polyFunctot");
  } else {
    func_2e2mu = tge_2e2mu->GetFunction("ratiofit");
    func_4mu   = tge_4mu->GetFunction("ratiofit");
    func_4e    = tge_4e->GetFunction("ratiofit");
  }

  func_2e2mu->SetLineColor(4);
  func_4mu->SetLineColor(2);
  func_4e->SetLineColor(3);
    
  TMultiGraph* mg_eff = new TMultiGraph("mg_eff","Total Efficiency "+sprocess+(is8TeV?" 8TeV":" 7TeV"));
  mg_eff->Add(tge_2e2mu); mg_eff->Add(tge_4mu); mg_eff->Add(tge_4e);

  TCanvas* c_eff = new TCanvas("c_eff","c_eff");
  c_eff->cd();

  mg_eff->Draw("AP");
  if (isTot) mg_eff->SetTitle("CMS Simulation, #sqrt{s} = "+(is8TeV?TString("8TeV"):TString("7TeV")));
  else mg_eff->SetTitle("Dijet ratio "+sprocess);
  mg_eff->GetXaxis()->SetTitleFont(42);
  mg_eff->GetXaxis()->SetTitleSize(0.045);
  mg_eff->GetYaxis()->SetTitleFont(42);
  mg_eff->GetYaxis()->SetTitleSize(0.045);
  mg_eff->GetXaxis()->SetTitle("m_{H} [GeV/c^{2}]");
  if (isTot) mg_eff->GetYaxis()->SetTitle("Efficiency");
  else mg_eff->GetYaxis()->SetTitle("Dijet ratio");

  mg_eff->GetYaxis()->SetRangeUser(0.,1.);


  TLegend* leg = new TLegend(0.35,0.20,0.80,0.45);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  
  TString legtitle = "H#rightarrow ZZ* #rightarrow 4l";
  leg->AddEntry((TObject*)0, legtitle.Data(), "");
  leg->AddEntry((TObject*)0, sprocessExt.Data(), "");
  leg->AddEntry(tge_2e2mu,"2e2#mu","PL");
  leg->AddEntry(tge_4mu,  "4#mu","PL");
  leg->AddEntry(tge_4e,  "4e","PL");

  leg->Draw("same");

  c_eff->SaveAs(baseDir+"/Efficiency_"+sprocess+(is8TeV?"_8TeV"+ratioOrNot+".eps":"_7TeV"+ratioOrNot+".eps"));
  c_eff->SaveAs(baseDir+"/Efficiency_"+sprocess+(is8TeV?"_8TeV"+ratioOrNot+".root":"_7TeV"+ratioOrNot+".root"));

  


}
