#include "TDirectory.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "THStack.h"
#include "TH2.h"
#include  <TF1.h>
#include "TChain.h"
#include "TLine.h"
#include "TCut.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"
#include "ZZAnalysis/AnalysisStep/interface/Category.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"
#include "TSpline.h"

#include <cmath>
#include "TMath.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooWorkspace.h"
#include "RooHist.h"
#include "RooChebychev.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "TPaveText.h"
#include "RooAddPdf.h"
#include "RooBreitWigner.h"
#include "RooFitResult.h"
#include "RooFFTConvPdf.h"
#include "RooAddition.h"
#include "RooMinuit.h"
#include "Math/MinimizerOptions.h"
#include <iomanip>
#include "RooAbsCollection.h"

using namespace RooFit;


#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>


int Wait() {
     cout << " Continue [<RET>|q]?  "; 
     char x;
     x = getchar();
     if ((x == 'q') || (x == 'Q')) return 1;
     return 0;
}

Double_t effSigma(TH1 *hist );
Double_t effSigma(RooAbsPdf *pdf, RooRealVar *obs, Int_t nbins);
float getFitEdge(float mass, float width, bool low);
float getFitEdgeHighMass(float mass, float width, bool low);

string schannel;
string scategory;
string ssample;




void fitZZShapeW(int ch, int cat,  int sample,
		     double rangeLow, double rangeHigh,
		     double fitValues[7], double fitErrors[7], double covQual[1]);

bool do7TeV;
bool doFFT;
// sample 1 = ggH , sample 2 = VBFH 
void all(int channels=0, int categ=-1, int sample = 0 )/*,bool doSfLepton=true)*/{


  if (channels == 0) schannel = "4mu";
  if (channels == 1) schannel = "4e";
  if (channels == 2) schannel = "2e2mu";

  if (categ < 0)  scategory = "All";
  if (categ == 0) scategory = "Untagged";
  if (categ == 1) scategory = "VBFtagged";
/*  if (categ == 2) scategory = "VBFtagged";
  if (categ == 3) scategory = "VHleptonictagged";
  if (categ == 4) scategory = "VHHadronictagged";
  if (categ == 5) scategory = "ttHtagged";
*/
  if (sample ==1) ssample = "qqZZ";
  if (sample ==2) ssample = "ggZZ";

  double fitValues[10];
  double fitErrors[10];
  double covQual[1];

  fitZZShapeW(channels,categ, sample,105., 165.,
		    fitValues,fitErrors,covQual);  
  
  cout << "pol1 value,error: " << fitValues[0] << " , " << fitErrors[0] << endl; 
  
  cout << "pol2 value,error: " << fitValues[1] << " , " << fitErrors[1] << endl; 
  
  cout << "pol3 value,error: " << fitValues[2] << " , " << fitErrors[2] << endl; 
  
  cout << "covQual of the fit: " << covQual[0] << endl;
  
  
  string filename = "bkg_shape_parametriztion_13TeV_" + ssample + "_" + schannel + "_" + scategory + "." + "yaml" ;
  
  ofstream outFile;
  outFile.open(filename);
  outFile <<"Sample   : "<< ssample << endl;
  outFile <<"Channel  : "<< schannel << endl;
  outFile <<"Category : "<< scategory <<endl;
  outFile <<"    chebPol1 : " <<"'"<<fitValues[0]<<"'"<<endl;
  outFile <<"    chebPol2 : " <<"'"<<fitValues[1]<<"'"<<endl;
  outFile <<"    chebPol3 : " <<"'"<<fitValues[2]<<"'"<<endl;
  outFile << endl;
  
}


void fitZZShapeW(int channels,int categ, int sample, 
		     double rangeLow, double rangeHigh,
		     double fitValues[7], double fitErrors[7], double covQual[1]){
 // ------ root settings ---------
  gROOT->Reset();  
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kFALSE);
  gStyle->SetPadGridY(kFALSE);
  //gStyle->SetOptStat("kKsSiourRmMen");
  gStyle->SetOptStat("iourme");
  //gStyle->SetOptStat("rme");
  //gStyle->SetOptStat("");
  gStyle->SetOptFit(11);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  // ------------------------------ 

  string scategory;
  string ssample;
  
  if (sample ==1) ssample = "qqZZ";
  if (sample ==2) ssample = "ggZZ";
   
  if (categ == 0  ) scategory = "Untagged";
  if (categ == 1  ) scategory = "VBFtagged";
  
/*  if (categ < 0)  scategory = "All";
  if (categ == 0) scategory = "Untagged";
  if (categ == 1) scategory = "1Jettagged";
  if (categ == 2) scategory = "VBFtagged";
  if (categ == 3) scategory = "VHleptonictagged";
  if (categ == 4) scategory = "VHHadronictagged";
  if (categ == 5) scategory = "ttHtagged";
*/
  ROOT::Math::MinimizerOptions::SetDefaultTolerance( 1.E-7);
  TFile * kfact = new TFile( "../Plotter/ggZZKfactors/Kfactor_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");  
  stringstream FileName[6];
  int howmany = 1;
  if(sample==1) FileName[0] << "root://lxcms03//data3/Higgs/160121/ZZTo4l/ZZ4lAnalysis.root";
  else if(sample==2) {
    FileName[0] << "root://lxcms03//data3/Higgs/150915/ggZZ4e/ZZ4lAnalysis.root";
    FileName[1] << "root://lxcms03//data3/Higgs/150915/ggZZ4mu/ZZ4lAnalysis.root";
    FileName[2] << "root://lxcms03//data3/Higgs/150915/ggZZ4tau/ZZ4lAnalysis.root";
    FileName[3] << "root://lxcms03//data3/Higgs/160121/ggZZ2e2mu/ZZ4lAnalysis.root";
    FileName[4] << "root://lxcms03//data3/Higgs/160121/ggZZ2mu2tau/ZZ4lAnalysis.root";
    FileName[5] << "root://lxcms03//data3/Higgs/160121/ggZZ2e2tau/ZZ4lAnalysis.root";
    howmany=6;
  }
  
  else {
    cout << "Wrong sample ." << endl;
    return;
  }
  
 
   

  TSpline3 * f1 = (TSpline3*)kfact->Get("sp_Kfactor");  
  
  TChain *ggTree = new TChain("ZZTree/candTree");

  for (int ifil = 0; ifil < howmany; ifil++) {
    cout << "Using " << FileName[ifil].str() << endl;
    ggTree->Add(FileName[ifil].str().c_str());
  } 

    // TTree* ggTree = (TTree*) ggFile->Get("ZZTree/candTree");

  float m4l;
  
  Short_t z1flav, z2flav; 
  float weight;
 
  Short_t nExtraLeptons;   
  float ZZPt,PVBF_VAJHU_old, PHJJ_VAJHU_old,GenHMass, GenZ1Phi,GenZ2Phi;
  Short_t nJets;
  Short_t nBTaggedJets;
  std::vector<float> * jetpt = 0;
  std::vector<float> * jeteta = 0;
  std::vector<float> * jetphi = 0;
  std::vector<float> * jetmass = 0;
  float jet30pt[10];
  float jet30eta[10];
  float jet30phi[10];
  float jet30mass[10];
  float Fisher;
  
  int  nentries = ggTree->GetEntries();
 
  //--- ggTree part
  ggTree->SetBranchAddress("ZZMass",&m4l);
  ggTree->SetBranchAddress("ZZMass",&m4l);
  ggTree->SetBranchAddress("Z1Flav",&z1flav);
  ggTree->SetBranchAddress("Z2Flav",&z2flav);
  ggTree->SetBranchAddress("genHEPMCweight",&weight);
  ggTree->SetBranchAddress("nExtraLep",&nExtraLeptons);
  ggTree->SetBranchAddress("nCleanedJets",&nJets);
  ggTree->SetBranchAddress("nCleanedJetsPt30BTagged",&nBTaggedJets);
  ggTree->SetBranchAddress("DiJetFisher",&Fisher);
  ggTree->SetBranchAddress("pvbf_VAJHU_highestPTJets",&PVBF_VAJHU_old);
  ggTree->SetBranchAddress("phjj_VAJHU_highestPTJets",&PHJJ_VAJHU_old);
  if(sample ==1){
  ggTree->SetBranchAddress("GenZ1Phi",&GenZ1Phi);
  ggTree->SetBranchAddress("GenZ2Phi",&GenZ2Phi);}
 
  ggTree->SetBranchAddress("JetPt",&jetpt);
  ggTree->SetBranchAddress("JetEta",&jeteta);
  ggTree->SetBranchAddress("JetPhi",&jetphi);
  ggTree->SetBranchAddress("JetMass",&jetmass);
  ggTree->SetBranchAddress("ZZPt",&ZZPt);

  ggTree->SetBranchAddress("GenHMass",&GenHMass);

  //--- rooFit part
  double xMin,xMax;
  // xInit = (double) massBin;
  xMin = rangeLow;
  xMax = rangeHigh ;
  cout << "Fit range: [" << xMin << " , " << xMax << "]" << endl;
  TH1F *hmass = new TH1F("hmass","hmass",200,xMin,xMax);
  //---------  
  
  RooRealVar x("mass","m_{4l}",140.,xMin,xMax,"GeV");
  RooRealVar w("myW","myW",1.0,0.,1000.);
 
 float newWeight = 0; 
 float GENdPhiZZ = 0;
 
  RooArgSet ntupleVarSet(x,w);
  RooDataSet dataset("mass4l","mass4l",ntupleVarSet,WeightVar("myW"));
  
 for(int k=0; k<nentries; k++){
    ggTree->GetEvent(k);

    int njet30 = 0;
    for (unsigned int ijet = 0; ijet < jetpt->size(); ijet++) { 
      if ( (*jetpt)[ijet] > 30. ) {
	jet30pt[njet30] = (*jetpt)[ijet];      
	jet30eta[njet30] = (*jeteta)[ijet];
	jet30phi[njet30] = (*jetphi)[ijet];
	jet30mass[njet30] = (*jetmass)[ijet];
	njet30++;
      }
    }  
//    int Cat = category(nExtraLeptons, ZZPt, m4l, njet30, nBTaggedJets, jet30pt, jet30eta, jet30phi,jet30mass, Fisher); 
    int Cat = categoryMor16(nJets, PVBF_VAJHU_old, PHJJ_VAJHU_old );
    if (categ >= 0 && categ != Cat ) continue;
    
    if(channels==0 && z1flav*z2flav != 28561) continue;
    if(channels==1 && z1flav*z2flav != 14641) continue;
    if (weight <= 0 ) cout << "Warning! Negative weight events" << endl;
    if(channels==2 && z1flav*z2flav != 20449) continue;
    
    if (sample ==1){     
    GENdPhiZZ = TMath::Abs(GenZ1Phi-GenZ2Phi); 
    if(GENdPhiZZ > TMath::Pi()) GENdPhiZZ = TMath::TwoPi() - GENdPhiZZ; 
    if (channels <= 1)(newWeight = 1.51583892176*(abs(GENdPhiZZ)>0.0&&abs(GENdPhiZZ)<=0.1)+1.49625666541*(abs(GENdPhiZZ)>0.1&&abs(GENdPhiZZ)<=0.2)+1.49552206191*(abs(GENdPhiZZ)>0.2&&abs(GENdPhiZZ)<=0.3)+1.48327315425*(abs(GENdPhiZZ)>0.3&&abs(GENdPhiZZ)<=0.4)+1.46558970113*(abs(GENdPhiZZ)>0.4&&abs(GENdPhiZZ)<=0.5)+1.49150088751*(abs(GENdPhiZZ)>0.5&&abs(GENdPhiZZ)<=0.6)+1.44118358045*(abs(GENdPhiZZ)>0.6&&abs(GENdPhiZZ)<=0.7)+1.44083060399*(abs(GENdPhiZZ)>0.7&&abs(GENdPhiZZ)<=0.8)+1.41433901912*(abs(GENdPhiZZ)>0.8&&abs(GENdPhiZZ)<=0.9)+1.42253421856*(abs(GENdPhiZZ)>0.9&&abs(GENdPhiZZ)<=1.0)+1.401037066*(abs(GENdPhiZZ)>1.0&&abs(GENdPhiZZ)<=1.1)+1.40853942881*(abs(GENdPhiZZ)>1.1&&abs(GENdPhiZZ)<=1.2)+1.38124774408*(abs(GENdPhiZZ)>1.2&&abs(GENdPhiZZ)<=1.3)+1.37055335743*(abs(GENdPhiZZ)>1.3&&abs(GENdPhiZZ)<=1.4)+1.347323316*(abs(GENdPhiZZ)>1.4&&abs(GENdPhiZZ)<=1.5)+1.34011343745*(abs(GENdPhiZZ)>1.5&&abs(GENdPhiZZ)<=1.6)+1.31266103651*(abs(GENdPhiZZ)>1.6&&abs(GENdPhiZZ)<=1.7)+1.29005506201*(abs(GENdPhiZZ)>1.7&&abs(GENdPhiZZ)<=1.8)+1.25532261479*(abs(GENdPhiZZ)>1.8&&abs(GENdPhiZZ)<=1.9)+1.25445564245*(abs(GENdPhiZZ)>1.9&&abs(GENdPhiZZ)<=2.0)+1.22404766442*(abs(GENdPhiZZ)>2.0&&abs(GENdPhiZZ)<=2.1)+1.17881678267*(abs(GENdPhiZZ)>2.1&&abs(GENdPhiZZ)<=2.2)+1.16262482714*(abs(GENdPhiZZ)>2.2&&abs(GENdPhiZZ)<=2.3)+1.10540114094*(abs(GENdPhiZZ)>2.3&&abs(GENdPhiZZ)<=2.4)+1.07474926569*(abs(GENdPhiZZ)>2.4&&abs(GENdPhiZZ)<=2.5)+1.02186459938*(abs(GENdPhiZZ)>2.5&&abs(GENdPhiZZ)<=2.6)+0.946334793286*(abs(GENdPhiZZ)>2.6&&abs(GENdPhiZZ)<=2.7)+0.857458082628*(abs(GENdPhiZZ)>2.7&&abs(GENdPhiZZ)<=2.8)+0.716607670482*(abs(GENdPhiZZ)>2.8&&abs(GENdPhiZZ)<=2.9)+1.13284178484*(abs(GENdPhiZZ)>2.9&&abs(GENdPhiZZ)<=3.1416));  
    
    if(channels ==2 )(newWeight = 1.51383448915*(abs(GENdPhiZZ)>0.0&&abs(GENdPhiZZ)<=0.1)+1.54173878018*(abs(GENdPhiZZ)>0.1&&abs(GENdPhiZZ)<=0.2)+1.49782963251*(abs(GENdPhiZZ)>0.2&&abs(GENdPhiZZ)<=0.3)+1.53495678292*(abs(GENdPhiZZ)>0.3&&abs(GENdPhiZZ)<=0.4)+1.47821703306*(abs(GENdPhiZZ)>0.4&&abs(GENdPhiZZ)<=0.5)+1.50433085929*(abs(GENdPhiZZ)>0.5&&abs(GENdPhiZZ)<=0.6)+1.52062624685*(abs(GENdPhiZZ)>0.6&&abs(GENdPhiZZ)<=0.7)+1.50701309003*(abs(GENdPhiZZ)>0.7&&abs(GENdPhiZZ)<=0.8)+1.49424315625*(abs(GENdPhiZZ)>0.8&&abs(GENdPhiZZ)<=0.9)+1.45053609615*(abs(GENdPhiZZ)>0.9&&abs(GENdPhiZZ)<=1.0)+1.46081252166*(abs(GENdPhiZZ)>1.0&&abs(GENdPhiZZ)<=1.1)+1.4716036222*(abs(GENdPhiZZ)>1.1&&abs(GENdPhiZZ)<=1.2)+1.4677000382*(abs(GENdPhiZZ)>1.2&&abs(GENdPhiZZ)<=1.3)+1.42240869064*(abs(GENdPhiZZ)>1.3&&abs(GENdPhiZZ)<=1.4)+1.39718402273*(abs(GENdPhiZZ)>1.4&&abs(GENdPhiZZ)<=1.5)+1.37559344752*(abs(GENdPhiZZ)>1.5&&abs(GENdPhiZZ)<=1.6)+1.39190131837*(abs(GENdPhiZZ)>1.6&&abs(GENdPhiZZ)<=1.7)+1.36856435056*(abs(GENdPhiZZ)>1.7&&abs(GENdPhiZZ)<=1.8)+1.31788480429*(abs(GENdPhiZZ)>1.8&&abs(GENdPhiZZ)<=1.9)+1.3140199508*(abs(GENdPhiZZ)>1.9&&abs(GENdPhiZZ)<=2.0)+1.27464174991*(abs(GENdPhiZZ)>2.0&&abs(GENdPhiZZ)<=2.1)+1.24234660682*(abs(GENdPhiZZ)>2.1&&abs(GENdPhiZZ)<=2.2)+1.24472740384*(abs(GENdPhiZZ)>2.2&&abs(GENdPhiZZ)<=2.3)+1.14625935167*(abs(GENdPhiZZ)>2.3&&abs(GENdPhiZZ)<=2.4)+1.10780499352*(abs(GENdPhiZZ)>2.4&&abs(GENdPhiZZ)<=2.5)+1.04205364674*(abs(GENdPhiZZ)>2.5&&abs(GENdPhiZZ)<=2.6)+0.973608545141*(abs(GENdPhiZZ)>2.6&&abs(GENdPhiZZ)<=2.7)+0.872169942668*(abs(GENdPhiZZ)>2.7&&abs(GENdPhiZZ)<=2.8)+0.734505279177*(abs(GENdPhiZZ)>2.8&&abs(GENdPhiZZ)<=2.9)+1.16315283723*(abs(GENdPhiZZ)>2.9&&abs(GENdPhiZZ)<=3.1416));
    }
    if (sample == 2)
    newWeight=f1->Eval(GenHMass);

    weight = weight * newWeight;
    ntupleVarSet.setRealValue("mass",m4l);
    ntupleVarSet.setRealValue("myW",weight);
    if(x.getVal()>xMin && x.getVal()<xMax)
      dataset.add(ntupleVarSet, weight);
    hmass->Fill(m4l);

  }
  //---------

  cout << "dataset n entries: " << dataset.sumEntries() << endl;


  TCanvas *c1 = new TCanvas("c1","c1",725,725);

  TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.35,0.95,0.97);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.02,0.95,0.35);
  pad2->Draw();


  //--- double CrystalBall
  //Chebyshev-Polynomial
  RooRealVar A1("A1","A1",-1,-3,3.);
  RooRealVar A2("A2","A2",0.5,-3.,3.);
  RooRealVar A3("A3","A3",0.,-3.,3.);
  RooChebychev model("model","model",x ,RooArgList(A1,A2,A3));
 
  RooFitResult *fitres = (RooFitResult*)model.fitTo(dataset,SumW2Error(1),Range(xMin,xMax),Strategy(2),NumCPU(8),Save(true));
  
  stringstream frameTitle;
  if(channels==0){frameTitle << "4#mu; "; }
  if(channels==1){frameTitle << "4e; ";}
  if(channels==2){frameTitle << "2e2#mu; ";}
  frameTitle << ssample << "; " << scategory;

  stringstream nameFileRoot;
  nameFileRoot << "fitM" << ssample << ".root";
  TFile *fileplot = TFile::Open(nameFileRoot.str().c_str(), "recreate");

  RooPlot* xframe = x.frame() ;
  xframe->SetTitle("");
  xframe->SetName("m4lplot");
  dataset.plotOn(xframe,DataError(RooAbsData::SumW2), MarkerStyle(kOpenCircle), MarkerSize(1.1) );
  int col;
  if(channels==0) col=kOrange+7;
  if(channels==1) col=kAzure+2;
  if(channels==2) col=kGreen+3;
  model.plotOn(xframe,LineColor(col));

  RooHist* hpull = xframe->pullHist();
  RooPlot* frame3 = x.frame(Title("Pull Distribution")) ;
  frame3->addPlotable(hpull,"P");




  // cosmetics
  TLegend *legend = new TLegend(0.70,0.15,0.85,0.25,NULL,"brNDC");
  legend->SetBorderSize(     0);
  legend->SetFillColor (     0);
  legend->SetTextAlign (    12);
  legend->SetTextFont  (    42);
  legend->SetTextSize  (0.03);

  TH1F *dummyPoints = new TH1F("dummyP","dummyP",1,0,1);
  TH1F *dummyLine = new TH1F("dummyL","dummyL",1,0,1);
  dummyPoints->SetMarkerStyle(kOpenCircle);
  dummyPoints->SetMarkerSize(1.1);
  dummyLine->SetLineColor(col);
  
  legend->AddEntry(dummyPoints, "Simulation", "pe");
  legend->AddEntry(dummyLine, "Parametric Model", "l");
  
  TPaveText *text = new TPaveText(0.15,0.90,0.77,0.98,"brNDC");
  text->AddText("CMS Simulation");
  text->SetBorderSize(0);
  text->SetFillStyle(0);
  text->SetTextAlign(12);
  text->SetTextFont(42);
  text->SetTextSize(0.03);

  TPaveText *titlet = new TPaveText(0.15,0.80,0.60,0.85,"brNDC");
  titlet->AddText(frameTitle.str().c_str());
  titlet->SetBorderSize(0);
  titlet->SetFillStyle(0);
  titlet->SetTextAlign(12);
  titlet->SetTextFont(132);
  titlet->SetTextSize(0.045);
  
  xframe->GetYaxis()->SetTitleOffset(1.5);

  cout << "EFF RMS = " << effSigma(hmass) << "    RMS = " << hmass->GetRMS() << endl;

  pad1->cd();
  stringstream nameFile, nameFileC, nameFilePng;
  nameFile << "fitM" << ssample << "_" << schannel<< "_category"<< scategory << ".pdf";
  nameFileC << "fitM" << ssample << "_" << schannel << "_category"<< scategory << ".C";
  nameFilePng << "fitM" << ssample << "_" << schannel << "_category"<< scategory << ".png";
  
  xframe->Draw(); 
  gPad->Update(); 
  // legend->Draw(); 
  text->Draw(); titlet->Draw();

  pad2->cd() ;
  frame3->Draw() ;
  frame3->SetMinimum(-3);
  frame3->SetMaximum(3);

  TLine *line1 = new TLine(105,0,165,0);
  line1->SetLineColor(kRed);
  line1->Draw();
  
  c1->Print (nameFile.str().c_str());
  c1->SaveAs(nameFileC.str().c_str());
  c1->SaveAs(nameFilePng.str().c_str());

  fileplot->cd();
  xframe->Write();
  // sigmat->Write();
  hmass->Write();
  fileplot->Close();

  if(fitValues!=0){
    fitValues[0] = A1.getVal();
    fitValues[1] = A2.getVal();
    fitValues[2] = A3.getVal();
  }  

  if(fitErrors!=0){
    fitErrors[0] = A1.getError();
    fitErrors[1] = A2.getError();
    fitErrors[2] = A3.getError();
  }

  covQual[0] = fitres->covQual();
  
}


float getFitEdge(float mass, float width, bool low) {
  double windowVal = max(width, float(1.));
  double lowside = (mass >= 275) ? 180. : 100.;
  double highside = (mass >= 650) ? 1500. : 800.;
  if (low) return std::max((mass - 20.*windowVal), lowside);
  else return std::min((mass + 15.*windowVal), highside);
}


float getFitEdgeHighMass(float mass, float width, bool low) {
  if (low) return std::max(200.,double(mass-5*width));
  else return std::min(1500.,double(mass+5*width));
}

//*************************************************************************************************
//Computes Eff Sigma
//*************************************************************************************************


Double_t effSigma(TH1 *hist )
{

  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
    return 0.;
  }
  
  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
    return 0.;
  }
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
  if(total < 100.) {
    cout << "effsigma: Too few entries " << total << endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;
  
  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }   
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;
  
  return widmin;
  
}


///-----------------------------------------------------------------------------
Double_t effSigma(RooAbsPdf *pdf, RooRealVar *obs, Int_t nbins)
{
  TH1 *hist = pdf->createHistogram(obs->GetName(), nbins);
  hist->Scale(nbins);

  return effSigma( hist);
}


void plotRegrVsNoRegr(int channel, int massBin) {
  stringstream filenom, filenoregr;
  filenom << "m4lplots/nominal/fitM" << massBin << "_channel" << channel << ".root";
  filenoregr << "m4lplots/noregr/fitM" << massBin << "_channel" << channel << ".root";

  int col;
  if(channel==0) col=kOrange+7;
  if(channel==1) col=kAzure+2;
  if(channel==2) col=kGreen+3;

  TCanvas *c1 = new TCanvas("c1","c1",750,750);

  TFile *tfilenom = TFile::Open(filenom.str().c_str());
  RooPlot *plotnom = (RooPlot*)tfilenom->Get("m4lplot");
  plotnom->SetMarkerStyle(kOpenSquare);
  plotnom->Draw();
  TPaveText *pavenom = (TPaveText*)tfilenom->Get("TPave");
  pavenom->SetTextColor(col);
  pavenom->Draw("same");

  TFile *tfilenoregr = TFile::Open(filenoregr.str().c_str());
  RooPlot *plotnoregr = (RooPlot*)tfilenoregr->Get("m4lplot");
  plotnoregr->Draw("same");
  TPaveText *pavenoregr = (TPaveText*)tfilenoregr->Get("TPave");
  pavenoregr->Draw("same");

  // cosmetics
  TLegend *legend = new TLegend(0.20,0.65,0.45,0.90,NULL,"brNDC");
  legend->SetBorderSize(     0);
  legend->SetFillColor (     0);
  legend->SetTextAlign (    12);
  legend->SetTextFont  (    42);
  legend->SetTextSize  (0.03);

  TH1F *dummyPointsNom = new TH1F("dummyPNom","dummyPNom",1,0,1);
  TH1F *dummyPointsNoRegr = new TH1F("dummyPNoregr","dummyPNoregr",1,0,1);
  TH1F *dummyLine = new TH1F("dummyL","dummyL",1,0,1);
  dummyPointsNoRegr->SetMarkerStyle(kFullCircle);
  dummyPointsNoRegr->SetMarkerSize(1.1);
  dummyPointsNom->SetMarkerStyle(kFullSquare);
  dummyPointsNom->SetMarkerColor(col);
  dummyPointsNom->SetLineColor(col);
  dummyPointsNom->SetMarkerSize(1.1);
  dummyLine->SetLineColor(col);
  
  legend->AddEntry(dummyPointsNoRegr, "Simulation (E_{std}-p comb.)", "pel");
  legend->AddEntry(dummyPointsNom, "Simulation (E_{regr}-p comb.)", "pel");

  legend->Draw();

  TPaveText *text = new TPaveText(0.15,0.90,0.77,0.98,"brNDC");
  text->AddText("CMS Simulation");
  text->SetBorderSize(0);
  text->SetFillStyle(0);
  text->SetTextAlign(12);
  text->SetTextFont(42);
  text->SetTextSize(0.03);

  text->Draw();

  stringstream frameTitle;
  if(channel==0){frameTitle << "4#mu, m_{H} = ";}
  if(channel==1){frameTitle << "4e, m_{H} = ";}
  if(channel==2){frameTitle << "2e2#mu, m_{H} = ";}
  frameTitle << massBin << " GeV";

  TPaveText *titlet = new TPaveText(0.15,0.80,0.60,0.85,"brNDC");
  titlet->AddText(frameTitle.str().c_str());
  titlet->SetBorderSize(0);
  titlet->SetFillStyle(0);
  titlet->SetTextAlign(12);
  titlet->SetTextFont(132);
  titlet->SetTextSize(0.045);

  titlet->Draw();

  c1->SaveAs("comp.pdf");

}

