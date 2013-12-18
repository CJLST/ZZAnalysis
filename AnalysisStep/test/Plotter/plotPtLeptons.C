bool debug=0;
bool drawgen=1;

void plotPtLeptons(int mass,int ich=2,bool fullsorted=0,TCanvas *c1=0x0){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TString filename ="root://lxcms02//data/Higgs/rootuplesIn/130613/PRODFSR_8TeV/";  
  
  string schannel;
  if (ich == 0) schannel = "4mu";
  if (ich == 1) schannel = "4e";
  if (ich == 2) schannel = "2e2mu";
  
  //filename.Append(schannel=="2e2mu"?"2mu2e":schannel);
  
  //filename.Form("/data3/2013/HZZ_root/130523_testJhuV3_testnewVH/PRODFSR_8TeV/ZZ4lAnalysis_jhuGenV3H126_allEvts.root");//10
  //filename.Form("/data3/2013/HZZ_root/130523_testJhuV3_testnewVH/PRODFSR_8TeV/ZZ4lAnalysis_H126_powheg15_allEvts.root");//11
  //filename.Append(Form("/HZZ4lTree_H%d.root",mass));
  filename.Append(Form("/ZZ4lAnalysis_H%d.root",mass));
  TString treename;treename.Form("ZZ%sTree/candTree",schannel.c_str());
  
  TFile *fMC = TFile::Open(filename.Data());
  TTree *tree = (TTree*)fMC->Get(treename.Data());
  Float_t genleppt[4],Genleppt[4],leppt[4];
  Int_t iBC;
  std::vector<float> *Z1Mass,*Z2Mass,*Z1Pt,*Z2Pt;
  std::vector<int> *ZZsel;
  std::vector<float>  *Leppt[4];//1=0, *Leppt2=0, *Leppt3=0, *Leppt4=0;
  TBranch *b_genleppt[4],*b_leppt[4],*b_Z1mass,*b_Z2mass,*b_Z1Pt,*b_Z2Pt,*b_iBC,*b_ZZsel;

  int col[4]={kBlue-5,kRed-5,kGreen-5,kGray+1};
  int colline[4]={kBlue,kRed+1,kGreen+3,kBlack};
  TH1F *hReco[4];
  TH1F *hGen[4];
  for(int j=0;j<4;j++){
    Leppt[j]=0;
    TString histname;histname.Form("P_{t}^{%d}_RECO",j+1);
    hReco[j]=new TH1F(histname.Data(),histname.Data(),100,4,79);
    hReco[j]->SetLineColor(col[j]);
    hReco[j]->SetLineStyle(2);
    hReco[j]->SetFillColor(col[j]);
    hReco[j]->SetFillStyle(3002);
    hReco[j]->GetXaxis()->SetTitle("p^{l}_{t} [GeV/c]");
    hReco[j]->GetYaxis()->SetTitle("a.u.");
    hReco[j]->SetLineWidth(2);

    histname.Form("P_{t}^{%d}_GEN",j+1);
    hGen[j]=new TH1F(histname.Data(),histname.Data(),100,4,104);
    hGen[j]->SetLineColor(colline[j]);
    hGen[j]->SetLineStyle(2);
    hGen[j]->SetFillColor(0);
    hGen[j]->SetFillStyle(0);
    hGen[j]->GetXaxis()->SetTitle("p^{l}_{t} [GeV/c]");
    hGen[j]->GetYaxis()->SetTitle("a.u.");
    hGen[j]->SetLineWidth(2);
  }
  
  for(int i=0;i<4;i++){
    TString varadd;varadd.Form("GenLep%dPt",i+1);
    tree->SetBranchAddress(varadd.Data(),&Genleppt[i],&b_genleppt[i]);
    varadd.Form("Lep%dPt",i+1);
    tree->SetBranchAddress(varadd.Data(),&Leppt[i],&b_leppt[i]);
  }
  //  tree->SetBranchAddress("Leppt1",&Leppt1,&b_leppt[0]);
  //  tree->SetBranchAddress("Leppt2",&Leppt2,&b_leppt[1]);
  //  tree->SetBranchAddress("Leppt3",&Leppt3,&b_leppt[2]);
  //  tree->SetBranchAddress("Leppt4",&Leppt4,&b_leppt[3]);
  tree->SetBranchAddress("iBC",&iBC,&b_iBC);
  tree->SetBranchAddress("ZZsel",&ZZsel,&b_ZZsel);
  tree->SetBranchAddress("Z1Mass",&Z1Mass,&b_Z1mass);
  tree->SetBranchAddress("Z2Mass",&Z2Mass,&b_Z2mass);
  tree->SetBranchAddress("Z1Pt",&Z1Pt,&b_Z1Pt);
  tree->SetBranchAddress("Z2Pt",&Z2Pt,&b_Z2Pt);

  int nentries = tree->GetEntries();
  for(int i=0;i<nentries;i++){
    tree->GetEntry(i);
    for(int j=0;j<4;j++){
      genleppt[j]=Genleppt[j];
    }
    if(genleppt[0]<genleppt[1]){
      if(debug)printf("inverting leptons 01\n");
      Float_t dummy = genleppt[0];genleppt[0]=genleppt[1];genleppt[1]=dummy;
    }
    if(genleppt[2]<genleppt[3]){
      if(debug)printf("inverting leptons 23\n");
      Float_t dummy = genleppt[2];genleppt[2]=genleppt[3];genleppt[3]=dummy;    
    }
    for(int j=0;j<4&&fullsorted;j++){
      bool sorted=false;
      while(!sorted){
	sorted=true;
	for(int j=0;j<3&&sorted;j++){	    
	  if(genleppt[j]<genleppt[j+1]){
	    double dummy2 = genleppt[j];genleppt[j]=genleppt[j+1];genleppt[j+1]=dummy2;    
	    sorted=false;
	  }
	}
      }
    }
    for(int j=0;j<4;j++)hGen[j]->Fill(genleppt[j]);
    
    if(iBC<0)continue;
    //if(ZZsel->at(iBC)<70)continue;
    if(ZZsel->at(iBC)<100)continue;
    if(debug)printf("entry %d/%d\n",i,nentries-1);
    for(int j=0;j<4;j++)leppt[j]=Leppt[j]->at(iBC);

    //Protections for Z2 tests
    //if(TMath::Abs(Z1Mass->at(iBC)-91.1876)>TMath::Abs(Z2Mass->at(iBC)-91.1876)){printf("ERROR 1 %f %f!!!!\n",Z1Mass->at(iBC),Z2Mass->at(iBC));return;}
    //if(TMath::Abs(leppt[0]+leppt[1]-Z1Pt->at(iBC))>0.01){printf("ERROR 2 %f %f %f %f %f %f!!!!\n",leppt[0],leppt[1],leppt[2],leppt[3],Z1Pt->at(iBC),Z2Pt->at(iBC));return;}
    //if(TMath::Abs(leppt[2]+leppt[3]-Z2Pt->at(iBC))>0.01){printf("ERROR 3 %f %f %f %f %f %f!!!!\n",leppt[0],leppt[1],leppt[2],leppt[3],Z1Pt->at(iBC),Z2Pt->at(iBC));return;}
    if(leppt[0]<leppt[1]){
      if(debug)printf("inverting leptons 01\n");
      Float_t dummy = leppt[0];leppt[0]=leppt[1];leppt[1]=dummy;
    }
    if(leppt[2]<leppt[3]){
      if(debug)printf("inverting leptons 23\n");
      Float_t dummy = leppt[2];leppt[2]=leppt[3];leppt[3]=dummy;    
    }
    //if(leppt[0]<10){printf("ERROR 6 %f %f %f %f %f %f!!!!\n",leppt[0],leppt[1],leppt[2],leppt[3],Z1Pt->at(iBC),Z2Pt->at(iBC));return;}
    //if(leppt[0]<20){printf("ERROR 7 %f %f %f %f %f %f!!!!\n",leppt[0],leppt[1],leppt[2],leppt[3],Z1Pt->at(iBC),Z2Pt->at(iBC));return;}
    //if(leppt[1]<10){printf("ERROR 8 %f %f %f %f %f %f!!!!\n",leppt[0],leppt[1],leppt[2],leppt[3],Z1Pt->at(iBC),Z2Pt->at(iBC));return;}
    if(fullsorted){
      bool sorted=false;
      while(!sorted){
	sorted=true;
	for(int j=0;j<3&&sorted;j++){
	  if(leppt[j]<leppt[j+1]){
	    dummy = leppt[j];leppt[j]=leppt[j+1];leppt[j+1]=dummy;    
	    sorted=false;
	  }
	}
      }
    }
    for(int j=0;j<4;j++)hReco[j]->Fill(leppt[j]);

  }

  if(debug)printf("plotting\n");
  TString cname;cname.Form("Leppt_%s_HMass%d",schannel.c_str(),mass);
  if(!c1)c1=new TCanvas(cname.Data(),cname.Data(),700,300);
  c1->SetTitle(cname.Data());
  c1->SetName(cname.Data());
  c1->cd();
  if(drawgen){
    hGen[3]->Draw();
    for(int j=2;j>=0;j--) hGen[j]->Draw("SAME");
    for(int j=3;j>=0;j--) hReco[j]->Draw("SAME");
  }
  else{
    hReco[3]->Draw();
    for(int j=2;j>=0;j--) hReco[j]->Draw("SAME");
  }
  TLegend *leg= new TLegend(0.54,0.4,0.88,0.85);
  leg->SetLineColor(0);
  leg->SetFillColor(0);  
  TString legtitle;legtitle.Form("H#rightarrow ZZ* #rightarrow %s",schannel.c_str());
  TString massstring;massstring.Form("m_{H} = %d GeV",mass);
  leg->AddEntry((TObject*)0, legtitle.Data(), "");
  leg->AddEntry((TObject*)0, massstring.Data(), "");
  TLegendEntry *before = leg->AddEntry("Before Analysis Selection","Before Analysis Selection","f");
  TLegendEntry *after  = leg->AddEntry("After Analysis Selection","After Analysis Selection","f");
  after->SetLineStyle(2);
  after->SetFillColor(kGray+1);
  after->SetLineColor(kGray+1);
  after->SetFillStyle(3002);
  before->SetLineStyle(2);
  before->SetFillColor(0);
  before->SetLineColor(kBlack);

  leg->Draw("SAME");
  if(fullsorted)cname.Append("_sorted");
  cname.Append(".png");
  c1->SaveAs(cname.Data());

  return;
}


void plotPtLeptons(){
  TCanvas *canv = new TCanvas("a","a",700,300);
  const int nPoints8TeV = 49;
  int masses8TeV[nPoints8TeV]   = {115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,135,140,145,150,160,170,180,190,200,220,250,275,300,325,350,375,400,425,450,475,500,525,550,575,600,650,700,750,800,850,900,950,1000};
  for(int im =0;im<nPoints8TeV;im++)
    for(int ich=0;ich<3;ich++)
      //for(int sorted=0;sorted<2;sorted++)
      plotPtLeptons(masses8TeV[im],ich,true,canv);
}
