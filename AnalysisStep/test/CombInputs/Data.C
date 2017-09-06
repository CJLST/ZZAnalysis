//
// Produce data trees for card maker.
//
// Usage: 
// root -q -b lib.C Data.C++
//

#include "TFile.h"
#include "TTree.h"
#include "ZZAnalysis/AnalysisStep/interface/Category.h"
#include "ZZAnalysis/AnalysisStep/interface/Discriminants.h"
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>




void all(int selAna =-10,  int channels=-1, int categ =-10){


  string schannel, scategory, sselAna;

  if (selAna == 0) sselAna = "Mor16";
  if (selAna == 1) sselAna = "ICHEP16";
  if (selAna == 2) sselAna = "Mor17";

  if (channels == 0) schannel = "4mu";
  if (channels == 1) schannel = "4e";
  if (channels == 2) schannel = "2e2mu";

  if (categ < 0 ) scategory = "ALL";
  if (categ == 0 && selAna == 0 ) scategory = "Untagged";
  if (categ == 1 && selAna == 0 ) scategory = "VBFtagged";

  if (categ == 0 && selAna == 1 ) scategory = "Untagged";
  if (categ == 1 && selAna == 1 ) scategory = "VBF1JetTagged";
  if (categ == 2 && selAna == 1 ) scategory = "VBF2JetTagged";
  if (categ == 3 && selAna == 1 ) scategory = "VHLeptTagged";
  if (categ == 4 && selAna == 1 ) scategory = "VHHadrTagged";
  if (categ == 5 && selAna == 1 ) scategory = "ttHTagged";

  //#FIXME: categories should be renamed, cf. prepareYields.C
  if (categ == 0 && selAna == 2 ) scategory = "UntaggedMor17";
  if (categ == 1 && selAna == 2 ) scategory = "VBF1JetTaggedMor17";
  if (categ == 2 && selAna == 2 ) scategory = "VBF2JetTaggedMor17";
  if (categ == 3 && selAna == 2 ) scategory = "VHLeptTaggedMor17";
  if (categ == 4 && selAna == 2 ) scategory = "VHHadrTaggedMor17";
  if (categ == 5 && selAna == 2 ) scategory = "ttHTaggedMor17";
  if (categ == 6 && selAna == 2 ) scategory = "VHMETTaggedMor17";



  bool useQGTagging = false;
  bool useVHMETTagged = true;


  Short_t z1flav, z2flav; 
  Double_t mass4l, kd;
  float m4l;
  Short_t ExtraZ;
  Short_t nExtraLeptons;
  Short_t nCleanedJets;  

  float ZZPt, p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, PHJ_VAJHU, p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, PWH_hadronic_VAJHU, PZH_hadronic_VAJHU, PFMET;
 
  float p_GG_SIG_ghg2_1_ghz1_1_JHUGen, p_QQB_BKG_MCFM;
  Short_t nJets;
  Short_t nBTaggedJets;
  std::vector<float> * JETQGLikeliHood = 0;
  std::vector<float> * jetpt = 0;
  std::vector<float> * jeteta = 0;
  std::vector<float> * jetphi = 0;
  std::vector<float> * jetmass = 0;
  float jetQGLL[100];
  float jetPHI[100];
  float jet30pt[10];
  float jet30eta[10];
  float jet30phi[10];
  float jet30mass[10];
  float Fisher;



  string FileName = "root://lxcms03//data3/Higgs/170222/AllData/ZZ4lAnalysis.root";
  cout << "Using " << FileName << endl;
  TFile* ggFile = TFile::Open(FileName.c_str());
  TTree* ggTree = (TTree*) ggFile->Get("ZZTree/candTree");
  int  nentries = ggTree->GetEntries();
    
  //--- ggTree part
  ggTree->SetBranchAddress("ZZMass",&m4l);
  ggTree->SetBranchAddress("Z1Flav",&z1flav);
  ggTree->SetBranchAddress("Z2Flav",&z2flav);
  ggTree->SetBranchAddress("nExtraLep",&nExtraLeptons);
  ggTree->SetBranchAddress("nCleanedJets",&nJets);
  ggTree->SetBranchAddress("nCleanedJetsPt30BTagged",&nBTaggedJets);
  ggTree->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal",&p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal);
  ggTree->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal",&p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal);

  ggTree->SetBranchAddress("DiJetFisher",&Fisher);
  ggTree->SetBranchAddress("nExtraZ",&ExtraZ);
  ggTree->SetBranchAddress("nCleanedJetsPt30",&nCleanedJets);
  ggTree->SetBranchAddress("JetQGLikelihood",&JETQGLikeliHood);
  ggTree->SetBranchAddress("p_JQCD_SIG_ghg2_1_JHUGen_JECNominal",&PHJ_VAJHU);
  ggTree->SetBranchAddress("p_JVBF_SIG_ghv1_1_JHUGen_JECNominal", &p_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
  ggTree->SetBranchAddress("pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal", &pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
  ggTree->SetBranchAddress("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal", &PWH_hadronic_VAJHU);
  ggTree->SetBranchAddress("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal",&PZH_hadronic_VAJHU);
  ggTree->SetBranchAddress("p_GG_SIG_ghg2_1_ghz1_1_JHUGen",&p_GG_SIG_ghg2_1_ghz1_1_JHUGen);
  ggTree->SetBranchAddress("p_QQB_BKG_MCFM",&p_QQB_BKG_MCFM);

  ggTree->SetBranchAddress("PFMET",&PFMET);
  ggTree->SetBranchAddress("JetPt",&jetpt);
  ggTree->SetBranchAddress("JetEta",&jeteta);
  ggTree->SetBranchAddress("JetPhi",&jetphi);
  ggTree->SetBranchAddress("JetMass",&jetmass);
  ggTree->SetBranchAddress("ZZPt",&ZZPt);

   
  stringstream output;
  output << "data_obs_"<< scategory<<"_"<< schannel <<".root";
  cout<< output.str()<<endl;
  TFile *newFile = new TFile(output.str().c_str(),"RECREATE");
  newFile->cd();
  TTree* newTree = new TTree("data_obs","data_obs");
  newTree->Branch("mass4l",&mass4l,"mass4l/D");
  newTree->Branch("kd",&kd,"kd/D");

  for(int k=0; k<nentries; k++){
    ggTree->GetEvent(k);
         
    int nj = 0;
    for (unsigned int nj = 0; nj < JETQGLikeliHood->size(); nj++) {
      jetQGLL[nj] = (*JETQGLikeliHood)[nj];
    }
    for(unsigned int kjet =0 ; kjet < jetphi->size(); kjet++){
      jetPHI[kjet] = (*jetphi)[kjet];
    } 
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
    int Cat = -10 ;
    if (selAna == 0) Cat = categoryMor16(nJets, p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal );	    
    if (selAna == 1) Cat = categoryIchep16(nExtraLeptons, ExtraZ, nCleanedJets, nBTaggedJets, jetQGLL, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, PHJ_VAJHU, p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, PWH_hadronic_VAJHU, PZH_hadronic_VAJHU, jetPHI, m4l, useQGTagging);
    if (selAna == 2) Cat = categoryMor17(nExtraLeptons, ExtraZ, nCleanedJets, nBTaggedJets, jetQGLL, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, PHJ_VAJHU, p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, PWH_hadronic_VAJHU, PZH_hadronic_VAJHU, jetPHI, m4l, PFMET, useVHMETTagged, useQGTagging);
 
    
    if (categ >= 0 && categ != Cat ) continue;
    
    if(channels==0 && z1flav*z2flav != 28561) continue;
    if(channels==1 && z1flav*z2flav != 14641) continue;
    if(channels==2 && z1flav*z2flav != 20449) continue;
    
    mass4l=m4l;
    kd = D_bkg_kin(p_GG_SIG_ghg2_1_ghz1_1_JHUGen,p_QQB_BKG_MCFM,z1flav*z2flav,m4l);
    newTree->Fill(); 
    
  }
  newFile->cd();
  newTree->Write("data_obs");
  newFile->Close();
  
}

void Data() {
  for (int FS=0; FS<3; ++FS){
    for (int Cat=0; Cat<7; ++Cat) {
      all(2,FS,Cat);
    }
  }
}
