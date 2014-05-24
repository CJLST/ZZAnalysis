{ 
  bool blind=false;
  bool zerobins=false;

  int rebin =20;
  //  float maxy = 21.5; //was 16.5
  //  float maxy = 10.;
  
  bool doHiMass = false;


  TString xlabel_M4l = "m_{4l} [GeV]";
  TString xlabel_Mll = "m_{ll} [GeV]";
  TString xlabel_MZ1 = "m_{Z1} [GeV]";
  TString xlabel_MZ2 = "m_{Z2} [GeV]";
  TString xlabel_Mllll = "m_{4l+#gamma} [GeV]";
  TString ylabel_M   = "GeV"; TString ylabel_P   = "GeV";
  TString xlabel_SIP = "SIP3D";
  TString xlabel_ISO = "R_{iso}";
  TString xlabel_CombRelIso2 = "R_{iso,i}+R_{iso,j}";


  TString outfile = "sum_"; 
  outfile += (long) LumiNormalization::getCME();
  outfile += "TeV.root";
  TFile* file= TFile::Open(outfile, "recreate");


  THStack*s;
  TH1F*d;

  //---------- mZZ
  THStack*s;
  TH1F*d;

  getSum( "ZZMass", s, d);
    
  string drawopt="";
  if (blind) drawopt="blind";
  //  drawopt +="m4l100";
  c = newCanvas("all_ZZMass70-800_10GeV");
  drawStack(s, d, 20, "", "", xlabel_M4l,ylabel_M,70,800,0,40); //was 22 at 7TeV

  c = newCanvas("all_ZZMass70-600_10GeV");
  drawStack(s, d, 20, "drawOverflow", "", xlabel_M4l,ylabel_M,70,600,0,40);//was 22 at 7TeV

  file->mkdir("ZZMass")->cd();  
  s->Write();
  d->Write();



  //---------- Zoom
  getSum( "ZZMass_zoom", s, d);
  
  THStack* sskip = skipBins(s, 1.5);
  TH1F* dskip = skipBins(d, 1.5);

  c = newCanvas("all_ZZMass_70-180_3GeV");
  string drawopt="";
  if (zerobins) drawopt="draw0Bins";
  drawStack(sskip, dskip, 6, drawopt, "", xlabel_M4l,ylabel_M,70.5,180,0,25); //was 10 at 7TeV
//   TLegend* legend = getLegend();
//   legend->SetY1NDC(0.7);
//   gPad->Modified();
//   legend->Modify();

//   c = newCanvas("all_ZZMass_100-170");
//   drawStack(sskip, dskip, 4, drawopt, "", xlabel_M4l,ylabel_M,100,170,0,7);

  file->mkdir("ZZMass_zoom")->cd();  
  s->Write();
  d->Write();

  //---------- Z1, Z2
  c = newCanvas("all_Z1Mass");
  getSum( "Z1", s, d);
  drawStack(s,d,5,"","",xlabel_MZ1,ylabel_M,40.,119.,0.,120.);//was 60 at 7TeV

  file->mkdir("Z1")->cd();  
  s->Write();
  d->Write();



  c = newCanvas("all_Z2Mass");
  getSum( "Z2", s, d);
  drawStack(s,d,5,"","",xlabel_MZ2,ylabel_M,12.5,119.,0.,45.); //was 21 at 7TeV

  file->mkdir("Z2")->cd();  
  s->Write();
  d->Write();

  //---------- High-mass

  if (doHiMass) {
    c = newCanvas("all_ZZMass_zz");
    getSum( "ZZMass_zz", s, d);
    drawStack(s, d, 40, drawopt, "", xlabel_M4l,ylabel_M,100,800,0,25); 

    c = newCanvas("all_Z1Mass_zz");
    getSum( "Z1_zz", s, d);
    drawStack(s,d,5,"","",xlabel_MZ1,ylabel_M,60.,119.,0.,60.); 

    c = newCanvas("all_Z2Mass_zz");
    getSum( "Z2_zz", s, d);
    drawStack(s,d,5,"","",xlabel_MZ2,ylabel_M,60.,119.,0.,20.);
  }
  
  //---------- Write out result  


   c = newCanvas("all_LDlow");
   getSum( "LD_lowmass", s, d);
   drawStack(s,d,2,"","","{D}_{bkg}^{kin}","",0,1,0,8); 
   file->mkdir("LD_lowmass")->cd();  
   s->Write();
   d->Write();
  
   if (false){
     c = newCanvas("all_LDhi");
     getSum( "LD_himass", s, d);
     drawStack(s,d,2,"","","{D}_{bkg}^{kin}","",0,1,0,16); 
     file->mkdir("LD_himass")->cd();  
     s->Write();
     d->Write();
   }

   if(false) {
     c = newCanvas("all_pseudoLD");
     getSum( "pseudoLD", s, d);
     drawStack(s,d,5,"","","pseudo-MELA","",0,1,0,5); 
     file->mkdir("pseudoLD")->cd();  
     s->Write();
     d->Write();
   }

  file->Close();
 
  return;
}
