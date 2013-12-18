void plot() {
  plot("4mu");
  plot("4e");
  plot("2e2mu");  
}



void plot(TString finalState) {
  
  

  // --------------------------------
  // input dir
  // --------------------------------


  //4mu
  //TString finalState = "4mu";
  //4e
  //TString finalState = "4e";
  //2e2mu
  //TString finalState = "2e2mu";
  

  TString phaseSpace = ""; // no CR plots
  //TString phaseSpace = "CREEEEos";
  //TString phaseSpace = "CRZLL";

  bool doHiMass = false;

  bool do2012 = true;

  TString epoch;
  TString XSectionFile;
  if (do2012) {
    TString inputDir="HZZ_plots/130613new/PRODFSR_8TeV";
    TString XSectionFile = "Xsection8TeV_v2.txt";
    //epoch = "25May2012"; //1616
    //epoch = "01Jun2012"; //2424
    //epoch = "08Jun2012"; //2968
    //epoch = "13Jun2012"; //3676
    //epoch = "TopUp24Jun2012";
    //epoch = "Moriond2013"; //12210
    epoch = Legacy2013; //19712

  } else {
    doHiMass = false; // himass is relevant only for 2012
    TString inputDir="HZZ_plots/130613new/PRODFSR";
    TString XSectionFile = "Xsection_v1.txt";
    epoch = "All2011";   //5051
    
  }

  TString lumiFile = "Luminosity.txt";
  

  bool blind=false;
  bool zerobins=false;
  
  bool noWeight = false; // do not scale step plots
  bool addSignal = true; // adding or not the Signal
  bool skipSmallBg = true;
  bool useDDBG = true;   // Replace backgrounds with data-driven estimates
   
  bool mergeVV = false;      // Merge di-bosons as VV
  bool splitZZ = false;      // Split ggZZ from qqZZ
  bool splitZZTau = false;   // Split tau from ZZ
  bool splitDYJets = true; // Split Zjets in b/noB
  
  
  
  // -----------------------------------------------------------
  
  if (! TString(gSystem->GetLibraries()).Contains("Number")) {
      gROOT->LoadMacro("$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/test/Plotter/root_lib/TFileServiceWrapper.cc+");
      gROOT->LoadMacro("$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/interface/Histograms.h");
      gROOT->LoadMacro("$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/test/Plotter/root_lib/Number.cc+");
      gROOT->LoadMacro("$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/test/Plotter/root_lib/XSecReader.cc+");
      gROOT->LoadMacro("$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/test/Plotter/root_lib/LumiNormalization.cc+");
      gROOT->LoadMacro("$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/test/Plotter/root_lib/HistoStack.cc+");
      gROOT->LoadMacro("macros.C");
      gROOT->LoadMacro("stackplot.C");
  }
  
  
  //--------------------------------
  // Set the style
  //--------------------------------

  gStyle->Reset();
  TStyle * style = getStyle("ZZ");
  //style->SetOptTitle(1);
  //style->SetPadGridX(false);
  //style->SetPadGridY(false);
  //style->SetOptStat(0);

  style->SetMarkerSize(0.8);
  //style->SetTitleOffset(0.85,"X");
  //style->SetPadBottomMargin(0.16);
  //style->SetPadRightMargin(0.05);
  style->cd();



  //-----------------------------------------
  // Get the normalization from MC xsection 
  //-----------------------------------------
  LumiNormalization lumiNorm(XSectionFile, lumiFile,
			     epoch, finalState, true,"");


  //-----------------------------------------
  // Samples
  //-----------------------------------------
   
  vector<TString> samples;
  
  if (epoch == "All2011"){
    if (finalState == "4mu") {
      samples.push_back("DoubleMu");
    } else if (finalState == "4e") {
      samples.push_back("DoubleEle");
    } else {
      samples.push_back("DoubleOr");
    } 
  }else if (epoch == "Moriond2013"){
    if (finalState == "4mu") {
      samples.push_back("DoubleMu_1963");
    } else if (finalState == "4e") {
      samples.push_back("DoubleEle_1963");
    } else {
      samples.push_back("DoubleOr_1963");
    } 
  }
  
//   if (epoch == "TopUp24Jun2012"){
//     if (finalState == "4mu") {
//       samples.push_back("DoubleMu_5261");
//     } else if (finalState == "4e") {
//       samples.push_back("DoubleEle_5261");
//     } else {
//       samples.push_back("DoubleOr_5261");
//     } 
//   } else if (epoch == "13Jun2012"){
//     if (finalState == "4mu") {
//       samples.push_back("DoubleMu_3676");
//     } else if (finalState == "4e") {
//       samples.push_back("DoubleEle_3676");
//     } else {
//       samples.push_back("DoubleOr_3676");
//     } 
//   } else if (epoch == "08Jun2012"){
//     if (finalState == "4mu") {
//       samples.push_back("DoubleMu_2968");
//     } else if (finalState == "4e") {
//       samples.push_back("DoubleEle_2968");
//     } else {
//       samples.push_back("DoubleOr_2968");
//     } 
//   } else if (epoch == "01Jun2012"){
//     if (finalState == "4mu") {
//       samples.push_back("DoubleMu_2420");
//     } else if (finalState == "4e") {
//       samples.push_back("DoubleEle_2420");
//     } else {
//       samples.push_back("DoubleOr_2420");
//     } 
//   } else if (epoch == "25May2012"){
//     if (finalState == "4mu") {
//       samples.push_back("DoubleMu_1616");
//     } else if (finalState == "4e") {
//       samples.push_back("DoubleEle_1616");
//     } else {
//       samples.push_back("DoubleOr_1616");
//     } 
//   }  else if (epoch == "All2011"){
//     if (finalState == "4mu") {
//       samples.push_back("DoubleMu_NewJSON");
//     } else if (finalState == "4e") {
//       samples.push_back("DoubleEle_NewJSON");
//     } else {
//       samples.push_back("DoubleOr_NewJSON");
//     } 
//   }  else if (epoch == "Complete2011"){
//     if (finalState == "4mu") {
//       samples.push_back("DoubleMu_4670");
//     } else if (finalState == "4e") {
//       samples.push_back("DoubleEle_4670");
//     } else {
//       samples.push_back("DoubleOr_4670");
//     } 
//   }   
  
  if (!skipSmallBg) {
    if (finalState == "4mu" ) {
      samples.push_back("QCDMuEnriched15to20");
      samples.push_back("QCDMuEnriched20to30");
      samples.push_back("QCDMuEnriched30to50");
      samples.push_back("QCDMuEnriched50to80");
      samples.push_back("QCDMuEnriched80to120");
      samples.push_back("QCDMuEnriched120to150");
      samples.push_back("QCDMuEnriched150");
    } else if (finalState == "4e" || finalState == "2e2mu") {     
      samples.push_back("QCDPt20to30EMEnriched");
      samples.push_back("QCD30to80EMEnriched");
      samples.push_back("QCD80to170EMEnriched");
      samples.push_back("QCD20to30BCtoE");
      samples.push_back("QCD30to80BCtoE");
      samples.push_back("QCD80to170BCtoE");
    }
  }
  
  //samples.push_back("TTTo2L2Nu2B");
  
  //samples.push_back("DYJetsToLLTuneZ2M10B");
  //samples.push_back("DYJetsToLLTuneZ2M50B");
  //samples.push_back("DYJetsToLLTuneZ2M10NoB");
  //samples.push_back("DYJetsToLLTuneZ2M50NoB");   
  
  if (!skipSmallBg) {
    
    samples.push_back("WJetsToLNu"); //    
    
    // Single top
    samples.push_back("Tschannel");       // was: TToBLNu_TuneZ2_s-channel
    samples.push_back("Tbarschannel");    // "
    samples.push_back("Ttchannel");       // was: TToBLNu_TuneZ2_t-channel
    samples.push_back("Tbartchannel");    // "
    samples.push_back("TtWchannelDR");    // was: TToBLNu_TuneZ2_tW-channel
    samples.push_back("TbartWchannelDR"); // "
    
    // WW
    samples.push_back("WWTo2L2Nu");  //
  }
  
  // WZ
  //samples.push_back("WZ"); // 
  
  
  // ZZ
  samples.push_back("ZZTo4e");
  samples.push_back("ZZTo4mu");
  samples.push_back("ZZTo4tau"); // FIXME: have a look
  samples.push_back("ZZTo2e2mu");
  samples.push_back("ZZTo2e2tau");
  samples.push_back("ZZTo2mu2tau");
  //     samples.push_back("ZZTo4L_tau");
  
  samples.push_back("ggZZ4l");
  samples.push_back("ggZZ2l2l");     
  
    if(addSignal) {
        samples.push_back("H126");
        samples.push_back("VBFH126");
        samples.push_back("WH126");
        samples.push_back("ZH126");
        samples.push_back("ttH126");
    
    //samples.push_back("H130");
    //      samples.push_back("H140");
    //      samples.push_back("H150");
    //      samples.push_back("H160");
    //      samples.push_back("H170");
    //      samples.push_back("H180");
    //      samples.push_back("H190");
    //      samples.push_back("H200");
    //      samples.push_back("H250");
    //      samples.push_back("H300");
    //      samples.push_back("H350");
    //      samples.push_back("H400");
    //      samples.push_back("H450");
    //      samples.push_back("H500");
    //      samples.push_back("H550");
    //      samples.push_back("H600");
   }
  

  //-----------------------------------------
  // Define groups of samples
  //-----------------------------------------

  HistoStack *stacker = new HistoStack();
   
  // Note that "data" is special!                              // Roberto 
  stacker->setGroup("data",   true, kBlack, 20, 0);            // 
  stacker->setGroup("QCD",    false, kOrange+6,      26, 4);   // was kRed-4 //kBlue-9
  stacker->setGroup("Z+jets",  false, kGreen-5,    28, 6);     // kAzure+2
  stacker->setGroup("Z",       false, kAzure+2,    28, 6);     // 63
  stacker->setGroup("Zbb/cc",    false, kSpring+2,  8, 2);	// kGreen-2 //  //kSpring+2 //kGreen -8
  stacker->setGroup("Z+gamma",   false, kBlue-9,   30, 2);	// kGreen-2
  stacker->setGroup("W+jets",  false, kBlue+2,     27, 5);	// kBlue+1
  stacker->setGroup("tt", false, kGray,  21, 3);               // kCyan  //kCyan-3
  stacker->setGroup("Single top", false, kCyan+2,  25, 3);     // kOrange +6;
  stacker->setGroup("VV", false, kViolet-4,    23, 1);           // ditto 
  //  stacker->setGroup("ZZ", false, kPink+2,    23, 1);           // kPink+9
  stacker->setGroup("ZZ", false, kAzure-9,    23, 1);           // kPink+9 //KMagenta-8
  //  stacker->setGroup("ZZ", false, kViolet-9,    23, 1);
  stacker->setGroup("ggZZ", false, kViolet,    23, 1);
  stacker->setGroup("ZZtau", false, kViolet+1,    23, 1);
  //   stacker->setGroup("ZZPo", false, kViolet,    23, 1);
  stacker->setGroup("WZ", false, kBlue-9,    22, 1); // Was kPink+8
  stacker->setGroup("WW", false, kPink+3,    21, 1);
  stacker->setGroup("h120",   true,  kOrange+10,    20, 50);
  stacker->setGroup("h125",   true,  kOrange+10,    20, 50);
  stacker->setGroup("h126",   true,  kOrange+10,    20, 50);
  stacker->setGroup("ggh126", true,  kOrange+10,    20, 50);
  stacker->setGroup("VBFh126",true,  kOrange+10,    20, 50);
  stacker->setGroup("h130",   true,  kYellow,    20, 50);
  stacker->setGroup("h140",   true,  kRed-4,    24, 51); // was:20
  stacker->setGroup("h150",   true,  kOrange+10,    20, 51);  
  stacker->setGroup("h160",   true,  kOrange, 21, 51);
  stacker->setGroup("h170",   true,  kOrange, 21, 51);
  stacker->setGroup("h180",   true,  kOrange, 21, 51);
  stacker->setGroup("h190",   true,  kOrange+7, 21, 51);
  stacker->setGroup("h200",   true,  kOrange+7, 21, 51);    
  stacker->setGroup("h250",   true,  kOrange-3,  22, 52);   
  stacker->setGroup("h300",   true,  kOrange-3,  22, 52); 
  stacker->setGroup("h350",   true,  kRed+1,  22, 52); //kOrange-3
  stacker->setGroup("h400",   true,  kOrange-3, 21, 51);
  stacker->setGroup("h450",   true,  kOrange-3,  22, 52);
  stacker->setGroup("h500",   true,  kOrange-3,  22, 52);
  stacker->setGroup("h550",   true,  kOrange-3,  22, 52);
  stacker->setGroup("h600",   true,  kOrange-3, 21, 51);
  
  //stacker->setGroup("Z+1jets",  false, kAzure+2,    27, 6);
  //stacker->setGroup("Z+2jets",  false, kAzure+3,    26, 6);
  //stacker->setGroup("Z+3jets",  false, kAzure+4,    25, 6);
  //stacker->setGroup("Z+4jets",  false, kAzure+5,    24, 6);
  //stacker->setGroup("Z+5jets",  false, kAzure+6,    23, 6);
  
  
  //--- There's no need to assign data samples to "data" groups, providing
  //--- that the sample name starts with "DoubleMu", "DoubleEle", "DoubleOr", or "data".
  
  stacker->assignToGroup("H120","h120");
  stacker->assignToGroup("H125","h125");


  stacker->assignToGroup("H126","h126");
  stacker->assignToGroup("ttH126","h126");
  stacker->assignToGroup("VBFH126","h126");
  stacker->assignToGroup("WH126","h126");
  stacker->assignToGroup("ZH126","h126");
/*
  //stacker->assignToGroup("H126","h126");
  stacker->assignToGroup("H126","h126");
  stacker->assignToGroup("ttH126","h126");
  stacker->assignToGroup("VBFH126","VBFh126");
  stacker->assignToGroup("WH126","VBFh126");
  stacker->assignToGroup("ZH126","VBFh126");
*/

  stacker->assignToGroup("H130","h130");
  stacker->assignToGroup("H140","h140");
  stacker->assignToGroup("H150","h150");
  stacker->assignToGroup("H160","h160");
  stacker->assignToGroup("H170","h170");
  stacker->assignToGroup("H180","h180");
  stacker->assignToGroup("H190","h190");
  stacker->assignToGroup("H200","h200");
  stacker->assignToGroup("H250","h250");
  stacker->assignToGroup("H300","h300");
  stacker->assignToGroup("H350","h350");
  stacker->assignToGroup("H400","h400");
  stacker->assignToGroup("H450","h450");
  stacker->assignToGroup("H500","h500");
  stacker->assignToGroup("H550","h550");
  stacker->assignToGroup("H600","h600");
  //    stacker->assignToGroup("VVjets", "VV");
  //    stacker->assignToGroup("ZZTo4L", "ZZ");
  if (mergeVV) {
    stacker->assignToGroup("ZZTo4e", "VV");
    stacker->assignToGroup("ZZTo4mu", "VV");
    stacker->assignToGroup("ZZTo4tau", "VV");
    stacker->assignToGroup("ZZTo2e2mu", "VV");
    stacker->assignToGroup("ZZTo2e2tau", "VV");
    stacker->assignToGroup("ZZTo2mu2tau", "VV");
    stacker->assignToGroup("ZZTo4L_tau", "VV");
    stacker->assignToGroup("ggZZ4l", "VV");
    stacker->assignToGroup("ggZZ2l2l", "VV");     
    stacker->assignToGroup("WZTo3LNu", "VV");
    stacker->assignToGroup("WZTo3LNuv2", "VV");
    stacker->assignToGroup("WWTo2L2Nu", "VV");
  } else {
    stacker->assignToGroup("ZZTo4e", "ZZ");
    stacker->assignToGroup("ZZTo4mu", "ZZ");
    stacker->assignToGroup("ZZTo2e2mu", "ZZ");
    if (splitZZTau) {
      stacker->assignToGroup("ZZTo4L_tau", "ZZtau");
      stacker->assignToGroup("ZZTo4tau", "ZZtau");       
      stacker->assignToGroup("ZZTo2e2tau", "ZZtau");
      stacker->assignToGroup("ZZTo2mu2tau", "ZZtau");
    } else {
      stacker->assignToGroup("ZZTo4L_tau", "ZZ");
      stacker->assignToGroup("ZZTo4tau", "ZZ");       
      stacker->assignToGroup("ZZTo2e2tau", "ZZ");
      stacker->assignToGroup("ZZTo2mu2tau", "ZZ");
    } 
    if (splitZZ) {
      stacker->assignToGroup("ggZZ4l", "ggZZ");
      stacker->assignToGroup("ggZZ2l2l", "ggZZ");
    } else {
      stacker->assignToGroup("ggZZ4l", "ZZ");
      stacker->assignToGroup("ggZZ2l2l", "ZZ");
    }
    stacker->assignToGroup("WZTo3LNu", "WZ");
    stacker->assignToGroup("WZTo3LNuv2", "WZ");
    stacker->assignToGroup("WWTo2L2Nu", "WW");
  }   
  stacker->assignToGroup("TTjets",  "tt");
  stacker->assignToGroup("TTTo2L2Nu2B",  "tt");
  stacker->assignToGroup("TToBLNu_TuneZ2_s-channel",  "Single top");
  stacker->assignToGroup("TToBLNu_TuneZ2_t-channel",  "Single top");
  stacker->assignToGroup("TToBLNu_TuneZ2_tW-channel",  "Single top");
  stacker->assignToGroup("Tschannel",  "Single top");
  stacker->assignToGroup("Tbarschannel",  "Single top");
  stacker->assignToGroup("Ttchannel",  "Single top");
  stacker->assignToGroup("Tbartchannel",  "Single top");
  stacker->assignToGroup("TtWchannelDR",  "Single top");
  stacker->assignToGroup("TbartWchannelDR",  "Single top");
  
  stacker->assignToGroup("DYJetsToLLTuneZ2M10","Z");
  stacker->assignToGroup("DYJetsToLLTuneZ2M50","Z");
   
  if (splitDYJets) {
    stacker->assignToGroup("DYJetsToLLTuneZ2M10B","Zbb/cc");
    stacker->assignToGroup("ZZTo4tau","Z+jets");
    stacker->assignToGroup("DYJetsToLLTuneZ2M50B","Zbb/cc");
    stacker->assignToGroup("ZZTo4tau","Z+jets");
  } else {
    stacker->assignToGroup("DYJetsToLLTuneZ2M10B","Z");
    stacker->assignToGroup("DYJetsToLLTuneZ2M10NoB","Z");
    stacker->assignToGroup("DYJetsToLLTuneZ2M50B","Z");
    stacker->assignToGroup("DYJetsToLLTuneZ2M50NoB","Z");
  }
 
  stacker->assignToGroup("ZGtoLLG",     "Z+gamma");
  stacker->assignToGroup("ZGToEEG",     "Z+gamma");
  stacker->assignToGroup("ZGToMuMuG",   "Z+gamma");
  stacker->assignToGroup("ZGToTauTauG", "Z+gamma");
  
  stacker->assignToGroup("WJetsToLNu",     "W+jets");
  
  stacker->assignToGroup("QCDMuEnriched15to20","QCD");
  stacker->assignToGroup("QCDMuEnriched20to30","QCD");
  stacker->assignToGroup("QCDMuEnriched30to50","QCD");
  stacker->assignToGroup("QCDMuEnriched50to80","QCD");
  stacker->assignToGroup("QCDMuEnriched80to120","QCD");
  stacker->assignToGroup("QCDMuEnriched120to150","QCD");
  stacker->assignToGroup("QCDMuEnriched150","QCD");
  stacker->assignToGroup("QCDPt20to30EMEnriched","QCD");
  stacker->assignToGroup("QCD30to80EMEnriched","QCD");
  stacker->assignToGroup("QCD80to170EMEnriched","QCD");
  stacker->assignToGroup("QCD20to30BCtoE","QCD");
  stacker->assignToGroup("QCD30to80BCtoE","QCD");
  stacker->assignToGroup("QCD80to170BCtoE","QCD");

  
  TString rootdir = TString("Plots") + finalState + "/";    
  TString rootdirCR = TString("Plots") + phaseSpace + "/";   

  for (int i = 0; i<samples.size();++i) {
    
    TString sName = samples[i];
    // Open the file
    float scaleFactor = 1.;
    
     
    TString fileName = inputDir +"/ZZ4lAnalysis_"+ sName +".root";
    //cout << fileName << endl;

    TFile *file = TFile::Open(fileName.Data());
    scaleFactor = lumiNorm.getScaleFactor(sName);


    //----------------------------------------------------------------------
    // Which plots
    //----------------------------------------------------------------------

    
    TH1F*  nEventComplete = (TH1F*) file->Get(rootdir + "nEventComplete");     
    nEventComplete->Sumw2();
    
    float Nevt_Gen = nEventComplete->GetBinContent(0);
    if (Nevt_Gen==0.) {
      Nevt_Gen=1.; // for backwards compatibility
      nEventComplete = (TH1F*) file->Get(rootdir + "nEventCompleteNorm");
    }
    
    // Do not scale data - FIXME: put Nevt_Gen =1 for data
    if (sName.Contains("Double",TString::kIgnoreCase) || sName.Contains("data",TString::kIgnoreCase)) Nevt_Gen=1.;
    
    //      TH1F*  nEvent = (TH1F*) file->Get(rootdir + "nEvent");
    //      nEvent->Sumw2();
    
    //cout << "Normalization factor: " << lumiNorm.getNormalizationFactor() << endl;
    cout << setprecision(8);
    cout << "--- sample name: " << sName << " Ngen: "<< Nevt_Gen << " scale factor: " << scaleFactor << " w: " << scaleFactor/Nevt_Gen << endl;    
    
    HCand* hCandM70  = new HCand(rootdir + "hCandZZCandM70",file);
    HCand* hCandM100 = new HCand(rootdir + "hCandZZCandM100",file);
    HCand* hCandM100_zz = new HCand(rootdir + "hCandZZCandM100_zz",file);

    if (phaseSpace != "") {
      TH1F*  hZ1Mass = (TH1F*) file->Get(rootdir + "hZ1Mass");
      HCand* hCandCR = new HCand(rootdirCR + "hCandCR",file);
    }

    HCand* hCandMela05   = new HCand(rootdir + "hCandZZCandM100_Mela05",file);
    HCand* hCandMRange = new HCand(rootdir + "hCandZZCandM121_131",file);
        
    //cout << rootdir << endl;
    
    //scaling the plots     
    if(!noWeight) {
      //       nEventComplete->Scale(scaleFactor/Nevt_Gen);
    }
    
    
    scaleFactor = scaleFactor/Nevt_Gen;
    
    //----------------------------------------------------------------------
    // Adding the plots to the stacker
    //----------------------------------------------------------------------
    
    stacker->add(sName, "nEventComplete", nEventComplete, scaleFactor);   


    HCand* thehCandM70 = hCandM70;
    HCand* thehCandM100 = hCandM100;

    // Use 
    bool plotForMela05 = false;
    bool plotForMassRange = false;
    int selection = 0;
    if (plotForMela05) {
      selection = 4;
      thehCandM70 = hCandMela05;
      thehCandM100 = hCandMela05;
    } else if (plotForMassRange) {
      selection = 3;
      thehCandM70 = hCandMRange;
      thehCandM100 = hCandMRange;
    }
    

    stacker->add(sName, "ZZMass_M70",  thehCandM70->hZZMass, scaleFactor);
    stacker->add(sName, "Z1Mass_M100", thehCandM100->hZ1Mass, scaleFactor);
    stacker->add(sName, "Z2Mass_M100", thehCandM100->hZ2Mass, scaleFactor);


    stacker->add(sName, "CosThetaStar", thehCandM100->hCosThetaStar, scaleFactor);
    stacker->add(sName, "CosTheta1", thehCandM100->hCosTheta1, scaleFactor);
    stacker->add(sName, "CosTheta2", thehCandM100->hCosTheta2, scaleFactor);
    stacker->add(sName, "Phi", thehCandM100->hPhi, scaleFactor);
    stacker->add(sName, "Phi1", thehCandM100->hPhi1, scaleFactor);
    //    stacker->add(sName, "LD", thehCandM100->hLD, scaleFactor);
    stacker->add(sName, "LD_lowmass", thehCandM100->hLD_lowmass, scaleFactor);
    stacker->add(sName, "LD_himass", thehCandM100->hLD_himass, scaleFactor);
    stacker->add(sName, "pseudoLD", thehCandM100->hPseudoLD, scaleFactor);



    
    stacker->add(sName, "ZZMass_M100_zz", hCandM100_zz->hZZMass, scaleFactor);
    stacker->add(sName, "Z1Mass_M100_zz", hCandM100_zz->hZ1Mass, scaleFactor);
    stacker->add(sName, "Z2Mass_M100_zz", hCandM100_zz->hZ2Mass, scaleFactor);

    if (phaseSpace != "") {
      stacker->add(sName, "hZ1Mass", hZ1Mass, scaleFactor);   
      stacker->add(sName, "Z1Mass_afterCR", hCandCR->hZ1Mass, scaleFactor);
      stacker->add(sName, "Z2Mass_afterCR", hCandCR->hZ2Mass, scaleFactor);
      stacker->add(sName, "ZaMass_afterCR", hCandCR->hZaMass, scaleFactor);
      stacker->add(sName, "ZbMass_afterCR", hCandCR->hZbMass, scaleFactor);
      stacker->add(sName, "ZZMass_afterCR", hCandCR->hZZMass, scaleFactor);
    }

    

  }

   
   TCanvas* c;
   THStack* s;
   TH2F* h2s;
   TH2F* h2b;
   TH1F* d;
   
   string xlabel_Mll = "m_{ll} [GeV]";
   string xlabel_MZ1 = "m_{Z1} [GeV]";
   string xlabel_MZ2 = "m_{Z2} [GeV]";
   string xlabel_Mllll = "m_{4l} [GeV]";
   string ylabel_M   = "GeV";
   string ylabel_P   = "GeV";
   string xlabel_SIP = "SIP3D";
   string xlabel_ISO = "R_{iso}";
   string xlabel_CombRelIso2 = "R_{iso,i}+R_{iso,j}";
   
   string xlabel_M4l;   
   if (finalState=="4mu") { 
     xlabel_M4l = "m_{4#mu} [GeV]"; //+#gamma
   } else if (finalState=="4e") {
     xlabel_M4l = "m_{4e} [GeV]"; //+#gamma
   } else {
     xlabel_M4l = "m_{2e2#mu} [GeV]"; //l+#gamma
   }
   
   //   TString hdir = "Rint:/" + finalState;
   gDirectory->cd("Rint:/");
   TDirectory* hdir = gDirectory->mkdir(finalState); //->cd();
   //   gDirectory->pwd();


   //----------------------------------------------------------------------
   // Step Plots
   //----------------------------------------------------------------------
    
   if(false){
   // Normalized # of events passing cuts
   c = newCanvas(finalState+"_StepPlot");
   c->SetLogy();
   //s = stacker->getStack("nEventComplete","-h120");
   s = stacker->getStack("nEventComplete","");
   d = stacker->getData("nEventComplete");

   THStack* s1 = selectStepPlotBins(s);
   TH1F*    d1 = selectStepPlotBins(d);

   drawStepPlot(s1,d1,"Cut","Events",2,8,1e-1,2e7);
   stepPlotYield(s, d);
   cout << endl;
   //stepPlotYield(s,d,2,10,3,2);
   stepPlotYield(s,d,8,10,6);

  
   if (false) {
     c = newCanvas(finalState+"_nEvent");
     c->SetLogy();
     s = stacker->getStack("nEvent");
     d = stacker->getData("nEvent");
     drawStepPlot(s,d,"Cut","Events",1,15);
     stepPlotYield(s, d,1,15);
   }
   }
   //return;

 
   if(false){


     //----------------------------------------------------------------------
     // Z1Plot
     //----------------------------------------------------------------------
     TString useGroups0="QCD;tt;Z+jets;Zbb/cc;Z+gamma;W+jets;Single top;";
     c = newCanvas(finalState+"_mZ1");
     s = stacker->getStack("hZ1Mass",useGroups0);
     
     d = stacker->getData("hZ1Mass");
     drawStack(s, d, 2, "", "", xlabel_MZ1,ylabel_M, 40, 120);
     
     
     //----------------------------------------------------------------------
     // HCand Plots CR
     //----------------------------------------------------------------------
     
     //TString useGroups="-h120;-h140;-h150;-h200;-h350";
     TString useGroups="QCD;tt;Z+jets;Zbb/cc;Z+gamma;W+jets;Single top;WW;WZ;ZZ;ggZZ";
     
     c = newCanvas(finalState+"_Z1Mass_"+phaseSpace);
     s = stacker->getStack("Z1Mass_afterCR",useGroups);
     d = stacker->getData("Z1Mass_afterCR");
     drawStack(s,d,5,"","",xlabel_MZ1,ylabel_M,40.,119.,0.,220.); //50.//180.
     
     
     c = newCanvas(finalState+"_Z2Mass_"+phaseSpace);
     s = stacker->getStack("Z2Mass_afterCR",useGroups);
     d = stacker->getData("Z2Mass_afterCR");
     drawStack(s,d,5,"","",xlabel_Mll,ylabel_M,0.,120.,0.,80.);  //25. //60.
     
     c = newCanvas(finalState+"_ZaMass_"+phaseSpace);
     s = stacker->getStack("ZaMass_afterCR",useGroups);
     d = stacker->getData("ZaMass_afterCR");
     drawStack(s,d,5,"","",xlabel_Mll,ylabel_M,0.,200.,0.,62.); 
     
     
     c = newCanvas(finalState+"_ZbMass_"+phaseSpace);
     s = stacker->getStack("ZbMass_afterCR",useGroups);
     d = stacker->getData("ZbMass_afterCR");
     drawStack(s,d,5,"","",xlabel_Mll,ylabel_M,0.,200.,0.,62.);
     
     
     c = newCanvas(finalState+"_ZZMass_"+phaseSpace);
     s = stacker->getStack("ZZMass_afterCR",useGroups);
     d = stacker->getData("ZZMass_afterCR");
     drawStack(s,d,20,"","",xlabel_Mllll,ylabel_M,100.,600.,0.,95.); //30.
     
   }
   
   //----------------------------------------------------------------------
   // Z1,Z2 Mass Plots
   if (true) {
     
     //--- 100-600 m4l range
     TString useGroups="Z+jets;ZZ;ggZZ;ZZtau;h126;h125";
     if (!useDDBG) useGroups = "";

     c = newCanvas(finalState+"_Z1Mass");
     THStack* s0 = stacker->getStack("Z1Mass_M100", useGroups);
     if (useDDBG) s = replaceDD(s0, lumiNorm.getLuminosity(), finalState, epoch, selection, "Z1");
     else s = s0;
     d = stacker->getData("Z1Mass_M100");
     storePlot(hdir,"Z1",s,d);
     drawStack(s,d,5,"","",xlabel_MZ1,ylabel_M,40.,119.,0.,150.); // 60 for 7TeV

     c = newCanvas(finalState+"_Z2Mass");
     THStack* s0 = stacker->getStack("Z2Mass_M100", useGroups);
     if (useDDBG) s = replaceDD(s0, lumiNorm.getLuminosity(), finalState, epoch, selection, "Z2");
     else s = s0;
     d = stacker->getData("Z2Mass_M100");
     storePlot(hdir,"Z2",s,d);
     drawStack(s,d,5,"","",xlabel_MZ2,ylabel_M,12.5,119.,0.,25.); //15 for 7TeV
   }
   
   if(doHiMass) {
     //--- 100-600 m4l range - ZZ cross section
     TString useGroups="Z+jets;ZZ;ggZZ;ZZtau";
     if (!useDDBG) useGroups = "";

     c = newCanvas(finalState+"_Z1Mass_zz");
     THStack* s0 = stacker->getStack("Z1Mass_M100_zz", useGroups);
     if (useDDBG) s = replaceDD(s0, lumiNorm.getLuminosity(), finalState, epoch, 2, "Z1");
     else s = s0;
     d = stacker->getData("Z1Mass_M100_zz");
     storePlot(hdir,"Z1_zz",s,d);
     drawStack(s,d,5,"","",xlabel_MZ1,ylabel_M,60.,119.,0.,35.); 

     c = newCanvas(finalState+"_Z2Mass_zz");
     THStack* s0 = stacker->getStack("Z2Mass_M100_zz", useGroups);
     if (useDDBG) s = replaceDD(s0, lumiNorm.getLuminosity(), finalState, epoch, 2, "Z2");
     else s = s0;
     d = stacker->getData("Z2Mass_M100_zz");
     storePlot(hdir,"Z2_zz",s,d);
     drawStack(s,d,5,"","",xlabel_MZ2,ylabel_M,60.,119.,0.,15.); 
     
   }


   // ----------------------------------------------------------------------
   // MELA 
   if (true) {
     TString useGroups="Z+jets;ZZ;ggZZ;ZZtau;h126;h125";
     if (selection==0 || selection==3) { // full range or plotForMassRange
       c = newCanvas(finalState+"_LD_lowmass");
       THStack* s0 = stacker->getStack("LD_lowmass", useGroups);
       if (useDDBG) s = replaceDD(s0, lumiNorm.getLuminosity(), finalState, epoch, selection, "LDlow");
       else s = s0;
       d = stacker->getData("LD_lowmass");
       storePlot(hdir,"LD_lowmass",s,d);
       drawStack(s,d,2,"","","D_{bkg}^{kin}","",0.,1.,0.,10.); //7 for 7TeV
     }
     
     if (selection==0) { // full 
       c = newCanvas(finalState+"_LD_himass");
       THStack* s0 = stacker->getStack("LD_himass","Z+jets;ZZ;ggZZ;ZZtau;h350");
       if (useDDBG) s = replaceDD(s0, lumiNorm.getLuminosity(), finalState, epoch, 0, "LDhi");
       else s = s0;
       d = stacker->getData("LD_himass");
       storePlot(hdir,"LD_himass",s,d);
       drawStack(s,d,2,"","","D_{bkg}^{kin}","",0.,1.,0.,14.); //10 for 7TeV
     }

     if (selection==3) { // plotForMassRange; we call replaceDD with sel=4 to use the normalization for mela>0.5
       c = newCanvas(finalState+"_PseudoLD");
       THStack* s0 = stacker->getStack("pseudoLD", useGroups);
       if (useDDBG) s = replaceDD(s0, lumiNorm.getLuminosity(), finalState, epoch, 4, "pseudoLD");
       else s = s0;
       d = stacker->getData("pseudoLD");
       storePlot(hdir,"pseudoLD",s,d);
       drawStack(s,d,2,"","","pseudo-MELA","",0.,1.,0.,14.); //10 for 7TeV
       
     }
     
   }


   //----------------------------------------------------------------------
   // ZZ Mass Plots
   if (true) {

     TString useGroups="Z+jets;ZZ;ggZZ;ZZtau;h120;h125;h126;h140;h150;h200;h350";
     if (!useDDBG) useGroups = "";


     //--- 70-600 range
     THStack* s0 = stacker->getStack("ZZMass_M70", useGroups);
     if (useDDBG) s = replaceDD(s0, lumiNorm.getLuminosity(), finalState, epoch, selection, "ZZ");
     else s = s0;
     d = stacker->getData("ZZMass_M70");     
     storePlot(hdir,"ZZMass",s,d);

     THStack* ss = s; //skipBins(s);
     TH1F* dd    = d; //skipBins(d);

     string drawopt="";
//     if (blind) drawopt="blind";
//     drawopt +="m4l100";
     
     if (plotForMela05) {
       c = newCanvas(finalState+"ZZMass_100-800");
       drawStack(ss, dd, 20, "", "", xlabel_M4l,ylabel_M,70,800,0,15);     
     } else {
       c = newCanvas(finalState+"ZZMass_70-800");
       drawStack(ss, dd, 20, "", "", xlabel_M4l,ylabel_M,70,800,0,25); //15 for 7TeV

       c = newCanvas(finalState+"ZZMass_70-600");
       drawStack(ss, dd, 20, "drawOverflow", "", xlabel_M4l,ylabel_M,70,600,0,25); //15 for 7TeV

     }

     //--- Zoom plot
     TString useGroups="Z+jets;ZZ;ggZZ;ZZtau;h126";
     THStack* s0 = stacker->getStack("ZZMass_M70", useGroups);
     if (useDDBG) s = replaceDD(s0, lumiNorm.getLuminosity(), finalState, epoch, selection, "ZZ");
     else s = s0;
     d = stacker->getData("ZZMass_M70");
     storePlot(hdir,"ZZMass_zoom",s,d);

     THStack* ss = skipBins(s, 1.5);
     TH1F* dd    = skipBins(d, 1.5);

     string drawopt="";
     if (zerobins) drawopt="draw0Bins";

     if (plotForMela05) {
       c = newCanvas(finalState+"ZZMass_100-180_3GeV");
       drawStack(ss, dd, 4, drawopt, "", xlabel_M4l,ylabel_M,70.5,180,0,8); 
     } else {
       c = newCanvas(finalState+"ZZMass_70-180_3GeV");
       drawStack(ss, dd, 6, drawopt, "", xlabel_M4l,ylabel_M,70.5,180,0,12); //7 for 7TeV
     } 
   }
   

   //----------------------------------------------------------------------
   // ZZ Mass Plots, high-mass range
   if (doHiMass) {
    
     //--- 100-600 range - ZZ
     TString useGroups="Z+jets;ZZ;ggZZ;ZZtau";
     THStack* s0 = stacker->getStack("ZZMass_M100_zz", useGroups);
     if (useDDBG) s = replaceDD(s0, lumiNorm.getLuminosity(), finalState, epoch, 2, "ZZ");
     else s = s0;
     d = stacker->getData("ZZMass_M100_zz");
     storePlot(hdir,"ZZMass_zz",s,d);

     THStack* ss = s; //skipBins(s);
     TH1F* dd    = d; //skipBins(d);

     c = newCanvas(finalState+"ZZMass_zz");
     string drawopt="";
     if (blind) drawopt="blind";
     //     drawopt +="m4l100";
     drawStack(ss, dd, 40, drawopt, "", xlabel_M4l,ylabel_M,100,800,0,18); 
     
     
   }

   
   
   delete stacker;
}


