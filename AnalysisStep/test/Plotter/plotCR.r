void plotCR() {
  
  

  // --------------------------------
  // input dir
  // --------------------------------

  //TString inputDir="/data3/HZZ_Histograms/Fall11MCJan16Data_Prod1/NewAnalysis";
  //TString inputDir="root://cmsphys05//data/b/botta/HZZ_Histograms/Fall11MCJan16Data_NewAnalyzer/Analyzer180512_STD";
  
  //The 6 different scenario
  //TString inputDir="/data1/b/botta/HZZ_root/270512/PRODSTD";
  //TString inputDir="/data1/b/botta/HZZ_root/220512/PRODNo12";
  //TString inputDir="/data1/b/botta/HZZ_root/220512/PRODLOWe";
  //TString inputDir="/data1/b/botta/HZZ_root/220512/PROD4Over4";
  
  //GreenLight 7 TeV [PRL JSON - JSON2011] - V5_2_0 
  //TString inputDir="/data1/b/botta/HZZ_root/270512/PRODSTD"; 
  //Preapproval 7 TeV [JSON2011] - V5_4_0 
  //TString inputDir="/data1/b/botta/HZZ_root/010612/PRODSTD"; //6/6
  //TString inputDir="/data1/b/botta/HZZ_root/030612/PRODSTD"; //4/4
  //Preapproval 8 TeV [JSONMay25] - V5_4_0 
  //TString inputDir="/data1/b/botta/HZZ_root/030612/PRODSTD_8TeV"; //4/4

  //Tree Prod 060612
  //TString inputDir="/data1/b/botta/HZZ_root/060612/PRODSTD_8TeV"; 
  //TString inputDir="/data1/b/botta/HZZ_root/060612/PRODSTD"; 
  //TString inputDir="/data1/b/botta/HZZ_root/060612/PRODFSR_8TeV"; 
  //TString inputDir="/data1/b/botta/HZZ_root/060612/PRODFSR"; 
  //TString inputDir="root://cmsphys05//data1/b/botta/HZZ_root/060612/PRODFSR"; 
  //TString inputDir="root://cmsphys05//data1/b/botta/HZZ_root/110612/PRODFSR"; 
  //TString inputDir="root://cmsphys05//data1/b/botta/HZZ_root/110612/PRODSTD"; 

  //Tree Prod 110612
  //TString inputDir="/data1/b/botta/HZZ_root/110612/PRODSTD_8TeV"; 
  //TString inputDir="/data1/b/botta/HZZ_root/110612/PRODSTD"; 
  //TString inputDir="/data1/b/botta/HZZ_root/110612/PRODFSR_8TeV"; 
  //TString inputDir="/data1/b/botta/HZZ_root/110612/PRODFSR"; 

  //Tree Prod 140612 - unblinding day
  //TString inputDir="/data1/b/botta/HZZ_root/140612/PRODFSR_8TeV"; 
  //TString inputDir="/data1/b/botta/HZZ_root/140612/PRODFSR"; 
  //TString inputDir="/data1/b/botta/HZZ_root/140612/PRODSTD_8TeV"; 
  //TString inputDir="/data1/b/botta/HZZ_root/140612/PRODSTD"; 

  //Tree Prod 200612 - second unblinding day
  //TString inputDir="/data1/b/botta/HZZ_root/200612/PRODFSR_8TeV"; 
  //TString inputDir="/data1/b/botta/HZZ_root/200612/PRODFSR"; 

  //Tree Prod 210612 - third unblinding day --> correct eScale
  //TString inputDir="root://cmsphys05//data1/b/botta/HZZ_root/210612/PRODFSR_8TeV"; 
  //TString inputDir="root://cmsphys05//data1/b/botta/HZZ_root/210612/PRODFSR"; 


  //Tree Prod 240612 - second unblinding day --> correct eScale
  //TString inputDir="root://cmsphys05//data1/b/botta/HZZ_root/240612/PRODFSR_8TeV"; 
  //TString inputDir="root://cmsphys05//data1/b/botta/HZZ_root/240612/PRODFSR"; 
  
  //===> ICHEP 

  //new set with correct LD
  //TString inputDir="root://lxcms00//data3/2012/HZZ_root/310812/PRODFSR_8TeV";
  //TString inputDir="root://lxcms00//data3/2012/HZZ_root/310812/PRODFSR";
  //and with muon pT lowered to 3
  //TString inputDir="root://lxcms00//data3/2012/HZZ_root/310812_mupt3/PRODFSR_8TeV";
  //TString inputDir="root://lxcms00//data3/2012/HZZ_root/310812_mupt3/PRODFSR";

  //===> after ICHEP starts here

  //hcp data and mc: no calibrations
  //TString inputDir="root://lxcms00//data3/2012/HZZ_root/131012_nocalib/PRODFSR_8TeV";
  //TString inputDir="root://lxcms00//data3/2012/HZZ_root/131012_nocalib/PRODFSR";
  //hcp data and mc: muon mom-scale + ele calib - no regression 
  //TString inputDir="root://lxcms00//data3/2012/HZZ_root/131012_noregr/PRODFSR_8TeV";
  //TString inputDir="root://lxcms00//data3/2012/HZZ_root/131012_noregr/PRODFSR";
  //hcp data and mc: muon mom-scale + ele calib + regression
  //TString inputDir="root://lxcms00//data3/2012/HZZ_root/131012/PRODFSR_8TeV";
  //TString inputDir="root://lxcms00//data3/2012/HZZ_root/151012/PRODFSR";

  //===> after first unblinding
 
  TString inputDir="root://lxcms00//data3/2012/HZZ_root/191012/PRODFSR";
  //TString inputDir="root://lxcms00//data3/2012/HZZ_root/191012/PRODFSR_8TeV";

  //4mu
  //TString finalState = "4mu";
  //4e
  //TString finalState = "4e";
  //2e2mu
  TString finalState = "2e2mu";
  
  //
  TString phaseSpace = "CREEMMss";
  //TString phaseSpace = "CRZLL";
  //TString phaseSpace = "CRZLLHiSIP";
  //TString phaseSpace = "CRZLLHiSIPKin";

  TString XSectionFile = "Xsection_v1.txt";
  //TString XSectionFile = "Xsection8TeV_v2.txt";
  TString lumiFile = "Luminosity.txt";
  
  //TString epoch = "Complete2011OldCalc"; //4670
  //TString epoch = "Complete2011"; //4968
  //TString epoch = "25May2012"; //1616
  //TString epoch = "01Jun2012"; //2424
  //TString epoch = "08Jun2012"; //2968
  //TString epoch = "13Jun2012"; //3676
  //TString epoch = "TopUp24Jun2012"; //5261
  //TString epoch = "almostRunC2012"; //6080


  TString epoch = "All2011"; //5051
  //TString epoch= "HCP2012"; //12100
  
  bool blind=false;
  
  bool noWeight = false; // do not scale step plots
  bool addSignal = true; // adding or not the Signal
  bool skipSmallBg = true;
  bool useDDBG = true;   // Replace backgrounds with data-driven estimates
   
  bool mergeVV = false;      // Merge di-bosons as VV
  bool splitZZ = false;      // Split ggZZ from qqZZ
  bool splitZZTau = false;   // Split tau from ZZ
  bool splitDYJets = true;   // Split Zjets in b/noB
  
  
  
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
      samples.push_back("DoubleMu_NewJSON");
    } else if (finalState == "4e") {
      samples.push_back("DoubleEle_NewJSON");
    } else {
      samples.push_back("DoubleOr_NewJSON");
    } 
  }else if (epoch == "HCP2012"){
    if (finalState == "4mu") {
      samples.push_back("DoubleMu_1210");
    } else if (finalState == "4e") {
      samples.push_back("DoubleEle_1210");
    } else {
      samples.push_back("DoubleOr_1210");
    } 
  }
  // else if (epoch == "TopUp24Jun2012"){
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
//   }else if (epoch == "08Jun2012"){
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
  
  samples.push_back("TTTo2L2Nu2B");
  
  samples.push_back("DYJetsToLLTuneZ2M10B");
  samples.push_back("DYJetsToLLTuneZ2M50B");
  samples.push_back("DYJetsToLLTuneZ2M10NoB");
  samples.push_back("DYJetsToLLTuneZ2M50NoB");   
  
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
  samples.push_back("WZ"); // 
  
  
  // ZZ
  
  samples.push_back("ZZTo4e");
  samples.push_back("ZZTo4mu");
  //samples.push_back("ZZTo4tau"); // FIXME: have a look
  samples.push_back("ZZTo2e2mu");
  samples.push_back("ZZTo2e2tau");
  samples.push_back("ZZTo2mu2tau");
  //     samples.push_back("ZZTo4L_tau");
  
  samples.push_back("ggZZ4l");
  samples.push_back("ggZZ2l2l");     
  
  if(addSignal) {
    //samples.push_back("H350");
    //    samples.push_back("H200");
    //    samples.push_back("H150");
    samples.push_back("H126");
    //samples.push_back("H120");
    
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
  stacker->setGroup("h125",   true,  kOrange-1,    20, 50);
  stacker->setGroup("h126",   true,  kOrange+10,    20, 50);
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
    stacker->assignToGroup("DYJetsToLLTuneZ2M10NoB","Z+jets");
    stacker->assignToGroup("DYJetsToLLTuneZ2M50B","Zbb/cc");
    stacker->assignToGroup("DYJetsToLLTuneZ2M50NoB","Z+jets");
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
    //     cout << fileName << endl;
    
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
    

    TH1F*  hZ1Mass = (TH1F*) file->Get(rootdir + "hZ1Mass_w"); 
    
    if(finalState=="4e"){
      TH1F*  hZ1MassEBEB = (TH1F*) file->Get(rootdir + "hZ1MassEBEB_w"); 
      TH1F*  hZ1MassEEEE = (TH1F*) file->Get(rootdir + "hZ1MassEEEE_w"); 
    }

    HCand* hCandFullSel = new HCand(rootdir + "hCandZZCandM100_w",file);
    HCand* hCandCR = new HCand(rootdirCR + "hCandCR_w",file);
    
        
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
    stacker->add(sName, "hZ1Mass", hZ1Mass, scaleFactor);   

    if(finalState=="4e"){
      stacker->add(sName, "hZ1MassEBEB", hZ1MassEBEB, scaleFactor);   
      stacker->add(sName, "hZ1MassEEEE", hZ1MassEEEE, scaleFactor);   
    }

    stacker->add(sName, "Z1Mass_afterSel", hCandFullSel->hZ1Mass, scaleFactor);
    stacker->add(sName, "Z2Mass_afterSel", hCandFullSel->hZ2Mass, scaleFactor);
    stacker->add(sName, "ZZMass_afterSel", hCandFullSel->hZZMass, scaleFactor);

    stacker->add(sName, "Z1Mass_afterCR", hCandCR->hZ1Mass, scaleFactor);
    stacker->add(sName, "Z2Mass_afterCR", hCandCR->hZ2Mass, scaleFactor);
//     stacker->add(sName, "ZaMass_afterCR", hCandCR->hZaMass, scaleFactor);
//     stacker->add(sName, "ZbMass_afterCR", hCandCR->hZbMass, scaleFactor);
    stacker->add(sName, "ZZMass_afterCR", hCandCR->hZZMass, scaleFactor);
    stacker->add(sName, "CombRelIso4_afterCR", hCandCR->hCombRelIso4, scaleFactor);
    stacker->add(sName, "SIP4_afterCR", hCandCR->hSIP4, scaleFactor);

    
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
     xlabel_M4l = "m_{4#mu} [GeV]";
   } else if (finalState=="4e") {
     xlabel_M4l = "m_{4e} [GeV]";
   } else {
     xlabel_M4l = "m_{2e2#mu} [GeV]";
   }
   
   TString hdir = "Rint:/" + finalState;
   gDirectory->cd("Rint:/");
   gDirectory->mkdir(finalState); //->cd();
   //   gDirectory->pwd();


   //----------------------------------------------------------------------
   // Step Plots
   //----------------------------------------------------------------------
    
   if(false){
   // Normalized # of events passing cuts
   c = newCanvas(finalState+"_StepPlot");
   c->SetLogy();
   //s = stacker->getStack("nEventComplete","h126");
   s = stacker->getStack("nEventComplete","");
   d = stacker->getData("nEventComplete");

   THStack* s1 = selectStepPlotBins(s);
   TH1F*    d1 = selectStepPlotBins(d);

   drawStepPlot(s1,d1,"Cut","Events",2,8,1e-1,2e8);
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
   drawStack(s, d, 1, "", "", xlabel_MZ1,ylabel_M, 40, 120);

   c = newCanvas(finalState+"_mZ1_zoom");
   s = stacker->getStack("hZ1Mass",useGroups0);

   d = stacker->getData("hZ1Mass");
   drawStack(s, d, 1, "", "", xlabel_MZ1,ylabel_M, 70, 110);
   

   if(finalState=="4e"){
   
     TString useGroups0="QCD;tt;Z+jets;Zbb/cc;Z+gamma;W+jets;Single top;";
     c = newCanvas(finalState+"_mZ1EBEB");
     s = stacker->getStack("hZ1MassEBEB",useGroups0);
     
     d = stacker->getData("hZ1MassEBEB");
     drawStack(s, d, 1, "", "", xlabel_MZ1,ylabel_M, 40, 120);

     c = newCanvas(finalState+"_mZ1EBEB_zoom");
     s = stacker->getStack("hZ1MassEBEB",useGroups0);
     
     d = stacker->getData("hZ1MassEBEB");
     drawStack(s, d, 1, "", "", xlabel_MZ1,ylabel_M, 70, 110);
     
     
     TString useGroups0="QCD;tt;Z+jets;Zbb/cc;Z+gamma;W+jets;Single top;";
     c = newCanvas(finalState+"_mZ1EEEE");
     s = stacker->getStack("hZ1MassEEEE",useGroups0);
     
     d = stacker->getData("hZ1MassEEEE");
     drawStack(s, d, 1, "", "", xlabel_MZ1,ylabel_M, 40, 120);


     TString useGroups0="QCD;tt;Z+jets;Zbb/cc;Z+gamma;W+jets;Single top;";
     c = newCanvas(finalState+"_mZ1EEEE_zoom");
     s = stacker->getStack("hZ1MassEEEE",useGroups0);
     
     d = stacker->getData("hZ1MassEEEE");
     drawStack(s, d, 1, "", "", xlabel_MZ1,ylabel_M, 70, 110);
     
   }

   }

   //----------------------------------------------------------------------
   // HCand Plots CR
   //----------------------------------------------------------------------

   //TString useGroups="-h120;-h140;-h150;-h200;-h350";
   TString useGroups="QCD;tt;Z+jets;Zbb/cc;Z+gamma;W+jets;Single top;WW;WZ;ZZ;ggZZ";

   c = newCanvas(finalState+"_Z1Mass_"+phaseSpace);
   s = stacker->getStack("Z1Mass_afterCR",useGroups);
   d = stacker->getData("Z1Mass_afterCR");
   drawStack(s,d,5,"","",xlabel_MZ1,ylabel_M,40.,119.,0.,1250.); //50.//180. //250.

   
   c = newCanvas(finalState+"_Z2Mass_"+phaseSpace);
   s = stacker->getStack("Z2Mass_afterCR",useGroups);
   d = stacker->getData("Z2Mass_afterCR");
   drawStack(s,d,5,"","",xlabel_Mll,ylabel_M,0.,120.,0.,550.);  //25. //60. //80. //250.

//    c = newCanvas(finalState+"_ZaMass_"+phaseSpace);
//    s = stacker->getStack("ZaMass_afterCR",useGroups);
//    d = stacker->getData("ZaMass_afterCR");
//    drawStack(s,d,5,"","",xlabel_Mll,ylabel_M,0.,200.,0.,62.); 

   
//    c = newCanvas(finalState+"_ZbMass_"+phaseSpace);
//    s = stacker->getStack("ZbMass_afterCR",useGroups);
//    d = stacker->getData("ZbMass_afterCR");
//    drawStack(s,d,5,"","",xlabel_Mll,ylabel_M,0.,200.,0.,62.);

   
   c = newCanvas(finalState+"_ZZMass_"+phaseSpace);
   s = stacker->getStack("ZZMass_afterCR",useGroups);
   d = stacker->getData("ZZMass_afterCR");
   drawStack(s,d,20,"","",xlabel_Mllll,ylabel_M,100.,600.,0.,450.); //30. //95. //150. //50.


   c = newCanvas(finalState+"_CombRelIso4_"+phaseSpace);
   s = stacker->getStack("CombRelIso4_afterCR",useGroups);
   d = stacker->getData("CombRelIso4_afterCR");
   drawStack(s,d,2,"","",xlabel_ISO,ylabel_M, 0.,5.,0.,180.); //30. //75. //150.

   c = newCanvas(finalState+"_SIP4_"+phaseSpace);
   s = stacker->getStack("SIP4_afterCR",useGroups);
   d = stacker->getData("SIP4_afterCR");
   drawStack(s,d,6,"","",xlabel_SIP,ylabel_M,0.,50.,0.,420.); //30. //170. //220.


   
  
  
   if (false) {

//      TH1F*    d4e_loose = new TH1F("d4e_loose", "d4e_loose", 2000, 0., 1000.);
//      TH1F*    d4mu_loose = new TH1F("d4mu_loose", "d4mu_loose", 2000, 0., 1000.);
//      TH1F*    d2e2mu_loose = new TH1F("d2e2mu_loose", "d2e2mu_loose", 2000, 0., 1000.);
     
//      TH1F*    d4e = new TH1F("d4e", "d4e", 2000, 0., 1000.);
//      TH1F*    d4mu = new TH1F("d4mu", "d4mu", 2000, 0., 1000.);
//      TH1F*    d2e2mu = new TH1F("d2e2mu", "d2e2mu", 2000, 0., 1000.);

//      TH1F*    d4e_zz = new TH1F("d4e_zz", "d4e_zz", 2000, 0., 1000.);
//      TH1F*    d4mu_zz = new TH1F("d4mu_zz", "d4mu_zz", 2000, 0., 1000.);
//      TH1F*    d2e2mu_zz = new TH1F("d2e2mu_zz", "d2e2mu_zz", 2000, 0., 1000.);


//      if (true) { // EPS
//        d4mu_zz  ->Fill(201.2); //A
//        d4mu     ->Fill(167.8); //B (lowmass)
//        d2e2mu_zz->Fill(162.9); //C
//        d4e      ->Fill(139.3); //D (lowmass)
//        d2e2mu_zz->Fill(207.1); //E
//        d4mu     ->Fill(144.9); //F (lowmass)
//        d2e2mu_zz->Fill(243.7); //G
//        d2e2mu_zz->Fill(257.9); //H
//        d4e      ->Fill(216.7); //I (lowmass)
//        d4mu_zz  ->Fill(238.5); //J
//        d2e2mu_zz->Fill(194.6); //K
//        d4mu     ->Fill(222.3); //L (lowmass)
//        d4mu     ->Fill(119.0); //M (lowmass)
//        d4e      ->Fill(125.7); //N (lowmass)
//        d2e2mu_zz->Fill(323.0); //O
//      }
     
//      if (true) { // LP
//        d4e_zz   ->Fill(190.2); //P
//        d4mu_zz  ->Fill(218.9); //Q
//        d4mu_zz  ->Fill(198.8); //R
//        d4mu_zz  ->Fill(308.6); //S
//        d4e_zz   ->Fill(361.8); //T
//        d4mu_zz  ->Fill(457.9); //U
//      }
     
//      if (true) { // Post LP upt to TS

//        d4mu  ->Fill(118.9); //V (lowmass)
//        d2e2mu_zz  ->Fill(315.2); //W
//        d2e2mu_zz  ->Fill(194.3); //X
//        d4mu       ->Fill(231.9); //Y (lowmass)
//        d4e_zz     ->Fill(293.5); //Z
//        d2e2mu_zz  ->Fill(391.1); //AA
//        d4e_zz     ->Fill(233.3); //AB
//        d2e2mu_zz  ->Fill(188.7); //AC       
//      } 

//      if (true) { // Post-TS
//        d2e2mu_zz  ->Fill(307.0); //AD
//        d2e2mu_zz  ->Fill(230.2); //AE
//        d2e2mu_zz  ->Fill(207.0); //AF
//        d2e2mu_zz  ->Fill(210.0); //AG
//        d2e2mu_zz  ->Fill(205.7); //AH
//        d2e2mu_zz  ->Fill(275.7); //AI
//        d4e        ->Fill(157.1); //AJ (lowmass)
//        d4mu_zz    ->Fill(193.9); //AK
//        //30Aug
//        d2e2mu_zz  ->Fill(323.9); //AL
//        d4e_zz     ->Fill(327.8); //AM
//        d4mu_zz    ->Fill(193.4); //AN
//        d2e2mu_zz  ->Fill(255.6); //AO
//        d4mu_zz    ->Fill(240.3); //AP
//        d2e2mu_zz  ->Fill(192.1); //AQ
//        d2e2mu_zz  ->Fill(309.5); //AR
//        d2e2mu     ->Fill(126.1); //AS
//        d4mu_zz    ->Fill(280.4); //AT
//        d4mu_zz    ->Fill(237.9); //AU
//        d4e_zz     ->Fill(371.0); //AV
//      }

//      // ADD Loose EVENTS (MZ1>50, MZ2>12)
//      if (true) {
//        d2e2mu_loose  ->Fill(142.4); //R-A
//        d4mu_loose    ->Fill(211.6); //R-B
//        d4mu_loose    ->Fill(118.8); //R-C
//        d2e2mu_loose  ->Fill(132.5); //R-D
//        d2e2mu_loose  ->Fill(181.5); //R-E
//        d2e2mu_loose  ->Fill(129.7); //R-F
//        d4e_loose     ->Fill(100.5); //R-G
//      }
     
//      d4mu->Add(d4mu_zz);
//      d4e->Add(d4e_zz);
//      d2e2mu->Add(d2e2mu_zz);

//      d4mu_loose->Add(d4mu);
//      d4e_loose->Add(d4e);
//      d2e2mu_loose->Add(d2e2mu);


//      TH1F* d;
//      TH1F* d_zz;
//      TH1F* d_loose;

//      if (finalState=="4mu") {
//        d = d4mu;
//        d_zz = d4mu_zz;
//        d_loose = d4mu_loose;
//      } else if (finalState=="4e") {
//        d = d4e;
//        d_zz = d4e_zz;
//        d_loose = d4e_loose;
//      } else if (finalState=="2e2mu") {
//        d = d2e2mu;
//        d_zz = d2e2mu_zz;
//        d_loose = d2e2mu_loose;
//      }


     // SELECTION RANGE

     TString useGroups="Z+jets;ZZ;ggZZ;ZZtau;h120;h140;h150;h200;h350";
     if (!useDDBG) useGroups = "";

     THStack* s0 = stacker->getStack("ZZMass_afterSel", useGroups);
     if (useDDBG) s = replaceDD(s0, lumiNorm.getLuminosity(), finalState, epoch, 0);
     else s = s0;
     d = stacker->getData("ZZMass_afterSel");
     gDirectory->cd(hdir);
     gDirectory->Add(s);
     if (d) gDirectory->Add(d);
     //c = newCanvas(finalState+"ZZMass");
     //c->SetLogy();
     //drawStack(s, d, 20, "", "", xlabel_M4l,ylabel_M,100,600,2.07e-2,4.7); // was: 2e-3, 2.06e-2
     //     drawStack(s, d, 30, "", "", xlabel_M4l,ylabel_M,100,600,2.2e-2,6.8);

     THStack* ss = skipBins(s);
     TH1F* dd = skipBins(d);

     c = newCanvas(finalState+"ZZMass_lin");
     string drawopt="";
     if (blind) drawopt="blind";
     drawopt +="m4l100";
     drawStack(ss, dd, 20, drawopt, "", xlabel_M4l,ylabel_M,100,600,0,15); 
     //     drawStack(ss, dd, 4, "", "", xlabel_M4l,ylabel_M,100,600,0,9);




     // For zoom plot
     TString useGroups="Z+jets;ZZ;ggZZ;ZZtau;h140;h120";

     THStack* s0 = stacker->getStack("ZZMass_afterSel", useGroups);
     if (useDDBG) s = replaceDD(s0, lumiNorm.getLuminosity(), finalState, epoch, 0);
     else s = s0;
     d = stacker->getData("ZZMass_afterSel");

     gDirectory->cd(hdir);
     gDirectory->mkdir("zoom")->cd();
     gDirectory->Add(s);
     if (d) gDirectory->Add(d);
     


     // THStack* s0 = stacker->getStack("ZZMass_afterMz_zz", useGroups);
//      if (useDDBG) s = replaceDD(s0, lumiNorm.getLuminosity(), finalState, 0);
//      else s = s0;
//      d = stacker->getData("ZZMass_afterMz_zz");
     

//      gDirectory->cd(hdir);
//      gDirectory->Add(s);
//      if (d) gDirectory->Add(d);

//      THStack* ss = skipBins(s);
//      TH1F* dd = skipBins(d);

// //      c = newCanvas(finalState+"Mass_High");
// //      c->SetLogy();
// //      drawStack(ss, dd, 20, "", "", xlabel_M4l,ylabel_M,100,600,2.07e-2,4.7);
// //      //drawStack(s, d, 30, "", "", xlabel_M4l,ylabel_M,100,600,2.2e-2,6.8);

//      c = newCanvas(finalState+"Mass_High_lin");
//      drawStack(ss, dd, 20, "", "", xlabel_M4l,ylabel_M,100,600,0,15);
//      //drawStack(s, d, 30, "", "", xlabel_M4l,ylabel_M,100,600,0,7);


//      //------------ Loose selection
//      THStack* s0 = stacker->getStack("ZZMass_afterMzLoose", useGroups);
//      if (useDDBG) s = replaceDD(s0, lumiNorm.getLuminosity(), finalState, 0);
//      else s = s0;
//      d = stacker->getData("ZZMass_afterMzLoose");

//      gDirectory->cd(hdir);
//      gDirectory->Add(s);
//      if (d) gDirectory->Add(d);
//      //     gDirectory->pwd();
//      //     gDirectory->ls();

//      THStack* ss = skipBins(s);
//      TH1F* dd = skipBins(d);

// //      c = newCanvas(finalState+"Mass_Low");
// //      c->SetLogy();
// //      drawStack(ss, dd, 20, "", "", xlabel_M4l,ylabel_M,100,600,2.07e-2,4.7); // was: 2e-3, 2.06e-2
// //      //     drawStack(s, d, 30, "", "", xlabel_M4l,ylabel_M,100,600,2.2e-2,6.8);

//      c = newCanvas(finalState+"Mass_Loose");
//      drawStack(ss, dd, 20, "", "", xlabel_M4l,ylabel_M,100,600,0,15);
//      //     drawStack(ss, dd, 4, "", "", xlabel_M4l,ylabel_M,100,600,0,9);


//      // For zoom plot
//      TString useGroups="Z+jets;ZZ;ggZZ;ZZtau;h140;h120";

//      THStack* s0 = stacker->getStack("ZZMass_afterMzLoose", useGroups);
//      if (useDDBG) s = replaceDD(s0, lumiNorm.getLuminosity(), finalState, 0);
//      else s = s0;
//      d = stacker->getData("ZZMass_afterMzLoose");


//      gDirectory->cd(hdir);
//      gDirectory->mkdir("zoom")->cd();
//      gDirectory->Add(s);
//      if (d) gDirectory->Add(d);
//      //     gDirectory->pwd();
// //      THStack* ss = skipBins(s);
// //      TH1F* dd = skipBins(d);

// //      c = newCanvas(finalState+"Mass_Loose_test");
// //      drawStack(ss, dd, 20, "", "", xlabel_M4l,ylabel_M,100,600,0,14);
     

   }

  delete stacker;
}


