#include <ZZAnalysis/AnalysisStep/interface/LeptonSFHelper.h>

#include <FWCore/MessageLogger/interface/MessageLogger.h>

using namespace std;

LeptonSFHelper::LeptonSFHelper(int year, std::string const &data_tag) :
  theYear(year),
  h_Ele_ID(nullptr),
  h_Ele_ID_Cracks(nullptr),
  h_Ele_Reco_lowPt(nullptr),
  h_Ele_Reco_midPt(nullptr),
  h_Ele_Reco_highPt(nullptr),
  h_Mu_SF(nullptr),
  h_Mu_Unc(nullptr) {

  TString basePath = Form("$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/");

  // -----ELECTRONS
  TString f_eleID, f_eleID_Cracks, f_eleReco_lowPt, f_eleReco_midPt, f_eleReco_highPt; // filenames
  
  if (year == 2016) {
    // 2016 preVFP Electrons
    if(data_tag.find("ULAPV") != std::string::npos) {
      f_eleID          = basePath+"ElectronSF_UL2016preVFP_nogap.root";
      f_eleID_Cracks   = basePath+"ElectronSF_UL2016preVFP_gap.root";

      f_eleReco_highPt = basePath+"egammaEffi_ptAbove20.txt_EGM2D_UL2016preVFP.root";
      f_eleReco_lowPt  = basePath+"egammaEffi_ptBelow20.txt_EGM2D_UL2016preVFP.root";
      f_eleReco_midPt  = "";      

    } else { // 2016 postVFP Electrons
      f_eleID          = basePath+"ElectronSF_UL2016postVFP_nogap.root";
      f_eleID_Cracks   = basePath+"ElectronSF_UL2016postVFP_gap.root";

      f_eleReco_highPt = basePath+"egammaEffi_ptAbove20.txt_EGM2D_UL2016postVFP.root";
      f_eleReco_lowPt  = basePath+"egammaEffi_ptBelow20.txt_EGM2D_UL2016postVFP.root"; 
      f_eleReco_midPt  = "";      
    }
  } else if (year == 2017) { // 2017 Electrons
    f_eleID            = basePath+"ElectronSF_UL2017_nogap.root";
    f_eleID_Cracks     = basePath+"ElectronSF_UL2017_gap.root";

    f_eleReco_highPt   = basePath+"egammaEffi_ptAbove20.txt_EGM2D_UL2017.root";
    f_eleReco_lowPt    = basePath+"egammaEffi_ptBelow20.txt_EGM2D_UL2017.root";
    f_eleReco_midPt    = "";      
    
  } else if (year == 2018) { // 2018 Electrons
    f_eleID            = basePath+"ElectronSF_UL2018_nogap.root";
    f_eleID_Cracks     = basePath+"ElectronSF_UL2018_gap.root";

    f_eleReco_highPt   = basePath+"egammaEffi_ptAbove20.txt_EGM2D_UL2018.root";
    f_eleReco_lowPt    = basePath+"egammaEffi_ptBelow20.txt_EGM2D_UL2018.root";
    f_eleReco_midPt    = "";      

  } else if (year == 2022) {
   if (data_tag.find("pre_EE") != std::string::npos) { // 2022 preEE
    //The most recent ID scale factors used in this context, calculated by Andro using the 2018UL MVA training, can be found at: /eos/user/a/anpetkov/SF2022eleID_EGMapproved/
     f_eleID           = basePath+"SF2022eleID_preEE.root";

     f_eleReco_highPt  = basePath+"egammaEffi_ptAbove75.txt_EGM2D_2022preEE.root";
     f_eleReco_midPt   = basePath+"egammaEffi_ptBelow75.txt_EGM2D_2022preEE.root";
     f_eleReco_lowPt   = basePath+"egammaEffi_ptBelow20.txt_EGM2D_2022preEE.root";

   } else { // 2022 postEE
     f_eleID           = basePath+"SF2022eleID_postEE.root";

     f_eleReco_highPt  = basePath+"egammaEffi_ptAbove75.txt_EGM2D_2022postEE.root";
     f_eleReco_midPt   = basePath+"egammaEffi_ptBelow75.txt_EGM2D_2022postEE.root";
     f_eleReco_lowPt   = basePath+"egammaEffi_ptBelow20.txt_EGM2D_2022postEE.root";
   }
   
  } else if (year >= 2023) { // FIXME: add 2024!
    if(data_tag.find("pre_BPix") != std::string::npos) { // 2023 preBPix
      //ID - for now using 2022 postEE
      std::cout<<"WARNING 2023 preBPix Electron ID SFs - for now using 2022postEE"<<std::endl;
      f_eleID          = basePath+"SF2022eleID_postEE.root";

      //RECO - SFs for Electrons in 2023PromptC from EG - https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammSFandSSRun3
      f_eleReco_highPt = basePath+"egammaEffi_ptAbove75.txt_EGM2D_2023preBPix.root";
      f_eleReco_midPt  = basePath+"egammaEffi_ptBelow75.txt_EGM2D_2023preBPix.root";
      f_eleReco_lowPt  = basePath+"egammaEffi_ptBelow20.txt_EGM2D_2023preBPix.root";

    } else { // 2023 postBPix
      //ID - for now using 2022 postEE
      std::cout<<"WARNING 2023 postBPix Electron ID SFs - for now using 2022postEE"<<std::endl;
      f_eleID          = basePath+"SF2022eleID_postEE.root";

      //RECO - SFs for Electrons in 2023PromptD from EG - https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammSFandSSRun3
      f_eleReco_highPt = basePath+"egammaEffi_ptAbove75.txt_EGM2D_2023postBPix.root";
      f_eleReco_midPt  = basePath+"egammaEffi_ptBelow75.txt_EGM2D_2023postBPix.root";
      f_eleReco_lowPt  = basePath+"egammaEffi_ptBelow20.txt_EGM2D_2023postBPix.root";
    }

  } else if (year<2016 or year>2024 ) {
    edm::LogError("LeptonSFHelper::") << "Ele SFs for " << year << " is not supported!";
    abort();
  }
  
  TFile* root_file = TFile::Open(f_eleID.Data(),"READ");
  h_Ele_ID = (TH2F*) root_file->Get("EGamma_SF2D")->Clone("h_Ele_ID");
  h_Ele_ID->SetDirectory(nullptr); // This is required to detach the clone from the file
  root_file->Close();
  
  if (f_eleID_Cracks != "") {
    root_file = TFile::Open(f_eleID_Cracks.Data(),"READ");
    h_Ele_ID_Cracks = (TH2F*) root_file->Get("EGamma_SF2D")->Clone("h_Ele_ID_Cracks");
    h_Ele_ID_Cracks->SetDirectory(nullptr);
    root_file->Close();
  }
  
  root_file = TFile::Open(f_eleReco_highPt.Data(),"READ");
  h_Ele_Reco_highPt = (TH2F*) root_file->Get("EGamma_SF2D")->Clone("h_Ele_Reco_highPt");
  h_Ele_Reco_highPt->SetDirectory(nullptr);
  root_file->Close();

  root_file = TFile::Open(f_eleReco_lowPt.Data(),"READ");
  h_Ele_Reco_lowPt = (TH2F*) root_file->Get("EGamma_SF2D")->Clone("h_Ele_Reco_lowPt");
  h_Ele_Reco_lowPt->SetDirectory(nullptr);
  root_file->Close();

  if (f_eleReco_midPt != "") {
    root_file = TFile::Open(f_eleReco_midPt.Data(),"READ");
    h_Ele_Reco_midPt = (TH2F*) root_file->Get("EGamma_SF2D")->Clone("h_Ele_Reco_midPt");
    h_Ele_Reco_midPt->SetDirectory(nullptr);
    root_file->Close();
  }


  // -----MUONS
  TString f_mu;

  if (year==2016) { // 2016 Muons
    f_mu = basePath+"final_HZZ_SF_2016UL_mupogsysts_newLoose.root";
  } else if (year==2017) { // 2017 Muons
    f_mu = basePath+"final_HZZ_SF_2017UL_mupogsysts_newLoose.root";
  } else if (year==2018) { // 2018 Muons
    f_mu = basePath+"final_HZZ_SF_2018UL_mupogsysts_newLoose.root";
  } else if (year==2022) { // 2022 Muons
    if(data_tag.find("pre_EE") != std::string::npos) { // 2022 Muons preEE
      f_mu = basePath+"final_HZZ_SF_Run3_2022_mupogsysts_newLoose_abseta3_fix_BCD_RMS.root";
    } else { // 2022 Muons postEE
      if (data_tag.find("MUON_ID_BYMVA") != std::string::npos) {
        f_mu = basePath + "mu_HZZ_2022_post_EE_MVA_ID.root"; // Muon MVA WP (2022postEE)
      } else {
        f_mu = basePath+"final_HZZ_SF_Run3_2022_mupogsysts_newLoose_abseta3_fix_EFG_RMS.root";
      }
    }
  } else if (year>=2023) { // 2023 Muons FIXME: add 2024
    // 2023 Muons preBPix/postBPix/ - root files taken from /afs/cern.ch/user/y/yujil/public/SF2023/
    if(data_tag.find("pre_BPix") != std::string::npos) { 
      f_mu = basePath+"final_HZZ_SF_2023C_RMS_mupogsysts.root";
    } else {
      f_mu = basePath+"final_HZZ_SF_2023D_RMS_mupogsysts.root";
    }
  } else if (year<2016 or year>2024 ) {
    edm::LogError("LeptonSFHelper::") << "Mu SFs for " << theYear << " is not supported!";
    abort();
  }

  root_file = TFile::Open(f_mu.Data(),"READ");
  h_Mu_SF  = (TH2D*)root_file->Get("FINAL")->Clone("h_Mu_SF");
  h_Mu_Unc = (TH2D*)root_file->Get("ERROR")->Clone("h_Mu_Unc");
  h_Mu_SF->SetDirectory(nullptr);
  h_Mu_Unc->SetDirectory(nullptr);
  root_file->Close();

  cout << "[LeptonSFHelper] SF maps opened from root files for " << year << " " << data_tag << endl;
}

LeptonSFHelper::~LeptonSFHelper() {}

pair<float, float> LeptonSFHelper::getSF(int flav, float pt, float eta, float SCeta, bool isCrack) const
{
   float RecoSF = 1.0;
   float SelSF = 1.0;
   float SF = 1.0;

   float RecoSF_Unc = 0.0;
   float SelSF_Unc = 0.0;
   float SFError = 0.0;
   
   //   cout << " flav = " << flav << " pt = " << pt << " eta = " << eta << " SCeta = " << SCeta << " isCrack = " << isCrack << endl;
   
   // Electron reconstruction SFs
   if(abs(flav) == 11) {
     if(pt < 20.) {
       RecoSF     = h_Ele_Reco_lowPt->GetBinContent(h_Ele_Reco_lowPt->GetXaxis()->FindBin(SCeta),h_Ele_Reco_lowPt->GetYaxis()->FindBin(15.));// FIXME: the histogram contains 1 pt bin only
       RecoSF_Unc = h_Ele_Reco_lowPt->GetBinError  (h_Ele_Reco_lowPt->GetXaxis()->FindBin(SCeta),h_Ele_Reco_lowPt->GetYaxis()->FindBin(15.));
     } else if(pt < 75. && h_Ele_Reco_midPt!= nullptr) {
       RecoSF     = h_Ele_Reco_midPt->GetBinContent(h_Ele_Reco_midPt->GetXaxis()->FindBin(SCeta),h_Ele_Reco_midPt->GetYaxis()->FindBin(std::min(pt,75.f)));
       RecoSF_Unc = h_Ele_Reco_midPt->GetBinError  (h_Ele_Reco_midPt->GetXaxis()->FindBin(SCeta),h_Ele_Reco_midPt->GetYaxis()->FindBin(std::min(pt,75.f)));
     } else {
       RecoSF     = h_Ele_Reco_highPt->GetBinContent(h_Ele_Reco_highPt->GetXaxis()->FindBin(SCeta),h_Ele_Reco_highPt->GetYaxis()->FindBin(std::min(pt,499.f)));
       RecoSF_Unc = h_Ele_Reco_highPt->GetBinError  (h_Ele_Reco_highPt->GetXaxis()->FindBin(SCeta),h_Ele_Reco_highPt->GetYaxis()->FindBin(std::min(pt,499.f)));
     }
     
     // Electron HZZ selection SF
     if (isCrack && h_Ele_ID_Cracks!=nullptr) {
       SelSF     = h_Ele_ID_Cracks->GetBinContent(h_Ele_ID_Cracks->FindFixBin(SCeta, std::min(pt,499.f)));
       SelSF_Unc = h_Ele_ID_Cracks->GetBinError  (h_Ele_ID_Cracks->FindFixBin(SCeta, std::min(pt,199.f)));
     } else {
       SelSF = h_Ele_ID->GetBinContent(h_Ele_ID->FindFixBin(SCeta, std::min(pt,499.f)));
       SelSF_Unc = h_Ele_ID->GetBinError  (h_Ele_ID->FindFixBin(SCeta, std::min(pt,499.f)));
     }

     SF = RecoSF*SelSF;   
     SFError = sqrt( RecoSF_Unc*RecoSF_Unc/(RecoSF*RecoSF) + SelSF_Unc*SelSF_Unc/(SelSF*SelSF) ); // assume full correlation between different electrons (and uncorrelated reco and sel uncertainties)
   }

   //Muon SF
   if(abs(flav) == 13 ) {
     //last bin contains the overflow
     SelSF = h_Mu_SF->GetBinContent(h_Mu_SF->GetXaxis()->FindBin(eta),h_Mu_SF->GetYaxis()->FindBin(std::min(pt,199.f)));
     SelSF_Unc = h_Mu_Unc->GetBinContent(h_Mu_Unc->GetXaxis()->FindBin(eta),h_Mu_Unc->GetYaxis()->FindBin(std::min(pt,199.f)));
     
     SF = SelSF;
     SFError = SelSF_Unc/SelSF; // assume full correlation between different muons (and uncorrelated reco and sel uncertainties)
   }

   return std::make_pair(SF, SFError);
}
