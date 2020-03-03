#include <ZZAnalysis/AnalysisStep/interface/LeptonSFHelper.h>

using namespace std;
using namespace edm;

LeptonSFHelper::LeptonSFHelper()
{
   // 2016 Electrons
   edm::FileInPath fipEleNotCracks_2016("ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ElectronSF_Legacy_2016_NoGap.root");
   root_file = TFile::Open(fipEleNotCracks_2016.fullPath().data(),"READ");
   h_Ele_notCracks_2016 = (TH2F*) root_file->Get("EGamma_SF2D")->Clone();
      
   edm::FileInPath fipEleCracks_2016("ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ElectronSF_Legacy_2016_Gap.root");
   root_file = TFile::Open(fipEleCracks_2016.fullPath().data(),"READ");
   h_Ele_Cracks_2016 = (TH2F*) root_file->Get("EGamma_SF2D")->Clone();
   
   edm::FileInPath fipEleReco_highPt_2016("ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/Ele_Reco_2016.root");
   root_file = TFile::Open(fipEleReco_highPt_2016.fullPath().data(),"READ");
   h_Ele_Reco_highPT_2016 = (TH2F*) root_file->Get("EGamma_SF2D")->Clone();
      
   edm::FileInPath fipEleReco_lowPt_2016("ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/Ele_Reco_LowEt_2016.root");
   root_file = TFile::Open(fipEleReco_lowPt_2016.fullPath().data(),"READ");
   h_Ele_Reco_lowPT_2016 = (TH2F*) root_file->Get("EGamma_SF2D")->Clone();
   
   // 2017 Electrons
   edm::FileInPath fipEleNotCracks_2017("ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ElectronSF_Legacy_2017_NoGap.root");
   root_file = TFile::Open(fipEleNotCracks_2017.fullPath().data(),"READ");
   h_Ele_notCracks_2017 = (TH2F*) root_file->Get("EGamma_SF2D")->Clone();
   
   edm::FileInPath fipEleCracks_2017("ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ElectronSF_Legacy_2017_Gap.root");
   root_file = TFile::Open(fipEleCracks_2017.fullPath().data(),"READ");
   h_Ele_Cracks_2017 = (TH2F*) root_file->Get("EGamma_SF2D")->Clone();
   
   edm::FileInPath fipEleReco_highPt_2017("ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/Ele_Reco_2017.root");
   root_file = TFile::Open(fipEleReco_highPt_2017.fullPath().data(),"READ");
   h_Ele_Reco_highPT_2017 = (TH2F*) root_file->Get("EGamma_SF2D")->Clone();
   
   edm::FileInPath fipEleReco_lowPt_2017("ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/Ele_Reco_LowEt_2017.root");
   root_file = TFile::Open(fipEleReco_lowPt_2017.fullPath().data(),"READ");
   h_Ele_Reco_lowPT_2017 = (TH2F*) root_file->Get("EGamma_SF2D")->Clone();
      

   // 2018 Electrons
   edm::FileInPath fipEleNotCracks_2018("ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ElectronSF_Legacy_2018_NoGap.root");
   root_file = TFile::Open(fipEleNotCracks_2018.fullPath().data(),"READ");
   h_Ele_notCracks_2018 = (TH2F*) root_file->Get("EGamma_SF2D")->Clone();
   
   edm::FileInPath fipEleCracks_2018("ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ElectronSF_Legacy_2018_Gap.root");
   root_file = TFile::Open(fipEleCracks_2018.fullPath().data(),"READ");
   h_Ele_Cracks_2018 = (TH2F*) root_file->Get("EGamma_SF2D")->Clone();
   
   edm::FileInPath fipEleReco_highPt_2018("ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/Ele_Reco_2018.root");
   root_file = TFile::Open(fipEleReco_highPt_2018.fullPath().data(),"READ");
   h_Ele_Reco_highPT_2018 = (TH2F*) root_file->Get("EGamma_SF2D")->Clone();
   
   edm::FileInPath fipEleReco_lowPt_2018("ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/Ele_Reco_LowEt_2018.root");
   root_file = TFile::Open(fipEleReco_lowPt_2018.fullPath().data(),"READ");
   h_Ele_Reco_lowPT_2018 = (TH2F*) root_file->Get("EGamma_SF2D")->Clone();
   
   // 2016 Muons
   edm::FileInPath fipMu_2016("/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/MuonSF/final_HZZ_SF_2016_legacy_mupogsysts.root"); //"ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ScaleFactors_mu_Moriond2017_v2.root");
   root_file = TFile::Open(fipMu_2016.fullPath().data(),"READ");
   h_Mu_SF_2016  = (TH2D*)root_file->Get("FINAL")->Clone();
   h_Mu_Unc_2016 = (TH2D*)root_file->Get("ERROR")->Clone();
   
   // 2017 Muons
   edm::FileInPath fipMu_2017("/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/MuonSF/final_HZZ_SF_2017_rereco_mupogsysts_3010.root"); //"ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ScaleFactors_mu_Moriond2018_final.root");
   root_file = TFile::Open(fipMu_2017.fullPath().data(),"READ");
   h_Mu_SF_2017  = (TH2D*)root_file->Get("FINAL")->Clone();
   h_Mu_Unc_2017 = (TH2D*)root_file->Get("ERROR")->Clone();
      
   // 2018 Muons
   edm::FileInPath fipMu_2018("/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/MuonSF/final_HZZ_SF_2018_rereco_mupogsysts_3010.root"); //"ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/final_HZZ_muon_SF_2018RunA2D_ER_2702.root");
   root_file = TFile::Open(fipMu_2018.fullPath().data(),"READ");
   h_Mu_SF_2018  = (TH2D*)root_file->Get("FINAL")->Clone();
   h_Mu_Unc_2018 = (TH2D*)root_file->Get("ERROR")->Clone();
   
   cout << "[LeptonSFHelper] SF maps opened from root files." << endl;
}

LeptonSFHelper::~LeptonSFHelper()
{
}

float LeptonSFHelper::getSF(int year, int flav, float pt, float eta, float SCeta, bool isCrack)
{
   float RecoSF = 1.0;
   float SelSF = 1.0;
   float SF = 1.0;
   
   //cout << "year = " << year << " flav = " << flav << " pt = " << pt << " eta = " << eta << " SCeta = " << SCeta << " isCrack = " << isCrack << endl;
   
   // Electron reconstruction SFs
   if(abs(flav) == 11) {
      if(year == 2016)
      {
         if(pt < 20.)
         {
            RecoSF = h_Ele_Reco_lowPT_2016->GetBinContent(h_Ele_Reco_lowPT_2016->GetXaxis()->FindBin(SCeta),h_Ele_Reco_lowPT_2016->GetYaxis()->FindBin(15.));// FIXME: the histogram contains 1 pt bin only
         }
         else
         {
            RecoSF = h_Ele_Reco_highPT_2016->GetBinContent(h_Ele_Reco_highPT_2016->GetXaxis()->FindBin(SCeta),h_Ele_Reco_highPT_2016->GetYaxis()->FindBin(std::min(pt,499.f)));
         }
      }
      else if(year == 2017)
      {
         if(pt < 20.)
         {
            RecoSF = h_Ele_Reco_lowPT_2017->GetBinContent(h_Ele_Reco_lowPT_2017->GetXaxis()->FindBin(SCeta),h_Ele_Reco_lowPT_2017->GetYaxis()->FindBin(15.));// FIXME: the histogram contains 1 pt bin only
         }
         else
         {
            RecoSF = h_Ele_Reco_highPT_2017->GetBinContent(h_Ele_Reco_highPT_2017->GetXaxis()->FindBin(SCeta),h_Ele_Reco_highPT_2017->GetYaxis()->FindBin(std::min(pt,499.f)));
         }
      }
      else if(year == 2018)
      {
         if(pt < 20.)
         {
            RecoSF = h_Ele_Reco_lowPT_2017->GetBinContent(h_Ele_Reco_lowPT_2017->GetXaxis()->FindBin(SCeta),h_Ele_Reco_lowPT_2017->GetYaxis()->FindBin(15.));// FIXME: the histogram contains 1 pt bin only
         }
         else
         {
            RecoSF = h_Ele_Reco_highPT_2017->GetBinContent(h_Ele_Reco_highPT_2017->GetXaxis()->FindBin(SCeta),h_Ele_Reco_highPT_2017->GetYaxis()->FindBin(std::min(pt,499.f)));
         }
      }
      else {
         edm::LogError("LeptonSFHelper::") << "Ele SFs for " << year << " is not supported!";
         abort();
      }
      
      // Electron HZZ selection SF
      if(year == 2016)
      {
         if(isCrack)
            {
               SelSF = h_Ele_Cracks_2016->GetBinContent(h_Ele_Cracks_2016->FindFixBin(SCeta, std::min(pt,199.f)));
            }
            else
            {
               SelSF = h_Ele_notCracks_2016->GetBinContent(h_Ele_notCracks_2016->FindFixBin(SCeta, std::min(pt,199.f)));
            }
         
      }
      
      else if(year == 2017)
      {
         if(isCrack)
            {
               SelSF = h_Ele_Cracks_2017->GetBinContent(h_Ele_Cracks_2017->FindFixBin(SCeta, std::min(pt,499.f)));
            }
            else
            {
               SelSF = h_Ele_notCracks_2017->GetBinContent(h_Ele_notCracks_2017->FindFixBin(SCeta, std::min(pt,499.f)));
            }
      }
      
      else if(year == 2018)
      {
         if(isCrack)
            {
               SelSF = h_Ele_Cracks_2018->GetBinContent(h_Ele_Cracks_2018->FindFixBin(SCeta, std::min(pt,499.f)));
            }
            else
            {
               SelSF = h_Ele_notCracks_2018->GetBinContent(h_Ele_notCracks_2018->FindFixBin(SCeta, std::min(pt,499.f)));
            }
      }
   
      
   else {
      edm::LogError("LeptonSFHelper::") << "Ele SFs for " << year << " is not supported!";
      abort();
      }
      
      SF = RecoSF*SelSF;
   }
   
   //Muon SF
   if(abs(flav) == 13 )
   {
      if(year == 2016)
      {
         SelSF = h_Mu_SF_2016->GetBinContent(h_Mu_SF_2016->GetXaxis()->FindBin(eta),h_Mu_SF_2016->GetYaxis()->FindBin(std::min(pt,199.f))); //last bin contains the overflow
      }
      else if(year == 2017)
      {
         SelSF = h_Mu_SF_2017->GetBinContent(h_Mu_SF_2017->GetXaxis()->FindBin(eta),h_Mu_SF_2017->GetYaxis()->FindBin(std::min(pt,199.f))); //last bin contains the overflow
      }
      else if(year == 2018)
      {
         SelSF = h_Mu_SF_2018->GetBinContent(h_Mu_SF_2018->GetXaxis()->FindBin(eta),h_Mu_SF_2018->GetYaxis()->FindBin(std::min(pt,199.f))); //last bin contains the overflow
      }
      else {
         edm::LogError("LeptonSFHelper::") << "Ele SFs for " << year << " is not supported!";
         abort();
      }
      
      SF = SelSF;
   }
   
    return SF;
}

float LeptonSFHelper::getSFError(int year, int flav, float pt, float eta, float SCeta, bool isCrack)
{
   float RecoSF = 1.0;
   float SelSF = 1.0;
   
   float RecoSF_Unc = 0.0;
   float SelSF_Unc = 0.0;
   float SFError = 0.0;
   
   if(abs(flav) == 11) {
      if(year == 2016)
      {
         if(pt < 20.)
         {
            RecoSF = h_Ele_Reco_lowPT_2016->GetBinContent(h_Ele_Reco_lowPT_2016->GetXaxis()->FindBin(SCeta),h_Ele_Reco_lowPT_2016->GetYaxis()->FindBin(15.));// FIXME: the histogram contains 1 pt bin only
            RecoSF_Unc = h_Ele_Reco_lowPT_2016->GetBinError(h_Ele_Reco_lowPT_2016->GetXaxis()->FindBin(SCeta),h_Ele_Reco_lowPT_2016->GetYaxis()->FindBin(15.));// FIXME: the histogram contains 1 pt bin only
         }
         else
         {
            RecoSF = h_Ele_Reco_highPT_2016->GetBinContent(h_Ele_Reco_highPT_2016->GetXaxis()->FindBin(SCeta),h_Ele_Reco_highPT_2016->GetYaxis()->FindBin(std::min(pt,499.f)));
            RecoSF_Unc = h_Ele_Reco_highPT_2016->GetBinError(h_Ele_Reco_highPT_2016->GetXaxis()->FindBin(SCeta),h_Ele_Reco_highPT_2016->GetYaxis()->FindBin(std::min(pt,499.f)));
         }
      }
      else if(year == 2017)
      {
         if(pt < 20.)
         {
            RecoSF = h_Ele_Reco_lowPT_2017->GetBinContent(h_Ele_Reco_lowPT_2017->GetXaxis()->FindBin(SCeta),h_Ele_Reco_lowPT_2017->GetYaxis()->FindBin(15.));// FIXME: the histogram contains 1 pt bin only
            RecoSF_Unc = h_Ele_Reco_lowPT_2017->GetBinError(h_Ele_Reco_lowPT_2017->GetXaxis()->FindBin(SCeta),h_Ele_Reco_lowPT_2017->GetYaxis()->FindBin(15.));// FIXME: the histogram contains 1 pt bin only
         }
         else
         {
            RecoSF = h_Ele_Reco_highPT_2017->GetBinContent(h_Ele_Reco_highPT_2017->GetXaxis()->FindBin(SCeta),h_Ele_Reco_highPT_2017->GetYaxis()->FindBin(std::min(pt,499.f)));
            RecoSF_Unc = h_Ele_Reco_highPT_2017->GetBinError(h_Ele_Reco_highPT_2017->GetXaxis()->FindBin(SCeta),h_Ele_Reco_highPT_2017->GetYaxis()->FindBin(std::min(pt,499.f)));
         }
      }
      else if(year == 2018)
      {
         if(pt < 20.)
         {
            RecoSF = h_Ele_Reco_lowPT_2017->GetBinContent(h_Ele_Reco_lowPT_2017->GetXaxis()->FindBin(SCeta),h_Ele_Reco_lowPT_2017->GetYaxis()->FindBin(15.));// FIXME: the histogram contains 1 pt bin only
            RecoSF_Unc = h_Ele_Reco_lowPT_2017->GetBinError(h_Ele_Reco_lowPT_2017->GetXaxis()->FindBin(SCeta),h_Ele_Reco_lowPT_2017->GetYaxis()->FindBin(15.));// FIXME: the histogram contains 1 pt bin only
         }
         else
         {
            RecoSF = h_Ele_Reco_highPT_2017->GetBinContent(h_Ele_Reco_highPT_2017->GetXaxis()->FindBin(SCeta),h_Ele_Reco_highPT_2017->GetYaxis()->FindBin(std::min(pt,499.f)));
            RecoSF_Unc = h_Ele_Reco_highPT_2017->GetBinError(h_Ele_Reco_highPT_2017->GetXaxis()->FindBin(SCeta),h_Ele_Reco_highPT_2017->GetYaxis()->FindBin(std::min(pt,499.f)));
         }
      }
      else {
         edm::LogError("LeptonSFHelper::") << "Ele SFs for " << year << " is not supported!";
         abort();
         }
      
      
      if(year == 2016)
      {
         if(isCrack)
         {
            SelSF = h_Ele_Cracks_2016->GetBinContent(h_Ele_Cracks_2016->FindFixBin(SCeta, std::min(pt,199.f)));
            SelSF_Unc = h_Ele_Cracks_2016->GetBinError(h_Ele_Cracks_2016->FindFixBin(SCeta, std::min(pt,199.f)));
         }
         else
         {
            SelSF = h_Ele_notCracks_2016->GetBinContent(h_Ele_notCracks_2016->FindFixBin(SCeta, std::min(pt,199.f)));
            SelSF_Unc = h_Ele_notCracks_2016->GetBinError(h_Ele_notCracks_2016->FindFixBin(SCeta, std::min(pt,199.f)));
         }
         
      }
      
      else if(year == 2017)
      {
         if(isCrack)
         {
            SelSF = h_Ele_Cracks_2017->GetBinContent(h_Ele_Cracks_2017->FindFixBin(SCeta, std::min(pt,499.f)));
            SelSF_Unc = h_Ele_Cracks_2017->GetBinError(h_Ele_Cracks_2017->FindFixBin(SCeta, std::min(pt,499.f)));
         }
         else
         {
            SelSF = h_Ele_notCracks_2017->GetBinContent(h_Ele_notCracks_2017->FindFixBin(SCeta, std::min(pt,499.f)));
            SelSF_Unc = h_Ele_notCracks_2017->GetBinError(h_Ele_notCracks_2017->FindFixBin(SCeta, std::min(pt,499.f)));
         }
      }
      
      else if(year == 2018)
      {
         if(isCrack)
         {
            SelSF = h_Ele_Cracks_2018->GetBinContent(h_Ele_Cracks_2018->FindFixBin(SCeta, std::min(pt,499.f)));
            SelSF_Unc = h_Ele_Cracks_2018->GetBinError(h_Ele_Cracks_2018->FindFixBin(SCeta, std::min(pt,499.f)));
         }
         else
         {
            SelSF = h_Ele_notCracks_2018->GetBinContent(h_Ele_notCracks_2018->FindFixBin(SCeta, std::min(pt,499.f)));
            SelSF_Unc = h_Ele_notCracks_2018->GetBinError(h_Ele_notCracks_2018->FindFixBin(SCeta, std::min(pt,499.f)));
         }
      }
      
      
      else {
         edm::LogError("LeptonSFHelper::") << "Ele SFs for " << year << " is not supported!";
         abort();
      }
      
      SFError = sqrt( RecoSF_Unc*RecoSF_Unc/(RecoSF*RecoSF) + SelSF_Unc*SelSF_Unc/(SelSF*SelSF) ); // assume full correlation between different electrons (and uncorrelated reco and sel uncertainties)
   }
   
   //Muon SF
   if(abs(flav) == 13 )
   {
      if(year == 2016)
      {
         SelSF = h_Mu_SF_2016->GetBinContent(h_Mu_SF_2016->GetXaxis()->FindBin(eta),h_Mu_SF_2016->GetYaxis()->FindBin(std::min(pt,199.f))); //last bin contains the overflow
         SelSF_Unc = h_Mu_Unc_2016->GetBinContent(h_Mu_Unc_2016->GetXaxis()->FindBin(eta),h_Mu_Unc_2016->GetYaxis()->FindBin(std::min(pt,199.f))); //last bin contains the overflow
      }
      else if(year == 2017)
      {
         SelSF = h_Mu_SF_2017->GetBinContent(h_Mu_SF_2017->GetXaxis()->FindBin(eta),h_Mu_SF_2017->GetYaxis()->FindBin(std::min(pt,199.f))); //last bin contains the overflow
         SelSF_Unc = h_Mu_Unc_2017->GetBinContent(h_Mu_Unc_2017->GetXaxis()->FindBin(eta),h_Mu_Unc_2017->GetYaxis()->FindBin(std::min(pt,199.f))); //last bin contains the overflow
      }
      else if(year == 2018)
      {
         SelSF = h_Mu_SF_2018->GetBinContent(h_Mu_SF_2018->GetXaxis()->FindBin(eta),h_Mu_SF_2018->GetYaxis()->FindBin(std::min(pt,199.f))); //last bin contains the overflow
         SelSF_Unc = h_Mu_Unc_2018->GetBinContent(h_Mu_Unc_2018->GetXaxis()->FindBin(eta),h_Mu_Unc_2018->GetYaxis()->FindBin(std::min(pt,199.f))); //last bin contains the overflow
      }
      else {
         edm::LogError("LeptonSFHelper::") << "Ele SFs for " << year << " is not supported!";
         abort();
      }
      
      SFError = SelSF_Unc/SelSF; // assume full correlation between different muons (and uncorrelated reco and sel uncertainties)
   }
   
   return SFError;
}

