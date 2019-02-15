#include <FWCore/ParameterSet/interface/FileInPath.h>

#include "../interface/BTaggingSFHelper.h"

#include <cmath>
#include "TString.h"
#include "TMath.h"
#include "TDirectory.h"

using namespace std;

BTaggingSFHelper::BTaggingSFHelper(std::string SFfilename, std::string effFileName)
{
    // Allow relative paths in python config file to be found in C++
    edm::FileInPath fip_sf(SFfilename);
    TDirectory* curdir = gDirectory;
    
    m_calib  = new BTagCalibration("CSVv2", fip_sf.fullPath().c_str()); // The CSVv2 designation doesn't matter if you are only reading, so don't bother to change it...
    BTagCalibrationReader* m_reader = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central");
    BTagCalibrationReader* m_reader_up = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "up");
    BTagCalibrationReader* m_reader_down = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "down");
    m_readers = std::vector<BTagCalibrationReader*>{ m_reader, m_reader_up, m_reader_down };
    for (BTagCalibrationReader*& mr:m_readers){
      mr->load(*m_calib, BTagEntry::FLAV_B, "comb");
      mr->load(*m_calib, BTagEntry::FLAV_C, "comb");
      mr->load(*m_calib, BTagEntry::FLAV_UDSG, "incl");
    }

    edm::FileInPath fip_eff(effFileName);
    m_fileEff = TFile::Open(fip_eff.fullPath().c_str(), "read");
    
    TString flavs[3] = {"b", "c", "udsg"};
    for(int flav=0; flav<3; flav++){
        TString name = Form("eff_%s_M_ALL",flavs[flav].Data());
        m_hEff[flav] = (TH1F*)m_fileEff->Get(name.Data());
    }
    curdir->cd(); // When m_fileEff is open, the current directory changes to m_fileEff in ROOT. Avoid this...
}

BTaggingSFHelper::~BTaggingSFHelper()
{
    if (m_fileEff && m_fileEff->IsOpen()) m_fileEff->Close(); // m_hEff[] are now invalid pointers!
    else if (m_fileEff) delete m_fileEff;
    for (BTagCalibrationReader*& mr:m_readers) delete mr;
    delete m_calib;
}

float BTaggingSFHelper::getSF(SFsyst syst, int jetFlavor, float pt, float eta)
{
    float SF = 1.0;
    
    BTagEntry::JetFlavor flav;
    if(abs(jetFlavor)==5){
        flav = BTagEntry::FLAV_B;
    }else if(abs(jetFlavor)==4){
        flav = BTagEntry::FLAV_C;
    }else{
        flav = BTagEntry::FLAV_UDSG;
    }

    constexpr float epsilon=1e-5;
    float myPt = pt;
    /*
    // Old code
    float MaxBJetPt = 669.9, MaxLJetPt = 999.9;  // value must be below the boundary
    float MaxJetEta = 2.5; // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
    */
	
    std::pair<float, float> eta_bounds = m_readers[central]->min_max_eta(flav, 0);
    if (eta_bounds.first==0.f) eta_bounds.first = -eta_bounds.second;
    if (eta_bounds.first>eta || eta>eta_bounds.second) return 1.;

    bool DoubleUncertainty = false;
    std::pair<float, float> pt_bounds = m_readers[central]->min_max_pt(flav, eta, 0); // Since this function uses eta, first determine eta bounds as above
    if (pt_bounds.first >= myPt){ myPt = pt_bounds.first+epsilon; DoubleUncertainty = true; }
    if (pt_bounds.second <= myPt){ myPt = pt_bounds.second-epsilon; DoubleUncertainty = true; }

    /*
    // Old code
    if((flav == BTagEntry::FLAV_B || flav == BTagEntry::FLAV_C) && pt>MaxBJetPt){
        myPt = MaxBJetPt;
        DoubleUncertainty = true;
    }
    if(flav == BTagEntry::FLAV_UDSG && pt>MaxLJetPt){
        myPt = MaxLJetPt;
        DoubleUncertainty = true;
    }
    */

    SF = m_readers[syst]->eval(flav, eta, myPt);
    
    if(DoubleUncertainty && syst!=central){
        float SFcentral = m_readers[central]->eval(flav, eta, myPt);
        SF = 2.f*(SF - SFcentral) + SFcentral;
    }

    // Old code
    //if (std::abs(eta) > MaxJetEta) SF = 1.; // Do not apply SF for jets with eta higher than the threshold
    return SF;
}


float BTaggingSFHelper::getEff(int jetFlavor, float pt, float eta)
{
    int flav;
    if(abs(jetFlavor)==5) flav = 0;
    else if(abs(jetFlavor)==4) flav = 1;
    else flav = 2;
    
    TH1F* h = m_hEff[flav];
    float aEta = TMath::Abs(eta);
    
    int binglobal = h->FindBin(pt, aEta);
    int binx, biny, binz;
    h->GetBinXYZ(binglobal, binx, biny, binz); // converts to x, y bins
    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();
    
    // under-overflows
    if (binx < 1) binx = 1;
    if (biny < 1) biny = 1;
    if (binx > nx) binx = nx;
    if (biny > ny) biny = ny;
    
    float eff = h->GetBinContent(binx, biny);
    
    // protection against wrongly measured efficiencies (low stat) --> reduce pT bin
    while(eff < 0.00000000001 && binx > 0){
        binx--;
        eff = h->GetBinContent(binx, biny);
    }
    
    return eff;
}




