#include <FWCore/ParameterSet/interface/FileInPath.h>

#include "../interface/BTaggingSFHelper.h"

#include "TString.h"
#include "TMath.h"

using namespace std;

BTaggingSFHelper::BTaggingSFHelper(std::string SFfilename, std::string effFileName)
{
    // Allow relative paths in python config file to be found in C++
    edm::FileInPath fip_sf(SFfilename);
    
    m_calib = new BTagCalibration("CSVv2", fip_sf.fullPath().c_str());
    
    m_reader    = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central");
    m_reader_up = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "up"     );
    m_reader_do = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "down"   );
    
    for(int i=0 ; i < 3; i++ ){
    m_readers[i][0] = m_reader;
    m_readers[i][1] = m_reader_up;
    m_readers[i][2] = m_reader_do;
    }
    
    for(int i=0 ; i < 3; i++ ){
        m_readers[i][0] = m_reader;
        m_readers[i][1] = m_reader_up;
        m_readers[i][2] = m_reader_do;
        if(i == 0){
            m_readers[i][0]->load(*m_calib, BTagEntry::FLAV_B, "comb");
            m_readers[i][1]->load(*m_calib, BTagEntry::FLAV_B, "comb");
            m_readers[i][2]->load(*m_calib, BTagEntry::FLAV_B, "comb");
        }
        else if(i == 1){
            m_readers[i][0]->load(*m_calib, BTagEntry::FLAV_C, "comb");
            m_readers[i][1]->load(*m_calib, BTagEntry::FLAV_C, "comb");
            m_readers[i][2]->load(*m_calib, BTagEntry::FLAV_C, "comb");
        }
        else if(i == 2){
            m_readers[i][0]->load(*m_calib, BTagEntry::FLAV_UDSG, "incl");
            m_readers[i][1]->load(*m_calib, BTagEntry::FLAV_UDSG, "incl");
            m_readers[i][2]->load(*m_calib, BTagEntry::FLAV_UDSG, "incl");
        }
    }
    
    edm::FileInPath fip_eff(effFileName);
    m_fileEff = new TFile(fip_eff.fullPath().c_str());
    
    TString flavs[3] = {"b", "c", "udsg"};
    for(int flav=0; flav<3; flav++){
        TString name = Form("eff_%s_M_ALL",flavs[flav].Data());
        m_hEff[flav] = (TH1F*)m_fileEff->Get(name.Data());
    }
}

BTaggingSFHelper::~BTaggingSFHelper()
{
    if (m_fileEff) delete m_fileEff;
}

float BTaggingSFHelper::getSF(SFsyst syst, int jetFlavor, float pt, float eta)
{
    float SF = 1.0;
    
    BTagEntry::JetFlavor flav;
    int myFlavIndex = -1; // indexes in the m_readers array
    int mySystIndex = (int) syst;
    
    if(abs(jetFlavor)==5){
        flav = BTagEntry::FLAV_B;
        myFlavIndex = 0;
    }else if(abs(jetFlavor)==4){
        flav = BTagEntry::FLAV_C;
        myFlavIndex = 1;
    }else{
        flav = BTagEntry::FLAV_UDSG;
        myFlavIndex = 2;
    }
    
    float myPt = pt;
    float MaxBJetPt = 669.9, MaxLJetPt = 999.9;  // value must be below the boundary
    bool DoubleUncertainty = false;
    if((myFlavIndex==0 || myFlavIndex==1) && pt>MaxBJetPt){
        myPt = MaxBJetPt;
        DoubleUncertainty = true;
    }
    if(myFlavIndex==2 && pt>MaxLJetPt){
        myPt = MaxLJetPt;
        DoubleUncertainty = true;
    }

    SF = m_readers[myFlavIndex][mySystIndex]->eval(flav, eta, myPt);
    
    if(DoubleUncertainty && syst!=central){
        float SFcentral = m_readers[myFlavIndex][central]->eval(flav, eta, myPt);
        SF = 2*(SF - SFcentral) + SFcentral;
    }
    
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




