/** \class LHECandidateFiller
 *
 *
 *  \author N. Amapane - Torino
 *  \author U. Sarica - JHU
 */


#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <FWCore/ParameterSet/interface/FileInPath.h>
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include <DataFormats/GeometryVector/interface/Point3DBase.h>

#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include <iomanip>
#include <iostream>
#include <string>
#include "TLorentzVector.h"
#include "TString.h"

using namespace std;


class LHECandidateFiller : public edm::EDProducer {
public:

  // Constructor
  explicit LHECandidateFiller(const edm::ParameterSet&);

  // Destructor
  ~LHECandidateFiller(){}

private:

  virtual void beginJob(){}
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){}

  virtual void add_PtEtaPhiMassId_Data(std::auto_ptr<pat::CompositeCandidate>& cand, string owner, int id, TLorentzVector mom, bool usePz);
  virtual void get_PtEtaPhiMassId_DataStrings(vector<string>& blist, string owner, bool usePz);

  bool isMC;

};


LHECandidateFiller::LHECandidateFiller(const edm::ParameterSet& iConfig) :
isMC(iConfig.getParameter<bool>("isMC"))
{
  consumesMany<LHEEventProduct>();
  produces<pat::CompositeCandidate>();
}


void LHECandidateFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  std::auto_ptr<pat::CompositeCandidate> cand(new pat::CompositeCandidate);

  if (isMC){
    edm::Handle<LHEEventProduct> lhe_evt;
    vector<edm::Handle<LHEEventProduct> > lhe_handles;
    iEvent.getManyByType(lhe_handles);
    if (lhe_handles.size()>0) lhe_evt = lhe_handles.front();

    if (lhe_evt.isValid()){
      // Setup
      const lhef::HEPEUP hepeup_ = lhe_evt->hepeup();
      const int nup_ = hepeup_.NUP;
      const std::vector<int> istup_ = hepeup_.ISTUP;
      const std::vector<std::pair<int, int>> mothup_ = hepeup_.MOTHUP;
      const std::vector<int> idup_ = hepeup_.IDUP;
      const std::vector<lhef::HEPEUP::FiveVector> pup_ = hepeup_.PUP;
      //

      // PDF scale
      float LHE_PDFScale = -1;
      if (lhe_evt->pdf()!=NULL) LHE_PDFScale = lhe_evt->pdf()->scalePDF;
      cand->addUserFloat("LHE_PDFScale", LHE_PDFScale);
      //

      // Mothers, daughters (+ associated particles)
      std::vector<std::pair<int, TLorentzVector>> pMother;
      std::vector<std::pair<int, TLorentzVector>> pDaughter;
      for (int ipart = 0; ipart<nup_; ipart++){
        if (istup_.at(ipart)==-1) pMother.push_back(std::pair<int, TLorentzVector>(idup_.at(ipart), TLorentzVector(pup_.at(ipart)[0], pup_.at(ipart)[1], pup_.at(ipart)[2], pup_.at(ipart)[3])));
        else if (istup_.at(ipart)==1){
          // TO BE IMPROVED DURING MIGRATION TO MELA V2:
          // - ADD MELACANDIDATE, LHEEVENT FROM LHEANALYZER
          // - DO SORTING, SELECTION TO FIND HIGGS/ZZ DAUGHTERS VS ASSOCIATED PARTICLES
          // - IMPLEMENT CONVERTLHE-STYLE COMPARATORS
          pDaughter.push_back(std::pair<int, TLorentzVector>(idup_.at(ipart), TLorentzVector(pup_.at(ipart)[0], pup_.at(ipart)[1], pup_.at(ipart)[2], pup_.at(ipart)[3])));
        }
      }
      //

      // Fill mothers
      int nMothers=pMother.size();
      cand->addUserInt("nLHEMothers", nMothers);
      for (int imot=0; imot<nMothers; imot++){
        TString strowner = Form("Mother%i", imot+1);
        string owner = strowner.Data();
        add_PtEtaPhiMassId_Data(cand, owner, pMother.at(imot).first, pMother.at(imot).second, true); // usePz=true since Pt=0
      }
      //

      // Fill daughters
      int nDaughters=pDaughter.size();
      cand->addUserInt("nLHEDaughters", nDaughters);
      for (int idau=0; idau<nDaughters; idau++){
        TString strowner = Form("Daughter%i", idau+1);
        string owner = strowner.Data();
        add_PtEtaPhiMassId_Data(cand, owner, pDaughter.at(idau).first, pDaughter.at(idau).second, false);
      }
      //

      // LHE weights
      // TO BE IMPROVED TO CONTAIN EVERYTHING
      float LHEweight_QCDscale_muR1_muF1=0;
      float LHEweight_QCDscale_muR1_muF2=0;
      float LHEweight_QCDscale_muR1_muF0p5=0;
      float LHEweight_QCDscale_muR2_muF1=0;
      float LHEweight_QCDscale_muR2_muF2=0;
      float LHEweight_QCDscale_muR2_muF0p5=0;
      float LHEweight_QCDscale_muR0p5_muF1=0;
      float LHEweight_QCDscale_muR0p5_muF2=0;
      float LHEweight_QCDscale_muR0p5_muF0p5=0;
      if (lhe_evt->weights().size()){
        if (lhe_evt->weights().size()>=9){
          LHEweight_QCDscale_muR1_muF1 = lhe_evt->weights().at(0).wgt / lhe_evt->originalXWGTUP(); // just for verification (should be 1)
          LHEweight_QCDscale_muR1_muF2 = lhe_evt->weights().at(1).wgt / lhe_evt->originalXWGTUP();
          LHEweight_QCDscale_muR1_muF0p5 = lhe_evt->weights().at(2).wgt / lhe_evt->originalXWGTUP();
          LHEweight_QCDscale_muR2_muF1 = lhe_evt->weights().at(3).wgt / lhe_evt->originalXWGTUP();
          LHEweight_QCDscale_muR2_muF2 = lhe_evt->weights().at(4).wgt / lhe_evt->originalXWGTUP();
          LHEweight_QCDscale_muR2_muF0p5 = lhe_evt->weights().at(5).wgt / lhe_evt->originalXWGTUP();
          LHEweight_QCDscale_muR0p5_muF1 = lhe_evt->weights().at(6).wgt / lhe_evt->originalXWGTUP();
          LHEweight_QCDscale_muR0p5_muF2 = lhe_evt->weights().at(7).wgt / lhe_evt->originalXWGTUP();
          LHEweight_QCDscale_muR0p5_muF0p5 = lhe_evt->weights().at(8).wgt / lhe_evt->originalXWGTUP();
        }
      }
      cand->addUserFloat("LHEweight_QCDscale_muR1_muF1", LHEweight_QCDscale_muR1_muF1);
      cand->addUserFloat("LHEweight_QCDscale_muR1_muF2", LHEweight_QCDscale_muR1_muF2);
      cand->addUserFloat("LHEweight_QCDscale_muR1_muF0p5", LHEweight_QCDscale_muR1_muF0p5);
      cand->addUserFloat("LHEweight_QCDscale_muR2_muF1", LHEweight_QCDscale_muR2_muF1);
      cand->addUserFloat("LHEweight_QCDscale_muR2_muF2", LHEweight_QCDscale_muR2_muF2);
      cand->addUserFloat("LHEweight_QCDscale_muR2_muF0p5", LHEweight_QCDscale_muR2_muF0p5);
      cand->addUserFloat("LHEweight_QCDscale_muR0p5_muF1", LHEweight_QCDscale_muR0p5_muF1);
      cand->addUserFloat("LHEweight_QCDscale_muR0p5_muF2", LHEweight_QCDscale_muR0p5_muF2);
      cand->addUserFloat("LHEweight_QCDscale_muR0p5_muF0p5", LHEweight_QCDscale_muR0p5_muF0p5);
      //
    }
  }
  iEvent.put(cand);
}



// TBC
void LHECandidateFiller::add_PtEtaPhiMassId_Data(std::auto_ptr<pat::CompositeCandidate>& cand, string owner, int id, TLorentzVector mom, bool usePz){
  vector<string> tmpBranchList;
  get_PtEtaPhiMassId_DataStrings(tmpBranchList, owner, usePz);
  for (unsigned int b=0; b<tmpBranchList.size(); b++){
    string branchname = tmpBranchList.at(b);
    if (branchname.find("Id")!=string::npos) cand->addUserInt(branchname, id);
    else if (branchname.find("Pt")!=string::npos) cand->addUserFloat(branchname, (float)mom.Pt());
    else if (branchname.find("Pz")!=string::npos) cand->addUserFloat(branchname, (float)mom.Z());
    else if (branchname.find("Eta")!=string::npos) cand->addUserFloat(branchname, (float)mom.Eta());
    else if (branchname.find("Phi")!=string::npos) cand->addUserFloat(branchname, (float)mom.Phi());
    else if (branchname.find("Mass")!=string::npos) cand->addUserFloat(branchname, (float)mom.M());
  }
}
void LHECandidateFiller::get_PtEtaPhiMassId_DataStrings(vector<string>& blist, string owner, bool usePz){
  string strGen = "LHE";
  vector<string> strtmp;

  strtmp.push_back("Pt");
  if (usePz) strtmp.push_back("Pz");
  else strtmp.push_back("Eta");
  strtmp.push_back("Phi");
  strtmp.push_back("Mass");
  strtmp.push_back("Id");

  for (unsigned int b=0; b<strtmp.size(); b++){
    string varname = strtmp.at(b);
    varname.insert(0, owner);
    varname.insert(0, strGen);
    blist.push_back(varname);
  }
}



#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(LHECandidateFiller);

