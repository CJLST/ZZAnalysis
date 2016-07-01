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

  edm::InputTag src_;
};


LHECandidateFiller::LHECandidateFiller(const edm::ParameterSet& iConfig) :
src_(iConfig.getParameter<edm::InputTag>("src"))
{
  produces<pat::CompositeCandidate>();
}


void LHECandidateFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  std::auto_ptr<pat::CompositeCandidate> result(new pat::CompositeCandidate);

  edm::Handle<LHEEventProduct> lhe_evt;
  iEvent.getByLabel(src_, lhe_evt);
  const lhef::HEPEUP hepeup_ = lhe_evt->hepeup();
  const int nup_ = hepeup_.NUP;
  const std::vector<int> istup_ = hepeup_.ISTUP;
  const std::vector<std::pair<int, int>> mothup_ = hepeup_.MOTHUP;
  const std::vector<int> idup_ = hepeup_.IDUP;
  const std::vector<lhef::HEPEUP::FiveVector> pup_ = hepeup_.PUP;


  // PDF scale
  float LHE_PDFScale = -1;
  if (lhe_evt->pdf()!=NULL) LHE_PDFScale = lhe_evt->pdf()->scalePDF;
  result.addUserFloat("LHE_PDFScale", LHE_PDFScale);
  //

  // Mothers, daughters (+ associated particles)
  std::vector<int> idMother;
  std::vector<TLorentzVector> pMother;
  std::vector<int> idDaughter;
  std::vector<TLorentzVector> pDaughter;
  for (int ipart = 0; ipart<nup_; ipart++){
    if (istup_.at(ipart)==-1){
      idMother.push_back(idup_.at(ipart));
      TLorentzVector tmpMom(pup_.at(ipart)[0], pup_.at(ipart)[1], pup_.at(ipart)[2], pup_.at(ipart)[3]);
      pMother.push_back(tmpMom);
    }
    else if(istup_.at(ipart)==1){
      idDaughter.push_back(idup_.at(ipart));
      TLorentzVector tmpMom(pup_.at(ipart)[0], pup_.at(ipart)[1], pup_.at(ipart)[2], pup_.at(ipart)[3]);
      pDaughter.push_back(tmpMom);
    }
  }
  int nMothers=idMother.size();
  int nDaughters=idDaughter.size();
  result.addUserInt("nLHEMothers", nMothers);
  result.addUserInt("nDaughters", nDaughters);
  for (int imot=0; imot<nMothers; imot++){
    result.addUserFloat();
  }



  // TBC

  // LHE weights
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
  result.addUserFloat("LHEweight_QCDscale_muR1_muF1", LHEweight_QCDscale_muR1_muF1);
  result.addUserFloat("LHEweight_QCDscale_muR1_muF2", LHEweight_QCDscale_muR1_muF2);
  result.addUserFloat("LHEweight_QCDscale_muR1_muF0p5", LHEweight_QCDscale_muR1_muF0p5);
  result.addUserFloat("LHEweight_QCDscale_muR2_muF1", LHEweight_QCDscale_muR2_muF1);
  result.addUserFloat("LHEweight_QCDscale_muR2_muF2", LHEweight_QCDscale_muR2_muF2);
  result.addUserFloat("LHEweight_QCDscale_muR2_muF0p5", LHEweight_QCDscale_muR2_muF0p5);
  result.addUserFloat("LHEweight_QCDscale_muR0p5_muF1", LHEweight_QCDscale_muR0p5_muF1);
  result.addUserFloat("LHEweight_QCDscale_muR0p5_muF2", LHEweight_QCDscale_muR0p5_muF2);
  result.addUserFloat("LHEweight_QCDscale_muR0p5_muF0p5", LHEweight_QCDscale_muR0p5_muF0p5);
  //

  iEvent.put(result);
}



// TBC
void LHECandidateFiller::addPtEtaPhiMassIdData(string owner, int id, TLorentzVector mom, bool usePz){
  vector<string> tmpBranchList;
  getPtEtaPhiMassIdDataStrings(owner, usePz);
  for (unsigned int b=0; b<tmpBranchList.size(); b++){
    string branchname = tmpBranchList.at(b);
    if(!addId || branchname.find("Id")==string::npos) reserveBranch(tmpBranchList.at(b), btype, doSetAddress);
    else{
      BaseTree::BranchTypes bInttype = BaseTree::bInt;
      if (btype==BaseTree::bVectorDouble) bInttype = BaseTree::bVectorInt;
      reserveBranch(tmpBranchList.at(b), bInttype, doSetAddress);
    }
  }
}
void LHECandidateFiller::getPtEtaPhiMIdBranches(vector<string>& blist, string owner, bool usePz){
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

