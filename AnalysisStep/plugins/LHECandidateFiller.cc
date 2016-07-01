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
#include "TH2F.h"
#include "TFile.h"
#include "TLorentzVector.h"

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
  using namespace edm;
  using namespace lhef;

  Handle<LHEEventProduct> lhe_evt;
  iEvent.getByLabel(src_, lhe_evt);
  const lhef::HEPEUP hepeup_ = lhe_evt->hepeup();
  const int nup_ = hepeup_.NUP;
  const std::vector<int> idup_ = hepeup_.IDUP;
  const std::vector<lhef::HEPEUP::FiveVector> pup_ = hepeup_.PUP;

  std::auto_ptr<pat::CompositeCandidate> result(new pat::CompositeCandidate);

  if (lhe_evt->pdf() != NULL) {
    std::cout << "PDF scale = " << std::setw(14) << std::fixed << lhe_evt->pdf()->scalePDF << std::endl;
    std::cout << "PDF 1 : id = " << std::setw(14) << std::fixed << lhe_evt->pdf()->id.first
      << " x = " << std::setw(14) << std::fixed << lhe_evt->pdf()->x.first
      << " xPDF = " << std::setw(14) << std::fixed << lhe_evt->pdf()->xPDF.first << std::endl;
    std::cout << "PDF 2 : id = " << std::setw(14) << std::fixed << lhe_evt->pdf()->id.second
      << " x = " << std::setw(14) << std::fixed << lhe_evt->pdf()->x.second
      << " xPDF = " << std::setw(14) << std::fixed << lhe_evt->pdf()->xPDF.second << std::endl;
  }
  std::cout << "Number of particles = " << nup_ << std::endl;
  for (unsigned int icount = 0; icount < (unsigned int)nup_; icount++) {
    std::cout << "# " << std::setw(14) << std::fixed << icount
      << std::setw(14) << std::fixed << idup_[icount]
      << std::setw(14) << std::fixed << (pup_[icount])[0]
      << std::setw(14) << std::fixed << (pup_[icount])[1]
      << std::setw(14) << std::fixed << (pup_[icount])[2]
      << std::setw(14) << std::fixed << (pup_[icount])[3]
      << std::setw(14) << std::fixed << (pup_[icount])[4]
      << std::endl;
  }
  if (lhe_evt->weights().size()) {
    std::cout << "weights:" << std::endl;
    for (size_t iwgt = 0; iwgt < lhe_evt->weights().size(); ++iwgt) {
      const LHEEventProduct::WGT& wgt = lhe_evt->weights().at(iwgt);
      std::cout << "\t" << wgt.id << ' '
        << std::scientific << wgt.wgt << std::endl;
    }
  }

  //result.addUserFloat("pbbh_VAJHU_dn", pbbh_VAJHU_dn);
  iEvent.put(result);
}



#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(LHECandidateFiller);

