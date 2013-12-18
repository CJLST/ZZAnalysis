
//
// $Id: RochesterMuonCorrector.cc,v 1.6 2013/01/31 11:08:08 namapane Exp $
//

/**
  \class    modules::RochesterMuonCorrectorT RochesterMuonCorrectorT.h "MuonAnalysis/MuonAssociators/interface/RochesterMuonCorrectorT.h"
  \brief    Applies Rochester corrections to muons            
  \author   Giovanni Petrucciani
  \version  $Id: RochesterMuonCorrector.cc,v 1.6 2013/01/31 11:08:08 namapane Exp $
*/


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "ZZAnalysis/AnalysisStep/interface/rochcor.h"
#include "ZZAnalysis/AnalysisStep/interface/rochcor2012.h"

namespace modules {

  template<typename T>
  class RochesterMuonCorrectorT : public edm::EDProducer {
    public:
      explicit RochesterMuonCorrectorT(const edm::ParameterSet & iConfig);
      virtual ~RochesterMuonCorrectorT() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
      /// Labels for input collections
      edm::InputTag src_;
      bool is55X;

      /// Rochester corrector
      rochcor corrector_;
      rochcor2012 corrector12_;
  };

} // namespace

template<typename T>
modules::RochesterMuonCorrectorT<T>::RochesterMuonCorrectorT(const edm::ParameterSet & iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src"))
{
  is55X=false;
  char * base = getenv("CMSSW_VERSION");
  if (base!=NULL) {
    std::string baseDir(base);
    if (baseDir.find("CMSSW_5_")!=std::string::npos)
      is55X=true;
    
  }
  else {
    printf("NO CMSSW Version found\n");
  }

  if (iConfig.getUntrackedParameter<bool>("fakeSmearing", false)) {
    corrector_.usefakeSmear=true;
    corrector12_.usefakeSmear=true;
  }

  produces<std::vector<T> >(); 
}

template<typename T>
void 
modules::RochesterMuonCorrectorT<T>::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;
    using namespace std;

    Handle<View<T> > src;
    iEvent.getByLabel(src_, src);

    unsigned int nsrc = src->size();
    auto_ptr<vector<T> > out(new vector<T>());
    out->reserve(nsrc);

    unsigned int run = iEvent.id().run(); 

    for (unsigned int i = 0; i < nsrc; ++i) {
        T mu = (*src)[i];
	TLorentzVector p4(mu.px(),mu.py(),mu.pz(),mu.energy());


        if (run <100 && !is55X) { //Monte Carlo 2011
	  corrector_.momcor_mc(p4, mu.charge(), 0.0, 0);
	}
	else if  (run <100 && is55X) {
	  corrector12_.momcor_mc(p4, mu.charge(), 0.0, 0);
	}
	else if (run>100&&run<=180252) { //2011 Data
            corrector_.momcor_data(p4, mu.charge(), 0.0, run <= 173692 ? 0 : 1);
	}	 
	else  if (run>190000) { //2012 Data
            corrector12_.momcor_data(p4, mu.charge(), 0.0, 0.0);
	}	 

	math::XYZTLorentzVector newP4(p4.Px(),p4.Py(),p4.Pz(),p4.Energy());
	mu.setP4(newP4);


        out->push_back(mu);


    }

    iEvent.put(out);
}


namespace modules {
    //typedef modules::RochesterMuonCorrectorT<reco::Muon>  RochesterMuonCorrector;
    typedef modules::RochesterMuonCorrectorT<pat::Muon>   RochesterPATMuonCorrector;
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace modules;
//DEFINE_FWK_MODULE(RochesterMuonCorrector);
DEFINE_FWK_MODULE(RochesterPATMuonCorrector);
