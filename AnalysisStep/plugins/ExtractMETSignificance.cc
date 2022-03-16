/*
** class  : ExtractMETSignificance
** author : L. Cadamuro (LLR)
** date   : 4 November 2016
** brief  : takes pat::MET in input and produces a standalone significance and covariance collection from it
**          used as a replacement of the previous standalone significance producer
*/

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <DataFormats/METReco/interface/CommonMETData.h>

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;

class ExtractMETSignificance : public edm::EDProducer {
    public: 
        /// Constructor
        explicit ExtractMETSignificance(const edm::ParameterSet&);
        /// Destructor
        ~ExtractMETSignificance(){};  

    private:
        virtual void beginJob(){};  
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob(){};

        edm::EDGetTokenT<View<pat::MET>> theMETTag;
};

ExtractMETSignificance::ExtractMETSignificance(const edm::ParameterSet& iConfig) :
theMETTag(consumes<View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET")))
{
    produces<double>("METSignificance");
    produces<math::Error<2>::type>("METCovariance");
}

void ExtractMETSignificance::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    Handle<View<pat::MET> > METHandle;
    iEvent.getByToken(theMETTag, METHandle);
    const pat::MET& patMET = (*METHandle)[0];
    
    const reco::METCovMatrix cov = patMET.getSignificanceMatrix();
    double sig = patMET.significance();

    //std::auto_ptr<double> significance (new double);
    std::unique_ptr<double> significance (new double);
	(*significance) = sig;

    //std::auto_ptr<math::Error<2>::type> covPtr(new math::Error<2>::type());
	std::unique_ptr<math::Error<2>::type> covPtr(new math::Error<2>::type());
	(*covPtr)(0,0) = cov(0,0);
    (*covPtr)(1,0) = cov(1,0);
    (*covPtr)(1,1) = cov(1,1);

    //iEvent.put( covPtr, "METCovariance" );
    //iEvent.put( significance, "METSignificance" );
	iEvent.put( std::move(covPtr), "METCovariance" );
	iEvent.put( std::move(significance), "METSignificance" );
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ExtractMETSignificance);
