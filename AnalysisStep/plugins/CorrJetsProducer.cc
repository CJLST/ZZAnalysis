// -*- C++ -*-
//
// Apply JEC to the pruned jet mass
//
// https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging#Recipes_to_apply_JEC_on_the_prun

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "PhysicsTools/PatAlgos/plugins/JetCorrFactorsProducer.h"

class CorrJetsProducer : public edm::EDProducer {
   public:
      explicit CorrJetsProducer(const edm::ParameterSet&);
      ~CorrJetsProducer();

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&) override;

      edm::EDGetTokenT<pat::JetCollection>         jetToken_;
      edm::EDGetTokenT<reco::VertexCollection>  vertexToken_;
      edm::EDGetTokenT<double>                          rho_;
      std::string                                   payload_;
      bool                                           isData_;
};

CorrJetsProducer::CorrJetsProducer(const edm::ParameterSet& iConfig):
    jetToken_   ( consumes<pat::JetCollection>     ( iConfig.getParameter<edm::InputTag>( "jets"    ) ) ),
    vertexToken_( consumes<reco::VertexCollection> ( iConfig.getParameter<edm::InputTag>( "vertex"  ) ) ),
    rho_        ( consumes<double>                 ( iConfig.getParameter<edm::InputTag>( "rho"     ) ) ),
    payload_    (                                    iConfig.getParameter<std::string>  ( "payload"   ) ),
    isData_     (                                    iConfig.getParameter<bool>         ( "isData"    ) )
{
    produces< pat::JetCollection >( "corrJets" );
}

CorrJetsProducer::~CorrJetsProducer() {}

void CorrJetsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;

    auto_ptr< pat::JetCollection > corrJets( new pat::JetCollection );

    Handle<pat::JetCollection>   jets;
    iEvent.getByToken(jetToken_, jets);
    size_t mult = jets->size();

    Handle<reco::VertexCollection>  vertices;
    iEvent.getByToken(vertexToken_, vertices);
    size_t nVtx = vertices->size();

    Handle< double >        rhoHandle;
    iEvent.getByToken(rho_, rhoHandle);
    double rho =           *rhoHandle;

    // Jet Energy Corrections
    ESHandle<JetCorrectorParametersCollection>       JetCorParColl;
    iSetup.get<JetCorrectionsRecord>().get(payload_, JetCorParColl); 
    vector<string> jecAK8PayloadNames;
    jecAK8PayloadNames.push_back("L2Relative");
    jecAK8PayloadNames.push_back("L3Absolute");
    if(isData_)
    jecAK8PayloadNames.push_back("L2L3Residual");
    vector<JetCorrectorParameters> vPar;
    for ( vector<std::string>::const_iterator payloadBegin = jecAK8PayloadNames.begin(), 
                                              payloadEnd   = jecAK8PayloadNames.end()  , 
                                              ipayload     = payloadBegin; 
                                              ipayload    != payloadEnd; 
                                            ++ipayload )  vPar.push_back( (*JetCorParColl)[*ipayload] );
    // Make the FactorizedJetCorrector
    shared_ptr<FactorizedJetCorrector> jecAK8 = shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) ); 
    jecAK8->setRho( rho  );
    jecAK8->setNPV( nVtx );
    for ( size_t i=0; i<mult; ++i ) {
        const pat::Jet& jet = (*jets)[i];
        jecAK8->setJetA  ( jet.jetArea()               );
        jecAK8->setJetPt ( jet.correctedP4(0).pt()     );
        jecAK8->setJetEta( jet.correctedP4(0).eta()    );
        jecAK8->setJetE  ( jet.correctedP4(0).energy() );
        float corrMass   = jet.userFloat("ak8PFJetsCHSPrunedMass") * jecAK8->getCorrection();
        pat::Jet* cloneJet = jet.clone();
        cloneJet->addUserFloat("ak8PFJetsCHSCorrPrunedMass", corrMass );
        corrJets->push_back( *cloneJet );
    } 
    iEvent.put(corrJets, "corrJets");
}

DEFINE_FWK_MODULE(CorrJetsProducer);
