// -*- C++ -*-
//
// Package:    ZJetsFilterMerger
// Class:      ZJetsFilterMerger
// 
/**\class ZJetsFilterMerger ZJetsFilterMerger.cc psi2s1s/ZJetsFilterMerger/src/ZJetsFilterMerger.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Roberto Covarelli - TORINO
//         Created:  Tue Nov 22 20:39:54 CST 2011
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <iostream>

using namespace edm;
using namespace std;

//
// class declaration
//

class ZJetsFilterMerger : public edm::EDFilter {
   public:
      explicit ZJetsFilterMerger(const edm::ParameterSet&);
      ~ZJetsFilterMerger();


   private:

      virtual bool filter(edm::Event&, const edm::EventSetup&) override;

      // ----------member data ---------------------------
    
       int option;             // option = 0 : pass all
                               // option = 1 : lheNb==0 & nGenStatus2bHad==0 (use with DYJetsLL)
                               // option = 2 : lheNb>0 (use with DYBJetsLL)
                               // option = 3 : nGenStatus2bHad==0 (use with DYJetsToLL_BGenFilter)
       edm::EDGetTokenT<LHEEventProduct> LHEtoken_;
       edm::EDGetTokenT<std::vector<reco::GenParticle> > Gentoken_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ZJetsFilterMerger::ZJetsFilterMerger(const edm::ParameterSet& iConfig):
option(iConfig.getUntrackedParameter("option", 0))
{
  LHEtoken_ = consumes<LHEEventProduct>(iConfig.getParameter< edm::InputTag > ("theLHESrc"));
  Gentoken_ = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter< edm::InputTag > ("theGenSrc"));

}


ZJetsFilterMerger::~ZJetsFilterMerger()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ZJetsFilterMerger::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   if (option == 0) return true;

   int nBLHE = 0;
   int nBGen = 0;

   edm::Handle< LHEEventProduct > EvtHandle ;
   iEvent.getByToken( LHEtoken_ , EvtHandle ) ;

   edm::Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByToken( Gentoken_ , genParticles);

   for (int i = 0; i < EvtHandle->hepeup().NUP; ++i) {
     if (EvtHandle->hepeup().ISTUP[i] != 1) {
       continue;
     }
     if (abs(EvtHandle->hepeup().IDUP[i]) == 5) nBLHE++;
   } 

   for(std::vector<reco::GenParticle>::const_iterator genParticle=genParticles->begin(); genParticle!=genParticles->end(); ++genParticle){
     if( ( abs(genParticle->pdgId())/100 == 5 || abs(genParticle->pdgId())/1000 == 5)) {
     cout << genParticle->pdgId() << " " << genParticle->status() << endl;
    }

// && genParticle->status()==2) nBGen++;
   }

   if (option == 1 && nBLHE == 0 && nBGen == 0) return true;   
   if (option == 2 && nBLHE > 0) return true;   
   if (option == 3 && nBGen > 0) return true;   
   return false;   
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
//define this as a plug-in
DEFINE_FWK_MODULE(ZJetsFilterMerger);
