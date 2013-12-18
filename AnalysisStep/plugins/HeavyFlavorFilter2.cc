// -*- C++ -*-
//
// Package:    HeavyFlavorFilter2
// Class:      HeavyFlavorFilter2
// 
/**\class HeavyFlavorFilter2 HeavyFlavorFilter2.cc ZbbAnalysis/AnalysisStep/src/HeavyFlavorFilter2.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco Musich ,Reviser: Stefano Casasso 40 2-A16,+41227671519,
//         Created:  Thu Apr 14 11:37:27 CEST 2011
// $Id: HeavyFlavorFilter2.cc,v 1.4 2012/09/11 12:12:45 namapane Exp $
//
//*****USAGE*****
//status2 (3): if True checks HF quarks with status 2 (3)
//hDaughterVeto: if True checks only HF quarks with strings or clusters (i.e. hadrons) as daughters
//zDaughterVeto: if True checks only HF quarks not having Z as daughter (i.e. not initial state partons)
//


// system include files
#include <memory>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// includes from Alberto's code
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Common/interface/GetProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"

//
// class declaration
//
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

//
// class declaration
//

class HeavyFlavorFilter2 : public edm::EDFilter {
   public:
      explicit HeavyFlavorFilter2(const edm::ParameterSet&);
      ~HeavyFlavorFilter2();

   private:
       virtual void beginJob() ;
       virtual bool filter(edm::Event&, const edm::EventSetup&);
       virtual void endJob() ;
  
      // ----------member data ---------------------------
  edm::InputTag genParticles_;
  bool status2_,status3_,hdau_veto_,zdau_veto_;
  double ptcut_;
  long wasRun;
  long wasAccept_b,wasAccept_c, wasAccept_tot;

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
HeavyFlavorFilter2::HeavyFlavorFilter2(const edm::ParameterSet& iConfig)
{
  genParticles_=iConfig.getParameter<edm::InputTag>("src");
  status3_=iConfig.getParameter<bool>("status3");
  status2_=iConfig.getParameter<bool>("status2");
  hdau_veto_=iConfig.getParameter<bool>("hDaughterVeto");
  zdau_veto_=iConfig.getParameter<bool>("zDaughterVeto");
  ptcut_=iConfig.getParameter<double>("ptcut");
  //now do what ever initialization is needed
  
}


HeavyFlavorFilter2::~HeavyFlavorFilter2()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
HeavyFlavorFilter2::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   wasRun++;
   using namespace edm;
   using namespace std;
   using namespace reco;
   Bool_t isb_,isc_;
   Bool_t isStatus2_,isStatus3_;
   //   Bool_t hasbdaughter_;
   //   Bool_t hascdaughter_;
   Bool_t hasZdaughter_;
   Bool_t hasHdaughter_;

   //Bool_t theFilter_b=false;
   //Bool_t theFilter_c=false;

   Int_t theFilter_b=0;
   Int_t theFilter_c=0;
   

   edm::Handle<GenParticleCollection> genParticlesCollection;
   iEvent.getByLabel(genParticles_, genParticlesCollection);

   for(GenParticleCollection::const_iterator genp = genParticlesCollection->begin();genp != genParticlesCollection->end(); ++ genp ) {  // loop over GEN particles
     
     isb_=false;
     isc_=false;
     isStatus2_=false;
     isStatus3_=false;
//      hasbdaughter_=false;
//      hascdaughter_=false;
     hasZdaughter_=false;
     hasHdaughter_=false;

     if (fabs(genp->pdgId())==4)isc_=true;
     if (fabs(genp->pdgId())==5)isb_=true;

     if (genp->status()==2)isStatus2_=true;
     if (genp->status()==3)isStatus3_=true;

     for( size_t i = 0; i < genp->numberOfDaughters(); i++) {
//        if( genp->daughter(i)->pdgId() == 5 ) hasbdaughter_ = true;
//        if( genp->daughter(i)->pdgId() == 4 ) hascdaughter_ = true;
       if( genp->daughter(i)->pdgId() == 23 ) hasZdaughter_ = true;
       if( genp->daughter(i)->pdgId() == 91 || genp->daughter(i)->pdgId() == 92 ) hasHdaughter_ = true;
     }
     
     if(isb_ || isc_){
     

       if (status2_ && status3_ && (isStatus2_ || isStatus3_)){
	if (zdau_veto_ && hdau_veto_){
	  if (hasZdaughter_ || hasHdaughter_==false)continue;
	  else if (genp->pt()>ptcut_){
	    if (isb_)theFilter_b++;
	    if (isc_)theFilter_c++;
	    //break;
	  }
	}

	else if (zdau_veto_){
	  if (hasZdaughter_)continue;
	  else if (genp->pt()>ptcut_){
	    if (isb_)theFilter_b++;
	    if (isc_)theFilter_c++;
	    //break;
	  }
	}

	else if (hdau_veto_){
	  if (hasHdaughter_==false)continue;
	  else if (genp->pt()>ptcut_){
	    if (isb_)theFilter_b++;
	    if (isc_)theFilter_c++;
	    //break;
	  }
	}

	else if (zdau_veto_==false && hdau_veto_==false){
	  if (genp->pt()>ptcut_){
	    if (isb_)theFilter_b++;
	    if (isc_)theFilter_c++;
	    //break;
	  }
	}
	
       }

       else if(status2_ && isStatus2_){
	if (hdau_veto_){
	  if (hasHdaughter_==false)continue;
	  else if (genp->pt()>ptcut_){
	    if (isb_)theFilter_b++;
	    if (isc_)theFilter_c++;
	    //break;
	  }
	}
	else {
	  if (genp->pt()>ptcut_){
	    if (isb_)theFilter_b++;
	    if (isc_)theFilter_c++;
	    //break;
	  }
	}
	
      }

       else if(status3_ && isStatus3_){
	if (zdau_veto_){
	  if (hasZdaughter_)continue;
	  else if (genp->pt()>ptcut_){
	    if (isb_)theFilter_b++;
	    if (isc_)theFilter_c++;
	    //break;
	  }
	}
	else {
	  if (genp->pt()>ptcut_){
	    if (isb_)theFilter_b++;
	    if (isc_)theFilter_c++;
	    //break;
	  }
	}
	
      }
      
     }//end isb or isc


   }//end loop over genParticles
   

   if(theFilter_b >=2) wasAccept_b++;
   if(theFilter_c >=2) wasAccept_c++;

   //   cout << "b " <<  theFilter_b << " c " << theFilter_c << endl;

   if( theFilter_b>=2 || theFilter_c >=2) {
     wasAccept_tot++;
     return true;
   }
   else return false;
   
}

// ------------ method called once each job just before starting event loop  ------------
void 
HeavyFlavorFilter2::beginJob()
{
  //edm::Service<TFileService> fileService;
  wasRun=0;
  wasAccept_b=0;
  wasAccept_c=0;
  wasAccept_tot=0;

  //h_genParticlesSize_ = fileService->make<TH1F>("h_genParticlesSize_","size of genParticles collection",2000,-0.5,1999.5); 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HeavyFlavorFilter2::endJob() {
  
  std::cout<<"***********************************"<<std::endl;
  std::cout<<"* HeavyFlavorFilter2::endJob() "<<std::endl;
  std::cout<<"***********************************"<<std::endl;
  std::cout<<"* Was run on       : "<<wasRun<<" events; accepted: " <<wasAccept_tot << std::endl;
  std::cout<<"*    for b quark: "<<std::endl;
  std::cout<<"*         Was accepted on  : "<<wasAccept_b<<" events"<<std::endl;
  std::cout<<"*         Filter efficiency: "<<(double(wasAccept_b)/wasRun)*100.<<"%"<<std::endl;
  std::cout<<"*    for c quark: "<<std::endl;
  std::cout<<"*         Was accepted on  : "<<wasAccept_c<<" events"<<std::endl;
  std::cout<<"*         Filter efficiency: "<<(double(wasAccept_c)/wasRun)*100.<<"%"<<std::endl;
  std::cout<<"***********************************"<<std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HeavyFlavorFilter2);
