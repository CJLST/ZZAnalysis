// Imported from CMGTools/Common/plugins/PATPFMuonEmbedder.cc rev. 1.1
// Original Id: PATPFMuonEmbedder.cc,v 1.1 2012/05/06 08:57:26 cbern Exp 
//
// Original Author:  Michail Bachtis
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#if CMSSW_VERSION>500
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "Math/GenVector/VectorUtil.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"

#include <iostream>
//
// class decleration


template <typename T>
class TunePMuonEmbedder : public edm::EDProducer {
   public:

  

  explicit TunePMuonEmbedder(const edm::ParameterSet& iConfig):
    src_(iConfig.getParameter<edm::InputTag>("src"))
     {
       produces<std::vector<T> >();
     }

  ~TunePMuonEmbedder() {}
   private:



  virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
  {
    using namespace edm;
    using namespace reco;

    std::auto_ptr<std::vector<T> > out(new std::vector<T>);
    
    Handle<std::vector<T> > cands;
    if(iEvent.getByLabel(src_,cands)) 
      for(unsigned int  i=0;i!=cands->size();++i){
	T muon = cands->at(i);
	reco::TrackRef track = (muon::tevOptimized(muon.globalTrack(),muon.innerTrack(),muon.tpfmsTrack(),muon.pickyTrack(), 200., 17., 40., 0.25)).first;
	math::XYZTLorentzVector p4(track->px(),track->py(),track->pz(),sqrt(track->p()*track->p()+0.1057*0.1057));
	muon.setP4(p4);
	out->push_back(muon);
      }
    iEvent.put(out);
  }
     

      // ----------member data ---------------------------
      edm::InputTag src_;

};



typedef TunePMuonEmbedder<reco::Muon> TunePRecoMuonEmbedder;
typedef TunePMuonEmbedder<pat::Muon> TunePPatMuonEmbedder;
 
DEFINE_FWK_MODULE(TunePRecoMuonEmbedder);
DEFINE_FWK_MODULE(TunePPatMuonEmbedder);

#endif
