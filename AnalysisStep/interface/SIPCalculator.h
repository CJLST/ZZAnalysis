#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/PatCandidates/interface/Muon.h"


class SIPCalculator {
public:
  SIPCalculator() {

  }

  ~SIPCalculator() {}
  
  void initialize(const edm::EventSetup& es) {
    es.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);

  }

  double calculate(const pat::Muon& muon,const reco::Vertex& PV) {
    if (muon.innerTrack().isNonnull() && PV.isValid()) {
      reco::TransientTrack tt = trackBuilder->build(muon.innerTrack());
      std::pair<bool,Measurement1D> result =IPTools::signedImpactParameter3D(tt,
         GlobalVector(muon.innerTrack()->px(),muon.innerTrack()->py(),muon.innerTrack()->pz()),
  	 PV);
      return fabs(result.second.value()/result.second.error());
    } else {
      return -99;
    }
  }


 private:
    edm::ESHandle<TransientTrackBuilder> trackBuilder;    

  

};
