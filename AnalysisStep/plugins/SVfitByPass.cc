// bypasses SVfit and saves the same userfloats
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>
#include <DataFormats/Candidate/interface/ShallowCloneCandidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <DataFormats/METReco/interface/CommonMETData.h>
//#include <TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h> //FRA 2018 data: Not used anymore, we use ClassicSvfit now
#include <TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h>
#include "Math/LorentzVector.h"

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;
using METUncertainty = pat::MET::METUncertainty;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

// ------------------------------------------------------------------

class SVfitBypass : public edm::EDProducer {
 public:
  /// Constructor
  explicit SVfitBypass(const edm::ParameterSet&);
    
  /// Destructor
  ~SVfitBypass(){
  };  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};
  
  edm::EDGetTokenT<View<reco::CompositeCandidate> > theCandidateTag;
  // std::vector <edm::EDGetTokenT<View<pat::MET> > > vtheMETTag; // polymorphism of view --> good for reco::PFMET and pat::MET! 
  edm::EDGetTokenT<View<pat::MET>> theMETTag; // polymorphism of view --> good for reco::PFMET and pat::MET! 
  edm::EDGetTokenT<double> theSigTag;
  edm::EDGetTokenT<math::Error<2>::type> theCovTag;

  edm::EDGetTokenT<double> theMETdxUPTag;
  edm::EDGetTokenT<double> theMETdyUPTag;
  edm::EDGetTokenT<double> theMETdxDOWNTag;
  edm::EDGetTokenT<double> theMETdyDOWNTag;
  edm::EDGetTokenT<double> theMETdxUPEESTag;
  edm::EDGetTokenT<double> theMETdyUPEESTag;
  edm::EDGetTokenT<double> theMETdxDOWNEESTag;
  edm::EDGetTokenT<double> theMETdyDOWNEESTag;

  //int sampleType;
  bool _usePairMET;
  //bool _useMVAMET;
};

// ------------------------------------------------------------------



SVfitBypass::SVfitBypass(const edm::ParameterSet& iConfig):
theCandidateTag(consumes<View<reco::CompositeCandidate> >(iConfig.getParameter<InputTag>("srcPairs"))),
theMETTag(consumes<View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
theSigTag(consumes<double>(iConfig.getParameter<edm::InputTag>("srcSig"))),
theCovTag(consumes<math::Error<2>::type>(iConfig.getParameter<edm::InputTag>("srcCov"))),
theMETdxUPTag(consumes<double>(iConfig.getParameter<edm::InputTag>("METdxUP"))),
theMETdyUPTag(consumes<double>(iConfig.getParameter<edm::InputTag>("METdyUP"))),
theMETdxDOWNTag(consumes<double>(iConfig.getParameter<edm::InputTag>("METdxDOWN"))),
theMETdyDOWNTag(consumes<double>(iConfig.getParameter<edm::InputTag>("METdyDOWN"))),
theMETdxUPEESTag(consumes<double>(iConfig.getParameter<edm::InputTag>("METdxUP_EES"))),
theMETdyUPEESTag(consumes<double>(iConfig.getParameter<edm::InputTag>("METdyUP_EES"))),
theMETdxDOWNEESTag(consumes<double>(iConfig.getParameter<edm::InputTag>("METdxDOWN_EES"))),
theMETdyDOWNEESTag(consumes<double>(iConfig.getParameter<edm::InputTag>("METdyDOWN_EES")))
{
  //theCandidateTag = iConfig.getParameter<InputTag>("srcPairs");
  //_useMVAMET = iConfig.getUntrackedParameter<bool>("useMVAMET");

  _usePairMET = iConfig.getParameter<bool>("usePairMET");

  // const std::vector<edm::InputTag>& inMET = iConfig.getParameter<std::vector<edm::InputTag> >("srcMET");
  // for (std::vector<edm::InputTag>::const_iterator it = inMET.begin(); it != inMET.end(); ++it)
  // {      
  //   // vtheMETTag.emplace_back(consumes<edm::View<reco::MET> >(*it) );
  //   vtheMETTag.emplace_back(consumes<edm::View<pat::MET> >(*it) );
  // }

  produces<pat::CompositeCandidateCollection>();
}  



void SVfitBypass::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  

  // Get lepton pairs
  Handle<View<reco::CompositeCandidate> > pairHandle;
  iEvent.getByToken(theCandidateTag, pairHandle);
  
  unsigned int elNumber = pairHandle->size();
  
  // Output collection
  std::unique_ptr<pat::CompositeCandidateCollection> result( new pat::CompositeCandidateCollection );

  // get event pat MET to be saved in output
  double METx = 0.;
  double METy = 0.; 
  double uncorrMETx = -999.;
  double uncorrMETy = -999.; 
  TMatrixD covMET(2, 2);
  float significance = -999.;
  Handle<View<pat::MET> > METHandle;

  // JES shift
  double METx_UP_JES   = 0.;
  double METy_UP_JES   = 0.;
  double METx_DOWN_JES = 0.;
  double METy_DOWN_JES = 0.;

  // TES shift
  double METx_UP_TES = 0.;
  double METy_UP_TES = 0.;
  double METx_DOWN_TES = 0.;
  double METy_DOWN_TES = 0.;

  // EES shift
  double METx_UP_EES = 0.;
  double METy_UP_EES = 0.;
  double METx_DOWN_EES = 0.;
  double METy_DOWN_EES = 0.;

  if (!_usePairMET)
  {
    // iEvent.getByToken(vtheMETTag.at(0), METHandle);
    iEvent.getByToken(theMETTag, METHandle);
    const pat::MET& patMET = (*METHandle)[0];
    METx = patMET.px();
    METy = patMET.py();
    Handle<double> significanceHandle;
    Handle<math::Error<2>::type> covHandle;
    iEvent.getByToken (theSigTag, significanceHandle);
    iEvent.getByToken (theCovTag, covHandle);
    covMET[0][0] = (*covHandle)(0,0);
    covMET[1][0] = (*covHandle)(1,0);
    covMET[0][1] = covMET[1][0]; // (1,0) is the only one saved
    covMET[1][1] = (*covHandle)(1,1);
    significance = (float) (*significanceHandle);

     // now get the MET shifted UP/DOWN by the JES
     LorentzVector patMET_UP_JES   = patMET.shiftedP4(METUncertainty::JetEnUp);
     METx_UP_JES = patMET_UP_JES.px();
     METy_UP_JES = patMET_UP_JES.py();
     LorentzVector patMET_DOWN_JES = patMET.shiftedP4(METUncertainty::JetEnDown);
     METx_DOWN_JES = patMET_DOWN_JES.px();
     METy_DOWN_JES = patMET_DOWN_JES.py();

     // TES shift of MET
     Handle<double> METxUPHandle;
     Handle<double> METyUPHandle;
     Handle<double> METxDOWNHandle;
     Handle<double> METyDOWNHandle;

     iEvent.getByToken (theMETdxUPTag, METxUPHandle);
     iEvent.getByToken (theMETdyUPTag, METyUPHandle);
     iEvent.getByToken (theMETdxDOWNTag, METxDOWNHandle);
     iEvent.getByToken (theMETdyDOWNTag, METyDOWNHandle);

     METx_UP_TES   = *METxUPHandle;
     METy_UP_TES   = *METyUPHandle;
     METx_DOWN_TES = *METxDOWNHandle;
     METy_DOWN_TES = *METyDOWNHandle;

     // EES shift of MET (E->tau ES)
     Handle<double> METxUPEESHandle;
     Handle<double> METyUPEESHandle;
     Handle<double> METxDOWNEESHandle;
     Handle<double> METyDOWNEESHandle;

     iEvent.getByToken (theMETdxUPEESTag, METxUPEESHandle);
     iEvent.getByToken (theMETdyUPEESTag, METyUPEESHandle);
     iEvent.getByToken (theMETdxDOWNEESTag, METxDOWNEESHandle);
     iEvent.getByToken (theMETdyDOWNEESTag, METyDOWNEESHandle);

     METx_UP_EES   = *METxUPEESHandle;
     METy_UP_EES   = *METyUPEESHandle;
     METx_DOWN_EES = *METxDOWNEESHandle;
     METy_DOWN_EES = *METyDOWNEESHandle;
  }

  // loop on all the pairs
  for (unsigned int i = 0; i < elNumber; ++i)
  {
    if (_usePairMET)
    {
      // iEvent.getByToken(vtheMETTag.at(i), METHandle);
      iEvent.getByToken(theMETTag, METHandle);
      //metNumber = METHandle->size();

      // const PFMET* pfMET = (PFMET*) &((*METHandle)[0]) ; // all this to transform the type of the pointer!
      // const pat::MET* patMET = &((*METHandle)[0]);
      const pat::MET* patMET = &((*METHandle)[i]);
      const reco::METCovMatrix& covMETbuf = patMET->getSignificanceMatrix();
      significance = (float) patMET->significance();

      METx = patMET->px();
      METy = patMET->py();

      uncorrMETx = ( patMET->hasUserFloat("uncorrPx") ) ? patMET->userFloat("uncorrPx") : -999;
      uncorrMETy = ( patMET->hasUserFloat("uncorrPy") ) ? patMET->userFloat("uncorrPy") : -999;

      covMET[0][0] = covMETbuf(0,0);
      covMET[1][0] = covMETbuf(1,0);
      covMET[0][1] = covMETbuf(0,1);
      covMET[1][1] = covMETbuf(1,1);

     // now get the MET shifted UP/DOWN by the JES
     LorentzVector patMET_UP_JES   = patMET->shiftedP4(METUncertainty::JetEnUp);
     METx_UP_JES = patMET_UP_JES.px();
     METy_UP_JES = patMET_UP_JES.py();
     LorentzVector patMET_DOWN_JES = patMET->shiftedP4(METUncertainty::JetEnDown);
     METx_DOWN_JES = patMET_DOWN_JES.px();
     METy_DOWN_JES = patMET_DOWN_JES.py();

     // TES shift of MET
     Handle<double> METxUPHandle;
     Handle<double> METyUPHandle;
     Handle<double> METxDOWNHandle;
     Handle<double> METyDOWNHandle;

     iEvent.getByToken (theMETdxUPTag, METxUPHandle);
     iEvent.getByToken (theMETdyUPTag, METyUPHandle);
     iEvent.getByToken (theMETdxDOWNTag, METxDOWNHandle);
     iEvent.getByToken (theMETdyDOWNTag, METyDOWNHandle);

     METx_UP_TES   = *METxUPHandle;
     METy_UP_TES   = *METyUPHandle;
     METx_DOWN_TES = *METxDOWNHandle;
     METy_DOWN_TES = *METyDOWNHandle;

     // EES shift of MET (E->tau ES)
     Handle<double> METxUPEESHandle;
     Handle<double> METyUPEESHandle;
     Handle<double> METxDOWNEESHandle;
     Handle<double> METyDOWNEESHandle;

     iEvent.getByToken (theMETdxUPEESTag, METxUPEESHandle);
     iEvent.getByToken (theMETdyUPEESTag, METyUPEESHandle);
     iEvent.getByToken (theMETdxDOWNEESTag, METxDOWNEESHandle);
     iEvent.getByToken (theMETdyDOWNEESTag, METyDOWNEESHandle);

     METx_UP_EES   = *METxUPEESHandle;
     METy_UP_EES   = *METyUPEESHandle;
     METx_DOWN_EES = *METxDOWNEESHandle;
     METy_DOWN_EES = *METyDOWNEESHandle;
    }

    // Get the pair and the two leptons composing it
    const CompositeCandidate& pairBuf = (*pairHandle)[i];
    pat::CompositeCandidate pair(pairBuf);

    float SVfitMass = -999.;
    // add user floats: SVfit mass, met properties, etc..  
    pair.addUserFloat("SVfitMass", SVfitMass);
    pair.addUserFloat("SVfitTransverseMass", -999);
    pair.addUserFloat("SVfit_pt",  -999);
    pair.addUserFloat("SVfit_eta", -999);
    pair.addUserFloat("SVfit_phi", -999);
    pair.addUserFloat("SVfitMassUnc", -999);
    pair.addUserFloat("SVfitTransverseMassUnc", -999);
    pair.addUserFloat("SVfit_ptUnc",  -999);
    pair.addUserFloat("SVfit_etaUnc",  -999);
    pair.addUserFloat("SVfit_phiUnc",  -999);
    pair.addUserFloat("SVfit_METRho", -999);
    pair.addUserFloat("SVfit_METPhi", -999);
    pair.addUserFloat("MEt_px", (float) METx);
    pair.addUserFloat("MEt_py", (float) METy);
    pair.addUserFloat("uncorrMEt_px", (float) uncorrMETx);
    pair.addUserFloat("uncorrMEt_py", (float) uncorrMETy);
    pair.addUserFloat("MEt_cov00", (float) covMET[0][0]);
    pair.addUserFloat("MEt_cov01", (float) covMET[0][1]);
    pair.addUserFloat("MEt_cov10", (float) covMET[1][0]);
    pair.addUserFloat("MEt_cov11", (float) covMET[1][1]);
    pair.addUserFloat("MEt_significance", significance);
    pair.addUserFloat("MEt_px_UP_JES", (float) METx_UP_JES);
    pair.addUserFloat("MEt_py_UP_JES", (float) METy_UP_JES);
    pair.addUserFloat("MEt_px_DOWN_JES", (float) METx_DOWN_JES);
    pair.addUserFloat("MEt_py_DOWN_JES", (float) METy_DOWN_JES);
    pair.addUserFloat("MEt_px_UP_TES", (float) METx_UP_TES);
    pair.addUserFloat("MEt_py_UP_TES", (float) METy_UP_TES);
    pair.addUserFloat("MEt_px_DOWN_TES", (float) METx_DOWN_TES);
    pair.addUserFloat("MEt_py_DOWN_TES", (float) METy_DOWN_TES);
    pair.addUserFloat("MEt_px_UP_EES", (float) METx_UP_EES);
    pair.addUserFloat("MEt_py_UP_EES", (float) METy_UP_EES);
    pair.addUserFloat("MEt_px_DOWN_EES", (float) METx_DOWN_EES);
    pair.addUserFloat("MEt_py_DOWN_EES", (float) METy_DOWN_EES);

    result->push_back(pair);
  }
  iEvent.put(std::move(result));
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(SVfitBypass);

