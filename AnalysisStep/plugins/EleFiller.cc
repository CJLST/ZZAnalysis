/** \class EleFiller
 *
 *  No description available.
 *
 *  $Date: 2013/05/24 15:42:42 $
 *  $Revision: 1.28 $
 *  \author N. Amapane (Torino)
 *  \author S. Bolognesi (JHU)
 *  \author C. Botta (CERN)
 *  \author S. Casasso (Torino)
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>

#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>

#include <vector>
#include <string>
#include "TRandom3.h"

using namespace edm;
using namespace std;
using namespace reco;



class EleFiller : public edm::EDProducer {
 public:
  /// Constructor
  explicit EleFiller(const edm::ParameterSet&);
    
  /// Destructor
  ~EleFiller();

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  edm::EDGetTokenT<pat::ElectronRefVector> electronToken;
  int sampleType;
  int setup;
  const StringCutObjectSelector<pat::Electron, true> cut;
  const CutSet<pat::Electron> flags;
  edm::EDGetTokenT<double> rhoToken;
  edm::EDGetTokenT<vector<Vertex> > vtxToken;
  EDGetTokenT<ValueMap<float> > BDTValueMapToken;
  string correctionFile;
  TRandom3 rgen_;
};


EleFiller::EleFiller(const edm::ParameterSet& iConfig) :
  electronToken(consumes<pat::ElectronRefVector>(iConfig.getParameter<edm::InputTag>("src"))),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  cut(iConfig.getParameter<std::string>("cut")),
  flags(iConfig.getParameter<ParameterSet>("flags")),
  correctionFile(iConfig.getParameter<std::string>("correctionFile")),
  rgen_(0)
{
  rhoToken = consumes<double>(LeptonIsoHelper::getEleRhoTag(sampleType, setup));
  vtxToken = consumes<vector<Vertex> >(edm::InputTag("goodPrimaryVertices"));

  BDTValueMapToken = consumes<ValueMap<float> >(iConfig.getParameter<InputTag>("mvaValuesMap"));

  produces<pat::ElectronCollection>();
	
}
EleFiller::~EleFiller(){
}


void
EleFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  

  // Get leptons and rho
  edm::Handle<pat::ElectronRefVector> electronHandle;
  iEvent.getByToken(electronToken, electronHandle);

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken, rhoHandle);
  double rho = *rhoHandle;

  edm::Handle<vector<Vertex> > vertices;
  iEvent.getByToken(vtxToken,vertices);

  Handle<ValueMap<float> > BDTValues;
  iEvent.getByToken(BDTValueMapToken, BDTValues);

  // Output collection
  auto result = std::make_unique<pat::ElectronCollection>();

  for (unsigned int i = 0; i< electronHandle->size(); ++i){

    //---Clone the pat::Electron
    pat::Electron l(*((*electronHandle)[i].get()));

    //--- PF ISO
    // for cone size R=0.4 :
    //float PFChargedHadIso   = l.chargedHadronIso();
    //float PFNeutralHadIso   = l.neutralHadronIso();
    //float PFPhotonIso       = l.photonIso();
    // for cone size R=0.3 :
    float PFChargedHadIso   = l.pfIsolationVariables().sumChargedHadronPt;
    float PFNeutralHadIso   = l.pfIsolationVariables().sumNeutralHadronEt;
    float PFPhotonIso       = l.pfIsolationVariables().sumPhotonEt;

    float SCeta = l.superCluster()->eta(); 
    float fSCeta = fabs(SCeta);

    float combRelIsoPF = LeptonIsoHelper::combRelIsoPF(sampleType, setup, rho, l);

    //--- SIP, dxy, dz
    float IP      = fabs(l.dB(pat::Electron::PV3D));
    float IPError = l.edB(pat::Electron::PV3D);
    float SIP     = IP/IPError;

    float dxy = 999.;
    float dz  = 999.;
    const Vertex* vertex = 0;
    if (vertices->size()>0) {
      vertex = &(vertices->front());
      dxy = fabs(l.gsfTrack()->dxy(vertex->position()));
      dz  = fabs(l.gsfTrack()->dz(vertex->position()));
    } 

	  
    // RunII BDT ID
    float BDT = 0.;
	
	//BDT running VID
	 BDT = (*BDTValues)[(*electronHandle)[i]];
	//cout << "BDT = " << (*BDTValues)[(*electronHandle)[i]] << endl;
    
    float pt = l.pt();

	  
    bool isBDT = false;

	if (setup==2016 || setup==2017 || setup==2018) 
	{
	  //WP taken from https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/ElectronIdentification/python/Identification/mvaElectronID_Fall17_iso_V2_cff.py#L21 
	 	isBDT         = (pt<=10 && ((fSCeta<0.8                  && BDT >  1.26402092475) ||
                                  (fSCeta>=0.8 && fSCeta<1.479 && BDT >  1.17808089508) ||
                                  (fSCeta>=1.479               && BDT >  1.33051972806)))
                 	 || (pt>10  && ((fSCeta<0.8                  && BDT >  2.36464785939) ||
                                  (fSCeta>=0.8 && fSCeta<1.479 && BDT >  2.07880614597) ||
                                  (fSCeta>=1.479               && BDT >  1.08080644615)));
	}
	else
	{
		std::cerr << "[ERROR] EleFiller: no BDT setup for: " << setup << " year!" << std::endl;
	}

    //-- Missing hit  
	 int missingHit;
	 missingHit = l.gsfTrack()->hitPattern().numberOfAllHits(HitPattern::MISSING_INNER_HITS);
	 
    //-- Flag for crack electrons (which use different efficiency SFs)
    bool isCrack = l.isGap(); 
    //--- Trigger matching
    int HLTMatch = 0; //FIXME
	 
	 
	 float scaleErr;
	 float smear_err_up;
	 
	 scaleErr= 1.;
	 smear_err_up = 1.;

	  
    //--- Embed user variables
    l.addUserFloat("PFChargedHadIso",PFChargedHadIso);
    l.addUserFloat("PFNeutralHadIso",PFNeutralHadIso);
    l.addUserFloat("PFPhotonIso",PFPhotonIso);
    l.addUserFloat("combRelIsoPF",combRelIsoPF);
    l.addUserFloat("SCeta",SCeta);
    l.addUserFloat("rho",rho);
    l.addUserFloat("SIP",SIP);
    l.addUserFloat("dxy",dxy);
    l.addUserFloat("dz",dz);
    l.addUserFloat("BDT",BDT);    
    l.addUserFloat("isBDT",isBDT);
    l.addUserFloat("isCrack",isCrack);
    l.addUserFloat("HLTMatch", HLTMatch);
    l.addUserFloat("missingHit", missingHit);
    l.addUserFloat("scale_unc",scaleErr);
    l.addUserFloat("smear_unc",smear_err_up);

    //--- MC parent code 
//     MCHistoryTools mch(iEvent);
//     if (mch.isMC()) {
//       int MCParentCode = 0;
//       //      int MCParentCode = mch.getParentCode(&l); //FIXME: does not work on cmg
//       l.addUserFloat("MCParentCode",MCParentCode);
//     }

    //--- Check selection cut. Being done here, flags are not available; but this way we 
    //    avoid wasting time on rejected leptons.
    if (!cut(l)) continue;

    //--- Embed flags (ie flags specified in the "flags" pset)
    for(CutSet<pat::Electron>::const_iterator flag = flags.begin(); flag != flags.end(); ++flag) {
      l.addUserFloat(flag->first,int((*(flag->second))(l)));
    }

    result->push_back(l);
  }
  iEvent.put(std::move(result));
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(EleFiller);

