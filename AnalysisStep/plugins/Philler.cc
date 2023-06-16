/** \class Philler
 *
 *  Produce collection of photons to be used for TLEs
 *
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/one/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"




#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;



class Philler : public edm::one::EDProducer<> {
 public:
  /// Constructor
  explicit Philler(const edm::ParameterSet&);

  /// Destructor
  ~Philler(){
  };

 private:
  virtual void beginJob(){};
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  edm::EDGetTokenT<edm::View<reco::Photon>> photonToken;
  edm::EDGetTokenT<pat::ElectronCollection> electronToken;

  bool recomputeBDT = false;
  int sampleType;
  int setup;
  const StringCutObjectSelector<pat::Photon, true> cut;
  const CutSet<pat::Photon> flags;
  edm::EDGetTokenT<double> rhoToken;
  edm::EDGetTokenT<vector<Vertex> > vtxToken;
  EDGetTokenT<ValueMap<float> > BDTValueMapToken;
  EDGetTokenT<ValueMap<float> > BDTValueMapToken2;

};


Philler::Philler(const edm::ParameterSet& iConfig) :
  photonToken(consumes<edm::View<reco::Photon>>(iConfig.getParameter<edm::InputTag>("src"))),
  electronToken(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("srcElectron"))),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  cut(iConfig.getParameter<std::string>("cut")),
  flags(iConfig.getParameter<ParameterSet>("flags")),
  //myMVATrig(0),
  BDTValueMapToken(consumes<ValueMap<float> >(iConfig.getParameter<InputTag>("mvaValuesMap"))),
  BDTValueMapToken2(consumes<ValueMap<float> >(iConfig.getParameter<InputTag>("mvaValuesMap2")))

{
  rhoToken = consumes<double>(LeptonIsoHelper::getEleRhoTag(sampleType, setup));
  vtxToken = consumes<vector<Vertex> >(edm::InputTag("goodPrimaryVertices"));
  produces<pat::PhotonCollection>();

}


void
Philler::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<edm::View<reco::Photon>> photonHandle;
  iEvent.getByToken(photonToken, photonHandle);

  edm::Handle<pat::ElectronCollection> electronHandle;
  iEvent.getByToken(electronToken, electronHandle);



  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken, rhoHandle);
  double rho = *rhoHandle;

  edm::Handle<vector<Vertex> > vertices;
  iEvent.getByToken(vtxToken,vertices);

  Handle<ValueMap<float> > BDTValues;
  iEvent.getByToken(BDTValueMapToken, BDTValues);

  Handle<ValueMap<float> > BDTValues2;
  iEvent.getByToken(BDTValueMapToken2, BDTValues2);


  // Output collection
  auto result = std::make_unique<pat::PhotonCollection>();

  for (unsigned int i = 0; i< photonHandle->size(); ++i){

    pat::Photon l(*(photonHandle->ptrAt(i)));
    const auto ph_SC = l.parentSuperCluster();

    float min_ele_dR = 999;
    bool has_matching_ele = false;
    float curr_dR = 999;
    int min_dR_index = -1;
    for(size_t i_ele = 0; i_ele < electronHandle->size(); ++i_ele) {
      const auto ele_parent_SC = electronHandle->at(i_ele).parentSuperCluster();

      curr_dR = reco::deltaR(l, electronHandle->at(i_ele));
      if(ph_SC==ele_parent_SC) {
         has_matching_ele = true;
         min_dR_index = i_ele;
         min_ele_dR = curr_dR;

      }
//      if(curr_dR < min_ele_dR) {
//        min_ele_dR = curr_dR;
//        min_dR_index = i_ele;
//      }
    }

    if(has_matching_ele) continue;
//    if(min_ele_dR < 0.1) continue;
    //---Clone the pat::Photon
//    pat::Photon l(*((*photonHandle)[i].get()));

    //--- PF ISO
    // for cone size R=0.3 :
//    float PFChargedHadIso   = l.pfIsolationVariables().sumChargedHadronPt;
//    float PFNeutralHadIso   = l.pfIsolationVariables().sumNeutralHadronEt;
//    float PFPhotonIso       = l.pfIsolationVariables().sumPhotonEt;

    //float fSCeta = fabs(l.eta());
    float fSCeta = fabs(l.superCluster()->eta());

    float combRelIsoPF = LeptonIsoHelper::combRelIsoPF(sampleType, setup, rho, l);

    //--- SIP, dxy, dz
    //float IP      = 999.; //fabs(l.dB(pat::Photon::PV3D));
    //float IPError = 999.; //l.edB(pat::Photon::PV3D);
    float SIP     = 0.; //IP/IPError;

    //float dxy = 999.;
    //float dz  = 999.;

    // RunII BDT ID
    float BDT = 0.;

    BDT =(*BDTValues)[photonHandle->ptrAt(i)]; // l.userFloat("PhotonMVAEstimatorRun2Spring15NonTrig25nsV1Values");

    // WP for Spring15-based BDT as proposed in https://indico.cern.ch/event/439325/session/1/contribution/21/attachments/1156760/1663207/slides_20150918.pdf
    bool isBDT = false;
    bool pass_WP_2 = false;
    // temporary ID https://twiki.cern.ch/twiki/bin/view/CMS/MultivariatePhotonIdentificationRun2

    // WP for TLE ID v1
    if(l.isEB() && fSCeta < 0.8 && BDT > -0.68)
        isBDT = true;// 0.374;
    if(l.isEB() && fSCeta >= 0.8 && BDT > -0.66)
        isBDT = true;// 0.374;
    if(l.isEE() && BDT > -0.6)
        isBDT = true;// 0.336;

    // tighter WP for TLE ID v1
    if(l.isEB() && fSCeta < 0.8 && BDT > -0.379)
        pass_WP_2 = true;
    if(l.isEB() && fSCeta >= 0.8 && BDT > -0.369)
        pass_WP_2 = true;
    if(l.isEE() && BDT > -0.294)
        pass_WP_2 = true;


    // Photon 90% MVA ID
    /*
    if(l.isEB() && BDT > 0.374)
        isBDT = true;// 0.374;
    if(l.isEE() && BDT > 0.336)
        isBDT = true;// 0.336;
    */

   if(l.gsfTrack().isNonnull() && l.gsfTrack().isAvailable())
	  {
// This is no longer supported
//		#if CMSSW_VERSION_MAJOR < 9
//	   l.gsfTrack()->hitPattern().numberOfHits(HitPattern::MISSING_INNER_HITS);
//		#else
	   l.gsfTrack()->hitPattern().numberOfAllHits(HitPattern::MISSING_INNER_HITS);
//		#endif
	  }

    //-- Flag for crack photons (which use different efficiency SFs)
    bool isCrack = (l.isEBEEGap())||(l.isEBGap())||(l.isEEGap());// l.isGap();
    //--- Trigger matching
    int HLTMatch = 0; //FIXME

    const edm::Ptr<reco::Photon> phRecoPtr(photonHandle->ptrAt(i));
    auto superCluster = phRecoPtr->superCluster();

    /*
    float pfSCfbrem = 999.;

    if(superCluster.isAvailable() && superCluster.isNonnull() ) {
        //LogWarning("") << "in if";
        if (superCluster->clustersSize() > 1) {
            //LogWarning("") << "in 2nd if";
            //CaloCluster_iterator first = superCluster->clustersBegin() ;
            //LogWarning("") << "after it";
            pfSCfbrem = 0.; //(superCluster->energy()-(*first)->energy()) / superCluster->energy();
        } else {
            pfSCfbrem = 0.;
        }
    }*/
    //--- Embed user variables
    //l.addUserFloat("PFChargedHadIso",PFChargedHadIso);
    //l.addUserFloat("PFNeutralHadIso",PFNeutralHadIso);
    //l.addUserFloat("PFPhotonIso",PFPhotonIso);
    l.addUserFloat("combRelIsoPF",combRelIsoPF);
    l.addUserFloat("rho",rho);
    l.addUserFloat("SIP",SIP);
    //l.addUserFloat("dxy",dxy);
    //l.addUserFloat("dz",dz);
    l.addUserFloat("BDT",BDT);
    l.addUserFloat("isBDT",isBDT);
    l.addUserFloat("isCrack",isCrack);
    l.addUserFloat("HLTMatch", HLTMatch);

    int ele_is_ID = -1;
    float ele_ID = -1;
    float ele_is_ISO = -1;
    float ele_is_SIP = -1;
    int ele_charge = -2;
    float ele_pt = -1;
    if(min_dR_index != -1) {
        pat::Electron ele((electronHandle->at(min_dR_index)));
        if(ele.hasUserFloat("isBDT")) ele_is_ID = ele.userFloat("isBDT");
        if(ele.hasUserFloat("BDT")) ele_ID = ele.userFloat("BDT");
        if(ele.hasUserFloat("combRelIsoPF")) ele_is_ISO = ele.userFloat("combRelIsoPF");
        if(ele.hasUserFloat("SIP")) ele_is_SIP = ele.userFloat("SIP");
        ele_charge = ele.charge();
        ele_pt = ele.pt();
    }
    l.addUserFloat("min_ele_dR", min_ele_dR);
    l.addUserFloat("ele_has_matching_ele", has_matching_ele);
    l.addUserFloat("ele_ID", ele_ID);
    l.addUserFloat("ele_is_ID", ele_is_ID);
    l.addUserFloat("ele_ISO", ele_is_ISO);
    l.addUserFloat("ele_SIP", ele_is_SIP);
    l.addUserFloat("ele_charge", ele_charge);
    l.addUserFloat("ele_pt", ele_pt);

    float BDT2 = (*BDTValues2)[photonHandle->ptrAt(i)];
    bool isBDT2 = false;
    if(l.isEB() && fSCeta < 0.8 && BDT2 > -0.68)
        isBDT2 = true;// 0.374;
    if(l.isEB() && fSCeta >= 0.8 && BDT2 > -0.66)
        isBDT2 = true;// 0.374;
    if(l.isEE() && BDT2 > -0.6)
        isBDT2 = true;// 0.336;
    l.addUserFloat("BDT2", BDT2);
    l.addUserFloat("isBDT2", isBDT2);
    l.addUserFloat("pass_WP_2", pass_WP_2);
    //l.addUserFloat("pfSCfbrem", pfSCfbrem);
    //l.addUserFloat("is_tle", true);
    //edm::LogPrint("") << "isBDT: " << isBDT;
    //edm::LogPrint("") << "In Philler, reading userFloat: isBDT" << l.userFloat("isBDT") << " eta " << l.eta() << "pt " << l.pt();

    //--- Check selection cuts.
    if (!cut(l)) continue;

    //--- Embed flags (ie flags specified in the "flags" pset)
    for(CutSet<pat::Photon>::const_iterator flag = flags.begin(); flag != flags.end(); ++flag) {
      l.addUserFloat(flag->first,int((*(flag->second))(l)));
    }

    result->push_back(l);
  }
  iEvent.put(std::move(result));
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(Philler);
