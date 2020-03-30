#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/Utilities/interface/EDMException.h>

#include <DataFormats/PatCandidates/interface/Jet.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <JetMETCorrections/Objects/interface/JetCorrectionsRecord.h>
#include <JetMETCorrections/Modules/interface/JetResolution.h>

#include <DataFormats/Math/interface/deltaR.h>

#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>
#include <ZZAnalysis/AnalysisStep/interface/BTaggingSFHelper.h>

#include "TRandom3.h"

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;


class JetFiller : public edm::EDProducer {
 public:
  /// Constructor
  explicit JetFiller(const edm::ParameterSet&);
    
  /// Destructor
  ~JetFiller(){};  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  edm::EDGetTokenT<edm::View<pat::Jet> > jetToken;
  int sampleType;
  int setup;
  const StringCutObjectSelector<pat::Jet, true> cut;
  bool isMC_;
  const std::string bTaggerName;
  float bTaggerThreshold;
  bool applyJEC_;
  const std::string jecType;
  bool applyJER_;
  const std::string jerType;
  std::string bTagSFFile;
  std::string bTagMCEffFile;
  BTaggingSFHelper bTagSFHelper;
  const CutSet<pat::Jet> flags;
  edm::EDGetTokenT<double> rhoToken;
  edm::EDGetTokenT<edm::ValueMap<float> > qgToken;
  edm::EDGetTokenT<edm::ValueMap<float> > axis2Token;
  edm::EDGetTokenT<edm::ValueMap<int> > multToken;
  edm::EDGetTokenT<edm::ValueMap<float> > ptDToken;

  JME::JetResolution resolution;
  JME::JetResolutionScaleFactor resolution_sf;

};


JetFiller::JetFiller(const edm::ParameterSet& iConfig) :
  jetToken(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("src"))),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  cut(iConfig.getParameter<std::string>("cut")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  bTaggerName(iConfig.getParameter<std::string>("bTaggerName")),
  bTaggerThreshold(iConfig.getParameter<double>("bTaggerThreshold")),
  applyJEC_(iConfig.getParameter<bool>("applyJEC")),
  jecType(iConfig.getParameter<std::string>("jecType")),
  applyJER_(iConfig.getParameter<bool>("applyJER")),
  jerType(iConfig.getParameter<std::string>("jerType")),
  bTagSFFile(iConfig.getParameter<std::string>("bTagSFFile")),
  bTagMCEffFile(iConfig.getParameter<std::string>("bTagMCEffFile")),
  bTagSFHelper(bTagSFFile,bTagMCEffFile),
  flags(iConfig.getParameter<edm::ParameterSet>("flags"))
{

  rhoToken = consumes<double>(LeptonIsoHelper::getEleRhoTag(sampleType, setup));

  qgToken = consumes<edm::ValueMap<float> >(edm::InputTag("QGTagger", "qgLikelihood"));
  axis2Token = consumes<edm::ValueMap<float> >(edm::InputTag("QGTagger", "axis2"));
  multToken = consumes<edm::ValueMap<int> >(edm::InputTag("QGTagger", "mult"));
  ptDToken = consumes<edm::ValueMap<float> >(edm::InputTag("QGTagger", "ptD"));

  produces<pat::JetCollection>();
}


void
JetFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //--- Get jets
  Handle<edm::View<pat::Jet> > jetHandle;
  iEvent.getByToken(jetToken, jetHandle);

  //--- rho
  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken, rhoHandle);
  double rho = *rhoHandle;

  //--- q/g tagger
  Handle<edm::ValueMap<float> > qgHandle; 
  iEvent.getByToken(qgToken, qgHandle);
  Handle<edm::ValueMap<float> > axis2Handle; 
  iEvent.getByToken(axis2Token, axis2Handle);
  Handle<edm::ValueMap<int> > multHandle; 
  iEvent.getByToken(multToken, multHandle);
  Handle<edm::ValueMap<float> > ptDHandle; 
  iEvent.getByToken(ptDToken, ptDHandle);

  //--- JEC uncertanties 
  ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get(jecType,JetCorParColl);
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty jecUnc(JetCorPar);

  //--- Output collection
  auto result = std::make_unique<pat::JetCollection>();

  for(auto jet = jetHandle->begin(); jet != jetHandle->end(); ++jet){

    pat::Jet j(*jet);
    double jpt = j.pt();
    double jeta = j.eta();
    double jabseta = fabs(jeta);
    double raw_jpt = j.correctedJet("Uncorrected").pt();

    // 20170220: using a float instead of a double  changes of the JER seed from 99494 to 99495, and changes post JER jet pT.
    // Note that while PAT::Candidate has this fucntion as double, we only save float accuracy in miniAOD anyway.
    double jphi = j.phi();


    //--- Retrieve the q/g likelihood
    edm::RefToBase<pat::Jet> jetRef(edm::Ref<edm::View<pat::Jet> >(jetHandle, jet - jetHandle->begin()));
    float qgLikelihood = (*qgHandle)[jetRef];
    float axis2 = (*axis2Handle)[jetRef];
    int mult = (*multHandle)[jetRef];
    float ptD = (*ptDHandle)[jetRef];


    //--- Get JEC uncertainties 
    jecUnc.setJetEta(jeta);
    jecUnc.setJetPt(jpt);
    float jes_unc = jecUnc.getUncertainty(true);
     
    float pt_jesup = jpt * (1.0 + jes_unc);
    float pt_jesdn = jpt * (1.0 - jes_unc);


    //--- loose jet ID, cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13TeVRun2016 
    float NHF  = j.neutralHadronEnergyFraction();
    float NEMF = j.neutralEmEnergyFraction();
    float CHF  = j.chargedHadronEnergyFraction();
    float CEMF = j.chargedEmEnergyFraction();
    int NumConst = j.chargedMultiplicity()+j.neutralMultiplicity();
    int NumNeutralParticles = j.neutralMultiplicity();
    float CHM  = j.chargedMultiplicity();
    //   float MUF  = j.muonEnergyFraction();

    bool JetID = true;
	  
	 if ( setup == 2016 )
	 { // Tight jet ID https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
	 	JetID      = ((NHF<0.90 && NEMF<0.90 && NumConst>1) && ((jabseta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || jabseta>2.4) && jabseta<=2.7) ||
      				 ( NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 && jabseta>2.7 && jabseta<=3.0 ) ||
						 ( NEMF<0.90 && NumNeutralParticles>10 && jabseta >3.0 );
	 }
	 
	 else if ( setup == 2017 )
	 {
	   // Tight jet ID https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017 without JetIDLepVeto
	 	JetID      = ((NHF<0.90 && NEMF<0.90 && NumConst>1) && ((jabseta<=2.4 && CHF>0 && CHM>0) || jabseta>2.4) && jabseta<=2.7) ||
      				 ( NEMF<0.99 && NEMF>0.02 && NumNeutralParticles>2 && jabseta>2.7 && jabseta<=3.0 ) ||
      				 ( NEMF<0.90 && NHF>0.02 && NumNeutralParticles>10 && jabseta>3.0 );
	 }
	 
	 else if ( setup == 2018) 
	 {
	   // Tight jet ID https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13TeVRun2018 without JetIDLepVeto 
	 	JetID      =   ( CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && jabseta<=2.6) ||
                     ( CHM>0 && NEMF<0.99 && NHF < 0.9 && jabseta>2.6 && jabseta<=2.7) ||
                     ( NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 && jabseta>2.7 && jabseta<=3.0 ) ||
                     ( NEMF<0.90 && NHF>0.2 && NumNeutralParticles>10 && jabseta>3.0 );
	 }
	 
	 else
	 {
     throw cms::Exception("JetID") << "Jet ID is not defined for the given setup (" << setup << ")!";
	 }
	 

    bool PUjetID = true;
   
    //Recommended tight PU JET ID https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
    if (applyJEC_) PUjetID = bool(j.userInt("pileupJetIdUpdated:fullId") & (1 << 0));
    else PUjetID = bool(j.userInt("pileupJetId:fullId") & (1 << 0));

    //--- b tagging and scaling factors
    float bTagger;
    bTagger = j.bDiscriminator(bTaggerName) + j.bDiscriminator((bTaggerName + "b")); //one should just sum for doing b tagging, the b and bb probabilities

    bool isBtagged = bTagger > bTaggerThreshold;
    bool isBtaggedWithSF   = isBtagged;
    bool isBtaggedWithSFUp = isBtagged;
    bool isBtaggedWithSFDn = isBtagged;
    if(isMC_){
      int flav = j.hadronFlavour();
      TRandom3 rand;
      rand.SetSeed(abs(static_cast<int>(sin(jphi)*100000)));
      float R = rand.Uniform();
      float SF   = bTagSFHelper.getSF(central,flav,jpt,jeta);
      float SFUp = bTagSFHelper.getSF(up     ,flav,jpt,jeta);
      float SFDn = bTagSFHelper.getSF(down   ,flav,jpt,jeta);
      float bTagMCEff = bTagSFHelper.getEff(flav,jpt,jeta);
      if(SF  <=1 && isBtagged && R<1.-SF  ) isBtaggedWithSF   = false;
      if(SFUp<=1 && isBtagged && R<1.-SFUp) isBtaggedWithSFUp = false;
      if(SFDn<=1 && isBtagged && R<1.-SFDn) isBtaggedWithSFDn = false;
      if(SF  >1 && !isBtagged && R<(1.-SF  )/(1.-1./bTagMCEff)) isBtaggedWithSF   = true;
      if(SFUp>1 && !isBtagged && R<(1.-SFUp)/(1.-1./bTagMCEff)) isBtaggedWithSFUp = true;
      if(SFDn>1 && !isBtagged && R<(1.-SFDn)/(1.-1./bTagMCEff)) isBtaggedWithSFDn = true;
       //for debug
//      if(jpt>30 && isBtagged)
//      {
//			cout<<"jet pT="<<jpt<<", eta="<<jeta<<endl;
//			cout<<"bTag name="<<bTaggerName<<endl;
//			cout<<"bTag score="<<bTagger<<endl;
//			cout<<" flav = "<<flav<<endl;
//			cout<<" SF   = "<<SF  <<endl;
//			cout<<" SFUp = "<<SFUp<<endl;
//			cout<<" SFDn = "<<SFDn<<endl;
//			cout<<" bTagMCEff = "<<bTagMCEff<<endl;
//			cout<<" R = "<<R<<endl;
//			cout<<" isBtagged         = "<<isBtagged        <<endl;
//			cout<<" isBtaggedWithSF   = "<<isBtaggedWithSF  <<endl;
//			cout<<" isBtaggedWithSFUp = "<<isBtaggedWithSFUp<<endl;
//			cout<<" isBtaggedWithSFDn = "<<isBtaggedWithSFDn<<endl;
//      }
		 
    }


    //--- JER

    j.addUserFloat("pt_JEC_noJER", jpt);
    float pt_jer   = -1.;
    float pt_jerup = -1.;
    float pt_jerdn = -1.;

    if(isMC_ && applyJER_){
      resolution = JME::JetResolution::get(iSetup, jerType+"_pt");
      resolution_sf = JME::JetResolutionScaleFactor::get(iSetup, jerType);


      JME::JetParameters res_parameters = {{JME::Binning::JetPt, jpt}, {JME::Binning::JetEta, jeta}, {JME::Binning::Rho, rho}};
      float res_pt  = resolution.getResolution(res_parameters);

      //JME::JetParameters sf_parameters = {{JME::Binning::JetEta, jeta}, {JME::Binning::Rho, rho}};
      JME::JetParameters sf_parameters = {{JME::Binning::JetPt, jpt}, {JME::Binning::JetEta, jeta}, {JME::Binning::Rho, rho}};
      float sf    = resolution_sf.getScaleFactor(sf_parameters);
      float sf_up = resolution_sf.getScaleFactor(sf_parameters, Variation::UP);
      float sf_dn = resolution_sf.getScaleFactor(sf_parameters, Variation::DOWN);

      //-- 'hybrid' method = scaling for well matched jets, smearing for the rest

      //- check matching to gen jets
      const reco::GenJet* genJet = j.genJet();
      bool matchedJet = genJet 
	&& ( reco::deltaR(jeta,jphi,genJet->eta(),genJet->phi()) < 0.2 ) 
	&& ( fabs(jpt-genJet->pt()) < 3*res_pt*jpt );

      if(matchedJet){
	//- apply scaling
	float gen_pt = genJet->pt();
	pt_jer   = max( 0., gen_pt + sf   *(jpt-gen_pt) );
	pt_jerup = max( 0., gen_pt + sf_up*(jpt-gen_pt) );
	pt_jerdn = max( 0., gen_pt + sf_dn*(jpt-gen_pt) );
      }else{
	//- apply smearing
	TRandom3 rand;
	rand.SetSeed(abs(static_cast<int>(sin(jphi)*100000)));
	float smear = rand.Gaus(0,1.);
	float sigma   = sqrt(sf   *sf   -1.) * res_pt*jpt;
	float sigmaup = sqrt(sf_up*sf_up-1.) * res_pt*jpt;
	float sigmadn = sqrt(sf_dn*sf_dn-1.) * res_pt*jpt;
	pt_jer   = max( 0., smear*sigma   + jpt );
	pt_jerup = max( 0., smear*sigmaup + jpt );
	pt_jerdn = max( 0., smear*sigmadn + jpt );
      } 

      j.setP4(reco::Particle::PolarLorentzVector(pt_jer, jeta, jphi, (pt_jer/jpt)*j.mass()));
    }

    //cout<<"jet pT="<<jpt<<", eta="<<jeta<<endl;
    //--- Embed user variables
    j.addUserFloat("qgLikelihood",qgLikelihood);
    j.addUserFloat("axis2",axis2);
    j.addUserFloat("mult",mult);
    j.addUserFloat("ptD",ptD);
    j.addUserFloat("jes_unc", jes_unc);
    j.addUserFloat("pt_jesup", pt_jesup);
    j.addUserFloat("pt_jesdn", pt_jesdn);
    j.addUserFloat("pt_jerup", pt_jerup);
    j.addUserFloat("pt_jerdn", pt_jerdn);
    j.addUserFloat("RawPt", raw_jpt);
    j.addUserFloat("JetID",JetID);
    j.addUserFloat("PUjetID",PUjetID);
    j.addUserFloat("bTagger",bTagger);
    j.addUserFloat("isBtagged",isBtagged);
    j.addUserFloat("isBtaggedWithSF",isBtaggedWithSF);
    j.addUserFloat("isBtaggedWithSF_Up",isBtaggedWithSFUp);
    j.addUserFloat("isBtaggedWithSF_Dn",isBtaggedWithSFDn);


    //--- Apply selection cuts
    if (!cut(j)) continue;


    //--- Embed flags (ie flags specified in the "flags" pset)
    for(CutSet<pat::Jet>::const_iterator flag = flags.begin(); flag != flags.end(); ++flag) {
      j.addUserFloat(flag->first,int((*(flag->second))(j)));
    }

    
    result->push_back(j);
  }

  //--- Reorder jets by pT
  std::sort(result->begin(),result->end(), [](const Jet& j1, const Jet& j2){return j1.pt()>j2.pt();});

  iEvent.put(std::move(result));
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(JetFiller);

