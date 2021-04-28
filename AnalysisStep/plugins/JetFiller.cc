#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/Utilities/interface/EDMException.h>
#include "FWCore/ParameterSet/interface/FileInPath.h"

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

  std::string jecUncFile_;
  std::vector<string> uncSources {};
  std::vector<JetCorrectionUncertainty*> splittedUncerts_;
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

  if (setup == 2016)
    {
      edm::FileInPath jecUncFile("ZZAnalysis/AnalysisStep/data/JEC/Regrouped_Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.txt");
      jecUncFile_ = jecUncFile.fullPath();
      uncSources.push_back("Total");
      uncSources.push_back("Absolute"); uncSources.push_back("Absolute_2016");
      uncSources.push_back("BBEC1"); uncSources.push_back("BBEC1_2016");
      uncSources.push_back("EC2"); uncSources.push_back("EC2_2016");
      uncSources.push_back("FlavorQCD");
      uncSources.push_back("HF"); uncSources.push_back("HF_2016");
      uncSources.push_back("RelativeBal");
      uncSources.push_back("RelativeSample_2016");
    }
  else if (setup == 2017)
    {
      edm::FileInPath jecUncFile("ZZAnalysis/AnalysisStep/data/JEC/Regrouped_Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.txt");
      jecUncFile_ = jecUncFile.fullPath();
      uncSources.push_back("Total");
      uncSources.push_back("Absolute"); uncSources.push_back("Absolute_2017");
      uncSources.push_back("BBEC1"); uncSources.push_back("BBEC1_2017");
      uncSources.push_back("EC2"); uncSources.push_back("EC2_2017");
      uncSources.push_back("FlavorQCD");
      uncSources.push_back("HF"); uncSources.push_back("HF_2017");
      uncSources.push_back("RelativeBal");
      uncSources.push_back("RelativeSample_2017");
    }
  else if (setup == 2018)
    {
      edm::FileInPath jecUncFile("ZZAnalysis/AnalysisStep/data/JEC/Regrouped_Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt");
      jecUncFile_ = jecUncFile.fullPath();
      uncSources.push_back("Total");
      uncSources.push_back("Absolute"); uncSources.push_back("Absolute_2018");
      uncSources.push_back("BBEC1"); uncSources.push_back("BBEC1_2018");
      uncSources.push_back("EC2"); uncSources.push_back("EC2_2018");
      uncSources.push_back("FlavorQCD");
      uncSources.push_back("HF"); uncSources.push_back("HF_2018");
      uncSources.push_back("RelativeBal");
      uncSources.push_back("RelativeSample_2018");
    }
  else cout << "jecUncFile NOT FOUND!";
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


  // JEC uncertainty (part 1) - No splitting
  // JetCorrectorParametersCollection refers to the JEC file read from db in MasterPy/ZZ4lAnalysis.py
  ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get(jecType,JetCorParColl);
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty jecUnc(JetCorPar);

  // JEC uncertainty (Part 2) - Splitting (May 2020)
  // Run 2 reduced set of uncertainties from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECUncertaintySources#Run_2_reduced_set_of_uncertainty
  // List of uncertainties: ['Absolute', 'Absolute_201*', 'BBEC1', 'BBEC1_201*', 'EC2', 'EC2_201*', 'FlavorQCD', 'HF', 'HF_201*', 'RelativeBal', 'RelativeSample_201*'] + 'Total'
  if(applyJEC_ && isMC_)
    {
      for (unsigned s_unc = 0; s_unc < uncSources.size(); s_unc++)
	{
	  JetCorrectorParameters corrParams = JetCorrectorParameters(jecUncFile_, uncSources[s_unc]);
	  splittedUncerts_.push_back(new JetCorrectionUncertainty(corrParams));
	}
    }

  //--- Output collection
  auto result = std::make_unique<pat::JetCollection>();
  int jet_number = 0;
  for(auto jet = jetHandle->begin(); jet != jetHandle->end(); ++jet){
    jet_number++;

    pat::Jet j(*jet);
    double jpt = j.pt();
    double jeta = j.eta();
    double jabseta = fabs(jeta);
    double raw_jpt = j.correctedJet("Uncorrected").pt();

    // 20170220: using a float instead of a double changes of the JER seed from 99494 to 99495, and changes post JER jet pT.
    // Note that while PAT::Candidate has this function as double, we only save float accuracy in miniAOD anyway.

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
    float jes_unc = jecUnc.getUncertainty(true); //It takes as argument "bool fDirection": true = up, false = dn; symmetric values
    float pt_jesup = jpt * (1.0 + jes_unc); // set the shifted pT up
    float pt_jesdn = jpt * (1.0 - jes_unc); // set the shifted pT dn
    //j.setP4(j.p4() * (1 + jes_unc)); // Checked that only pt and energy/mass are affected by JEC variations considering the p4() components, not eta/phi

    vector<float> jes_unc_split {};
    vector<float> pt_jesup_split {};
    vector<float> pt_jesdn_split {};
    float singleContr_jes_unc = 0;

    if(applyJEC_ && isMC_)
      {
	for (unsigned s_unc = 0; s_unc < uncSources.size(); s_unc++)
	  {
	    singleContr_jes_unc = 0;
	    splittedUncerts_[s_unc]->setJetEta(jeta);
	    splittedUncerts_[s_unc]->setJetPt(jpt);
	    singleContr_jes_unc = splittedUncerts_[s_unc]->getUncertainty(true); //It takes as argument "bool fDirection": true = up, false = dn; symmetric values
	    jes_unc_split.push_back(singleContr_jes_unc);
	    pt_jesup_split.push_back( jpt * (1.0 + singleContr_jes_unc));
	    pt_jesdn_split.push_back( jpt * (1.0 - singleContr_jes_unc));
	  }
      }
    else
      {
	for(unsigned s_unc = 0; s_unc < uncSources.size(); s_unc++)
	  {
	    jes_unc_split.push_back(-999.);
	    pt_jesup_split.push_back(-999.);
	    pt_jesdn_split.push_back(-999.);
	  }
      }
    // For DEBUG: JEC splitting
  //   cout << "=============================" << endl;
  //   cout << "Previous total JES UNCERTAINTY = " << jes_unc << endl;
  //   cout << "----- NUMBER of CONSIDERED UNCERTAINTIES = " << uncSources.size() << " -----" << endl;
  //   cout << "JEC unc sources considered:\n";
  //   for (unsigned s_unc = 0; s_unc < uncSources.size(); s_unc++)
  //     {
	// cout << s_unc << " " << uncSources[s_unc] << '\t' << jes_unc_split[s_unc] << "," << endl;
  //     }
  //
  //   cout << "----- JET NUMERO " << jet_number << " -----" << endl;
  //   cout << "pT JesPt_UP = " << pt_jesup << endl;
  //   cout << "pT JesPt_DN = " << pt_jesdn << endl;
  //
  //   for (unsigned s_unc = 0; s_unc < uncSources.size(); s_unc++)
  //     {
	// cout << s_unc << " pT JesPt_UP SPLIT = " << pt_jesup_split[s_unc] << endl;
	// cout << s_unc << " pT JesPt_DN SPLIT = " << pt_jesdn_split[s_unc] << endl;
	// cout << "------------------------------------------------" << endl;
  //     }


    //--- Loose jet ID, cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13TeVRun2016
    float NHF  = j.neutralHadronEnergyFraction();
    float NEMF = j.neutralEmEnergyFraction();
    float CHF  = j.chargedHadronEnergyFraction();
    float CEMF = j.chargedEmEnergyFraction();
    int NumConst = j.chargedMultiplicity()+j.neutralMultiplicity();
    int NumNeutralParticles = j.neutralMultiplicity();
    float CHM  = j.chargedMultiplicity();
    //   float MUF  = j.muonEnergyFraction();

    bool JetID = true;
	  
    // if ( setup == 2016 )
    // { // Tight jet ID https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
    //   JetID    = ((NHF<0.90 && NEMF<0.90 && NumConst>1) && ((jabseta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || jabseta>2.4) && jabseta<=2.7) ||
    //              ( NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 && jabseta>2.7 && jabseta<=3.0 ) ||
    //              ( NEMF<0.90 && NumNeutralParticles>10 && jabseta >3.0 );
    // }

    // else if ( setup == 2017 )
    // {
    //  // Tight jet ID https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017 without JetIDLepVeto
    //   JetID    = ((NHF<0.90 && NEMF<0.90 && NumConst>1) && ((jabseta<=2.4 && CHF>0 && CHM>0) || jabseta>2.4) && jabseta<=2.7) ||
    //              ( NEMF<0.99 && NEMF>0.02 && NumNeutralParticles>2 && jabseta>2.7 && jabseta<=3.0 ) ||
    //              ( NEMF<0.90 && NHF>0.02 && NumNeutralParticles>10 && jabseta>3.0 );
    // }

    // else if ( setup == 2018) 
    // {
    //  // Tight jet ID https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13TeVRun2018 without JetIDLepVeto 
    //   JetID    = ( CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && jabseta<=2.6) ||
    //              ( CHM>0 && NEMF<0.99 && NHF < 0.9 && jabseta>2.6 && jabseta<=2.7) ||
    //              ( NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 && jabseta>2.7 && jabseta<=3.0 ) ||
    //              ( NEMF<0.90 && NHF>0.2 && NumNeutralParticles>10 && jabseta>3.0 );
    // }

    // else
    // {
    //  throw cms::Exception("JetID") << "Jet ID is not defined for the given setup (" << setup << ")!";
    // }
	 
    //--- JetID for UL, only one WP provided. Cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13TeVUL
    //--- Lepton veto not included, as we perform our own cleaning and selection between leps and jets
    if ( setup == 2016 )
    {
      JetID    =  (jabseta<=2.4 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF<0.9) || 
                  (jabseta>2.4 && jabseta<=2.7 && NEMF<0.99 && NHF < 0.9) ||
                  (jabseta>2.7 && jabseta<=3.0 && NEMF>0.0 && NEMF<0.99 && NHF<0.9 && NumNeutralParticle>1) ||
                  (jabseta>3.0 && NEMF<0.90 && NHF>0.2 && NumNeutralParticle>10);
    }

    else if ( (setup == 2017) || (setup == 2018) )
    {
      JetID    = (jabseta<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF<0.9) || 
                 (jabseta>2.6 && jabseta<=2.7 && CEMF<0.8 && CHM>0 && NEMF<0.99 && NHF < 0.9) ||
                 (jabseta>2.7 && jabseta<=3.0 && NEMF>0.01 && NEMF<0.99 && NumNeutralParticle>1) ||
                 (jabseta>3.0 && NEMF<0.90 && NHF>0.2 && NumNeutralParticle>10);
    }

    else
    {
     throw cms::Exception("JetID") << "Jet ID is not defined for the given setup (" << setup << ")!";
    } 

    bool PUjetID = false;
    float PUjetID_score = 0;

    //Recommended tight PU JET ID https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
    //Recommended tight WP for jet pileup ID taken from https://github.com/alefisico/cmssw/blob/PUID_102X/RecoJets/JetProducers/python/PileupJetIDCutParams_cfi.py
    // Computed using 2016 81X training and used for all three years
    //4 Eta Categories  0-2.5   2.5-2.75   2.75-3.0   3.0-5.0
    //Tight Id
    //Pt010_Tight    = cms.vdouble( 0.69, -0.35, -0.26, -0.21),
    //Pt1020_Tight   = cms.vdouble( 0.69, -0.35, -0.26, -0.21),
    //Pt2030_Tight   = cms.vdouble( 0.69, -0.35, -0.26, -0.21),
    //Pt3050_Tight   = cms.vdouble( 0.86, -0.10, -0.05, -0.01),

    if ( applyJEC_ && ( setup == 2017 || setup == 2018 || (setup == 2016 && (!(isMC_)) )))
    {
      PUjetID_score = j.userFloat("pileupJetIdUpdated:fullDiscriminant");
    }
    else
    {
      PUjetID_score = j.userFloat("pileupJetId:fullDiscriminant");
    }

    // WP applied explicitely
    //cout << "PU ID score = " << PUjetID_score << " PT = " << jpt << " ETA = " << jabseta << endl;
    // PUjetID       = ( 
		  //    (jpt > 0 && jpt <= 10 && jabseta > 0 && jabseta<=2.5) && PUjetID_score > 0.69  ||
		  //    (jpt > 10 && jpt <= 20 && jabseta > 0 && jabseta<=2.5) && PUjetID_score > 0.69 ||
		  //    (jpt > 20 && jpt <= 30 && jabseta > 0 && jabseta<=2.5) && PUjetID_score > 0.69 ||
		  //    (jpt > 30 && jpt <= 50 && jabseta > 0 && jabseta<=2.5) && PUjetID_score > 0.86 ||

		  //    (jpt > 0 && jpt <= 10 && jabseta > 2.5 && jabseta<=2.75) && PUjetID_score > -0.35  ||
		  //    (jpt > 10 && jpt <= 20 && jabseta > 2.5 && jabseta<=2.75) && PUjetID_score > -0.35 ||
		  //    (jpt > 20 && jpt <= 30 && jabseta > 2.5 && jabseta<=2.75) && PUjetID_score > -0.35 ||
		  //    (jpt > 30 && jpt <= 50 && jabseta > 2.5 && jabseta<=2.75) && PUjetID_score > -0.10 ||

		  //    (jpt > 0 && jpt <= 10 && jabseta > 2.75 && jabseta<=3.0) && PUjetID_score > -0.26  ||
		  //    (jpt > 10 && jpt <= 20 && jabseta > 2.75 && jabseta<=3.0) && PUjetID_score > -0.26 ||
		  //    (jpt > 20 && jpt <= 30 && jabseta > 2.75 && jabseta<=3.0) && PUjetID_score > -0.26 ||
		  //    (jpt > 30 && jpt <= 50 && jabseta > 2.75 && jabseta<=3.0) && PUjetID_score > -0.05 ||

		  //    (jpt > 0 && jpt <= 10 && jabseta > 3.0 && jabseta<=5.0) && PUjetID_score > -0.21  ||
		  //    (jpt > 10 && jpt <= 20 && jabseta > 3.0 && jabseta<=5.0) && PUjetID_score > -0.21 ||
		  //    (jpt > 20 && jpt <= 30 && jabseta > 3.0 && jabseta<=5.0) && PUjetID_score > -0.21 ||
		  //    (jpt > 30 && jpt <= 50 && jabseta > 3.0 && jabseta<=5.0) && PUjetID_score > -0.01 
		  //   );

    //--- UL Jet PU ID. Changed WPs and pT regions for the tranining. Cf. https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetIDUL#Recommendations_for_2018_UL_data
    //--- 4 Eta Categories: 0-2.5   2.5-2.75   2.75-3.0   3.0-5.0
    //--- Tight Jet PU ID:
    //--- Pt1020_Tight = cms.vdouble( 0.77, 0.38, -0.31, -0.21),
    //--- Pt2030_Tight = cms.vdouble( 0.90, 0.60, -0.12, -0.13),
    //--- Pt3040_Tight = cms.vdouble( 0.96, 0.82, 0.20, 0.09),
    //--- Pt4050_Tight = cms.vdouble( 0.98, 0.92, 0.47, 0.29),
    PUjetID = ( 
         (jpt > 10 && jpt <= 20 && jabseta > 0 && jabseta<=2.5) && PUjetID_score > 0.77 ||
         (jpt > 20 && jpt <= 30 && jabseta > 0 && jabseta<=2.5) && PUjetID_score > 0.90 ||
         (jpt > 30 && jpt <= 40 && jabseta > 0 && jabseta<=2.5) && PUjetID_score > 0.96 ||
         (jpt > 40 && jpt <= 50 && jabseta > 0 && jabseta<=2.5) && PUjetID_score > 0.98 ||

         (jpt > 10 && jpt <= 20 && jabseta > 2.5 && jabseta<=2.75) && PUjetID_score > 0.38 ||
         (jpt > 20 && jpt <= 30 && jabseta > 2.5 && jabseta<=2.75) && PUjetID_score > 0.60 ||
         (jpt > 30 && jpt <= 40 && jabseta > 2.5 && jabseta<=2.75) && PUjetID_score > 0.82 ||
         (jpt > 40 && jpt <= 50 && jabseta > 2.5 && jabseta<=2.75) && PUjetID_score > 0.92 ||

         (jpt > 10 && jpt <= 20 && jabseta > 2.75 && jabseta<=3.0) && PUjetID_score > -0.31 ||
         (jpt > 20 && jpt <= 30 && jabseta > 2.75 && jabseta<=3.0) && PUjetID_score > -0.12 ||
         (jpt > 30 && jpt <= 40 && jabseta > 2.75 && jabseta<=3.0) && PUjetID_score > 0.20  ||
         (jpt > 40 && jpt <= 50 && jabseta > 2.75 && jabseta<=3.0) && PUjetID_score > 0.47  ||

         (jpt > 10 && jpt <= 20 && jabseta > 3.0 && jabseta<=5.0) && PUjetID_score > -0.21 ||
         (jpt > 20 && jpt <= 30 && jabseta > 3.0 && jabseta<=5.0) && PUjetID_score > -0.13 ||
         (jpt > 30 && jpt <= 40 && jabseta > 3.0 && jabseta<=5.0) && PUjetID_score > 0.09  ||
         (jpt > 40 && jpt <= 50 && jabseta > 3.0 && jabseta<=5.0) && PUjetID_score > 0.29 
      );
    //if (PUjetID) cout << "PUjetID!" << endl;
    
    /*if ( applyJEC_ && ( setup == 2017 || setup == 2018 || (setup == 2016 && (!(isMC_)) )))
      {
        PUjetID_score = j.userFloat("pileupJetIdUpdated:fullDiscriminant"); 
        PUjetID = bool(j.userInt("pileupJetIdUpdated:fullId") & (1 << 0));
      }
    //if (applyJEC_) PUjetID = bool(j.userInt("pileupJetIdUpdated:fullId") & (1 << 0)); // MODIFIED according to update jet collection for 2016 MC (needed to include DeepCSV)                  
      else
      {
        PUjetID_score = j.userFloat("pileupJetId:fullDiscriminant"); 
        PUjetID = bool(j.userInt("pileupJetId:fullId") & (1 << 0));
      }
    */

//--- b tagging and scaling factors
    float bTagger;
    bTagger = j.bDiscriminator(bTaggerName) + j.bDiscriminator((bTaggerName + "b")); //one should just sum for doing b tagging, the b and bb probabilities
    //cout << "b tag is = " << bTagger << endl;

    // Check of tagger labels stored in the MiniAOD and recognized by the bDiscriminator
    //#const std::vector<std::pair<std::string, float> > & getPairDiscri() const;
    //auto& pd = j.getPairDiscri();
    //for (size_t pd_obj = 0; pd_obj < pd.size(); ++pd_obj)
    //  {
    //cout << pd_obj << "  Discriminator: " << pd.at(pd_obj).first << " \t " << pd.at(pd_obj).second  << endl;
    // }

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
    // Adding variables to take into account single contributions to the total JEC uncertainty:
    // UNCERTAINTY...
    j.addUserFloat("jes_unc_split_Total", jes_unc_split[0]);
    j.addUserFloat("jes_unc_split_Abs", jes_unc_split[1]);
    j.addUserFloat("jes_unc_split_Abs_year", jes_unc_split[2]);
    j.addUserFloat("jes_unc_split_BBEC1", jes_unc_split[3]);
    j.addUserFloat("jes_unc_split_BBEC1_year", jes_unc_split[4]);
    j.addUserFloat("jes_unc_split_EC2", jes_unc_split[5]);
    j.addUserFloat("jes_unc_split_EC2_year", jes_unc_split[6]);
    j.addUserFloat("jes_unc_split_FlavQCD", jes_unc_split[7]);
    j.addUserFloat("jes_unc_split_HF", jes_unc_split[8]);
    j.addUserFloat("jes_unc_split_HF_year", jes_unc_split[9]);
    j.addUserFloat("jes_unc_split_RelBal", jes_unc_split[10]);
    j.addUserFloat("jes_unc_split_RelSample_year", jes_unc_split[11]);
    // ... and pT UP/DN VARIATIONS
    j.addUserFloat("pt_jesup_split_Total", pt_jesup_split[0]);
    j.addUserFloat("pt_jesdn_split_Total", pt_jesdn_split[0]);
    j.addUserFloat("pt_jesup_split_Abs", pt_jesup_split[1]);
    j.addUserFloat("pt_jesdn_split_Abs", pt_jesdn_split[1]);
    j.addUserFloat("pt_jesup_split_Abs_year", pt_jesup_split[2]);
    j.addUserFloat("pt_jesdn_split_Abs_year", pt_jesdn_split[2]);
    j.addUserFloat("pt_jesup_split_BBEC1", pt_jesup_split[3]);
    j.addUserFloat("pt_jesdn_split_BBEC1", pt_jesdn_split[3]);
    j.addUserFloat("pt_jesup_split_BBEC1_year", pt_jesup_split[4]);
    j.addUserFloat("pt_jesdn_split_BBEC1_year", pt_jesdn_split[4]);
    j.addUserFloat("pt_jesup_split_EC2", pt_jesup_split[5]);
    j.addUserFloat("pt_jesdn_split_EC2", pt_jesdn_split[5]);
    j.addUserFloat("pt_jesup_split_EC2_year", pt_jesup_split[6]);
    j.addUserFloat("pt_jesdn_split_EC2_year", pt_jesdn_split[6]);
    j.addUserFloat("pt_jesup_split_FlavQCD", pt_jesup_split[7]);
    j.addUserFloat("pt_jesdn_split_FlavQCD", pt_jesdn_split[7]);
    j.addUserFloat("pt_jesup_split_HF", pt_jesup_split[8]);
    j.addUserFloat("pt_jesdn_split_HF", pt_jesdn_split[8]);
    j.addUserFloat("pt_jesup_split_HF_year", pt_jesup_split[9]);
    j.addUserFloat("pt_jesdn_split_HF_year", pt_jesdn_split[9]);
    j.addUserFloat("pt_jesup_split_RelBal", pt_jesup_split[10]);
    j.addUserFloat("pt_jesdn_split_RelBal", pt_jesdn_split[10]);
    j.addUserFloat("pt_jesup_split_RelSample_year", pt_jesup_split[11]);
    j.addUserFloat("pt_jesdn_split_RelSample_year", pt_jesdn_split[11]);
    //////////////////////////////////////////////////////////////////////////////////////////////
    j.addUserFloat("pt_jerup", pt_jerup);
    j.addUserFloat("pt_jerdn", pt_jerdn);
    j.addUserFloat("RawPt", raw_jpt);
    j.addUserFloat("JetID",JetID);
    j.addUserFloat("PUjetID",PUjetID);
    j.addUserFloat("PUjetID_score",PUjetID_score);
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


    result->push_back(j);}


  //--- Reorder jets by pT
  std::sort(result->begin(),result->end(), [](const Jet& j1, const Jet& j2){return j1.pt()>j2.pt();});

  iEvent.put(std::move(result));
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(JetFiller);
