/** \class ZZjjCandidateFiller
 *
 *
 *  \author R. Covarelli - Torino
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Jet.h>

#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>
#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>
#include <ZZMatrixElement/MELA/interface/Mela.h>
//#include <ZZAnalysis/AnalysisStep/interface/ZZMassErrors.h>
//#include <ZZAnalysis/AnalysisStep/interface/MCHistoryTools.h>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>
#include <ZZAnalysis/AnalysisStep/interface/CompositeCandMassResolution.h>
#include <FWCore/ParameterSet/interface/FileInPath.h>

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include <DataFormats/GeometryVector/interface/Point3DBase.h>
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include <DataFormats/PatCandidates/interface/Jet.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <JetMETCorrections/Objects/interface/JetCorrectionsRecord.h>
#include <JetMETCorrections/Modules/interface/JetResolution.h>

#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h>
#include <RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h>
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>

#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <ZZAnalysis/AnalysisStep/interface/Fisher.h>
#include <ZZAnalysis/AnalysisStep/interface/Comparators.h>
#include <ZZAnalysis/AnalysisStep/interface/utils.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>
#include <ZZAnalysis/AnalysisStep/interface/JetCleaner.h>
#include <KinZfitter/KinZfitter/interface/KinZfitter.h>

#include "TH2F.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include <string>

using namespace zzanalysis;

// bool doKinFitJJ = true;

class ZZjjCandidateFiller : public edm::EDProducer {
public:
  /// Constructor
  explicit ZZjjCandidateFiller(const edm::ParameterSet&);
    
  /// Destructor
  virtual ~ZZjjCandidateFiller();

private:
  typedef map<const reco::Candidate*, const pat::PFParticle*> FSRToLepMap;

  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  void getPairMass(const reco::Candidate* lp, const reco::Candidate* lm, FSRToLepMap& photons, float& mass, int& ID);

  edm::EDGetTokenT<edm::View<reco::CompositeCandidate> > candidateToken;
  const CutSet<pat::CompositeCandidate> preBestCandSelection;
  const CutSet<pat::CompositeCandidate> cuts;
  int sampleType;
  int setup;
  float superMelaMass;
  Mela* mela;
  bool embedDaughterFloats;
  bool isMerged;
  // bool ZRolesByMass;
  reco::CompositeCandidate::role_collection rolesZ1Z2;  
  reco::CompositeCandidate::role_collection rolesZ2Z1;
  bool isMC;
  bool recomputeIsoForFSR;
  TH2F* corrSigmaMu;
  TH2F* corrSigmaEle;
  Comparators::ComparatorTypes bestCandType;
  KinZfitter *kinZfitter;
  edm::EDGetTokenT<double> rhoForMuToken;
  edm::EDGetTokenT<double> rhoForEleToken;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetToken;
  edm::EDGetTokenT<vector<reco::MET> > METToken;
  edm::EDGetTokenT<edm::View<reco::Candidate> > softLeptonToken;
  edm::EDGetTokenT<edm::View<reco::CompositeCandidate> > ZCandToken;
  std::string candidateLabel;
};


static int SetupToSqrts(int setup) {
  if (setup==2011) return 7;
  else if (setup==2012) return 8;
  else if (setup==2015) return 13;
  else if (setup==2016) return 13;
  else return 0;
}


ZZjjCandidateFiller::ZZjjCandidateFiller(const edm::ParameterSet& iConfig) :
  candidateToken(consumes<edm::View<reco::CompositeCandidate> >(iConfig.getParameter<edm::InputTag>("src"))),
  preBestCandSelection(iConfig.getParameter<edm::ParameterSet>("bestCandAmong")),
  cuts(iConfig.getParameter<edm::ParameterSet>("flags")),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  superMelaMass(iConfig.getParameter<double>("superMelaMass")),
  embedDaughterFloats(iConfig.getUntrackedParameter<bool>("embedDaughterFloats",true)),
  // ZRolesByMass(iConfig.getParameter<bool>("ZRolesByMass")),
  isMerged(iConfig.getParameter<bool>("isMerged")),
  isMC(iConfig.getParameter<bool>("isMC")),
  recomputeIsoForFSR(iConfig.getParameter<bool>("recomputeIsoForFSR")),
  corrSigmaMu(0),
  corrSigmaEle(0),
  kinZfitter(0)
{
  mela = new Mela(SetupToSqrts(setup), superMelaMass, TVar::SILENT);
  mela->setCandidateDecayMode(TVar::CandidateDecay_ZZ);

  produces<pat::CompositeCandidateCollection>();
  rhoForMuToken = consumes<double>(LeptonIsoHelper::getMuRhoTag(sampleType, setup));
  rhoForEleToken = consumes<double>(LeptonIsoHelper::getEleRhoTag(sampleType, setup));
  jetToken = consumes<edm::View<pat::Jet> >(edm::InputTag("cleanJets"));
  METToken = consumes<vector<reco::MET> >(edm::InputTag("slimmedMETs"));
  softLeptonToken = consumes<edm::View<reco::Candidate> >(edm::InputTag("softLeptons"));
  ZCandToken = consumes<edm::View<reco::CompositeCandidate> >(edm::InputTag("ZCand"));
  
  rolesZ1Z2 = {"Z1", "Z2"};
  rolesZ2Z1 = {"Z2", "Z1"};
  
  edm::FileInPath fip("ZZAnalysis/AnalysisStep/data/ebeOverallCorrections.Legacy2013.v0.root");
  std::string ebePath=fip.fullPath();

  if (setup != 2015) {// FIXME:  EbE corrections to be updated for Run II
    // EbE corrections
    TFile* fCorrSigma = new TFile(ebePath.data()); // FIXME: is leaked
    std::string sigmaCorrType = (isMC?"mc":"reco");
    std::string sigmaCorrYear = ""; 
    if (setup==2011) sigmaCorrYear = "42x";
    else if (setup==2012) sigmaCorrYear = "53x";

    corrSigmaMu=  (TH2F*)fCorrSigma->Get(("mu_"+sigmaCorrType+sigmaCorrYear).data()); 
    corrSigmaEle= (TH2F*)fCorrSigma->Get(("el_"+sigmaCorrType+sigmaCorrYear).data());
  }
  
  string cmp=iConfig.getParameter<string>("bestCandComparator");
  if      (cmp=="byBestZ1bestZ2") bestCandType=Comparators::byBestZ1bestZ2;
  else if (cmp=="byBestZqq")      bestCandType=Comparators::byBestZqq;
  else if (cmp=="byBestKD")       bestCandType=Comparators::byBestKD;
  else if (cmp=="byBestKD_VH")    bestCandType=Comparators::byBestKD_VH;
  else if (cmp=="byBestPsig")    bestCandType=Comparators::byBestPsig;
  else if (cmp=="byMHWindow")    bestCandType=Comparators::byMHWindow;
  else abort();

  //-- kinematic refitter
  kinZfitter = new KinZfitter(!isMC);

  candidateLabel = iConfig.getParameter<edm::InputTag>("src").label();
}

ZZjjCandidateFiller::~ZZjjCandidateFiller(){
  delete kinZfitter;
  delete mela;
}

void ZZjjCandidateFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){  
  using namespace edm;
  using namespace std;
  using namespace reco;

  std::auto_ptr<pat::CompositeCandidateCollection> result(new pat::CompositeCandidateCollection);

  // const float ZmassValue = 91.1876;

  double rhoForMu, rhoForEle;
  {
    edm::Handle<double> rhoHandle;
    iEvent.getByToken(rhoForMuToken, rhoHandle);
    rhoForMu = *rhoHandle;
    iEvent.getByToken(rhoForEleToken, rhoHandle);
    rhoForEle = *rhoHandle;
  }

  //--- JEC uncertanties for fat jets
  ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get("AK8PFchs",JetCorParColl);
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty jecUnc(JetCorPar);
  
  JME::JetResolution resolution_pt, resolution_phi;

  resolution_pt = JME::JetResolution::get(iSetup, "AK4PFchs_pt"); 
  resolution_phi = JME::JetResolution::get(iSetup, "AK4PFchs_phi");  

  // Get LLLL candidates (both resolved and merged jets
  Handle<edm::View<CompositeCandidate> > LLLLCands;
  iEvent.getByToken(candidateToken, LLLLCands);

  // Get jets
  Handle<edm::View<pat::Jet> > CleanJets;
  iEvent.getByToken(jetToken, CleanJets);

  // Get MET
  Handle<vector<reco::MET> > pfmetcoll;
  iEvent.getByToken(METToken, pfmetcoll);
  math::XYZTLorentzVector pfmet;
  if(pfmetcoll.isValid()) pfmet = pfmetcoll->front().p4(); // standard MET is pfmet.pt();

  // Get leptons (in order to store extra leptons)
  Handle<View<reco::Candidate> > softleptoncoll;
  iEvent.getByToken(softLeptonToken, softleptoncoll);
  vector<reco::CandidatePtr> goodisoleptonPtrs;
  for( View<reco::Candidate>::const_iterator lep = softleptoncoll->begin(); lep != softleptoncoll->end(); ++ lep ){ 
    if((bool)userdatahelpers::getUserFloat(&*lep,"isGood") 
       //       && (bool)userdatahelpers::getUserFloat(&*lep,"isIsoFSRUncorr") // with old FSR strategy
       && (bool)userdatahelpers::getUserFloat(&*lep,"passCombRelIsoPFFSRCorr") // with new FSR strategy
       ){
      const reco::CandidatePtr lepPtr(softleptoncoll,lep-softleptoncoll->begin());
      goodisoleptonPtrs.push_back(lepPtr);
    }
  }

  // Get Z Candidates
  // Handle<View<CompositeCandidate> > ZCands;
  // iEvent.getByToken(ZCandToken, ZCands);

  // Get processID
//   edm::Handle<GenEventInfoProduct> gen;
//   iEvent.getByLabel( "generator", gen );
//   int processID = gen->signalProcessID();

  // to calculate mass resolution
  // CompositeCandMassResolution errorBuilder;       
  // errorBuilder.init(iSetup);

  vector<int> bestCandIdx(preBestCandSelection.size(),-1); 
  vector<float> maxPtSum(preBestCandSelection.size(),-1); 
  vector< vector<int> > preSelCands(preBestCandSelection.size());  

  //----------------------------------------------------------------------
  //--- Loop over input candidates
  // if (LLLLCands->size()>1) cout << "Size before " << candidateLabel << " = " << LLLLCands->size() << endl;
  for (View<CompositeCandidate>::const_iterator cand = LLLLCands->begin(); cand != LLLLCands->end(); ++cand) {

    int icand = distance(LLLLCands->begin(), cand);

    pat::CompositeCandidate myCand(*cand);

    if (embedDaughterFloats){
      userdatahelpers::embedDaughterData(myCand);
    }

    //--- Set id of the Z1 and "Z1"/"Z2" labels. This allows to call e.g. aHiggs->daughter("Z1"). 
    // if ZRolesByMass is true, 'Z1' and iZ1 refer to the  Z closest to mZ; otherwise Z1 = daughter(0). The latter is used for control regions.
    const reco::CompositeCandidate::role_collection* ZRoles = &rolesZ1Z2;
    // int iZ1 = 0;
    // int iZ2 = 1;
    /* if (ZRolesByMass) {
      if(std::abs(myCand.daughter(0)->mass()-ZmassValue)>=std::abs(myCand.daughter(1)->mass()-ZmassValue)){
      swap(iZ1,iZ2);
      ZRoles = &rolesZ2Z1;
      }
      }*/

    // ZRoles = &rolesZ2Z1;
    myCand.setRoles(*ZRoles);
    myCand.applyRoles();

    //--- Z pointers
    const reco::Candidate* Z1= myCand.daughter(0);
    const reco::Candidate* Z2= myCand.daughter(1);
    vector<const reco::Candidate*> Zs ={ Z1, Z2 }; // in the original order

    //--- Lepton pointers in the original order
    const reco::Candidate* Z1J1 = Z1->daughter(0);
    const reco::Candidate* Z1J2 = Z1->daughter(1);
    const reco::Candidate* Z2L1 = Z2->daughter(0);
    const reco::Candidate* Z2L2 = Z2->daughter(1);
    vector<const reco::Candidate*> ZZLeps ={ Z1J1, Z1J2, Z2L1, Z2L2 }; // array, in the original order

    // Create corresponding array of fourmomenta; will add FSR (below)
    vector<math::XYZTLorentzVector> pij(4);
    std::transform(ZZLeps.begin(), ZZLeps.end(), pij.begin(), [](const reco::Candidate* c){return c->p4(); });

    //--- Collect FSR photons and map them to the corresponding leptons
    FSRToLepMap FSRMap;
    // for (unsigned iZ=0; iZ<2; ++iZ) {
    for (unsigned ifsr=2; ifsr<Z2->numberOfDaughters(); ++ifsr) {
      const pat::PFParticle* fsr = static_cast<const pat::PFParticle*>(Zs[1]->daughter(ifsr));
      int ilep = 2+fsr->userFloat("leptIdx");
      FSRMap[ZZLeps[ilep]]= fsr;
      pij[ilep]+=fsr->p4();
    }
    // }

    //--- Lepton four-vectors in the original order; with FSR added
    math::XYZTLorentzVector p11 = pij[0];
    math::XYZTLorentzVector p12 = pij[1];
    math::XYZTLorentzVector p21 = pij[2];
    math::XYZTLorentzVector p22 = pij[3];
    // int id11 = Z1J1->pdgId();
    // int id12 = Z1J2->pdgId();
    int id21 = Z2L1->pdgId();
    int id22 = Z2L2->pdgId();
    int candChannel = id21*id22;
    float rho = 0.;

    // Recompute isolation for all four leptons if FSR is present
    //  for (int zIdx=0; zIdx<2; ++zIdx) {
    float worstMuIso=0;
    float worstEleIso=0;
    //cout << "LLLLCands ==== " << endl;
    for (int dauIdx=0; dauIdx<2; ++dauIdx) {
      const reco::Candidate* z = myCand.daughter(1);
      const reco::Candidate* d = z->daughter(dauIdx);
      // cout << "LLLLCands " << d->isMuon() << endl;
      float combRelIsoPFCorr = 0;
      rho = ((d->isMuon()) ? rhoForMu : rhoForEle);
      if (recomputeIsoForFSR) {  //FIXME: will recompute iso for individual leptons in the new scheme
        float fsrCorr = 0; // The correction to PFPhotonIso
        for (FSRToLepMap::const_iterator ifsr=FSRMap.begin(); ifsr!=FSRMap.end(); ++ifsr) {
          double dR = ROOT::Math::VectorUtil::DeltaR(ifsr->second->p4(), d->momentum());
          // Check if the photon is in the lepton's iso cone and not vetoed
          if (dR<0.3 && ((d->isMuon() && dR > 0.01) ||
            (d->isElectron() && (fabs((static_cast<const pat::Electron*>(d->masterClone().get()))->superCluster()->eta()) < 1.479 || dR > 0.08)))) {
            fsrCorr += ifsr->second->pt();
          }
        }
        combRelIsoPFCorr =  LeptonIsoHelper::combRelIsoPF(sampleType, setup, rho, d, fsrCorr);
        string base;
        stringstream str;
        str << "d1." << "d" << dauIdx << ".";
        str >> base;
        myCand.addUserFloat(base+"combRelIsoPFFSRCorr", combRelIsoPFCorr);
        myCand.addUserFloat(base+"passCombRelIsoPFFSRCorr", combRelIsoPFCorr < LeptonIsoHelper::isoCut(d)); // FIXME: not the most elegant solution; hard coded right now to see how things evolve about lepton isolation requirements. 
      }
      else {
        combRelIsoPFCorr = userdatahelpers::getUserFloat(d, "combRelIsoPFFSRCorr");
      }
      if (d->isMuon()) worstMuIso  = max(worstMuIso, combRelIsoPFCorr);
      else             worstEleIso = max(worstEleIso, combRelIsoPFCorr);
    }
    string base = "d1.";
    myCand.addUserFloat(base+"worstMuIso", worstMuIso);
    myCand.addUserFloat(base+"worstEleIso", worstEleIso);
    // }

    //--- Sign-ordered leptons and leptopn four-vectors (without FSR), to be used to compute mZa, mZb, mZalpha, mZbeta
    const reco::Candidate* Z1Lp(Z1J1);
    const reco::Candidate* Z1Lm(Z1J2);
    const reco::Candidate* Z2Lp(Z2L1);
    const reco::Candidate* Z2Lm(Z2L2);

    // Sort leptons for OS Z candidates; no sorting for the same-sign collections used for CRs
    // if (Z1Lp->charge() < 0 && Z1Lp->charge()*Z1Lm->charge()<0) {
    //  swap(Z1Lp,Z1Lm);    
    // }
    if (Z2Lp->charge() < 0 && Z2Lp->charge()*Z2Lm->charge()<0) {
      swap(Z2Lp, Z2Lm);
    }

    math::XYZTLorentzVector p1p(Z1Lp->p4());
    math::XYZTLorentzVector p1m(Z1Lm->p4());
    math::XYZTLorentzVector p2p(Z2Lp->p4());
    math::XYZTLorentzVector p2m(Z2Lm->p4());

    bool passSmartMLL = true;
    /* if (((ZaID==-121||ZaID==-169) && std::abs(mZa-ZmassValue)<std::abs(mZ1-ZmassValue) && mZb<12) ||
       ((ZalphaID==-121||ZalphaID==-169) && std::abs(mZalpha-ZmassValue)<std::abs(mZ1-ZmassValue) && mZbeta<12)) passSmartMLL = false;  */


    //--- QCD suppression cut
    vector<const reco::Candidate*> lep;
    lep.push_back(Z1Lm);
    lep.push_back(Z1Lp);
    lep.push_back(Z2Lm);
    lep.push_back(Z2Lp);

    //--- worst SIP value
    vector<double> SIPS ={ myCand.userFloat("d1.d0.SIP"), myCand.userFloat("d1.d1.SIP") };
    sort(SIPS.begin(), SIPS.end());
    float SIP4 = SIPS[1];

    //--- Sorted pTs
    vector<double> ptS;
    //ptS.push_back(Z1Lm->pt());
    //ptS.push_back(Z1Lp->pt());
    ptS.push_back(Z2Lm->pt());
    ptS.push_back(Z2Lp->pt());
    sort(ptS.begin(), ptS.end());
   
    //----------------------------------------------------------------------
    //--- Embed angular information and probabilities to build discriminants
    float costheta1=0, costheta2=0, phi=0, costhetastar=0, phistar1=0;
    float xi=0, xistar=0;

    // detaJJ, Mjj and Fisher. These are per-event variables in the SR, but not necessarily in the CR as we clean jets also 
    // for loose-but-not-tight leptons.
    
    float DiJetMass  = -99;
    float DiJetDEta  = -99;
    float DiJetFisher  = -99;

    unsigned int nCandidates=0; // Should equal 3 after the loop below
    for (int jecnum = 0; jecnum < 3; jecnum++){

      // Lepton TLorentzVectors, including FSR 
      SimpleParticleCollection_t daughters;
      
       // multiplier: +1 means JEC up, -1 means JEC down
      double jecnum_multiplier = 0;
      if (jecnum==1) jecnum_multiplier = 1.;
      else if (jecnum==2) jecnum_multiplier = -1.;

       // JEC corrected masses
      float newZ1Mass = 0.;
      float newZZMass = 0.;
      float jratio=0., jratio1=0., jratio2=0.;
      const reco::Candidate* z = myCand.daughter(0);
      if (isMerged) {
	const pat::Jet* thejet = dynamic_cast <const pat::Jet*> (z->masterClone().get());
	jecUnc.setJetEta(thejet->eta());
	jecUnc.setJetPt(thejet->pt());
	float jec_unc = jecUnc.getUncertainty(true);
        jratio = 1. + jecnum_multiplier * jec_unc;
	jecUnc.setJetEta(Z1J1->eta());
	jecUnc.setJetPt(Z1J1->pt());
	float jec_unc1 = jecUnc.getUncertainty(true);
	jratio1 = 1. + jecnum_multiplier * jec_unc1;
	jecUnc.setJetEta(Z1J2->eta());
	jecUnc.setJetPt(Z1J2->pt());
	float jec_unc2 = jecUnc.getUncertainty(true);
        jratio2 = 1. + jecnum_multiplier * jec_unc2;
	newZ1Mass = jratio*myCand.userFloat("d0.ak8PFJetsCHSCorrPrunedMass");
	newZZMass = (jratio*z->p4() + myCand.daughter(1)->p4()).mass();
	
      } else {	  
	const reco::Candidate* d1 = z->daughter(0);
	const pat::Jet* thejet1 = dynamic_cast <const pat::Jet*> (d1->masterClone().get());
	float jec_unc1 = thejet1->userFloat("jec_unc");
	jratio1 = 1. + jecnum_multiplier * jec_unc1;
	const reco::Candidate* d2 = z->daughter(1);
	const pat::Jet* thejet2 = dynamic_cast <const pat::Jet*> (d2->masterClone().get());
	float jec_unc2 = thejet2->userFloat("jec_unc");
	jratio2 = 1. + jecnum_multiplier * jec_unc2;
	newZ1Mass = (jratio1*d1->p4() + jratio2*d2->p4()).mass();
	newZZMass = (jratio1*d1->p4() + jratio2*d2->p4() + myCand.daughter(1)->p4()).mass();
      }
      if (jecnum == 1) {
	myCand.addUserFloat("Z1Mass_JecUp", newZ1Mass);
	myCand.addUserFloat("ZZMass_JecUp", newZZMass);
      } else {
	myCand.addUserFloat("Z1Mass_JecDown", newZ1Mass);
	myCand.addUserFloat("ZZMass_JecDown", newZZMass);
      } 

      /* if (isMerged) {
	daughters.push_back(SimpleParticle_t(0, TLorentzVector(p11.x()*jratio, p11.y()*jratio, p11.z()*jratio, p11.t()*jratio)));
	daughters.push_back(SimpleParticle_t(0, TLorentzVector(p12.x()*jratio, p12.y()*jratio, p12.z()*jratio, p12.t()*jratio)));
	daughters.push_back(SimpleParticle_t(id21, TLorentzVector(p21.x(), p21.y(), p21.z(), p21.t())));
	daughters.push_back(SimpleParticle_t(id22, TLorentzVector(p22.x(), p22.y(), p22.z(), p22.t())));
	} else { */
      daughters.push_back(SimpleParticle_t(0, TLorentzVector(p11.x()*jratio1, p11.y()*jratio1, p11.z()*jratio1, p11.t()*jratio1)));
      daughters.push_back(SimpleParticle_t(0, TLorentzVector(p12.x()*jratio2, p12.y()*jratio2, p12.z()*jratio2, p12.t()*jratio2)));
      daughters.push_back(SimpleParticle_t(id21, TLorentzVector(p21.x(), p21.y(), p21.z(), p21.t())));
      daughters.push_back(SimpleParticle_t(id22, TLorentzVector(p22.x(), p22.y(), p22.z(), p22.t())));
      // } 
      
      //--- Compute angles, better done here
      if (jecnum == 0) {
	TUtil::computeAngles(
			     daughters.at(0).second, daughters.at(0).first,
			     daughters.at(1).second, daughters.at(1).first,
			     daughters.at(2).second, daughters.at(2).first,
			     daughters.at(3).second, daughters.at(3).first,
			     costhetastar, costheta1, costheta2, phi, phistar1
			     );
	//--- compute higgs azimuthal angles, xi
	TLorentzVector Z14vec = daughters.at(0).second + daughters.at(1).second;
	TLorentzVector higgs = Z14vec + daughters.at(2).second + daughters.at(3).second;
	TVector3 Xaxis(1, 0, 0);
	xi = higgs.Phi();
	// boost Z1 into rest frame of higgs
	// xistar is the angle between Z decay plane and x-axis 
	Z14vec.Boost(-higgs.BoostVector());
 	xistar = Z14vec.Phi();
      }
      
      SimpleParticleCollection_t associated;
      
      vector<const pat::Jet*> cleanedJetsPt30Jec;
      vector<float> jec_ratio;
    
      // MELA variables
      for (edm::View<pat::Jet>::const_iterator jet = CleanJets->begin(); jet != CleanJets->end(); ++jet){
        // calculate JEC uncertainty up/down
        float jec_unc = jet->userFloat("jec_unc");
        float ratio = 1. + jecnum_multiplier * jec_unc;
        float newPt = jet->pt() * ratio;
        // apply pt>30GeV cut
        if (newPt<=30.0) continue;
        // additional jets cleaning 
        if (!jetCleaner::isGood(myCand, *jet)) continue;
        // remove jets belonging to the candidate
        bool belongs = false;
        if (isMerged) {
          for (int dauIdx=0; dauIdx<2; ++dauIdx) {       
            const reco::Candidate* z = myCand.daughter(1);
            const reco::Candidate* d = z->daughter(dauIdx);
            double dR = ROOT::Math::VectorUtil::DeltaR(jet->p4(), d->momentum());
            if (dR < 0.4) belongs = true;
          }
          const reco::Candidate* z2 = myCand.daughter(0);
          double dR2 = ROOT::Math::VectorUtil::DeltaR(jet->p4(), z2->momentum());
          // cout << "signal merged jet pt = " << z2->pt() << " eta = " << z2->eta() << std::endl;
          if (dR2 < 0.8) belongs = true;
        } else {
          for (int theCand=0; theCand<2; ++theCand) {
            for (int dauIdx=0; dauIdx<2; ++dauIdx) {
              const reco::Candidate* z = myCand.daughter(theCand);
              const reco::Candidate* d = z->daughter(dauIdx);
              double dR = ROOT::Math::VectorUtil::DeltaR(jet->p4(), d->momentum());
              // if (theCand == 0) cout << "resolved merged jet pt = " << d->pt() << " eta = " << d->eta() << std::endl;
              if (dR < 0.4) belongs = true;
            }
          }
        }
        if (belongs) continue;
        // store jets and up/down ratio
        // cout << "other jet pt = " << jet->pt() << " eta = " << jet->eta() << std::endl;
        cleanedJetsPt30Jec.push_back(&*jet);
        jec_ratio.push_back(ratio);
      }
      if (jecnum==0 && cleanedJetsPt30Jec.size()>1){
        const pat::Jet& jet1 = *(cleanedJetsPt30Jec.at(0));
        const pat::Jet& jet2 = *(cleanedJetsPt30Jec.at(1));
        DiJetDEta = jet1.eta()-jet2.eta();
        DiJetMass = (jet1.p4()+jet2.p4()).M();
        DiJetFisher = fisher(DiJetMass, DiJetDEta);
      }
      for (unsigned int ijet = 0; ijet<cleanedJetsPt30Jec.size(); ijet++){
        TLorentzVector jet(
          cleanedJetsPt30Jec[ijet]->p4().x()*jec_ratio.at(ijet),
          cleanedJetsPt30Jec[ijet]->p4().y()*jec_ratio.at(ijet),
          cleanedJetsPt30Jec[ijet]->p4().z()*jec_ratio.at(ijet),
          cleanedJetsPt30Jec[ijet]->p4().t()*jec_ratio.at(ijet)
          );
        associated.push_back(SimpleParticle_t(0, jet));
      }
      mela->setInputEvent(&daughters, &associated, 0, 0); nCandidates++;
    }

    mela->setCurrentCandidateFromIndex(0);

    // float p0plus_VAJHU=0;
    // mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
    // mela->computeP(p0plus_VAJHU, true);
    float p0minus_VAJHU=0;
    mela->setProcess(TVar::H0minus, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p0minus_VAJHU, true);
    float p0hplus_VAJHU=0;
    mela->setProcess(TVar::H0hplus, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p0hplus_VAJHU, true);
    // float p2bplus_VAJHU=0;
    // mela->setProcess(TVar::H2_g5, TVar::JHUGen, TVar::ZZGG);
    // mela->computeP(p2bplus_VAJHU, true);
    float p2_VAJHU=0;
    mela->setProcess(TVar::H2_g1g5, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p2_VAJHU, true);

    // float pqqZJJ_VAMCFM=0;
    // mela->setProcess(TVar::bkgZJets, TVar::MCFM, TVar::JJQCD);
    // mela->computeP(pqqZJJ_VAMCFM, true);
    float bkg_VAMCFM=0;
    mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
    mela->computeP(bkg_VAMCFM, true);
    float ggzz_VAMCFM=0;
    mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
    mela->computeP(ggzz_VAMCFM, true);
    // float p0plus_VAMCFM=0;
    // mela->setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    // mela->computeP(p0plus_VAMCFM, true);
    float ggzz_p0plus_VAMCFM=0;
    mela->setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela->computeP(ggzz_p0plus_VAMCFM, true);
    float Dgg10_VAMCFM=0;
    mela->computeD_gg(TVar::MCFM, TVar::D_gg10, Dgg10_VAMCFM);

    float pvbf_VAJHU_highestPTJets=-1;
    float phjj_VAJHU_highestPTJets=-1;
    float pvbf_VAJHU_highestPTJets_up=-1;
    float phjj_VAJHU_highestPTJets_up=-1;
    float pvbf_VAJHU_highestPTJets_dn=-1;
    float phjj_VAJHU_highestPTJets_dn=-1;
    float pqqZJJ_VAMCFM=-1;
    float p0plus_VAJHU=-1;
    float pqqZJJ_VAMCFM_up=-1;
    float p0plus_VAJHU_up=-1;
    float pqqZJJ_VAMCFM_dn=-1;
    float p0plus_VAJHU_dn=-1;
    float p0plus_VAMCFM=-1;
    float p2bplus_VAJHU=-1;
    float p0plus_VAMCFM_up=-1;
    float p2bplus_VAJHU_up=-1;
    float p0plus_VAMCFM_dn=-1;
    float p2bplus_VAJHU_dn=-1;
    
    // Do these loops at the end to avoid switching particles off first and then on again
    for (unsigned int jecnum=0; jecnum<nCandidates; jecnum++){
      mela->setCurrentCandidateFromIndex(jecnum);
      MELACandidate* melaCand = mela->getCurrentCandidate();

      if (melaCand!=0){
        unsigned int nGoodJets=melaCand->getNAssociatedJets();
        bool hasAtLeastOneJet = (nGoodJets>0);
        //bool hasAtLeastTwoJets = (nGoodJets>1);
	float pqqZJJ_VAMCFM_temp=-1;
	float p0plus_VAJHU_temp=-1;
	float p0plus_VAMCFM_temp=-1;
	float p2bplus_VAJHU_temp=-1;
        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
	mela->computeP(p0plus_VAJHU_temp, true);
        mela->setProcess(TVar::H2_g5, TVar::JHUGen, TVar::ZZGG);
	mela->computeP(p2bplus_VAJHU_temp, true);
        mela->setProcess(TVar::bkgZJets, TVar::MCFM, TVar::JJQCD);
	mela->computeP(pqqZJJ_VAMCFM_temp, true);
        mela->setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
	mela->computeP(p0plus_VAMCFM_temp, true);

        if (hasAtLeastOneJet){
          float phjj_VAJHU_highestPTJets_temp = -1;
          float pvbf_VAJHU_highestPTJets_temp = -1;
      
          for (unsigned int firstjet = 0; firstjet < nGoodJets; firstjet++){
            for (unsigned int secondjet = 1; secondjet < nGoodJets; secondjet++){
              if (secondjet<=firstjet) continue;
              for (unsigned int disableJet=0; disableJet<nGoodJets; disableJet++){
                bool flag=false;
                if (disableJet==firstjet || disableJet==secondjet) flag=true;
                melaCand->getAssociatedJet(disableJet)->setSelected(flag);
              }

              float pvbf_temp = -1;
              mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
              mela->computeProdP(pvbf_temp, true);
              float phjj_temp = -1;
              mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJQCD);
              mela->computeProdP(phjj_temp, true);
            
              if (firstjet == 0 && secondjet == 1){
                phjj_VAJHU_highestPTJets_temp = phjj_temp;
                pvbf_VAJHU_highestPTJets_temp = pvbf_temp;
              
              }
            }
          
          }

          for (unsigned int disableJet=0; disableJet<nGoodJets; disableJet++) melaCand->getAssociatedJet(disableJet)->setSelected(true); // Turn everything back on

          if (jecnum == 0){
            phjj_VAJHU_highestPTJets = phjj_VAJHU_highestPTJets_temp;
            pvbf_VAJHU_highestPTJets = pvbf_VAJHU_highestPTJets_temp;
	    p0plus_VAJHU = p0plus_VAJHU_temp;
	    p0plus_VAMCFM = p0plus_VAMCFM_temp;
	    p2bplus_VAJHU = p2bplus_VAJHU_temp;
	    pqqZJJ_VAMCFM = pqqZJJ_VAMCFM_temp;
          }
          else if (jecnum == 1){
            phjj_VAJHU_highestPTJets_up = phjj_VAJHU_highestPTJets_temp;
            pvbf_VAJHU_highestPTJets_up = pvbf_VAJHU_highestPTJets_temp;
	    p0plus_VAJHU_up = p0plus_VAJHU_temp;
	    p0plus_VAMCFM_up = p0plus_VAMCFM_temp;
	    p2bplus_VAJHU_up = p2bplus_VAJHU_temp;
	    pqqZJJ_VAMCFM_up = pqqZJJ_VAMCFM_temp;
          }
          else if (jecnum == 2){
            phjj_VAJHU_highestPTJets_dn = phjj_VAJHU_highestPTJets_temp;
            pvbf_VAJHU_highestPTJets_dn = pvbf_VAJHU_highestPTJets_temp;
	    p0plus_VAJHU_dn = p0plus_VAJHU_temp;
	    p0plus_VAMCFM_dn = p0plus_VAMCFM_temp;
	    p2bplus_VAJHU_dn = p2bplus_VAJHU_temp;
	    pqqZJJ_VAMCFM_dn = pqqZJJ_VAMCFM_temp; 
          }
        } // if hasAtLeastOneJet
      } // End if melaCand!=0
    } // for jecnum = 0 to 2
    
    //----------------------------------------------------------------------
    //--- kinematic refitting using Z mass constraint (not for merged!)

    float ZZMassRefit = -1.;
    float Z1MassRefit = -1.;
    float ZZMassUnrefitErr = -1.;
 
    if(!isMerged){
      
      vector<TLorentzVector> selectedLeptons;
      vector<TLorentzVector> selectedJets;
      // std::map<unsigned int, TLorentzVector> selectedFsrMap;
      
      for(unsigned ilep=0; ilep<4; ilep++){
	
	reco::Candidate* oneLep = (reco::Candidate*)ZZLeps[ilep];
	if (oneLep->hasMasterClone()) oneLep = (reco::Candidate*)oneLep->masterClone().get();
	TLorentzVector p4;
	p4.SetPxPyPzE(oneLep->px(),oneLep->py(),oneLep->pz(),oneLep->energy());
	
	if(FSRMap.find(ZZLeps[ilep])!=FSRMap.end()){
	  pat::PFParticle fsr = *(FSRMap[ZZLeps[ilep]]);
	  TLorentzVector p4fsr;
	  p4fsr.SetPxPyPzE(fsr.px(),fsr.py(),fsr.pz(),fsr.energy());
	  p4 += p4fsr;
	}
	if (ilep<2) selectedJets.push_back(p4);
	else selectedLeptons.push_back(p4);
	
      }
      
      kinZfitter->Setup2L2Q(selectedLeptons,selectedJets,resolution_pt,resolution_phi,rho);
      kinZfitter->KinRefitZlepZhad();
      
      // To get refit mZZ
      ZZMassRefit = kinZfitter->GetRefitMZZ2L2Q();
      // To get refit hadronic mZ (mjj)
      Z1MassRefit = kinZfitter->GetRefitMZhad();
      ZZMassUnrefitErr = -1.;
      
      // four 4-vectors after refitting order by Z1_1,Z1_2,Z2_1,Z2_2
      //vector<TLorentzVector> p4 = kinZfitter->GetRefitP4s(); 
    }
      
    mela->resetInputEvent();

    //----------------------------------------------------------------------
    //--- 4l vertex fits (experimental)

    //CandConstraintFit::fit(&myCand, iSetup);
    /* if (doVtxFit && abs(id11)==13 && abs(id12)==13 && abs(id21)==13 && abs(id22)==13) {
      
      edm::ESHandle<TransientTrackBuilder> theTTBuilder;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);

      //Creating a KinematicParticleFactory
      KinematicParticleFactoryFromTransientTrack factory;
      
      const ParticleMass muon_mass = 0.1056583;
      const ParticleMass electron_mass = 0.0005;
      float muon_sigma = 0.0000000001;
      float electron_sigma = 0.0000000001;
      
      //initial chi2 and ndf before kinematic fits. The chi2 of the reconstruction is not considered 
      float chi;
      float ndof;
  
      vector<RefCountedKinematicParticle> Particles;
      //vector<TransientTrack> t_tks;
    
      for (unsigned k = 0; k < myCand.numberOfDaughters(); ++k ) {
        const reco::Candidate* Z = myCand.daughter(k);
        for (unsigned l = 0; l < Z->numberOfDaughters(); ++l ) {
          chi = 0.; ndof = 0.;

          const reco::Candidate* lepton= Z->daughter(l);
          
          if (lepton->isGlobalMuon() || lepton->isTrackerMuon()){
            TransientTrack tt = theTTBuilder->build(lepton->get<TrackRef>());
            Particles.push_back(factory.particle (tt,muon_mass,chi,ndof,muon_sigma));
            //t_tks.push_back(tt);
          }
          else if (lepton->isElectron()){
            TransientTrack tt = theTTBuilder->build(lepton->get<GsfTrackRef>());
            Particles.push_back(factory.particle (tt,electron_mass,chi,ndof,electron_sigma));
            //t_tks.push_back(tt);
          }
        }
      }
      
      //cout << "Number of particle for constrain fitter= " << Particles.size()<< endl;
      
      if (Particles.size()>=4){
        KinematicParticleVertexFitter fitter; 
        RefCountedKinematicTree myTree = fitter.fit(Particles); 
        
        if ( !myTree->isEmpty()) {
          //accessing the tree components
          myTree->movePointerToTheTop();
          
          RefCountedKinematicParticle allLeptonsCand     = myTree->currentParticle();
          RefCountedKinematicVertex allLeptonsVertex     = myTree->currentDecayVertex();
          
          // if(dbg) cout << "m(" << myCand->numberOfDaughters() << "l): " << allLeptonsCand->currentState().mass() << " +- "
          //                   << sqrt(allLeptonsCand->currentState().kinematicParametersError().matrix()(6,6) ) << endl;
          
          reco::Vertex constrainedVertex(reco::Vertex::Point(allLeptonsVertex->position()),
                                         allLeptonsVertex->error().matrix_new(), 
                                         allLeptonsVertex->chiSquared(), 
                                         allLeptonsVertex->degreesOfFreedom(),0);       
          
          //      if(dbg) cout << "kinematicFit vertex, ndof, chi2, prob: " 
          //                   << allLeptonsVertex->position() << " , " 
          //                   << allLeptonsVertex->degreesOfFreedom() << " , "
          //                   << allLeptonsVertex->chiSquared()   << " , "
          //                   << TMath::Prob(allLeptonsVertex->chiSquared(),allLeptonsVertex->degreesOfFreedom()) << endl;
          
          //myCand->addUserData("ConstrainedCandVtx",constrainedVertex);
          myCand.addUserFloat("CFitM",allLeptonsCand->currentState().mass());
          myCand.addUserFloat("CFitSigmaM",allLeptonsCand->currentState().kinematicParametersError().matrix()(6,6));
          myCand.addUserFloat("CFitNdof",allLeptonsVertex->degreesOfFreedom());
          myCand.addUserFloat("CFitChi2",allLeptonsVertex->chiSquared());

        } else {
          cout << " ERROR CandConstraintFit: KinematicParticleVertexFitter failed " << endl;
        } 
      }
      } */
  
    //----------------------------------------------------------------------
    //--- Embed variables
    myCand.addUserFloat("candChannel",    candChannel);
    myCand.addUserFloat("SIP4",           SIP4);
    myCand.addUserFloat("pt1",            ptS[1]); // leading-pT
    myCand.addUserFloat("pt2",            ptS[0]); // sub-leading pT
    // myCand.addUserFloat("mZa",            mZa);
    // myCand.addUserFloat("mZb",            mZb);
    // myCand.addUserFloat("ZaID",           ZaID);
    // myCand.addUserFloat("ZbID",           ZbID);
    // myCand.addUserFloat("mZalpha",        mZalpha);
    // myCand.addUserFloat("mZbeta",         mZbeta);
    // myCand.addUserFloat("ZalphaID",       ZalphaID);
    // myCand.addUserFloat("ZbetaID",        ZbetaID);
    myCand.addUserFloat("mLL4",           Z2->mass()); // smallest mass of any AF/OS pair
    myCand.addUserFloat("mLL6",           Z2->mass());   // smallest mass of any AF/AS pair
    myCand.addUserFloat("passSmartMLL",   passSmartMLL);
    myCand.addUserFloat("costheta1",      costheta1);
    myCand.addUserFloat("costheta2",      costheta2);
    myCand.addUserFloat("phi",            phi);
    myCand.addUserFloat("costhetastar",   costhetastar);
    myCand.addUserFloat("phistar1",       phistar1);
    myCand.addUserFloat("xistar",         xistar);  //azimuthal angle of higgs in rest frame of higgs
    myCand.addUserFloat("xi",             xi);      //azimuthal angle of higgs in lab frame

    myCand.addUserFloat("m4l",            (Z1Lm->p4()+Z1Lp->p4()+Z2Lm->p4()+Z2Lp->p4()).mass()); // mass without FSR

    // if(!isMerged) {
    myCand.addUserFloat("ZZMassRefit"   , ZZMassRefit);
    myCand.addUserFloat("Z1MassRefit", Z1MassRefit);
    myCand.addUserFloat("ZZMassUnrefitErr", ZZMassUnrefitErr);
      // }

    // Jet quantities
    myCand.addUserFloat("DiJetMass", DiJetMass);
    myCand.addUserFloat("DiJetDEta", DiJetDEta);
    myCand.addUserFloat("DiJetFisher", DiJetFisher);

//Mela v2
    myCand.addUserFloat("p0plus_VAJHU", p0plus_VAJHU);
    myCand.addUserFloat("p0minus_VAJHU", p0minus_VAJHU);
    myCand.addUserFloat("p0hplus_VAJHU", p0hplus_VAJHU);
    myCand.addUserFloat("p2bplus_VAJHU", p2bplus_VAJHU);
    myCand.addUserFloat("p2_VAJHU", p2_VAJHU);
    myCand.addUserFloat("phjj_VAJHU_highestPTJets", phjj_VAJHU_highestPTJets);
    myCand.addUserFloat("pvbf_VAJHU_highestPTJets", pvbf_VAJHU_highestPTJets);
    myCand.addUserFloat("phjj_VAJHU_highestPTJets_up", phjj_VAJHU_highestPTJets_up);
    myCand.addUserFloat("p0plus_VAJHU_up", p0plus_VAJHU_up);
    myCand.addUserFloat("p0plus_VAMCFM_up", p0plus_VAMCFM_up);
    myCand.addUserFloat("p2bplus_VAJHU_up", p2bplus_VAJHU_up);
    myCand.addUserFloat("pqqZJJ_VAMCFM_up", pqqZJJ_VAMCFM_up);
    myCand.addUserFloat("pvbf_VAJHU_highestPTJets_up", pvbf_VAJHU_highestPTJets_up);
    myCand.addUserFloat("phjj_VAJHU_highestPTJets_dn", phjj_VAJHU_highestPTJets_dn);
    myCand.addUserFloat("p0plus_VAJHU_dn", p0plus_VAJHU_dn);
    myCand.addUserFloat("p0plus_VAMCFM_dn", p0plus_VAMCFM_dn);
    myCand.addUserFloat("p2bplus_VAJHU_dn", p2bplus_VAJHU_dn);
    myCand.addUserFloat("pqqZJJ_VAMCFM_dn", pqqZJJ_VAMCFM_dn);
    myCand.addUserFloat("pvbf_VAJHU_highestPTJets_dn", pvbf_VAJHU_highestPTJets_dn);

    myCand.addUserFloat("pqqZJJ_VAMCFM", pqqZJJ_VAMCFM);
    myCand.addUserFloat("bkg_VAMCFM", bkg_VAMCFM);
    myCand.addUserFloat("p0plus_VAMCFM", p0plus_VAMCFM);
    myCand.addUserFloat("ggzz_VAMCFM", ggzz_VAMCFM);
    myCand.addUserFloat("ggzz_p0plus_VAMCFM", ggzz_p0plus_VAMCFM);
    myCand.addUserFloat("Dgg10_VAMCFM", Dgg10_VAMCFM);

    // VH
    // myCand.addUserFloat("pzh_VAJHU",pzh_VAJHU);

    //--- MC matching. To be revised, cf. MuFiller, EleFiller
//     if (isMC) {
//       int refID = 25; // FIXME: handle ZZ (sigId = 23)
//       bool MC_isRight = (myCand.userFloat("d0.d0.MCParentCode")==refID &&
//                       myCand.userFloat("d0.d1.MCParentCode")==refID &&
//                       myCand.userFloat("d1.d0.MCParentCode")==refID &&
//                       myCand.userFloat("d1.d1.MCParentCode")==refID);
//       bool MC_isRightPair = false; //FIXME to be 

//       myCand.addUserFloat("MC_isRight",     MC_isRight);
//       myCand.addUserFloat("MC_isRightPair", MC_isRightPair);
//     }
 


    //----------------------------------------------------------------------    
    //--- Check if candedate passes the "bestCandAmong" selections (2011 PRL logic)
    int iCRname=0;
    for(CutSet<pat::CompositeCandidate>::const_iterator bca = preBestCandSelection.begin(); bca != preBestCandSelection.end(); ++bca){
      int preBestCandResult= int((*(bca->second))(myCand));     
      
      bool okSubjets = true;
      if (Z1->numberOfDaughters() < 2) okSubjets = false;
    // In miniAOD,daughters 1 and 2 are the soft-drop subjets if these exist, otherwise they are two random constituents. Check here that subjets are really there.
      if (isMerged && okSubjets) {
	const pat::Jet* testjet = dynamic_cast <const pat::Jet*> (Z1->masterClone().get());
	if (testjet->hasSubjets("SoftDrop") ) {
	  const pat::JetPtrCollection subJets = testjet->subjets("SoftDrop");
	  if (subJets.size() != 2) okSubjets = false;
	  if (okSubjets && fabs(Z1J1->pt() - subJets.at(0)->pt()) > 0.002) okSubjets = false; 
	}
      }  
      
      if (preBestCandResult && okSubjets){
        // Fill preSelCands matrix
        preSelCands[iCRname].push_back(icand);
      }
      iCRname++;
    }
    
    result->push_back(myCand);

  } // End of loop over input candidates

  // if (LLLLCands->size()>1) cout << "Size of candidates " << candidateLabel << " = " << (int)result->size() << endl; 
  //--- For each of the bestCandAmong preselections, find the best candidate and store its index (bestCandIdx)
  Comparators::BestCandComparator myComp(*result, bestCandType);
  for (int iCRname=0; iCRname<(int)preSelCands.size(); ++iCRname) {
    if (preSelCands[iCRname].size() > 0) {
      bestCandIdx[iCRname] = *std::min_element( preSelCands[iCRname].begin(), preSelCands[iCRname].end(), myComp);
    }
  }

  //--- Embed best candidate flag (must be done in a separate loop)
  for (int i = 0; i< (int)result->size(); ++i) {
    pat::CompositeCandidate& myCand = (*result)[i];

    int iCRname=0;
    for(CutSet<pat::CompositeCandidate>::const_iterator bca = preBestCandSelection.begin(); bca != preBestCandSelection.end(); ++bca){
      bool isBestCand = (i==bestCandIdx[iCRname]);
      myCand.addUserFloat(bca->first,isBestCand);
      iCRname++;
    }

    //--- Embed flags (ie cuts specified in the "flags" pset).
    //    We do this here so that isBestCand is available within the cuts.
    for(CutSet<pat::CompositeCandidate>::const_iterator cut = cuts.begin(); cut != cuts.end(); ++cut) {
      myCand.addUserFloat(cut->first,int((*(cut->second))(myCand)));
    }
  }
  
  iEvent.put(result);

}


  
void 
ZZjjCandidateFiller::getPairMass(const reco::Candidate* lp, const reco::Candidate* lm, ZZjjCandidateFiller::FSRToLepMap& photons, float& mass, int& ID) {
  math::XYZTLorentzVector llp4 = lp->p4()+lm->p4();
  auto lpp = photons.find(lp);
  auto lmp = photons.find(lm);
  if (lpp!=photons.end()) llp4+=lpp->second->p4();
  if (lmp!=photons.end()) llp4+=lmp->second->p4();
  mass=llp4.mass();
  ID=lp->pdgId()*lm->pdgId();
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ZZjjCandidateFiller);

