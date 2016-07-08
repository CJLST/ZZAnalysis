/** \class ZZCandidateFiller
 *
 *
 *  \author N. Amapane - Torino
 *  \author C. Botta - Torino
 *  \author G. Ortona - LLR
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

#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h>
#include <RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h>
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>

#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
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

bool doVtxFit = false;

class ZZCandidateFiller : public edm::EDProducer {
public:
  /// Constructor
  explicit ZZCandidateFiller(const edm::ParameterSet&);

  /// Destructor
  ~ZZCandidateFiller();

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
  bool ZRolesByMass;
  reco::CompositeCandidate::role_collection rolesZ1Z2;
  reco::CompositeCandidate::role_collection rolesZ2Z1;
  bool isMC;
  bool doKinFit;
  TH2F* corrSigmaMu;
  TH2F* corrSigmaEle;
  Comparators::ComparatorTypes bestCandType;
  KinZfitter *kinZfitter;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetToken;
  edm::EDGetTokenT<pat::METCollection> metToken;
  edm::EDGetTokenT<edm::View<reco::Candidate> > softLeptonToken;
  edm::EDGetTokenT<edm::View<reco::CompositeCandidate> > ZCandToken;
};


ZZCandidateFiller::ZZCandidateFiller(const edm::ParameterSet& iConfig) :
  candidateToken(consumes<edm::View<reco::CompositeCandidate> >(iConfig.getParameter<edm::InputTag>("src"))),
  preBestCandSelection(iConfig.getParameter<edm::ParameterSet>("bestCandAmong")),
  cuts(iConfig.getParameter<edm::ParameterSet>("flags")),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  superMelaMass(iConfig.getParameter<double>("superMelaMass")),
  embedDaughterFloats(iConfig.getUntrackedParameter<bool>("embedDaughterFloats",true)),
  ZRolesByMass(iConfig.getParameter<bool>("ZRolesByMass")),
  isMC(iConfig.getParameter<bool>("isMC")),
  doKinFit(iConfig.getParameter<bool>("doKinFit")),
  corrSigmaMu(0),
  corrSigmaEle(0),
  kinZfitter(0)
{
  produces<pat::CompositeCandidateCollection>();

  mela = new Mela(SetupToSqrts(setup), superMelaMass, TVar::SILENT);
  mela->setCandidateDecayMode(TVar::CandidateDecay_ZZ);

  jetToken = consumes<edm::View<pat::Jet> >(edm::InputTag("cleanJets"));
  metToken = consumes<pat::METCollection>(edm::InputTag("slimmedMETs"));
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
  else if (cmp=="byBestKD")       bestCandType=Comparators::byBestKD;
  else if (cmp=="byBestKD_VH")    bestCandType=Comparators::byBestKD_VH;
  else if (cmp=="byBestPsig")    bestCandType=Comparators::byBestPsig;
  else if (cmp=="byMHWindow")    bestCandType=Comparators::byMHWindow;
  else abort();

  //-- kinematic refitter
  kinZfitter = new KinZfitter(!isMC);

}

ZZCandidateFiller::~ZZCandidateFiller(){
  delete kinZfitter;
  delete mela;
}


void ZZCandidateFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  using namespace std;
  using namespace reco;

  std::auto_ptr<pat::CompositeCandidateCollection> result(new pat::CompositeCandidateCollection);

  const float ZmassValue = PDGHelpers::Zmass;

  // Get LLLL candidates
  Handle<edm::View<CompositeCandidate> > LLLLCands;
  iEvent.getByToken(candidateToken, LLLLCands);

  // Get jets
  Handle<edm::View<pat::Jet> > CleanJets;
  iEvent.getByToken(jetToken, CleanJets);

  // Get MET
  float PFMET = 0.;
  float PFMETPhi = 0.;
  Handle<pat::METCollection> metHandle;
  iEvent.getByToken(metToken, metHandle);
  if(metHandle.isValid()){
    PFMET = metHandle->front().pt();
    PFMETPhi = metHandle->front().phi();
  }

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
  Handle<View<CompositeCandidate> > ZCands;
  iEvent.getByToken(ZCandToken, ZCands);

  // Get processID
//   edm::Handle<GenEventInfoProduct> gen;
//   iEvent.getByLabel( "generator", gen );
//   int processID = gen->signalProcessID();

  // to calculate mass resolution
  CompositeCandMassResolution errorBuilder;
  errorBuilder.init(iSetup);

  vector<int> bestCandIdx(preBestCandSelection.size(),-1);
  vector<float> maxPtSum(preBestCandSelection.size(),-1);
  vector< vector<int> > preSelCands(preBestCandSelection.size());
  //----------------------------------------------------------------------
  //--- Loop over input candidates
  for( View<CompositeCandidate>::const_iterator cand = LLLLCands->begin(); cand != LLLLCands->end(); ++ cand ) {
    int icand = distance(LLLLCands->begin(),cand);

    pat::CompositeCandidate myCand(*cand);

    if (embedDaughterFloats){
      userdatahelpers::embedDaughterData(myCand);
    }

    //--- Set id of the Z1 and "Z1"/"Z2" labels. This allows to call e.g. aHiggs->daughter("Z1").
    // if ZRolesByMass is true, 'Z1' and iZ1 refer to the  Z closest to mZ; otherwise Z1 = daughter(0). The latter is used for control regions.
    const reco::CompositeCandidate::role_collection* ZRoles = &rolesZ1Z2;
    int iZ1 = 0;
    int iZ2 = 1;
    if (ZRolesByMass) {
      if(std::abs(myCand.daughter(0)->mass()-ZmassValue)>=std::abs(myCand.daughter(1)->mass()-ZmassValue)){
        swap(iZ1,iZ2);
        ZRoles = &rolesZ2Z1;
      }
    }
    myCand.setRoles(*ZRoles);
    myCand.applyRoles();

    //--- Z pointers
    const reco::Candidate* Z1= myCand.daughter(iZ1);
    const reco::Candidate* Z2= myCand.daughter(iZ2);
    vector<const reco::Candidate*> Zs = {Z1, Z2}; // in the original order

    //--- Lepton pointers in the original order
    const reco::Candidate* Z1L1= Z1->daughter(0);
    const reco::Candidate* Z1L2= Z1->daughter(1);
    const reco::Candidate* Z2L1= Z2->daughter(0);
    const reco::Candidate* Z2L2= Z2->daughter(1);
    vector<const reco::Candidate*> ZZLeps = {Z1L1,Z1L2,Z2L1,Z2L2}; // array, in the original order

    // Create corresponding array of fourmomenta; will add FSR (below)
    vector<math::XYZTLorentzVector> pij(4);
    std::transform(ZZLeps.begin(), ZZLeps.end(),pij.begin(), [](const reco::Candidate* c){return c->p4();});

    //--- Collect FSR photons and map them to the corresponding leptons
    FSRToLepMap FSRMap;
    for (unsigned iZ=0; iZ<2; ++iZ) {
      for (unsigned ifsr=2; ifsr<Zs[iZ]->numberOfDaughters(); ++ifsr) {
    const pat::PFParticle* fsr = static_cast<const pat::PFParticle*>(Zs[iZ]->daughter(ifsr));
    int ilep = iZ*2+fsr->userFloat("leptIdx");
    FSRMap[ZZLeps[ilep]]= fsr;
    pij[ilep]+=fsr->p4();
      }
    }

    //--- Lepton four-vectors in the original order; with FSR added
    math::XYZTLorentzVector p11 = pij[0];
    math::XYZTLorentzVector p12 = pij[1];
    math::XYZTLorentzVector p21 = pij[2];
    math::XYZTLorentzVector p22 = pij[3];
    int id11 = Z1L1->pdgId();
    int id12 = Z1L2->pdgId();
    int id21 = Z2L1->pdgId();
    int id22 = Z2L2->pdgId();
    int candChannel = id11*id12*id21*id22;

    if((id11 == 22 && id12 == 22) || (id21 == 22 && id22 == 22)) LogError("Z with 2 tle") << "Found a Z candidate made up of 2 trackless electrons";

    if(id11 == 22) id11 = -1 * id12;
    if(id12 == 22) id12 = -1 * id11;
    if(id21 == 22) id21 = -1 * id22;
    if(id22 == 22) id22 = -1 * id21;


    // Compute worst-lepton isolation
    for (int zIdx=0; zIdx<2; ++zIdx) {
      float worstMuIso=0;
      float worstEleIso=0;
      for (int dauIdx=0; dauIdx<2; ++dauIdx) {
    const reco::Candidate* z = myCand.daughter(zIdx);
    const reco::Candidate* d = z->daughter(dauIdx);
    float combRelIsoPFCorr = userdatahelpers::getUserFloat(d,"combRelIsoPFFSRCorr");
    if (d->isMuon()) worstMuIso  = max(worstMuIso,  combRelIsoPFCorr);
    else             worstEleIso = max(worstEleIso, combRelIsoPFCorr);
      }
      string base = (zIdx==0?"d0.":"d1.");
      myCand.addUserFloat(base+"worstMuIso",worstMuIso);
      myCand.addUserFloat(base+"worstEleIso",worstEleIso);
    }


    //----------------------------------------------------------------------
    //--- Alternative lepton pairings: "smart cut" and QCD suppression and

    //--- Sign-ordered leptons and leptopn four-vectors (without FSR), to be used to compute mZa, mZb, mZalpha, mZbeta
    const reco::Candidate* Z1Lp(Z1L1);
    const reco::Candidate* Z1Lm(Z1L2);
    const reco::Candidate* Z2Lp(Z2L1);
    const reco::Candidate* Z2Lm(Z2L2);

    // Sort leptons for OS Z candidates; no sorting for the same-sign collections used for CRs
    if (Z1Lp->charge() < 0 && Z1Lp->charge()*Z1Lm->charge()<0) {
      swap(Z1Lp,Z1Lm);
    }
    if (Z2Lp->charge() < 0 && Z2Lp->charge()*Z2Lm->charge()<0) {
      swap(Z2Lp,Z2Lm);
    }

    math::XYZTLorentzVector p1p(Z1Lp->p4());
    math::XYZTLorentzVector p1m(Z1Lm->p4());
    math::XYZTLorentzVector p2p(Z2Lp->p4());
    math::XYZTLorentzVector p2m(Z2Lm->p4());

    // Build the other SF/OS combination
    float mZ1= Z1->mass();
    float mZa, mZb;
    int ZaID, ZbID;
    getPairMass(Z1Lp,Z2Lm,FSRMap,mZa,ZaID);
    getPairMass(Z1Lm,Z2Lp,FSRMap,mZb,ZbID);

    // For same-sign CRs, the Z2 leptons are same sign, so we need to check also the other combination.
    float mZalpha, mZbeta;
    int ZalphaID, ZbetaID;
    getPairMass(Z1Lp,Z2Lp,FSRMap,mZalpha,ZalphaID);
    getPairMass(Z1Lm,Z2Lm,FSRMap,mZbeta,ZbetaID);

    // Sort (mZa,mZb and) (mZalpha,mZbeta) so that a and alpha are the ones closest to mZ
    if (std::abs(mZa-ZmassValue)>=std::abs(mZb-ZmassValue)) {
      swap(mZa,mZb);
      swap(ZaID,ZbID);
    }
    if (std::abs(mZalpha-ZmassValue)>=std::abs(mZbeta-ZmassValue)) {
      swap(mZalpha,mZbeta);
      swap(ZalphaID,ZbetaID);
    }

    // "smart cut" mll logic: veto the candidate if by swapping leptons we find a better Z1 and the Z2 is below 12 GeV.
    // To handle same-sign CRs, we have to check both alternate pairings, and consider those that have a SF/OS Z1.
    bool passSmartMLL = true;
    if (((ZaID==-121||ZaID==-169) && std::abs(mZa-ZmassValue)<std::abs(mZ1-ZmassValue) && mZb<12) ||
        ((ZalphaID==-121||ZalphaID==-169) && std::abs(mZalpha-ZmassValue)<std::abs(mZ1-ZmassValue) && mZbeta<12)) passSmartMLL = false;


    //--- QCD suppression cut
    vector<const reco::Candidate*> lep;
    lep.push_back(Z1Lm);
    lep.push_back(Z1Lp);
    lep.push_back(Z2Lm);
    lep.push_back(Z2Lp);

    float mll6 = 9999;
    float mll4 = 9999;
    for (int i=0;i<4;++i) {
      for (int j=i+1;j<4;++j) {
        float mll = (lep[i]->p4()+lep[j]->p4()).mass();
        mll6 = min(mll, mll6);
        if (lep[i]->charge()*lep[j]->charge()<0) { //OS
          mll4 = min (mll,mll4);
        }
      }
    }


    //--- worst SIP value
    vector<double> SIPS = {myCand.userFloat("d0.d0.SIP"), myCand.userFloat("d0.d1.SIP"), myCand.userFloat("d1.d0.SIP"), myCand.userFloat("d1.d1.SIP")};
    sort(SIPS.begin(),SIPS.end());
    float SIP4 = SIPS[3];

    //--- Sorted pTs
    vector<double> ptS;
    ptS.push_back(Z1Lm->pt());
    ptS.push_back(Z1Lp->pt());
    ptS.push_back(Z2Lm->pt());
    ptS.push_back(Z2Lp->pt());
    sort(ptS.begin(),ptS.end());

    //--- Mass and Lepton uncertainties
    std::vector<double> errs;
    float massError = errorBuilder.getMassResolutionWithComponents(myCand, errs);
    int offset =0;
    float sigma[2][3] = {{0,0,0}, {0,0,0}};

    myCand.addUserFloat("massError",      massError);
    myCand.addUserFloat("massError11",    errs[0]);
    sigma[0][0] = errs[0];
    myCand.addUserFloat("massError12",    errs[1]);
    sigma[0][1] = errs[1];
    if (myCand.daughter(0)->numberOfDaughters()==3){
      myCand.addUserFloat("massError13",    errs[2]);
      sigma[0][2] = errs[2];
      offset = 1;
    }
    myCand.addUserFloat("massError21",    errs[2+offset]);
    sigma[1][0] = errs[2+offset];
    myCand.addUserFloat("massError22",    errs[3+offset]);
    sigma[1][1] = errs[3+offset];
    if (myCand.daughter(1)->numberOfDaughters()==3){
      myCand.addUserFloat("massError23",    errs[4+offset]);
      sigma[1][2]=errs[4+offset];
    }

    float massErrorCorr=0;
    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
        const reco::Candidate* l=cand->daughter(i)->daughter(j);
        const TH2F* h;
        if (l->isMuon()) h = corrSigmaMu;
        else             h = corrSigmaEle;
    float ebecorr=1.;
    if (h!=0) {
      int ptBin  = min(max(1,h->GetXaxis()->FindBin(l->pt())), h->GetNbinsX());
      int etaBin = min(max(1,h->GetYaxis()->FindBin(fabs(l->eta()))), h->GetNbinsY());
      ebecorr = h->GetBinContent(ptBin, etaBin);
    }
        massErrorCorr+= (sigma[i][j]*ebecorr)*(sigma[i][j]*ebecorr);
      }
    }
    massErrorCorr += (sigma[0][2])*(sigma[0][2]);
    massErrorCorr += (sigma[1][2])*(sigma[1][2]);
    massErrorCorr = sqrt(massErrorCorr);
    myCand.addUserFloat("massErrorCorr",      massErrorCorr);


    //--- store good isolated leptons that are not involved in the current ZZ candidate
    int nExtraLep = 0;
    SimpleParticleCollection_t associatedLeptons;
    for (vector<reco::CandidatePtr>::const_iterator lepPtr = goodisoleptonPtrs.begin(); lepPtr != goodisoleptonPtrs.end(); ++lepPtr){
      const reco::Candidate* lep = lepPtr->get();
      if (
        reco::deltaR(lep->p4(), Z1L1->p4()) > 0.02 &&
        reco::deltaR(lep->p4(), Z1L2->p4()) > 0.02 &&
        reco::deltaR(lep->p4(), Z2L1->p4()) > 0.02 &&
        reco::deltaR(lep->p4(), Z2L2->p4()) > 0.02
        ){
        myCand.addUserCand("ExtraLep"+to_string(nExtraLep), *lepPtr);

        SimpleParticle_t theLepton(
        lep->pdgId(),
        TLorentzVector(lep->p4().x(), lep->p4().y(), lep->p4().z(), lep->p4().t())
        );
        bool inserted=false;
        for (SimpleParticleCollection_t::iterator ielo=associatedLeptons.begin(); ielo<associatedLeptons.end(); ielo++){
          if (lep->pt()>(*ielo).second.Pt()){
            inserted=true;
            associatedLeptons.insert(ielo, theLepton);
            break;
          }
        }
        if (!inserted) associatedLeptons.push_back(theLepton);
        nExtraLep++;
      }
    }
    myCand.addUserFloat("nExtraLep",nExtraLep);

    // Leptonically decaying WH
    if (nExtraLep>=1){
      // Take leading-pT lepton to compute fake neutrino
      int nuid = -associatedLeptons.at(0).first + (associatedLeptons.at(0).first>0 ? -1 : +1);

      // Take neutrino momentum from the MET, using a W mass constraint to solve for the z component
      float a = associatedLeptons.at(0).second.X();
      float b = associatedLeptons.at(0).second.Y();
      float c = associatedLeptons.at(0).second.Z();
      float f = associatedLeptons.at(0).second.T();
      TLorentzVector myLep(a, b, c, f);
      float x = PFMET*cos(PFMETPhi);
      float y = PFMET*sin(PFMETPhi);
      float m = PDGHelpers::Wmass;
      float delta = pow(c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y), 2) - 4*(-4*c*c + 4*f*f)*(-a*a*a*a - 2*a*a*b*b - b*b*b*b - 2*a*a*c*c - 2*b*b*c*c - c*c*c*c + 2*a*a*f*f + 2*b*b*f*f + 2*c*c*f*f - f*f*f*f - 2*a*a*m*m - 2*b*b*m*m - 2*c*c*m*m + 2*f*f*m*m - m*m*m*m - 4*a*a*a*x - 4*a*b*b*x - 4*a*c*c*x + 4*a*f*f*x - 4*a*m*m*x - 4*a*a*x*x + 4*f*f*x*x - 4*a*a*b*y - 4*b*b*b*y - 4*b*c*c*y + 4*b*f*f*y - 4*b*m*m*y - 8*a*b*x*y - 4*b*b*y*y + 4*f*f*y*y);

      if (delta>=0.){
        float z1 = (-c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y) - sqrt(delta)) / (2*(-4*c*c + 4*f*f));
        float z2 = (-c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y) + sqrt(delta)) / (2*(-4*c*c + 4*f*f));
        TLorentzVector myNu0(x, y, z1, sqrt(x*x+y*y+z1*z1));
        TLorentzVector myNu1(x, y, z2, sqrt(x*x+y*y+z2*z2));
        associatedLeptons.push_back(SimpleParticle_t(nuid, myNu0));
        associatedLeptons.push_back(SimpleParticle_t(nuid, myNu1));
      }
      else{
        TLorentzVector myNu(x, y, 0, TMath::Sqrt(x*x+y*y));
        associatedLeptons.push_back(SimpleParticle_t(nuid, myNu));
      }
    }


    //--- store Z candidates whose leptons are not involved in the current ZZ candidate
    int nExtraZ = 0;
    vector<const CompositeCandidate*> extraZs;
    for( View<CompositeCandidate>::const_iterator zcand = ZCands->begin(); zcand != ZCands->end(); ++ zcand ) {
      if( reco::deltaR( zcand->daughter(0)->p4(), Z1L1->p4() ) > 0.02 &&
      reco::deltaR( zcand->daughter(0)->p4(), Z1L2->p4() ) > 0.02 &&
      reco::deltaR( zcand->daughter(0)->p4(), Z2L1->p4() ) > 0.02 &&
      reco::deltaR( zcand->daughter(0)->p4(), Z2L2->p4() ) > 0.02 &&
      reco::deltaR( zcand->daughter(1)->p4(), Z1L1->p4() ) > 0.02 &&
      reco::deltaR( zcand->daughter(1)->p4(), Z1L2->p4() ) > 0.02 &&
      reco::deltaR( zcand->daughter(1)->p4(), Z2L1->p4() ) > 0.02 &&
      reco::deltaR( zcand->daughter(1)->p4(), Z2L2->p4() ) > 0.02    ){
    const reco::CandidatePtr myZCand(ZCands,zcand-ZCands->begin());
    if((bool)userdatahelpers::getUserFloat(&*myZCand,"GoodIsoLeptons")){
      nExtraZ++;
      extraZs.push_back(&*zcand);
      myCand.addUserCand("assocZ"+to_string(nExtraZ),myZCand);
    }
      }
    }
    myCand.addUserFloat("nExtraZ",nExtraZ);


    /**********************/
    /**********************/
    /***** BEGIN MELA *****/
    /**********************/
    /**********************/

    // Lepton TLorentzVectors, including FSR
    SimpleParticleCollection_t daughters;
    daughters.push_back(SimpleParticle_t(0, TLorentzVector(p11.x(), p11.y(), p11.z(), p11.t())));
    daughters.push_back(SimpleParticle_t(0, TLorentzVector(p12.x(), p12.y(), p12.z(), p12.t())));
    daughters.push_back(SimpleParticle_t(id21, TLorentzVector(p21.x(), p21.y(), p21.z(), p21.t())));
    daughters.push_back(SimpleParticle_t(id22, TLorentzVector(p22.x(), p22.y(), p22.z(), p22.t())));

    //--- Compute angles, better done here
    float costheta1=0, costheta2=0, phi=0, costhetastar=0, phistar1=0;
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
    float xi = higgs.Phi();
    // boost Z1 into rest frame of higgs
    // xistar is the angle between Z decay plane and x-axis
    Z14vec.Boost(-higgs.BoostVector());
    float xistar = Z14vec.Phi();
    // detaJJ, Mjj and Fisher. These are per-event variables in the SR, but not necessarily in the CR as we clean jets also
    // for loose-but-not-tight leptons.
    float DiJetMass  = -99;
    float DiJetDEta  = -99;
    float DiJetFisher  = -99;

    unsigned int nCandidates=0; // Should equal 3 after the loop below
    for (int jecnum = 0; jecnum < 3; jecnum++){
      SimpleParticleCollection_t associated;

      // multiplier: +1 means JEC up, -1 means JEC down
      double jecnum_multiplier = 0;
      if (jecnum==1) jecnum_multiplier = 1.;
      else if (jecnum==2) jecnum_multiplier = -1.;

      vector<const pat::Jet*> cleanedJetsPt30Jec;
      vector<float> jec_ratio;
      for (edm::View<pat::Jet>::const_iterator jet = CleanJets->begin(); jet != CleanJets->end(); ++jet){
        // calculate JEC uncertainty up/down
        float jec_unc = jet->userFloat("jec_unc");
        float ratio = 1. + jecnum_multiplier * jec_unc;
        float newPt = jet->pt() * ratio;
        // apply pt>30GeV cut
        if (newPt<=30.) continue;
        // additional jets cleaning for loose leptons belonging to this candidate (for CRs only;
        // does nothing for the SR as jets are already cleaned with all tight isolated leptons )
        if (!jetCleaner::isGood(myCand, *jet)) continue;
        // store jets and up/down ratio
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
      for (unsigned int ilep=0; ilep<associatedLeptons.size(); ilep++) associated.push_back(associatedLeptons.at(ilep));
      mela->setInputEvent(&daughters, &associated, 0, 0);
      for (unsigned int ijet = 0; ijet<cleanedJetsPt30Jec.size(); ijet++){
        TLorentzVector jet(
          cleanedJetsPt30Jec[ijet]->p4().x()*jec_ratio.at(ijet),
          cleanedJetsPt30Jec[ijet]->p4().y()*jec_ratio.at(ijet),
          cleanedJetsPt30Jec[ijet]->p4().z()*jec_ratio.at(ijet),
          cleanedJetsPt30Jec[ijet]->p4().t()*jec_ratio.at(ijet)
          );
        SimpleParticleCollection_t stableTopDaughters; // Just a collection with one jet as the top
        stableTopDaughters.push_back(SimpleParticle_t(0, jet));
        mela->appendTopCandidate(&stableTopDaughters);
      }
      nCandidates++;
    }

    mela->setCurrentCandidateFromIndex(0);

    /****************************/
    /***** SPIN-0 JHUGEN ME *****/
    /****************************/
    /***** PURE MEs *****/
    // 0+m
    float p0plus_VAJHU=0;
    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p0plus_VAJHU, true);
    // 0+L1
    float p0_g1prime2_VAJHU=0;
    mela->setProcess(TVar::H0_g1prime2, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p0_g1prime2_VAJHU, true);
    // 0+h
    float p0hplus_VAJHU=0;
    mela->setProcess(TVar::H0hplus, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p0hplus_VAJHU, true);
    // 0-
    float p0minus_VAJHU=0;
    mela->setProcess(TVar::H0minus, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p0minus_VAJHU, true);
    // 0+hzgs
    float p0hplus_zgs_VAJHU=0;
    mela->setProcess(TVar::H0_Zgs, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p0hplus_zgs_VAJHU, true);
    // 0+L1zgs
    float p0_g1prime2_zgs_VAJHU=0;
    mela->setProcess(TVar::H0_Zgsg1prime2, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p0_g1prime2_zgs_VAJHU, true);
    // 0-zgs
    float p0minus_zgs_VAJHU=0;
    mela->setProcess(TVar::H0_Zgs_PS, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p0minus_zgs_VAJHU, true);
    // 0+hgsgs
    float p0hplus_gsgs_VAJHU=0;
    mela->setProcess(TVar::H0_gsgs, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p0hplus_gsgs_VAJHU, true);
    // 0-gsgs
    float p0minus_gsgs_VAJHU=0;
    mela->setProcess(TVar::H0_gsgs_PS, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p0minus_gsgs_VAJHU, true);

    /***** MIXTURE MEs *****/
    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    // 0+m -- 0+h
    float pg1g2_VAJHU=0;
    (mela->selfDHggcoupl)[0][0]=1.;
    (mela->selfDHzzcoupl)[0][0][0]=1.;
    (mela->selfDHzzcoupl)[0][1][0]=1.;
    mela->computeP(pg1g2_VAJHU, true);
    pg1g2_VAJHU -= p0plus_VAJHU+p0hplus_VAJHU;
    // 0+m -- 0+h (phase(0+h)=pi/2)
    float pg1g2_pi2_VAJHU=0;
    (mela->selfDHggcoupl)[0][0]=1.;
    (mela->selfDHzzcoupl)[0][0][0]=1.;
    (mela->selfDHzzcoupl)[0][1][1]=1.;
    mela->computeP(pg1g2_pi2_VAJHU, true);
    pg1g2_pi2_VAJHU -= p0plus_VAJHU+p0hplus_VAJHU;
    // 0+m -- 0-
    float pg1g4_VAJHU=0;
    (mela->selfDHggcoupl)[0][0]=1.;
    (mela->selfDHzzcoupl)[0][0][0]=1.;
    (mela->selfDHzzcoupl)[0][3][0]=1.;
    mela->computeP(pg1g4_VAJHU, true);
    pg1g4_VAJHU -= p0plus_VAJHU+p0minus_VAJHU;
    // 0+m -- 0- (phase(0-)=pi/2)
    float pg1g4_pi2_VAJHU=0;
    (mela->selfDHggcoupl)[0][0]=1.;
    (mela->selfDHzzcoupl)[0][0][0]=1.;
    (mela->selfDHzzcoupl)[0][3][1]=1.;
    mela->computeP(pg1g4_pi2_VAJHU, true);
    pg1g4_pi2_VAJHU -= p0plus_VAJHU+p0minus_VAJHU;
    // 0+m -- 0+L1 (no phase=pi/2 equivalent)
    float pg1g1prime2_VAJHU=0;
    (mela->selfDHggcoupl)[0][0]=1.;
    (mela->selfDHzzcoupl)[0][0][0]=1.;
    (mela->selfDHzzcoupl)[0][11][0]=1.;
    mela->computeP(pg1g1prime2_VAJHU, true);
    pg1g1prime2_VAJHU -= p0plus_VAJHU+p0_g1prime2_VAJHU;
    // 0+m ZZ -- 0+L1 Zg*
    float p0plus_zz_g1prime2_zgs_VAJHU=0;
    (mela->selfDHggcoupl)[0][0]=1.;
    (mela->selfDHzzcoupl)[0][0][0]=1.;
    (mela->selfDHzzcoupl)[0][30][0]=1.;
    mela->computeP(p0plus_zz_g1prime2_zgs_VAJHU, true);
    p0plus_zz_g1prime2_zgs_VAJHU -= p0plus_VAJHU+p0_g1prime2_zgs_VAJHU;
    // 0+m ZZ -- 0+L1 Zg* (phase(0+L1 Zgs*)=pi/2)
    float p0plus_zz_g1prime2_zgs_pi2_VAJHU=0;
    (mela->selfDHggcoupl)[0][0]=1.;
    (mela->selfDHzzcoupl)[0][0][0]=1.;
    (mela->selfDHzzcoupl)[0][30][1]=1.;
    mela->computeP(p0plus_zz_g1prime2_zgs_pi2_VAJHU, true);
    p0plus_zz_g1prime2_zgs_pi2_VAJHU -= p0plus_VAJHU+p0_g1prime2_zgs_VAJHU;
    // 0+m ZZ -- 0+h Zg*
    float p0plus_zz_0hplus_zgs_VAJHU=0;
    (mela->selfDHggcoupl)[0][0]=1.;
    (mela->selfDHzzcoupl)[0][0][0]=1.;
    (mela->selfDHzzcoupl)[0][4][0]=1.;
    mela->computeP(p0plus_zz_0hplus_zgs_VAJHU, true);
    p0plus_zz_0hplus_zgs_VAJHU -= p0plus_VAJHU+p0hplus_zgs_VAJHU;
    // 0+m ZZ -- 0- Zg*
    float p0plus_zz_0minus_zgs_VAJHU=0;
    (mela->selfDHggcoupl)[0][0]=1.;
    (mela->selfDHzzcoupl)[0][0][0]=1.;
    (mela->selfDHzzcoupl)[0][6][0]=1.;
    mela->computeP(p0plus_zz_0minus_zgs_VAJHU, true);
    p0plus_zz_0minus_zgs_VAJHU -= p0plus_VAJHU+p0minus_zgs_VAJHU;
    // 0+m ZZ -- 0+h g*g*
    float p0plus_zz_0hplus_gsgs_VAJHU=0;
    (mela->selfDHggcoupl)[0][0]=1.;
    (mela->selfDHzzcoupl)[0][0][0]=1.;
    (mela->selfDHzzcoupl)[0][7][0]=1.;
    mela->computeP(p0plus_zz_0hplus_gsgs_VAJHU, true);
    p0plus_zz_0hplus_gsgs_VAJHU -= p0plus_VAJHU+p0hplus_gsgs_VAJHU;
    // 0+m ZZ -- 0- g*g*
    float p0plus_zz_0minus_gsgs_VAJHU=0;
    (mela->selfDHggcoupl)[0][0]=1.;
    (mela->selfDHzzcoupl)[0][0][0]=1.;
    (mela->selfDHzzcoupl)[0][9][0]=1.;
    mela->computeP(p0plus_zz_0minus_gsgs_VAJHU, true);
    p0plus_zz_0minus_gsgs_VAJHU -= p0plus_VAJHU+p0minus_gsgs_VAJHU;

    /****************************/
    /***** SPIN-1 JHUGEN ME *****/
    /****************************/
    // qqb 1- VV
    float p1_VAJHU=0;
    mela->setProcess(TVar::H1minus, TVar::JHUGen, TVar::ZZQQB);
    mela->computeP(p1_VAJHU, true);
    // 1- VV
    float p1_prodIndep_VAJHU=0;
    mela->setProcess(TVar::H1minus, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela->computeP(p1_prodIndep_VAJHU, true);
    // qqb 1+ VV
    float p1plus_VAJHU=0;
    mela->setProcess(TVar::H1plus, TVar::JHUGen, TVar::ZZQQB);
    mela->computeP(p1plus_VAJHU, true);
    // 1+ VV
    float p1plus_prodIndep_VAJHU=0;
    mela->setProcess(TVar::H1plus, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela->computeP(p1plus_prodIndep_VAJHU, true);

    /****************************/
    /***** SPIN-2 JHUGEN ME *****/
    /****************************/
    // Should find a nicer way -- U. Sarica
    // gg 2+m VV (g1=1, g5=1)
    float p2plus_gg_VAJHU=0;
    mela->setProcess(TVar::H2_g1g5, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p2plus_gg_VAJHU, true);
    // qqb 2+m VV (g1=1, g5=1)
    float p2plus_qqb_VAJHU=0;
    mela->setProcess(TVar::H2_g1g5, TVar::JHUGen, TVar::ZZQQB);
    mela->computeP(p2plus_qqb_VAJHU, true);
    // 2+m VV (g1=1, g5=1)
    float p2plus_prodIndep_VAJHU=0;
    mela->setProcess(TVar::H2_g1g5, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela->computeP(p2plus_prodIndep_VAJHU, true);

    // gg 2+h2 VV (g2=1)
    float p2h2plus_gg_VAJHU=0;
    mela->setProcess(TVar::H2_g2, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p2h2plus_gg_VAJHU, true);
    // qqb 2+h2 VV (g2=1)
    float p2h2plus_qqb_VAJHU=0;
    mela->setProcess(TVar::H2_g2, TVar::JHUGen, TVar::ZZQQB);
    mela->computeP(p2h2plus_qqb_VAJHU, true);
    // 2+h2 VV (g2=1)
    float p2h2plus_prodIndep_VAJHU=0;
    mela->setProcess(TVar::H2_g2, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela->computeP(p2h2plus_prodIndep_VAJHU, true);

    // gg 2+h3 VV (g3=1)
    float p2h3plus_gg_VAJHU=0;
    mela->setProcess(TVar::H2_g3, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p2h3plus_gg_VAJHU, true);
    // qqb 2+h3 VV (g3=1)
    float p2h3plus_qqb_VAJHU=0;
    mela->setProcess(TVar::H2_g3, TVar::JHUGen, TVar::ZZQQB);
    mela->computeP(p2h3plus_qqb_VAJHU, true);
    // 2+h3 VV (g3=1)
    float p2h3plus_prodIndep_VAJHU=0;
    mela->setProcess(TVar::H2_g3, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela->computeP(p2h3plus_prodIndep_VAJHU, true);

    // gg 2+h4 VV (g4=1)
    float p2h4plus_gg_VAJHU=0;
    mela->setProcess(TVar::H2_g4, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p2h4plus_gg_VAJHU, true);
    // qqb 2+h4 VV (g4=1)
    float p2h4plus_qqb_VAJHU=0;
    mela->setProcess(TVar::H2_g4, TVar::JHUGen, TVar::ZZQQB);
    mela->computeP(p2h4plus_qqb_VAJHU, true);
    // 2+h4 VV (g4=1)
    float p2h4plus_prodIndep_VAJHU=0;
    mela->setProcess(TVar::H2_g4, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela->computeP(p2h4plus_prodIndep_VAJHU, true);

    // gg 2+b VV (g5=1)
    float p2bplus_gg_VAJHU=0;
    mela->setProcess(TVar::H2_g5, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p2bplus_gg_VAJHU, true);
    // qqb 2+b VV (g5=1)
    float p2bplus_qqb_VAJHU=0;
    mela->setProcess(TVar::H2_g5, TVar::JHUGen, TVar::ZZQQB);
    mela->computeP(p2bplus_qqb_VAJHU, true);
    // 2+b VV (g5=1)
    float p2bplus_prodIndep_VAJHU=0;
    mela->setProcess(TVar::H2_g5, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela->computeP(p2bplus_prodIndep_VAJHU, true);

    // gg 2+h6 VV (g6=1)
    float p2h6plus_gg_VAJHU=0;
    mela->setProcess(TVar::H2_g6, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p2h6plus_gg_VAJHU, true);
    // qqb 2+h6 VV (g6=1)
    float p2h6plus_qqb_VAJHU=0;
    mela->setProcess(TVar::H2_g6, TVar::JHUGen, TVar::ZZQQB);
    mela->computeP(p2h6plus_qqb_VAJHU, true);
    // 2+h6 VV (g6=1)
    float p2h6plus_prodIndep_VAJHU=0;
    mela->setProcess(TVar::H2_g6, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela->computeP(p2h6plus_prodIndep_VAJHU, true);

    // gg 2+h7 VV (g7=1)
    float p2h7plus_gg_VAJHU=0;
    mela->setProcess(TVar::H2_g7, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p2h7plus_gg_VAJHU, true);
    // qqb 2+h7 VV (g7=1)
    float p2h7plus_qqb_VAJHU=0;
    mela->setProcess(TVar::H2_g7, TVar::JHUGen, TVar::ZZQQB);
    mela->computeP(p2h7plus_qqb_VAJHU, true);
    // 2+h7 VV (g7=1)
    float p2h7plus_prodIndep_VAJHU=0;
    mela->setProcess(TVar::H2_g7, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela->computeP(p2h7plus_prodIndep_VAJHU, true);

    // gg 2-h8 VV (g8=1)
    float p2hminus_gg_VAJHU=0;
    mela->setProcess(TVar::H2_g8, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p2hminus_gg_VAJHU, true);
    // qqb 2-h8 VV (g8=1)
    float p2hminus_qqb_VAJHU=0;
    mela->setProcess(TVar::H2_g8, TVar::JHUGen, TVar::ZZQQB);
    mela->computeP(p2hminus_qqb_VAJHU, true);
    // 2-h8 VV (g8=1)
    float p2hminus_prodIndep_VAJHU=0;
    mela->setProcess(TVar::H2_g8, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela->computeP(p2hminus_prodIndep_VAJHU, true);

    // gg 2-h9 VV (g9=1)
    float p2h9minus_gg_VAJHU=0;
    mela->setProcess(TVar::H2_g9, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p2h9minus_gg_VAJHU, true);
    // qqb 2-h9 VV (g9=1)
    float p2h9minus_qqb_VAJHU=0;
    mela->setProcess(TVar::H2_g9, TVar::JHUGen, TVar::ZZQQB);
    mela->computeP(p2h9minus_qqb_VAJHU, true);
    // 2-h9 VV (g9=1)
    float p2h9minus_prodIndep_VAJHU=0;
    mela->setProcess(TVar::H2_g9, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela->computeP(p2h9minus_prodIndep_VAJHU, true);

    // gg 2-h10 VV (g10=1)
    float p2h10minus_gg_VAJHU=0;
    mela->setProcess(TVar::H2_g10, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p2h10minus_gg_VAJHU, true);
    // qqb 2-h10 VV (g10=1)
    float p2h10minus_qqb_VAJHU=0;
    mela->setProcess(TVar::H2_g10, TVar::JHUGen, TVar::ZZQQB);
    mela->computeP(p2h10minus_qqb_VAJHU, true);
    // 2-h10 VV (g10=1)
    float p2h10minus_prodIndep_VAJHU=0;
    mela->setProcess(TVar::H2_g10, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela->computeP(p2h10minus_prodIndep_VAJHU, true);

    /****************************************/
    /***** SIGNAL OR BACKGROUND MCFM ME *****/
    /****************************************/
    float p0plus_VAMCFM=0;
    mela->setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela->computeP(p0plus_VAMCFM, true);
    float ggzz_p0plus_VAMCFM=0;
    mela->setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela->computeP(ggzz_p0plus_VAMCFM, true);
    float ggzz_VAMCFM=0;
    mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
    mela->computeP(ggzz_VAMCFM, true);
    float bkg_VAMCFM=0;
    mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
    mela->computeP(bkg_VAMCFM, true);
    float bkg_prodIndep_VAMCFM=0;
    mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZINDEPENDENT);
    mela->computeP(bkg_prodIndep_VAMCFM, true);
    float pZJJ_VAMCFM=0;
    mela->setProcess(TVar::bkgZJets, TVar::MCFM, TVar::JJQCD);
    mela->computeP(pZJJ_VAMCFM, true);
    float Dgg10_VAMCFM=0;
    mela->computeD_gg(TVar::MCFM, TVar::D_gg10, Dgg10_VAMCFM);

    /*********************/
    /***** SuperMELA *****/
    /*********************/
    // Signal
    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
    float p0plus_m4l=0; mela->computePM4l(TVar::SMSyst_None, p0plus_m4l);
    float p0plus_m4l_ScaleUp=0; mela->computePM4l(TVar::SMSyst_ScaleUp, p0plus_m4l_ScaleUp);
    float p0plus_m4l_ScaleDown=0; mela->computePM4l(TVar::SMSyst_ScaleDown, p0plus_m4l_ScaleDown);
    float p0plus_m4l_ResUp=0; mela->computePM4l(TVar::SMSyst_ResUp, p0plus_m4l_ResUp);
    float p0plus_m4l_ResDown=0; mela->computePM4l(TVar::SMSyst_ResDown, p0plus_m4l_ResDown);
    // Background
    mela->setProcess(TVar::bkgZZ, TVar::JHUGen, TVar::ZZGG);
    float bkg_m4l=0; mela->computePM4l(TVar::SMSyst_None, bkg_m4l);
    float bkg_m4l_ScaleUp=0; mela->computePM4l(TVar::SMSyst_ScaleUp, bkg_m4l_ScaleUp);
    float bkg_m4l_ScaleDown=0; mela->computePM4l(TVar::SMSyst_ScaleDown, bkg_m4l_ScaleDown);
    float bkg_m4l_ResUp=0; mela->computePM4l(TVar::SMSyst_ResUp, bkg_m4l_ResUp);
    float bkg_m4l_ResDown=0; mela->computePM4l(TVar::SMSyst_ResDown, bkg_m4l_ResDown);

    /*********************************************/
    /***** Probabilities for H + 1/2 leptons *****/
    /*********************************************/
    float pwh_leptonic_VAJHU=-1.;
    float pzh_leptonic_VAJHU=-1.;

    TLorentzVector highestWHMENeutrino(0, 0, 0, 0); // Keep track of the neutrino used
    if (nExtraLep>0){
      MELACandidate* melaCand = mela->getCurrentCandidate();

      if (melaCand!=0){
        /***** WH *****/
        int nNeutrinos = melaCand->getNAssociatedNeutrinos();
        for (int inu=0; inu<nNeutrinos; inu++){
          for (int disableNu=0; disableNu<nNeutrinos; disableNu++) melaCand->getAssociatedNeutrino(disableNu)->setSelected(disableNu==inu); // Disable all neutrinos other than index==inu
          // Compute WH MEs
          float pwh_leptonic_VAJHU_tmp;
          mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Lep_WH);
          mela->computeProdP(pwh_leptonic_VAJHU_tmp, true);
          // Select the MEs to record
          // For now, pick neutrino based on best SM ME, will revise later -- U. Sarica
          if (pwh_leptonic_VAJHU_tmp>pwh_leptonic_VAJHU){
            pwh_leptonic_VAJHU = pwh_leptonic_VAJHU_tmp;
            highestWHMENeutrino = melaCand->getAssociatedNeutrino(inu)->p4;
          }
        }
        for (int disableNu=0; disableNu<nNeutrinos; disableNu++) melaCand->getAssociatedNeutrino(disableNu)->setSelected(true); // Return every nu selection back to true
        /***** End WH *****/
        /****************/
        /***** ZH *****/
        unsigned int nSortedVs = melaCand->getNSortedVs(); // Be careful, sortedV==0,1 are (guaranteed to be) the ZZ daughters! One needs to start any loop from index=2.
        const unsigned int iSortedVstart=2;
        double dZmass=20000;
        int chosenZ=-1;
        // Choose the Z by mass closest to mZ (~equivalent to ordering by best SM ME but would be equally valid for BSM MEs as well)
        for (unsigned int iV=iSortedVstart; iV<nSortedVs; iV++){
          MELAParticle* associatedV = melaCand->getSortedV(iV);
          if (!PDGHelpers::isAZBoson(associatedV->id)) continue;
          if (!PDGHelpers::isALepton(associatedV->getDaughter(0)->id)) continue;
          if (fabs(associatedV->m()-PDGHelpers::Zmass)<dZmass){ dZmass=associatedV->m()-PDGHelpers::Zmass; chosenZ=(int)iV; }
        }
        if(chosenZ>=0){
          // Disable every associated Z boson and its daughters unless it is the chosen one
          for (unsigned int disableV=iSortedVstart; disableV<nSortedVs; disableV++){
            bool flag=(((int)disableV)==chosenZ);
            MELAParticle* einV = melaCand->getSortedV(disableV);
            if (PDGHelpers::isAZBoson(einV->id)){
              einV->setSelected(flag);
              for (int iVj=0; iVj<einV->getNDaughters(); iVj++) einV->getDaughter(iVj)->setSelected(flag);
            }
          }
          // Compute ZH MEs
          mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Lep_ZH);
          mela->computeProdP(pzh_leptonic_VAJHU, true);
          // Re-enable every associated Z boson and its daughters, should be exactly the same loop as used in disabling them except for the setSelected flag.
          for (unsigned int disableV=iSortedVstart; disableV<nSortedVs; disableV++){
            MELAParticle* einV = melaCand->getSortedV(disableV);
            if (PDGHelpers::isAZBoson(einV->id)){
              einV->setSelected(true);
              for (int iVj=0; iVj<einV->getNDaughters(); iVj++) einV->getDaughter(iVj)->setSelected(true);
            }
          }
        }
        /***** End ZH *****/
      } // End if melaCand!=0
    } // End if nExtraLep>0

    /******************************************/
    /***** Probabilities for H + 1/2 jets *****/
    /******************************************/
    float pvbf_VAJHU_highestPTJets=-1;
    float phjj_VAJHU_highestPTJets=-1;
    float pvbf_VAJHU_highestPTJets_up=-1;
    float phjj_VAJHU_highestPTJets_up=-1;
    float pvbf_VAJHU_highestPTJets_dn=-1;
    float phjj_VAJHU_highestPTJets_dn=-1;

    float pvbf_VAJHU_bestDjet=-1;
    float phjj_VAJHU_bestDjet=-1;
    float pvbf_VAJHU_bestDjet_up=-1;
    float phjj_VAJHU_bestDjet_up=-1;
    float pvbf_VAJHU_bestDjet_dn=-1;
    float phjj_VAJHU_bestDjet_dn=-1;

    float pAux_vbf_VAJHU = 1;
    float phj_VAJHU = -1;
    float pwh_hadronic_VAJHU = -1;
    float pzh_hadronic_VAJHU = -1;
    float ptth_VAJHU = -1;
    float pbbh_VAJHU = -1;
    float pAux_vbf_VAJHU_up = 1;
    float phj_VAJHU_up = -1;
    float pwh_hadronic_VAJHU_up = -1;
    float pzh_hadronic_VAJHU_up = -1;
    float ptth_VAJHU_up = -1;
    float pbbh_VAJHU_up = -1;
    float pAux_vbf_VAJHU_dn = 1;
    float phj_VAJHU_dn = -1;
    float pwh_hadronic_VAJHU_dn = -1;
    float pzh_hadronic_VAJHU_dn = -1;
    float ptth_VAJHU_dn = -1;
    float pbbh_VAJHU_dn = -1;

    // Do these loops at the end to avoid switching particles off first and then on again
    for (unsigned int jecnum=0; jecnum<nCandidates; jecnum++){
      mela->setCurrentCandidateFromIndex(jecnum);
      MELACandidate* melaCand = mela->getCurrentCandidate();

      if (melaCand!=0){
        unsigned int nGoodJets=melaCand->getNAssociatedJets();
        bool hasAtLeastOneJet = (nGoodJets>0);
        bool hasAtLeastTwoJets = (nGoodJets>1);
        if (hasAtLeastOneJet){
          float phjj_VAJHU_highestPTJets_temp = -1;
          float pvbf_VAJHU_highestPTJets_temp = -1;

          float djet_max = -1;
          float phjj_VAJHU_bestDjet_temp = -1;
          float pvbf_VAJHU_bestDjet_temp = -1;
          float pAux_vbf_VAJHU_temp = 1;
          float phj_VAJHU_temp = -1;
          float pwh_hadronic_VAJHU_temp = -1;
          float pzh_hadronic_VAJHU_temp = -1;
          float ptth_VAJHU_temp = -1;
          float pbbh_VAJHU_temp = -1;

          for (unsigned int firstjet = 0; firstjet < nGoodJets; firstjet++){ // Loop over first jet
            for (unsigned int secondjet = 1; secondjet < nGoodJets; secondjet++){ // Loop over second jet
              if (secondjet<=firstjet) continue;

              // Disable jets and tops
              for (unsigned int disableJet=0; disableJet<nGoodJets; disableJet++) melaCand->getAssociatedJet(disableJet)->setSelected((disableJet==firstjet || disableJet==secondjet)); // Disable the other jets
              unsigned int nDisabledStableTops=0;
              for (int itop=0; itop<melaCand->getNAssociatedTops(); itop++){
                MELATopCandidate* einTop = melaCand->getAssociatedTop(itop);
                if (einTop->getNDaughters()==3) einTop->setSelected(false); // All unstable tops are disabled in the loop for jets (where "jet"=="stable top") since we are looping over jecnum
                else{
                  einTop->setSelected((nDisabledStableTops==firstjet || nDisabledStableTops==secondjet)); // Disable the other stable tops
                  nDisabledStableTops++;
                }
              }

              float pvbf_temp = -1;
              mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
              mela->computeProdP(pvbf_temp, true);
              float phjj_temp = -1;
              mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJQCD);
              mela->computeProdP(phjj_temp, true);

              double djet_temp = pvbf_temp / phjj_temp;
              if (djet_temp > djet_max){
                phjj_VAJHU_bestDjet_temp = phjj_temp;
                pvbf_VAJHU_bestDjet_temp = pvbf_temp;
                djet_max = djet_temp;
              }

              if (firstjet == 0 && secondjet == 1){ // If leading/subleading-pT jets
                phjj_VAJHU_highestPTJets_temp = phjj_temp;
                pvbf_VAJHU_highestPTJets_temp = pvbf_temp;

                float pwh_temp = -1;
                float pzh_temp = -1;
                float ptth_temp = -1;
                float pbbh_temp = -1;
                mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_WH);
                mela->computeProdP(pwh_temp, true);
                mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_ZH);
                mela->computeProdP(pzh_temp, true);
                mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ttH);
                mela->computeProdP(ptth_temp, true);
                mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::bbH);
                mela->computeProdP(pbbh_temp, true);
                pwh_hadronic_VAJHU_temp = pwh_temp;
                pzh_hadronic_VAJHU_temp = pzh_temp;
                ptth_VAJHU_temp = ptth_temp;
                pbbh_VAJHU_temp = pbbh_temp;
              }
            } // End loop over second jet

            if (!hasAtLeastTwoJets){ // Compute H + 1 jet
              for (unsigned int disableJet=0; disableJet<nGoodJets; disableJet++) melaCand->getAssociatedJet(disableJet)->setSelected((disableJet==firstjet)); // Disable everything except the single jet

              float phj_temp = -1;
              float pjvbf_temp = -1;
              float pAux_vbf_temp = -1;
              mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JQCD);
              mela->computeProdP(phj_temp, true);
              mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
              mela->computeProdP(pjvbf_temp, true); // Un-integrated ME
              mela->getPAux(pAux_vbf_temp); // = Integrated / un-integrated

              djet_max = pjvbf_temp*pAux_vbf_temp / phj_temp;
              phj_VAJHU_temp = phj_temp;
              pvbf_VAJHU_highestPTJets_temp = pjvbf_temp;
              pAux_vbf_VAJHU_temp = pAux_vbf_temp;
              pvbf_VAJHU_bestDjet_temp = pvbf_VAJHU_highestPTJets_temp;
            }
          } // End loop over first jet

          for (unsigned int disableJet=0; disableJet<nGoodJets; disableJet++) melaCand->getAssociatedJet(disableJet)->setSelected(true); // Turn all jets back on
          for (int itop=0; itop<melaCand->getNAssociatedTops(); itop++) melaCand->getAssociatedTop(itop)->setSelected(true); // Turn all tops back on

          if (jecnum == 0){
            phjj_VAJHU_highestPTJets = phjj_VAJHU_highestPTJets_temp;
            pvbf_VAJHU_highestPTJets = pvbf_VAJHU_highestPTJets_temp;
            phjj_VAJHU_bestDjet = phjj_VAJHU_bestDjet_temp;
            pvbf_VAJHU_bestDjet = pvbf_VAJHU_bestDjet_temp;
            pAux_vbf_VAJHU = pAux_vbf_VAJHU_temp;
            phj_VAJHU = phj_VAJHU_temp;
            pwh_hadronic_VAJHU = pwh_hadronic_VAJHU_temp;
            pzh_hadronic_VAJHU = pzh_hadronic_VAJHU_temp;
            ptth_VAJHU = ptth_VAJHU_temp;
            pbbh_VAJHU = pbbh_VAJHU_temp;
          }
          else if (jecnum == 1){
            phjj_VAJHU_highestPTJets_up = phjj_VAJHU_highestPTJets_temp;
            pvbf_VAJHU_highestPTJets_up = pvbf_VAJHU_highestPTJets_temp;
            phjj_VAJHU_bestDjet_up = phjj_VAJHU_bestDjet_temp;
            pvbf_VAJHU_bestDjet_up = pvbf_VAJHU_bestDjet_temp;
            pAux_vbf_VAJHU_up = pAux_vbf_VAJHU_temp;
            phj_VAJHU_up = phj_VAJHU_temp;
            pwh_hadronic_VAJHU_up = pwh_hadronic_VAJHU_temp;
            pzh_hadronic_VAJHU_up = pzh_hadronic_VAJHU_temp;
            ptth_VAJHU_up = ptth_VAJHU_temp;
            pbbh_VAJHU_up = pbbh_VAJHU_temp;
          }
          else if (jecnum == 2){
            phjj_VAJHU_highestPTJets_dn = phjj_VAJHU_highestPTJets_temp;
            pvbf_VAJHU_highestPTJets_dn = pvbf_VAJHU_highestPTJets_temp;
            phjj_VAJHU_bestDjet_dn = phjj_VAJHU_bestDjet_temp;
            pvbf_VAJHU_bestDjet_dn = pvbf_VAJHU_bestDjet_temp;
            pAux_vbf_VAJHU_dn = pAux_vbf_VAJHU_temp;
            phj_VAJHU_dn = phj_VAJHU_temp;
            pwh_hadronic_VAJHU_dn = pwh_hadronic_VAJHU_temp;
            pzh_hadronic_VAJHU_dn = pzh_hadronic_VAJHU_temp;
            ptth_VAJHU_dn = ptth_VAJHU_temp;
            pbbh_VAJHU_dn = pbbh_VAJHU_temp;
          }
        } // End if hasAtLeastOneJet
      } // End if melaCand!=0
    } // End for jecnum = 0 to 2

    // IMPORTANT: Reset input events at the end all calculations!
    mela->resetInputEvent();

    /********************/
    /********************/
    /***** END MELA *****/
    /********************/
    /********************/

    //----------------------------------------------------------------------

    //----------------------------------------------------------------------
    //--- kinematic refitting using Z mass constraint

    float ZZMassRefit = 0.;
    float ZZMassRefitErr = 0.;
    float ZZMassUnrefitErr = 0.;

    if(doKinFit){

      vector<reco::Candidate *> selectedLeptons;
      std::map<unsigned int, TLorentzVector> selectedFsrMap;

      for(unsigned ilep=0; ilep<4; ilep++){

    selectedLeptons.push_back((reco::Candidate*)(ZZLeps[ilep]->masterClone().get()));

    if(FSRMap.find(ZZLeps[ilep])!=FSRMap.end()){
      pat::PFParticle fsr = *(FSRMap[ZZLeps[ilep]]);
      TLorentzVector p4;
      p4.SetPxPyPzE(fsr.px(),fsr.py(),fsr.pz(),fsr.energy());
      selectedFsrMap[ilep] = p4;
    }

      }

      kinZfitter->Setup(selectedLeptons, selectedFsrMap);
      kinZfitter->KinRefitZ();

      ZZMassRefit = kinZfitter->GetRefitM4l();
      ZZMassRefitErr = kinZfitter->GetRefitM4lErrFullCov();
      ZZMassUnrefitErr = kinZfitter->GetM4lErr();

      // four 4-vectors after refitting order by Z1_1,Z1_2,Z2_1,Z2_2
      //vector<TLorentzVector> p4 = kinZfitter->GetRefitP4s();

    }


    //----------------------------------------------------------------------
    //--- 4l vertex fits (experimental)

    //CandConstraintFit::fit(&myCand, iSetup);
    if (doVtxFit && abs(id11)==13 && abs(id12)==13 && abs(id21)==13 && abs(id22)==13) {

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
    }

    //----------------------------------------------------------------------
    //--- Embed variables
    myCand.addUserFloat("candChannel",    candChannel);
    myCand.addUserFloat("SIP4",           SIP4);
    myCand.addUserFloat("pt1",            ptS[3]); // leading-pT
    myCand.addUserFloat("pt2",            ptS[2]); // sub-leading pT
    myCand.addUserFloat("mZa",            mZa);
    myCand.addUserFloat("mZb",            mZb);
    myCand.addUserFloat("ZaID",           ZaID);
    myCand.addUserFloat("ZbID",           ZbID);
    myCand.addUserFloat("mZalpha",        mZalpha);
    myCand.addUserFloat("mZbeta",         mZbeta);
    myCand.addUserFloat("ZalphaID",       ZalphaID);
    myCand.addUserFloat("ZbetaID",        ZbetaID);
    myCand.addUserFloat("mLL4",           mll4); // smallest mass of any AF/OS pair
    myCand.addUserFloat("mLL6",           mll6);   // smallest mass of any AF/AS pair
    myCand.addUserFloat("passSmartMLL",   passSmartMLL);
    myCand.addUserFloat("costheta1",      costheta1);
    myCand.addUserFloat("costheta2",      costheta2);
    myCand.addUserFloat("phi",            phi);
    myCand.addUserFloat("costhetastar",   costhetastar);
    myCand.addUserFloat("phistar1",       phistar1);
    myCand.addUserFloat("xistar",         xistar);  //azimuthal angle of higgs in rest frame of higgs
    myCand.addUserFloat("xi",             xi);      //azimuthal angle of higgs in lab frame

    myCand.addUserFloat("m4l",            (Z1Lm->p4()+Z1Lp->p4()+Z2Lm->p4()+Z2Lp->p4()).mass()); // mass without FSR
    if(doKinFit) {
      myCand.addUserFloat("ZZMassRefit"   , ZZMassRefit);
      myCand.addUserFloat("ZZMassRefitErr", ZZMassRefitErr);
      myCand.addUserFloat("ZZMassUnrefitErr", ZZMassUnrefitErr);
    }

    // Jet quantities
    myCand.addUserFloat("DiJetMass", DiJetMass);
    myCand.addUserFloat("DiJetDEta", DiJetDEta);
    myCand.addUserFloat("DiJetFisher", DiJetFisher);

    // MELA probabilities
    // JHUGen
    myCand.addUserFloat("p0plus_VAJHU",   p0plus_VAJHU);
    myCand.addUserFloat("p0_g1prime2_VAJHU", p0_g1prime2_VAJHU);
    myCand.addUserFloat("p0hplus_VAJHU", p0hplus_VAJHU);
    myCand.addUserFloat("p0minus_VAJHU", p0minus_VAJHU);

    myCand.addUserFloat("p0_g1prime2_zgs_VAJHU", p0_g1prime2_zgs_VAJHU);
    myCand.addUserFloat("p0hplus_zgs_VAJHU", p0hplus_zgs_VAJHU);
    myCand.addUserFloat("p0minus_zgs_VAJHU", p0minus_zgs_VAJHU);
    myCand.addUserFloat("p0hplus_gsgs_VAJHU", p0hplus_gsgs_VAJHU);
    myCand.addUserFloat("p0minus_gsgs_VAJHU", p0minus_gsgs_VAJHU);

    myCand.addUserFloat("pg1g1prime2_VAJHU", pg1g1prime2_VAJHU);
    myCand.addUserFloat("pg1g2_VAJHU", pg1g2_VAJHU);
    myCand.addUserFloat("pg1g2_pi2_VAJHU", pg1g2_pi2_VAJHU);
    myCand.addUserFloat("pg1g4_VAJHU", pg1g4_VAJHU);
    myCand.addUserFloat("pg1g4_pi2_VAJHU", pg1g4_pi2_VAJHU);

    myCand.addUserFloat("p0plus_zz_g1prime2_zgs_VAJHU", p0plus_zz_g1prime2_zgs_VAJHU);
    myCand.addUserFloat("p0plus_zz_g1prime2_zgs_pi2_VAJHU", p0plus_zz_g1prime2_zgs_pi2_VAJHU);
    myCand.addUserFloat("p0plus_zz_0hplus_zgs_VAJHU", p0plus_zz_0hplus_zgs_VAJHU);
    myCand.addUserFloat("p0plus_zz_0minus_zgs_VAJHU", p0plus_zz_0minus_zgs_VAJHU);
    myCand.addUserFloat("p0plus_zz_0hplus_gsgs_VAJHU", p0plus_zz_0hplus_gsgs_VAJHU);
    myCand.addUserFloat("p0plus_zz_0minus_gsgs_VAJHU", p0plus_zz_0minus_gsgs_VAJHU);

    myCand.addUserFloat("p1_VAJHU", p1_VAJHU);
    myCand.addUserFloat("p1_prodIndep_VAJHU", p1_prodIndep_VAJHU);
    myCand.addUserFloat("p1plus_VAJHU", p1plus_VAJHU);
    myCand.addUserFloat("p1plus_prodIndep_VAJHU", p1plus_prodIndep_VAJHU);

    myCand.addUserFloat("p2plus_gg_VAJHU", p2plus_gg_VAJHU); // 2+mg1g5
    myCand.addUserFloat("p2plus_qqb_VAJHU", p2plus_qqb_VAJHU); // 2+mg1g5
    myCand.addUserFloat("p2plus_prodIndep_VAJHU", p2plus_prodIndep_VAJHU); // 2+mg1g5
    myCand.addUserFloat("p2h2plus_gg_VAJHU", p2h2plus_gg_VAJHU); // 2+h2
    myCand.addUserFloat("p2h2plus_qqb_VAJHU", p2h2plus_qqb_VAJHU); // 2+h2
    myCand.addUserFloat("p2h2plus_prodIndep_VAJHU", p2h2plus_prodIndep_VAJHU); // 2+h2
    myCand.addUserFloat("p2h3plus_gg_VAJHU", p2h3plus_gg_VAJHU); // 2+h3
    myCand.addUserFloat("p2h3plus_qqb_VAJHU", p2h3plus_qqb_VAJHU); // 2+h3
    myCand.addUserFloat("p2h3plus_prodIndep_VAJHU", p2h3plus_prodIndep_VAJHU); // 2+h3
    myCand.addUserFloat("p2h4plus_gg_VAJHU", p2h4plus_gg_VAJHU); // 2+h4
    myCand.addUserFloat("p2h4plus_qqb_VAJHU", p2h4plus_qqb_VAJHU); // 2+h4
    myCand.addUserFloat("p2h4plus_prodIndep_VAJHU", p2h4plus_prodIndep_VAJHU); // 2+h4
    myCand.addUserFloat("p2bplus_gg_VAJHU", p2bplus_gg_VAJHU); // 2+h5
    myCand.addUserFloat("p2bplus_qqb_VAJHU", p2bplus_qqb_VAJHU); // 2+h5
    myCand.addUserFloat("p2bplus_prodIndep_VAJHU", p2bplus_prodIndep_VAJHU); // 2+h5
    myCand.addUserFloat("p2h6plus_gg_VAJHU", p2h6plus_gg_VAJHU); // 2+h6
    myCand.addUserFloat("p2h6plus_qqb_VAJHU", p2h6plus_qqb_VAJHU); // 2+h6
    myCand.addUserFloat("p2h6plus_prodIndep_VAJHU", p2h6plus_prodIndep_VAJHU); // 2+h6
    myCand.addUserFloat("p2h7plus_gg_VAJHU", p2h7plus_gg_VAJHU); // 2+h7
    myCand.addUserFloat("p2h7plus_qqb_VAJHU", p2h7plus_qqb_VAJHU); // 2+h7
    myCand.addUserFloat("p2h7plus_prodIndep_VAJHU", p2h7plus_prodIndep_VAJHU); // 2+h7
    myCand.addUserFloat("p2hminus_gg_VAJHU", p2hminus_gg_VAJHU); // 2-h8
    myCand.addUserFloat("p2hminus_qqb_VAJHU", p2hminus_qqb_VAJHU); // 2-h8
    myCand.addUserFloat("p2hminus_prodIndep_VAJHU", p2hminus_prodIndep_VAJHU); // 2-h8
    myCand.addUserFloat("p2h9minus_gg_VAJHU", p2h9minus_gg_VAJHU); // 2-h9
    myCand.addUserFloat("p2h9minus_qqb_VAJHU", p2h9minus_qqb_VAJHU); // 2-h9
    myCand.addUserFloat("p2h9minus_prodIndep_VAJHU", p2h9minus_prodIndep_VAJHU); // 2-h9
    myCand.addUserFloat("p2h10minus_gg_VAJHU", p2h10minus_gg_VAJHU); // 2-h10
    myCand.addUserFloat("p2h10minus_qqb_VAJHU", p2h10minus_qqb_VAJHU); // 2-h10
    myCand.addUserFloat("p2h10minus_prodIndep_VAJHU", p2h10minus_prodIndep_VAJHU); // 2-h10

    // MCFM
    myCand.addUserFloat("p0plus_VAMCFM", p0plus_VAMCFM);
    myCand.addUserFloat("ggzz_VAMCFM", ggzz_VAMCFM);
    myCand.addUserFloat("ggzz_p0plus_VAMCFM", ggzz_p0plus_VAMCFM);
    myCand.addUserFloat("bkg_VAMCFM", bkg_VAMCFM);
    myCand.addUserFloat("bkg_prodIndep_VAMCFM", bkg_prodIndep_VAMCFM);
    myCand.addUserFloat("pZJJ_VAMCFM", pZJJ_VAMCFM);
    myCand.addUserFloat("Dgg10_VAMCFM", Dgg10_VAMCFM);

    // SuperMELA
    myCand.addUserFloat("p0plus_m4l", p0plus_m4l);
    myCand.addUserFloat("p0plus_m4l_ScaleUp", p0plus_m4l_ScaleUp);
    myCand.addUserFloat("p0plus_m4l_ScaleDown", p0plus_m4l_ScaleDown);
    myCand.addUserFloat("p0plus_m4l_ResUp", p0plus_m4l_ResUp);
    myCand.addUserFloat("p0plus_m4l_ResDown", p0plus_m4l_ResDown);
    myCand.addUserFloat("bkg_m4l", bkg_m4l);
    myCand.addUserFloat("bkg_m4l_ScaleUp", bkg_m4l_ScaleUp);
    myCand.addUserFloat("bkg_m4l_ScaleDown", bkg_m4l_ScaleDown);
    myCand.addUserFloat("bkg_m4l_ResUp", bkg_m4l_ResUp);
    myCand.addUserFloat("bkg_m4l_ResDown", bkg_m4l_ResDown);

    // JHUGen production
    myCand.addUserFloat("pwh_leptonic_VAJHU", pwh_leptonic_VAJHU);
    myCand.addUserFloat("pzh_leptonic_VAJHU", pzh_leptonic_VAJHU);

    myCand.addUserFloat("phjj_VAJHU_highestPTJets", phjj_VAJHU_highestPTJets);
    myCand.addUserFloat("pvbf_VAJHU_highestPTJets", pvbf_VAJHU_highestPTJets);
    myCand.addUserFloat("phjj_VAJHU_highestPTJets_up", phjj_VAJHU_highestPTJets_up);
    myCand.addUserFloat("pvbf_VAJHU_highestPTJets_up", pvbf_VAJHU_highestPTJets_up);
    myCand.addUserFloat("phjj_VAJHU_highestPTJets_dn", phjj_VAJHU_highestPTJets_dn);
    myCand.addUserFloat("pvbf_VAJHU_highestPTJets_dn", pvbf_VAJHU_highestPTJets_dn);
    myCand.addUserFloat("phjj_VAJHU_bestDjet", phjj_VAJHU_bestDjet);
    myCand.addUserFloat("pvbf_VAJHU_bestDjet", pvbf_VAJHU_bestDjet);
    myCand.addUserFloat("phjj_VAJHU_bestDjet_up", phjj_VAJHU_bestDjet_up);
    myCand.addUserFloat("pvbf_VAJHU_bestDjet_up", pvbf_VAJHU_bestDjet_up);
    myCand.addUserFloat("phjj_VAJHU_bestDjet_dn", phjj_VAJHU_bestDjet_dn);
    myCand.addUserFloat("pvbf_VAJHU_bestDjet_dn", pvbf_VAJHU_bestDjet_dn);

    myCand.addUserFloat("pAux_vbf_VAJHU", pAux_vbf_VAJHU);
    myCand.addUserFloat("pAux_vbf_VAJHU_up", pAux_vbf_VAJHU_up);
    myCand.addUserFloat("pAux_vbf_VAJHU_dn", pAux_vbf_VAJHU_dn);

    myCand.addUserFloat("phj_VAJHU", phj_VAJHU);
    myCand.addUserFloat("phj_VAJHU_up", phj_VAJHU_up);
    myCand.addUserFloat("phj_VAJHU_dn", phj_VAJHU_dn);

    myCand.addUserFloat("pwh_hadronic_VAJHU", pwh_hadronic_VAJHU);
    myCand.addUserFloat("pwh_hadronic_VAJHU_up", pwh_hadronic_VAJHU_up);
    myCand.addUserFloat("pwh_hadronic_VAJHU_dn", pwh_hadronic_VAJHU_dn);

    myCand.addUserFloat("pzh_hadronic_VAJHU", pzh_hadronic_VAJHU);
    myCand.addUserFloat("pzh_hadronic_VAJHU_up", pzh_hadronic_VAJHU_up);
    myCand.addUserFloat("pzh_hadronic_VAJHU_dn", pzh_hadronic_VAJHU_dn);

    myCand.addUserFloat("ptth_VAJHU", ptth_VAJHU);
    myCand.addUserFloat("ptth_VAJHU_up", ptth_VAJHU_up);
    myCand.addUserFloat("ptth_VAJHU_dn", ptth_VAJHU_dn);

    myCand.addUserFloat("pbbh_VAJHU", pbbh_VAJHU);
    myCand.addUserFloat("pbbh_VAJHU_up", pbbh_VAJHU_up);
    myCand.addUserFloat("pbbh_VAJHU_dn", pbbh_VAJHU_dn);


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

      if (preBestCandResult){
        // Fill preSelCands matrix
        preSelCands[iCRname].push_back(icand);
      }
      iCRname++;
    }
    result->push_back(myCand);

  } // End of loop over input candidates


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
ZZCandidateFiller::getPairMass(const reco::Candidate* lp, const reco::Candidate* lm, ZZCandidateFiller::FSRToLepMap& photons, float& mass, int& ID) {
  math::XYZTLorentzVector llp4 = lp->p4()+lm->p4();
  auto lpp = photons.find(lp);
  auto lmp = photons.find(lm);
  if (lpp!=photons.end()) llp4+=lpp->second->p4();
  if (lmp!=photons.end()) llp4+=lmp->second->p4();
  mass=llp4.mass();
  ID=lp->pdgId()*lm->pdgId();
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ZZCandidateFiller);

