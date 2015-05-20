/** \class ZZCandidateFiller
 *
 *
 *  $Date: 2013/12/12 16:19:40 $
 *  $Revision: 1.72 $
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
#include <ZZMatrixElement/MELA/src/computeAngles.h>
#include <ZZMatrixElement/MEMCalculators/interface/MEMCalculators.h>
//#include <ZZAnalysis/AnalysisStep/interface/ZZMassErrors.h>
//#include <ZZAnalysis/AnalysisStep/interface/MCHistoryTools.h>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>
#include <ZZAnalysis/AnalysisStep/interface/CompositeCandMassResolution.h>
#include <ZZAnalysis/AnalysisStep/interface/DiLeptonKinFitter.h>
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
#include <DataFormats/METReco/interface/PFMET.h>
#include <ZZAnalysis/AnalysisStep/interface/Fisher.h>
#include <ZZAnalysis/AnalysisStep/interface/Comparators.h>
#include <ZZAnalysis/AnalysisStep/interface/utils.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>


#include "TH2F.h"
#include "TFile.h"

#include <string>

using namespace zzanalysis;

bool doKinFit = false;
bool doVtxFit = false;
bool doMEKD = true;

class ZZCandidateFiller : public edm::EDProducer {
public:
  /// Constructor
  explicit ZZCandidateFiller(const edm::ParameterSet&);
    
  /// Destructor
  ~ZZCandidateFiller(){};  

private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  void getPairMass(const reco::Candidate* lp, const reco::Candidate* lm, map<const reco::Candidate*, math::XYZTLorentzVector>& photons, float& mass, int& ID);

  edm::InputTag theCandidateTag;
  const CutSet<pat::CompositeCandidate> preBestCandSelection;
  const CutSet<pat::CompositeCandidate> cuts;
  int sampleType;
  int setup;
  float superMelaMass;
  MEMs combinedMEM;
  Mela* myMela;
  bool embedDaughterFloats;
  bool ZRolesByMass;
  reco::CompositeCandidate::role_collection rolesZ1Z2;  
  reco::CompositeCandidate::role_collection rolesZ2Z1;
  bool isMC;
  TH2F* corrSigmaMu;
  TH2F* corrSigmaEle;
  Comparators::ComparatorTypes bestCandType;
};


static int SetupToSqrts(int setup) {
  if (setup==2011) return 7;
  else if (setup==2012) return 8;
  else if (setup==2015) return 13;
  else return 0;
}


ZZCandidateFiller::ZZCandidateFiller(const edm::ParameterSet& iConfig) :
  theCandidateTag(iConfig.getParameter<edm::InputTag>("src")),
  preBestCandSelection(iConfig.getParameter<edm::ParameterSet>("bestCandAmong")),
  cuts(iConfig.getParameter<edm::ParameterSet>("flags")),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  superMelaMass(iConfig.getParameter<double>("superMelaMass")),
  combinedMEM(SetupToSqrts(setup),superMelaMass,"CTEQ6L"),
  embedDaughterFloats(iConfig.getUntrackedParameter<bool>("embedDaughterFloats",true)),
  ZRolesByMass(iConfig.getParameter<bool>("ZRolesByMass")),
  isMC(iConfig.getParameter<bool>("isMC")),
  corrSigmaMu(0),
  corrSigmaEle(0)
{
  produces<pat::CompositeCandidateCollection>();
  
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
  else if (cmp=="byBestKD")           bestCandType=Comparators::byBestKD;
  else if (cmp=="byBestKD_VH")           bestCandType=Comparators::byBestKD_VH;
  else abort();

  //-- Non-MEM discriminants
  myMela = combinedMEM.m_MELA;
}


void ZZCandidateFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){  
  using namespace edm;
  using namespace std;
  using namespace reco;

  std::auto_ptr<pat::CompositeCandidateCollection> result(new pat::CompositeCandidateCollection);

  const float ZmassValue = 91.1876;

  double rhoForMu, rhoForEle;
  {
    edm::Handle<double> rhoHandle;
    iEvent.getByLabel(LeptonIsoHelper::getMuRhoTag(sampleType, setup), rhoHandle);
    rhoForMu = *rhoHandle;
    iEvent.getByLabel(LeptonIsoHelper::getEleRhoTag(sampleType, setup), rhoHandle);
    rhoForEle = *rhoHandle;
  }

  // Get LLLL candidates
  Handle<View<CompositeCandidate> > LLLLCands;
  iEvent.getByLabel(theCandidateTag, LLLLCands);

  // Get jets
  Handle<edm::View<pat::Jet> > CleanJets;
  iEvent.getByLabel("cleanJets", CleanJets);
  vector<const pat::Jet*> cleanedJetsPt30;
  for(edm::View<pat::Jet>::const_iterator jet = CleanJets->begin(); jet != CleanJets->end(); ++jet){
    if(jet->pt()>30) cleanedJetsPt30.push_back(&*jet);
  }

  // Get MET
  Handle<vector<reco::MET> > pfmetcoll;
  iEvent.getByLabel("slimmedMETs", pfmetcoll);
  math::XYZTLorentzVector pfmet;
  if(pfmetcoll.isValid()) pfmet = pfmetcoll->front().p4(); // standard MET is pfmet.pt();

  // Get leptons (in order to store extra leptons)
  Handle<View<reco::Candidate> > softleptoncoll;
  iEvent.getByLabel("softLeptons", softleptoncoll);
  vector<reco::CandidatePtr> goodisoleptonPtrs;
  for( View<reco::Candidate>::const_iterator lep = softleptoncoll->begin(); lep != softleptoncoll->end(); ++ lep ){ 
    if((bool)userdatahelpers::getUserFloat(&*lep,"isGood") && (bool)userdatahelpers::getUserFloat(&*lep,"isIsoFSRUncorr")){
      const reco::CandidatePtr lepPtr(softleptoncoll,lep-softleptoncoll->begin());
      goodisoleptonPtrs.push_back(lepPtr);
    }
  }

  // Get Z Candidates
  Handle<View<CompositeCandidate> > ZCands;
  iEvent.getByLabel("ZCand", ZCands);

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
    //--- Lepton pointers in the original order
    const reco::Candidate* Z1L1= Z1->daughter(0);
    const reco::Candidate* Z1L2= Z1->daughter(1);
    const reco::Candidate* Z2L1= Z2->daughter(0);
    const reco::Candidate* Z2L2= Z2->daughter(1);

    //--- Lepton four-vectors in the original order; with FSR added (below)
    math::XYZTLorentzVector p11 = Z1L1->p4();
    math::XYZTLorentzVector p12 = Z1L2->p4();
    math::XYZTLorentzVector p21 = Z2L1->p4();
    math::XYZTLorentzVector p22 = Z2L2->p4();
    int id11 = Z1L1->pdgId();
    int id12 = Z1L2->pdgId();
    int id21 = Z2L1->pdgId();
    int id22 = Z2L2->pdgId();
    int candChannel = id11*id12*id21*id22;


    //--- Pick FSR photons; add their fourmomentum to p11...p22
    vector<math::XYZTLorentzVector> photons;
    int d0FSR = (static_cast<const pat::CompositeCandidate*>(Z1->masterClone().get()))->userFloat("dauWithFSR");
    if (d0FSR>=0) {
      if (Z1->numberOfDaughters()!=3) cout << "ERROR: ZZCandidateFiller: problem in FSR" << endl;
      const math::XYZTLorentzVector& fsr = Z1->daughter(2)->p4();
      photons.push_back(fsr);
      if (d0FSR==0) {	
        p11 = p11 + fsr;
      } else if (d0FSR==1){
        p12 = p12 + fsr;
      }
    }
    int d1FSR = (static_cast<const pat::CompositeCandidate*>(Z2->masterClone().get()))->userFloat("dauWithFSR");
    if (d1FSR>=0) {
      if (Z2->numberOfDaughters()!=3) cout << "ERROR: ZZCandidateFiller: problem in FSR" << endl;
      const math::XYZTLorentzVector& fsr = Z2->daughter(2)->p4();
      photons.push_back(fsr);
      if (d1FSR==0) {
        p21 = p21 + fsr;
      } else if (d1FSR==1){
        p22 = p22 + fsr;
      }
    }

    // Recompute isolation for all four leptons if FSR is present
    for (int zIdx=0; zIdx<2; ++zIdx) {
      float worstMuIso=0;
      float worstEleIso=0;
      for (int dauIdx=0; dauIdx<2; ++dauIdx) { 
	const reco::Candidate* z = myCand.daughter(zIdx);
	const reco::Candidate* d = z->daughter(dauIdx);
	float fsrCorr = 0; // The correction to PFPhotonIso
	for (vector<math::XYZTLorentzVector>::const_iterator ifsr=photons.begin(); ifsr!=photons.end(); ++ifsr) {
	  double dR = ROOT::Math::VectorUtil::DeltaR(*ifsr,d->momentum());
	  // Check if the photon is in the lepton's iso cone and not vetoed
	  if (dR<0.4 && ((d->isMuon() && dR > 0.01) ||
			 (d->isElectron() && (fabs((static_cast<const pat::Electron*>(d->masterClone().get()))->superCluster()->eta()) < 1.479 || dR > 0.08)))) {
	    fsrCorr += ifsr->pt();
	  }
	}
	float rho = ((d->isMuon())?rhoForMu:rhoForEle);
	float combRelIsoPFCorr =  LeptonIsoHelper::combRelIsoPF(sampleType, setup, rho, d, fsrCorr);
	if (d->isMuon()) worstMuIso  = max(worstMuIso,  combRelIsoPFCorr);
	else             worstEleIso = max(worstEleIso, combRelIsoPFCorr);
	string base;
	stringstream str;
	str << "d" << zIdx << "." << "d" << dauIdx << ".";
	str >> base;
	myCand.addUserFloat(base+"combRelIsoPFFSRCorr",combRelIsoPFCorr);
	myCand.addUserFloat(base+"passCombRelIsoPFFSRCorr",combRelIsoPFCorr < (d->isMuon()?0.4:0.5)); // FIXME: not the most elegant solution; hard coded right now to see how things evolve about lepton isolation requirements. 
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
    vector<const reco::Candidate*> ZLeps = {Z1L1,Z1L2,Z2L1,Z2L2};
    map<const reco::Candidate*, math::XYZTLorentzVector> FSR;
    if (d0FSR>=0) FSR[ZLeps[d0FSR]]   = Z1->daughter(2)->p4();
    if (d1FSR>=0) FSR[ZLeps[d1FSR+2]] = Z2->daughter(2)->p4();

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
    getPairMass(Z1Lp,Z2Lm,FSR,mZa,ZaID);
    getPairMass(Z1Lm,Z2Lp,FSR,mZb,ZbID);

    // For same-sign CRs, the Z2 leptons are same sign, so we need to check also the other combination. 
    float mZalpha, mZbeta;
    int ZalphaID, ZbetaID;
    getPairMass(Z1Lp,Z2Lp,FSR,mZalpha,ZalphaID);
    getPairMass(Z1Lm,Z2Lm,FSR,mZbeta,ZbetaID);

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
    for( vector<reco::CandidatePtr>::const_iterator lepPtr = goodisoleptonPtrs.begin(); lepPtr != goodisoleptonPtrs.end(); ++ lepPtr ) {
      const reco::Candidate* lep = lepPtr->get();
      if( reco::deltaR( lep->p4(), Z1L1->p4() ) > 0.02 &&
	  reco::deltaR( lep->p4(), Z1L2->p4() ) > 0.02 &&
	  reco::deltaR( lep->p4(), Z2L1->p4() ) > 0.02 &&
	  reco::deltaR( lep->p4(), Z2L2->p4() ) > 0.02 ){
	nExtraLep++;
	myCand.addUserCand("ExtraLep"+to_string(nExtraLep),*lepPtr);
      }
    }
    myCand.addUserFloat("nExtraLep",nExtraLep);

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
	if((bool)userdatahelpers::getUserFloat(&*myZCand,"GoodLeptons")){
	  nExtraZ++;
	  extraZs.push_back(&*zcand);
	  myCand.addUserCand("assocZ"+to_string(nExtraZ),myZCand);
	}
      }
    }
    myCand.addUserFloat("nExtraZ",nExtraZ);



    //----------------------------------------------------------------------
    //--- Embed angular information and probabilities to build discriminants

    // Lepton TLorentzVectors, including FSR
    TLorentzVector pL11(p11.x(),p11.y(),p11.z(),p11.t());
    TLorentzVector pL12(p12.x(),p12.y(),p12.z(),p12.t());
    TLorentzVector pL21(p21.x(),p21.y(),p21.z(),p21.t());
    TLorentzVector pL22(p22.x(),p22.y(),p22.z(),p22.t());

    std::vector<TLorentzVector> partP;
    partP.push_back(pL11);
    partP.push_back(pL12);
    partP.push_back(pL21);
    partP.push_back(pL22);

    std::vector<int> partId; 
    partId.push_back(id11);
    partId.push_back(id12);
    partId.push_back(id21);
    partId.push_back(id22);

    //--- compute angles
    float costheta1=0, costheta2=0, phi=0, costhetastar=0, phistar1=0;
    mela::computeAngles(pL11,id11,pL12,id12,pL21,id21,pL22,id22,costhetastar,costheta1,costheta2,phi,phistar1);
  
    //--- compute higgs azimuthal angles, xi
    TLorentzVector higgs = pL11+pL12+pL21+pL22;
    TLorentzVector Z14vec = pL11+pL12;
    TVector3 Xaxis(1,0,0);
    float xi = higgs.Phi();
    // boost Z1 into rest frame of higgs
    // xistar is the angle between Z decay plane and x-axis 
    Z14vec.Boost(-higgs.BoostVector());
    float xistar = Z14vec.Phi();
    // higgs.Boost(-higgs.BoostVector());


    double mekd_ld=-1,mekd_pseudold=-1,mekd_gravld=-1,mekd_A=0,mekd_B=0,mekd_AP=0,mekd_BP=0,mekd_AG=0,mekd_BG=0;
    if(doMEKD){ //FIXME: do we still need this?
      int prodid=id11*id21*id12*id22;
      if(prodid==28561 || prodid==14641 || prodid==20449){
        combinedMEM.computeKD(MEMNames::kSMHiggs,MEMNames::kqqZZ,MEMNames::kMEKD,partP,partId, mekd_ld,mekd_A,mekd_B);
        combinedMEM.computeKD(MEMNames::k0minus,MEMNames::kSMHiggs,MEMNames::kMEKD,partP,partId,mekd_pseudold,mekd_AP,mekd_BP);
        combinedMEM.computeKD(MEMNames::k2mplus_gg,MEMNames::kSMHiggs,MEMNames::kMEKD,partP,partId,mekd_gravld,mekd_AG,mekd_BG);
      }
    }  

    // probabilities
    double p0plus_VAJHU,p0minus_VAJHU,p0plus_VAMCFM,p0hplus_VAJHU;
    double p1_VAJHU,p1plus_VAJHU,p2_VAJHU,p2qqb_VAJHU; 
    double bkg_VAMCFM,bkg_prodIndep_VAMCFM;
    double ggzz_VAMCFM,ggzz_c1_VAMCFM,ggzz_c5_VAMCFM,ggzz_ci_VAMCFM;
    double ggzz_p0plus_VAMCFM;

    // exotic spin-2 models
    double p2hplus_VAJHU, p2hminus_VAJHU,p2bplus_VAJHU;
    double p2hplus_qqb_VAJHU,p2hplus_prodIndep_VAJHU,p2hminus_qqb_VAJHU,p2hminus_prodIndep_VAJHU,p2bplus_qqb_VAJHU,p2bplus_prodIndep_VAJHU;
    double p2h2plus_gg_VAJHU,p2h2plus_qqbar_VAJHU,p2h2plus_prodIndep_VAJHU,p2h3plus_gg_VAJHU,p2h3plus_qqbar_VAJHU,p2h3plus_prodIndep_VAJHU;
    double p2h6plus_gg_VAJHU,p2h6plus_qqbar_VAJHU,p2h6plus_prodIndep_VAJHU,p2h7plus_gg_VAJHU,p2h7plus_qqbar_VAJHU,p2h7plus_prodIndep_VAJHU;
    double p2h9minus_gg_VAJHU,p2h9minus_qqbar_VAJHU,p2h9minus_prodIndep_VAJHU,p2h10minus_gg_VAJHU,p2h10minus_qqbar_VAJHU,p2h10minus_prodIndep_VAJHU;

    // production independent probabilites
    double p1_prodIndep_VAJHU,p1plus_prodIndep_VAJHU,p2_prodIndep_VAJHU;
    double p0plus_m4l,bkg_m4l; //supermela
    double p0plus_m4l_ScaleUp,bkg_m4l_ScaleUp,p0plus_m4l_ScaleDown,bkg_m4l_ScaleDown,p0plus_m4l_ResUp,bkg_m4l_ResUp,p0plus_m4l_ResDown,bkg_m4l_ResDown; // supermela uncertainties

    double p0_g1prime2_VAJHU, Dgg10_VAMCFM, pzzzg_VAJHU, pzzgg_VAJHU, p0Zgs_VAJHU, p0gsgs_VAJHU, pzzzg_PS_VAJHU, pzzgg_PS_VAJHU, p0Zgs_PS_VAJHU, p0gsgs_PS_VAJHU, p0Zgs_g1prime2_VAJHU;

    combinedMEM.computeME(MEMNames::kSMHiggs     ,MEMNames::kJHUGen    ,partP,partId,p0plus_VAJHU);   // higgs, vector algebra, JHUgen
    combinedMEM.computeME(MEMNames::k0minus      ,MEMNames::kJHUGen    ,partP,partId,p0minus_VAJHU);  // pseudoscalar, vector algebra, JHUgen
    combinedMEM.computeME(MEMNames::kSMHiggs     ,MEMNames::kMCFM      ,partP,partId,p0plus_VAMCFM);  // higgs, vector algebra, MCFM
    combinedMEM.computeME(MEMNames::k0hplus      ,MEMNames::kJHUGen    ,partP,partId,p0hplus_VAJHU);  // 0h+ (high dimensional operator), vector algebra, JHUgen
    combinedMEM.computeME(MEMNames::k1minus      ,MEMNames::kJHUGen    ,partP,partId,p1_VAJHU);       // zprime, vector algebra, JHUgen
    combinedMEM.computeME(MEMNames::k1minus_prodIndep      ,MEMNames::kJHUGen    ,partP,partId,p1_prodIndep_VAJHU);           // zprime, vector algebra, JHUgen
    combinedMEM.computeME(MEMNames::k1plus       ,MEMNames::kJHUGen    ,partP,partId,p1plus_VAJHU);   // 1+ (axial vector), vector algebra, JHUgen
    combinedMEM.computeME(MEMNames::k1plus_prodIndep       ,MEMNames::kJHUGen    ,partP,partId,p1plus_prodIndep_VAJHU);   // 1+ (axial vector), vector algebra, JHUgen

    combinedMEM.computeME(MEMNames::k2mplus_gg   ,MEMNames::kJHUGen    ,partP,partId,p2_VAJHU);       // graviton, vector algebra, JHUgen
    combinedMEM.computeME(MEMNames::k2mplus_prodIndep   ,MEMNames::kJHUGen    ,partP,partId,p2_prodIndep_VAJHU);              // graviton, vector algebra, JHUgen
    combinedMEM.computeME(MEMNames::k2mplus_qqbar,MEMNames::kJHUGen    ,partP,partId,p2qqb_VAJHU);    // graviton produced by qqbar vector algebra, JHUgen
    combinedMEM.computeME(MEMNames::k2hplus,      MEMNames::kJHUGen    ,partP,partId,p2hplus_VAJHU);    // graviton produced by qqbar vector algebra, JHUgen
    combinedMEM.computeME(MEMNames::k2hminus,     MEMNames::kJHUGen    ,partP,partId,p2hminus_VAJHU);    // graviton produced by qqbar vector algebra, JHUgen
    combinedMEM.computeME(MEMNames::k2bplus,      MEMNames::kJHUGen    ,partP,partId,p2bplus_VAJHU);    // graviton produced by qqbar vector algebra, JHUgen

    combinedMEM.computeME(MEMNames::k2hplus_qqbar,MEMNames::kJHUGen    ,partP,partId,           p2hplus_qqb_VAJHU );    // graviton produced by qqbar vector algebra, JHUgen
    combinedMEM.computeME(MEMNames::k2hplus_prodIndep,MEMNames::kJHUGen  ,partP,partId, p2hplus_prodIndep_VAJHU   );    // graviton produced by qqbar vector algebra, JHUgen
    combinedMEM.computeME(MEMNames::k2hminus_qqbar,MEMNames::kJHUGen    ,partP,partId,  p2hminus_qqb_VAJHU        );    // graviton produced by qqbar vector algebra, JHUgen
    combinedMEM.computeME(MEMNames::k2hminus_prodIndep,MEMNames::kJHUGen ,partP,partId, p2hminus_prodIndep_VAJHU  );    // graviton produced by qqbar vector algebra, JHUgen
    combinedMEM.computeME(MEMNames::k2bplus_qqbar,MEMNames::kJHUGen    ,partP,partId,           p2bplus_qqb_VAJHU );    // graviton produced by qqbar vector algebra, JHUgen
    combinedMEM.computeME(MEMNames::k2bplus_prodIndep,MEMNames::kJHUGen  ,partP,partId, p2bplus_prodIndep_VAJHU   );    // graviton produced by qqbar vector algebra, JHUgen

    combinedMEM.computeME(MEMNames::k2h2plus_gg         ,MEMNames::kJHUGen,partP,partId,p2h2plus_gg_VAJHU         );
    combinedMEM.computeME(MEMNames::k2h2plus_qqbar      ,MEMNames::kJHUGen,partP,partId,p2h2plus_qqbar_VAJHU      );
    combinedMEM.computeME(MEMNames::k2h2plus_prodIndep  ,MEMNames::kJHUGen,partP,partId,p2h2plus_prodIndep_VAJHU  );
    combinedMEM.computeME(MEMNames::k2h3plus_gg         ,MEMNames::kJHUGen,partP,partId,p2h3plus_gg_VAJHU         );
    combinedMEM.computeME(MEMNames::k2h3plus_qqbar      ,MEMNames::kJHUGen,partP,partId,p2h3plus_qqbar_VAJHU      );
    combinedMEM.computeME(MEMNames::k2h3plus_prodIndep  ,MEMNames::kJHUGen,partP,partId,p2h3plus_prodIndep_VAJHU  );
    combinedMEM.computeME(MEMNames::k2h6plus_gg         ,MEMNames::kJHUGen,partP,partId,p2h6plus_gg_VAJHU         );
    combinedMEM.computeME(MEMNames::k2h6plus_qqbar      ,MEMNames::kJHUGen,partP,partId,p2h6plus_qqbar_VAJHU      );
    combinedMEM.computeME(MEMNames::k2h6plus_prodIndep  ,MEMNames::kJHUGen,partP,partId,p2h6plus_prodIndep_VAJHU  );
    combinedMEM.computeME(MEMNames::k2h7plus_gg         ,MEMNames::kJHUGen,partP,partId,p2h7plus_gg_VAJHU         );
    combinedMEM.computeME(MEMNames::k2h7plus_qqbar      ,MEMNames::kJHUGen,partP,partId,p2h7plus_qqbar_VAJHU      );
    combinedMEM.computeME(MEMNames::k2h7plus_prodIndep  ,MEMNames::kJHUGen,partP,partId,p2h7plus_prodIndep_VAJHU  );
    combinedMEM.computeME(MEMNames::k2h9minus_gg        ,MEMNames::kJHUGen,partP,partId,p2h9minus_gg_VAJHU        );
    combinedMEM.computeME(MEMNames::k2h9minus_qqbar     ,MEMNames::kJHUGen,partP,partId,p2h9minus_qqbar_VAJHU     );
    combinedMEM.computeME(MEMNames::k2h9minus_prodIndep ,MEMNames::kJHUGen,partP,partId,p2h9minus_prodIndep_VAJHU );
    combinedMEM.computeME(MEMNames::k2h10minus_gg       ,MEMNames::kJHUGen,partP,partId,p2h10minus_gg_VAJHU       );
    combinedMEM.computeME(MEMNames::k2h10minus_qqbar    ,MEMNames::kJHUGen,partP,partId,p2h10minus_qqbar_VAJHU    );
    combinedMEM.computeME(MEMNames::k2h10minus_prodIndep,MEMNames::kJHUGen,partP,partId,p2h10minus_prodIndep_VAJHU);

    combinedMEM.computeME(MEMNames::kqqZZ_prodIndep ,MEMNames::kMCFM      ,partP,partId,bkg_prodIndep_VAMCFM);     // background, vector algebra, MCFM
    combinedMEM.computeME(MEMNames::kqqZZ        ,MEMNames::kMCFM      ,partP,partId,bkg_VAMCFM);     // background, vector algebra, MCFM
    combinedMEM.computeME(MEMNames::k0_g1prime2,MEMNames::kJHUGen,partP,partId,p0_g1prime2_VAJHU);     // background, vector algebra, MCFM
    combinedMEM.computeME(MEMNames::kggHZZ_10,MEMNames::kMCFM,partP,partId,Dgg10_VAMCFM);     // background, vector algebra, MCFM
    //combinedMEM.computeME(MEMNames::kggZZ        ,MEMNames::kMCFM      ,partP,partId,ggzz_VAMCFM);    // background, vector algebra, MCFM for ggzz

    combinedMEM.computeME(MEMNames::k0_Zgs                       ,MEMNames::kJHUGen    ,partP,partId,p0Zgs_VAJHU);  // SM Higgs to Zgamma star 
    combinedMEM.computeME(MEMNames::k0_gsgs      ,MEMNames::kJHUGen    ,partP,partId,p0gsgs_VAJHU);  // SM Higgs to gamma star gamma star 
    combinedMEM.computeME(MEMNames::k0_Zgs_PS                    ,MEMNames::kJHUGen    ,partP,partId,p0Zgs_PS_VAJHU);  // SM Higgs to Zgamma star 
    combinedMEM.computeME(MEMNames::k0_gsgs_PS      ,MEMNames::kJHUGen    ,partP,partId,p0gsgs_PS_VAJHU);  // SM Higgs to gamma star gamma star
    combinedMEM.computeME(MEMNames::k0_Zgs_g1prime2, MEMNames::kJHUGen, partP, partId, p0Zgs_g1prime2_VAJHU);  // SM Higgs to Zgamma star Lambda1 

    vector<complex<double> > *coupling = new vector<complex<double> >;
    vector<complex<double> > *couplingprod= new vector<complex<double> >;
    combinedMEM.computeME(MEMNames::kggZZ, MEMNames::kMCFM, partP, partId, ggzz_VAMCFM);
    combinedMEM.computeME(MEMNames::kggZZ_SMHiggs, MEMNames::kMCFM, partP, partId, ggzz_p0plus_VAMCFM);

    complex<double> coup(1.,0.); 
    coupling->push_back(coup);
    combinedMEM.computeME(MEMNames::kggZZ_SMHiggs, MEMNames::kJHUGen, partP, partId, couplingprod, coupling, ggzz_c1_VAMCFM);
    coupling->clear();
    coup=5.;
    coupling->push_back(coup);
    combinedMEM.computeME(MEMNames::kggZZ_SMHiggs, MEMNames::kJHUGen, partP, partId, couplingprod, coupling, ggzz_c5_VAMCFM);
    coupling->clear();
    complex<double> coup2(0.,1.);
    coupling->push_back(coup2);
    combinedMEM.computeME(MEMNames::kggZZ_SMHiggs, MEMNames::kJHUGen, partP, partId, couplingprod, coupling, ggzz_ci_VAMCFM);

    // Supermela signal
    // m4l probability as in datacards
    combinedMEM.computePm4l(partP,partId, MEMNames::kNone,      p0plus_m4l,           bkg_m4l);
    // probabilities for systematics
    combinedMEM.computePm4l(partP,partId, MEMNames::kScaleUp,   p0plus_m4l_ScaleUp,   bkg_m4l_ScaleUp); 
    combinedMEM.computePm4l(partP,partId, MEMNames::kScaleDown, p0plus_m4l_ScaleDown, bkg_m4l_ScaleDown);
    combinedMEM.computePm4l(partP,partId, MEMNames::kResolUp,   p0plus_m4l_ResUp,     bkg_m4l_ResUp);
    combinedMEM.computePm4l(partP,partId, MEMNames::kResolDown, p0plus_m4l_ResDown,   bkg_m4l_ResDown);


    // spinMELA
    double pg1g4_mela, pg1g4_VAJHU, pg1g2_mela, pg1g2_VAJHU, pg1g4_pi2_VAJHU, pg1g2_pi2_VAJHU, pg1g1prime2_VAJHU, pzzzg_g1prime2_VAJHU, pzzzg_g1prime2_pi2_VAJHU;
    combinedMEM.computeME_Interference(MEMNames::kg1g4,     MEMNames::kAnalytical, partP, partId, pg1g4_mela);
    combinedMEM.computeME_Interference(MEMNames::kg1g4,     MEMNames::kJHUGen,     partP, partId, pg1g4_VAJHU);
    combinedMEM.computeME_Interference(MEMNames::kg1g4_pi_2,MEMNames::kJHUGen,     partP, partId, pg1g4_pi2_VAJHU);
    combinedMEM.computeME_Interference(MEMNames::kg1g2_pi_2,MEMNames::kJHUGen,     partP, partId, pg1g2_pi2_VAJHU);
    combinedMEM.computeME_Interference(MEMNames::kg1g2,     MEMNames::kAnalytical, partP, partId, pg1g2_mela);
    combinedMEM.computeME_Interference(MEMNames::kg1g2,     MEMNames::kJHUGen,     partP, partId, pg1g2_VAJHU);     
    combinedMEM.computeME_Interference(MEMNames::k_g1g1prime2,     MEMNames::kJHUGen,     partP, partId, pg1g1prime2_VAJHU);     

    combinedMEM.computeME_Interference(MEMNames::kzzzg,MEMNames::kJHUGen,partP, partId,pzzzg_VAJHU);
    combinedMEM.computeME_Interference(MEMNames::kzzgg,MEMNames::kJHUGen,partP, partId,pzzgg_VAJHU);
    combinedMEM.computeME_Interference(MEMNames::kzzzg_PS,MEMNames::kJHUGen,partP, partId,pzzzg_PS_VAJHU);
    combinedMEM.computeME_Interference(MEMNames::kzzgg_PS,MEMNames::kJHUGen,partP, partId,pzzgg_PS_VAJHU);
    combinedMEM.computeME_Interference(MEMNames::kzzzg_g1prime2, MEMNames::kJHUGen, partP, partId, pzzzg_g1prime2_VAJHU);
    combinedMEM.computeME_Interference(MEMNames::kzzzg_g1prime2_pi_2, MEMNames::kJHUGen, partP, partId, pzzzg_g1prime2_pi2_VAJHU);


    

    //----------------------------------------------------------------------
    //--- Add production probabilities
//    TLorentzVector nullFourVector(0, 0, 0, 0);

    double phjj_VAJHU_old = -1.;
    double pvbf_VAJHU_old = -1.;
    double phjj_VAJHU_old_up = -1.;
    double pvbf_VAJHU_old_up = -1.;
    double phjj_VAJHU_old_dn = -1.;
    double pvbf_VAJHU_old_dn = -1.;
    double phjj_VAJHU_new = -1.;
    double pvbf_VAJHU_new = -1.;
    double phjj_VAJHU_new_up = -1.;
    double pvbf_VAJHU_new_up = -1.;
    double phjj_VAJHU_new_dn = -1.;
    double pvbf_VAJHU_new_dn = -1.;

    double pAux_vbf_VAJHU = 1;
    double pAux_vbf_VAJHU_up = 1;
    double pAux_vbf_VAJHU_dn = 1;

    double phj_VAJHU = -1.;
    double phj_VAJHU_up = -1.;
    double phj_VAJHU_dn = -1.;

    double pwh_hadronic_VAJHU = -1.;
    double pwh_hadronic_VAJHU_up = -1.;
    double pwh_hadronic_VAJHU_dn = -1.;

    double pzh_hadronic_VAJHU = -1.;
    double pzh_hadronic_VAJHU_up = -1.;
    double pzh_hadronic_VAJHU_dn = -1.;

    double ptth_VAJHU = -1.;
    double ptth_VAJHU_up = -1.;
    double ptth_VAJHU_dn = -1.;

    double pbbh_VAJHU = -1.;
    double pbbh_VAJHU_up = -1.;
    double pbbh_VAJHU_dn = -1.;

    bool hasAtLeastOneJet = (cleanedJetsPt30.size() > 0);
    bool hasAtLeastTwoJets = (cleanedJetsPt30.size() > 1);
    if (hasAtLeastOneJet) {
      std::vector<TLorentzVector> partPprod;
      partPprod.push_back(pL11);
      partPprod.push_back(pL12);
      partPprod.push_back(pL21);
      partPprod.push_back(pL22);

      std::vector<int> partIdprod;
      partIdprod.push_back(id11);
      partIdprod.push_back(id12);
      partIdprod.push_back(id21);
      partIdprod.push_back(id22);
      partIdprod.push_back(0); // If you know the jet flavors, even better, put the pdg ids for vhMELA
      partIdprod.push_back(0); // If you know the jet flavors, even better, put the pdg ids for vhMELA

      for (int jecnum = 0; jecnum < 3; jecnum++){
        TLorentzVector higgs_undec = pL11+pL12+pL21+pL22;

        TLorentzVector jet1, jet2;
        double djet_max = -1;
        double phjj_VAJHU_old_temp = -1;
        double pvbf_VAJHU_old_temp = -1;
        double phjj_VAJHU_new_temp = -1;
        double pvbf_VAJHU_new_temp = -1;
        double pAux_vbf_VAJHU_temp = 1;

        double phj_VAJHU_temp = -1;
        double pwh_hadronic_VAJHU_temp = -1;
        double pzh_hadronic_VAJHU_temp = -1;
        double ptth_VAJHU_temp = -1;
        double pbbh_VAJHU_temp = -1;

/*
        double jecnum_multiplier = 0;
        if (jecnum==1) jecnum_multiplier = 1.;
        else if (jecnum==2) jecnum_multiplier = -1.;
*/
        for (unsigned int firstjet = 0; firstjet < cleanedJetsPt30.size(); firstjet++){
          Float_t ratio1 = 1. /*+ jecnum_multiplier * cleanedJetsPt30[firstjet]->uncOnFourVectorScale()*/;
          jet1.SetXYZT(cleanedJetsPt30[firstjet]->p4().x(), cleanedJetsPt30[firstjet]->p4().y(), cleanedJetsPt30[firstjet]->p4().z(), cleanedJetsPt30[firstjet]->p4().t());
          jet1 *= ratio1;
          partPprod.push_back(jet1);

          for (unsigned int secondjet = 1; secondjet < cleanedJetsPt30.size(); secondjet++){
            if (secondjet <= firstjet) continue;
            Float_t ratio2 = 1. /*+ jecnum_multiplier * cleanedJetsPt30[secondjet]->uncOnFourVectorScale()*/;


            jet2.SetXYZT(cleanedJetsPt30[secondjet]->p4().x(), cleanedJetsPt30[secondjet]->p4().y(), cleanedJetsPt30[secondjet]->p4().z(), cleanedJetsPt30[secondjet]->p4().t());
            jet2 *= ratio2;
            partPprod.push_back(jet2);

            double phjj_temp = -1;
            double pvbf_temp = -1;
            combinedMEM.computeME(MEMNames::kJJ_SMHiggs_GG, MEMNames::kJHUGen, partPprod, partIdprod, phjj_temp);
            combinedMEM.computeME(MEMNames::kJJ_SMHiggs_VBF, MEMNames::kJHUGen, partPprod, partIdprod, pvbf_temp);
            double djet_temp = pvbf_temp / (pvbf_temp + phjj_temp);
            if (djet_temp > djet_max){
              phjj_VAJHU_new_temp = phjj_temp;
              pvbf_VAJHU_new_temp = pvbf_temp;
              djet_max = djet_temp;
            }

            if (firstjet == 0 && secondjet == 1){
              phjj_VAJHU_old_temp = phjj_temp;
              pvbf_VAJHU_old_temp = pvbf_temp;

              float pwh_temp = -1;
              float pzh_temp = -1;
              float ptth_temp = -1;
              float pbbh_temp = -1;
              myMela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::WH);
              myMela->computeProdP(jet1, 0, jet2, 0, higgs_undec, pwh_temp);
              myMela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZH);
              myMela->computeProdP(jet1, 0, jet2, 0, higgs_undec, pzh_temp);
              myMela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ttH);
              myMela->computeProdP(jet1, 0, jet2, 0, higgs_undec, ptth_temp);
              myMela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::bbH);
              myMela->computeProdP(jet1, 0, jet2, 0, higgs_undec, pbbh_temp);
              pwh_hadronic_VAJHU_temp = (double)pwh_temp;
              pzh_hadronic_VAJHU_temp = (double)pzh_temp;
              ptth_VAJHU_temp = (double)ptth_temp;
              pbbh_VAJHU_temp = (double)pbbh_temp;
            }

            partPprod.pop_back();
          }
          if (!hasAtLeastTwoJets){ // Compute H + 1 jet directly through Mela
            jet2.SetXYZT(0, 0, 0, 0);

            float phj_temp = -1;
            float pjvbf_temp = -1;
            float pAux_vbf_temp = -1;

            myMela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JH);
            myMela->computeProdP(jet1, 0, jet2, 0, higgs_undec, phj_temp);

            myMela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
            myMela->computeProdP(jet1, 0, jet2, 0, higgs_undec, pjvbf_temp); // Un-integrated ME
            myMela->get_PAux(pAux_vbf_temp); // = Integrated / un-integrated

            djet_max = (double)(pjvbf_temp*pAux_vbf_temp / (pjvbf_temp*pAux_vbf_temp + phj_temp));
            phj_VAJHU_temp = (double)phj_temp;
            pvbf_VAJHU_old_temp = (double)pjvbf_temp;
            pAux_vbf_VAJHU_temp = (double)pAux_vbf_temp;
            pvbf_VAJHU_new_temp = pvbf_VAJHU_old_temp;
          }

          partPprod.pop_back();
        }
        if (jecnum == 0){
          phjj_VAJHU_old = phjj_VAJHU_old_temp;
          pvbf_VAJHU_old = pvbf_VAJHU_old_temp;
          phjj_VAJHU_new = phjj_VAJHU_new_temp;
          pvbf_VAJHU_new = pvbf_VAJHU_new_temp;

          pAux_vbf_VAJHU = pAux_vbf_VAJHU_temp;

          phj_VAJHU = phj_VAJHU_temp;
          pwh_hadronic_VAJHU = pwh_hadronic_VAJHU_temp;
          pzh_hadronic_VAJHU = pzh_hadronic_VAJHU_temp;
          ptth_VAJHU = ptth_VAJHU_temp;
          pbbh_VAJHU = pbbh_VAJHU_temp;
        }
        if (jecnum == 1){
          phjj_VAJHU_old_up = phjj_VAJHU_old_temp;
          pvbf_VAJHU_old_up = pvbf_VAJHU_old_temp;
          phjj_VAJHU_new_up = phjj_VAJHU_new_temp;
          pvbf_VAJHU_new_up = pvbf_VAJHU_new_temp;

          pAux_vbf_VAJHU_up = pAux_vbf_VAJHU_temp;

          phj_VAJHU_up = phj_VAJHU_temp;
          pwh_hadronic_VAJHU_up = pwh_hadronic_VAJHU_temp;
          pzh_hadronic_VAJHU_up = pzh_hadronic_VAJHU_temp;
          ptth_VAJHU_up = ptth_VAJHU_temp;
          pbbh_VAJHU_up = pbbh_VAJHU_temp;
        }
        if (jecnum == 2){
          phjj_VAJHU_old_dn = phjj_VAJHU_old_temp;
          pvbf_VAJHU_old_dn = pvbf_VAJHU_old_temp;
          phjj_VAJHU_new_dn = phjj_VAJHU_new_temp;
          pvbf_VAJHU_new_dn = pvbf_VAJHU_new_temp;

          pAux_vbf_VAJHU_dn = pAux_vbf_VAJHU_temp;

          phj_VAJHU_dn = phj_VAJHU_temp;
          pwh_hadronic_VAJHU_dn = pwh_hadronic_VAJHU_temp;
          pzh_hadronic_VAJHU_dn = pzh_hadronic_VAJHU_temp;
          ptth_VAJHU_dn = ptth_VAJHU_temp;
          pbbh_VAJHU_dn = pbbh_VAJHU_temp;
        }
      }


    }

    //----------------------------------------------------------------------
      
    // Use extraLeps, extraZs, cleanedJetsPt30, and pfmet to compute additional discriminants

    float pzh_VAJHU=-1.;
    // ZH discriminant
    if (extraZs.size()>=1) {
      myMela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZH);
      double selfDHvvcoupl[SIZE_HVV_VBF][2] = { { 0 } };
      for (vector<const CompositeCandidate*>::const_iterator iZ = extraZs.begin(); iZ!=extraZs.end(); ++iZ){
	//TLorentzVector (math::XYZTLorentzVector)
	const reco::Candidate* Zal1 = (*iZ)->daughter(0);
	const reco::Candidate* Zal2 = (*iZ)->daughter(1);

	TLorentzVector aZLeps[2] = {tlv(Zal1->p4()), tlv(Zal2->p4())};
	int aZLepsId[2] = {Zal1->pdgId(),Zal2->pdgId()};
	float pzh_VAJHU_tmp;
	myMela->computeProdP(aZLeps, partP.data(), aZLepsId, partId.data(), false, selfDHvvcoupl, pzh_VAJHU_tmp); 
	pzh_VAJHU = std::max(pzh_VAJHU_tmp,pzh_VAJHU);
      }
    } 
    
    // ...
  
    
    //----------------------------------------------------------------------
    //--- Z1 kinematic refit (experimental)
    float m4lRef=0;
    float chi2LepZ1Ref=0;
    //float chi2ProbLepZ1Ref=0;
    float mZ1Ref=0;
    if (doKinFit && abs(id11)==13 && abs(id12)==13 && abs(id21)==13 && abs(id22)==13 ) {
      
      // Get the Z1 leptons 4 four-vectors refitted with the mass constraint (the fsr photon is already taken into account)
      DiLeptonKinFitter* fitter_dilep = new DiLeptonKinFitter( "fitter_dilep", "fitter_dilep", 91.1876);
      std::pair<TLorentzVector,TLorentzVector> lepZ1Ref = fitter_dilep->fit(Z1);

      chi2LepZ1Ref = fitter_dilep->getS();
      //chi2ProbLepZ1Ref = TMath::Prob(fitter_dilep->getS(), fitter_dilep->getNDF());
    
      math::XYZTLorentzVector pRef11(lepZ1Ref.first.Px(),lepZ1Ref.first.Py(),lepZ1Ref.first.Pz(),lepZ1Ref.first.E());
      math::XYZTLorentzVector pRef12(lepZ1Ref.second.Px(),lepZ1Ref.second.Py(),lepZ1Ref.second.Pz(),lepZ1Ref.second.E());
      
      m4lRef  = (pRef11+pRef12+p21+p22).mass();
      mZ1Ref  = (pRef11+pRef12).mass();

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

    if(doMEKD){
      myCand.addUserFloat("MEKD_LD", mekd_ld);
      myCand.addUserFloat("MEKD_PseudoLD", mekd_pseudold);
      myCand.addUserFloat("MEKD_GravLD", mekd_gravld);
    }
    myCand.addUserFloat("m4l",            (Z1Lm->p4()+Z1Lp->p4()+Z2Lm->p4()+Z2Lp->p4()).mass()); // mass without FSR
    if (doKinFit && abs(id11)==13 && abs(id12)==13 && abs(id21)==13 && abs(id22)==13 ) {
      myCand.addUserFloat("m4l_Z1Fit",  m4lRef ); // mass from Z1 refitted (FSR not considered in the refit procedure)
      myCand.addUserFloat("chi2_Z1Fit", chi2LepZ1Ref);
      myCand.addUserFloat("mZ1Fit",  mZ1Ref);
    }
    // add probabilities
    myCand.addUserFloat("p0plus_VAJHU",   p0plus_VAJHU);   // higgs, vector algebra, JHUgen
    myCand.addUserFloat("p0minus_VAJHU",  p0minus_VAJHU);  // pseudoscalar, vector algebra, JHUgen
    myCand.addUserFloat("p0plus_VAMCFM",  p0plus_VAMCFM);  // higgs, vector algebra, MCFM
    myCand.addUserFloat("p0hplus_VAJHU",  p0hplus_VAJHU);// 0h+ (high dimensional operator), vector algebra, JHUgen
    myCand.addUserFloat("p1_VAJHU",       p1_VAJHU);       // zprime, vector algebra, JHUgen,
    myCand.addUserFloat("p1plus_VAJHU",   p1plus_VAJHU);// 1+ (axial vector), vector algebra, JHUgen,
    myCand.addUserFloat("p1_prodIndep_VAJHU",       p1_prodIndep_VAJHU);       // zprime, vector algebra, JHUgen,
    myCand.addUserFloat("p1plus_prodIndep_VAJHU",   p1plus_prodIndep_VAJHU);// 1+ (axial vector), vector algebra, JHUgen,
    myCand.addUserFloat("p2_VAJHU",       p2_VAJHU);       // graviton, vector algebra, JHUgen,
    myCand.addUserFloat("p2qqb_VAJHU",    p2qqb_VAJHU);       // graviton, vector algebra, JHUgen,
    //backgrounds
    myCand.addUserFloat("bkg_VAMCFM",     bkg_VAMCFM);     // background, vector algebra, MCFM
    //myCand.addUserFloat("ggzz_VAMCFM",    ggzz_VAMCFM);    //  background, vector algebra, MCFM for ggzz   
    myCand.addUserFloat("ggzz_VAMCFM",  ggzz_VAMCFM); //background, vector algebra, MCFM for ggzz
    myCand.addUserFloat("ggzz_p0plus_VAMCFM",  ggzz_p0plus_VAMCFM); //background, vector algebra, MCFM for ggzz
    myCand.addUserFloat("ggzz_c1_VAMCFM",  ggzz_c1_VAMCFM); //higgs + background + interference, vector algebra, MCFM for ggzz
    myCand.addUserFloat("ggzz_c5_VAMCFM",  ggzz_c5_VAMCFM); //higgs (w/ 5 times coupling) + background + interference, vector algebra, MCFM for ggzz
    myCand.addUserFloat("ggzz_ci_VAMCFM",  ggzz_ci_VAMCFM); //higgs (w/ complex coupling) + background + interference, vector algebra, MCFM for ggzz
    myCand.addUserFloat("p2_prodIndep_VAJHU",    p2_prodIndep_VAJHU);       // graviton, vector algebra, JHUgen,
    myCand.addUserFloat("p2hplus_VAJHU",    p2hplus_VAJHU);       // graviton, vector algebra, JHUgen,
    myCand.addUserFloat("p2hminus_VAJHU",    p2hminus_VAJHU);       // graviton, vector algebra, JHUgen,
    myCand.addUserFloat("p2bplus_VAJHU",    p2bplus_VAJHU);       // graviton, vector algebra, JHUgen,
    myCand.addUserFloat("bkg_prodIndep_VAMCFM",     bkg_prodIndep_VAMCFM);     // background, vector algebra, MCFM

    myCand.addUserFloat("p2hplus_qqb_VAJHU"         ,   p2hplus_qqb_VAJHU);
    myCand.addUserFloat("p2hplus_prodIndep_VAJHU"   ,   p2hplus_prodIndep_VAJHU);
    myCand.addUserFloat("p2hminus_qqb_VAJHU"        ,   p2hminus_qqb_VAJHU);
    myCand.addUserFloat("p2hminus_prodIndep_VAJHU"  ,   p2hminus_prodIndep_VAJHU);
    myCand.addUserFloat("p2bplus_qqb_VAJHU"         ,   p2bplus_qqb_VAJHU);
    myCand.addUserFloat("p2bplus_prodIndep_VAJHU"   ,   p2bplus_prodIndep_VAJHU);

    myCand.addUserFloat("p2h2plus_gg_VAJHU"         ,   p2h2plus_gg_VAJHU);
    myCand.addUserFloat("p2h2plus_qqbar_VAJHU"      ,   p2h2plus_qqbar_VAJHU);
    myCand.addUserFloat("p2h2plus_prodIndep_VAJHU"  ,   p2h2plus_prodIndep_VAJHU);
    myCand.addUserFloat("p2h3plus_gg_VAJHU"         ,   p2h3plus_gg_VAJHU);
    myCand.addUserFloat("p2h3plus_qqbar_VAJHU"      ,   p2h3plus_qqbar_VAJHU);
    myCand.addUserFloat("p2h3plus_prodIndep_VAJHU"  ,   p2h3plus_prodIndep_VAJHU);
    myCand.addUserFloat("p2h6plus_gg_VAJHU"         ,   p2h6plus_gg_VAJHU);
    myCand.addUserFloat("p2h6plus_qqbar_VAJHU"      ,   p2h6plus_qqbar_VAJHU);
    myCand.addUserFloat("p2h6plus_prodIndep_VAJHU"  ,   p2h6plus_prodIndep_VAJHU);
    myCand.addUserFloat("p2h7plus_gg_VAJHU"         ,   p2h7plus_gg_VAJHU);
    myCand.addUserFloat("p2h7plus_qqbar_VAJHU"      ,   p2h7plus_qqbar_VAJHU);
    myCand.addUserFloat("p2h7plus_prodIndep_VAJHU"  ,   p2h7plus_prodIndep_VAJHU);
    myCand.addUserFloat("p2h9minus_gg_VAJHU"        ,   p2h9minus_gg_VAJHU);
    myCand.addUserFloat("p2h9minus_qqbar_VAJHU"     ,   p2h9minus_qqbar_VAJHU);
    myCand.addUserFloat("p2h9minus_prodIndep_VAJHU" ,   p2h9minus_prodIndep_VAJHU);
    myCand.addUserFloat("p2h10minus_gg_VAJHU"       ,   p2h10minus_gg_VAJHU);       
    myCand.addUserFloat("p2h10minus_qqbar_VAJHU"    ,   p2h10minus_qqbar_VAJHU);  
    myCand.addUserFloat("p2h10minus_prodIndep_VAJHU",   p2h10minus_prodIndep_VAJHU);

    //pt/rapidity
    //myCand.addUserFloat("p0_pt",          p0_pt);          // multiplicative probability for signal pt
    //myCand.addUserFloat("p0_y",           p0_y);           // multiplicative probability for signal y
    //myCand.addUserFloat("bkg_pt",         bkg_pt);         // multiplicative probability for bkg pt
    //myCand.addUserFloat("bkg_y",          bkg_y);          // multiplicative probability for bkg y
    
    // supermela
    myCand.addUserFloat("p0plus_m4l",     p0plus_m4l);  // signal m4l probability as in datacards
    myCand.addUserFloat("bkg_m4l",        bkg_m4l);     // backgroun m4l probability as in datacards
    myCand.addUserFloat("p0plus_m4l_ScaleUp",p0plus_m4l_ScaleUp);// signal m4l probability for systematics
    myCand.addUserFloat("bkg_m4l_ScaleUp",bkg_m4l_ScaleUp);// backgroun m4l probability for systematics
    myCand.addUserFloat("p0plus_m4l_ScaleDown",p0plus_m4l_ScaleDown);// signal m4l probability for systematics
    myCand.addUserFloat("bkg_m4l_ScaleDown",bkg_m4l_ScaleDown);// backgroun m4l probability for systematics
    myCand.addUserFloat("p0plus_m4l_ResUp",p0plus_m4l_ResUp);// signal m4l probability for systematics
    myCand.addUserFloat("bkg_m4l_ResUp",bkg_m4l_ResUp);// backgroun m4l probability for systematics
    myCand.addUserFloat("p0plus_m4l_ResDown",p0plus_m4l_ResDown);// signal m4l probability for systematics
    myCand.addUserFloat("bkg_m4l_ResDown",bkg_m4l_ResDown);// backgroun m4l probability for systematics
  
    // spinMELA
    myCand.addUserFloat("pg1g4_mela",      pg1g4_mela);
    myCand.addUserFloat("pg1g4_VAJHU",     pg1g4_VAJHU);
    myCand.addUserFloat("pg1g4_pi2_VAJHU", pg1g4_pi2_VAJHU);
    myCand.addUserFloat("pg1g2_pi2_VAJHU", pg1g2_pi2_VAJHU);
    myCand.addUserFloat("pg1g2_mela",      pg1g2_mela);
    myCand.addUserFloat("pg1g2_VAJHU",     pg1g2_VAJHU);

    myCand.addUserFloat("p0_g1prime2_VAJHU",p0_g1prime2_VAJHU);
    myCand.addUserFloat("pg1g1prime2_VAJHU",pg1g1prime2_VAJHU);
    myCand.addUserFloat("Dgg10_VAMCFM",Dgg10_VAMCFM);

    myCand.addUserFloat("pzzzg_VAJHU",pzzzg_VAJHU);
    myCand.addUserFloat("pzzgg_VAJHU",pzzgg_VAJHU);
    myCand.addUserFloat("p0Zgs_VAJHU",p0Zgs_VAJHU);
    myCand.addUserFloat("p0gsgs_VAJHU",p0gsgs_VAJHU);

    myCand.addUserFloat("pzzzg_PS_VAJHU",pzzzg_PS_VAJHU);
    myCand.addUserFloat("pzzgg_PS_VAJHU",pzzgg_PS_VAJHU);
    myCand.addUserFloat("p0Zgs_PS_VAJHU",p0Zgs_PS_VAJHU);
    myCand.addUserFloat("p0gsgs_PS_VAJHU",p0gsgs_PS_VAJHU);

    myCand.addUserFloat("p0Zgs_g1prime2_VAJHU", p0Zgs_g1prime2_VAJHU);
    myCand.addUserFloat("pzzzg_g1prime2_VAJHU", pzzzg_g1prime2_VAJHU);
    myCand.addUserFloat("pzzzg_g1prime2_pi2_VAJHU", pzzzg_g1prime2_pi2_VAJHU);

    // Production MELA
    myCand.addUserFloat("phjj_VAJHU_old", phjj_VAJHU_old);
    myCand.addUserFloat("pvbf_VAJHU_old", pvbf_VAJHU_old);
    myCand.addUserFloat("phjj_VAJHU_old_up", phjj_VAJHU_old_up);
    myCand.addUserFloat("pvbf_VAJHU_old_up", pvbf_VAJHU_old_up);
    myCand.addUserFloat("phjj_VAJHU_old_dn", phjj_VAJHU_old_dn);
    myCand.addUserFloat("pvbf_VAJHU_old_dn", pvbf_VAJHU_old_dn);
    myCand.addUserFloat("phjj_VAJHU_new", phjj_VAJHU_new);
    myCand.addUserFloat("pvbf_VAJHU_new", pvbf_VAJHU_new);
    myCand.addUserFloat("phjj_VAJHU_new_up", phjj_VAJHU_new_up);
    myCand.addUserFloat("pvbf_VAJHU_new_up", pvbf_VAJHU_new_up);
    myCand.addUserFloat("phjj_VAJHU_new_dn", phjj_VAJHU_new_dn);
    myCand.addUserFloat("pvbf_VAJHU_new_dn", pvbf_VAJHU_new_dn);

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


    // VH
    myCand.addUserFloat("pzh_VAJHU",pzh_VAJHU);

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
        //      cout << "Pass best cand presel for " << iCRname << " " << myCand.mass() << " " << Z1->mass() << " " << Z2->mass() << " " << Z2->daughter(0)->pt() << " " << Z2->daughter(1)->pt() << " " <<  endl;
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
      
//       // For debug purposes
//       cout << "preSelCands[iCRname].size() = " << preSelCands[iCRname].size() << endl;
//       // for (int i=0; i<(int)preSelCands[iCRname].size(); ++i){
//       for (vector<int>::const_iterator iCand = preSelCands[iCRname].begin(); iCand<preSelCands[iCRname].end(); ++iCand) {
//      const reco::Candidate* dau0 = (*result)[*iCand].daughter(0)->masterClone().get();
//      const reco::Candidate* dau1 = (*result)[*iCand].daughter(1)->masterClone().get();
//              cout << "  [ZZCandidateFiller] candidate " << *iCand << ": mass daughter 0 = " << dau0->mass() << ", mass daughter 1 = " << dau1->mass() 
//           << ", ptSum dau0 = "<< dau0->daughter(0)->masterClone().get()->pt() + dau0->daughter(1)->masterClone().get()->pt()
//           << ", ptSum dau1 = "<< dau1->daughter(0)->masterClone().get()->pt() + dau1->daughter(1)->masterClone().get()->pt() << endl;          
//       }      
//       cout << "[ZZCandidateFiller] was chosen candidate with index = " << *std::min_element( preSelCands[iCRname].begin(), preSelCands[iCRname].end(), myComp) << endl;

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
ZZCandidateFiller::getPairMass(const reco::Candidate* lp, const reco::Candidate* lm, map<const reco::Candidate*, math::XYZTLorentzVector>& photons, float& mass, int& ID) {

//   cout << "GPM " << (long) lp << " " << (long) lm << photons.size();
//   for  (map<const reco::Candidate*, math::XYZTLorentzVector>::iterator i = photons.begin(); i!=photons.end(); ++i) {
//     cout << " (" << (long) i->first << " " << i->second << ") ";
//   }
  

  math::XYZTLorentzVector llp4 = lp->p4()+lm->p4();
//   cout << llp4.mass();
  auto lpp = photons.find(lp);
  auto lmp = photons.find(lm);
  if (lpp!=photons.end()) llp4+=lpp->second;
  if (lmp!=photons.end()) llp4+=lmp->second;
//   cout << " " << mass <<  endl;
  mass=llp4.mass();
  ID=lp->pdgId()*lm->pdgId();
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ZZCandidateFiller);

