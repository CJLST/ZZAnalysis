/** \class ZZCandidateFiller
 *
 *  No description available.
 *
 *  $Date: 2013/12/12 16:19:40 $
 *  $Revision: 1.72 $
 *  \author N. Amapane - Torino
 *  \author C. Botta - Torino
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>

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
#include <AnalysisDataFormats/CMGTools/interface/PFJet.h>
#include "ZZAnalysis/AnalysisStep/interface/VBFCandidateJetSelector.h"
#include <ZZAnalysis/AnalysisStep/interface/Fisher.h>



#include "TH2F.h"
#include "TFile.h"

#include <string>

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

  edm::InputTag theCandidateTag;
  const CutSet<pat::CompositeCandidate> preBestCandSelection;
  const CutSet<pat::CompositeCandidate> cuts;
  int sampleType;
  float superMelaMass;
  Mela mela;
  MEMs combinedMEM;
  bool embedDaughterFloats;
  bool ZRolesByMass;
  reco::CompositeCandidate::role_collection rolesZ1Z2;  
  reco::CompositeCandidate::role_collection rolesZ2Z1;
  bool isMC;
  TH2F* corrSigmaMu;
  TH2F* corrSigmaEle;
};




ZZCandidateFiller::ZZCandidateFiller(const edm::ParameterSet& iConfig) :
  theCandidateTag(iConfig.getParameter<edm::InputTag>("src")),
  preBestCandSelection(iConfig.getParameter<edm::ParameterSet>("bestCandAmong")),
  cuts(iConfig.getParameter<edm::ParameterSet>("flags")),
  sampleType(iConfig.getParameter<int>("sampleType")),
  superMelaMass(iConfig.getParameter<double>("superMelaMass")),
  //FIXME: should use LEPTON_SETUP instead of sampleType for mela and combinedMEM. This will be an issue for samples rescaled to different sqrts (none at present)
  mela((sampleType==2011)?7:8,superMelaMass),
  combinedMEM((sampleType==2011)?7.:8,superMelaMass,"CTEQ6L"),
  embedDaughterFloats(iConfig.getUntrackedParameter<bool>("embedDaughterFloats",true)),
  ZRolesByMass(iConfig.getParameter<bool>("ZRolesByMass")),
  isMC(iConfig.getParameter<bool>("isMC"))
{
  produces<pat::CompositeCandidateCollection>();
  
  rolesZ1Z2.push_back("Z1");
  rolesZ1Z2.push_back("Z2");
  rolesZ2Z1.push_back("Z2");
  rolesZ2Z1.push_back("Z1");
  
  string ebePath;
//   if (sampleType==2011) { 
//     edm::FileInPath fip("ZZAnalysis/AnalysisStep/data/ebeOverallCorrections.HCP2012.v1.root");
//     ebePath=fip.fullPath();
//   } else {
//     edm::FileInPath fip("ZZAnalysis/AnalysisStep/data/ebeOverallCorrections.Legacy2013.v0.root");
//     ebePath=fip.fullPath();
//   }
  
  edm::FileInPath fip("ZZAnalysis/AnalysisStep/data/ebeOverallCorrections.Legacy2013.v0.root");
  ebePath=fip.fullPath();

  
  TFile* fCorrSigma = new TFile(ebePath.data()); // FIXME: is leaked
  
  std::string sigmaCorrType = (isMC?"mc":"reco");
  std::string sigmaCorrYear = ((sampleType==2011)?"42x":"53x"); //FIXME: should use LEPTON_SETUP instead. This will be an issue for samples rescaled to different sqrts (none at present)

  corrSigmaMu=  (TH2F*)fCorrSigma->Get(("mu_"+sigmaCorrType+sigmaCorrYear).data()); 
  corrSigmaEle= (TH2F*)fCorrSigma->Get(("el_"+sigmaCorrType+sigmaCorrYear).data());
}


void ZZCandidateFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){  
  using namespace edm;
  using namespace std;
  using namespace reco;

  std::auto_ptr<pat::CompositeCandidateCollection> result(new pat::CompositeCandidateCollection);

  const float ZmassValue = 91.1876;

  // Get LLLL candidates
  Handle<View<CompositeCandidate> > LLLLCands;
  iEvent.getByLabel(theCandidateTag, LLLLCands);

  // Get jets
  Handle<edm::View<cmg::PFJet> > pfjetscoll;
  iEvent.getByLabel("cmgPFJetSel", pfjetscoll);

  // Get processID
//   edm::Handle<GenEventInfoProduct> gen;
//   iEvent.getByLabel( "generator", gen );
//   int processID = gen->signalProcessID();

  // to calculate mass resolution
  CompositeCandMassResolution errorBuilder;       
  errorBuilder.init(iSetup);

  float minDZMass=999;
  int bestZFlavour=0;
  // Look for the best Z in the event. In 2e2mu, there one isBestZ in the Zee and one in the Zmm collections!
  for( View<CompositeCandidate>::const_iterator cand = LLLLCands->begin(); cand != LLLLCands->end(); ++ cand ) {
    for (int i=0;i<2;++i) {
      if (userdatahelpers::getUserFloat(cand->daughter(i),"isBestZ")) {
	float dZMass=fabs(cand->daughter(i)->mass()-ZmassValue);
	if (dZMass<minDZMass){
	  minDZMass=dZMass;
	  bestZFlavour=cand->daughter(i)->daughter(0)->pdgId();
	}
      }
    }
  }

  vector<int> bestCandIdx(preBestCandSelection.size(),-1); 
  vector<float> maxPtSum(preBestCandSelection.size(),-1); 
  for( View<CompositeCandidate>::const_iterator cand = LLLLCands->begin(); cand != LLLLCands->end(); ++ cand ) {
    int i = distance(LLLLCands->begin(),cand);

    pat::CompositeCandidate myCand(*cand);

    if (embedDaughterFloats){  
      userdatahelpers::embedDaughterData(myCand);
    }
    

    const reco::CompositeCandidate::role_collection* ZRoles = &rolesZ1Z2;
    int iZ1 = 0; // index of the Z closest to mZ. This will be different to rolesZ1Z2
    int iZ2 = 1;
    if(fabs(myCand.daughter(0)->mass()-ZmassValue)>=fabs(myCand.daughter(1)->mass()-ZmassValue)){
      swap(iZ1,iZ2);
      if (ZRolesByMass) { 
	ZRoles = &rolesZ2Z1;
      }
    }

    //--- Set "Z1"/"Z2" labels. This allows to call e.g. aHiggs->daughter("Z1"). 
    // if ZRolesByMass is true, 'Z1' becoes Z closest to mZ; otherwise Z1 = daughter(0). The latter is used for control regions.
    myCand.setRoles(*ZRoles);
    myCand.applyRoles();

    //--- Embed shortcut variables
    vector<double> combRelIso;
    combRelIso.push_back(myCand.userFloat("d0.d0.combRelIso"));
    combRelIso.push_back(myCand.userFloat("d0.d1.combRelIso"));
    combRelIso.push_back(myCand.userFloat("d1.d0.combRelIso"));
    combRelIso.push_back(myCand.userFloat("d1.d1.combRelIso"));
    sort(combRelIso.begin(),combRelIso.end());
    float iso34 = combRelIso[2]+combRelIso[3];    // Pairwise isolation (for cut on all pairs)
  
    vector<double> SIPS;
    SIPS.push_back(myCand.userFloat("d0.d0.SIP"));
    SIPS.push_back(myCand.userFloat("d0.d1.SIP"));
    SIPS.push_back(myCand.userFloat("d1.d0.SIP"));
    SIPS.push_back(myCand.userFloat("d1.d1.SIP"));
    sort(SIPS.begin(),SIPS.end());
    float SIP4 = SIPS[3]; // worst SIP value


    const reco::Candidate* Z1= myCand.daughter(iZ1);
    const reco::Candidate* Z2= myCand.daughter(iZ2);
    const reco::Candidate* Z1Lp= Z1->daughter(0);	
    const reco::Candidate* Z1Lm= Z1->daughter(1);
    const reco::Candidate* Z2Lp= Z2->daughter(0);
    const reco::Candidate* Z2Lm= Z2->daughter(1);

    // Get the 4 four-vectors and add FSR (to be used to compute angles)
    math::XYZTLorentzVector p11 = Z1Lp->p4();
    math::XYZTLorentzVector p12 = Z1Lm->p4();
    math::XYZTLorentzVector p21 = Z2Lp->p4();
    math::XYZTLorentzVector p22 = Z2Lm->p4();
    int id11 = Z1Lp->pdgId();
    int id12 = Z1Lm->pdgId();
    int id21 = Z2Lp->pdgId();
    int id22 = Z2Lm->pdgId();

    int d0FSR = (static_cast<const pat::CompositeCandidate*>(Z1->masterClone().get()))->userFloat("dauWithFSR");
    if (d0FSR>=0) {
      if (Z1->numberOfDaughters()!=3) cout << "ERROR: ZZCandidateFiller: problem in FSR" << endl;
      if (d0FSR==0) {
 	p11 = p11 + Z1->daughter(2)->p4();
      } else if (d0FSR==1){
 	p12 = p12 + Z1->daughter(2)->p4();      
      }
    }

    int d1FSR = (static_cast<const pat::CompositeCandidate*>(Z2->masterClone().get()))->userFloat("dauWithFSR");
    if (d1FSR>=0) {
      if (Z2->numberOfDaughters()!=3) cout << "ERROR: ZZCandidateFiller: problem in FSR" << endl;
      if (d1FSR==0) {
 	p21 = p21 + Z2->daughter(2)->p4();
      } else if (d1FSR==1){
 	p22 = p22 + Z2->daughter(2)->p4();      
      }
    }

    TLorentzVector pL11(p11.x(),p11.y(),p11.z(),p11.t());
    TLorentzVector pL12(p12.x(),p12.y(),p12.z(),p12.t());
    TLorentzVector pL21(p21.x(),p21.y(),p21.z(),p21.t());
    TLorentzVector pL22(p22.x(),p22.y(),p22.z(),p22.t());

    // Now sort pointers according to the p/m name; do nothing for the same-sign collections used for CRs
    if (Z1Lp->charge() < 0 && Z1Lp->charge()*Z1Lm->charge()<0) swap(Z1Lp,Z1Lm);    
    if (Z2Lp->charge() < 0 && Z2Lp->charge()*Z2Lm->charge()<0) swap(Z2Lp,Z2Lm);



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
    

    //--- Embed angular information and KD
    float costheta1=0, costheta2=0, phi=0, costhetastar=0, phistar1=0;

    // higgs azimuthal angles
    float xi=0,xistar=0;

    //float ld=0, psig=0, pbkg=0;
    double mekd_ld=0,mekd_pseudold=0,mekd_gravld=0,mekd_A=0,mekd_B=0,mekd_AP=0,mekd_BP=0,mekd_AG=0,mekd_BG=0;

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

    ///--- compute angles
    mela::computeAngles(pL11,id11,pL12,id12,pL21,id21,pL22,id22,costhetastar,costheta1,costheta2,phi,phistar1);
  
    ///--- compute higgs azimuthal angles, xi
    TLorentzVector higgs = pL11+pL12+pL21+pL22;
    TLorentzVector Z14vec = pL11+pL12;
    TVector3 Xaxis(1,0,0);
    xi = higgs.Phi();

    // boost Z1 into rest frame of higgs
    // xistar is the angle between Z decay plane and x-axis 
    Z14vec.Boost(-higgs.BoostVector());
    xistar = Z14vec.Phi();
    // higgs.Boost(-higgs.BoostVector());

    if(doMEKD){ 
      int prodid=id11*id21*id12*id22;
      if(prodid==28561 || prodid==14641 || prodid==20449){
 	combinedMEM.computeKD(MEMNames::kSMHiggs,MEMNames::kqqZZ,MEMNames::kMEKD,partP,partId, mekd_ld,mekd_A,mekd_B);
 	combinedMEM.computeKD(MEMNames::k0minus,MEMNames::kSMHiggs,MEMNames::kMEKD,partP,partId,mekd_pseudold,mekd_AP,mekd_BP);
 	combinedMEM.computeKD(MEMNames::k2mplus_gg,MEMNames::kSMHiggs,MEMNames::kMEKD,partP,partId,mekd_gravld,mekd_AG,mekd_BG);
      }else{
	mekd_ld=-1;
	mekd_pseudold=-1;
	mekd_gravld=-1;
      }
    }
  
    //--- other kinematic discriminants
//     float m1 = Z1->mass();
//     float m2 = Z2->mass();
//     float mzz = myCand.mass();

    // probabilities
    //double p0plus_melaNorm,p0plus_mela,p0minus_mela,p0hplus_mela;
    double p0plus_VAJHU,p0minus_VAJHU,p0plus_VAMCFM,p0hplus_VAJHU;
    //double p1_mela,p1plus_mela;
    double p1_VAJHU,p1plus_VAJHU,p2_VAJHU,p2qqb_VAJHU; 
    double bkg_VAMCFM,bkg_prodIndep_VAMCFM; //bkg_VAMCFMNorm;
    double ggzz_VAMCFM,ggzz_c1_VAMCFM,ggzz_c5_VAMCFM,ggzz_ci_VAMCFM;
		double ggzz_p0plus_VAMCFM;
    //double p2_mela, p2qqb_mela;
    //double bkg_mela;
    //float pt4l, Y4l;
    //float p0_pt,p0_y,bkg_pt,bkg_y;   // rapidity/pt

    // exotic spin-2 models
    double p2hplus_VAJHU, p2hminus_VAJHU,p2bplus_VAJHU;
    //double p2hplus_mela, p2hminus_mela,p2bplus_mela;

    // production independent probabilites
    //double p1_prodIndep_mela,p1plus_prodIndep_mela,p2_prodIndep_mela;
    double p1_prodIndep_VAJHU,p1plus_prodIndep_VAJHU,p2_prodIndep_VAJHU;

    double p0plus_m4l,bkg_m4l; //supermela
    double p0plus_m4l_ScaleUp,bkg_m4l_ScaleUp,p0plus_m4l_ScaleDown,bkg_m4l_ScaleDown,p0plus_m4l_ResUp,bkg_m4l_ResUp,p0plus_m4l_ResDown,bkg_m4l_ResDown; // supermela uncertainties

    // No longer used - AJW
    //int flavor; flavor=0;
    //if(abs(id11)==abs(id21) && abs(id11)==11 )// 4e
    //  flavor = 1;
    //if(abs(id11)==abs(id21) && abs(id11)==13 )// 4mu
    //  flavor = 2;
    //if(abs(id11)!=abs(id21))                  // 2e2mu
    //  flavor = 3;
    
//     pt4l = myCand.p4().pt();    
//     Y4l  = myCand.p4().Rapidity();    
    
    /*
    combinedMEM.computeME(MEMNames::kSMHiggs     ,MEMNames::kMELA_HCP  ,partP,partId,p0plus_melaNorm);// higgs, analytic distribution, normalized
    combinedMEM.computeME(MEMNames::kSMHiggs     ,MEMNames::kAnalytical,partP,partId,p0plus_mela);    // higgs, analytic distribution
    combinedMEM.computeME(MEMNames::k0minus      ,MEMNames::kAnalytical,partP,partId,p0minus_mela);   // pseudoscalar, analytic distribution	   
    combinedMEM.computeME(MEMNames::k0hplus      ,MEMNames::kAnalytical,partP,partId,p0hplus_mela);   // 0h+, analytic distribution		   
    combinedMEM.computeME(MEMNames::k1minus_prodIndep      ,MEMNames::kAnalytical,partP,partId,p1_prodIndep_mela);        // zprime, analytic distribution			    
    combinedMEM.computeME(MEMNames::k1minus      ,MEMNames::kAnalytical,partP,partId,p1_mela);        // zprime, analytic distribution			    
    combinedMEM.computeME(MEMNames::k1plus       ,MEMNames::kAnalytical,partP,partId,p1plus_mela);    // 1+, analytic distribution				    
    combinedMEM.computeME(MEMNames::k1plus_prodIndep       ,MEMNames::kAnalytical,partP,partId,p1plus_prodIndep_mela);    // 1+, analytic distribution				    
    combinedMEM.computeME(MEMNames::k2mplus_gg   ,MEMNames::kAnalytical,partP,partId,p2_mela );	      // graviton, analytic distribution			    
    combinedMEM.computeME(MEMNames::k2mplus_prodIndep   ,MEMNames::kAnalytical,partP,partId,p2_prodIndep_mela );	      // graviton, analytic distribution			    
    combinedMEM.computeME(MEMNames::k2mplus_qqbar,MEMNames::kAnalytical,partP,partId,p2qqb_mela);     // graviton produced by qqbar vector algebra, analytical,   
    combinedMEM.computeME(MEMNames::k2hplus,MEMNames::kAnalytical,partP,partId,p2hplus_mela);     // graviton produced by qqbar vector algebra, analytical,   
    combinedMEM.computeME(MEMNames::k2hminus,MEMNames::kAnalytical,partP,partId,p2hminus_mela);     // graviton produced by qqbar vector algebra, analytical,   
    combinedMEM.computeME(MEMNames::k2bplus,MEMNames::kAnalytical,partP,partId,p2bplus_mela);     // graviton produced by qqbar vector algebra, analytical,   
    combinedMEM.computeME(MEMNames::kqqZZ        ,MEMNames::kAnalytical,partP,partId,bkg_mela);       // background,  analytic distribution
    */

    combinedMEM.computeME(MEMNames::kSMHiggs     ,MEMNames::kJHUGen    ,partP,partId,p0plus_VAJHU);   // higgs, vector algebra, JHUgen	   
    combinedMEM.computeME(MEMNames::k0minus      ,MEMNames::kJHUGen    ,partP,partId,p0minus_VAJHU);  // pseudoscalar, vector algebra, JHUgen   
    combinedMEM.computeME(MEMNames::kSMHiggs     ,MEMNames::kMCFM      ,partP,partId,p0plus_VAMCFM);  // higgs, vector algebra, MCFM            
    combinedMEM.computeME(MEMNames::k0hplus      ,MEMNames::kJHUGen    ,partP,partId,p0hplus_VAJHU);  // 0h+ (high dimensional operator), vector algebra, JHUgen
    combinedMEM.computeME(MEMNames::k1minus      ,MEMNames::kJHUGen    ,partP,partId,p1_VAJHU);	      // zprime, vector algebra, JHUgen,			    
    combinedMEM.computeME(MEMNames::k1minus_prodIndep      ,MEMNames::kJHUGen    ,partP,partId,p1_prodIndep_VAJHU);	      // zprime, vector algebra, JHUgen,			    
    combinedMEM.computeME(MEMNames::k1plus       ,MEMNames::kJHUGen    ,partP,partId,p1plus_VAJHU);   // 1+ (axial vector), vector algebra, JHUgen,		    
    combinedMEM.computeME(MEMNames::k1plus_prodIndep       ,MEMNames::kJHUGen    ,partP,partId,p1plus_prodIndep_VAJHU);   // 1+ (axial vector), vector algebra, JHUgen,		    

    combinedMEM.computeME(MEMNames::k2mplus_gg   ,MEMNames::kJHUGen    ,partP,partId,p2_VAJHU);	      // graviton, vector algebra, JHUgen,			    
    combinedMEM.computeME(MEMNames::k2mplus_prodIndep   ,MEMNames::kJHUGen    ,partP,partId,p2_prodIndep_VAJHU);	      // graviton, vector algebra, JHUgen,			    
    combinedMEM.computeME(MEMNames::k2mplus_qqbar,MEMNames::kJHUGen    ,partP,partId,p2qqb_VAJHU);    // graviton produced by qqbar vector algebra, JHUgen,       
    combinedMEM.computeME(MEMNames::k2hplus,MEMNames::kJHUGen    ,partP,partId,p2hplus_VAJHU);    // graviton produced by qqbar vector algebra, JHUgen,       
    combinedMEM.computeME(MEMNames::k2hminus,MEMNames::kJHUGen    ,partP,partId,p2hminus_VAJHU);    // graviton produced by qqbar vector algebra, JHUgen,       
    combinedMEM.computeME(MEMNames::k2bplus,MEMNames::kJHUGen    ,partP,partId,p2bplus_VAJHU);    // graviton produced by qqbar vector algebra, JHUgen,       
    combinedMEM.computeME(MEMNames::kqqZZ_prodIndep ,MEMNames::kMCFM      ,partP,partId,bkg_prodIndep_VAMCFM);     // background, vector algebra, MCFM
    combinedMEM.computeME(MEMNames::kqqZZ        ,MEMNames::kMCFM      ,partP,partId,bkg_VAMCFM);     // background, vector algebra, MCFM
    //combinedMEM.computeME(MEMNames::kggZZ        ,MEMNames::kMCFM      ,partP,partId,ggzz_VAMCFM);    // background, vector algebra, MCFM for ggzz

    //MEs for ggZZ + signal interference
//    double coupling[2]={0.0,0.0};
//    combinedMEM.computeME(MEMNames::kSMHiggs, MEMNames::kMCFM, partP, partId, coupling, ggzz_VAMCFM);
//    coupling[0]=1.0;
//    combinedMEM.computeME(MEMNames::kSMHiggs, MEMNames::kMCFM, partP, partId, coupling, ggzz_c1_VAMCFM);
//    coupling[0]=5.0;
//    combinedMEM.computeME(MEMNames::kSMHiggs, MEMNames::kMCFM, partP, partId, coupling, ggzz_c5_VAMCFM);
//    coupling[0]=0.0;coupling[1]=1.0;
//    combinedMEM.computeME(MEMNames::kSMHiggs, MEMNames::kMCFM, partP, partId, coupling, ggzz_ci_VAMCFM);    
		vector<complex<double>> *coupling = new vector<complex<double> >;
		vector<complex<double>> *couplingprod= new vector<complex<double> >;
    combinedMEM.computeME(MEMNames::kggZZ, MEMNames::kMCFM, partP, partId, ggzz_VAMCFM);
    combinedMEM.computeME(MEMNames::kggZZ_SMHiggs, MEMNames::kMCFM, partP, partId, ggzz_p0plus_VAMCFM);
		
		complex<double> coup(1.,0.); 
		coupling->push_back(coup);
    combinedMEM.computeME(MEMNames::kggZZ_SMHiggs, MEMNames::kJHUGen, partP, partId, couplingprod,coupling, ggzz_c1_VAMCFM);
		coupling->clear();
    coup=5.;
    coupling->push_back(coup);
    combinedMEM.computeME(MEMNames::kggZZ_SMHiggs, MEMNames::kJHUGen, partP, partId, couplingprod,coupling, ggzz_c5_VAMCFM);
		coupling->clear();
    coup.real(0.);
    coup.imag(1.);
    coupling->push_back(coup);
    combinedMEM.computeME(MEMNames::kggZZ_SMHiggs, MEMNames::kJHUGen, partP, partId, couplingprod,coupling, ggzz_ci_VAMCFM);  
    // save sapce by not storing KDs that we can recreate 

    //double MorindKD_mela,MorindKD_pseudo,MorindKD_highdim,MorindKD_vec,MorindKD_psvec,MorindKD_gggrav,MorindKD_qqgrav;
    //double ME_A,ME_B;
    // Moriond KDs 
    //combinedMEM.computeKD(MEMNames::kSMHiggs, MEMNames::kJHUGen, MEMNames::kqqZZ, MEMNames::kMCFM, &MEMs::probRatio, MorindKD_mela ,ME_A,ME_B);// => recreate as p0plus_VAJHU/(p0plus_VAJHU + bkg_VAMCFMNorm) 
    //combinedMEM.computeKD(MEMNames::kSMHiggs, MEMNames::kJHUGen, MEMNames::k0minus, MEMNames::kJHUGen, &MEMs::probRatio, MorindKD_pseudo ,ME_A,ME_B);// => recreate as  p0plus_VAJHU/(p0plus_VAJHU + p0minus_VAJHU) 
    //combinedMEM.computeKD(MEMNames::kSMHiggs, MEMNames::kJHUGen, MEMNames::k0hplus, MEMNames::kJHUGen, &MEMs::probRatio, MorindKD_highdim ,ME_A,ME_B);// => recreate as p0plus_VAJHU/(p0plus_VAJHU + p0hplus_VAJHU) 
    //combinedMEM.computeKD(MEMNames::kSMHiggs, MEMNames::kJHUGen, MEMNames::k1minus, MEMNames::kJHUGen, &MEMs::probRatio, MorindKD_psvec ,ME_A,ME_B);// => recreate as p0plus_VAJHU/(p0plus_VAJHU + p1_VAJHU)
    //combinedMEM.computeKD(MEMNames::kSMHiggs, MEMNames::kJHUGen, MEMNames::k1plus,MEMNames::kJHUGen, &MEMs::probRatio,  MorindKD_vec ,ME_A,ME_B);// => recreate as p0plus_VAJHU/(p0plus_VAJHU + p1plus_VAJHU)
    //combinedMEM.computeKD(MEMNames::kSMHiggs, MEMNames::kJHUGen, MEMNames::k2mplus_gg, MEMNames::kJHUGen, &MEMs::probRatio,  MorindKD_gggrav ,ME_A,ME_B);// => recreate as  p0plus_VAJHU/(p0plus_VAJHU + p2_VAJHU)
    //combinedMEM.computeKD(MEMNames::kSMHiggs, MEMNames::kJHUGen, MEMNames::k2mplus_qqbar, MEMNames::kJHUGen, &MEMs::probRatio, MorindKD_qqgrav ,ME_A,ME_B);// => recreate as  p0plus_VAJHU/(p0plus_VAJHU + p2qqb_VAJHU)

    //double HCPKD_mela,HCPKD_grav ,HCPKD_pseudo;
    //HCP KDs 
    //combinedMEM.computeKD(kSMHiggs,kMELA_HCP,kqqZZ,kMELA_HCP,&MEMs::probRatio,HCPKD_mela  ,ME_A,ME_B);  // => recreate as p0plus_melaNorm/(p0plus_melaNorm + bkg_mela)
    //combinedMEM.computeKD(kSMHiggs,kAnalytical,k0minus,kAnalytical,&MEMs::probRatio,HCPKD_pseudo, ME_A,ME_B);// => recreate as p0plus_mela/(p0plus_mela + p0minus_mela)
    //combinedMEM.computeKD(kSMHiggs,kAnalytical,k2mplus_gg,kAnalytical,&MEMs::probRatio,HCPKD_grav  ,ME_A,ME_B);// => recreate as p0plus_mela/(p0plus_mela + p2_mela)
    
    //check syncronization
    //     std::cout << "melaHCP"  << HCPKD_mela << " : " << p0plus_melaNorm/(p0plus_melaNorm + bkg_mela) << " : " <<  ld  << std::endl;
    //     std::cout << "pseudoHCP"<< HCPKD_pseudo << " : " << p0plus_mela/(p0plus_mela + p0minus_mela)  <<  " : " <<  psld  <<std::endl;
    //     std::cout << "pseudoHCP"<< HCPKD_grav << " : " << p0plus_mela/(p0plus_mela + p2_mela)  <<  " : " <<  gravld  <<std::endl;
    //     std::cout << "MorindKD_mela " <<  MorindKD_mela<< " : " <<  p0plus_VAJHU/(p0plus_VAJHU + bkg_VAMCFMNorm)   << std::endl;
    //     std::cout << "MorindKD_pseudo " << MorindKD_pseudo << " : " <<  p0plus_VAJHU/(p0plus_VAJHU + p0minus_VAJHU) << std::endl;
    //     std::cout << "MorindKD_highdim " << MorindKD_highdim << " : " << p0plus_VAJHU/(p0plus_VAJHU + p0hplus_VAJHU)  << std::endl;
    //     std::cout << "MorindKD_vec " << MorindKD_vec << " : " << p0plus_VAJHU/(p0plus_VAJHU + p1plus_VAJHU)<< std::endl;
    //     std::cout << "MorindKD_psvec " << MorindKD_psvec << " : " << p0plus_VAJHU/(p0plus_VAJHU + p1_VAJHU)<< std::endl;
    //     std::cout << "MorindKD_gggrav " <<  MorindKD_gggrav << " : " << p0plus_VAJHU/(p0plus_VAJHU + p2_VAJHU)<< std::endl;
    //     std::cout << "MorindKD_qqgrav " << MorindKD_qqgrav << " : " << p0plus_VAJHU/(p0plus_VAJHU + p2qqb_VAJHU)<< std::endl;
    //     std::cout << std::endl;


    //supermela signal
    // m4l probability as in datacards
    combinedMEM.computePm4l(partP,partId,
			    MEMNames::kNone,
			    p0plus_m4l,
			    bkg_m4l);
    // m4l probability for systematics
    combinedMEM.computePm4l(partP,partId,
			    MEMNames::kScaleUp,
			    p0plus_m4l_ScaleUp, 
			    bkg_m4l_ScaleUp); 
    // m4l probability for systematics
    combinedMEM.computePm4l(partP,partId,
			    MEMNames::kScaleDown,
			    p0plus_m4l_ScaleDown,
			    bkg_m4l_ScaleDown);
    // m4l probability for systematics
    combinedMEM.computePm4l(partP,partId,
			    MEMNames::kResolUp,
			    p0plus_m4l_ResUp,
			    bkg_m4l_ResUp);
    // m4l probability for systematics
    combinedMEM.computePm4l(partP,partId,
			    MEMNames::kResolDown,
			    p0plus_m4l_ResDown,
			    bkg_m4l_ResDown);


    // spinMELA
    double pg1g4_mela, pg1g4_VAJHU, pg1g2_mela, pg1g2_VAJHU, pg1g4_pi2_VAJHU, pg1g2_pi2_VAJHU;
    combinedMEM.computeME_Interference(MEMNames::kg1g4,     MEMNames::kAnalytical, partP, partId, pg1g4_mela);
    combinedMEM.computeME_Interference(MEMNames::kg1g4,     MEMNames::kJHUGen,     partP, partId, pg1g4_VAJHU);
    combinedMEM.computeME_Interference(MEMNames::kg1g4_pi_2,MEMNames::kJHUGen,     partP, partId, pg1g4_pi2_VAJHU);
    combinedMEM.computeME_Interference(MEMNames::kg1g2_pi_2,MEMNames::kJHUGen,     partP, partId, pg1g2_pi2_VAJHU);
    combinedMEM.computeME_Interference(MEMNames::kg1g2,     MEMNames::kAnalytical, partP, partId, pg1g2_mela);
    combinedMEM.computeME_Interference(MEMNames::kg1g2,     MEMNames::kJHUGen,     partP, partId, pg1g2_VAJHU);     

    // VBF jets
    vector<const cmg::PFJet*> cleanedJets;
    VBFCandidateJetSelector myVBFCandidateJetSelector;
    cleanedJets = myVBFCandidateJetSelector.cleanJets(myCand,pfjetscoll,sampleType); //   //FIXME: should use LEPTON_SETUP instead of sampleType for mela and combinedMEM. This will be an issue for samples rescaled to different sqrts (none at present)
    vector<const cmg::PFJet*> cleanedJetsPt30;
    for (unsigned int i=0; i < cleanedJets.size(); ++i){
      const cmg::PFJet& myjet = *(cleanedJets.at(i));  
      if (myjet.pt()>30) cleanedJetsPt30.push_back(&myjet);
    }
    float detajj =-99.f;
    float mjj  =-99.f;
    float VD   =-99.f;
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

    if (cleanedJetsPt30.size()>=2) {
      detajj = cleanedJetsPt30[0]->eta()-cleanedJetsPt30[1]->eta();
      mjj = (cleanedJetsPt30[0]->p4()+cleanedJetsPt30[1]->p4()).M();
      VD = fisher(mjj,detajj);

      //Add dijet production probabilities
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
      partIdprod.push_back(1);
      partIdprod.push_back(1);

      if(cleanedJetsPt30.size()==2){
	for(int jecnum=0;jecnum<3;jecnum++){
	  Float_t ratio1=1.;
	  Float_t ratio2=1.;
	  if(jecnum==1){
	    ratio1=1 + cleanedJetsPt30[0]->uncOnFourVectorScale();
	    ratio2=1 + cleanedJetsPt30[1]->uncOnFourVectorScale();
	  }
	  if(jecnum==2){
	    ratio1=1 - cleanedJetsPt30[0]->uncOnFourVectorScale();
	    ratio2=1 - cleanedJetsPt30[1]->uncOnFourVectorScale();
	  }
	  TLorentzVector jet1(cleanedJetsPt30[0]->p4().x(),cleanedJetsPt30[0]->p4().y(),cleanedJetsPt30[0]->p4().z(),cleanedJetsPt30[0]->p4().t());
	  TLorentzVector jet2(cleanedJetsPt30[1]->p4().x(),cleanedJetsPt30[1]->p4().y(),cleanedJetsPt30[1]->p4().z(),cleanedJetsPt30[1]->p4().t());
	  jet1.SetPtEtaPhiM(jet1.Pt()*ratio1,jet1.Eta(),jet1.Phi(),jet1.M()*ratio1);
	  jet2.SetPtEtaPhiM(jet2.Pt()*ratio2,jet2.Eta(),jet2.Phi(),jet2.M()*ratio2);
	  partPprod.push_back(jet1);
	  partPprod.push_back(jet2);
	  double phjj_temp,pvbf_temp;
	  combinedMEM.computeME(MEMNames::kJJ_SMHiggs_GG,MEMNames::kJHUGen,partPprod,partIdprod,phjj_temp);
	  combinedMEM.computeME(MEMNames::kJJ_SMHiggs_VBF,MEMNames::kJHUGen,partPprod,partIdprod,pvbf_temp);
	  if(jecnum==0){
	    phjj_VAJHU_old=phjj_temp;
	    pvbf_VAJHU_old=pvbf_temp;
	    phjj_VAJHU_new=phjj_VAJHU_old;
	    pvbf_VAJHU_new=pvbf_VAJHU_old;
	  }
	  if(jecnum==1){
	    phjj_VAJHU_old_up=phjj_temp;
	    pvbf_VAJHU_old_up=pvbf_temp;
	    phjj_VAJHU_new_up=phjj_VAJHU_old_up;
	    pvbf_VAJHU_new_up=pvbf_VAJHU_old_up;
	  }
	  if(jecnum==2){
	    phjj_VAJHU_old_dn=phjj_temp;
	    pvbf_VAJHU_old_dn=pvbf_temp;
	    phjj_VAJHU_new_dn=phjj_VAJHU_old_dn;
	    pvbf_VAJHU_new_dn=pvbf_VAJHU_old_dn;
	  }
	  partPprod.pop_back();
	  partPprod.pop_back();
	}
      }
      if(cleanedJetsPt30.size()>2){
	TLorentzVector jet1,jet2;
	double djet_max=-1.;
	double phjj_temp=-1.;
	double pvbf_temp=-1.;
	double phjj_VAJHU_old_temp =-1.;
	double pvbf_VAJHU_old_temp =-1.;
	double phjj_VAJHU_new_temp =-1.;
	double pvbf_VAJHU_new_temp =-1.;
	for(int jecnum=0;jecnum<3;jecnum++){
	  Float_t ratio1=1.;
	  Float_t ratio2=1.;
	  if(jecnum==1){
	    ratio1=1 + cleanedJetsPt30[0]->uncOnFourVectorScale();
	    ratio2=1 + cleanedJetsPt30[1]->uncOnFourVectorScale();
	  }
	  if(jecnum==2){
	    ratio1=1 - cleanedJetsPt30[0]->uncOnFourVectorScale();
	    ratio2=1 - cleanedJetsPt30[1]->uncOnFourVectorScale();
	  }
	  jet1.SetXYZT(cleanedJetsPt30[0]->p4().x(),cleanedJetsPt30[0]->p4().y(),cleanedJetsPt30[0]->p4().z(),cleanedJetsPt30[0]->p4().t());
	  jet2.SetXYZT(cleanedJetsPt30[1]->p4().x(),cleanedJetsPt30[1]->p4().y(),cleanedJetsPt30[1]->p4().z(),cleanedJetsPt30[1]->p4().t());
	  jet1.SetPtEtaPhiM(jet1.Pt()*ratio1,jet1.Eta(),jet1.Phi(),jet1.M()*ratio1);
	  jet2.SetPtEtaPhiM(jet2.Pt()*ratio2,jet2.Eta(),jet2.Phi(),jet2.M()*ratio2);
	  partPprod.push_back(jet1);
	  partPprod.push_back(jet2);
	  combinedMEM.computeME(MEMNames::kJJ_SMHiggs_GG,MEMNames::kJHUGen,partPprod,partIdprod,phjj_VAJHU_old_temp);
	  combinedMEM.computeME(MEMNames::kJJ_SMHiggs_VBF,MEMNames::kJHUGen,partPprod,partIdprod,pvbf_VAJHU_old_temp);
	  partPprod.pop_back();
	  partPprod.pop_back();
	  for(unsigned int firstjet=0;firstjet<cleanedJetsPt30.size();firstjet++){
	    for(unsigned int secondjet=0;secondjet<cleanedJetsPt30.size();secondjet++){
	      if(secondjet<=firstjet) continue;
	      ratio1=1;
	      ratio2=1;
	      if(jecnum==1){
		ratio1=1 + cleanedJetsPt30[firstjet]->uncOnFourVectorScale();
		ratio2=1 + cleanedJetsPt30[secondjet]->uncOnFourVectorScale();
	      }
	      if(jecnum==2){
		ratio1=1 - cleanedJetsPt30[firstjet]->uncOnFourVectorScale();
		ratio2=1 - cleanedJetsPt30[secondjet]->uncOnFourVectorScale();
	      }
	      jet1.SetXYZT(cleanedJetsPt30[firstjet]->p4().x(),cleanedJetsPt30[firstjet]->p4().y(),cleanedJetsPt30[firstjet]->p4().z(),cleanedJetsPt30[firstjet]->p4().t());
	      jet2.SetXYZT(cleanedJetsPt30[secondjet]->p4().x(),cleanedJetsPt30[secondjet]->p4().y(),cleanedJetsPt30[secondjet]->p4().z(),cleanedJetsPt30[secondjet]->p4().t());
	      jet1.SetPtEtaPhiM(jet1.Pt()*ratio1,jet1.Eta(),jet1.Phi(),jet1.M()*ratio1);
	      jet2.SetPtEtaPhiM(jet2.Pt()*ratio2,jet2.Eta(),jet2.Phi(),jet2.M()*ratio2);
	      partPprod.push_back(jet1);
	      partPprod.push_back(jet2);
	      phjj_temp=-1.;
	      pvbf_temp=-1.;
	      combinedMEM.computeME(MEMNames::kJJ_SMHiggs_GG,MEMNames::kJHUGen,partPprod,partIdprod,phjj_temp);
	      combinedMEM.computeME(MEMNames::kJJ_SMHiggs_VBF,MEMNames::kJHUGen,partPprod,partIdprod,pvbf_temp);
	      double djet_temp= pvbf_temp/(pvbf_temp+phjj_temp);
	      if(djet_temp > djet_max){
		phjj_VAJHU_new_temp=phjj_temp;
		pvbf_VAJHU_new_temp=pvbf_temp;
		djet_max=djet_temp;
	      }
	      partPprod.pop_back();
	      partPprod.pop_back();
	    }
	  }
	  //HERE
	  if(jecnum==0){
	    phjj_VAJHU_old=phjj_VAJHU_old_temp;
	    pvbf_VAJHU_old=pvbf_VAJHU_old_temp;
	    phjj_VAJHU_new=phjj_VAJHU_new_temp;
	    pvbf_VAJHU_new=pvbf_VAJHU_new_temp;
	  }
	  if(jecnum==1){
	    phjj_VAJHU_old_up=phjj_VAJHU_old_temp;
	    pvbf_VAJHU_old_up=pvbf_VAJHU_old_temp;
	    phjj_VAJHU_new_up=phjj_VAJHU_new_temp;
	    pvbf_VAJHU_new_up=pvbf_VAJHU_new_temp;
	  }
	  if(jecnum==2){
	    phjj_VAJHU_old_dn=phjj_VAJHU_old_temp;
	    pvbf_VAJHU_old_dn=pvbf_VAJHU_old_temp;
	    phjj_VAJHU_new_dn=phjj_VAJHU_new_temp;
	    pvbf_VAJHU_new_dn=pvbf_VAJHU_new_temp;
	  }
	}
      }
    }

    // Old-style
    float mZa = (Z1Lp->p4()+Z2Lm->p4()).mass();
    float mZb = (Z1Lm->p4()+Z2Lp->p4()).mass();
    if (mZa<mZb) swap(mZa,mZb); //sorting

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
    
  
    // Sorted pTs
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
  
    //      float massError2=0;
    float massErrorCorr=0;
    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
	const reco::Candidate* l=cand->daughter(i)->daughter(j);
	const TH2F* h;
	if (l->isMuon()) h = corrSigmaMu;
	else             h = corrSigmaEle;
	int ptBin  = min(max(1,h->GetXaxis()->FindBin(l->pt())), h->GetNbinsX());
	int etaBin = min(max(1,h->GetYaxis()->FindBin(fabs(l->eta()))), h->GetNbinsY());
	float ebecorr = h->GetBinContent(ptBin, etaBin);
	//	  cout << "ebe: " << l->isMuon() << " " << l->isElectron() << " " << ebecorr << endl;
	
	//	  massError2+=sqr(sigma[i][j]);
	massErrorCorr+= (sigma[i][j]*ebecorr)*(sigma[i][j]*ebecorr);
      }
    }
    
    
    massErrorCorr += (sigma[0][2])*(sigma[0][2]);
    massErrorCorr += (sigma[1][2])*(sigma[1][2]);
    massErrorCorr = sqrt(massErrorCorr);
    myCand.addUserFloat("massErrorCorr",      massErrorCorr);
  
    
    //--- 4l vertex fits
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
	  // 		       << sqrt(allLeptonsCand->currentState().kinematicParametersError().matrix()(6,6) ) << endl;
	  
	  reco::Vertex constrainedVertex(reco::Vertex::Point(allLeptonsVertex->position()),
					 allLeptonsVertex->error().matrix_new(), 
					 allLeptonsVertex->chiSquared(), 
					 allLeptonsVertex->degreesOfFreedom(),0);	
	  
	  // 	  if(dbg) cout << "kinematicFit vertex, ndof, chi2, prob: " 
	  // 		       << allLeptonsVertex->position() << " , " 
	  // 		       << allLeptonsVertex->degreesOfFreedom() << " , "
	  // 		       << allLeptonsVertex->chiSquared()   << " , "
	  // 		       << TMath::Prob(allLeptonsVertex->chiSquared(),allLeptonsVertex->degreesOfFreedom()) << endl;
	  
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
  
    myCand.addUserFloat("SIP4",           SIP4);
    myCand.addUserFloat("pt1",            ptS[3]); // leading-pT
    myCand.addUserFloat("pt2",            ptS[2]); // sub-leading pT
    myCand.addUserFloat("iso34",          iso34);
    myCand.addUserFloat("mZa",            mZa);
    myCand.addUserFloat("mZb",            mZb);
    myCand.addUserFloat("mLL4",           mll4); // smallest mass of any AF/OS pair
    myCand.addUserFloat("mLL6",           mll6);   // smallest mass of any AF/AS pair
    myCand.addUserFloat("costheta1",      costheta1);
    myCand.addUserFloat("costheta2",      costheta2);
    myCand.addUserFloat("phi",            phi);
    myCand.addUserFloat("costhetastar",   costhetastar);
    myCand.addUserFloat("phistar1",       phistar1);
    myCand.addUserFloat("xistar",         xistar);  //azimuthal angle of higgs in rest frame of higgs
    myCand.addUserFloat("xi",             xi);      //azimuthal angle of higgs in lab frame

    //    myCand.addUserFloat("phistar2",       phistar2);
    // outdated ... AJW
    //myCand.addUserFloat("LD",             ld);
    //myCand.addUserFloat("PSig",           psig);
    //myCand.addUserFloat("PBkg",           pbkg);
    //myCand.addUserFloat("pseudoLD",       psld);
    //myCand.addUserFloat("PSigPS",         psigPS);
    //myCand.addUserFloat("spin2PMLD",      gravld);
    //myCand.addUserFloat("PSig2PM",        psig2PM);
    if(doMEKD){
      myCand.addUserFloat("MEKD_LD", mekd_ld);
      myCand.addUserFloat("MEKD_PseudoLD", mekd_pseudold);
      myCand.addUserFloat("MEKD_GravLD", mekd_gravld);
    }
    myCand.addUserFloat("m4l",            (Z1Lm->p4()+Z1Lp->p4()+Z2Lm->p4()+Z2Lp->p4()).mass()); // mass without FSR
    if (doKinFit && abs(id11)==13 && abs(id12)==13 && abs(id21)==13 && abs(id22)==13 ) {
      myCand.addUserFloat("m4lRef",  m4lRef ); // mass from Z1 refitted (FSR not considered in the refit procedure)
      myCand.addUserFloat("chi2Fit", chi2LepZ1Ref);
      myCand.addUserFloat("mZ1Ref",  mZ1Ref);
    }
    // add probabilities

    /*   KDs for analytical calculations
    myCand.addUserFloat("p0plus_melaNorm",p0plus_melaNorm);// higgs, analytic distribution, normalized as for normal MELA distribution
    myCand.addUserFloat("p0plus_mela",    p0plus_mela);    // higgs, analytic distribution
    myCand.addUserFloat("p0minus_mela",   p0minus_mela);   // pseudoscalar, analytic distribution
    myCand.addUserFloat("p0hplus_mela",   p0hplus_mela);// 0h+, analytic distribution
    myCand.addUserFloat("p1_mela",        p1_mela);        // zprime, analytic distribution
    myCand.addUserFloat("p1plus_mela",    p1plus_mela);// 1+, analytic distribution 
    myCand.addUserFloat("p2_mela",        p2_mela);        // graviton, analytic distribution
    myCand.addUserFloat("p2qqb_mela",     p2qqb_mela);// graviton produced by qqbar vector algebra, analytical,
    myCand.addUserFloat("bkg_mela",       bkg_mela);       // background,  analytic distribution
    myCand.addUserFloat("p1_prodIndep_mela",        p1_prodIndep_mela);        // zprime, analytic distribution
    myCand.addUserFloat("p1plus_prodIndep_mela",    p1plus_prodIndep_mela);// 1+, analytic distribution 
    myCand.addUserFloat("p2_prodIndep_mela",        p2_prodIndep_mela);        // graviton, analytic distribution
    myCand.addUserFloat("p2hplus_mela",    p2hplus_mela);       // graviton, vector algebra, analytical,
    myCand.addUserFloat("p2hminus_mela",    p2hminus_mela);       // graviton, vector algebra, analytical,
    myCand.addUserFloat("p2bplus_mela",    p2bplus_mela);       // graviton, vector algebra, analytical,
    */

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

    // VBF
    myCand.addUserFloat("detajj",  detajj);
    myCand.addUserFloat("mjj",     mjj);
    myCand.addUserFloat("VD",      VD);
    myCand.addUserFloat("phjj_VAJHU_old",phjj_VAJHU_old);
    myCand.addUserFloat("pvbf_VAJHU_old",pvbf_VAJHU_old);
    myCand.addUserFloat("phjj_VAJHU_old_up",phjj_VAJHU_old_up);
    myCand.addUserFloat("pvbf_VAJHU_old_up",pvbf_VAJHU_old_up);
    myCand.addUserFloat("phjj_VAJHU_old_dn",phjj_VAJHU_old_dn);
    myCand.addUserFloat("pvbf_VAJHU_old_dn",pvbf_VAJHU_old_dn);
    myCand.addUserFloat("phjj_VAJHU_new",phjj_VAJHU_new);
    myCand.addUserFloat("pvbf_VAJHU_new",pvbf_VAJHU_new);
    myCand.addUserFloat("phjj_VAJHU_new_up",phjj_VAJHU_new_up);
    myCand.addUserFloat("pvbf_VAJHU_new_up",pvbf_VAJHU_new_up);
    myCand.addUserFloat("phjj_VAJHU_new_dn",phjj_VAJHU_new_dn);
    myCand.addUserFloat("pvbf_VAJHU_new_dn",pvbf_VAJHU_new_dn);


    //--- MC matching
    if (isMC) {
      int refID = 25; // FIXME: handle ZZ (sigId = 23)
      bool MC_isRight = (myCand.userFloat("d0.d0.MCParentCode")==refID &&
			 myCand.userFloat("d0.d1.MCParentCode")==refID &&
			 myCand.userFloat("d1.d0.MCParentCode")==refID &&
			 myCand.userFloat("d1.d1.MCParentCode")==refID);
      bool MC_isRightPair = false; //FIXME to be 

      myCand.addUserFloat("MC_isRight",     MC_isRight);
      myCand.addUserFloat("MC_isRightPair", MC_isRightPair);
    }
  

    
    
    //--- Select the "best candidate" among those passing the "bestCandAmong" selection (2011 PRL logic)
    int iCRname=0;
    for(CutSet<pat::CompositeCandidate>::const_iterator bca = preBestCandSelection.begin(); bca != preBestCandSelection.end(); ++bca){
      int preBestCandResult= int((*(bca->second))(myCand));

      if (preBestCandResult){
	// For the SR, we want the Z1 to be the best among ee and mumu
	bool hasBestZ = true;
	if ((bca->first)=="isBestCand" &&  Z1->daughter(0)->pdgId()!=bestZFlavour) hasBestZ = false;

	// Look for the daughter with largesr-pT leptons
	float ptSum = Z2->daughter(0)->pt()+Z2->daughter(1)->pt();
	if (hasBestZ && ptSum > maxPtSum[iCRname]){
	  maxPtSum[iCRname] = ptSum;
	  bestCandIdx[iCRname] = i;
	}
      }
      iCRname++;
    }
    result->push_back(myCand);
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
  

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ZZCandidateFiller);

