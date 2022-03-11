/** \class TauFiller
 *
 *  No description available.
 *
 *  $Date: 2013/05/24 15:42:42 $
 *  $Revision: 1.28 $
 *  \author G. Ortona (LLR)
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
//#include <DataFormats/TauReco/interface/PFTauDiscriminator.h>

//#include "DataFormats/VertexReco/interface/Vertex.h"
#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>
#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>
//#include "BDTId.h"

#include <vector>
#include <string>

#include <TFile.h>             // TFile
#include <TH1.h>               // TH1
#include <TGraphAsymmErrors.h> // TGraphAsymmErrors

using namespace edm;
using namespace std;
using namespace reco;

//bool recomputeBDT = false;

class TauFiller : public edm::EDProducer {
 public:
  /// Constructor
  explicit TauFiller(const edm::ParameterSet&);
    
  /// Destructor
  ~TauFiller();

  //ByIsolationMVA3oldDMwoLTraw
 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  edm::EDGetTokenT<pat::TauRefVector> theCandidateTag;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > theGenTag ;
  edm::EDGetTokenT<vector<Vertex> > theVtxTag ;
  const std::string theDiscriminatorTag;
  const StringCutObjectSelector<pat::Tau, true> cut;
  const CutSet<pat::Tau> flags;
  const bool ApplyTESCentralCorr; // shift the central TES value
  const std::string theTESYear;
  TFile* TESFile;
  TFile* EESFile;
  TH1* TESh1;
  TGraphAsymmErrors* EESgr;

  vector<string> tauIntDiscrims_; // tau discrims to be added as userInt
  vector<string> tauFloatDiscrims_; // tau discrims to be added as userFloats
};


TauFiller::TauFiller(const edm::ParameterSet& iConfig) :
  theCandidateTag(consumes<pat::TauRefVector>(iConfig.getParameter<InputTag>("src"))),
  theGenTag(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genCollection"))),
  theVtxTag(consumes<vector<Vertex>>(iConfig.getParameter<edm::InputTag>("vtxCollection"))),
  theDiscriminatorTag(iConfig.getParameter<std::string>("discriminator")),
  cut(iConfig.getParameter<std::string>("cut")),
  flags(iConfig.getParameter<ParameterSet>("flags")), 
  ApplyTESCentralCorr(iConfig.getParameter<bool>("ApplyTESCentralCorr")),
  theTESYear(iConfig.getParameter<std::string>("year"))

{
  produces<pat::TauCollection>();

  tauIntDiscrims_ = 
  {
    "decayModeFinding", // it is decayModeFindingOldDMs
    "decayModeFindingNewDMs",
    
    "byLooseCombinedIsolationDeltaBetaCorr3Hits",
    "byMediumCombinedIsolationDeltaBetaCorr3Hits",
    "byTightCombinedIsolationDeltaBetaCorr3Hits",
    
    //"byVLooseIsolationMVArun2v1DBoldDMwLT",
    //"byLooseIsolationMVArun2v1DBoldDMwLT",
    //"byMediumIsolationMVArun2v1DBoldDMwLT",
    //"byTightIsolationMVArun2v1DBoldDMwLT",
    //"byVTightIsolationMVArun2v1DBoldDMwLT",

    //"byVLooseIsolationMVArun2v1DBnewDMwLT",
    //"byLooseIsolationMVArun2v1DBnewDMwLT",
    //"byMediumIsolationMVArun2v1DBnewDMwLT",
    //"byTightIsolationMVArun2v1DBnewDMwLT",
    //"byVTightIsolationMVArun2v1DBnewDMwLT",

    //"byLooseIsolationMVArun2v1DBdR03oldDMwLT",
    //"byMediumIsolationMVArun2v1DBdR03oldDMwLT",
    //"byTightIsolationMVArun2v1DBdR03oldDMwLT",
    //"byVTightIsolationMVArun2v1DBdR03oldDMwLT",
    
    "againstMuonLoose3",
    "againstMuonTight3",

    "againstElectronVLooseMVA6",
    "againstElectronLooseMVA6",
    "againstElectronMediumMVA6",
    "againstElectronTightMVA6",
    "againstElectronVTightMVA6",

    "byVVLooseIsolationMVArun2017v2DBoldDMwLT2017", //FRA syncApr2018
    "byVLooseIsolationMVArun2017v2DBoldDMwLT2017", //FRA syncApr2018
    "byLooseIsolationMVArun2017v2DBoldDMwLT2017",  //FRA syncApr2018
    "byMediumIsolationMVArun2017v2DBoldDMwLT2017", //FRA syncApr2018
    "byTightIsolationMVArun2017v2DBoldDMwLT2017",  //FRA syncApr2018
    "byVTightIsolationMVArun2017v2DBoldDMwLT2017", //FRA syncApr2018
    "byVVTightIsolationMVArun2017v2DBoldDMwLT2017", //FRA syncApr2018
    
    //"byVLooseIsolationMVArun2017v1DBoldDMwLT2017", //FRA syncApr2018
    //"byLooseIsolationMVArun2017v1DBoldDMwLT2017",  //FRA syncApr2018
    //"byMediumIsolationMVArun2017v1DBoldDMwLT2017", //FRA syncApr2018
    //"byTightIsolationMVArun2017v1DBoldDMwLT2017",  //FRA syncApr2018
    //"byVTightIsolationMVArun2017v1DBoldDMwLT2017", //FRA syncApr2018
    
    //"byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017", //FRA syncApr2018
    //"byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017",  //FRA syncApr2018
    //"byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017", //FRA syncApr2018
    //"byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017",  //FRA syncApr2018
    //"byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017", //FRA syncApr2018
    
    "byVVVLooseDeepTau2017v2p1VSjet",
    "byVVLooseDeepTau2017v2p1VSjet", 
    "byVLooseDeepTau2017v2p1VSjet",  
    "byLooseDeepTau2017v2p1VSjet",   
    "byMediumDeepTau2017v2p1VSjet",  
    "byTightDeepTau2017v2p1VSjet",   
    "byVTightDeepTau2017v2p1VSjet",  
    "byVVTightDeepTau2017v2p1VSjet", 
    
    "byVVVLooseDeepTau2017v2p1VSe",  
    "byVVLooseDeepTau2017v2p1VSe", 
    "byVLooseDeepTau2017v2p1VSe",   
    "byLooseDeepTau2017v2p1VSe",   
    "byMediumDeepTau2017v2p1VSe",   
    "byTightDeepTau2017v2p1VSe",   
    "byVTightDeepTau2017v2p1VSe",   
    "byVVTightDeepTau2017v2p1VSe",   
    
    "byVLooseDeepTau2017v2p1VSmu",
    "byLooseDeepTau2017v2p1VSmu", 
    "byMediumDeepTau2017v2p1VSmu",
    "byTightDeepTau2017v2p1VSmu", 

  };

  tauFloatDiscrims_ =
  {
    "byCombinedIsolationDeltaBetaCorrRaw3Hits",
    "byIsolationMVArun2v1DBoldDMwLTraw",
    "byPhotonPtSumOutsideSignalCone",
    "footprintCorrection",
    "neutralIsoPtSumWeight",
    "photonPtSumOutsideSignalCone",
    "chargedIsoPtSum",
    "neutralIsoPtSum",
    "puCorrPtSum",
    //"byIsolationMVArun2017v1DBoldDMwLTraw2017",      //FRA syncApr2018
    "byIsolationMVArun2017v2DBoldDMwLTraw2017",      //FRA syncApr2018
    //"byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017", //FRA syncApr2018
    "byDeepTau2017v2p1VSjetraw",  
    "byDeepTau2017v2p1VSeraw",  
    "byDeepTau2017v2p1VSmuraw",
   };

  // TES input files
  edm::FileInPath TESFileName ("TauPOG/TauIDSFs/data/TauES_dm_DeepTau2017v2p1VSjet_"+theTESYear+".root");
  TESFile  = new TFile(TESFileName.fullPath().data());

  // TES input histos
  TH1::AddDirectory(false);
  TESh1 = dynamic_cast<TH1*>((const_cast<TFile*>(TESFile))->Get("tes"));
  
  // Davide: Use Legacy EES for the time being as sugested by tau POG

  std::string theEESYear = theTESYear;//"";
  
  //if      (theTESYear == "UL2016_preVFP" || theTESYear == "UL2016_postVFP") theEESYear = "2016Legacy";
  //else if (theTESYear == "UL2017")                                          theEESYear = "2017ReReco";
  //else if (theTESYear == "UL2018")                                          theEESYear = "2018ReReco";
  //else                                                                      theEESYear = "2016Legacy";
  
  // EES input file
  edm::FileInPath EESFileName("TauPOG/TauIDSFs/data/TauFES_eta-dm_DeepTau2017v2p1VSe_"+theEESYear+".root");
  EESFile = new TFile(EESFileName.fullPath().data());

  // EES input histos
  EESgr = dynamic_cast<TGraphAsymmErrors*>((const_cast<TFile*>(EESFile))->Get("fes"));
}

TauFiller::~TauFiller()
{
  delete TESFile;
  delete TESh1;
  delete EESFile;
  delete EESgr;
}

using LorentzVectorE = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double>>;

// Return original Tau Mother genParticle:
// given a genParticle, go back in the chain of the gen particles until
// you find the last tau (i.e. the tau that has as mother not a tau)
const reco::GenParticle getMother (const reco::GenParticle genP)
{
  //std::cout << "  * first genPart - pdg: " << genP.pdgId() << " - isPrompt: " << genP.statusFlags().isPrompt() << " - momentum: " << genP.momentum() << endl;
  reco::GenParticleRef genM = genP.motherRef(0);
  assert(genM.isNonnull() && genM.isAvailable());  // sanity
  if (std::abs(genP.pdgId())==15 && std::abs(genM->pdgId())!=15)
  {
    //std::cout << "    returning daughter  - pdg: " << genP.pdgId() << " - isPrompt: " << genP.statusFlags().isPrompt() << " - momentum: " << genP.momentum() << endl;
    return genP;
  }
  else
  {
    //std::cout << "    retrying with mother  - pdg: " << genM->pdgId() << " - isPrompt: " << genM->statusFlags().isPrompt() << " - momentum: " << genM->momentum() << endl;
    return getMother(*genM);
  }
}

void
TauFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  

//read one PFTauDiscriminator (set discriminatorSrc_ in to an edm::InputTag before)

  // Get leptons and discriminators
  edm::Handle<pat::TauRefVector> tauHandle;
  iEvent.getByToken(theCandidateTag, tauHandle);
    
  edm::Handle<vector<Vertex> >  vertexs;
  iEvent.getByToken(theVtxTag, vertexs);

  edm::Handle<edm::View<reco::GenParticle> > genHandle;
  iEvent.getByToken(theGenTag, genHandle);

  // TES corrections
  Int_t binDM0  = TESh1->GetXaxis()->FindBin((int)0);
  Int_t binDM1  = TESh1->GetXaxis()->FindBin(1);
  Int_t binDM10 = TESh1->GetXaxis()->FindBin(10);
  Int_t binDM11 = TESh1->GetXaxis()->FindBin(11);

  double Shift1Pr    = TESh1->GetBinContent(binDM0);
  double Shift1PrPi0 = TESh1->GetBinContent(binDM1);
  double Shift3Pr    = TESh1->GetBinContent(binDM10);
  double Shift3PrPi0 = TESh1->GetBinContent(binDM11);

  // EES corrections
  double EFakeShift1PrB    = EESgr->GetY()[0]; // barrel DM 0
  double EFakeShift1PrPi0B = EESgr->GetY()[1]; // barrel DM 1
  double EFakeShift1PrE    = EESgr->GetY()[2]; // endcap DM 0
  double EFakeShift1PrPi0E = EESgr->GetY()[3]; // endcap DM 1

  // Output collection
  //auto_ptr<pat::TauCollection> result( new pat::TauCollection() );
  std::unique_ptr<pat::TauCollection> result( new pat::TauCollection() );

  for (unsigned int itau = 0; itau < tauHandle->size(); ++itau){
    //cout << "- TauFiller itau: " << itau << endl;

    //---Clone the pat::Tau
    pat::Tau l(*((*tauHandle)[itau].get()));
    
    double shiftP = 1.;
    double shiftMass = 1.;
    bool isTESShifted = false;
    bool isTauMatched = false;
    bool isEESShifted = false;

    // Check if the original Tau Mother genParticle isPrompt
    bool isTauPrompt = false;
    if ( l.genJet())
    {
      // Get the constituents of the genJet
      std::vector<const reco::GenParticle*> genParts = l.genJet()->getGenConstituents();

      // Get the original tau mother starting from one of the constituents of the genJet,
      // it does not matther which one, so we use element '0' --> genParts[0]
      const reco::GenParticle returned = getMother(*(genParts[0]));
      //std::cout << "  Returned - pdg: " << returned.pdgId() << " - isPrompt: " << returned.statusFlags().isPrompt() << " - momentum: " << returned.momentum() << endl;

      // Check if the original tau mother isPrompt or not
      isTauPrompt = returned.statusFlags().isPrompt();
    }

    if ( l.genJet() && deltaR(l.p4(), l.genJet()->p4()) < 0.3 && l.genJet()->pt() > 15. && ((std::abs(l.genJet()->pt()-l.pt())/l.genJet()->pt()) < 1.0) && isTauPrompt && ApplyTESCentralCorr)
    {
      isTauMatched = true;
      isTESShifted = true;
      //cout << "---- gen get pt: " << l.genJet()->pt() << endl;
      if (l.decayMode()==0)       // 1prong
      {
        shiftP    = Shift1Pr;
        shiftMass = 1.;
      }
      else if (l.decayMode()==1)  // 1prong+pi0
      {
        shiftP    = Shift1PrPi0;
        shiftMass = Shift1PrPi0;
      }
      else if (l.decayMode()==10) // 3prong
      {
        shiftP    = Shift3Pr;
        shiftMass = Shift3Pr;
      }
      else if (l.decayMode()==11) // 3prong+pi0
      {
        shiftP    = Shift3PrPi0;
        shiftMass = Shift3PrPi0;
      }
      else  // these are not real taus and will be rejected --> we don't care about the shift and just put 1
      {
        isTESShifted = false;
        isTauMatched = false;
        shiftP    = 1.;
        shiftMass = 1.;
      }
      //if(l.decayMode()>=1 && l.decayMode()<=10){
      //  shiftP = Shift;
      //  shiftMass = Shift;
      //}
      //else if(l.decayMode()==0){
      //  shiftP = Shift;
      //  shiftMass = 1.;
      //}
    }

    int genmatch = 6; // 6 = fake
    if (isTauMatched)
    {
      genmatch = 5;
    }
    else // !isTauMatched
    {
      //https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2016#MC_Matching
      //https://github.com/KIT-CMS/Artus/blob/dictchanges/KappaAnalysis/src/Utility/GeneratorInfo.cc#L77-L165
      if (genHandle.isValid())
      {
        GenParticle closest = (GenParticle) (*genHandle)[0] ;
        double closestDR = 999;

        for (unsigned int iGen = 0; iGen < genHandle->size(); iGen++)
        {
          const GenParticle& genP = (*genHandle)[iGen];
          int genP_pdgId = std::abs(genP.pdgId());
          double deltaPTpT = std::abs( genP.pt() - l.pt() ) / genP.pt() ;
          if (genP.pt() > 8. && (genP_pdgId == 11 || genP_pdgId == 13) && (genP.statusFlags().isPrompt() || genP.statusFlags().isDirectPromptTauDecayProduct()) && deltaPTpT < 0.5)
          {
            double tmpDR = deltaR(l.p4(), genP.p4());
            if (tmpDR < closestDR)
            {
              closest = genP;
              closestDR = tmpDR;
            }
          }
        }

        if (closestDR < 0.3)
        {
          int pdgId = std::abs(closest.pdgId());
          if      (pdgId == 11 && closest.pt() > 8. && closest.statusFlags().isPrompt()) genmatch = 1;
          else if (pdgId == 13 && closest.pt() > 8. && closest.statusFlags().isPrompt()) genmatch = 2;
          else if (pdgId == 11 && closest.pt() > 8. && closest.statusFlags().isDirectPromptTauDecayProduct()) genmatch = 3;
          else if (pdgId == 13 && closest.pt() > 8. && closest.statusFlags().isDirectPromptTauDecayProduct()) genmatch = 4;
        }
      } // end genHandle.isValid()
    } // end !isTauMatched

    //E->tau ES
    if(ApplyTESCentralCorr)
    {
      if ((genmatch == 1 || genmatch == 3) && l.decayMode()==0)
      {
        shiftP = EFakeShift1PrB; // 1prong
        if (fabs(l.eta())> 1.5)
          shiftP = EFakeShift1PrE;
        shiftMass = 1.;
        isEESShifted = true;
      }
      if ((genmatch == 1 || genmatch == 3) && l.decayMode()==1)
      {
        shiftP    = EFakeShift1PrPi0B; // 1prong+pi0
	      shiftMass = EFakeShift1PrPi0B;
	      if (fabs(l.eta())> 1.5)
        {
	        shiftP    = EFakeShift1PrPi0E;
	        shiftMass = EFakeShift1PrPi0E;
	      }
	      isEESShifted = true;
      }
    }

    double pxS_Nominal = l.px()*shiftP;
    double pyS_Nominal = l.py()*shiftP;
    double pzS_Nominal = l.pz()*shiftP;
    double massS_Nominal = l.mass()*shiftMass;
    double enS_Nominal = TMath::Sqrt(pxS_Nominal*pxS_Nominal + pyS_Nominal*pyS_Nominal + pzS_Nominal*pzS_Nominal + massS_Nominal*massS_Nominal);
    math::XYZTLorentzVectorD p4S_Nominal( pxS_Nominal, pyS_Nominal, pzS_Nominal, enS_Nominal );


    //Up and Down variations: NominalTESUncertainty from python cfg
    //const float udShift[2] = {1.03, 0.97}; // 0: UP, 1: DOWN
    //const double udShift[2] = {1. + (NominalTESUncertainty/100.), 1. - (NominalTESUncertainty/100.)}; // 0: UP, 1: DOWN // unused

    float TESshiftDM0  = 1;
    float TESshiftDM1  = 1;
    float TESshiftDM10 = 1;
    float TESshiftDM11 = 1;

    double errDM0  = TESh1->GetBinError(binDM0);
    double errDM1  = TESh1->GetBinError(binDM1);
    double errDM10 = TESh1->GetBinError(binDM10);
    double errDM11 = TESh1->GetBinError(binDM11);

    TESshiftDM0  = errDM0 ;
    TESshiftDM1  = errDM1 ;
    TESshiftDM10 = errDM10;
    TESshiftDM11 = errDM11;

    float udshiftP[2] = {1., 1.};
    float udshiftMass[2] = {1., 1.};

    if(isTauMatched){

      if (l.decayMode()==0)       // 1prong
      {
        udshiftP[0]    =  Shift1Pr + TESshiftDM0; //udShift[0]; // up
        udshiftP[1]    =  Shift1Pr - TESshiftDM0; //udShift[1]; // down
        udshiftMass[0] = udshiftMass[1] = 1.; // no mass shift for pi0
        l.addUserFloat("TESshiftDM0",TESshiftDM0);
      }
      else if (l.decayMode()==1)  // 1prong+pi0
      {
        udshiftP[0]    = Shift1PrPi0 + TESshiftDM1; //udShift[0]; // up
        udshiftP[1]    = Shift1PrPi0 - TESshiftDM1; //udShift[1]; // down
        udshiftMass[0] = Shift1PrPi0 + TESshiftDM1; //udShift[0]; // up
        udshiftMass[1] = Shift1PrPi0 - TESshiftDM1; //udShift[1]; // down
        l.addUserFloat("TESshiftDM1",TESshiftDM1);
      }
      else if (l.decayMode()==10) // 3prong
      {
        udshiftP[0]    = Shift3Pr + TESshiftDM10; //udShift[0]; // up
        udshiftP[1]    = Shift3Pr - TESshiftDM10; //udShift[1]; // down
        udshiftMass[0] = Shift3Pr + TESshiftDM10; //udShift[0]; // up
        udshiftMass[1] = Shift3Pr - TESshiftDM10; //udShift[1]; // down
        l.addUserFloat("TESshiftDM10",TESshiftDM10);
      }
      else if (l.decayMode()==11) // 3prong
      {
        udshiftP[0]    = Shift3PrPi0 + TESshiftDM11; //udShift[0]; // up
        udshiftP[1]    = Shift3PrPi0 - TESshiftDM11; //udShift[1]; // down
        udshiftMass[0] = Shift3PrPi0 + TESshiftDM11; //udShift[0]; // up
        udshiftMass[1] = Shift3PrPi0 - TESshiftDM11; //udShift[1]; // down
        l.addUserFloat("TESshiftDM11",TESshiftDM11);
      }
      else  // these are not real taus and will be rejected --> we don't care about the shift and just put 1
      {
        isTESShifted = false;
        udshiftP[0]    = udshiftP[1]    = 1.;
        udshiftMass[0] = udshiftMass[1] = 1.;
      }
      //if(l.decayMode()>=1 && l.decayMode()<=10){
      //  udshiftP[0] = udShift[0]; // up
      //  udshiftP[1] = udShift[1]; // down
      //  udshiftMass[0] = udShift[0]; // up
      //  udshiftMass[1] = udShift[1]; // down
      //}
      //else if(l.decayMode()==0){
      //  udshiftP[0] = udShift[0]; // up
      //  udshiftP[1] = udShift[1]; // down
      //  udshiftMass[0] = udshiftMass[1] = 1.; // no mass shift for pi0
      //}
      //else isTESShifted = false;
    }

    if (ApplyTESCentralCorr && isTESShifted)
    {
      // up shift
      double pxS = l.px()*udshiftP[0];
      double pyS = l.py()*udshiftP[0];
      double pzS = l.pz()*udshiftP[0];
      double massS = l.mass()*udshiftMass[0];
      double enS = TMath::Sqrt(pxS*pxS + pyS*pyS + pzS*pzS + massS*massS);
      math::XYZTLorentzVectorD p4SUP( pxS, pyS, pzS, enS );
      // set userfloats
      l.addUserFloat("px_TauUp",p4SUP.px());
      l.addUserFloat("py_TauUp",p4SUP.py());
      l.addUserFloat("pz_TauUp",p4SUP.pz());
      l.addUserFloat("e_TauUp",p4SUP.energy());
      l.addUserFloat("m_TauUp",p4SUP.mass());

      // down shift
      pxS = l.px()*udshiftP[1];
      pyS = l.py()*udshiftP[1];
      pzS = l.pz()*udshiftP[1];
      massS = l.mass()*udshiftMass[1];
      enS = TMath::Sqrt(pxS*pxS + pyS*pyS + pzS*pzS + massS*massS);
      math::XYZTLorentzVectorD p4SDOWN( pxS, pyS, pzS, enS );
      // set userfloats
      l.addUserFloat("px_TauDown",p4SDOWN.px());
      l.addUserFloat("py_TauDown",p4SDOWN.py());
      l.addUserFloat("pz_TauDown",p4SDOWN.pz());
      l.addUserFloat("e_TauDown",p4SDOWN.energy());
      l.addUserFloat("m_TauDown",p4SDOWN.mass());
    }

    // Up Down shifts for e->tau fES
    float EESshiftDM0Bup   = EESgr->GetErrorYhigh(0); // barrel DM 0
    float EESshiftDM0Bdown = EESgr->GetErrorYlow (0);

    float EESshiftDM1Bup   = EESgr->GetErrorYhigh(1); // barrel DM 1
    float EESshiftDM1Bdown = EESgr->GetErrorYlow (1);

    float EESshiftDM0Eup   = EESgr->GetErrorYhigh(2); // endcap DM 0
    float EESshiftDM0Edown = EESgr->GetErrorYlow (2);

    float EESshiftDM1Eup   = EESgr->GetErrorYhigh(3); // endcap DM 1
    float EESshiftDM1Edown = EESgr->GetErrorYlow (3);

    float udEFakeshiftP[2] = {1., 1.};
    float udEFakeshiftMass[2] = {1., 1.};

    if ((genmatch == 1 || genmatch == 3) &&l.decayMode()==0)
    {
      udEFakeshiftP[0]    =  EFakeShift1PrB + EESshiftDM0Bup; // up
      udEFakeshiftP[1]    =  EFakeShift1PrB - EESshiftDM0Bdown; // down
      if(fabs(l.eta())>1.5)
      {
        udEFakeshiftP[0]    =  EFakeShift1PrE + EESshiftDM0Eup; // up
        udEFakeshiftP[1]    =  EFakeShift1PrE - EESshiftDM0Edown; // down
      }
      udEFakeshiftMass[0] = udEFakeshiftMass[1] = 1.; // no mass shift for pi0
      if(fabs(l.eta())<=1.5)
      {
        l.addUserFloat("EESshiftDM0up",EESshiftDM0Bup);
        l.addUserFloat("EESshiftDM0dw",EESshiftDM0Bdown);
      }
      else
      {
        l.addUserFloat("EESshiftDM0up",EESshiftDM0Eup);
        l.addUserFloat("EESshiftDM0dw",EESshiftDM0Edown);
      }
    }
    if ((genmatch == 1 || genmatch == 3) &&l.decayMode()==1)
    {
      udEFakeshiftP[0]    = EFakeShift1PrPi0B + EESshiftDM1Bup;   // up
      udEFakeshiftP[1]    = EFakeShift1PrPi0B - EESshiftDM1Bdown; // down
      udEFakeshiftMass[0] = EFakeShift1PrPi0B + EESshiftDM1Bup;   // up
      udEFakeshiftMass[1] = EFakeShift1PrPi0B - EESshiftDM1Bdown; // down
      if(fabs(l.eta())>1.5)
      {
        udEFakeshiftP[0]    = EFakeShift1PrPi0E + EESshiftDM1Eup;   // up
        udEFakeshiftP[1]    = EFakeShift1PrPi0E - EESshiftDM1Edown; // down
        udEFakeshiftMass[0] = EFakeShift1PrPi0E + EESshiftDM1Eup;   // up
        udEFakeshiftMass[1] = EFakeShift1PrPi0E - EESshiftDM1Edown; // down
      }
      if(fabs(l.eta())<=1.5)
      {
        l.addUserFloat("EESshiftDM1up",EESshiftDM1Bup);
        l.addUserFloat("EESshiftDM1dw",EESshiftDM1Bdown);
      }
      else
      {
        l.addUserFloat("EESshiftDM1up",EESshiftDM1Eup);
        l.addUserFloat("EESshiftDM1dw",EESshiftDM1Edown);
      }
    }

    if (ApplyTESCentralCorr && isEESShifted)
    {
      // up shift
      double pxEFakeS = l.px()*udEFakeshiftP[0];
      double pyEFakeS = l.py()*udEFakeshiftP[0];
      double pzEFakeS = l.pz()*udEFakeshiftP[0];
      double massEFakeS = l.mass()*udEFakeshiftMass[0];
      double enEFakeS = TMath::Sqrt(pxEFakeS*pxEFakeS + pyEFakeS*pyEFakeS + pzEFakeS*pzEFakeS + massEFakeS*massEFakeS);
      math::XYZTLorentzVectorD p4EFakeSUP( pxEFakeS, pyEFakeS, pzEFakeS, enEFakeS );
      // set userfloats
      l.addUserFloat("px_EleUp",p4EFakeSUP.px());
      l.addUserFloat("py_EleUp",p4EFakeSUP.py());
      l.addUserFloat("pz_EleUp",p4EFakeSUP.pz());
      l.addUserFloat("e_EleUp",p4EFakeSUP.energy());
      l.addUserFloat("m_EleUp",p4EFakeSUP.mass());

      // down shift
      pxEFakeS = l.px()*udEFakeshiftP[1];
      pyEFakeS = l.py()*udEFakeshiftP[1];
      pzEFakeS = l.pz()*udEFakeshiftP[1];
      massEFakeS = l.mass()*udEFakeshiftMass[1];
      enEFakeS = TMath::Sqrt(pxEFakeS*pxEFakeS + pyEFakeS*pyEFakeS + pzEFakeS*pzEFakeS + massEFakeS*massEFakeS);
      math::XYZTLorentzVectorD p4EFakeSDOWN( pxEFakeS, pyEFakeS, pzEFakeS, enEFakeS );
      // set userfloats
      l.addUserFloat("px_EleDown",p4EFakeSDOWN.px());
      l.addUserFloat("py_EleDown",p4EFakeSDOWN.py());
      l.addUserFloat("pz_EleDown",p4EFakeSDOWN.pz());
      l.addUserFloat("e_EleDown",p4EFakeSDOWN.energy());
      l.addUserFloat("m_EleDown",p4EFakeSDOWN.mass());
    }

    // MES: muon->tauh fake ES
    // no central correction needed
    // up/down uncertainties are set at 1%
    if ( genmatch == 2 || genmatch == 4 )
    {
      float MESshift = 0.01;
      l.addUserFloat("MESshiftup",MESshift);
      l.addUserFloat("MESshiftdw",MESshift);
    }

    //--- PF ISO
    float PFChargedHadIso   = l.chargedHadronIso();
    float PFNeutralHadIso   = l.neutralHadronIso();
    float PFPhotonIso       = l.photonIso();

    float combRelIsoPF = LeptonIsoHelper::combRelIsoPF(l);

    int numChargedParticlesSignalCone = l.signalChargedHadrCands().size();
    int numNeutralHadronsSignalCone = l.signalNeutrHadrCands().size();
    int numPhotonsSignalCone = l.signalGammaCands().size();
    int numParticlesSignalCone = l.signalCands().size();
    int numChargedParticlesIsoCone = l.isolationChargedHadrCands().size();
    int numNeutralHadronsIsoCone = l.isolationNeutrHadrCands().size();
    int numPhotonsIsoCone = l.isolationGammaCands().size();
    int numParticlesIsoCone = l.isolationCands().size();
    float leadChargedParticlePt=l.leadCand()->pt();
    float trackRefPt = (l.leadChargedHadrCand().isNonnull() ? l.leadChargedHadrCand()->pt() : 0.);

    
    //Decay mode
    //int decayMode = -1;
    //int A = l.signalPFChargedHadrCands().size();
    //int B = l.signalPFGammaCands().size();
    //if(A==1&&B==0)decayMode=1;
    //else if(A==1&&B>0)decayMode=2;
    //else if (A==3&&B==0)decayMode=3;
    float tauid = (l.isTauIDAvailable(theDiscriminatorTag) ? l.tauID(theDiscriminatorTag) : -999);
    //printf("A, B, tau %d %d %f \n",A,B,tauid);

    //if(decayMode<0&&tauid==0)edm::LogWarning("TauFiller: Unrecognized decay mode");
    /*

    //--- SIP, dxy, dz
    float IP      = fabs(l.dB(pat::Electron::PV3D));
    float IPError = l.edB(pat::Electron::PV3D);
    float SIP     = IP/IPError;
    */
    
    float dxy = 999.;
    float dz  = 999.;
    if (vertexs->size()>0) {
      //dxy = l.dxy();
      //const Vertex* vertex = &(vertexs->front());          
      //dz = l.vertex().z() - vertex[0].z();

      pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(l.leadChargedHadrCand().get());
      dz=packedLeadTauCand->dz();
      dxy=packedLeadTauCand->dxy();

      //For some reasons, the reference secondaryVertex() is empty EVEN if hasSecondaryVertex is true
      //To be asked to miniAOD people
      //if(l.hasSecondaryVertex()) {
        //dz  = l.secondaryVertex().get()->z()-vertex->z();      
      //}
    } 
   
    //--- Embed user variables
    l.addUserInt("isTESShifted",isTESShifted);
    l.addUserInt("isEESShifted",isEESShifted);
    l.addUserInt("isTauMatched",isTauMatched);
    l.addUserFloat("HPSDiscriminator",tauid); 
    l.addUserFloat("decayMode",l.decayMode()); 
    l.addUserFloat("dxy",dxy); 
    l.addUserFloat("dz",dz); 
    l.addUserFloat("genmatch",genmatch);
    l.addUserFloat("PFChargedHadIso",PFChargedHadIso); 
    l.addUserFloat("PFNeutralHadIso",PFNeutralHadIso); 
    l.addUserFloat("PFPhotonIso",PFPhotonIso); 
    l.addUserFloat("combRelIsoPF",combRelIsoPF); 
    l.addUserInt("numChargedParticlesSignalCone",numChargedParticlesSignalCone);
    l.addUserInt("numNeutralHadronsSignalCone",numNeutralHadronsSignalCone);
    l.addUserInt("numPhotonsSignalCone",numPhotonsSignalCone);
    l.addUserInt("numParticlesSignalCone",numParticlesSignalCone);
    l.addUserInt("numChargedParticlesIsoCone",numChargedParticlesIsoCone);
    l.addUserInt("numNeutralHadronsIsoCone",numNeutralHadronsIsoCone);
    l.addUserInt("numPhotonsIsoCone",numPhotonsIsoCone);
    l.addUserInt("numParticlesIsoCone",numParticlesIsoCone);
    l.addUserFloat("leadChargedParticlePt",leadChargedParticlePt);
    l.addUserFloat("trackRefPt",trackRefPt); 
    
    l.addUserFloat("passCombRelIsoPFFSRCorr",true);//Fake to pass bareZCand cut

    // fill all userfloats
    for (unsigned int iuf = 0; iuf < tauFloatDiscrims_.size(); iuf++)
    {
      string ID = tauFloatDiscrims_.at(iuf);
      l.addUserFloat (ID.c_str(), l.isTauIDAvailable(ID.c_str()) ? l.tauID (ID.c_str()) : -999.);
    }

    // fill all userints
    for (unsigned int iui = 0; iui < tauIntDiscrims_.size(); iui++)
    {
      string ID = tauIntDiscrims_.at(iui);
      int ui = -999;
      if (l.isTauIDAvailable(ID.c_str()))
      {
        ui = ( (l.tauID (ID.c_str()) > 0.5) ? 1 : 0);
      }
      l.addUserInt (ID.c_str(), ui);
    }


    //--- MC parent code 
    const reco::GenParticle* genL= l.genParticleRef().get();
    float px=0,py=0,pz=0,E=0,fromH=0;
    float pxHad=0, pyHad=0, pzHad=0, EHad=0; // hadronic gen tau
    int status=99999, id=99999;

    if(genL){
      px =genL->p4().Px();
      py =genL->p4().Py();
      pz =genL->p4().Pz();
      E =genL->p4().E();
      status =genL->status();
      id =genL->pdgId();

      //cout << "Tau filler: " << i << " [px, id] = " << l.px() << " , " << l.pdgId() << " | (px, py, pz, e) " << px << " " << py << " " << pz << " " << E << " | ID: " << genL->pdgId() << " | status: " << genL->status() << endl;

      // build hadronic gen tau (all visible sons)
      for (unsigned int iDau = 0; iDau < genL->numberOfDaughters(); iDau++)
      {
          const Candidate * Dau = genL->daughter( iDau );
          int dauId = Dau->pdgId();
          if (abs(dauId) != 12 && abs(dauId) != 14 && abs(dauId) != 16)
          {
              pxHad += Dau->p4().Px();
              pyHad += Dau->p4().Py();
              pzHad += Dau->p4().Pz();
              EHad += Dau->p4().E();
          }
      }
      
      //search if it comes from H
      Handle<edm::View<reco::GenParticle> > genHandle;
      iEvent.getByToken(theGenTag, genHandle);
      for(unsigned int ipruned = 0; ipruned< genHandle->size(); ++ipruned){
        int pdgmot = (&(*genHandle)[ipruned])->pdgId();
        if(abs(pdgmot)==25){
          if(userdatahelpers::isAncestor(&(*genHandle)[ipruned],genL)){
            fromH=1;
            break;
          }
        }
      }
    }
    l.addUserFloat("genPx",px);
    l.addUserFloat("genPy",py);
    l.addUserFloat("genPz",pz);
    l.addUserFloat("genE",E);
    l.addUserInt("status", status);
    l.addUserInt("id", id);

    l.addUserFloat("fromH",fromH);

    l.addUserFloat("genHadPx",px);
    l.addUserFloat("genHadPy",py);
    l.addUserFloat("genHadPz",pz);
    l.addUserFloat("genHadE",E);

    //     MCHistoryTools mch(iEvent);
    //     if (mch.isMC()) {
    //       int MCParentCode = 0;
    //       //      int MCParentCode = mch.getParentCode(&l); //FIXME: does not work on cmg
    //       l.addUserFloat("MCParentCode",MCParentCode);
    //     }
    
    //cout<<"  unshifted (px,py,pz)  : " << l.px() << "," << l.py() << "," << l.pz() << endl;
    //cout<<"  unshifted (pt,eta,phi): " << l.pt() << "," << l.eta() << "," << l.phi() << endl;
    //cout<<"  genmatch: " << genmatch << endl;
    //cout<<"  DM      : " << l.decayMode() << endl;

    // apply the actual shift of the central value here
    l.setP4( p4S_Nominal );

    //cout<<"  shifted (px,py,pz)  : " << l.px() << "," << l.py() << "," << l.pz() << endl;
    //cout<<"  shifted (pt,eta,phi): " << l.pt() << "," << l.eta() << "," << l.phi() << endl;

    //--- Check selection cut. Being done here, flags are not available; but this way we 
    //    avoid wasting time on rejected leptons.
    if (!cut(l)) continue;

    //--- Embed flags (ie flags specified in the "flags" pset)
    for(CutSet<pat::Tau>::const_iterator flag = flags.begin(); flag != flags.end(); ++flag) {
      l.addUserFloat(flag->first,int((*(flag->second))(l)));
    }
    //cout << " ---------> DOPO pattau - p4: " << LorentzVectorE(l.p4()) << endl;
    result->push_back(l);
  }
  //iEvent.put(result);
  iEvent.put(std::move(result));
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(TauFiller);

