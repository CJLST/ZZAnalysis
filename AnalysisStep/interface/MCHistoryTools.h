#ifndef MCHistoryTools_h
#define MCHistoryTools_h

/** \class MCHistoryTools
 *
 *  Some handy tools to examine genMatches of PAT::Object and correct some flaws in the default matching.
 *
 *  \author N. Amapane - Torino
 *  \author C. Botta   - CERN
 */

#include <FWCore/Framework/interface/Event.h>
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include <DataFormats/HepMCCandidate/interface/GenParticleFwd.h>
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>

#include <vector>
#include <string>

class MCHistoryTools {
 public:
  /// Constructor
  MCHistoryTools(const edm::Event & event, std::string sampleName="");
  MCHistoryTools(const edm::Event & event, std::string sampleName, edm::Handle<edm::View<reco::Candidate> > & genParticles, edm::Handle<GenEventInfoProduct> & genInfo);

  /// Destructor
  virtual ~MCHistoryTools();
  
  bool isMC() {return ismc;}
  
  int genFinalState();

  // flavour of the associated V decay (0 if n/a)
  int genAssociatedFS();

  /// ZZInAcceptance = 60, 120 on MZ21, MZ2 (mass from the MC Z)
  /// ZZ4lInEtaAcceptance = eta cut on the 4 leptons (from the direct Z daughter)
  /// ZZ4lInEtaPtAcceptance = eta, pt cuts on the 4 leptons (from the direct Z daughter)
  void genAcceptance(bool& gen_ZZ4lInEtaAcceptance, bool& gen_ZZ4lInEtaPtAcceptance);

  const reco::Candidate * genH() {init(); return theGenH;}

  const std::vector<const reco::Candidate *>& genZs() {init(); return theGenZ;}

  /// Z1 is defined as the Z closest to nominal mass.
//   const reco::Candidate * genZ1() {init(); return theGenZ1;}
//   const reco::Candidate * genZ2() {init(); return theGenZ2;}
  
  /// The leptons coming from Zs (in no specific order) 
  const std::vector<const reco::Candidate *>& genZLeps() {init(); return theGenLeps;}
    
  /// The leptons coming from ZZ or HZZ, sorted according to the reco-level criteria. Equals to genZLeps() if there are !=4 leptons from ZZ.
  const std::vector<const reco::Candidate *>& sortedGenZZLeps() {init(); return theSortedGenLepts;}
  
  /// The leptons coming from associated W, Z, or ttbar
  const std::vector<const reco::Candidate *>& genAssociatedLeps() {init(); return theAssociatedLeps;}

  /// The FSR photons
  const std::vector<const reco::Candidate *>& genFSR() {init(); return theGenFSR;}

  /// Find the actual lepton parent (first parent in MC history with a different pdgID)
  const reco::GenParticle* getParent(const reco::GenParticle* genLep);

  /// Same as the above, but try recovery in case default PAT matching fails.
  /// (This happens rather often for electron due to brems/FSR, since the default matching handles this poorly).
  /// The recovery consists in matching between a selected list of genleptons, cf. getMatch().
  const reco::GenParticle* getParent(const pat::Electron* lep, const std::vector<const reco::Candidate *>& gen4lep);

  /// Manual matching with closest same-flavour gen lepton (of any status). 
  /// This was tested to work great when the provided candidates are e.g. only the signal ones.
  const reco::GenParticle* getMatch(const pat::Electron* lep, const std::vector<const reco::Candidate *>& gen4lep);  

  /// Return the code of the particle's parent: 25 for H->Z->l; 23 for Z->l; +-15 for tau->l if genlep is e,mu.
  int getParentCode(const reco::GenParticle* genLep);

  /// Same as above, but if no match is found, search for a match within gen4lep
  /// Cf. getParent(lep, gen4lep) for details.
  int getParentCode(const pat::Electron* lep, const std::vector<const reco::Candidate *>& gen4lep);

  /// Shortcut for the above, using genZLeps() for the recovery matching.
  int getParentCode(const pat::Electron* lep) {return getParentCode(lep, genZLeps());}

  unsigned int getProcessID() {return processID;}
  float gethepMCweight() {return hepMCweight;}
    

  static int fsrMatch(const reco::Candidate* recoFSR, const std::vector<const reco::Candidate*>& genFSRs);

 private:
  edm::Handle<edm::View<reco::Candidate> > particles;
  bool ismc;
  unsigned int processID;
  float hepMCweight;
  
  bool isInit;
  const reco::Candidate * theGenH; 
  std::vector<const reco::Candidate *> theGenZ;
  std::vector<const reco::Candidate *> theAssociatedV;
  std::vector<const reco::Candidate *> theGenLeps;
  std::vector<const reco::Candidate *> theSortedGenLepts;
  std::vector<const reco::Candidate *> theAssociatedLeps;
  std::vector<const reco::Candidate *> theGenFSR;
  
  void init();

};
#endif

