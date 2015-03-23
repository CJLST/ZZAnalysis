#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>

#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>


using namespace std;
using namespace reco;

void userdatahelpers::embedDaughterData(pat::CompositeCandidate& cand) {

  for (unsigned i = 0; i<cand.numberOfDaughters(); ++i) {
    const reco::Candidate* d = cand.daughter(i)->masterClone().get();

    // We need the concrete object to access the method userFloat(). 
    // (A more general solution would be to creat a StringObjectFunction on the fly for each 
    // entry in userFloatNames(). That's maybe too time consuming (to be checked))
    if (const pat::CompositeCandidate* cc = dynamic_cast<const pat::CompositeCandidate*>(d)) {
      embedDaughterData(cand, i, cc);      
    } else if (const pat::Muon* mu = dynamic_cast<const pat::Muon*>(d)) {
      embedDaughterData(cand, i, mu);
    } else if (const pat::Electron* ele = dynamic_cast<const pat::Electron*>(d)) {
      embedDaughterData(cand, i, ele);
    } else {
      cout << "DaughterDataEmbedder: Unsupported daughter type" << endl;
    }
  }
}


float userdatahelpers::getUserFloat(const reco::Candidate* c, const char* name){
  //if(c->hasMasterClone()) c = c->masterClone().get();
  if (const pat::Muon* mu = dynamic_cast<const pat::Muon*>(c)) {
    return mu->userFloat(name);
  } else if (const pat::Electron* ele = dynamic_cast<const pat::Electron*>(c)) {
    return ele->userFloat(name);
  } else if (const pat::CompositeCandidate* cc = dynamic_cast<const pat::CompositeCandidate*>(c)) {
    return cc->userFloat(name);
  }
  return 0;
}


const PhotonPtrVector*  
userdatahelpers::getUserPhotons(const reco::Candidate* c){
  if(c->hasMasterClone())  c = c->masterClone().get();
  if (abs(c->pdgId())==13) {
    const pat::Muon* mu = static_cast<const pat::Muon*>(c);
    if (mu->hasUserData("FSRCandidates")){
      return mu->userData<PhotonPtrVector>("FSRCandidates");
    } else return 0;
  } else if (abs(c->pdgId())==11) {
    const pat::Electron* ele = static_cast<const pat::Electron*>(c);
    if (ele->hasUserData("FSRCandidates")){
      return ele->userData<PhotonPtrVector>("FSRCandidates");
    } else return 0;
  } else {
    cout << "ERROR: userdatahelpers::getUserPhotons" << endl;
    abort();
  }
}


void 
userdatahelpers::getSortedLeptons(const pat::CompositeCandidate& cand, vector<const Candidate*>& leptons, vector<string>& labels, bool is4l) {

  if (is4l) { // Regular 4 lepton SR/CR
    // Pointers to Z and leptons
    const Candidate* Z1   = cand.daughter("Z1");
    const Candidate* Z2   = cand.daughter("Z2");
    const Candidate* Z1Lp = Z1->daughter(0);
    const Candidate* Z1Ln = Z1->daughter(1);
    const Candidate* Z2Lp = Z2->daughter(0);
    const Candidate* Z2Ln = Z2->daughter(1);
    // corresponding prefixes for UserFloats
    string Z1Label = "d0.";
    string Z2Label = "d1.";
    if (cand.daughter("Z1")==cand.daughter(1)) swap(Z1Label,Z2Label);
    string Z1LpLabel = Z1Label + "d0.";
    string Z1LnLabel = Z1Label + "d1.";
    string Z2LpLabel = Z2Label + "d0.";
    string Z2LnLabel = Z2Label + "d1.";

    // Set charge according to the label; do nothing for the same-sign collections used for CRs
    if (Z1Lp->charge() < 0 && Z1Lp->charge()*Z1Ln->charge()<0) {
      swap(Z1Lp,Z1Ln);
      swap(Z1LpLabel,Z1LnLabel);
    }
    if (Z2Lp->charge() < 0 && Z2Lp->charge()*Z2Ln->charge()<0) {
      swap(Z2Lp,Z2Ln);
      swap(Z2LpLabel,Z2LnLabel);
    }
      
    // Put the four daughter leptons in a vector, ordered in a standard way
    leptons.resize(4);
    leptons[0]=Z1Lp;
    leptons[1]=Z1Ln;
    leptons[2]=Z2Lp;
    leptons[3]=Z2Ln;
    labels.resize(4);
    labels[0] = Z1LpLabel;
    labels[1] = Z1LnLabel;
    labels[2] = Z2LpLabel;
    labels[3] = Z2LnLabel;

  } else { // Z+l
    const Candidate* Z1 = cand.daughter(0); // the Z
    const reco::Candidate* Z1Lp = Z1->daughter(0);
    const reco::Candidate* Z1Ln = Z1->daughter(1);
    string Z1LpLabel = "d0.d0.";
    string Z1LnLabel = "d0.d1.";
    if (Z1Lp->charge() < 0 && Z1Lp->charge()*Z1Ln->charge()<0) {
      swap(Z1Lp,Z1Ln);
      swap(Z1LpLabel,Z1LnLabel);
    }
    
    leptons.resize(3);
    leptons[0]=Z1Lp;
    leptons[1]=Z1Ln;
    leptons[2] = cand.daughter(1);
    labels.resize(3);
    labels[0]= Z1LpLabel;
    labels[1]=Z1LnLabel;
    labels[2]="d1.";
  }
}
