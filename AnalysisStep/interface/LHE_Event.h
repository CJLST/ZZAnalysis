#ifndef EVENTBASE_H
#define EVENTBASE_H

#include <vector>
#include <ZZMatrixElement/MELA/interface/MELACandidate.h>
#include "TLorentzVector.h"

class LHE_Event{
public:

  // Constructors

  LHE_Event(){}
  ~LHE_Event(){ wipeAll(); }

  // Member functions

  void constructVVCandidates(int isZZ=1, int fstype=0);
  void addVVCandidateMother(MELAParticle* mother);
  void addVVCandidateAppendages();

  int getNZZCandidates() const{ return ZZcandidates.size(); }
  int getNLeptons() const{ return leptons.size(); }
  int getNNeutrinos() const{ return neutrinos.size(); }
  int getNPhotons() const{ return photons.size(); }
  int getNJets() const{ return jets.size(); }
  int getNIntermediates() const{ return intermediates.size(); }
  int getNParticles() const{ return particles.size(); }

  MELACandidate* getZZCandidate(int index) const;
  MELAParticle* getLepton(int index) const;
  MELAParticle* getNeutrino(int index) const;
  MELAParticle* getPhoton(int index) const;
  MELAParticle* getJet(int index) const;
  MELAParticle* getIntermediate(int index) const;
  MELAParticle* getParticle(int index) const;

  void addParticle(MELAParticle* myParticle){ particles.push_back(myParticle); }
  void addIntermediate(MELAParticle* myParticle){ intermediates.push_back(myParticle); }
  void addLepton(MELAParticle* myParticle, bool genuineParticle=true){ leptons.push_back(myParticle); if (genuineParticle) addParticle(myParticle); }
  void addNeutrino(MELAParticle* myParticle, bool genuineParticle=true){ neutrinos.push_back(myParticle); if (genuineParticle) addParticle(myParticle); }
  void addPhoton(MELAParticle* myParticle, bool genuineParticle=true){ photons.push_back(myParticle); if (genuineParticle) addParticle(myParticle); }
  void addJet(MELAParticle* myParticle, bool genuineParticle=true){ jets.push_back(myParticle); if (genuineParticle) addParticle(myParticle); }

protected:
  std::vector<MELAParticle*> particles;
  std::vector<MELAParticle*> intermediates;
  std::vector<MELAParticle*> leptons;
  std::vector<MELAParticle*> neutrinos;
  std::vector<MELAParticle*> photons;
  std::vector<MELAParticle*> jets;
  std::vector<MELACandidate*> ZZcandidates;

  void addZZCandidate(MELACandidate* myParticle); // Protected to avoid adding external ZZCandidates and DELETING THEM TWICE!

  template<typename ParticleType> void wipeArray(std::vector<ParticleType*>& particleArray, bool doDelete=true){ if (doDelete){ for (unsigned int i=0; i<particleArray.size(); i++){ ParticleType* delpar = particleArray.at(i); delete delpar; } } particleArray.clear(); };
  void wipeAll(){ leptons.clear(); neutrinos.clear(); photons.clear(); jets.clear(); wipeArray(ZZcandidates, true); wipeArray(intermediates, false); wipeArray(particles, false); };
};


#endif
