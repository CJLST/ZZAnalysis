#include <iostream>
#include <ZZAnalysis/AnalysisStep/interface/LHE_Event.h>

using namespace PDGHelpers;

void LHE_Event::addZZCandidate(MELACandidate* myParticle){
  bool isIdentical = (getNZZCandidates()>0);
  for (int cc=0; cc<getNZZCandidates(); cc++){
    MELACandidate* testCand = ZZcandidates.at(cc);

    for (int i=0; i<2; i++){
      MELAParticle* testV = testCand->getSortedV(i);
      MELAParticle* partV = myParticle->getSortedV(i);
      if (partV==0 || testV==0){
        if (partV==0 && testV==0) continue; // Check the next intermediate V if there is any.
        else { isIdentical=false; break; } // They are definitely not the same.
      }
      if (testV->getNDaughters() != partV->getNDaughters()) { isIdentical=false;  break; }; // Again, they cannot be the same.
      for (int j=0; j<testV->getNDaughters(); j++){
        MELAParticle* testD = testV->getDaughter(j);
        MELAParticle* partD = partV->getDaughter(j);
        isIdentical = isIdentical && (testD==partD);
        if (!isIdentical) break;
      }
      if (!isIdentical) break;
    }
    if (!isIdentical) break;
  }
  if (!isIdentical) ZZcandidates.push_back(myParticle);
  else { delete myParticle; myParticle=0; }
}

void LHE_Event::constructVVCandidates(int isZZ, int fstype){
  /*
  fstype  / ZZ==1 / WW==0  / Yukawa==2 / Zgam=3 / gamgam=4 / Z+nj=5
  fstype=0: 4l    / lnulnu / 2l        / 2l     / gam      / 2l
  fstype=1: 4q    / 4q     / 2q        / 2q     / -        / 2q
  fstype=2: 2l2q  / lnu2q  / -         / -      / -        / -
  fstype=3: 2l2nu / -      / -         / -      / -        / -
  fstype=4: 2q2nu / -      / -         / -      / -        / -
  fstype=5: 4nu   / -      / -         / 2nu    / -        / 2nu
  fstype=-1: Any
  fstype=-2: 2l2X
  fstype=-3: 2nu2X
  fstype=-4: 2q2X
  */

  if (
    (isZZ<=0 && fstype>2)
    ||
    (isZZ==1 && fstype>5)
    ||
    (isZZ==2 && fstype>1)
    ||
    (isZZ==3 && (fstype>1 && fstype!=5))
    ||
    (isZZ==4 && fstype>0)
    ||
    (isZZ==5 && (fstype>1 && fstype!=5))
    ||
    isZZ>5
    ||
    (fstype<-1 && isZZ>1)
    ||
    fstype<-4
    ){
    if (isZZ<0) std::cerr << "No " << "undecayed" << " candidate with final state " << fstype << " is possible!" << std::endl;
    else if (isZZ==0) std::cerr << "No " << "WW" << " candidate with final state " << fstype << " is possible!" << std::endl;
    else if (isZZ==1) std::cerr << "No " << "ZZ" << " candidate with final state " << fstype << " is possible!" << std::endl;
    else if (isZZ==2) std::cerr << "No " << "f-fbar" << " candidate with final state " << fstype << " is possible!" << std::endl;
    else if (isZZ==3) std::cerr << "No " << "Zgamma" << " candidate with final state " << fstype << " is possible!" << std::endl;
    else if (isZZ==4) std::cerr << "No " << "gammagamma" << " candidate with final state " << fstype << " is possible!" << std::endl;
    else if (isZZ==5) std::cerr << "No " << "Z+(n)jets" << " candidate with final state " << fstype << " is possible!" << std::endl;
    else if (isZZ>5) std::cerr << "Unknown candidate with final state " << fstype << "!" << std::endl;
    else if (fstype<-4) std::cerr << "Unknown final state " << fstype << "!" << std::endl;
    return;
  }

  std::vector<MELAParticle*> lepMinusPlus[3][2]; // l-, l+
  std::vector<MELAParticle*> lepNuNubar[3][2]; // nu, nub
  std::vector<MELAParticle*> quarkAntiquark[7][2]; // q, qb

  for (std::vector<MELAParticle*>::iterator it = leptons.begin(); it<leptons.end(); it++){ // Leptons
    int iFirst=0;
    int iSecond=0;

    if (abs((*it)->id)==11) iFirst = 0;
    else if (abs((*it)->id)==13) iFirst = 1;
    else if (abs((*it)->id)==15) iFirst = 2;
    if ((*it)->id<0) iSecond=1;
    lepMinusPlus[iFirst][iSecond].push_back(*it);
  }
  for (std::vector<MELAParticle*>::iterator it = neutrinos.begin(); it<neutrinos.end(); it++){ // Neutrinos
    int iFirst=0;
    int iSecond=0;

    if (abs((*it)->id)==12) iFirst = 0;
    else if (abs((*it)->id)==14) iFirst = 1;
    else if (abs((*it)->id)==16) iFirst = 2;
    if ((*it)->id<0) iSecond=1;
    lepNuNubar[iFirst][iSecond].push_back(*it);
  }
  for (std::vector<MELAParticle*>::iterator it = jets.begin(); it<jets.end(); it++){ // Jets
    int iFirst=abs((*it)->id); // Yes, 0-6, 0 being unknown
    if (PDGHelpers::isAGluon(iFirst)) continue;
    int iSecond=0;

    if ((*it)->id<0) iSecond=1;
    quarkAntiquark[iFirst][iSecond].push_back(*it);
  }

  std::vector<MELAParticle*> tmpVhandle;

  if (isZZ==1 || isZZ==3 || isZZ==5){ // ZZ

    if (fstype<0 || (isZZ==1 && (fstype==0 || fstype==2 || fstype==3)) || (isZZ==3 && fstype==0)){ // Z->2l
      for (int c=0; c<3; c++){
        for (unsigned int i=0; i<lepMinusPlus[c][0].size(); i++){
          for (unsigned int j=0; j<lepMinusPlus[c][1].size(); j++){
            TLorentzVector pV = lepMinusPlus[c][0].at(i)->p4+lepMinusPlus[c][1].at(j)->p4;
            MELAParticle* V = new MELAParticle(23, pV);
            V->addDaughter(lepMinusPlus[c][0].at(i));
            V->addDaughter(lepMinusPlus[c][1].at(j));
            tmpVhandle.push_back(V);
          }
        }
      }
    }
    if (fstype<0 || (isZZ==1 && (fstype==3 || fstype==4 || fstype==5)) || (isZZ==3 && fstype==5)){ // Z->2nu
      for (int c=0; c<3; c++){
        for (unsigned int i=0; i<lepNuNubar[c][0].size(); i++){
          for (unsigned int j=0; j<lepNuNubar[c][1].size(); j++){
            TLorentzVector pV = lepNuNubar[c][0].at(i)->p4+lepNuNubar[c][1].at(j)->p4;
            MELAParticle* V = new MELAParticle(23, pV);
            V->addDaughter(lepNuNubar[c][0].at(i));
            V->addDaughter(lepNuNubar[c][1].at(j));
            tmpVhandle.push_back(V);
          }
        }
      }
    }
    if (fstype<0 || (isZZ==1 && (fstype==1 || fstype==2 || fstype==4)) || (isZZ==3 && fstype==1)){ // Z->2q
      for (int c=1; c<7; c++){
        for (unsigned int i=0; i<quarkAntiquark[c][0].size(); i++){
          for (unsigned int j=0; j<quarkAntiquark[c][1].size(); j++){
            TLorentzVector pV = quarkAntiquark[c][0].at(i)->p4+quarkAntiquark[c][1].at(j)->p4;
            MELAParticle* V = new MELAParticle(23, pV);
            V->addDaughter(quarkAntiquark[c][0].at(i));
            V->addDaughter(quarkAntiquark[c][1].at(j));
            tmpVhandle.push_back(V);
          }
        }
      }
    }

  }
  else if(isZZ==0){ // WW

    if (fstype<0 || fstype==0 || fstype==2){ // W->lnu
      for (int c=0; c<3; c++){
        for (unsigned int i=0; i<lepMinusPlus[c][1].size(); i++){
          for (unsigned int j=0; j<lepNuNubar[c][0].size(); j++){
            TLorentzVector pV = lepMinusPlus[c][1].at(i)->p4+lepNuNubar[c][0].at(j)->p4;
            MELAParticle* V = new MELAParticle(24, pV);
            V->addDaughter(lepMinusPlus[c][1].at(i));
            V->addDaughter(lepNuNubar[c][0].at(j));
            tmpVhandle.push_back(V);
          }
        }
        for (unsigned int i=0; i<lepMinusPlus[c][0].size(); i++){
          for (unsigned int j=0; j<lepNuNubar[c][1].size(); j++){
            TLorentzVector pV = lepMinusPlus[c][0].at(i)->p4+lepNuNubar[c][1].at(j)->p4;
            MELAParticle* V = new MELAParticle(-24, pV);
            V->addDaughter(lepMinusPlus[c][0].at(i));
            V->addDaughter(lepNuNubar[c][1].at(j));
            tmpVhandle.push_back(V);
          }
        }
      }
    }
    if (fstype<0 || fstype==1 || fstype==2){ // W->2q
      for (int c=1; c<7; c++){
        for (int d=1; d<7; d++){
          if (d==c) continue;
          for (unsigned int i=0; i<quarkAntiquark[c][0].size(); i++){
            for (unsigned int j=0; j<quarkAntiquark[d][1].size(); j++){
              int totalcharge = quarkAntiquark[c][0].at(i)->charge() + quarkAntiquark[d][1].at(j)->charge();
              if (abs(totalcharge)!=1) continue;

              TLorentzVector pV = quarkAntiquark[c][0].at(i)->p4+quarkAntiquark[d][1].at(j)->p4;
              MELAParticle* V = new MELAParticle(24*totalcharge, pV);
              V->addDaughter(quarkAntiquark[c][0].at(i));
              V->addDaughter(quarkAntiquark[d][1].at(j));
              tmpVhandle.push_back(V);
            }
          }
        }
      }
    }
  }
  else if (isZZ==2){ // H->f fbar

    if (fstype<0 || fstype==0){ // H->2l
      for (int c=0; c<3; c++){
        for (unsigned int i=0; i<lepMinusPlus[c][0].size(); i++){
          for (unsigned int j=0; j<lepMinusPlus[c][1].size(); j++){
            MELAParticle* F1 = lepMinusPlus[c][0].at(i);
            MELAParticle* F2 = lepMinusPlus[c][1].at(j);

            TLorentzVector pH = F1->p4+F2->p4;
            MELACandidate* cand = new MELACandidate(25, pH, true);
            cand->addDaughter(F1);
            cand->addDaughter(F2);

            double defaultHVVmass = HVVmass;
            setHVVmass(Zeromass);
            cand->sortDaughters();
            setHVVmass(defaultHVVmass);
            addZZCandidate(cand);
          }
        }
      }
    }
    if (fstype<0 || fstype==1){ // H->2q
      for (int c=1; c<7; c++){
        for (unsigned int i=0; i<quarkAntiquark[c][0].size(); i++){
          for (unsigned int j=0; j<quarkAntiquark[c][1].size(); j++){
            MELAParticle* F1 = quarkAntiquark[c][0].at(i);
            MELAParticle* F2 = quarkAntiquark[c][1].at(j);

            TLorentzVector pH = F1->p4+F2->p4;
            MELACandidate* cand = new MELACandidate(25, pH, true);
            cand->addDaughter(F1);
            cand->addDaughter(F2);

            double defaultHVVmass = HVVmass;
            setHVVmass(Zeromass);
            cand->sortDaughters();
            setHVVmass(defaultHVVmass);
            addZZCandidate(cand);
          }
        }
      }
    }

  }
  else if (isZZ==5){ // Z->f fbar

    if (fstype<0 || fstype==0){ // Z->2l
      for (int c=0; c<3; c++){
        for (unsigned int i=0; i<lepMinusPlus[c][0].size(); i++){
          for (unsigned int j=0; j<lepMinusPlus[c][1].size(); j++){
            MELAParticle* F1 = lepMinusPlus[c][0].at(i);
            MELAParticle* F2 = lepMinusPlus[c][1].at(j);

            TLorentzVector pCand = F1->p4+F2->p4;
            MELACandidate* cand = new MELACandidate(23, pCand, true);
            cand->addDaughter(F1);
            cand->addDaughter(F2);

            double defaultHVVmass = HVVmass;
            setHVVmass(Zeromass);
            cand->sortDaughters();
            setHVVmass(defaultHVVmass);
            addZZCandidate(cand);
          }
        }
      }
    }
    if (fstype<0 || fstype==1){ // Z->2q
      for (int c=1; c<7; c++){
        for (unsigned int i=0; i<quarkAntiquark[c][0].size(); i++){
          for (unsigned int j=0; j<quarkAntiquark[c][1].size(); j++){
            MELAParticle* F1 = quarkAntiquark[c][0].at(i);
            MELAParticle* F2 = quarkAntiquark[c][1].at(j);

            TLorentzVector pCand = F1->p4+F2->p4;
            MELACandidate* cand = new MELACandidate(23, pCand, true);
            cand->addDaughter(F1);
            cand->addDaughter(F2);

            double defaultHVVmass = HVVmass;
            setHVVmass(Zeromass);
            cand->sortDaughters();
            setHVVmass(defaultHVVmass);
            addZZCandidate(cand);
          }
        }
      }
    }
    if (fstype<0 || fstype==5){ // Z->2nu
      for (int c=0; c<3; c++){
        for (unsigned int i=0; i<lepNuNubar[c][0].size(); i++){
          for (unsigned int j=0; j<lepNuNubar[c][1].size(); j++){
            MELAParticle* F1 = lepNuNubar[c][0].at(i);
            MELAParticle* F2 = lepNuNubar[c][1].at(j);

            TLorentzVector pCand = F1->p4+F2->p4;
            MELACandidate* cand = new MELACandidate(23, pCand, true);
            cand->addDaughter(F1);
            cand->addDaughter(F2);

            double defaultHVVmass = HVVmass;
            setHVVmass(Zeromass);
            cand->sortDaughters();
            setHVVmass(defaultHVVmass);
            addZZCandidate(cand);
          }
        }
      }
    }

  }
  else{ // Undecayed
    for (std::vector<MELAParticle*>::iterator it = intermediates.begin(); it<intermediates.end(); it++){ // Add directly
      if (isAHiggs((*it)->id)){
        TLorentzVector pH = (*it)->p4;
        MELACandidate* cand = new MELACandidate(25, pH, true);

        double defaultHVVmass = HVVmass;
        setHVVmass(Zeromass);
        cand->sortDaughters();
        setHVVmass(defaultHVVmass);
        addZZCandidate(cand);
      }
    }
  }

  if (debugVars::debugFlag) std::cout << "Number of V/ZZ before sorting photons: " << tmpVhandle.size() << " " << getNZZCandidates() << std::endl;

  if (isZZ==3 || isZZ==4){
    for (unsigned int i=0; i<photons.size(); i++){ // Copy the photons
      TLorentzVector pV = photons.at(i)->p4;
      MELAParticle* V = new MELAParticle(photons.at(i)->id, pV);
      V->addDaughter(photons.at(i)); // Photon is its own daughter!
      tmpVhandle.push_back(V);
    }
  }

  if (debugVars::debugFlag) std::cout << "Number of V/ZZ after sorting photons: " << tmpVhandle.size() << " " << getNZZCandidates() << std::endl;

  if (
    ((fstype<0 || fstype==1 || fstype==2 || fstype==4) && (isZZ==0 || isZZ==1)) // W/Z->2j reco.-level
    ||
    ((fstype<0 || fstype==1) && isZZ==2) // H->2j reco.-level
    ||
    ((fstype<0 || fstype==1) && isZZ==3) // H->Zgam with Z->2j
    ||
    ((fstype<0 || fstype==1) && isZZ==5) // Z+2jets with Z->2j
    ){
    for (unsigned int i=0; i<quarkAntiquark[0][0].size(); i++){
      if (quarkAntiquark[0][0].at(i)->id!=0) continue;
      for (unsigned int j=i+1; j<quarkAntiquark[0][0].size(); j++){
        if (quarkAntiquark[0][0].at(j)->id!=0) continue;
        if (isZZ==0 || isZZ==1 || isZZ==3){
          TLorentzVector pV = quarkAntiquark[0][0].at(i)->p4+quarkAntiquark[0][0].at(j)->p4;
          MELAParticle* V = new MELAParticle(0, pV);
          V->addDaughter(quarkAntiquark[0][0].at(i));
          V->addDaughter(quarkAntiquark[0][0].at(j));
          tmpVhandle.push_back(V);
        }
        else if (isZZ==2){
          MELAParticle* F1 = quarkAntiquark[0][0].at(i);
          MELAParticle* F2 = quarkAntiquark[0][0].at(j);

          TLorentzVector pH = F1->p4+F2->p4;
          MELACandidate* cand = new MELACandidate(25, pH, true);
          cand->addDaughter(F1);
          cand->addDaughter(F2);

          double defaultHVVmass = HVVmass;
          setHVVmass(Zeromass);
          cand->sortDaughters();
          setHVVmass(defaultHVVmass);
          addZZCandidate(cand);
        }
        else if (isZZ==5){
          MELAParticle* F1 = quarkAntiquark[0][0].at(i);
          MELAParticle* F2 = quarkAntiquark[0][0].at(j);

          TLorentzVector pCand = F1->p4+F2->p4;
          MELACandidate* cand = new MELACandidate(23, pCand, true);
          cand->addDaughter(F1);
          cand->addDaughter(F2);

          double defaultHVVmass = HVVmass;
          setHVVmass(Zeromass);
          cand->sortDaughters();
          setHVVmass(defaultHVVmass);
          addZZCandidate(cand);
        }
      }
    }
    if (debugVars::debugFlag) std::cout << "Number of V/ZZ after sorting reco. jets: " << tmpVhandle.size() << " " << getNZZCandidates() << std::endl;
  }



  for (unsigned int i=0; i<tmpVhandle.size(); i++){
    for (unsigned int j=i; j<tmpVhandle.size(); j++){
      if (tmpVhandle.at(i)==tmpVhandle.at(j)) continue;
      if ((tmpVhandle.at(i)->charge()+tmpVhandle.at(j)->charge())!=0) continue;

      MELAParticle* Vi1 = tmpVhandle.at(i)->getDaughter(0);
      MELAParticle* Vi2 = tmpVhandle.at(i)->getDaughter(1);
      MELAParticle* Vj1 = tmpVhandle.at(j)->getDaughter(0);
      MELAParticle* Vj2 = tmpVhandle.at(j)->getDaughter(1);
      /*
      std::cout << "11: " << Vi1->id << '\t' << Vi1->x() << '\t' << Vi1->y() << '\t' << Vi1->z() << '\t' << Vi1->t() << '\t' << std::endl;
      std::cout << "12: " << Vi2->id << '\t' << Vi2->x() << '\t' << Vi2->y() << '\t' << Vi2->z() << '\t' << Vi2->t() << '\t' << std::endl;
      std::cout << "21: " << Vj1->id << '\t' << Vj1->x() << '\t' << Vj1->y() << '\t' << Vj1->z() << '\t' << Vj1->t() << '\t' << std::endl;
      std::cout << "22: " << Vj2->id << '\t' << Vj2->x() << '\t' << Vj2->y() << '\t' << Vj2->z() << '\t' << Vj2->t() << '\t' << std::endl;
      */
      if (Vi1==Vj1 || (Vi2==Vj2 && Vi2 != 0)) continue;
      bool createCandidate=true;
      if (isZZ<=1 && fstype<-1){
        unsigned int partcounter=0;
        if (fstype==-2){ // Count leptons
          if (Vi1!=0 && PDGHelpers::isALepton(Vi1->id)) partcounter++;
          if (Vi2!=0 && PDGHelpers::isALepton(Vi2->id)) partcounter++;
          if (Vj1!=0 && PDGHelpers::isALepton(Vj1->id)) partcounter++;
          if (Vj2!=0 && PDGHelpers::isALepton(Vj2->id)) partcounter++;
        }
        else if (fstype==-3){ // Count neutrinos
          if (Vi1!=0 && PDGHelpers::isANeutrino(Vi1->id)) partcounter++;
          if (Vi2!=0 && PDGHelpers::isANeutrino(Vi2->id)) partcounter++;
          if (Vj1!=0 && PDGHelpers::isANeutrino(Vj1->id)) partcounter++;
          if (Vj2!=0 && PDGHelpers::isANeutrino(Vj2->id)) partcounter++;
        }
        else if (fstype==-4){ // Count jets
          if (Vi1!=0 && PDGHelpers::isAJet(Vi1->id)) partcounter++;
          if (Vi2!=0 && PDGHelpers::isAJet(Vi2->id)) partcounter++;
          if (Vj1!=0 && PDGHelpers::isAJet(Vj1->id)) partcounter++;
          if (Vj2!=0 && PDGHelpers::isAJet(Vj2->id)) partcounter++;
        }
        if (partcounter<2) createCandidate=false;
      }
      if (!createCandidate) continue;

      if (debugVars::debugFlag){
        if (Vi1!=0) std::cout << "Vi1 not zero. Id: " << Vi1->id << std::endl;
        if (Vi2!=0) std::cout << "Vi2 not zero. Id: " << Vi2->id << std::endl;
        if (Vj1!=0) std::cout << "Vj1 not zero. Id: " << Vj1->id << std::endl;
        if (Vj2!=0) std::cout << "Vj2 not zero. Id: " << Vj2->id << std::endl;
      }

      TLorentzVector pH(0, 0, 0, 0);
      if (Vi1!=0) pH = pH + Vi1->p4;
      if (Vi2!=0) pH = pH + Vi2->p4;
      if (Vj1!=0) pH = pH + Vj1->p4;
      if (Vj2!=0) pH = pH + Vj2->p4;
      MELACandidate* cand = new MELACandidate(25, pH, true);

      if (Vi1!=0) cand->addDaughter(Vi1);
      if (Vi2!=0) cand->addDaughter(Vi2);
      if (Vj1!=0) cand->addDaughter(Vj1);
      if (Vj2!=0) cand->addDaughter(Vj2);

      double defaultHVVmass = HVVmass;
      if (isZZ==0) setHVVmass(Wmass);
      else if (isZZ!=4) setHVVmass(Zmass);
      else setHVVmass(Zeromass);

      if (debugVars::debugFlag) std::cout << "Sorting daughters..." << std::endl;
      cand->sortDaughters();
      if (debugVars::debugFlag) std::cout << "Sorted daughters successfully!" << std::endl;
      setHVVmass(defaultHVVmass);

      addZZCandidate(cand);

      if (debugVars::debugFlag) std::cout << "Added candidate for V" << i << " V" << j << std::endl;
    }
  }

  for (unsigned int i=0; i<tmpVhandle.size(); i++) delete tmpVhandle.at(i);
  if (debugVars::debugFlag) std::cout << "tmpVhandle deletion step is done." << std::endl;
  tmpVhandle.clear();
}

MELACandidate* LHE_Event::getZZCandidate(int index)const{
  if ((int)ZZcandidates.size()>index) return ZZcandidates.at(index);
  else return 0;
}
MELAParticle* LHE_Event::getLepton(int index)const{
  if ((int)leptons.size()>index) return leptons.at(index);
  else return 0;
}
MELAParticle* LHE_Event::getNeutrino(int index)const{
  if ((int)neutrinos.size()>index) return neutrinos.at(index);
  else return 0;
}
MELAParticle* LHE_Event::getPhoton(int index)const{
  if ((int)photons.size()>index) return photons.at(index);
  else return 0;
}
MELAParticle* LHE_Event::getJet(int index)const{
  if ((int)jets.size()>index) return jets.at(index);
  else return 0;
}
MELAParticle* LHE_Event::getIntermediate(int index)const{
  if ((int)intermediates.size()>index) return intermediates.at(index);
  else return 0;
}
MELAParticle* LHE_Event::getParticle(int index)const{
  if ((int)particles.size()>index) return particles.at(index);
  else return 0;
}

void LHE_Event::addVVCandidateMother(MELAParticle* mother){
  for (std::vector<MELACandidate*>::iterator it = ZZcandidates.begin(); it<ZZcandidates.end(); it++) (*it)->addMother(mother);
}
void LHE_Event::addVVCandidateAppendages(){
  for (std::vector<MELACandidate*>::iterator it = ZZcandidates.begin(); it<ZZcandidates.end(); it++){
    for (std::vector<MELAParticle*>::iterator iL = leptons.begin(); iL<leptons.end(); iL++){ if ((*iL)->passSelection) (*it)->addAssociatedLeptons(*iL); }
    for (std::vector<MELAParticle*>::iterator iL = neutrinos.begin(); iL<neutrinos.end(); iL++){ if ((*iL)->passSelection) (*it)->addAssociatedNeutrinos(*iL); }
    for (std::vector<MELAParticle*>::iterator iP = photons.begin(); iP<photons.end(); iP++){ if ((*iP)->passSelection) (*it)->addAssociatedPhotons(*iP); }
    for (std::vector<MELAParticle*>::iterator iJ = jets.begin(); iJ<jets.end(); iJ++){ if ((*iJ)->passSelection) (*it)->addAssociatedJets(*iJ); }
    (*it)->addAssociatedVs();
  }
}





