/** \file
 *
 *
 *  \author A. Tarabini - LLR
 */

#include <ZZAnalysis/AnalysisStep/interface/GenTools.h>
#include <DataFormats/HepMCCandidate/interface/GenParticleFwd.h>
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>
#include <DataFormats/PatCandidates/interface/PackedGenParticle.h>

#include "ZZAnalysis/AnalysisStep/interface/HZZ4LGENAna.h"

#include <DataFormats/PatCandidates/interface/Jet.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <JetMETCorrections/Objects/interface/JetCorrectionsRecord.h>
#include <JetMETCorrections/Modules/interface/JetResolution.h>

#include <MelaAnalytics/CandidateLOCaster/interface/MELACandidateRecaster.h>

#include <tuple>


GenTools::GenTools(edm::Handle<reco::GenParticleCollection> _pruned, edm::Handle<edm::View<pat::PackedGenParticle> > _packed, edm::Handle<edm::View<reco::GenJet> > _genJ){
  pruned = _pruned;
  packed = _packed;
  genJ = _genJ;
}

GenTools::~GenTools(){}


void GenTools::init(){
  reco::GenParticleCollection::const_iterator genPart;
  int j = -1;
  int nGENLeptons=0;

  theLepts.clear();
  theLeptsId.clear();
  Lepts.clear();
  LeptsStatus.clear();
  LeptsId.clear();
  LeptsMom.clear();
  LeptsMomMom.clear();
  Lepts_RelIso.clear();
  theJets_pt30_eta4p7.clear();
  theJets_pt30_eta2p5.clear();
  Lep_Hindex.clear();
  theHiggs.clear();
  theZs.clear();
  theZsMom.clear();
  theZsDaughters.clear();
  theJets_pt30_eta4p7.clear();
  theJets_pt30_eta2p5.clear();
  for (int i=0; i<4; ++i) {Lep_Hindex_tmp[i]=-1;};//position of Higgs candidate leptons in lep_p4: 0 = Z1 lead, 1 = Z1 sub, 2 = Z2 lead, 3 = Z3 sub

  for(genPart = pruned->begin(); genPart != pruned->end(); genPart++) {
      j++;

      if (abs(genPart->pdgId())==11  || abs(genPart->pdgId())==13 || abs(genPart->pdgId())==15) {

          if (!(genPart->status()==1 || abs(genPart->pdgId())==15)) continue;
          if (!(genAna.MotherID(&pruned->at(j))==23 || genAna.MotherID(&pruned->at(j))==443 || genAna.MotherID(&pruned->at(j))==553 || abs(genAna.MotherID(&pruned->at(j)))==24) ) continue;

          nGENLeptons++;
          // Collect FSR photons
          TLorentzVector lep_dressed;
          lep_dressed.SetPtEtaPhiE(genPart->pt(),genPart->eta(),genPart->phi(),genPart->energy());
          set<int> gen_fsrset;
          for(size_t k=0; k<packed->size();k++){
              if( (*packed)[k].status() != 1) continue; // stable particles only
              if( (*packed)[k].pdgId() != 22) continue; // only photons
              double this_dR_lgamma = deltaR(genPart->eta(), genPart->phi(), (*packed)[k].eta(), (*packed)[k].phi());
              bool idmatch=false;
              if ((*packed)[k].mother(0)->pdgId()==genPart->pdgId() ) idmatch=true;
              const reco::Candidate * mother = (*packed)[k].mother(0);
              for(size_t m=0;m<mother->numberOfMothers();m++) {
                  if ( (*packed)[k].mother(m)->pdgId() == genPart->pdgId() ) idmatch=true;
              }
              if (!idmatch) continue;
              if(this_dR_lgamma<0.3) {
                  gen_fsrset.insert(k);
                  TLorentzVector gamma;
                  gamma.SetPtEtaPhiE((*packed)[k].pt(),(*packed)[k].eta(),(*packed)[k].phi(),(*packed)[k].energy());
                  lep_dressed = lep_dressed+gamma;
              }
          } // Dressed leptons loop

          Lepts.push_back(lep_dressed);
          LeptsStatus.push_back(genPart->status());
          LeptsId.push_back(genPart->pdgId());
          LeptsMom.push_back(genAna.MotherID(&pruned->at(j)));
          LeptsMomMom.push_back(genAna.MotherMotherID(&pruned->at(j)));
          // if(whereBreak == "lepts"){return;};


          TLorentzVector thisLep;
          thisLep.SetPtEtaPhiM(lep_dressed.Pt(),lep_dressed.Eta(),lep_dressed.Phi(),lep_dressed.M());
          // GEN iso calculation
          double this_GENiso=0.0;
          double this_GENneutraliso=0.0;
          double this_GENchargediso=0.0;
          for(size_t k=0; k<packed->size();k++){
              if( (*packed)[k].status() != 1 ) continue; // stable particles only
              if (abs((*packed)[k].pdgId())==12 || abs((*packed)[k].pdgId())==14 || abs((*packed)[k].pdgId())==16) continue; // exclude neutrinos
              if ((abs((*packed)[k].pdgId())==11 || abs((*packed)[k].pdgId())==13)) continue; // exclude leptons
              if (gen_fsrset.find(k)!=gen_fsrset.end()) continue; // exclude particles which were selected as fsr photons
              double this_dRvL = deltaR(thisLep.Eta(), thisLep.Phi(), (*packed)[k].eta(), (*packed)[k].phi());
              if(this_dRvL<0.3) {
                  this_GENiso = this_GENiso + (*packed)[k].pt();
                  if ((*packed)[k].charge()==0) this_GENneutraliso = this_GENneutraliso + (*packed)[k].pt();
                  if ((*packed)[k].charge()!=0) this_GENchargediso = this_GENchargediso + (*packed)[k].pt();
              }
          } // GEN iso loop
          this_GENiso = this_GENiso/thisLep.Pt();
          Lepts_RelIso.push_back(this_GENiso);
          // END GEN iso calculation
      } // leptons
      if (genPart->pdgId()==25) {
          v.SetPtEtaPhiM(genPart->pt(), genPart->eta(), genPart->phi(), genPart->mass());
          theHiggs.push_back(v);
      }
      if ((genPart->pdgId()==23 || genPart->pdgId()==443 || genPart->pdgId()==553) && (genPart->status()>=20 && genPart->status()<30) ) {
          const reco::Candidate *Zdau0=genPart->daughter(0);
          int ZdauId = fabs(Zdau0->pdgId());
          if (fabs(Zdau0->pdgId())==23) {
              int ndau = genPart->numberOfDaughters();
              for (int d=0; d<ndau; d++) {
                  const reco::Candidate *Zdau=genPart->daughter(d);
                  // if (verbose) cout<<"ZDau "<<d<<" id "<<fabs(Zdau->pdgId())<<endl;
                  if (fabs(Zdau->pdgId())<17) {
                      ZdauId = fabs(Zdau->pdgId());
                      break;
                  }
              }
          }
          if (Zdau0) theZsDaughters.push_back(ZdauId);
          theZsMom.push_back(genAna.MotherID(&pruned->at(j)));
          v.SetPtEtaPhiM(genPart->pt(), genPart->eta(), genPart->phi(), genPart->mass());
          theZs.push_back(v);
          // if (Zdau0) _GENZ_DaughtersId.push_back(ZdauId);
          // _GENZ_MomId.push_back(genAna.MotherID(&pruned->at(j)));
          // _GENZ_pt.push_back(genPart->pt());
          // _GENZ_eta.push_back(genPart->eta());
          // _GENZ_phi.push_back(genPart->phi());
          // _GENZ_mass.push_back(genPart->mass());
      }

      // if (abs(genPart->pdgId())>500 && abs(genPart->pdgId())<600 && genPart->status()==2) {
      //     _nGENStatus2bHad+=1;
      // }
  }

  if (Lepts.size()>=4) {

      unsigned int L1_nocuts=99; unsigned int L2_nocuts=99; unsigned int L3_nocuts=99; unsigned int L4_nocuts=99;
      bool passedFiducialSelectionNoCuts = mZ1_mZ2_ext(L1_nocuts, L2_nocuts, L3_nocuts, L4_nocuts, false);
      if (passedFiducialSelectionNoCuts) {
          TLorentzVector Z1_1, Z1_2, Z2_1, Z2_2;
          Z1_1.SetPtEtaPhiM(Lepts.at(L1_nocuts).Pt(),Lepts.at(L1_nocuts).Eta(),Lepts.at(L1_nocuts).Phi(),Lepts.at(L1_nocuts).M());
          Z1_2.SetPtEtaPhiM(Lepts.at(L2_nocuts).Pt(),Lepts.at(L2_nocuts).Eta(),Lepts.at(L2_nocuts).Phi(),Lepts.at(L2_nocuts).M());
          Z2_1.SetPtEtaPhiM(Lepts.at(L3_nocuts).Pt(),Lepts.at(L3_nocuts).Eta(),Lepts.at(L3_nocuts).Phi(),Lepts.at(L3_nocuts).M());
          Z2_2.SetPtEtaPhiM(Lepts.at(L4_nocuts).Pt(),Lepts.at(L4_nocuts).Eta(),Lepts.at(L4_nocuts).Phi(),Lepts.at(L4_nocuts).M());
          // _GENdPhiZZ = deltaPhi((Z1_1+Z1_2).Phi(),(Z2_1+Z2_2).Phi());
          // _GENmassZZ = (Z1_1+Z1_2+Z2_1+Z2_2).M();
          // _GENpTZZ = (Z1_1+Z1_2+Z2_1+Z2_2).Pt();
      }
  }

  /////// DO THE FIDUCIAL VOLUME CALCULATION //////////////
  passedFiducial=false;
  int nFiducialLeptons = 0;
  int nFiducialPtLead=0;
  int nFiducialPtSublead=0;

  for (unsigned int i=0; i<Lepts.size(); ++i) {

      TLorentzVector thisLep;
      thisLep.SetPtEtaPhiM(Lepts.at(i).Pt(),Lepts.at(i).Eta(),Lepts.at(i).Phi(),Lepts.at(i).M());

      if ( ( (abs(LeptsId.at(i)) == 13 && thisLep.Pt() > 5.0 && abs(thisLep.Eta()) < 2.4)
             || (abs(LeptsId.at(i)) == 11 && thisLep.Pt() > 7.0 && abs(thisLep.Eta()) < 2.5) )
           && Lepts_RelIso.at(i)<0.35) {
          nFiducialLeptons++;
          if (thisLep.Pt()>20) nFiducialPtLead++;
          if (thisLep.Pt()>10) nFiducialPtSublead++;
      }
  }
  if (nFiducialLeptons>=4 && nFiducialPtLead>=1 && nFiducialPtSublead>=2) {

      // START FIDUCIAL EVENT TOPOLOGY CUTS
      unsigned int L1=99; unsigned int L2=99; unsigned int L3=99; unsigned int L4=99;
      passedFiducial = mZ1_mZ2_ext(L1, L2, L3, L4, true);

      Lep_Hindex_tmp[0] = L1; Lep_Hindex_tmp[1] = L2; Lep_Hindex_tmp[2] = L3; Lep_Hindex_tmp[3] = L4;
      for(int i=0; i<4; i++){
      Lep_Hindex.push_back(Lep_Hindex_tmp[i]);
      }

      if (passedFiducial) {

          TLorentzVector LS3_Z1_1, LS3_Z1_2, LS3_Z2_1, LS3_Z2_2;
          LS3_Z1_1.SetPtEtaPhiM(Lepts[L1].Pt(),Lepts[L1].Eta(),Lepts[L1].Phi(),Lepts[L1].M());
          LS3_Z1_2.SetPtEtaPhiM(Lepts[L2].Pt(),Lepts[L2].Eta(),Lepts[L2].Phi(),Lepts[L2].M());
          LS3_Z2_1.SetPtEtaPhiM(Lepts[L3].Pt(),Lepts[L3].Eta(),Lepts[L3].Phi(),Lepts[L3].M());
          LS3_Z2_2.SetPtEtaPhiM(Lepts[L4].Pt(),Lepts[L4].Eta(),Lepts[L4].Phi(),Lepts[L4].M());

          theLepts = {LS3_Z1_1,LS3_Z1_2,LS3_Z2_1,LS3_Z2_2};
          theLeptsId = {LeptsId[L1],LeptsId[L2],LeptsId[L3],LeptsId[L4]};
      }

      bool passedMassOS = true; bool passedElMuDeltaR = true; bool passedDeltaR = true;
      unsigned int N=Lepts.size();
      for(unsigned int i = 0; i<N; i++) {
          for(unsigned int j = i+1; j<N; j++) {

              // only consider the leptons from Z1 and Z2
              if (!(i==L1 || i==L2 || i==L3 || i==L4)) continue;
              if (!(j==L1 || j==L2 || j==L3 || j==L4)) continue;

              TLorentzVector li, lj;
              li.SetPtEtaPhiM(Lepts[i].Pt(),Lepts[i].Eta(),Lepts[i].Phi(),Lepts[i].M());
              lj.SetPtEtaPhiM(Lepts[j].Pt(),Lepts[j].Eta(),Lepts[j].Phi(),Lepts[j].M());

              TLorentzVector mll = li+lj;

              if(LeptsId[i]*LeptsId[j]<0) {
                  if(mll.M()<=4) { passedMassOS = false; break; }
              }

              if(abs(LeptsId[i]) != abs(LeptsId[j])) {
                  double deltaR = li.DeltaR(lj);
                  if(deltaR<=0.02) { passedElMuDeltaR = false; break; }
              }
              double deltaRll = li.DeltaR(lj);
              if(deltaRll<=0.02) { passedDeltaR = false; break; }
          }
      }
      if(passedMassOS==false || passedElMuDeltaR==false || passedDeltaR==false) passedFiducial=false;
      // if(whereBreak == "fidsel") {return;}

      if (passedFiducial) {

          // DO GEN JETS
          edm::View<reco::GenJet>::const_iterator genjet;
          for(genjet = genJ->begin(); genjet != genJ->end(); genjet++) {

              double pt = genjet->pt();  double eta = genjet->eta();
              if (pt<30.0 || abs(eta)>4.7) continue;

              bool inDR_pt30_eta4p7 = false;
              unsigned int N=Lepts.size();
              for(unsigned int i = 0; i<N; i++) {
                  //if (GENlep_status[i]!=1) continue;
                  if (!(abs(LeptsId[i])==11 || abs(LeptsId[i])==13)) continue;
                  TLorentzVector genlep;
                  genlep.SetPtEtaPhiM(Lepts[i].Pt(),Lepts[i].Eta(),Lepts[i].Phi(),Lepts[i].M());
                  double dR = deltaR(genlep.Eta(), genlep.Phi(), genjet->eta(),genjet->phi());
                  if(dR<0.4) {
                      inDR_pt30_eta4p7=true;
                  }
              }

              // count number of gen jets which no gen leptons are inside its cone
              if (!inDR_pt30_eta4p7) {
                  v.SetPtEtaPhiM(genjet->pt(), genjet->eta(), genjet->phi(), genjet->mass());
                  theJets_pt30_eta4p7.push_back(v);
                  if (abs(genjet->eta())<2.5) {
                      theJets_pt30_eta2p5.push_back(v);
                  }
              }

          }// loop over gen jets
      } //passedFiducial
  } // 4 fiducial leptons
}


bool GenTools::mZ1_mZ2_ext(unsigned int& L1, unsigned int& L2, unsigned int& L3, unsigned int& L4, bool makeCuts){

    double offshell = 999.0; bool findZ1 = false; bool passZ1 = false;

    L1 = 0; L2 = 0;

    unsigned int N = Lepts.size();

    for(unsigned int i=0; i<N; i++){
        for(unsigned int j=i+1; j<N; j++){


            if((LeptsId[i]+LeptsId[j])!=0) continue;

            TLorentzVector li, lj;
            li.SetPtEtaPhiM(Lepts[i].Pt(),Lepts[i].Eta(),Lepts[i].Phi(),Lepts[i].M());
            lj.SetPtEtaPhiM(Lepts[j].Pt(),Lepts[j].Eta(),Lepts[j].Phi(),Lepts[j].M());


            if (makeCuts) {
                if ( abs(LeptsId[i]) == 13 && (li.Pt() < 5.0 || abs(li.Eta()) > 2.4)) continue;
                if ( abs(LeptsId[i]) == 11 && (li.Pt() < 7.0 || abs(li.Eta()) > 2.5)) continue;
                if ( Lepts_RelIso[i]>0.35) continue;

                if ( abs(LeptsId[j]) == 13 && (lj.Pt() < 5.0 || abs(lj.Eta()) > 2.4)) continue;
                if ( abs(LeptsId[j]) == 11 && (lj.Pt() < 7.0 || abs(lj.Eta()) > 2.5)) continue;
                if ( Lepts_RelIso[j]>0.35) continue;
            }

            TLorentzVector mll = li+lj;

            if(abs(mll.M()-91.1876)<offshell){
                double mZ1 = mll.M();
                L1 = i; L2 = j; findZ1 = true; offshell = abs(mZ1-91.1876);
            }
        }
    }

    TLorentzVector l1, l2;
    l1.SetPtEtaPhiM(Lepts[L1].Pt(),Lepts[L1].Eta(),Lepts[L1].Phi(),Lepts[L1].M());
    l2.SetPtEtaPhiM(Lepts[L2].Pt(),Lepts[L2].Eta(),Lepts[L2].Phi(),Lepts[L2].M());
    TLorentzVector ml1l2 = l1+l2;

    if(ml1l2.M()>40 && ml1l2.M()<120 && findZ1) passZ1 = true;
    if (!makeCuts) passZ1 = true;

    double pTL34 = 0.0; bool findZ2 = false;
    //bool m4lwindow=false; double window_lo=70.0; double window_hi=140.0;

    //cout<<"findZ2"<<endl;
    for(unsigned int i=0; i<N; i++){
        if(i==L1 || i==L2) continue; // can not be the lep from Z1
        for(unsigned int j=i+1; j<N; j++){
            if(j==L1 || j==L2) continue; // can not be the lep from Z1
            if((LeptsId[i]+LeptsId[j])!=0) continue;

            TLorentzVector li, lj;
            li.SetPtEtaPhiM(Lepts[i].Pt(),Lepts[i].Eta(),Lepts[i].Phi(),Lepts[i].M());
            lj.SetPtEtaPhiM(Lepts[j].Pt(),Lepts[j].Eta(),Lepts[j].Phi(),Lepts[j].M());
            TLorentzVector Z2 = li+lj;

            if (makeCuts) {
                if ( abs(LeptsId[i]) == 13 && (li.Pt() < 5.0 || abs(li.Eta()) > 2.4)) continue;
                if ( abs(LeptsId[i]) == 11 && (li.Pt() < 7.0 || abs(li.Eta()) > 2.5)) continue;
                if ( Lepts_RelIso[i]>0.35) continue;

                if ( abs(LeptsId[j]) == 13 && (lj.Pt() < 5.0 || abs(lj.Eta()) > 2.4)) continue;
                if ( abs(LeptsId[j]) == 11 && (lj.Pt() < 7.0 || abs(lj.Eta()) > 2.5)) continue;
                if ( Lepts_RelIso[j]>0.35) continue;
            }

            if ( (li.Pt()+lj.Pt())>=pTL34 ) {
                double mZ2 = Z2.M();
                // if (verbose) cout<<"_GEN mZ2: "<<mZ2<<endl;
                if( (mZ2>12 && mZ2<120) || (!makeCuts) ) {
                    L3 = i; L4 = j; findZ2 = true;
                    pTL34 = li.Pt()+lj.Pt();
                    // if (verbose) cout<<"is the new _GEN cand"<<endl;
                    //if (m4l>window_lo && m4l<window_hi) m4lwindow=true;
                } else {
                    // still assign L3 and L4 to this pair if we don't have a passing Z2 yet
                    if (findZ2 == false) {L3 = i; L4 = j;}
                    //cout<<"is not new _GEN cand"<<endl;
                }
            }

        } // lj
    } // li

    if(passZ1 && findZ2) return true;
    else return false;
}
