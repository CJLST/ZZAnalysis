/** \file
 *
 *  Retrieve information for boson decays to leptons from the genparticle collection, in the following way:
 *  genH is the generated H (if more than one is present, the last one with daughter = Z is taken)
 *  genZs are the Zs (required to have parent!=Z to avoid double counting)
 *  genLeps are the leptons, required to have either parent=Z or (for generators where the Z 
 * is not explicit in the MC history) status = 3. beware of how FSR is described in the MC history...
 *
 *
 *  $Date: 2013/10/25 15:26:57 $
 *  $Revision: 1.19 $
 *  \author N. Amapane - CERN
 */

#include <ZZAnalysis/AnalysisStep/interface/MCHistoryTools.h>

#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>

#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>

using namespace std;
using namespace reco;
using namespace edm;

namespace {
  bool dbg = false;
}

MCHistoryTools::MCHistoryTools(const edm::Event & event) :
  ismc(false),
  processID(0),
  hepMCweight(1),
  isInit(false),
  theGenH(0) {
  //  if(event.getByLabel("genParticles", particles)){  // genParticles are not available in cmgTuple, only in PAT
  if(event.getByLabel("genParticlesPruned", particles)){
    ismc=true;
    
    edm::Handle<GenEventInfoProduct> gen;
    event.getByLabel( "generator", gen );
    processID = gen->signalProcessID();

//   Process IDs for current samples (Fall11/Summer12) 
//   Generally corresopond to MSUB for Pythia samples, cf. for example: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/Configuration/GenProduction/python/EightTeV/WH_ZH_TTH_HToZZTo4L_M_115_TuneZ2star_8TeV_pythia6_tauola_cff.py?view=markup
//   
//   0-5,35 = DYjets; also WZJets in Fall11
//   0,1,2 = ZZJetsTo4L (MadGraph)
//   661 = GluGluToZZ (gg2zz)
//   10011 = GluGluToHToZZTo4L_M-*_8TeV-powheg-pythia6; GluGluToHToZZTo4L_M-*_mll1_7TeV-powheg-pythia6; VBF_HToZZTo4L_M-*_8TeV-powheg-pythia6, VBF_ToHToZZTo4L_M-*_7TeV-powheg-pythia6
//   11113 = ZZ2e2mu
//   11115 = ZZ2e2tau
//   11315 = ZZ2mu2tau
//   11111 = ZZTo4e
//   11313 = ZZTo4mu
//   11515 = ZZTo4tau
//   322200 = TT
//   23     = WZ (in Summer12)
//   24  = ZH production
//   26  = WH production
//   121 = gg to ttH
//   122 = qq to ttH

//   100 = JHU samples
//   9999 = HH
  
//   take the MC weight
      GenEventInfoProduct  genInfo = *(gen.product());
      hepMCweight = genInfo.weight();
      
  }
}

MCHistoryTools::~MCHistoryTools(){}


// Find the actual lepton parent (first parent in MC history with a different pdgID)
const reco::GenParticle* MCHistoryTools::getParent(const reco::GenParticle* genLep) {
  if(genLep==0) return 0; 
  int flavor = genLep->pdgId();
  
  while (genLep->mother()!=0 && genLep->mother()->pdgId() == flavor) {
    //cout  << " getparent " << genLep->mother()->pdgId();
    genLep = (const GenParticle*) genLep->mother();
  }
  //cout  << " getparent:ret: " << genLep->mother()->pdgId() << endl;
  return (const GenParticle*) genLep->mother();
}

// Same as the above, but try recovery in case default PAT matching fails.
// (This happens rather often for electron due to brems/FSR, since the default matching handles this poorly).
// The recovery consists in matching between a selected list of genleptons, cf. getMatch().
const reco::GenParticle* MCHistoryTools::getParent(const pat::Electron* lep, const vector<const Candidate *>& gen4lep) {
  const reco::GenParticle* parent= getParent((lep->genParticleRef()).get());
  if (parent==0) {
    parent=getParent(getMatch(lep, gen4lep));
  }
  return parent;
}

// Manual matching with closest same-flavour gen lepton (of any status). 
// This was test to work great when the provided candidates are e.g. only the signal ones.
//FIXME: check and charge!
const reco::GenParticle* MCHistoryTools::getMatch(const pat::Electron* lep, const vector<const Candidate *>& gen4lep) {
  float mindeltaR=9999;
  const Candidate * genmatch = 0;
  int lcode = lep->pdgId();
  for (vector<const Candidate *>::const_iterator iglep = gen4lep.begin(); iglep!=gen4lep.end(); ++iglep) {
    if ((*iglep)->pdgId() == lcode) {
      float dR = deltaR(*lep,*(*iglep));
      if (dR<mindeltaR && dR < 0.1) {
	mindeltaR=dR;
	genmatch = *iglep;
      }
    }
  }
  if(dbg)  cout << "Recovery match: dR: " << mindeltaR << " " << lep->pt() << " " << ((genmatch!=0)?(genmatch->pt()):-1) << endl;
  return (const GenParticle*) genmatch;
}

//Return the code of the particle's parent: 25 for H->Z->l; 23 for Z->l; +-15 for tau->l if genlep is e,mu.
int MCHistoryTools::getParentCode(const reco::GenParticle* genLep) {
  int parentId = 0;
  
  const Candidate* particle = getParent(genLep);
  
  if (particle) {
    parentId = particle->pdgId();
    //     cout << "getParentCode1 : " << parentId;
    if (parentId == 23 && particle->mother()!=0) {
      if (particle->mother()->pdgId() == 25) parentId = 25;
    }
  }
  //   cout <<  " getParentCode1: result " <<  parentId << endl;
  return parentId;
}

// Same as above, but if no match is found, search for a match within gen4lep
// Cf. getParent(lep, gen4lep) for details.
int MCHistoryTools::getParentCode(const pat::Electron* lep, const vector<const Candidate *>& gen4lep) {
  int parentId = getParentCode((lep->genParticleRef()).get());
  
  if (parentId==0) { // Recover for bad matching to status 1 electrons due to brems
    const Candidate * particle = getParent(lep, gen4lep);
    if (particle) {
      parentId = particle->pdgId();
      //       cout << "getParentCode2 : " << parentId;
      if (parentId == 23 && particle->mother()!=0) {
	if (particle->mother()->pdgId() == 25) parentId = 25;
      }
    }
  }
  //   cout <<  " getParentCode2: result " <<  parentId << endl;
  return parentId;
}


void
MCHistoryTools::init() {
  if (isInit) return;
  if (!ismc) return;  

//   int nMu=0;
//   int nEle=0;
//   int nTau=0;

  for( View<Candidate>::const_iterator p = particles->begin(); p != particles->end(); ++ p ) {
    int id = abs(p->pdgId());


    if (dbg){
      if (id==13 || id==11 || id ==23 || id==25) {
	cout << "Genpart: id " << id << " pt " << p->pt() << " eta " << p->eta() << " phi " << p->phi()
	     << " status " << p->status()
	     << " parent id " << (p->mother()!=0?p->mother()->pdgId():0)
	  // << " vertex " << p->vertex()
	     << endl;
      }
    }

    //--- H
    if (id==25) {
      if (theGenH==0 || (p->daughter(0)->pdgId()==25)) { // Handle HH samples - genH will be the H decaying to ZZ in this case
	theGenH = &*p;
      }

    //--- Z
    } else if (id==23 && p->mother()!=0 && p->mother()->pdgId()!=23) {// avoid double counting
      if (processID==24 && p->mother()->pdgId()!=25) { //This is an associated Z
	theAssociatedV.push_back(&*p);
      } else { //ZZ or H->ZZ
	theGenZ.push_back(&*p);
      }

    // --- H in JHU samples      
    } else if (id==39 && processID==100){ 
      theGenH = &*p;  
    }
    
    // Lepton (as first daughter of a Z, or status = 3 for ggZZ as the Z are not present in the MC History of ggZZ Summer 12 samples; note that 
    // some special samples (HH) don't use status=3, although that would pick the original leptons in most samples (which could be interesting at some 
    // point to have the leptons in WZ and TT samples)
    else if ((id== 13 || id==11 || id==15) && ((p->mother()!=0 && p->mother()->pdgId()==23) || (processID==661 && p->status()==3))) {
      if (processID==24 && p->mother()!=0 && p->mother()->mother()!=0 && p->mother()->mother()->pdgId()!=25) continue; // ZH: skip leptons from associated Z
      theGenLeps.push_back(&*p);
    }
  } // end loop on particles

  if (theGenLeps.size()!=theGenZ.size()*2) {
    if (processID==661) {// ggZZ samples in Summer12 miss Zs in the MC history.
      // We could build Zs here, at least for 2e2mu where this is unproblematic
    } else if (processID==24 || processID==26 || processID==121 || processID==122) { 
      // For 2012, VH/ttH samples are inclusive in Z decays.
    } else if  (processID==0 || processID==1 || processID==2){
      // For ZZJetsTo4L (MadGraph) samples which contain events with 4 leptons and only 1 Z.
    } else {
      cout << "ERROR: MCHistoryTools::init: unexpected genparticle content for processID= " << processID << " : " << theGenLeps.size() << " " << theGenZ.size() << endl;
      abort();
    }
  }
  isInit = true;


  if (dbg) {
    cout << "MCHistoryTools: "  << processID << " " << genFinalState() << " " << (theGenH==0) << " " << theGenZ.size() << " " << theGenLeps.size() // << " " << nMu << " " << nEle << " " << nTau 
	 << endl;
  }

}

void
MCHistoryTools::genAcceptance(bool& gen_ZZInAcceptance, bool& gen_ZZ4lInEtaAcceptance, bool& gen_ZZ4lInEtaPtAcceptance, bool& gen_m4l_180){
  if (!ismc) return;  
  init();

  float gen_mZ1 = -1.;
  float gen_mZ2 = -1.;

  gen_ZZInAcceptance = false;
  gen_ZZ4lInEtaAcceptance = false;
  gen_ZZ4lInEtaPtAcceptance = false;
  gen_m4l_180 = false;
  int gen_Z1_flavour =0;
  int gen_Z2_flavour =0;
  float gen_4leptonsMass =-1.;

  const float ZmassValue = 91.1876;

  if (theGenZ.size()==2) {
    gen_mZ1 = theGenZ[0]->p4().mass(); // FIXME should take the 2 gen l with status 1!
    gen_mZ2 = theGenZ[1]->p4().mass();
    gen_Z1_flavour = abs(theGenZ[0]->daughter(0)->pdgId());
    gen_Z2_flavour = abs(theGenZ[1]->daughter(0)->pdgId());

    if ( fabs(ZmassValue - gen_mZ2) < fabs(ZmassValue - gen_mZ1)) {
      swap(gen_mZ1,gen_mZ2);
      swap(gen_Z1_flavour, gen_Z2_flavour);
    }

    if (gen_mZ1>60. && gen_mZ1<120. && gen_mZ2>60. && gen_mZ2<120.) {
      gen_ZZInAcceptance = true;
    }


    if ( gen_Z1_flavour<15 && gen_Z2_flavour<15) {    
      gen_ZZ4lInEtaAcceptance = true;
      gen_ZZ4lInEtaPtAcceptance = true;


      if (theGenLeps.size()>=4) {
	for (int i=0; i<4; ++i){	  
	  //FIXME should take the 2 gen l with status 1!
	  if ((abs(theGenLeps[i]->pdgId()) == 11 && !(theGenLeps[i]->pt() > 7. && fabs(theGenLeps[i]->eta()) < 2.5)) ||
	      (abs(theGenLeps[i]->pdgId()) == 13 && !(theGenLeps[i]->pt() > 5. && fabs(theGenLeps[i]->eta()) < 2.4))) { 
	    gen_ZZ4lInEtaPtAcceptance = false;
	  }
	  if ((abs(theGenLeps[i]->pdgId()) == 11 && !(fabs(theGenLeps[i]->eta()) < 2.5)) ||
	      (abs(theGenLeps[i]->pdgId()) == 13 && !(fabs(theGenLeps[i]->eta()) < 2.4))) { 
	    gen_ZZ4lInEtaAcceptance = false;
	  }
	}
      }
    }
  }

  if (theGenLeps.size()==4){
    gen_4leptonsMass = (theGenLeps[0]->p4()+theGenLeps[1]->p4()+theGenLeps[2]->p4()+theGenLeps[3]->p4()).mass();
    if (gen_4leptonsMass >= 180){
      gen_m4l_180 = true;
    }
  }
}

int
MCHistoryTools::genFinalState(){
  if (!ismc) return -1;  
  init();
  
  if (theGenH!=0 && theGenZ.size()!=2) {
    // cout << "ERROR: MCHistoryTools: genH!=0 but genZ.size()==" << theGenZ.size() << endl;
    if (abs(theGenH->daughter(0)->pdgId())<10) {
      // This can happen due to a known problem. cf https://hypernews.cern.ch/HyperNews/CMS/get/generators/1405.html
      // abort();
      return BUGGY;
    }
  }

  int gen_finalState = NONE;  
  if (theGenZ.size()==2){
    int gen_Z1_flavour = abs(theGenZ[0]->daughter(0)->pdgId());
    int gen_Z2_flavour = abs(theGenZ[1]->daughter(0)->pdgId());

    if (gen_Z1_flavour == 11 && gen_Z2_flavour == 11) {
      gen_finalState = EEEE;
    } else if (gen_Z1_flavour == 13 && gen_Z2_flavour == 13) {
      gen_finalState = MMMM;
    } else if ((gen_Z1_flavour == 11 && gen_Z2_flavour == 13) || 
	       (gen_Z1_flavour == 13 && gen_Z2_flavour == 11)) {
      gen_finalState = EEMM;
    } else if (gen_Z1_flavour == 15 || gen_Z2_flavour == 15) {
      gen_finalState = LLTT;
    } else if (processID==24 || processID==26 || processID==121 || processID==122) { // ZH 8TeV samples are inclusive in Z decays.
      return NONE;
    } else {
      cout << "ERROR: MCHistoryTools: processID: " << processID << " Z flavour= " << gen_Z1_flavour << " " << gen_Z2_flavour << endl;
      abort();
    }
  } else if (theGenZ.size()==0 && theGenLeps.size()==4 && processID==661) {
    // Handle samples where Zs are not explicit in the MC history
    int nele=0;
    int nmu=0;
    int ntau=0;
    for (int i=0; i<4; ++i){    
      int id = abs(theGenLeps[i]->pdgId());
      if (id==11) nele++;
      if (id==13) nmu++;
      if (id==15) ntau++;
    }
    if (nele==4)  gen_finalState = EEEE;
    else if (nmu==4)  gen_finalState = MMMM;
    else if (nmu==2 && nele==2) gen_finalState = EEMM;
    else if (ntau==2 || ntau==4) gen_finalState = LLTT;
    else {
      cout << "ERROR: MCHistoryTools: leptons: " << nele << " " << nmu << " " << ntau << endl;
      abort();
    }
  } 
  return gen_finalState;
}

int
MCHistoryTools::genAssociatedFS(){
  if (!ismc) return -1;  
  init();

  int id=0;
  if (processID==24){
    if (theAssociatedV.size()!=1) {
      cout << "ERROR: ZH with " << theAssociatedV.size() << " associated V" << endl;  
    } else {
      id = abs(theAssociatedV.front()->daughter(0)->pdgId());
    }
  }
  
  return id;
}


/*
  vector<double> genPt(4);

  float gen_m4l = -1.; // defined for ZZ for the time being
  float gen_mZ1 = -1.;
  float gen_mZ2 = -1.;
  float gen_Z1_flavour =0;
  float gen_Z2_flavour =0;
  bool gen_ZZInAcceptance = false;
  bool gen_ZZ4lInEtaPtAcceptance = false;
  bool gen_ZZ4lInEtaAcceptance = false;
  int gen_finalState = -1;
  float genWeight = 1.;

    
    if (isZZ && genZ.size()!=2) {
      cout << "ERROR:# Z= " << genZ.size() <<endl;
      abort();
    }
		
    // ZZ or H->ZZ
    if (genZ.size()==2) {
      
      gen_m4l = (genZ[0]->p4()+ // FIXME should take the 4 gen l with status 1!
		 genZ[1]->p4()).mass(); // checked to be the same of theGenH->p4().mass().
      gen_mZ1 = genZ[0]->p4().mass(); // FIXME should take the 2 gen l with status 1!
      gen_mZ2 = genZ[1]->p4().mass();
      gen_Z1_flavour = abs(genZ[0]->daughter(0)->pdgId());
      gen_Z2_flavour = abs(genZ[1]->daughter(0)->pdgId());
      
      // Sort Z1/Z2 according to standard definition
      if ( fabs(ZmassValue - gen_mZ2) < fabs(ZmassValue - gen_mZ1)) {
	swap(gen_mZ1,gen_mZ2);
	swap(gen_Z1_flavour, gen_Z2_flavour);
      }
      hmZ1vsmZ2->Fill(gen_mZ1,gen_mZ2);


      // GEN acceptance
      if (gen_mZ1>60. && gen_mZ1<120. && gen_mZ2>60. && gen_mZ2<120.) {
	gen_ZZInAcceptance = true;
      }
      
      if ( gen_Z1_flavour<15 && gen_Z2_flavour<15) {
	if (gen4lep.size()!=4) {
	  cout << "ERROR:# gen4l = " << gen4lep.size() <<endl;
	  abort();
	}
	
	
	gen_ZZ4lInEtaPtAcceptance = true;
	gen_ZZ4lInEtaAcceptance = true;
	for (int i=0; i<4; ++i){
	  
	  genPt[i] = gen4lep[i]->pt();
	  
	  //FIXME should take the 2 gen l with status 1!
	  if ((abs(gen4lep[i]->pdgId()) == 11 && !(gen4lep[i]->pt() > 7. && fabs(gen4lep[i]->eta()) < 2.5)) ||
	      (abs(gen4lep[i]->pdgId()) == 13 && !(gen4lep[i]->pt() > 5. && fabs(gen4lep[i]->eta()) < 2.4))) { 
	    gen_ZZ4lInEtaPtAcceptance = false;
	  }
	  if ((abs(gen4lep[i]->pdgId()) == 11 && !(fabs(gen4lep[i]->eta()) < 2.5)) ||
	      (abs(gen4lep[i]->pdgId()) == 13 && !(fabs(gen4lep[i]->eta()) < 2.4))) { 
	    gen_ZZ4lInEtaAcceptance = false;
	  }
	}   
	
	sort(genPt.begin(),genPt.end());	
	hgenPt1->Fill(genPt[0]);
	hgenPt2->Fill(genPt[1]);
	hgenPt3->Fill(genPt[2]);
	hgenPt4->Fill(genPt[3]);
	
      }
      
       if (gen_Z1_flavour == 11 && gen_Z2_flavour == 11) {
	gen_finalState = EEEE;
	++gen_ZZ4e;
	if (gen_ZZ4lInEtaAcceptance) ++gen_ZZ4e_EtaAcceptance;
	if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ4e_LeptonAcceptance;
	if (gen_ZZInAcceptance) {
	  ++gen_ZZ4e_MassAcceptance;
	  if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ4e_MassPtAcceptance;
	}
      } else if (gen_Z1_flavour == 13 && gen_Z2_flavour == 13) {
	gen_finalState = MMMM;
	++gen_ZZ4mu;
	if (gen_ZZ4lInEtaAcceptance) ++gen_ZZ4mu_EtaAcceptance;
	if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ4mu_LeptonAcceptance;
	if (gen_ZZInAcceptance) {
	  ++gen_ZZ4mu_MassAcceptance;
	  if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ4mu_MassPtAcceptance;
	}
      } else if ((gen_Z1_flavour == 11 && gen_Z2_flavour == 13) || 
		 (gen_Z1_flavour == 13 && gen_Z2_flavour == 11)) {
	gen_finalState = EEMM;
	++gen_ZZ2mu2e;
	if (gen_ZZ4lInEtaAcceptance) ++gen_ZZ2mu2e_EtaAcceptance;
	if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ2mu2e_LeptonAcceptance;
	if (gen_ZZInAcceptance) {
	  ++gen_ZZ2mu2e_MassAcceptance;
	  if (gen_ZZ4lInEtaPtAcceptance) ++gen_ZZ2mu2e_MassPtAcceptance;
	}
      } else if (gen_Z1_flavour == 15 || gen_Z2_flavour == 15) {
	gen_finalState = 10;
	++gen_ZZ2l2tau;
      } else {
	// something wrong?
	cout << "ERROR: Z flavour= " << gen_Z1_flavour << " " << gen_Z2_flavour << endl;
	abort();
      }
    }
  }
}
*/
