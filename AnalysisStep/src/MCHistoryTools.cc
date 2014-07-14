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

#include <boost/algorithm/string/predicate.hpp>

using namespace std;
using namespace reco;
using namespace edm;

namespace {
  bool dbg = false;
}

MCHistoryTools::MCHistoryTools(const edm::Event & event, string sampleName) :
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
//   0,1,2 = ZZJetsTo4L (MadGraph); ZZZ (Madgraph)
//   66  =  WZZ_8TeV-aMCatNLO-herwig
//   661 = GluGluToZZ (gg2zz);  phantom samples
//   10011 = GluGluToHToZZTo4L_M-*_8TeV-powheg-pythia6; GluGluToHToZZTo4L_M-*_mll1_7TeV-powheg-pythia6; VBF_HToZZTo4L_M-*_8TeV-powheg-pythia6, VBF_ToHToZZTo4L_M-*_7TeV-powheg-pythia6
//   11113 = ZZ2e2mu
//   11115 = ZZ2e2tau
//   11315 = ZZ2mu2tau
//   11111 = ZZTo4e
//   11313 = ZZTo4mu
//   11515 = ZZTo4tau
//   10131 = ggTo2e2mu_BSMHContinInterf-MCFM67
//   10132 = ggTo2e2mu_Contin-MCFM67
//   322200 = TT
//   23     = WZ (in Summer12)
//   24  = ZH production
//   26  = WH production
//   121 = gg to ttH
//   122 = qq to ttH
//   9999 = HH, powheg15jhuGenV3 samples

//   100 = old JHU samples AND gg2zz samples with off-shell Higgs (ggTo2l2l_Continuum, ggTo2l2l_H, ggTo2l2l_ContinuumInterfH, same with ggTo4l)
    // We override the processID in this case, based on sampleName
    if (processID == 100) {
      // FIXME: fix gg2ZZ samples
      if (boost::starts_with(sampleName,"ZHiggs")) processID=900024;
      if (boost::starts_with(sampleName,"WHiggs")) processID=900026; 
      if (boost::starts_with(sampleName,"ggTo"))   processID=900661; 
    }

    if (processID == 0) {
      if (boost::starts_with(sampleName,"ZZZJets")) processID=900101;      
    }
    

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
      if ((processID==24 || processID==900024) && p->mother()->pdgId()!=25) { //This is an associated Z
	theAssociatedV.push_back(&*p);
      } else { //ZZ or H->ZZ
	theGenZ.push_back(&*p);
      }

    // --- H in some JHU samples has id=39 instead of 25
    } else if (id==39 && processID==100){ 
      theGenH = &*p;  
    }
    
    // Lepton (as first daughter of a Z, or status = 3 for ggZZ and phantom as the Zs are not present 
    // in the MC History of ggZZ Summer 12 samples, or leptons may not come from the Z
    else if ((id== 13 || id==11 || id==15) && ((p->mother()!=0 && p->mother()->pdgId()==23) || ((processID==661||processID==900661) && p->status()==3))) { 

      // ZH: skip leptons from associated Z
      if ((processID==24||processID==900024) && p->mother()!=0 && p->mother()->mother()!=0 && p->mother()->mother()->pdgId()!=25) {
	// Can save these to a separate list of leptons from associated prod
	continue;
      }
      
      theGenLeps.push_back(&*p);
    }
  } // end loop on particles


  bool isOK=true;
  // Check consistency of what we have collected (FIXME: remove Z stuff)
  if (theGenLeps.size()!=theGenZ.size()*2) {
    if (processID==24 || processID==26 || processID==121 || processID==122) {
      // For 2012, VH/ttH samples are inclusive in Z decays, assume everything is fine
    } else if (processID==661 || processID==900661 || processID==0 || processID==1 || processID==2 || processID==66 || processID==900101) {
      // Samples which miss Zs in the MC history, assume everything is fine      
    } else {
      isOK = false;
    }
  }

 
  if (!isOK) {    
      cout << "ERROR: MCHistoryTools::init: unexpected genparticle content for processID= " << processID << " : " << theGenLeps.size() << " " << theGenZ.size() << endl;
      abort();    
  }
  

  // Sort leptons, as done for the signal, for cases where we have 4.
  if (theGenLeps.size()==4) {
    const float ZmassValue = 91.1876;  
    float minDZMass=1E36;
    float iZ11=-1, iZ12=-1, iZ21=-1, iZ22=-1;
    
    // Find Z1 as closest-to-mZ l+l- combination
    for (int i=0; i<4; ++i) {
      for (int j=i+1; j<4; ++j) {
	if (theGenLeps[i]->pdgId()+theGenLeps[j]->pdgId()==0) { // Same flavour, opposite sign
	  float dZMass = std::abs((theGenLeps[i]->p4()+theGenLeps[j]->p4()).mass()-ZmassValue);
	  if (dZMass<minDZMass){
	    minDZMass=dZMass;
	    iZ11=i;
	    iZ12=j;
	  }
	}
      }
    }    

    // Z2 is from remaining 2 leptons
    if (iZ11!=-1 && iZ12!=-1){
      for (int i=0; i<4; ++i) {
	if (i!=iZ11 && i!=iZ12) {
	  if (iZ21==-1) iZ21=i;
	  else iZ22=i;
	}
      }
    }

    if (iZ22==-1 || theGenLeps[iZ21]->pdgId()+theGenLeps[iZ22]->pdgId()!=0) { //Test remaining conditions: Z2 is found and SF, OS
      cout << "MCHistoryTools: Cannot sort leptons ";
      for (int i=0; i<4; ++i) cout << theGenLeps[i]->pdgId() << " ";
      cout << iZ11 << " " << iZ12 << " " << iZ21 << " " << iZ22 << endl;
      abort();
    }
    
    // Sort leptons by sign
    if (theGenLeps[iZ11]->pdgId() < 0 ) {
      swap(iZ11,iZ12);
    }
    if (theGenLeps[iZ21]->pdgId() < 0 ) {
      swap(iZ21,iZ22);
    }

    theSortedGenLepts.push_back(theGenLeps[iZ11]);
    theSortedGenLepts.push_back(theGenLeps[iZ12]);
    theSortedGenLepts.push_back(theGenLeps[iZ21]);
    theSortedGenLepts.push_back(theGenLeps[iZ22]);

//     cout << " Gen Lepton sorting: " << sampleName << " " 
// 	 << (theSortedGenLepts[0]->p4()+theSortedGenLepts[1]->p4()).mass() << " " 
// 	 << (theSortedGenLepts[2]->p4()+theSortedGenLepts[3]->p4()).mass() << " | "
// 	 << (theSortedGenLepts[0]->p4()+theSortedGenLepts[3]->p4()).mass() << " "  
// 	 << (theSortedGenLepts[2]->p4()+theSortedGenLepts[1]->p4()).mass() << " "
// 	 << theSortedGenLepts[0]->pdgId() << " " 
// 	 << theSortedGenLepts[1]->pdgId() << " " 
// 	 << theSortedGenLepts[2]->pdgId() << " " 
// 	 << theSortedGenLepts[3]->pdgId() << " " 
// 	 << endl;
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

  gen_ZZInAcceptance = false; // obsolete

  gen_ZZ4lInEtaAcceptance = false;
  gen_ZZ4lInEtaPtAcceptance = false;
  gen_m4l_180 = false;
  int gen_Z1_flavour =0;
  int gen_Z2_flavour =0;
  float gen_4leptonsMass =-1.;

  const float ZmassValue = 91.1876;

  // This is done using gen Z. Obsolete, to be removed!
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
  }
  

  int nlInEtaAcceptance = 0;
  int nlInEtaPtAcceptance = 0;

  for (unsigned int i=0; i<theGenLeps.size(); ++i){	  
    //FIXME should take the 2 gen l with status 1!
    if ((abs(theGenLeps[i]->pdgId()) == 11 && !(theGenLeps[i]->pt() > 7. && fabs(theGenLeps[i]->eta()) < 2.5)) ||
	(abs(theGenLeps[i]->pdgId()) == 13 && !(theGenLeps[i]->pt() > 5. && fabs(theGenLeps[i]->eta()) < 2.4))) { 
      ++nlInEtaPtAcceptance;
    }
    if ((abs(theGenLeps[i]->pdgId()) == 11 && !(fabs(theGenLeps[i]->eta()) < 2.5)) ||
	(abs(theGenLeps[i]->pdgId()) == 13 && !(fabs(theGenLeps[i]->eta()) < 2.4))) { 
      ++nlInEtaAcceptance;
    }
  }

  if (nlInEtaPtAcceptance>=4) gen_ZZ4lInEtaPtAcceptance = true;
  if (nlInEtaAcceptance>=4) gen_ZZ4lInEtaAcceptance = true;


  // FIXME: still needed?
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

    // FIXME this does not make much sense now that we re-pair Zs in the MC history.
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
    } else if (processID==661) { // in phantom samples, some of the 4 leptons may not be listed as coming from the 2 Z. Proceed to the following logic based on leptons.
      
    } else {
      cout << "ERROR: MCHistoryTools: processID: " << processID << " Z flavour= " << gen_Z1_flavour << " " << gen_Z2_flavour << endl;
      abort();
    }
  } else if (theGenZ.size()==0 && theGenLeps.size()==4 && (processID==661 || processID==900661)) {
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
