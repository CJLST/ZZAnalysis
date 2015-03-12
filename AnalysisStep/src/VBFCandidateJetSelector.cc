#include "../interface/VBFCandidateJetSelector.h"
#include <iostream>

bool debug = false;

std::vector<const pat::Jet*> 
VBFCandidateJetSelector::cleanJets(std::vector<const reco::Candidate*> leptons, 
				   edm::Handle<edm::View<pat::Jet> > jets, int year)
{
  for(edm::View<pat::Jet>::const_iterator j = jets->begin(); j != jets->end(); ++j) {

    float jeta=fabs(j->eta());
                
          //            bool looseJetID = (j->getSelection("cuts_looseJetId") > 0) ;
          /*
          bool looseJetID = (j->component(5).fraction() < 0.99 && 
                             j->component(4).fraction() < 0.99 && 
                             j->nConstituents() > 1 && 
                             ( j->component(1).fraction() > 0 || jeta > 2.4 )  &&
                             ( ( j->component(1).number() + j->component(2).number() + j->component(3).number() ) > 0 || jeta > 2.4 ) &&
                             ( j->component(2).fraction() < 0.99 || jeta > 2.4 ) );
          */    
          bool looseJetID = true;            
          bool passPU = true; //j->passPuJetId("full", PileupJetIdentifier::kLoose);
          
                //HARD CODED implementation of JetMET V00-03-04 WPs - for synch only
                //#4 Eta Categories  0-2.5 2.5-2.75 2.75-3.0 3.0-5.0
                //Pt010_Loose    = cms.vdouble(-0.95,-0.96,-0.94,-0.95),
                //Pt1020_Loose   = cms.vdouble(-0.95,-0.96,-0.94,-0.95),
                //Pt2030_Loose   = cms.vdouble(-0.63,-0.60,-0.55,-0.45),
                //Pt3050_Loose   = cms.vdouble(-0.63,-0.60,-0.55,-0.45), 
                float jpt=j->pt();
                float jpumva=0.;
                //if(year==2012)jpumva=j->puMva("full53x");
                //else jpumva=j->puMva("full");
                jpumva=j->userFloat("pileupJetId:fullDiscriminant");
                if(jpt>20){
                  if(jeta>3.){
                    if(jpumva<=-0.45)passPU=false;
                  }else if(jeta>2.75){
                    if(jpumva<=-0.55)passPU=false;
                  }else if(jeta>2.5){
                    if(jpumva<=-0.6)passPU=false;
                  }else if(jpumva<=-0.63)passPU=false;
                }else{
                  if(jeta>3.){
                    if(jpumva<=-0.95)passPU=false;
                  }else if(jeta>2.75){
                    if(jpumva<=-0.94)passPU=false;
                  }else if(jeta>2.5){
                    if(jpumva<=-0.96)passPU=false;
                  }else if(jpumva<=-0.95)passPU=false;
                }

    bool isDeltaR = true;
                
    if ( looseJetID && passPU ) {

      for(unsigned int l=0; l<leptons.size(); ++l){
        const reco::Candidate* myLep = leptons.at(l); 
        if (reco::deltaR(j->p4(),myLep->p4()) < 0.4)
          isDeltaR = false;
      }

      if (isDeltaR && j->pt()>20 && jeta<4.7)
        selected_.push_back(&(*j)) ;
    }

    if (debug) std::cout << "VBF: " << j->pt() << " " << looseJetID << " " << passPU << " " << isDeltaR << " " << jpumva << std::endl;

  }

  return selected_;

}
