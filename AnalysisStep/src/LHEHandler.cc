#include <ZZAnalysis/AnalysisStep/interface/LHEHandler.h>

typedef std::vector<int> vectorInt;
typedef std::vector<std::pair<int, int>> vectorIntPair;

using namespace std;
using namespace PDGHelpers;


LHEHandler::LHEHandler(edm::Handle<LHEEventProduct>* lhe_evt_, int VVMode_, int VVDecayMode_, bool doKinematics_) :
VVMode(VVMode_),
VVDecayMode(VVDecayMode_),
doKinematics(doKinematics_),
genEvent(0),
genCand(0)
{
  setHandle(lhe_evt_);
  extract();
}
LHEHandler::LHEHandler(int VVMode_, int VVDecayMode_, bool doKinematics_) :
VVMode(VVMode_),
VVDecayMode(VVDecayMode_),
doKinematics(doKinematics_),
genEvent(0),
genCand(0)
{
  clear();
}
LHEHandler::~LHEHandler(){ clear(); }


void LHEHandler::setHandle(edm::Handle<LHEEventProduct>* lhe_evt_){ clear(); lhe_evt=lhe_evt; }
void LHEHandler::clear(){
  lhe_evt=0;

  if (genEvent!=0) delete genEvent;
  genCand=0;

  for (unsigned int p=0; p<particleList.size(); p++){
    MELAParticle* tmpPart = (MELAParticle*)particleList.at(p);
    if (tmpPart!=0) delete tmpPart;
  }
  particleList.clear();

  LHEWeight.clear();
  PDFScale=0;
}

MELACandidate* LHEHandler::getBestCandidate(){ return genCand; }
float LHEHandler::getLHEWeight(unsigned int whichWeight, float defaultValue){
  if (whichWeight<LHEWeight.size()) return LHEWeight.at(whichWeight);
  else return defaultValue;
}
float LHEHandler::getPDFScale(){ return PDFScale; }



void LHEHandler::extract(){
  if (lhe_evt!=0){
    if (lhe_evt->isValid()){

      readEvent();

      if (doKinematics){

        genEvent = new LHE_Event();
        vectorInt hasGenHiggs;

        for (unsigned int p=0; p<particleList.size(); p++){
          MELAParticle* genPart = particleList.at(p); // Has mother info from LHE reading
          if (isAHiggs(genPart->id)) hasGenHiggs.push_back(p);
          if (genPart->genStatus==1){
            if (isALepton(genPart->id)) genEvent->addLepton(genPart);
            else if (isANeutrino(genPart->id)) genEvent->addNeutrino(genPart);
            else if (isAPhoton(genPart->id)) genEvent->addPhoton(genPart);
            else if (isAGluon(genPart->id) || isAQuark(genPart->id)) genEvent->addJet(genPart);
          }
        }

        genEvent->constructVVCandidates(VVMode, VVDecayMode);
        for (unsigned int p=0; p<particleList.size(); p++){
          MELAParticle* genPart = particleList.at(p);
          if (genPart->genStatus==-1) genEvent->addVVCandidateMother(genPart);
        }
        genEvent->addVVCandidateAppendages();

        genCand=0;
        if (hasGenHiggs.size()>0){
          for (unsigned int gk=0; gk<hasGenHiggs.size(); gk++){
            MELACandidate* tmpCand = matchAHiggsToParticle(*genEvent, particleList.at(hasGenHiggs.at(gk)));
            if (tmpCand!=0){
              if (genCand==0) genCand=tmpCand;
              else genCand = candComparator(genCand, tmpCand, VVMode);
            }
          }
        }
        if (genCand==0) genCand = candidateSelector(*genEvent, VVMode);

      }
      else{ genCand=0; genEvent=0; }

    }
  }
}

void LHEHandler::readEvent(){
  // Particles
  if (doKinematics){
    const lhef::HEPEUP hepeup_ = (*lhe_evt)->hepeup();
    const int nup = hepeup_.NUP;
    const vectorInt istup = hepeup_.ISTUP;
    const vectorIntPair mothup = hepeup_.MOTHUP;
    const vectorInt idup = hepeup_.IDUP;
    const std::vector<lhef::HEPEUP::FiveVector> pup = hepeup_.PUP;
    vectorInt motherIDs_first;
    vectorInt motherIDs_second;
    for (int ipart = 0; ipart<nup; ipart++){
      SimpleParticle_t simplePart(idup.at(ipart), TLorentzVector(pup.at(ipart)[0], pup.at(ipart)[1], pup.at(ipart)[2], pup.at(ipart)[3]));
      MELAParticle* onePart = new MELAParticle(simplePart.first, simplePart.second);
      onePart->setGenStatus(istup.at(ipart));
      particleList.push_back(onePart);

      motherIDs_first.push_back(mothup.at(ipart).first);
      motherIDs_second.push_back(mothup.at(ipart).second);
    }
    // Link the mothers
    for (int a=0; a<nup; a++){
      if (motherIDs_first.at(a)>0) particleList.at(a)->addMother(particleList.at(motherIDs_first.at(a)-1));
      if (motherIDs_second.at(a)>0 && motherIDs_first.at(a)!=motherIDs_second.at(a)) particleList.at(a)->addMother(particleList.at(motherIDs_second.at(a)-1)); // Refrain from adding the same mothers twice
    }
  }
  //

  // PDF scale
  PDFScale = -1;
  if ((*lhe_evt)->pdf()!=NULL) PDFScale = (*lhe_evt)->pdf()->scalePDF;
  //

  // LHE weights (Maximum 9 or size of weights array)
  for (int iw=0; iw<min(9, (int)((*lhe_evt)->weights().size())); iw++) LHEWeight.push_back((*lhe_evt)->weights().at(iw).wgt / (*lhe_evt)->originalXWGTUP());
  //
}


MELACandidate* LHEHandler::matchAHiggsToParticle(LHE_Event& ev, MELAParticle* genH){
  MELACandidate* cand=0;
  for (int t=0; t<ev.getNZZCandidates(); t++){
    MELACandidate* tmpCand = ev.getZZCandidate(t);
    // FIX ME
    double genhmassquant = genH->m()+genH->pt()+fabs(genH->z());
    double massdiff = fabs(genhmassquant-tmpCand->m()-tmpCand->pt()-fabs(tmpCand->z()));
    double massratio = 0;
    if (genhmassquant>0) massratio = massdiff / genhmassquant;
    if (massratio<0.001){
      if (cand==0) cand = tmpCand;
      else{
        TLorentzVector vGen = genH->p4;
        TLorentzVector vTmp = tmpCand->p4;
        TLorentzVector vCur = cand->p4;

        double dot_tmp = vTmp.Dot(vGen);
        double dot_curr = vCur.Dot(vGen);
        if (fabs(dot_tmp-vGen.M2())<fabs(dot_curr - vGen.M2())) cand = tmpCand;
      }
    }
  }
  return cand;
}

MELACandidate* LHEHandler::candidateSelector(LHE_Event& ev, int isZZ){
  MELACandidate* cand=0;
  for (int t=0; t<ev.getNZZCandidates(); t++){
    MELACandidate* tmpCand = ev.getZZCandidate(t);
    //if (!tmpCand->passSelection) continue;
    if (cand==0) cand=tmpCand;
    else cand = candComparator(cand, tmpCand, isZZ);
  }
  return cand;
}

MELACandidate* LHEHandler::candComparator(MELACandidate* cand1, MELACandidate* cand2, int isZZ){
  MELACandidate* theChosenOne=0;

  double defaultHVVmass = PDGHelpers::HVVmass;
  if (isZZ==0) PDGHelpers::setHVVmass(PDGHelpers::Wmass);
  else PDGHelpers::setHVVmass(PDGHelpers::Zmass);

  double diffmass1 = fabs(cand1->getSortedV(0)->m()-PDGHelpers::HVVmass);
  double diffmass2 = fabs(cand2->getSortedV(0)->m()-PDGHelpers::HVVmass);
  double Z2scsumpt_cand1=0, Z2scsumpt_cand2=0;
  MELAParticle* c11 = cand1->getSortedV(1)->getDaughter(0);
  MELAParticle* c12 = cand1->getSortedV(1)->getDaughter(1);
  MELAParticle* c21 = cand2->getSortedV(1)->getDaughter(0);
  MELAParticle* c22 = cand2->getSortedV(1)->getDaughter(1);
  if (c11!=0) Z2scsumpt_cand1 += c11->pt();
  if (c12!=0) Z2scsumpt_cand1 += c12->pt();
  if (c21!=0) Z2scsumpt_cand2 += c21->pt();
  if (c22!=0) Z2scsumpt_cand2 += c22->pt();
  if (
    (diffmass1>diffmass2)
    ||
    (diffmass1==diffmass2 && Z2scsumpt_cand2>Z2scsumpt_cand1)
    ) theChosenOne = cand2;
  else theChosenOne = cand1;

  PDGHelpers::setHVVmass(defaultHVVmass);
  return theChosenOne;
}


