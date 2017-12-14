#include <ZZAnalysis/AnalysisStep/interface/LHEHandler.h>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include "TLorentzVector.h"
#include "TString.h"

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
  setHandle(lhe_evt_); // Also calls clear()
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


void LHEHandler::setHandle(edm::Handle<LHEEventProduct>* lhe_evt_){ clear(); lhe_evt=lhe_evt_; }
void LHEHandler::clear(){
  lhe_evt=0;

  if (genEvent!=0){ delete genEvent; genEvent=0; }
  genCand=0;

  for (unsigned int p=0; p<particleList.size(); p++){
    MELAParticle* tmpPart = (MELAParticle*)particleList.at(p);
    if (tmpPart!=0) delete tmpPart;
  }
  particleList.clear();

  LHEWeight.clear();
  LHEWeight_PDFVariationUpDn.clear();
  LHEWeight_AsMZUpDn.clear();
  PDFid.clear();
  PDFScale=0;
  LHEOriginalWeight=1;
}

MELACandidate* LHEHandler::getBestCandidate(){ return genCand; }
float const& LHEHandler::getLHEOriginalWeight() const{
  return this->LHEOriginalWeight;
}
float LHEHandler::getLHEWeight(unsigned int whichWeight, float defaultValue) const{
  if (whichWeight<LHEWeight.size()) return LHEWeight.at(whichWeight);
  else return defaultValue;
}
float LHEHandler::getLHEWeight_PDFVariationUpDn(int whichUpDn, float defaultValue) const{
  if (whichUpDn>0 && LHEWeight_PDFVariationUpDn.size()>1) return LHEWeight_PDFVariationUpDn.at(1);
  else if (whichUpDn<0 && LHEWeight_PDFVariationUpDn.size()>0) return LHEWeight_PDFVariationUpDn.at(0);
  else return defaultValue;
}
float LHEHandler::getLHEWeigh_AsMZUpDn(int whichUpDn, float defaultValue) const{
  if (whichUpDn>0 && LHEWeight_AsMZUpDn.size()>1) return LHEWeight_AsMZUpDn.at(1);
  else if (whichUpDn<0 && LHEWeight_AsMZUpDn.size()>0) return LHEWeight_AsMZUpDn.at(0);
  else return defaultValue;
}
float const& LHEHandler::getPDFScale() const{ return PDFScale; }



void LHEHandler::extract(){
  if (lhe_evt==0) cerr << "LHEHandler::extract: lhe_evt==0" << endl;

  if (lhe_evt!=0){
    if (!lhe_evt->isValid()) cerr << "LHEHandler::extract: lhe_evt invalid!" << endl;

    if (lhe_evt->isValid()){

      readEvent();

      if (doKinematics){

        genEvent = new LHE_Event();
        vectorInt hasGenHiggs;

        {
          int p=0;
          for (MELAParticle* genPart:particleList){
            if (isAHiggs(genPart->id)) hasGenHiggs.push_back(p);
            if (genPart->genStatus==1){
              if (isALepton(genPart->id)) genEvent->addLepton(genPart);
              else if (isANeutrino(genPart->id)) genEvent->addNeutrino(genPart);
              else if (isAPhoton(genPart->id)) genEvent->addPhoton(genPart);
              else if (isAGluon(genPart->id) || isAQuark(genPart->id)) genEvent->addJet(genPart);
            }
          }
          p++;
        }

        genEvent->constructVVCandidates(VVMode, VVDecayMode);
        for (MELAParticle* genPart:particleList){ if (genPart->genStatus==-1) genEvent->addVVCandidateMother(genPart); }
        genEvent->addVVCandidateAppendages();

        genCand=nullptr;
        if (hasGenHiggs.size()>0){
          for (int iH:hasGenHiggs){
            MELACandidate* tmpCand = matchAHiggsToParticle(*genEvent, particleList.at(iH));
            if (tmpCand){
              if (!genCand) genCand=tmpCand;
              else genCand = candComparator(genCand, tmpCand, VVMode);
            }
          }
        }
        if (!genCand) genCand = candidateSelector(*genEvent, VVMode);

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
  if ((*lhe_evt)->pdf()!=NULL){
    PDFScale = (*lhe_evt)->pdf()->scalePDF;
    PDFid.push_back((*lhe_evt)->pdf()->id.first);
    if (PDFid.at(0) != (*lhe_evt)->pdf()->id.second) PDFid.push_back((*lhe_evt)->pdf()->id.second);
  }
  //

  // LHE weights (Maximum 9 or size of weights array for QCD variations, 2 for PDF variations and 2 for alphas(mZ) variations)
  LHEOriginalWeight = (*lhe_evt)->originalXWGTUP();
  vector<float> tmpWgtArray;
  vector<float> LHEPDFVariationWgt;
  vector<float> LHEPDFAlphaSMZWgt;
  for (int iw=0; iw<(int)((*lhe_evt)->weights().size()); iw++){
    int wgtid=atoi((*lhe_evt)->weights().at(iw).id.c_str());
    float wgtval=(*lhe_evt)->weights().at(iw).wgt / LHEOriginalWeight;
    //cout << "PDF id = " << PDFid.at(0) << " " << wgtid << " -> " << wgtval << endl;

    if (wgtid<2000) LHEWeight.push_back(wgtval);
    else if (wgtid<3000) tmpWgtArray.push_back(wgtval); // Add PDF replicas and alphas(mZ) variations from the same pdf
  }

  // Handle LO samples 
  if (tmpWgtArray.size()==0) {
    for (int iw=0; iw<(int)((*lhe_evt)->weights().size()); iw++) {
      int wgtid=atoi((*lhe_evt)->weights().at(iw).id.c_str());
      if (wgtid>=10 && wgtid<=110) {
        float wgtval=(*lhe_evt)->weights().at(iw).wgt / LHEOriginalWeight;
        tmpWgtArray.push_back(wgtval);
      }
    }
  }

  if (!tmpWgtArray.empty()){
    if (tmpWgtArray.size()>=100){
      for (unsigned int iwgt=0; iwgt<100; iwgt++) addByLowestInAbs(tmpWgtArray.at(iwgt), LHEPDFVariationWgt); // Since weights could be (-) or (+), use LowestInAbs
      if (tmpWgtArray.size()>100){ for (unsigned int iwgt=100; iwgt<tmpWgtArray.size(); iwgt++) LHEPDFAlphaSMZWgt.push_back(tmpWgtArray.at(iwgt)); } // Dn, Up
    }
  }
  // Find the proper PDF and alphas(mZ) variations
  if (!LHEWeight.empty()){
    float centralWeight=1;
    centralWeight = LHEWeight.at(0);
    LHEWeight_PDFVariationUpDn.push_back(findNearestOneSigma(centralWeight, 1, LHEPDFVariationWgt));
    LHEWeight_PDFVariationUpDn.push_back(findNearestOneSigma(centralWeight, -1, LHEPDFVariationWgt));
    if (LHEPDFAlphaSMZWgt.size()>1){
      float asdn = LHEPDFAlphaSMZWgt.at(0);
      float asup = LHEPDFAlphaSMZWgt.at(1);
      // Rescale alphas(mZ) variations from 0.118+-0.001 to 0.118+-0.0015
      LHEWeight_AsMZUpDn.push_back(centralWeight + (asup-centralWeight)*1.5);
      LHEWeight_AsMZUpDn.push_back(centralWeight + (asdn-centralWeight)*1.5);
    }
  }
  //
}


MELACandidate* LHEHandler::matchAHiggsToParticle(LHE_Event& ev, MELAParticle const* genH){
  MELACandidate* cand=0;
  for (int t=0; t<ev.getNZZCandidates(); t++){
    MELACandidate* tmpCand = ev.getZZCandidate(t);

    double dotproduct = sqrt(genH->p4.Vect().Dot(tmpCand->p4.Vect()) + genH->t()*tmpCand->t());
    double genhdotproduct = sqrt(genH->p4.Vect().Dot(genH->p4.Vect()) + genH->t()*genH->t());
    double massdiff = fabs(genhdotproduct-dotproduct);
    double massratio = 0;
    if (genhdotproduct>0) massratio = massdiff / genhdotproduct;
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

  double HVVmass = PDGHelpers::Zeromass;
  if (isZZ==0) HVVmass = PDGHelpers::Wmass;
  else if (isZZ==1 || isZZ==3) HVVmass = PDGHelpers::Zmass;

  double diffmass1 = fabs(cand1->getSortedV(0)->m()-HVVmass);
  double diffmass2 = fabs(cand2->getSortedV(0)->m()-HVVmass);
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
  return theChosenOne;
}

void LHEHandler::addByLowestInAbs(float val, std::vector<float>& valArray){
  bool inserted = false;
  for (std::vector<float>::iterator it = valArray.begin(); it<valArray.end(); it++){
    if (fabs(*it)>=fabs(val)){
      inserted=true;
      valArray.insert(it, val);
      break;
    }
  }
  if (!inserted) valArray.push_back(val);
}

float LHEHandler::findNearestOneSigma(float ref, int lowhigh, std::vector<float> const& wgt_array){
  int nrep = wgt_array.size();
  int pos_low=-1, pos_high=nrep;
  for (int irep=0; irep<nrep; irep++){ // Assume ordered from low to high
    float tmp = fabs(wgt_array.at(irep));
    if (fabs(tmp)<fabs(ref)) pos_low=irep;
    else if (fabs(tmp)>fabs(ref)){
      pos_high=irep;
      break;
    }
  }
  if (lowhigh==-1 && nrep>0){ // Low
    float ninst = pos_low+1;
    if (ninst==0) return ref;
    int nprogress = (ninst*0.6827+0.5); // truncate down
    int taken = ninst-nprogress;
    return wgt_array.at(taken);
  }
  else if (lowhigh==1 && nrep>0){
    float ninst = nrep-pos_high;
    if (ninst==0) return ref;
    int nprogress = (ninst*0.6827+0.5); // truncate down
    int taken = nprogress+pos_high-1;
    return wgt_array.at(taken);
  }
  else return ref;
}

