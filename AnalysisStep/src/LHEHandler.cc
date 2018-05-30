#include <ZZAnalysis/AnalysisStep/interface/LHEHandler.h>
#include <MelaAnalytics/EventContainer/interface/HiggsComparators.h>
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
using namespace HiggsComparators;

const static bool useNNPDF30 = true;

LHEHandler::LHEHandler(edm::Handle<LHEEventProduct>* lhe_evt_, int VVMode_, int VVDecayMode_, bool doKinematics_, int year_) :
VVMode(VVMode_),
VVDecayMode(VVDecayMode_),
doKinematics(doKinematics_),
year(year_),
genEvent(0),
genCand(0)
{
  setHandle(lhe_evt_); // Also calls clear()
  extract();
}
LHEHandler::LHEHandler(int VVMode_, int VVDecayMode_, bool doKinematics_, int year_) :
VVMode(VVMode_),
VVDecayMode(VVDecayMode_),
doKinematics(doKinematics_),
year(year_),
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
  powhegOriginalWeight=1;
  defaultNLOweight=1;
}

MELACandidate* LHEHandler::getBestCandidate(){ return genCand; }
float const& LHEHandler::getLHEOriginalWeight() const{
  return this->LHEOriginalWeight;
}
float const& LHEHandler::getPowhegOriginalWeight() const{
  return this->powhegOriginalWeight;
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
float LHEHandler::reweightNNLOtoNLO() const{
  if (year == 2017)
    return defaultNLOweight; //note this is already divided by originalPowhegWeight
  throw cms::Exception("LHEWeights") << "Shouldn't be calling this function for " << year;
}



void LHEHandler::extract(){
  if (lhe_evt==0) cerr << "LHEHandler::extract: lhe_evt==0" << endl;

  if (lhe_evt!=0){
    if (!lhe_evt->isValid()) cerr << "LHEHandler::extract: lhe_evt invalid!" << endl;

    if (lhe_evt->isValid()){

      readEvent();

      if (doKinematics){

        genEvent = new MELAEvent();
        vectorInt hasGenHiggs;

        {
          int p=0;
          for (MELAParticle* genPart:particleList){
            if (isAHiggs(genPart->id)){
              hasGenHiggs.push_back(p);
              if (VVMode==-1 && (genPart->genStatus==1 || genPart->genStatus==2)) genEvent->addIntermediate(genPart);
            }
            if (genPart->genStatus==1){
              if (isALepton(genPart->id)) genEvent->addLepton(genPart);
              else if (isANeutrino(genPart->id)) genEvent->addNeutrino(genPart);
              else if (isAPhoton(genPart->id)) genEvent->addPhoton(genPart);
              else if (isAGluon(genPart->id) || isAQuark(genPart->id)) genEvent->addJet(genPart);
            }
            p++;
          }
        }

        genEvent->constructVVCandidates(VVMode, VVDecayMode);
        for (MELAParticle* genPart:particleList){ if (genPart->genStatus==-1) genEvent->addVVCandidateMother(genPart); }
        genEvent->addVVCandidateAppendages();

        genCand=nullptr;
        if (!hasGenHiggs.empty()){
          for (int iH:hasGenHiggs){
            MELACandidate* tmpCand = matchAHiggsToParticle(*genEvent, particleList.at(iH));
            if (tmpCand){
              if (!genCand) genCand=tmpCand;
              else genCand = candComparator(genCand, tmpCand, HiggsComparators::BestZ1ThenZ2ScSumPt, VVMode);
            }
          }
        }
        if (!genCand) genCand = candidateSelector(*genEvent, HiggsComparators::BestZ1ThenZ2ScSumPt, VVMode);

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
  vector<float> LHEPDFAlphaSMZWgt;
  vector<float> LHEPDFVariationWgt;
  bool founddefaultNLOweight = false;
  bool foundpowhegOriginalWeight = false;
  AlternateWeightsType weightstype = unknown;
  powhegOriginalWeight = LHEOriginalWeight;
  //first find the main powheg weight (should be one of the first)
  for (const auto& weight : (*lhe_evt)->weights()) {
    if (weight.id == "1001") {
      powhegOriginalWeight = weight.wgt;
      foundpowhegOriginalWeight = true;
      break;
    }
  }

  for (const auto& weight : (*lhe_evt)->weights()) {
    int wgtid;
    try {
      wgtid = stoi(weight.id.c_str());
    } catch (std::invalid_argument& e) {
      continue;  //we don't use non-numerical indices, but they exist in some 2016 MC samples
    }
    float wgtval=weight.wgt / powhegOriginalWeight;
    //cout << "PDF id = " << PDFid.at(0) << " " << wgtid << " -> " << wgtval << endl;
    if (year == 2016) {
      if (wgtid<2000) LHEWeight.push_back(wgtval);
      else if (wgtid<3000) LHEPDFVariationWgt.push_back(wgtval); // Add PDF replicas and alphas(mZ) variations from the same pdf
    } else if (year == 2017) {
      //Madgraph 0 offset
      if (weightstype == unknown && wgtid == 1) {weightstype = madgraph_0offset; LHEWeight.push_back(wgtval);}
      else if (weightstype == madgraph_0offset && 2 <= wgtid && wgtid <= 9) LHEWeight.push_back(wgtval);
      else if (weightstype == madgraph_0offset && 10 <= wgtid && wgtid <= 120) {/*do nothing, these are the NNLO variations*/}
      else if (weightstype == madgraph_0offset && wgtid == 121) {founddefaultNLOweight = true; defaultNLOweight = wgtval;}
      else if (weightstype == madgraph_0offset && 122 <= wgtid && wgtid <= 223) LHEPDFVariationWgt.push_back(wgtval);
      else if (weightstype == madgraph_0offset && 224 <= wgtid && wgtid <= 1080) {/*do nothing, these are other various weights*/}

      //QCD variations for all the other weightstypes
      else if (weightstype == unknown && 1001 <= wgtid && wgtid <= 1009) LHEWeight.push_back(wgtval);

      //Madgraph 1000 offset
      else if (wgtid == 1010 && weightstype == unknown) weightstype = madgraph_1000offset;
      else if (weightstype == madgraph_1000offset && 1011 <= wgtid && wgtid <= 1120) {/*do nothing, these are the NNLO variations*/}
      else if (weightstype == madgraph_1000offset && wgtid == 1121) {founddefaultNLOweight = true; defaultNLOweight = wgtval;}
      else if (weightstype == madgraph_1000offset && 1122 <= wgtid && wgtid <= 1223) LHEPDFVariationWgt.push_back(wgtval);
      else if (weightstype == madgraph_1000offset && 1224 <= wgtid && wgtid <= 2080) {/*do nothing, these are other various weights*/}

      //powheg
      else if (weightstype == unknown && wgtid == 2000) {weightstype = powheg; /*but do nothing, this is an NNLO variation*/}
      else if (weightstype == powheg && 1500 <= wgtid && wgtid <= 1602 && useNNPDF30) {
        if (wgtid == 1500) {founddefaultNLOweight = true; defaultNLOweight = wgtval;}
        else LHEPDFVariationWgt.push_back(wgtval);
      }
      else if (weightstype == powheg && 1500 <= wgtid && wgtid <= 1602) {/*do nothing, these are the NLO pdf for NNPDF30 and variations*/}
      else if (weightstype == powheg && wgtid == 1700)                  {/*do nothing, this is the NNLO pdf for NNPDF30*/}
      else if (weightstype == powheg && (wgtid == 1800 || wgtid == 1850 || wgtid == 1900 || wgtid == 1950)) {/*do nothing, these are LO pdfs*/}
      else if (weightstype == powheg && 2001 <= wgtid && wgtid <= 2111) {/*do nothing, these are more NNLO variations*/}
      else if (weightstype == powheg && wgtid >= 4000)                  {/*do nothing, these are other various weights*/}
      else if (weightstype == powheg && 3000 <= wgtid && wgtid <= 3102 && useNNPDF30) {/*do nothing, these are the NLO pdf for NNPDF31 and variations*/}
      else if (weightstype == powheg && wgtid == 3000) {founddefaultNLOweight = true; defaultNLOweight = wgtval;}
      else if (weightstype == powheg && 3001 <= wgtid && wgtid <= 3102) LHEPDFVariationWgt.push_back(wgtval);

      else throw cms::Exception("LHEWeights") << "Don't know what to do with alternate weight id = " << wgtid << "(weightstype == " << weightstype << ")";
    } else {
      throw cms::Exception("LHEWeights") << "Unknown year " << year;
    }
  }

  if (year == 2017 && !(*lhe_evt)->weights().empty() && !(foundpowhegOriginalWeight && founddefaultNLOweight && LHEWeight.size() == 9 && LHEPDFVariationWgt.size() == 102)) {
    throw cms::Exception("LHEWeights")
            << "For 2017 MC, expect to find either\n"
            << " - no alternate LHE weights, or\n"
            << " - all of the following:\n"
            << "   - muR and muF variations (1001-1009, found " << LHEWeight.size() << " of them, " << (foundpowhegOriginalWeight ? "" : "not ") << "including 1001)\n"
            << "   - the default NLO PDF weight (3000, " << (founddefaultNLOweight ? "found" : "didn't find") << " it)\n"
            << "   - NLO PDF weight variations (3001-3102, found " << LHEPDFVariationWgt.size() << " of them)";
  }

  if (year == 2017 && weightstype == madgraph_1000offset) { //but not 0offset!  0offset does it the same way as powheg
    LHEWeight = {LHEWeight[0], LHEWeight[3], LHEWeight[6],  //note LHEWeight[8] is always defined here,
                 LHEWeight[1], LHEWeight[4], LHEWeight[7],  //because weightstype != unknown implies that there are weights present
                 LHEWeight[2], LHEWeight[5], LHEWeight[8]}; //so if LHEweight.size() != 9 you already got an exception
  }                                                         //in the previous block

  // Handle LO samples 
  if (year == 2016 && LHEPDFVariationWgt.size()==0) {
    for (int iw=0; iw<(int)((*lhe_evt)->weights().size()); iw++) {
      int wgtid=atoi((*lhe_evt)->weights().at(iw).id.c_str());
      if (wgtid>=10 && wgtid<=110) {
        float wgtval=(*lhe_evt)->weights().at(iw).wgt / LHEOriginalWeight;
        LHEPDFVariationWgt.push_back(wgtval);
      }
    }
  }

  if (LHEPDFVariationWgt.size() > 100) {
    auto firstalphasweight = LHEPDFVariationWgt.begin() + 100; //= iterator to LHEPDFVariationWgt[100] = 101st entry
    LHEPDFAlphaSMZWgt.assign(firstalphasweight, LHEPDFVariationWgt.end());
    LHEPDFVariationWgt.erase(firstalphasweight, LHEPDFVariationWgt.end());
  }


  // Find the proper PDF and alphas(mZ) variations
  if (!LHEWeight.empty()){
    switch (year) {
      case 2016: {
        std::sort(LHEPDFVariationWgt.begin(), LHEPDFVariationWgt.end(), compareAbs);
        float centralWeight = LHEWeight.at(0);
        LHEWeight_PDFVariationUpDn.push_back(findNearestOneSigma(centralWeight, 1, LHEPDFVariationWgt));
        LHEWeight_PDFVariationUpDn.push_back(findNearestOneSigma(centralWeight, -1, LHEPDFVariationWgt));
        if (LHEPDFAlphaSMZWgt.size()>1){
          float asdn = LHEPDFAlphaSMZWgt.at(0);
          float asup = LHEPDFAlphaSMZWgt.at(1);
          // Rescale alphas(mZ) variations from 0.118+-0.001 to 0.118+-0.0015
          LHEWeight_AsMZUpDn.push_back(centralWeight + (asup-centralWeight)*1.5);
          LHEWeight_AsMZUpDn.push_back(centralWeight + (asdn-centralWeight)*1.5);
        }
        break;
      }
      case 2017: {
        //https://arxiv.org/pdf/1706.00428v2.pdf page 88
        //we are working with NNPDF31_nlo_hessian_pdfas
        //see the routine starting on line 101 of PDFSet.cc in LHAPDF-6.2.1
        //based on NNPDF31_nlo_hessian_pdfas.info, which confirms that errorType() is symmhessian+as
        float centralWeight = defaultNLOweight;
        if (centralWeight == 0) {
          LHEWeight_PDFVariationUpDn = {0, 0};
          LHEWeight_AsMZUpDn = {0, 0};
          edm::LogWarning warning("ZeroWeight");
          warning << "default NLO PDF weight is 0\nIncoming particle id and pz:";
          for (const auto& p : particleList)
            if (p->genStatus == -1)
              warning << "\n" << p->id << " " << p->z();
          break;
        }

        float errorsquared = 0;
        for (const auto& wt : LHEPDFVariationWgt) {
          float difference = wt - centralWeight;
          errorsquared += difference*difference;
        }
        float error = sqrt(errorsquared);
        LHEWeight_PDFVariationUpDn = {(centralWeight + error) / reweightNNLOtoNLO(), (centralWeight - error) / reweightNNLOtoNLO()};

        if (LHEPDFAlphaSMZWgt.size()>1){
          float asdn = LHEPDFAlphaSMZWgt.at(0);
          float asup = LHEPDFAlphaSMZWgt.at(1);
          // Rescale alphas(mZ) variations from 0.118+-0.002 to 0.118+-0.0015
          //                           Note this number ^ is different than 2016!
          LHEWeight_AsMZUpDn.push_back((centralWeight + (asup-centralWeight)*0.75) / reweightNNLOtoNLO());
          LHEWeight_AsMZUpDn.push_back((centralWeight + (asdn-centralWeight)*0.75) / reweightNNLOtoNLO());
        }
        break;
      }
      default: {
        throw cms::Exception("LHEWeights") << "Unknown year " << year;
      }
    }
  }
  //
}


bool LHEHandler::compareAbs(float val1, float val2) {
  return std::abs(val1) < std::abs(val2);
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

