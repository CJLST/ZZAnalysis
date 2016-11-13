#include <ZZAnalysis/AnalysisStep/interface/MELAComputation.h>

using namespace std;


MELAComputation::MELAComputation(MELAHypothesis* targetP_) :
  targetP(targetP_),
  pME(0.),
  pAux(1.),
  cMEAvg(1.)
{
  setOption(targetP->getOption());
  resetMaximizationCache();
}
MELAComputation::~MELAComputation(){
  addededP.clear();
  subtractedP.clear();
  multipliedP.clear();
  dividedP.clear();
  maximize_num.clear();
  maximize_denom.clear();
}

void MELAComputation::addContingencies(vector<MELAHypothesis*>& allHypos){
  addContingency(allHypos, opt->addedAliases, addedP);
  addContingency(allHypos, opt->subtractedAliases, subtractedP);
  addContingency(allHypos, opt->multipliedAliases, multipliedP);
  addContingency(allHypos, opt->dividedAliases, dividedP);

  addContingency(allHypos, opt->maximizationNumerators, maximize_num);
  addContingency(allHypos, opt->maximizationDenominators, maximize_denom);
}
void MELAComputation::addContingency(vector<MELAHypothesis*>& allHypos, vector<string>& source, vector<MELAHypothesis*>& dest){
  for (unsigned int is=0; is<source.size(); is++){
    // Match by alias, not by name!
    string aliasToMatch = source.at(is);

    for (unsigned int ih=0; ih<allHypos.size(); ih++){
      if (allHypos.at(ih)->getOption()->isAliased() && aliasToMatch==allHypos.at(ih)->getOption()->getAlias()){
        dest.push_back(allHypos.at(ih));
        break;
      }
    }

  }
}

Bool_t MELAComputation::testMaximizationCache(){
  Float_t testCache = ((maximize_num.size()+maximize_denom.size())>0 ? 1. : -1.);
  // Always use type UseME for such comparisons, the others don't make much sense
  for (unsigned int ip=0; ip<maximize_num.size(); ip++){ testCache *= maximize_num.at(ip)->getVal(MELAHypothesis::UseME); }
  for (unsigned int ip=0; ip<maximize_denom.size(); ip++){ testCache /= maximize_denom.at(ip)->getVal(MELAHypothesis::UseME); }
  if (testCache>=maximizationCachedVal){
    maximizationCachedVal = testCache;
    return true;
  }
  else return false;
}
void MELAComputation::update(){
  if (contUpdate && testMaximizationCache()){
    pME = extractVal(MELAHypothesis::UseME);
    pAux = extractVal(MELAHypothesis::UsePAux);
    cMEAvg = extractVal(MELAHypothesis::UsePConstant);
    contUpdate = false;
  }
}
void MELAComputation::forceUpdate(){
  if (testMaximizationCache()){
    pME = extractVal(MELAHypothesis::UseME);
    pAux = extractVal(MELAHypothesis::UsePAux);
    cMEAvg = extractVal(MELAHypothesis::UsePConstant);
  }
}

Float_t MELAComputation::extractVal(MELAHypothesis::METype valtype){
  Float_t tmp = targetP->getVal(valtype);
  for (unsigned int ime=0; ime<addedP.size(); ime++){ if (addedP.at(ime)!=0) tmp += addedP.at(ime)->getVal(valtype); }
  for (unsigned int ime=0; ime<subtractedP.size(); ime++){ if (subtractedP.at(ime)!=0) tmp -= subtractedP.at(ime)->getVal(valtype); }

  Float_t factor=1;
  for (unsigned int ime=0; ime<multipliedP.size(); ime++){ if (subtractedP.at(ime)!=0) factor *= multipliedP.at(ime)->getVal(valtype); }
  for (unsigned int ime=0; ime<dividedP.size(); ime++){ if (dividedP.at(ime)!=0) factor /= dividedP.at(ime)->getVal(valtype); }

  if (factor==factor) tmp *= factor;
  else tmp=0;

  return tmp;
}
Float_t MELAComputation::getVal(MELAHypothesis::METype valtype) const{
  switch (valtype){
  case MELAHypothesis::UseME:
    return pME;
  case MELAHypothesis::UsePAux:
    return pAux;
  case MELAHypothesis::UsePConstant:
    return cMEAvg;
  default:
    return 0;
  }
}

void MELAComputation::resetMaximizationCache(){
  contUpdate=true;
  maximizationCachedVal=-1.;
}
void MELAComputation::reset(){
  resetMaximizationCache();
  targetP->reset();
  pME = targetP->getVal(MELAHypothesis::UseME);
  pAux = targetP->getVal(MELAHypothesis::UsePAux);
  cMEAvg = targetP->getVal(MELAHypothesis::UsePConstant);
}



