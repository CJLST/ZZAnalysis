/** \class MELAComputation
*
*
*  \author N. Amapane - Torino
*  \author U. Sarica - JHU
*
*
*  Description:
*
*  This class is just a formula of different MELA hypotheses.
*  Such a formula implementation is needed in order to extract interference contributions to the probabilities,
*  or in case there are additional divisors or multiplier to the value that goes into the branch.
*  This class also decides on which ME instance to pick, in case a maximization of a ME ratio is requested
*  (eg. BestDVBFHJJ or best SM JHUGen Lep_WH).
*
*/
#ifndef MELACOMPUTATION_H
#define MELACOMPUTATION_H

#include <ZZAnalysis/AnalysisStep/interface/MELAHypothesis.h>


class MELAComputation{

protected:

  MELAHypothesis* targetP;
  MELAOptionParser* opt;
  Bool_t contUpdate;
  Float_t pME;
  Float_t pAux;
  Float_t cMEAvg;
  Float_t maximizationCachedVal;

  std::vector<MELAHypothesis*> addedP;
  std::vector<MELAHypothesis*> subtractedP;
  std::vector<MELAHypothesis*> multipliedP;
  std::vector<MELAHypothesis*> dividedP;
  // The calculation done is (targetP-subtractedP)*multipliedP/dividedP

  std::vector<MELAHypothesis*> maximize_num;
  std::vector<MELAHypothesis*> maximize_denom;
  // If the options specify maximization of num/denom, keep track of the value and only update when needed.


  void addContingency(std::vector<MELAHypothesis*>& allHypos, std::vector<std::string>& source, std::vector<MELAHypothesis*>& dest);
  Bool_t testMaximizationCache(); // Used in update()
  Float_t extractVal(MELAHypothesis::METype valtype); // Computation after all contingencies are added

public:

  MELAOptionParser* getOption(){ return opt; }

  void addContingencies(std::vector<MELAHypothesis*>& allHypos); // Fill addedP, subtractedP, dividedP, multipliedP etc. vectors based on the specifications in the options

  void update();
  void forceUpdate();
  Float_t getVal(MELAHypothesis::METype valtype) const; // Extract pME, pAux or cMEAvg
  MELAHypothesis* getHypothesis(){ return targetP; } // Return the hypothesis, most likely to group in a vector
  std::string getName() const{ if (opt!=0) return opt->getName(); else return ""; }
  std::string getAlias() const{ if (opt!=0) return opt->getAlias(); else return ""; }
  std::string getCluster() const{ if (opt!=0) return opt->getCluster(); else return ""; }

  MELAComputation(MELAHypothesis* targetP_);
  virtual ~MELAComputation();
  
  void reset();

};


#endif
