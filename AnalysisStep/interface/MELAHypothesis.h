/** \class MELAHypothesis
*
*
*  \author N. Amapane - Torino
*  \author U. Sarica - JHU
*
*
*  Description:
*
*  This class is the main class to compute any MELA probability.
*  It is meant to be completely blind to anything else might be happening outside itself.
*
*/
#ifndef MELAHYPOTHESIS_H
#define MELAHYPOTHESIS_H

#include <ZZMatrixElement/MELA/interface/Mela.h>
#include "MELAOptionParser.h"


class MELAHypothesis{

protected:

  Mela* mela;
  MELAOptionParser* opt;
  Bool_t optIsOwned;

  float pME;
  float pAux;
  float cMEAvg;

public:

  enum METype{
    UseME,
    UsePAux,
    UsePConstant
  };

  Float_t getVal(METype valtype) const;
  MELAOptionParser* getOption(){ return opt; }

  MELAHypothesis(
    Mela* mela_,
    MELAOptionParser* opt_
    );
  MELAHypothesis(
    Mela* mela_,
    std::string stropt
    );
  virtual ~MELAHypothesis(){ if (optIsOwned) delete opt; opt=0; }

  void computeP(MELACandidate* cand); // Wrapper
  void computeP(unsigned int index); // Wrapper
  void computeP(); // Main function
  void computePM4l(MELACandidate* cand); // Wrapper
  void computePM4l(unsigned int index); // Wrapper
  void computePM4l(); // Main function
  void reset();

};


#endif
