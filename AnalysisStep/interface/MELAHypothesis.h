/** \class MELABranch
*
*
*  \author N. Amapane - Torino
*  \author U. Sarica - JHU
*/
#ifndef MELAHYPOTHESIS_H
#define MELAHYPOTHESIS_H

#include <ZZMatrixElement/MELA/interface/Mela.h>
#include <ZZAnalysis/AnalysisStep/interface/MELAOptionParser.h>


class MELAHypothesis{

protected:

  Mela* mela;
  MELAOptionParser* opt;

  float pME;
  float pAux;
  float cMEAvg;

  void reset();
  void setVal();

public:

  enum METype{
    UseME,
    UsePAux,
    UsePConstant
  };

  float getVal(METype valtype);

  MELAHypothesis(
    Mela* mela_,
    MELAOptionParser* opt_
    );
  virtual ~MELAHypothesis(){}

  void computeP(MELACandidate* cand); // Wrapper
  void computeP(unsigned int index); // Wrapper
  void computeP(); // Main function
  void computePM4l(MELACandidate* cand); // Wrapper
  void computePM4l(unsigned int index); // Wrapper
  void computePM4l(); // Main function

};


#endif
