/** \class MELABranch
*
*
*  \author N. Amapane - Torino
*  \author U. Sarica - JHU
*/
#ifndef MELAOPTIONPARSER_H
#define MELAOPTIONPARSER_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <cassert>
#include "TString.h"
#include <ZZMatrixElement/MELA/interface/TVar.hh>


class MELAOptionParser{

protected:

  std::vector<std::string> rawOptions;
  std::string couplingsString;

public:

  TVar::Process proc;
  TVar::Production prod;
  TVar::MatrixElement ME;
  TVar::SuperMelaSyst superSyst;
  Bool_t includePAux;
  UInt_t isGenProb;
  Float_t hmass;
  Float_t h2mass;
  Float_t hwidth;
  Float_t h2width;

  SpinZeroCouplings coupl_H;
  SpinOneCouplings coupl_Zp;
  SpinTwoCouplings coupl_X;

  MELAOptionParser(std::string stropts);
  ~MELAOptionParser(){}

  void analyze();
  void splitOption(std::string rawoption, std::string& wish, std::string& value, char delimiter=':');
  void splitOptionRecursive(std::string rawoption, std::vector<std::string>& splitoptions, char delimiter=';');
  void interpretOption(std::string wish, std::string value);

  Bool_t isGen(){ return Bool_t(isGenProb>0); }
  Bool_t hasPAux(){ return includePAux; }

protected:

  Bool_t checkListVariable(std::vector<std::string>& list, std::string var);

  void setProcess(std::string opt);
  void setProduction(std::string opt);
  void setME(std::string opt);
  void setSuperMelaSyst(std::string opt);
  void extractCoupling(std::string opt);

};

#endif

