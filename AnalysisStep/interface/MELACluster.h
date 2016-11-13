/** \class MELACluster
*
*
*  \author N. Amapane - Torino
*  \author U. Sarica - JHU
*
*
*  Description:
*
*  This class contains the matrix element computations
*  that belong to a specific choice of candidates.
*  
*
*/
#ifndef MELACLUSTER_H
#define MELACLUSTER_H

#include <ZZAnalysis/AnalysisStep/interface/MELAComputation.h>


class MELACluster{

protected:

  std::string name;
  std::vector<MELAComputation*> computers;

public:

  MELACluster(std::string name_);
  virtual ~MELACluster();

  std::string getName(){ return name; }
  void addComputation(MELAComputation* comp);
  void update();
  void forceUpdate();
  void reset();

};


#endif
