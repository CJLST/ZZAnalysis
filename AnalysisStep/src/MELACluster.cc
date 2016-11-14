#include <ZZAnalysis/AnalysisStep/interface/MELACluster.h>

using namespace std;


MELACluster::MELACluster(string name_) : name(name_){}
MELACluster::~MELACluster(){}
void MELACluster::addComputation(MELAComputation* comp){ computers.push_back(comp); }
void MELACluster::computeAll(){
  // Re-compute all related hypotheses
  for (unsigned int ic=0; ic<computers.size(); ic++){
    // Avoid re-computing MEs twice (could happen through copy-computations)
    if (!computers.at(ic)->getOption()->isCopy()) computers.at(ic)->getHypothesis()->computeP();
  }
}
void MELACluster::update(){ for (unsigned int ic=0; ic<computers.size(); ic++) computers.at(ic)->update(); }
void MELACluster::forceUpdate(){ for (unsigned int ic=0; ic<computers.size(); ic++) computers.at(ic)->forceUpdate(); }
void MELACluster::reset(){ for (unsigned int ic=0; ic<computers.size(); ic++) computers.at(ic)->reset(); }

