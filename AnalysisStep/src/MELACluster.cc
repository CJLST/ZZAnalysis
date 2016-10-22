#include <ZZAnalysis/AnalysisStep/interface/MELACluster.h>

using namespace std;


MELACluster::MELACluster(string name_) : name(name_){}
MELACluster::~MELACluster(){}
void MELACluster::addComputation(MELAComputation* comp){ computers.push_back(comp); }
void MELACluster::update(){ for (unsigned int ic=0; ic<computers.size(); ic++) computers.at(ic)->update(); }
void MELACluster::reset(){ for (unsigned int ic=0; ic<computers.size(); ic++) computers.at(ic)->reset(); }
