/** \class LHEHandler
 *
 *
 *  \author N. Amapane - Torino
 *  \author U. Sarica - JHU
 */
#ifndef LHEHADNLER_H
#define LHEHANDLER_H

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <FWCore/ParameterSet/interface/FileInPath.h>
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include <DataFormats/GeometryVector/interface/Point3DBase.h>

#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include <ZZMatrixElement/MELA/interface/Mela.h>
#include <ZZAnalysis/AnalysisStep/interface/LHE_Event.h>

#include <iomanip>
#include <iostream>
#include <string>
#include "TLorentzVector.h"
#include "TString.h"


class LHEHandler{
public:

  LHEHandler(int VVMode_, int VVDecayMode_, bool doKinematics_);
  LHEHandler(edm::Handle<LHEEventProduct>* lhe_evt_, int VVMode_, int VVDecayMode_, bool doKinematics_);
  virtual ~LHEHandler();
  
  void setHandle(edm::Handle<LHEEventProduct>* lhe_evt_);
  void extract();
  void clear();

  MELACandidate* getBestCandidate();
  float getLHEWeight(unsigned int whichWeight, float defaultValue=1);
  float getPDFScale();

protected:

  // VVMode and VVDecayMode: See comment lines within LHE_Event::constructVVCandidates
  int VVMode;
  int VVDecayMode;
  bool doKinematics;

  edm::Handle<LHEEventProduct>* lhe_evt;
  vector<MELAParticle*> particleList;
  LHE_Event* genEvent;
  MELACandidate* genCand;
  vector<float> LHEWeight;
  float PDFScale;

  void readEvent();
  //
  MELACandidate* matchAHiggsToParticle(LHE_Event& ev, MELAParticle* genH);
  MELACandidate* candidateSelector(LHE_Event& ev, int isZZ);
  MELACandidate* candComparator(MELACandidate* cand1, MELACandidate* cand2, int isZZ);

};


#endif
