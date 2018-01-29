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
#include <MelaAnalytics/EventContainer/interface/MELAEvent.h>

#include <vector>


class LHEHandler{
public:

  LHEHandler(int VVMode_, int VVDecayMode_, bool doKinematics_, int year_);
  LHEHandler(edm::Handle<LHEEventProduct>* lhe_evt_, int VVMode_, int VVDecayMode_, bool doKinematics_, int year_);
  virtual ~LHEHandler();
  
  void setHandle(edm::Handle<LHEEventProduct>* lhe_evt_);
  void extract();
  void clear();

  MELACandidate* getBestCandidate();
  float const& getLHEOriginalWeight() const; // Weight written in the <event> block, supposed to = genhepmcweight if no Pythia reweighting is done
  float getLHEWeight(unsigned int whichWeight, float defaultValue=1) const; // = {Weights written in LHE weight variations} / getLHEOriginalWeight()
  float getLHEWeight_PDFVariationUpDn(int whichUpDn, float defaultValue=1) const; // = {Weights written in LHE weight variations} / getLHEOriginalWeight()
  float getLHEWeigh_AsMZUpDn(int whichUpDn, float defaultValue=1) const; // = {Weights written in LHE weight variations} / getLHEOriginalWeight()
  float const& getPDFScale() const;
  float reweightNNLOtoNLO() const;

  static bool compareAbs(float val1, float val2);
  static float findNearestOneSigma(float ref, int lowhigh, std::vector<float> const& wgt_array);

protected:

  // VVMode and VVDecayMode: See comment lines within MELAEvent::constructVVCandidates
  const int VVMode;
  const int VVDecayMode;
  const bool doKinematics;
  const int year;

  edm::Handle<LHEEventProduct>* lhe_evt;
  vector<MELAParticle*> particleList;
  MELAEvent* genEvent;
  MELACandidate* genCand;

  int defaultNLOweight;
  float LHEOriginalWeight;
  vector<float> LHEWeight;
  vector<float> LHEWeight_PDFVariationUpDn;
  vector<float> LHEWeight_AsMZUpDn;
  vector<int> PDFid;
  float PDFScale;

  void readEvent();

};


#endif
