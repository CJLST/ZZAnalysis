#ifndef ZZ4lConfigHelper_h
#define ZZ4lConfigHelper_h

/** \class ZZ4lConfigHelper
 *
 *  No description available.
 *
 *  $Date: 2013/01/30 21:54:21 $
 *  $Revision: 1.7 $
 *  \author N. Amapane - CERN
 */

#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <DataFormats/Common/interface/TriggerResults.h>
#include <FWCore/Common/interface/TriggerNames.h>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>

class ZZ4lConfigHelper {
 public:
  /// Constructor
  ZZ4lConfigHelper(const edm::ParameterSet& pset);

  /// Destructor
  virtual ~ZZ4lConfigHelper(){}

  /// Pass skim
  bool passSkim(const edm::Event & event)  { short bw=0; return passSkim(event,bw); }
  
  /// Pass skim (set bit in trigworld)
  bool passSkim(const edm::Event & event, short& trigworld);
  
  /// Pass trigger requests
  bool passTrigger(const edm::Event & event) { short bw=0; return passTrigger(event,bw); }

  /// Pass trigger requests (and set bits in trigworld)
  bool passTrigger(const edm::Event & event, short& trigworld);


  /// Pass MC filters specified in the card "MCFilterPath"
  bool passMCFilter(const edm::Event & event);
  
  /// Pass the specified filter
  bool passFilter(const edm::Event & event, const std::string& filterPath, bool fromHLT=false);

  bool isMC() {return isMC_;};

  Channel channel() {return theChannel;}

  /// Running condition to be emulated (2011 or 2012)
  int setup() {return theSetup;};

  /// Type of MC sample (2011, 2012); can be different from setup()!
  int sampleType() {return theSampleType;}

  std::string PD;

 private:
  Channel theChannel;
  bool isMC_;
  int theSetup;
  int theSampleType;
  std::vector<std::string> skimPaths;
  std::string MCFilter;
  edm::EventID cachedEvtId;  
  edm::Handle<edm::TriggerResults> triggerResults;
  edm::Handle<edm::TriggerResults> triggerResultsHLT;
  const edm::TriggerNames* triggerNames;
  const edm::TriggerNames* triggerNamesHLT;

  void eventInit(const edm::Event & event);

};
#endif

