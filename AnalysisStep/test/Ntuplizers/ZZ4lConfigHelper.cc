#include "ZZ4lConfigHelper.h"
#include "ZZ4lConfigHelper.h"
#include <ZZAnalysis/AnalysisStep/interface/bitops.h>

using namespace std;
using namespace edm;

ZZ4lConfigHelper::ZZ4lConfigHelper(const ParameterSet& pset) :
  PD(pset.getParameter<std::string>("PD")),
  isMC_(pset.getUntrackedParameter<bool>("isMC")),
  theSetup(pset.getParameter<int>("setup")),
  theSampleType(pset.getParameter<int>("sampleType")),
  skimPaths(pset.getParameter<std::vector<std::string> >("skimPaths")),
  MCFilter(pset.getParameter<std::string>("MCFilterPath")),
  anyTrigger(true) // take the OR of all trigger paths
  
{
  string channel = pset.getUntrackedParameter<string>("channel");
  theChannel = finalState(channel);
  
  // Check for inconsistent configurations
  if ( ( theSampleType!=2011 && theSampleType!=2012 && theSampleType!=2015 && theSampleType!=2016 && theSampleType!=2017 && theSampleType!=2018) ||
       ( theSetup!=2011 && theSetup!=2012 && theSetup!=2015 && theSetup!=2016 && theSetup!=2017 && theSetup!=2018) ||
       ( theSampleType!=theSetup ) // No sample rescaling supported as of now.
       // We may add exception for MC only when needed.
       ) {
    cout << "ERROR: ZZ4lConfigHelper: inconsistent setup: sampleType=" << theSampleType << ", setup=" << theSetup << ", isMC=" <<isMC_ << endl;
    abort();
  }
  
  
  if ((isMC_&&PD!="") || (!isMC_ && (PD!="DoubleEle" && PD!="DoubleMu" && PD!="MuEG" && PD!="DoubleEG" && PD!="DoubleMuon" && PD!="MuonEG" && PD!="SingleElectron" && PD!="SingleMuon" && PD!="EGamma"))) {
    cout << "ERROR: ZZ4lConfigHelper: isMC: " << isMC_ << " PD: " << PD << endl;
    abort();
  }    

  if (!isMC_&&MCFilter!="") {
    cout << "ERROR: ZZ4lConfigHelper: MCFilter= " << MCFilter << " when isMC=0" 
	 << endl;
    abort();
  }    
  
}

bool 
ZZ4lConfigHelper::passMCFilter(const edm::Event & event, edm::Handle<edm::TriggerResults> & trigRes){
  if (MCFilter=="") return true;
  return passFilter(event, trigRes, MCFilter);
}

bool 
ZZ4lConfigHelper::passSkim(const edm::Event & event, edm::Handle<edm::TriggerResults> & trigRes, short& trigworld){
  bool evtPassSkim = false;
  if (skimPaths.size()==0) {
    evtPassSkim=true;
  } else {
    for (vector<string>::const_iterator name = skimPaths.begin(); name!= skimPaths.end(); ++name) {
      if (passFilter(event, trigRes, *name)) {
	evtPassSkim = true; 
	break;
      }
    }
  }
  if (evtPassSkim) set_bit_16(trigworld,15);
  return evtPassSkim;
}

bool 
ZZ4lConfigHelper::passTrigger(const edm::Event & event, edm::Handle<edm::TriggerResults> & trigRes, short& trigworld){

  bool passDiMu  = passFilter(event, trigRes, "triggerDiMu");
  bool passDiEle = passFilter(event, trigRes, "triggerDiEle");
  bool passMuEle  = passFilter(event, trigRes, "triggerMuEle"); 
  if ((!isMC_) && theSetup == 2011 ) { // follow changes in trigger menu in data 2011 (see wiki)
    int irun=event.id().run(); 
    if (irun>=175973) {
      passMuEle  = passFilter(event, trigRes, "triggerMuEle3");
    } else if (irun>=167914) {
      passMuEle  = passFilter(event, trigRes, "triggerMuEle2");
    }
  }
  bool passTriEle = false;
  if (theSetup == 2012 || theSetup >= 2015) {
    passTriEle = passFilter(event, trigRes, "triggerTriEle");
  }
  bool passTriMu = false;
  bool passSingleEle = false;
  if (theSetup >= 2015) {
    passTriMu = passFilter(event, trigRes, "triggerTriMu");
    passSingleEle = passFilter(event, trigRes, "triggerSingleEle");
  }
  bool passSingleMu = false;
  if (theSetup >= 2016) {
    passSingleMu = passFilter(event, trigRes, "triggerSingleMu");
  }


  bool evtPassTrigger = false;

  // Check all triggers together if anyTrigger is specified (or for CRs)
  if (anyTrigger || theChannel==ZLL || theChannel==ZL || theChannel==ZZ) {
      if ((PD=="" && (passDiEle || passDiMu || passMuEle || passTriEle || passTriMu || passSingleEle || passSingleMu)) ||
	  ((PD=="DoubleEle"||PD=="DoubleEG" ||PD=="EGamma" ) && (passDiEle || passTriEle)) ||
	  ((PD=="DoubleMu" ||PD=="DoubleMuon") && (passDiMu || passTriMu) && !passDiEle && !passTriEle) ||
	  ((PD=="MuEG"     ||PD=="MuonEG"    ) && passMuEle && !passDiMu && !passTriMu && !passDiEle && !passTriEle) ||
	  ((PD=="SingleElectron" || PD=="EGamma") && passSingleEle && !passMuEle && !passDiMu && !passTriMu && !passDiEle && !passTriEle) ||
	  (PD=="SingleMuon" && passSingleMu && !passSingleEle && !passMuEle && !passDiMu && !passTriMu && !passDiEle && !passTriEle)
	  ) {
	evtPassTrigger = true;
      } 
  }
  
  else {
    cout << "[ERROR]: ZZ4lConfigHelper: unexpected config " << theChannel << endl;
    abort();
  }

  
  if (evtPassTrigger) set_bit_16(trigworld,0);
  if (passDiMu) set_bit_16(trigworld,1);
  if (passDiEle) set_bit_16(trigworld,2);
  if (passMuEle) set_bit_16(trigworld,3);
  if (passTriEle) set_bit_16(trigworld,4);
  if (passTriMu) set_bit_16(trigworld,5);
  if (passSingleEle) set_bit_16(trigworld,6);
  if (passSingleMu) set_bit_16(trigworld,7);
  

  return evtPassTrigger;
}

bool
ZZ4lConfigHelper::passMETTrigger(const edm::Event & event, edm::Handle<edm::TriggerResults> & trigRes){

  bool passMETFilter  = true;
	
//  if (theSetup >= 2017) {
//    passMETFilter = passFilter(event, trigRes, "triggerMETFilters");
//  }
  return passMETFilter;
}


bool
ZZ4lConfigHelper::passFilter(const edm::Event & event, edm::Handle<edm::TriggerResults> & trigRes, const string& filterPath) {

  //if (event.id()==cachedEvtId) return;
  //cachedEvtId = event.id();

  // Initialize trigger results table
  if (trigRes.isValid()) {
    triggerResults = trigRes;
    triggerNames = &(event.triggerNames(*triggerResults));
  } else {
    cout << "ERROR: failed to get TriggerResults" << endl;
  }

  //  for (unsigned i=0; i<triggerNames->size(); i++) cout << triggerNames->triggerName(i) << endl;
  unsigned i =  triggerNames->triggerIndex(filterPath);
  
  if (i== triggerNames->size()){
    cout << "ERROR: ZZ4lConfigHelper::isTriggerBit: path does not exist! " << filterPath << endl;
    abort();
  }
  //  cout << " Trigger result for " << filterPath << " : accept=" << triggerResults->accept(i) << endl;
  return triggerResults->accept(i);

}




