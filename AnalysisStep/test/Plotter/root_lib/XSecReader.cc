
/*
 *  See header file for a description of this class.
 *
 *  $Date: 2013/06/28 16:19:33 $
 *  $Revision: 1.7 $
 *  \author G. Cerminara - NEU Boston & INFN Torino
 */

#include "XSecReader.h"
#include "TString.h"
#include "TRegexp.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

namespace {
  bool dbg = false;
};



using namespace std;

XSecReader::XSecReader(const TString& xsecFileName, const TString& lumiFileName) {
  ifstream xsecFile(xsecFileName.Data());
  string line;
  // Read the xsection file and store numbers in a map (per sample)

  TRegexp empty1("^ *#");
  TRegexp empty2("^ *$");

  int iline = 0;
  while (getline(xsecFile,line)) {
    ++iline;
    if( line == "" || line[0] == '#' ) continue; // Skip comments and empty lines
    TString tline(line);
    if (tline.Contains(empty1) || tline.Contains(empty2)) continue;
    stringstream linestr;
    linestr << line;
    TString finalState, sampleName;
    double nEvents;
    double xsec, br;
    linestr >> finalState >> sampleName >> nEvents >> xsec >> br;
    double weight = xsec * br / nEvents;
    if (dbg) {
      cout << "sample: " << sampleName << " final state: " << finalState << " #ev: " << nEvents
	   << " xsec: " << xsec << " br: " << br << " weight: " << weight << endl;    
    }
    bool set = false;
    if(finalState == "fourEle" || finalState == "4e" || finalState == "all") {
      theWeightMapFourE[sampleName] = weight;
      theInitialNEvFourE[sampleName] = nEvents;
      set = true;
    }  
    if(finalState == "fourMu"|| finalState == "4mu" || finalState == "all") {
      theWeightMapFourMu[sampleName] = weight;
      theInitialNEvFourMu[sampleName] = nEvents;
      set = true;
    } 
    if(finalState == "twoMutwoE"|| finalState == "2e2mu" || finalState == "all") {
      theWeightMap2Mu2E[sampleName] = weight;
      theInitialNEv2Mu2E[sampleName] = nEvents;
      set = true;
    }
    if (!set) {
      cerr << "[XSecReader]*** Error: final state in the map is invalid: \"" << finalState << "\" for sample " << sampleName << " in file " << xsecFileName << endl;
      cout << "input line " << iline << " : ->" << line << "<-" << endl;
    }
  }
  xsecFile.close();

  ifstream lumiFile(lumiFileName.Data());

  // Read the lumi map and store numbers in a map (per sample)
  while (getline(lumiFile,line)) {
    if( line == "" || line[0] == '#' ) continue; // Skip comments and empty lines
    stringstream linestr;
    linestr << line;
    TString finalState;
    TString epoch;
    bool dqApplied;
    double lumi = -1;
    linestr >> finalState >> epoch >> dqApplied >> lumi;
    if (dbg) {
      cout << "final state: " << finalState << " epoch: " << epoch << " DQ applied: " <<
	dqApplied << " lumi: " << lumi << endl;
    }
    if(dqApplied) {
      theLumiMapDQ[epoch][finalState] = lumi;
    } else {
      theLumiMapNoDQ[epoch][finalState] = lumi;
    }
  }
  lumiFile.close();

}


XSecReader::~XSecReader(){}


double XSecReader::getWeight(const TString& sampleName, const TString& epoch,
			     const TString& finalState, const bool dqApplied) const {
  const map<TString, double> *theWeightMap = 0;
  if(finalState == "fourEle" || finalState == "4e") {
    theWeightMap = &theWeightMapFourE;
  } else if(finalState == "fourMu" || finalState == "4mu" || finalState == "all") {
    theWeightMap = &theWeightMapFourMu;
  } else if(finalState == "twoMutwoE"|| finalState == "2e2mu") {
    theWeightMap = &theWeightMap2Mu2E;    
  }else {
    cerr << "[XSecReader]*** Error: final state is invalid: " << finalState << endl;
  }

  if(theWeightMap->find(sampleName) == theWeightMap->end()) cout << "111" << endl;
  if(theLumiMapDQ.find(epoch) == theLumiMapDQ.end()) cout << "112" << endl;
  if(theLumiMapDQ.find(epoch)->second.find(finalState) == theLumiMapDQ.find(epoch)->second.end()) cout << "113" << endl;

  if(theWeightMap->find(sampleName) == theWeightMap->end() ||
     theLumiMapDQ.find(epoch) == theLumiMapDQ.end() ||
     theLumiMapDQ.find(epoch)->second.find(finalState) == theLumiMapDQ.find(epoch)->second.end()) {
//     if(theWeightMap->find(sampleName) == theWeightMap->end()) 
//       cout << "[XSecReader]***ERROR: invalid sample name: " << sampleName << endl;
    cout << "[XSecReader]***ERROR: invalid labels; sampleName: " << sampleName
	 << " epoch: " << epoch
	 << " finalState: " << finalState
	 << " dqApplied: " << dqApplied << endl;
    return -1;
  }
  double lumi = -1;
  if(dqApplied) {
    lumi = ((theLumiMapDQ.find(epoch)->second).find(finalState))->second;
  } else {
    lumi = ((theLumiMapNoDQ.find(epoch)->second).find(finalState))->second;
  }
//   cout << "[XSecReader] lumi: " << lumi
//        << " weight: " << theWeightMap->find(sampleName)->second << endl;

  return lumi*theWeightMap->find(sampleName)->second;
}


double XSecReader::getInitNEv(const TString& sampleName, const TString& finalState) const {
  const map<TString, double> *theNEvMap = 0;
  if(finalState == "fourEle" || finalState == "4e") {
    theNEvMap = &theInitialNEvFourE;
  } else if(finalState == "fourMu" || finalState == "4mu" || finalState == "all") {
    theNEvMap = &theInitialNEvFourMu;
  }  else if(finalState == "twoMutwoE"|| finalState == "2e2mu") {
    theNEvMap = &theInitialNEv2Mu2E;
  }else {
    cerr << "[XSecReader]*** Error: final state is invalid: " << finalState << endl;
  }
  
  map<TString, double>::const_iterator elem = theNEvMap->find(sampleName);
  
  if(elem == theNEvMap->end()) {
    cout << "[XSecReader]***ERROR: invalid labels; sampleName: " << sampleName << endl;
    return -1;
  } else {
    return elem->second;
  }
}


double XSecReader::getLuminosity(const TString& epoch, const TString& finalState,
				 const bool dqApplied) const {
  double ret = -1;
  if(dqApplied) {
    ret = theLumiMapDQ.find(epoch)->second.find(finalState)->second;  
  } else {
    ret = theLumiMapNoDQ.find(epoch)->second.find(finalState)->second;
  }
  return ret;
}
