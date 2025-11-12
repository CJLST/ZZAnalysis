#ifndef DATASET_H
#define DATASET_H

#include <TChain.h>
#include <TString.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TCollection.h>

#include <iostream>

/**
 * Dataset
 * ----------
 * Utility that mimics the interface of a TFile containing the Runs, Events and 
 * (optionally) AllEvents trees, setting up TChains from all Chunks of a 
 * given sample in the specified folder.
 *
 */
class Dataset  {
public:
  /// Construct from path and sample name, eg: "ggH125".
  /// Regex-style wildcards allowed, e.g.: ".*2022." = all 2022 data (all PDs, all eras)
  /// delayOpen: speeds up initialization by skipping file opening and reading headers when they are added.
  ///            This result in much faster initialization, but:
  ///            - File errors (corrupted files, missing trees etc) will be found at the first access 
  ///            - GetEntriesFast() will return TChain::kBigNumber, regardless of the real number of entries,
  ///              until the chain is opened for the first time or GetEntries() is called.
  /// verbosityLevel: 0 = only err; 2=debug
  explicit Dataset(const TString& path,
		   const TString& sampleName,
		   bool delayOpen=false,
		   int verbosityLevel=0);

  ~Dataset();

  /// (Re)scan and add all Chunk*/ZZ4lAnalysis.root under baseDir (non-recursive).
  void AddChunks(const TString& baseDir, bool delayOpen);

  /// Mimic TFile::Get: return the chain if user asks the tree by name
  TObject* Get(const char* name);

  /// Mimic direct tree access in interactive root and PyROOT
  TChain* Runs;
  TChain* Events;
  TChain* AllEvents;
  
private:
  TString sampleName_;
  int verbosityLevel_;

};

#endif
