#include "ZZAnalysis/AnalysisStep/interface/Dataset.h"
#include <RtypesCore.h>
#include <regex>
#include <stdexcept>
#include <TFile.h>

using namespace std;

Dataset::Dataset(const TString& path, const TString& sampleName, bool delayOpen, int verbosityLevel)
  : Runs(nullptr),
    Events(nullptr),
    AllEvents(nullptr),
    sampleName_(sampleName),
    verbosityLevel_(verbosityLevel) {
  AddChunks(path, delayOpen);
}

Dataset::~Dataset() {
  delete Runs;
  delete Events;
  delete AllEvents;
}

void Dataset::AddChunks(const TString& path, bool delayOpen) {
  TSystemDirectory basedir(path, path);
  TList* subdirs = basedir.GetListOfFiles();
  if (!subdirs) {
    cerr << "[Dataset] Cannot list directory: " << path << endl;
    return;
  }

  regex chunkRe(sampleName_+"_Chunk[0-9]+$");
  
  TIter next(subdirs);
  int fileCount = 0;
  while (TObject* o = next()) {
    const char* name = o->GetName();

    Long64_t mode = 0;
    if (delayOpen) mode = TChain::kBigNumber;
    
    if (!name || name[0] == '.') continue; // skip . .. and hidden files
    if (regex_match(name, chunkRe)) {
      TString filePath = TString::Format("%s/%s/ZZ4lAnalysis.root",path.Data(), name);
	
      if (!gSystem->AccessPathName(filePath)) {
	// Check which trees are contained in the first file
	if (Events==nullptr) {
	  auto f = TFile::Open(filePath);
	  if (f->Get("Runs") != nullptr and f->Get("Runs") != nullptr) {
	    Runs = new TChain("Runs");
	    Events = new TChain("Events");
	  } else {
	    throw runtime_error("File "+filePath+" has no Runs/Events tree");
	  }
	  if (f->Get("AllEvents") != nullptr) { // optional tree
	    AllEvents = new TChain("AllEvents");
	  }
	  f->Close();
	}

	// Note: -1 forces reading the tree's header to read the number of events; otherwise
	int n1 = Runs->Add(filePath,mode);
	int n2 = Events->Add(filePath,mode);
	int n3 = -1;
	if (AllEvents!=nullptr) n3 = AllEvents->Add(filePath, mode);
	if (verbosityLevel_>=2) cout << "Adding " << filePath << endl;
	if (n1==0 || n2==0 || n3 ==0) {
	  throw runtime_error(TString("Tree ") +(n1==0?"Runs ":"")+(n2==0?"Events ":"")+(n3==0?"AllEvents ":"")+"not found in "+filePath);
	}
	fileCount++;
      } else {
	throw runtime_error("Missing file: "+filePath);
      }
    }
  }

  if (fileCount==0) {
    throw runtime_error("No valid tree found in  "+path+ " for datset "+sampleName_);
  } else {
    cout << "Dataset " << sampleName_ << ": opened " << fileCount << " files" << endl;
  }
}


TObject* Dataset::Get(const char* name) {
  TString req(name);
  if (req == "Events") return Events;
  else if (req == "AllEvents") return AllEvents;
  else return nullptr;
}
