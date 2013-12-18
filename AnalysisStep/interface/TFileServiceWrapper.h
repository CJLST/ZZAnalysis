#ifndef TFileServiceWrapper_h
#define TFileServiceWrapper_h

/** \class TFileServiceWrapper
 *
 *  A wrapper around TFileService, that can be replaced with something working in bare root.
 *
 *  $Date: 2012/06/19 22:41:13 $
 *  $Revision: 1.1 $
 *  \author N. Amapane - CERN
 */

#include <TH1F.h>
#include <vector>
#include <FWCore/ServiceRegistry/interface/Service.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>

class TFileServiceWrapper {
 public:
  /// Constructor
  TFileServiceWrapper(){}

  /// Destructor
  virtual ~TFileServiceWrapper(){}
  
  // Operations
  TH1F* makeTH1F(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup){
    edm::Service<TFileService> fs;
    return fs->make<TH1F>(name, title, nbinsx, xlow, xup);
  }
};
#endif

