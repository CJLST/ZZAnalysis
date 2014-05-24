#ifndef XSecReader_H
#define XSecReader_H

/** \class XSecReader
 *  Reads the Xsection.txt and the Luminosity.txt files which 
 *  store the Xsections for the different MC samples and the luminosity
 *  for the different data samples. 
 *
 *  $Date: 2012/05/14 16:04:23 $
 *  $Revision: 1.1 $
 *  \author G. Cerminara - NEU Boston & INFN Torino
 */

#include <map>

class TString;


class XSecReader {
public:
  /// Constructor
  XSecReader(const TString& xsecFileName, const TString& lumiFileName);

  /// Destructor
  virtual ~XSecReader();
  
  // Operations
  // Get the weight = lumi*xsec*BR/# initial events 
  double getWeight(const TString& sampleName, const TString& epoch,
		   const TString& finalState, const bool dqApplied) const;
  // Get the weight = lumi*xsec*BR/# initial events 
  double getInitNEv(const TString& sampleName, const TString& finalState) const;
  
  double getLuminosity(const TString& epoch, const TString& finalState, const bool dqApplied) const;


  
protected:

private:
  std::map<TString, double> theWeightMapFourMu;
  std::map<TString, double> theWeightMapFourE;
  std::map<TString, double> theWeightMap2Mu2E;
  std::map<TString, double> theInitialNEvFourE;
  std::map<TString, double> theInitialNEvFourMu;
  std::map<TString, double> theInitialNEv2Mu2E;

  std::map<TString, std::map<TString, double> > theLumiMapDQ;
  std::map<TString, std::map<TString, double> > theLumiMapNoDQ;

};
#endif

