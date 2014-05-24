#ifndef LumiNormalization_H
#define LumiNormalization_H

/** \class LumiNormalization
 *  Class to get the normalization of each sample.
 *  The relative normalization is get from MC xsections while the overall normalization
 *  can be derived from the number of events under the Z peak.
 *  The samples used in this luminosity normalization must be added explicitly throught the add method.
 *  Look to the XSection.txt and Luminosity.txt files for naming conventions.
 *
 *  $Date: 2012/06/20 22:34:01 $
 *  $Revision: 1.2 $
 *  \author G. Cerminara - NEU Boston & INFN Torino
 */

#include "TString.h"
#include "Number.h"
#include <vector>

class XSecReader;

class LumiNormalization {
public:
  /// Constructor
  // Create a tool for normalization. You need to specify the XSection file and the lumi file
  // and the elements to look for the right lumi in the map.
  // also the inputDir for the histo to be used for data/MC normalization in required (che be a null string)
  LumiNormalization(const TString& xsecFileName, const TString& lumiFileName,
		    const TString& epoch, const TString& finalState, bool dqApplied,
		    const TString& inputDir);

  /// Destructor
  virtual ~LumiNormalization();

  // Operations
  // Add data to the count of events to be used for data/MC normalization
  // the histo used to get the count is zmass_restricted_1
  // the default file for this histo is: histos_zzanalysis_data_FINALSTATE.root
  void addData();


  // Add data to the count of events to be used for data/MC normalization
  // the histo used to get the count is zmass_restricted_1
  // This method allows to specify the name of the file which contains this histo
  void addData(const TString& filename);


  // Add MC sample to the count of events to be used for data/MC normalization
  // the histo used to get the count is zmass_restricted_1
  // the default file for this histo is: theInputDir+"/histos_zzanalysis_mc_"+sampleName+"_"+theFinalState+".root"
  void addMC(const TString& sampleName);


  // Add MC sample to the count of events to be used for data/MC normalization
  // the histo used to get the count is zmass_restricted_1
  // This method allows to specify the name of the file which contains this histo
  void addMC(const TString& sampleName, const TString& filename);



  // Add MC samples to the count of events to be used for data/MC normalization
  // the histo used to get the count is zmass_restricted_1
  // the default file for this histo is: theInputDir+"/histos_zzanalysis_mc_"+sampleName+"_"+theFinalState+".root"
  void addMC(const std::vector<TString>& sampleNames);



  // Get the relative normalization of MC/data under the Z peak
  double getNormalizationFactor() const;


  // Get the overall scale factor (relative normalization of MC/data under the Z peak + theo. xsec)
  double getScaleFactor(const TString& sampleName) const;


  // add an additional scale facto
  void setAdditionalScale(double scale);


  // switch off the normalization of MC/data under the Z peak
  void normalizeToZPeak(bool doNormalize);


  // set the name of the file histo to be used for the normalization of MC/data under the Z peak
  void setNormalizHistoName(const TString& histoname);


  // Get initial # of events (as in XSection.txt scaled for lumi e xsec)
  Number getInitialNEv(const TString& sampleName) const;

  
  static double getLuminosity();
  
  static int getCME() {
    return theCME;
  }
      

protected:

private:
  double nData;
  double nMC;
  TString theInputDir;

  TString theEpoch;
  TString theFinalState;
  TString theNormHistoName;
  bool isDqApplied;
  double additionalScale;
  bool normalizeToZ;

  XSecReader *xsecRead;
  static double theLumi;
  static int theCME;
  
};
#endif

