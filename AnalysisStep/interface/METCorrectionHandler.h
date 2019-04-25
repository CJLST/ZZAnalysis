#ifndef METCORRECTIONHANDLER_H
#define METCORRECTIONHANDLER_H

#include <vector>
#include <unordered_map>
#include <utility>
#include <ZZAnalysis/AnalysisStep/interface/METObject.h>
#include <ZZAnalysis/AnalysisStep/interface/SystematicVariations.h>
#include "TString.h"
#include "TVar.hh"


class METCorrectionHandler{
protected:
  TVar::VerbosityLevel verbosity;
  bool applyCorrection;
  int theDataYear;
  TString theDataPeriod;
  std::vector<float> lumilist;
  std::vector<std::vector<std::pair<float, float>>> values_data_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, std::vector<std::vector<std::pair<float, float>>>> values_MC_map;

  void setDataPeriod(TString const& s);
  void readFile(TString const& strinput);

public:
  METCorrectionHandler(TString const& str_data_period);
  ~METCorrectionHandler();

  bool setup();
  void reset();

  bool hasMETCorrection() const{ return applyCorrection; }
  void correctMET(float const& genMET, float const& genMETPhi, METObject* obj, bool useFastSim) const;

  void printParameters() const;

  void setVerbosity(TVar::VerbosityLevel v){ verbosity=v; }

  std::vector<TString> getValidDataPeriods() const;
  static float getIntegratedLuminosity(TString const& period);
  static void autoExpandEnvironmentVariables(TString& path);

};



#endif
