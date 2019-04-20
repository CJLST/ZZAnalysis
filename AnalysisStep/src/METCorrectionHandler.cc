#include <cassert>
#include <string>
#include <regex>
#include <fstream>
#include <sstream>
#include <ZZAnalysis/AnalysisStep/interface/METCorrectionHandler.h>
#include "MELAStreamHelpers.hh"
#include "TDirectory.h"
#include "TRandom3.h"


using namespace std;
using namespace MELAStreamHelpers;


METCorrectionHandler::METCorrectionHandler(TString const& str_data_period) :
  verbosity(TVar::ERROR),
  applyCorrection(false)
{
  setDataPeriod(str_data_period);
  setup();
}

METCorrectionHandler::~METCorrectionHandler(){ this->reset(); }

std::vector<TString> METCorrectionHandler::getValidDataPeriods() const{
  std::vector<TString> res;
  if (theDataYear == 2016) res = std::vector<TString>{ "2016B", "2016C", "2016D", "2016E", "2016F", "2016G", "2016H" };
  else if (theDataYear == 2017) res = std::vector<TString>{ "2017B", "2017C", "2017D", "2017E", "2017F" };
  else if (theDataYear == 2018) res = std::vector<TString>{ "2018A", "2018B", "2018C", "2018D" };
  else{
    MELAerr << "SampleHelpers::getValidDataPeriods: Data periods for year " << theDataYear << " are undefined." << endl;
    assert(0);
  }
  return res;
}
float METCorrectionHandler::getIntegratedLuminosity(TString const& period){
  // To install brilcalc, do
  // export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda/bin:$PATH
  // pip install brilws --user --upgrade (pip install brilws --user if running for the first time)
  // and then run the different brilcalc lumi ... commands
  // Using brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt
  if (period == "2016") return 35.882515397;
  else if (period == "2016B") return 5.711130443;
  else if (period == "2016C") return 2.572903492;
  else if (period == "2016D") return 4.242291558;
  else if (period == "2016E") return 4.025228139;
  else if (period == "2016F") return 3.104509131;
  else if (period == "2016G") return 7.575824256;
  else if (period == "2016H") return 8.650628378;
  // Using brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt
  else if (period == "2017") return 41.529152052;
  else if (period == "2017B") return 4.793969901;
  else if (period == "2017C") return 9.632746391;
  else if (period == "2017D") return 4.247792713;
  else if (period == "2017E") return 9.314581018;
  else if (period == "2017F") return 13.540062029;
  //// Using brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PREAPPROVED.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt
  // Using brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt
  else if (period == "2018") return 59.740565209;
  else if (period == "2018A") return 14.027614284;
  else if (period == "2018B") return 7.066552173;
  else if (period == "2018C") return 6.898816878;
  else if (period == "2018D") return 31.747581874;
  else if (period == "2018_HEMaffected") return 38.662770627;
  else{
    MELAerr << "SampleHelpers::getIntegratedLuminosity(" << period << "): Period is not defined." << endl;
    assert(0);
    return -1;
  }
}
void METCorrectionHandler::setDataPeriod(TString const& s){
  theDataPeriod = s;
  theDataYear = -1;
  if (theDataPeriod.Contains("2016")) theDataYear = 2016;
  else if (theDataPeriod.Contains("2017")) theDataYear = 2017;
  else if (theDataPeriod.Contains("2018")) theDataYear = 2018;
  else{
    MELAerr << "SampleHelpers::setDataPeriod: Could not recognize the data period string " << s << " to assign the data year." << endl;
    assert(0);
  }
}
void METCorrectionHandler::autoExpandEnvironmentVariables(TString& path){
  std::string strpath=path.Data();
  static std::regex env("\\$\\{([^}]+)\\}");
  std::smatch match;
  while (std::regex_search(strpath, match, env)){
    const char* s = getenv(match[1].str().c_str());
    const std::string var(s == NULL ? "" : s);
    strpath.replace(match[0].first, match[0].second, var);
  }
  path=strpath;
}

bool METCorrectionHandler::setup(){
  using namespace SystematicsHelpers;

  bool res = true;
  TDirectory* curdir = gDirectory;

  this->reset();

  TString _pkgpath_="${CMSSW_BASE}/src/ZZAnalysis/AnalysisStep/data/";
  autoExpandEnvironmentVariables(_pkgpath_);
  TString strdatacore = Form("%s%s/%i/", _pkgpath_.Data(), "METSFs", theDataYear);
  if (theDataPeriod == ""){
    MELAerr << "METCorrectionHandler::setup: Need to have a valid data period set." << endl;
    assert(0);
  }

  std::vector<TString> strValidPeriods = getValidDataPeriods();
  const size_t nValidPeriods = strValidPeriods.size();

  if (theDataYear == 2016){
    applyCorrection = true;
    // Use reserve since the push is done at the reading stage.
    lumilist.reserve(nValidPeriods);
    values_data_map.reserve(nValidPeriods);
    values_MC_map[sNominal] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[sNominal].reserve(nValidPeriods);
    values_MC_map[eJECDn] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJECDn].reserve(nValidPeriods);
    values_MC_map[eJECUp] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJECUp].reserve(nValidPeriods);
    values_MC_map[eJERDn] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJERDn].reserve(nValidPeriods);
    values_MC_map[eJERUp] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJERUp].reserve(nValidPeriods);
    values_MC_map[ePUDn] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[ePUDn].reserve(nValidPeriods);
    values_MC_map[ePUUp] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[ePUUp].reserve(nValidPeriods);
    for (TString const& period:strValidPeriods){
      if (theDataPeriod == "2016" || theDataPeriod == period){
        if (verbosity >= TVar::DEBUG) MELAout << "METCorrectionHandler::setup: Adding data period " << period << "..." << endl;
        if (period != "2016B") this->readFile(strdatacore+Form("fitparameters_Run%s-17Jul2018_MC.txt", period.Data()));
        else this->readFile(strdatacore+Form("fitparameters_Run%s-17Jul2018_ver2_MC.txt", period.Data()));
        lumilist.push_back(getIntegratedLuminosity(period));
      }
    }
    // Accumulate integrated luminosity and divide by the last one
    for (unsigned int il=1; il<lumilist.size(); il++) lumilist.at(il) += lumilist.at(il-1);
    for (unsigned int il=0; il<lumilist.size(); il++) lumilist.at(il) /= lumilist.back();
  }
  else if (theDataYear == 2017){
    applyCorrection = true;
    // Use reserve since the push is done at the reading stage.
    lumilist.reserve(nValidPeriods);
    values_data_map.reserve(nValidPeriods);
    values_MC_map[sNominal] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[sNominal].reserve(nValidPeriods);
    values_MC_map[eJECDn] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJECDn].reserve(nValidPeriods);
    values_MC_map[eJECUp] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJECUp].reserve(nValidPeriods);
    values_MC_map[eJERDn] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJERDn].reserve(nValidPeriods);
    values_MC_map[eJERUp] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJERUp].reserve(nValidPeriods);
    values_MC_map[ePUDn] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[ePUDn].reserve(nValidPeriods);
    values_MC_map[ePUUp] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[ePUUp].reserve(nValidPeriods);
    for (TString const& period:strValidPeriods){
      if (theDataPeriod == "2017" || theDataPeriod == period){
        if (verbosity >= TVar::DEBUG) MELAout << "METCorrectionHandler::setup: Adding data period " << period << "..." << endl;
        this->readFile(strdatacore+Form("fitparameters_Run%s-31Mar2018-v1_MC.txt", period.Data())); lumilist.push_back(getIntegratedLuminosity(period));
      }
    }
    if (theDataPeriod == "2017F-09May2018"){
      if (verbosity >= TVar::DEBUG) MELAout << "METCorrectionHandler::setup: Adding data period Run2017F-09May2018..." << endl;
      this->readFile(strdatacore+"fitparameters_Run2017F-09May2018-v1_MC.txt"); lumilist.push_back(getIntegratedLuminosity("2017F"));
    }
    // Accumulate integrated luminosity and divide by the last one
    for (unsigned int il=1; il<lumilist.size(); il++) lumilist.at(il) += lumilist.at(il-1);
    for (unsigned int il=0; il<lumilist.size(); il++) lumilist.at(il) /= lumilist.back();
  }
  else if (theDataYear == 2018){
    applyCorrection = true;
    // Use reserve since the push is done at the reading stage.
    lumilist.reserve(nValidPeriods);
    values_data_map.reserve(nValidPeriods);
    values_MC_map[sNominal] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[sNominal].reserve(nValidPeriods);
    values_MC_map[eJECDn] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJECDn].reserve(nValidPeriods);
    values_MC_map[eJECUp] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJECUp].reserve(nValidPeriods);
    values_MC_map[eJERDn] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJERDn].reserve(nValidPeriods);
    values_MC_map[eJERUp] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[eJERUp].reserve(nValidPeriods);
    values_MC_map[ePUDn] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[ePUDn].reserve(nValidPeriods);
    values_MC_map[ePUUp] = std::vector<std::vector<std::pair<float, float>>>(); values_MC_map[ePUUp].reserve(nValidPeriods);
    for (TString const& period:strValidPeriods){
      if (theDataPeriod == "2018" || theDataPeriod == period){
        if (verbosity >= TVar::DEBUG) MELAout << "METCorrectionHandler::setup: Adding data period " << period << "..." << endl;
        if (period != "2018D") this->readFile(strdatacore+Form("fitparameters_Run%s-17Sep2018_MC.txt", period.Data()));
        else this->readFile(strdatacore+Form("fitparameters_Run%s-PromptReco_MC.txt", period.Data()));
        lumilist.push_back(getIntegratedLuminosity(period));
      }
    }
    // Accumulate integrated luminosity and divide by the last one
    for (unsigned int il=1; il<lumilist.size(); il++) lumilist.at(il) += lumilist.at(il-1);
    for (unsigned int il=0; il<lumilist.size(); il++) lumilist.at(il) /= lumilist.back();
  }

  curdir->cd();

  return res;
}
void METCorrectionHandler::reset(){
  applyCorrection = false;
  lumilist.clear();
  values_data_map.clear();
  values_MC_map.clear();
}

void METCorrectionHandler::readFile(TString const& strinput){
  using namespace SystematicsHelpers;

  // Notice that the systematics list does not contain eMETUp and eMETDn. This is because MET up and down variations involve varying the nominal values by the given errors.
  std::vector<SystematicVariationTypes> const validSysts={
    sNominal, eJECUp, eJECDn, eJERUp, eJERDn, ePUUp, ePUDn
  };
  const size_t nCopies = validSysts.size()+1;
  vector<pair<float, float>> value_error_S1; value_error_S1.reserve(nCopies);
  vector<pair<float, float>> value_error_S2; value_error_S2.reserve(nCopies);
  vector<pair<float, float>> value_error_S3; value_error_S3.reserve(nCopies);
  vector<pair<float, float>> value_error_F1; value_error_F1.reserve(nCopies);
  vector<pair<float, float>> value_error_F2; value_error_F2.reserve(nCopies);

  ifstream fin;
  fin.open(strinput.Data(), std::ios_base::in);
  if (fin.good()){
    string str_in;
    while (!fin.eof()){
      getline(fin, str_in);
      vector<pair<float, float>>* input_par_pair=nullptr;
      if (str_in.find("Sigma 1")!=string::npos) input_par_pair = &value_error_S1;
      else if (str_in.find("Sigma 2")!=string::npos) input_par_pair = &value_error_S2;
      else if (str_in.find("Sigma 3")!=string::npos) input_par_pair = &value_error_S3;
      else if (str_in.find("Frac 1")!=string::npos) input_par_pair = &value_error_F1;
      else if (str_in.find("Frac 2")!=string::npos) input_par_pair = &value_error_F2;
      if (input_par_pair){
        string dummy;
        float val, err;
        stringstream ss(str_in);
        //MELAout << str_in << endl;
        ss >> dummy >> dummy >> dummy >> val >> dummy >> err;
        input_par_pair->emplace_back(val, err);
      }
    }
  }
  fin.close();

  if (value_error_S1.size()!=nCopies || value_error_S2.size()!=nCopies || value_error_S3.size()!=nCopies || value_error_F1.size()!=nCopies || value_error_F2.size()!=nCopies){
    MELAerr << "METCorrectionHandler::readFile(" << strinput << "): Some values do not have size " << nCopies << endl;
    assert(0);
  }

  values_data_map.push_back(std::vector<std::pair<float, float>>());
  values_data_map.back().reserve(5);
  values_data_map.back().emplace_back(value_error_S1.front());
  values_data_map.back().emplace_back(value_error_S2.front());
  values_data_map.back().emplace_back(value_error_S3.front());
  values_data_map.back().emplace_back(value_error_F1.front());
  values_data_map.back().emplace_back(value_error_F2.front());

  size_t whichElement = 1;
  for (auto& syst:validSysts){
    auto* it = &(values_MC_map.find(syst)->second); it->push_back(std::vector<std::pair<float, float>>());
    it->back().reserve(5);
    it->back().emplace_back(value_error_S1.at(whichElement));
    it->back().emplace_back(value_error_S2.at(whichElement));
    it->back().emplace_back(value_error_S3.at(whichElement));
    it->back().emplace_back(value_error_F1.at(whichElement));
    it->back().emplace_back(value_error_F2.at(whichElement));
    whichElement++;
  }
}

void METCorrectionHandler::correctMET(float const& genMET, float const& genMETPhi, METObject* obj, bool /*useFastSim*/) const{
  using namespace SystematicsHelpers;

  if (!applyCorrection) return;
  if (!obj) return;

  std::vector<SystematicVariationTypes> const validSysts={
    sNominal, eMETDn, eMETUp, eJECDn, eJECUp, eJERDn, eJERUp, ePUDn, ePUUp
  };
  const TLorentzVector v_genMET(genMET*cos(genMETPhi), genMET*sin(genMETPhi), 0, genMET);

  TRandom3 rand;
  rand.SetSeed(std::abs(static_cast<int>(sin(genMETPhi)*100000)));
  float era_x = rand.Uniform();
  float frac_x = rand.Uniform();
  int i_era = -1;
  int i_frac = -1;
  for (auto const& lumi_era:lumilist){
    i_era++;
    if (era_x<=lumi_era) break;
  }
  if (i_era<0){
    MELAerr << "METCorrectionHandler::correctMET: i_era = " << i_era << " < 0." << endl;
    assert(0);
  }
  else if (i_era>=(int) values_data_map.size()){
    MELAerr << "METCorrectionHandler::correctMET: i_era = " << i_era << " >= size of eras = " << values_data_map.size() << endl;
    assert(0);
  }
  std::vector<std::pair<float, float>> const& data_values = values_data_map.at(i_era);
  float const& F1 = data_values.at(3).first;
  float const& F2 = data_values.at(4).first;
  if (frac_x<=F1) i_frac = 0;
  else if (frac_x<=(1.f - F1)*F2+F1) i_frac = 1;
  else i_frac = 2;

  if (verbosity>=TVar::DEBUG) MELAout << "METCorrectionHandler::correctMET: The i_era and i_frac random numbers: " << i_era << ", " << i_frac << endl;

  bool correctUpperError=false;
  float scale_nominal=1.f;
  float* met_value=nullptr;
  float* phi_value=nullptr;
  for (auto& syst:validSysts){
    auto systEff = syst;
    switch (syst){
    case sNominal:
      met_value = &(obj->extras.met);
      phi_value = &(obj->extras.phi);
      break;
    case eMETDn:
      systEff = sNominal;
      met_value = &(obj->extras.met_METdn);
      phi_value = &(obj->extras.phi_METdn);
      break;
    case eMETUp:
      systEff = sNominal;
      met_value = &(obj->extras.met_METup);
      phi_value = &(obj->extras.phi_METup);
      break;
    case eJECDn:
      met_value = &(obj->extras.met_JECdn);
      phi_value = &(obj->extras.phi_JECdn);
      break;
    case eJECUp:
      met_value = &(obj->extras.met_JECup);
      phi_value = &(obj->extras.phi_JECup);
      break;
    case eJERDn:
      met_value = &(obj->extras.met_JERdn);
      phi_value = &(obj->extras.phi_JERdn);
      break;
    case eJERUp:
      met_value = &(obj->extras.met_JERup);
      phi_value = &(obj->extras.phi_JERup);
      break;
    case ePUDn:
      met_value = &(obj->extras.met_PUdn);
      phi_value = &(obj->extras.phi_PUdn);
      break;
    case ePUUp:
      met_value = &(obj->extras.met_PUup);
      phi_value = &(obj->extras.phi_PUup);
      break;
    default:
      MELAout << "METCorrectionHandler::correctMET: Unsupported systematic " << syst << " in switch statement." << endl;
      assert(0);
      break;
    }

    TLorentzVector v_MET(
      *met_value * cos(*phi_value),
      *met_value * sin(*phi_value),
      0,
      *met_value
    );
    if (verbosity>=TVar::DEBUG) MELAout << "METCorrectionHandler::correctMET: The old MET vector: " << v_MET << endl;

    v_MET = v_MET - v_genMET;
    //v_MET = v_MET - TLorentzVector(genMET*cos(*phi_value), genMET*sin(*phi_value), 0, genMET);
    if (verbosity>=TVar::DEBUG) MELAout << "METCorrectionHandler::correctMET: The MET difference vector: " << v_MET << endl;

    auto const& mc_values = values_MC_map.find(systEff)->second;
    float data_sigma = data_values.at(i_frac).first;
    float mc_sigma = mc_values.at(i_era).at(i_frac).first;
    if (syst == eMETDn){
      data_sigma -= data_values.at(i_frac).second;
      mc_sigma += mc_values.at(i_era).at(i_frac).second;
    }
    else if (syst == eMETUp){
      data_sigma += data_values.at(i_frac).second;
      mc_sigma -= mc_values.at(i_era).at(i_frac).second;
    }

    float scale = data_sigma/mc_sigma;
    if (syst == sNominal){
      scale_nominal = scale;
      correctUpperError = (scale_nominal<1.f);
    }
    // When correcting the METup value to 1, require that the deviation of the nominal scale factor is <7 sigma. Otherwise, leave the METup scale as it is.
    if (correctUpperError && syst == eMETUp && scale<1.f){
      float scale_diff = fabs(scale - scale_nominal);
      float scale_diff_nominal = fabs(1.f - scale_nominal);
      if (scale_diff>0.f && scale_diff_nominal/scale_diff<7.f) scale = 1.f;
    }

    if (verbosity>=TVar::DEBUG) MELAout << "METCorrectionHandler::correctMET: The MET scale for systematic " << syst << ": " << data_sigma << " / " << mc_sigma << " = " << scale << endl;

    v_MET *= scale;
    v_MET += v_genMET;
    *met_value = v_MET.Pt();
    *phi_value = v_MET.Phi();

    if (verbosity>=TVar::DEBUG) MELAout << "METCorrectionHandler::correctMET: The new MET vector: " << v_MET << endl;
  }
}

void METCorrectionHandler::printParameters() const{
  using namespace SystematicsHelpers;

  if (!applyCorrection) return;

  MELAout << "Begin METCorrectionHandler::printParameters" << endl;
  MELAout << "Normalized lumi list: " << lumilist << endl;
  unsigned int i_era = 0;
  for (auto const& ed:values_data_map){
    MELAout << "Era " << i_era << " values in data: " << ed << endl;
    i_era++;
  }

  std::vector<SystematicVariationTypes> const validSysts={
    sNominal, eJECDn, eJECUp, eJERDn, eJERUp, ePUDn, ePUUp
  };
  for (auto const& syst:validSysts){
    i_era=0;
    for (auto const& ed:values_MC_map.find(syst)->second){
      MELAout << "Era " << i_era << " values in MC[" << syst << "]: " << ed << endl;
      i_era++;
    }
  }
  MELAout << "End METCorrectionHandler::printParameters" << endl;
}
