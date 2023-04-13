/** \file
 *
 *  \author T. Sculac
 */

#include <ZZAnalysis/AnalysisStep/interface/PhotonIDHelper.h>

#include <iostream>
#include <map>

using namespace std;
using namespace edm;
using namespace pat;
using namespace reco;

bool PhotonIDHelper::isCutBasedID_Loose(int year, const pat::Photon& photon){
   
   bool isID = false;
   
   // [FIXME] Removing photon ID for the moment since we do not use it in the analysis. Be careful if you want to use it, because the code will run locally but if you try to submit the jobs on Condor they will fail becua replaces "-" with "_". I have contacted EGamma to provide a fix for this so maybe check if a new version of photon ID that does not have this problem exists before using it.
   
   if (year == 2016) isID = photon.photonID("cutBasedPhotonID_Fall17_94X_V2_loose");
   if (year == 2017) isID = photon.photonID("cutBasedPhotonID_Fall17_94X_V2_loose");
   if (year == 2018) isID = photon.photonID("cutBasedPhotonID_Fall17_94X_V2_loose");
   
   return isID;
}

bool PhotonIDHelper::isCutBasedID_Medium(int year, const pat::Photon& photon){
   
   bool isID = false;
   
   // [FIXME] Removing photon ID for the moment since we do not use it in the analysis. Be careful if you want to use it, because the code will run locally but if you try to submit the jobs on Condor they will fail becua replaces "-" with "_". I have contacted EGamma to provide a fix for this so maybe check if a new version of photon ID that does not have this problem exists before using it.
   
   if (year == 2016) isID = photon.photonID("cutBasedPhotonID_Fall17_94X_V2_medium");
   if (year == 2017) isID = photon.photonID("cutBasedPhotonID_Fall17_94X_V2_medium");
   if (year == 2018) isID = photon.photonID("cutBasedPhotonID_Fall17_94X_V2_medium");
   
   return isID;
}

bool PhotonIDHelper::isCutBasedID_Tight(int year, const pat::Photon& photon){
   
   bool isID = false;
   
   // [FIXME] Removing photon ID for the moment since we do not use it in the analysis. Be careful if you want to use it, because the code will run locally but if you try to submit the jobs on Condor they will fail becua replaces "-" with "_". I have contacted EGamma to provide a fix for this so maybe check if a new version of photon ID that does not have this problem exists before using it.
   
   if (year == 2016) isID = photon.photonID("cutBasedPhotonID_Fall17_94X_V2_tight");
   if (year == 2017) isID = photon.photonID("cutBasedPhotonID_Fall17_94X_V2_tight");
   if (year == 2018) isID = photon.photonID("cutBasedPhotonID_Fall17_94X_V2_tight");
   
   return isID;
}