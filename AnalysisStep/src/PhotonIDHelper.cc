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
   
   if (year == 2016) isID = photon.photonID("cutBasedPhotonID_Fall17_94X_V2_loose");
   if (year == 2017) isID = photon.photonID("cutBasedPhotonID_Fall17_94X_V2_loose");
   if (year == 2018) isID = photon.photonID("cutBasedPhotonID_Fall17_94X_V2_loose");
   
   return isID;
}

bool PhotonIDHelper::isCutBasedID_Medium(int year, const pat::Photon& photon){
   
   bool isID = false;
      
   if (year == 2016) isID = photon.photonID("cutBasedPhotonID_Fall17_94X_V2_medium");
   if (year == 2017) isID = photon.photonID("cutBasedPhotonID_Fall17_94X_V2_medium");
   if (year == 2018) isID = photon.photonID("cutBasedPhotonID_Fall17_94X_V2_medium");
   
   return isID;
}

bool PhotonIDHelper::isCutBasedID_Tight(int year, const pat::Photon& photon){
   
   bool isID = false;
      
   if (year == 2016) isID = photon.photonID("cutBasedPhotonID_Fall17_94X_V2_tight");
   if (year == 2017) isID = photon.photonID("cutBasedPhotonID_Fall17_94X_V2_tight");
   if (year == 2018) isID = photon.photonID("cutBasedPhotonID_Fall17_94X_V2_tight");
   
   return isID;
}