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
   
   if (year == 2016) isID = photon.photonID("cutBasedPhotonID-Fall17-94X-V2-loose");
   if (year == 2017) isID = photon.photonID("cutBasedPhotonID-Fall17-94X-V2-loose");
   if (year == 2018) isID = photon.photonID("cutBasedPhotonID-Fall17-94X-V2-loose");
   
   return isID;
}
