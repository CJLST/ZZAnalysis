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
   
   // [FIXME] Removing photon ID for the moment since we do not use it in the analysis. Be careful if you want to use it, because the code will run locally but if you try to submit the jobs on Condor they will fail becuase dumpPython() replaces "-" with "_". I have contacted EGamma to provide a fix for this so maybe check if a new version of photon ID that does not have this problem exists before using it.
   
   if (year == 2016) isID = false; //photon.photonID("cutBasedPhotonID-Fall17-94X-V2-loose");
   if (year == 2017) isID = false; //photon.photonID("cutBasedPhotonID-Fall17-94X-V2-loose");
   if (year == 2018) isID = false; //photon.photonID("cutBasedPhotonID-Fall17-94X-V2-loose");
   
   return isID;
}
