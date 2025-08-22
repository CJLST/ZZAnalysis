#!/bin/bash 
#
# Instructions:
# cmsenv
# wget -O ${TMPDIR}/checkout.csh https://raw.githubusercontent.com/CJLST/ZZAnalysis/Run3/checkout_13X.csh
# cd $CMSSW_BASE/src
# chmod u+x ${TMPDIR}/checkout.csh
# ${TMPDIR}/checkout.csh

############## For CMSSW_13_3_3 (2022 data/MC) or CMSSW_14_1_6 (2023 data/MC)

#exit when any command fails
set -e

git cms-init

# New Jet PU ID: dedicated training for each year
#git cms-addpkg  RecoJets/JetProducers

# STXS Categorisation: now directly implemented in CMSSW
#git cms-addpkg GeneratorInterface/RivetInterface
#git cms-addpkg SimDataFormats/HTXS

# Updated for UL. See: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018
git clone -b run3ID https://github.com/swagata87/EgammaPostRecoTools.git  EgammaUser/EgammaPostRecoTools
#FIXME: hack to run in 13X; recipe should be fixed
sed -i s/"12 : \[0,1,6\]"/"12 : [0,1,6],13 : [0,1,2,3]"/ EgammaUser/EgammaPostRecoTools/python/EgammaPostRecoTools.py

#### Please do not add any custom (non-CMSSW) package before this line ####
#ZZAnalysis
git clone https://github.com/CJLST/ZZAnalysis.git ZZAnalysis
(cd ZZAnalysis; git checkout Run3)

# Muon MVA
#git cms-addpkg CondFormats/EgammaObjects
#git cms-addpkg CommonTools/MVAUtils
#git clone https://github.com/bonanomi/MuonMVAReader.git MuonMVAReader
#(cd MuonMVAReader; git checkout 3d53269)

#Common LHE tools (private FW update, based on v1.4.2)
git clone https://github.com/pianonick411/CommonLHETools
(cd CommonLHETools; git checkout -b from-v135 v1.3.5)

#MELA
git clone https://github.com/JHUGen/JHUGenMELA.git JHUGenMELA
(cd JHUGenMELA; git checkout -b from-v243 v2.4.3; ./setup.sh)

#MELA Analytics
git clone https://github.com/MELALabs/MelaAnalytics.git
(cd MelaAnalytics; git checkout -b from-v23 v2.3; ./setup.sh)

#Move MELA libraries to the proper place, so that we can avoid its silly
#env settings 
mkdir -p ${CMSSW_BASE}/lib/${SCRAM_ARCH}
ln -s ${CMSSW_BASE}/src/JHUGenMELA/MELA/data/*/*.so \
      ${CMSSW_BASE}/src/MelaAnalytics/CandidateLOCaster/lib/*.so \
      ${CMSSW_BASE}/src/MelaAnalytics/GenericMEComputer/lib/*.so \
      ${CMSSW_BASE}/src/MelaAnalytics/EventContainer/lib/*.so \
      ${CMSSW_BASE}/lib/${SCRAM_ARCH}

#kinematic refitting
git clone https://github.com/ferrico/KinZfitter.git KinZfitter
(cd KinZfitter ; git checkout -b from-f781891 f781891)
sed -i '/SimTracker\/Records/d' KinZfitter/HelperFunction/BuildFile.xml
sed -i '/SimTracker\/Records/d' KinZfitter/KinZfitter/BuildFile.xml
sed -i '/#include "RooMinuit.h"/d' KinZfitter/KinZfitter/interface/KinZfitter.h

#Fix some memory issues (#47289 queued CMSSW_15_0_X))
git cms-addpkg PhysicsTools/NanoAODTools
git fetch https://github.com/namapane/cmssw.git NAT-dev2:namapane_NAT-dev2
git cherry-pick aa9ecbd04d6 98f8692142f

if [[ $CMSSW_VERSION == CMSSW_13_3_* ]]; then
 #Pick the fix from #43536 (haddNano.py); in release since 13_0_18, 14_0_2, 14_1_0; it was not backported to 13_3_X
 git cms-addpkg PhysicsTools/NanoAOD
 git cms-cherry-pick-pr 43536 CMSSW_13_0_X
   
 #Pick more fixes in NanoAODTools (support for "S" branch types); in release since 14_0_0, 14_1_0; queued in 13_3_X (X>3)
 git fetch https://github.com/namapane/cmssw.git NAT-dev:namapane_NAT-dev
 git cherry-pick 3e73ca4c2f8
fi
   
#get nanoAODTools modules
git clone https://github.com/cms-cat/nanoAOD-tools-modules.git PhysicsTools/NATModules
(cd PhysicsTools/NATModules; git checkout -b from-906b58d 906b58d)


#CommonLHETools requires the MELA env to be set for compilation
(eval `MelaAnalytics/setup.sh env`; cd CommonLHETools; scram b -j4)

#Now we can compile everything
scram b -j4

#Set up protection for user pulling but missing that this recipe has been updated
python3 ZZAnalysis/AnalysisStep/python/validateCheckout.py
