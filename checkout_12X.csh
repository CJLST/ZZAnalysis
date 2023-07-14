#!/bin/bash 
#
# Instructions:
# wget -O ${TMPDIR}/checkout_12X.csh https://raw.githubusercontent.com/CJLST/ZZAnalysis/Run3/checkout_12X.csh
# cd $CMSSW_BASE/src
# cmsenv
# chmod u+x ${TMPDIR}/checkout_12X.csh
# ${TMPDIR}/checkout_12X.csh

############## For CMSSW_12_6_5

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
#git cms-addpkg EgammaAnalysis/ElectronTools
#(rm -rf EgammaAnalysis/ElectronTools/data;git clone https://github.com/jainshilpi/EgammaAnalysis-ElectronTools.git -b ULSSfiles_correctScaleSysMC EgammaAnalysis/ElectronTools/data;)

#### Please do not add any custom (non-CMSSW) package before this line ####
#ZZAnalysis
git clone https://github.com/CJLST/ZZAnalysis.git ZZAnalysis
(cd ZZAnalysis; git checkout Run3)

# Muon MVA
#git cms-addpkg CondFormats/EgammaObjects
#git cms-addpkg CommonTools/MVAUtils
#git clone https://github.com/bonanomi/MuonMVAReader.git MuonMVAReader
#(cd MuonMVAReader; git checkout 3d53269)

#Common LHE tools
git clone https://github.com/usarica/CommonLHETools.git
(cd CommonLHETools; git checkout -b from-v142 v1.4.2)

#MELA
git clone https://github.com/JHUGen/JHUGenMELA.git JHUGenMELA
(cd JHUGenMELA; git checkout -b from-v238 v2.3.8; ./setup.sh)

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

#kinematic refitting (obsolete?)
git clone https://github.com/mhl0116/KinZfitter-1.git KinZfitter
(cd KinZfitter ; git checkout -b from-27daebb 27daebb)

#NanoAODTools
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
(cd PhysicsTools/NanoAODTools ; git checkout -b from-c32f055 c32f055)

#CommonLHETools requires the MELA env to be set for compilation
(eval `MelaAnalytics/setup.sh env`; cd CommonLHETools; scram b -j4)

#Now we can compile everything
scram b -j4

