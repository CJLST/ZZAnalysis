#!/bin/tcsh -fe
#
# Instructions:
# wget -O /tmp/checkout_80X.csh https://raw.githubusercontent.com/CJLST/ZZAnalysis/miniAOD_80X/checkout_80X.csh
# cd $CMSSW_BASE/src
# cmsenv
# source /tmp/checkout_80X.csh


############## For CMSSW_8_0_21
git cms-init
# Preliminary electron scale and smearing corrections according to https://twiki.cern.ch/twiki/bin/view/CMS/EGMSmearer
git cms-merge-topic -u shervin86:Moriond2017_JEC_energyScales
git cms-addpkg EgammaAnalysis/ElectronTools
(cd EgammaAnalysis/ElectronTools/data ; git clone -b master https://github.com/ECALELFS/ScalesSmearings.git)


#### Please do not add any custom (non-CMSSW) package before this line ####

#ZZAnalysis
git clone https://github.com/CJLST/ZZAnalysis.git ZZAnalysis
(cd ZZAnalysis; git checkout miniAOD_80X)

#effective areas (to be updated)
git clone -n https://github.com/latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools
(cd Muon/MuonAnalysisTools ; git checkout master -- interface/MuonEffectiveArea.h )

#Note we did not eventually adopt the latest EA update, we kept V00-00-30-00 of UserCode/EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h .  c0db796 corresponds to it.
git clone -n https://github.com/cms-analysis/EgammaAnalysis-ElectronTools EGamma/EGammaAnalysisTools
(cd EGamma/EGammaAnalysisTools; git checkout c0db796 -- interface/ElectronEffectiveArea.h) 

#MuScleFit: probably tbf
#git clone https://github.com/scasasso/usercode MuScleFit

#MELA
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
(cd ZZMatrixElement; git checkout -b from-v203 v2.0.3)
# replace ZZMatrixElement/MELA/setup.sh -j 8)
pushd ${CMSSW_BASE}/src/ZZMatrixElement/MELA/fortran/
make all
mv libjhugenmela.so ../data/${SCRAM_ARCH}/
popd

#kinematic refitting
git clone https://github.com/VBF-HZZ/KinZfitter.git
(cd KinZfitter ; git checkout -b from-dd5f616 dd5f616)

#muon momentum scale corrections (76X)
git clone https://github.com/bachtis/Analysis.git -b KaMuCa_V4 KaMuCa 

#Jet energy corrections (CMGTools)
#(mkdir -p CMGTools/Common; cd CMGTools/Common ; wget https://raw.githubusercontent.com/CERN-PH-CMG/cmg-cmssw/a875832047532c5469aa9795751f0363cd5d9244/CMGTools/Common/plugins/JetEnergyCorrector.h)

