#!/bin/sh
#
# Instructions:
# wget -P /tmp https://github.com/hengne/ZZAnalysis/blob/miniAOD/checkout_70X.csh
# cd $CMSSW_BASE/src
# cmsenv
# source /tmp/checkout_70X.sh


############## test version for CMSSW_7_3_3

#init a new release
release=CMSSW_7_3_3
export SCRAM_ARCH=slc6_amd64_gcc481
alias cmsenv='eval `scramv1 runtime -sh`'
alias cmsrel='scramv1 project CMSSW'

scram project $release
cd $release/src
cmsenv


#electron MVA ID
#git cms-merge-topic HuguesBrun:trigElecIdInCommonIsoSelection720 

#ZZAnalysis
git clone https://github.com/CJLST/ZZAnalysis.git ZZAnalysis

cd ZZAnalysis
git checkout from-miniAOD_hengne_dev_733
git checkout 7c06aa1 -- AnalysisStep/plugins/EleFiller.cc
cd -

#effective areas (to be updated)
git clone -n https://github.com/latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools
cd Muon/MuonAnalysisTools 
git checkout master -- interface/MuonEffectiveArea.h 
cd -

#Note we did not eventually adopt the latest EA update, we kept V00-00-30-00 of UserCode/EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h .  c0db796 corresponds to it.
git clone -n https://github.com/cms-analysis/EgammaAnalysis-ElectronTools EGamma/EGammaAnalysisTools
cd EGamma/EGammaAnalysisTools
git checkout c0db796 -- interface/ElectronEffectiveArea.h
cd -


#MuScleFit: probably tbf
#git clone https://github.com/scasasso/usercode MuScleFit

#MELA
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
cd ZZMatrixElement
git checkout -b from-d375a59 d375a59
cd -

#MELA dependencies
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit

#photon ISO information for FSR 
git clone -n https://github.com/VBF-HZZ/UFHZZAnalysisRun2
cd UFHZZAnalysisRun2 
git checkout origin/csa14 FSRPhotons #This does not set the correct branch, but picks the right one anyway
cd -

#q/g tagging (here dirtily backported from 74X, to be simplified in the future)
git cms-addpkg RecoJets/JetProducers
git cms-addpkg RecoJets/JetAlgorithms
wget https://raw.githubusercontent.com/cms-sw/cmssw/CMSSW_7_4_X/RecoJets/JetProducers/interface/QGTagger.h -O RecoJets/JetProducers/interface/QGTagger.h
wget https://raw.githubusercontent.com/cms-sw/cmssw/CMSSW_7_4_X/RecoJets/JetProducers/plugins/QGTagger.cc -O RecoJets/JetProducers/plugins/QGTagger.cc
wget https://raw.githubusercontent.com/cms-sw/cmssw/CMSSW_7_4_X/RecoJets/JetProducers/plugins/BuildFile.xml -O RecoJets/JetProducers/plugins/BuildFile.xml
wget https://raw.githubusercontent.com/cms-sw/cmssw/CMSSW_7_4_X/RecoJets/JetProducers/python/QGTagger_cfi.py -O RecoJets/JetProducers/python/QGTagger_cfi.py
wget https://raw.githubusercontent.com/cms-sw/cmssw/CMSSW_7_4_X/RecoJets/JetAlgorithms/src/QGLikelihoodCalculator.cc -O RecoJets/JetAlgorithms/src/QGLikelihoodCalculator.cc

#Jet energy corrections (CMGTools)
#(mkdir -p CMGTools/Common; cd CMGTools/Common ; wget https://raw.githubusercontent.com/CERN-PH-CMG/cmg-cmssw/a875832047532c5469aa9795751f0363cd5d9244/CMGTools/Common/plugins/JetEnergyCorrector.h)

#Not needed anymore (dependency removed)
#git clone https://github.com/HZZ4l/CombinationPy.git HZZ4L_Combination/CombinationPy

# Not needed, for the time being
#git clone https://github.com/msnowball/HCSaW Higgs/Higgs_CS_and_Width
#cd Higgs/Higgs_CS_and_Width
#git filter-branch --subdirectory-filter Higgs_CS_and_Width
#cd -
