#!/bin/tcsh -fe
#
# Instructions:
# wget -P /tmp https://raw.githubusercontent.com/CJLST/ZZAnalysis/miniAOD_74X/checkout_70X.csh
# cd $CMSSW_BASE/src
# cmsenv
# source /tmp/checkout_70X.csh


############## For miniAOD/CMSSW_7_4_5

# init git env.
git cms-init

#electron MVA ID (still with Phys14 training)
#git cms-merge-topic sregnard:Phys14ElectronMvaIdFor745

# put EID and PUPPI together
git cms-merge-topic VirginiaCMS:merge-puppi-Phys14ElectronMvaIdFor745

#ZZAnalysis
git clone https://github.com/CJLST/ZZAnalysis.git ZZAnalysis
(cd ZZAnalysis; git checkout miniAOD_74X_puppi)

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
(cd ZZMatrixElement ; git checkout -b from-V00-02-00 V00-02-00)
(cd ZZMatrixElement/MELA/data; mkdir slc6_amd64_gcc491; cp slc6_amd64_gcc481/* slc6_amd64_gcc491/)

#MELA dependencies
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
(cd HiggsAnalysis/CombinedLimit; git checkout 74x-root6)

#photon ISO information for FSR 
git clone -n https://github.com/VBF-HZZ/UFHZZAnalysisRun2
(cd UFHZZAnalysisRun2 ; git checkout origin/csa14 FSRPhotons) #This does not set the correct branch, but picks the right one anyway


#Jet energy corrections (CMGTools)
#(mkdir -p CMGTools/Common; cd CMGTools/Common ; wget https://raw.githubusercontent.com/CERN-PH-CMG/cmg-cmssw/a875832047532c5469aa9795751f0363cd5d9244/CMGTools/Common/plugins/JetEnergyCorrector.h)

#Not needed anymore (dependency removed)
#git clone https://github.com/HZZ4l/CombinationPy.git HZZ4L_Combination/CombinationPy

# Not needed, for the time being
#git clone https://github.com/msnowball/HCSaW Higgs/Higgs_CS_and_Width
#cd Higgs/Higgs_CS_and_Width
#git filter-branch --subdirectory-filter Higgs_CS_and_Width
#cd -
