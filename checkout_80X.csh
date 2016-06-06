#!/bin/tcsh -fe
#
# Instructions:
# wget -O /tmp/checkout_80X.csh https://raw.githubusercontent.com/CJLST/ZZAnalysis/miniAOD_80X/checkout_80X.csh
# cd $CMSSW_BASE/src
# cmsenv
# source /tmp/checkout_80X.csh


############## For CMSSW_8_0_6

#electron momentum scale corrections (76X).
git cms-merge-topic -u matteosan1:smearer_76X

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
(cd ZZMatrixElement ; git checkout -b from-c6d45de c6d45de)

#photon ISO information for FSR 
git clone -n https://github.com/VBF-HZZ/UFHZZAnalysisRun2
(cd UFHZZAnalysisRun2 ; git checkout origin/csa14 FSRPhotons) #This does not set the correct branch, but picks the right one anyway

#kinematic refitting
git clone https://github.com/tocheng/KinZfitter.git
(cd KinZfitter ; git checkout -b from-v1.0 v1.0)

#muon momentum scale corrections (76X)
git clone https://github.com/bachtis/Analysis.git -b KaMuCa_V2 KaMuCa 

#hack the KinZfitter to use the corrected muon pT error
sed -i 's/reco::Muon/pat::Muon/g' KinZfitter/HelperFunction/interface/HelperFunction.h
sed -i 's/reco::Muon/pat::Muon/g' KinZfitter/HelperFunction/src/HelperFunction.cc
sed -i 's/double pterr = mu->muonBestTrack()->ptError();/double pterr = mu->userFloat("correctedPtError");/g' KinZfitter/HelperFunction/src/HelperFunction.cc

#Jet energy corrections (CMGTools)
#(mkdir -p CMGTools/Common; cd CMGTools/Common ; wget https://raw.githubusercontent.com/CERN-PH-CMG/cmg-cmssw/a875832047532c5469aa9795751f0363cd5d9244/CMGTools/Common/plugins/JetEnergyCorrector.h)

# Not needed, for the time being
#git clone https://github.com/msnowball/HCSaW Higgs/Higgs_CS_and_Width
#cd Higgs/Higgs_CS_and_Width
#git filter-branch --subdirectory-filter Higgs_CS_and_Width
#cd -
