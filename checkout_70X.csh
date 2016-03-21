#!/bin/tcsh -fe
#
# Instructions:
# wget -pP /tmp https://raw.githubusercontent.com/CJLST/ZZAnalysis/miniAOD_76X/checkout_70X.csh
# cd $CMSSW_BASE/src
# cmsenv
# source /tmp/checkout_70X.csh


############## For CMSSW_7_6_3_patch2

#electron momentum scale corrections (76X).
git cms-merge-topic -u matteosan1:smearer_76X

#ZZAnalysis
git clone -b 2l2q https://github.com/CJLST/ZZAnalysis.git ZZAnalysis
(cd ZZAnalysis; git checkout 2l2q)

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
(cd ZZMatrixElement ; git checkout -b from-V00-02-01-patch1 V00-02-01-patch1)

#temporary patch for MELA
cp /afs/cern.ch/user/c/covarell/public/forUlash/MEMCalculators.h ZZMatrixElement/MEMCalculators/interface/
cp /afs/cern.ch/user/c/covarell/public/forUlash/MEMCalculators.cpp ZZMatrixElement/MEMCalculators/src/

#MELA dependencies
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
(cd HiggsAnalysis/CombinedLimit; git checkout 74x-root6)

#photon ISO information for FSR 
git clone -n https://github.com/VBF-HZZ/UFHZZAnalysisRun2
(cd UFHZZAnalysisRun2 ; git checkout origin/csa14 FSRPhotons) #This does not set the correct branch, but picks the right one anyway

#kinematic refitting
git clone https://github.com/covarell/KinZfitter.git
(cd KinZfitter ; git checkout -b from-v1.0 v1.0)

#muon momentum scale corrections (76X)
git clone https://github.com/bachtis/Analysis.git -b KaMuCa_V2 KaMuCa 

#hack the KinZfitter to use the corrected muon pT error
sed -i 's/reco::Muon/pat::Muon/g' KinZfitter/HelperFunction/interface/HelperFunction.h
sed -i 's/reco::Muon/pat::Muon/g' KinZfitter/HelperFunction/src/HelperFunction.cc
sed -i 's/double pterr = mu->muonBestTrack()->ptError();/double pterr = mu->userFloat("correctedPtError");/g' KinZfitter/HelperFunction/src/HelperFunction.cc


#Jet energy corrections (CMGTools)
#(mkdir -p CMGTools/Common; cd CMGTools/Common ; wget https://raw.githubusercontent.com/CERN-PH-CMG/cmg-cmssw/a875832047532c5469aa9795751f0363cd5d9244/CMGTools/Common/plugins/JetEnergyCorrector.h)

#Not needed anymore (dependency removed)
#git clone https://github.com/HZZ4l/CombinationPy.git HZZ4L_Combination/CombinationPy

# Not needed, for the time being
#git clone https://github.com/msnowball/HCSaW Higgs/Higgs_CS_and_Width
#cd Higgs/Higgs_CS_and_Width
#git filter-branch --subdirectory-filter Higgs_CS_and_Width
#cd -
