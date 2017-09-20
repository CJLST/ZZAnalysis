#!/bin/tcsh -fe
#
# Instructions:
# wget -O /tmp/checkout_80X.csh https://raw.githubusercontent.com/CJLST/ZZAnalysis/miniAOD_80X/checkout_80X.csh
# cd $CMSSW_BASE/src
# cmsenv
# source /tmp/checkout_90X.csh


############## For CMSSW_9_2_10
git cms-init
# Preliminary electron scale and smearing corrections according to https://twiki.cern.ch/twiki/bin/view/CMS/EGMSmearer
#FIXME this includes some changes to EgammaAnalysis/ElectronTools; in particular, regressionWeights_cfi that is not present in cmssw proper (?)
#git cms-merge-topic rafaellopesdesa:EgammaAnalysis80_EGMSmearer_Moriond17_23Jan

git cms-addpkg EgammaAnalysis/ElectronTools
#FIXME the following scale/smearing are for Moriond 2017, to be updated  
#(cd EgammaAnalysis/ElectronTools/data ; git clone -b master https://github.com/ECALELFS/ScalesSmearings.git ; cd ScalesSmearings ; git checkout tags/Moriond17_23Jan_v1)

#MET recipe according to https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription#Instructions_for_8_0_X_X_26_patc
#FIXME git cms-merge-topic cms-met:METRecipe_8020 -u
#FIXME git cms-merge-topic cms-met:METRecipe_80X_part2 -u


#### Please do not add any custom (non-CMSSW) package before this line ####

#ZZAnalysis
git clone https://github.com/CJLST/ZZAnalysis.git ZZAnalysis
(cd ZZAnalysis; git checkout miniAOD_80X)

#effective areas (to be updated)
#FIXME git clone -n https://github.com/latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools
#FIXME (cd Muon/MuonAnalysisTools ; git checkout master -- interface/MuonEffectiveArea.h )

#Note we did not eventually adopt the latest EA update, we kept V00-00-30-00 of UserCode/EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h .  c0db796 corresponds to it.
#FIXME git clone -n https://github.com/cms-analysis/EgammaAnalysis-ElectronTools EGamma/EGammaAnalysisTools
#FIXME (cd EGamma/EGammaAnalysisTools; git checkout c0db796 -- interface/ElectronEffectiveArea.h) 

#MuScleFit: probably tbf
#git clone https://github.com/scasasso/usercode MuScleFit

# Higgs Combination Package, Needed for the Double Crystall Ball function. 
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
(cd HiggsAnalysis/CombinedLimit; git checkout 81x-root606-integration)

#MELA
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
(cd ZZMatrixElement; git checkout -b from-v210p1 v2.1.0.p1)
# replace ZZMatrixElement/MELA/setup.sh -j 8)
pushd ${CMSSW_BASE}/src/ZZMatrixElement/MELA/fortran/
make all
mv libjhugenmela.so ../data/${SCRAM_ARCH}/
popd

#kinematic refitting
git clone https://github.com/mhl0116/KinZfitter-1.git KinZfitter
(cd KinZfitter ; git checkout -b from-ee2d8ef ee2d8ef)

#muon momentum scale corrections (76X)
git clone https://github.com/bachtis/Analysis.git -b KaMuCa_V4 KaMuCa 


