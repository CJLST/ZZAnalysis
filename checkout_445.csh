#!/bin/tcsh -fe

############## For V5_15_0/CMSSW_4_4_5 POST-LEGACY 

setenv CVSROOT ":ext:${USER}@lxplus5.cern.ch:/afs/cern.ch/user/c/cvscmssw/public/CMSSW"

scram setup /afs/cern.ch/cms/slc5_amd64_gcc462/external/git-toolfile/1.0/etc/scram.d/git.xml
rehash

git clone https://github.com/CJLST/ZZAnalysis.git ZZAnalysis

git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement ; (cd ZZMatrixElement; git checkout -b from-V00-00-26 V00-00-26)

git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit; 
cd HiggsAnalysis/CombinedLimit
git pull origin master
git checkout -b from-V01-13-02 HiggsAnalysis-CombinedLimit-V01-13-02
git checkout a8b620f33381f7bf345210d17ee68a50d8e9acb8 src/HZZ2L2QRooPdfs.cc
#git checkout 9f6fb24b65d8f93eaf9db5f5bee12df3c51783c2 src/HZZ4LRooPdfs.cc 
git checkout c842d93c26ecd7040ace3b5b056d00bf169f926c  src/HZZ4LRooPdfs.cc
git checkout a8b620f33381f7bf345210d17ee68a50d8e9acb8 interface/HZZ2L2QRooPdfs.h
git checkout aaed00c3a568d333fc77ec164980f93a4a1181ef interface/HZZ4LRooPdfs.h
cd -
wget www.cern.ch/amapane/H4l/CMSSW/444/HiggsAnalysis/LinkDef.h; mv LinkDef.h HiggsAnalysis/CombinedLimit/src

#mkdir Higgs; (cd Higgs; cvs co -r V00-03-01 -d Higgs_CS_and_Width UserCode/Snowball/Higgs/Higgs_CS_and_Width)
git clone https://github.com/UFLHEP/HCSaW.git Higgs; (cd Higgs; git checkout -b from-bf8904f224 bf8904f224 ; rm -r Higgs_CS_and_Width_Fermiophobic Higgs_CS_and_Width_SM4)
mkdir MuScleFit; (cd MuScleFit; cvs co -r muscle_v4_2_0 -d Calibration UserCode/scasasso/MuScleFit/Calibration)
mkdir  HZZ4L_Combination; (cd HZZ4L_Combination; cvs co -r bonato_supermela_20121107 -d CombinationPy UserCode/HZZ4L_Combination/CombinationPy)
mkdir -p  Muon/MuonAnalysisTools; (cd  Muon/MuonAnalysisTools; cvs co -r 1.7 -d interface UserCode/sixie/Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h)

#For mu ghost cleaning
cvs co -r U09-01-05-01 DataFormats/MuonReco

# Egamma official recipe for corrections/smearing
cvs co -r V08-11-10-02 RecoEgamma/EgammaTools
cvs co -r V00-00-05-44X EgammaAnalysis/ElectronTools 
cd EgammaAnalysis/ElectronTools/data
cat download.url | xargs wget
cd -
mkdir -p EGamma/EGammaAnalysisTools; (cd EGamma/EGammaAnalysisTools; cvs co -r V00-00-30-BP42X -d interface UserCode/EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h)
cvs co -r V06-23-01      CondFormats/DataRecord
cvs co -r V01-02-13      CondFormats/EcalObjects
cvs co -r V00-04-01      CondFormats/EgammaObjects
cvs co -r V05-08-24      RecoEcal/EgammaCoreTools

mkdir AnalysisDataFormats; (cd AnalysisDataFormats; cvs co -r lucieg_Ap8 -d CMGTools UserCode/CMG/AnalysisDataFormats/CMGTools)
mkdir CMGTools; (cd CMGTools; cvs co -r V00-02-10 -d External UserCode/CMG/CMGTools/External; cvs co -r cbern_Apr15 -d RootTools UserCode/CMG/CMGTools/RootTools; cvs co -d Production UserCode/CMG/CMGTools/Production)
cvs co -r V02-05-11 DataFormats/CaloRecHit
cvs co -r lhx_12JAN2012_v1 DataFormats/METReco
cvs co -r V06-04-39 DataFormats/PatCandidates
cvs co -r V08-07-52 PhysicsTools/PatAlgos
cvs co -r V03-09-18-03 PhysicsTools/PatUtils
cvs co -r V08-03-16 PhysicsTools/Utilities

