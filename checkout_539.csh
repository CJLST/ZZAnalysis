#!/bin/tcsh -fe

############## For V5_15_0/CMSSW_5_3_9 POST-LEGACY

setenv CVSROOT ":ext:${USER}@lxplus.cern.ch:/afs/cern.ch/user/c/cvscmssw/public/CMSSW"

if ( ! -X git ) then 
    scram setup /afs/cern.ch/cms/slc5_amd64_gcc462/external/git-toolfile/1.0/etc/scram.d/git.xml 
    rehash
endif

git clone https://github.com/CJLST/ZZAnalysis.git ZZAnalysis
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement ; (cd ZZMatrixElement; git checkout -b from-V00-00-27 V00-00-27)

# cvs co -r V02-06-00 HiggsAnalysis/CombinedLimit
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit; (cd HiggsAnalysis/CombinedLimit; git pull origin master; git checkout -b from-V02-06-00 HiggsAnalysis-CombinedLimit-V02-06-00)

mkdir MuScleFit; (cd MuScleFit; cvs co -r muscle_v4_2_0 -d Calibration UserCode/scasasso/MuScleFit/Calibration)
#mkdir Higgs; (cd Higgs; cvs co -r V00-03-01 -d Higgs_CS_and_Width UserCode/Snowball/Higgs/Higgs_CS_and_Width)
git clone https://github.com/UFLHEP/HCSaW.git Higgs; (cd Higgs/Higgs_CS_and_Width; git checkout -b from-bf8904f224 bf8904f224 )
mkdir  HZZ4L_Combination; (cd HZZ4L_Combination; cvs co -r bonato_supermela_20121107 -d CombinationPy UserCode/HZZ4L_Combination/CombinationPy)
mkdir -p  Muon/MuonAnalysisTools; (cd  Muon/MuonAnalysisTools; cvs co -r 1.7 -d interface UserCode/sixie/Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h)
mkdir -p CMGTools/Common; (cd CMGTools/Common;  cvs co -r lucieg_Ap18 -d plugins UserCode/CMG/CMGTools/Common/plugins/JetEnergyCorrector.h)

# Egamma official recipe for corrections/smearing
cvs co -r V09-00-01 RecoEgamma/EgammaTools
mkdir -p EGamma/EGammaAnalysisTools; (cd EGamma/EGammaAnalysisTools; cvs co -r V00-00-30-00 -d interface UserCode/EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h)
#cvs co -r V00-00-30-01 -d EGamma/EGammaAnalysisTools UserCode/EGamma/EGammaAnalysisTools # New effective areas
cvs co -r V00-00-08 EgammaAnalysis/ElectronTools 
cd EgammaAnalysis/ElectronTools/data
cat download.url | xargs wget
cd -

mkdir AnalysisDataFormats; (cd AnalysisDataFormats; cvs co -r lucieg_Ap8 -d CMGTools UserCode/CMG/AnalysisDataFormats/CMGTools)
 (cd CMGTools; cvs co -r V00-03-04 -d External UserCode/CMG/CMGTools/External)
 (cd CMGTools; cvs co -r cbern_Apr15 -d RootTools UserCode/CMG/CMGTools/RootTools)
 (cd CMGTools; cvs co -d Production UserCode/CMG/CMGTools/Production)
cvs co -r V02-05-11 DataFormats/CaloRecHit
cvs co -r V06-05-06-10 DataFormats/PatCandidates
cvs co -r V00-02-14 DataFormats/StdDictionaries
cvs co -r V08-09-56 PhysicsTools/PatAlgos
cvs co -r V03-09-28 PhysicsTools/PatUtils
cvs co -r V03-03-12-02 RecoMET/METProducers


if ( $SCRAM_ARCH == "slc6_amd64_gcc472" || $SCRAM_ARCH == "slc6_amd64_gcc481") then
cat <<EOF >! tmp.txt
<Flags CXXFLAGS="-Wno-delete-non-virtual-dtor -Wno-error=unused-but-set-variable -Wno-error=unused-variable"/> 
EOF
cat tmp.txt >> CMGTools/RootTools/BuildFile.xml
cat tmp.txt >> RecoEgamma/EgammaTools/BuildFile.xml
cat tmp.txt >> ZZAnalysis/AnalysisStep/BuildFile.xml

endif
