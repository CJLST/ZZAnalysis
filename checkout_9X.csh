#!/bin/tcsh -fe
#
# Instructions:
# wget -O ${TMPDIR}/checkout_9X.csh https://raw.githubusercontent.com/CJLST/ZZAnalysis/miniAOD_80X/checkout_9X.csh
# cd $CMSSW_BASE/src
# cmsenv
# source ${TMPDIR}/checkout_9X.csh


############## For CMSSW_9_4_2
git cms-init

#Preliminary electron scale and smearing corrections according to https://twiki.cern.ch/twiki/bin/view/CMS/EGMSmearer
git cms-merge-topic cms-egamma:EGM_94X_v1
(cd EgammaAnalysis/ElectronTools/data ; git clone https://github.com/ECALELFS/ScalesSmearings.git ; cd ScalesSmearings ; git checkout Run2017_17Nov2017_v1)

#Electron MVA ID in 94X according to https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariateElectronIdentificationRun2#Recipes_for_regular_users
# MVA ID V2 now, not yet available as part of official recepie
git cms-merge-topic guitargeek:ElectronID_MVA2017_V2_HZZ_940pre3
rm -rf RecoEgamma/ElectronIdentification/data #Delete old BDT weights so we can clone new ones
git clone -b ElectronID_MVA2017_V2 https://github.com/guitargeek/RecoEgamma-ElectronIdentification RecoEgamma/ElectronIdentification/data/
#Hack to make our run_cfg.py job script work with new implementation of ID
sed -i "s@-float('Inf')@-999999.@g" RecoEgamma/ElectronIdentification/python/Identification/mvaElectronID_Fall17_iso_V1_cff.py
sed -i "s@float('Inf')@999999.@g" RecoEgamma/ElectronIdentification/python/Identification/mvaElectronID_Fall17_iso_V1_cff.py
sed -i "s@-float('Inf')@-999999.@g" RecoEgamma/ElectronIdentification/python/Identification/mvaElectronID_Fall17_iso_V2_cff.py
sed -i "s@float('Inf')@999999.@g" RecoEgamma/ElectronIdentification/python/Identification/mvaElectronID_Fall17_iso_V2_cff.py
sed -i "s@-float('Inf')@-999999.@g" RecoEgamma/ElectronIdentification/python/Identification/mvaElectronID_Fall17_noIso_V1_cff.py
sed -i "s@float('Inf')@999999.@g" RecoEgamma/ElectronIdentification/python/Identification/mvaElectronID_Fall17_noIso_V1_cff.py
sed -i "s@-float('Inf')@-999999.@g" RecoEgamma/ElectronIdentification/python/Identification/mvaElectronID_Fall17_noIso_V2_cff.py
sed -i "s@float('Inf')@999999.@g" RecoEgamma/ElectronIdentification/python/Identification/mvaElectronID_Fall17_noIso_V2_cff.py

#Simplified template cross section
git cms-addpkg GeneratorInterface/RivetInterface
#The tool is currently broken so it needs this hack to work
sed -i 's@#include "FWCore/Framework/interface/stream/EDProducer.h"@#include "FWCore/Framework/interface/EDProducer.h"@g' GeneratorInterface/RivetInterface/plugins/HTXSRivetProducer.cc
sed -i 's@class HTXSRivetProducer : public edm::stream::EDProducer<> {@class HTXSRivetProducer : public edm::EDProducer {@g' GeneratorInterface/RivetInterface/plugins/HTXSRivetProducer.cc

#### Please do not add any custom (non-CMSSW) package before this line ####

#ZZAnalysis
git clone https://github.com/CJLST/ZZAnalysis.git ZZAnalysis
(cd ZZAnalysis; git checkout miniAOD_80X)

#MuScleFit: probably tbf
#git clone https://github.com/scasasso/usercode MuScleFit

# Higgs Combination Package, Needed for the Double Crystall Ball function. 
git clone -n https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
(cd HiggsAnalysis/CombinedLimit; git checkout origin/81x-root606-integration -- interface/HZZ2L2QRooPdfs.h interface/HZZ4LRooPdfs.h src/HZZ2L2QRooPdfs.cc src/HZZ4LRooPdfs.cc BuildFile.xml)
cat <<EOF > HiggsAnalysis/CombinedLimit/src/classes.h
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"
EOF
cat <<EOF > HiggsAnalysis/CombinedLimit/src/classes_def.xml
<lcgdict>
        <class name="RooSigPlusInt" />
        <class name="RooVBFZZPdf" />
        <class name="RooVBFZZPdf_v2" />
        <class name="RooaDoubleCBxBW" />
        <class name="RooggZZPdf" />
        <class name="RooggZZPdf_v2" />
        <class name="RooqqZZPdf" />
        <class name="RooqqZZPdf_v2" />
        <class name="RooDoubleCB" />
        <class name="RooCPSHighMassGGH" />
        <class name="RooCPSHighMassGGHNoInterf" />
        <class name="RooCPSHighMassVBF" />
        <class name="RooCPSHighMassVBFNoInterf" />
        <class name="RooBWHighMassGGH" />
        <class name="RooTsallis" />
        <class name="RooRelBW" />
        <class name="RooRelBW1" />
        <class name="RooRelBWHighMass" />
        <class name="RooRelBWUF" />
        <class name="RooRelBWUFParam" />
        <class name="RooRelBWUFParamWidth" />
        <class name="RooRelBWUF_SM4" />
        <class name="RooTwoETwoMuMassRes" />
        <class name="RooTwoETwoMuMassShapePdf2" />
        <class name="RooFourEMassRes" />
        <class name="RooFourEMassShapePdf2" />
        <class name="RooFourMuMassRes" />
        <class name="RooFourMuMassShapePdf2" />
        <class name="Roo4lMasses2D" />
        <class name="Roo4lMasses2D_Bkg" />
        <class name="Roo4lMasses2D_BkgGGZZ" />
        <class name="RooBetaFunc_v2" />
        <class name="RooLevelledExp" />
        <class name="Triangle" />
        <class name="RooFermi" />
        <class name="RooCB" />
</lcgdict>
EOF


#MELA Analytics
git clone https://github.com/usarica/MelaAnalytics.git

#MELA
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
(cd ZZMatrixElement; git checkout -b from-v214 v2.1.4)
# replace ZZMatrixElement/MELA/setup.sh -j 8
pushd ${CMSSW_BASE}/src/ZZMatrixElement/MELA/COLLIER/
  pkgname="collier-1.2"
  pkgdir="COLLIER-1.2"
  tarname=$pkgname".tar.gz"
  tarweb="https://www.hepforge.org/archive/collier/"$tarname
  libname="libcollier.so"
  tmpdir="colliertmp"
  wget $tarweb
  mkdir $tmpdir
  tar -xvzf $tarname -C $tmpdir
  rm $tarname
  mv $tmpdir"/"$pkgdir"/src/"* ./
  rm -rf $tmpdir
  make
  mv $libname "../data/"$SCRAM_ARCH"/"$libname
popd
pushd ${CMSSW_BASE}/src/ZZMatrixElement/MELA/fortran/
  make all
  mv libjhugenmela.so ../data/${SCRAM_ARCH}/
popd

#kinematic refitting
git clone https://github.com/mhl0116/KinZfitter-1.git KinZfitter
(cd KinZfitter ; git checkout -b from-ee2d8ef ee2d8ef)

#muon momentum scale corrections (76X)
git clone https://github.com/bachtis/Analysis.git -b KaMuCa_V4 KaMuCa 


