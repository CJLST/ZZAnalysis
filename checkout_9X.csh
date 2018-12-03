#!/bin/tcsh -fe
#
# Instructions:
# wget -O ${TMPDIR}/checkout_9X.csh https://raw.githubusercontent.com/CJLST/ZZAnalysis/miniAOD_80X/checkout_9X.csh
# cd $CMSSW_BASE/src
# cmsenv
# chmod u+x ${TMPDIR}/checkout_9X.csh
# ${TMPDIR}/checkout_9X.csh

############## For CMSSW_9_4_9
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


#MET corrections according to https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription#Instructions_for_9_4_X_X_9_for_2
git cms-merge-topic cms-met:METFixEE2017_949_v2

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
git clone -b 94x https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit

#MELA Analytics
git clone https://github.com/usarica/MelaAnalytics.git
(cd MelaAnalytics; git checkout -b from-v11 v1.1)

#Common LHE tools
git clone https://github.com/usarica/CommonLHETools.git
(cd CommonLHETools; git checkout -b from-v121 v1.2.1)

#MELA
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
(cd ZZMatrixElement; git checkout -b from-v216 v2.1.6)
# replace ZZMatrixElement/MELA/setup.sh -j 8
(                                                                 \
  cd ${CMSSW_BASE}/src/ZZMatrixElement/MELA/COLLIER/             ;\
  set pkgname="collier-1.2"                                      ;\
  set pkgdir="COLLIER-1.2"                                       ;\
  set tarname=$pkgname".tar.gz"                                  ;\
  set tarweb="https://www.hepforge.org/archive/collier/"$tarname ;\
  set libname="libcollier.so"                                    ;\
  set tmpdir="colliertmp"                                        ;\
  wget $tarweb                                                   ;\
  mkdir $tmpdir                                                  ;\
  tar -xvzf $tarname -C $tmpdir                                  ;\
  rm $tarname                                                    ;\
  mv $tmpdir"/"$pkgdir"/src/"* ./                                ;\
  rm -rf $tmpdir                                                 ;\
  make                                                           ;\
  mv $libname "../data/"$SCRAM_ARCH"/"$libname                   ;\
)
(                                                                 \
  cd ${CMSSW_BASE}/src/ZZMatrixElement/MELA/fortran/             ;\
  make all                                                       ;\
  mv libjhugenmela.so ../data/${SCRAM_ARCH}/                     ;\
)

#kinematic refitting
git clone https://github.com/mhl0116/KinZfitter-1.git KinZfitter
(cd KinZfitter ; git checkout -b from-27daebb 27daebb)

#muon momentum scale corrections (76X)
git clone https://github.com/bachtis/Analysis.git -b KaMuCa_V4 KaMuCa


