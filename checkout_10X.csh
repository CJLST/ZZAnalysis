#!/bin/tcsh -fe
#
# Instructions:
# wget -O ${TMPDIR}/checkout_10X.csh https://raw.githubusercontent.com/CJLST/ZZAnalysis/miniAOD_80X/checkout_10X.csh
# cd $CMSSW_BASE/src
# cmsenv
# chmod u+x ${TMPDIR}/checkout_10X.csh
# ${TMPDIR}/checkout_10X.csh

############## For CMSSW_10_2_5_patch1
git cms-init

#Preliminary electron scale and smearing corrections according to https://twiki.cern.ch/twiki/bin/view/CMS/EGMSmearer
# git cms-merge-topic cms-egamma:EGM_94X_v1
# (cd EgammaAnalysis/ElectronTools/data ; git clone https://github.com/ECALELFS/ScalesSmearings.git ; cd ScalesSmearings ; git checkout Run2017_17Nov2017_v1)

#MET corrections according to https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription#Instructions_for_9_4_X_X_0_for_M
# git cms-merge-topic cms-met:METRecipe94x

#Simplified template cross section
# git cms-addpkg GeneratorInterface/RivetInterface

#### Please do not add any custom (non-CMSSW) package before this line ####

#ZZAnalysis
git clone https://github.com/CJLST/ZZAnalysis.git ZZAnalysis
(cd ZZAnalysis; git checkout miniAOD_80X)

#MELA Analytics
git clone https://github.com/usarica/MelaAnalytics.git
(cd MelaAnalytics; git checkout -b from-v11 v1.1)

#Common LHE tools
git clone https://github.com/usarica/CommonLHETools.git
(cd CommonLHETools; git checkout -b from-v121 v1.2.1)

#MELA
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
(cd ZZMatrixElement; git checkout -b from-v217b2 v2.1.7b2)
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


