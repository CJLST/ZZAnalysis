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
git cms-merge-topic cms-egamma:EgammaID_1023
git cms-merge-topic cms-egamma:EgammaPostRecoTools
git clone git@github.com:cms-egamma/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data
(cd EgammaAnalysis/ElectronTools/data ; git checkout ScalesSmearing2018_Dev)
git cms-merge-topic cms-egamma:EgammaPostRecoTools_dev

#MET corrections according to https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription#Instructions_for_9_4_X_X_0_for_M
# git cms-merge-topic cms-met:METRecipe94x

#Simplified template cross section
git cms-addpkg GeneratorInterface/RivetInterface
git cms-addpkg SimDataFormats/HTXS
git remote add amarini https://github.com/amarini/cmssw.git
git fetch amarini 
git checkout amarini/htxs_stage1p1_cmssw949_v2 GeneratorInterface/RivetInterface
git checkout amarini/htxs_stage1p1_cmssw949_v2 SimDataFormats/HTXS

#### Please do not add any custom (non-CMSSW) package before this line ####
#ZZAnalysis
git clone https://github.com/CJLST/ZZAnalysis.git ZZAnalysis
(cd ZZAnalysis; git checkout miniAOD_80X)

#MELA Analytics
git clone https://github.com/usarica/MelaAnalytics.git
(cd MelaAnalytics; git checkout -b from-v11 v1.1)

#Common LHE tools
git clone https://github.com/usarica/CommonLHETools.git
(cd CommonLHETools; git checkout -b from-v124 v1.2.4)

#MELA
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
(cd ZZMatrixElement; git checkout -b from-v218 v2.1.8)
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


