#!/bin/tcsh -fe
#
# Instructions:
# wget -O ${TMPDIR}/checkout_10X.csh https://raw.githubusercontent.com/CJLST/ZZAnalysis/miniAOD_80X/checkout_10X.csh
# cd $CMSSW_BASE/src
# cmsenv
# chmod u+x ${TMPDIR}/checkout_10X.csh
# ${TMPDIR}/checkout_10X.csh

############## For CMSSW_10_2_18
git cms-init

#Preliminary electron scale and smearing corrections according to https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPostRecoRecipes#2018_Preliminary_Energy_Correcti
#We need the ElectronTools package to calculate smear and scale uncertainties so just download the ScaleAndSmearing files manualy 
git cms-merge-topic cms-egamma:EgammaPostRecoTools
git cms-merge-topic cms-egamma:PhotonIDValueMapSpeedup1029
git cms-merge-topic cms-egamma:slava77-btvDictFix_10210
git cms-addpkg EgammaAnalysis/ElectronTools
(rm -rf EgammaAnalysis/ElectronTools/data;git clone https://github.com/cms-data/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data;)

# New Jet PU ID: dedicated training for each year
git cms-addpkg  RecoJets/JetProducers

# STXS Categorisation: now directly implemented in CMSSW
git cms-addpkg GeneratorInterface/RivetInterface
git cms-addpkg SimDataFormats/HTXS

# 2016 and 2018 retraining for electron BDT
git cms-merge-topic mkovac:Electron_XGBoost_MVA_2016_and_2018_CMSSW_10_2_15

#MET corrections according to https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription#Instructions_for_9_4_X_X_0_for_M
git cms-merge-topic cms-met:METFixEE2017_949_v2_backport_to_102X

#Simplified template cross section
## need at least the first two lines
# git cms-addpkg GeneratorInterface/RivetInterface
# git cms-addpkg SimDataFormats/HTXS
# git remote add bonanomi https://github.com/bonanomi/cmssw.git
# git fetch bonanomi 
# git checkout bonanomi/hstxs1p2_CMSSW_10_2_X GeneratorInterface/RivetInterface
# git checkout bonanomi/hstxs1p2_CMSSW_10_2_X SimDataFormats/HTXS

#### Please do not add any custom (non-CMSSW) package before this line ####
#ZZAnalysis
git clone https://github.com/CJLST/ZZAnalysis.git ZZAnalysis
(cd ZZAnalysis; git checkout Run2_CutBased_BTag16)

# Muon MVA
git clone https://github.com/bonanomi/MuonMVAReader.git MuonMVAReader

#MELA Analytics
git clone https://github.com/MELALabs/MelaAnalytics.git

#Common LHE tools
#git clone https://github.com/usarica/CommonLHETools.git
#(cd CommonLHETools; git checkout -b from-v134 v1.3.4)
#--->Hack for scram update preventing shell scripts in BuildFiles
git clone https://github.com/namapane/CommonLHETools.git
(cd CommonLHETools; git checkout v134-fixScram)


#MELA
#git clone https://github.com/JHUGen/JHUGenMELA.git JHUGenMELA
#(cd JHUGenMELA; git checkout -b from-v231 v2.3.1)
#--->Hack for scram update preventing shell scripts in BuildFiles
git clone https://github.com/namapane/JHUGenMELA.git JHUGenMELA
(cd JHUGenMELA; git checkout v231-fixScram)

# replace JHUGenMELA/MELA/setup.sh -j 8
(                                                                 \
  cd ${CMSSW_BASE}/src/JHUGenMELA/MELA/COLLIER/                  ;\
  set pkgname="collier-1.2.0"                                    ;\
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
  cd ${CMSSW_BASE}/src/JHUGenMELA/MELA/fortran/                  ;\
  make all                                                       ;\
  mv libjhugenmela.so ../data/${SCRAM_ARCH}/                     ;\
)
(                                                                 \
  cd ${CMSSW_BASE}/src/JHUGenMELA/MELA/                          ;\
  ./downloadNNPDF.sh                                             ;\
)

#download MCFM lib (cannot be done in BuildFile.xml any longer)
$CMSSW_BASE/src/JHUGenMELA/MELA/data/retrieve.csh $SCRAM_ARCH mcfm_707

#kinematic refitting
git clone https://github.com/mhl0116/KinZfitter-1.git KinZfitter
(cd KinZfitter ; git checkout -b from-27daebb 27daebb)




