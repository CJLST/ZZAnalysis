#!/bin/tcsh -fe
#
# Instructions:
# wget -O ${TMPDIR}/checkout_10X.csh https://raw.githubusercontent.com/acappati/ZZAnalysis/Run2UL_22_nano/checkout_10X.csh
# cd $CMSSW_BASE/src
# cmsenv
# chmod u+x ${TMPDIR}/checkout_10X.csh
# ${TMPDIR}/checkout_10X.csh

############## For CMSSW_10_6_26
git cms-init

# New Jet PU ID: dedicated training for each year
git cms-addpkg  RecoJets/JetProducers

# STXS Categorisation: now directly implemented in CMSSW
git cms-addpkg GeneratorInterface/RivetInterface
git cms-addpkg SimDataFormats/HTXS

# Updated for UL. See: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018
git cms-addpkg RecoEgamma/EgammaTools
git clone https://github.com/cms-egamma/EgammaPostRecoTools.git
mv EgammaPostRecoTools/python/EgammaPostRecoTools.py RecoEgamma/EgammaTools/python/.
git cms-addpkg EgammaAnalysis/ElectronTools
(rm -rf EgammaAnalysis/ElectronTools/data;git clone https://github.com/jainshilpi/EgammaAnalysis-ElectronTools.git -b ULSSfiles_correctScaleSysMC EgammaAnalysis/ElectronTools/data;)

#### Please do not add any custom (non-CMSSW) package before this line ####
#ZZAnalysis
git clone https://github.com/acappati/ZZAnalysis.git ZZAnalysis
(cd ZZAnalysis; git checkout Run2UL_22_nano)

# Muon MVA
git cms-addpkg CondFormats/EgammaObjects
git cms-addpkg CommonTools/MVAUtils
git clone https://github.com/bonanomi/MuonMVAReader.git MuonMVAReader
(cd MuonMVAReader; git checkout 3d53269)

#MELA Analytics
git clone https://github.com/MELALabs/MelaAnalytics.git
(cd MelaAnalytics; git checkout -b from-v22 v2.2)

#Common LHE tools
git clone https://github.com/usarica/CommonLHETools.git
(cd CommonLHETools; git checkout -b from-v141 v1.4.1)

#MELA
git clone https://github.com/JHUGen/JHUGenMELA.git JHUGenMELA
(cd JHUGenMELA; git checkout -b from-v235 v2.3.5)

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

#NanoAODTools 
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git NanoAODTools
#We may have to switch to our fork at some point since the central version is unmaintained
#git clone https://github.com/CJLST/nanoAOD-tools.git PhysicsTools/NanoAODTools
#(cd PhysicsTools/NanoAODTools ; git checkout -b from-ded19d8 ded19d8)
