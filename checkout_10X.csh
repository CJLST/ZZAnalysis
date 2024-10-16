#!/bin/tcsh -fe
#
# Instructions:
# cmssw-el7
# cmsrel CMSSW_10_6_30
# cd CMSSW_10_6_30/src
# cmsenv
# wget -O ${TMPDIR}/checkout_10X.csh https://raw.githubusercontent.com/CJLST/ZZAnalysis/Run2UL_22/checkout.csh
# chmod u+x ${TMPDIR}/checkout.csh
# ${TMPDIR}/checkout.csh

############## For CMSSW_10_6_30
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
git clone https://github.com/CJLST/ZZAnalysis.git ZZAnalysis
(cd ZZAnalysis; git checkout Run2UL_22)

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
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
(cd PhysicsTools/NanoAODTools ; git checkout -b from-c32f055 c32f055)

echo
echo "***Note: Please use ZZAnalysis/start_el7.sh to open a slc7 container in new shells, in order to have access to the Condor queues.***"
