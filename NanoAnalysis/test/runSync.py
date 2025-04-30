#!/usr/bin/env python3
'''Script to prepare sync files

Run with: 
./runSync.py
'''

from __future__ import print_function
from ZZAnalysis.NanoAnalysis.tools import setConf, getConf, insertAfter

# Check that the checkout recipe has been properly updated 
from ZZAnalysis.AnalysisStep.validateCheckout import validateCheckout 
if not validateCheckout() :
    exit(1)

#SampleToRun = "MCsync_Rereco"
#SampleToRun = "MCsync_UL"
#SampleToRun = "Data2022"
SampleToRun = "MC2022"
#SampleToRun = "ggh125_UL"

### Customize processing variables
#setConf("DEBUG", True)
setConf("runMELA", True)
setConf("bestCandByMELA", False)
setConf("APPLYMUCORR", False)
setConf("APPLYELECORR", False)

#setConf("SYNCMODE", True) # Force muon resolution correction with fixed +1 sigma smearing
setConf("APPLYJETCORR", False)

json = None #replace this if needed

################################################################################
if SampleToRun == "MC2022" :
    # 2022 MC sample
    setConf("SAMPLENAME", "sync")
    setConf("DATA_TAG", "post_EE")
    setConf("XSEC", 1.)
    setConf("LEPTON_SETUP", 2022)
    setConf("IsMC", True)
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("fileNames",[
        "/store/mc/Run3Summer22EENanoAODv12/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg2-JHUGenV752-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/2540000/25c8f5ff-9de0-4a0c-9e2f-757332ad392f.root", #ggH 13.2 kevts
        "/store/mc/Run3Summer22EENanoAODv12/ZHto2Zto4L_M125_TuneCP5_13p6TeV_powheg2-minlo-HZJ-JHUGenV752-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/2520000/60da2336-6355-438e-ad2b-c8c4d83e50fe.root" #ZH 32.6 kevts 
        ])



#####################################################################
### This import should be done AFTER all customizations (setConf calls)
from ZZAnalysis.NanoAnalysis.nanoZZ4lAnalysis import *
######################################################################

### Tweak postprocessor parameters as necessary
p.prefetch=True # Prefetch remote files
p.longTermCache=True # keep prefetched files (useful for rerunning tests several times)

p.haddFileName = "ZZ4lSync.root"

### Select specific events to debug
#p.cut = "event==1583673"

### Print out detailed candidate information for debug purposes
#from ZZAnalysis.NanoAnalysis.dumpEvents import dumpEvents
#p.cut = None # Remove preselction
#insertAfter(p.modules,"lepFiller",dumpEvents(level=-1),getConf("NANOVERSION", 11)) 

#p.branchsel=None #Read all branches
#p.outputbranchsel=None #Output all branches
from PhysicsTools.NanoAODTools.postprocessing.framework.branchselection import BranchSelection
p.outputbranchsel=BranchSelection(['drop *',
                                   'keep run',
                                   'keep event',
                                   'keep luminosityBlock',
#                  'keep Flag*',
                                   'keep Electron_pt',
                                   'keep Electron_uncorrected_pt',
                                   'keep Electron_eta',
                                   'keep Electron_phi',
                                   'keep Electron_charge',
                                   'keep Electron_pdgId',
                                   'keep Electron_fsrPhotonIdx',
                                   'keep Electron_sip3d',
                                   'keep Electron_mvaHZZIso',
                                   'keep Muon_pt',
                                   'keep Muon_uncorrected_pt',
                                   'keep Muon_eta',
                                   'keep Muon_phi',
                                   'keep Muon_charge',
                                   'keep Muon_pdgId',
                                   'keep Muon_fsrPhotonIdx',
                                   'keep Muon_sip3d',
#                                   'keep Lepton*',
#                                   'drop Lepton_ZZ*',
                                   'keep FsrPhoton*',
                                   'drop FsrPhoton_mass',
                                   'drop FsrPhoton_genFsrIdx',
                                   'keep ZZCand*',
                                   'drop ZZCand_rapidity',
                                   'keep puWeight',
                                   'drop ZZCand_nExtra*'
                                   ])


#replace JSON
p.json = json

#from ZZAnalysis.NanoAnalysis.dumpEvents import dumpEvents
#insertAfter(p.modules,"lepFiller",dumpEvents(level=-1,nanoVersion=12))

### Run the postprocessor
p.run()
