#!/usr/bin/env python
###
# Example for running the analysis locally, after customizing variables.
# Run with: 
# python runLocal.py
###
from __future__ import print_function
from ZZAnalysis.NanoAnalysis.tools import setConf, getConf, insertAfter

#SampleToRun = "Data2022"
SampleToRun = "MCsync_UL"

### Customize processing variables
#setConf("runMELA", False)
#setConf("bestCandByMELA", False)

## Select specific events to debug
#setConf("preselection","run==316239  && luminosityBlock==226 && event==284613817")

## Force filling K factors and weights (default: all off)
#setConf("APPLY_K_NNLOQCD_ZZGG", 1) # 0:None; 1: NNLO/LO; 2: NNLO/NLO; 3: NLO/LO
#setConf("APPLY_K_NNLOQCD_ZZQQB", True)
#setConf("APPLY_K_NNLOEW_ZZQQB", True)
#setConf("APPLY_QCD_GGF_UNCERT", True)

setConf("PROCESS_CR", True)
setConf("DEBUG", False)
setConf("SYNCMODE", True) # Force muon resolution correction with fixed +1 sigma smearing


################################################################################
if SampleToRun == "Data2022" :
    # 2022 data sample from /MuonEG/Run2022D-PromptNanoAODv10_v1-v1/NANOAOD
    setConf("IsMC", False)
    setConf("LEPTON_SETUP", 2022)
    setConf("PD", "any")
    setConf("SAMPLENAME", "test")
    setConf("TRIGPASSTHROUGH", True)
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("fileNames",[
        "/store/data/Run2022D/MuonEG/NANOAOD/PromptNanoAODv10_v2-v1/50000/68f42f42-3274-46ec-b23d-bfadc13012c2.root",
        ])


################################################################################
elif SampleToRun == "ggh125_UL" : ### 2018 UL test sample
    setConf("SAMPLENAME", "ggH125")
    setConf("XSEC", 48.58*0.0002745)
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("fileNames",[
        "/store/mc/RunIISummer20UL18NanoAODv2/WplusH_HToZZTo4L_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/270000/3B6A5CB5-2B7C-924D-85B4-FC3B0C1F4909.root",
        ])

################################################################################
elif SampleToRun == "MCsync_UL" :
    # Custom-reprocessed Rereco nanoAOD file with updated FSR and electron MVA,
    # no packing for genparticle p3; 26000 events
    # corresponding to:/store/mc/RunIISummer20UL17MiniAODv2/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v2/130000/3E4E8D55-3993-2B43-AF3B-7AB45BBE0BDA.root
    # Note that electron mva variable is electron_mvaHZZIso as for nanoAOD v10
    setConf("SAMPLENAME", "ggH125")
    setConf("XSEC", 48.58*0.0002745)
    setConf("LEPTON_SETUP", 2017)
    setConf("store","")
    setConf("fileNames",["/eos/user/n/namapane/H4lnano/ggH125_2017UL_fixedFSR.root"])


################################################################################
elif SampleToRun == "MCsync_Rereco" :
     # Custom-reprocessed Rereco nanoAOD file with updated FSR,
     # corresponding to:/store/mc/RunIIAutumn18NanoAODv7/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/BA6D7F40-ED5E-7D4E-AB14-CE8A9C5DE7EC.root
     # NOTE: by default, electron cuts are for UL in lepFiller.py!
    setConf("SAMPLENAME", "ggH125")
    setConf("XSEC", 48.58*0.0002745)

    setConf("store","")
    setConf("fileNames",["/eos/user/n/namapane/H4lnano/ggH125_fixedFSR.root"])


#####################################################################
### This import should be done AFTER all customizations (setConf calls)
from ZZAnalysis.NanoAnalysis.nanoZZ4lAnalysis import *
######################################################################

### Tweak postprocessor parameters as necessary
p.prefetch=False # Skip local prefetching
#p.longTermCache=True # keep prefetched file (useful for rerunning tests several times)
if len(p.inputFiles) == 1 :
    p.haddFileName = None # Skip final hadd
#p.maxEntries = 1000

### Print out detailed candidate information for debug purposes
#p.cut = None # Remove preselction
#from ZZAnalysis.NanoAnalysis.dumpEvents import dumpEvents
#insertAfter(p.modules,"lepFiller",dumpEvents(level=-1)) 

#p.branchsel=None #Read all branches
#p.outputbranchsel=None #Output all branches


### Run the postprocessor
p.run()
