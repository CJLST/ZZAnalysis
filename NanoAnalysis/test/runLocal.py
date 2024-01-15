#!/usr/bin/env python3
###
# Example for running the analysis locally, after customizing variables.
# Run with: 
# python runLocal.py
###
from __future__ import print_function
from ZZAnalysis.NanoAnalysis.tools import setConf, getConf, insertAfter

#SampleToRun = "MCsync_Rereco"
#SampleToRun = "MCsync_UL"
#SampleToRun = "Data2022"
SampleToRun = "MC2022"

### Customize processing variables
#setConf("runMELA", False)
#setConf("bestCandByMELA", False)
setConf("APPLYMUCORR", False) #NOTE: mu corrections removed for comparision with mini since they are not deterministic on nanoAOD.


## Force filling K factors and weights (default: all off)
#setConf("APPLY_K_NNLOQCD_ZZGG", 1) # 0:None; 1: NNLO/LO; 2: NNLO/NLO; 3: NLO/LO
#setConf("APPLY_K_NNLOQCD_ZZQQB", True)
#setConf("APPLY_K_NNLOEW_ZZQQB", True)
#setConf("APPLY_QCD_GGF_UNCERT", True)

setConf("PROCESS_CR", True)
setConf("DEBUG", False)
setConf("SYNCMODE", True) # Force muon resolution correction with fixed +1 sigma smearing
#setConf("ADD_ALLEVENTS", True) # Add extra tree of gen info for all events

json = None #replace this if needed


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
    setConf("LEPTON_SETUP", 2018)
    setConf("DATA_TAG", "UL")
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("fileNames",[
        "/store/mc/RunIISummer20UL18NanoAODv2/WplusH_HToZZTo4L_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/270000/3B6A5CB5-2B7C-924D-85B4-FC3B0C1F4909.root",
        ])

################################################################################
elif SampleToRun == "MCsync_UL" :
    # Custom-reprocessed Rereco nanoAOD file with updated FSR and electron MVA,
    # no packing for genparticle p3; 26000 events
    # corresponding to:/store/mc/RunIISummer20UL17MiniAODv2/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v2/130000/3E4E8D55-3993-2B43-AF3B-7AB45BBE0BDA.root
    setConf("SAMPLENAME", "ggH125")
    setConf("XSEC", 48.58*0.0002745)
    setConf("LEPTON_SETUP", 2017)
    setConf("NANOVERSION", 10) # variable defined as per nanoAOD v10 (notably electron_mvaHZZIso)
    setConf("DATA_TAG", "UL")
    setConf("store","")
    setConf("fileNames",["/eos/user/n/namapane/H4lnano/ggH125_2017UL_fixedFSR.root"])
#    setConf("fileNames",["/eos/user/n/namapane/H4lnano/ggH125_2017UL_fixedFSR_nopacking.root"]) # with no packing of muon eta, phi, mass


################################################################################
elif SampleToRun == "MCsync_Rereco" :
     # Custom-reprocessed Rereco nanoAOD file with updated FSR,
     # corresponding to:/store/mc/RunIIAutumn18NanoAODv7/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/BA6D7F40-ED5E-7D4E-AB14-CE8A9C5DE7EC.root
    setConf("APPLYMUCORR", True)
    setConf("SAMPLENAME", "ggH125")
    setConf("XSEC", 48.58*0.0002745)
    setConf("NANOVERSION", 9)
    setConf("store","")
    setConf("fileNames",["/eos/user/n/namapane/H4lnano/ggH125_fixedFSR.root"])


################################################################################
elif SampleToRun == "MC2022" :
    # 2022 MC sample
    #/store/mc/Run3Summer22EEMiniAODv3/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg2-JHUGenV752-pythia8/MINIAODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/2540000/0b3804ca-2ae7-46df-bfab-49f92f76b047.root
    setConf("SAMPLENAME", "ggH125")
    setConf("XSEC", 52.23*0.0002745)
    setConf("LEPTON_SETUP", 2022)
    setConf("IsMC", True)
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("fileNames",[
        "/store/mc/Run3Summer22EENanoAODv12/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg2-JHUGenV752-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/2540000/25c8f5ff-9de0-4a0c-9e2f-757332ad392f.root",
#        "/store/mc/Run3Summer22EENanoAODv12/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg2-JHUGenV752-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/2530000/8f306f2b-1284-41b8-a98f-744267f64b9c.root",
        ])
#    json = {"1": [[1245, 1245],[1306, 1306],[1410, 1410],[1692, 1692],[1903, 1903],[1910, 1910],[1915, 1915],[1927, 1927],[1939, 1939],[1940, 1940],[1944, 1944],[1945, 1945],[1956, 1956],[1960, 1960],[1965, 1965],[1967, 1967],[1968, 1968],[1969, 1969],[2104, 2104]]}



#####################################################################
### This import should be done AFTER all customizations (setConf calls)
from ZZAnalysis.NanoAnalysis.nanoZZ4lAnalysis import *
######################################################################

### Tweak postprocessor parameters as necessary
p.prefetch=True # Prefetch remote files
p.longTermCache=True # keep prefetched files (useful for rerunning tests several times)
if len(p.inputFiles) == 1 :
    p.haddFileName = None # Skip final hadd
#p.maxEntries = 1000

### Select specific events to debug
#p.cut = "run==316239  && luminosityBlock==226 && event==284613817"

### Print out detailed candidate information for debug purposes
#from ZZAnalysis.NanoAnalysis.dumpEvents import dumpEvents
#p.cut = None # Remove preselction
#insertAfter(p.modules,"lepFiller",dumpEvents(level=-1),getConf("NANOVERSION", 11)) 

#p.branchsel=None #Read all branches
#p.outputbranchsel=None #Output all branches

#replace JSON
p.json = json

### Run the postprocessor
p.run()
