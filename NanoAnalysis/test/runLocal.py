#!/usr/bin/env python3
###
# Example for running the analysis locally, after customizing variables.
# Run with: 
# python runLocal.py
###
from __future__ import print_function
from ZZAnalysis.NanoAnalysis.tools import setConf, getConf, insertAfter

# Check that the checkout recipe has been properly updated 
from ZZAnalysis.AnalysisStep.validateCheckout import validateCheckout 
if not validateCheckout() :
    exit(1)

#SampleToRun = "MCsync_2018Rereco" # for mini vs nano sync
#SampleToRun = "MCsync_2017UL" # for mini vs nano sync
#SampleToRun = "Data2022"
SampleToRun = "MC2022EE"
#SampleToRun = "MC2023postBPix"
#SampleToRun = "MELA_Test"
#SampleToRun = "ggh125_2018UL"
#SampleToRun = "forNanoDoc" # To prepare variable lists with inspectNanoFile.py
#SampleToRun = "Data2024"
#SampleToRun = "MC2024"

### Customize processing variables.
#setConf("runMELA", False)
#setConf("bestCandByMELA", False)
#setConf("APPLYMUCORR", False)
#setConf("APPLYELECORR", False)

## Force filling K factors and weights (default: all off)
#setConf("APPLY_K_NNLOQCD_ZZGG", 1) # 0:None; 1: NNLO/LO; 2: NNLO/NLO; 3: NLO/LO
#setConf("APPLY_K_NNLOQCD_ZZQQB", True)
#setConf("APPLY_K_NNLOEW_ZZQQB", True)
#setConf("APPLY_QCD_GGF_UNCERT", True)

setConf("PROCESS_CR", True)
setConf("PROCESS_ZL", True)
setConf("DEBUG", False)
setConf("SYNCMODE", True) # Force muon resolution correction with fixed +1 sigma smearing
#setConf("ADD_ALLEVENTS", True) # Add extra tree of gen info for all events
#setConf("FILTER_EVENTS", 'Z') # Store all events which contain a good Z candidate
#setConf("FILTER_EVENTS", '3L_20_10') # for trigger studies
#setConf("FILTER_EVENTS", 'NoFilter') # don't skip events with no candidates
#setConf("TRIGPASSTHROUGH", True) #don't skip events failing triggers
#setConf("APPLYJETCORR", False)
#setConf("CANDSTOSTORE",'AllWithRelaxedMuId')

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
elif SampleToRun == "Data2024" :
    setConf("DATA_TAG","")
    setConf("PD","MuEG")
    setConf("PROCESS_CR",True)
    setConf("PROCESS_ZL",True)
    setConf("FILTER_EVENTS","Z")
    setConf("NANOVERSION",15)
    setConf("LEPTON_SETUP",2024)
    setConf("APPLYELECORR",False)
    setConf("APPLYMUCORR",True)
    setConf("APPLYJETCORR",True)
    setConf("IsMC",False)
    setConf("SAMPLENAME","MuonEG2024Cv1")
    setConf("XSEC",-1.0)
    setConf("fileNames",['root://cms-xrd-global.cern.ch//store/data/Run2024C/MuonEG/NANOAOD/MINIv6NANOv15-v1/2530000/694ee48f-c22e-4ae3-83f8-19f98edb481b.root',
                         #'root://cms-xrd-global.cern.ch//store/data/Run2024C/MuonEG/NANOAOD/MINIv6NANOv15-v1/2530000/d5d27820-494a-4dba-aadb-d1d210a90a33.root'
                         ])

################################################################################
elif SampleToRun == "ggh125_2018UL" : ### 2018 UL test sample
    setConf("SAMPLENAME", "ggH125")
    setConf("XSEC", 48.58*0.0002745)
    setConf("LEPTON_SETUP", 2018)
    setConf("NANOVERSION", 9)    
    setConf("DATA_TAG", "UL")
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("fileNames",[
        "/store/mc/RunIISummer20UL18NanoAODv2/WplusH_HToZZTo4L_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/270000/3B6A5CB5-2B7C-924D-85B4-FC3B0C1F4909.root",
        ])

################################################################################
elif SampleToRun == "MCsync_2017UL" :
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
elif SampleToRun == "MCsync_2018Rereco" :
     # Custom-reprocessed Rereco nanoAOD file with updated FSR,
     # corresponding to:/store/mc/RunIIAutumn18NanoAODv7/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/BA6D7F40-ED5E-7D4E-AB14-CE8A9C5DE7EC.root
    setConf("APPLYMUCORR", True)
    setConf("SAMPLENAME", "ggH125")
    setConf("XSEC", 48.58*0.0002745)
    setConf("NANOVERSION", 9)
    setConf("store","")
    setConf("fileNames",["/eos/user/n/namapane/H4lnano/ggH125_fixedFSR.root"])


################################################################################
elif SampleToRun == "MC2022EE" :
    # 2022 MC sample
    setConf("SAMPLENAME", "ggH125")
    setConf("DATA_TAG", "")
    setConf("XSEC", 52.23*0.0002745)
    setConf("LEPTON_SETUP", 2022)
    setConf("IsMC", True)
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("APPLY_QCD_GGF_UNCERT", True) # for ggH
#   setConf("MUON_ID_BYMVA", True)
    setConf("fileNames",[
        "/store/mc/Run3Summer22EENanoAODv12/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg2-JHUGenV752-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/2540000/25c8f5ff-9de0-4a0c-9e2f-757332ad392f.root", # 13158 events
#        "/store/mc/Run3Summer22EENanoAODv12/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg2-JHUGenV752-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/2530000/8f306f2b-1284-41b8-a98f-744267f64b9c.root",
        ])
#    json = {"1": [[1245, 1245],[1306, 1306],[1410, 1410],[1692, 1692],[1903, 1903],[1910, 1910],[1915, 1915],[1927, 1927],[1939, 1939],[1940, 1940],[1944, 1944],[1945, 1945],[1956, 1956],[1960, 1960],[1965, 1965],[1967, 1967],[1968, 1968],[1969, 1969],[2104, 2104]]}

################################################################################
elif SampleToRun == "MC2023postBPix" :
    # 2023 MC sample
    setConf("SAMPLENAME", "ggH125")
    setConf("DATA_TAG", "post_BPix")
    setConf("XSEC", 52.23*0.0002745)
    setConf("LEPTON_SETUP", 2023)
    setConf("IsMC", True)
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("APPLY_QCD_GGF_UNCERT", True) # for ggH
    setConf("fileNames",[
        "/store/mc/Run3Summer23BPixNanoAODv12/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg-jhugen-pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_postBPix_v6-v2/50000/4daf03b4-93c9-4d35-bf53-738bb6cba90e.root", # 12 kevts
        ])

################################################################################
elif SampleToRun == "MC2024" :
    setConf("NANOVERSION", 15)
    setConf("DATA_TAG", "")
    setConf("XSEC", 0.5677*0.0002745)
    setConf("LEPTON_SETUP", 2024)
    setConf("IsMC", True)
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("fileNames",[
        "/store/mc/RunIII2024Summer24NanoAODv15/WminusH-Hto2Zto4L_Par-M-125_TuneCP5_13p6TeV_powhegMINLO-jhugen-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/120000/cbb0799a-f7e4-49cb-b5ba-71693538169b.root", # 11520 events
        ])


################################################################################
elif SampleToRun == "forNanoDoc" :
    # Create a file with a complete set of variables to feed to inspectNanoFile to generate variable documentation
    setConf("SAMPLENAME", "ggH125")
    setConf("DATA_TAG", "post_EE")
    setConf("XSEC", 52.23*0.0002745)
    setConf("LEPTON_SETUP", 2022)
    setConf("IsMC", True)
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("runMELA", True)
    setConf("APPLYMUCORR", True)
    setConf("APPLYELECORR", True)
    setConf("APPLYJETCORR", True)
    # setConf("APPLY_K_NNLOQCD_ZZGG", 1) # requires mcHistoryTools before weightFiller when AllEvents=true, which is not needed in practical cases
    # setConf("APPLY_K_NNLOQCD_ZZQQB", True) # ditto
    setConf("APPLY_K_NNLOEW_ZZQQB", True)
    setConf("APPLY_QCD_GGF_UNCERT", True)
    setConf("PROCESS_CR", True)
    setConf("PROCESS_ZL", True)
    setConf("ADD_ALLEVENTS", True)
    setConf("fileNames",["/store/mc/Run3Summer22EENanoAODv12/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg2-JHUGenV752-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/2540000/25c8f5ff-9de0-4a0c-9e2f-757332ad392f.root"])


################################################################################
elif SampleToRun == "MELA_Test" : 
    setConf("SAMPLENAME", "ggH125")
    setConf("LEPTON_SETUP", 2022)  
    setConf("XSEC", 290.58626*0.0002745)
    setConf("IsMC", True)
    setConf("ADD_ALLEVENTS", True)
    setConf("NANOVERSION", 15)
    setConf("store", "")
    setConf("fileNames", ["/eos/user/n/nipinto/old_CMSSW_13_3_3/src/ggH_test.root"]) # private reprocessing to add LHE mothers/daughters as in v15
    

################################################################################
### Tweak postprocessor parameters as necessary

def customizeProcessForLocal(p) :
    p.prefetch=True # Prefetch remote files
    p.longTermCache=True # keep prefetched files (useful for rerunning tests several times)
    if len(p.inputFiles) == 1 :
        p.haddFileName = None # Skip final hadd

    p.json = json # replace JSON
        
    ### Run only on the first N events in the file
    #p.maxEntries = 10000

    ### Select specific events to debug
    #p.cut = "run==316239  && luminosityBlock==226 && event==284613817"

    ### Print out detailed candidate information for debug purposes
    #from ZZAnalysis.NanoAnalysis.dumpEvents import dumpEvents
    #p.cut = None # Remove preselction
    #insertAfter(p.modules,"lepFiller",dumpEvents(level=-1),getConf("NANOVERSION", 11)) 

    ### Dump MC and LHE history for selected events
    #from ZZAnalysis.NanoAnalysis.mcHistoryDump import mcHistoryDump
    #p.modules.append(mcHistoryDump( printGen=True, printLHE=True))

    ### Read all branches
    #p.branchsel=None 
    #p.outputbranchsel=None #Output all branches

    
setConf("customizations", customizeProcessForLocal, append=True) 


#####################################################################
### This import should be done AFTER all configuration (setConf calls)
from ZZAnalysis.NanoAnalysis.nanoZZ4lAnalysis import *
######################################################################

### Run the postprocessor
p.run()
