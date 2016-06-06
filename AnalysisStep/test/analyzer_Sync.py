
LEPTON_SETUP = 2016
PD = ""
MCFILTER = ""
ELECORRTYPE = "None" # "None", "Moriond", "Paper", or "RunII"
ELEREGRESSION = "None" # "None", "Moriond", "PaperNoComb", or "Paper" 
APPLYMUCORR = True #FIXME: False causes a crash in the Z fitter, use MUCORRTYPE="None" for the time being (to be fixed)
MUCORRTYPE = "None" # "None", "MC_76X_13TeV". See above.
BUNCH_SPACING = 25
#FSRMODE = "Legacy" # Legacy or Run II

CRSync = False #add CR paths

#KEEPLOOSECOMB = True # Do not skip loose lepton ZZ combinations (for debugging

#For DATA: 
#IsMC = False
#PD = "DoubleEle"

# Get absolute path
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------

execfile(PyFilePath + "analyzer.py")

process.source.inputCommands = cms.untracked.vstring("keep *", "drop LHERunInfoProduct_*_*_*", "drop LHEEventProduct_*_*_*")

### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------

process.source.fileNames = cms.untracked.vstring(

    ## Spring16 MiniAODv1 files for synchronization
    '/store/mc/RunIISpring16MiniAODv1/VBF_HToZZTo4L_M190_13TeV_powheg2_JHUgenV6_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1/00000/28ADE6D5-021F-E611-B1A4-00145E5521B9.root',
    '/store/mc/RunIISpring16MiniAODv1/WminusH_HToZZTo4L_M150_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1/00000/C6DF34D4-CF20-E611-A8EE-782BCB27B958.root'

    )

process.calibratedPatElectrons.isSynchronization = cms.bool(True)
process.calibratedMuons.isSynchronization = cms.bool(True)

process.maxEvents.input = -1
#process.source.skipEvents = cms.untracked.uint32(5750)

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


### ----------------------------------------------------------------------
### Analyzer for Plots
### ----------------------------------------------------------------------


#process.source.eventsToProcess = cms.untracked.VEventRange("1:72285")

# Debug
process.dumpUserData =  cms.EDAnalyzer("dumpUserData",
     dumpTrigger = cms.untracked.bool(True),
#     muonSrc = cms.InputTag("slimmedMuons"),
#     electronSrc = cms.InputTag("slimmedElectrons"),
     muonSrc = cms.InputTag("appendPhotons:muons"), 
     electronSrc = cms.InputTag("appendPhotons:electrons"),
     candidateSrcs = cms.PSet(
        Z     = cms.InputTag("ZCand"),                                  
        ZZ  = cms.InputTag("ZZCand"),
#        ZLL  = cms.InputTag("ZLLCand"),
#        ZL  = cms.InputTag("ZlCand"),
     ),
     jetSrc = cms.InputTag("cleanJets"),
)

# Create lepton sync file
#process.PlotsZZ.dumpForSync = True;
#process.p = cms.EndPath( process.PlotsZZ)

# Keep all events in the tree, even if no candidate is selected
#process.ZZTree.skipEmptyEvents = False


# Also process CRs
if (CRSync) :
    process.CRPath = cms.Path(process.CR)
    process.CRtrees = cms.EndPath(process.CRZLLTree + process.CRZLTree)

# replace the paths in analyzer.py
process.trees = cms.EndPath(process.ZZTree)

#Dump reconstructed variables
#process.appendPhotons.debug = cms.untracked.bool(True)
#process.dump = cms.Path(process.dumpUserData)

#Print MC history
#process.mch = cms.EndPath(process.printTree)
