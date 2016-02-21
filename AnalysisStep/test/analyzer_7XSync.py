
LEPTON_SETUP = 2015
PD = ""
MCFILTER = ""
ELECORRTYPE = "RunII" # "None", "Moriond", "Paper", or "RunII"
ELEREGRESSION = "None" # "None", "Moriond", "PaperNoComb", or "Paper" 
APPLYMUCORR = True
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

    ## Fall15 MiniAODv1 files for a first synchronization on MC
    '/store/mc/RunIIFall15MiniAODv1/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/E2490ECF-CBA7-E511-9B19-001E67398458.root',
    '/store/mc/RunIIFall15MiniAODv1/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/282C35FB-68A3-E511-A0C4-0CC47A4C8E5E.root',
    '/store/mc/RunIIFall15MiniAODv1/WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/20000/E2DA5AA7-C5AC-E511-97E0-0CC47A4C8E98.root',

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


#process.source.eventsToProcess = cms.untracked.VEventRange("1:122490", "1:1343", "1:177684")

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
