from past.builtins import execfile
### Uncomment options to replace default values
LEPTON_SETUP = 2022
#ELECORRTYPE = "None" # "None" to switch off
#ELEREGRESSION = "None" # "None" to switch off
APPLYMUCORR = False  # Switch off muon scale corrections
APPLYJEC = False     # Switch off JEC
APPLYJER = False     # Switch off JER
RECORRECTMET = False # Switch off MET corr
#KINREFIT = True    # control KinZFitter (very slow)
PROCESS_CR = False   # False = Skip CR paths and trees
#ADDLOOSEELE = True  # Run paths for loose electrons
#APPLYTRIG = False    # Skip events failing required triggers. They are stored with sel<0 if set to False
#KEEPLOOSECOMB = True # Do not skip loose lepton ZZ combinations (for debugging)
ADDZTREE = False # Add tree for Z analysis
ADDLHEKINEMATICS = False  #
FAILED_TREE_LEVEL = False # To print candTree_failed, if you don't want to save it comment this line

PD = ""
MCFILTER = ""

### For DATA:
#IsMC = False
#PD = "DoubleMu"
#DATA_TAG = "ReReco" # Change to "PromptReco" for Run2018 period D

# Get absolute path
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------

execfile(PyFilePath + "analyzer.py")
#execfile(PyFilePath + "prod/pyFragments/RecoProbabilities.py") #FIXME2022: fails as 13.6 TeV is unknown to SuperDijetMela::SetupResolutionModel: Model MELADifermionResolutionModel

if not IsMC:
	process.source.inputCommands = cms.untracked.vstring("keep *", "drop LHERunInfoProduct_*_*_*", "drop LHEEventProduct_*_*_*") ###FIXME In 9X this removes all collections for MC

### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------

process.source.fileNames = cms.untracked.vstring(
	"/store/mc/Run3Summer22EEMiniAODv3/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg2-JHUGenV752-pythia8/MINIAODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/2540000/0b3804ca-2ae7-46df-bfab-49f92f76b047.root" #ggH, 43129 evts
### 

    )

process.calibratedMuons.isSynchronization = cms.bool(True)

### Events to be processed/picked/skipped
process.maxEvents.input = -1
#process.source.skipEvents = cms.untracked.uint32(5750)
#process.source.eventsToProcess = cms.untracked.VEventRange("1:1711:848227")


# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100


### ----------------------------------------------------------------------
### Debug options
### ----------------------------------------------------------------------

process.dumpUserData =  cms.EDAnalyzer("dumpUserData",
     dumpTrigger = cms.untracked.bool(True),
     muonSrcs = cms.PSet(
#       slimmedMuons = cms.InputTag("slimmedMuons"),
#	calibratedMuons = cms.InputTag("calibratedMuons"),				       
        muons = cms.InputTag("appendPhotons:muons"),
     ),
     electronSrcs = cms.PSet(
#       slimmedElectron = cms.InputTag("slimmedElectrons"),
        electrons = cms.InputTag("appendPhotons:electrons"),
#        RSE = cms.InputTag("appendPhotons:looseElectrons"),
#        TLE = cms.InputTag("appendPhotons:electronstle"), #These are actually photons, should add a photonSrcs section for them.
     ),
     candidateSrcs = cms.PSet(
        Z     = cms.InputTag("ZCand"),
#        ZRSE     = cms.InputTag("ZCandlooseEle"),
#        ZTLE     = cms.InputTag("ZCandtle"),
        ZZ  = cms.InputTag("ZZCand"),
#        ZZRSE     = cms.InputTag("ZZCandlooseEle"),
#        ZZTLE     = cms.InputTag("ZZCandtle"),
#        ZLL  = cms.InputTag("ZLLCand"),
#        ZL  = cms.InputTag("ZlCand"),
     ),
     jetSrc = cms.InputTag("cleanJets"),
)

### Dump reconstructed variables
#process.appendPhotons.debug = cms.untracked.bool(True)
#process.fsrPhotons.debug = cms.untracked.bool(True)
#process.dump = cms.Path(process.dumpUserData)

### Print MC history
#process.mch = cms.EndPath(process.printTree)

### Create lepton sync file
#process.PlotsZZ.dumpForSync = True;
#process.p = cms.EndPath( process.PlotsZZ)

### Keep all events in the tree, even if no candidate is selected
#process.ZZTree.skipEmptyEvents = False

### Replace the paths in analyzer.py
#process.trees = cms.EndPath(process.ZZTree)

### Monitor memory usage
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)
