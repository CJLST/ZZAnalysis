### Uncomment options to replace default values
LEPTON_SETUP = 2018
#ELECORRTYPE = "None" # "None" to switch off
#ELEREGRESSION = "None" # "None" to switch off
#APPLYMUCORR = False  # Switch off muon scale corrections
#APPLYJEC = False     # Switch off JEC
#APPLYJER = False     # Switch off JER
#RECORRECTMET = False # Switch off MET corr
#KINREFIT = True    # control KinZFitter (very slow)
PROCESS_CR = False   # False = Skip CR paths and trees
#ADDLOOSEELE = True  # Run paths for loose electrons
#APPLYTRIG = False    # Skip events failing required triggers. They are stored with sel<0 if set to False
#KEEPLOOSECOMB = True # Do not skip loose lepton ZZ combinations (for debugging)
ADDZTREE = False # Add tree for Z analysis
ADDLHEKINEMATICS = True  #
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
execfile(PyFilePath + "prod/pyFragments/RecoProbabilities.py")

if not IsMC:
	process.source.inputCommands = cms.untracked.vstring("keep *", "drop LHERunInfoProduct_*_*_*", "drop LHEEventProduct_*_*_*") ###FIXME In 9X this removes all collections for MC

### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------

process.source.fileNames = cms.untracked.vstring(


	
### UL - 2018 sync files
'/store/mc/RunIISummer20UL18MiniAOD/VBF_HToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/30000/D5965022-AF5A-9941-88ED-3C8CBB3E47E5.root',
'/store/mc/RunIISummer20UL18MiniAOD/WplusH_HToZZTo4L_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/110000/A49F09E0-0169-4247-BB63-BB9BFF40C41B.root',

#"/store/mc/RunIISummer20UL18MiniAODv2/ttH_HToZZ_4LFilter_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/260000/C0DA604A-0549-7844-B58B-846A0F0D9EDD.root",
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
