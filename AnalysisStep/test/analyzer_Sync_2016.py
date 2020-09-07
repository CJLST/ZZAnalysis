
#DATA_TAG = "ReReco" # Change to PromptReco for Run2016 period H
LEPTON_SETUP = 2016  # current default = 2017 = Moriond2017
APPLYMUCORR = True  # Switch off muon scale corrections
APPLYJEC = True     #
APPLYJER = True     #
RECORRECTMET = True #
#KINREFIT = True    # control KinZFitter (very slow)
PROCESS_CR = True   # Uncomment to run CR paths and trees
#ADDLOOSEELE = True  # Run paths for loose electrons
#APPLYTRIG = False    # hack for samples missing correct triggers - use with caution
#KEEPLOOSECOMB = True # Do not skip loose lepton ZZ combinations (for debugging)
ADDZTREE = True      # Add tree for Z analysis

PD = ""
MCFILTER = ""

#For DATA: 
#IsMC = False
#PD = "DoubleMu"

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
### LEGACY PAPER - 2016 sync files
'/store/mc/RunIISummer16MiniAODv2/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV709_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/20000/8A6DC1B7-D1D3-E711-8D88-002590DE6E32.root',
'/store/mc/RunIISummer16MiniAODv2/ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUGenV709_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/30000/CCDA7461-FDD7-E711-AFA4-008CFAF2224C.root',                                                                                                                                                                                                        
'/store/mc/RunIISummer16MiniAODv2/ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV709_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/20000/622DA168-2DD2-E711-AE8A-E0071B7A5650.root'
)

#process.calibratedPatElectrons.isSynchronization = cms.bool(True) #process.calibratedPatElectrons.isSynchronization = cms.bool(True) # Not needed anymore since new EGamma smearing is event deterministic
#process.calibratedMuons.isSynchronization = cms.bool(True)

process.maxEvents.input = -1
#process.source.skipEvents = cms.untracked.uint32(5750)

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


### ----------------------------------------------------------------------
### Analyzer for Plots
### ----------------------------------------------------------------------


#process.source.eventsToProcess = cms.untracked.VEventRange("1:8670")

# Debug
process.dumpUserData =  cms.EDAnalyzer("dumpUserData",
     dumpTrigger = cms.untracked.bool(True),
     muonSrcs = cms.PSet(
#       slimmedMuons = cms.InputTag("slimmedMuons"),
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

# Create lepton sync file
#process.PlotsZZ.dumpForSync = True;
#process.p = cms.EndPath( process.PlotsZZ)

# Keep all events in the tree, even if no candidate is selected
#process.ZZTree.skipEmptyEvents = False

# replace the paths in analyzer.py
#process.trees = cms.EndPath(process.ZZTree)

#Dump reconstructed variables
#process.appendPhotons.debug = cms.untracked.bool(True)
#process.fsrPhotons.debug = cms.untracked.bool(True)
#process.dump = cms.Path(process.dumpUserData)

#Print MC history
#process.mch = cms.EndPath(process.printTree)


#Monitor memory usage
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)
