
LEPTON_SETUP = 2017
PD = ""
MCFILTER = ""
ELECORRTYPE = "RunII" # "None", "Moriond", "Paper", or "RunII"
ELEREGRESSION = "None" # "None", "Moriond", "PaperNoComb", or "Paper" 
APPLYMUCORR = True
BUNCH_SPACING = 25

#For DATA: 
#IsMC = False
#PD = "DoubleEle"

# Get absolute path
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------

execfile(PyFilePath + "ntuplizerZCand.py")

process.TFileService.fileName = cms.string('ZCandNtuple.root')

process.source.inputCommands = cms.untracked.vstring("keep *", "drop LHERunInfoProduct_*_*_*", "drop LHEEventProduct_*_*_*")

### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------

process.source.fileNames = cms.untracked.vstring(
    
    '/store/mc/RunIISummer16MiniAODv2/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/221CC46F-2FC6-E611-8FFC-0CC47A1E0488.root',
    '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-herwigpp_30M/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/043B3BC8-8BB7-E611-8CCE-0090FAA573E0.root',

    )

process.calibratedPatElectrons.isSynchronization = cms.bool(True)
process.calibratedMuons.isSynchronization = cms.bool(True)

process.maxEvents.input = -1
#process.source.skipEvents = cms.untracked.uint32(5750)

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#process.source.eventsToProcess = cms.untracked.VEventRange("1:122490", "1:1343", "1:177684")

# Debug
process.dumpUserData =  cms.EDAnalyzer("dumpUserData",
     dumpTrigger = cms.untracked.bool(True),
     muonSrcs =  cms.PSet(
        muons = cms.InputTag("appendPhotons:muons"),
     ),
     electronSrcs = cms.PSet(
        electrons = cms.InputTag("appendPhotons:electrons"),
     ),
     candidateSrcs = cms.PSet(
        Z = cms.InputTag("ZCand"),
     ),
     jetSrc = cms.InputTag("cleanJets"),
)

# Keep all events in the tree, even if no candidate is selected
#process.ZTree.skipEmptyEvents = False


# replace the paths in ntuplizerZCand.py
#process.trees = cms.EndPath(process.ZTree)

#Dump reconstructed variables
#process.appendPhotons.debug = cms.untracked.bool(True)
#process.dump = cms.Path(process.dumpUserData)

#Print MC history
#process.mch = cms.EndPath(process.printTree)
