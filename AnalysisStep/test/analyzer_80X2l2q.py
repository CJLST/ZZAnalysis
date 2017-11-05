
LEPTON_SETUP = 2016
PD = ""
MCFILTER = ""
ELECORRTYPE = "RunII" # "None", "Moriond", "Paper", or "RunII"
ELEREGRESSION = "Paper" # "None", "Moriond", "PaperNoComb", or "Paper" 
APPLYMUCORR = True
BUNCH_SPACING = 25
#FSRMODE = "Legacy" # Legacy or Run II

CRSync = False #add CR paths

#KEEPLOOSECOMB = True # Do not skip loose lepton ZZ combinations (for debugging

#For DATA: 
IsMC = False
PD = "DoubleEle"

# Get absolute path
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------

execfile(PyFilePath + "analyzer2l2q.py")

process.source.inputCommands = cms.untracked.vstring("keep *", "drop LHERunInfoProduct_*_*_*")

### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------

process.source.fileNames = cms.untracked.vstring(
    '/store/mc/RunIISummer16MiniAODv2/GluGluHToZZTo2L2Q_M1000_13TeV_powheg2_JHUgenV698_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/0ECDBF89-22D5-E611-BFC1-90B11C27F101.root',
    '/store/mc/RunIISummer16MiniAODv2/GluGluHToZZTo2L2Q_M1000_13TeV_powheg2_JHUgenV698_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/109961C9-17D6-E611-A8CA-0CC47A7C3636.root',
    '/store/mc/RunIISummer16MiniAODv2/GluGluHToZZTo2L2Q_M1000_13TeV_powheg2_JHUgenV698_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/3665EE1D-17D6-E611-96D5-0025905A612E.root'
    )

process.calibratedPatElectrons.isSynchronization = cms.bool(True)
process.calibratedMuons.isSynchronization = cms.bool(True)
# process.zJetsFilter.option = 1

process.maxEvents.input = 3000
#process.source.skipEvents = cms.untracked.uint32(5750)

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100


### ----------------------------------------------------------------------
### Analyzer for Plots
### ----------------------------------------------------------------------


#process.source.eventsToProcess = cms.untracked.VEventRange("1:122490", "1:1343", "1:177684")

# Debug
process.dumpUserData =  cms.EDAnalyzer("dumpUserData",
                                       dump = cms.untracked.bool(True),
                                       dumpTrigger = cms.untracked.bool(True),
                                       #     muonSrc = cms.InputTag("slimmedMuons"),
                                       #     electronSrc = cms.InputTag("slimmedElectrons"),
                                       muonSrc = cms.InputTag("appendPhotons:muons"), 
                                       electronSrc = cms.InputTag("appendPhotons:electrons"),
                                       candidateSrcs = cms.PSet(Z     = cms.InputTag("ZCand"),                                  
                                                                ZZ  = cms.InputTag("ZZCand"),
                                                                ZZfat  = cms.InputTag("ZZCandFat"),
                                                                # Zb     = cms.InputTag("bareZCand"),
                                                                # Zjjb  = cms.InputTag("bareZjjCand"),     
                                                                # ZZb  = cms.InputTag("bareZZCand"),
                                                                # ZL  = cms.InputTag("ZlCand"),
                                                                ),
                                       #   jetSrc = cms.InputTag("cleanJets"),
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
#process.ZZTree.skipEmptyEvents = False

#Dump reconstructed variables
#process.appendPhotons.debug = cms.untracked.bool(True)
#process.dump = cms.Path(process.dumpUserData)

#Print MC history
#process.mch = cms.EndPath(process.printTree)
