
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
#IsMC = False
#PD = "DoubleEle"

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
     #'/store/mc/RunIIFall15MiniAODv1/BulkGravToZZToZlepZhad_narrow_M-800_13TeV-madgraph/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/20000/14E8B66A-E5B0-E511-8CDC-B083FED177B1.root'
    '/store/mc/RunIISummer16MiniAODv2/BulkGravToZZToZlepZhad_narrow_M-1000_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/12EF1D5C-A4B6-E611-9596-00266CFCC618.root',
    '/store/mc/RunIISummer16MiniAODv2/BulkGravToZZToZlepZhad_narrow_M-1000_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/42CB5991-A4B6-E611-8243-02163E013F94.root',
    '/store/mc/RunIISummer16MiniAODv2/BulkGravToZZToZlepZhad_narrow_M-1000_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/46D85D63-A4B6-E611-A2E3-848F69FD0EAB.root',
    '/store/mc/RunIISummer16MiniAODv2/BulkGravToZZToZlepZhad_narrow_M-1000_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/4CE3CF95-A3B6-E611-9A9D-D4AE526A0C89.root',
    '/store/mc/RunIISummer16MiniAODv2/BulkGravToZZToZlepZhad_narrow_M-1000_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/5EE83C52-A4B6-E611-BD85-7CD30ABD2EE8.root',
    )

process.calibratedPatElectrons.isSynchronization = cms.bool(True)
process.calibratedMuons.isSynchronization = cms.bool(True)
# process.zJetsFilter.option = 1

process.maxEvents.input = 1000
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
process.ZZTree.skipEmptyEvents = False

#Dump reconstructed variables
#process.appendPhotons.debug = cms.untracked.bool(True)
#process.dump = cms.Path(process.dumpUserData)

#Print MC history
#process.mch = cms.EndPath(process.printTree)
