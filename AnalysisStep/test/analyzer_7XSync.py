
LEPTON_SETUP = 2015
PD = ""
MCFILTER = ""
ELECORRTYPE = "RunII" # "None", "Moriond", "Paper", or "RunII"
ELEREGRESSION = "None" # "None", "Moriond", "PaperNoComb", or "Paper" 
APPLYMUCORR = False
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

process.calibratedElectrons.isSynchronization = cms.bool(True) #This replaces the smearing with a shift of 1 sigma

### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------

process.source.fileNames = cms.untracked.vstring(
        
#    'file:/afs/cern.ch/user/g/gortona/work/public/skimmedSamples/DYskim_Spring15.root' # The old CR MC sync file
   
    ## Old Spring15 sync files
#     '/store/mc/RunIISpring15DR74/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/68791C0A-3013-E511-88FD-D4AE5269F5FF.root',
#     '/store/mc/RunIISpring15DR74/WplusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/60000/04BD6860-9F08-E511-8A80-842B2B1858FB.root',
#     '/store/mc/RunIISpring15DR74/WminusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/70000/4A9FED55-DF0C-E511-A4B2-3417EBE6471D.root',
#     '/store/mc/RunIISpring15DR74/ZH_HToZZ_4LFilter_M125_13TeV_powheg-minlo-HZJ_JHUgen_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/104B7067-0C02-E511-8FFB-0030487D07BA.root',


    
    ## a ggH miniAODv2 file
    # '/store/mc/RunIISpring15MiniAODv2/GluGluHToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/2E642290-C56D-E511-AD68-008CFA110C88.root', 


    ## MiniAODv2 sync files
    '/store/mc/RunIISpring15MiniAODv2/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/30000/3E964C5D-1D6E-E511-8B9A-0050560207C5.root',
    '/store/mc/RunIISpring15MiniAODv2/WminusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/30000/D8CA6B54-056F-E511-BB1A-02163E014CE3.root',
    '/store/mc/RunIISpring15MiniAODv2/WplusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/D22BEE88-C26D-E511-B330-002590A81EF0.root',
    '/store/mc/RunIISpring15MiniAODv2/ttH_HToZZ_4LFilter_M125_13TeV_powheg_JHUgen_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/80000/84F62DD7-1475-E511-9F59-009C02AB98A6.root'

    )


### FIXME: DEBUGGING ONLY, turns off smearing if used with special tag
#process.calibratedPatElectrons.synchronization = True
#process.calibratedMuons.fakeSmearing = cms.bool(True)

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
