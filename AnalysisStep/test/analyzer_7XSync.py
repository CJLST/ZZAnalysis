### ----------------------------------------------------------------------
###
### HCP Synchronization file for 533.
###
### FIXME: the path includes the rochester correction (which is a passtrough for 53X right now);
### the syncronization is supposed to be done with no escale applied
###
###----------------------------------------------------------------------

LEPTON_SETUP = 2012
PD = ""
MCFILTER = ""
ELECORRTYPE   = "None" # "None", "Moriond", or "Paper"
ELEREGRESSION = "None" # "None", "Moriond", "PaperNoComb", or "Paper" 
APPLYMUCORR = False

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


### ----------------------------------------------------------------------
  ### Replace parameters
### ----------------------------------------------------------------------
process.source.fileNames = cms.untracked.vstring(
    #'/store/cmst3/user/gpetrucc/miniAOD/v1/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_PU_S14_PAT.root'
    #'/store/cmst3/user/gpetrucc/miniAOD/v1/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_PU_S14_PAT_big.root' #S14 = 50ns scenario, GT: PLS170_V6AN1
    '/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/148E558C-946F-E411-AFA7-7845C4FC3A52.root'
    #"/store/mc/Spring14miniaod/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/0881ABEB-2709-E411-9E42-00145EDD7581.root" # 1st file from the central 20bx25 sample, GT: PLS170_V7AN1
    )


### FIXME: DEBUGGING ONLY, turns off smearing if used with special tag
#process.calibratedPatElectrons.synchronization = True
#process.calibratedMuons.fakeSmearing = cms.bool(True)

#process.appendPhotons.debug = cms.untracked.bool(True)

process.maxEvents.input = 5000

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000


### ----------------------------------------------------------------------
### Analyzer for Plots
### ----------------------------------------------------------------------


#process.source.eventsToProcess = cms.untracked.VEventRange("1:122490", "1:1343", "1:177684")
#process.source.eventsToProcess = cms.untracked.VEventRange("1:30602")
#process.source.eventsToProcess = cms.untracked.VEventRange("191062:303:330091192")

# Debug
process.dumpUserData =  cms.EDAnalyzer("dumpUserData",
     dumpTrigger = cms.untracked.bool(True),
     muonSrc = cms.InputTag("slimmedMuons"),
     electronSrc = cms.InputTag("slimmedElectrons"),
#     muonSrc = cms.InputTag("appendPhotons:muons"), 
#     electronSrc = cms.InputTag("appendPhotons:electrons"),
     candidateSrcs = cms.PSet(
        Zmm   = cms.InputTag("MMCand"),
        Zee   = cms.InputTag("EECand"),
#        Z     = cms.InputTag("ZCand"),                                  
        MMMM  = cms.InputTag("MMMMCand"),
        EEEE  = cms.InputTag("EEEECand"),
        EEMM  = cms.InputTag("EEMMCand"),
     )
)


#Print MC history
#process.p = cms.EndPath(process.printTree)

#Dump reconstructed variables
#process.dump = cms.Path(process.dumpUserData)

#process.Plots4mu.dumpMC   = cms.untracked.bool(True)
#process.Plots4e.dumpMC    = cms.untracked.bool(True)
#process.Plots2e2mu.dumpMC = cms.untracked.bool(True)



# replace the paths in analyzer.py
#process.p = cms.EndPath( process.Plots4mu + process.Plots4e + process.Plots2e2mu )
process.trees = cms.EndPath(process.ZZ4muTree * process.ZZ4eTree * process.ZZ2e2muTree )

