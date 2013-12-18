### ----------------------------------------------------------------------
###
### HCP Synchronization file for 533.
###
### FIXME: the path includes the rochester correction (which is a passtrough for 53X right now);
### the syncronization is supposed to be done with no escale applied
###
###----------------------------------------------------------------------

LEPTON_SETUP = 2011
PD = ""
MCFILTER = ""
ELECORRTYPE   = "Paper" # "None", "Moriond", or "Paper"
ELEREGRESSION = "Paper" # "None", "Moriond", "PaperNoComb", or "Paper" 
APPLYMUCORR = True

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
#    'root://lxcms00//data3/2012/HZZ_cmgTuple/synchHCP2/H120_44X_V590.root'
    'root://lxcms00//data3/2013/HZZ_cmgTuple/BE445_H1207TeV.root' #445 V5_15_0
)


process.maxEvents.input = -1

### FIXME: DEBUGGING ONLY, turns off smearing if used with special tag
process.calibratedPatElectrons.synchronization = True
#process.calibratedPatElectrons.smearingRatio = 1
#process.calibratedMuons.fakeSmearing = cms.untracked.bool(True)
process.calibratedMuons.fakeSmearing = cms.bool(True)

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
#     muonSrc = cms.InputTag("patMuonsWithTrigger"),
#   electronSrc = cms.InputTag("patElectronsWithTrigger::PAT"),
     muonSrc = cms.InputTag("appendPhotons:muons"), 
     electronSrc = cms.InputTag("appendPhotons:electrons"),
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

process.Plots4mu.dumpMC   = cms.untracked.bool(True)
process.Plots4e.dumpMC    = cms.untracked.bool(True)
process.Plots2e2mu.dumpMC = cms.untracked.bool(True)

# replace the paths in analyzer.py
process.p = cms.EndPath( process.Plots4mu + process.Plots4e + process.Plots2e2mu )
process.trees = cms.EndPath(process.ZZ4muTree * process.ZZ4eTree * process.ZZ2e2muTree )

