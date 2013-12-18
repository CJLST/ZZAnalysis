### ----------------------------------------------------------------------
###
### Example analayzer
###
###----------------------------------------------------------------------

LEPTON_SETUP = 2012
PD = ""
MCFILTER = ""


# Get absolute path
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------

execfile(PyFilePath + "MasterPy/ZZ4lAnalysis.py")         # 2012 reference analysis


### ----------------------------------------------------------------------
  ### Replace parameters
### ----------------------------------------------------------------------
process.source.fileNames = cms.untracked.vstring(
#    '/store/cmst3/user/cmgtools/CMG/ZZTo4mu_8TeV-powheg-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM/V5/PAT_CMG_V5_4_0/cmgTuple_1.root'        # qqZZ 
#    '/store/cmst3/user/cmgtools/CMG/GluGluToZZTo4L_8TeV-gg2zz-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM/V5/PAT_CMG_V5_4_0/cmgTuple_28.root' # ggZZ
#    'root://cmsphys05//data/b/botta/V5_4_0/patTuple_H126Summer12.root' #Summer12 H126 for FSR synch 
       'root://cmsphys05//data/b/botta/V5_4_0/cmgTuple_H126Summer12.root' #Summer12 H126 for FSR synch
    )


process.maxEvents.input = -1


### ----------------------------------------------------------------------
### Output root file (monitoring histograms)
### ----------------------------------------------------------------------
process.TFileService=cms.Service('TFileService',
                                fileName=cms.string('ZZ4lAnalysis.root')
                                )


### ----------------------------------------------------------------------
### Additional Path 
### ----------------------------------------------------------------------
process.CRPath = cms.Path(process.CR)
#process.CRPath = cms.Path(process.bareZCand + process.ZCand +  process.ZlCand)
#process.HF = cms.Path(process.heavyflavorfilter)

### ----------------------------------------------------------------------
### Analyzer for Plots
### ----------------------------------------------------------------------


PlotSetup = cms.EDAnalyzer("ZZ4lAnalyzer",
                           channel = cms.untracked.string('aChannel'),
                           candCollection = cms.untracked.string('aCand'),
                           isMC = cms.untracked.bool(IsMC),
                           sampleType = cms.int32(SAMPLE_TYPE),
                           setup = cms.int32(LEPTON_SETUP),
                           skimPaths = cms.vstring(SkimPaths),
                           PD = cms.string(PD),
                           MCFilterPath = cms.string(MCFILTER),
                           skipEmptyEvents = cms.bool(True),
                           )

process.Plots4mu = PlotSetup.clone()
process.Plots4mu.channel           = 'MMMM'
process.Plots4mu.candCollection    = 'MMMMCand'

process.Plots4e = PlotSetup.clone()
process.Plots4e.channel            = 'EEEE'
process.Plots4e.candCollection     = 'EEEECand'

process.Plots2e2mu = PlotSetup.clone()
process.Plots2e2mu.channel         = 'EEMM'
process.Plots2e2mu.candCollection  = 'EEMMCand'


# # All events together
# process.PlotsZZ    = cms.EDAnalyzer("ZZ4lAnalyzer",
#                                     channel = cms.untracked.string('ZZ'),
#                                     candCollection = cms.untracked.string('ZZCand'),
#                                     isMC = cms.untracked.bool(IsMC),
#                                     setup = cms.int32(LEPTON_SETUP),
#                                     skimPaths = cms.vstring(SkimPaths),
#                                     PD = cms.string(""),
#                                     MCFilterPath = cms.string(""),
#                                     )

### ----------------------------------------------------------------------
### Analyzer for Trees
### ----------------------------------------------------------------------


TreeSetup = cms.EDAnalyzer("HZZ4lNtupleMaker",
                                   channel = cms.untracked.string('aChannel'),
                                   CandCollection = cms.untracked.string('aCand'),
                                   fileName = cms.untracked.string('candTree'),
                                   isMC = cms.untracked.bool(IsMC),
                                   sampleType = cms.int32(SAMPLE_TYPE),
                                   setup = cms.int32(LEPTON_SETUP),
                                   skimPaths = cms.vstring(SkimPaths),
                                   PD = cms.string(PD),
                                   MCFilterPath = cms.string(MCFILTER),
                                   skipEmptyEvents = cms.bool(True),
                                   )

process.ZZ4muTree = TreeSetup.clone()
process.ZZ4muTree.channel = 'MMMM'
process.ZZ4muTree.CandCollection = 'MMMMCand'

process.ZZ4eTree = TreeSetup.clone()
process.ZZ4eTree.channel = 'EEEE'
process.ZZ4eTree.CandCollection = 'EEEECand'

process.ZZ2e2muTree = TreeSetup.clone()
process.ZZ2e2muTree.channel = 'EEMM'
process.ZZ2e2muTree.CandCollection = 'EEMMCand'

#process.source.eventsToProcess = cms.untracked.VEventRange("1:122490", "1:1343", "1:177684")

# Debug
process.dumpUserData =  cms.EDAnalyzer("dumpUserData",
     dumpTrigger = cms.untracked.bool(True),
     muonSrc = cms.InputTag("appendPhotons:muons"), 
     electronSrc = cms.InputTag("appendPhotons:electrons"),
     candidateSrcs = cms.PSet(
        Zmm   = cms.InputTag("MMCand"),
        Zee   = cms.InputTag("EECand"),
        Z     = cms.InputTag("ZCand"),                                          
        MMMM  = cms.InputTag("MMMMCand"),
        EEEE  = cms.InputTag("EEEECand"),
        EEMM  = cms.InputTag("EEMMCand"),
     )
)


#Print MC history
#process.p = cms.EndPath(process.printTree)

#Dump reconstructed variables
#process.dump = cms.Path(process.dumpUserData)

process.p = cms.EndPath( process.Plots4mu + process.Plots4e + process.Plots2e2mu )
process.trees = cms.EndPath(process.ZZ4muTree * process.ZZ4eTree * process.ZZ2e2muTree )

