from ZZAnalysis.AnalysisStep.defaults import *

### ----------------------------------------------------------------------
###
### Example ntuplizer for lepton variables at Z candidate level
###
###----------------------------------------------------------------------

# Set defaults for variables used in this file (in case they are not defined by a caller script)
declareDefault("PD", "", globals()) # "" for MC, "DoubleEle", "DoubleMu", or "MuEG" for data 
declareDefault("MCFILTER", "", globals())
declareDefault("XSEC", 1, globals())

    
# Get absolute path
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------
execfile(PyFilePath + "MasterPy/ZZ4lAnalysis.py")


### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------
process.maxEvents.input = -1
#process.options.wantSummary = False


### ----------------------------------------------------------------------
### Output root file (monitoring histograms)
### ----------------------------------------------------------------------
process.TFileService=cms.Service('TFileService',
                                fileName=cms.string('ZZ4lAnalysis.root')
                                )


### ----------------------------------------------------------------------
### Ntuplizer
### ----------------------------------------------------------------------
process.ZTree = cms.EDAnalyzer("ZNtupleMaker",
                               channel = cms.untracked.string('ZZ'),
                               CandCollection = cms.untracked.string('ZCand'),
                               fileName = cms.untracked.string('candTree'),
                               isMC = cms.untracked.bool(IsMC),
                               sampleType = cms.int32(SAMPLE_TYPE),
                               setup = cms.int32(LEPTON_SETUP),
                               skimPaths = cms.vstring(SkimPaths),
                               PD = cms.string(PD),
                               MCFilterPath = cms.string(MCFILTER),
                               skipEmptyEvents = cms.bool(True),
                               sampleName = cms.string(SAMPLENAME),
                               xsec = cms.double(XSEC)
                               )


### ----------------------------------------------------------------------
### Debug
### ----------------------------------------------------------------------

# Define candidates to be dumped
process.ZFiltered = cms.EDFilter("PATCompositeCandidateRefSelector",
                                 src = cms.InputTag("ZCand"),
                                 cut = cms.string("userFloat('GoodLeptons')")
                                 )
# Select only events with one such candidate
process.ZSelection= cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("ZFiltered"),
                                 minNumber = cms.uint32(1)
                                 )

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
     )
)


### ----------------------------------------------------------------------
### Paths
### ----------------------------------------------------------------------

# overwrite master py: no need to build ZZ candidates
process.Candidates = cms.Path(
       process.muons             +
       process.electrons         + process.cleanSoftElectrons +
       process.fsrPhotons        + process.boostedFsrPhotons +
       process.appendPhotons     +
       process.softLeptons       +
       process.cleanJets         +
       process.bareZCand         + process.ZCand
    )

# ntuplize
process.trees = cms.EndPath(process.ZTree)

# Dump reconstructed variables
#process.dump = cms.Path(process.ZFiltered + process.ZSelection + process.dumpUserData)

    

