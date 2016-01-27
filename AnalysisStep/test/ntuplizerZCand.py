### ----------------------------------------------------------------------
###
### Example analyzer for Z->ll ntuples
###
###----------------------------------------------------------------------


try:
    IsMC
except NameError:
    IsMC = True

try:
    LEPTON_SETUP
except NameError:
    LEPTON_SETUP = 2015

try:
    PD
except NameError:
    PD = ""             # "" for MC, "DoubleEle", "DoubleMu", or "MuEG" for data 

try:
    MCFILTER
except NameError:
    MCFILTER = ""

try: 
    XSEC
except NameError:
    XSEC=1
    
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
process.maxEvents.input = -1
#process.options.wantSummary = False


### ----------------------------------------------------------------------
### Output root file (monitoring histograms)
### ----------------------------------------------------------------------
process.TFileService=cms.Service('TFileService',
                                fileName=cms.string('ZZ4lAnalysis.root')
                                )

### ----------------------------------------------------------------------
### Analyzer for Trees
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
     muonSrc = cms.InputTag("appendPhotons:muons"), 
     electronSrc = cms.InputTag("appendPhotons:electrons"),
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

    

