### ----------------------------------------------------------------------
###
### Simple analyzer to inspect user data attached to candidates.
###
###----------------------------------------------------------------------
from past.builtins import execfile

LEPTON_SETUP = 2012
PD = ""
MCFILTER = ""
ELECORRTYPE   = "Paper" # "None", "Moriond", or "Paper"
ELEREGRESSION = "Paper" # "None", "Moriond", "PaperNoComb", or "Paper" 
APPLYMUCORR = True


### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------
execfile("MasterPy/ZZ4lAnalysis.py")


### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------
process.source.fileNames = cms.untracked.vstring(
    'root://lxcms00//data3/2013/HZZ_cmgTuple/BE539_H1258TeV.root' #533 V5_15_0 version
#    'root://lxcms00//data3/2013/HZZ_cmgTuple/V5150_VBFH1258TeV.root' #533 V5_15_0 VBF file
    )
 
process.maxEvents.input = 100
# process.options.wantSummary = False


### ----------------------------------------------------------------------
### Add control region collections (not included in default sequence)
### ----------------------------------------------------------------------
process.CRPath = cms.Path(process.CR)


#process.source.eventsToProcess = cms.untracked.VEventRange("1:21715")


process.dumpUserData =  cms.EDAnalyzer("dumpUserData",                                   
#     muonSrc = cms.InputTag("patMuonsWithTrigger"),
#     electronSrc = cms.InputTag("patElectronsWithTrigger::PAT"),
     muonSrc = cms.InputTag("appendPhotons:muons"),
     electronSrc = cms.InputTag("appendPhotons:electrons"),
     candidateSrcs = cms.PSet(
#        Z   = cms.InputTag("ZCand"),
        ZZ  = cms.InputTag("ZZCand"),
#        LL    = cms.InputTag("LLCand"),
#        ZLL   =cms.InputTag("ZLLCand"),    # Starting point for all CRs
#        LLLL  = cms.InputTag("LLLLCand"), # For RFC studies
     )
)

process.analysis = cms.EndPath(process.dumpUserData)





