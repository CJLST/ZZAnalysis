### ----------------------------------------------------------------------
###
### An example to illustrate how one can implement event selection in the .py 
###
###----------------------------------------------------------------------
from past.builtins import execfile

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------
execfile("MasterPy/ZZ4lAnalysis.py")


### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------
process.source.fileNames = cms.untracked.vstring(
       '/store/cmst3/user/cmgtools/CMG/GluGluToHToZZTo4L_M-130_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V4_0_2/patTuple_PF2PAT_1_1_YgI.root'
)

process.maxEvents.input = 100
#process.options.wantSummary = False


### ----------------------------------------------------------------------
### Add control region collections (not included in default sequence)
### ----------------------------------------------------------------------
process.CRPath = cms.Path(process.CR)


### ----------------------------------------------------------------------
### An example of a selection. (this is the PRL selection, which in fact is
### already available in ZZ4lAnalyzer_cfg.py)
### ----------------------------------------------------------------------

PRESEL         = ("userFloat('isBestCand') &&" +
                  "userFloat('mZa')>12")
ISO            = ("userFloat('iso34')<0.35")
SIP            = ("userFloat('SIP4')<4")
KIN            = ("daughter('Z1').mass>60 && daughter('Z1').mass<120 &&" +
                  "daughter('Z2').mass>12 && daughter('Z2').mass<120 &&" +
                  "mass>100."
                  )
FULLSEL = ( PRESEL + "&&" +  ISO + "&&" + SIP+ "&&" +  KIN )


### Create a filtered collection of candidates satisfying FULLSEL
process.EEMMFiltered = cms.EDFilter("PATCompositeCandidateRefSelector",
    src = cms.InputTag("EEMMCand"),
    cut = cms.string(FULLSEL)
)

### Select only events with one such candidate
process.selection = cms.EDFilter("CandViewCountFilter",
                                src = cms.InputTag("EEMMFiltered"),
                                minNumber = cms.uint32(1)
                            )

# Ideally, you will pass the filtered collecton EEMMFiltered to your analyzer
# as shown below (this is fake as this analyzer is dummy)
process.dummy =  cms.EDAnalyzer("dumpUserData",
                                src = cms.InputTag("EEMMFiltered") 
                            )

# Run  process.dummy only on events passing process.selection
process.sel = cms.Path(process.EEMMFiltered + process.selection +
                       process.dummy)







