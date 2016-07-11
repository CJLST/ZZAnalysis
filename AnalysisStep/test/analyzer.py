from ZZAnalysis.AnalysisStep.defaults import *
from ZZAnalysis.AnalysisStep.couplings import *

### ----------------------------------------------------------------------
###
### Example analyzer
###
###----------------------------------------------------------------------

# Set defaults for variables used in this file (in case they are not defined by a caller script)
declareDefault("PD", "", globals()) # "" for MC, "DoubleEle", "DoubleMu", or "MuEG" for data
declareDefault("MCFILTER", "", globals())
declareDefault("XSEC", 1, globals())
declareDefault("PROCESS_CR", False, globals())

# Couplings for reweighting
declareDefault("REWEIGHTING_TYPE", "none", globals())
declareDefault("SPIN", 0, globals())
couplings = Couplings()
for coupling in couplings.allnames():
    declareDefault(coupling, 0, globals())
    couplings[coupling] = globals()[coupling]

# LHE info
#  VVDECAYMODE\VVMODE  / ZZ==1 / WW==0  / Yukawa==2 / Zgam=3 / gamgam=4 / Z+nj=5
#                     0: 4l    / lnulnu / 2l        / 2l     / gam      / 2l
#                     1: 4q    / 4q     / 2q        / 2q     / -        / 2q
#                     2: 2l2q  / lnu2q  / -         / -      / -        / -
#                     3: 2l2nu / -      / -         / -      / -        / -
#                     4: 2q2nu / -      / -         / -      / -        / -
#                     5: 4nu   / -      / -         / 2nu    / -        / 2nu
#                    -1: [ Any                                                 ]
#                    -2: [ 2l2X         ]
#                    -3: [ 2nu2X        ]
#                    -4: [ 2q2X         ]
declareDefault("VVMODE", 0, globals())
declareDefault("VVDECAYMODE", 0, globals())
declareDefault("ADDLHEKINEMATICS", False, globals())
declareDefault("SAMPLEPRODUCTIONID", "", globals()) # Reserve for reweighting, not yet used
declareDefault("SAMPLEPROCESSID", "", globals()) # Reserve for reweighting, not yet used

# K factors
declareDefault("APPLY_K_NNLOQCD_ZZGG", 0, globals()) # 0: Do not; 1: NNLO/LO; 2: NNLO/NLO; 3: NLO/LO
declareDefault("APPLY_K_NNLOQCD_ZZQQB", False, globals())
declareDefault("APPLY_K_NLOEW_ZZQQB", False, globals())

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
### Analyzers for Plots
### ----------------------------------------------------------------------

# All events together
process.PlotsZZ    = cms.EDAnalyzer("ZZ4lAnalyzer",
                                    channel = cms.untracked.string('ZZ'),
                                    candCollection = cms.untracked.string('ZZCand'),
                                    isMC = cms.untracked.bool(IsMC),
                                    sampleType = cms.int32(SAMPLE_TYPE),
                                    setup = cms.int32(LEPTON_SETUP),
                                    skimPaths = cms.vstring(SkimPaths),
                                    PD = cms.string(PD),
                                    MCFilterPath = cms.string(MCFILTER),
                                    sampleName = cms.string(SAMPLENAME),
                                    dumpForSync = cms.untracked.bool(False),
                                    )

### Control Region Plots

# PlotCRSetup    = cms.EDAnalyzer("ZZ4lAnalyzerCR",
#                                 channel = cms.untracked.string('aChannel'),
#                                 candCollection = cms.untracked.string('aCand'),
#                                 isMC = cms.untracked.bool(IsMC),
#                                 sampleType = cms.int32(SAMPLE_TYPE),
#                                 setup = cms.int32(LEPTON_SETUP),
#                                 skimPaths = cms.vstring(SkimPaths),
#                                 PD = cms.string(PD),
#                                 MCFilterPath = cms.string(MCFILTER),
#                                 )

# process.PlotsCRZLL = PlotCRSetup.clone()
# process.PlotsCRZLL.channel               = "ZLL"
# process.PlotsCRZLL.candCollection        = 'ZLLCand'


#Count events with at least 1 Z
# process.ZFiltered = cms.EDFilter("PATCompositeCandidateRefSelector",
#     src = cms.InputTag("ZCand"),
#     cut = cms.string("userFloat('GoodLeptons')")
# )
# process.sStep4 = cms.EDFilter("CandViewCountFilter",
#                               src = cms.InputTag("ZFiltered"),
#                               minNumber = cms.uint32(1)
#                               )
#process.step4 = cms.Path(process.SkimSequence + process.ZFiltered + process.sStep4 )


#Count events with at least 1 ZZ
# process.ZZFiltered = cms.EDFilter("PATCompositeCandidateRefSelector",
#     src = cms.InputTag("ZZCand"),
#     cut = cms.string("userFloat('GoodLeptons')")
# )
# process.sStep5 = cms.EDFilter("CandViewCountFilter",
#                                 src = cms.InputTag("ZZFiltered"),
#                                 minNumber = cms.uint32(1)
#                             )
#process.step5 = cms.Path(process.SkimSequence + process.ZZFiltered + process.sStep5 )


### ----------------------------------------------------------------------
### Analyzer for Trees
### ----------------------------------------------------------------------

TreeSetup = cms.EDAnalyzer("HZZ4lNtupleMaker",
                           channel = cms.untracked.string('aChannel'),
                           CandCollection = cms.untracked.string('ZZCand'),
                           fileName = cms.untracked.string('candTree'),
                           isMC = cms.untracked.bool(IsMC),
                           sampleType = cms.int32(SAMPLE_TYPE),
                           setup = cms.int32(LEPTON_SETUP),
                           skimPaths = cms.vstring(SkimPaths),
                           PD = cms.string(PD),
                           MCFilterPath = cms.string(MCFILTER),
                           skipEmptyEvents = cms.bool(True),
                           sampleName = cms.string(SAMPLENAME),
                           superMelaMass = cms.double(SUPERMELA_MASS),
                           xsec = cms.double(XSEC),
                           spin = cms.int32(int(SPIN)),
                           HVVcouplings_real = cms.vdouble(*couplings.getcouplings(spin=0, WW=False, imag=False)),
                           HVVcouplings_imag = cms.vdouble(*couplings.getcouplings(spin=0, WW=False, imag=True)),
                           ZVVcouplings_real = cms.vdouble(*couplings.getcouplings(spin=1, imag=False)),
                           ZVVcouplings_imag = cms.vdouble(*couplings.getcouplings(spin=1, imag=True)),
                           GVVcouplings_real = cms.vdouble(*couplings.getcouplings(spin=2, gg=False, imag=False)),
                           GVVcouplings_imag = cms.vdouble(*couplings.getcouplings(spin=2, gg=False, imag=True)),
                           Gggcouplings_real = cms.vdouble(*couplings.getcouplings(spin=2, gg=True, imag=False)),
                           Gggcouplings_imag = cms.vdouble(*couplings.getcouplings(spin=2, gg=True, imag=True)),
                           reweightingtype = cms.string(REWEIGHTING_TYPE),
                           VVMode = cms.int32(int(VVMODE)),
                           VVDecayMode = cms.int32(int(VVDECAYMODE)),
                           AddLHEKinematics = cms.bool(ADDLHEKINEMATICS),
                           sampleProductionId = cms.string(SAMPLEPRODUCTIONID),
                           sampleProcessId = cms.string(SAMPLEPROCESSID),
                           Apply_K_NNLOQCD_ZZGG = cms.int32(int(APPLY_K_NNLOQCD_ZZGG)),
                           Apply_K_NNLOQCD_ZZQQB = cms.bool(APPLY_K_NNLOQCD_ZZQQB),
                           Apply_K_NLOEW_ZZQQB = cms.bool(APPLY_K_NLOEW_ZZQQB),
                           )

### Signal region
process.ZZTree = TreeSetup.clone()
process.ZZTree.channel = 'ZZ'

### Trees for control regions
process.CRZLLTree = TreeSetup.clone()
process.CRZLLTree.channel = 'ZLL'
process.CRZLLTree.CandCollection = 'ZLLCand'

### Trilepton CR, for fake rate
process.CRZLTree = TreeSetup.clone()
process.CRZLTree.channel = 'ZL'
process.CRZLTree.CandCollection = 'ZlCand'

### Loose electron candidates, signal region
process.ZZTreelooseEle = TreeSetup.clone()
process.ZZTreelooseEle.channel = 'ZZ'
process.ZZTreelooseEle.is_loose_ele_selection = cms.bool(True)
process.ZZTreelooseEle.CandCollection = 'ZZCandlooseEle'

#### Loose electron control regions
process.CRZLLTreelooseEle = TreeSetup.clone()
process.CRZLLTreelooseEle.channel = 'ZLL'
process.CRZLLTreelooseEle.CandCollection = 'ZLLCandlooseEle'

#### Loose electron Trilepton CR, for fake rate
process.CRZLTreelooseEle = TreeSetup.clone()
process.CRZLTreelooseEle.channel = 'ZL'
process.CRZLTreelooseEle.CandCollection = 'ZlCandlooseEle'


### TLE Signal region
process.ZZTreetle = TreeSetup.clone()
process.ZZTreetle.channel = 'ZZ'
process.ZZTreetle.is_loose_ele_selection = cms.bool(True)
process.ZZTreetle.CandCollection = 'ZZCandtle'

### TLE Trees for control regions only
process.CRZLLTreetle = TreeSetup.clone()
process.CRZLLTreetle.channel = 'ZLL'
process.CRZLLTreetle.CandCollection = 'ZLLCandtle'
#process.CRZLLTreetle.is_loose_ele_selection = cms.bool(True)
#process.CRZLLTreetle.CandCollection_regular = cms.untracked.string('ZLLCand')

### TLE Trilepton CR, for fake rate
process.CRZLTreetle = TreeSetup.clone()
process.CRZLTreetle.channel = 'ZL'
process.CRZLTreetle.CandCollection = 'ZlCandtle'
#process.CRZLTreetle.is_loose_ele_selection = cms.bool(True)
#process.CRZLTreetle.CandCollection_regular = cms.untracked.string('ZlCand')



# Debug
#Define candidates to be dumped
process.ZZFiltered = cms.EDFilter("PATCompositeCandidateRefSelector",
                                  src = cms.InputTag("ZZCand"),
                                  cut = cms.string("userFloat('isBestCand')")
                                  )
### Select only events with one such candidate
process.ZZSelection= cms.EDFilter("CandViewCountFilter",
                                  src = cms.InputTag("ZZFiltered"),
                                  minNumber = cms.uint32(1)
                                  )

### Select CR events
process.CRFiltered = cms.EDFilter("PATCompositeCandidateRefSelector",
                                  src = cms.InputTag("ZLLCand"),
                                  cut = cms.string("((userFloat('isBestCRZLLss')&&userFloat('CRZLLss')))||(userFloat('isBestCRZLLos_2P2F')&&userFloat('CRZLLos_2P2F'))||(userFloat('isBestCRZLLos_3P1F')&&userFloat('CRZLLos_3P1F'))")
                                  )
### Select only events with one such candidate
process.CRSelection= cms.EDFilter("CandViewCountFilter",
                                  src = cms.InputTag("CRFiltered"),
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
        Z     = cms.InputTag("ZCand"),
        ZZ  = cms.InputTag("ZZCand"),
        ZLL   =cms.InputTag("ZLLCand"),    # Starting point for all CRs
     )
)

if (PROCESS_CR or not IsMC):
    process.CRPath = cms.Path(process.CR)
    if (not IsMC):
        process.dump = cms.Path(process.ZZFiltered + process.ZZSelection + process.dumpUserData)
        process.dumpCR = cms.Path(process.CRFiltered + process.CRSelection + process.dumpUserData)
    process.trees = cms.EndPath( process.ZZTree + process.CRZLLTree + process.CRZLTree)
else:
#    process.CRPath = cms.Path(process.CRZl) #still needed by the plotter
    process.trees = cms.EndPath(process.ZZTree)

process.plots = cms.EndPath(process.PlotsZZ)


if (ADDLOOSEELE) :
    if (PROCESS_CR or not IsMC):
        process.CRPath += process.CRlooseEle
        process.trees += cms.Sequence( process.ZZTreelooseEle + process.CRZLLTreelooseEle + process.CRZLTreelooseEle)
        process.CRPath += process.CRtle
        process.trees += cms.Sequence( process.ZZTreetle + process.CRZLLTreetle + process.CRZLTreetle)
    else:
        process.trees += cms.Sequence(process.ZZTreelooseEle)
        process.trees += cms.Sequence(process.ZZTreetle)
