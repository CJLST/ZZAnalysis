from ZZAnalysis.AnalysisStep.defaults import *
from ZZAnalysis.AnalysisStep.miscenums import *

### ----------------------------------------------------------------------
###
### Example analyzer
###
###----------------------------------------------------------------------

# Set defaults for variables used in this file (in case they are not defined by a caller script)
declareDefault("PD", "", globals()) # "" for MC, "DoubleEle", "DoubleMu", or "MuEG" for data
declareDefault("MCFILTER", "", globals())
declareDefault("XSEC", 1, globals())
declareDefault("GENXSEC", 1, globals())
declareDefault("GENBR", 1, globals())
declareDefault("PROCESS_CR", False, globals())
declareDefault("ADDZTREE", False, globals())

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
declareDefault("VVMODE", 1, globals())
declareDefault("VVDECAYMODE", 0, globals())
declareDefault("ADDLHEKINEMATICS", False, globals())

# K factors
declareDefault("APPLY_K_NNLOQCD_ZZGG", 0, globals()) # 0: Do not; 1: NNLO/LO; 2: NNLO/NLO; 3: NLO/LO
declareDefault("APPLY_K_NNLOQCD_ZZQQB", False, globals())
declareDefault("APPLY_K_NLOEW_ZZQQB", False, globals())

#failed events
declareDefault("SKIP_EMPTY_EVENTS", True, globals())
declareDefault("FAILED_TREE_LEVEL", 0, globals())

#ggF uncertainties for HTXS
declareDefault("APPLY_QCD_GGF_UNCERT", False, globals() )

if FAILED_TREE_LEVEL and not SKIP_EMPTY_EVENTS:
    raise ValueError(
                     "Inconsistent options: FAILED_TREE_LEVEL={}, SKIP_EMPTY_EVENTS={}\n"
                     "If you want to write a failed tree, set SKIP_EMPTY_EVENTS=True"
                     .format(FAILED_TREE_LEVEL, SKIP_EMPTY_EVENTS)
                    )

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
                           metSrc = metTag,
                           applyTrigger = cms.bool(APPLYTRIG), #Skip events failing required triggers. They are stored with sel<0 if set to false
                           applyTrigEff = cms.bool(False), #Add trigger efficiency as a weight, for samples where the trigger cannot be applied (obsoltete)
                           skipEmptyEvents = cms.bool(SKIP_EMPTY_EVENTS),
                           failedTreeLevel = cms.int32(FAILED_TREE_LEVEL),
                           sampleName = cms.string(SAMPLENAME),
									GenXSEC = cms.double(GENXSEC),
									GenBR = cms.double(GENBR),
                           dataTag=cms.string(DATA_TAG), #added for recognizing UL16 pre/post VFP

                           # MELA parameters
                           superMelaMass = cms.double(SUPERMELA_MASS),

                           # Reco MEs to pick from the candidate
                           recoProbabilities = cms.vstring(),

                           # LHE info. parameters
                           lheProbabilities = cms.vstring(),
                           xsec = cms.double(XSEC),
                           VVMode = cms.int32(VVMODE),
                           VVDecayMode = cms.int32(VVDECAYMODE),
                           AddLHEKinematics = cms.bool(ADDLHEKINEMATICS),
                           Apply_K_NNLOQCD_ZZGG = cms.int32(APPLY_K_NNLOQCD_ZZGG),
                           Apply_K_NNLOQCD_ZZQQB = cms.bool(APPLY_K_NNLOQCD_ZZQQB),
                           Apply_K_NLOEW_ZZQQB = cms.bool(APPLY_K_NLOEW_ZZQQB),
									Apply_QCD_GGF_UNCERT = cms.bool(APPLY_QCD_GGF_UNCERT),
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

#### Loose electron control region where Z1 is from tight real RSE and the regular leptons are fake
process.CRZLLTreeZ1RSE = TreeSetup.clone()
process.CRZLLTreeZ1RSE.channel = 'ZLL'
process.CRZLLTreeZ1RSE.CandCollection = 'ZLLCandZ1RSE'



#### Loose electron Trilepton CR, for fake rate
process.CRZLTreelooseEle = TreeSetup.clone()
process.CRZLTreelooseEle.channel = 'ZL'
process.CRZLTreelooseEle.CandCollection = 'ZlCandlooseEle'


#### TLE Signal region
#process.ZZTreetle = TreeSetup.clone()
#process.ZZTreetle.channel = 'ZZ'
#process.ZZTreetle.is_loose_ele_selection = cms.bool(True)
#process.ZZTreetle.CandCollection = 'ZZCandtle'
#
#### TLE Trees for control regions only
#process.CRZLLTreetle = TreeSetup.clone()
#process.CRZLLTreetle.channel = 'ZLL'
#process.CRZLLTreetle.CandCollection = 'ZLLCandtle'
##process.CRZLLTreetle.is_loose_ele_selection = cms.bool(True)
##process.CRZLLTreetle.CandCollection_regular = cms.untracked.string('ZLLCand')
#
#### TLE Trilepton CR, for fake rate
#process.CRZLTreetle = TreeSetup.clone()
#process.CRZLTreetle.channel = 'ZL'
#process.CRZLTreetle.CandCollection = 'ZlCandtle'
##process.CRZLTreetle.is_loose_ele_selection = cms.bool(True)
##process.CRZLTreetle.CandCollection_regular = cms.untracked.string('ZlCand')


### ----------------------------------------------------------------------
### Z tree
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
                               metSrc = metTag,
                               skipEmptyEvents = cms.bool(True),
                               sampleName = cms.string(SAMPLENAME),
                               xsec = cms.double(XSEC)
                               )



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
     ),
    jetSrc = cms.InputTag("cleanJets"),
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
        process.trees += cms.Sequence( process.ZZTreelooseEle + process.CRZLLTreelooseEle + process.CRZLTreelooseEle + process.CRZLLTreeZ1RSE)
        #process.CRPath += process.CRtle
        #process.trees += cms.Sequence( process.ZZTreetle + process.CRZLLTreetle + process.CRZLTreetle)
    else:
        process.trees += cms.Sequence(process.ZZTreelooseEle)
        #process.trees += cms.Sequence(process.ZZTreetle)

if (ADDZTREE) :
     process.trees += cms.Sequence(process.ZTree)
