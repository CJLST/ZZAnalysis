### ----------------------------------------------------------------------
###
### Example analayzer
###
###----------------------------------------------------------------------


try:
    IsMC
except NameError:
    IsMC = True

try:
    LEPTON_SETUP
except NameError:
    LEPTON_SETUP = 2012 # define the set of effective areas, rho corrections, etc.

try:
    PD
except NameError:
    PD = ""             # "" for MC, "DoubleEle", "DoubleMu", or "MuEG" for data 

try:
    MCFILTER
except NameError:
    MCFILTER = ""


# Get absolute path
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------

#execfile(PyFilePath + "MasterPy/ZZ4lAnalysisPRL2011.py")         # 2011 reference analysis
execfile(PyFilePath + "MasterPy/ZZ4lAnalysis.py")         # 2012 reference analysis


### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------
process.source.fileNames = cms.untracked.vstring(

#        'root://cmsphys05//data/b/botta/V5_2_0/cmgTuple_H120Fall11_noSmearing.root' #Fall11 H120 for May, 21 synch exercise
#        'root://cmsphys05//data/b/botta/V5_4_0/cmgTuple_H120Fall11_noSmearing.root' #Fall11 H120 for FSR synch
         'root://cmsphys05//data/b/botta/V5_4_0/cmgTuple_H126Summer12.root' #Summer12 H126 for FSR synch        
    )


# from CMGTools.Production.datasetToSource import *
# process.source = datasetToSource(
#    'cmgtools',
#    '/GluGluToHToZZTo4L_M-120_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/V5/PAT_CMG_V5_2_0/',
#    'patTuple.*.root'
#    )

process.maxEvents.input = -1
#process.options.wantSummary = False


#Add my own cuts
#process.EEMMCand.flags +=
# SIPcut = cms.string("userFloat('SIP4')<4.")


### ----------------------------------------------------------------------
### Output root file (monitoring histograms)
### ----------------------------------------------------------------------
process.TFileService=cms.Service('TFileService',
                                fileName=cms.string('ZZ4lAnalysis.root')
                                )


### ----------------------------------------------------------------------
### Analyzers for Plots
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
                           sampleName = cms.string(SAMPLENAME),
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

### Control Region Plots

PlotCRSetup    = cms.EDAnalyzer("ZZ4lAnalyzerCR",
                                channel = cms.untracked.string('aChannel'),
                                candCollection = cms.untracked.string('aCand'),
                                isMC = cms.untracked.bool(IsMC),
                                sampleType = cms.int32(SAMPLE_TYPE),
                                setup = cms.int32(LEPTON_SETUP),
                                skimPaths = cms.vstring(SkimPaths),
                                PD = cms.string(PD),
                                MCFilterPath = cms.string(MCFILTER),
                                )

process.PlotsCRZLL = PlotCRSetup.clone()
process.PlotsCRZLL.channel               = "ZLL"
process.PlotsCRZLL.candCollection        = 'ZLLCand'

# process.PlotsCRZMM = PlotCRSetup.clone()
# process.PlotsCRZMM.channel               = "EEMM"
# process.PlotsCRZMM.candCollection        = 'CRZMM'

# process.PlotsCRZEE = PlotCRSetup.clone()
# process.PlotsCRZEE.channel               = "EEMM"
# process.PlotsCRZEE.candCollection        = 'CRZEE'

# process.PlotsCRZLLHiSIP = PlotCRSetup.clone()
# process.PlotsCRZLLHiSIP.channel          = "EEMM"
# process.PlotsCRZLLHiSIP.candCollection   = 'CRZLLHiSIP'

# process.PlotsCRZLLHiSIPMM = PlotCRSetup.clone()
# process.PlotsCRZLLHiSIPMM.channel        = "EEMM"
# process.PlotsCRZLLHiSIPMM.candCollection = 'CRZLLHiSIPMM'

# process.PlotsCRZLLHiSIPKin = PlotCRSetup.clone()
# process.PlotsCRZLLHiSIPKin.channel       = "EEMM"
# process.PlotsCRZLLHiSIPKin.candCollection= 'CRZLLHiSIPKin'

# process.PlotsCRMMEEss = PlotCRSetup.clone()
# process.PlotsCRMMEEss.channel            = "EEMM"
# process.PlotsCRMMEEss.candCollection     = 'CRMMEEss'

# process.PlotsCREEMMss = PlotCRSetup.clone()
# process.PlotsCREEMMss.channel            = "EEMM"
# process.PlotsCREEMMss.candCollection     = 'CREEMMss'

# process.PlotsCRMMMMss = PlotCRSetup.clone()
# process.PlotsCRMMMMss.channel            = "MMMM"
# process.PlotsCRMMMMss.candCollection     = 'CRMMMMss'

# process.PlotsCREEEEss = PlotCRSetup.clone()
# process.PlotsCREEEEss.channel            = "EEEE"
# process.PlotsCREEEEss.candCollection     = 'CREEEEss'

# process.PlotsCRMMEEos = PlotCRSetup.clone()
# process.PlotsCRMMEEos.channel            = "EEMM"
# process.PlotsCRMMEEos.candCollection     = 'CRMMEEos'

# process.PlotsCREEMMos = PlotCRSetup.clone()
# process.PlotsCREEMMos.channel            = "EEMM"
# process.PlotsCREEMMos.candCollection     = 'CREEMMos'

# process.PlotsCRMMMMos = PlotCRSetup.clone()
# process.PlotsCRMMMMos.channel            = "MMMM"
# process.PlotsCRMMMMos.candCollection     = 'CRMMMMos'

# process.PlotsCREEEEos = PlotCRSetup.clone()
# process.PlotsCREEEEos.channel            = "EEEE"
# process.PlotsCREEEEos.candCollection     = 'CREEEEos'


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
                                   CandCollection = cms.untracked.string('aCand'),
                                   fileName = cms.untracked.string('candTree'),
                                   isMC = cms.untracked.bool(IsMC),
                                   sampleType = cms.int32(SAMPLE_TYPE),
                                   setup = cms.int32(LEPTON_SETUP),
                                   skimPaths = cms.vstring(SkimPaths),
                                   PD = cms.string(PD),
                                   MCFilterPath = cms.string(MCFILTER),
                                   skipEmptyEvents = cms.bool(True),
                                   onlyBestCandidate =  cms.bool(False),
                                   sampleName = cms.string(SAMPLENAME),
                                   superMelaMass = cms.double(SUPERMELA_MASS),
                                   )

process.ZZ4muTree = TreeSetup.clone()
process.ZZ4muTree.channel = 'MMMM'
process.ZZ4muTree.CandCollection = 'MMMMCand'
process.ZZ4muTree.onlyBestCandidate = True

process.ZZ4eTree = TreeSetup.clone()
process.ZZ4eTree.channel = 'EEEE'
process.ZZ4eTree.CandCollection = 'EEEECand'
process.ZZ4muTree.onlyBestCandidate = True

process.ZZ2e2muTree = TreeSetup.clone()
process.ZZ2e2muTree.channel = 'EEMM'
process.ZZ2e2muTree.CandCollection = 'EEMMCand'
process.ZZ4muTree.onlyBestCandidate = True

### Trees for control regions only
TreeSetup.fileName = 'crTree' #No idea why we should change this...

process.CRZLLTree = TreeSetup.clone()
process.CRZLLTree.channel = 'ZLL'
process.CRZLLTree.CandCollection = 'ZLLCand'
process.CRZLLTree.onlyBestCandidate = True


# Trilepton CR, for fake rate
process.CRZLTree = TreeSetup.clone()
process.CRZLTree.channel = 'ZL'
process.CRZLTree.CandCollection = 'ZlCand'


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


process.dumpUserData =  cms.EDAnalyzer("dumpUserData",
     dumpTrigger = cms.untracked.bool(True),
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

if (not IsMC):
    process.CRPath = cms.Path(process.CR)
    process.dump = cms.Path(process.ZZFiltered + process.ZZSelection + process.dumpUserData)
    process.p = cms.EndPath( process.Plots4mu + process.Plots4e + process.Plots2e2mu + process.PlotsCRZLL )
    process.trees = cms.EndPath( process.ZZ4muTree * process.ZZ4eTree * process.ZZ2e2muTree * process.CRZLLTree + process.CRZLTree)
else:
    process.CRPath = cms.Path(process.CRZl) #still needed by the plotter
    process.p = cms.EndPath( process.Plots4mu + process.Plots4e + process.Plots2e2mu)
    process.trees = cms.EndPath( process.ZZ4muTree * process.ZZ4eTree * process.ZZ2e2muTree)
    

