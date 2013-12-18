### ----------------------------------------------------------------------
###
###  ZZ4lAnalyzer
###
###----------------------------------------------------------------------

### ----------------------------------------------------------------------
### Flags that need to be setted:
### ----------------------------------------------------------------------

IsMC = True 

PATVERSION = "V5"

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------
#execfile("MasterPy/ZZ4lAnalysisPRL2011.py") # 2011 PRL analysis
execfile("MasterPy/ZZ4lAnalysis.py")         # 2012 reference analysis

### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------
process.source.fileNames = cms.untracked.vstring(
    #'file:/data3/HZZ_Pattuple/CMG/V4_0_2/patTuple_PF2PAT_1_1_YgI.root'
    #'/store/cmst3/user/cmgtools/CMG/GluGluToHToZZTo4L_M-130_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V4_0_2/patTuple_PF2PAT_10_1_Ryi.root'
    '/store/cmst3/user/cmgtools/CMG/GluGluToHToZZTo4L_M-130_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/V5/PAT_CMG_V5_2_0/patTuple_1.root'
   #'root://lxcms00//data3/HZZ_Pattuple/CMG/V5_1_0/patTuple.root' # Reference file for May, 10 sync excercise
    #'root://cmsphys05//data/b/botta/V5_2_0/CMGTools/CMSSW_4_4_4/src/CMGTools/Common/prod/patTuple.root' # DY file for the V5_2_0 test
   )

process.maxEvents.input = 100
#process.options.wantSummary = False

process.CRPath = cms.Path(process.CR)

### ----------------------------------------------------------------------
### Output root file (monitoring histograms)
### ----------------------------------------------------------------------

process.TFileService=cms.Service('TFileService',
                                 fileName=cms.string('ZZ4lAnalysisTree.root')
                                 )

process.ZZ4muTree = cms.EDAnalyzer("HZZ4lNtupleMaker",
                                   channel = cms.untracked.string('MMMM'),
                                   CandCollection = cms.untracked.string('MMMMCand'),
                                   fileName = cms.untracked.string('candTree'),
                                   isMC = cms.untracked.bool(IsMC),
                                   setup = cms.int32(LEPTON_SETUP),
                                   skimPaths = cms.vstring(SkimPaths),
                                   PD = cms.string(""),
                                   MCFilterPath = cms.string(""),                                    
                                   )

process.ZZ4eTree = cms.EDAnalyzer("HZZ4lNtupleMaker",
                                  channel = cms.untracked.string('EEEE'),
                                  CandCollection = cms.untracked.string('EEEECand'),
                                  fileName = cms.untracked.string('candTree'),
                                  isMC = cms.untracked.bool(IsMC),
                                  setup = cms.int32(LEPTON_SETUP),
                                  skimPaths = cms.vstring(SkimPaths),
                                  PD = cms.string(""),
                                  MCFilterPath = cms.string(""),                                    
                                  )

process.ZZ2e2muTree = cms.EDAnalyzer("HZZ4lNtupleMaker",
                                     channel = cms.untracked.string('EEMM'),
                                     CandCollection = cms.untracked.string('EEMMCand'),
                                     fileName = cms.untracked.string('candTree'),
                                     setup = cms.int32(LEPTON_SETUP),
                                     isMC = cms.untracked.bool(IsMC),
                                     skimPaths = cms.vstring(SkimPaths),
                                     PD = cms.string(""),
                                     MCFilterPath = cms.string(""),                                    
                                     )

process.CRMMEEssTree = cms.EDAnalyzer("HZZ4lNtupleMaker",
                                      channel = cms.untracked.string('EEMM'),
                                      CandCollection = cms.untracked.string('CRMMEEss'),
                                      fileName = cms.untracked.string('crTree'),
                                      setup = cms.int32(LEPTON_SETUP),
                                      isMC = cms.untracked.bool(IsMC),
                                      skimPaths = cms.vstring(SkimPaths),                                     
                                      PD = cms.string(""),
                                      MCFilterPath = cms.string(""),
                                      )



process.CREEMMssTree = cms.EDAnalyzer("HZZ4lNtupleMaker",
                                      channel = cms.untracked.string('EEMM'),
                                      CandCollection = cms.untracked.string('CREEMMss'),
                                      fileName = cms.untracked.string('crTree'),
                                      setup = cms.int32(LEPTON_SETUP),
                                      isMC = cms.untracked.bool(IsMC),
                                      skimPaths = cms.vstring(SkimPaths),                                     
                                      PD = cms.string(""),
                                      MCFilterPath = cms.string(""),
                                      )



process.CRMMMMssTree = cms.EDAnalyzer("HZZ4lNtupleMaker",
                                      channel = cms.untracked.string('MMMM'),
                                      CandCollection = cms.untracked.string('CRMMMMss'),
                                      fileName = cms.untracked.string('crTree'),
                                      setup = cms.int32(LEPTON_SETUP),
                                      isMC = cms.untracked.bool(IsMC),
                                      skimPaths = cms.vstring(SkimPaths),                                     
                                      PD = cms.string(""),
                                      MCFilterPath = cms.string(""),
                                      )



process.CREEEEssTree = cms.EDAnalyzer("HZZ4lNtupleMaker",
                                      channel = cms.untracked.string('EEEE'),
                                      CandCollection = cms.untracked.string('CREEEEss'),
                                      fileName = cms.untracked.string('crTree'),
                                      setup = cms.int32(LEPTON_SETUP),
                                      isMC = cms.untracked.bool(IsMC),
                                      skimPaths = cms.vstring(SkimPaths),                                     
                                      PD = cms.string(""),
                                      MCFilterPath = cms.string(""),
                                      )



process.CRMMEEosTree = cms.EDAnalyzer("HZZ4lNtupleMaker",
                                      channel = cms.untracked.string('EEMM'),
                                      CandCollection = cms.untracked.string('CRMMEEos'),
                                      fileName = cms.untracked.string('crTree'),
                                      setup = cms.int32(LEPTON_SETUP),
                                      isMC = cms.untracked.bool(IsMC),
                                      skimPaths = cms.vstring(SkimPaths),                                     
                                      PD = cms.string(""),
                                      MCFilterPath = cms.string(""),
                                      )



process.CREEMMosTree = cms.EDAnalyzer("HZZ4lNtupleMaker",
                                      channel = cms.untracked.string('EEMM'),
                                      CandCollection = cms.untracked.string('CREEMMos'),
                                      fileName = cms.untracked.string('crTree'),
                                      setup = cms.int32(LEPTON_SETUP),
                                      isMC = cms.untracked.bool(IsMC),
                                      skimPaths = cms.vstring(SkimPaths),                                     
                                      PD = cms.string(""),
                                      MCFilterPath = cms.string(""),
                                      )


process.CRMMMMosTree = cms.EDAnalyzer("HZZ4lNtupleMaker",
                                      channel = cms.untracked.string('MMMM'),
                                      CandCollection = cms.untracked.string('CRMMMMos'),
                                      fileName = cms.untracked.string('crTree'),
                                      setup = cms.int32(LEPTON_SETUP),
                                      isMC = cms.untracked.bool(IsMC),
                                      skimPaths = cms.vstring(SkimPaths),                                     
                                      PD = cms.string(""),
                                      MCFilterPath = cms.string(""),
                                      )



process.CREEEEosTree = cms.EDAnalyzer("HZZ4lNtupleMaker",
                                      channel = cms.untracked.string('EEEE'),
                                      CandCollection = cms.untracked.string('CREEEEos'),
                                      fileName = cms.untracked.string('crTree'),
                                      setup = cms.int32(LEPTON_SETUP),
                                      isMC = cms.untracked.bool(IsMC),
                                      skimPaths = cms.vstring(SkimPaths),                                     
                                      PD = cms.string(""),
                                      MCFilterPath = cms.string(""),
                                      )

process.CRZLLHiSIPTree = cms.EDAnalyzer("HZZ4lNtupleMaker",
                                      channel = cms.untracked.string('MMMM'),
                                      CandCollection = cms.untracked.string('CRZLLHiSIP'),
                                      fileName = cms.untracked.string('crTree'),
                                      setup = cms.int32(LEPTON_SETUP),
                                      isMC = cms.untracked.bool(IsMC),
                                      skimPaths = cms.vstring(SkimPaths),                                     
                                      PD = cms.string(""),
                                      MCFilterPath = cms.string(""),
                                      )


process.CRZLLTree = cms.EDAnalyzer("HZZ4lNtupleMaker",
                                      channel = cms.untracked.string('MMMM'),
                                      CandCollection = cms.untracked.string('CRZLL'),
                                      fileName = cms.untracked.string('crTree'),
                                      setup = cms.int32(LEPTON_SETUP),
                                      isMC = cms.untracked.bool(IsMC),
                                      skimPaths = cms.vstring(SkimPaths),                                     
                                      PD = cms.string(""),
                                      MCFilterPath = cms.string(""),
                                      )


                                       

process.trees = cms.EndPath( process.ZZ4muTree * process.ZZ4eTree * process.ZZ2e2muTree * process.CRMMMMssTree * process.CREEEEssTree * process.CREEMMssTree * process.CRMMEEssTree * process.CRMMMMosTree * process.CREEEEosTree * process.CREEMMosTree * process.CRMMEEosTree * process.CRZLLHiSIPTree * process.CRZLLTree )

