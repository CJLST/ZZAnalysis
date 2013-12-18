### ----------------------------------------------------------------------
###
### Simple analyzer to inspect user data attached to candidates.
###
###----------------------------------------------------------------------

#ELECORRTYPE = "Summer12_DR53X_HCP2012"
APPLYELEREGRESSION = False
APPLYMUCORR = False

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------
execfile("MasterPy/ZZ4lAnalysis.py")


### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------
process.source.fileNames = cms.untracked.vstring(
#    'file:/data3/HZZ_Pattuple/CMG/V4_0_2/patTuple_PF2PAT_1_1_YgI.root'
#    'root://cmsphys05//data/b/botta/V5_2_0/CMGTools/CMSSW_4_4_4/src/CMGTools/Common/prod/cmgTuple.root' # DY file for the V5_2_0 test
#    '/store/cmst3/user/cmgtools/CMG/GluGluToHToZZTo4L_M-130_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/V5/PAT_CMG_V5_2_0/cmgTuple_1.root'
#    'root://cmsphys05//data/b/botta/V5_4_0/cmgTuple_H120Fall11_noSmearing.root' #Fall11 H120 for FSR synch
#    'root://lxcms00//data3/2012/HZZ_cmgTuple/synchHCP1/H125_53X.root' #HCP sync file
    'root://lxcms00//data3/2012/HZZ_cmgTuple/synchHCP1/H125_53X_pt3ele.root' #Version with looser ele pT filter
#    'root://lxcms00//data3/2012/HZZ_cmgTuple/synchHCP1/VBFH125_53X.root'  #HCP VBF sync file

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
#        Zmm   = cms.InputTag("MMCand"),
#        Zee   = cms.InputTag("EECand"),
#        Zll   = cms.InputTag("ZCand"),
#        LL    = cms.InputTag("LLCand"),
        MMMM  = cms.InputTag("MMMMCand"),
        EEEE  = cms.InputTag("EEEECand"),
        EEMM  = cms.InputTag("EEMMCand"),
#        ZLL   =cms.InputTag("ZLLCand"),    # Starting point for all CRs
#        LLLL  = cms.InputTag("LLLLCand"), # For RFC studies
     )
)

process.analysis = cms.EndPath(process.dumpUserData)





