import FWCore.ParameterSet.Config as cms
from ZZAnalysis.AnalysisStep.defaults import *
import os, sys

process = cms.Process("ZZ")

### ----------------------------------------------------------------------
### Flags that need to be set
### ----------------------------------------------------------------------

# Set defaults for variables used in this file (in case they are not defined by a caller script)

declareDefault("IsMC", True, globals())

# Set of effective areas, rho corrections, etc. (can be 2011, 2012, 2015 or 2016)
declareDefault("LEPTON_SETUP", 2016, globals())

# Flag that reflects the actual sqrts of the sample (can be 2011, 2012, 2015 or 2016)
# Can differ from SAMPLE_TYPE for samples that are rescaled to a different sqrts.
declareDefault("SAMPLE_TYPE", LEPTON_SETUP, globals())

# Control global tag to be used (for data or rereco)
declareDefault("DATA_TAG", "ReReco", globals())

#Optional name of the sample/dataset being analyzed
declareDefault("SAMPLENAME", "", globals())

#Type of electron scale correction/smearing: "None", "RunII"
declareDefault("ELECORRTYPE", "RunII", globals())

#Apply muon scale correction
declareDefault("APPLYMUCORR", True, globals())

#Reapply JEC
declareDefault("APPLYJEC", True, globals())

#Apply JER
declareDefault("APPLYJER", True, globals())

#Recorrect MET
declareDefault("RECORRECTMET", True, globals())

#FSR mode
declareDefault("FSRMODE", "RunII", globals())

#Bunch spacing (can be 25 or 50)
declareDefault("BUNCH_SPACING", 25, globals())

#Mass used for SuperMELA
declareDefault("SUPERMELA_MASS", 125, globals())

#Selection flow strategy
declareDefault("SELSETUP", "allCutsAtOncePlusSmart", globals())

#Best candidate comparator (see interface/Comparators.h)
declareDefault("BESTCANDCOMPARATOR", "byBestKD", globals())

# Set to True to make candidates with the full combinatorial of loose leptons (for debug; much slower)
declareDefault("KEEPLOOSECOMB", False, globals())

# Activate the Z kinematic refit (very slow)
declareDefault("KINREFIT", False, globals())

# Activate paths for loose electron categories
declareDefault("ADDLOOSEELE", False, globals())

# Activate trigger paths in MC; note that for 2016, only reHLT samples have the correct triggers!!!
declareDefault("APPLYTRIG", True, globals())

# CMSSW version 8X or 9X
CMSSW_VERSION = os.environ['CMSSW_VERSION']
CMSSWVERSION = int(CMSSW_VERSION.split("_")[1])

if SELSETUP=="Legacy" and not BESTCANDCOMPARATOR=="byBestZ1bestZ2":
    print "WARNING: In ZZ4lAnalysis.py the SELSETUP=\"Legacy\" flag is meant to reproduce the Legacy results, ignoring the setting of the BESTCANDCOMPARATOR: ",BESTCANDCOMPARATOR
    BESTCANDCOMPARATOR = "byBestZ1bestZ2"


# The isolation cuts for electrons and muons. FIXME: there is an hardcoded instance of these values in src/LeptonIsoHelper.cc !!
ELEISOCUT = 99999. # [FIXME] Remove isolation cuts from the code completely
MUISOCUT = 0.35

### ----------------------------------------------------------------------
### Set the GT
### ----------------------------------------------------------------------
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag


if (SAMPLE_TYPE == 2011) :
    if IsMC:
        process.GlobalTag.globaltag = 'START44_V13::All'
    else:
#        process.GlobalTag.globaltag = 'GR_R_44_V5::All'
        process.GlobalTag.globaltag = 'GR_R_44_V15C::All'
#        process.GlobalTag.connect = "sqlite_file:/afs/cern.ch/user/a/alcaprod/public/Alca/GlobalTag/GR_R_44_V15C.db"

elif (SAMPLE_TYPE == 2012) :
    if IsMC:
#        process.GlobalTag.globaltag = 'START53_V7G::All'   #for 53X MC
#        process.GlobalTag.globaltag = 'START53_V23::All'   #for 53X MC, updated JEC on top of pattuples produced with START53_V7G
        process.GlobalTag.globaltag = 'PLS170_V6AN1::All'
    else:
#        process.GlobalTag.globaltag = 'FT_53_V21_AN3::All' # For 22Jan2013
#        process.GlobalTag.globaltag = 'FT_53_V21_AN4::All' # For 22Jan2013, updated JEC on top of pattuples produced with FT_53_V21_AN3
        process.GlobalTag.globaltag = 'GR_70_V2_AN1::All'

elif (SAMPLE_TYPE == 2015) :
    if IsMC:
        process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
    else:
        process.GlobalTag.globaltag = '76X_dataRun2_v15'

### ICHEP 2016, superseeded
# elif (SAMPLE_TYPE == 2016):
#     from Configuration.AlCa.GlobalTag import GlobalTag
#     if IsMC:
#         process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_miniAODv2_v1', '')
#     else:
#         process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_ICHEP16JEC_v0', '')

elif (SAMPLE_TYPE == 2016):
    if IsMC:
        process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mcRun2_asymptotic_v3', '')
    else:
        process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v10', '')

elif (SAMPLE_TYPE == 2017):
    if IsMC:
        process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v17', '')
    else:
        process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v11', '')

elif (SAMPLE_TYPE == 2018):
    if IsMC:
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v18', '')
    else:
        process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v4', '')



print '\t',process.GlobalTag.globaltag

### ----------------------------------------------------------------------
### Standard stuff
### ----------------------------------------------------------------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


### ----------------------------------------------------------------------
### Source
### ----------------------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      '/store/cmst3/user/cmgtools/CMG/GluGluToHToZZTo4L_M-130_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/V5/PAT_CMG_V5_2_0/patTuple_1.root'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


### ----------------------------------------------------------------------
### Trigger bit Requests
### ----------------------------------------------------------------------
import HLTrigger.HLTfilters.hltHighLevel_cfi

process.hltFilterDiMu  = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterDiEle = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterMuEle = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterMuEle2 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterMuEle3 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterTriEle = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterTriMu  = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterSingleEle = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterSingleMu  = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterDiMu.TriggerResultsTag  = cms.InputTag("TriggerResults","","HLT")
process.hltFilterDiEle.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterMuEle.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterMuEle2.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterMuEle3.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterTriEle.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterTriMu.TriggerResultsTag  = cms.InputTag("TriggerResults","","HLT")
process.hltFilterSingleEle.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterSingleMu.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterDiMu.throw  = cms.bool(False) #FIXME: beware of this!
process.hltFilterDiEle.throw = cms.bool(False) #FIXME: beware of this!
process.hltFilterMuEle.throw = cms.bool(False) #FIXME: beware of this!
process.hltFilterMuEle2.throw = cms.bool(False) #FIXME: beware of this!
process.hltFilterMuEle3.throw = cms.bool(False) #FIXME: beware of this!
process.hltFilterTriEle.throw = cms.bool(False) #FIXME: beware of this!
process.hltFilterTriMu.throw  = cms.bool(False) #FIXME: beware of this!
process.hltFilterSingleEle.throw = cms.bool(False) #FIXME: beware of this!
process.hltFilterSingleMu.throw  = cms.bool(False) #FIXME: beware of this!

### 2016 triggers - final
if (LEPTON_SETUP == 2016):
   process.hltFilterDiEle.HLTPaths = ["HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*"]
   process.hltFilterDiMu.HLTPaths = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*"]
   process.hltFilterMuEle.HLTPaths = ["HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v*","HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v*"]
   process.hltFilterTriEle.HLTPaths = ["HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*"]
   process.hltFilterTriMu.HLTPaths = ["HLT_TripleMu_12_10_5_v*"]
   process.hltFilterSingleEle.HLTPaths = ["HLT_Ele25_eta2p1_WPTight_Gsf_v*","HLT_Ele27_WPTight_Gsf_v*","HLT_Ele27_eta2p1_WPLoose_Gsf_v*", "HLT_Ele32_eta2p1_WPTight_Gsf_v*"]
   process.hltFilterSingleMu.HLTPaths = ["HLT_IsoMu20_v*","HLT_IsoTkMu20_v*","HLT_IsoMu22_v*","HLT_IsoTkMu22_v*","HLT_IsoMu24_v*","HLT_IsoTkMu24_v*"]

   process.triggerTriEle = cms.Path(process.hltFilterTriEle)
   process.triggerTriMu  = cms.Path(process.hltFilterTriMu )
   process.triggerSingleEle = cms.Path(process.hltFilterSingleEle)
   process.triggerSingleMu  = cms.Path(process.hltFilterSingleMu )

### 2017 triggers - final
elif (LEPTON_SETUP == 2017):
   process.hltFilterDiEle.HLTPaths = ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*","HLT_DoubleEle33_CaloIdL_MW_v*"]
   process.hltFilterDiMu.HLTPaths = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*"]
   process.hltFilterMuEle.HLTPaths = ["HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v*","HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v*","HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v*"]
   process.hltFilterTriEle.HLTPaths = ["HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*"]
   process.hltFilterTriMu.HLTPaths = ["HLT_TripleMu_10_5_5_DZ_v*","HLT_TripleMu_12_10_5_v*"]
   process.hltFilterSingleEle.HLTPaths = ["HLT_Ele35_WPTight_Gsf_v*","HLT_Ele38_WPTight_Gsf_v*","HLT_Ele40_WPTight_Gsf_v*"]
   process.hltFilterSingleMu.HLTPaths = ["HLT_IsoMu27_v*"]

   process.triggerTriEle = cms.Path(process.hltFilterTriEle)
   process.triggerTriMu  = cms.Path(process.hltFilterTriMu )
   process.triggerSingleEle = cms.Path(process.hltFilterSingleEle)
   process.triggerSingleMu  = cms.Path(process.hltFilterSingleMu )

### 2018 triggers - FIXME: to be updated (26/6/18)
elif (LEPTON_SETUP == 2018):
   process.hltFilterDiEle.HLTPaths = ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*","HLT_DoubleEle25_CaloIdL_MW_v*"]
   process.hltFilterDiMu.HLTPaths = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*"]
   process.hltFilterMuEle.HLTPaths = ["HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v*","HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v*"]
   process.hltFilterTriEle.HLTPaths = [""]
   process.hltFilterTriMu.HLTPaths = ["HLT_TripleMu_10_5_5_DZ_v*","HLT_TripleMu_12_10_5_v*"]
   process.hltFilterSingleEle.HLTPaths = ["HLT_Ele32_WPTight_Gsf_v*"]
   process.hltFilterSingleMu.HLTPaths = ["HLT_IsoMu24_v*"]

   process.triggerTriEle = cms.Path(process.hltFilterTriEle)
   process.triggerTriMu  = cms.Path(process.hltFilterTriMu )
   process.triggerSingleEle = cms.Path(process.hltFilterSingleEle)
   process.triggerSingleMu  = cms.Path(process.hltFilterSingleMu )



process.triggerDiMu   = cms.Path(process.hltFilterDiMu)
process.triggerDiEle  = cms.Path(process.hltFilterDiEle)
process.triggerMuEle  = cms.Path(process.hltFilterMuEle)


### ----------------------------------------------------------------------
### MET FILTERS
### ----------------------------------------------------------------------
#process.METFilters  = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
#process.METFilters.TriggerResultsTag  = cms.InputTag("TriggerResults","","RECO")
#if (IsMC):
#   process.METFilters.TriggerResultsTag  = cms.InputTag("TriggerResults","","PAT")
#
#if (LEPTON_SETUP == 2017):#MET Filters available in miniAOD as described here https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
#   if (IsMC):
#      process.METFilters.HLTPaths = ["Flag_goodVertices","Flag_globalSuperTightHalo2016Filter","Flag_HBHENoiseFilter","Flag_HBHENoiseIsoFilter","Flag_EcalDeadCellTriggerPrimitiveFilter","Flag_BadPFMuonFilter","Flag_BadChargedCandidateFilter"]
#   else:
#      process.METFilters.HLTPaths = ["Flag_goodVertices","Flag_globalSuperTightHalo2016Filter","Flag_HBHENoiseFilter","Flag_HBHENoiseIsoFilter","Flag_EcalDeadCellTriggerPrimitiveFilter","Flag_BadPFMuonFilter","Flag_BadChargedCandidateFilter","Flag_eeBadScFilter"]
#
#process.triggerMETFilters = cms.Path(process.METFilters)

### ----------------------------------------------------------------------
### MC Filters and tools
### ----------------------------------------------------------------------

process.heavyflavorfilter = cms.EDFilter('HeavyFlavorFilter2',
#                                 src= cms.InputTag("genParticles"), # genParticles available only in PAT
                                 src= cms.InputTag("prunedGenParticles"),
                                 status2 = cms.bool(True),
                                 status3 = cms.bool(False),
                                 hDaughterVeto = cms.bool(False),
                                 zDaughterVeto = cms.bool(True),
                                 ptcut=cms.double(0)
                                 )


process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.drawTree = cms.EDAnalyzer("ParticleTreeDrawer",
                                   src = cms.InputTag("prunedGenParticles"),
                                   printP4 = cms.untracked.bool(False),
                                   printPtEtaPhi = cms.untracked.bool(False),
                                   printVertex = cms.untracked.bool(False),
                                   printStatus = cms.untracked.bool(True),
                                   printIndex = cms.untracked.bool(False) )


process.printTree = cms.EDAnalyzer("ParticleListDrawer",
                                   maxEventsToPrint = cms.untracked.int32(-1),
                                   printVertex = cms.untracked.bool(False),
                                   src = cms.InputTag("prunedGenParticles")
                                   )


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   calibratedPatElectrons = cms.PSet(
                                                       initialSeed = cms.untracked.uint32(1),
                                                       engineName = cms.untracked.string('TRandom3')
                                                       ),
                                                   )

# FIXME Add total kinematics filter for MC

process.goodPrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlineSlimmedPrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(True),
)


### ----------------------------------------------------------------------
### HTXS categorisation
### ----------------------------------------------------------------------
if(IsMC):
   process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
   process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
															  inputPruned = cms.InputTag("prunedGenParticles"),
															  inputPacked = cms.InputTag("packedGenParticles"),
															  )
   process.myGenerator = cms.EDProducer("GenParticles2HepMCConverter",
													 genParticles = cms.InputTag("mergedGenParticles"),
													 genEventInfo = cms.InputTag("generator"),
													 signalParticlePdgIds = cms.vint32(25), ## for the Higgs analysis
													 )
   process.rivetProducerHTXS = cms.EDProducer('HTXSRivetProducer',
															 HepMCCollection = cms.InputTag('myGenerator','unsmeared'),
															 LHERunInfo = cms.InputTag('externalLHEProducer'),
															 ProductionMode = cms.string('AUTO'),
															 )
   process.htxs = cms.Path(process.mergedGenParticles*process.myGenerator*process.rivetProducerHTXS)

### ----------------------------------------------------------------------
### ----------------------------------------------------------------------
### Loose lepton selection + cleaning + embeddding of user data
### ----------------------------------------------------------------------
### ----------------------------------------------------------------------

SIP =  "userFloat('SIP')<4"
GOODLEPTON = "userFloat('ID') && " + SIP  # Lepton passing tight ID + SIP [ISO is asked AFTER FSR!!!]

if (LEPTON_SETUP >= 2016) : # Run II
    TIGHTMUON = "userFloat('isPFMuon') || (userFloat('isTrackerHighPtMuon') && pt>200)"
else:
    TIGHTMUON = "userFloat('isPFMuon')" # Run I



#------- MUONS -------

#--- Mu e-scale corrections (Rochester)
# process.calibratedMuons =  cms.EDProducer("RochesterPATMuonCorrector",
#                                                src = cms.InputTag("patMuonsWithTrigger")
#                                               )

#--- Mu e-scale corrections (MuScleFit)
#process.calibratedMuons = cms.EDProducer("MuScleFitPATMuonCorrector",
#                         src = cms.InputTag("slimmedMuons"),
#                         debug = cms.bool(False),
#                         identifier = cms.string("Summer12_DR53X_smearReReco"),
#                         applySmearing = cms.bool(IsMC),
#                         fakeSmearing = cms.bool(False)
#                         )

#--- Mu e-scale corrections (KalmanMuonCalibrator, 2015)
process.calibratedMuons = cms.EDProducer("KalmanPATMuonCorrector",
                                         src = cms.InputTag("slimmedMuons"),
                                         identifier = cms.string(""),
                                         isMC = cms.bool(IsMC),
                                         isSynchronization = cms.bool(False),
                                         )

#--- Set correct identifier for muon corrections
if LEPTON_SETUP == 2016: # Rochester corrections for 2016 data
      process.calibratedMuons = cms.EDProducer("RochesterPATMuonCorrector",
                                            src = cms.InputTag("slimmedMuons"),
                                            identifier = cms.string("RoccoR2016"),
                                            isMC = cms.bool(IsMC),
                                            isSynchronization = cms.bool(False),
                                            )
elif LEPTON_SETUP == 2017:# Rochester corrections for 2017 data
     process.calibratedMuons = cms.EDProducer("RochesterPATMuonCorrector",
                                         src = cms.InputTag("slimmedMuons"),
                                         identifier = cms.string("RoccoR2017"),
                                         isMC = cms.bool(IsMC),
                                         isSynchronization = cms.bool(False),
                                         )
elif LEPTON_SETUP == 2018:# Rochester corrections for 2018 data
     process.calibratedMuons = cms.EDProducer("RochesterPATMuonCorrector",
                                         src = cms.InputTag("slimmedMuons"),
                                         identifier = cms.string("RoccoR2018"),
                                         isMC = cms.bool(IsMC),
                                         isSynchronization = cms.bool(False),
                                         )

else:
    if APPLYMUCORR:
        print "APPLYMUCORR not configured for LEPTON_SETUP =", LEPTON_SETUP
        sys.exit()

#--- Mu Ghost cleaning
process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
                                   src = cms.InputTag("calibratedMuons"),
                                   preselection = cms.string("track.isNonnull"),
                                   passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
                                   fractionOfSharedSegments = cms.double(0.499))


process.bareSoftMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("cleanedMu"),
    cut = cms.string("pt>5 && abs(eta)<2.4 && (isGlobalMuon || (isTrackerMuon && numberOfMatches>0)) && muonBestTrackType!=2")
#    Lowering pT cuts
#    cut = cms.string("(isGlobalMuon || (isTrackerMuon && numberOfMatches>0)) &&" +
#                     "pt>3 && p>3.5 && abs(eta)<2.4")
)


# MC matching. As the genParticles are no more available in cmg, we re-match with genParticlesPruned.
process.muonMatch = cms.EDProducer("MCMatcher", # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                   src     = cms.InputTag("softMuons"), # RECO objects to match
                                   matched = cms.InputTag("prunedGenParticles"),   # mc-truth particle collection
                                   mcPdgId     = cms.vint32(13), # one or more PDG ID (13 = muon); absolute values (see below)
                                   checkCharge = cms.bool(True), # True = require RECO and MC objects to have the same charge
                                   mcStatus = cms.vint32(1),     # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                                   maxDeltaR = cms.double(0.5),  # Minimum deltaR for the match
                                   maxDPtRel = cms.double(0.5),  # Minimum deltaPt/Pt for the match
                                   resolveAmbiguities = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
                                   resolveByMatchQuality = cms.bool(False), # False = just match input in order; True = pick lowest deltaR pair first
                                   )

process.softMuons = cms.EDProducer("MuFiller",
    src = cms.InputTag("bareSoftMuons"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    cut = cms.string("userFloat('dxy')<0.5 && userFloat('dz')<1."),
    TriggerResults = cms.InputTag('TriggerResults','','HLT'),
    flags = cms.PSet(
        ID = cms.string(TIGHTMUON), # tight muon ID
        isSIP = cms.string(SIP),
        isGood = cms.string(GOODLEPTON),
        isIsoFSRUncorr  = cms.string("userFloat('combRelIsoPF')<" + str(MUISOCUT)),
#       Note: passCombRelIsoPFFSRCorr is currently set in LeptonPhotonMatcher for new FSR strategy; in ZZCandidateFiller for the old one
    )
)


if APPLYMUCORR :
    process.muons =  cms.Sequence(process.calibratedMuons + process.cleanedMu + process.bareSoftMuons + process.softMuons)
else:
    process.cleanedMu.src = cms.InputTag("slimmedMuons")
    process.muons =  cms.Sequence(process.cleanedMu + process.bareSoftMuons + process.softMuons)




#------- ELECTRONS -------

#--- Run2 electron momentum scale and resolution corrections

process.selectedSlimmedElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string("pt>5 && abs(eta)<2.5")
)

if (LEPTON_SETUP == 2016):
   from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
   setupEgammaPostRecoSeq(process,
                          runEnergyCorrections=True, #corrections by default are fine so no need to re-run
                          era='2016-Legacy')

if (LEPTON_SETUP == 2017):
   from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
   setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=True, #corrections by default are fine so no need to re-run
                       era='2017-Nov17ReReco')

## Preliminary Moriond 19 corrections from https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPostRecoRecipes#2018_Preliminary_Energy_Correcti
if (LEPTON_SETUP == 2018):
    process.load("RecoEgamma.EgammaTools.calibratedEgammas_cff")
    process.calibratedPatElectrons.correctionFile = "EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2018_Step2Closure_CoarseEtaR9Gain"

#--- Set up electron ID (VID framework)
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format to be DataFormat.MiniAOD, as appropriate
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff','RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff']
# add them to the VID producer
for idmod in my_id_modules:
	setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
# and don't forget to add the producer 'process.egmGsfElectronIDSequence' to the path, i.e. process.electrons


process.bareSoftElectrons = cms.EDFilter("PATElectronRefSelector",
   src = cms.InputTag("calibratedPatElectrons"),
   cut = cms.string("pt>7 && abs(eta)<2.5")
   )

process.softElectrons = cms.EDProducer("EleFiller",
   src    = cms.InputTag("bareSoftElectrons"),
   sampleType = cms.int32(SAMPLE_TYPE),
   setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
   cut = cms.string("userFloat('dxy')<0.5 && userFloat('dz')<1"),
   flags = cms.PSet(
        ID = cms.string("userFloat('isBDT')"),
        isSIP = cms.string(SIP),
        isGood = cms.string(GOODLEPTON),
        isIsoFSRUncorr  = cms.string("userFloat('combRelIsoPF')<"+str(ELEISOCUT))
#       Note: passCombRelIsoPFFSRCorr is currently set in LeptonPhotonMatcher for new FSR strategy; in ZZCandidateFiller for the old one
        ),
        mvaValuesMap = cms.InputTag(""),
   	correctionFile = cms.string(""),
   )

process.softElectrons.correctionFile = process.calibratedPatElectrons.correctionFile

if (LEPTON_SETUP == 2016):
   process.softElectrons.mvaValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV2RawValues")
#94X BDT with ID and Isolation
if (LEPTON_SETUP == 2017):
   process.softElectrons.mvaValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV2RawValues")

#2018 version of BDT with ID and Isolation
if (LEPTON_SETUP == 2018):
   process.softElectrons.mvaValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV2RawValues")

#process.electrons = cms.Sequence(process.selectedSlimmedElectrons + process.calibratedPatElectrons + process.egmGsfElectronIDSequence + process.bareSoftElectrons + process.softElectrons) # (use this version when running VID)
#process.electrons = cms.Sequence(process.selectedSlimmedElectrons + process.calibratedPatElectrons + process.bareSoftElectrons + process.softElectrons) # (use this version without VID)

# Handle special cases
if (ELECORRTYPE == "None" ):   # No correction at all. Skip correction modules.
   process.bareSoftElectrons.src = cms.InputTag('selectedSlimmedElectrons')
   process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('selectedSlimmedElectrons')
   process.electronMVAValueMapProducer.srcMiniAOD =  cms.InputTag('selectedSlimmedElectrons')

   process.electrons = cms.Sequence(process.selectedSlimmedElectrons + process.egmGsfElectronIDSequence + process.bareSoftElectrons + process.softElectrons) # (use this version when running VID)

elif ELECORRTYPE == "RunII" :
   process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('calibratedPatElectrons')
   process.electronMVAValueMapProducer.srcMiniAOD= cms.InputTag('calibratedPatElectrons')

if (LEPTON_SETUP == 2016):
   process.bareSoftElectrons.src = cms.InputTag('selectedSlimmedElectrons')
   process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('selectedSlimmedElectrons')
   process.electronMVAValueMapProducer.srcMiniAOD =  cms.InputTag('selectedSlimmedElectrons')
   process.electrons = cms.Sequence(process.selectedSlimmedElectrons + process.egmGsfElectronIDSequence + process.bareSoftElectrons + process.softElectrons) # (use this version when running VID)

if (LEPTON_SETUP == 2017): 
   process.bareSoftElectrons.src = cms.InputTag('selectedSlimmedElectrons')
   process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('selectedSlimmedElectrons')
   process.electronMVAValueMapProducer.srcMiniAOD =  cms.InputTag('selectedSlimmedElectrons')
   process.electrons = cms.Sequence(process.selectedSlimmedElectrons + process.egmGsfElectronIDSequence + process.bareSoftElectrons + process.softElectrons) # (use this version when running VID)

if (LEPTON_SETUP == 2018):
	process.selectedSlimmedElectrons.src = cms.InputTag("slimmedElectrons")
	process.electrons = cms.Sequence(process.selectedSlimmedElectrons + process.calibratedPatElectrons + process.egmGsfElectronIDSequence + process.bareSoftElectrons + process.softElectrons)


#--- TrackLess Electrons
process.bareSoftPhotons = cms.EDFilter("PATPhotonRefSelector",
   src = cms.InputTag("slimmedPhotons"),
   cut = cms.string("pt>7 && abs(eta)<2.5")
   )

process.softPhotons = cms.EDProducer("Philler",
   src    = cms.InputTag("bareSoftPhotons"),
   srcElectron = cms.InputTag("softElectrons"),
   mvaValuesMap = cms.InputTag("photonMVAValueMapProducer:TLEMVAEstimatorRun2Fall15V1Values"),
   mvaValuesMap2 = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig25nsV2Values"),
   sampleType = cms.int32(SAMPLE_TYPE),
   setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
   cut = cms.string("1"), # dxy, dz not applied
   flags = cms.PSet(
        ID = cms.string("userFloat('isBDT')"),
        isSIP = cms.string(SIP),
        isGood = cms.string(GOODLEPTON),
        pass_lepton_ID = cms.string("userFloat('isBDT')"),
        pass_lepton_SIP = cms.string(SIP),
#        isIsoFSRUncorr  = cms.string("userFloat('combRelIsoPF')<"+ELEISOCUT), #TLE isolation is not corrected for FSR gammas.
        ),
   )


process.electronMatch = cms.EDProducer("MCMatcher",       # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                       src         = cms.InputTag("bareSoftElectrons"), # RECO objects to match
                                       matched     = cms.InputTag("prunedGenParticles"), # mc-truth particle collection
                                       mcPdgId     = cms.vint32(11),               # one or more PDG ID (11 = electron); absolute values (see below)
                                       checkCharge = cms.bool(True),               # True = require RECO and MC objects to have the same charge
                                       mcStatus    = cms.vint32(1),                # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                                       maxDeltaR   = cms.double(0.5),              # Minimum deltaR for the match
                                       maxDPtRel   = cms.double(0.5),              # Minimum deltaPt/Pt for the match
                                       resolveAmbiguities    = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
                                       resolveByMatchQuality = cms.bool(False),    # False = just match input in order; True = pick lowest deltaR pair first
                                       )


### ----------------------------------------------------------------------
### Lepton Cleaning (clean electrons collection from muons)
### ----------------------------------------------------------------------

process.cleanSoftElectrons = cms.EDProducer("PATElectronCleaner",
    # pat electron input source
    src = cms.InputTag("softElectrons"),
    # preselection (any string-based cut for pat::Electron)
    preselection = cms.string(''),
    # overlap checking configurables
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
           src       = cms.InputTag("softMuons"), # Start from loose lepton def
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string("userFloat('isGood')"),
           deltaR              = cms.double(0.05),
           checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), # overlaps don't cause the electron to be discared
        )
    ),
    # finalCut (any string-based cut for pat::Electron)
    finalCut = cms.string(''),
)



### ----------------------------------------------------------------------
### Search for FSR candidates
### ----------------------------------------------------------------------

# Create a photon collection; cfg extracted from "UFHZZAnalysisRun2.FSRPhotons.fsrPhotons_cff"
process.fsrPhotons = cms.EDProducer("PhotonFiller",
    electronSrc = cms.InputTag("cleanSoftElectrons"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    photonSel = cms.string(FSRMODE)  # "skip", "passThrough", "Legacy", "RunII"
)

import PhysicsTools.PatAlgos.producersLayer1.pfParticleProducer_cfi
process.boostedFsrPhotons = PhysicsTools.PatAlgos.producersLayer1.pfParticleProducer_cfi.patPFParticles.clone(
    pfCandidateSource = 'fsrPhotons'
)

process.appendPhotons = cms.EDProducer("LeptonPhotonMatcher",
    muonSrc = cms.InputTag("softMuons"),
    electronSrc = cms.InputTag("cleanSoftElectrons"),
    photonSrc = cms.InputTag("boostedFsrPhotons"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    photonSel = cms.string(FSRMODE),  # "skip", "passThrough", "Legacy", "RunII"
    muon_iso_cut = cms.double(MUISOCUT),
    electron_iso_cut = cms.double(ELEISOCUT),
    debug = cms.untracked.bool(False),
    )

if ADDLOOSEELE:
    process.appendPhotons.looseElectronSrc = cms.InputTag("cleanSoftLooseElectrons")
#    process.appendPhotons.tleSrc = cms.InputTag("softPhotons")
#    process.appendPhotons.TLEMinPt = cms.double(25.)

# All leptons, any F/C.
# CAVEAT: merging creates copies of the objects, so that CandViewShallowCloneCombiner is not able to find
# overlaps between merged collections and the original ones.
process.softLeptons = cms.EDProducer("CandViewMerger",
#    src = cms.VInputTag(cms.InputTag("softMuons"), cms.InputTag("cleanSoftElectrons"))
    src = cms.VInputTag(cms.InputTag("appendPhotons:muons"), cms.InputTag("appendPhotons:electrons"))
)




### ----------------------------------------------------------------------
### ----------------------------------------------------------------------
### BUILD CANDIDATES
### ----------------------------------------------------------------------
### ----------------------------------------------------------------------



### ----------------------------------------------------------------------
### Dileptons: combine/merge leptons into intermediate (bare) collections;
###            Embed additional user variables into final collections
### ----------------------------------------------------------------------
TWOGOODLEPTONS = ("userFloat('d0.isGood') && userFloat('d1.isGood')") # Z made of 2 good leptons (ISO not yet applied)

### NOTE: Isolation cut has been moved to ZZ candidates as we now correct for FSR of all four photons.
### Because if this, isBestZ flags are no longer correct; BESTZ_AMONG is set to "" for safety
# ZISO           = ("( (abs(daughter(0).pdgId)==11 && userFloat('d0.combRelIsoPFFSRCorr')<0.5) || (abs(daughter(0).pdgId)==13 && userFloat('d0.combRelIsoPFFSRCorr')<0.4) ) && ( (abs(daughter(1).pdgId)==11 && userFloat('d1.combRelIsoPFFSRCorr')<0.5) || (abs(daughter(1).pdgId)==13 && userFloat('d1.combRelIsoPFFSRCorr')<0.4) )") #ISO after FSR
# ZLEPTONSEL     = TWOGOODLEPTONS + "&&" + ZISO
# BESTZ_AMONG = ( ZLEPTONSEL ) # "Best Z" chosen among those with 2 leptons with ID, SIP, ISO
# BESTZ_AMONG = ("")

ZLEPTONSEL     = TWOGOODLEPTONS # Note: this is without ISO

Z1PRESEL    = (ZLEPTONSEL + " && mass > 40 && mass < 120") # Note: this is without ISO

BESTZ_AMONG = ( Z1PRESEL + "&& userFloat('d0.passCombRelIsoPFFSRCorr') && userFloat('d1.passCombRelIsoPFFSRCorr')" )

TWOGOODISOLEPTONS = ( TWOGOODLEPTONS + "&& userFloat('d0.passCombRelIsoPFFSRCorr') && userFloat('d1.passCombRelIsoPFFSRCorr')" )

# Cut to filter out unneeded ll combinations as upstream as possible
if KEEPLOOSECOMB:
    KEEPLOOSECOMB_CUT = 'mass > 0 && abs(daughter(0).pdgId())==abs(daughter(1).pdgId())' # Propagate also combinations of loose leptons (for debugging); just require same-flavour
else:
    if FSRMODE == "RunII" : # Just keep combinations of tight leptons (passing ID, SIP and ISO)
        KEEPLOOSECOMB_CUT = "mass > 0 && abs(daughter(0).pdgId())==abs(daughter(1).pdgId()) && daughter(0).masterClone.userFloat('isGood') && daughter(1).masterClone.userFloat('isGood') && daughter(0).masterClone.userFloat('passCombRelIsoPFFSRCorr') &&  daughter(1).masterClone.userFloat('passCombRelIsoPFFSRCorr')"
    else :
        print "KEEPLOOSECOMB == False && FSRMODE =! RunII", FSRMODE, "is no longer supported"
        sys.exit()

### ----------------------------------------------------------------------
### Dileptons (Z->ee, Z->mm)
### ----------------------------------------------------------------------

# l+l- (SFOS, both e and mu)
process.bareZCand = cms.EDProducer("PATCandViewShallowCloneCombiner",
    decay = cms.string('softLeptons@+ softLeptons@-'),
    cut = cms.string(KEEPLOOSECOMB_CUT), # see below
    checkCharge = cms.bool(True)
)


if KEEPLOOSECOMB:
    process.bareZCand.cut = cms.string('mass > 0 && abs(daughter(0).pdgId())==abs(daughter(1).pdgId())') # Propagate also combinations of loose leptons (for debugging)
else:
    if FSRMODE == "RunII" : # Just keep combinations of tight leptons (passing ID, SIP and ISO)
        process.bareZCand.cut = cms.string("mass > 0 && abs(daughter(0).pdgId())==abs(daughter(1).pdgId()) && daughter(0).masterClone.userFloat('isGood') && daughter(1).masterClone.userFloat('isGood') && daughter(0).masterClone.userFloat('passCombRelIsoPFFSRCorr') &&  daughter(1).masterClone.userFloat('passCombRelIsoPFFSRCorr')")
    else : # Just keep combinations of tight leptons (passing ID and SIP; iso cannot be required at this point, with the legacy FSR logic)
        process.bareZCand.cut = cms.string("mass > 0 && abs(daughter(0).pdgId())==abs(daughter(1).pdgId()) && daughter(0).masterClone.userFloat('isGood') && daughter(1).masterClone.userFloat('isGood')")

process.ZCand = cms.EDProducer("ZCandidateFiller",
    src = cms.InputTag("bareZCand"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    bestZAmong = cms.string(BESTZ_AMONG),
    FSRMode = cms.string(FSRMODE), # "skip", "Legacy", "RunII"
    flags = cms.PSet(
        GoodLeptons = cms.string(ZLEPTONSEL),
        GoodIsoLeptons = cms.string(TWOGOODISOLEPTONS),
        Z1Presel = cms.string(Z1PRESEL),
    )
)


# ll, same flavour/any charge, for control regions only
process.bareLLCand = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('softLeptons softLeptons'),
    #cut = cms.string('deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)>0.02'), # protect against ghosts
    cut = cms.string('deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)>0.02 && abs(daughter(0).pdgId())==abs(daughter(1).pdgId())'), # protect against ghosts && same flavour
    checkCharge = cms.bool(False)
)
process.LLCand = cms.EDProducer("ZCandidateFiller",
    src = cms.InputTag("bareLLCand"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    bestZAmong = cms.string(BESTZ_AMONG),
    FSRMode = cms.string(FSRMODE), # "skip", "Legacy", "RunII"
    flags = cms.PSet(
        GoodLeptons = cms.string(ZLEPTONSEL),
        Z1Presel = cms.string(Z1PRESEL),
    )
)


### ----------------------------------------------------------------------
### TriLeptons (for fake rate)
### ----------------------------------------------------------------------
Z_PLUS_LEP_MIJ=("sqrt(pow(daughter(0).daughter({0}).energy+daughter(1).energy, 2) - " +
                "     pow(daughter(0).daughter({0}).px    +daughter(1).px, 2) -" +
                "     pow(daughter(0).daughter({0}).py    +daughter(1).py, 2) -" +
                "     pow(daughter(0).daughter({0}).pz    +daughter(1).pz, 2))")

process.ZlCand = cms.EDProducer("PATCandViewShallowCloneCombiner",
    decay = cms.string('ZCand softLeptons'),
    cut = cms.string("deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).eta, daughter(1).phi)>0.02 &&" + # Ghost suppression
                     "deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).eta, daughter(1).phi)>0.02 &&" +
                     ("(daughter(0).daughter(0).charge == daughter(1).charge || %s > 4) && " % ( Z_PLUS_LEP_MIJ.format(0))) +       # mLL>4 for the OS pair (Giovanni's impl)
                     ("(daughter(0).daughter(1).charge == daughter(1).charge || %s > 4) && " % ( Z_PLUS_LEP_MIJ.format(1))) +
                     "daughter(0).masterClone.userFloat('isBestZ') &&" + # This includes the Z1 isolation requirement
                     "daughter(0).masterClone.userFloat('Z1Presel')"
                     ),
    checkCharge = cms.bool(False)
)


### ----------------------------------------------------------------------
### Quadrileptons: combine/merge leptons into intermediate (bare) collections;
###                Embed additional user variables into final collections
### ----------------------------------------------------------------------

FOURGOODLEPTONS    =  ("userFloat('d0.GoodLeptons') && userFloat('d1.GoodLeptons')" +
                       "&& userFloat('d0.worstEleIso') <" + str(ELEISOCUT) +
                       "&& userFloat('d1.worstEleIso') <" + str(ELEISOCUT) +
                       "&& userFloat('d0.worstMuIso') <" + str(MUISOCUT) +
                       "&& userFloat('d1.worstMuIso') <" + str(MUISOCUT)
                       ) #ZZ made of 4 tight leptons passing SIP and ISO


Z1MASS            = "daughter('Z1').mass>40 && daughter('Z1').mass<120"
Z2MASS            = "daughter('Z2').mass>4  && daughter('Z2').mass<120" # (was > 4 in Synch) to deal with m12 cut at gen level
#MLL3On4_12        = "userFloat('mZa')>12" # mll>12 on 3/4 pairs;
#MLLALLCOMB        = "userFloat('mLL6')>4" # mll>4 on 6/6 AF/AS pairs;
MLLALLCOMB        = "userFloat('mLL4')>4" # mll>4 on 4/4 AF/OS pairs;
SMARTMALLCOMB     = "userFloat('passSmartMLL')" # Require swapped-lepton Z2' to be >12 IF Z1' is SF/OS and closer to 91.1876 than mZ1
PT20_10           = ("userFloat('pt1')>20 && userFloat('pt2')>10") #20/10 on any of the 4 leptons
M4l100            = "mass>100"




if SELSETUP=="Legacy": # Default Configuration (Legacy paper): cut on selected best candidate
    print "SELSETUP=Legacy no longer supported after 8fb591a because we now set combRelIsoPFFSRCorr at the level of ZZ; will need to re-introduce it at  the Z level so that ZZCand.bestZAmong (ie isBestZ flag) works again as before"
    sys.exit()
    HASBESTZ          = "daughter('Z1').masterClone.userFloat('isBestZ')"
    BESTCAND_AMONG = (FOURGOODLEPTONS + "&&" +
                      HASBESTZ        + "&&" +
                      Z1MASS
                      )

    SR           = (FOURGOODLEPTONS + "&&" +
                           HASBESTZ        + "&&" +
                           Z1MASS          + "&&" +
                           Z2MASS          + "&&" +
                           MLLALLCOMB      + "&&" +
                           PT20_10         + "&&" +
                           "mass>70"        + "&&" +
                           "daughter('Z2').mass>12")


elif SELSETUP=="allCutsAtOnce": # Select best candidate among those passing all cuts

    BESTCAND_AMONG = (FOURGOODLEPTONS + "&&" +
                      Z1MASS          + "&&" +
                      Z2MASS          + "&&" +
                      MLLALLCOMB      + "&&" +
                      PT20_10         + "&&" +
                      "mass>70"       + "&&" +
                      "daughter('Z2').mass>12"
                      )

    SR = BESTCAND_AMONG



elif SELSETUP=="allCutsAtOnceButMZ2": # Select best candidate among those passing all cuts except mZ2, leave only mZ2 cut at the end

    BESTCAND_AMONG = (FOURGOODLEPTONS + "&&" +
                      Z1MASS          + "&&" +
                      Z2MASS          + "&&" +
                      MLLALLCOMB      + "&&" +
                      PT20_10         + "&&" +
                      "mass>70"
                      )

    SR           = (BESTCAND_AMONG + "&&" +
                           "daughter('Z2').mass>12")



elif SELSETUP=="allCutsAtOncePlusMZb": # Apply also cut on mZb, -> expected to reduce the background but also signal efficiency

    BESTCAND_AMONG = (FOURGOODLEPTONS + "&&" +
                      Z1MASS          + "&&" +
                      Z2MASS          + "&&" +
                      MLLALLCOMB      + "&&" +
                      PT20_10         + "&&" +
                      "mass>70"       + "&&" +
                      "userFloat('mZb')>12" + "&&" +
                      "daughter('Z2').mass>12"
                      )

    SR = BESTCAND_AMONG

elif SELSETUP=="allCutsAtOncePlusSmart": # Apply smarter mZb cut

    BESTCAND_AMONG = (FOURGOODLEPTONS + "&&" +
                      Z1MASS          + "&&" +
                      Z2MASS          + "&&" +
                      MLLALLCOMB      + "&&" +
                      PT20_10         + "&&" +
                      "mass>70"       + "&&" +
                      SMARTMALLCOMB + "&&" +
                      "daughter('Z2').mass>12"
                      )

    SR = BESTCAND_AMONG


else:
    print "Please choose one of the following string for SELSETUP: 'Legacy', 'allCutsAtOnceButMZ2', 'allCutsAtOnce', 'allCutsAtOncePlusMZb', 'allCutsAtOncePlusSmart'"
    sys.exit()


FULLSEL            = (SR      + "&&" +
                      M4l100)


# Ghost Suppression cut
NOGHOST4l = ("deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).daughter(0).eta, daughter(1).daughter(0).phi)>0.02 &&"+
             "deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).daughter(1).eta, daughter(1).daughter(1).phi)>0.02 &&"+
             "deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).daughter(0).eta, daughter(1).daughter(0).phi)>0.02 &&"+
             "deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).daughter(1).eta, daughter(1).daughter(1).phi)>0.02") # protect against ghosts


# Preselection of 4l candidates
LLLLPRESEL = NOGHOST4l # Just suppress candidates with overlapping leptons

#LLLLPRESEL = (NOGHOST4l + "&&" +
#              FOURGOODLEPTONS) # Only retain candidates with 4 tight leptons (ID, iso, SIP)


# ZZ Candidates
if FSRMODE == "Legacy":
    print "\nERROR: FSRMODE=Legacy is no longer supported. Aborting...\n"
    exit(1)


process.bareZZCand= cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('ZCand ZCand'),
    cut = cms.string(LLLLPRESEL),
    checkCharge = cms.bool(True)
)
process.ZZCand = cms.EDProducer("ZZCandidateFiller",
    src = cms.InputTag("bareZZCand"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP),
    superMelaMass = cms.double(SUPERMELA_MASS),
    isMC = cms.bool(IsMC),
    bestCandAmong = cms.PSet(isBestCand = cms.string(BESTCAND_AMONG)),
    bestCandComparator = cms.string(BESTCANDCOMPARATOR),
    ZRolesByMass = cms.bool(True),
    doKinFit = cms.bool(KINREFIT),
    flags = cms.PSet(
        GoodLeptons =  cms.string(FOURGOODLEPTONS),
        Z2Mass  = cms.string(Z2MASS),
        MAllComb = cms.string(MLLALLCOMB),
        SR = cms.string(SR),
        FullSel70 = cms.string(SR), #Obsolete, use "SR"
        FullSel = cms.string(FULLSEL),
    ),
    recoProbabilities = cms.vstring(),
    #These are actually no longer needed after we dropped the Legacy FSR algorithm
    muon_iso_cut = cms.double(MUISOCUT),
    electron_iso_cut = cms.double(ELEISOCUT),
)



### ----------------------------------------------------------------------
### Z+LL Control regions.
### By constryction, daughter(0) is the Z, daughter(1) the LL pair.
### By setting ZRolesByMass = False, daugther('Z1') is set as daughter(0),
### so we can use the same cuts as for the SR.
### ----------------------------------------------------------------------

Z2MM = "abs(daughter(1).daughter(0).pdgId()*daughter(1).daughter(1).pdgId())==169" #Z2 = mumu
Z2EE = "abs(daughter(1).daughter(0).pdgId()*daughter(1).daughter(1).pdgId())==121" #Z2 = ee
Z2LL_SS = "daughter(1).daughter(0).pdgId()==daughter(1).daughter(1).pdgId()"       #Z2 = same-sign, same-flavour
Z2LL_OS = "(daughter(1).daughter(0).pdgId + daughter(1).daughter(1).pdgId) == 0"   #Z2 = l+l-
Z2MM_OS = "daughter(1).daughter(0).pdgId()*daughter(1).daughter(1).pdgId()==-169" #Z2 = mu+mu-
Z2MM_SS = "daughter(1).daughter(0).pdgId()*daughter(1).daughter(1).pdgId()== 169" #Z2 = mu-mu-/mu+mu+
Z2EE_OS = "daughter(1).daughter(0).pdgId()*daughter(1).daughter(1).pdgId()==-121" #Z2 = e+e-
Z2EE_SS = "daughter(1).daughter(0).pdgId()*daughter(1).daughter(1).pdgId()== 121" #Z2 = e-e-/e+e+
Z2ID    = "userFloat('d1.d0.ID')     && userFloat('d1.d1.ID')"                    #ID on LL leptons
Z2SIP   = "userFloat('d1.d0.SIP')< 4 && userFloat('d1.d1.SIP')< 4"                #SIP on LL leptons
CR_Z2MASS = "daughter(1).mass>4  && daughter(1).mass<120"                        #Mass on LL; cut at 4


# Define cuts for selection of the candidates among which the best one is chosen.
CR_BESTCANDBASE = ("userFloat('d0.Z1Presel') && userFloat('d0.worstEleIso') <" + str(ELEISOCUT) +
                   "&& userFloat('d0.worstMuIso') <" + str(MUISOCUT) ) # To be revised

CR_BESTCANDBASE_AA   = ("userFloat('d0.Z1Presel') && userFloat('d0.worstEleIso') <" + str(ELEISOCUT) +
                        "&& userFloat('d0.worstMuIso') <" + str(MUISOCUT) + "&&" +
                        Z2SIP) # base for AA CR: # Z1 with tight leptons passing SIP and ISO, mass cuts; SIP on Z2


CR_BESTZLLss = ""
if SELSETUP == "Legacy":
    # No longer supported; cf comment for "Legacy" in the SR.
    CR_BESTZLLss = CR_BESTCANDBASE_AA + "&&" + "userFloat('d0.isBestZ') &&"+ Z2LL_SS
elif SELSETUP == "allCutsAtOnceButMZ2":
    CR_BESTZLLss = CR_BESTCANDBASE_AA + "&&" + Z2LL_SS + "&&" +CR_Z2MASS + "&&" + MLLALLCOMB + "&&" + PT20_10 + "&& mass>70"
elif SELSETUP == "allCutsAtOnce":
    CR_BESTZLLss = CR_BESTCANDBASE_AA + "&&" + Z2LL_SS + "&&" +CR_Z2MASS + "&&" + MLLALLCOMB + "&&" + PT20_10 + "&&" + "mass>70" + "&&" + "daughter(1).mass>12"
elif SELSETUP == "allCutsAtOncePlusMZb":
    CR_BESTZLLss = CR_BESTCANDBASE_AA + "&&" + Z2LL_SS + "&&" +CR_Z2MASS + "&&" + MLLALLCOMB + "&&" + PT20_10 + "&&" + "mass>70" + "&&" + "daughter(1).mass>12" + "&&" + "userFloat('mZb')>12"
elif SELSETUP == "allCutsAtOncePlusSmart":
    CR_BESTZLLss = CR_BESTCANDBASE_AA + "&&" + Z2LL_SS + "&&" +CR_Z2MASS + "&&" + MLLALLCOMB + "&&" + PT20_10 + "&&" + "mass>70" + "&&" + "daughter(1).mass>12" + "&&" + SMARTMALLCOMB


# Base for the selection cut applied on the best candidate. This almost fully (except for M4l100) overlaps with the cuts defined above, except for startegies where the best candidate is chosen at the beginning (Legacy, allCutsAtOnceButMZ2).
CR_BASESEL = (CR_Z2MASS + "&&" +              # mass cuts on LL
              MLLALLCOMB + "&&" +             # mass cut on all lepton pairs
              PT20_10    + "&&" +             # pT> 20/10 over all 4 l
              "daughter(1).mass>12 &&" +      # mZ2 >12
              "mass>70" )                     # m4l cut

##### CR based on Z+2 opposite sign leptons that pass the loose selection #####

# check weather the 2 leptons from Z2 pass the tight selection
PASSD0 = "(userFloat('d1.d0.isGood') && userFloat('d1.d0.passCombRelIsoPFFSRCorr'))" # FIXME, passCombRelIsoPFFSRCorr result is hard coded
PASSD1 = "(userFloat('d1.d1.isGood') && userFloat('d1.d1.passCombRelIsoPFFSRCorr'))" # FIXME, passCombRelIsoPFFSRCorr result is hard coded
# ... and fill some useful variable needed for the CR logic
FAILD0 = "!" + PASSD0
FAILD1 = "!" + PASSD1
BOTHFAIL = FAILD0 + "&&" + FAILD1
BOTHPASS = PASSD0 + "&&" + PASSD1
PASSD0_XOR_PASSD1 = "((" + PASSD0 + "&&" + FAILD1 + ") || (" + PASSD1 + "&&" + FAILD0 + "))"
PASSD0_OR_PASSD1  = "(" + PASSD0 + "||" + PASSD1 + ")"


CR_BESTZLLos = (CR_BESTCANDBASE_AA    + "&&" +
                CR_BASESEL            + "&&" +
                Z2LL_OS               + "&&" +
                SMARTMALLCOMB         )

# CR 3P1F
CR_BESTZLLos_3P1F = (CR_BESTZLLos + "&&" + PASSD0_OR_PASSD1)
CR_ZLLosSEL_3P1F  = (CR_BESTZLLos + "&&" + PASSD0_XOR_PASSD1) # Is the CR_BESTZLLos request redundant?


# CR 2P2F
CR_BESTZLLos_2P2F   = (CR_BESTZLLos)
CR_ZLLosSEL_2P2F    = (CR_BESTZLLos + "&&" + BOTHFAIL)  # Is the CR_BESTZLLos request redundant?

################################################################################

#Prefix for Z2 daughters, to be used in the "cut" string
#D1D0 = "daughter(1).daughter(0).masterClone"
#D1D1 = "daughter(1).daughter(1).masterClone"

# Z (OSSF,both e/mu) + LL (any F/C, with no ID/iso); this is the starting point for control regions
process.bareZLLCand= cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('ZCand LLCand'),
    cut = cms.string(NOGHOST4l),
#                     "&& !(" + Z2LL_OS + "&&" + D1D0+".userFloat('isGood')&&"+D1D1+".userFloat('isGood')&&"+D1D0+".userFloat('passCombRelIsoPFFSRCorr')&&"+D1D1+".userFloat('passCombRelIsoPFFSRCorr'))"), # Reject OS tight combinations. Note that BOTHFAIL cannot be used here ("d0.d1.xxx" is available only in and after ZZCandidateFiller.
    checkCharge = cms.bool(False)
)
process.ZLLCand = cms.EDProducer("ZZCandidateFiller",
    src = cms.InputTag("bareZLLCand"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP),
    superMelaMass = cms.double(SUPERMELA_MASS),
    isMC = cms.bool(IsMC),
    bestCandComparator = cms.string(BESTCANDCOMPARATOR),
    bestCandAmong = cms.PSet(
      isBestCand    = cms.string("0"), #do not set SR best cand flag
      isBestCRZLLss = cms.string(CR_BESTZLLss),
      isBestCRZLLos_2P2F = cms.string(CR_BESTZLLos_2P2F),
      isBestCRZLLos_3P1F = cms.string(CR_BESTZLLos_3P1F)

    ),
    ZRolesByMass = cms.bool(False),  # daughter('Z1') = daughter(0)
    doKinFit = cms.bool(KINREFIT),
    flags = cms.PSet(
      SR = cms.string(SR),
      CRZLLss = cms.string(CR_BASESEL),             #combine with proper isBestCRZLLss for AA ss/os CRss
      CRZLLos_2P2F = cms.string(CR_ZLLosSEL_2P2F),
      CRZLLos_3P1F = cms.string(CR_ZLLosSEL_3P1F),
    ),
    recoProbabilities = cms.vstring(),
    #These are actually no longer needed after we dropped the Legacy FSR algorithm
    muon_iso_cut = cms.double(MUISOCUT),
    electron_iso_cut = cms.double(ELEISOCUT),
)



### ----------------------------------------------------------------------
### Jets
### ----------------------------------------------------------------------

from RecoJets.JetProducers.PileupJetIDParams_cfi import full_80x_chs
from RecoJets.JetProducers.PileupJetIDCutParams_cfi import full_80x_chs_wp

process.load("CondCore.CondDB.CondDB_cfi")



process.load("RecoJets.JetProducers.PileupJetID_cfi")
process.pileupJetIdUpdated = process.pileupJetId.clone(
    jets=cms.InputTag("slimmedJets"),
    inputIsCorrected=False,
    applyJec=True,
    vertexes=cms.InputTag("offlineSlimmedPrimaryVertices")
    )

### q/g likelihood
qgDatabaseVersion = 'cmssw8020_v2'
process.QGPoolDBESSource = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(messageLevel = cms.untracked.int32(1)),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
        cms.PSet(
            record = cms.string('QGLikelihoodRcd'),
            tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_AK4PFchs'),
            label  = cms.untracked.string('QGL_AK4PFchs')
        ),
      ),
      connect = cms.string('sqlite_fip:ZZAnalysis/AnalysisStep/data/QGTagging/QGL_'+qgDatabaseVersion+'.db')
)
process.es_prefer_qg = cms.ESPrefer('PoolDBESSource','QGPoolDBESSource')
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets = cms.InputTag( 'slimmedJets' )
process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')


theBTagger=""
theBTaggerThr=0
theBTagSFFile=""
theBTagMCEffFile=""
if (LEPTON_SETUP == 2016): #CSVv2M, from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
   theBTagger="pfCombinedInclusiveSecondaryVertexV2BJetTags"
   theBTaggerThr=0.8484
   theBTagSFFile="ZZAnalysis/AnalysisStep/data/BTagging/CSVv2Moriond17_2017_1_26_BtoH.csv" #Preliminary Moriond17 SF, from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
   theBTagMCEffFile="ZZAnalysis/AnalysisStep/data/BTagging/bTagEfficiencies_80X_ICHEP.root"
elif (LEPTON_SETUP == 2017): #DeepCSV, from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
   #FIXME: Add 2018 correctly!
   theBTagger="pfDeepCSVJetTags:probb"
   theBTaggerThr=0.4941
   theBTagSFFile="ZZAnalysis/AnalysisStep/data/BTagging/DeepCSV_94XSF_V3_B_F.csv"
   theBTagMCEffFile="ZZAnalysis/AnalysisStep/data/BTagging/bTagEfficiencies_94X_Moriond18_v1.root"
elif (LEPTON_SETUP == 2018): #DeepCSV, from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
   theBTagger="pfDeepCSVJetTags:probb"
   theBTaggerThr=0.4184
   theBTagSFFile="ZZAnalysis/AnalysisStep/data/BTagging/DeepCSV_102XSF_V1.csv"
   theBTagMCEffFile="ZZAnalysis/AnalysisStep/data/BTagging/bTagEfficiencies_102X_Moriond19.root" 
else:
   sys.exit("ZZ4lAnalysis.py: Need to define the btagging for the new setup!")

process.dressedJets = cms.EDProducer("JetFiller",
    src = cms.InputTag("slimmedJets"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP),
    cut = cms.string("pt>20 && abs(eta)<4.7 && userFloat('JetID') && userFloat('PUjetID')"),
    isMC = cms.bool(IsMC),
    bTaggerName = cms.string(theBTagger),
    bTaggerThreshold = cms.double(theBTaggerThr),
    applyJEC = cms.bool(APPLYJEC),
    jecType = cms.string("AK4PFchs"),
    applyJER = cms.bool(APPLYJER),
    jerType = cms.string("AK4PFchs"),
    bTagSFFile = cms.string(theBTagSFFile),
    bTagMCEffFile = cms.string(theBTagMCEffFile),
    flags = cms.PSet()
    )


### Load JEC
#Taken from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC#Recommended_for_MC

if (APPLYJEC and SAMPLE_TYPE == 2016):
    if IsMC:
        process.jec = cms.ESSource("PoolDBESSource",
            DBParameters = cms.PSet(
                messageLevel = cms.untracked.int32(1)
                ),
            timetype = cms.string('runnumber'),
            toGet = cms.VPSet(
                cms.PSet(
                    record = cms.string('JetCorrectionsRecord'),
                    tag    = cms.string('JetCorrectorParametersCollection_Summer16_23Sep2016V4_MC_AK4PFchs'), #for 80X/Moriond17
                    label  = cms.untracked.string('AK4PFchs')
                    ),
                ),
             connect = cms.string('sqlite_fip:ZZAnalysis/AnalysisStep/data/JEC/Summer16_23Sep2016V4_MC.db'), #for 80X/Moriond17
            )
    else:
        process.jec = cms.ESSource("PoolDBESSource",
            DBParameters = cms.PSet(
                messageLevel = cms.untracked.int32(1)
                ),
            timetype = cms.string('runnumber'),
            toGet = cms.VPSet(
                cms.PSet(
                    record = cms.string('JetCorrectionsRecord'),
                    tag    = cms.string('JetCorrectorParametersCollection_Summer16_23Sep2016AllV4_DATA_AK4PFchs'), #for 80X/Moriond17
                    label  = cms.untracked.string('AK4PFchs')
                    ),
                ),
            connect = cms.string('sqlite_fip:ZZAnalysis/AnalysisStep/data/JEC/Summer16_23Sep2016AllV4_DATA.db'), #for 80X/ICHEP16
            )

    ## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

    ### reapply JEC
    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors
    process.patJetCorrFactorsReapplyJEC = updatedPatJetCorrFactors.clone(
        src = cms.InputTag("slimmedJets"),
        levels = ['L1FastJet','L2Relative','L3Absolute'],
        payload = 'AK4PFchs' )
    if not IsMC:
        process.patJetCorrFactorsReapplyJEC.levels = ['L1FastJet','L2Relative','L3Absolute','L2L3Residual']

    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets
    process.patJetsReapplyJEC = updatedPatJets.clone(
        jetSource = cms.InputTag("slimmedJets"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
        )

    ### add pileup id and discriminant to patJetsReapplyJEC
    process.patJetsReapplyJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
    process.patJetsReapplyJEC.userData.userInts.src += ['pileupJetIdUpdated:fullId']

    ### Replace inputs in QGTagger and dressedJets
    process.QGTagger.srcJets = cms.InputTag( 'patJetsReapplyJEC')
    process.dressedJets.src = cms.InputTag('patJetsReapplyJEC')

if (APPLYJEC and SAMPLE_TYPE == 2017):
    if IsMC:
        process.jec = cms.ESSource("PoolDBESSource",
            DBParameters = cms.PSet(
                messageLevel = cms.untracked.int32(1)
                ),
            timetype = cms.string('runnumber'),
            toGet = cms.VPSet(
                cms.PSet(
                    record = cms.string('JetCorrectionsRecord'),
#                    tag    = cms.string('JetCorrectorParametersCollection_Fall17_17Nov2017_V6_MC_AK4PFchs'), #for 94X/Moriond18
                     tag    = cms.string('JetCorrectorParametersCollection_Fall17_17Nov2017_V8_MC_AK4PFchs'), #for 94X/MET fix
                    label  = cms.untracked.string('AK4PFchs')
                    ),
                ),
#             connect = cms.string('sqlite_fip:ZZAnalysis/AnalysisStep/data/JEC/Fall17_17Nov2017_V6_MC.db'), #for 94X/Moriond18
              connect = cms.string('sqlite_fip:ZZAnalysis/AnalysisStep/data/JEC/Fall17_17Nov2017_V8_MC.db'), #for 94X/MET fix
            )
    else:
        process.jec = cms.ESSource("PoolDBESSource",
            DBParameters = cms.PSet(
                messageLevel = cms.untracked.int32(1)
                ),
            timetype = cms.string('runnumber'),
            toGet = cms.VPSet(
                cms.PSet(
                    record = cms.string('JetCorrectionsRecord'),
                    tag    = cms.string('JetCorrectorParametersCollection_Fall17_17Nov2017BCDEF_V6_DATA_AK4PFchs'), #for 94X/Moriond18
                    label  = cms.untracked.string('AK4PFchs')
                    ),
                ),
            connect = cms.string('sqlite_fip:ZZAnalysis/AnalysisStep/data/JEC/Fall17_17Nov2017BCDEF_V6_DATA.db'),
            )

    ## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

    ### reapply JEC
    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors
    process.patJetCorrFactorsReapplyJEC = updatedPatJetCorrFactors.clone(
        src = cms.InputTag("slimmedJets"),
        levels = ['L1FastJet','L2Relative','L3Absolute'],
        payload = 'AK4PFchs' )
    if not IsMC:
        process.patJetCorrFactorsReapplyJEC.levels = ['L1FastJet','L2Relative','L3Absolute','L2L3Residual']

    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets
    process.patJetsReapplyJEC = updatedPatJets.clone(
        jetSource = cms.InputTag("slimmedJets"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
        )

    ### add pileup id and discriminant to patJetsReapplyJEC
    process.patJetsReapplyJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
    process.patJetsReapplyJEC.userData.userInts.src += ['pileupJetIdUpdated:fullId']

    ### Replace inputs in QGTagger and dressedJets
    process.QGTagger.srcJets = cms.InputTag( 'patJetsReapplyJEC')
    process.dressedJets.src = cms.InputTag('patJetsReapplyJEC')


if (APPLYJEC and SAMPLE_TYPE == 2018):
    if IsMC:
        process.jec = cms.ESSource("PoolDBESSource",
            DBParameters = cms.PSet(
                messageLevel = cms.untracked.int32(1)
                ),
            timetype = cms.string('runnumber'),
            toGet = cms.VPSet(
                cms.PSet(
                    record = cms.string('JetCorrectionsRecord'),
                     tag    = cms.string('JetCorrectorParametersCollection_Autumn18_V8_MC_AK4PFchs'), #for 10_2_X MC
                    label  = cms.untracked.string('AK4PFchs')
                    ),
                ),
              connect = cms.string('sqlite_fip:ZZAnalysis/AnalysisStep/data/JEC/Autumn18_V8_MC.db'), #for 10_2_X MC
            )
    else: # [FIXME] Preliminary recommendation is to use no JEC for data, change this once it is available
        process.jec = cms.ESSource("PoolDBESSource",
            DBParameters = cms.PSet(
                messageLevel = cms.untracked.int32(1)
                ),
            timetype = cms.string('runnumber'),
            toGet = cms.VPSet(
                cms.PSet(
                    record = cms.string('JetCorrectionsRecord'),
                    tag    = cms.string('JetCorrectorParametersCollection_Autumn18_RunABCD_V8_DATA_AK4PFchs'),
                    label  = cms.untracked.string('AK4PFchs')
                    ),
                ),
            connect = cms.string('sqlite_fip:ZZAnalysis/AnalysisStep/data/JEC/Autumn18_RunABCD_V8_DATA.db'),
            )

    ## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

    ### reapply JEC
    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors
    process.patJetCorrFactorsReapplyJEC = updatedPatJetCorrFactors.clone(
        src = cms.InputTag("slimmedJets"),
        levels = ['L1FastJet','L2Relative','L3Absolute'],
        payload = 'AK4PFchs' )
    if not IsMC:
        process.patJetCorrFactorsReapplyJEC.levels = ['L1FastJet','L2Relative','L3Absolute','L2L3Residual']

    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets
    process.patJetsReapplyJEC = updatedPatJets.clone(
        jetSource = cms.InputTag("slimmedJets"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
        )

    ### add pileup id and discriminant to patJetsReapplyJEC
    process.patJetsReapplyJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
    process.patJetsReapplyJEC.userData.userInts.src += ['pileupJetIdUpdated:fullId']

    ### Replace inputs in QGTagger and dressedJets
    process.QGTagger.srcJets = cms.InputTag( 'patJetsReapplyJEC')
    process.dressedJets.src = cms.InputTag('patJetsReapplyJEC')


# JER from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyResolution#Accessing_factors_from_database
if (APPLYJER and SAMPLE_TYPE == 2017):
   process.load('Configuration.StandardSequences.Services_cff')
   process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
   from CondCore.DBCommon.CondDBSetup_cfi import *
   process.jer = cms.ESSource("PoolDBESSource",
                                 CondDBSetup,
                                 connect = cms.string('sqlite_fip:ZZAnalysis/AnalysisStep/data/JER/Autumn18_V1_MC.db'),#[FIXME] Include correct file
                                 toGet = cms.VPSet(
                                                   cms.PSet(
                                                            record = cms.string('JetResolutionRcd'),
                                                            tag    = cms.string('JR_Autumn18_V1_MC_PtResolution_AK4PFchs'),
                                                            label  = cms.untracked.string('AK4PFchs_pt')
                                                            ),
                                                   cms.PSet(
                                                            record = cms.string('JetResolutionRcd'),
                                                            tag    = cms.string('JR_Autumn18_V1_MC_PhiResolution_AK4PFchs'),
                                                            label  = cms.untracked.string('AK4PFchs_phi')
                                                            ),
                                                   cms.PSet(
                                                            record = cms.string('JetResolutionScaleFactorRcd'),
                                                            tag    = cms.string('JR_Autumn18_V1_MC_SF_AK4PFchs'),
                                                            label  = cms.untracked.string('AK4PFchs')
                                                            )
                                                   )
                                 )
   process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')

if (APPLYJER and SAMPLE_TYPE == 2017):
   process.load('Configuration.StandardSequences.Services_cff')
   process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
   from CondCore.DBCommon.CondDBSetup_cfi import *
   process.jer = cms.ESSource("PoolDBESSource",
                                 CondDBSetup,
                                 connect = cms.string('sqlite_fip:ZZAnalysis/AnalysisStep/data/JER/Autumn18_V1_MC.db'),#[FIXME] Include correct file
                                 toGet = cms.VPSet(
                                                   cms.PSet(
                                                            record = cms.string('JetResolutionRcd'),
                                                            tag    = cms.string('JR_Autumn18_V1_MC_PtResolution_AK4PFchs'),
                                                            label  = cms.untracked.string('AK4PFchs_pt')
                                                            ),
                                                   cms.PSet(
                                                            record = cms.string('JetResolutionRcd'),
                                                            tag    = cms.string('JR_Autumn18_V1_MC_PhiResolution_AK4PFchs'),
                                                            label  = cms.untracked.string('AK4PFchs_phi')
                                                            ),
                                                   cms.PSet(
                                                            record = cms.string('JetResolutionScaleFactorRcd'),
                                                            tag    = cms.string('JR_Autumn18_V1_MC_SF_AK4PFchs'),
                                                            label  = cms.untracked.string('AK4PFchs')
                                                            )
                                                   )
                                 )
   process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')


if (APPLYJER and SAMPLE_TYPE == 2018):
    process.load('Configuration.StandardSequences.Services_cff')
    process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
    from CondCore.DBCommon.CondDBSetup_cfi import *
    process.jer = cms.ESSource("PoolDBESSource",
                               CondDBSetup,
                               connect = cms.string('sqlite_fip:ZZAnalysis/AnalysisStep/data/JER/Autumn18_V1_MC.db'),
                               toGet = cms.VPSet(
                                                 cms.PSet(
                                                          record = cms.string('JetResolutionRcd'),
                                                          tag    = cms.string('JR_Autumn18_V1_MC_PtResolution_AK4PFchs'),
                                                          label  = cms.untracked.string('AK4PFchs_pt')
                                                          ),
                                                 cms.PSet(
                                                          record = cms.string('JetResolutionRcd'),
                                                          tag    = cms.string('JR_Autumn18_V1_MC_PhiResolution_AK4PFchs'),
                                                          label  = cms.untracked.string('AK4PFchs_phi')
                                                          ),
                                                 cms.PSet(
                                                          record = cms.string('JetResolutionScaleFactorRcd'),
                                                          tag    = cms.string('JR_Autumn18_V1_MC_SF_AK4PFchs'),
                                                          label  = cms.untracked.string('AK4PFchs')
                                                          )
                                                 )
                               )
    process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')

### Clean jets wrt. good (preFSR-)isolated leptons
process.cleanJets = cms.EDProducer("JetsWithLeptonsRemover",
                                   Jets      = cms.InputTag("dressedJets"),
                                   Muons     = cms.InputTag("appendPhotons:muons"),
                                   Electrons = cms.InputTag("appendPhotons:electrons"),
                                   Diboson   = cms.InputTag(""),
                                   JetPreselection      = cms.string(""),
                                   MuonPreselection     = cms.string("userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')"),
                                   ElectronPreselection = cms.string("userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')"),
                                   DiBosonPreselection  = cms.string(""),
                                   MatchingType = cms.string("byDeltaR"),
                                   cleanFSRFromLeptons = cms.bool(True),
                                   DebugPlots = cms.untracked.bool(False),
                                   DebugPrintOuts = cms.untracked.bool(False)
                                   )

if FSRMODE=="Legacy" :
    process.cleanJets.MuonPreselection     = "userFloat('isGood') && userFloat('isIsoFSRUncorr')"
    process.cleanJets.ElectronPreselection = "userFloat('isGood') && userFloat('isIsoFSRUncorr')"
    process.cleanJets.cleanFSRFromLeptons = False


### ----------------------------------------------------------------------
### Missing ET
### ----------------------------------------------------------------------

metTag = cms.InputTag("slimmedMETs")

if (RECORRECTMET and SAMPLE_TYPE == 2016):
   
    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    runMetCorAndUncFromMiniAOD(process,
                               isData=(not IsMC),
                               )
    metTag = cms.InputTag("slimmedMETs","","ZZ")

    ### somehow MET recorrection gets this lost again...
    process.patJetsReapplyJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
    process.patJetsReapplyJEC.userData.userInts.src += ['pileupJetIdUpdated:fullId']

### Recorrect MET, cf. https://indico.cern.ch/event/759372/contributions/3149378/attachments/1721436/2779341/metreport.pdf slide 10
###                and https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription#Instructions_for_9_4_X_X_9_for_2
if (RECORRECTMET and SAMPLE_TYPE == 2017):

    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

    runMetCorAndUncFromMiniAOD(process,
                               isData=(not IsMC),
                               fixEE2017 = True,
                               fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
                               postfix = "ModifiedMET"
                               )
    metTag = cms.InputTag("slimmedMETsModifiedMET","","ZZ")

    ### somehow MET recorrection gets this lost again...
    process.patJetsReapplyJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
    process.patJetsReapplyJEC.userData.userInts.src += ['pileupJetIdUpdated:fullId']

if (RECORRECTMET and SAMPLE_TYPE == 2018):
    
    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    runMetCorAndUncFromMiniAOD(process,
                               isData=(not IsMC),
                               )
    metTag = cms.InputTag("slimmedMETs","","ZZ")

    ### somehow MET recorrection gets this lost again...
    process.patJetsReapplyJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
    process.patJetsReapplyJEC.userData.userInts.src += ['pileupJetIdUpdated:fullId']


### ----------------------------------------------------------------------
### Paths
### ----------------------------------------------------------------------

process.preSkimCounter = cms.EDProducer("EventCountProducer")
process.PVfilter =  cms.Path(process.preSkimCounter+process.goodPrimaryVertices)

if APPLYJEC:
    process.Jets = cms.Path(process.pileupJetIdUpdated + process.patJetCorrFactorsReapplyJEC + process.patJetsReapplyJEC + process.QGTagger + process.dressedJets )
else:
    process.Jets = cms.Path( process.QGTagger + process.dressedJets )

if (RECORRECTMET and SAMPLE_TYPE == 2016):
    if IsMC:
        process.MET = cms.Path(process.fullPatMetSequence)
    else:
        process.MET = cms.Path(process.fullPatMetSequence)

if (RECORRECTMET and SAMPLE_TYPE == 2017):
    if IsMC:
        process.MET = cms.Path(process.fullPatMetSequenceModifiedMET)
    else:
        process.MET = cms.Path(process.fullPatMetSequenceModifiedMET)

if (RECORRECTMET and SAMPLE_TYPE == 2018):
    if IsMC:
        process.MET = cms.Path(process.fullPatMetSequence)
    else:
        process.MET = cms.Path(process.fullPatMetSequence)



### ----------------------------------------------------------------------
### Filters
### ----------------------------------------------------------------------
### Create filter for events with one candidate in the SR
process.ZZCandSR = cms.EDFilter("PATCompositeCandidateRefSelector",
    src = cms.InputTag("ZZCand"),
    cut = cms.string(SR) # That is, all candidates with m4l>70
 )

process.ZZCandFilter = cms.EDFilter("CandViewCountFilter",
                                src = cms.InputTag("ZZCandSR"),
                                minNumber = cms.uint32(1)
                            )

# Prepare lepton collections
process.Candidates = cms.Path(
       process.muons             +
       process.electrons         + process.cleanSoftElectrons +
       process.fsrPhotons        + process.boostedFsrPhotons +
       process.appendPhotons     +
       process.softLeptons       +
       process.cleanJets         +
# Build 4-lepton candidates
       process.bareZCand         + process.ZCand     +
       process.bareZZCand        + process.ZZCand
    )

# Optional sequence to build control regions. To get it, add
#process.CRPath = cms.Path(process.CRZl) # only trilepton
#OR
#process.CRPath = cms.Path(process.CR)   # trilep+4lep CRs

process.CRZl = cms.Sequence(
       process.bareZCand         + process.ZCand     +
       process.ZlCand
   )

process.CR = cms.Sequence(
       process.bareZCand         + process.ZCand     +
       process.bareLLCand        + process.LLCand    +
       process.bareZLLCand       + process.ZLLCand   +
       process.ZlCand
   )


### Skim, triggers and MC filters (Only store filter result, no filter is applied)

# 2012 skim, Reimplementation by Giovanni
#process.load("ZZAnalysis.AnalysisStep.Skim2012_cfg")
#process.SkimSequence = cms.Sequence(process.HZZSkim2012)
#process.Skim = cms.Path(process.SkimSequence)

#SkimPaths = cms.vstring('Skim')
SkimPaths = cms.vstring('PVfilter') #Do not apply skim, just require a good PV

# process.HF = cms.Path(process.heavyflavorfilter)

# FIXME total kin filter?


if (ADDLOOSEELE) :
    import os
    execfile(os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/MasterPy/LooseEle.py")
#    execfile(os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/MasterPy/TracklessEle.py")


process.GlobalTag.toGet = cms.VPSet( (cms.PSet(
        connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
        label = cms.untracked.string('electron_eb_ECALonly'),
        record = cms.string('GBRDWrapperRcd'),
        tag = cms.string('GEDelectron_EBCorrection_80X_EGM_v4')
    ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBCorrection_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EBUncertainty_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBUncertainty_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EECorrection_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EECorrection_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EEUncertainty_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EEUncertainty_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBCorrection_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBCorrection_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBUncertainty_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBUncertainty_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EECorrection_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EECorrection_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EEUncertainty_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EEUncertainty_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBCorrection_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBCorrection_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBUncertainty_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBUncertainty_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EECorrection_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EECorrection_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EEUncertainty_80X_EGM_v4')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EEUncertainty_80X_EGM_v4')
        ),
    )
)
