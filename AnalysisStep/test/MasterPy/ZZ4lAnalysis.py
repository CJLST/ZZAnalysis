import FWCore.ParameterSet.Config as cms
from ZZAnalysis.AnalysisStep.defaults import *

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

#Optional name of the sample/dataset being analyzed
declareDefault("SAMPLENAME", "", globals())

#Type of electron scale correction/smearing: "None", "RunII"
declareDefault("ELECORRTYPE", "None", globals())

#Apply electron escale regression
declareDefault("ELEREGRESSION", "None", globals())

#Apply muon scale correction
declareDefault("APPLYMUCORR", False, globals())

#muon scale correction identifier, "None" for pass-through, "MC_76X_13TeV" for 2015
declareDefault("MUCORRTYPE", "None", globals())

#Reapply JEC
declareDefault("APPLYJEC", True, globals())

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
declareDefault("KINREFIT", True, globals())


if SELSETUP=="Legacy" and not BESTCANDCOMPARATOR=="byBestZ1bestZ2":
    print "WARNING: In ZZ4lAnalysis.py the SELSETUP=\"Legacy\" flag is meant to reproduce the Legacy results, ignoring the setting of the BESTCANDCOMPARATOR: ",BESTCANDCOMPARATOR
    BESTCANDCOMPARATOR = "byBestZ1bestZ2"


# The isolation cuts for electrons and muons. FIXME: there is an hardcoded instance of these values in src/LeptonIsoHelper.cc !!
ELEISOCUT = "0.35"
MUISOCUT = "0.35"

### ----------------------------------------------------------------------
### Set the GT
### ----------------------------------------------------------------------
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

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
    from Configuration.AlCa.GlobalTag import GlobalTag
    if IsMC:
        process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
    else:
        process.GlobalTag.globaltag = '76X_dataRun2_v15'

else:
    from Configuration.AlCa.GlobalTag import GlobalTag
    if IsMC:
        process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
    else:
        process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

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

# MuEG

if (LEPTON_SETUP == 2011):
    # DoubleEle
    process.hltFilterDiEle.HLTPaths = ["HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*","HLT_TripleEle10_CaloIdL_TrkIdVL_v*"] # to run on DATA/MC 2011
    if (IsMC):
        # DoubleMu
        process.hltFilterDiMu.HLTPaths = ["HLT_Mu17_Mu8_v*"] # to run on MC 2011    
        process.hltFilterMuEle.HLTPaths = ["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v*","HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*"] # to run on MC 2011

    else :
        process.hltFilterDiMu.HLTPaths = ["HLT_DoubleMu7_v*", "HLT_Mu13_Mu8_v*", "HLT_Mu17_Mu8_v*"] # to run on DATA 2011 NB: Emulation is needed [FIXME]
        process.hltFilterMuEle.HLTPaths = ["HLT_Mu17_Ele8_CaloIdL_v*", "HLT_Mu8_Ele17_CaloIdL_v*"] # For run 1-167913
        process.hltFilterMuEle2.HLTPaths = ["HLT_Mu17_Ele8_CaloIdL_v*", "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*"] # run 167914-175972
        process.hltFilterMuEle3.HLTPaths = ["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v*","HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*"] # : run 175973-180252
        process.triggerMuEle2  = cms.Path(process.hltFilterMuEle2)
        process.triggerMuEle3  = cms.Path(process.hltFilterMuEle3)

elif (LEPTON_SETUP == 2012):
    # DoubleEle
    process.hltFilterDiEle.HLTPaths = ["HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"] # to run on DATA/MC 2012
    process.hltFilterDiMu.HLTPaths = ["HLT_Mu17_Mu8_v*", "HLT_Mu17_TkMu8_v*"] # to run on DATA/MC 2012
    process.hltFilterMuEle.HLTPaths = ["HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*","HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"] # to run on DATA/MC 2012
    process.hltFilterTriEle.HLTPaths = ["HLT_Ele15_Ele8_Ele5_CaloIdL_TrkIdVL_v*"]
    process.triggerTriEle  = cms.Path(process.hltFilterTriEle)

elif (LEPTON_SETUP == 2015):
    #FIXME: For now, the MC paths are the ones used in the RunIISpring15DR74 MC samples for 25ns,7e33 conditions.
    process.hltFilterDiEle.HLTPaths = ["HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"]
    #process.hltFilterDiEle.HLTPaths = ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"] #for 25ns,14e33 (not necessary in 2015 data)
    process.hltFilterDiMu.HLTPaths = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*"]
    process.hltFilterMuEle.HLTPaths = ["HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v*","HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v*"]
    #process.hltFilterMuEle.HLTPaths = ["HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v*","HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v*"] #for 25ns,14e33 (not necessary in 2015 data)
    process.hltFilterTriEle.HLTPaths = ["HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*"]
    process.hltFilterTriMu.HLTPaths = ["HLT_TripleMu_12_10_5_v*"]
    process.hltFilterSingleEle.HLTPaths = ["HLT_Ele23_WPLoose_Gsf_v*"]
    process.triggerTriEle = cms.Path(process.hltFilterTriEle)
    process.triggerTriMu  = cms.Path(process.hltFilterTriMu )
    process.triggerSingleEle = cms.Path(process.hltFilterSingleEle)

elif (LEPTON_SETUP == 2016):
    if (IsMC):
	#At the moment, the HLT paths are not present in the "tranche 1" background MC and MiniAODv1 signal MC. They will be added in "tranche 2"/"tranche 3" and MiniAODv2. 
	process.hltFilterDiEle.HLTPaths = ["*"]
        process.hltFilterDiMu.HLTPaths = ["*"]
        process.hltFilterMuEle.HLTPaths = ["*"]
        process.hltFilterTriEle.HLTPaths = ["*"]
        process.hltFilterTriMu.HLTPaths = ["*"]
        process.hltFilterSingleEle.HLTPaths = ["*"]
        process.hltFilterSingleMu.HLTPaths = ["*"]
    
    else:
        process.hltFilterDiEle.HLTPaths = ["HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*"]
        process.hltFilterDiMu.HLTPaths = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*"]
        process.hltFilterMuEle.HLTPaths = ["HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v*","HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v*"]
        process.hltFilterTriEle.HLTPaths = ["HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*"]
        process.hltFilterTriMu.HLTPaths = ["HLT_TripleMu_12_10_5_v*"]
        process.hltFilterSingleEle.HLTPaths = ["HLT_Ele25_eta2p1_WPTight_Gsf_v*","HLT_Ele27_WPTight_Gsf_v*","HLT_Ele27_eta2p1_WPLoose_Gsf_v*"]
        process.hltFilterSingleMu.HLTPaths = ["HLT_IsoMu20_v*","HLT_IsoTkMu20_v*","HLT_IsoMu22_v*","HLT_IsoTkMu22_v*"]

    process.triggerTriEle = cms.Path(process.hltFilterTriEle)
    process.triggerTriMu  = cms.Path(process.hltFilterTriMu )
    process.triggerSingleEle = cms.Path(process.hltFilterSingleEle)
    process.triggerSingleMu  = cms.Path(process.hltFilterSingleMu )

process.triggerDiMu   = cms.Path(process.hltFilterDiMu)
process.triggerDiEle  = cms.Path(process.hltFilterDiEle)
process.triggerMuEle  = cms.Path(process.hltFilterMuEle)


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
### ----------------------------------------------------------------------
### Loose lepton selection + cleaning + embeddding of user data
### ----------------------------------------------------------------------
### ----------------------------------------------------------------------

SIP =  "userFloat('SIP')<4"
GOODLEPTON = "userFloat('ID') && " + SIP  # Lepton passing tight ID + SIP [ISO is asked AFTER FSR!!!]

if LEPTON_SETUP == 2016:
    TIGHTMUON = "userFloat('isPFMuon') || (userFloat('isTrackerHighPtMuon') && pt>200)"
else:
    TIGHTMUON = "userFloat('isPFMuon')"



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
                                         identifier = cms.string(MUCORRTYPE),
                                         isMC = cms.bool(IsMC),
                                         isSynchronization = cms.bool(False),
                                         )

#--- Set correct identifier for muon corrections
if LEPTON_SETUP == 2011: # (MuScleFit)
    if IsMC:
        process.calibratedMuons.identifier = cms.string("Fall11_START42")
    else:
        process.calibratedMuons.identifier = cms.string("Data2011_42X")
elif LEPTON_SETUP == 2012: # (MuScleFit)
    if IsMC:
        process.calibratedMuons.identifier = cms.string("Summer12_DR53X_smearReReco")
    else:
        process.calibratedMuons.identifier = cms.string("Data2012_53X_ReReco")
else: # (KalmanMuonCalibrator, 2015)
    if IsMC:
        process.calibratedMuons.identifier = cms.string("MC_76X_13TeV")
    else:
        process.calibratedMuons.identifier = cms.string("DATA_76X_13TeV")


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
    flags = cms.PSet(
        ID = cms.string(TIGHTMUON), # tight muon ID
        isSIP = cms.string(SIP),
        isGood = cms.string(GOODLEPTON),
        isIsoFSRUncorr  = cms.string("userFloat('combRelIsoPF')<" + MUISOCUT),
#       Note: passCombRelIsoPFFSRCorr is currently set in LeptonPhotonMatcher for new FSR strategy; in ZZCandidateFiller for the old one
    )
)


if APPLYMUCORR :
    process.muons =  cms.Sequence(process.calibratedMuons + process.cleanedMu + process.bareSoftMuons + process.softMuons)
else:
    process.cleanedMu.src = cms.InputTag("slimmedMuons")
    process.muons =  cms.Sequence(process.cleanedMu + process.bareSoftMuons + process.softMuons)
    



#------- ELECTRONS -------

##--- Electron regression+calibrarion must be applied after BDT is recomputed
## NOTE patElectronsWithRegression->eleRegressionEnergy;  calibratedElectrons-> calibratedPatElectrons
## Default: NEW ECAL regression + NEW calibration + NEW combination
#process.load('EgammaAnalysis.ElectronTools.electronRegressionEnergyProducer_cfi')
#process.eleRegressionEnergy.inputElectronsTag = cms.InputTag('slimmedElectrons')
#process.eleRegressionEnergy.energyRegressionType = 2 ## 1: ECAL regression w/o subclusters 2 (default): ECAL regression w/ subclusters)
##process.eleRegressionEnergy.vertexCollection = cms.InputTag('goodPrimaryVertices')

#process.load("EgammaAnalysis.ElectronTools.calibratedPatElectrons_cfi")
#process.calibratedPatElectrons.correctionsType = 2 # 1 = old regression, 2 = new regression, 3 = no regression, 0 = nothing
#process.calibratedPatElectrons.combinationType = 3
#process.calibratedPatElectrons.lumiRatio = cms.double(1.0)
#process.calibratedPatElectrons.isMC    = IsMC
#process.calibratedPatElectrons.synchronization = cms.bool(False)

#if (LEPTON_SETUP == 2011):
#   process.eleRegressionEnergy.rhoCollection = cms.InputTag('kt6PFJetsForIso:rho')
#   if (IsMC):
#       process.calibratedPatElectrons.inputDataset = "Fall11"
#   else :
#       process.calibratedPatElectrons.inputDataset = "Jan16ReReco"
#elif (LEPTON_SETUP == 2012) :
#   if (IsMC):
#       process.calibratedPatElectrons.inputDataset = "Summer12_LegacyPaper"
#   else :
#       process.calibratedPatElectrons.inputDataset = "22Jan2013ReReco"
#else :
#   if (IsMC):
#       process.calibratedPatElectrons.inputDataset = ""
#   else :
#       process.calibratedPatElectrons.inputDataset = ""


#--- Run2 electron momentum scale and resolution corrections

process.selectedSlimmedElectrons = cms.EDFilter("PATElectronSelector", 
    ## this protects against a crash in electron calibration
    ## due to electrons with eta > 2.5
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string("pt>5 && abs(eta)<2.5")
)

process.calibratedPatElectrons = cms.EDProducer("CalibratedPatElectronProducerRun2",
    electrons = cms.InputTag('selectedSlimmedElectrons'),
    gbrForestName = cms.string("gedelectron_p4combination_25ns"),
    isMC = cms.bool(IsMC),
    isSynchronization = cms.bool(False),
    correctionFile = cms.string("EgammaAnalysis/ElectronTools/data/76X_16DecRereco_2015")
)

if (BUNCH_SPACING == 50):
    process.calibratedPatElectrons.grbForestName = cms.string("gedelectron_p4combination_50ns")


#--- Set up electron ID (VID framework)
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format to be DataFormat.MiniAOD, as appropriate
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_V1_cff']
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
        isIsoFSRUncorr  = cms.string("userFloat('combRelIsoPF')<"+ELEISOCUT)
#       Note: passCombRelIsoPFFSRCorr is currently set in LeptonPhotonMatcher for new FSR strategy; in ZZCandidateFiller for the old one
        ),
   mvaValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16V1Values"), # (when running VID)
   )

process.electrons = cms.Sequence(process.selectedSlimmedElectrons + process.calibratedPatElectrons + process.egmGsfElectronIDSequence + process.bareSoftElectrons + process.softElectrons) # (use this version when running VID)
#process.electrons = cms.Sequence(process.selectedSlimmedElectrons + process.calibratedPatElectrons + process.bareSoftElectrons + process.softElectrons) # (use this version without VID)

# Handle special cases
if ELEREGRESSION == "None" and (ELECORRTYPE == "None" or BUNCH_SPACING == 50) :   # No correction at all. Skip correction modules.
    process.bareSoftElectrons.src = cms.InputTag('slimmedElectrons')
    process.electrons = cms.Sequence(process.egmGsfElectronIDSequence + process.bareSoftElectrons + process.softElectrons)

#elif ELEREGRESSION == "None" and ELECORRTYPE == "RunII" :
#    process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('calibratedPatElectrons') # (when running VID)
#    process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('calibratedPatElectrons') # (when running VID)

#elif ELEREGRESSION == "Moriond" and ELECORRTYPE == "Moriond" : # Moriond corrections: OLD ECAL regression + OLD calibration + OLD combination 
#    process.electrons = cms.Sequence(process.eleRegressionEnergy + process.calibratedPatElectrons + process.bareSoftElectrons + process.softElectrons)
#    if (LEPTON_SETUP == 2011):
#        process.eleRegressionEnergy.regressionInputFile = cms.string("EgammaAnalysis/ElectronTools/data/eleEnergyReg2011Weights_V1.root")
#    else :
#        process.eleRegressionEnergy.regressionInputFile = cms.string("EgammaAnalysis/ElectronTools/data/eleEnergyReg2012Weights_V1.root")
#    process.eleRegressionEnergy.energyRegressionType = 1
#    process.calibratedPatElectrons.correctionsType   = 1
#    process.calibratedPatElectrons.combinationType   = 1
#
#elif ELEREGRESSION == "Paper" and ELECORRTYPE == "None" : # NEW ECAL regression + NO calibration + NO combination
#    process.electrons = cms.Sequence(process.eleRegressionEnergy + process.calibratedPatElectrons + process.bareSoftElectrons + process.softElectrons)
#    process.eleRegressionEnergy.energyRegressionType = 2
#    process.calibratedPatElectrons.correctionsType   = 0
#    process.calibratedPatElectrons.combinationType   = 0
#
#elif ELEREGRESSION == "Paper" and ELECORRTYPE == "PaperNoComb" : # NEW ECAL regression + NEW calibration + NO combination
#    process.electrons = cms.Sequence(process.eleRegressionEnergy + process.calibratedPatElectrons + process.bareSoftElectrons + process.softElectrons)
#    process.eleRegressionEnergy.energyRegressionType = 2
#    process.calibratedPatElectrons.correctionsType   = 1
#    process.calibratedPatElectrons.combinationType   = 0
    

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
#    looseElectronSrc = cms.InputTag("cleanSoftLooseElectrons"),
    photonSrc = cms.InputTag("boostedFsrPhotons"),
#    tleSrc = cms.InputTag("softPhotons"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    photonSel = cms.string(FSRMODE),  # "skip", "passThrough", "Legacy", "RunII"
    muon_iso_cut = cms.double(MUISOCUT),
    electron_iso_cut = cms.double(ELEISOCUT),
    )



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

### ----------------------------------------------------------------------
### Dileptons (Z->ee, Z->mm)
### ----------------------------------------------------------------------

# l+l- (SFOS, both e and mu)
process.bareZCand = cms.EDProducer("PATCandViewShallowCloneCombiner",
    decay = cms.string('softLeptons@+ softLeptons@-'),
    cut = cms.string('0'), # see below
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


# ll, any combination of flavour/charge, for control regions only
process.bareLLCand = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('softLeptons softLeptons'),
    cut = cms.string('deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)>0.02'), # protect against ghosts
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
                       "&& userFloat('d0.worstEleIso') <" + ELEISOCUT +
                       "&& userFloat('d1.worstEleIso') <" + ELEISOCUT +
                       "&& userFloat('d0.worstMuIso') <" + MUISOCUT +
                       "&& userFloat('d1.worstMuIso') <" + MUISOCUT
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
    #These are actually no longer needed after we dropped the Legacy FSR algorithm
    muon_iso_cut = cms.double(MUISOCUT),
    electron_iso_cut = cms.double(ELEISOCUT),
)

# LHECandidate
process.LHECand = cms.EDProducer(
    "LHECandidateFiller",
    isMC = cms.bool(IsMC)
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
CR_BESTCANDBASE = ("userFloat('d0.Z1Presel') && userFloat('d0.worstEleIso') <" + ELEISOCUT +
                   "&& userFloat('d0.worstMuIso') <" + MUISOCUT ) # To be revised

CR_BESTCANDBASE_AA   = ("userFloat('d0.Z1Presel') && userFloat('d0.worstEleIso') <" + ELEISOCUT +
                        "&& userFloat('d0.worstMuIso') <" + MUISOCUT + "&&" +
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
    #These are actually no longer needed after we dropped the Legacy FSR algorithm
    muon_iso_cut = cms.double(MUISOCUT),
    electron_iso_cut = cms.double(ELEISOCUT),
)



### ----------------------------------------------------------------------
### Jets
### ----------------------------------------------------------------------

# q/g likelihood
process.load("CondCore.CondDB.CondDB_cfi")

qgDatabaseVersion = 'v2b'
QGPoolDBESSource = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(messageLevel = cms.untracked.int32(1)),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(),
      connect = cms.string('frontier://FrontierProd/CMS_COND_PAT_000'),
)
for type in ['AK4PFchs']:
  QGPoolDBESSource.toGet.extend(cms.VPSet(cms.PSet(
    record = cms.string('QGLikelihoodRcd'),
    tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_'+type),
    label  = cms.untracked.string('QGL_'+type)
  )))
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets = cms.InputTag( 'slimmedJets' )
process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')

process.dressedJets = cms.EDProducer("JetFiller",
    src = cms.InputTag("slimmedJets"),
    cut = cms.string("pt>20 && abs(eta)<4.7"),
    bTaggerName = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
    jecType = cms.string("AK4PFchs"),
    flags = cms.PSet(
        isBtagged = cms.string("userFloat('bTagger')>0.800"),
        )
    )

### Load JEC
if APPLYJEC: 
    if IsMC: 
        process.jec = cms.ESSource("PoolDBESSource",
            DBParameters = cms.PSet(
                messageLevel = cms.untracked.int32(1)
                ),
            timetype = cms.string('runnumber'),
            toGet = cms.VPSet(
                cms.PSet(
                    record = cms.string('JetCorrectionsRecord'),
#                    tag    = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_MC_AK4PFchs'), #for 76X
                    tag    = cms.string('JetCorrectorParametersCollection_Spring16_25nsV1_MC_AK4PFchs'), #for 80X
                    label  = cms.untracked.string('AK4PFchs')
                    ),
                ), 
#            connect = cms.string('sqlite_fip:ZZAnalysis/AnalysisStep/data/JEC/Fall15_25nsV2_MC.db'), #for 76X
             connect = cms.string('sqlite_fip:ZZAnalysis/AnalysisStep/data/JEC/Spring16_25nsV1_MC.db'), #for 80X
            
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
                    tag    = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_DATA_AK4PFchs'), #for 76X FIXME: update for 80X
                    label  = cms.untracked.string('AK4PFchs')
                    ),
                ), 
            connect = cms.string('sqlite_fip:ZZAnalysis/AnalysisStep/data/JEC/Fall15_25nsV2_DATA.db'), #for 76X FIXME: update for 80X   
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

    ### Replace inputs in QGTagger and dressedJets
    process.QGTagger.srcJets = cms.InputTag( 'patJetsReapplyJEC')
    process.dressedJets.src = cms.InputTag('patJetsReapplyJEC')



# Clean jets wrt. good (preFSR-)isolated leptons
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
                                   DebugPlots = cms.untracked.bool(False)
                                   )

if FSRMODE=="Legacy" :
    process.cleanJets.MuonPreselection     = "userFloat('isGood') && userFloat('isIsoFSRUncorr')"
    process.cleanJets.ElectronPreselection = "userFloat('isGood') && userFloat('isIsoFSRUncorr')"
    process.cleanJets.cleanFSRFromLeptons = False


### ----------------------------------------------------------------------
### Paths
### ----------------------------------------------------------------------

process.preSkimCounter = cms.EDProducer("EventCountProducer")
process.PVfilter =  cms.Path(process.preSkimCounter+process.goodPrimaryVertices)

if APPLYJEC:
    process.Jets = cms.Path( process.patJetCorrFactorsReapplyJEC + process.patJetsReapplyJEC + process.QGTagger + process.dressedJets )
else:
    process.Jets = cms.Path( process.QGTagger + process.dressedJets )

# LHE info
# Prepare lepton collections
process.LHECandidates = cms.Path(process.LHECand)

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

### 2011 HZZ Skim
#process.afterSkimCounter = cms.EDProducer("EventCountProducer")
#process.load("ZZAnalysis.AnalysisStep.HZZSkim_cfg")
#process.skim = cms.Path(process.skim2011 + process.afterSkimCounter) # the 2011 skim

### 2012 skim.
#FIXME this is the version from /afs/cern.ch/user/p/psilva/public/HZZSkim/PDWG_HZZSkim_cff.py
#      which is buggy!!
#process.load("ZZAnalysis.AnalysisStep.PDWG_HZZSkim_cff") 
#SkimPaths = cms.vstring('HZZ4ePath', 'HZZ2e2mPath', 'HZZ2m2ePath', 'HZZ4mPath', 'HZZem2ePath', 'HZZem2mPath')

# Reimplementation by Giovanni
#process.load("ZZAnalysis.AnalysisStep.Skim2012_cfg")
#process.SkimSequence = cms.Sequence(process.HZZSkim2012)
#process.Skim = cms.Path(process.SkimSequence)


#SkimPaths = cms.vstring('Skim')
SkimPaths = cms.vstring('PVfilter') #Do not apply skim 

# process.HF = cms.Path(process.heavyflavorfilter)

# FIXME total kin filter?



