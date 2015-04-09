import FWCore.ParameterSet.Config as cms
process = cms.Process("ZZ")

### ----------------------------------------------------------------------
### Flags that need to be setted
### ----------------------------------------------------------------------

try:
    IsMC
except NameError:
    IsMC = True

#Set of effective areas, rho corrections, etc. (can be 2011, 2012 or 2015)
try:
    LEPTON_SETUP
except NameError:
    LEPTON_SETUP = 2015

#Flag that reflects the actual sqrts of the sample (can be 2011, 2012 or 2015)
try:
    SAMPLE_TYPE
except NameError:
    SAMPLE_TYPE = LEPTON_SETUP # LEPTON_SETUP can differ from SAMPLE_TYPE for samples that are rescaled to a different sqrts.

#Optional name of the sample/dataset being analyzed
try:
    SAMPLENAME
except NameError:
    SAMPLENAME = ""

#Type of electron scale correction/smearing
try:
    ELECORRTYPE
except NameError:
    ELECORRTYPE = "Paper"

#Apply electron escale regression
try:
    ELEREGRESSION
except NameError:
    ELEREGRESSION = "Paper"
    
#Apply muon scale correction
try:
    APPLYMUCORR
except NameError:
    APPLYMUCORR = True

#Mass used for SuperMELA
try:
    SUPERMELA_MASS
except NameError:
    SUPERMELA_MASS = 125.6

#Selection flow strategy
try:
    SELSETUP
except NameError:
    SELSETUP = "allCutsAtOncePlusSmart"

#Best candidate comparator (see interface/Comparators.h)
try:
    BESTCANDCOMPARATOR
except NameError:
    BESTCANDCOMPARATOR = "byBestZ1bestZ2"

if SELSETUP=="Legacy" and not BESTCANDCOMPARATOR=="byBestZ1bestZ2":
    print "WARNING: In ZZ4lAnalysis.py the SELSETUP=\"Legacy\" flag is meant to reproduce the Legacy results, ignoring the setting of the BESTCANDCOMPARATOR: ",BESTCANDCOMPARATOR
    BESTCANDCOMPARATOR = "byBestZ1bestZ2"


# Set to True to make candidates with the full combinatorial of loose leptons (for debug; much slower)
KEEPLOOSELEP = False

# The isolation cuts for electrons and muons
ELEISOCUT = "0.5"
MUISOCUT = "0.4"



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
        
else: 
    if IsMC:
        process.GlobalTag.globaltag = 'PHYS14_25_V1' #FIXME: for Phys14-25ns only
    else:
        process.GlobalTag.globaltag = '' #FIXME: not yet available for RunII

print process.GlobalTag.globaltag

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
process.hltFilterDiMu.TriggerResultsTag  = cms.InputTag("TriggerResults","","HLT")
process.hltFilterDiEle.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterMuEle.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterMuEle2.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterMuEle3.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterTriEle.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterDiMu.throw  = cms.bool(False) #FIXME: beware of this!
process.hltFilterDiEle.throw = cms.bool(False) #FIXME: beware of this!
process.hltFilterMuEle.throw = cms.bool(False) #FIXME: beware of this!
process.hltFilterMuEle2.throw = cms.bool(False) #FIXME: beware of this!
process.hltFilterMuEle3.throw = cms.bool(False) #FIXME: beware of this!
process.hltFilterTriEle.throw = cms.bool(False) #FIXME: beware of this!

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
    #FIXME: This is the menu used in the PHYS14 samples, i.e. these paths are temporary and will NOT be used in Run II.
    process.hltFilterDiEle.HLTPaths = ["HLT_Ele23_Ele12_CaloId_TrackId_Iso_v*"]
    process.hltFilterDiMu.HLTPaths = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*"]
    process.hltFilterMuEle.HLTPaths = ["HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v*","HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v*"]
    process.hltFilterTriEle.HLTPaths = ["HLT_Ele17_Ele12_Ele10_CaloId_TrackId_v*"]
    process.triggerTriEle  = cms.Path(process.hltFilterTriEle)

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

# # Mu e-scale corrections (Rochester)
# process.calibratedMuons =  cms.EDProducer("RochesterPATMuonCorrector",
#                                                src = cms.InputTag("patMuonsWithTrigger")
#                                               )


# Mu e-scale corrections (MuScleFit)
process.calibratedMuons = cms.EDProducer("MuScleFitPATMuonCorrector", 
                         src = cms.InputTag("slimmedMuons"), 
                         debug = cms.bool(False), 
                         identifier = cms.string("Summer12_DR53X_smearReReco"), 
                         applySmearing = cms.bool(IsMC), 
                         fakeSmearing = cms.bool(False)
                         )

# Set correct identifier for MuScleFit muon corrections
if LEPTON_SETUP == 2011:
    if IsMC:
        process.calibratedMuons.identifier = cms.string("Fall11_START42")
    else:
        process.calibratedMuons.identifier = cms.string("Data2011_42X")
elif LEPTON_SETUP == 2012:
    if IsMC:
        process.calibratedMuons.identifier = cms.string("Summer12_DR53X_smearReReco")
    else:
        process.calibratedMuons.identifier = cms.string("Data2012_53X_ReReco")
else:
    if IsMC:
        process.calibratedMuons.identifier = cms.string("") #FIXME: not yet available for RunII
    else:
        process.calibratedMuons.identifier = cms.string("") #FIXME: not yet available for RunII
        

### Mu Ghost cleaning
process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
                                   src = cms.InputTag("calibratedMuons"),
                                   preselection = cms.string("track.isNonnull"),
                                   passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
                                   fractionOfSharedSegments = cms.double(0.499))

if APPLYMUCORR == False :
    process.cleanedMu.src = cms.InputTag("slimmedMuons")



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
        ID = cms.string("userFloat('isPFMuon')" ), # PF ID
        isSIP = cms.string(SIP),
        isGood = cms.string(GOODLEPTON),
        isIsoFSRUncorr  = cms.string("userFloat('combRelIsoPF')<" + MUISOCUT),
#        isGLB = cms.string("isGlobalMuon"),
#        isTk = cms.string("isTrackerMuon"),
    )
)


if APPLYMUCORR :
    process.muons =  cms.Sequence(process.calibratedMuons + process.cleanedMu + process.bareSoftMuons + process.softMuons)
else:
    process.cleanedMu.src = src = cms.InputTag("slimmedMuons")
    process.muons =  cms.Sequence(process.cleanedMu + process.bareSoftMuons + process.softMuons)
    



#------- ELECTRONS -------
#--- Electron regression+calibrarion must be applied after BDT is recomputed
# NOTE patElectronsWithRegression->eleRegressionEnergy;  calibratedElectrons-> calibratedPatElectrons
# Default: NEW ECAL regression + NEW calibration + NEW combination
process.load('EgammaAnalysis.ElectronTools.electronRegressionEnergyProducer_cfi')
process.eleRegressionEnergy.inputElectronsTag = cms.InputTag('slimmedElectrons')
process.eleRegressionEnergy.energyRegressionType = 2 ## 1: ECAL regression w/o subclusters 2 (default): ECAL regression w/ subclusters)
#process.eleRegressionEnergy.vertexCollection = cms.InputTag('goodPrimaryVertices')
process.load("EgammaAnalysis.ElectronTools.calibratedPatElectrons_cfi")
process.calibratedPatElectrons.correctionsType = 2 # 1 = old regression, 2 = new regression, 3 = no regression, 0 = nothing
process.calibratedPatElectrons.combinationType = 3
process.calibratedPatElectrons.lumiRatio = cms.double(1.0)
process.calibratedPatElectrons.isMC    = IsMC
process.calibratedPatElectrons.synchronization = cms.bool(False)

if (LEPTON_SETUP == 2011):
   process.eleRegressionEnergy.rhoCollection = cms.InputTag('kt6PFJetsForIso:rho')
   if (IsMC):
       process.calibratedPatElectrons.inputDataset = "Fall11"
   else :
       process.calibratedPatElectrons.inputDataset = "Jan16ReReco"
elif (LEPTON_SETUP == 2012) :
   if (IsMC):
       process.calibratedPatElectrons.inputDataset = "Summer12_LegacyPaper"
   else :
       process.calibratedPatElectrons.inputDataset = "22Jan2013ReReco"
else :
   if (IsMC):
       process.calibratedPatElectrons.inputDataset = "" #FIXME: not yet available for RunII
   else :
       process.calibratedPatElectrons.inputDataset = "" #FIXME: not yet available for RunII

    

process.bareSoftElectrons = cms.EDFilter("PATElectronRefSelector",
   src = cms.InputTag("calibratedPatElectrons"),
   cut = cms.string("pt>7 && abs(eta)<2.5")
   )

process.softElectrons = cms.EDProducer("EleFiller",
   src    = cms.InputTag("bareSoftElectrons"),
   sampleType = cms.int32(SAMPLE_TYPE),          
   setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
   cut = cms.string("userFloat('dxy')<0.5 && userFloat('dz')<1 && userFloat('missingHit')<=1"),
   flags = cms.PSet(
        ID = cms.string("userFloat('isBDT')"), # BDT MVA ID
        isSIP = cms.string(SIP),
        isGood = cms.string(GOODLEPTON),
        isIsoFSRUncorr  = cms.string("userFloat('combRelIsoPF')<"+ELEISOCUT)
        )
   )


process.electrons = cms.Sequence(process.eleRegressionEnergy + process.calibratedPatElectrons + process.bareSoftElectrons + process.softElectrons)

# Handle special cases
if ELEREGRESSION == "None" and ELECORRTYPE == "None" :   # No correction at all. Skip correction modules.
    process.bareSoftElectrons.src = cms.InputTag('slimmedElectrons')
    process.electrons = cms.Sequence(process.bareSoftElectrons + process.softElectrons)

elif ELEREGRESSION == "Moriond" and ELECORRTYPE == "Moriond" : # Moriond corrections: OLD ECAL regression + OLD calibration + OLD combination 
    if (LEPTON_SETUP == 2011):
        process.eleRegressionEnergy.regressionInputFile = cms.string("EgammaAnalysis/ElectronTools/data/eleEnergyReg2011Weights_V1.root")
    else :
        process.eleRegressionEnergy.regressionInputFile = cms.string("EgammaAnalysis/ElectronTools/data/eleEnergyReg2012Weights_V1.root")
    process.eleRegressionEnergy.energyRegressionType = 1
    process.calibratedPatElectrons.correctionsType   = 1
    process.calibratedPatElectrons.combinationType   = 1

elif ELEREGRESSION == "Paper" and ELECORRTYPE == "None" : # NEW ECAL regression + NO calibration + NO combination
    process.eleRegressionEnergy.energyRegressionType = 2
    process.calibratedPatElectrons.correctionsType   = 0
    process.calibratedPatElectrons.combinationType   = 0

elif ELEREGRESSION == "Paper" and ELECORRTYPE == "PaperNoComb" : # NEW ECAL regression + NEW calibration + NO combination
    process.eleRegressionEnergy.energyRegressionType = 2
    process.calibratedPatElectrons.correctionsType   = 1
    process.calibratedPatElectrons.combinationType   = 0
    

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
           preselection        = cms.string("(isGlobalMuon || userFloat('isPFMuon'))&&userFloat('SIP')<4"), #
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
process.fsrPhotons = cms.EDProducer("FSRPhotonProducer",
    srcCands = cms.InputTag("packedPFCandidates"),
    muons = cms.InputTag("slimmedMuons"), 
    ptThresh = cms.double(2.0),
    extractMuonFSR = cms.bool(False),
)
import PhysicsTools.PatAlgos.producersLayer1.pfParticleProducer_cfi 
process.boostedFsrPhotons = PhysicsTools.PatAlgos.producersLayer1.pfParticleProducer_cfi.patPFParticles.clone(
    pfCandidateSource = 'fsrPhotons'
)

process.appendPhotons = cms.EDProducer("LeptonPhotonMatcher",
    muonSrc = cms.InputTag("softMuons"),
    electronSrc = cms.InputTag("cleanSoftElectrons"),
    photonSrc = cms.InputTag("boostedFsrPhotons"),
    matchFSR = cms.bool(True)
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
BESTZ_AMONG = ("")

ZLEPTONSEL     = TWOGOODLEPTONS # Note: this is without ISO

Z1PRESEL    = (ZLEPTONSEL + " && mass > 40 && mass < 120") # Note: this is without ISO


### ----------------------------------------------------------------------
### Dileptons (Z->ee, Z->mm)
### ----------------------------------------------------------------------

# l+l- (SFOS, both e and mu)
process.bareZCand = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('softLeptons@+ softLeptons@-'),
    cut = cms.string('mass > 0 && abs(daughter(0).pdgId())==abs(daughter(1).pdgId())'),
    checkCharge = cms.bool(True)
)

#if KEEPLOOSELEP:
#    process.bareZCand.cut = cms.string('mass > 0 && abs(daughter(0).pdgId())==abs(daughter(1).pdgId())')
#else:
#    process.bareZCand.cut = cms.string('mass > 0 && abs(daughter(0).pdgId())==abs(daughter(1).pdgId()) &&' + TWOGOODLEPTONS)
    
process.ZCand = cms.EDProducer("ZCandidateFiller",
    src = cms.InputTag("bareZCand"),
    sampleType = cms.int32(SAMPLE_TYPE),                     
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    bestZAmong = cms.string(BESTZ_AMONG),
    flags = cms.PSet(
        GoodLeptons = cms.string(ZLEPTONSEL),
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
    flags = cms.PSet(
        GoodLeptons = cms.string(ZLEPTONSEL),
        Z1Presel = cms.string(Z1PRESEL),
    )
)


### ----------------------------------------------------------------------
### TriLeptons (for fake rate)
### ----------------------------------------------------------------------
#FIXME: to be reviseed since ISOLATION is no longer included in isBestZ and Z1Presel!
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
                     "daughter(0).masterClone.userFloat('isBestZ') &&" +
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
Z2MASS            = "daughter('Z2').mass>4  && daughter('Z2').mass<120" # (was > 4 in Synch) to deal with m12 cut at gen level #FIXME
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
    
    FULLSEL70           = (FOURGOODLEPTONS + "&&" +
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

    FULLSEL70 = BESTCAND_AMONG



elif SELSETUP=="allCutsAtOnceButMZ2": # Select best candidate among those passing all cuts except mZ2, leave only mZ2 cut at the end

    BESTCAND_AMONG = (FOURGOODLEPTONS + "&&" +
                      Z1MASS          + "&&" +
                      Z2MASS          + "&&" +
                      MLLALLCOMB      + "&&" +
                      PT20_10         + "&&" +
                      "mass>70"                        
                      )

    FULLSEL70           = (BESTCAND_AMONG + "&&" +
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

    FULLSEL70 = BESTCAND_AMONG

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

    FULLSEL70 = BESTCAND_AMONG


else:
    print "Please choose one of the following string for SELSETUP: 'Legacy', 'allCutsAtOnceButMZ2', 'allCutsAtOnce', 'allCutsAtOncePlusMZb', 'allCutsAtOncePlusSmart'"
    sys.exit()


FULLSEL            = (FULLSEL70      + "&&" +
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
    flags = cms.PSet(
        GoodLeptons =  cms.string(FOURGOODLEPTONS),
        Z2Mass  = cms.string(Z2MASS),
        MAllComb = cms.string(MLLALLCOMB),
        FullSel70 = cms.string(FULLSEL70),
        FullSel = cms.string(FULLSEL),
    )
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
              M4l100 )                        # m4l>100 

##### CR based on Z+2 opposite sign leptons that pass the loose selection #####

# check weather the 2 leptons from Z2 pass the tight selection
PASSD0 = "(userFloat('d1.d0.isGood') && userFloat('d1.d0.passCombRelIsoPFFSRCorr'))" # FIXME, passCombRelIsoPFFSRCorr result is hard coded
PASSD1 = "(userFloat('d1.d1.isGood') && userFloat('d1.d1.passCombRelIsoPFFSRCorr'))" # FIXME, passCombRelIsoPFFSRCorr result is hard coded
# ... and fill some useful variable needed for the CR logic
FAILD0 = "!" + PASSD0
FAILD1 = "!" + PASSD1
BOTHFAIL = FAILD0 + "&&" + FAILD1
PASSD0_XOR_PASSD1 = "((" + PASSD0 + "&&" + FAILD1 + ") || (" + PASSD1 + "&&" + FAILD0 + "))"
PASSD0_OR_PASSD1  = "(" + PASSD0 + "||" + PASSD1 + ")"

CR_BESTZLLos = (CR_BESTCANDBASE_AA    + "&&" +  
                Z2LL_OS               + "&&" +  
                CR_Z2MASS             + "&&" +
                MLLALLCOMB            + "&&" +
                PT20_10               + "&&" + 
                "mass > 70 &&"               +
                "daughter(1).mass>12" + "&&" +
                SMARTMALLCOMB         )

# CR 3P1F
CR_BESTZLLos_3P1F = (CR_BESTZLLos + "&&" + PASSD0_OR_PASSD1)                 
CR_ZLLosSEL_3P1F  = (CR_BASESEL    + "&&" + PASSD0_XOR_PASSD1)


# CR 2P2F
CR_BESTZLLos_2P2F   = (CR_BESTZLLos)
CR_ZLLosSEL_2P2F    = (CR_BASESEL + "&&" + BOTHFAIL)

################################################################################




# Z (OSSF,both e/mu) + LL (any F/C, with no ID/iso); this is the starting point for control regions
process.bareZLLCand= cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('ZCand LLCand'),
    cut = cms.string(NOGHOST4l), #Note that LLLLPRESEL cannot be used here
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
      isBestCRZLL = cms.string(CR_BESTCANDBASE+ "&&" + # FIXME: To be revised!!!
                         "userFloat('d0.isBestZ') &&"+
                         Z2ID
                         ),
      isBestCRZLLss = cms.string(CR_BESTZLLss),
      isBestCRZLLos_2P2F = cms.string(CR_BESTZLLos_2P2F),
      isBestCRZLLos_3P1F = cms.string(CR_BESTZLLos_3P1F)

    ),
    ZRolesByMass = cms.bool(False),  # daughter('Z1') = daughter(0)
    flags = cms.PSet(
      SR = cms.string(BESTCAND_AMONG),
      CRZLL =  cms.string(CR_BASESEL),             # with isBestCRZLL flag = no SIP CR
      CRZLLss = cms.string(CR_BASESEL),             #combine with proper isBestCRZLLss for AA ss/os CRss    
      CRZLLos_2P2F = cms.string(CR_ZLLosSEL_2P2F),        
      CRZLLos_3P1F = cms.string(CR_ZLLosSEL_3P1F),        
    )
)



### ----------------------------------------------------------------------
### Recorrect jets
### ----------------------------------------------------------------------

UPDATE_JETS = False #FIXME: deactivated for now; what should be done for RunII ?

if (UPDATE_JETS and LEPTON_SETUP==2012) :
#    from CMGTools.Common.miscProducers.cmgPFJetCorrector_cfi import cmgPFJetCorrector
    process.cmgPFJetSel = cms.EDProducer( "PFJetCorrector",
                                        # make sure your jet and rho collections are compatible
                                        src = cms.InputTag( 'slimmedJets' ),
                                        rho = cms.InputTag( 'kt6PFJets:rho:RECO' ),
                                        vertices = cms.InputTag('offlinePrimaryVertices'),
                                        payload = cms.string('AK5PF'),
                                        levels = cms.vstring('L1FastJet','L2Relative','L3Absolute'),
                                        sort = cms.bool(True),
                                        verbose = cms.untracked.bool( False ))    
    if not IsMC:
        process.cmgPFJetSel.levels = ['L1FastJet','L2Relative','L3Absolute','L2L3Residual']

#        print 'Correction levels', process.cmgPFJetSel.levels


### ----------------------------------------------------------------------
### Paths
### ----------------------------------------------------------------------

process.preSkimCounter = cms.EDProducer("EventCountProducer")
process.PVfilter =  cms.Path(process.preSkimCounter+process.goodPrimaryVertices)


# Prepare lepton collections
process.Candidates = cms.Path(
       process.muons             +
       process.electrons         + process.cleanSoftElectrons +
       process.fsrPhotons        + process.boostedFsrPhotons +
       process.appendPhotons     +
       process.softLeptons       +
# Build 4-lepton candidates
       process.bareZCand         + process.ZCand     +  
       process.bareZZCand        + process.ZZCand
    )

if (UPDATE_JETS and LEPTON_SETUP==2012) :
    process.Candidates.insert(0, process.cmgPFJetSel)

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



