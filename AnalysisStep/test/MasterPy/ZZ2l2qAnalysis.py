import FWCore.ParameterSet.Config as cms
process = cms.Process("ZZ")

### ----------------------------------------------------------------------
### Flags that need to be set
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
    ELECORRTYPE = "RunII"

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

#Bunch spacing (can be 25 or 50)
try:
    BUNCH_SPACING
except NameError:
    BUNCH_SPACING = 25

#Mass used for SuperMELA
try:
    SUPERMELA_MASS
except NameError:
    SUPERMELA_MASS = 125

#Selection flow strategy
try:
    SELSETUP
except NameError:
    SELSETUP = "allCutsAtOncePlusSmart"

#Best candidate comparator (see interface/Comparators.h)
try:
    BESTCANDCOMPARATOR
except NameError:
    BESTCANDCOMPARATOR = "byBestZqq"

if SELSETUP=="Legacy" and not BESTCANDCOMPARATOR=="byBestZ1bestZ2":
    print "WARNING: In ZZ4lAnalysis.py the SELSETUP=\"Legacy\" flag is meant to reproduce the Legacy results, ignoring the setting of the BESTCANDCOMPARATOR: ",BESTCANDCOMPARATOR
    BESTCANDCOMPARATOR = "byBestZ1bestZ2"


# Set to True to make candidates with the full combinatorial of loose leptons (for debug; much slower)
try:
    KEEPLOOSECOMB
except NameError:
    KEEPLOOSECOMB = False

# The isolation cuts for electrons and muons. FIXME: there is an hardcoded instance of these values in src/LeptonIsoHelper.cc !!
ELEISOCUT = "0.35"
MUISOCUT = "0.35"


# Which FSR mode to use: "Legacy" or "RunII"
try:
    FSRMODE
except NameError:
    FSRMODE = "RunII"


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
    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
    if IsMC:
        process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'
    else:
        process.GlobalTag.globaltag = '76X_dataRun2_v15'

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
# process.hltFilterMuEle = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
# process.hltFilterMuEle2 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
# process.hltFilterMuEle3 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
# process.hltFilterTriEle = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
# process.hltFilterTriMu  = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
# process.hltFilterSingleEle = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterDiMu.TriggerResultsTag  = cms.InputTag("TriggerResults","","HLT")
process.hltFilterDiEle.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
# process.hltFilterMuEle.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
# process.hltFilterMuEle2.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
# process.hltFilterMuEle3.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
# process.hltFilterTriEle.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
# process.hltFilterTriMu.TriggerResultsTag  = cms.InputTag("TriggerResults","","HLT")
# process.hltFilterSingleEle.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterDiMu.throw  = cms.bool(False) #FIXME: beware of this!
process.hltFilterDiEle.throw = cms.bool(False) #FIXME: beware of this!
# process.hltFilterMuEle.throw = cms.bool(False) #FIXME: beware of this!
# process.hltFilterMuEle2.throw = cms.bool(False) #FIXME: beware of this!
# process.hltFilterMuEle3.throw = cms.bool(False) #FIXME: beware of this!
# process.hltFilterTriEle.throw = cms.bool(False) #FIXME: beware of this!
# process.hltFilterTriMu.throw  = cms.bool(False) #FIXME: beware of this!
# process.hltFilterSingleEle.throw = cms.bool(False) #FIXME: beware of this!

# MuEG

if (LEPTON_SETUP == 2011):
    # DoubleEle
    process.hltFilterDiEle.HLTPaths = ["HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*","HLT_TripleEle10_CaloIdL_TrkIdVL_v*"] # to run on DATA/MC 2011
    if (IsMC):
        # DoubleMu
        process.hltFilterDiMu.HLTPaths = ["HLT_Mu17_Mu8_v*"] # to run on MC 2011    
  #      process.hltFilterMuEle.HLTPaths = ["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v*","HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*"] # to run on MC 2011

    else :
        process.hltFilterDiMu.HLTPaths = ["HLT_DoubleMu7_v*", "HLT_Mu13_Mu8_v*", "HLT_Mu17_Mu8_v*"] # to run on DATA 2011 NB: Emulation is needed [FIXME]
  #     process.hltFilterMuEle.HLTPaths = ["HLT_Mu17_Ele8_CaloIdL_v*", "HLT_Mu8_Ele17_CaloIdL_v*"] # For run 1-167913
  #      process.hltFilterMuEle2.HLTPaths = ["HLT_Mu17_Ele8_CaloIdL_v*", "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*"] # run 167914-175972
  #      process.hltFilterMuEle3.HLTPaths = ["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v*","HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*"] # : run 175973-180252
  #      process.triggerMuEle2  = cms.Path(process.hltFilterMuEle2)
  #      process.triggerMuEle3  = cms.Path(process.hltFilterMuEle3)

elif (LEPTON_SETUP == 2012):
    # DoubleEle
    process.hltFilterDiEle.HLTPaths = ["HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"] # to run on DATA/MC 2012
    process.hltFilterDiMu.HLTPaths = ["HLT_Mu17_Mu8_v*", "HLT_Mu17_TkMu8_v*"] # to run on DATA/MC 2012
  #  process.hltFilterMuEle.HLTPaths = ["HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*","HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"] # to run on DATA/MC 2012
  #  process.hltFilterTriEle.HLTPaths = ["HLT_Ele15_Ele8_Ele5_CaloIdL_TrkIdVL_v*"]
  #  process.triggerTriEle  = cms.Path(process.hltFilterTriEle)

elif (LEPTON_SETUP == 2015):
    #FIXME: For now, the MC paths are the ones used in the RunIISpring15DR74 MC samples for 25ns,7e33 conditions.
    process.hltFilterDiEle.HLTPaths = ["HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"]
    #process.hltFilterDiEle.HLTPaths = ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"] #for 25ns,14e33 (not necessary in 2015 data)
    process.hltFilterDiMu.HLTPaths = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*"]
  #  process.hltFilterMuEle.HLTPaths = ["HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v*","HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v*"]
    #process.hltFilterMuEle.HLTPaths = ["HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v*","HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v*"] #for 25ns,14e33 (not necessary in 2015 data)
  #  process.hltFilterTriEle.HLTPaths = ["HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*"]
  #  process.hltFilterTriMu.HLTPaths = ["HLT_TripleMu_12_10_5_v*"]
  #  process.hltFilterSingleEle.HLTPaths = ["HLT_Ele23_WPLoose_Gsf_v*"]
  #  process.triggerTriEle = cms.Path(process.hltFilterTriEle)
  #  process.triggerTriMu  = cms.Path(process.hltFilterTriMu )
  #  process.triggerSingleEle = cms.Path(process.hltFilterSingleEle)

process.triggerDiMu   = cms.Path(process.hltFilterDiMu)
process.triggerDiEle  = cms.Path(process.hltFilterDiEle)
# process.triggerMuEle  = cms.Path(process.hltFilterMuEle)


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
                                         identifier = cms.string("MC_76X_13TeV"),
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
    cut = cms.string("pt>10 && abs(eta)<2.4 && (isGlobalMuon || (isTrackerMuon && numberOfMatches>0)) && muonBestTrackType!=2")
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
    cut = cms.string("pt>10 && abs(eta)<2.5")
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


##--- Set up electron ID (VID framework)
#from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
## turn on VID producer, indicate data format to be DataFormat.MiniAOD, as appropriate
#dataFormat = DataFormat.MiniAOD
#switchOnVIDElectronIdProducer(process, dataFormat)
## define which IDs we want to produce
#my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff']
## add them to the VID producer
#for idmod in my_id_modules:
#    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
## and don't forget to add the producer 'process.egmGsfElectronIDSequence' to the path, i.e. process.electrons


process.bareSoftElectrons = cms.EDFilter("PATElectronRefSelector",
   src = cms.InputTag("calibratedPatElectrons"),
   cut = cms.string("pt>10 && abs(eta)<2.5")
   )

process.softElectrons = cms.EDProducer("EleFiller",
   src    = cms.InputTag("bareSoftElectrons"),
   sampleType = cms.int32(SAMPLE_TYPE),          
   setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
   cut = cms.string("userFloat('dxy')<0.5 && userFloat('dz')<1"),# removed cut because variable is in Spring15 ID  && userFloat('missingHit')<=1"),
   flags = cms.PSet(
        ID = cms.string("userFloat('isBDT')"),
        isSIP = cms.string(SIP),
        isGood = cms.string(GOODLEPTON),
        isIsoFSRUncorr  = cms.string("userFloat('combRelIsoPF')<"+ELEISOCUT)
#       Note: passCombRelIsoPFFSRCorr is currently set in LeptonPhotonMatcher for new FSR strategy; in ZZCandidateFiller for the old one
        ),
   #mvaValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"), # (when running VID)
   )


#process.electrons = cms.Sequence(process.selectedSlimmedElectrons + process.calibratedPatElectrons + process.egmGsfElectronIDSequence + process.bareSoftElectrons + process.softElectrons) # (use this version when running VID)
process.electrons = cms.Sequence(process.selectedSlimmedElectrons + process.calibratedPatElectrons + process.bareSoftElectrons + process.softElectrons)

# Handle special cases
if ELEREGRESSION == "None" and (ELECORRTYPE == "None" or BUNCH_SPACING == 50) :   # No correction at all. Skip correction modules.
    process.bareSoftElectrons.src = cms.InputTag('slimmedElectrons')
    process.electrons = cms.Sequence(process.bareSoftElectrons + process.softElectrons)

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
    sampleType = cms.int32(SAMPLE_TYPE),                     
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    photonSel = cms.string(FSRMODE)  # "skip", "passThrough", "Legacy", "RunII"
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

Z2PRESEL    = ("mass > 40 && mass < 150") # Note: this is without ISO

BESTZ2_AMONG = ( Z2PRESEL ) ## && userFloat('d1.passCombRelIsoPFFSRCorr')" )

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
        Z1Presel = cms.string(Z1PRESEL),
    )
)


# ll, any combination of flavour/charge, for control regions only
process.bareZjjCand = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('cleanJets cleanJets'),
    cut = cms.string('pt>100 && mass>40 && mass <150'), # protect against ghosts
    checkCharge = cms.bool(False)
)
process.ZjjCand = cms.EDProducer("ZCandidateFiller",
    src = cms.InputTag("bareZjjCand"),
    sampleType = cms.int32(SAMPLE_TYPE),  
    embedDaughterFloats = cms.untracked.bool(True),                   
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    bestZAmong = cms.string(BESTZ2_AMONG),
    FSRMode = cms.string("skip"), # "skip", "Legacy", "RunII"
    flags = cms.PSet(
        GoodLeptons = cms.string('0'),
        Z1Presel = cms.string(Z2PRESEL),
    )
)


### ----------------------------------------------------------------------
### TriLeptons (for fake rate)
### ----------------------------------------------------------------------
#FIXME: to be reviseed since ISOLATION is no longer included in isBestZ and Z1Presel!
# Z_PLUS_LEP_MIJ=("sqrt(pow(daughter(0).daughter({0}).energy+daughter(1).energy, 2) - " +
#                "     pow(daughter(0).daughter({0}).px    +daughter(1).px, 2) -" +  
#                 "     pow(daughter(0).daughter({0}).py    +daughter(1).py, 2) -" +  
#                 "     pow(daughter(0).daughter({0}).pz    +daughter(1).pz, 2))")

# process.ZlCand = cms.EDProducer("PATCandViewShallowCloneCombiner",
#     decay = cms.string('ZCand softLeptons'),
#     cut = cms.string("deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).eta, daughter(1).phi)>0.02 &&" + # Ghost suppression
#                      "deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).eta, daughter(1).phi)>0.02 &&" +
#                      ("(daughter(0).daughter(0).charge == daughter(1).charge || %s > 4) && " % ( Z_PLUS_LEP_MIJ.format(0))) +       # mLL>4 for the OS pair (Giovanni's impl)
#                      ("(daughter(0).daughter(1).charge == daughter(1).charge || %s > 4) && " % ( Z_PLUS_LEP_MIJ.format(1))) +
#                      "daughter(0).masterClone.userFloat('isBestZ') &&" +
#                      "daughter(0).masterClone.userFloat('Z1Presel')"
#                      ),
#     checkCharge = cms.bool(False)
# )


### ----------------------------------------------------------------------
### Quadrileptons: combine/merge leptons into intermediate (bare) collections;
###                Embed additional user variables into final collections
### ----------------------------------------------------------------------

FOURGOODLEPTONS    =  ("userFloat('d1.GoodLeptons')" +
                #        "&& userFloat('d0.worstEleIso') <" + ELEISOCUT +
                       "&& userFloat('d1.worstEleIso') <" + ELEISOCUT +
                #        "&& userFloat('d0.worstMuIso') <" + MUISOCUT +
                       "&& userFloat('d1.worstMuIso') <" + MUISOCUT
                       ) #ZZ made of 2 tight leptons passing SIP and ISO


Z1MASS            = "daughter('Z1').mass>70 && daughter('Z1').mass<110"
Z2MASS            = "daughter('Z2').mass>55 && daughter('Z2').mass<120 && daughter('Z2').pt>100" # (was > 4 in Synch) to deal with m12 cut at gen level #FIXME
#MLL3On4_12        = "userFloat('mZa')>12" # mll>12 on 3/4 pairs; 
#MLLALLCOMB        = "userFloat('mLL6')>4" # mll>4 on 6/6 AF/AS pairs;
MLLALLCOMB        = "userFloat('mLL4')>4" # mll>4 on 4/4 AF/OS pairs;
SMARTMALLCOMB     = "userFloat('passSmartMLL')" # Require swapped-lepton Z2' to be >12 IF Z1' is SF/OS and closer to 91.1876 than mZ1
PT20_10           = ("userFloat('pt1')>40 && userFloat('pt2')>20") #20/10 on any of the 4 leptons
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
                    #  MLLALLCOMB      + "&&" +
                      PT20_10         + "&&" +
                      "mass>300"       
                   #   SMARTMALLCOMB + "&&" +
                   #   "daughter('Z2').mass>12"
                      )

    SR = BESTCAND_AMONG


else:
    print "Please choose one of the following string for SELSETUP: 'Legacy', 'allCutsAtOnceButMZ2', 'allCutsAtOnce', 'allCutsAtOncePlusMZb', 'allCutsAtOncePlusSmart'"
    sys.exit()


FULLSEL            = (SR      + "&&" +
                      M4l100)


# Ghost Suppression cut
NOGHOST4l = ("deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).daughter(0).eta, daughter(1).daughter(0).phi)>0.02 && "+
             "deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).daughter(1).eta, daughter(1).daughter(1).phi)>0.02 && "+
             "deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).daughter(0).eta, daughter(1).daughter(0).phi)>0.02 && "+
             "deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).daughter(1).eta, daughter(1).daughter(1).phi)>0.02 && "+  # protect against ghosts
             "deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(0).daughter(1).eta, daughter(0).daughter(1).phi)>0.4 && "+ # remove merged jets
             "abs(daughter(0).daughter(0).eta)<2.5 && "+
             "abs(daughter(0).daughter(1).eta)<2.5 && " + # jets only in tracker region
             "daughter(0).daughter(0).mass<60 && "+
             "daughter(0).daughter(1).mass<60" ) # remove merged jets


# Preselection of 4l candidates
LLLLPRESEL = NOGHOST4l # Just suppress candidates with overlapping leptons

#LLLLPRESEL = (NOGHOST4l + "&&" +
#              FOURGOODLEPTONS) # Only retain candidates with 4 tight leptons (ID, iso, SIP)


# ZZ Candidates

if FSRMODE == "Legacy":
    RECOMPUTEISOFORFSR = True
elif FSRMODE == "RunII":
    RECOMPUTEISOFORFSR = False



process.bareZZCand= cms.EDProducer("PATCandViewShallowCloneCombiner",
    decay = cms.string('ZjjCand ZCand'),
    cut = cms.string(LLLLPRESEL),
    checkCharge = cms.bool(True)
)

process.ZZCand = cms.EDProducer("ZZjjCandidateFiller",
    src = cms.InputTag("bareZZCand"),
    sampleType = cms.int32(SAMPLE_TYPE),                    
    setup = cms.int32(LEPTON_SETUP),
    superMelaMass = cms.double(SUPERMELA_MASS),
    isMC = cms.bool(IsMC),
    bestCandAmong = cms.PSet(isBestCand = cms.string(BESTCAND_AMONG)),
    bestCandComparator = cms.string(BESTCANDCOMPARATOR),
    ZRolesByMass = cms.bool(False),
    recomputeIsoForFSR = cms.bool(RECOMPUTEISOFORFSR),
    flags = cms.PSet(
        GoodLeptons =  cms.string(FOURGOODLEPTONS),
        Z2Mass  = cms.string(Z2MASS),
        MAllComb = cms.string(MLLALLCOMB),
        SR = cms.string(SR),
        FullSel70 = cms.string(SR), #Obsolete, use "SR"
        FullSel = cms.string(FULLSEL),
    )
)



### ----------------------------------------------------------------------
### Jets
### ----------------------------------------------------------------------


### Load JEC
process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *

if IsMC: 
    process.jec = cms.ESSource("PoolDBESSource",
        DBParameters = cms.PSet(
            messageLevel = cms.untracked.int32(1)
        ),
        timetype = cms.string('runnumber'),
        toGet = cms.VPSet(
        cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_MC_AK4PFchs'),
            label  = cms.untracked.string('AK4PFchs')
            ),
        ), 
        connect = cms.string('sqlite_fip:ZZAnalysis/AnalysisStep/data/JEC/Fall15_25nsV2_MC.db'),
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
            tag    = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_DATA_AK4PFchs'),
            label  = cms.untracked.string('AK4PFchs')
            ),
        ), 
        connect = cms.string('sqlite_fip:ZZAnalysis/AnalysisStep/data/JEC/Fall15_25nsV2_DATA.db'),
    )
## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

### reapply JEC
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
process.patJetCorrFactorsReapplyJEC = patJetCorrFactorsUpdated.clone(
  src = cms.InputTag("slimmedJets"),
  levels = ['L1FastJet','L2Relative','L3Absolute'],
  payload = 'AK4PFchs' )
if not IsMC:
    process.patJetCorrFactorsReapplyJEC.levels = ['L1FastJet','L2Relative','L3Absolute','L2L3Residual']

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
process.patJetsReapplyJEC = patJetsUpdated.clone(
  jetSource = cms.InputTag("slimmedJets"),
  jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
  )


# embed q/g likelihood
process.load('RecoJets.JetProducers.QGTagger_cfi')
#process.QGTagger.srcJets = cms.InputTag( 'slimmedJets' )
process.QGTagger.srcJets = cms.InputTag( 'patJetsReapplyJEC' )
process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')

process.dressedJets = cms.EDProducer("JetFiller",
    #src = cms.InputTag("slimmedJets"),
    src = cms.InputTag("patJetsReapplyJEC"),
    cut = cms.string("pt>30 && abs(eta)<4.7"),
    bTaggerName = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
    jecType = cms.string("AK4PFchs"),
    flags = cms.PSet(
        isBtagged = cms.string("userFloat('bTagger')>0.800"),
    )
)


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

#process.Jets = cms.Path( process.QGTagger + process.dressedJets )
# reapply JEC
process.Jets = cms.Path( process.patJetCorrFactorsReapplyJEC + process.patJetsReapplyJEC + process.QGTagger + process.dressedJets )

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
       process.bareZjjCand       + process.ZjjCand     +  
       process.bareZZCand        + process.ZZCand
    )

# Optional sequence to build control regions. To get it, add
#process.CRPath = cms.Path(process.CRZl) # only trilepton
#OR
#process.CRPath = cms.Path(process.CR)   # trilep+4lep CRs

# process.CRZl = cms.Sequence(
#        process.bareZCand         + process.ZCand     +  
#        process.ZlCand            
#    )

# process.CR = cms.Sequence(
#        process.bareZCand         + process.ZCand     +  
#        process.bareLLCand        + process.LLCand    +
#        process.bareZLLCand       + process.ZLLCand   +
#        process.ZlCand            
#    )

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



