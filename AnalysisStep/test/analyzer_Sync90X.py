
#DATA_TAG = "ReReco" # Change to PromptReco for Run2016 period H
LEPTON_SETUP = 2017  # current default = 2017 = Moriond2017
ELECORRTYPE = "None" # "None" to switch off
ELEREGRESSION = "None" # "None" to switch off
APPLYMUCORR = True  # Switch off muon scale corrections
APPLYJEC = True     #
APPLYJER = True     #
RECORRECTMET = True #
#KINREFIT = False    # control KinZFitter (very slow)
PROCESS_CR = False   # Uncomment to run CR paths and trees
#ADDLOOSEELE = True  # Run paths for loose electrons
#APPLYTRIG = False    # hack for samples missing correct triggers - use with caution
#KEEPLOOSECOMB = True # Do not skip loose lepton ZZ combinations (for debugging)
ADDZTREE = True      # Add tree for Z analysis

PD = ""
MCFILTER = ""

#For DATA: 
IsMC = False
PD = "DoubleMu"

# Get absolute path
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------

execfile(PyFilePath + "analyzer.py")
execfile(PyFilePath + "prod/pyFragments/RecoProbabilities.py")

if not IsMC:
	process.source.inputCommands = cms.untracked.vstring("keep *", "drop LHERunInfoProduct_*_*_*", "drop LHEEventProduct_*_*_*") ###FIXME In 9X this removes all collections for MC

### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------

process.source.fileNames = cms.untracked.vstring(


### ICHEP 16 sync files
    ## reHLT miniAODv2 files for sync
#    '/store/mc/RunIISummer16MiniAODv2/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0062F021-8DC0-E611-9C48-90B11CBCFFEA.root',
#    '/store/mc/RunIISpring16MiniAODv2/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/9AEF77A8-4B3B-E611-97F0-44A842240F8D.root',
                                                 ## '/store/mc/RunIISpring16MiniAODv2/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/20000/FE43CA0A-993A-E611-91FB-0242AC130002.root'
    ## High-mass reHLT sync file
#'/store/mc/RunIISpring16MiniAODv2/Graviton2PBToZZTo4L_width0p1_M-2000_13TeV-JHUgenV6-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/20000/EC5DACC2-7539-E611-9499-001C23BED42C.root'


### Moriond 17 sync files
    ## Low mass
#    '/store/mc/RunIISummer16MiniAODv2/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/22F32262-3FC5-E611-B373-D4AE526DEDB7.root',
#    '/store/mc/RunIISummer16MiniAODv2/WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/387FA719-E6CC-E611-A1F0-FA163E7D6032.root',
#    '/store/mc/RunIISummer16MiniAODv2/ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/50DE4DA2-1EC1-E611-9A3C-002590E7E010.root'

    ## High mass
#    '/store/mc/RunIISummer16MiniAODv2/Graviton2PMToZZTo4L_M750_13TeV_JHUgenV6_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/9AA3CE64-F3CB-E611-BB50-0025907B5048.root'

    ### 920 relval
#    '/store/relval/CMSSW_9_2_0/DoubleMuon/MINIAOD/91X_dataRun2_relval_v6_RelVal_doubMu2016B-v1/10000/363F2659-453C-E711-8B1C-0CC47A7C345E.root'

    ### Run2017C sample data file
	 '/store/data/Run2017B/DoubleMuon/MINIAOD/17Nov2017-v1/30000/0852E0CB-E7D7-E711-B2DA-0025905C3DCE.root',
							
	 ### 2017MC sample file
#	 '/store/mc/RunIIFall17MiniAOD/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v2/60000/267F8236-BFDA-E711-B4BF-008CFAC9422C.root',
    )

#process.calibratedPatElectrons.isSynchronization = cms.bool(True)
process.calibratedMuons.isSynchronization = cms.bool(True)

process.maxEvents.input = -1
#process.source.skipEvents = cms.untracked.uint32(5750)

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


### ----------------------------------------------------------------------
### Analyzer for Plots
### ----------------------------------------------------------------------


#process.source.eventsToProcess = cms.untracked.VEventRange("1:8670")

# Debug
process.dumpUserData =  cms.EDAnalyzer("dumpUserData",
     dumpTrigger = cms.untracked.bool(True),
     muonSrcs = cms.PSet(
#       slimmedMuons = cms.InputTag("slimmedMuons"),
        muons = cms.InputTag("appendPhotons:muons"),
     ),
     electronSrcs = cms.PSet(
#       slimmedElectron = cms.InputTag("slimmedElectrons"),
        electrons = cms.InputTag("appendPhotons:electrons"),
#        RSE = cms.InputTag("appendPhotons:looseElectrons"),
#        TLE = cms.InputTag("appendPhotons:electronstle"), #These are actually photons, should add a photonSrcs section for them.
     ),
     candidateSrcs = cms.PSet(
        Z     = cms.InputTag("ZCand"),
#        ZRSE     = cms.InputTag("ZCandlooseEle"),
#        ZTLE     = cms.InputTag("ZCandtle"),
        ZZ  = cms.InputTag("ZZCand"),
#        ZZRSE     = cms.InputTag("ZZCandlooseEle"),
#        ZZTLE     = cms.InputTag("ZZCandtle"),
#        ZLL  = cms.InputTag("ZLLCand"),
#        ZL  = cms.InputTag("ZlCand"),
     ),
     jetSrc = cms.InputTag("cleanJets"),
)

# Create lepton sync file
#process.PlotsZZ.dumpForSync = True;
#process.p = cms.EndPath( process.PlotsZZ)

# Keep all events in the tree, even if no candidate is selected
#process.ZZTree.skipEmptyEvents = False

# replace the paths in analyzer.py
#process.trees = cms.EndPath(process.ZZTree)

#Dump reconstructed variables
#process.appendPhotons.debug = cms.untracked.bool(True)
#process.fsrPhotons.debug = cms.untracked.bool(True)
#process.dump = cms.Path(process.dumpUserData)

#Print MC history
#process.mch = cms.EndPath(process.printTree)


#Monitor memory usage
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)
