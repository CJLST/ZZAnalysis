import FWCore.ParameterSet.Config as cms

process = cms.Process("llrTree")
# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load('Configuration.StandardSequences.GeometryDB_cff')
#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")


#----------------------------------------------------------
# Setup 'analysis'  options
#----------------------------------------------------------------------

from ZZAnalysis.AnalysisStep.Options import *
# you can add as much options in this file and make cfg highly flexible
# and avoid making many cfg files with just few different parameters
#-----------------------------------------------------------------------------------------------------
#EXAMPLE; cmsRun zzPATSkim.py maxEvents=10 isMC=0 dataset=Jan16ReReco globaltag=FT_R_42_V24
#HELP: python ZZ4lAnalyzer_cfg.py --help  (display list of all the options)
#------------------------------------------------------------------------------------------------------
options=getAnalysisStepOptions()  

options.maxEvents=-1
options.isMC=0
options.goodlumi=''
options.leptonSetup=2012
#After parsing, options cannot be changed
options.parseArguments()

print "cmsRun with options: "
print "===================="
print options


### ----------------------------------------------------------------------
### Set the GT
### ----------------------------------------------------------------------
my_global_tag=''
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#if bool(options.isMC) and options.globaltag=='NOT_PROVIDED':
  #my_global_tag='START44_V12'
#elif options.globaltag=='NOT_PROVIDED':
  #print "Not MC and no global tag provided. Using default"
  #my_global_tag='GR_P_V32'
#else:
  #my_global_tag=options.globaltag

#print options.globaltag
#process.GlobalTag.globaltag = my_global_tag+'::All'


if (options.leptonSetup == 2011) :
    if IsMC:
        process.GlobalTag.globaltag = 'START44_V13::All'
    else: 
#        process.GlobalTag.globaltag = 'GR_R_44_V5::All'
        process.GlobalTag.globaltag = 'GR_R_44_V15C::All'
        process.GlobalTag.connect = "sqlite_file:/afs/cern.ch/user/a/alcaprod/public/Alca/GlobalTag/GR_R_44_V15C.db"

else : 
    if options.isMC:
        process.GlobalTag.globaltag = 'START53_V10::All'   #for 53X MC
    else:
        process.GlobalTag.globaltag = 'GR_P_V41_AN1::All' # PrompReco(2012C-v2)


print "GLOBAL TAG: " + str(process.GlobalTag.globaltag)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
  #'/store/cernproduction/hzz4l/CMG/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/V5_B/PAT_CMG_V5_10_0/cmgTuple_26.root'
  '/store/cernproduction/hzz4l/CMG/DoubleElectron/Run2012B-13Jul2012-v1/AOD/V5/PAT_CMG_V5_10_0/cmgTuple_1.root'
  
),
		        lumisToProcess = cms.untracked.VLuminosityBlockRange()
)


if len(options.files) > 0 :
    process.source.fileNames =  cms.untracked.vstring(options.files)
#print process.source.fileNames

if not bool(options.isMC):
  if options.goodlumi!='':
    print "JSON: " + options.goodlumi
    import FWCore.PythonUtilities.LumiList  as LumiList
    myLumis = LumiList.LumiList(filename = options.goodlumi).getCMSSWString().split(',')
    #process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
    process.source.lumisToProcess.extend(myLumis)
  else :
      print "***WARNING: Running on data with no JSON goodlumi selection."

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

filesuffix = str(options.tag)
if options.tag=='':
    filesuffix = ""

process.TFileService = cms.Service('TFileService',
                                   fileName=cms.string(options.filePrepend+"llrNtuple"+filesuffix+".root" ),
                                   closeFileFast = cms.untracked.bool(True)
                                   )




process.goodPrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(True),
)

#--- Electron regression. Must be applied after BDT is recomputed
process.patElectronsWithRegression = cms.EDProducer("RegressionEnergyPatElectronProducer",
   debug = cms.untracked.bool(False),
   inputPatElectronsTag = cms.InputTag('patElectronsWithTrigger'),
   regressionInputFile = cms.string("EGamma/EGammaAnalysisTools/data/eleEnergyRegWeights_V1.root"),
   energyRegressionType = cms.uint32(1),
   rhoCollection = cms.InputTag('kt6PFJets:rho:RECO'),
   vertexCollection = cms.InputTag('goodPrimaryVertices')
   )



# Electron calibration. Default parameters are for calibration only.
process.calibratedElectrons = cms.EDProducer("CalibratedPatElectronProducer",
   inputPatElectronsTag = cms.InputTag("patElectronsWithRegression"),
   inputDataset = cms.string(options.dataset),
   applyCorrections=cms.int32(1),
   isMC    = cms.bool(bool(options.isMC)),
   isAOD   = cms.bool(True),
   updateEnergyError = cms.bool(True),
   debug   = cms.bool(False),
   )
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   calibratedElectrons = cms.PSet(
                                                       initialSeed = cms.untracked.uint32(1),
                                                       engineName = cms.untracked.string('TRandom3')
                                                       ),
                                                   )

process.bareSoftElectrons = cms.EDFilter("PATElectronRefSelector",
   src = cms.InputTag("calibratedElectrons"),
   cut = cms.string("pt>7 && abs(eta)<2.5 &&" +
                    "gsfTrack.trackerExpectedHitsInner.numberOfHits<=1"
                    )
   )
 
### ----------------------------------------------------------------------
### Lepton Cleaning (clean electrons collection from muons)
### ----------------------------------------------------------------------

process.cleanElectrons = cms.EDProducer("PATElectronCleaner",
    # pat electron input source
    src = cms.InputTag("bareSoftElectrons"),
    # preselection (any string-based cut for pat::Electron)
    preselection = cms.string(''),
    # overlap checking configurables
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
           src       = cms.InputTag("patMuonsWithTrigger"), # Start from loose lepton def
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string("(isGlobalMuon || userFloat('isPFMuon'))"), #
           deltaR              = cms.double(0.05),  
           checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), # overlaps don't cause the electron to be discared
        )
    ),
    # finalCut (any string-based cut for pat::Electron)
    finalCut = cms.string(''),
)

process.electrons = cms.Sequence(process.patElectronsWithRegression + process.calibratedElectrons + process.bareSoftElectrons + process.cleanElectrons)

APPLYELEREGRESSION = True
#APPLYELEREGRESSION = False
ELECORRTYPE = options.dataset
ELECORRTYPEDUMMY = ELECORRTYPE

# Handle special cases
if APPLYELEREGRESSION == False and ELECORRTYPE == "None" :   #No correction at all
    process.bareSoftElectrons.src = cms.InputTag('patElectronsWithTrigger')
    process.electrons = cms.Sequence(process.bareSoftElectrons)

elif APPLYELEREGRESSION == True and ELECORRTYPE == "None" : # Regression only
    process.calibratedElectrons.applyCorrections = 10
    # The damn module wants a meaningful value of inputDataset even if irrelevant!
    if (bool(options.isMC)) :
        process.calibratedElectrons.inputDataset = cms.string(ELECORRTYPEDUMMY) 
    else :
        process.calibratedElectrons.inputDataset = cms.string("Prompt") 

elif APPLYELEREGRESSION == False and ELECORRTYPE != "None" : #Calibration only
        process.calibratedElectrons.applyCorrections = 999
        process.calibratedElectrons.inputPatElectronsTag =  cms.InputTag('patElectronsWithTrigger')
        process.electrons = cms.Sequence(process.calibratedElectrons + process.bareSoftElectrons + process.cleanElectrons)

LEPTON_SETUP = options.leptonSetup
process.produceTree = cms.EDAnalyzer("NtupleProducer",
                                     MuonTag         = cms.InputTag('patMuonsWithTrigger'),
			       EleTag          = cms.InputTag('cleanElectrons'),
                                     PhotonTag       = cms.InputTag('cmgPhotonSel'),
                                     JetTag          = cms.InputTag('cmgPFJetSel'),
                                     SCTag           = cms.InputTag('SelectedEtSuperCluster'),
                                     VerticesTag     = cms.InputTag('goodPrimaryVertices'),
                                     HLTTag          = cms.InputTag('TriggerResults','','HLT'),
                                     TriggerEventTag = cms.InputTag('hltTriggerSummaryAOD::HLT'),
                                     HLTFilters  = cms.VInputTag( #HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL
                                                      'hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter::HLT',
                                                      'hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter::HLT',
                                                      'hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ::HLT',
                                                      #HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50
                                                      'hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter::HLT',
                                                      'hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter::HLT',
                                                      #HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50
                                                      'hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter::HLT',
                                                       #HLT_Ele15_Ele8_Ele5_CaloIdL_TrkIdVL
                                                      'hltL1NonIsoHLT3LegEleIdTripleElectronEt15Et8Et5EleIdDphiFilter::HLT'),
                                                   #   'hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4PMMassFilter::HLT'),
                                     MCTag           = cms.InputTag('genParticlesPruned'),
                                     isMC	= cms.bool(bool(options.isMC)),
                                     lepton_setup = cms.int32(LEPTON_SETUP)
                                     )




print "Setting path for LLR Trees production."
#process.p = cms.Path(process.goodPrimaryVertices+process.cleanElectrons+process.calibratedElectrons+process.produceTree)
process.p = cms.Path(process.goodPrimaryVertices + process.electrons + process.produceTree)

