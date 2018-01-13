
#DATA_TAG = "ReReco" # Change to PromptReco for Run2016 period H
LEPTON_SETUP = 2017  # current default = 2017 = Moriond2017
ELECORRTYPE = "None" # "None" to switch off
ELEREGRESSION = "None" # "None" to switch off
APPLYMUCORR = True  # Switch off muon scale corrections
APPLYJEC = False     #
APPLYJER = False     #
RECORRECTMET = False #
#KINREFIT = False    # control KinZFitter (very slow)
PROCESS_CR = False   # Uncomment to run CR paths and trees
#ADDLOOSEELE = True  # Run paths for loose electrons
APPLYTRIG = False    # hack for samples missing correct triggers - use with caution
#KEEPLOOSECOMB = True # Do not skip loose lepton ZZ combinations (for debugging)
ADDZTREE = True      # Add tree for Z analysis

PD = ""
MCFILTER = ""


# Get absolute path
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"


### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------

execfile(PyFilePath + "ntuplizerLeptons.py")

process.TFileService.fileName = cms.string('LeptonNtuple.root')

#process.source.inputCommands = cms.untracked.vstring("keep *", "drop LHERunInfoProduct_*_*_*", "drop LHEEventProduct_*_*_*")


### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------

process.source.fileNames = cms.untracked.vstring(

    ## Signal, Fall15 MiniAODv2
	 #'/store/mc/RunIIFall15MiniAODv2/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/024F30C4-96B9-E511-91A0-24BE05C616C1.root',
     
    ## Background, Fall15 MiniAODv2
    '/store/mc/RunIISpring16MiniAODv1/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/60218150-0F01-E611-BBEC-0CC47A78A41C.root',

    )

#process.calibratedPatElectrons.isSynchronization = cms.bool(True)
process.calibratedMuons.isSynchronization = cms.bool(True)

process.maxEvents.input = -1
#process.source.skipEvents = cms.untracked.uint32(5750)

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#process.source.eventsToProcess = cms.untracked.VEventRange("1:122490", "1:1343", "1:177684")


## put back preselection
#process.bareSoftMuons.cut = cms.string("")
#process.softMuons.cut = cms.string("")
#process.bareSoftElectrons.cut = cms.string("")
#process.softElectrons.cut = cms.string("")


# replace the paths in ntuplizerLeptons.py
#process.trees = cms.EndPath(process.LeptonTree)


