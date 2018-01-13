from ZZAnalysis.AnalysisStep.defaults import *

### ----------------------------------------------------------------------
###
### Example ntuplizer for lepton studies
###
###----------------------------------------------------------------------

# Set defaults for variables used in this file (in case they are not defined by a caller script)
declareDefault("PD", "", globals()) # "" for MC, "DoubleEle", "DoubleMu", or "MuEG" for data 
declareDefault("MCFILTER", "", globals())
declareDefault("XSEC", 1, globals())
declareDefault("GENXSEC", 1, globals())
declareDefault("GENBR", 1, globals())
LEPTON_SETUP = 2017  # current default = 2017 = Moriond2017
ELECORRTYPE = "None" # "None" to switch off
ELEREGRESSION = "None" # "None" to switch off
APPLYMUCORR = True  # Switch off muon scale corrections
APPLYJEC = False     # 
APPLYJER = False     #
RECORRECTMET = False #


# Get absolute path
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------
execfile(PyFilePath + "MasterPy/ZZ4lAnalysis.py")


### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------
process.maxEvents.input = -1

## remove any preselection (to be re-applied in the plotter if needed)
process.bareSoftMuons.cut = cms.string("")
process.softMuons.cut = cms.string("")
process.bareSoftElectrons.cut = cms.string("")
process.softElectrons.cut = cms.string("")


### ----------------------------------------------------------------------
### Output root file (monitoring histograms)
### ----------------------------------------------------------------------
process.TFileService=cms.Service('TFileService',
                                fileName=cms.string('ZZ4lAnalysis.root')
                                )


### ----------------------------------------------------------------------
### Ntuplizer
### ----------------------------------------------------------------------
process.LeptonTree = cms.EDAnalyzer('LeptonNtupleMaker',
     fileName = cms.untracked.string('lepTree'),
     genParticleSrc = cms.InputTag("prunedGenParticles"),
     electronSrc = cms.InputTag("softElectrons"),
     muonSrc = cms.InputTag("softMuons"),
)


### ----------------------------------------------------------------------
### Paths
### ----------------------------------------------------------------------

# overwrite master py: only the beginning of the lepton flow matters
process.triggerTriEle = cms.Path()
process.triggerTriMu = cms.Path()
process.triggerSingleEle = cms.Path()
process.triggerDiMu = cms.Path()
process.triggerDiEle = cms.Path()
process.triggerMuEle = cms.Path()
process.Jets = cms.Path()
process.Candidates = cms.Path(process.muons + process.electrons)

# ntuplize
process.trees = cms.EndPath(process.LeptonTree)
