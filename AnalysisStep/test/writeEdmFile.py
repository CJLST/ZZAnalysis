### writeEdmFile.py -----------------------------------------------------------
###
### Example to write candidates and ancillary info in an output edm file
### 
### $Id: writeEdmFile.py,v 1.1 2012/10/16 20:20:49 namapane Exp $
###----------------------------------------------------------------------------


IsMC = True
LEPTON_SETUP = 2012
PD = ""
ELECORRTYPE = "Summer12_DR53X_HCP2012"
APPLYELEREGRESSION = False
APPLYELECALIB = False
APPLYMUCORR = False


# Get absolute path
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------

execfile(PyFilePath + "MasterPy/ZZ4lAnalysis.py")         # 2012 reference analysis


### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------
process.source.fileNames = cms.untracked.vstring(
    'root://lxcms00//data3/2012/HZZ_cmgTuple/synchHCP2/H125_53X_V5100.root' # V5_10_0 version
    )

process.maxEvents.input = -1

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


### ----------------------------------------------------------------------
### Additional Paths
### ----------------------------------------------------------------------
process.CRPath = cms.Path(process.CR)



### ----------------------------------------------------------------------
### Output module
### ----------------------------------------------------------------------

from ZZAnalysis.AnalysisStep.ZZ4lEventContent import ZZ4lEventContent

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('test.root'),
    outputCommands = ZZ4lEventContent.outputCommands,
    #SelectEvents = cms.untracked.PSet(
    #SelectEvents = cms.vstring('path')
    #)
    )

process.outp = cms.EndPath(process.out)
