#----------------------------------------------------------------------
#
# Configuration for data, DoubleMu stream
#
#----------------------------------------------------------------------

LEPTON_SETUP = 2015

ELECORRTYPE = "RunII"
ELEREGRESSION = "None"
APPLYMUCORR = True

FSRMODE = "RunII"

# Load deafult job config
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"
execfile(PyFilePath + "analyzer2l2q.py")        

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#process.appendPhotons.photonSel="skip"

#process.appendPhotons.photonSel = "passThrough"
#process.FSR  = cms.EDAnalyzer("FSRAnalyzer",)
#process.FSRPath = cms.EndPath(process.FSR)
