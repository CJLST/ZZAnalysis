#----------------------------------------------------------------------
#
# Configuration for data, DoubleMu stream
#
#----------------------------------------------------------------------


# Load deafult job config
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"
execfile(PyFilePath + "analyzer.py")        

print 'Loading trackless_e_analyzer.py'
execfile(PyFilePath + "trackless_e_analyzer.py")        

print 'Loading loose_electron_analyzer.py'
execfile(PyFilePath + "loose_ele_analyzer.py")        



# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#process.appendPhotons.photonSel="skip"

#process.appendPhotons.photonSel = "passThrough"
#process.FSR  = cms.EDAnalyzer("FSRAnalyzer",)
#process.FSRPath = cms.EndPath(process.FSR)
