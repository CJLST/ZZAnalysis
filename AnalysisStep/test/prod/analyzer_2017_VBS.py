#----------------------------------------------------------------------
#
# Configuration for ZZjj VBS analysis 
#
#----------------------------------------------------------------------


# Load deafult job config
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"
execfile(PyFilePath + "analyzer.py")        


# use SMP ZZ arbitrator (does not alter RSE/TLE)
process.ZLLCand.bestCandComparator = cms.string('byBestZ1bestZ2')
process.ZZCand.bestCandComparator = cms.string('byBestZ1bestZ2')

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#process.appendPhotons.photonSel="skip"

#process.appendPhotons.photonSel = "passThrough"
#process.FSR  = cms.EDAnalyzer("FSRAnalyzer",)
#process.FSRPath = cms.EndPath(process.FSR)
