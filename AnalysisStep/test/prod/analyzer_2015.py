#----------------------------------------------------------------------
#
# Configuration for data, DoubleMu stream
#
#----------------------------------------------------------------------

LEPTON_SETUP = 2015

ELECORRTYPE   = "None"
ELEREGRESSION = "None"
APPLYMUCORR = False

# Load deafult job config
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"
execfile(PyFilePath + "analyzer.py")        

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
