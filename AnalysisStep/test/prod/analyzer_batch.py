#----------------------------------------------------------------------
#
# Configuration for batch jobs
#
#----------------------------------------------------------------------
from past.builtins import execfile

#FIXME: switch off all corrections until they become available
APPLYMUCORR = False  # Switch off muon scale corrections
APPLYJEC = False     # Switch off JEC
APPLYJER = False     # Switch off JER
RECORRECTMET = False # Switch off MET corr

# Load deafult job config
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"
execfile(PyFilePath + "analyzer.py")

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
