LEPTON_SETUP = 2012

IsMC = False
PD = "MuEG"
MCFILTER = ""
ELECORRTYPE = "2012Jul13ReReco"
APPLYELEREGRESSION = True
APPLYELECALIB = True
APPLYMUCORR = True

import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"
execfile(PyFilePath + "analyzer.py")        
execfile(PyFilePath + "prod/json_2012.py")      



process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000




process.source.fileNames = cms.untracked.vstring(
    'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/MuEG/Run2012B-13Jul2012-v1/AOD/V5_B/PAT_CMG_V5_10_0/cmgTuple_1036.root',
    'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/MuEG/Run2012B-13Jul2012-v1/AOD/V5_B/PAT_CMG_V5_10_0/cmgTuple_770.root',

    )


