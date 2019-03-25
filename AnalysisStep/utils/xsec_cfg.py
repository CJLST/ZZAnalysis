# run like
# cmsRun xsec_cfg.py datasetName=/WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM maxEvents=100000

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

from ZZAnalysis.AnalysisStep.eostools import listFiles

def get_max_files(DAS_name, max_files) :
  result = []
  file_names = listFiles(DAS_name, "dbs")
  xrd_prefix = 'root://cms-xrd-global.cern.ch/'
#  xrd_prefix = 'root://xrootd-cms.infn.it/'

  for i in range(min(len(file_names), max_files)) :
     if file_names[i]:
        result.append(file_names[i])

  return result




options = VarParsing ('analysis')

options.register('datasetName',
		'',
		VarParsing.multiplicity.singleton,
		VarParsing.varType.string,
		"DAS-style name of a primary dataset, e.g. /ZZTo4L_13TeV-sherpa/RunIISummer16MiniAODv2-PUMoriond17_v2_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM")

#options.register('maxEvents',
#		100000,
#		VarParsing.multiplicity.singleton,
#		VarParsing.varType.int,
#		"number of MC events to run over, more means less statistical uncertainty on xsection")



options.parseArguments()
process = cms.Process('XSec')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents),
)

process.load('FWCore.MessageService.MessageLogger_cfi')
# How often do you want to see the "Begin processing the .." 
process.MessageLogger.cerr.FwkReport.reportEvery = 1000



process.source = cms.Source ("PoolSource",
    fileNames = cms.untracked.vstring(*get_max_files(options.datasetName, 50)),
)
process.xsec = cms.EDAnalyzer("GenXSecAnalyzer")

process.ana = cms.Path(process.xsec)
process.schedule = cms.Schedule(process.ana)
