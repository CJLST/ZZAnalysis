
"""
This configuration is just an interface to the analyzer.py configuration
where all the sequences, paths are defined. This config can be used with 
command line options and is sutable for easy T3Cooker production at T3
clusters. 
"""

import FWCore.ParameterSet.Config as cms

from ZZAnalysis.AnalysisStep.Options import *
# and avoid making many cfg files with just few different parameters
#-----------------------------------------------------------------------------------------------------
#EXAMPLE; cmsRun zzPATSkim.py maxEvents=10 isMC=0 dataset=Jan16ReReco globaltag=FT_R_42_V24
#HELP: python ZZ4lAnalyzer_cfg.py --help  (display list of all the options)
#------------------------------------------------------------------------------------------------------
options = getAnalysisOptions()

options.maxEvents = -1
options.isMC = 0
options.goodlumi = ''
options.setup = 2012
options.tune = ''
#After parsing, options cannot be changed
options.parseArguments()

print "cmsRun with options: "
print "===================="
print options

# Tune can be:
# for DATA: "DoubleEle", "DoubleMu","MuEG"
# for MC:   "", "MCAllEvents", "DYJets_B", "DYJets_NoB"
#
# setup can be 2011 or 2012 
#
# elecorrtype can be:
# for DATA 2012: "13JulReReco", "06AugReReco", "24AugReReco", "Prompt2012C" -> pass "2012Jul13ReReco" to the module
# for DATA 2011: "Jan16ReReco"
# for MC:  "Summer12_53X" -> pass "Summer12_DR53X_HCP2012" to the module
#          "Summer12_52X" (not supported)
#          "Fall11"

tune  = options.tune
setup = options.setup
elecorrtype = options.elecorrtype

# Set global parameters
IsMC = bool(options.isMC)
IsMC = True
if tune == "DoubleEle" or tune == "DoubleMu" or tune == "MuEG" :
    IsMC = False

PD = ""
MCFILTER = ""
if tune != "" :
    if tune == "DYJets_B" or  tune == "DYJets_NoB" :
        MCFILTER = "HF" # customization is done below
    elif tune == "DoubleMu" or tune == "DoubleEle" or tune == "MuEG" :
        PD = tune
    elif tune == "MCAllEvents" :
        PD = "" # customization is done below
    else :
        raise ValueError, "invalid tune", tune

#print "tune=%(tune)s, IsMC=%(IsMC)d, PD=%(PD)s, MCFILTER=%(MCFILTER)s, elecorrtype=%(elecorrtype)s" %locals()

if elecorrtype == "13JulReReco" or elecorrtype == "06AugReReco" or elecorrtype == "24AugReReco" or elecorrtype == "Prompt2012C" :
    elecorrInputDataset = "2012Jul13ReReco"
elif elecorrtype == "Summer12_53X" :
    elecorrInputDataset = "Summer12_DR53X_HCP2012"
else :
    elecorrInputDataset = elecorrtype
 
#------------------------------------------------------------------------------------------------------
# read the configuration analyzer.py with global variables defined in namespace
#   
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"
# Read CFG file so that it is customized with the above globals
namespace = {'IsMC':IsMC, 'PD':PD, 'MCFILTER':MCFILTER, 'ELECORRTYPE':elecorrInputDataset}
print "Namespace: ",
print namespace
execfile(PyFilePath +"analyzer.py",namespace) #read the config file
process = namespace.get('process') 
#------------------------------------------------------------------------------------------------------
    
if tune == "MCAllEvents" :
    process.ZZ4muTree.skipEmptyEvents = False
    process.ZZ4eTree.skipEmptyEvents = False
    process.ZZ2e2muTree.skipEmptyEvents = False

if tune == "DYJets_B" :
    process.HF = cms.Path(process.heavyflavorfilter)
elif tune == "DYJets_NoB" :
    process.HF = cms.Path(~process.heavyflavorfilter)
    


#------------------------------------------------------------------------------------------------------
# setup JSON in case of data 
#
json = options.goodlumi
json = "USE_DEFAULT_JSON" #temporary 
process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
if IsMC == False:
    if json!='' and json!='USE_DEFAULT_JSON':
        print "JSON: " + options.goodlumi
        import FWCore.PythonUtilities.LumiList  as LumiList
        myLumis = LumiList.LumiList(filename = json).getCMSSWString().split(',')
        process.source.lumisToProcess.extend(myLumis)
    elif json=='USE_DEFAULT_JSON':
        if tune == "DoubleEle" or tune == "DoubleMu" or tune == "MuEG" : 
	  if setup==2011 :
	      execfile(PyFilePath + "prod/json_2011.py", namespace) 
	  elif setup==2012 :
	      execfile(PyFilePath + "prod/json_2012.py", namespace) 
	  else :
	      raise ValueError, "invalid setup", setup
    else :
        print "***WARNING: Running on data with no JSON goodlumi selection."
#print process.source.lumisToProcess
#print len(process.source.lumisToProcess)

#------------------------------------------------------------------------------------------------------

    
#------------------------------------------------------------------------------------------------------
#setup input and output file names
#
process.maxEvents.input = cms.untracked.int32(options.maxEvents) 
#print process.maxEvents.input 

if len(options.files) > 0 :
    process.source.fileNames =  cms.untracked.vstring(options.files)
print "Input files: "    
print process.source.fileNames 
    
filesuffix = str(options.tag)
if options.tag=='':
    filesuffix = ""
process.TFileService.fileName = cms.string(options.filePrepend+"ZZ4lAnalysis"+filesuffix+".root")
print "Output path: "+str(process.TFileService.fileName)
#------------------------------------------------------------------------------------------------------






          






