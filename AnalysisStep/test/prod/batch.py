#!/bin/env python

# Adapted from /CMGTools/Production/scripts/crisBatch.py

import sys
import imp
import copy
import os
import shutil
import pickle
import math
from CMGTools.Production.batchmanager import BatchManager
from CMGTools.Production.datasetToSource import *

def chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]

def split(comps):
    # import pdb; pdb.set_trace()
    splitComps = []
    for comp in comps:
        chunkSize = 1
        numJobs = 0
### Interpret splitFactor = # of jobs:
#        if hasattr( comp, 'splitFactor') and comp.splitFactor>1:
#            chunkSize = len(comp.files) / comp.splitFactor
#             if len(comp.files) % comp.splitFactor:
#                 chunkSize += 1
### Interpret splitFactor = #files/job
        if hasattr( comp, 'splitFactor') and comp.splitFactor<len(comp.files):
            chunkSize = comp.splitFactor
            for ichunk, chunk in enumerate( chunks( comp.files, chunkSize)):
                numJobs=numJobs+1
                newComp = copy.deepcopy(comp)
                newComp.files = chunk
                newComp.samplename = comp.name
                newComp.name = '{name}_Chunk{index}'.format(name=newComp.name,
                                                       index=ichunk)
                splitComps.append( newComp )
        else:
            numJobs=1
            splitComps.append( comp )
        print comp.name, ": files=", len(comp.files), ' chunkSize=', chunkSize, "jobs=", numJobs
    return splitComps


def batchScriptCERN( index, remoteDir=''):
   '''prepare the LSF version of the batch script, to run on LSF'''
#   print "INDEX", index
#   print "remotedir", remoteDir
   script = """#!/bin/tcsh
#BSUB -q 8nh
#BSUB -o job_%J.txt
#ulimit -v 3000000
limit
cat /proc/cpuinfo
cat /proc/meminfo
cd $CMSSW_BASE/src
cmsenv
cd -
echo 'Environment:'
echo
env
echo
echo 'Copying' ${LS_SUBCWD} to ${PWD} 
cp -rf $LS_SUBCWD .
echo '...done'
echo
echo Workdir content:
ls -l
echo
cd `find . -type d | grep /`
pwd
echo 'Running at:' `date`
cmsRun run_cfg.py >& log.txt
set cmsRunStatus=$?
echo 'cmsRun done at: ' `date` with exit status: $cmsRunStatus
if ( $cmsRunStatus != 0 ) echo $cmsRunStatus > exitStatus.txt
gzip log.txt
echo
echo 'ls: '
pwd
ls -l
echo
echo 'Sending the job directory back...'
cp *.root *.txt *.gz $LS_SUBCWD
if ( -z ZZ4lAnalysis.root ) then
 echo 'Empty file:  ZZ4lAnalysis.root'
 if ( -s ZZ4lAnalysis.root ) then
   echo Retrying...
   sleep 10
   cp *.root $LS_SUBCWD
 endif
endif
echo 'destination dir: ls: '
cd $LS_SUBCWD
pwd
ls -l
setenv ROOT_HIST 0
if ( -s ZZ4lAnalysis.root ) then
 root -q -b '${CMSSW_BASE}/src/ZZAnalysis/AnalysisStep/test/prod/rootFileIntegrity.r(\"ZZ4lAnalysis.root\")'
else
 echo moving empty file
 mv ZZ4lAnalysis.root ZZ4lAnalysis.root.empty
endif

echo '...done at' `date`
exit $cmsRunStatus
""" 
   return script


def batchScriptLocal( index ):
   '''prepare a local version of the batch script, to run using nohup'''

   script = """#!/bin/bash
echo 'running'
cmsRun run_cfg.py
""" 
   return script


def tuneProcess(process):
    return;
            
class MyBatchManager( BatchManager ):
   '''Batch manager specific to cmsRun processes.''' 
         
   def PrepareJobUser(self, jobDir, value ):
       '''Prepare one job. This function is called by the base class.'''
#       print value
#       print splitComponents[value]

       # import pdb; pdb.set_trace()
       #prepare the batch script
       scriptFileName = jobDir+'/batchScript.sh'
       scriptFile = open(scriptFileName,'w')
       # storeDir = self.remoteOutputDir_.replace('/castor/cern.ch/cms','')
       mode = self.RunningMode(options.batch)
       if mode == 'LXPLUS':
           scriptFile.write( batchScriptCERN( value ) )
       elif mode == 'LOCAL':
           scriptFile.write( batchScriptLocal( value ) ) 
       scriptFile.close()
       os.system('chmod +x %s' % scriptFileName)

       # Tune can be:
       # for DATA: "DoubleEle", "DoubleMu","MuEG"
       # for MC:   "", "MCAllEvents", "DYJets_B", "DYJets_NoB"
       #
       # setup can be 2011 or 2012 

       SAMPLENAME  = splitComponents[value].samplename
       tune  = splitComponents[value].tune
       setup = splitComponents[value].setup

       
       # Set global parameters
       IsMC = True
       PD = ""
       MCFILTER = ""
       SUPERMELA_MASS = 125.6
       if tune != "" :
           tunes = tune.split(",")
           if ("DoubleMu" in tunes) or ("DoubleEle" in tunes) or ("MuEG" in tunes):
               if len(tunes)!=1 : # Do not allow other options for data
                   raise ValueError, "invalid tune", tune
               IsMC = False
               PD = tune
               
           if ("DYJets_B" in tunes) or  ("DYJets_NoB" in tunes) :
               MCFILTER = "HF" # Note: further customization of the filter module is done below
           if ("MH126" in tunes) :
               SUPERMELA_MASS = 126
           # customization for "MCAllEvents" is done below

           #FIXME: should check tunes for consistency

       print SAMPLENAME, tune, IsMC, PD, MCFILTER, SUPERMELA_MASS

       # Read CFG file so that it is customized with the above globals
       namespace = {'IsMC':IsMC, 'PD':PD, 'MCFILTER':MCFILTER, 'SUPERMELA_MASS':SUPERMELA_MASS, 'SAMPLENAME':SAMPLENAME}
       execfile(cfgFileName,namespace)
#       handle = open(cfgFileName, 'r')
#       cfo = imp.load_source("pycfg", cfgFileName, handle)
#       process=cfo.process
#       handle.close()                  

       process = namespace.get('process') 
       process.source.fileNames = splitComponents[value].files

       if "MCAllEvents" in tune:
           process.ZZ4muTree.skipEmptyEvents = False
           process.ZZ4eTree.skipEmptyEvents = False
           process.ZZ2e2muTree.skipEmptyEvents = False

       #tuneProcess( process )
       cfgFile = open(jobDir+'/run_cfg.py','w')
#       cfgFile.write('import FWCore.ParameterSet.Config as cms\n\n')
       # cfgFile.write('import os,sys\n')
       cfgFile.write( process.dumpPython() )
       cfgFile.write( '\n' )
       del namespace
       del process

       # Tune process
       # JSON for data
       if IsMC==False :
           if setup==2011 :
               cfgFile.write( open('json_2011.py').read() ) 
           elif setup==2012 :
               cfgFile.write( open('json_2012.py').read() )
           else :
                raise ValueError, "invalid setup", setup
        
       if "DYJets_B" in tune :
            cfgFile.write( 'process.HF = cms.Path(process.heavyflavorfilter)\n\n' )
       elif "DYJets_NoB" in tune :
            cfgFile.write( 'process.HF = cms.Path(~process.heavyflavorfilter)\n\n' )
       cfgFile.close()


class Component(object):

    def __init__(self, name, user, dataset, pattern, splitFactor, tune, setup ):
        self.name = name
        print "checking "+self.name
        self.files = datasetToSource( user, dataset, pattern).fileNames # , True for readCache (?)
        self.splitFactor = splitFactor
        self.tune = tune
        self.setup = setup
        
      
if __name__ == '__main__':
    batchManager = MyBatchManager()
    batchManager.parser_.usage="""
    %prog [options] <cfgFile>

    Run Colin's python analysis system on the batch.
    Job splitting is determined by your configuration file.
    """

    options, args = batchManager.ParseOptions()

    cfgFileName = args[0]

    handle = open(cfgFileName, 'r')
    cfo = imp.load_source("pycfg", cfgFileName, handle)
    setup = cfo.LEPTON_SETUP
    components = [ Component(na, us, da, pattern, sp, tune, setup) for na, us,da,pattern,sp,tune in cfo.samples ]
    handle.close()


    splitComponents = split( components )
    listOfValues = range(0, len(splitComponents))
    listOfNames = [comp.name for comp in splitComponents]

    batchManager.PrepareJobs( listOfValues, listOfNames )

    waitingTime = 0.05
    batchManager.SubmitJobs( waitingTime )

