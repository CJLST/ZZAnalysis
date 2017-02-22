#!/bin/env python

import sys
import imp
import copy
import os
import shutil
import pickle
import math
import pprint
import subprocess
from datetime import date
from optparse import OptionParser
from ZZAnalysis.AnalysisStep.eostools import *
from ZZAnalysis.AnalysisStep.readSampleInfo import *


def chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]

def split(comps):
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
        if hasattr( comp, 'splitFactor') : #and comp.splitFactor<len(comp.files)
            chunkSize = comp.splitFactor
            for ichunk, chunk in enumerate( chunks( comp.files, chunkSize)):
                numJobs=numJobs+1
                newComp = copy.deepcopy(comp)
                newComp.files = chunk
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
if ( ! $?LS_SUBCWD ) then
  #running interactively
  set LS_SUBCWD=`pwd`
  cd `mktemp -d`
endif
limit
cat /proc/cpuinfo
cat /proc/meminfo
cd $CMSSW_BASE/src
cmsenv
setenv X509_USER_PROXY ~/x509up_u${uid}
cd -
echo 'Environment:'
echo
env
echo
echo 'Copying' ${LS_SUBCWD} to ${PWD} 
cp -rf $LS_SUBCWD/ .
echo '...done'
echo
#echo Workdir content:
#ls -l
echo
cd `find . -type d | grep /`
pwd
echo 'Running at:' `date`
cmsRun run_cfg.py |& grep -v -e 'MINUIT WARNING' -e 'Second derivative zero' -e 'Negative diagonal element' -e 'added to diagonal of error matrix' > log.txt
set cmsRunStatus=$?
echo 'cmsRun done at: ' `date` with exit status: $cmsRunStatus
if ( $cmsRunStatus != 0 ) echo $cmsRunStatus > exitStatus.txt
gzip log.txt
echo
#echo 'ls: '
#pwd
#ls -l
#echo
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
#pwd
#ls -l
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

            
class MyBatchManager:
    '''Batch manager specific to cmsRun processes.''' 

    def __init__(self):        
        # define options and arguments ====================================
        self.parser_ = OptionParser()
        self.parser_.add_option("-o", "--output-dir", dest="outputDir",
                                help="Name of the local output directory for your jobs. This directory will be created automatically.",
                                default=None)

        self.parser_.add_option("-f", "--force", action="store_true",
                                dest="force", default=False,
                                help="Don't ask any questions, just over-write")

        self.parser_.add_option("-n", "--negate", action="store_true",
                                dest="negate", default=False,
                                help="create jobs, but does not submit the jobs.")

        self.parser_.add_option("-b", "--batch", dest="batch",
                                help="batch command. default is: 'bsub -q 8nh < batchScript.sh'. You can also use 'nohup < ./batchScript.sh &' to run locally.",
                                default="bsub -q 8nh < ./batchScript.sh")

        self.parser_.add_option("-p", "--pdf", dest="PDFstep",
                                help="Step of PDF systematic uncertainty evaluation. It could be 1 or 2.",
                                default=0)
       
        self.parser_.add_option("-s", "--secondary-input-dir", dest="secondaryInputDir",
                                help="Name of the local input directory for your PDF jobs",
                                default=None)

        self.parser_.add_option("-i", "--input", dest="cfgFileName",
                                help="input cfg",
                                default="analyzer_2015.py")

        self.parser_.add_option("-d", "--debug", action="store_true",
                                dest="verbose",default =False,
                                help="Activate verbose output",)


        (self.options_,self.args_) = self.parser_.parse_args()


        if ( len(self.args_)!=1) :
            print "Please specify sample.csv file to be used.\n"
            sys.exit(1)
             
        csvfile = self.args_[0]

        if (os.path.exists(csvfile) == False ):
            print "File", csvfile, "does not exist; please specify a valid one.\n"
            sys.exit(1)

        # Handle output directory
        outputDir = self.options_.outputDir
        if outputDir==None:
#             today = date.today()
#             outputDir = 'OutCmsBatch_%s' % today.strftime("%d%h%y_%H%M")
            gitrevision = subprocess.check_output(['git', "rev-parse", "--short", "HEAD"]) #revision of the git area where the command is exectuted
            outputDir = "PROD_" + csvfile.replace('.csv','') + "_"+gitrevision.rstrip()
            print 'output directory not specified, using %s' % outputDir
        self.outputDir_ = os.path.abspath(outputDir)
        self.workingDir = str(self.outputDir_)
        if( os.path.isdir(self.outputDir_) == True ):
            input = ''
            if not self.options_.force:
                while input != 'y' and input != 'n':
                    input = raw_input( 'The directory ' + self.outputDir_ + ' exists. Are you sure you want to continue? its contents will be overwritten [y/n] ' )
            if input == 'n':
                sys.exit(1)
            else:
                os.system( 'rm -rf ' + self.outputDir_)
        print 'Job folder: %s' % self.outputDir_
        self.mkdir( self.outputDir_ )
        with open(os.path.join(self.outputDir_, ".gitignore"), "w") as f:
            f.write("**/*\n")

        #FIXME check this???
        if not self.options_.secondaryInputDir == None:
            self.secondaryInputDir_ = os.path.abspath(self.options_.secondaryInputDir) + '/AAAOK/'
        else: 
            self.secondaryInputDir_ = None


    def mkdir( self, dirname ):
        try:
            os.makedirs(dirname)
        except OSError:
            pass
    

       
    def PrepareJobs(self, listOfValues, listOfDirNames=None):
        print 'PREPARING JOBS ======== '
        self.listOfJobs_ = []

        if listOfDirNames is None:
	    print 'here'
            for value in listOfValues:
                self.PrepareJob( value )
        else:
	    print 'yay!'
            for value, name in zip( listOfValues, listOfDirNames):
		print 'value=', value
		print 'name=', name
                self.PrepareJob( value, name )
        if batchManager.options_.verbose:
            print "list of jobs:"
            pp = pprint.PrettyPrinter(indent=4)
            pp.pprint( self.listOfJobs_)



    def PrepareJob( self, value, dirname=None):
       '''Prepare a job for a given value.
       
       calls PrepareJobUser, which should be overloaded by the user.
       '''
       print 'PrepareJob : %s' % value 
       dname = dirname
       if dname  is None:
           dname = 'Job_{value}'.format( value=value )
       jobDir = '/'.join( [self.outputDir_, dname])
       print '\t',jobDir 
       self.mkdir( jobDir )
       self.listOfJobs_.append( jobDir )
       if not self.secondaryInputDir_ == None: self.inputPDFDir_ = '/'.join( [self.secondaryInputDir_, dname])

       print 'dirname=',dirname
       print 'self.outputDir_=',self.outputDir_
       chunk_number = dirname[dirname.find('Chunk') + len('Chunk')]
        
       print 'chunk=', chunk_number
       chunk_number = float(chunk_number)
       if chunk_number == 0 :
           self.PrepareJobUserTemplate( jobDir, value )
       self.PrepareJobUserFromTemplate(jobDir, value)

    def PrepareJobUserFromTemplate(self, jobDir, value ):
        scriptFileName = jobDir+'/batchScript.sh'
        scriptFile = open(scriptFileName,'w')
        scriptFile.write( batchScriptCERN( value ) )
        scriptFile.close()
        os.system('chmod +x %s' % scriptFileName)
        template_name = splitComponents[value].samplename + 'run_template_cfg.py'

#	working_dir = os.path.dirname(self.outputDir_)


	template_file_name = '%s/%s'%(self.outputDir_, template_name) #splitComponents[value].samplename + '_run_template_cfg.py' 
#        shutil.copyfile(template_file_name, '%s/run_cfg.py'%jobDir)  
        new_job_path = '%s/run_cfg.py'%jobDir
	files = splitComponents[value].files 
	files = ["'%s'"%f for f in files]
	files = ', '.join(files)
	with open(new_job_path, 'w') as new_job_cfg :
	    with open(template_file_name) as f:
	        for line in f :
		    if line.find('REPLACE')  > 1 :
			actual_source_string = "fileNames = cms.untracked.vstring(%s),\n"%files 
        	        line = actual_source_string
		    new_job_cfg.write(line)
 
    def PrepareJobUserTemplate(self, jobDir, value ):
       '''Prepare one job. This function is called by the base class.'''
#       print value
#       print splitComponents[value]

       #prepare the batch script
       scriptFileName = jobDir+'/batchScript.sh'
       scriptFile = open(scriptFileName,'w')
       scriptFile.write( batchScriptCERN( value ) )
       scriptFile.close()
       os.system('chmod +x %s' % scriptFileName)
       
       print '\t',splitComponents[value].pyFragments

       variables = splitComponents[value].variables
       pyFragments = splitComponents[value].pyFragments
       setup = splitComponents[value].setup
       
       variables['IsMC'] = True
       if 'PD' in variables and not variables['PD'] == '': variables['IsMC'] = False
       #else: variables['PD'] = ""

       if batchManager.options_.verbose:
           print 'value ',value
           print 'scv ',splitComponents[value].name
       
       variables['SAMPLENAME'] = splitComponents[value].samplename
       variables['XSEC'] = splitComponents[value].xsec 
       #variables = {'IsMC':IsMC, 'PD':PD, 'MCFILTER':MCFILTER, 'SUPERMELA_MASS':SUPERMELA_MASS, 'SAMPLENAME':SAMPLENAME, 'XSEC':XSEC, 'SKIM_REQUIRED':SKIM_REQUIRED}

       print "\tParameters: ", variables

       execfile(cfgFileName,variables)
       
       process = variables.get('process') 
       process.source = splitComponents[value].source
       process.source.fileNames = cms.untracked.vstring('REPLACE_WITH_FILES') # splitComponents[value].files

       for fragment in pyFragments:
           execfile('pyFragments/{0:s}'.format(fragment),variables)  

       # PDF step 1 case: create also a snippet to be used later in step 2 phase
       if splitComponents[value].pdfstep == 1:
           cfgSnippetPDFStep2 = open(jobDir+'/inputForPDFstep2.py','w')
           cfgSnippetPDFStep2.write('process.source.fileNames = ["file:{0:s}/{1:s}"]\n'.format(self.outputDir_+'/AAAOK'+jobDir.replace(self.outputDir_,''), process.weightout.fileName.value()))
           cfgSnippetPDFStep2.write('process.source.secondaryFileNames = [')
           for item in splitComponents[value].files: cfgSnippetPDFStep2.write("'%s',\n" % item)
           cfgSnippetPDFStep2.write(']')
           cfgSnippetPDFStep2.write( '\n' )
           cfgSnippetPDFStep2.close()


       template_name = variables['SAMPLENAME'] + 'run_template_cfg.py'
       print 'Save template as=', '%s/%s'%(self.outputDir_, template_name)
       print 'jobDir=',jobDir
       cfgFile = open('%s/%s'%(self.outputDir_, template_name),'w')
       cfgFile.write( process.dumpPython() )
       cfgFile.write( '\n' )

       del process

       if splitComponents[value].pdfstep == 2:
           cfgSnippetPDFStep2 = open(self.inputPDFDir_+'/inputForPDFstep2.py','r')
           shutil.copyfileobj(cfgSnippetPDFStep2,cfgFile)
           cfgSnippetPDFStep2.close()
           

       

        
       cfgFile.close()


class Component(object):

    def __init__(self, name, prefix, dataset, pattern, splitFactor, variables, pyFragments, xsec, BR, setup, pdfstep):
        self.name = name
        self.samplename = name
        print "checking "+self.name
        self.source = datasetToSource( prefix, dataset, pattern)
        self.files = self.source.fileNames
        print len(self.files)
        if len(self.files)==0 :
            sys.exit("ERROR: no input files found")
        self.splitFactor = int(splitFactor)
        self.variables = variables
        self.pyFragments = pyFragments
        self.xsec = float(xsec)*float(BR)
        self.setup = setup
        self.pdfstep = int(pdfstep)
        if self.pdfstep <0 or self.pdfstep>2:
            print "Unknown PDF step", pdfstep
            sys.exit(1)
        
        
      
if __name__ == '__main__':
    batchManager = MyBatchManager()
    
    cfgFileName = batchManager.options_.cfgFileName #"analyzer_2015.py" # This is the python job config. FIXME make it configurable.
    sampleCSV  = batchManager.args_[0]            # This is the csv file with samples to be analyzed./    

    handle = open(cfgFileName, 'r')
    cfo = imp.load_source("pycfg", cfgFileName, handle)
    setup = cfo.LEPTON_SETUP

    components = []
    sampleDB = readSampleDB(sampleCSV)
    for sample, settings in sampleDB.iteritems():
        if settings['execute']:
            pdfstep = batchManager.options_.PDFstep
            components.append(Component(sample, settings['prefix'], settings['dataset'], settings['pattern'], settings['splitLevel'], settings['::variables'],settings['::pyFragments'],settings['crossSection'], settings['BR'], setup, pdfstep)) #FIXME-RB not bool(settings['pdf']))) #settings['pdf'] used here as full sel, without cuts.
    
    handle.close()


    splitComponents = split( components )

    listOfValues = range(0, len(splitComponents))
    listOfNames = [comp.name for comp in splitComponents]

    batchManager.PrepareJobs( listOfValues, listOfNames )

    waitingTime = 0.05
#FIXME to be implemented; should check batchManager.options_.negate; can simply call resubmit.csh
#batchManager.SubmitJobs( waitingTime )

