#!/usr/bin/env python
"""
Some tools to handle accesing dataset information on dbs and files on eos
stolen from https://github.com/CERN-PH-CMG/cmg-cmssw/blob/CMGTools-from-CMSSW_7_4_0/CMGTools/Production/python/eostools.py
and from CMGTools.Production.datasetToSource
"""

import sys
import os
import re
import FWCore.ParameterSet.Config as cms
import subprocess

class cmsFileManip:
    """A class to interact with files/directories"""
    def runCommand( self, cmd, setShell=True ) :
        myCommand = subprocess.Popen( cmd,
                                      stdin=subprocess.PIPE,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE,shell=setShell )
        ( out, err ) = myCommand.communicate()
        if myCommand.returncode != 0:
            print >> sys.stderr, "Command (%s) failed with return code: %d" % ( cmd, myCommand.returncode )
            print >> sys.stderr, err

        return out,err,myCommand.returncode



def setCAFPath():
    """Hack to get the CAF scripts on the PYTHONPATH"""
    caf = '/afs/cern.ch/cms/caf/python'

    if caf not in sys.path:
        sys.path.append(caf)
setCAFPath()


def runXRDCommand(path, cmd, *args):
    """Run an xrd command.
    """

    lfn = eosToLFN(path)
    command = ['xrd', 'eoscms', cmd, lfn]
    command.extend(args)
    runner = cmsFileManip()
    # print ' '.join(command)
    return runner.runCommand(command, False) # fails when shell=True



def eosToLFN( path ):
    """Converts a EOS PFN to an LFN.
    Just strip out /eos/cms from path.
    If this string is not found, return path.
    ??? Shouldn't we raise an exception instead?"""
    return path.replace('root://eoscms.cern.ch/', '').replace('/eos/cms','')



def lfnToPFN( path, tfcProt = 'rfio'):
    """Converts an LFN to a PFN. For example:
    /store/cmst3/user/cbern/CMG/TauPlusX/Run2011A-03Oct2011-v1/AOD/V2/PAT_CMG_V2_4_0/H2TAUTAU_Nov21
    ->
    root://eoscms//eos/cms/store/cmst3/user/cbern/CMG/TauPlusX/Run2011A-03Oct2011-v1/AOD/V2/PAT_CMG_V2_4_0/H2TAUTAU_Nov21?svcClass=cmst3&stageHost=castorcms
    This function only checks path, and does not access the storage system.
    If the path is in /store/cmst3, it assumes that the CMST3 svcClass is to be used.
    Otherwise, is uses the default one.
    """

    if path.startswith("/store/"):
        path = path.replace("/store/","root://eoscms.cern.ch//eos/cms/store/")
    if path.startswith("/pnfs/psi.ch/cms/trivcat/"):
        path = path.replace("/pnfs/psi.ch/cms/trivcat/","root://t3se01.psi.ch//")
    #print "path to cmsFile():", path
    entity = cmsIO.cmsFile( path, tfcProt )
#    tokens = cmsIO.splitPFN(entity.pfn)
    pfn = '%s://%s//%s/' % (entity.protocol,entity.host,entity.path)

    pfn = entity.pfn
    if tfcProt == 'rfio' and \
        entity.path.startswith("/eos/cms/") and \
                str(entity.stat()).startswith("Error 3011: Unable to stat"):

            pfn.replace("/eos/cms","/castor/cern.ch/cms")
            pfn.replace("eoscms","castorcms")
    return pfn

def runDBS(dataset, instance = 'prod/global'):
    cmd = '"file dataset='+dataset +' instance=%s"'%instance

    command = ['/cvmfs/cms.cern.ch/common/dasgoclient' , '--limit=0', '--query', cmd]
    runner = cmsFileManip()
    run_command = ' '.join(command)
    return runner.runCommand(run_command)


def listFiles(sample, path, rec = False, full_info = False):
    """Provides a list of files with different methods according to path. Valid paths: 'list', 'dbs', 'dbs-USER', a local filesystem path, an eos path
    """
    result = []

    # -- list from a local file --
    if path=="list" :
        with open(sample) as f:
            result = f.readlines()
        return result

    # -- listing from dbs --
    elif path=="dbs" :
        files, _, _ =runDBS(sample)
        for line in files.split('\n'):
#            print line
#            result.append("root://cms-xrd-global.cern.ch//"+line)
            result.append(line)
        return result

    # -- listing from user dbs --
    elif path=="dbs-USER" :
	print 'Querying USER db'
	files, _, _ =runDBS(sample, 'prod/phys03')
        for line in files.split('\n'):
            result.append(line)
        return result


    # -- listing from path=local dir -- untested!
    elif path=="local" :
        if os.path.isdir(sample):
           if not rec:
             # not recursive
             return [ 'file:'+'/'.join([path,file]) for file in os.listdir(sample)]
           else:
             # recursive, directories are put in the list first,
             # followed by the list of all files in the directory tree
             allFiles = []
             for root,dirs,files in os.walk(sample):
                 result.extend( [ '/'.join([root,dir]) for dir in dirs] )
                 allFiles.extend( [ 'file:'+'/'.join([root,file]) for file in files] )
             result.extend(allFiles)
             return result

    # -- listing from EOS (path = eos path prefix, normally "/eos/cms/")
    else  :
        cmd = 'dirlist'
        if rec:
            cmd = 'dirlistrec'

        print path+sample, cmd
        files, _, _ = runXRDCommand(path+sample, cmd)

        for line in files.split('\n'):
            tokens = [t for t in line.split() if t]
            if tokens:
               #convert to an LFN
               if full_info:
                   result.append( tokens)
               else:
                   result.append( eosToLFN(tokens[4]) )
        return result

# Fill a cms.source with list of files. Arguments:
# Prefix: where the list of files is obtained from. Can be dbs, dbs-USER, eos, list, or a local path
# dataset: the dataset name if prefix is "dbs", the local path if "local", the file with a list if "list"
# fileprefix: path to be prepended to each file. Normally empty; other options:
# /eos/cms/ for files available on eos
# root://cms-xrd-global.cern.ch/ to force global redirector

def datasetToSource( prefix, dataset, fileprefix=''):

#    print user, dataset, fileprefix
    recursive=True #FIXME: this is needed for central production, but care is needed if other stuff is present in the EOS path

    data=listFiles(dataset, prefix, recursive)

    rootre = re.compile('.*root')
    files = [fileprefix+f for f in data if rootre.match(f)]

#    print files

    source = cms.Source(
	"PoolSource",
	noEventSort = cms.untracked.bool(True),
	duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
	fileNames = cms.untracked.vstring(),
        secondaryFileNames = cms.untracked.vstring(),
        )

    source.fileNames.extend( files )

    return source
