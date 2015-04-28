#!/usr/bin/env python
"""
Some tools stolen from https://github.com/CERN-PH-CMG/cmg-cmssw/blob/CMGTools-from-CMSSW_7_4_0/CMGTools/Production/python/eostools.py
and from CMGTools.Production.datasetToSource
"""

import sys
import os
import FWCore.ParameterSet.Config as cms

def setCAFPath():
    """Hack to get the CAF scripts on the PYTHONPATH"""
    caf = '/afs/cern.ch/cms/caf/python'

    if caf not in sys.path:
        sys.path.append(caf)
setCAFPath()
import cmsIO


def runXRDCommand(path, cmd, *args):
    """Run an xrd command.
    """
    
    lfn = eosToLFN(path)
    #print "lfn:", lfn, cmd
    tokens = cmsIO.splitPFN(lfnToPFN(lfn))
    
    command = ['xrd', tokens[1], cmd, tokens[2]]
    command.extend(args)
    runner = cmsIO.cmsFileManip()
    # print ' '.join(command)
    return runner.runCommand(command)



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



def listFiles(path, rec = False, full_info = False):
    """Provides a list of the specified directory
    """
    # -- listing on the local filesystem --
    if os.path.isdir( path ):
        if not rec:
            # not recursive
            return [ '/'.join([path,file]) for file in os.listdir( path )]
        else:
            # recursive, directories are put in the list first,
            # followed by the list of all files in the directory tree
            result = []
            allFiles = []
            for root,dirs,files in os.walk(path):
                result.extend( [ '/'.join([root,dir]) for dir in dirs] )
                allFiles.extend( [ '/'.join([root,file]) for file in files] )
            result.extend(allFiles)
            return result
    # -- listing on EOS --
    cmd = 'dirlist'
    if rec:
        cmd = 'dirlistrec'
    files, _, _ = runXRDCommand(path, cmd)
    result = []
    for line in files.split('\n'):
        tokens = [t for t in line.split() if t]
        if tokens:
            #convert to an LFN
            # result.append(tuple(tokens))
            #COLIN need same interface for eos and local fs
            if full_info:
                result.append( tokens)
            else:
                result.append( tokens[4] )
    return result


def datasetToSource( user, dataset, pattern='.*root', readCache=False):

#    print user, dataset, pattern
#    data = createDataset(user, dataset, pattern, readCache)
    data=listFiles(dataset)

#    print data
    
    source = cms.Source(
	"PoolSource",
	noEventSort = cms.untracked.bool(True),
	duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
	fileNames = cms.untracked.vstring(),
        secondaryFileNames = cms.untracked.vstring(),
        )
    
    source.fileNames.extend( data )

    return source
