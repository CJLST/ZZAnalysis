#!/bin/env python
###
# Read a csv production file for MINIAOD and output the corresponding one for NANO.
###

from __future__ import print_function
import sys
# import imp
# import copy
# import os
# import shutil
# import pickle
# import math
# import pprint
# import subprocess
#from datetime import date
from optparse import OptionParser
from ZZAnalysis.AnalysisStep.eostools import *
from ZZAnalysis.AnalysisStep.readSampleInfo import *



if __name__ == '__main__' :
    parser = OptionParser()
    (options, args) = parser.parse_args()
    
    csvfile = args[0] 
#    sampleDB = readSampleDB(csvfile)
#    for sample, settings in sampleDB.iteritems() :
#        print(sample, settings)


    infoFile                  = open(csvfile, "r")
#    defaults                  = {}
    header                    = None
    datasetIdx = -1
    prefixIdx = -1
    patternIdx = -1
    
    for line in infoFile:
        line                    = line.strip()
        if not line:              continue
        if line.startswith("#"):
            print(line)
        else :
            data = line.split(",")
            if header:            
                if len(data) != len(header):
                    raise ValueError, "Inconsistent number of columns in data '" + line + "', expected header = " + str(header)

                dataset = data[datasetIdx]
                instance = 'prod/global'
                cmd = '"child dataset='+ dataset +' instance=%s"'%instance
                command = ['/cvmfs/cms.cern.ch/common/dasgoclient' , '--limit=0', '--query', cmd]
                run_command = ' '.join(command)
                runner = cmsFileManip()
                cmdout, _, _ = runner.runCommand(run_command)
                result = []
                for line in cmdout.split('\n'):
                    if line != "" : result.append(line)
                if len(result) != 1 :
                    print("ERROR: das output", cmdout, "contains >1 child") # FIXME extend script to ask for choice
                    exit(1)                    
                data[datasetIdx] = result[0]

                # Now find if nanoAODs are accessible at CERN
                cmd = '"site dataset='+ data[datasetIdx] +' instance=%s"'%instance
                command = ['/cvmfs/cms.cern.ch/common/dasgoclient' , '--limit=0', '--query', cmd]
                run_command = ' '.join(command)
                cmdout, _, _ = runner.runCommand(run_command)
                for line in cmdout.split('\n'):
                    if "CERN" in line :
                        print("NOTE:", data[0], "available at", line)

#                data[prefixIdx] = "dbs"
                data[patternIdx] = "root://cms-xrd-global.cern.ch/"
                
                print(",".join(data))

            # Read header information (first line in file)
            else:
                header = []
                for idatum,datum in enumerate(data):
                    if not datum: break
                    if "=" in datum:
                        (datum, default)  = map(str.strip, datum.split("="))
                    if datum=="dataset" : datasetIdx = idatum
                    if datum=="prefix"  : prefixIdx = idatum
                    if datum=="pattern" : patternIdx = idatum
                    header.append(datum)
                
                print(line)
