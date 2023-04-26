#!/bin/env python
# Add data files for different PDs, after archiving the Chunk folders.
# Assumes that haddChunks has already been run to merge chunks.
#
# Usage:
# haddData.py [folder]
# [folder] = location of  (default: "AAAAOK")
# Merged data are placed in a folder Data[year]/
# in the CURRENT directory (not uder [folder]!)

#
from __future__ import print_function
import sys
import os
import glob

if len(sys.argv) > 1 :
    jobpath = str(sys.argv[1])
else: 
    jobpath = "AAAOK"

# Check if the specified file is a nanoAOD file.
def checkNano(file):
    from ROOT import TFile
    rf = TFile.Open(file)
    keys = rf.GetListOfKeys()
    if (keys.FindObject("Events") != None and keys.FindObject("Runs") != None and  keys.FindObject("LuminosityBlocks") != None) :
        isNano = True
    else :
        isNano= False
    rf.Close()
    return isNano


if (os.path.exists(jobpath) == False ):
    print("Folder", jobpath, "does not exist; please specify a valid one.\n")
    sys.exit(1)


# Archive Chunk folders, if any
chunks = glob.glob(jobpath+"/*_Chunk*")
if len(chunks) :
    print("Archiving Chunks...")
    os.system('cd '+jobpath+'; mkdir -p Chunks; mv *_Chunk* Chunks')


dataYears = ["2016","2017","2018","2022","2023"]

mergedYears = []
isNano = False

for year in dataYears :
    rootfiles=glob.glob(jobpath+'/*'+year+"*/ZZ4lAnalysis.root")
    if len(rootfiles) :
        dest = "Data"+year
        os.system('mkdir -p '+dest)
        mergedYears.append(dest)
        if len(rootfiles) == 1 : # only one file, move it
            haddCmd = 'cp ' + rootfiles[0] + " " + dest
        else : # hadd files
            if checkNano(rootfiles[0]) :
                isNano = True
                haddCmd = 'haddnano.py ' + dest + '/ZZ4lAnalysis.root ' + ' '.join(rootfiles)
            else :
                haddCmd = 'hadd -ff '  + dest + '/ZZ4lAnalysis.root ' + ' '.join(rootfiles)
        print ("hadding files:", haddCmd)
        os.system(haddCmd)

if len(mergedYears) :
    # ... still have merge different years into a single file 
    dest = "AllData"
    os.system('mkdir -p '+dest)    
    if isNano :
        haddCmd = 'haddnano.py ' + dest + '/ZZ4lAnalysis.root Data*/ZZ4lAnalysis.root'
    else:
        haddCmd = 'hadd -ff '  + dest + '/ZZ4lAnalysis.root Data*/ZZ4lAnalysis.root'
    print ("hadding all years:", haddCmd)
    os.system(haddCmd)
    
