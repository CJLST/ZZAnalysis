#!/usr/bin/env python3
#
# Example of steering file for secondary tree processing in order to add variables in friend trees.
# For use in batch processing, use:
# batch_Condor.py samples_secondary.csv -i secondaryProd.py


from ZZAnalysis.NanoAnalysis.tools import setConf, getConf

NANOVERSION = getConf("NANOVERSION", 12)
LEPTON_SETUP = getConf("LEPTON_SETUP", 2022)
fileNames = getConf("fileNames", ["ZZ4lAnalysis.root",])

# In batch processing, config fragments can be specified in the pyFragment field;
# in interactive production they have to be called explicitly.
#import prod.pyFragments.DefaultProbs 

melaSettings = getConf("probabilities", None)
from ZZAnalysis.NanoAnalysis.initializeMELA import * 
mela = initializeMELA(True, LEPTON_SETUP)
from ZZAnalysis.NanoAnalysis.RecoProbFiller import *

sequence = [RecoProbFiller(mela, NANOVERSION, melaSettings)]

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
p = PostProcessor(".", fileNames,
                  prefetch=False, longTermCache=False,
                  cut=None,
                  branchsel=None,
                  outputbranchsel=None, # select branches to be written out
                  jsonInput=None,
                  modules=sequence,
                  noOut=False,
                  #haddFileName="ZZ4lAnalysis_ext.root",
                  maxEntries=0,
                  firstEntry=0,
                  friend=True,
#                  postfix="_ext",
                  provenance = False
                  )

p.run()

# Add a soft link to input file, so that the original and the friend tree can be found in the same place. 
import os, re
for file in fileNames :
    file = re.sub('^file:','', file)
    os.symlink(file, os.path.basename(file))
