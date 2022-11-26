#!/usr/bin/env python
##
# Print a list of selected candidates from nanoAOD trees.
##
from __future__ import print_function
import math
import argparse
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

parser = argparse.ArgumentParser(description='Print selected 4l events in the specified nanoAOD files')

parser.add_argument('file', help='filename to read', nargs='?', default="ggH125_fixedFSR_Skim.root")

args = parser.parse_args()

nanoFile = args.file


t = ROOT.TChain("Events")
t.Add(nanoFile)

iEntry = 0
while t.GetEntry(iEntry):
    iEntry+=1

    eventId='{}:{}:{}'.format(t.run,t.luminosityBlock,t.event)
    
    iZZ = t.bestCandIdx
    if iZZ < 0 : continue # no candidate passes selection

    print('{}:{:.2f}:{:.2f}:{:.2f}:{:.3f}:{:.3f}:{:.3f}'.format(eventId,
                                                                t.ZZCand_mass[iZZ],
                                                                t.ZZCand_Z1mass[iZZ],
                                                                t.ZZCand_Z2mass[iZZ],
                                                                t.ZZCand_KD[iZZ],
                                                                t.overallEventWeight,
                                                                t.ZZCand_dataMCWeight[iZZ],
                                                                ))
