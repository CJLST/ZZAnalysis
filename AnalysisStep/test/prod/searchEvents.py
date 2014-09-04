#! /usr/bin/env python


##################################
## R. Bellan (UNITO) - Sep 2014 ##
##################################

import sys, os, commands, math, re, string

import ROOT
from ROOT import gSystem, gPad
gSystem.Load("libFWCoreFWLite")
from ROOT import AutoLibraryLoader
AutoLibraryLoader.enable()
gSystem.Load("libDataFormatsFWLite")

from DataFormats.FWLite import Events, Handle, Lumis

def Green(st):
    return '\033[0;32m'+str(st)+'\033[00m'


inputdir = '/data3/2014/HZZ_out/140519/PRODFSR_8TeV/'
sample   = 'DYJetsToLLTuneZ2M50-NoB'
wantedEvents = [70756397,55091083,26404299]






def findEvents(wantedEvents, sample, inputdir):
    failure, output = commands.getstatusoutput('ls {0:s}/{1:s}_Chunk*/*.root | grep Chunk'.format(inputdir,sample))
    files = output.split()
    checkTree('ZZ4muTree/candTree', wantedEvents, files)
    checkTree('ZZ4eTree/candTree', wantedEvents, files)
    checkTree('ZZ2e2muTree/candTree', wantedEvents, files)

def checkTree(treename,wantedEvents, files):

    tree = ROOT.TChain(treename)
    for f in files:
        tree.Add(f)
    entries = tree.GetEntries()


    cfgfiles = []

    for jentry in xrange(entries):           
        # get the next tree in the chain and verify
        ientry = tree.LoadTree(jentry)
        if ientry < 0: break
    # copy next entry into memory and verify
        nb = tree.GetEntry(jentry)
        if nb<=0: continue
        # use the values directly from the tree
        #    print mychain.MHT, weight
        if tree.EventNumber in wantedEvents:
            print Green("Found {0:.0f} in local tree {1:s} of file {2:s}".format(tree.EventNumber, treename, tree.GetCurrentFile().GetName()))
            cfg = tree.GetCurrentFile().GetName().replace('ZZ4lAnalysis.root','run_cfg.py')
            print cfg
            if cfg not in cfgfiles:
                cfgfiles.append(cfg)
        
    for cfg in cfgfiles:
        print cfg
        infoFile                  = open(cfg, "r")
        edmRootFiles = []
        for line in infoFile:
            line                    = line.strip()
            if not line:              continue
            #if not line.startswith("'") '.root' not in line: continue
            if '.root' not in line: continue
            if not line.startswith("'") and 'fileNames' not in line: continue
            print line.split("'")[1::2][0]
            edmRootFile = 'root://eoscms//eos/cms'+line.split("'")[1::2][0]

            events = Events(edmRootFile)
            for event in events:
                if event.eventAuxiliary().event() in wantedEvents:
                    print Green("Found {0:.0f} in edm file {1:s}".format(event.eventAuxiliary().event(), edmRootFile)) 


findEvents(wantedEvents, sample, inputdir)
