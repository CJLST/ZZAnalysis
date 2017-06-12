#! /usr/bin/env python


##################################
## R. Bellan (UNITO) - Sep 2014 ##
##################################
##
## N.A. 2017 - extended to produce a csv file for batch.py in order to rerun on selected events of a given tree.

inputdir = '/data3/Higgs/170222/AAAOK/' # path of the production Chunks
samples   = ['MuonEG2016B',    'MuonEG2016C',    'MuonEG2016D',    'MuonEG2016E',    'MuonEG2016F',    'MuonEG2016G',    'MuonEG2016Hv2',    'MuonEG2016Hv3',    
             'DoubleEG2016B',  'DoubleEG2016C',  'DoubleEG2016D',  'DoubleEG2016E',  'DoubleEG2016F',  'DoubleEG2016G',  'DoubleEG2016Hv2',  'DoubleEG2016Hv3',  
             'DoubleMu2016B',  'DoubleMu2016C',  'DoubleMu2016D',  'DoubleMu2016E',  'DoubleMu2016F',  'DoubleMu2016G',  'DoubleMu2016Hv2',  'DoubleMu2016Hv3',  
             'SingleEle2016B', 'SingleEle2016C', 'SingleEle2016D', 'SingleEle2016E', 'SingleEle2016F', 'SingleEle2016G', 'SingleEle2016Hv2', 'SingleEle2016Hv3', 
             'SingleMuon2016B','SingleMuon2016C','SingleMuon2016D','SingleMuon2016E','SingleMuon2016F','SingleMuon2016G','SingleMuon2016Hv2','SingleMuon2016Hv3',
             ] #which samples to look at

getAllEvents = True        # Take all events from the trees (passing the criteria sepcifed below)
wantedEvents = [312810956] # othewise, search for event n. in this list. (Note that only event N is used, lumi and run are ignored...)

prepareCSV = True          # Create a searchEvents.csv file in order to reprocess events 
verifyEdmFile = False      # open file to check that event is present (useful for jobs with many input files/job)
debug = False

unblind = True            # Unblinded (affects data only)
range = 'low'              # full,low,high mass region (only for blinded data)


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



def findEvents(wantedEvents, sample, inputdir):
    failure, output = commands.getstatusoutput('ls {0:s}/{1:s}_Chunk*/*.root | grep Chunk'.format(inputdir,sample))
    files = output.split()
    return checkTree(sample, 'ZZTree/candTree', wantedEvents, files)

def checkTree(sample, treename,wantedEvents, files):

    tree = ROOT.TChain(treename)
    for f in files:
        tree.Add(f)
    entries = tree.GetEntries()

    cfgfiles = []
    eventIds = []
    LFNs = []

    for jentry in xrange(entries):           
        # get the next tree in the chain and verify
        ientry = tree.LoadTree(jentry)
        if ientry < 0: break

        # copy next entry into memory and verify
        nb = tree.GetEntry(jentry)
        if nb<=0: continue

        if getAllEvents :
            # Select events to be considered            
            if tree.ZZsel>=90 :
                mass4l = tree.ZZMass
                if tree.RunNumber>100 and unblind==False :
                ## blind Higgs peak in data
                    isLowUnblind  = (mass4l>=70 and mass4l<=110)
                    isHighUnblind = (mass4l>=150 and mass4l<=300)
                    if range=='low' and not isLowUnblind : continue
                    elif range=='high' and not isHighUnblind : continue
                    elif range=='full' and not (isLowUnblind or isHighUnblind) : continue

                eventId = "{0:d}:{1:d}:{2:d}".format(tree.RunNumber, tree.LumiNumber, tree.EventNumber)
                eventIds.append(eventId)
                print Green("Found {0:s} in local tree {1:s} of file {2:s}".format(eventId, treename, tree.GetCurrentFile().GetName()))

                cfg = tree.GetCurrentFile().GetName().replace('ZZ4lAnalysis.root','run_cfg.py')
                if debug : print cfg
                if cfg not in cfgfiles:
                    cfgfiles.append(cfg)                    

        else : 
            if tree.EventNumber in wantedEvents:
                cfg = tree.GetCurrentFile().GetName().replace('ZZ4lAnalysis.root','run_cfg.py')
                if debug : print cfg
                if cfg not in cfgfiles:
                    cfgfiles.append(cfg)
        
    for cfg in cfgfiles:
        if debug : print cfg
        infoFile                  = open(cfg, "r")
        edmRootFiles = []
        for line in infoFile:
            line                    = line.strip()
            if not line:              continue
            if '.root' not in line: continue
            if not line.startswith("'") and 'fileNames' not in line: continue
            LFN = line.split("'")[1::2]

            if verifyEdmFile : 
                edmRootFile = 'root://eoscms//eos/cms'+LFN[0] #FIXME: must loop over all files...
                events = Events(edmRootFile)
                for event in events:
                    if event.eventAuxiliary().event() in wantedEvents:
                        print Green("Found {0:.0f} in edm file {1:s}".format(event.eventAuxiliary().event(), edmRootFile))
                        LFNs.append(LFN)
            else:                
                print Green(LFN)
                LFNs+=LFN
            #   print Green("candidate{0:d}:{1:d}:{2:d} in edm file {1:s}".format(event.eventAuxiliary().run(), event.eventAuxiliary().lumi(), event.eventAuxiliary().event(), edmRootFile))
    if len(eventIds) == 0 :
        return 0
    f = open('pyFragments/searchEvents_'+sample+'.py', 'w')
    f.write('process.source.fileNames = cms.untracked.vstring()\n')
    f.write('process.source.fileNames.extend(['+',\n'.join('\'%s\'' % val for val in LFNs)+ '])\n')
    f.write('process.source.eventsToProcess = cms.untracked.VEventRange('+ ',\n'.join('\'%s\'' % val for val in eventIds)+')\n')
    f.close()
    f = open('searchEvents.csv', 'a')
    PD = ''
    if   'MuonEG' in sample: PD = 'MuEG'
    elif 'DoubleEG' in sample: PD = 'DoubleEle'
    elif 'DoubleMu' in sample: PD = 'DoubleMuon'
    elif 'SingleEle' in sample: PD = 'SingleElectron'
    elif 'SingleMu' in sample: PD = 'SingleMuon'

    addVariables=';SKIP_EMPTY_EVENTS=False' # Add SKIP_EMPTY_EVENTS to avoid problems with changes in trigger requirements (events may have to be picked from different PDs...)
    addVariables+=';KINREFIT=True'
    if 'Hv2' in sample or 'Hv3' in sample : addVariables += ';DATA_TAG=PromptReco' # FIXME: needed at the time of Moriond17
    f.write(sample+',,-1,,,,source.fileNames,,20,PD='+PD+';PROCESS_CR=True;ADDLOOSEELE=False'+addVariables+',json_2016.py;RecoProbabilities.py;searchEvents_'+sample+'.py, # '+str(len(eventIds))+' events \n')
    f.close()
    return len(eventIds)


if prepareCSV :
    f = open('searchEvents.csv', 'w')
    f.write('identifier,process,crossSection=-1,BR=1,execute=True,dataset,prefix,pattern,splitLevel,::variables,::pyFragments,comment\n')
    f.close()


total = 0
for sample in samples :
    found = findEvents(wantedEvents, sample, inputdir)
    print 'Sample:', sample, 'Found: ', found, 'events'
    total += found
print 'Total: ' , total, 'events'
