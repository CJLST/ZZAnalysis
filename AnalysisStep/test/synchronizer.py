#!/usr/bin/env python

import ROOT
import math
import optparse
import os, sys
from syncUtils import *


# define function for parsing options
def parseOptions():

    usage = ('usage: %prog [options] datasetList\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    parser.add_option('-i', '--input', dest='inFile', type='string', default="ZZ4lAnalysis.root",    help='input file')
    parser.add_option('-n', '--noOutput', dest='noOutput', action='store_true', default=False, help='do not write sync file in output')
    parser.add_option('-o', '--output', dest='outFile', type='string', default="eventlist.txt",    help='output sync file')
    parser.add_option('-f', '--finalState', dest='finalState', type='string', default="all",    help='final states: all, 4e, 4mu, 2e2mu')


    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    if not "." in opt.outFile:
        print "Please use an extension for the output file (e.g. \".txt\")"
        sys.exit()
    
            

def loop():
    

    inFileName = opt.inFile
    outFileName = opt.outFile
    finalState = opt.finalState
        
    print "Processing file: ",inFileName,"..."

    cands = []
    totCounter = 0
    chanCounter = {}
    for aChan in ["4mu","4e","2e2mu"]:

        chanCounter[aChan] = 0

        if finalState!="all" and aChan!=finalState: continue
        
        tree = ROOT.TChain("ZZ"+aChan+"Tree/candTree")
        tree.Add(inFileName)
        tree.SetBranchStatus("*",0)

        # Variables we are interested in for the sync
        tree.SetBranchStatus("iBC",1)
        tree.SetBranchStatus("ZZsel",1)        
        tree.SetBranchStatus("RunNumber",1)
        tree.SetBranchStatus("LumiNumber",1)
        tree.SetBranchStatus("EventNumber",1)            
        tree.SetBranchStatus("ZZMass",1)
        tree.SetBranchStatus("Z1Mass",1)
        tree.SetBranchStatus("Z2Mass",1)
        tree.SetBranchStatus("ZZMassErr",1)
        tree.SetBranchStatus("ZZMassErrCorr",1)
        tree.SetBranchStatus("p0plus_VAJHU",1)
        tree.SetBranchStatus("p0minus_VAJHU",1)
        tree.SetBranchStatus("p0hplus_VAJHU",1)
        tree.SetBranchStatus("p1plus_VAJHU",1)
        tree.SetBranchStatus("p1_VAJHU",1)
        tree.SetBranchStatus("p2_VAJHU",1)
        tree.SetBranchStatus("p2qqb_VAJHU",1)            
        tree.SetBranchStatus("bkg_VAMCFM",1)
        tree.SetBranchStatus("ZZPt",1)
        tree.SetBranchStatus("JetPt",1)
        tree.SetBranchStatus("DiJetMass",1)
        tree.SetBranchStatus("DiJetDEta",1)            

        iEntry=0
        while tree.GetEntry(iEntry):

            # print "   Inspecting entry n. ",iEntry,"..."
            iEntry+=1
            ZZsel       = tree.ZZsel
            iBC         = tree.iBC
            run         = tree.RunNumber
            lumi        = tree.LumiNumber
            event       = tree.EventNumber

            theEvent = Event(iBC,run,lumi,event)

            if iBC>=0 and ZZsel[iBC]>=90:
                totCounter += 1
                chanCounter[aChan] += 1
                mass4l        = tree.ZZMass[iBC]
                mZ1           = tree.Z1Mass[iBC]
                mZ2           = tree.Z2Mass[iBC]
                massErrRaw    = tree.ZZMassErr[iBC]
                massErrCorr   = tree.ZZMassErrCorr[iBC]
                p0plus_VAJHU  = tree.p0plus_VAJHU[iBC]
                p0minus_VAJHU = tree.p0minus_VAJHU[iBC]
                p0hplus_VAJHU = tree.p0hplus_VAJHU[iBC]
                p1plus_VAJHU  = tree.p1plus_VAJHU[iBC] 
                p1_VAJHU      = tree.p1_VAJHU[iBC]     
                p2_VAJHU      = tree.p2_VAJHU[iBC]     
                p2qqb_VAJHU   = tree.p2qqb_VAJHU[iBC]              
                bkg_VAMCFM    = tree.bkg_VAMCFM[iBC]                    
                pt4l          = tree.ZZPt[iBC]
                jetpt         = tree.JetPt
                njets30       = 0
                jet1pt        = 0.
                jet2pt        = 0.
                mjj           = tree.DiJetMass
                detajj        = tree.DiJetDEta
                fishjj        = -1.                

                jets30 = []
                
                for i in range(len(jetpt)):                    
                    if jetpt[i]>30.:
                        jets30.append(jetpt[i])                    
                    
                theKDs = KDs(p0plus_VAJHU,p0minus_VAJHU,p0hplus_VAJHU,p1plus_VAJHU,p1_VAJHU,p2_VAJHU,p2qqb_VAJHU,bkg_VAMCFM)
                theCand = Candidate(theEvent,mass4l,mZ1,mZ2,massErrRaw,massErrCorr,pt4l,jets30,mjj,detajj,theKDs)
                cands.append(theCand)



    # Sort candidates on a event number basis
    sortedCands = sorted(cands, key=lambda cand: float(cand.eventInfo.event)) 

    if not opt.noOutput:
        # Print in sync format
        outFile      = open(outFileName.replace(".","_"+opt.finalState+"."),"w")
        line = ""
        
        for aCand in sortedCands:
            line += aCand.eventInfo.printOut()
            line += ":"
            line += aCand.printOut()
            line += "\n"                

        
        outFile.write(line)
        outFile.close()
        
        print "Output written in file: ",outFileName,"\n"

    counterStr = str(totCounter) + "/" + str(chanCounter["4e"]) + "/" + str(chanCounter["4mu"]) + "/" + str(chanCounter["2e2mu"])
    print "\n## Selected events all/4e/4mu/2e2mu : "+counterStr
    

        

if __name__ == "__main__":

    parseOptions()
    loop()
