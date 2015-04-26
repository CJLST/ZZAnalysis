#!/usr/bin/env python

import ROOT
import math
import optparse
import os, sys
from syncUtils import *

# See interface/FinalStates.h for the convention
CRdict = {"CRZLLss":21, "CRZLLos_2P2F":22, "CRZLLos_3P1F":23}

# define function for parsing options
def parseOptions():

    usage = ('usage: %prog [options] datasetList\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    parser.add_option('-i', '--input', dest='inFile', type='string', default="ZZ4lAnalysis.root",    help='input file')
    parser.add_option('-o', '--output', dest='outFile', type='string', default="eventlist.txt",    help='output sync file')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()


def test_bit(mask, iBit):

    return (mask >> iBit) & 1
            

def loop():
    

    inFileName = opt.inFile
    outFileNameStr = opt.outFile
        
    print "Processing file: ",inFileName,"..."

    crCands = {}
    totCounter = 0
    crCounter = {}
    for aCR in CRdict.keys():
        crCounter[aCR] = 0
        crCands[aCR]   = []
    
    tree = ROOT.TChain("CRZLLTree/candTree")
    tree.Add(inFileName)
    tree.SetBranchStatus("*",0)

    # Variables we are interested in for the sync
    tree.SetBranchStatus("iBC",1)
    tree.SetBranchStatus("ZZsel",1)
    tree.SetBranchStatus("CRflag",1)    
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
    tree.SetBranchStatus("nExtraLep",1)
    tree.SetBranchStatus("nCleanedJetsPt30BTagged",1)
    tree.SetBranchStatus("JetPt",1)
    tree.SetBranchStatus("JetEta",1)
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
        CRflag      = tree.CRflag


        theEvent = Event(iBC,run,lumi,event)

        for iCand in range(len(CRflag)):

            for aCR in CRdict.keys():
            
                if test_bit(CRflag[iCand],CRdict[aCR]):
                    
                    crCounter[aCR] += 1                    
                    mass4l        = tree.ZZMass[iCand]
                    mZ1           = tree.Z1Mass[iCand]
                    mZ2           = tree.Z2Mass[iCand]
                    massErrRaw    = tree.ZZMassErr[iCand]
                    massErrCorr   = tree.ZZMassErrCorr[iCand]
                    p0plus_VAJHU  = tree.p0plus_VAJHU[iCand]
                    p0minus_VAJHU = tree.p0minus_VAJHU[iCand]
                    p0hplus_VAJHU = tree.p0hplus_VAJHU[iCand]
                    p1plus_VAJHU  = tree.p1plus_VAJHU[iCand] 
                    p1_VAJHU      = tree.p1_VAJHU[iCand]     
                    p2_VAJHU      = tree.p2_VAJHU[iCand]     
                    p2qqb_VAJHU   = tree.p2qqb_VAJHU[iCand]              
                    bkg_VAMCFM    = tree.bkg_VAMCFM[iCand]                    
                    pt4l          = tree.ZZPt[iCand]
                    nExtraLep     = tree.nExtraLep[iCand]
                    jetpt         = tree.JetPt
                    jeteta        = tree.JetEta
                    njets30Btag   = tree.nCleanedJetsPt30BTagged
                    mjj           = tree.DiJetMass
                    detajj        = tree.DiJetDEta

                    jets30pt = []
                    jets30eta = []

                    for i in range(len(jetpt)):                    
                        if jetpt[i]>30.:
                            jets30pt.append(jetpt[i])
                            jets30eta.append(jeteta[i])
                            
                    theKDs = KDs(p0plus_VAJHU,p0minus_VAJHU,p0hplus_VAJHU,p1plus_VAJHU,p1_VAJHU,p2_VAJHU,p2qqb_VAJHU,bkg_VAMCFM)
                    theCand = Candidate(theEvent,mass4l,mZ1,mZ2,massErrRaw,massErrCorr,pt4l,nExtraLep,jets30pt,jets30eta,njets30Btag,mjj,detajj,theKDs)
                    crCands[aCR].append(theCand)


    crSortedCands = crCands
    for aCR in CRdict.keys():
        # Sort candidates on a event number basis
        crSortedCands[aCR] = sorted(crCands[aCR], key=lambda cand: float(cand.eventInfo.event))

        # Print in sync format
        outFileName = outFileNameStr.replace(".txt","_"+aCR+".txt")
        outFile = open(outFileName,"w")        
        line = ""
        for aCand in crSortedCands[aCR]:
            line += aCand.eventInfo.printOut()
            line += ":"
            line += aCand.printOut()
            line += "\n"
            
        outFile.write(line)
        outFile.close()        

        
        print "\n## Events in CR "+aCR+": ",crCounter[aCR],"events"
        print "Output written in file: ",outFileName




        

if __name__ == "__main__":

    parseOptions()
    loop()
