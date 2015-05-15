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

    hfile = ROOT.TFile(inFileName)
    hCounters = hfile.Get("ZZTree/Counters")

    for aChan in ["4mu","4e","2e2mu"]:

        chanCounter[aChan] = 0

        if finalState!="all" and aChan!=finalState: continue
        
#        tree = ROOT.TChain("ZZ"+aChan+"Tree/candTree")
        tree = ROOT.TChain("ZZTree/candTree")
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
        tree.SetBranchStatus("Z1Flav",1)
        tree.SetBranchStatus("Z2Flav",1)
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
        tree.SetBranchStatus("p0plus_m4l",1)
        tree.SetBranchStatus("bkg_m4l",1)
        tree.SetBranchStatus("Dgg10_VAMCFM",1)
        tree.SetBranchStatus("pvbf_VAJHU_new",1)
        tree.SetBranchStatus("phjj_VAJHU_new",1)        
        tree.SetBranchStatus("ZZPt",1)
        tree.SetBranchStatus("nExtraLep",1)
        tree.SetBranchStatus("nCleanedJetsPt30BTagged",1)
        tree.SetBranchStatus("JetPt",1)
        tree.SetBranchStatus("JetEta",1)
        tree.SetBranchStatus("JetPhi",1)
        tree.SetBranchStatus("JetMass",1)
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

            if iBC>=0 and ZZsel>=90 :
                ZZflav      = tree.Z1Flav*tree.Z2Flav;
                if  (aChan=="4e" and ZZflav!=14641) or (aChan=="4mu" and ZZflav!=28561) or (aChan=="2e2mu" and ZZflav!=20449) : continue

                totCounter += 1
                chanCounter[aChan] += 1
                mass4l        = tree.ZZMass
                mZ1           = tree.Z1Mass
                mZ2           = tree.Z2Mass
                massErrRaw    = tree.ZZMassErr
                massErrCorr   = tree.ZZMassErrCorr
                p0plus_VAJHU  = tree.p0plus_VAJHU
                p0minus_VAJHU = tree.p0minus_VAJHU
                p0hplus_VAJHU = tree.p0hplus_VAJHU
                p1plus_VAJHU  = tree.p1plus_VAJHU 
                p1_VAJHU      = tree.p1_VAJHU     
                p2_VAJHU      = tree.p2_VAJHU     
                p2qqb_VAJHU   = tree.p2qqb_VAJHU              
                bkg_VAMCFM    = tree.bkg_VAMCFM
                p0plus_m4l    = tree.p0plus_m4l
                bkg_m4l       = tree.bkg_m4l
                Dgg10_VAMCFM  = tree.Dgg10_VAMCFM
                pvbf_VAJHU    = tree.pvbf_VAJHU_new
                phjj_VAJHU    = tree.phjj_VAJHU_new
                pt4l          = tree.ZZPt
                nExtraLep     = tree.nExtraLep
                jetpt         = tree.JetPt
                jeteta        = tree.JetEta
                jetphi        = tree.JetPhi
                jetmass       = tree.JetMass
                njets30Btag   = tree.nCleanedJetsPt30BTagged
                mjj           = tree.DiJetMass
                detajj        = tree.DiJetDEta

                jets30pt = []
                jets30eta = []
                jets30phi = []
                jets30mass = []

                for i in range(len(jetpt)): 
                    if jetpt[i]>30.:
                        jets30pt.append(jetpt[i])
                        jets30eta.append(jeteta[i])
                        jets30phi.append(jetphi[i])
                        jets30mass.append(jetmass[i])
                    
                theKDs = KDs(p0plus_VAJHU,p0minus_VAJHU,p0hplus_VAJHU,p1plus_VAJHU,p1_VAJHU,p2_VAJHU,p2qqb_VAJHU,bkg_VAMCFM,p0plus_m4l,bkg_m4l,Dgg10_VAMCFM,pvbf_VAJHU,phjj_VAJHU)
                theCand = Candidate(theEvent,mass4l,mZ1,mZ2,massErrRaw,massErrCorr,pt4l,nExtraLep,jets30pt,jets30eta,jets30phi,jets30mass,njets30Btag,mjj,detajj,theKDs)
                cands.append(theCand)



    # Sort candidates on a event number basis
    sortedCands = sorted(cands, key=lambda cand: float(cand.eventInfo.event)) 

    if not opt.noOutput:
        # Print in sync format
        outFileName = outFileName.replace(".","_"+opt.finalState+".")
        outFile = open(outFileName,"w")
        line = ""
        
        for aCand in sortedCands:
            line += aCand.eventInfo.printOut()
            line += ":"
            line += aCand.printOut()
            line += "\n"                

        
        outFile.write(line)
        outFile.close()
        
        print "Output written in file: ",outFileName,"\n"

        print "## Total/4e/4mu/2e2mu/2l2tau : ", int(hCounters.GetBinContent(1)),  "/",  int(hCounters.GetBinContent(3)),  "/",  int(hCounters.GetBinContent(2)),  "/",  int(hCounters.GetBinContent(4)), "/", int (hCounters.GetBinContent(5))

    counterStr = str(totCounter) + "/" + str(chanCounter["4e"]) + "/" + str(chanCounter["4mu"]) + "/" + str(chanCounter["2e2mu"])
    print "\n## Selected events all/4e/4mu/2e2mu : "+counterStr
    

        

if __name__ == "__main__":

    parseOptions()
    loop()
