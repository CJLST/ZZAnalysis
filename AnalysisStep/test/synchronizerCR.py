#!/usr/bin/env python

import ROOT
import math
import optparse
import os, sys
from syncUtils import *
from operator import attrgetter

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
    tree.SetBranchStatus("CRflag",1)
    tree.SetBranchStatus("RunNumber",1)
    tree.SetBranchStatus("LumiNumber",1)
    tree.SetBranchStatus("EventNumber",1)
    tree.SetBranchStatus("genHEPMCweight",1)
    tree.SetBranchStatus("PUWeight",1)
    tree.SetBranchStatus("dataMCWeight",1)
    tree.SetBranchStatus("ZZMass",1)
    tree.SetBranchStatus("Z1Mass",1)
    tree.SetBranchStatus("Z2Mass",1)
    tree.SetBranchStatus("Z1Flav",1)
    tree.SetBranchStatus("Z2Flav",1)
    tree.SetBranchStatus("ZZMassErr",1)
    tree.SetBranchStatus("ZZMassErrCorr",1)
    tree.SetBranchStatus("ZZMassRefit",1)
    tree.SetBranchStatus("ZZMassRefitErr",1)
    tree.SetBranchStatus("p_GG_SIG_ghg2_1_ghz1_1_JHUGen",1)
    tree.SetBranchStatus("p_GG_SIG_ghg2_1_ghz4_1_JHUGen",1)
    tree.SetBranchStatus("p_GG_SIG_ghg2_1_ghz2_1_JHUGen",1)
    tree.SetBranchStatus("p_QQB_SIG_ZPqqLR_1_gZPz2_1_JHUGen",1)
    tree.SetBranchStatus("p_QQB_SIG_ZPqqLR_1_gZPz1_1_JHUGen",1)
    tree.SetBranchStatus("p_GG_SIG_gXg1_1_gXz1_1_gXz5_1_JHUGen",1)
    tree.SetBranchStatus("p_QQB_SIG_XqqLR_1_gXz1_1_gXz5_1_JHUGen",1)
    tree.SetBranchStatus("p_QQB_BKG_MCFM",1)
    tree.SetBranchStatus("p_m4l_SIG",1)
    tree.SetBranchStatus("p_m4l_BKG",1)
    tree.SetBranchStatus("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal",1)
    tree.SetBranchStatus("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal",1)
    tree.SetBranchStatus("p_JQCD_SIG_ghg2_1_JHUGen_JECNominal",1)
    tree.SetBranchStatus("p_JVBF_SIG_ghv1_1_JHUGen_JECNominal",1)
    tree.SetBranchStatus("pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal",1)
    tree.SetBranchStatus("p_HadWH_SIG_ghv1_1_JHUGen_JECNominal",1)
    tree.SetBranchStatus("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal",1)
    tree.SetBranchStatus("ZZPt",1)
    tree.SetBranchStatus("nExtraLep",1)
    tree.SetBranchStatus("nExtraZ",1)
    tree.SetBranchStatus("nCleanedJetsPt30",1)
    tree.SetBranchStatus("nCleanedJetsPt30BTagged",1)
    tree.SetBranchStatus("JetPt",1)
    tree.SetBranchStatus("JetEta",1)
    tree.SetBranchStatus("JetPhi",1)
    tree.SetBranchStatus("JetMass",1)
    tree.SetBranchStatus("JetQGLikelihood",1)
    tree.SetBranchStatus("DiJetMass",1)
    tree.SetBranchStatus("DiJetDEta",1)

    iEntry=0
    while tree.GetEntry(iEntry):

        # print "   Inspecting entry n. ",iEntry,"..."
        iEntry+=1
        run         = tree.RunNumber
        lumi        = tree.LumiNumber
        event       = tree.EventNumber
        CRflag      = tree.CRflag


        theEvent = Event(run,lumi,event)

        #for iCand in range(len(CRflag)):
        #if event != storedEvent:
        if (True) : #I don't want to reindent everything
            for aCR in CRdict.keys():

                if test_bit(CRflag,CRdict[aCR]):

                    crCounter[aCR] += 1
                    mass4l        = tree.ZZMass
                    mZ1           = tree.Z1Mass
                    mZ2           = tree.Z2Mass
                    massErrRaw    = tree.ZZMassErr
                    massErrCorr   = tree.ZZMassErrCorr
                    m4lRefit      = tree.ZZMassRefit
                    m4lRefitErr   = tree.ZZMassRefitErr
                    p_GG_SIG_ghg2_1_ghz1_1_JHUGen  = tree.p_GG_SIG_ghg2_1_ghz1_1_JHUGen
                    p_GG_SIG_ghg2_1_ghz4_1_JHUGen = tree.p_GG_SIG_ghg2_1_ghz4_1_JHUGen
                    p_GG_SIG_ghg2_1_ghz2_1_JHUGen = tree.p_GG_SIG_ghg2_1_ghz2_1_JHUGen
                    p_QQB_SIG_ZPqqLR_1_gZPz2_1_JHUGen  = tree.p_QQB_SIG_ZPqqLR_1_gZPz2_1_JHUGen
                    p_QQB_SIG_ZPqqLR_1_gZPz1_1_JHUGen      = tree.p_QQB_SIG_ZPqqLR_1_gZPz1_1_JHUGen
                    p_GG_SIG_gXg1_1_gXz1_1_gXz5_1_JHUGen      = tree.p_GG_SIG_gXg1_1_gXz1_1_gXz5_1_JHUGen
                    p_QQB_SIG_XqqLR_1_gXz1_1_gXz5_1_JHUGen   = tree.p_QQB_SIG_XqqLR_1_gXz1_1_gXz5_1_JHUGen
                    p_QQB_BKG_MCFM    = tree.p_QQB_BKG_MCFM
                    p_m4l_SIG    = tree.p_m4l_SIG
                    p_m4l_BKG       = tree.p_m4l_BKG
                    p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal    = tree.JJVBF_SIG_ghv1_1_JHUGen_JECNominal
                    p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal    = tree.JJQCD_SIG_ghg2_1_JHUGen_JECNominal
                    p_JQCD_SIG_ghg2_1_JHUGen_JECNominal     = tree.p_JQCD_SIG_ghg2_1_JHUGen_JECNominal
                    p_JVBF_SIG_ghv1_1_JHUGen_JECNominal     = tree.p_JVBF_SIG_ghv1_1_JHUGen_JECNominal
                    pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal     = tree.pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal
                    p_HadWH_SIG_ghv1_1_JHUGen_JECNominal = tree.p_HadWH_SIG_ghv1_1_JHUGen_JECNominal
                    p_HadZH_SIG_ghz1_1_JHUGen_JECNominal = tree.p_HadZH_SIG_ghz1_1_JHUGen_JECNominal
                    pt4l          = tree.ZZPt
                    nExtraLep     = tree.nExtraLep
                    nExtraZ       = tree.nExtraZ
                    jetpt         = tree.JetPt
                    jeteta        = tree.JetEta
                    jetphi        = tree.JetPhi
                    jetmass       = tree.JetMass
                    jetQGLikelihood = tree.JetQGLikelihood
                    njets30       = tree.nCleanedJetsPt30
                    njets30Btag   = tree.nCleanedJetsPt30BTagged
                    mjj           = tree.DiJetMass
                    detajj        = tree.DiJetDEta
                    weight        = sign(tree.genHEPMCweight) * tree.PUWeight * tree.dataMCWeight

                    jets30pt = []
                    jets30eta = []
                    jets30phi = []
                    jets30mass = []
                    jets30QGLikelihood = []

                    for i in range(len(jetpt)):
                        if jetpt[i]>30.:
                            jets30pt.append(jetpt[i])
                            jets30eta.append(jeteta[i])
                            jets30phi.append(jetphi[i])
                            jets30mass.append(jetmass[i])
                            jets30QGLikelihood.append(jetQGLikelihood[i])

                    theKDs = KDs(p_GG_SIG_ghg2_1_ghz1_1_JHUGen,p_GG_SIG_ghg2_1_ghz4_1_JHUGen,p_GG_SIG_ghg2_1_ghz2_1_JHUGen,p_QQB_SIG_ZPqqLR_1_gZPz2_1_JHUGen,p_QQB_SIG_ZPqqLR_1_gZPz1_1_JHUGen,p_GG_SIG_gXg1_1_gXz1_1_gXz5_1_JHUGen,p_QQB_SIG_XqqLR_1_gXz1_1_gXz5_1_JHUGen,p_QQB_BKG_MCFM,p_m4l_SIG,p_m4l_BKG,p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,p_HadWH_SIG_ghv1_1_JHUGen_JECNominal,p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,njets30,jets30QGLikelihood,jets30phi,tree.Z1Flav*tree.Z2Flav,mass4l)
                    theCand = Candidate(theEvent,mass4l,mZ1,mZ2,massErrRaw,massErrCorr,m4lRefit,m4lRefitErr,pt4l,nExtraLep,nExtraZ,jets30pt,jets30eta,jets30phi,jets30mass,njets30,njets30Btag,mjj,detajj,theKDs,weight,jets30QGLikelihood,p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,p_HadWH_SIG_ghv1_1_JHUGen_JECNominal,p_HadZH_SIG_ghz1_1_JHUGen_JECNominal)
                    crCands[aCR].append(theCand)


    crSortedCands = crCands
    for aCR in CRdict.keys():
        # Sort candidates on a run / lumisection / event number basis
        crSortedCands[aCR] = sorted(crCands[aCR], key=attrgetter('eventInfo.run', 'eventInfo.lumi', 'eventInfo.event'))

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
