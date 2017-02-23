#!/usr/bin/env python

import ROOT
import math
import optparse
import os, sys
from syncUtils import *
from operator import attrgetter

# define function for parsing options
def parseOptions():

    usage = ('usage: %prog [options] datasetList\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    parser.add_option('-i', '--input', dest='inFile', type='string', default="ZZ4lAnalysis.root",    help='input file')
    parser.add_option('-m', '--synchMode', dest='synchMode', type='string', default='HZZ', help='produce sync file in HZZ/VBS format')
    parser.add_option('-n', '--noOutput', dest='noOutput', action='store_true', default=False, help='do not write sync file in output')
    parser.add_option('-o', '--output', dest='outFile', type='string', default="eventlist.txt",    help='output sync file')
    parser.add_option('-f', '--finalState', dest='finalState', type='string', default="all",    help='final states: all, 4e, 4mu, 2e2mu')
    parser.add_option('-l', '--long', dest='longOutput', action='store_true', default=True,    help='long output')
    parser.add_option('-u', '--unblind', dest='unblind', action='store_true', default=False,    help='Unblinded (affects data only)')
    parser.add_option('-r', '--range', dest='range', type='string', default='full',    help='full,low,high mass region (only for blinded data)')
    parser.add_option('-s', '--selection', dest='selection', type='string', default="REG",    help='select REG/TLE/RSE trees')


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
    selection = opt.selection

    is_loose_ele = False
    if selection == 'REG' :
        tree_name = "ZZTree/candTree"
    elif selection == 'TLE' :
        tree_name = "ZZTreetle/candTree"
        is_loose_ele = True
    elif selection == 'RSE' :
        tree_name = "ZZTreelooseEle/candTree"
        is_loose_ele = True

    print "Processing file: ",inFileName,"..."

    cands = []
    totCounter = 0
    chanCounter = {}

    hfile = ROOT.TFile(inFileName)
    hCounters = hfile.Get("ZZTree/Counters")

    list_final_states = ["4mu","4e","2e2mu"] if not is_loose_ele else ["2m1e1Loose", "3e1Loose"]


    for aChan in list_final_states :

        chanCounter[aChan] = 0

        if finalState!="all" and aChan!=finalState: continue

        isMC = False
#        tree = ROOT.TChain("ZZ"+aChan+"Tree/candTree")
        tree = ROOT.TChain(tree_name)
        tree.Add(inFileName)
        #tree.SetBranchStatus("*",0)

        ## Variables we are interested in for the sync
        #tree.SetBranchStatus("ZZsel",1)
        #tree.SetBranchStatus("RunNumber",1)
        #tree.SetBranchStatus("LumiNumber",1)
        #tree.SetBranchStatus("EventNumber",1)
        #if tree.GetBranch("genHEPMCweight") :
        #    isMC = True
        #    tree.SetBranchStatus("genHEPMCweight",1)
        #    tree.SetBranchStatus("PUWeight",1)
        #    tree.SetBranchStatus("dataMCWeight",1)
        #tree.SetBranchStatus("ZZMass",1)
        #tree.SetBranchStatus("Z1Mass",1)
        #tree.SetBranchStatus("Z2Mass",1)
        #tree.SetBranchStatus("Z1Flav",1)
        #tree.SetBranchStatus("Z2Flav",1)
        #tree.SetBranchStatus("ZZMassErr",1)
        #tree.SetBranchStatus("ZZMassErrCorr",1)
        #if opt.longOutput:
        #    tree.SetBranchStatus("ZZMassRefit",1)
        #    tree.SetBranchStatus("ZZMassRefitErr",1)
        ##tree.SetBranchStatus("p_GG_SIG_ghg2_1_ghz1_1_JHUGen",1)
        ##tree.SetBranchStatus("p_GG_SIG_ghg2_1_ghz4_1_JHUGen",1)
        ##tree.SetBranchStatus("p_GG_SIG_ghg2_1_ghz2_1_JHUGen",1)
        ##tree.SetBranchStatus("p_QQB_SIG_ZPqqLR_1_gZPz2_1_JHUGen",1)
        ##tree.SetBranchStatus("p_QQB_SIG_ZPqqLR_1_gZPz1_1_JHUGen",1)
        ##tree.SetBranchStatus("p_GG_SIG_gXg1_1_gXz1_1_gXz5_1_JHUGen",1)
        ##tree.SetBranchStatus("p_QQB_SIG_XqqLR_1_gXz1_1_gXz5_1_JHUGen",1)
        ##tree.SetBranchStatus("p_QQB_BKG_MCFM",1)
        ##tree.SetBranchStatus("p_m4l_SIG",1)
        ##tree.SetBranchStatus("p_m4l_BKG",1)
        ##tree.SetBranchStatus("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal",1)
        ##tree.SetBranchStatus("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal",1)
        ##tree.SetBranchStatus("p_JQCD_SIG_ghg2_1_JHUGen_JECNominal",1)
        ##tree.SetBranchStatus("p_JVBF_SIG_ghv1_1_JHUGen_JECNominal",1)
        ##tree.SetBranchStatus("pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal",1)
        ##tree.SetBranchStatus("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal",1)
        ##tree.SetBranchStatus("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal",1)
        #tree.SetBranchStatus("ZZPt",1)
        #tree.SetBranchStatus("nExtraLep",1)
        #tree.SetBranchStatus("nExtraZ",1)
        #tree.SetBranchStatus("nCleanedJetsPt30",1)
        #tree.SetBranchStatus("nCleanedJetsPt30BTagged",1)
        #tree.SetBranchStatus("JetPt",1)
        #tree.SetBranchStatus("JetEta",1)
        #tree.SetBranchStatus("JetPhi",1)
        #tree.SetBranchStatus("JetMass",1)
        #tree.SetBranchStatus("JetQGLikelihood",1)
        #tree.SetBranchStatus("DiJetMass",1)
        #tree.SetBranchStatus("DiJetDEta",1)
        #tree.SetBranchStatus("LepLepId",1)
        #tree.SetBranchStatus("LepPt",1)

        iEntry=0
        while tree.GetEntry(iEntry):

            # print "   Inspecting entry n. ",iEntry,"..."
            iEntry+=1
            ZZsel       = tree.ZZsel
            run         = tree.RunNumber
            lumi        = tree.LumiNumber
            event       = tree.EventNumber



            pass_selection = False


            if selection == 'REG' :
                if ZZsel>=90 : pass_selection = True

            if selection == 'RSE' :
                if ZZsel>=120 : pass_selection = True
            if selection == 'TLE' :
                TLE_index = -1
#                print tree.LepLepId[0]
                for i, ID in enumerate(tree.LepLepId) :
                    if abs(ID) == 22 :
                        TLE_index = i
                if TLE_index < 0 :
                    print 'Did not find TLE'
                    continue
                if ZZsel>=120 and tree.TLE_dR_Z < 1.6 and tree.LepPt[TLE_index] >= 30. : pass_selection = True

	    if opt.synchMode == 'VBS' and ZZsel<120 : pass_selection = False

            if pass_selection :

                if not is_loose_ele :
                    ZZflav        = tree.Z1Flav*tree.Z2Flav
                    if (aChan=="4e" and ZZflav != 14641) or (aChan=="4mu" and ZZflav!=28561) or (aChan=="2e2mu" and ZZflav!=20449) : continue
                else :
#                    print 'Event is ', abs(tree.Z1Flav*tree.Z2Flav)
                    ZZflav = abs(tree.Z1Flav*tree.Z2Flav)
                    if not ((aChan=="2m1e1Loose" and (ZZflav==169*242 or ZZflav==169*121)) or (aChan=="3e1Loose" and (ZZflav==121*242 or ZZflav==121*121))) : continue
                mass4l        = tree.ZZMass

                if run>100 and opt.unblind==False :
                ## blind Higgs peak in data
                    isLowUnblind  = (mass4l>=70 and mass4l<=110)
                    isHighUnblind = (mass4l>=150 and mass4l<=500)
                    if opt.range=='low' and not isLowUnblind : continue
                    elif opt.range=='high' and not isHighUnblind : continue
                    elif opt.range=='full' and not (isLowUnblind or isHighUnblind) : continue

                totCounter += 1
                chanCounter[aChan] += 1

                theCand = Candidate(tree, opt)
                cands.append(theCand)



    # Sort candidates on a run / lumisection / event number basis
    sortedCands = sorted(cands, key=attrgetter('run', 'lumi', 'event'))

    if not opt.noOutput:
        # Print in sync format
        ext = opt.finalState
        if selection != "REG" :
            ext += selection
        if opt.unblind==False and opt.range!='all' :
            ext+='_'+opt.range
        outFileName = outFileName.replace(".","_"+ext+".")
        outFile = open(outFileName,"w")
        line = ""

        for aCand in sortedCands:
            line += aCand.printOut(opt)
            line += "\n"


        outFile.write(line)
        outFile.close()

        print "Output written in file: ",outFileName,"\n"
        counterStr=str(int(hCounters.GetBinContent(1)))+"/"+str(int(hCounters.GetBinContent(3)))+"/"+str(int(hCounters.GetBinContent(2)))+"/"+str(int(hCounters.GetBinContent(4)))+"/"+str(int(hCounters.GetBinContent(5)))
        print "## Total/4e/4mu/2e2mu/2l2tau : "+counterStr
    if not is_loose_ele :
        counterStr = str(totCounter) + "/" + str(chanCounter["4e"]) + "/" + str(chanCounter["4mu"]) + "/" + str(chanCounter["2e2mu"])
        print "\n## Selected events all/4e/4mu/2e2mu : "+counterStr
    else :
        counterStr = str(totCounter) + "/" + str(chanCounter["2m1e1Loose"]) + "/" + str(chanCounter["3e1Loose"])
        print "\n## Selected events all/2m1e1Loose/3e1Loose : "+counterStr



if __name__ == "__main__" :
    parseOptions()
    loop()
