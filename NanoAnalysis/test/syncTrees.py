#!/usr/bin/env python3
###
# compare events in nanoAOD trees and print out differences.
##
import math
from ROOT import *
PyConfig.IgnoreCommandLineOptions = True
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

import sys
if len(sys.argv) !=3 :
    print ("Please specify the 2 root files to compare")
    
#file1 = "25c8f5ff-9de0-4a0c-9e2f-757332ad392f_Skim_ref.root"
#file2 = "25c8f5ff-9de0-4a0c-9e2f-757332ad392f_Skim.root"

file1 = sys.argv[1]
file2 = sys.argv[2]

region = 'SR'
#region = '3P1F'
#region = '2P2F'
#region = 'SS'

def selectFinalState(aCand) :
    return True; # all candidates
#    return aCand.Z1flav*aCand.Z2flav == 28561 #4mu only


checkExtra = False
verbose = 0


compareZZVars=["mass",
               "massPreFSR",
               "Z1mass",
               "Z2mass",
               "Z1flav",
               "Z2flav",
               "KD",
               "dataMCWeight",
               ]

if checkExtra:
    compareZZVars.extend(["nExtraLep",
                          "nExtraZ",
                          ])

#FIXME to be added
compareEvtVars = ["dataMCWeight",
                  "ggH_NNLOPS_weight,",
                  "KFactor_QCD_qqZZ_M",
                  "KFactor_QCD_ggZZ_Nominal",
                  ]


compareWeights = False
compareKD = True
compareExtra = True and region == 'SR' # FIXME to be implemented for CRs
massTolerance = 0.05 # in GeV; account for rounding due to data packing

# definitions for CRs in mini
CRdict = {"SS":21, "2P2F":22, "3P1F":23}
    
tree1 = TChain("Events")
tree1.Add(file1)

tree2 = TChain("Events")
tree2.Add(file2)

#mini=1, nano=2

iEntry1=0
nMatch=0
nDiffer=0
missing_1=0
missing_2=0
found2=[False]*tree2.GetEntries() # True = event in nano tree has already been found in mini tree, or can be skipped (no candidate, etc)
maxM4lDiff = 0.

lastfound = -1

h_m4lDiff = TH1F("h_m4lDiff","m4l_1-m4l_2", 2000, -1, 1)


def printLeps(tree, cand) :
    lepPts = list(tree.Electron_pt) + list(tree.Muon_pt)
    lepIds = list(tree.Electron_pdgId) + list(tree.Muon_pdgId)
    lepEtas = list(tree.Electron_eta) + list(tree.Muon_eta)
    lepPhis = list(tree.Electron_phi) + list(tree.Muon_phi)
    lepIdxs = [cand.Z1l1Idx,cand.Z1l2Idx,cand.Z2l1Idx,cand.Z2l2Idx]
    for i, lIdx in enumerate(lepIdxs) :
        print("   ", lepIds[lIdx], '{:.5f} {:.3f} {:.3f}'.format(lepPts[lIdx], lepEtas[lIdx], lepPhis[lIdx]), end=" ")
    print()

while tree1.GetEntry(iEntry1):
#    print("tree1: "+str(tree1.RunNumber)+":"+str(tree1.LumiNumber)+":"+str(tree1.EventNumber))
    iEntry1+=1

    if verbose>=2 and iEntry1%100==0 :
        print ("...", iEntry1)

    if region == 'SR' :     iBC1 = tree2.bestCandIdx
    elif region == 'SS' :   iBC1 = tree2.ZLLbestSSIdx
    elif region == '2P2F' : iBC1 = tree2.ZLLbest2P2FIdx
    elif region == '3P1F' : iBC1 = tree2.ZLLbest3P1FIdx
    
    if iBC1<0 : # or not tree1.HLT_passZZ4l: ### FIXME
        continue

    found = False
    iEntry2 = max(0,lastfound-3) # Assume events are approximately ordered in both files to speed up things
    # iEntry2 = 0 # Re-start from the beginning, if files are in random order (very slow)
    while tree2.GetEntry(iEntry2) :
        thisEntry2 = iEntry2
        iEntry2 += 1
        if found2[thisEntry2] : continue # was alredy found: skip
#        print(" tree2: "+str(tree2.run)+":"+str(tree2.luminosityBlock)+":"+str(tree2.event))
        
        if region == 'SR' :     iBC2 = tree2.bestCandIdx
        elif region == 'SS' :   iBC2 = tree2.ZLLbestSSIdx
        elif region == '2P2F' : iBC2 = tree2.ZLLbest2P2FIdx
        elif region == '3P1F' : iBC2 = tree2.ZLLbest3P1FIdx
        
        if iBC2 < 0 : # or not tree2.HLT_passZZ4l: ### FIXME
            found2[thisEntry2] = True # Setting this as found since it should not be checked further
            continue
            
        if tree1.run==tree2.run and tree1.luminosityBlock==tree2.luminosityBlock and tree1.event==tree2.event :
            found2[thisEntry2] = True
            found = True
            lastfound=thisEntry2;
            break
    
    if found : #event has a candidate in both trees
        ZZs_1 = Collection(tree1, 'ZZCand')
        ZZs_2 = Collection(tree2, 'ZZCand')
        theZZ1 = ZZs_1[iBC1]
        theZZ2 = ZZs_2[iBC2]

        differ = False
        if selectFinalState(theZZ1) :
            for var in compareZZVars :
                val1 = eval("theZZ1."+var)
                val2 = eval("theZZ2."+var)
                if  val1 != val2 :
                    if differ == False :
                        differ = True
                        print(f"{tree1.run}:{tree1.luminosityBlock}:{tree1.event}: candidate differences:")
                        print(f"  mass: {theZZ1.mass} {theZZ2.mass}")
                        print(f"  Z1mass: {theZZ1.Z1mass} {theZZ2.Z1mass}")
                        print(f"  Z2mass: {theZZ1.Z2mass} {theZZ2.Z2mass}")
                    if var not in ["mass", "Z1mass", "Z2mass"] :
                        print(f"  {var}: {val1} {val2}")
                    if var=="mass" : # may be a difference in FSR?
                        if theZZ1.massPreFSR == theZZ2.massPreFSR : 
                            print("  Note: same massPreFSR: {theZZ1.massPreFSR}. Probably a difference in FSR")
                            # Fixme add debug
                        
        if differ :
            nDiffer+=1
            if region=='SR': #FIXME
                print(f"  Leptons in {file1} :")
                printLeps(tree1, theZZ1)
                print(f"  Leptons in {file2} :")
                printLeps(tree2, theZZ2)

            m4lDiff=theZZ1.mass-theZZ2.mass
            h_m4lDiff.Fill(m4lDiff)
            maxM4lDiff=max(abs(m4lDiff),maxM4lDiff)

        nMatch+=1

    else :
        print(f"Missing in tree2: {tree1.run}:{tree1.luminosityBlock}:{tree1.event}")
        missing_2 +=1
        if verbose > 0 :
            printLeps(tree1, "   ")

for iEntry2,found in enumerate(found2):
    if not found:
        tree2.GetEntry(iEntry2)
        print(f"Missing in tree1: {tree2.run}:{tree2.luminosityBlock}:{tree2.event}")
        missing_1 +=1
        if verbose > 0 :
            printLeps(tree2, "   ")

print("Matches in", region, ":", nMatch, '(+'+str(missing_1)+',-'+str(missing_2)+')')
print(f" With candidate differences: {nDiffer}")
print(f" max m4l diff: {maxM4lDiff}")

# c1 = TCanvas("M4ldiff","M4ldiff")
# c1.SetLogy()
# h_m4lDiff.GetXaxis().SetRangeUser(-0.1, 0.1)
# h_m4lDiff.GetXaxis().SetTitle("m_{mini} - m_{nano} [GeV]")
# h_m4lDiff.Draw()
# c1.Print("m4ldiff_"+flavour+".png")
