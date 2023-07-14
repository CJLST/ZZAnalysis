#!/usr/bin/env python
###
# compare events in a nanoAOD and a miniAOD CJLST file and print out differences.
##
from __future__ import print_function
import math
from ROOT import *

region = 'SR'
#region = '3P1F'
#region = '2P2F'
#region = 'SS'

verbose = 1

#cjlstFile   = "../../AnalysisStep/test/ZZ4lAnalysis_sync_ggH_102X.root" 
#nanoFile = "ggH125_Rereco_fixedFSR.root"

# 2017 UL ggH, 26000 events
cjlstFile   = "../../AnalysisStep/test/ZZ4lAnalysis_sync_ggH2017_106X_nomuscale.root" 
nanoFile = "ggH125_2017UL_fixedFSR_nomucorr_full.root"
#nanoFile = "ggH125_2017UL_fixedFSR_Skim.root"

#2022 MC ggH, 13k events
#cjlstFile = "../../AnalysisStep/test/ZZ4lAnalysis_sync_SRCR.root"
#nanoFile = "e5e2fe04-7fb1-43ed-a81a-3acded81f0e7_Skim.root"



compareWeights = False
compareKD = True
compareExtra = True and region == 'SR' # FIXME to be implemented for CRs
massTolerance = 0.05 # in GeV; account for rounding due to data packing

if region == 'SR' :
    treeMiniName = 'ZZTree/candTree'
    nanoPrefix = 'treeNano.ZZCand_'
else :
    treeMiniName = 'CRZLLTree/candTree'
    nanoPrefix = 'treeNano.ZLLCand_'

# definitions for CRs in mini
CRdict = {"SS":21, "2P2F":22, "3P1F":23}
def test_bit(mask, iBit):
    return (mask >> iBit) & 1
    
treeMini = TChain(treeMiniName)
treeMini.Add(cjlstFile)

treeNano = TChain("Events")
treeNano.Add(nanoFile)



iEntryMini=0
nMatch=0;
foundNano=[False]*treeNano.GetEntries() # True = event in nano tree has already been found in mini tree, or can be skipped (no candidate, etc)
maxM4lDiff = 0.

lastfound = -1

h_m4lDiff = TH1F("h_m4lDiff","h_m4lDiff", 2000, -1, 1)

def printLeps_mini(treeMini, prefix="") :
    print(prefix, end="")
    for i in range(4): print(treeMini.LepLepId[i], '{:.3f} {:.3f}'.format(treeMini.LepPt[i], treeMini.LepEta[i]), end=" ")
    print()

def printLeps_nano(treeNano, prefix) :
    print(prefix, end="")
    lepPts = list(treeNano.Electron_pt) + list(treeNano.Muon_pt)
    lepIds = list(treeNano.Electron_pdgId) + list(treeNano.Muon_pdgId)
    lepEtas = list(treeNano.Electron_eta) + list(treeNano.Muon_eta)
    nanoLepIdxs = [eval(nanoPrefix+'Z1l1Idx[iBC]'),eval(nanoPrefix+'Z1l2Idx[iBC]'), eval(nanoPrefix+'Z2l1Idx[iBC]'),eval(nanoPrefix+'Z2l2Idx[iBC]')]
    case1 = False
    case2 = False
    if verbose>0 : 
        for i, lIdx in enumerate(nanoLepIdxs) :
            print(lepIds[lIdx], '{:.3f} {:.3f}'.format(lepPts[lIdx], lepEtas[lIdx]), end=" ")        
            laId = abs(lepIds[lIdx])
            lPt = lepPts[lIdx]
            laEta = abs(lepEtas[lIdx])

            if (laId == 11 and laEta>2.49) or laEta > 2.39 :
                case1 = True
                print ("*", end="")
            if laId == 11 and lPt > 9. and lPt < 11.:
                case2 = True
                print ("**", end="")
        print()
        if case1 : print ("   * : close to eta bound")
        if case2 : print ("   ** : electron close to 10 GeV pt threshold")


while treeMini.GetEntry(iEntryMini):

    iEntryMini+=1

    if verbose>=2 and iEntryMini%100==0 :
        print ("...", iEntryMini)
    
    if region == 'SR' :
        if treeMini.ZZsel<0 : continue
    else :
        if treeMini.ZZsel>=90: continue # Do not consider CRs in events with a SR candidate, as prescribed
        if not test_bit(treeMini.CRflag,CRdict[region]) : continue

    found = False
    iEntryNano = lastfound # Assume events are ordered in both files to speed up things
    end = treeNano.GetEntries()-1
    while iEntryNano < end :
        iEntryNano += 1
#        if foundNano[iEntryNano] : continue # was alredy found: skip
        treeNano.GetEntry(iEntryNano)
        
        if region == 'SR' :     iBC = treeNano.bestCandIdx
        elif region == 'SS' :   iBC = treeNano.ZLLbestSSIdx
        elif region == '2P2F' : iBC = treeNano.ZLLbest2P2FIdx
        elif region == '3P1F' : iBC = treeNano.ZLLbest3P1FIdx


        if iBC < 0 or not treeNano.HLT_passZZ4l: # no candidate passes the selection
            # in this event, or the event does not pass the required triggers
            # (for samples processed with TRIGPASSTHROUGH=True)
            foundNano[iEntryNano] = True
            continue
            
        if treeMini.RunNumber==treeNano.run and treeMini.LumiNumber==treeNano.luminosityBlock and treeMini.EventNumber==treeNano.event :
            foundNano[iEntryNano] = True
            found = True
            lastfound=iEntryNano;
            break

        # FIXME: handle case where events are not ordered: at end ol loop (iEntryNano=treeNano.GetEntries()-1), set iEntryNano=-1, end=lastfound
    
    if found :
        t2_ZZMass = eval(nanoPrefix+'mass[iBC]')
        t2_ZZMassPreFSR=eval(nanoPrefix+'massPreFSR[iBC]')
        t2_Z1Mass=eval(nanoPrefix+'Z1mass[iBC]')
        t2_Z2Mass=eval(nanoPrefix+'Z2mass[iBC]')
        t2_Z1flav=eval(nanoPrefix+'Z1flav[iBC]')
        t2_Z2flav=eval(nanoPrefix+'Z2flav[iBC]')

        if compareExtra :
            t2_nExtraLep=eval(nanoPrefix+'nExtraLep[iBC]')
            t2_nExtraZ=eval(nanoPrefix+'nExtraZ[iBC]')
        
#        t2_dataMC=eval(nanoPrefix+'dataMCWeight[iBC]')
        t2_hasFSR = abs(t2_ZZMass-t2_ZZMassPreFSR)>0.02

        t1_hasFSR = treeMini.fsrPt.size()>0
        ps1=treeMini.p_GG_SIG_ghg2_1_ghz1_1_JHUGen
        pb1=treeMini.p_QQB_BKG_MCFM
        KD_mini = ps1/(ps1+pb1)

        h_m4lDiff.Fill(treeMini.ZZMass-t2_ZZMass)

        m4lDiff=abs(treeMini.ZZMass-t2_ZZMass)
        maxM4lDiff=max(m4lDiff,maxM4lDiff)
 
        if t2_hasFSR != t1_hasFSR:
            print("FSR Diff: "+str(treeMini.RunNumber)+":"+str(treeMini.LumiNumber)+":"+str(treeMini.EventNumber),
                  "mini:", '{:.2f} {:.2f} {:.2f}'.format(treeMini.ZZMass,treeMini.Z1Mass,treeMini.Z2Mass), "hasFSR:",t1_hasFSR ,
                  "nano:", '{:.2f} {:.2f} {:.2f}'.format(t2_ZZMass, t2_Z1Mass, t2_Z2Mass), "hasFSR:", t2_hasFSR)
            

        elif m4lDiff>massTolerance or abs(treeMini.Z1Mass-t2_Z1Mass)>massTolerance or abs(treeMini.Z2Mass-t2_Z2Mass)>massTolerance:
            # check for FSR
            
            print("Mass Diff: "+str(treeMini.RunNumber)+":"+str(treeMini.LumiNumber)+":"+str(treeMini.EventNumber), m4lDiff)
            print("   mini:", '{:.2f} {:.2f} {:.2f} {:.4f} {:.4f} {:.2f}'.format(treeMini.ZZMass,treeMini.Z1Mass,treeMini.Z2Mass, ps1, pb1, KD_mini), "hasFSR:",t1_hasFSR , "\n"
                  "   nano:", '{:.2f} {:.2f} {:.2f} {:.2f}'.format(t2_ZZMass, t2_Z1Mass, t2_Z2Mass, eval(nanoPrefix+'KD[iBC]')), "hasFSR:", t2_hasFSR)
            if region=='SR': #FIXME
                printLeps_mini(treeMini, "   mini: ")                
                printLeps_nano(treeNano, "   nano: ")

        if compareExtra :
            if(t2_nExtraLep!=treeMini.nExtraLep or t2_nExtraZ!=treeMini.nExtraZ) :
                print("nExtra diff:"+str(treeMini.RunNumber)+":"+str(treeMini.LumiNumber)+":"+str(treeMini.EventNumber),
                      treeMini.nExtraLep, t2_nExtraLep, treeMini.nExtraZ, t2_nExtraZ)
        if compareKD:
            KD_nano = eval(nanoPrefix+'KD[iBC]')
            if abs(1.-KD_nano/KD_mini) > 0.005 :
                print("KD diff: "+str(treeMini.RunNumber)+":"+str(treeMini.LumiNumber)+":"+str(treeMini.EventNumber),KD_mini, KD_nano)


#treeMini.Z2Mass

#        h_RWdataMC.Fill(t2_dataMC/treeMini.dataMCWeight)
#        if abs((treeMini.PUWeight-treeNano.puWeight)/treeMini.PUWeight)>0.5 or abs((treeMini.dataMCWeight-t2_dataMC)/treeMini.dataMCWeight)>0.05 :
#            print("W Diff: "+str(treeMini.RunNumber)+":"+str(treeMini.LumiNumber)+":"+str(treeMini.EventNumber), treeMini.PUWeight, treeNano.puWeight, treeMini.dataMCWeight, t2_dataMC)

#        h_ggH_NNLOPS_weight.Fill(treeMini.ggH_NNLOPS_weight,treeNano.ggH_NNLOPS_Weight)

#        h_K_QCDGG.Fill(treeMini.KFactor_QCD_qqZZ_M,treeNano.KFactor_QCD_qqZZ_M)
#        h_K_QCDQQ.Fill(treeMini.KFactor_QCD_ggZZ_Nominal,treeNano.KFactor_QCD_ggZZ_Nominal)

        nMatch+=1
    else :
        print("Missing in nano: "+str(treeMini.RunNumber)+":"+str(treeMini.LumiNumber)+":"+str(treeMini.EventNumber), treeMini.ZZsel)
        print('   {:.2f} {:.2f} {:.2f}'.format(treeMini.ZZMass,treeMini.Z1Mass,treeMini.Z2Mass), "hasFSR:",t1_hasFSR)
        printLeps_mini(treeMini, "   ")

for iEntryNano,found in enumerate(foundNano):
    if not found:
        treeNano.GetEntry(iEntryNano)
        print("Missing in mini: "+str(treeNano.run)+":"+str(treeNano.luminosityBlock)+":"+str(treeNano.event))
        printLeps_nano(treeNano, "   ")

print("Matches in", region, ":", nMatch)
print("max m4l diff:", maxM4lDiff)
