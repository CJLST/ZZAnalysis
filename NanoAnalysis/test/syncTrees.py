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


cjlstFile   = "../../AnalysisStep/test/ZZ4lAnalysis.root" 
nanoFile = "nanoZZ4lAnalysis.root"

compareWeights = False
compareKD = False
massTolerance = 0.1 # in GeV; account for rounding due to data packing
massTolerance = 0.5 #FIXME

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
foundNano=[False]*treeNano.GetEntries()

while treeMini.GetEntry(iEntryMini):

    iEntryMini+=1
    if region == 'SR' :
        if treeMini.ZZsel<0 : continue
    else :
        if treeMini.ZZsel>=90: continue # Do not consider CRs in events with a SR candidate, as prescribed
        if not test_bit(treeMini.CRflag,CRdict[region]) : continue

    iEntryNano=0
    found = False
    while treeNano.GetEntry(iEntryNano):
        iEntryNano+=1
        if foundNano[iEntryNano-1] : continue # was alredy found: skip

        if region == 'SR' :     iBC = treeNano.bestCandIdx
        elif region == 'SS' :   iBC = treeNano.ZLLbestSSIdx
        elif region == '2P2F' : iBC = treeNano.ZLLbest2P2FIdx
        elif region == '3P1F' : iBC = treeNano.ZLLbest3P1FIdx


        if iBC < 0 : # no best candidate in this event
            foundNano[iEntryNano-1] = True
            continue
            
        if treeMini.RunNumber==treeNano.run and treeMini.LumiNumber==treeNano.luminosityBlock and treeMini.EventNumber==treeNano.event :
            foundNano[iEntryNano-1] = True
            found = True
            break

    if found :
        t2_ZZMass = eval(nanoPrefix+'mass[iBC]')
        t2_ZZMassPreFSR=eval(nanoPrefix+'massPreFSR[iBC]')
        t2_Z1Mass=eval(nanoPrefix+'Z1mass[iBC]')
        t2_Z2Mass=eval(nanoPrefix+'Z2mass[iBC]')
        t2_Z1flav=eval(nanoPrefix+'Z1flav[iBC]')
        t2_Z2flav=eval(nanoPrefix+'Z2flav[iBC]')
#        t2_dataMC=eval(nanoPrefix+'dataMCWeight[iBC]')
        t2_hasFSR = abs(t2_ZZMass-t2_ZZMassPreFSR)>0.02

        t1_hasFSR = treeMini.fsrPt.size()>0
        ps1=treeMini.p_GG_SIG_ghg2_1_ghz1_1_JHUGen
        pb1=treeMini.p_QQB_BKG_MCFM
        KD_mini = ps1/(ps1+pb1)

 
        if t2_hasFSR != t1_hasFSR:
            print("FSR Diff: "+str(treeMini.RunNumber)+":"+str(treeMini.LumiNumber)+":"+str(treeMini.EventNumber),
                  "mini:", '{:.2f} {:.2f} {:.2f}'.format(treeMini.ZZMass,treeMini.Z1Mass,treeMini.Z2Mass), "hasFSR:",t1_hasFSR ,
                  "nano:", '{:.2f} {:.2f} {:.2f}'.format(t2_ZZMass, t2_Z1Mass, t2_Z2Mass), "hasFSR:", t2_hasFSR)
            

        elif abs(treeMini.ZZMass-t2_ZZMass)>massTolerance or abs(treeMini.Z1Mass-t2_Z1Mass)>massTolerance or abs(treeMini.Z2Mass-t2_Z2Mass)>massTolerance:
            # check for FSR
            
            print("Mass Diff: "+str(treeMini.RunNumber)+":"+str(treeMini.LumiNumber)+":"+str(treeMini.EventNumber),
                  "mini:", '{:.2f} {:.2f} {:.2f} {:.4f} {:.4f} {:.2f}'.format(treeMini.ZZMass,treeMini.Z1Mass,treeMini.Z2Mass, ps1, pb1, KD_mini), "hasFSR:",t1_hasFSR ,
                  "nano:", '{:.2f} {:.2f} {:.2f} {:.2f}'.format(t2_ZZMass, t2_Z1Mass, t2_Z2Mass, eval(nanoPrefix+'KD[iBC]')), "hasFSR:", t2_hasFSR)
            if region=='SR': #FIXME
                print ("   mini ", end ="")
                for i in range(4): print("lep"+str(i), treeMini.LepPt[i],treeMini.LepLepId[i], end=" ")
                nanoLepIdxs = [eval(nanoPrefix+'Z1l1Idx[iBC]'),eval(nanoPrefix+'Z1l2Idx[iBC]'), eval(nanoPrefix+'Z2l1Idx[iBC]'),eval(nanoPrefix+'Z2l2Idx[iBC]')]
                lepPts = list(treeNano.Muon_pt) + list(treeNano.Electron_pt)
                lepIds = list(treeNano.Muon_pdgId) + list(treeNano.Electron_pdgId)
                print ("\n   nano ", end="")
                for i, lIdx in enumerate(nanoLepIdxs) :
                    print(lepPts[lIdx], lepIds[lIdx], end=" ")
                    print()
                    print(lepPts)


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

for iEntryNano,found in enumerate(foundNano):
    if not found:
        treeNano.GetEntry(iEntryNano)
        print("Missing in mini: "+str(treeNano.run)+":"+str(treeNano.luminosityBlock)+":"+str(treeNano.event))

print("Matches in", region, ":", nMatch)
