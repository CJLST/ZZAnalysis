### 
# Add extra objects/variables to the best ZZ and CR candidates in the event, for categorization
# -extra lepts
# -extra Zs
#
# TODO:
# -generalize to add the same quantities for CRs (ZLLCand) (move actual code into a function to be called for each candidate)
# -add other categorization variables (MELA discriminants besides KD, etc)
# -move here the computation of data/MC scale factors here from ZZFiller,
# so that it is not called for the whole combinatorial
#
###

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
#from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR

class ZZExtraFiller(Module):
    def __init__(self, region):
        print("***ZZExtraFiller", flush=True)
        self.region = region

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("ZZCand_nExtraLep", "I", lenVar="nZZCand")
        self.out.branch("ZZCand_nExtraZ", "I", lenVar="nZZCand")

    def analyze(self, event):
        ZZs = Collection(event, 'ZZCand')
        Zs = Collection(event, 'ZCand')
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        leps = list(electrons) + list(muons)
        nlep=len(leps)

        candIdx = event.bestCandIdx
        if candIdx <0: return True
        theZZ = ZZs[candIdx]        
        theZZleps = [theZZ.Z1l1Idx, theZZ.Z1l2Idx, theZZ.Z2l1Idx, theZZ.Z2l2Idx]

        # Extra leps
        extraLeps = []
        for i in range(nlep) :
 #           print(i, leps[i].pt, leps[i].ZZFullSel)
            if i in theZZleps : continue
            if leps[i].ZZFullSel : extraLeps.append(i)

#        print("---ZZLeps:", theZZleps)
#        print(extraLeps)

        # Extra Zs
        extraZs = []
        for iZ, Z in enumerate(Zs) :
            if Z.l1Idx in theZZleps or Z.l2Idx in theZZleps : continue
            extraZs.append(iZ)

        # Store the result as a collection variable, even if it is filled only for the best candidate
        nExtraLepsA = [-1]*event.nZZCand
        nExtraZsA = [-1]*event.nZZCand
        nExtraLepsA[candIdx] = len(extraLeps)
        nExtraZsA[candIdx] = len(extraZs)
        
        self.out.fillBranch("ZZCand_nExtraLep", nExtraLepsA)
        self.out.fillBranch("ZZCand_nExtraZ", nExtraZsA)


        return True
