### 
# Add extra objects/variables to the best ZZ and CR candidates in the event, for categorization
# -extra lepts
# -extra Zs
#
###

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
#from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR

class ZZExtraFiller(Module):
    def __init__(self, region):
        self.region = region

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("ZZCand_nExtraLep", "O", lenVar="nZZCand") 
        self.out.branch("ZZCand_nExtraZ", "O", lenVar="nZZCand") 

    def analyze(self, event):
        ZZs = Collection(event, 'ZZCand')
        Zs = Collection(event, 'ZCand')
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        leps = list(muons)+ list(electrons)
        nlep=len(leps)

        candIdx = event.bestCandIdx
        if candIdx <0: return True
        theZZ = ZZs[candIdx]        
        theZZleps = [theZZ.Z1l1Idx, theZZ.Z1l2Idx, theZZ.Z2l1Idx, theZZ.Z2l2Idx]

        # Extra leps
        extraLeps = []
        for i in range(nlep) :
            if i in theZZleps : continue
            if leps[i].ZZFullSel : extraLeps.append(i)

        # Extra Zs
        extraZs = []
        for iZ, Z in enumerate(Zs) :
            if Z.l1Idx in theZZleps or Z.l2Idx in theZZleps : continue
            extraZs.append[iZ]

        # Store the result as a collection variable, even if it is filled only for the best candidate
        nExtraLepsA = [-1]*event.nZZCand
        nExtraZsA = [-1]*event.nZZCand
        nExtraLepsA[candIdx] = len(extraLeps)
        nExtraZsA[candIdx] = len(extraZs)
        
        self.out.fillBranch("ZZCand_nExtraLep", nExtraLepsA)
        self.out.fillBranch("ZZCand_nExtraZ", nExtraZsA)


        return True
