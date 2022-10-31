### 
# -Jet-Lepton cross-cleaning
# -Additional JES, JEC
# FIXME: to be implemented.
###

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR

class jetFiller(Module):
    def __init__(self, era):
        self.era = era

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Jet_lepVeto", "O", lenVar="nJet") # vetoed by one of the candidate's leptons
        self.out.branch("Jet_uncorrectedPt", "I", lenVar="nJet") # Original pT
        self.out.branch("Jet_Pt", "I", lenVar="nJet") # after correction        
        # up/down variations...

    def analyze(self, event):
        jets = Collection(event, 'Jet')
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        leps = list(muons)+ list(electrons)
        nlep=len(leps)

## Apply JES, JER...


## Jet-lepton cleaning with best candidate's leptons and FSR..


## b-tagging..


        return True
