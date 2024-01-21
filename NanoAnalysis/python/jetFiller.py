### 
# -Jet-Lepton cross-cleaning. Should be called after JES/JEC modules.
# TODO:
# - to be implemented
# - Add b-tagging info (?)
###

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.HeppyCore.utils.deltar import deltaR

class jetFiller(Module):
    def __init__(self):
        print("***jetFiller", flush=True)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Jet_ZZMask", "O", lenVar="nJet") # jet is vetoed by one of the candidate's leptons
        self.out.branch("JetLeadingIdx", "S") # index of leading jet after cleaning
        self.out.branch("JetSubleadingIdx", "S") #index of subleading jet after cleaning
        self.out.branch("nCleanedJetsPt30", "B") #number of jets above 30 GeV
        self.out.branch("nCleanedJetsPt30_jesUp", "B")
        self.out.branch("nCleanedJetsPt30_jesDn", "B")
        # up/down variations...

    def analyze(self, event):
        jets = Collection(event, 'Jet')
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        leps = list(muons)+ list(electrons)
        nlep=len(leps)

        ## Jet-lepton cleaning with best candidate's leptons and FSR..
        mask = [False]*event.nJet
        leadingJetIdx = -1
        subleadingJetIdx =-1
        nCleanedJetsPt30=0
        nCleanedJetsPt30_jesDn=0
        nCleanedJetsPt30_jesUp=0
        # TO BE IMPLEMENTED
        # According to the analysis recipe at https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsZZ4lRunIILegacy#Jets, "[jets] must be cleaned with a DeltaR>0.4 cut wrt all tight leptons in the event passing the SIP and isolation cut computed after FSR correction, as well as with all FSR collected photons attached to these leptons."
        # Note: the current implementation on miniAODs (https://github.com/CJLST/ZZAnalysis/blob/Run2UL_22_nano/AnalysisStep/plugins/JetsWithLeptonsRemover.cc) probably does something different than this. To be reviewed.


        self.out.fillBranch("Jet_ZZMask", mask)
        self.out.fillBranch("JetLeadingIdx", leadingJetIdx)
        self.out.fillBranch("JetSubleadingIdx", subleadingJetIdx)
        self.out.fillBranch("nCleanedJetsPt30", nCleanedJetsPt30)
        self.out.fillBranch("nCleanedJetsPt30_jesUp", nCleanedJetsPt30_jesUp)
        self.out.fillBranch("nCleanedJetsPt30_jesDn", nCleanedJetsPt30_jesDn)
        
        return True
