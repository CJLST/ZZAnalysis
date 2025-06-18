### 
# Jet-Lepton cross-cleaning. Should be called after JES/JEC modules.
# Cleaning is based on the fraction of jet pt that is carried by all leptons passing full sel
# and their FSR photons in dR < 0.4. The jet is masked if the fraction is > 0.5.
# TODO:
# - Add b-tagging info (?)
# - take into account jet veto maps for nCleanedJetsPt30, JetLeadingIdx, etc?
###

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.HeppyCore.utils.deltar import deltaR

class jetFiller(Module):
    def __init__(self):
        print("***jetFiller", flush=True)
        self.dRMin = 0.4 # dR between lepton (or FSR) and jet to assume overlap
        self.EFthreshold = 0.5 # threshold of leptons pt over jet pt to veto the jet

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Jet_ZZMask", "O", lenVar="nJet", title="jet is vetoed by selected leptons or FSR photons")
        self.out.branch("Jet_ZZLepEF", "F", lenVar="nJet", title="Fraction of jet pt carried by the vetoing leptons or FSR photons", limitedPrecision=6),
        self.out.branch("JetLeadingIdx", "S", title="index of leading jet after cleaning")
        self.out.branch("JetSubleadingIdx", "S", title="index of subleading jet after cleaning")
        self.out.branch("nCleanedJetsPt30", "B", title="number of cleaned jets above 30 GeV")
        self.out.branch("nCleanedJetsPt30_jesUp", "B", title="number of cleaned jets, up JES variation")
        self.out.branch("nCleanedJetsPt30_jesDn", "B", title="number of cleaned jets, down JES variation")
        # up/down variations...

    def analyze(self, event):
        jets = Collection(event, 'Jet')
        electrons = Collection(event, "Electron")
        fsrs = Collection(event, "FsrPhoton")
        muons = Collection(event, "Muon")
        leps = list(muons)+ list(electrons)
        nlep=len(leps)

        ## Jet-lepton cleaning with best candidate's leptons and FSR..
        mask = [False]*event.nJet
        jet_lepPtF = [0.]*event.nJet
        leadingJetIdx = -1
        subleadingJetIdx =-1
        leadingJetPt = 0.
        subleadingJetPt = 0.
        nCleanedJetsPt30=0
        nCleanedJetsPt30_jesDn=0
        nCleanedJetsPt30_jesUp=0

        # According to the analysis recipe at https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsZZ4lRunIILegacy#Jets, "[jets] must be cleaned with a DeltaR>0.4 cut wrt all tight leptons in the event passing the SIP and isolation cut computed after FSR correction, as well as with all FSR collected photons attached to these leptons."
        # Note: the current implementation on miniAODs (https://github.com/CJLST/ZZAnalysis/blob/Run2UL_22_nano/AnalysisStep/plugins/JetsWithLeptonsRemover.cc) probably does something different than this. To be reviewed.
        for ij, jet in enumerate(jets) :            
            for lep in leps :
                if not lep.ZZFullSel : continue
                if deltaR(lep.eta, lep.phi, jet.eta, jet.phi) < self.dRMin :
                    jet_lepPtF[ij] += lep.pt
                if lep.fsrPhotonIdx >=0 :
                    fsr = fsrs[lep.fsrPhotonIdx]
                    if deltaR(fsr.eta, fsr.phi, jet.eta, jet.phi) < self.dRMin :
                        jet_lepPtF[ij] += fsr.pt
            jet_lepPtF[ij] /= jet.pt
            if jet_lepPtF[ij] > self.EFthreshold :                
                mask[ij] = True
            else : 
                # Decide pt cut based on jet eta:
                if 2.5 < abs(jet.eta) < 3.0:
                    pt_cut = 50.0
                else:
                    pt_cut = 30.0
                if jet.pt > pt_cut : nCleanedJetsPt30 += 1
                #FIXME: add jesUp, jesDn

                # NOTE: The 50 GeV pT cut in the forward eta region (2.5 < |eta| < 3.0)
                # was introduced following JME recommendations
                # Ref: https://gitlab.cern.ch/cms-jetmet/coordination/coordination/-/issues/113
                # TODO: This should be removed once a JSON-level fix is implemented.
                
                # Note: we cannot rely on the fact that the jet collection is sorted by pt since JES can change this.  
                if jet.pt > leadingJetPt:
                    subleadingJetPt = leadingJetPt
                    subleadingJetIdx = leadingJetIdx
                    leadingJetPt = jet.pt
                    leadingJetIdx = ij
                elif jet.pt > subleadingJetPt :
                    subleadingJetPt = jet.pt
                    subleadingJetIdx = ij
        
        self.out.fillBranch("Jet_ZZMask", mask)
        self.out.fillBranch("Jet_ZZLepEF", jet_lepPtF)
        self.out.fillBranch("JetLeadingIdx", leadingJetIdx)
        self.out.fillBranch("JetSubleadingIdx", subleadingJetIdx)
        self.out.fillBranch("nCleanedJetsPt30", nCleanedJetsPt30)
        self.out.fillBranch("nCleanedJetsPt30_jesUp", nCleanedJetsPt30_jesUp)
        self.out.fillBranch("nCleanedJetsPt30_jesDn", nCleanedJetsPt30_jesDn)
        
        return True
