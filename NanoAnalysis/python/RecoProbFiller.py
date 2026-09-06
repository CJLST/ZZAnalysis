import copy
import math
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.MELAProbHelper import MELAProbHelper
from ZZAnalysis.NanoAnalysis.tools import branchCollection
import Mela
from operator import itemgetter


class RecoProbFiller(Module):
    """Calculates proabilities with Reco-level information. 
    MELA = Pointer to MELA passed from nanoZZ4lAnalysis.py 
    MELASettings = dictionary with settings for the probabilities to be computed
    processCR = if True, also fill ZLLCand with the same probabilities as ZZCand
    """
    
    def __init__(self, MELA, NANOVERSION, MELASettings=None, processCR=False):
        self.MELA = MELA
        self.NANOVERSION = NANOVERSION
        self.MELASettings = MELASettings
        self.processCR = processCR
        self.EFthreshold = 0.5 # threshold of leptons pt over jet pt to veto the jet, to be reapplied to the varied jets. See the nominal one in jetFiller.py
        self.probHelpers = {"ZZCand": MELAProbHelper(self.MELA, self.MELASettings, "Reco", candColl="ZZCand")}
        if self.processCR:
            self.probHelpers["ZLLCand"] = MELAProbHelper(self.MELA, self.MELASettings, "Reco", candColl="ZLLCand")
        print("***RecoProbFiller: set for: ", self.probHelpers["ZZCand"].names,
              "processCR:", self.processCR, flush=True)
        self.filler = {}
        self.jetVariations = []

    def setJetVariations(self, jetCorrector):
        """Configure the JES/JER kinematic points produced upstream.

        The jet-correction module owns the exact JES labels, so deriving the
        list from that module keeps the MELA output synchronized with either
        the split-11 or total-JES configuration.  Data have no shifted jets.
        """
        hasJetDependentProbabilities = any(
            probHelper._isJetDependent(prob)
            for probHelper in self.probHelpers.values()
            for prob in probHelper.sortedSettings
        )
        if not getattr(jetCorrector, "is_mc", False) or not hasJetDependentProbabilities:
            self.jetVariations = []
        elif isinstance(jetCorrector.evaluator_JES, dict):
            self.jetVariations = [
                f"{label}_{direction}"
                for label in jetCorrector.evaluator_JES
                for direction in ("ScaleUp", "ScaleDn")
            ]
        else:
            self.jetVariations = ["scaleUp", "scaleDn"]
        if getattr(jetCorrector, "is_mc", False) and hasJetDependentProbabilities:
            self.jetVariations.extend(("smearUp", "smearDn"))

        # Jet-varied MELA probabilities belong only to the signal-region
        # ZZCand collection.  Keep nominal probabilities for the ZLL control
        # region, but never book or compute JES/JER probability variations.
        for candColl, probHelper in self.probHelpers.items():
            probHelper.setJetVariations(
                self.jetVariations if candColl == "ZZCand" else []
            )

        print("***RecoProbFiller: jet-varied MELA points:", self.jetVariations, flush=True)


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.filler["ZZCand"] = branchCollection(wrappedOutputTree, lenVar="nZZCand")
        self.bookAngles("ZZCand", self.filler["ZZCand"])
        if self.processCR:
            self.filler["ZLLCand"] = branchCollection(wrappedOutputTree, lenVar="nZLLCand")
            self.bookAngles("ZLLCand", self.filler["ZLLCand"])
        
        if self.MELASettings is not None:
            for probHelper in self.probHelpers.values():
                probHelper.bookProbs(wrappedOutputTree)
            
    def bookAngles(self, collName, filler) :
        filler.branch(collName + "_costheta1", "F", limitedPrecision=12)
        filler.branch(collName + "_costheta2", "F", limitedPrecision=12)
        filler.branch(collName + "_Phi", "F", limitedPrecision=12)
        filler.branch(collName + "_costhetastar", "F", limitedPrecision=12)
        filler.branch(collName + "_Phi1", "F", limitedPrecision=12)
                
    def _buildMELAInputs(self, cands, event, leps, fsrPhotons, jets):
        jets_idx = [i for i in (event.JetLeadingIdx, event.JetSubleadingIdx) if i >= 0]
        jets_MELA = Mela.SimpleParticleCollection_t()
        for idx in jets_idx:
            jet = jets[idx]
            p4 = jet.p4()
            jets_MELA.add_particle(Mela.SimpleParticle_t(0, p4.Px(), p4.Py(), p4.Pz(), p4.E()))

        candsDaughters = [Mela.SimpleParticleCollection_t() for _ in cands]
        candsAssociated = [copy.deepcopy(jets_MELA) for _ in cands]

        for iCand, aCand in enumerate(cands):
            theCandLepIdxs = [aCand.Z1l1Idx, aCand.Z1l2Idx, aCand.Z2l1Idx, aCand.Z2l2Idx]
            theCandLeps = [leps[i] for i in theCandLepIdxs]
            dressedLepsp4 = [self.getDressedP4(lep=l, fsrPhotons=fsrPhotons) for l in theCandLeps]
            daughters = candsDaughters[iCand]
            associated = candsAssociated[iCand]
            daughters.add_particle(Mela.SimpleParticle_t(theCandLeps[0].pdgId, dressedLepsp4[0].Px(), dressedLepsp4[0].Py(), dressedLepsp4[0].Pz(), dressedLepsp4[0].E()))
            daughters.add_particle(Mela.SimpleParticle_t(theCandLeps[1].pdgId, dressedLepsp4[1].Px(), dressedLepsp4[1].Py(), dressedLepsp4[1].Pz(), dressedLepsp4[1].E()))
            daughters.add_particle(Mela.SimpleParticle_t(theCandLeps[2].pdgId, dressedLepsp4[2].Px(), dressedLepsp4[2].Py(), dressedLepsp4[2].Pz(), dressedLepsp4[2].E()))
            daughters.add_particle(Mela.SimpleParticle_t(theCandLeps[3].pdgId, dressedLepsp4[3].Px(), dressedLepsp4[3].Py(), dressedLepsp4[3].Pz(), dressedLepsp4[3].E()))

            extralep_idx = [i for i in (aCand.extraLep1Idx, aCand.extraLep2Idx) if i>=0]
            for idx in extralep_idx:
                lep = leps[idx]
                p4 = lep.p4()
                associated.add_particle(Mela.SimpleParticle_t(lep.pdgId, p4.Px(), p4.Py(), p4.Pz(), p4.E()))

        return candsDaughters, candsAssociated


    def _selectVariedJets(self, jets, variation):
        """Return shifted leading/subleading jets after the nominal selection.

        Jet-lepton cleaning is reevaluated because ``Jet_ZZLepEF`` has the
        nominal jet pt in its denominator.  Nothing from this temporary
        selection is written to the output tree.
        """
        selected = []
        pt_name = variation + "_pt"
        mass_name = variation + "_mass"
        for idx, jet in enumerate(jets):
            if jet.jetId != 6:
                continue
            pt = getattr(jet, pt_name)
            leptonPt = jet.ZZLepEF * jet.pt
            overlapsLeptons = pt <= 0. or leptonPt / pt > self.EFthreshold
            if overlapsLeptons or pt <= jet.ptThreshold:
                continue
            mass = getattr(jet, mass_name)
            selected.append((pt, idx, mass))
        selected.sort(key=itemgetter(0), reverse=True)
        return selected[:2]

    def _buildVariedAssociated(self, cands, leps, jets):
        """Build per-candidate associated objects for every jet variation."""
        associatedVariations = {}
        selectedJetCounts = {}
        associatedLeptons = {}
        for iCand, aCand in enumerate(cands):
            associatedLeptons[iCand] = []
            for idx in (aCand.extraLep1Idx, aCand.extraLep2Idx):
                if idx < 0:
                    continue
                lep = leps[idx]
                p4 = lep.p4()
                associatedLeptons[iCand].append(Mela.SimpleParticle_t(
                    lep.pdgId, p4.Px(), p4.Py(), p4.Pz(), p4.E()
                ))
        for variation in self.jetVariations:
            selectedJets = self._selectVariedJets(jets, variation)
            selectedJetCounts[variation] = len(selectedJets)
            varied = []
            for iCand, aCand in enumerate(cands):
                associated = Mela.SimpleParticleCollection_t()
                for pt, idx, mass in selectedJets:
                    jet = jets[idx]
                    px = pt * math.cos(jet.phi)
                    py = pt * math.sin(jet.phi)
                    pz = pt * math.sinh(jet.eta)
                    energy = math.sqrt(max(mass * mass + px * px + py * py + pz * pz, 0.))
                    associated.add_particle(Mela.SimpleParticle_t(0, px, py, pz, energy))
                for lep in associatedLeptons[iCand]:
                    associated.add_particle(lep)
                varied.append(associated)
            associatedVariations[variation] = varied
        return associatedVariations, selectedJetCounts

    def analyze(self, event):
        leps = Collection(event, 'Lepton')
        fsrPhotons = Collection(event, "FsrPhoton")
        jets = Collection(event, 'Jet')
        selectedJetCountNominal = sum(
            idx >= 0 for idx in (event.JetLeadingIdx, event.JetSubleadingIdx)
        )

        for collName, probHelper in self.probHelpers.items():
            cands = Collection(event, collName)
            theFiller = self.filler[collName]
            candsDaughters, candsAssociated = self._buildMELAInputs(cands, event, leps, fsrPhotons, jets)
            candsAssociatedVariations, selectedJetCounts = self._buildVariedAssociated(cands, leps, jets)
            for c in candsDaughters:
                self.MELA.setInputEvent(c, None, None, False)
                qH, mZ1, mZ2, helcosthetaZ1, helcosthetaZ2, helPhi, costhetastar, phistarZ1 = self.MELA.computeDecayAngles() 
                theFiller.appendValue(collName + "_costheta1", helcosthetaZ1)
                theFiller.appendValue(collName + "_costheta2", helcosthetaZ2)
                theFiller.appendValue(collName + "_Phi", helPhi)
                theFiller.appendValue(collName + "_costhetastar", costhetastar)
                theFiller.appendValue(collName + "_Phi1", phistarZ1)
            theFiller.fillBranches(cands)
            if self.MELASettings is not None:
                probHelper.fillProbs(
                    candsDaughters,
                    candsAssociated,
                    None,
                    candAssociatedVariations=candsAssociatedVariations,
                    selectedJetCounts=selectedJetCounts,
                    selectedJetCountNominal=selectedJetCountNominal,
                )
        return True

    def getDressedP4(self, lep, fsrPhotons):
        '''Returns the dressed 4-momentum including FSR photon if available'''
        p4 = lep.p4()
        if hasattr(lep, 'fsrPhotonIdx') and lep.fsrPhotonIdx >= 0:
            p4 += fsrPhotons[lep.fsrPhotonIdx].p4()
        return p4
