### Analyze MC Truth
# -Optionally dump MC history
# -Add ZZ gen valriables:
# -GenZZFinalState: product of IDs of the four gen ZZ leptons
# -GenZZ_*Idx: ZZ lepton indices in the GenPart collection (FIXME: unsorted)
# -FsrPhoton_genFsrIdx: index of the closest gen FSR from Z->l (e, mu)
# 
###

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR
from ROOT import TLorentzVector

ZMASS = 91.1876
MIN_MZ1 = 40
MAX_MZ1 = 120
MIN_MZ2 = 12
MAX_MZ2 = 120

class genFiller(Module):
    def __init__(self, dump=False):
        print("-----> INIT GEN FILLER <-----")
        self.writeHistFile = False

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("nDressedLeptons", "I")
        self.out.branch("DressedLeptons_pt", "F", lenVar="nDressedLeptons")
        self.out.branch("GenRelIso", "F", lenVar="nDressedLeptons")

    # Find particle's real mother, (first parent in MC history with a different pdgID)
    def Mother(self, part, gen) :
        idxMother= part.genPartIdxMother
        while idxMother>=0 and gen[idxMother].pdgId == part.pdgId:
            idxMother = gen[idxMother].genPartIdxMother
        idMother=0
        if idxMother >=0 : idMother = gen[idxMother].pdgId
        return idxMother, idMother

    # Return the ID of the leptons's parent: 25 for H->Z->l; 23 for Z->l; +-15 for tau->l if genlep is e,mu.
    def getParentID(self, part, gen) :
        pIdx, pID = self.Mother(part, gen)
        if pIdx < 0 : return 0
        ppIdx = gen[pIdx].genPartIdxMother
        if pID == 23 and ppIdx>=0 and gen[ppIdx].pdgId == 25 :
            pID = 25
        return pID

    def dressLeptons(self, genpart, packedpart):
        lep_dressed = TLorentzVector()
        lep_dressed.SetPtEtaPhiM(genpart.pt, genpart.eta, genpart.phi, genpart.mass)
        for pp in packedpart :
            if pp.status != 1: continue
            if pp.pdgId != 22: continue
            dR_lgamma = deltaR(genpart.eta, genpart.phi, pp.eta, pp.phi)
            idMatch = False
            if pp.genPartIdxMother < 0: continue
            if (packedpart[pp.genPartIdxMother].pdgId == genpart.pdgId): idMatch = True

            if not idMatch: continue

            if dR_lgamma < 0.3:
                lep_dressed += pp.p4()

        return lep_dressed

    def computeGenIso(self, current_lepton, packedpart):
        genIso = 0.0
        for pp in packedpart :
            if pp.status != 1: continue
            if ((abs(pp.pdgId) != 12) or (abs(pp.pdgId) != 14) or (abs(pp.pdgId) != 16)): continue 
            if ((abs(pp.pdgId) == 11) or (abs(pp.pdgId) == 13)): continue
            # TODO: include if (gen_fsrset.find(k)!=gen_fsrset.end()) continue;
            dRvL = deltaR(current_lepton.Eta(), current_lepton.Phi(), pp.eta, pp.phi)
            if dRvL<0.3:
                genIso += pp.pt
        return genIso

    def GenHiggsCounter(self, nGENHiggs):
        nGENHiggs += 1
        # TODO : Create H cand from gp

    def buildLLPair(self, l1, l2):
        l_a = TLorentzVector()
        l_b = TLorentzVector()
        l_a.SetPtEtaPhiM(l1.Pt(), l1.Eta(), l1.Phi(), l1.M())
        l_b.SetPtEtaPhiM(l2.Pt(), l2.Eta(), l2.Phi(), l2.M())

        return l_a, l_b

    def buildZ1Mass(self, Leptons, LeptonsId, makeCuts):
        offshell = 999.0
        findZ1 = False
        idx_l1 = 0
        idx_l2 = 0
        for i, l1 in enumerate(Leptons):
            for j, l2 in enumerate(Leptons):
                if j <= i : continue
                if ((LeptonsId[i] + LeptonsId[j]) != 0) : continue
                l_i, l_j = self.buildLLPair(l1, l2)

                # TODO: Implement pt, eta, Iso cuts
                ## if (makeCuts and self.checkCuts()) : continue

                mll = TLorentzVector()
                mll = l_i + l_j
                if(abs(mll.M() - ZMASS) < offshell) :
                    mZ1 = mll.M()
                    idx_l1 = i
                    idx_l2 = j
                    findZ1 = True
                    offshell = abs(mZ1 - ZMASS)

        return offshell, findZ1, idx_l1, idx_l2

    def buildZ2Mass(self, Leptons, LeptonsId, idx_l1, idx_l2, makeCuts):
        pT_l3l4 = 0.0
        idx_l3 = 0
        idx_l4 = 0
        findZ2 = False
        for i, l1 in enumerate(Leptons):
            if ((i == idx_l1) or (i == idx_l2)) : continue
            for j, l2 in enumerate(Leptons):
                if j <= i : continue
                if ((j == idx_l1) or (j == idx_l2)) : continue
                if ((LeptonsId[i] + LeptonsId[j]) != 0) : continue
                l_i, l_j = self.buildLLPair(l1, l2)

                Z2 = TLorentzVector()
                Z2 = l_i + l_j

                # TODO: Implement pt, eta, Iso cuts
                ## if (makeCuts and self.checkCuts()) : continue

                if (l_i.Pt() + l_j.Pt() >= pT_l3l4):
                    mass_Z2 = Z2.M()
                    if ((mass_Z2>MIN_MZ2 and mass_Z2<MAX_MZ2) or (not makeCuts)) :
                        findZ2 = True
                        idx_l3 = i
                        idx_l4 = j
                        pT_l3l4 = l_i.Pt() + l_j.Pt()
                    else :
                        if not findZ2:
                            idx_l3 = i
                            idx_l4 = j
        return findZ2, idx_l3, idx_l4

    def buildZMasses(self, Leptons, LeptonsId, makeCuts = False):
        passFidSel = False
        passZ1 = False

        offshell, findZ1, idx_l1, idx_l2 = self.buildZ1Mass(Leptons, LeptonsId, makeCuts)

        z1_l1 = Leptons[idx_l1]
        z1_l2 = Leptons[idx_l2]
        l1, l2 = self.buildLLPair(z1_l1, z1_l2)
        ml1l2 = TLorentzVector()
        ml1l2 = l1 + l2

        if (ml1l2.M()>MIN_MZ1 and ml1l2.M()<MAX_MZ1 and findZ1) : passZ1 = True
        if not makeCuts : passZ1 = True

        findZ2, idx_l3, idx_l4 = self.buildZ2Mass(Leptons, LeptonsId, idx_l1, idx_l2, makeCuts)

        z_leps_idx = [idx_l1, idx_l2, idx_l3, idx_l4]

        if (passZ1 and findZ2) : passFidSel = True

        return passFidSel, z_leps_idx

    def buildZCands(self, Leptons, LeptonsId):
        passFidSel_noCuts, z_leps_idx = self.buildZMasses(Leptons, LeptonsId)
        if passFidSel_noCuts:
            l1 = Leptons[z_leps_idx[0]]; l2 = Leptons[z_leps_idx[1]]
            l3 = Leptons[z_leps_idx[2]]; l4 = Leptons[z_leps_idx[3]]
            Z1_l1, Z1_l2 = self.buildLLPair(l1, l2)
            Z2_l1, Z2_l2 = self.buildLLPair(l3, l4)
        return Z1_l1, Z1_l2, Z2_l1, Z2_l2

    def init_collections(self):
        Leptons = []
        LeptonsId = []
        LeptonsReco = []

        return Leptons, LeptonsId, LeptonsReco

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        genpart=Collection(event,"GenPart")

        dressedLeptons = [-1]*len(genpart)
        Lepts_RelIso   = [-1]*len(genpart)

        nGENHiggs = 0.0

        Leptons, LeptonsId, LeptonsReco = self.init_collections()

        for i, gp in enumerate(genpart) :
            if ((abs(gp.pdgId) == 11) or (abs(gp.pdgId) == 13) or (abs(gp.pdgId) == 15)) :
                if (not((gp.status == 1) or (abs(gp.pdgId) == 15))): continue
                # TODO: Check conditions on motherID
                # if (!(genAna.MotherID(&genParticles->at(j))==23 || genAna.MotherID(&genParticles->at(j))==443 || genAna.MotherID(&genParticles->at(j))==553 || abs(genAna.MotherID(&genParticles->at(j)))==24) ) continue;

                # Dress leptons
                # PackedGenParticles in miniAOD is GenPart.status == 1
                lep_dressed = self.dressLeptons(gp, genpart)
                Leptons.append(lep_dressed)
                LeptonsId.append(gp.pdgId)
                # LeptonsReco # TODO : Figure out (&genParticles->at(j));

                current_lepton = lep_dressed
                genIso = self.computeGenIso(current_lepton, genpart)
                genIso = genIso / current_lepton.Pt()

                if (gp.pdgId == 25): self.GenHiggsCounter(nGENHiggs)

                Lepts_RelIso[i] = genIso
                dressedLeptons[i] = lep_dressed.Pt()

        if(len(Leptons)>=4) :
            Z1_l1, Z1_l2, Z2_l1, Z2_l2 = self.buildZCands(Leptons, LeptonsId)

        self.out.fillBranch("nDressedLeptons", len(dressedLeptons))
        self.out.fillBranch("DressedLeptons_pt", dressedLeptons)
        self.out.fillBranch("GenRelIso", Lepts_RelIso)
        return True

