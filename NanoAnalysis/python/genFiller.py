### Gen Filler
# -Define fiducial selection at gen-level
# -Save bool for events that pass fiducial sel:
# -Save properties of Gen Z cands and corresponding leptons
# 
###

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR
from ROOT import TLorentzVector

from ZZAnalysis.NanoAnalysis.tools import Mother

ZMASS = 91.1876
MIN_MZ1 = 40
MAX_MZ1 = 120
MIN_MZ2 = 12
MAX_MZ2 = 120

class genFiller(Module):
    '''
        Module that builds gen-level Z and ZZ candidates
        that satisfy the fiducial selection for the HZZ analysis.
    '''
    def __init__(self, dump=False):
        self.printGenHist = dump # print MC history

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("nFidDressedLeps", "I")
        self.out.branch("FidDressedLeps_pt", "F", lenVar="nFidDressedLeps")
        self.out.branch("FidDressedLeps_RelIso", "F", lenVar="nFidDressedLeps")
        self.out.branch("FidZZ_Z1l1Idx", "I") # Indices in the GenPart
        self.out.branch("FidZZ_Z1l2Idx", "I")
        self.out.branch("FidZZ_Z2l1Idx", "I")
        self.out.branch("FidZZ_Z2l2Idx", "I")
        self.out.branch("FidZZ_mass", "F")
        self.out.branch("FidZZ_pt", "F")
        self.out.branch("FidZZ_eta", "F")
        self.out.branch("FidZZ_phi", "F")
        self.out.branch("FidZZ_rapidity", "F")
        self.out.branch("FidZ1_mass", "F")
        self.out.branch("FidZ2_mass", "F")
        self.out.branch("passedFiducial", "B")

    def dressLeptons(self, genpart, packedpart):
        '''
            Dress leptons from gammas with dR(l,gamma)<0.3
            Save the dressed leptons collection and the
            indices of the gammas used to dress the leptons
        '''
        fsr_gamma_idx = []
        lep_dressed = TLorentzVector()
        lep_dressed.SetPtEtaPhiM(genpart.pt, genpart.eta, genpart.phi, genpart.mass)
        for idx, pp in enumerate(packedpart) :
            if pp.status != 1: continue
            if pp.pdgId != 22: continue
            dR_lgamma = deltaR(genpart.eta, genpart.phi, pp.eta, pp.phi)
            idMatch = False
            if pp.genPartIdxMother < 0: continue
            if (packedpart[pp.genPartIdxMother].pdgId == genpart.pdgId): idMatch = True

            if not idMatch: continue

            if dR_lgamma < 0.3:
                fsr_gamma_idx.append(idx)
                lep_dressed += pp.p4()

        return lep_dressed, fsr_gamma_idx

    def computeGenIso(self, current_lepton, packedpart, fsr_gamma_idx):
        '''
            Compute isolation for gen-level leptons.
            Exclude photons used for FSR recovery from the computation.
        '''
        genIso = 0.0
        for idx, pp in enumerate(packedpart) :
            if pp.status != 1: continue
            if ((abs(pp.pdgId) != 12) or (abs(pp.pdgId) != 14) or (abs(pp.pdgId) != 16)): continue 
            if ((abs(pp.pdgId) == 11) or (abs(pp.pdgId) == 13)): continue
            if (idx in fsr_gamma_idx): continue
            dRvL = deltaR(current_lepton.Eta(), current_lepton.Phi(), pp.eta, pp.phi)
            if dRvL<0.3:
                genIso += pp.pt
        genIso = genIso / current_lepton.Pt()
        return genIso

    def GenHiggsCounter(self, nGENHiggs):
        nGENHiggs += 1
        # TODO : Create H cand from gp

    def buildLLPair(self, l1, l2):
        '''
            Util function that returns a pair of TLorentzVectors. 
        '''
        l_a = TLorentzVector()
        l_b = TLorentzVector()
        l_a.SetPtEtaPhiM(l1.Pt(), l1.Eta(), l1.Phi(), l1.M())
        l_b.SetPtEtaPhiM(l2.Pt(), l2.Eta(), l2.Phi(), l2.M())

        return l_a, l_b

    def unzipLeptons(self, LeptonsCollection):
        '''
            Util function that unzips a LeptonsCollection.
            LeptonsCollection contains Leptons, LeptonsId and
            Lepts_RelIso.
        '''
        Leptons, LeptonsId, Lepts_RelIso = LeptonsCollection
        return Leptons, LeptonsId, Lepts_RelIso

    def checkCuts(self, LeptonsCollection, idx_1, idx_2):
        '''
            Util function that checks if two leptons satisfy the
            quality (kinematics and Iso) cuts used to define the
            fiducial volume.
        '''
        passCuts = True
        Leptons, LeptonsId, Lepts_RelIso = self.unzipLeptons(LeptonsCollection)
        id_1 = LeptonsId[idx_1]; id_2 = LeptonsId[idx_2]
        l_1 = Leptons[idx_1]; l_2 = Leptons[idx_2]
        iso_1 = Lepts_RelIso[idx_1]; iso_2 = Lepts_RelIso[idx_2]

        if (abs(id_1) == 13 and (l_1.Pt() < 5.0 or abs(l_1.Eta()) > 2.4)) :
            passCuts = False
        if (abs(id_1) == 11 and (l_1.Pt() < 7.0 or abs(l_1.Eta()) > 2.5)) :
            passCuts = False
        if (iso_1 > 0.35) :
            passCuts = False

        if (abs(id_2) == 13 and (l_2.Pt() < 5.0 or abs(l_2.Eta()) > 2.4)) :
            passCuts = False
        if (abs(id_2) == 11 and (l_2.Pt() < 7.0 or abs(l_2.Eta()) > 2.5)) :
            passCuts = False
        if (iso_2 > 0.35) :
            passCuts = False

        return passCuts

    def buildZ1Mass(self, LeptonsCollection, makeCuts):
        '''
            Build the Z1 candidate (neareast to true Z mass).
            Returns the difference in mass between Z1 and Z,
            a bool that specifies if a Z1 is present in the event
            and the indices of the two leptons used to build Z1.
        '''
        Leptons, LeptonsId, Lepts_RelIso = self.unzipLeptons(LeptonsCollection)
        offshell = 999.0
        findZ1 = False
        idx_l1 = 0
        idx_l2 = 0
        for i, l1 in enumerate(Leptons):
            for j, l2 in enumerate(Leptons):
                if j <= i : continue
                if ((LeptonsId[i] + LeptonsId[j]) != 0) : continue
                l_i, l_j = self.buildLLPair(l1, l2)

                passCuts = self.checkCuts(LeptonsCollection, i, j)
                if makeCuts :
                    if not passCuts : continue

                mll = TLorentzVector()
                mll = l_i + l_j
                if(abs(mll.M() - ZMASS) < offshell) :
                    mZ1 = mll.M()
                    idx_l1 = i
                    idx_l2 = j
                    findZ1 = True
                    offshell = abs(mZ1 - ZMASS)

        return offshell, findZ1, idx_l1, idx_l2

    def buildZ2Mass(self, LeptonsCollection, idx_l1, idx_l2, makeCuts):
        '''
            Build the Z2 candidate (skip leptons used for Z1).
            Returns a bool that specifies if a Z2 is present in the event
            and the indices of the two leptons used to build Z2.
        '''
        Leptons, LeptonsId, Lepts_RelIso = self.unzipLeptons(LeptonsCollection)
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

                passCuts = self.checkCuts(LeptonsCollection, i, j)
                if makeCuts :
                    if not passCuts : continue

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

    def buildZMasses(self, LeptonsCollection, makeCuts = False):
        '''
            Builds the Z candidates and checks if the event passes
            the fiducial selection.
            The function returns the (sorted) indices of the four leptons
            used to define Z1 and Z2 and a bool that specifies if the
            event passes the fiducial selection.
        '''
        Leptons, LeptonsId, Lepts_RelIso = self.unzipLeptons(LeptonsCollection)
        passFidSel = False
        passZ1 = False

        offshell, findZ1, idx_l1, idx_l2 = self.buildZ1Mass(LeptonsCollection, makeCuts)

        z1_l1 = Leptons[idx_l1]
        z1_l2 = Leptons[idx_l2]
        l1, l2 = self.buildLLPair(z1_l1, z1_l2)
        ml1l2 = TLorentzVector()
        ml1l2 = l1 + l2

        if (ml1l2.M()>MIN_MZ1 and ml1l2.M()<MAX_MZ1 and findZ1) : passZ1 = True
        if not makeCuts : passZ1 = True

        findZ2, idx_l3, idx_l4 = self.buildZ2Mass(LeptonsCollection, idx_l1, idx_l2, makeCuts)

        z_leps_idx = [idx_l1, idx_l2, idx_l3, idx_l4]

        if (passZ1 and findZ2) : passFidSel = True

        return passFidSel, z_leps_idx

    def getZCands(self, Leptons, z_leps_idx):
        '''
            Util function that returns the leptons that
            compose the Z1 and Z2 candidates.
        '''
        l1 = Leptons[z_leps_idx[0]]; l2 = Leptons[z_leps_idx[1]]
        l3 = Leptons[z_leps_idx[2]]; l4 = Leptons[z_leps_idx[3]]
        Z1_l1, Z1_l2 = self.buildLLPair(l1, l2)
        Z2_l1, Z2_l2 = self.buildLLPair(l3, l4)

        return Z1_l1, Z1_l2, Z2_l1, Z2_l2

    def getZIndex(self, LeptonsId, z_leps_idx):
        '''
            Util function that returns the IDs of the leptons that
            compose the Z1 and Z2 candidates.
        '''
        # idx_1 = LeptonsId[z_leps_idx[0]]; idx_2 = LeptonsId[z_leps_idx[1]]
        # idx_3 = LeptonsId[z_leps_idx[2]]; idx_4 = LeptonsId[z_leps_idx[3]]
        idx_1 = z_leps_idx[0]; idx_2 = z_leps_idx[1]
        idx_3 = z_leps_idx[2]; idx_4 = z_leps_idx[3]

        return idx_1, idx_2, idx_3, idx_4

    def getExtraLeps(self, LeptonsCollection, passFidSel, z_leps_idx):
        '''
            Function used to count the number of extra leptons in the event.
            Retruns the collection of extra leptons and of their IDs.
        '''
        Leptons, LeptonsId, Lepts_RelIso = self.unzipLeptons(LeptonsCollection)

        ExtraLeps = [-1]*len(Leptons)
        ExtraLepsId = [-1]*len(Leptons)

        if passFidSel:
            for idx, l_i in enumerate(Leptons):
                if((idx not in z_leps_idx) and (Lepts_RelIso[idx]<0.35)) :
                    ExtraLeps[idx] = Leptons[idx]
                    ExtraLepsId[idx] = LeptonsId[idx]

        return ExtraLeps, ExtraLepsId

    def buildZCands(self, LeptonsCollection, passFidSel, z_idx):
        '''
            Util function that combines `getZCands` and `getZIndex`
            and returns the collection of sorted leptons that make
            the two Z cands in the event and their IDs.
            For events that do not pass the fiducial selection,
            dummy collections filled with -1 are returned.
        '''
        Leptons, LeptonsId, Lepts_RelIso = self.unzipLeptons(LeptonsCollection)

        if passFidSel:
            Z1_l1, Z1_l2, Z2_l1, Z2_l2 = self.getZCands(Leptons, z_idx)
            idx_1, idx_2, idx_3, idx_4 = self.getZIndex(LeptonsId, z_idx)
            ZCands = [Z1_l1, Z1_l2, Z2_l1, Z2_l2]
            ZIdx   = [idx_1, idx_2, idx_3, idx_4]
            return ZCands, ZIdx
        else:
            return [-1]*len(Leptons), [-1]*len(Leptons)

    def countFiducialLeps(self, LeptonsCollection):
        '''
            Function that returns the number of leptons that
            pass the pT, eta and Iso criteria required by the HZZ analysis.
            The function also returns the nr of leading (pT>20) and
            subleading (pT>10) leptons that satisfy the pT, eta and Iso criteria.
        '''
        Leptons, LeptonsId, Lepts_RelIso = self.unzipLeptons(LeptonsCollection)

        nFidLeps = 0; nFidPtLead = 0; nFidPtSubLead = 0
        for l_i, idx_i, iso_i in zip(Leptons, LeptonsId, Lepts_RelIso):
            thisLep = TLorentzVector()
            thisLep.SetPtEtaPhiM(l_i.Pt(), l_i.Eta(), l_i.Phi(), l_i.M())
            if(((abs(idx_i) == 13 and thisLep.Pt() > 5.0 and abs(thisLep.Eta()) < 2.4)
                or (abs(idx_i) == 11 and thisLep.Pt() > 7.0 and abs(thisLep.Eta()) < 2.5))
                and iso_i < 0.35) :
                nFidLeps += 1
                if thisLep.Pt() > 20 : nFidPtLead += 1
                if thisLep.Pt() > 10 : nFidPtSubLead += 1 

        return nFidLeps, nFidPtLead, nFidPtSubLead

    def checkEventTopology(self, LeptonsCollection, zFid_leps_idx):
        '''
            Util function that checks topology criteria for leptons and
            returns booleans that specify if mll(OS)>4 GeV, and dR(l,l)>0.02
            for all the pairs of leptons in the event.
        '''
        Leptons, LeptonsId, Lepts_RelIso = self.unzipLeptons(LeptonsCollection)
        passedMassOS = True; passedElMuDeltaR = True; passedDeltaR = True;

        for i, l1 in enumerate(Leptons):
            for j, l2 in enumerate(Leptons):
                if j <= i : continue
                if i not in zFid_leps_idx: continue
                if j not in zFid_leps_idx: continue
                l_i, l_j = self.buildLLPair(l1, l2)
                mll = TLorentzVector()
                mll = l_i + l_j

                if(LeptonsId[i]*LeptonsId[j]<0):
                    if(mll.M()<=4):
                        passedMassOS = False
                        break

                deltaR_ll = deltaR(l_i.Eta(), l_i.Phi(), l_j.Eta(), l_j.Phi())

                if(abs(LeptonsId[i]) != abs(LeptonsId[j])):
                    if(deltaR_ll<=0.02):
                        passedElMuDeltaR = False
                        break

                if deltaR_ll <= 0.02:
                    passedDeltaR = False
                    break

        return passedMassOS, passedElMuDeltaR, passedDeltaR

    def init_collections(self):
        '''
            Util function to return empty lists for the
            Leptons information needed.
        '''
        Leptons = []
        LeptonsId = []
        LeptonsReco = []

        return Leptons, LeptonsId, LeptonsReco

    def fill_HCand_branches(self, ZCands_fidSel, ZIdx_fidSel):
        '''
            Function that fills branches for gen-level
            Higgs boson candindate.
        '''
        if ZCands_fidSel[0] == -1:
            z1mass = -1
            z2mass = -1
            zzmass = -1
            z1l1idx = -1
            z1l2idx = -1
            z2l1idx = -1
            z2l2idx = -1
            zzrapidity = -1
            zzpt = -1
            zzeta = -1
            zzphi = -1
        else:
            gen_H_cand = ZCands_fidSel[0]+ZCands_fidSel[1]+ZCands_fidSel[2]+ZCands_fidSel[3]
            zzmass = (gen_H_cand).M()
            z1l1idx = ZIdx_fidSel[0]
            z1l2idx = ZIdx_fidSel[1]
            z2l1idx = ZIdx_fidSel[2]
            z2l2idx = ZIdx_fidSel[3]
            zzrapidity = (gen_H_cand).Rapidity()
            zzpt = (gen_H_cand).Pt()
            zzeta = (gen_H_cand).Eta()
            zzphi = (gen_H_cand).Phi()
            z1mass = (ZCands_fidSel[0]+ZCands_fidSel[1]).M()
            z2mass = (ZCands_fidSel[2]+ZCands_fidSel[3]).M()

        self.out.fillBranch("FidZZ_mass", zzmass)
        self.out.fillBranch("FidZZ_Z1l1Idx", z1l1idx) #FIXME: to be sorted with standard criteria
        self.out.fillBranch("FidZZ_Z1l2Idx", z1l2idx)
        self.out.fillBranch("FidZZ_Z2l1Idx", z2l1idx)
        self.out.fillBranch("FidZZ_Z2l2Idx", z2l2idx)
        self.out.fillBranch("FidZZ_pt", zzpt)
        self.out.fillBranch("FidZZ_eta", zzeta)
        self.out.fillBranch("FidZZ_phi", zzphi)
        self.out.fillBranch("FidZZ_rapidity", zzrapidity)
        self.out.fillBranch("FidZ1_mass", z1mass)
        self.out.fillBranch("FidZ2_mass", z2mass)

    def analyze(self, event):
        '''
            Process event and return True (go to next module)
            or False (fail, go to next event).
            Stores to branches the relevant gen-level observables.
        '''

        genpart=Collection(event,"GenPart")

        dressedLeptons = [-1]*len(genpart)
        Lepts_RelIso   = [-1]*len(genpart)

        nGENHiggs = 0.0

        Leptons, LeptonsId, LeptonsReco = self.init_collections()

        for i, gp in enumerate(genpart) :
            if ((abs(gp.pdgId) == 11) or (abs(gp.pdgId) == 13) or (abs(gp.pdgId) == 15)) :
                if (not((gp.status == 1) or (abs(gp.pdgId) == 15))): continue
                mom_idx, mom_id = Mother(gp, genpart)
                if (not((mom_id==23) or (mom_id==443) or (mom_id==553) or (abs(mom_id)==24))): continue

                # Dress leptons
                # PackedGenParticles in miniAOD is GenPart.status == 1
                lep_dressed, fsr_gamma_idx = self.dressLeptons(gp, genpart)
                Leptons.append(lep_dressed)
                LeptonsId.append(gp.pdgId)

                current_lepton = lep_dressed
                genIso = self.computeGenIso(current_lepton, genpart, fsr_gamma_idx)

                if (gp.pdgId == 25): self.GenHiggsCounter(nGENHiggs)
                # TODO: Add GENZ variables (if needed for the analysis)

                Lepts_RelIso[i] = genIso
                dressedLeptons[i] = lep_dressed.Pt()

        LeptonsCollection = [Leptons, LeptonsId, Lepts_RelIso]

        if(len(Leptons)>=4):
            passFidSel_noCut, z_leps_idx = self.buildZMasses(LeptonsCollection)
            ZCands_noCut, ZIdx_noCut = self.buildZCands(LeptonsCollection, passFidSel_noCut, z_leps_idx)

        nFidLeps, nFidPtLead, nFidPtSubLead = self.countFiducialLeps(LeptonsCollection)
        passFidSel = False

        if(nFidLeps >= 4 and nFidPtLead >= 1 and nFidPtSubLead >= 2) :
            passFidSel, zFid_leps_idx = self.buildZMasses(LeptonsCollection, makeCuts = True)
            ZCands_fidSel, ZIdx_fidSel = self.buildZCands(LeptonsCollection, passFidSel, zFid_leps_idx)
            ExtraLep_fidSel, ExtraLepIdx_fidSel = self.getExtraLeps(LeptonsCollection, passFidSel, zFid_leps_idx)

            passedMassOS, passedElMuDeltaR, passedDeltaR = self.checkEventTopology(LeptonsCollection, zFid_leps_idx)
            if((passedMassOS == False) or (passedElMuDeltaR == False) or (passedDeltaR == False)): passFidSel = False

            self.fill_HCand_branches(ZCands_fidSel, ZIdx_fidSel)

            # TODO: Add GenJets
            # TODO: Add MELA

        self.out.fillBranch("nFidDressedLeps", len(dressedLeptons))
        self.out.fillBranch("FidDressedLeps_pt", dressedLeptons)
        self.out.fillBranch("FidDressedLeps_RelIso", Lepts_RelIso)
        self.out.fillBranch("passedFiducial", passFidSel)

        return True

