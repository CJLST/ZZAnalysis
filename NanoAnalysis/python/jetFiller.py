### 
# Jet-Lepton cross-cleaning; pre-compute jet counts, index of leading/subleading jets, etc.
# Should be called after JES/JEC modules.
# TODO:
# - take into account jet veto maps for nCleanedJetsPt30, JetLeadingIdx, etc?
###

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.HeppyCore.utils.deltar import deltaR
from array import array

class jetFiller(Module):
    def __init__(self, year, cleanCRs=True) :
        ''' Add jet variables; apply jet cleaning; store count of cleaned jets, etc. as a reference.
        Parameters:
          cleanCRs: True = consider also failing leptons in CRs for the cleaning
        '''
        print("***jetFiller", flush=True)
        self.dRMin = 0.4 # dR between lepton (or FSR) and jet to assume overlap

        self.cleanCRs = cleanCRs
        self.year = year

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Jet_ZZId", "O", lenVar="nJet", title="jet passes ID criteria")
        self.out.branch("Jet_ZZMask", "b", lenVar="nJet", title="jet is vetoed by selected leptons or FSR photons. Bit array: bit 0=any P lep, 1=F in ZL, 2=F in SS CR, 3=F in 2P2F, 4=F in 3P1F") # Note: b=UChar_t 
        self.out.branch("Jet_ZZLepEF", "F", lenVar="nJet", title="Fraction of jet pt carried by the vetoing leptons or FSR photons", limitedPrecision=6)
        self.out.branch("Jet_ptThreshold", "F", lenVar="nJet", title="pT threshold applied to the jet for counting nCleanedJetsPt30, etc.")
        self.out.branch("JetLeadingIdx", "S", title="index of leading jet after cleaning")
        self.out.branch("JetSubleadingIdx", "S", title="index of subleading jet after cleaning")
        self.out.branch("nCleanedJetsPt30", "B", title="number of cleaned jets above 30 GeV")
        # self.out.branch("nCleanedJetsPt30_jesUp", "B", title="number of cleaned jets, up JES variation")
        # self.out.branch("nCleanedJetsPt30_jesDn", "B", title="number of cleaned jets, down JES variation")
        self.out.branch("nCleanedJetsPt30BTagged", "B", title="number of cleaned jets above 30 GeV passing the b-tagging requirement")
        self.out.branch("nCleanedJetsPt30BTagged_bTagSF", "B", title="number of cleaned jets above 30 GeV passing the b-tagging requirement with SF applied")
        # self.out.branch("nCleanedJetsPt30BTagged_bTagSFUp_correlated", "B", title="number of cleaned jets above 30 GeV passing the b-tagging requirement with SF up variation applied (variation correlated across years)")
        # self.out.branch("nCleanedJetsPt30BTagged_bTagSFUp_uncorrelated", "B", title="number of cleaned jets above 30 GeV passing the b-tagging requirement with SF up variation applied (variation uncorrelated across years)")
        # self.out.branch("nCleanedJetsPt30BTagged_bTagSFDn_correlated", "B", title="number of cleaned jets above 30 GeV passing the b-tagging requirement with SF down variation applied (variation correlated across years)")
        # self.out.branch("nCleanedJetsPt30BTagged_bTagSFDn_uncorrelated", "B", title="number of cleaned jets above 30 GeV passing the b-tagging requirement with SF down variation applied (variation uncorrelated across years)")

    def getJetAttr(self, jet, name, default=False):
        try:
            return getattr(jet, name)
        except RuntimeError as err:
            if "Unknown branch Jet_%s" % name in str(err):
                return default
            raise

    def analyze(self, event):
        jets = Collection(event, 'Jet')
        leps = Collection(event, "Lepton")
        fsrs = Collection(event, "FsrPhoton")

        ### Mark SR/CRs to which failing leptons belong. Set as a bit mask with:
        # bit 0=P (FullSel); 1=F in ZL CR; 2=F in SS CR; 3=F in 2P2F CR; 4=F in 3P1F CR
        for l in leps: l.regions = 1 if l.ZZFullSel else 0 # Set bit 0
        if self.cleanCRs and event.bestCandIdx<0 : # In SR events, no CR candidates are built
            if event.ZLCand_lepIdx>=0 : # ZL CR (no ZZ/ZLL candidate can be present in this case)
                l = leps[event.ZLCand_lepIdx]
                if not l.ZZFullSel :
                    l.regions |= 2 # Set bit 1 = F in ZL CR
            else: # ZLL CRs (may contain up to one candidate per CR)
                ZLLs = Collection(event, "ZLLCand")
                for iCR, ZLLIdx in enumerate([event.ZLLbestSSIdx, event.ZLLbest2P2FIdx, event.ZLLbest3P1FIdx]) :
                    if ZLLIdx >= 0 :
                        CRCand = ZLLs[ZLLIdx]
                        for lepIdx in [CRCand.Z2l1Idx, CRCand.Z2l2Idx] : # Z2 leptons
                            l = leps[lepIdx]
                            if not l.ZZFullSel :
                                l.regions |= 1<<(iCR+2) # Set bits 2 to 4 = F in ZLL CRs

        ZZId = [False]*event.nJet
        mask = array('B', event.nJet*[0]) # Note: UChar_t is uppercase 'B' in python array!
        jet_lepPtF = [0.]*event.nJet
        jet_ptThreshold = [30.]*event.nJet
        leadingJetIdx = -1
        subleadingJetIdx =-1
        leadingJetPt = 0.
        subleadingJetPt = 0.
        nCleanedJetsPt30 = 0
        # nCleanedJetsPt30_jesDn = 0
        # nCleanedJetsPt30_jesUp = 0
        nCleanedJetsPt30BTagged = 0
        nCleanedJetsPt30BTagged_bTagSF = 0
        # nCleanedJetsPt30BTagged_bTagSFUp_correlated = 0
        # nCleanedJetsPt30BTagged_bTagSFUp_uncorrelated = 0
        # nCleanedJetsPt30BTagged_bTagSFDn_correlated = 0
        # nCleanedJetsPt30BTagged_bTagSFDn_uncorrelated = 0

        for ij, jet in enumerate(jets) :
            # jet ID: tight ID and tightLepVeto ID, cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13p6TeV#nanoAOD_Flags
            # Note: jetId is a bitmask; in v15 samples the value is jetId==6, but in v9 samples additional bits were set,
            # so we are using a mask as recommended.
            # The cut at |eta|<=5 is explicitly added as it is not implemented consistently in the different jetID recipes.
            if (jet.jetId & 6) == 6 and abs(jet.eta) <= 5.0 :
                ZZId[ij] = True
                
            # Mitigation for HF spikes, cf https://indico.cern.ch/event/1624984/#30-jet-performance-calibration, slide 14.
            # Summary of applied pT cuts following JME recommendations:
            # Jets with abs(eta) < 2.5 -> pt_cut = 30
            # Jets with 2.5 <= abs(eta) < 3.0 -> pt_cut = 50 in 2022 to 2024 (*), Ref: https://gitlab.cern.ch/cms-jetmet/coordination/coordination/-/issues/113
            # Jets with abs(eta) >= 3.0 in 2022/2023 -> pt_cut = 50
            # Jets with abs(eta) >= 3.0 in other years -> pt_cut = 30
            # (*) For 2025 and 26 it should be rechecked, cf. slides linked above.
            if self.year >= 2022 and self.year <= 2024 :
                if 2.5 <= abs(jet.eta) < 3.0 or (abs(jet.eta) >= 3.0 and self.year in [2022, 2023]):
                    jet_ptThreshold[ij] = 50.

            # Mask jet if it overlaps with a P (FullSel) lepton, or with a FSR photon attached to a P lepton; also,
            # for CR events, mask jet if it overlaps with a F lepton belonging to a CR candidate, or his FSR photon.
            # The bitmask value can be used to keep track of which CR originated the veto, in case >1 CR candidates
            # are present in the same event, or as a bool for "inclusive" cleaning (i.e. mask with all F leptons of
            # all CR candidates in the event considered together).
            # Also, compute the fraction of the pT of leptons+FSR in the code; for this value, which is no longer used
            # for masking, the F leptons of all CR candidates in the event are considered together, for simplicity.
            for lep in leps :
                if lep.regions==0 : continue # does not pass ZZFullSel, and is not a F lepton in CRs: skip
                if deltaR(lep.eta, lep.phi, jet.eta, jet.phi) < self.dRMin :
                    jet_lepPtF[ij] += lep.pt
                    mask[ij] |= lep.regions
                # Consider FSR photon attached to the lepton; it may overlap with jet even if the lepton doesn't
                if lep.fsrPhotonIdx >=0 : 
                    fsr = fsrs[lep.fsrPhotonIdx]
                    if deltaR(fsr.eta, fsr.phi, jet.eta, jet.phi) < self.dRMin :
                        jet_lepPtF[ij] += fsr.pt
                        mask[ij] |= lep.regions # Use the regions of the corresponding lepton

            jet_lepPtF[ij] /= jet.pt

            # Store pre-computed values of nJets and leading/subleading jet; using "inclusive" cleaning in CRs.
            if (not mask[ij]) and ZZId[ij] and jet.pt > jet_ptThreshold[ij] :
                nCleanedJetsPt30 += 1
                #FIXME: add jesUp, jesDn
                isBtagged = self.getJetAttr(jet, "isBtagged", False)
                isBtaggedSF = self.getJetAttr(jet, "isBtaggedwithSF", isBtagged)
                # isBtaggedSFUp_correlated = self.getJetAttr(jet, "isBtaggedwithSF_up_correlated", isBtaggedSF)
                # isBtaggedSFUp_uncorrelated = self.getJetAttr(jet, "isBtaggedwithSF_up_uncorrelated", isBtaggedSF)
                # isBtaggedSFDn_correlated = self.getJetAttr(jet, "isBtaggedwithSF_down_correlated", isBtaggedSF)
                # isBtaggedSFDn_uncorrelated = self.getJetAttr(jet, "isBtaggedwithSF_down_uncorrelated", isBtaggedSF)

                if isBtagged:
                    nCleanedJetsPt30BTagged += 1
                if isBtaggedSF:
                    nCleanedJetsPt30BTagged_bTagSF += 1
                # if isBtaggedSFUp_correlated:
                #     nCleanedJetsPt30BTagged_bTagSFUp_correlated += 1
                # if isBtaggedSFUp_uncorrelated:
                #     nCleanedJetsPt30BTagged_bTagSFUp_uncorrelated += 1
                # if isBtaggedSFDn_correlated:
                #     nCleanedJetsPt30BTagged_bTagSFDn_correlated += 1
                # if isBtaggedSFDn_uncorrelated:
                #     nCleanedJetsPt30BTagged_bTagSFDn_uncorrelated += 1

                # IDX of Leading/subleading jets passing all selections (including pT).
                # Note: we cannot rely on the fact that the jet collection is sorted by pt since JES can change this.  
                if jet.pt > leadingJetPt:
                    subleadingJetPt = leadingJetPt
                    subleadingJetIdx = leadingJetIdx
                    leadingJetPt = jet.pt
                    leadingJetIdx = ij
                elif jet.pt > subleadingJetPt :
                    subleadingJetPt = jet.pt
                    subleadingJetIdx = ij
        
        self.out.fillBranch("Jet_ZZId", ZZId)
        self.out.fillBranch("Jet_ZZMask", mask)
        self.out.fillBranch("Jet_ZZLepEF", jet_lepPtF)
        self.out.fillBranch("Jet_ptThreshold", jet_ptThreshold)
        self.out.fillBranch("JetLeadingIdx", leadingJetIdx)
        self.out.fillBranch("JetSubleadingIdx", subleadingJetIdx)
        self.out.fillBranch("nCleanedJetsPt30", nCleanedJetsPt30)
        # self.out.fillBranch("nCleanedJetsPt30_jesUp", nCleanedJetsPt30_jesUp)
        # self.out.fillBranch("nCleanedJetsPt30_jesDn", nCleanedJetsPt30_jesDn)
        self.out.fillBranch("nCleanedJetsPt30BTagged", nCleanedJetsPt30BTagged)
        self.out.fillBranch("nCleanedJetsPt30BTagged_bTagSF", nCleanedJetsPt30BTagged_bTagSF)
        # self.out.fillBranch("nCleanedJetsPt30BTagged_bTagSFUp_correlated", nCleanedJetsPt30BTagged_bTagSFUp_correlated)
        # self.out.fillBranch("nCleanedJetsPt30BTagged_bTagSFUp_uncorrelated", nCleanedJetsPt30BTagged_bTagSFUp_uncorrelated)
        # self.out.fillBranch("nCleanedJetsPt30BTagged_bTagSFDn_correlated", nCleanedJetsPt30BTagged_bTagSFDn_correlated)
        # self.out.fillBranch("nCleanedJetsPt30BTagged_bTagSFDn_uncorrelated", nCleanedJetsPt30BTagged_bTagSFDn_uncorrelated)

        return True

