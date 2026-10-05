### 
# Jet-Lepton cross-cleaning. Should be called after JES/JEC modules.
# Cleaning is based on the fraction of jet pt that is carried by all leptons passing full sel
# and their FSR photons in dR < 0.4. The jet is masked if the fraction is > 0.5.
# TODO:
# - take into account jet veto maps for nCleanedJetsPt30, JetLeadingIdx, etc?
###

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.HeppyCore.utils.deltar import deltaR

class jetFiller(Module):
    def __init__(self, year):
        print("***jetFiller", flush=True)
        self.dRMin = 0.4 # dR between lepton (or FSR) and jet to assume overlap
        self.EFthreshold = 0.5 
        #############################################################################
        # NOTE: This is the threshold of leptons pt over jet pt to veto the jet. It #
        # is also used in RecoProbFiller to reapply the veto to the varied jets, and#
        # should be kept consistent with this one.                                  #
        #############################################################################

        self.year = year

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Jet_ZZMask", "O", lenVar="nJet", title="jet is vetoed by selected leptons or FSR photons")
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
        electrons = Collection(event, "Electron")
        fsrs = Collection(event, "FsrPhoton")
        muons = Collection(event, "Muon")
        leps = list(muons)+ list(electrons)
        nlep=len(leps)

        ## Jet-lepton cleaning with best candidate's leptons and FSR..
        mask = [False]*event.nJet
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

        # According to the analysis recipe at https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsZZ4lRunIILegacy#Jets, "[jets] must be cleaned with a DeltaR>0.4 cut wrt all tight leptons in the event passing the SIP and isolation cut computed after FSR correction, as well as with all FSR collected photons attached to these leptons."
        # Note: the current implementation on miniAODs (https://github.com/CJLST/ZZAnalysis/blob/Run2UL_22_nano/AnalysisStep/plugins/JetsWithLeptonsRemover.cc) probably does something different than this. To be reviewed.
        for ij, jet in enumerate(jets) :
            if self.year >=2022 : #Larger pT threshold at high eta, only for Run3
                if 2.5 <= abs(jet.eta) < 3.0 or (abs(jet.eta) >= 3.0 and self.year in [2022, 2023]):
                    jet_ptThreshold[ij] = 50.

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
                # Consider only jets passing  tight ID and tightLepVeto ID, cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13p6TeV#nanoAOD_Flags
                # FIXME: To be checked/updated for Run2 
                if jet.jetId == 6 :
                    #############################################################################
                    # NOTE: RecoProbFiller also uses the jet ID requirement, and should be kept #
                    #consistent with this one.                                                  #
                    #############################################################################
                    #Summary of applied pT cuts following JME recommendations:
                    #Jets with abs(eta) < 2.5 -> pt_cut = 30
                    #Jets with 2.5 <= abs(eta) < 3.0 -> pt_cut = 50 , Ref: https://gitlab.cern.ch/cms-jetmet/coordination/coordination/-/issues/113
                    #Jets with abs(eta) >= 3.0 in 2022/2023 -> pt_cut = 50
                    #Jets with abs(eta) >= 3.0 in other years -> pt_cut = 30
                    #See slide 14: https://indico.cern.ch/event/1615783/contributions/6811120/attachments/3186812/5672346/20251204_JetMET_PerformanceRun3.pdf

                    if jet.pt > jet_ptThreshold[ij]:
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

