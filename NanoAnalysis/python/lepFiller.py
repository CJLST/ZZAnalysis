###
# -re-associate FSR to leptons according to H4l selection
# -compute FSR-corrected isolation
# -set variables for relaxed and full ID
# No filtering is applied.
###
from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.HeppyCore.utils.deltar import deltaR

def eleBDTCut(ele, era, preUL, nanoVersion=10) :
    pt = ele.pt
    fSCeta = abs(ele.eta + ele.deltaEtaSC)
    if era == 2017 or era == 2018 :
        if nanoVersion < 10:
            BDT = ele.mvaFall17V2Iso     # Using 2017 WP and training (ElectronMVAEstimatorRun2Fall17IsoV2Values) since this is the only one available in Run2 UL nanoAODs.
        else :
            BDT = ele.mvaHZZIso
        if preUL : # pre-UL WP for Run II (miniAOD branch: Run2_CutBased_BTag16)
            return (pt<=10. and     ((fSCeta<0.8                   and BDT > 0.85216885148) or \
                                     (fSCeta>=0.8 and fSCeta<1.479 and BDT > 0.82684550976) or \
                                     (fSCeta>=1.479                and BDT > 0.86937630022))) \
                    or (pt>10. and  ((fSCeta<0.8                   and BDT > 0.98248928759) or \
                                     (fSCeta>=0.8 and fSCeta<1.479 and BDT > 0.96919224579) or \
                                     (fSCeta>=1.479                and BDT > 0.79349796445)))

        else: # UL WP (miniAOD branch Run2_CutBased_UL)
            return (pt<=10. and     ((fSCeta<0.8                   and BDT > 0.9128577458) or \
                                     (fSCeta>=0.8 and fSCeta<1.479 and BDT > 0.9056792368) or \
                                     (fSCeta>=1.479                and BDT > 0.9439440575))) \
                    or (pt>10. and  ((fSCeta<0.8                   and BDT > 0.1559788054) or \
                                     (fSCeta>=0.8 and fSCeta<1.479 and BDT > 0.0273863727) or \
                                     (fSCeta>=1.479                and BDT > -0.5532483665)))

    elif era >= 2022 :
        BDT = ele.mvaHZZIso
        #2018 UL tuning (ElectronMVAEstimatorRun2Summer18ULIdIsoValues)
        return (pt<=10. and     ((fSCeta<0.8                   and BDT > 0.9044286167) or \
                                 (fSCeta>=0.8 and fSCeta<1.479 and BDT > 0.9094166886) or \
                                 (fSCeta>=1.479                and BDT > 0.9443653660))) \
                or (pt>10. and  ((fSCeta<0.8                   and BDT > 0.1968600840) or \
                                 (fSCeta>=0.8 and fSCeta<1.479 and BDT > 0.0759172100) or \
                                 (fSCeta>=1.479                and BDT > -0.5169136775)))



class lepFiller(Module):
    def __init__(self, cuts, era):
        print("***lepFiller: era:", era, flush=True)
        self.writeHistFile=False
        self.cuts = cuts
        self.passEleBDT = cuts["passEleBDT"]
        self.eleRelaxedId = cuts["eleRelaxedId"]
        self.eleRelaxedIdNoSIP = cuts["eleRelaxedIdNoSIP"]
        self.eleFullId = cuts["eleFullId"]
        self.eleFullIdNoSIP = cuts["eleFullIdNoSIP"]
        self.muRelaxedId = cuts["muRelaxedId"]
        self.muRelaxedIdNoSIP = cuts["muRelaxedIdNoSIP"]
        self.muFullId = cuts["muFullId"]
        self.muFullIdNoSIP = cuts["muFullIdNoSIP"]
        self.era = era

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree        
        
        self.out.branch("FsrPhoton_mass", "F", lenVar="nFsrPhoton") # Hack so that photon.p4() works
        self.out.branch("FsrPhoton_dROverEt2", "F", lenVar="nFsrPhoton") # Overwrite existing value

        self.out.branch("Electron_passBDT", "O", lenVar="nElectron")
        self.out.branch("Electron_ZZFullSel", "O", lenVar="nElectron")   # pass full SR selection (for electrons it is the same as FullID)
        self.out.branch("Electron_ZZFullId", "O", lenVar="nElectron")    # pass full ID selection
        self.out.branch("Electron_ZZFullSelNoSIP", "O", lenVar="nElectron") # pass full ID without SIP (base for CR SIP  method)
        self.out.branch("Electron_ZZRelaxedId", "O", lenVar="nElectron") # pass relaxed ID including SIP (base for SS and OS CRs)
        self.out.branch("Electron_ZZRelaxedIdNoSIP", "O", lenVar="nElectron") # pass relaxed ID without SIP (widest subset requirement for CRs)
        self.out.branch("Electron_fsrPhotonIdx", "I", lenVar="nElectron") # Overwrite existing value
        self.out.branch("Electron_pfRelIso03FsrCorr", "F", lenVar="nElectron")
        self.out.branch("Electron_passIso", "O", lenVar="nElectron")

        self.out.branch("Muon_ZZFullSel", "O", lenVar="nMuon")   # pass full SR selection (FullID + isolation)
        self.out.branch("Muon_ZZFullId", "O", lenVar="nMuon")    # pass full ID selection
        self.out.branch("Muon_ZZFullSelNoSIP", "O", lenVar="nMuon") # pass full ID + iso, without SIP (base for CR SIP  method)
        self.out.branch("Muon_ZZRelaxedId", "O", lenVar="nMuon")  # pass relaxed ID including SIP (base for SS and OS CRs)
        self.out.branch("Muon_ZZRelaxedIdNoSIP", "O", lenVar="nMuon")  # pass relaxed ID without SIP (widest subest requirement for CRs)
        self.out.branch("Muon_fsrPhotonIdx", "I", lenVar="nMuon") # Overwrite existing value
        self.out.branch("Muon_pfRelIso03FsrCorr", "F", lenVar="nMuon")
        self.out.branch("Muon_passIso", "O", lenVar="nMuon")

#    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
#        pass

    # Recompute lepton isolation subtracting energy of the given FSR photons
    def isoFsrCorr(self, l, selectedFSR) : 
        combRelIsoPFFSRCorr = l.pfRelIso03_all
        for f in selectedFSR :
            dR = deltaR(l.eta, l.phi, f.eta, f.phi)
            if dR >0.01 and dR < 0.3 :
                combRelIsoPFFSRCorr = max(0., l.pfRelIso03_all-f.pt/l.pt)
        return combRelIsoPFFSRCorr


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")

        # IDs (no iso)
        eleBDT = list(self.passEleBDT(e, self.era) for e in electrons)
        eleRelaxedId = list(self.eleRelaxedId(e) for e in electrons)
        muRelaxedId = list(self.muRelaxedId(m) for m in muons)
        eleRelaxedIdNoSIP = list(self.eleRelaxedIdNoSIP(e) for e in electrons)
        muRelaxedIdNoSIP = list(self.muRelaxedIdNoSIP(m) for m in muons)
        eleFullId = list(self.eleFullId(e, self.era) for e in electrons)
        muFullId = list(self.muFullId(m, self.era) for m in muons)
        eleFullIdNoSIP = list(self.eleFullIdNoSIP(e, self.era) for e in electrons)
        muFullIdNoSIP = list(self.muFullIdNoSIP(m, self.era) for m in muons)

        # Skip events that do not contain enough tight leptons (note: for CRs, this should be modified)
#FIXME: a filter here could speed up things, but must be handled properly
#        if (len(eleTight)+len(muTight)<4) : return False

        
        # Re-match FSR with Loose leptons, as default associati1on may be wrong due to different cuts, masking/stealing.
        # Store new FSR idxs, dRET2, and commpute FSR-corrected isolation
        fsrPhoton_myMuonIdx=[-1]*len(fsrPhotons)
        fsrPhoton_myElectronIdx=[-1]*len(fsrPhotons)
        fsrPhoton_mydROverEt2=[-1.]*len(fsrPhotons)
        muFsrPhotonIdx=[-1]*len(muons)
        eleFsrPhotonIdx=[-1]*len(electrons)

        muPhotondREt2 = [1e9]*len(muons)
        elePhotondREt2 = [1e9]*len(electrons)
        
        for ifsr, fsr in enumerate (fsrPhotons):
            if (fsr.pt < 2. or abs(fsr.eta) > 2.4 or fsr.relIso03 > self.cuts["fsr_Iso"] ) : continue
            dRmin = 0.5 # max dR for association
            closestMu=-1
            closestEle=-1
            for ilep, lep in enumerate(muons):
                if muRelaxedId[ilep] :
                    dR =  deltaR(lep.eta, lep.phi, fsr.eta, fsr.phi)
                    if dR < dRmin and dR > 0.001 and dR/fsr.pt/fsr.pt < self.cuts["fsr_dRET2"]:
                        dRmin = dR
                        closestMu=ilep

            for ilep, lep in enumerate(electrons):
                if eleRelaxedId[ilep] :
                    dR =  deltaR(lep.eta, lep.phi, fsr.eta, fsr.phi)
                    if dR < dRmin and dR > 0.001 and dR/fsr.pt/fsr.pt < self.cuts["fsr_dRET2"]:
                        dRmin = dR
                        closestMu=-1
                        closestEle=ilep

            if (closestMu>=0 or closestEle>=0) :
                dREt2 = dRmin/fsr.pt/fsr.pt
                fsrPhoton_mydROverEt2[ifsr] = dREt2
                if (closestMu>=0) :
                    fsrPhoton_myMuonIdx[ifsr] = closestMu
                    if dREt2<muPhotondREt2[closestMu] : #keep only the lowest-dRET2 for each lepton
                        muPhotondREt2[closestMu] = dREt2
                        muFsrPhotonIdx[closestMu] = ifsr
                if (closestEle>=0) :
                    fsrPhoton_myElectronIdx[ifsr] = closestEle
                    if dREt2<elePhotondREt2[closestEle] :
                        elePhotondREt2[closestEle] = dREt2
                        eleFsrPhotonIdx[closestEle] = ifsr


        # Recompute isolation removing all selected FSRs
        ele_isoFsrCorr = [-1.]*len(electrons)
        mu_isoFsrCorr = [-1.]*len(muons)
        ele_passIso = [False]*len(electrons)
        mu_passIso = [False]*len(muons)
        eleFullSel = [False]*len(electrons)
        eleFullSelNoSIP = [False]*len(electrons)
        muFullSel = [False]*len(muons)
        muFullSelNoSIP = [False]*len(muons)

        selectedFSR = []
        for ifsr in muFsrPhotonIdx + eleFsrPhotonIdx :
            if ifsr>=0 : selectedFSR.append(fsrPhotons[ifsr])
            
        for ilep,lep in enumerate(electrons) :
            ele_isoFsrCorr[ilep] = self.isoFsrCorr(lep, selectedFSR) 
            ele_passIso[ilep] = True # no iso cut for electrons 
            eleFullSel[ilep] = eleFullId[ilep] and ele_passIso[ilep]
            eleFullSelNoSIP[ilep] = eleFullIdNoSIP[ilep] and ele_passIso[ilep]

        for ilep,lep in enumerate(muons) :
            mu_isoFsrCorr[ilep] = self.isoFsrCorr(lep, selectedFSR)
            mu_passIso[ilep] = mu_isoFsrCorr[ilep] < self.cuts["relIso"]
            muFullSel[ilep] = muFullId[ilep] and mu_passIso[ilep]
            muFullSelNoSIP[ilep] = muFullIdNoSIP[ilep] and mu_passIso[ilep]

        fsrM = [0.]*len(fsrPhotons)
        self.out.fillBranch("FsrPhoton_mass", fsrM)
        self.out.fillBranch("FsrPhoton_dROverEt2", fsrPhoton_mydROverEt2)

        self.out.fillBranch("Electron_passBDT", eleBDT)
        self.out.fillBranch("Electron_ZZFullSel", eleFullSel)
        self.out.fillBranch("Electron_ZZFullId", eleFullId)
        self.out.fillBranch("Electron_ZZRelaxedId", eleRelaxedId)
        self.out.fillBranch("Electron_ZZRelaxedIdNoSIP", eleRelaxedIdNoSIP)
        self.out.fillBranch("Electron_ZZFullSelNoSIP", eleFullSelNoSIP)
        self.out.fillBranch("Electron_fsrPhotonIdx", eleFsrPhotonIdx)
        self.out.fillBranch("Electron_pfRelIso03FsrCorr", ele_isoFsrCorr)
        self.out.fillBranch("Electron_passIso", ele_passIso)

        self.out.fillBranch("Muon_ZZFullSel", muFullSel)
        self.out.fillBranch("Muon_ZZFullId", muFullId)
        self.out.fillBranch("Muon_ZZRelaxedId", muRelaxedId)
        self.out.fillBranch("Muon_ZZRelaxedIdNoSIP", muRelaxedIdNoSIP)
        self.out.fillBranch("Muon_ZZFullSelNoSIP", muFullSelNoSIP)
        self.out.fillBranch("Muon_fsrPhotonIdx", muFsrPhotonIdx)
        self.out.fillBranch("Muon_pfRelIso03FsrCorr", mu_isoFsrCorr)
        self.out.fillBranch("Muon_passIso", mu_passIso)

        return True
