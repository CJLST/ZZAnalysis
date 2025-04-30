from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.HeppyCore.utils.deltar import deltaR
from ZZAnalysis.NanoAnalysis.tools import getLeptons


class ZZIDStudies (Module):
    def __init__(self):
        """Add variables for specific ID studies, to be able to quickly draw ROC curves.
        """

        # ZZCand flags for muon ID studies. Each Flag will be set to true for a candidate if all of its muons pass the specified ID.
        self.muonIDs=[dict(name="ZZFullSel", sel=lambda l : l.ZZFullId and l.passIso), # Standard ZZ selection; this is used for setting default bestCandIdx
                      dict(name="ZZRelaxedIDOnly",sel=lambda l : l.pt>5 and abs(l.eta)<2.4 and (l.isGlobal or (l.isTracker and l.nStations>0))),# ZZ relaxed mu ID without dxy, dz, SIP, isolation cuts (for optimization). Note: this is looser than nanoAOD presel.
                      dict(name="ZZFullIDOnly",   sel=lambda l : l.pt>5 and abs(l.eta)<2.4 and (l.isGlobal or (l.isTracker and l.nStations>0)) and (l.isPFcand or (l.highPtId>0 and l.pt>200.))),# ZZ full ID without dxy, dz, SIP, and isolation cuts (for optimization)
                      dict(name="looseId", sel=lambda l : l.looseId),   # POG CutBasedIdLoose
                      dict(name="mediumId", sel=lambda l : l.mediumId), # POG CutBasedIdMedium
                      dict(name="mediumPromptId", sel=lambda l : l.mediumPromptId), # POG CutBasedIdMediumPrompt (=mediumId + tighter dxy, dz cuts)
                      dict(name="tightId", sel=lambda l : l.tightId), # POG CutBasedIdTight
                      dict(name="highPtId", sel=lambda l : l.highPtId>0), # >0 = POG tracker high pT; 2 = global high pT, which includes the former
                      dict(name="isPFcand", sel=lambda l : l.isPFcand),
                      dict(name="isGlobal", sel=lambda l : l.isGlobal),  # Note: this is looser than nanoAOD presel.
                      dict(name="isTracker", sel=lambda l : l.isTracker),# Note: this is looser than nanoAOD presel.
                      dict(name="isTrackerArb", sel=lambda l : l.isTracker and l.nStations>0), # Arbitrated tracker muon. Note: this is looser than nanoAOD presel.
                      dict(name="inTimeMuon", sel=lambda l : l.inTimeMuon), # 
                      ]

        # ZZCand variables storeing the worst value of a given quantity among all muons of of a candidate, for cut optimization studies.
        # Worst is intended as lowest value (as for an MVA), unless the variable's name starts with "max".
        self.muonIDVars=[dict(name="maxdxy", sel=lambda l : abs(l.dxy)),
                         dict(name="maxdz", sel=lambda l : abs(l.dz)),
                         dict(name="maxsip3d", sel=lambda l : abs(l.sip3d)),
                         dict(name="maxip3d", sel=lambda l : abs(l.ip3d)),
                         dict(name="maxpfRelIso03FsrCorr", sel=lambda l : l.pfRelIso03FsrCorr), # FSR-corrected iso, DR=0.3
                         dict(name="maxpfRelIso03_all", sel=lambda l : l.pfRelIso03_all),
                         dict(name="maxpfRelIso04_all", sel=lambda l : l.pfRelIso04_all),
                         dict(name="maxminiPFRelIso_all", sel=lambda l : l.miniPFRelIso_all), # miniIso
                         dict(name="mvaLowPt", sel=lambda l : l.mvaLowPt), # additional presel (l.looseId and l.sip3d<4. and l.dxy<0.5 and l.dz < 1) is required, cf: https://github.com/cms-sw/cmssw/blob/90f498af750cf4271c0a988fef352b0698012a40/PhysicsTools/PatAlgos/plugins/PATMuonProducer.cc#L762-L764
                         # dict(name="promptMVA", sel=lambda l : l.mvaTTH), # should add H4l preselection for consistencty with mvaLowPt; this is looser than the original recommendation (https://twiki.cern.ch/twiki/bin/viewauth/CMS/LeptonMVA). Was retrained and renamed "promptMVA" in v14
                         # dict(name="mvaMuID", sel=lambda l : l.mvaMuID), # muon MVA from 22-001. Note: Was retrained in v14; using H4l preselection for consistency, see above
                         ]

        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        for ID in self.muonIDs :
            self.out.branch("ZZCand_mu"+ID["name"], "O", lenVar="nZZCand", title=f'True if all muons of the cand pass {ID["name"]}')
        for var in self.muonIDVars :
            self.out.branch("ZZCand_mu"+var["name"], "F", lenVar="nZZCand", title=f'Worst value of {var["name"].removeprefix("max")} among all muons of the cand', limitedPrecision=16)

        self.out.branch("ZExtraMu1Idx", "S", title="Index of leading extra muon in Z events, for data/MC studies")
        self.out.branch("ZExtraMu2Idx", "S", title="Index of subleading extra muon in Z events, for data/MC studies")

        
    def analyze(self, event):

        ZZCand_passID    = [[] for il in range(len(self.muonIDs))]
        ZZCand_worstVar  = [[] for il in range(len(self.muonIDVars))]

        ZZs = Collection(event, 'ZZCand')

        # Set flags for IDs passed by all muons of candidate
        for iZZ, ZZ in enumerate(ZZs) :
            zzleps = getLeptons(ZZ, event)
            for iID, ID in enumerate(self.muonIDs) :
                passId = True
                for ilep in range(4):
                    lep = zzleps[ilep]
                    if (abs(lep.pdgId)==13 and not ID["sel"](lep)) or \
                       (abs(lep.pdgId)==11 and not lep.ZZFullSel) : # Protection in case a looser preselection for electrons was used
                        passId = False
                        continue
                ZZCand_passID[iID].append(passId)

            # Set worst value of selection variable among all candidate's muons
            for iVar, var in enumerate(self.muonIDVars):
                worstVar = 99999.
                if var["name"].startswith("max") : worstVar = -99999.
                for ilep in range(4) :
                    if abs(zzleps[ilep].pdgId)==11 : continue # Consider only muons
                    else :
                        if var["name"].startswith("max") :
                            worstVar = max(worstVar, var["sel"](zzleps[ilep]))                                    
                        else :
                            worstVar = min(worstVar, var["sel"](zzleps[ilep]))
                ZZCand_worstVar[iVar].append(worstVar)
 
        for iID, ID in enumerate(self.muonIDs) :
            self.out.fillBranch("ZZCand_mu"+ID["name"], ZZCand_passID[iID])
        for iVar, var in enumerate(self.muonIDVars) :
            self.out.fillBranch("ZZCand_mu"+var["name"], ZZCand_worstVar[iVar])


        ### Search for an additional muons in Z events, to be used for data-MC studies.
        # 
        # This is similar to the ZLCand CR that is used to determine the fake rate, but:
        # -the additional lepton is not required to pass  dxy, dz, sip requirements (only pt, eta)
        # -no QCD suppression (the mLL>4 cut on all OS pairs) cut is applied
        # -do not discard events where >1 extra lepton is present (but store the leading
        # and subleading extra leptons so that it is still possible to consider only events with exactly 1)

        ZExtraMu1Idx = ZExtraMu2Idx = -1
        ZExtraMu1Pt = ZExtraMu2Pt = -1.
        if event.bestZIdx >= 0 and event.nMuon+event.nElectron>2 :
            Zs = Collection(event, 'ZCand')
            theZ = Zs[event.bestZIdx]
            if theZ.mass > 40 and theZ.mass < 120:
                leps = Collection(event, 'Lepton')
                Zl1 = leps[theZ.l1Idx]
                Zl2 = leps[theZ.l2Idx]
                for i,aL in enumerate(leps):
                    # Search for additional muons, with ghost suppression DR cut
                    if i != theZ.l1Idx and i!= theZ.l2Idx and abs(aL.pdgId)==13 and aL.pt>5. and abs(aL.eta) < 2.4 and \
                       deltaR(aL.eta, aL.phi, Zl1.eta, Zl1.phi) > 0.02 and \
                       deltaR(aL.eta, aL.phi, Zl2.eta, Zl2.phi) > 0.02 :
                       if aL.pt > ZExtraMu1Pt :
                           ZExtraMu2Pt, ZExtraMu2Idx = ZExtraMu1Pt, ZExtraMu1Idx
                           ZExtraMu1Pt = aL.pt
                           ZExtraMu1Idx = i
                       elif aL.pt > ZExtraMu2Pt :
                           ZExtraMu2Pt = aL.pt
                           ZExtraMu2Idx = i

        self.out.fillBranch("ZExtraMu1Idx", ZExtraMu1Idx)
        self.out.fillBranch("ZExtraMu2Idx", ZExtraMu2Idx)
        
        return True
