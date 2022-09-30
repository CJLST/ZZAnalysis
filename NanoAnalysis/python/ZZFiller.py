###
# Build Z and ZZ candidates.
# ZZ candidates are stored as a nanoAOD collection, so that several candidates per event can be saved (for eg. ID optimization studies).
# The index of the best ZZ candidate in each event is stored as ZZBestCand_Idx
# Filter: only events with at least one ZZ candidate are kept
###
from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import *

from functools import cmp_to_key
from ctypes import CDLL, c_float, c_int, c_bool
from JHUGenMELA.MELA.mela import Mela, SimpleParticle_t, SimpleParticleCollection_t, TVar

class StoreOption:
    BestCandOnly, AllCands, AllWithRelaxedMuId = range(0,3)


class ZZFiller(Module):

    def __init__(self, runMELA, bestCandByMELA, isMC, year):
        self.writeHistFile = False
        self.isMC = isMC
        self.year = year
        self.runMELA = runMELA or bestCandByMELA
        if bestCandByMELA :
            self.bestCandCmp = self.bestCandByDbkgKin
        else:
            self.bestCandCmp = self.bestCandByZ1Z2
        self.DEBUG = False
        self.ZmassValue = 91.1876;

        # Which candidates should be stored: BestCandOnly, AllCands, AllWithRelaxedMuId (for muon ID studies)
        self.candsToStore = StoreOption.BestCandOnly

        # Lepton requirements for the best candidate: full ID + iso
        self.bestCandPresel =  (lambda l : l.ZZFullId and l.passIso)

        # Pre-selection of leptons to build Z and ZZ candidates; normally it is the full ID + iso,
        # but can be relaxed to store loose candidates for ID studies (see below)
        self.leptonPresel = (lambda l : l.ZZFullId and l.passIso)
        self.muonIDs=[]
        self.muonIDVars=[]

        # Use relaxed muon preselection for ID studies: fully relax muon ID. Electron ID is unchanged.
        if self.candsToStore == StoreOption.AllWithRelaxedMuId :
            self.leptonPresel = (lambda l : (abs(l.pdgId)==13 and l.pt>5 and abs(l.eta) < 2.4) or (abs(l.pdgId)==11 and l.ZZFullId))

            # Add flags for muon ID studies. Each Flag will be set to true for a candidate if all of its muons pass the specified ID.
            self.muonIDs=[dict(name="ZZFullSel", sel=lambda l : l.ZZFullId and l.passIso), # Standard ZZ selection; this is used for setting default bestCandIdx
                          dict(name="ZZRelaxedIDNoSIP",sel=lambda l : l.pt>5 and abs(l.eta)<2.4 and (l.isGlobal or (l.isTracker and l.nStations>0))),# ZZ relaxed mu ID without dxy, dz, SIP cuts (for optimization). Note: this is looser than nanoAOD presel.
                          dict(name="ZZFullIDNoSIP",   sel=lambda l : l.pt>5 and abs(l.eta)<2.4 and (l.isGlobal or (l.isTracker and l.nStations>0)) and (l.isPFcand or (l.highPtId>0 and l.pt>200.))),# ZZ full ID without dxy, dz, SIP, and isolation cuts (for optimization)                      
                          dict(name="looseId", sel=lambda l : l.looseId),   # POG CutBasedIdLoose
                          dict(name="mediumId", sel=lambda l : l.mediumId), # POG CutBasedIdMedium
                          dict(name="mediumPromptId", sel=lambda l : l.mediumPromptId), # POG CutBasedIdMediumPrompt (=mediumId + tighter dxy, dz cuts)
                          dict(name="tightId", sel=lambda l : l.tightId), # POG CutBasedIdTight
                          dict(name="highPtId", sel=lambda l : l.highPtId>0), # >0 = POG tracker high pT; 2 = global high pT, which includes the former
                          dict(name="isPFcand", sel=lambda l : l.isPFcand),
                          dict(name="isGlobal", sel=lambda l : l.isGlobal),  # Note: this is looser than nanoAOD presel.
                          dict(name="isTracker", sel=lambda l : l.isTracker),# Note: this is looser than nanoAOD presel.
                          dict(name="isTrackerArb", sel=lambda l : l.isTracker and l.nStations>0), # Arbitrated tracker muon. Note: this is looser than nanoAOD presel.
                          ]

            # Add variable to store the worst value of a given quantity among the 4 leptons of a candidate, for optimization studies.
            # Worst is intended as lowest value (as for an MVA), unless the variable's name starts with "max".
            self.muonIDVars=[dict(name="maxsip3d", sel=lambda l : l.sip3d if (abs(l.dxy)<0.5 and abs(l.dz) < 1) else 999.), # dxy, dz cuts included with SIP
                             dict(name="maxpfRelIso03FsrCorr", sel=lambda l : l.pfRelIso03FsrCorr),
                             dict(name="mvaLowPt", sel=lambda l : l.mvaLowPt if (l.looseId and l.sip3d<4. and l.dxy<0.5 and l.dz < 1) else -2.), # additional presel required, cf: https://cmssdt.cern.ch/dxr/CMSSW/source/PhysicsTools/PatAlgos/plugins/PATMuonProducer.cc#1027-1046
                             ]


        # Data-MC SFs. 
        # NanoAODTools provides a module based on LeptonEfficiencyCorrector.cc, but that does not seem to be flexible enough for us:
        # https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/common/lepSFProducer.py
        self.lib = CDLL('libZZAnalysisAnalysisStep.so') # used for also for MELA c-constants
        self.lepSFHelper = self.lib.get_LeptonSFHelper() # Note for 2016 UL samples: requires passing bool preVFP
        self.lib.LeptonSFHelper_getSF.restype = c_float
        self.lib.LeptonSFHelper_getSFError.restype = c_float

        if self.runMELA :
            self.mela = Mela(13, 125, TVar.ERROR)
            self.mela.setCandidateDecayMode(TVar.CandidateDecay_ZZ)
            self.lib.D_bkg_kin.restype = c_float        


        # Example of adding control histograms (requires self.writeHistFile = True)
        # def beginJob(self,histFile=None, histDirName=None):
        #    Module.beginJob(self, histFile, histDirName+"_ZZFiller")
        #    self.histFile=None # Hack to prevent histFile to be closed before other modules write their histograms
        #    self.h_ZZMass = ROOT.TH1F('ZZMass','ZZMass',130,70,200)
        #    self.addObject(self.h_ZZMass)


    def endJob(self):
        self.lib.del_LeptonSFHelper(self.lepSFHelper)


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree


        self.out.branch("nZZCand", "I")
        self.out.branch("ZZCand_mass", "F", lenVar="nZZCand")
        self.out.branch("ZZCand_massPreFSR", "F", lenVar="nZZCand")
        self.out.branch("ZZCand_Z1mass", "F", lenVar="nZZCand")
        self.out.branch("ZZCand_Z1flav", "I", lenVar="nZZCand")
        self.out.branch("ZZCand_Z2mass", "F", lenVar="nZZCand")
        self.out.branch("ZZCand_Z2flav", "I", lenVar="nZZCand")
        self.out.branch("ZZCand_KD", "F", lenVar="nZZCand")
        self.out.branch("ZZCand_Z2sumpt", "F", lenVar="nZZCand")
        # Note: lepton indices are numbered for leps=list(muons)+list(electrons) and run up to nlep=len(leps);
        # no special ordering of l1, l2 is applied
        self.out.branch("ZZCand_Z1l1Idx", "I", lenVar="nZZCand") 
        self.out.branch("ZZCand_Z1l2Idx", "I", lenVar="nZZCand")
        self.out.branch("ZZCand_Z2l1Idx", "I", lenVar="nZZCand")
        self.out.branch("ZZCand_Z2l2Idx", "I", lenVar="nZZCand")
        for ID in self.muonIDs :
            self.out.branch("ZZCand_mu"+ID["name"], "O", lenVar="nZZCand")
        for var in self.muonIDVars :
            self.out.branch("ZZCand_mu"+var["name"], "F", lenVar="nZZCand")
        if self.isMC :
            self.out.branch("ZZCand_dataMCWeight", "F", lenVar="nZZCand")

        self.out.branch("ZZBestCand_Idx", "I")


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        eventId='{}:{}:{}'.format(event.run,event.luminosityBlock,event.event)
        if self.DEBUG : print ('Event '+eventId)
 
        # Collections
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")
        leps = list(muons)+ list(electrons)
        nlep=len(leps)

        ### Z combinatorial over selected leps (after FSR-corrected ISO cut for muons)
        Zs = []
        for i,l1 in enumerate(leps):            
            if self.leptonPresel(l1) :
                for j in range(i+1,nlep):
                    l2 = leps[j]
                    if self.leptonPresel(l2) and l1.pdgId == -l2.pdgId : #OS,SF
                        aZ = self.ZCand(i, j, leps, fsrPhotons)
                        zmass = aZ.M
                        if self.DEBUG: print('Z={:.4g} pt1={:.3g} pt2={:.3g}'.format(zmass, l1.pt, l2.pt), l1.fsrPhotonIdx,  l2.fsrPhotonIdx)
                        if (zmass>12. and zmass<120.):
                            Zs.append(aZ)

# At some point, we may want to store Z candidates in the event, in a dedicated ZFiller module; like:
#                         Z_d1[iZ]=aZ.l1Idx
#                         Z_d2[iZ]=aZ.l2Idx
#                         Z_fsr1Idx[iZ]=aZ.fsr1Idx
#                         Z_fsr2Idx[iZ]=aZ.fsr2Idx
#                         Z_pt[iZ]=aZ.p4.pt()
#                         Z_eta[iZ]=aZ.p4.eta()
#                         Z_phi[iZ]=aZ.p4.phi()
#                         Z_mass[iZ]=aZ.p4.M()
#                         Z_finalState[iZ]=aZ.finalState()


        ### Build all ZZ combinations passing the ZZ selection
        ZZs = []
        bestCandIdx = -1
        if len(Zs) >= 2:
            for iZ,aZ in enumerate(Zs):
                for jZ in range(iZ+1, len(Zs)):
                    # check that these Zs are mutually exclusive (not sharing the same lepton) 
                    if Zs[iZ].l1Idx==Zs[jZ].l1Idx or Zs[iZ].l2Idx==Zs[jZ].l2Idx or Zs[iZ].l2Idx==Zs[jZ].l1Idx or Zs[iZ].l2Idx==Zs[jZ].l2Idx: continue
                    
                    # set Z1 and Z2
                    Z1Idx, Z2Idx = jZ, iZ
                    if abs(Zs[iZ].M-self.ZmassValue) < abs(Zs[jZ].M-self.ZmassValue):
                        Z1Idx, Z2Idx = Z2Idx, Z1Idx

                    Z1=Zs[Z1Idx]
                    Z2=Zs[Z2Idx]

                    # Z1 mass cut
                    if Z1.M <= 40. : continue

                    zzleps = [Z1.l1, Z1.l2, Z2.l1, Z2.l2]
                    lepPts = []
                    # QCD suppression on all OS pairs, regardelss of flavour
                    passQCD = True    # QCD suppression on all OS pairs, regardelss of flavour
                    passDeltaR = True # DeltaR>0.02 cut among all leptons to protect against split tracks
                    for k in range(4):
                        lepPts.append(zzleps[k].pt)
                        for l in range (k+1,4):
                            if zzleps[k].pdgId*zzleps[l].pdgId < 0 and (zzleps[k].p4()+zzleps[l].p4()).M()<=4.:
                                passQCD = False
                                break
                            if deltaR(zzleps[k].eta, zzleps[k].phi, zzleps[l].eta, zzleps[l].phi) <= 0.02 :
                                passDeltaR = False
                                break

                    if not (passQCD and passDeltaR) : continue

                    # trigger acceptance cuts (20,10 GeV)
                    lepPts.sort()
                    if not (lepPts[3]>20. and lepPts[2]>10.) : continue

                    #"Smart cut" on alternate pairings for same-sign candidates
                    if abs(Z1.l1.pdgId) == abs(Z2.l1.pdgId):
                        mZa, mZb = 0., 0.
                        if Z1.l1.pdgId == -Z2.l1.pdgId:
                            mZa=(Z1.l1DressedP4+Z2.l1DressedP4).M()
                            mZb=(Z1.l2DressedP4+Z2.l2DressedP4).M()
                        elif Z1.l1.pdgId == -Z2.l2.pdgId:
                            mZa=(Z1.l1DressedP4+Z2.l2DressedP4).M()
                            mZb=(Z1.l2DressedP4+Z2.l1DressedP4).M()
                        if (abs(mZa-self.ZmassValue)>abs(mZb-self.ZmassValue)) : mZa, mZb = mZb, mZa
                        if (abs(mZa-self.ZmassValue)<abs(Z1.M-self.ZmassValue)) and mZb < 12.: continue

                    #Compute D_bkg^kin
                    p_GG_SIG_ghg2_1_ghz1_1_JHUGen = 0.
                    p_QQB_BKG_MCFM = 1.
                    if self.runMELA:
                        daughters = SimpleParticleCollection_t()
                        daughters.push_back(SimpleParticle_t(Z1.l1.pdgId, Z1.l1DressedP4))
                        daughters.push_back(SimpleParticle_t(Z1.l2.pdgId, Z1.l2DressedP4))
                        daughters.push_back(SimpleParticle_t(Z2.l1.pdgId, Z2.l1DressedP4))
                        daughters.push_back(SimpleParticle_t(Z2.l2.pdgId, Z2.l2DressedP4))
                        self.mela.setInputEvent(daughters, 0, 0, 0)
                        self.mela.setProcess(TVar.HSMHiggs, TVar.JHUGen, TVar.ZZGG)
                        p_GG_SIG_ghg2_1_ghz1_1_JHUGen = self.mela.computeP(True)
                        self.mela.setProcess(TVar.bkgZZ, TVar.MCFM, TVar.ZZQQB)
                        p_QQB_BKG_MCFM = self.mela.computeP(True)
                        self.mela.resetInputEvent()
                    ZZ = self.ZZCand(Z1, Z2, p_GG_SIG_ghg2_1_ghz1_1_JHUGen, p_QQB_BKG_MCFM)

                    # If criteria for best cand are statisfied, check if this is the best cand and set bestCandIdx
                    ZZ.passBestCandPresel = True
                    for ilep in range(4):
                        if not self.bestCandPresel(zzleps[ilep]) : ZZ.passBestCandPresel = False
                    if ZZ.passBestCandPresel :
                        if bestCandIdx<0 or self.bestCandCmp(ZZ,ZZs[bestCandIdx]) < 0:
                            bestCandIdx = len(ZZs)

                    # Set flags for IDs passed by all leptons of candidate (muon only, for the time being)
                    passId = [True]*len(self.muonIDs)
                    for iID, ID in enumerate(self.muonIDs) :
                        for ilep in range(4):
                            if abs(zzleps[ilep].pdgId)==13 and not ID["sel"](zzleps[ilep]) :
                                passId[iID] = False
                                continue
                    ZZ.passId = passId

                    # Set worst value of specified selection variables (muon only, for the time being)
                    worstVar = [99.]*len(self.muonIDVars)
                    for iVar, var in enumerate(self.muonIDVars):
                        if var["name"].startswith("max") : worstVar[iVar] = -1.
                        for iilep in range(4) :
                            if abs(zzleps[iilep].pdgId)==11 : continue
                            else :
                                if var["name"].startswith("max") :
                                    worstVar[iVar] = max(worstVar[iVar], var["sel"](zzleps[iilep]))                                    
                                else :
                                    worstVar[iVar] = min(worstVar[iVar], var["sel"](zzleps[iilep]))
                    ZZ.worstVar = worstVar


                    if self.DEBUG: print("ZZ:", len(ZZs),  (Z1.p4+Z2.p4).M(), ZZ.Z1.M, ZZ.Z2.M, ZZ.Z2.sumpt(), ZZ.finalState(), p_GG_SIG_ghg2_1_ghz1_1_JHUGen, p_QQB_BKG_MCFM, ZZ.KD, ZZ.passBestCandPresel)
                    ZZs.append(ZZ)


            if len(ZZs) == 0 : return False
            
            if self.DEBUG: print("bestCand:", bestCandIdx)

            if self.candsToStore == StoreOption.BestCandOnly : # keep only the best cand, as single element of the collection
                ZZs = [ZZs[bestCandIdx]]
                bestCandIdx = 0

            ZZCand_mass = [0.]*len(ZZs)
            ZZCand_massPreFSR = [0.]*len(ZZs)
            ZZCand_Z1mass = [0.]*len(ZZs)
            ZZCand_Z1flav = [0.]*len(ZZs)
            ZZCand_Z2mass = [0.]*len(ZZs)
            ZZCand_Z2flav = [0.]*len(ZZs)
            ZZCand_Z1l1Idx = [-1]*len(ZZs)
            ZZCand_Z1l2Idx = [-1]*len(ZZs)
            ZZCand_Z2l1Idx = [-1]*len(ZZs)
            ZZCand_Z2l2Idx = [-1]*len(ZZs)
            ZZCand_KD = [0.]*len(ZZs)
            ZZCand_Z2sumpt = [0.]*len(ZZs)
            ZZCand_wDataMC = [0.]*len(ZZs)
            ZZCand_passID    = [[] for il in range(len(self.muonIDs))]
            ZZCand_worstVar  = [[] for il in range(len(self.muonIDVars))]
            
            for iZZ, ZZ in enumerate(ZZs) :
                ZZCand_mass[iZZ] = ZZ.p4.M()
                ZZCand_massPreFSR[iZZ] = ZZ.massPreFSR()
                ZZCand_Z1mass[iZZ] = ZZ.Z1.M
                ZZCand_Z1flav[iZZ] = ZZ.Z1.finalState()
                ZZCand_Z2mass[iZZ] = ZZ.Z2.M
                ZZCand_Z2flav[iZZ] = ZZ.Z2.finalState()
                ZZCand_Z1l1Idx[iZZ] = ZZ.Z1.l1Idx
                ZZCand_Z1l2Idx[iZZ] = ZZ.Z1.l2Idx
                ZZCand_Z2l1Idx[iZZ] = ZZ.Z2.l1Idx
                ZZCand_Z2l2Idx[iZZ] = ZZ.Z2.l2Idx
                ZZCand_KD[iZZ] = ZZ.KD
                ZZCand_Z2sumpt[iZZ] = ZZ.Z2.sumpt()
                if self.isMC: ZZCand_wDataMC[iZZ] =  self.getDataMCWeight(ZZ.leps())
                for iID, ID in enumerate(self.muonIDs) :
                    ZZCand_passID[iID].append(ZZ.passId[iID])
                for iVar, var in enumerate(self.muonIDVars) :
                    ZZCand_worstVar[iVar].append(ZZ.worstVar[iVar])
                

            self.out.fillBranch("nZZCand", len(ZZs))
            self.out.fillBranch("ZZCand_mass", ZZCand_mass)
            self.out.fillBranch("ZZCand_massPreFSR", ZZCand_massPreFSR)
            self.out.fillBranch("ZZCand_Z1mass", ZZCand_Z1mass)
            self.out.fillBranch("ZZCand_Z1flav", ZZCand_Z1flav)
            self.out.fillBranch("ZZCand_Z2mass", ZZCand_Z2mass)
            self.out.fillBranch("ZZCand_Z2flav", ZZCand_Z2flav)
            self.out.fillBranch("ZZCand_KD", ZZCand_KD)
            self.out.fillBranch("ZZCand_Z2sumpt", ZZCand_Z2sumpt)
            self.out.fillBranch("ZZCand_Z1l1Idx", ZZCand_Z1l1Idx)
            self.out.fillBranch("ZZCand_Z1l2Idx", ZZCand_Z1l2Idx)
            self.out.fillBranch("ZZCand_Z2l1Idx", ZZCand_Z2l1Idx)
            self.out.fillBranch("ZZCand_Z2l2Idx", ZZCand_Z2l2Idx)
            for iID, ID in enumerate(self.muonIDs) :
                self.out.fillBranch("ZZCand_mu"+ID["name"], ZZCand_passID[iID])
            for iVar, var in enumerate(self.muonIDVars) :
                self.out.fillBranch("ZZCand_mu"+var["name"], ZZCand_worstVar[iVar])

            if self.isMC:
                self.out.fillBranch("ZZCand_dataMCWeight", ZZCand_wDataMC)

            self.out.fillBranch("ZZBestCand_Idx", bestCandIdx)


            ### Fill control plot (example)
            # self.h_ZZMass.Fill(ZZCand_mass[bestCandIdx])

            return True

        return False # No candidate found, do not write the event.


    # Temporary class to store information on a ZZ candidate.
    # NOTE: We may need to move Zs into the Event as persistent objects, and to be re-used for CRs. They would be built by a separate module in that case.
    class ZCand: 
        def __init__(self, l1Idx, l2Idx, leps, fsrPhotons):
            # FIXME: we may want to set a default order for i, j (eg 1=-, 2=+)
            self.l1Idx = l1Idx
            self.l2Idx = l2Idx
            self.l1 = leps[l1Idx]
            self.l2 = leps[l2Idx]
            self.fsr1Idx = self.l1.fsrPhotonIdx
            self.fsr2Idx = self.l2.fsrPhotonIdx
    
            self.l1DressedP4 = self.l1.p4()
            self.l2DressedP4 = self.l2.p4()
            if self.fsr1Idx>=0 : self.l1DressedP4 += fsrPhotons[self.fsr1Idx].p4()
            if self.fsr2Idx>=0 : self.l2DressedP4 += fsrPhotons[self.fsr2Idx].p4()
    
            self.p4 = self.l1DressedP4 + self.l2DressedP4
    
            self.M = self.p4.M() # cache the mass as it is used often
    
        def sumpt(self) : # sum of lepton pTs, used to sort candidates
            return self.l1.pt+self.l2.pt
    
        def finalState(self) :
            return self.l1.pdgId*self.l2.pdgId
    
    
    # Temporary class to store information on a ZZ candidate.
    class ZZCand:
        def __init__(self, Z1, Z2, p_GG_SIG_ghg2_1_ghz1_1_JHUGen=0., p_QQB_BKG_MCFM=1.):
            self.Z1 = Z1
            self.Z2 = Z2
            self.p4 = Z1.p4+Z2.p4
            self.p_GG_SIG_ghg2_1_ghz1_1_JHUGen = p_GG_SIG_ghg2_1_ghz1_1_JHUGen
            self.p_QQB_BKG_MCFM = p_QQB_BKG_MCFM
            self.KD = p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(p_GG_SIG_ghg2_1_ghz1_1_JHUGen+p_QQB_BKG_MCFM) # without c-constants, for candidate sorting

        def finalState(self) :
            return self.Z1.finalState()*self.Z2.finalState()

        def massPreFSR(self) :
            return (self.Z1.l1.p4()+self.Z1.l2.p4()+self.Z2.l1.p4()+self.Z2.l2.p4()).M()

        def leps(self) :
            return([self.Z1.l1, self.Z1.l2, self.Z2.l1, self.Z2.l2])


    ### Compute lepton efficiency scale factor
    def getDataMCWeight(self, leps) :
       dataMCWeight = 1.
       for lep in leps:           
           myLepID = abs(lep.pdgId)
           mySCeta = lep.eta
           isCrack = False
           if myLepID==11 :
               mySCeta = lep.eta + lep.deltaEtaSC
               # FIXME: isGap() is not available in nanoAODs, and cannot be recomputed easily based on eta, phi. We thus use the non-gap SFs for all electrons.

           # Deal with very rare cases when SCeta is out of 2.5 bounds
           mySCeta = min(mySCeta,2.49)
           mySCeta = max(mySCeta,-2.49)

           SF = self.lib.LeptonSFHelper_getSF(self.lepSFHelper, c_int(self.year), c_int(myLepID), c_float(lep.pt), c_float(lep.eta), c_float(mySCeta), c_bool(isCrack))
#           SF_Unc = lepSFHelper.getSFError(year,myLepID,lep.pt,lep.eta, mySCeta, isCrack)
           dataMCWeight *= SF

       return dataMCWeight


    ### Comparators to select the best candidate in the event. Return -1 if a is better than b, +1 otherwise
    # Choose by abs(MZ1-MZ), or sum(PT) if same Z1
    def bestCandByZ1Z2(self,a,b): 
        if abs(a.Z1.M-b.Z1.M) < 1e-4 : # same Z1: choose the candidate with highest-pT Z2 leptons
            print ("cmp by sumpt")
            if a.Z2.sumpt() > b.Z2.sumpt() :
                return -1
            else :
                return 1
        else : # choose based on Z1 masses
            if abs(a.Z1.M-self.ZmassValue) < abs(b.Z1.M-self.ZmassValue) :
                return -1 
            else :
                return 1

    # Choose by DbkgKin
    def bestCandByDbkgKin(self,a,b): 
        if abs((a.p4).M() - (b.p4).M())<1e-4 and a.finalState()==b.finalState() and (a.finalState() == 28561 or a.finalState()==14641) : # Equivalent: same masss (tolerance 100 keV) and same FS -> same leptons and FSR
            return self.bestCandByZ1Z2(a,b)
        if a.KD > b.KD : return -1 # choose by best dbkgkin
        else: return 1
        

