###
# Build Z and ZZ candidates.
# ZZ candidates are stored as a nanoAOD collection, so that several candidates per event can be saved (for eg. ID optimization studies).
# The index of the best ZZ candidate in each event is stored as bestCandIdx
# Filter: only events with at least one ZZ candidate are kept
###
from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.HeppyCore.utils.deltar import deltaR
from ROOT import LeptonSFHelper

from functools import cmp_to_key
from ROOT import Mela, SimpleParticle_t, SimpleParticleCollection_t, TVar, TLorentzVector
from ctypes import c_float

class StoreOption:
    # Define which candidates should be stored:
    # BestCandSRCR = only the best SR candidate and the best candidate for each of the ZLL CRs are stored. This is normally the default for production jobs
    # BestCandOnly = only the best SR candidate in the event is retained (useful for development/debugging)
    # AllSRCands   = keep all SR candidates passing the full selection cuts (including permutation of leptons); no CR
    # AllWithRelaxedMuId = keep any SR candidate that can be made, even if leptons don't pass ID cuts (for ID cut optimization).
    BestCandSRCR, BestCandOnly, AllCands, AllWithRelaxedMuId = range(0,4)


class ZZFiller(Module):

    def __init__(self, runMELA, bestCandByMELA, isMC, year, processCR, debug=False):
        print("***ZZFiller: isMC:", isMC, "year:", year, flush=True)
        self.writeHistFile = False
        self.isMC = isMC
        self.year = year
        self.runMELA = runMELA or bestCandByMELA
        if bestCandByMELA :
            self.bestCandCmp = self.bestCandByDbkgKin
        else:
            self.bestCandCmp = self.bestCandByZ1Z2

        self.addSSCR = processCR
        self.addOSCR = processCR    
        self.addSIPCR = processCR

        self.DEBUG = debug
        self.ZmassValue = 91.1876;

        self.candsToStore = StoreOption.BestCandOnly # Only store the best candidate for the SR


        # Pre-selection of leptons to build Z and LL candidates.
        # Normally it is the full ID + iso if only the SR is considered, or the relaxed ID if CRs are also filled,
        # but can be relaxed further to store loose candidates for ID studies (see below)
        if self.addSIPCR or self.addOSCR or self.addSSCR :
            self.leptonPresel = (lambda l : l.ZZRelaxedIdNoSIP) # minimal selection good for all CRs: no SIP, no ID, no iso
        else : # SR only
            self.leptonPresel = (lambda l : l.ZZFullSel)

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
        self.lepSFHelper = LeptonSFHelper(False) # FIXME for 2016 UL samples: requires passing bool preVFP

        if self.runMELA :
            sqrts=13.;
            if year>=2022 :
                sqrts=13.6
            self.mela = Mela(sqrts, 125, TVar.ERROR)
            self.mela.setCandidateDecayMode(TVar.CandidateDecay_ZZ)
            print("", flush=True) # avoids MELA init messages to mix with job output 

        # Example of adding control histograms (requires self.writeHistFile = True)
        # def beginJob(self,histFile=None, histDirName=None):
        #    Module.beginJob(self, histFile, histDirName+"_ZZFiller")
        #    self.histFile=None # Hack to prevent histFile to be closed before other modules write their histograms
        #    self.h_ZZMass = ROOT.TH1F('ZZMass','ZZMass',130,70,200)
        #    self.addObject(self.h_ZZMass)


    def endJob(self):
         print("", flush=True)


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch("nZCand", "I", title="Z candidates passing the full H4l selection")
        self.out.branch("ZCand_mass", "F", lenVar="nZCand", title="mass")
#         self.out.branch("ZCand_pt", "F", lenVar="nZCand")
#         self.out.branch("ZCand_eta", "F", lenVar="nZCand")
#         self.out.branch("ZCand_phi", "F", lenVar="nZCand")
        self.out.branch("ZCand_flav", "F", lenVar="nZCand", title="Product of the pdgIds of the 2 daughters")
        self.out.branch("ZCand_l1Idx", "S", lenVar="nZCand", title="index of 1st daughter in Electron+Muon merged collection")
        self.out.branch("ZCand_l2Idx", "S", lenVar="nZCand", title="index of 2nd daughter in Electron+Muon merged collection")
        self.out.branch("ZCand_fsr1Idx", "S", lenVar="nZCand", title="index of FSR associated to l1 (-1 if none)")
        self.out.branch("ZCand_fsr2Idx", "S", lenVar="nZCand", title="index of FSR associated to l2 (-1 if none)") 
        self.out.branch("bestZIdx", "S", title="Best Z in the event (mass closest to mZ)")

        self.out.branch("nZZCand", "I", title="ZZ candidates passing the full H4l selection")
        self.out.branch("ZZCand_mass", "F", lenVar="nZZCand", title="mass")
#         self.out.branch("ZZCand_pt", "F", lenVar="nZZCand")
#         self.out.branch("ZZCand_eta", "F", lenVar="nZZCand")
#         self.out.branch("ZZCand_phi", "F", lenVar="nZZCand")
        self.out.branch("ZZCand_massPreFSR", "F", lenVar="nZZCand", title="mass without FSR photons")
        self.out.branch("ZZCand_Z1mass", "F", lenVar="nZZCand", title="Z1 mass")
        self.out.branch("ZZCand_Z1flav", "I", lenVar="nZZCand", title="Product of the pdgIds of the 2 Z1 daughters")
        self.out.branch("ZZCand_Z2mass", "F", lenVar="nZZCand", title="Z2 mass")
        self.out.branch("ZZCand_Z2flav", "I", lenVar="nZZCand", title="Product of the pdgIds of the 2 Z2 daughters")
        self.out.branch("ZZCand_KD", "F", lenVar="nZZCand", title="Kinematic discriminant for the choice of best candidate")
        self.out.branch("ZZCand_Z2sumpt", "F", lenVar="nZZCand", title="sum of Z2 daughter pts (used in the choice of best candidate)")
        # Note: lepton indices are numbered for leps=list(electrons)+list(muons) and run up to nlep=len(leps);
        # no special ordering of l1, l2 is applied
        self.out.branch("ZZCand_Z1l1Idx", "S", lenVar="nZZCand", title="Index of 1st Z1 daughter in the Electron+Muon merged collection")
        self.out.branch("ZZCand_Z1l2Idx", "S", lenVar="nZZCand", title="Index of 2nd Z1 daughter in the Electron+Muon merged collection")
        self.out.branch("ZZCand_Z2l1Idx", "S", lenVar="nZZCand", title="Index of 1st Z2 daughter in the Electron+Muon merged collection")
        self.out.branch("ZZCand_Z2l2Idx", "S", lenVar="nZZCand", title="Index of 2nd Z2 daughter in the Electron+Muon merged collection")
        for ID in self.muonIDs :
            self.out.branch("ZZCand_mu"+ID["name"], "O", lenVar="nZZCand")
        for var in self.muonIDVars :
            self.out.branch("ZZCand_mu"+var["name"], "F", lenVar="nZZCand")
        if self.isMC :
            self.out.branch("ZZCand_dataMCWeight", "F", lenVar="nZZCand", title="data/MC efficiency correction weight")
        self.out.branch("bestCandIdx", "S", title="Seleced ZZ candidate in the event")


        if self.addSSCR or self. addOSCR or self.addSIPCR :
            self.out.branch("nZLLCand", "I", title="Z+LL control region candidates")
            self.out.branch("ZLLCand_mass", "F", lenVar="nZLLCand")
            self.out.branch("ZLLCand_massPreFSR", "F", lenVar="nZLLCand")
#             self.out.branch("ZLLCand_pt", "F", lenVar="nZLLCand")
#             self.out.branch("ZLLCand_eta", "F", lenVar="nZLLCand")
#             self.out.branch("ZLLCand_phi", "F", lenVar="nZLLCand")
            self.out.branch("ZLLCand_Z1mass", "F", lenVar="nZLLCand")
            self.out.branch("ZLLCand_Z1flav", "I", lenVar="nZLLCand")
            self.out.branch("ZLLCand_Z2mass", "F", lenVar="nZLLCand")
            self.out.branch("ZLLCand_Z2flav", "S", lenVar="nZLLCand")
            self.out.branch("ZLLCand_Z1l1Idx", "S", lenVar="nZLLCand") 
            self.out.branch("ZLLCand_Z1l2Idx", "S", lenVar="nZLLCand")
            self.out.branch("ZLLCand_Z2l1Idx", "S", lenVar="nZLLCand")
            self.out.branch("ZLLCand_Z2l2Idx", "S", lenVar="nZLLCand")
            self.out.branch("ZLLCand_KD", "F", lenVar="nZLLCand")
            self.out.branch("ZLLbestSSIdx", "S", title="best candidate for the SS CR")
            self.out.branch("ZLLbest2P2FIdx", "S", title="best candidate for the 2P2F CR")
            self.out.branch("ZLLbest3P1FIdx", "S", title="best candidate for the 3P1F CR")
            self.out.branch("ZLLbestSIPCRIdx", "S", title="best candidate for the SIP CR")



    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        eventId='{}:{}:{}'.format(event.run,event.luminosityBlock,event.event)
        if self.DEBUG : print ('Event '+eventId)
 
        # Collections
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")
        leps = list(electrons)+list(muons)
        nlep=len(leps)

        ### Z combinatorial over selected leps (after FSR-corrected ISO cut for muons)
        Zs = []
        for i,l1 in enumerate(leps):
            if self.leptonPresel(l1) :
                for j in range(i+1,nlep):
                    l2 = leps[j]
                    nPassLep = 0
                    isOSSF = False
                    isSR   = False
                    is1FCR = False
                    is2FCR = False
                    isSSCR = False
                    isSIPCR = False
                    if self.leptonPresel(l2):
                        if l1.pdgId == -l2.pdgId : #OS,SF: candidate for signal region 
                            isOSSF = True
                            if l1.ZZRelaxedId and l2.ZZRelaxedId:
                                nPassLep = int(l1.ZZFullSel) + int(l2.ZZFullSel) # 0,1,2 leptons passing full sel
                                if nPassLep == 2 : isSR = True # SR (to be used to set the best candidate)
                                if self.addOSCR :
                                    if nPassLep == 0 : is2FCR = True
                                    elif nPassLep == 1 : is1FCR = True
                        elif l1.pdgId == l2.pdgId : # SS control regions
                            if self.addSSCR and l1.ZZRelaxedId and l2.ZZRelaxedId : isSSCR = True
                            if self.addSIPCR and l1.ZZFullSelNoSIP and l2.ZZFullSelNoSIP : isSIPCR = True
                            if not (isSSCR or isSIPCR) : continue
                        else:
                            continue
                        
                        # Set a default order for OS leptons in the Z candidate: 1=l+, 2=l-
                        idx1, idx2 = i, j
                        if (l1.pdgId*l2.pdgId<0 and l1.pdgId>0) :
                            idx1, idx2 = idx2, idx1
                        
                        aZ = self.ZCand(idx1, idx2, leps, fsrPhotons)
                        aZ.isOSSF = isOSSF
                        aZ.isSR   = isSR
                        aZ.is1FCR = is1FCR
                        aZ.is2FCR = is2FCR
                        aZ.isSSCR = isSSCR
                        aZ.isSIPCR = isSIPCR
                        
                        zmass = aZ.M
                        if self.DEBUG: print('Z={:.4g} pt1={:.3g} pt2={:.3g} fsr1={} fsr2={} SR={} 1F={} 2F={} SS={} SSSIP={}'.format(zmass, l1.pt, l2.pt, l1.fsrPhotonIdx,  l2.fsrPhotonIdx, isSR, is1FCR, is2FCR, isSSCR, isSIPCR))
                        if (zmass>12. and zmass<120.):
                            Zs.append(aZ)


        ### Build ZZ and ZLL combinations passing the ZZ selection
        if len(Zs) >= 2:
            ZZs = []
            bestCandIdx = -1

            ZLLs = []
            ZLLsTemp = []
            best3P1FCRIdx = -1
            best2P2FCRIdx = -1
            bestSSCRIdx = -1
            bestSIPCRIdx = -1

            ### Signal region
            for iZ,aZ in enumerate(Zs):
                for jZ in range(iZ+1, len(Zs)):                    
                    if Zs[iZ].isOSSF and Zs[jZ].isOSSF : # 2 OSSSF Zs
                        ZZ = self.makeCand(Zs[iZ],Zs[jZ])
                        if ZZ == None: continue
                        ZZs.append(ZZ)

                        if self.DEBUG : print("ZZ:", len(ZZs)-1, ZZ.p4.M(), ZZ.Z1.M, ZZ.Z2.M, ZZ.Z2.sumpt(), ZZ.finalState(), ZZ.p_GG_SIG_ghg2_1_ghz1_1_JHUGen, ZZ.p_QQB_BKG_MCFM, ZZ.KD, (ZZ.Z1.isSR and ZZ.Z2.isSR))

                        #Search for the the best cand in the SR (ie among those passing full ID cuts)
                        if ZZ.Z1.isSR and ZZ.Z2.isSR :
                            if bestCandIdx<0 or self.bestCandCmp(ZZ,ZZs[bestCandIdx]) < 0: bestCandIdx = len(ZZs)-1

            if self.DEBUG : print("bestCand:", bestCandIdx)


            ### ZLL combinations for control regions, made of 1 good Z + 1 ll pair;
            ### these are considered only if no SR candidate is present in the event
            if bestCandIdx < 0 and (self.addSSCR or self. addOSCR or self.addSIPCR) :
                for iZ1,Z1 in enumerate(Zs):
                    if Z1.isSR : 
                        for iZ2,Z2 in enumerate(Zs):
                            if Z2.is1FCR or Z2.is2FCR or Z2.isSSCR or Z2.isSIPCR :
                                ZLL = self.makeCand(Z1, Z2, sortZsByMass=False, fillIDVars=False)
                                if ZLL == None: continue
                                if ZLL.Z2.is2FCR  and (best2P2FCRIdx<0 or self.bestCandCmp(ZLL,ZLLsTemp[best2P2FCRIdx]) < 0) : best2P2FCRIdx = len(ZLLsTemp)
                                if ZLL.Z2.is1FCR  and (best3P1FCRIdx<0 or self.bestCandCmp(ZLL,ZLLsTemp[best3P1FCRIdx]) < 0) : best3P1FCRIdx = len(ZLLsTemp)
                                if ZLL.Z2.isSSCR  and (bestSSCRIdx<0 or self.bestCandCmp(ZLL,ZLLsTemp[bestSSCRIdx]) < 0) : bestSSCRIdx = len(ZLLsTemp)
                                if ZLL.Z2.isSIPCR and (bestSIPCRIdx<0 or self.bestCandCmp(ZLL,ZLLsTemp[bestSIPCRIdx]) < 0) : bestSIPCRIdx = len(ZLLsTemp)
                                ZLLsTemp.append(ZLL)
                                
                # Only one candidate per method should be selected per event: check if
                if best2P2FCRIdx >= 0 and best3P1FCRIdx >= 0 :
                    print ("WARNING: event", eventId, "has CR candidates in both 2P2F amd 3P1F regions")   #FIXME choose the best among the two

                # Store only ZLL candidates that belong to at least 1 CR
                for iZLL, ZLL in enumerate(ZLLsTemp) :
                    select = False
                    if iZLL == best2P2FCRIdx :
                        best2P2FCRIdx = len(ZLLs)
                        select = True
                    if iZLL == best3P1FCRIdx :
                        best3P1FCRIdx = len(ZLLs)
                        select = True
                    if iZLL == bestSSCRIdx :
                        bestSSCRIdx = len(ZLLs)
                        select = True
                    if iZLL == bestSIPCRIdx :
                        bestSIPCRIdx = len(ZLLs)
                        select = True
                    if select : ZLLs.append(ZLL)
                    if self.DEBUG: print("ZLL:", iZLL, ZLL.p4.M(), ZLL.Z1.M, ZLL.Z2.M, ZLL.Z2.sumpt(), ZLL.finalState(), ZLL.p_GG_SIG_ghg2_1_ghz1_1_JHUGen, ZLL.p_QQB_BKG_MCFM, ZLL.KD,
                                         "2P2F:", int(iZLL==best2P2FCRIdx), "3P1F:", int(iZLL==best3P1FCRIdx), "SS:", int(iZLL == bestSSCRIdx), "SIP:", int(iZLL == bestSIPCRIdx))
                    
            if self.candsToStore == StoreOption.BestCandOnly : # keep only the best cand as single element of the ZZ collection
                if bestCandIdx >= 0 :
                    ZZs = [ZZs[bestCandIdx]]
                    bestCandIdx = 0
                else :
                    ZZs = []

            ### Skip events with no candidates
            if len(ZZs) == 0 and len(ZLLs) == 0: return False

            ### Now fill the variables to be stored as output
            # Fill Zs - we are only interested in signal ones
            bestZIdx = -1
            ZCand_mass = []
            ZCand_flav = []
            ZCand_l1Idx = []
            ZCand_l2Idx = []
            ZCand_fsr1Idx = []
            ZCand_fsr2Idx = []

            for iZ, aZ in enumerate(Zs) :
                if aZ.isSR :
                    if bestZIdx<0 or abs(aZ.M-self.ZmassValue)<abs(ZCand_mass[bestZIdx]-self.ZmassValue) :
                        bestZIdx = len(ZCand_mass)
                    ZCand_mass.append(aZ.p4.M())
#                     ZCand_pt.append(aZ.p4.pt())
#                     ZCand_eta.append(aZ.p4.eta())
#                     ZCand_phi.append(aZ.p4.phi())
                    ZCand_flav.append(aZ.finalState())
                    ZCand_l1Idx.append(aZ.l1Idx)
                    ZCand_l2Idx.append(aZ.l2Idx)
                    ZCand_fsr1Idx.append(aZ.fsr1Idx)
                    ZCand_fsr2Idx.append(aZ.fsr2Idx)

            self.out.fillBranch("nZCand", len(ZZs))
            self.out.fillBranch("ZCand_mass", ZCand_mass)
#             self.out.fillBranch("ZCand_pt", ZCand_pt)
#             self.out.fillBranch("ZCand_eta", ZCand_eta)
#             self.out.fillBranch("ZCand_phi", ZCand_phi)
            self.out.fillBranch("ZCand_flav", ZCand_flav)
            self.out.fillBranch("ZCand_l1Idx", ZCand_l1Idx)
            self.out.fillBranch("ZCand_l2Idx", ZCand_l2Idx)
            self.out.fillBranch("ZCand_fsr1Idx", ZCand_fsr1Idx)
            self.out.fillBranch("ZCand_fsr2Idx", ZCand_fsr2Idx)
            self.out.fillBranch("bestZIdx", bestZIdx)

            ### Fill ZZ candidates
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
                if self.year < 2022 : #FIXME
                    if self.isMC: ZZCand_wDataMC[iZZ] =  self.getDataMCWeight(ZZ.leps())
                else :
                    ZZCand_wDataMC[iZZ] = 1.
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

            self.out.fillBranch("bestCandIdx", bestCandIdx)


            if self.addSSCR or self. addOSCR or self.addSIPCR :
                ZLLCand_mass   = [0.]*len(ZLLs)
                ZLLCand_massPreFSR = [0.]*len(ZLLs)
#                 ZLLCand_pt     = [0.]*len(ZLLs)
#                 ZLLCand_eta    = [0.]*len(ZLLs)
#                 ZLLCand_phi    = [0.]*len(ZLLs)
                ZLLCand_Z1mass = [0.]*len(ZLLs)
                ZLLCand_Z1flav = [0.]*len(ZLLs)
                ZLLCand_Z2mass = [0.]*len(ZLLs)
                ZLLCand_Z2flav = [0.]*len(ZLLs)
                ZLLCand_Z1l1Idx = [-1]*len(ZLLs)
                ZLLCand_Z1l2Idx = [-1]*len(ZLLs)
                ZLLCand_Z2l1Idx = [-1]*len(ZLLs)
                ZLLCand_Z2l2Idx = [-1]*len(ZLLs)
                ZLLCand_KD     = [0.]*len(ZLLs)

                for iZLL, ZLL in enumerate(ZLLs) :
                    ZLLCand_mass[iZLL] = ZLL.p4.M()
                    ZLLCand_massPreFSR[iZLL] = ZLL.massPreFSR()
#                     ZLLCand_pt[iZLL] = ZLL.p4.pt()
#                     ZLLCand_eta[iZLL] = ZLL.p4.eta()
#                     ZLLCand_phi[iZLL] = ZLL.p4.phi()
                    ZLLCand_Z1mass[iZLL] = ZLL.Z1.M
                    ZLLCand_Z1flav[iZLL] = ZLL.Z1.finalState()
                    ZLLCand_Z2mass[iZLL] = ZLL.Z2.M
                    ZLLCand_Z2flav[iZLL] = ZLL.Z2.finalState()
                    ZLLCand_Z1l1Idx[iZLL] = ZLL.Z1.l1Idx
                    ZLLCand_Z1l2Idx[iZLL] = ZLL.Z1.l2Idx
                    ZLLCand_Z2l1Idx[iZLL] = ZLL.Z2.l1Idx
                    ZLLCand_Z2l2Idx[iZLL] = ZLL.Z2.l2Idx
                    ZLLCand_KD[iZLL] = ZLL.KD

                self.out.fillBranch("nZLLCand", len(ZLLs))
                self.out.fillBranch("ZLLCand_mass",   ZLLCand_mass)
                self.out.fillBranch("ZLLCand_massPreFSR",   ZLLCand_massPreFSR)
#                 self.out.fillBranch("ZLLCand_pt",     ZLLCand_pt)
#                 self.out.fillBranch("ZLLCand_eta",    ZLLCand_eta)
#                 self.out.fillBranch("ZLLCand_phi",    ZLLCand_phi)
                self.out.fillBranch("ZLLCand_Z1mass", ZLLCand_Z1mass)
                self.out.fillBranch("ZLLCand_Z1flav", ZLLCand_Z1flav)
                self.out.fillBranch("ZLLCand_Z2mass", ZLLCand_Z2mass)
                self.out.fillBranch("ZLLCand_Z2flav", ZLLCand_Z2flav)
                self.out.fillBranch("ZLLCand_Z1l1Idx", ZLLCand_Z1l1Idx)
                self.out.fillBranch("ZLLCand_Z1l2Idx", ZLLCand_Z1l2Idx)
                self.out.fillBranch("ZLLCand_Z2l1Idx", ZLLCand_Z2l1Idx)
                self.out.fillBranch("ZLLCand_Z2l2Idx", ZLLCand_Z2l2Idx)
                self.out.fillBranch("ZLLCand_KD",     ZLLCand_KD)
                if self.addSSCR :
                    self.out.fillBranch("ZLLbestSSIdx",  bestSSCRIdx)
                if self.addOSCR :
                    self.out.fillBranch("ZLLbest2P2FIdx", best2P2FCRIdx)
                    self.out.fillBranch("ZLLbest3P1FIdx", best3P1FCRIdx)
                if self.addSIPCR :
                    self.out.fillBranch("ZLLbestSIPCRIdx", bestSIPCRIdx)


            ### Fill control plot (example)
            # self.h_ZZMass.Fill(ZZCand_mass[bestCandIdx])

            return True

        return False # No candidate found, do not write the event.


    # Temporary class to store information on a ZZ candidate.
    # NOTE: We may need to move Zs into the Event as persistent objects, and to be re-used for CRs. They would be built by a separate module in that case.
    class ZCand: 
        def __init__(self, l1Idx, l2Idx, leps, fsrPhotons):
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

           SF = self.lepSFHelper.getSF(self.year, myLepID, lep.pt, lep.eta, mySCeta, isCrack)
#           SF_Unc = self.lepSFHelper.getSFError(year, myLepID, lep.pt, lep.eta, mySCeta, isCrack)
           dataMCWeight *= SF

       return dataMCWeight


    ### Comparators to select the best candidate in the event. Return -1 if a is better than b, +1 otherwise
    # Choose by abs(MZ1-MZ), or sum(PT) if same Z1
    def bestCandByZ1Z2(self,a,b): 
        if abs(a.Z1.M-b.Z1.M) < 1e-4 : # same Z1: choose the candidate with highest-pT Z2 leptons
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
        

    # Build a ZZ object from given Za, Zb pair, if it passes selection cuts; None is returned otherwise.
    # All relevant candidate variables are computed for the candidate. Options:
    #  sortZsByMass : set Z1 and Z2 according to closest-Mz criteria (for SR); otherwise, specified order is
    #                 kept (useful for CRs)
    def makeCand(self, Za, Zb, sortZsByMass=True, fillIDVars=True) :
        # check that these Zs are mutually exclusive (not sharing the same lepton) 
        if Za.l1Idx==Zb.l1Idx or Za.l2Idx==Zb.l2Idx or Za.l2Idx==Zb.l1Idx or Za.l2Idx==Zb.l2Idx: return None
        
        # set Z1 and Z2
        Z1, Z2 = Za, Zb
        if sortZsByMass and abs(Zb.M-self.ZmassValue) < abs(Za.M-self.ZmassValue):
            Z1, Z2 = Z2, Z1

        # Z1 mass cut
        if Z1.M <= 40. : return None

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

        if not (passQCD and passDeltaR) : return None

        # trigger acceptance cuts (20,10 GeV)
        lepPts.sort()
        if not (lepPts[3]>20. and lepPts[2]>10.) : return None

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
            if (abs(mZa-self.ZmassValue)<abs(Z1.M-self.ZmassValue)) and mZb < 12.: return None

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
            res = c_float(0.)
            self.mela.computeP(res,True)
            p_GG_SIG_ghg2_1_ghz1_1_JHUGen = res.value
            self.mela.setProcess(TVar.bkgZZ, TVar.MCFM, TVar.ZZQQB)
            self.mela.computeP(res,True)
            p_QQB_BKG_MCFM = res.value
            self.mela.resetInputEvent()
        if (p_GG_SIG_ghg2_1_ghz1_1_JHUGen+p_QQB_BKG_MCFM == 0.) :
            print ("ERROR", p_GG_SIG_ghg2_1_ghz1_1_JHUGen, p_QQB_BKG_MCFM)
            p_QQB_BKG_MCFM = 1. # FIXME: fix for error with message: "TUtil::CheckPartonMomFraction: At least one of the parton momentum fractions is greater than 1."
        ZZ = self.ZZCand(Z1, Z2, p_GG_SIG_ghg2_1_ghz1_1_JHUGen, p_QQB_BKG_MCFM)

        # Set flags for IDs passed by all leptons of candidate (muon only for the time being), which are 
        # used for studies on ID tuning
        if fillIDVars:
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

        return ZZ
