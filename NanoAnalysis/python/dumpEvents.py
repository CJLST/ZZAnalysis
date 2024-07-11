###
# Print variables, for debug purposes.
# level
# -1: lepton information (must be run after LepFiller)
# >=1 final candidates (must be run after ZZFiller)
###
from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
#from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR

### Dump reconstructed quantities, for debug purposes

charge={1:"+", -1:"-"}

class dumpEvents(Module):
    def __init__(self, level=1, nanoVersion=11):
        self.writeHistFile = False
        self.level=level
        self.nanoVersion=nanoVersion
        self.printTrigger = True
        
    def analyze(self, event):

        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")

        eventId='{}:{}:{}'.format(event.run,event.luminosityBlock,event.event)
        print(eventId, "-----------")

        if self.level >= 3 or self.level == -1: #full detail
            for i,lep in enumerate(muons):
               print(eventId+" ", end="")
               print('#{} mu{} pt={:.3g} eta={:.3g} phi={:.3g} dxy={:.3g} dz={:.3g} SIP={:.3g}'.format(i, charge[lep.charge],
                                                                                                       lep.pt,
                                                                                                       lep.eta,
                                                                                                       lep.phi,
                                                                                                       abs(lep.dxy),
                                                                                                       abs(lep.dz),
                                                                                                       lep.sip3d,
                                                                                                       ),
                           end="")
               print(' GLB={} TK={} PF={}'.format(int(lep.isGlobal),
                                                  int(lep.isTracker),
                                                  int(lep.isPFcand),
                                                  #int(lep.highPtId>0),# 1=trk,2=global ;  isHighPtId={}
                                                  ),
                     end="")
               print(' combRelIsoPF={:.3g} combRelIsoPFFSRCorr={:.3g}'.format(lep.pfRelIso03_all, lep.pfRelIso03FsrCorr),
                     end="")
               print(' relaxedId={} FullId={} FullSel={}'.format(lep.ZZRelaxedId,lep.ZZFullId,lep.ZZFullSel),
                     end="")
               if (lep.fsrPhotonIdx>=0):
                   fsr=fsrPhotons[lep.fsrPhotonIdx]
                   
                   print(' FSR: pt={:.3g} eta={:.3g}, phi={:.3g}, dREt2={:.3g}, Iso={:.3g}, isTrue={}'.format(fsr.pt, fsr.eta, fsr.phi, fsr.dROverEt2, fsr.relIso03, (None if event.run>1 else (fsr.genFsrIdx>=0))))
               else:
                   print()
   
            for i,lep in enumerate(electrons):
#  #            if not eleLooseId(lep) and not DEBUG : continue # Print only leptons passing loose ID
               print(eventId+" ", end="")
               print('#{} e{}  pt={:.3g} eta={:.3g} phi={:.3g} dxy={:.3g} dz={:.3g} SIP={:.3g}'.format(i, charge[lep.charge],
                                                                                                       lep.pt,
                                                                                                       lep.eta,
                                                                                                       lep.phi,
                                                                                                       abs(lep.dxy),
                                                                                                       abs(lep.dz),
                                                                                                       lep.sip3d),
                     end="")

               passBDT=lep.passBDT
               if self.nanoVersion < 10:
                   BDT = lep.mvaFall17V2Iso
               else: 
                   BDT = lep.mvaHZZIso
               print(' scEta={:.3g} BDT={:.3g} passBDT={} relaxedId={} FullId={} FullSel={}'.format(lep.eta+lep.deltaEtaSC, BDT, passBDT, lep.ZZRelaxedId, lep.ZZFullId, lep.ZZFullSel), end="")
               if (lep.fsrPhotonIdx>=0):
                   fsr=fsrPhotons[lep.fsrPhotonIdx]
                   print(' FSR: pt={:.3g} eta={:.3g}, phi={:.3g}, dREt2={:.3g}, Iso={:.3g}, isTrue={}'.format(fsr.pt, fsr.eta, fsr.phi, fsr.dROverEt2, fsr.relIso03, (None if event.run>1 else (fsr.genFsrIdx>=0))))
               else:
                   print()
   

        if self.level >=1 : # final candidates
            # Sync-style printout
            leps = list(muons)+ list(electrons)
#            FIXME: obsolete data format, must be updated
            ZZFlav = (leps[event.Z1_l1Idx].pdgId*leps[event.Z2_l1Idx].pdgId)**2
            print ('{}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{}:{:.6g}:{:.3f}:{:.3f}:{:.3f}:{:.3f}'.format(
                eventId,
                event.ZZ_mass,
                event.Z1_mass,
                event.Z2_mass,
                event.ZZ_massPreFSR,
                ZZFlav,
                event.ZZ_Dbkgkin,
                event.Generator_weight,
                event.puWeight,
                event.ZZ_dataMCWeight,
                event.totalWeight,
                ))

 
        if self.printTrigger :
            print("Passing triggers:") 
            branchList = event._tree.GetListOfBranches()
            for b in branchList :
                if "HLT" in b.GetName() :
                    if eval("event."+b.GetName()) :
                        print(b.GetName())

            print("TrigObjs:")
            tos =  Collection(event, "TrigObj")
            for to in tos:
                # filterBits:
                # for Electron: 0 => CaloIdL_TrackIdL_IsoVL, 1 => 1e (WPTight), 2 => 1e (WPLoose), 3 => OverlapFilter PFTau, 4 => 2e, 5 => 1e-1mu, 6 => 1e-1tau, 7 => 3e, 8 => 2e-1mu, 9 => 1e-2mu, 10 => 1e (32_L1DoubleEG_AND_L1SingleEGOr), 11 => 1e (CaloIdVT_GsfTrkIdT), 12 => 1e (PFJet), 13 => 1e (Photon175_OR_Photon200)
                # for Muon: 0 => TrkIsoVVL, 1 => Iso, 2 => OverlapFilter PFTau, 3 => 1mu, 4 => 2mu, 5 => 1mu-1e, 6 => 1mu-1tau, 7 => 3mu, 8 => 2mu-1e, 9 => 1mu-2e, 10 => 1mu (Mu50), 11 => 1mu (Mu100), 12 => 1mu-1photon
                print(to.id, to.l2pt, to.pt, to.eta, bin(to.filterBits))
        return True

