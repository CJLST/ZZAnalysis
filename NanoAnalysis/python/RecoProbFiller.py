from __future__ import print_function
import copy
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from  ZZAnalysis.NanoAnalysis.initializeMELA import check_enum
from ZZAnalysis.NanoAnalysis.ZZExtraFiller import *
from ZZAnalysis.NanoAnalysis.MELAProbHelper import MELAProbHelper
import Mela


class RecoProbFiller(Module):
    """Calculates proabilities with Reco-level information. 
    MELA = Pointer to MELA passed from nanoZZ4lAnalysis.py 
    MELASettings = dictionary with settings for the probabilities to be computed
    """
    
    def __init__(self, MELA, NANOVERSION, MELASettings = None):
        self.MELA = MELA
        self.NANOVERSION = NANOVERSION
        self.sortedSettings = []
        self.MELASettings = MELASettings
        self.ProbHelper = MELAProbHelper(self.MELA, self.MELASettings, "Reco")
        print("***RecoProbFiller: set for: ", self.ProbHelper.names, flush=True)


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.MELASettings != None: 
            self.ProbHelper.bookProbs(wrappedOutputTree)
        

    def analyze(self, event):
        if event.nZZCand==0: return True
        cands = Collection(event, 'ZZCand')
        leps = Collection(event, 'Lepton')
        fsrPhotons = Collection(event, "FsrPhoton")
        jets = Collection(event, 'Jet')

        # Add leading and sub-leadig jets, if present
        jets_idx = [i for i in (event.JetLeadingIdx, event.JetSubleadingIdx) if i >= 0]
        jets_MELA = Mela.SimpleParticleCollection_t()
        for idx in jets_idx:
            jet = jets[idx]
            p4 = jet.p4()
            jets_MELA.add_particle(Mela.SimpleParticle_t(0, p4.Px(), p4.Py(), p4.Pz(), p4.E()))

        # Cache the daughters for each candidate
        candsDaughters = [Mela.SimpleParticleCollection_t()]*len(cands)
        # The jets are the same for every candidate
        candsAssociated = [copy.deepcopy(jets_MELA)]*len(cands)

        for iCand, aCand in enumerate(cands):
            theCandLepIdxs = [aCand.Z1l1Idx, aCand.Z1l2Idx, aCand.Z2l1Idx, aCand.Z2l2Idx]  
            theCandLeps = [leps[i] for i in theCandLepIdxs]
            dressedLepsp4 = [ZZExtraFiller.getDressedP4(self = None, lep = l, fsrPhotons=fsrPhotons) for l in theCandLeps]
            # dressedLepsp4 = theCandLeps
            daughters = candsDaughters[iCand]
            associated = candsAssociated[iCand]
            daughters.add_particle(Mela.SimpleParticle_t(theCandLeps[0].pdgId, dressedLepsp4[0].Px(), dressedLepsp4[0].Py(), dressedLepsp4[0].Pz(), dressedLepsp4[0].E()))
            daughters.add_particle(Mela.SimpleParticle_t(theCandLeps[1].pdgId, dressedLepsp4[1].Px(), dressedLepsp4[1].Py(), dressedLepsp4[1].Pz(), dressedLepsp4[1].E()))
            daughters.add_particle(Mela.SimpleParticle_t(theCandLeps[2].pdgId, dressedLepsp4[2].Px(), dressedLepsp4[2].Py(), dressedLepsp4[2].Pz(), dressedLepsp4[2].E()))
            daughters.add_particle(Mela.SimpleParticle_t(theCandLeps[3].pdgId, dressedLepsp4[3].Px(), dressedLepsp4[3].Py(), dressedLepsp4[3].Pz(), dressedLepsp4[3].E()))

            extralep_idx = [i for i in (event.ZZCand_extraLep1Idx[iCand], event.ZZCand_extraLep2Idx[iCand]) if i >= 0]
            for idx in extralep_idx:
                lep = leps[idx]
                p4 = lep.p4()
                associated.add_particle(Mela.SimpleParticle_t(lep.pdgId, p4.Px(), p4.Py(), p4.Pz(), p4.E()))

        self.ProbHelper.fillProbs(candsDaughters, candsAssociated, None)

        return True
    
