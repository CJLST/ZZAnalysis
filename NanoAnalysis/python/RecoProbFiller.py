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


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.MELASettings != None: 
            # print("***RecoProbFiller: Booking Probs")
            self.ProbHelper.bookProbs(wrappedOutputTree)

                        

        
    def analyze(self, event):
        if event.nZZCand==0: return True
        cands = Collection(event, 'ZZCand')
        leps = Collection(event, 'Lepton')
        fsrPhotons = Collection(event, "FsrPhoton")
        # Cache the daughters for each candidate
        candsDaughters = [Mela.SimpleParticleCollection_t()]*len(cands)
        candsAssociated = [None]*len(cands) # FIXME to be added
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
        
            self.ProbHelper.fillProbs(candsDaughters, candsAssociated, None)

        return True
    
