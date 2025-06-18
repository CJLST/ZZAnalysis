### 
# Compute the gen-level angles. 
#
#
###

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

#from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR
from ctypes import c_float
import Mela

class genExtraFiller(Module):
    def __init__(self, MELA):
        print("***genExtraFiller")
        self.MELA = MELA

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.bookExtra("GenPart")
        

    def bookExtra(self, collName) :
        theLenVar="n"+collName
        # Book MELA angle branches
        self.out.branch("GenZZ" + "_costheta1", "F",  limitedPrecision=12)
        self.out.branch("GenZZ" + "_costheta2", "F",  limitedPrecision=12)
        self.out.branch("GenZZ" + "_Phi", "F", limitedPrecision=12)
        self.out.branch("GenZZ" + "_costhetastar", "F", limitedPrecision=12)
        self.out.branch("GenZZ" + "_Phi1", "F",  limitedPrecision=12)

    def analyze(self, event) :
        self.addExtra('GenPart', event)
        

        return True

    def addExtra(self, collName, event) :
        genpart = Collection(event, collName)
        print("Z1l1idx: ", event.GenZZ_Z1l1Idx)
       
    
        genLeps = [genpart[event.GenZZ_Z1l1Idx], genpart[event.GenZZ_Z1l2Idx], genpart[event.GenZZ_Z2l1Idx], genpart[event.GenZZ_Z2l2Idx]]

 
        # MELA angle arrays
        helcosthetaZ1s = [-999.] * len(genpart)
        helcosthetaZ2s = [-999.] * len(genpart)
        helPhis = [-999.] * len(genpart)
        costhetastars = [-999.] * len(genpart)
        phistarZ1s = [-999.] * len(genpart)
        mZ1s = [-999.] * len(genpart)
        mZ2s = [-999.] * len(genpart)
            

        if self.MELA != None: 

            daughters = Mela.SimpleParticleCollection_t()
            daughters.add_particle(Mela.SimpleParticle_t(genLeps[0].pdgId, genLeps[0].pt, genLeps[0].eta, genLeps[0].phi, genLeps[0].mass))
            # print("daughter 1: ", [genLeps[0].pdgId, genLeps[0].pt, genLeps[0].eta, genLeps[0].phi, getMass(genLeps[0].pdgId)])
            
            daughters.add_particle(Mela.SimpleParticle_t(genLeps[1].pdgId, genLeps[1].pt, genLeps[1].eta, genLeps[1].phi, genLeps[1].mass))
            # print("daughter 2: ", [genLeps[1].pdgId, genLeps[1].pt, genLeps[1].eta, genLeps[1].phi, getMass(genLeps[1].pdgId)])
            daughters.add_particle(Mela.SimpleParticle_t(genLeps[2].pdgId, genLeps[2].pt, genLeps[2].eta, genLeps[2].phi, genLeps[2].mass))
            # print("daughter 3: ", [genLeps[2].pdgId, genLeps[2].pt, genLeps[2].eta, genLeps[2].phi, getMass(genLeps[2].pdgId)])
            daughters.add_particle(Mela.SimpleParticle_t(genLeps[3].pdgId, genLeps[3].pt, genLeps[3].eta, genLeps[3].phi, genLeps[3].mass))
            # print("daughter 4: ", [genLeps[3].pdgId, genLeps[3].pt, genLeps[3].eta, genLeps[3].phi, getMass(genLeps[3].pdgId)])
            

            self.MELA.setInputEvent(daughters, None, None, 0)

            qH, mZ1, mZ2, helcosthetaZ1, helcosthetaZ2, helPhi, costhetastar, phistarZ1 = self.MELA.computeDecayAngles() 

            self.MELA.resetInputEvent()

            # Store angles
            helcosthetaZ1s = helcosthetaZ1
            helcosthetaZ2s = helcosthetaZ2
            helPhis = helPhi
            costhetastars = costhetastar
            phistarZ1s = phistarZ1
            mZ1s = mZ1
            mZ2s = mZ2  

        # Fill MELA angle branches
        self.out.fillBranch("GenZZ_costheta1", helcosthetaZ1s)
        self.out.fillBranch("GenZZ_costheta2", helcosthetaZ2s)
        self.out.fillBranch("GenZZ_Phi", helPhis)
        self.out.fillBranch("GenZZ_costhetastar", costhetastars)
        self.out.fillBranch("GenZZ_Phi1", phistarZ1s)


def getMass(pdgId):
    if abs(pdgId) == 11:
        return 0.000511
    elif abs(pdgId) == 13:
        return 0.1057
    elif abs(pdgId) == 15:
        return 1.777
    


        
