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
        # Book MELA angle branches
        self.out.branch("GenZZ_qH", "F", limitedPrecision=12, title="The mass of the Higgs candidate as reconstructed by the 4-leptons at Gen-level.")
        self.out.branch("GenZZ_mZ1", "F",  limitedPrecision=12, title="The mass of the first decay particle as reconstructed by 2 of the Gen-level leptons.")
        self.out.branch("GenZZ_mZ2", "F",  limitedPrecision=12, title="The mass of the second decay particle as reconstructed by 2 of the Gen-level leptons.")
        self.out.branch("GenZZ_costheta1", "F",  limitedPrecision=12, title="In the Higgs' rest frame, theta_1 is the angle between the momentum of Z1 and the momentum of one of its decay products.")
        self.out.branch("GenZZ_costheta2", "F",  limitedPrecision=12, title="In the Higgs' rest frame, theta_2 is the angle between the momentum of Z2 and the momentum of one of its decay products.")
        self.out.branch("GenZZ_Phi", "F", limitedPrecision=12, title="In the Higgs' rest frame, phi is the angle between the planes formed by the decay products of the two Z bosons.")
        self.out.branch("GenZZ_costhetastar", "F", limitedPrecision=12, title="In the Higgs' rest frame, theta_star is the angle between the beamline and the momentum of one of the Higgs' decay products.")
        self.out.branch("GenZZ_Phi1", "F",  limitedPrecision=12, title="In the Higgs' rest frame, phi_1 is the angle between the decay plane of Z1 and the beamline.")

    def analyze(self, event) :
        self.addExtra('GenPart', event)
        

        return True

    def addExtra(self, collName, event) :
        genpart = Collection(event, collName)
        print("Z1l1idx: ", event.GenZZ_Z1l1Idx)
       
    
        genLeps = [genpart[event.GenZZ_Z1l1Idx], genpart[event.GenZZ_Z1l2Idx], genpart[event.GenZZ_Z2l1Idx], genpart[event.GenZZ_Z2l2Idx]]

 
        # MELA angles
        qH = -999. 
        helcosthetaZ1 = -999.
        helcosthetaZ2 = -999.
        helPhi = -999.
        costhetastar = -999.
        phistarZ1 =-999.
        mZ1 = -999.
        mZ2 = -999.
            

        if self.MELA != None: 

            daughters = Mela.SimpleParticleCollection_t()
            daughters.add_particle(Mela.SimpleParticle_t(genLeps[0].pdgId, genLeps[0].pt, genLeps[0].eta, genLeps[0].phi, genLeps[0].mass))
            daughters.add_particle(Mela.SimpleParticle_t(genLeps[1].pdgId, genLeps[1].pt, genLeps[1].eta, genLeps[1].phi, genLeps[1].mass))
            daughters.add_particle(Mela.SimpleParticle_t(genLeps[2].pdgId, genLeps[2].pt, genLeps[2].eta, genLeps[2].phi, genLeps[2].mass))
            daughters.add_particle(Mela.SimpleParticle_t(genLeps[3].pdgId, genLeps[3].pt, genLeps[3].eta, genLeps[3].phi, genLeps[3].mass))
            

            self.MELA.setInputEvent(daughters, None, None, 0)

            qH, mZ1, mZ2, helcosthetaZ1, helcosthetaZ2, helPhi, costhetastar, phistarZ1 = self.MELA.computeDecayAngles() 

            self.MELA.resetInputEvent()

            

        # Fill MELA angle branches
        self.out.fillBranch("GenZZ_qH", qH)
        self.out.fillBranch("GenZZ_mZ1", mZ1)
        self.out.fillBranch("GenZZ_mZ2", mZ2)
        self.out.fillBranch("GenZZ_costheta1", helcosthetaZ1)
        self.out.fillBranch("GenZZ_costheta2", helcosthetaZ2)
        self.out.fillBranch("GenZZ_Phi", helPhi)
        self.out.fillBranch("GenZZ_costhetastar", costhetastar)
        self.out.fillBranch("GenZZ_Phi1", phistarZ1)



    


        
