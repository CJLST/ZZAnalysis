from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.HeppyCore.utils.deltar import deltaR
import Mela


class genAngProbFiller(Module):
    """Calculates angles and proabilities with LHE-level information. 
    MELA = Pointer to MELA passed from nanoZZ4lAnalysis.py 
    """
    
    def __init__(self, MELA, settingsDict = None):
        print("***genAngProbFiller", flush=True)
        self.MELA = MELA
        self.MELAsettings = settingsDict
            
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("LHEMela_qH", "F")
        self.out.branch("LHEMela_mZ1", "F")
        self.out.branch("LHEMela_mZ2", "F")
        self.out.branch("LHEMela_costheta1", "F")
        self.out.branch("LHEMela_costheta2", "F")
        self.out.branch("LHEMela_Phi", "F")
        self.out.branch("LHEMela_costhetastar", "F")
        self.out.branch("LHEMela_Phi1", "F")
        self.out.branch("LHEMela_nativeProb", "F")
        self.out.branch("LHEMela_nativeProbprop", "F")
        
    def analyze(self, event):
        LHEPart = Collection(event, 'LHEPart')
        LHEMothers = filter(lambda p: p.MELAStatus==1, LHEPart)
        LHEDaughters = filter(lambda p: p.MELAStatus==2, LHEPart)
        LHEAssociated = filter(lambda p: p.MELAStatus==3, LHEPart)
        mothers = Mela.SimpleParticleCollection_t()
        daughters = Mela.SimpleParticleCollection_t()
        associated = Mela.SimpleParticleCollection_t()
        

        for i, mp in enumerate(LHEMothers): 
            temp_particle = Mela.SimpleParticle_t(mp.pdgId, mp.pt, mp.eta, mp.phi, mp.mass, True)
            mothers.add_particle(temp_particle)
        
        for i, dp in enumerate(LHEDaughters): 
            temp_particle = Mela.SimpleParticle_t(dp.pdgId, dp.pt, dp.eta, dp.phi, dp.mass, True)
            daughters.add_particle(temp_particle)
        
        for i, ap in enumerate(LHEAssociated): 
            temp_particle = Mela.SimpleParticle_t(ap.pdgId, ap.pt, ap.eta, ap.phi, ap.mass, True)
            associated.add_particle(temp_particle)
        
        
        self.MELA.setInputEvent(daughters, associated, mothers, 1)
        qH, mZ1, mZ2, costheta1, costheta2, Phi, costhetastar, Phi1 = self.MELA.computeDecayAngles()
        self.out.fillBranch("LHEMela_costheta1", costheta1)
        self.out.fillBranch("LHEMela_costheta2", costheta2)
        self.out.fillBranch("LHEMela_Phi", Phi)
        self.out.fillBranch("LHEMela_Phi1", Phi1)
        self.out.fillBranch("LHEMela_costhetastar", costhetastar)
        
        if self.MELAsettings != None: 
            self.MELA.differentiate_HWW_HZZ = self.MELAsettings["separatewwzz"]
            self.MELA.setProcess(self.MELAsettings["process"], self.MELAsettings["matrixelement"], self.MELAsettings["production"])
            for coupling, vals, in self.MELAsettings["couplings"].items():
                setattr(self.MELA, coupling, vals)

            # ME = getattr(Mela, 'MatrixElement')
            # PROC = getattr(Mela, 'Process')
            # ghz1 = getattr(self.MELA, 'ghz1')
            # ghg2 = getattr(self.MELA, 'ghg2')

            # print("Matrix Element: ", ME, "process: ", PROC, "ghz1: ", ghz1, "ghg2: ", ghg2)

            if self.MELAsettings["prod"] and self.MELAsettings["dec"]: 
                nativeprob = self.MELA.computeProdDecP(False)
            elif self.MELAsettings["prod"]: 
                nativeprob = self.MELA.computeProdP(False)
            elif self.MELAsettings["dec"]:
                nativeprob = self.MELA.computeP(False)
                # print(nativeprob)
            else:
                raise KeyError("Need to specify either production, decay, or computeprop!")
            
            if self.MELAsettings["computeprop"] and (self.MELAsettings["prod"] or self.MELAsettings["dec"]):
                    nativeprob_prop = self.MELA.getXPropagator(self.MELAsettings["propscheme"])
                    self.out.fillBranch("LHEMela_nativeProbprop", nativeprob_prop)
            elif self.MELAsettings["computeprop"]:
                nativeprob = self.MELA.getXPropagator(self.MELAsettings["propscheme"])
            
            self.out.fillBranch("LHEMela_nativeProb", nativeprob)
            
        

        self.MELA.resetInputEvent()
        
       
        
        
        return True
    




    # def analyze(self, event):
    #     LHEPart = Collection(event, 'LHEPart')
    #     for i, lp in enumerate(LHEPart):
    #         MELA_Stati = []
    #         if lp.status == -1: 
    #             MELA_Stati.append(1)
    #         elif lp.status == 1: 
    #             parentIdx = lp.firstMotherIdx
    #             parentPdg = LHEPart[parentIdx].pdgId

    #             gParentIdx = LHEPart[parentIdx].firstMotherIdx 
    #             if gParentIdx == -1: 
    #                 gParentPdg = 0
    #             else: 
    #                 gParentPdg = LHEPart[gParentIdx].pdgId 
    #             if parentPdg == 23 and gParentPdg == 25:
    #                 MELA_Stati.append(2)
    #             else: 
    #                 MELA_Stati.append(3)
    #         else: 
    #             MELA_Stati.append(-1)
    #         self.out.fillBranch("LHEPart_MELAStatus", MELA_Stati)
    #     return True