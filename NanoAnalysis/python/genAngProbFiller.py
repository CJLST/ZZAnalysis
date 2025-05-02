from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.HeppyCore.utils.deltar import deltaR
import Mela


class genAngProbFiller(Module):
    """Calculates angles and proabilities with LHE-level information. 
    MELA = Pointer to MELA passed from nanoZZ4lAnalysis.py 
    """
    
    def __init__(self, MELA, NANOVERSION, settingsDict = None):
        print("***genAngProbFiller", flush=True)
        self.MELA = MELA
        self.MELAsettings = settingsDict
        self.NANOVERSION = NANOVERSION
        print("THIS IS NANOVERSION: ", self.NANOVERSION)
            
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
        
        mothers = Mela.SimpleParticleCollection_t()
        daughters = Mela.SimpleParticleCollection_t()
        associated = Mela.SimpleParticleCollection_t()


        
        ## Only run if mother-daughter associations are available, i.e. nanoAODv15 or newer. 
        if self.NANOVERSION >= 15: 
            LHEMothers = filter(lambda p: p.MELAStatus==1, LHEPart)
            LHEDaughters = filter(lambda p: p.MELAStatus==2, LHEPart)
            LHEAssociated = filter(lambda p: p.MELAStatus==3, LHEPart)
            for i, mp in enumerate(LHEMothers): 
                temp_particle = Mela.SimpleParticle_t(mp.pdgId, mp.pt, mp.eta, mp.phi, mp.mass, True)
                mothers.add_particle(temp_particle)
            
            for i, dp in enumerate(LHEDaughters): 
                temp_particle = Mela.SimpleParticle_t(dp.pdgId, dp.pt, dp.eta, dp.phi, dp.mass, True)
                daughters.add_particle(temp_particle)
            
            for i, ap in enumerate(LHEAssociated): 
                temp_particle = Mela.SimpleParticle_t(ap.pdgId, ap.pt, ap.eta, ap.phi, ap.mass, True)
                associated.add_particle(temp_particle)
        else: 
            # print("genAngProbFiller: NANOAODv14 or older, using workaround")
            for i, lp in enumerate(LHEPart): 
                temp_particle = Mela.SimpleParticle_t(lp.pdgId, lp.pt, lp.eta, lp.phi, lp.mass, True)
                if lp.status == -1: 
                    mothers.add_particle(temp_particle)
                elif lp.status == 1: 
                    if abs(lp.pdgId) in [11, 13, 15] and i >= len(LHEPart) - 4:
                        daughters.add_particle(temp_particle)
                    
                elif lp.status == 2: 
                    if lp.pdgId == 25: 
                        higgs = temp_particle
                        hMass = lp.mass
                    else:
                        continue
        
                
        #Check if selected 4-leps match the higgs: 
        if abs(hMass - daughters.MTotal()) < 0.01:
            # self.MELA.setInputEvent(daughters, associated, mothers, 1)
            self.MELA.setInputEvent(daughters, None, None, 0)
            qH, mZ1, mZ2, costheta1, costheta2, Phi, costhetastar, Phi1 = self.MELA.computeDecayAngles()
            self.out.fillBranch("LHEMela_costheta1", costheta1)
            self.out.fillBranch("LHEMela_costheta2", costheta2)
            self.out.fillBranch("LHEMela_Phi", Phi)
            self.out.fillBranch("LHEMela_Phi1", Phi1)
            self.out.fillBranch("LHEMela_costhetastar", costhetastar)
        else: 
            print(hMass - daughters.MTotal())
            
        
        # self.MELA.setInputEvent(daughters, associated, mothers, 1)
        # qH, mZ1, mZ2, costheta1, costheta2, Phi, costhetastar, Phi1 = self.MELA.computeDecayAngles()
        # self.out.fillBranch("LHEMela_costheta1", costheta1)
        # self.out.fillBranch("LHEMela_costheta2", costheta2)
        # self.out.fillBranch("LHEMela_Phi", Phi)
        # self.out.fillBranch("LHEMela_Phi1", Phi1)
        # self.out.fillBranch("LHEMela_costhetastar", costhetastar)

        
        if self.MELAsettings != None: 
            
            self.MELA.differentiate_HWW_HZZ = self.MELAsettings["separatewwzz"]
            self.MELA.setProcess(self.MELAsettings["process"], self.MELAsettings["matrixelement"], self.MELAsettings["production"])

            # Needs to be changed:
            # for coupling, vals, in self.MELAsettings["couplings"].items():
            #     setattr(self.MELA, coupling, vals)

            #Need to loop over all dictionaries in the "probabilities" list passed by the pyFragment. 
            for i, prob in enumerate(self.MELAsettings["probabilities"]):

                for coupling, vals, in prob:
                    if coupling != 'name': 
                        setattr(self.MELA, coupling, vals)


                

                eventSide = ""
                if self.MELAsettings["prod"] and self.MELAsettings["dec"]: 
                    probability = self.MELA.computeProdDecP(False)
                    eventSide = "ProdDec"
                elif self.MELAsettings["prod"]: 
                    probability = self.MELA.computeProdP(False)
                    eventSide = "Prod"
                elif self.MELAsettings["dec"]:
                    probability = self.MELA.computeP(False)
                    eventSide = "Dec"
                    # print(nativeprob)
                else:
                    raise KeyError("Need to specify either production, decay, or computeprop!")
                
                
                
                if self.MELAsettings["computeprop"] and (self.MELAsettings["prod"] or self.MELAsettings["dec"]):
                        proabilityprop = self.MELA.getXPropagator(self.MELAsettings["propscheme"])
                        self.out.branch("LHEMela_"+prob["name"]+eventSide, "F")
                        self.out.fillBranch("LHEMela_"+prob["name"]+eventSide, proabilityprop)
                elif self.MELAsettings["computeprop"]:
                    probability = self.MELA.getXPropagator(self.MELAsettings["propscheme"])
                
                # self.out.fillBranch("LHEMela_nativeProb", nativeprob)
                self.out.branch("LHEMela_"+prob["name"]+eventSide, "F")
                self.out.fillBranch("LHEMela_"+prob["name"]+eventSide, probability)
            
        

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