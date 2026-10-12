from __future__ import print_function
import copy
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.HeppyCore.utils.deltar import deltaR
from ZZAnalysis.NanoAnalysis.initializeMELA import check_enum
from ZZAnalysis.NanoAnalysis.MELAProbHelper import MELAProbHelper
import Mela


class LHEAngProbFiller(Module):
    """Calculates angles and proabilities with LHE-level information. 
    MELA = Pointer to MELA passed from nanoZZ4lAnalysis.py 
    """
    
    def __init__(self, MELA, NANOVERSION, MELASettings = None):
        self.MELA = MELA
        self.NANOVERSION = NANOVERSION
        self.MELASettings = MELASettings
        self.sortedSettings = []
        self.ProbHelper = MELAProbHelper(self.MELA, self.MELASettings, "LHE")
        print("***LHEAngProbFiller: set for: ", self.ProbHelper.names, flush=True)

       
            
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("LHEMela_costheta1", "F", limitedPrecision=16, title="In the Higgs' rest frame, theta_1 is the angle between the momentum of Z1 and the momentum of one of its decay products.")
        self.out.branch("LHEMela_costheta2", "F", limitedPrecision=16, title="In the Higgs' rest frame, theta_2 is the angle between the momentum of Z2 and the momentum of one of its decay products.")
        self.out.branch("LHEMela_Phi", "F", limitedPrecision=16, title="In the Higgs' rest frame, phi is the angle between the planes formed by the decay products of the two Z bosons.")
        self.out.branch("LHEMela_costhetastar", "F", limitedPrecision=16, title="In the Higgs' rest frame, theta_star is the angle between the beamline and the momentum of one of the Higgs' decay products.")
        self.out.branch("LHEMela_Phi1", "F", limitedPrecision=16, title="In the Higgs' rest frame, phi_1 is the angle between the decay plane of Z1 and the beamline.")
        if self.MELASettings != None: 
            self.ProbHelper.bookProbs(wrappedOutputTree)


        
    def analyze(self, event):
        LHEPart = Collection(event, 'LHEPart')
        
        mothers = Mela.SimpleParticleCollection_t()
        daughters = Mela.SimpleParticleCollection_t()
        associated = Mela.SimpleParticleCollection_t()
        hMass = -999.


        ## Only run if mother-daughter associations are available, i.e. nanoAODv15 or newer. 
        if self.NANOVERSION >= 15: 
            LHEMothers = filter(lambda p: p.MELAStatus==1, LHEPart)
            LHEDaughters = filter(lambda p: p.MELAStatus==2, LHEPart)
            LHEAssociated = filter(lambda p: p.MELAStatus==3, LHEPart)

            status2 = filter(lambda p: p.status==2, LHEPart)
            for i, mp in enumerate(LHEMothers): 
                temp_particle = Mela.SimpleParticle_t(mp.pdgId, mp.pt, mp.eta, mp.phi, mp.mass, True)
                mothers.add_particle(temp_particle)

            for i, hp in enumerate(status2): 
                if hp.pdgId == 25: 
                        hMass = hp.mass

            for i, dp in enumerate(LHEDaughters): 
                temp_particle = Mela.SimpleParticle_t(dp.pdgId, dp.pt, dp.eta, dp.phi, dp.mass, True)
                daughters.add_particle(temp_particle)

            nDaughters = len(daughters.toList())
            for i, ap in enumerate(LHEAssociated): 
                temp_particle = Mela.SimpleParticle_t(ap.pdgId, ap.pt, ap.eta, ap.phi, ap.mass, True)
                if nDaughters == 0 and hMass < 0. and abs(ap.pdgId) in [11, 13, 15] : # For ZZ bkg samples, compute angles from the 4 leptons from Z decays
                    daughters.add_particle(temp_particle)
                else:
                    associated.add_particle(temp_particle)
                        
        elif self.NANOVERSION < 15: 
            # NANOAODv14 or older, using workaround"
            statusOneLepsIdx = []
            statusOneAssocIdx = []
            for i, lp in enumerate(LHEPart): 
                if lp.status == -1: 
                    mothers.add_particle(Mela.SimpleParticle_t(lp.pdgId, lp.pt, lp.eta, lp.phi, lp.mass, True))
                elif lp.status == 1: 
                    ## Build a list of all status one leptons and associated and handle collections at the end. 
                    if abs(lp.pdgId) in [11, 13, 15]:
                        statusOneLepsIdx.append(i)
                    else: 
                        statusOneAssocIdx.append(i)
                elif lp.status == 2 and lp.pdgId == 25: 
                        hMass = lp.mass
            
            ## Add last four leptons as daughters, remaining leptons are associated. 
            for i in statusOneLepsIdx[-4:]: 
                lp = LHEPart[i]
                temp_particle = Mela.SimpleParticle_t(lp.pdgId, lp.pt, lp.eta, lp.phi, lp.mass, True)
                daughters.add_particle(temp_particle)
            del statusOneLepsIdx[-4:]

            for i in statusOneLepsIdx:
                lp = LHEPart[i]
                temp_particle = Mela.SimpleParticle_t(lp.pdgId, lp.pt, lp.eta, lp.phi, lp.mass, True)
                associated.add_particle(temp_particle)
            
            for i in statusOneAssocIdx: 
                lp = LHEPart[i]
                temp_particle = Mela.SimpleParticle_t(lp.pdgId, lp.pt, lp.eta, lp.phi, lp.mass, True)
                associated.add_particle(temp_particle)

                
        else: 
            print("**LHEAngProbFiller: No version of NANOAOD specified!")
            
        #Check if selected 4-leps match the higgs: 
        if len(daughters.toList()) == 4:
            self.MELA.setInputEvent(daughters, None, None, 0)
            qH, mZ1, mZ2, costheta1, costheta2, Phi, costhetastar, Phi1 = self.MELA.computeDecayAngles()
            if (hMass > 0.) and (abs(hMass - daughters.MTotal()) > 0.01) : # Don't trigger this when there is no Higgs present. Some samples, such as ggF produced with MCFM-JHUGen, produce the four lepton system without the intermediate higgs. 
                print(f"WARNING: LHEAngProbFiller: The invariant mass of the four LHE-leptons, {daughters.MTotal()}, is  different from the mass of the LHE-Higgs {hMass}! Expected a difference of less than 0.01, obtained a difference of ", abs(hMass - daughters.MTotal()))
        else:
            qH, mZ1, mZ2, costheta1, costheta2, Phi, costhetastar, Phi1 = 0.,0.,0.,-999.,-999.,-999.,-999.,-999.
            print(f"WARNING: LHEAngProbFiller: {len(daughters.toList())} LHE-leptons were selected for this event (4 expected)!", len(associated.toList()))

        self.out.fillBranch("LHEMela_costheta1", costheta1)
        self.out.fillBranch("LHEMela_costheta2", costheta2)
        self.out.fillBranch("LHEMela_Phi", Phi)
        self.out.fillBranch("LHEMela_Phi1", Phi1)
        self.out.fillBranch("LHEMela_costhetastar", costhetastar)

        self.MELA.resetInputEvent()

        self.ProbHelper.fillProbs([daughters], [associated], [mothers])
        
        return True
    
