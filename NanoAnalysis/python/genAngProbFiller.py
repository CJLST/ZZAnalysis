from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.HeppyCore.utils.deltar import deltaR
import ZZAnalysis.NanoAnalysis.initializeMELA
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
        for i, prob in enumerate(self.MELAsettings): 
            self.out.branch(f"LHEMela_{prob['Name']}", "F")
    
        
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
        elif self.NANOVERSION < 15: 
            # print("genAngProbFiller: NANOAODv14 or older, using workaround")
            for i, lp in enumerate(LHEPart): 
                temp_particle = Mela.SimpleParticle_t(lp.pdgId, lp.pt, lp.eta, lp.phi, lp.mass, True)
                if lp.status == -1: 
                    mothers.add_particle(temp_particle)
                elif lp.status == 1: 
                    if abs(lp.pdgId) in [11, 13] and i >= len(LHEPart) - 4:
                        daughters.add_particle(temp_particle)
                    elif i < len(LHEPart) - 4: 
                        associated.add_particle(temp_particle)
                elif lp.status == 2: 
                    if lp.pdgId == 25: 
                        higgs = temp_particle
                        hMass = lp.mass
                    else:
                        continue
                else: 
                    continue
        else: 
            print("**genAngProbFiller: No version of NANOAOD specified!")
        
                
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
            self.MELA.resetInputEvent()
        else: 
            print("Selected 4 leptons too different from H-mass in LHE: ", hMass - daughters.MTotal())
            self.MELA.resetInputEvent()

        if self.MELAsettings != None: 
            self.MELA.setInputEvent(daughters, associated, mothers, 1)
            setupInputs = {
                "Name": "Default_SHOULDBERENAMED",
                "Process": None, 
                "MatrixElement": None, 
                "Production":None, 
                "Prod": None, 
                "Dec": None, 
                "Couplings": None, 
                "isgen": None, 
                "computeprop": None, 
                "propscheme": Mela.ResonancePropagatorScheme.FixedWidth, 
                "decaymode": Mela.CandidateDecayMode.CandidateDecay_ZZ,
                "separatewwzz":False,
                "useconstant":False,
                "match_mX":False,
                "lepton_interference":Mela.LeptonInterference.DefaultLeptonInterf,
                "ispm4l": None,
                "dividep": None 
            }
            

            for p, prob in enumerate(self.MELAsettings):
                ### Only calculate LHE-level probabilities in this function. 
                if prob["isgen"]:  
                    ### Parse MELA settings for the desired probability and overwrite setup for all given parameters:
                    for key, val in prob.items(): setupInputs[key] = prob[key]

                    ### Define everything 
                    MELA_Name = setupInputs["Name"]
                    MELA_Process = initializeMELA.check_enum(setupInputs["Process"], Mela.Process)
                    MELA_MatrixElement = initializeMELA.check_enum(setupInputs["MatrixElement"], Mela.MatrixElement)
                    MELA_Production = initializeMELA.check_enum(setupInputs["Production"], Mela.Production)
                    MELA_prod = setupInputs["Prod"]
                    MELA_dec = setupInputs["Dec"]
                    MELA_computeprop = setupInputs["computeprop"]
                    MELA_propscheme = initializeMELA.check_enum(setupInputs["propscheme"], Mela.ResonancePropagatorScheme)
                    MELA_decaymode = initializeMELA.check_enum(setupInputs["decaymode"], Mela.CandidateDecayMode)
                    MELA_separatewwzz = setupInputs["separatewwzz"]
                    MELA_useconstant = setupInputs["useconstant"]
                    MELA_matchMx = setupInputs["match_mX"]
                    MELA_leptoninterference = initializeMELA.check_enum(setupInputs["lepton_interference"], Mela.LeptonInterference)
                    MELA_ispm4l = setupInputs["ispm4l"]
                    MELA_divideP = setupInputs["dividep"]


                    ### Configure MELA for the event. 
                    self.MELA.setProcess(MELA_Process, MELA_MatrixElement, MELA_Production)
                    self.MELA.differentiate_HWW_HZZ = MELA_separatewwzz
                    self.MELA.setCandidateDecayMode(MELA_decaymode)
                    self.MELA.setMelaLeptonInterference(MELA_leptoninterference)

                    if MELA_matchMx: 
                        self.MELA.setMelaHiggsMassWidth(daughters.MTotal(), 0.00001, 0)
                        self.MELA.setMelaHiggsMassWidth(daughters.MTotal(), 0.00001, 1)

                    for coupl, coupl_val in setupInputs["Couplings"]: 
                        setattr(self.MELA, coupl, coupl_val)

        
                    if MELA_prod and MELA_dec: 
                        probability = self.MELA.computeProdDecP(MELA_useconstant)
                    elif MELA_prod: 
                        probability = self.MELA.computeProdP(MELA_useconstant)
                    elif MELA_dec:
                        probability = self.MELA.computeP(MELA_useconstant)
                    elif MELA_ispm4l: 
                        probability = self.MELA.computePM4L(Mela.SuperMelaSyst.SMSyst_None)
                        probability_ScaleUp = self.MELA.computePM4l(Mela.SuperMelaSyst.SMSyst_ScaleUp)
                        probability_ScaleDown = self.MELA.computePM4l(Mela.SuperMelaSyst.SMSyst_ScaleDown)
                        probability_SystUp = self.MELA.computePM4l(Mela.SuperMelaSyst.SMSyst_ResUp)
                        probability_SystDown = self.MELA.computePM4l(Mela.SuperMelaSyst.SMSyst_ResDown)

                        self.out.branch("LHEMela_"+setupInputs["Name"]+"_ScaleUp", "F")
                        self.out.fillBranch("LHEMela_"+setupInputs["Name"]+"_ScaleUp", probability_ScaleUp)

                        self.out.branch("LHEMela_"+setupInputs["Name"]+"_ScaleDown", "F")
                        self.out.fillBranch("LHEMela_"+setupInputs["Name"]+"_ScaleDown", probability_ScaleDown)

                        self.out.branch("LHEMela_"+setupInputs["Name"]+"_SystUp", "F")
                        self.out.fillBranch("LHEMela_"+setupInputs["Name"]+"_SystUp", probability_SystUp)

                        self.out.branch("LHEMela_"+setupInputs["Name"]+"_SystDown", "F")
                        self.out.fillBranch("LHEMela_"+setupInputs["Name"]+"_SystDown", probability_SystDown)


                    else:
                        raise KeyError("Need to specify either production, decay, pm4l, or computeprop!")
                    
                    if MELA_computeprop and (MELA_prod or MELA_dec): 
                        probabilityprop = self.MELA.getXPropagator(MELA_propscheme)
                        self.out.branch("LHEMela_"+setupInputs["Name"]+"_prop", "F")
                        self.out.fillBranch("LHEMela_"+setupInputs["Name"]+"_prop", probabilityprop)
                    elif MELA_computeprop: 
                        probability = self.MELA.getXPropagator(MELA_propscheme)
                    

                    #TODO: Figure out a way to ensure that denominator probability is always computed before the probability to be normalized. 
                    if MELA_divideP is not None:
                        if MELA_divideP == MELA_Name: 
                            denominator = probability
                        elif MELA_divideP != MELA_Name: 
                            probability /= denominator


                
                
                
                
                
                    # self.out.fillBranch("LHEMela_nativeProb", nativeprob)
                    # self.out.branch("LHEMela_"+prob["name"], "F")
                    self.out.fillBranch("LHEMela_"+prob["name"], probability)
                    self.MELA.resetInputEvent()
        
       
        
        
        return True
    



