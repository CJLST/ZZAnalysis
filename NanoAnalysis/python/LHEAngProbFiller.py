from __future__ import print_function
import copy
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.HeppyCore.utils.deltar import deltaR
from  ZZAnalysis.NanoAnalysis.initializeMELA import check_enum
import Mela


class LHEAngProbFiller(Module):
    """Calculates angles and proabilities with LHE-level information. 
    MELA = Pointer to MELA passed from nanoZZ4lAnalysis.py 
    """
    
    def __init__(self, MELA, NANOVERSION, MELASettings = None):
        self.MELA = MELA
        self.NANOVERSION = NANOVERSION
        self.sortedSettings = []

        if MELASettings != None:
            defaults = {
                "Name": "Default_You_Should_Rename_This",
                "Process": None, 
                "MatrixElement": None, 
                "Production":None, 
                "Prod": None, 
                "Dec": None, 
                "Couplings": None, 
                "isgen": None, 
                "computeprop": None, 
                "propscheme": "FixedWidth", 
                "decaymode": "CandidateDecay_ZZ",
                "separatewwzz":False,
                "useconstant":False,
                "match_mX":False,
                "lepton_interference":"DefaultLeptonInterf",
                "ispm4l": None,
                "dividep": None
            }
            
            for prob in MELASettings:
                if (not prob["isgen"]): continue # Only calculate LHE-level probabilities in this module.

                ### Merge specific settings with defaults and check for unsupported values
                fprob=copy.deepcopy(defaults)
                for key, value in prob.items():
                    if key in defaults:
                        fprob[key] = value
                    else:
                        raise(ValueError(f"LHEAngProbFiller: unknown parameter {key} in {prob['Name']}"))

                # Add branch name so it does not need to be remade within loops
                fprob["branchname"] = f"LHEMela_P_{prob['Name']}"

                ### Sort the MELASettings dictionary so that all probabilities with divideP are last 
                if (fprob["dividep"]==None):
                    self.sortedSettings.insert(0,fprob)
                else:
                    self.sortedSettings.append(fprob)

            ### Add index of probability to be used for dividep.
            names = [d["Name"] for d in self.sortedSettings]
            for prob in self.sortedSettings:
                dp = prob["dividep"]
                prob["dividep_idx"] = (-1 if dp == None else names.index(dp))

            print(f"***LHEAngProbFiller: probs: {names}", flush=True)
            
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("LHEMela_qH", "F", title="The mass of the Higgs candidate as reconstructed by the 4-leptons at LHE-level.")
        self.out.branch("LHEMela_mZ1", "F", title="The mass of the first decay particle as reconstructed by 2 of the LHE-level leptons.")
        self.out.branch("LHEMela_mZ2", "F", title="The mass of the second decay particle as reconstructed by 2 of the LHE-level leptons.")
        self.out.branch("LHEMela_costheta1", "F", limitedPrecision=16, title="In the Higgs' rest frame, theta_1 is the angle between the momentum of Z1 and the momentum of one of its decay products.")
        self.out.branch("LHEMela_costheta2", "F", limitedPrecision=16, title="In the Higgs' rest frame, theta_2 is the angle between the momentum of Z2 and the momentum of one of its decay products.")
        self.out.branch("LHEMela_Phi", "F", limitedPrecision=16, title="In the Higgs' rest frame, phi is the angle between the planes formed by the decay products of the two Z bosons.")
        self.out.branch("LHEMela_costhetastar", "F", limitedPrecision=16, title="In the Higgs' rest frame, theta_star is the angle between the beamline and the momentum of one of the Higgs' decay products.")
        self.out.branch("LHEMela_Phi1", "F", limitedPrecision=16, title="In the Higgs' rest frame, phi_1 is the angle between the decay plane of Z1 and the beamline.")
        if len(self.sortedSettings) !=0 :
            for i, prob in enumerate(self.sortedSettings):
                self.out.branch(prob["branchname"], "F", limitedPrecision=16, title="User-defined LHE-level probability")
                if prob.get("ispm4l", False): 
                    self.out.branch(prob["branchname"]+"_ScaleUp", "F", limitedPrecision=16, title="User-defined LHE-level m4l probability with Scale uncertainties up")
                    self.out.branch(prob["branchname"]+"_ScaleDown", "F", limitedPrecision=16, title="User-defined LHE-level m4l probability with Scale uncertainties down")
                    self.out.branch(prob["branchname"]+"_SystUp", "F", limitedPrecision=16, title="User-defined LHE-level m4l probability with Systematic uncertainties up")
                    self.out.branch(prob["branchname"]+"_SystDown", "F", limitedPrecision=16, title="User-defined LHE-level m4l probability with Systematic uncertainties down")
                if prob["computeprop"]: 
                    self.out.branch(prob["branchname"]+"_prop", "F", limitedPrecision=16, title="User-defined LHE-level probability with non-default propagator scheme")

        
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

            status2 = filter(lambda p: p.status==2, LHEPart)
            for i, mp in enumerate(LHEMothers): 
                temp_particle = Mela.SimpleParticle_t(mp.pdgId, mp.pt, mp.eta, mp.phi, mp.mass, True)
                mothers.add_particle(temp_particle)
            
            for i, dp in enumerate(LHEDaughters): 
                temp_particle = Mela.SimpleParticle_t(dp.pdgId, dp.pt, dp.eta, dp.phi, dp.mass, True)
                daughters.add_particle(temp_particle)
            
            for i, ap in enumerate(LHEAssociated): 
                temp_particle = Mela.SimpleParticle_t(ap.pdgId, ap.pt, ap.eta, ap.phi, ap.mass, True)
                associated.add_particle(temp_particle)
            
            for i, hp in enumerate(status2): 
                temp_particle = Mela.SimpleParticle_t(hp.pdgId, hp.pt, hp.eta, hp.phi, hp.mass, True)
                if hp.pdgId == 25: 
                        higgs = temp_particle
                        hMass = hp.mass

        elif self.NANOVERSION < 15: 
            # print("genAngProbFiller: NANOAODv14 or older, using workaround")
            for i, lp in enumerate(LHEPart): 
                temp_particle = Mela.SimpleParticle_t(lp.pdgId, lp.pt, lp.eta, lp.phi, lp.mass, True)
                if lp.status == -1: 
                    mothers.add_particle(temp_particle)
                elif lp.status == 1: 
                    if i >= len(LHEPart) - 4: 
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
        if abs(hMass - daughters.MTotal()) < 0.01 and len(daughters.toList()) == 4:
            # self.MELA.setInputEvent(daughters, associated, mothers, 1)
            self.MELA.setInputEvent(daughters, None, None, 0)
            qH, mZ1, mZ2, costheta1, costheta2, Phi, costhetastar, Phi1 = self.MELA.computeDecayAngles()
            self.out.fillBranch("LHEMela_costheta1", costheta1)
            self.out.fillBranch("LHEMela_costheta2", costheta2)
            self.out.fillBranch("LHEMela_Phi", Phi)
            self.out.fillBranch("LHEMela_Phi1", Phi1)
            self.out.fillBranch("LHEMela_costhetastar", costhetastar)
        else: 
            if len(daughters.toList()) != 4: 
                print(f"WARNING: LHEAngProbFiller: {len(daughters.toList())} LHE-leptons were selected for this event (4 expected)!")
            elif abs(hMass - daughters.MTotal()) < 0.01: 
                print(f"WARNING: LHEAngProbFiller: The invariant mass of the four LHE-leptons, {daughters.MTotal()}, is too different from the mass of the LHE-Higgs {hMass}! Expected a difference of less than 0.01, obtained a difference of ", hMass - daughters.MTotal())

        self.MELA.resetInputEvent()

        if len(self.sortedSettings) != 0 and len(daughters.toList()) == 4:
            vprob = [0.]*len(self.sortedSettings)
            for iprob, prob in enumerate(self.sortedSettings):

                ### Reset the event and the default probability settings per probability to be calculated. 
                self.MELA.setInputEvent(daughters, associated, mothers, 1)

                ### Define everything 
                MELA_Name = prob["Name"]
                MELA_Process = check_enum(prob["Process"], Mela.Process)
                MELA_MatrixElement = check_enum(prob["MatrixElement"], Mela.MatrixElement)
                MELA_Production = check_enum(prob["Production"], Mela.Production)
                MELA_prod = prob["Prod"]
                MELA_dec = prob["Dec"]
                MELA_computeprop = prob["computeprop"]
                MELA_propscheme = check_enum(prob["propscheme"], Mela.ResonancePropagatorScheme)
                MELA_separatewwzz = prob["separatewwzz"]
                MELA_useconstant = prob["useconstant"]
                MELA_matchMx = prob["match_mX"]
                MELA_leptoninterference = check_enum(prob["lepton_interference"], Mela.LeptonInterference)
                MELA_ispm4l = prob["ispm4l"]
                MELA_divideP_idx = prob["dividep_idx"]
                MELA_branchname = prob["branchname"]

                ### Configure MELA for the event. 
                self.MELA.setProcess(MELA_Process, MELA_MatrixElement, MELA_Production)
                self.MELA.differentiate_HWW_HZZ = MELA_separatewwzz
                self.MELA.setMelaLeptonInterference(MELA_leptoninterference)
                

                if MELA_matchMx: 
                    self.MELA.setMelaHiggsMassWidth(daughters.MTotal(), 0.00001, 0)
                    self.MELA.setMelaHiggsMassWidth(daughters.MTotal(), 0.00001, 1)
                
                for coupl, coupl_val in prob["Couplings"].items(): 
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

                    self.out.fillBranch(MELA_branchname+"_ScaleUp", probability_ScaleUp)
                    self.out.fillBranch(MELA_branchname+"_ScaleDown", probability_ScaleDown)
                    self.out.fillBranch(MELA_branchname+"_SystUp", probability_SystUp)
                    self.out.fillBranch(MELA_branchname+"_SystDown", probability_SystDown)
                elif MELA_computeprop==False :
                    raise KeyError(f"LHEAngProbFiller: need to specify either (production and/or decay) or pm4l or computeprop for {MELA_Name}")
                
                if MELA_computeprop and not MELA_ispm4l:
                    if (MELA_prod or MELA_dec): 
                        probabilityprop = self.MELA.getXPropagator(MELA_propscheme)
                        self.out.fillBranch(MELA_branchname+"_prop", probabilityprop)
                    elif MELA_computeprop: 
                        probability = self.MELA.getXPropagator(MELA_propscheme)
                
                # Handle divideP 
                vprob[iprob] = probability
                if MELA_divideP_idx != -1:
                    den = vprob[MELA_divideP_idx]
                    if den != 0: 
                        probability /= den
                    else: 
                        print("**LHEAngProbFiller: Protecting against division by 0!")
                

                self.out.fillBranch(MELA_branchname, probability)
                self.MELA.resetInputEvent()
        
       
        
        
        return True
    
