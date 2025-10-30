from __future__ import print_function
import copy
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from  ZZAnalysis.NanoAnalysis.initializeMELA import check_enum
from ZZAnalysis.NanoAnalysis.ZZExtraFiller import *
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
                if (prob["isgen"]): continue # Only calculate Reco-level probabilities in this module.

                ### Merge specific settings with defaults and check for unsupported values
                fprob=copy.deepcopy(defaults)
                for key, value in prob.items():
                    if key in defaults:
                        fprob[key] = value
                    else:
                        raise(ValueError(f"RecoProbFiller: unknown parameter {key} in {prob['Name']}"))

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

            print(f"***RecoProbFiller: probs: {names}", flush=True)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if len(self.sortedSettings) !=0 :
            for p, prob in enumerate(self.sortedSettings):
                self.out.branch(f"ZZCand_P_{prob['Name']}", "F", lenVar = "nZZCand", limitedPrecision=16, title="User-defined Reco-level probability")
                if prob.get("ispm4l", False):
                    self.out.branch("ZZCand_P_"+prob["Name"]+"_ScaleUp", "F", lenVar = "nZZCand", limitedPrecision=16, title="User-defined Reco-level m4l probability with Scale uncertainties up")
                    self.out.branch("ZZCand_P_"+prob["Name"]+"_ScaleDown", "F", lenVar = "nZZCand", limitedPrecision=16, title="User-defined Reco-level m4l probability with Scale uncertainties down")
                    self.out.branch("ZZCand_P_"+prob["Name"]+"_SystUp", "F", lenVar = "nZZCand", limitedPrecision=16, title="User-defined Reco-level m4l probability with Systematic uncertainties up")
                    self.out.branch("ZZCand_P_"+prob["Name"]+"_SystDown", "F", lenVar = "nZZCand", limitedPrecision=16, title="User-defined Reco-level m4l probability with Systematic uncertainties down")
                if prob["computeprop"]: 
                    self.out.branch("ZZCand_P_"+prob["Name"]+"_prop", "F", lenVar = "nZZCand", limitedPrecision=16, title="User-defined Reco-level probability with non-default propagator scheme")
                        

        
    def analyze(self, event):
        if len(self.sortedSettings) != 0 :
            cands = Collection(event, 'ZZCand')
            if len (cands) == 0 :
                return True
            electrons = Collection(event, "Electron")
            muons = Collection(event, "Muon")
            self.leps = list(electrons) + list(muons)
            fsrPhotons = Collection(event, "FsrPhoton")
            
            mothers = Mela.SimpleParticleCollection_t()
            daughters = Mela.SimpleParticleCollection_t()
            associated = Mela.SimpleParticleCollection_t()

            for iprob, prob in enumerate(self.sortedSettings):
                vprob = [0.]*len(self.sortedSettings) # to be used to retrieve denominators for probs whith dividep
                self.MELA.setInputEvent(daughters, None, None, 0)
                
                ### Parse MELA settings for the desired probability
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

                ### Configure MELA for the event. 
                self.MELA.setProcess(MELA_Process, MELA_MatrixElement, MELA_Production)
                self.MELA.differentiate_HWW_HZZ = MELA_separatewwzz
                self.MELA.setMelaLeptonInterference(MELA_leptoninterference)                    

                if MELA_matchMx: 
                    self.MELA.setMelaHiggsMassWidth(daughters.MTotal(), 0.00001, 0)
                    self.MELA.setMelaHiggsMassWidth(daughters.MTotal(), 0.00001, 1)
                    
                for coupl, coupl_val in prob["Couplings"].items(): 
                    setattr(self.MELA, coupl, coupl_val) # FIXME: This needs to be reset at the beginning of the loop!

                # Define an array to fill with a probability for each candidate
                probVec = [-999.]*len(cands)
                    
                # Do so again for special cases where additional output is needed
                if MELA_ispm4l: 
                    probVec_ScaleUp = [-999.]*len(cands)
                    probVec_ScaleDown = [-999.]*len(cands)
                    probVec_SystUp = [-999.]*len(cands)
                    probVec_SystDown = [-999.]*len(cands)
                    
                if MELA_computeprop and (MELA_prod or MELA_dec): 
                    probPropVec = [-999.]*len(cands)
                    
                # Lifted from ZZExtraFiller 
                for iCand, aCand in enumerate(cands):
                    theCandLepIdxs = [aCand.Z1l1Idx, aCand.Z1l2Idx, aCand.Z2l1Idx, aCand.Z2l2Idx]
                        
                    theCandLeps = [self.leps[i] for i in theCandLepIdxs]
                    dressedLepsp4 = [ZZExtraFiller.getDressedP4(self = None, lep = l, fsrPhotons=fsrPhotons) for l in theCandLeps]
                    daughters.add_particle(Mela.SimpleParticle_t(theCandLeps[0].pdgId, dressedLepsp4[0].Px(), dressedLepsp4[0].Py(), dressedLepsp4[0].Pz(), dressedLepsp4[0].E()))
                    daughters.add_particle(Mela.SimpleParticle_t(theCandLeps[1].pdgId, dressedLepsp4[1].Px(), dressedLepsp4[1].Py(), dressedLepsp4[1].Pz(), dressedLepsp4[1].E()))
                    daughters.add_particle(Mela.SimpleParticle_t(theCandLeps[2].pdgId, dressedLepsp4[2].Px(), dressedLepsp4[2].Py(), dressedLepsp4[2].Pz(), dressedLepsp4[2].E()))
                    daughters.add_particle(Mela.SimpleParticle_t(theCandLeps[3].pdgId, dressedLepsp4[3].Px(), dressedLepsp4[3].Py(), dressedLepsp4[3].Pz(), dressedLepsp4[3].E()))
                    
                    ### Reset the event and the default probability settings per probability to be calculated. 
                    if MELA_prod and MELA_dec: 
                        probVec[iCand] = self.MELA.computeProdDecP(MELA_useconstant)
                    elif MELA_prod: 
                        probVec[iCand] = self.MELA.computeProdP(MELA_useconstant)
                    elif MELA_dec:
                        probVec[iCand] = self.MELA.computeP(MELA_useconstant)
                    elif MELA_ispm4l: 
                        probVec[iCand] = self.MELA.computePM4L(Mela.SuperMelaSyst.SMSyst_None)
                        probVec_ScaleUp[iCand] = self.MELA.computePM4l(Mela.SuperMelaSyst.SMSyst_ScaleUp)
                        probVec_ScaleDown[iCand] = self.MELA.computePM4l(Mela.SuperMelaSyst.SMSyst_ScaleDown)
                        probVec_SystUp[iCand] = self.MELA.computePM4l(Mela.SuperMelaSyst.SMSyst_ResUp)
                        probVec_SystDown[iCand] = self.MELA.computePM4l(Mela.SuperMelaSyst.SMSyst_ResDown)
                    else:
                        raise KeyError(f"RecoProbFiller: {MELA_Name}: Need to specify either production, decay, pm4l, or computeprop!") # FIXME: then the next 2 cases should come before this else? 
                        
                    if MELA_computeprop and (MELA_prod or MELA_dec): 
                        probPropVec[iCand] = self.MELA.getXPropagator(MELA_propscheme)
                    elif MELA_computeprop: 
                        probVec[iCand] = self.MELA.getXPropagator(MELA_propscheme)

                    # Handle divideP
                    vprob[iprob] = probVec[iCand]
                    if MELA_divideP_idx != -1:
                        den = vprob[MELA_divideP_idx] #Because of prob sorting, this has already been computed
                        if den != 0: 
                            probVec[iCand] /= den
                        else: 
                            print("**RecoProbFiller: Protecting against division by 0.")
                            
                    ### Reset the input event on a per-candidate basis. 
                    self.MELA.resetInputEvent()

                    # Write out all branches for this probability
                    self.out.fillBranch("ZZCand_P_"+prob["Name"], probVec)

                    if MELA_ispm4l: 
                        self.out.fillBranch("ZZCand_P_"+prob["Name"]+"_ScaleUp", probVec_ScaleUp)
                        self.out.fillBranch("ZZCand_P_"+prob["Name"]+"_ScaleDown", probVec_ScaleDown)
                        self.out.fillBranch("ZZCand_P_"+prob["Name"]+"_SystUp", probVec_SystUp)
                        self.out.fillBranch("ZZCand_P_"+prob["Name"]+"_SystDown", probVec_SystDown)
                    
                    if MELA_computeprop and (MELA_prod or MELA_dec): 
                        self.out.fillBranch("ZZCand_P_"+prob["Name"]+"_prop", probPropVec)
        
        return True
    
