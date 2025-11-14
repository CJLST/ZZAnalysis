from __future__ import print_function
import copy
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from  ZZAnalysis.NanoAnalysis.initializeMELA import check_enum
from ZZAnalysis.NanoAnalysis.ZZExtraFiller import *
import Mela

class MELAProbHelper(Module): 
    """Class for handling computation of probabilities with MELA.
    MELA = MELA Object passed from nanoZZ4lAnalysis.py 
    MELASettings = dictionary with settings for the probabilities to be computed
    isGen = The level of information with which the probability is calculated. Possible values are True (for LHE-level) and False (for reco-level).  
    """

    def __init__(self, MELA, MELASettings, isGen):
        self.MELA = MELA
        self.sortedSettings = []
        self.isGen = isGen

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
        

        ### Split into reco and LHE probs
        for prob in MELASettings:
            if (prob["isgen"] != self.isGen): continue 

            ### Merge specific settings with defaults and check for unsupported values
            fprob=copy.deepcopy(defaults)
            for key, value in prob.items():
                if key in defaults:
                    fprob[key] = value
                else:
                    raise(ValueError(f"MELAProbHelper: unknown parameter {key} in {prob['Name']}"))

            # Add branch name so it does not need to be remade within loops
            if self.isGen: 
                fprob["branchname"] = f"LHEMela_P_{prob['Name']}"
            else: 
                fprob["branchname"] = f"ZZCand_P_{prob['Name']}"
            
            ### Sort the MELASettings dictionary so that all probabilities with divideP are last 
            if (fprob["dividep"]==None):
                self.sortedSettings.insert(0,fprob)
            else: 
                self.sortedSettings.append(fprob)

        ### Add index of probability to be used for dividep.
        names = [d["Name"] for d in self.sortedSettings]
        for prob in self.sortedSettings:
            dp = prob["dividep"]
            prob["dividep_eval"] = (None if dp == None else f"aCand.P_{dp}") # string to be evaluated to extract denominator
            # prob["dividep_idx"] = (-1 if dp == None else names.index(dp)) # prob index, more efficient but would require keeping probs for all cands

        print(f"***MELAProbHelper: probs: {names}", flush=True)
        
    def bookProbs(self, wrappedOutputTree): 
        #this needs a per-module lenVar. 
        self.out = wrappedOutputTree
        if len(self.sortedSettings) !=0 :
            for p, prob in enumerate(self.sortedSettings):
                print(prob["branchname"])
                if self.isGen: 
                    self.out.branch(prob["branchname"], "F", lenVar = None, limitedPrecision=16, title="User-defined LHE-level probability")
                    if prob.get("ispm4l", False):
                        self.out.branch(prob["branchname"]+"_ScaleUp", "F", lenVar = None, limitedPrecision=16, title="User-defined Reco-level m4l probability with Scale uncertainties up")
                        self.out.branch(prob["branchname"]+"_ScaleDown", "F", lenVar = None, limitedPrecision=16, title="User-defined Reco-level m4l probability with Scale uncertainties down")
                        self.out.branch(prob["branchname"]+"_SystUp", "F", lenVar = None, limitedPrecision=16, title="User-defined Reco-level m4l probability with Systematic uncertainties up")
                        self.out.branch(prob["branchname"]+"_SystDown", "F", lenVar = None, limitedPrecision=16, title="User-defined Reco-level m4l probability with Systematic uncertainties down")
                    if prob["computeprop"]: 
                        self.out.branch(prob["branchname"]+"_prop", "F", limitedPrecision=16, title="User-defined weight to translate from POWHEG complex propagator scheme to JHUGen Breit-Wigner scheme")

                else: 
                    self.out.branch(prob["branchname"], "F", lenVar = "nZZCand", limitedPrecision=16, title="User-defined LHE-level probability")
                    if prob.get("ispm4l", False):
                        self.out.branch(prob["branchname"]+"_ScaleUp", "F", lenVar = "nZZCand", limitedPrecision=16, title="User-defined Reco-level m4l probability with Scale uncertainties up")
                        self.out.branch(prob["branchname"]+"_ScaleDown", "F", lenVar = "nZZCand", limitedPrecision=16, title="User-defined Reco-level m4l probability with Scale uncertainties down")
                        self.out.branch(prob["branchname"]+"_SystUp", "F", lenVar = "nZZCand", limitedPrecision=16, title="User-defined Reco-level m4l probability with Systematic uncertainties up")
                        self.out.branch(prob["branchname"]+"_SystDown", "F", lenVar = "nZZCand", limitedPrecision=16, title="User-defined Reco-level m4l probability with Systematic uncertainties down")
                    if prob["computeprop"]: 
                        self.out.branch(prob["branchname"]+"_prop", "F", lenVar = "nZZCand", limitedPrecision=16, title="User-defined weight to translate from POWHEG complex propagator scheme to JHUGen Breit-Wigner scheme")
    
    def fillProbs(self, candDaughters, candAssociated, candMothers): 
        if len(self.sortedSettings) == 0: return True
        else: 
            for iprob, prob in enumerate(self.sortedSettings):
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
                MELA_divideP_eval = prob["dividep_eval"]
                MELA_branchname = prob["branchname"]

                ### Configure MELA for the event. 
                self.MELA.setProcess(MELA_Process, MELA_MatrixElement, MELA_Production)
                self.MELA.differentiate_HWW_HZZ = MELA_separatewwzz
                self.MELA.setMelaLeptonInterference(MELA_leptoninterference)                    

                

                for coupl, coupl_val in prob["Couplings"].items(): 
                    setattr(self.MELA, coupl, coupl_val)

                # Define arrays to fill with the probabilites for each candidate
                probVec = [-999.]*len(candDaughters)
                if MELA_ispm4l: 
                    probVec_ScaleUp = [-999.]*len(candDaughters)
                    probVec_ScaleDown = [-999.]*len(candDaughters)
                    probVec_SystUp = [-999.]*len(candDaughters)
                    probVec_SystDown = [-999.]*len(candDaughters)
                if MELA_computeprop and (MELA_prod or MELA_dec): 
                    probPropVec = [-999.]*len(candDaughters)

                # Compute prob for each candidate
                for iCand, aCand in enumerate(candDaughters):
                    vprob = [0.]*len(self.sortedSettings) # to be used to retrieve denominators for probs whith dividep for the current cand

                    if MELA_matchMx: 
                        self.MELA.setMelaHiggsMassWidth(candDaughters[iCand].MTotal(), 0.00001, 0)
                        self.MELA.setMelaHiggsMassWidth(candDaughters[iCand].MTotal(), 0.00001, 1)

                    ### Reset the event and the default probability settings per probability to be calculated. 

                    if self.isGen:
                        self.MELA.setInputEvent(candDaughters[iCand], candAssociated[iCand], candMothers[iCand], 1)
                    else: 
                        self.MELA.setInputEvent(candDaughters[iCand], candAssociated[iCand], None, 0)

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
                    elif MELA_computeprop==False :
                        raise KeyError(f"MELAProbHelper: need to specify either (production and/or decay) or pm4l or computeprop for {MELA_Name}")

                    if MELA_computeprop and not MELA_ispm4l:
                        if (MELA_prod or MELA_dec): 
                            probPropVec[iCand] = self.MELA.getXPropagator(MELA_propscheme)
                        elif MELA_computeprop: 
                            probVec[iCand] = self.MELA.getXPropagator(MELA_propscheme)

                    # Handle divideP
                    vprob[iprob] = probVec[iCand]
                    if MELA_divideP_eval != None:
                        den = eval(MELA_divideP_eval) #Because of prob sorting, this has already been stored
                        if den != 0: 
                            probVec[iCand] /= den
                        else: 
                            print("**MELAProbHelper: Protecting against division by 0.")

                # Write out all branches for this probability
                if self.isGen: 
                    self.out.fillBranch(MELA_branchname, probVec[0])
                    if MELA_ispm4l: 
                        self.out.fillBranch(MELA_branchname+"_ScaleUp", probVec_ScaleUp[0])
                        self.out.fillBranch(MELA_branchname+"_ScaleDown", probVec_ScaleDown[0])
                        self.out.fillBranch(MELA_branchname+"_SystUp", probVec_SystUp[0])
                        self.out.fillBranch(MELA_branchname+"_SystDown", probVec_SystDown[0])
                    if MELA_computeprop and (MELA_prod or MELA_dec): 
                        self.out.fillBranch(MELA_branchname+"_prop", probPropVec[0])
                else: 
                    self.out.fillBranch(MELA_branchname, probVec)

                    if MELA_ispm4l: 
                        self.out.fillBranch(MELA_branchname+"_ScaleUp", probVec_ScaleUp)
                        self.out.fillBranch(MELA_branchname+"_ScaleDown", probVec_ScaleDown)
                        self.out.fillBranch(MELA_branchname+"_SystUp", probVec_SystUp)
                        self.out.fillBranch(MELA_branchname+"_SystDown", probVec_SystDown)
                            
                    if MELA_computeprop and (MELA_prod or MELA_dec): 
                        self.out.fillBranch(MELA_branchname+"_prop", probPropVec)
                
                
        return True
        

                        

