from __future__ import print_function
import copy
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from  ZZAnalysis.NanoAnalysis.initializeMELA import check_enum
from ZZAnalysis.NanoAnalysis.ZZExtraFiller import *
import Mela

class MELAProbHelper(): 
    """Class for handling computation of probabilities with MELA.
    MELA = MELA Object passed from nanoZZ4lAnalysis.py 
    MELASettings = dictionary with settings for the probabilities to be computed
    ModuleContext = The level of information with which the probability is calculated. Possible values are "LHE" (for LHE-level) and "Reco" (for reco-level). These are passed from the module calling the computation. 
    A probability with context defined as "Any" will be computed at lhe and reco level.  
    """

    def __init__(self, MELA, MELASettings, ModuleContext):
        self.MELA = MELA
        self.sortedSettings = []
        self.ModuleContext = ModuleContext
        self.names = []
        
        if ModuleContext not in ["LHE", "Reco"] :
            raise valueError("MELAProbHelper: invalid ModuleContext", ModuleContext)

        if MELASettings != None:
            defaults = {
                "Name": "Default_You_Should_Rename_This",
                "Process": None, 
                "MatrixElement": None, 
                "Production":None, 
                "Prod": None, 
                "Dec": None, 
                "Couplings": None, 
                "context": None, 
                "computeprop": None, 
                "propscheme": "FixedWidth", 
                "decaymode": "CandidateDecay_ZZ",
                "separatewwzz":False,
                "useconstant":False,
                "addPAux":False,
                "match_mX":False,
                "lepton_interference":"DefaultLeptonInterf",
                "ispm4l": None,
                'addPmavjj': False,
                'addPmavjj_true': False,
                "dividep": None
            }
        

            ### Split into reco and LHE probs
            for prob in MELASettings:
                if prob["context"] not in ["LHE", "Reco", "Any"] :
                    raise valueError(f"MELAProbHelper: invalid context {prob['context']} for {prob['Name']}")
                    
                if (prob["context"] != ModuleContext) and (prob["context"] != "Any"): continue 

                ### Merge specific settings with defaults and check for unsupported values
                fprob=copy.deepcopy(defaults)
                for key, value in prob.items():
                    if key in defaults:
                        fprob[key] = value
                    else:
                        raise(ValueError(f"MELAProbHelper: unknown parameter {key} in {prob['Name']}"))

                # Add branch name so it does not need to be remade within loops
                if self.ModuleContext == "LHE": 
                    fprob["branchname"] = f"LHEMela_P_{prob['Name']}"
                else: 
                    fprob["branchname"] = f"ZZCand_P_{prob['Name']}"
                
                ### Sort the MELASettings dictionary so that all probabilities with divideP are last 
                if (fprob["dividep"]==None):
                    self.sortedSettings.insert(0,fprob)
                else: 
                    self.sortedSettings.append(fprob)

            ### Add index of probability to be used for dividep.
            self.names = [d["Name"] for d in self.sortedSettings]
            for prob in self.sortedSettings:
                dp = prob["dividep"]
                if dp == None :
                    prob["dividep_idx"] = -1
                elif ModuleContext != "LHE" :
                        raise(ValueError(f"MELAProbHelper: for {prob['Name']}: 'dividep' is supported only for 'context=LHE'"))                    
                else:
                    prob["dividep_idx"] = self.names.index(dp) # Index of probability to be used as denominator.
        
    def bookProbs(self, wrappedOutputTree): 
        #this needs a per-module lenVar. 
        self.out = wrappedOutputTree
        # print("***MELAHELPER: ", len(self.sortedSettings))
        if len(self.sortedSettings) !=0 :
            for p, prob in enumerate(self.sortedSettings):
                # print(prob["branchname"])
                if self.ModuleContext == "LHE":
                    lenVar_ = None
                else:
                    lenVar_ = "nZZCand"
                    
                self.out.branch(prob["branchname"], "F", lenVar=lenVar_, limitedPrecision=16, title=f"User-defined {self.ModuleContext}-level probability")
                if prob.get("ispm4l", False):
                    self.out.branch(prob["branchname"]+"_ScaleUp", "F", lenVar=lenVar_, limitedPrecision=16, title="User-defined {self.ModuleContext}-level m4l probability with Scale uncertainties up")
                    self.out.branch(prob["branchname"]+"_ScaleDown", "F", lenVar=lenVar_, limitedPrecision=16, title="User-defined {self.ModuleContext}-level m4l probability with Scale uncertainties down")
                    self.out.branch(prob["branchname"]+"_SystUp", "F", lenVar=lenVar_, limitedPrecision=16, title="User-defined {self.ModuleContext}-level m4l probability with Systematic uncertainties up")
                    self.out.branch(prob["branchname"]+"_SystDown", "F", lenVar=lenVar_, limitedPrecision=16, title="User-defined {self.ModuleContext}-level m4l probability with Systematic uncertainties down")
                if prob["computeprop"]: 
                    self.out.branch(prob["branchname"]+"_prop", "F", lenVar=lenVar_, limitedPrecision=16, title="User-defined weight to translate from POWHEG complex propagator scheme to JHUGen Breit-Wigner scheme")
                if prob["addPAux"]:
                    self.out.branch(prob["branchname"]+"_aux", "F", lenVar=lenVar_, limitedPrecision=16, title="User-defined auxiliary probability")
                if prob["addPmavjj"]:
                    self.out.branch(prob["branchname"]+"_mavjj", "F", lenVar=lenVar_, title="User-defined mavjj probability")
                if prob["addPmavjj_true"]:
                    self.out.branch(prob["branchname"]+"_mvajj_true", "F", lenVar=lenVar_, title="User-defined mvajj_true probability")
    def fillProbs(self, candDaughters, candAssociated, candMothers): 
        if len(self.sortedSettings) == 0: return True
        else: 
            vprob = [0.]*len(self.sortedSettings) # to be used to retrieve denominators for probs whith dividep for the current cand
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
                MELA_addPmavjj = prob["addPmavjj"]
                MELA_addPmavjj_true = prob["addPmavjj_true"]
                MELA_divideP_idx = prob["dividep_idx"]
                MELA_branchname = prob["branchname"]
                MELA_addPAux = prob["addPAux"]

                ### Configure MELA for the event. 
                self.MELA.setProcess(MELA_Process, MELA_MatrixElement, MELA_Production)
                self.MELA.differentiate_HWW_HZZ = MELA_separatewwzz
                self.MELA.setMelaLeptonInterference(MELA_leptoninterference)

                # Define arrays to fill with the probabilites for each candidate
                probVec = [-999.]*len(candDaughters)
                if MELA_ispm4l: 
                    probVec_ScaleUp = [-999.]*len(candDaughters)
                    probVec_ScaleDown = [-999.]*len(candDaughters)
                    probVec_SystUp = [-999.]*len(candDaughters)
                    probVec_SystDown = [-999.]*len(candDaughters)
                if MELA_computeprop: 
                    probPropVec = [-999.]*len(candDaughters)

                if MELA_addPAux:
                    probVec_PAux = [-999.]*len(candDaughters)

                if MELA_addPmavjj:
                    probVec_mavjj = [-999.]*len(candDaughters)
                if MELA_addPmavjj_true:
                    probVec_mvajj_true = [-999.]*len(candDaughters)

                # Compute prob for each candidate
                for iCand, aCand in enumerate(candDaughters):
                    
                    for coupl, coupl_val in prob["Couplings"].items(): 
                        setattr(self.MELA, coupl, coupl_val)

                    if MELA_matchMx: 
                        self.MELA.setMelaHiggsMassWidth(candDaughters[iCand].MTotal(), 0.00001, 0)
                        self.MELA.setMelaHiggsMassWidth(candDaughters[iCand].MTotal(), 0.00001, 1)

                    ### Reset the event and the default probability settings per probability to be calculated. 

                    if self.ModuleContext == "LHE":
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
                        probVec[iCand] = self.MELA.computePM4l(Mela.SuperMelaSyst.SMSyst_None)
                        probVec_ScaleUp[iCand] = self.MELA.computePM4l(Mela.SuperMelaSyst.SMSyst_ScaleUp)
                        probVec_ScaleDown[iCand] = self.MELA.computePM4l(Mela.SuperMelaSyst.SMSyst_ScaleDown)
                        probVec_SystUp[iCand] = self.MELA.computePM4l(Mela.SuperMelaSyst.SMSyst_ResUp)
                        probVec_SystDown[iCand] = self.MELA.computePM4l(Mela.SuperMelaSyst.SMSyst_ResDown)
                    elif MELA_computeprop==False :
                        raise KeyError(f"MELAProbHelper: need to specify either (production and/or decay) or pm4l or computeprop for {MELA_Name}")

                    if MELA_computeprop:
                        probPropVec[iCand] = self.MELA.getXPropagator(MELA_propscheme)


                    if MELA_addPAux:
                        probVec_PAux[iCand] = self.MELA.getPAux()

                    if MELA_addPmavjj:
                        probVec_mavjj[iCand] = self.MELA.computeDijetConvBW(False)
                    if MELA_addPmavjj_true:
                        probVec_mvajj_true[iCand] = self.MELA.computeDijetConvBW(True)

                    # Handle divideP
                    if self.ModuleContext == "LHE": 
                        vprob[iprob] = probVec[0] # Note: iCand is always 0 for LHE
                     
                        if MELA_divideP_idx != -1:
                            den = vprob[MELA_divideP_idx] 
                            if den != 0: 
                                probVec[0] /= den
                            else: 
                                print("***MELAProbHelper: Protecting against division by 0.")

                # Write out all branches for this probability
                if self.ModuleContext == "LHE": # Branches are not collections in this case
                    probVec=probVec[0]
                    if MELA_ispm4l: 
                        probVec_ScaleUp   = probVec_ScaleUp[0]
                        probVec_ScaleDown = probVec_ScaleDown[0]
                        probVec_SystUp    = probVec_SystUp[0]
                        probVec_SystDown  = probVec_SystDown[0]
                    if MELA_computeprop and (MELA_prod or MELA_dec): 
                        probPropVec = probPropVec[0]
                    if MELA_addPAux:
                        probVec_PAux = probVec_PAux[0]
                    if MELA_addPmavjj:
                        probVec_mavjj = probVec_mavjj[0]
                    if MELA_addPmavjj_true:
                        probVec_mvajj_true = probVec_mvajj_true[0]

                self.out.fillBranch(MELA_branchname, probVec)

                if MELA_ispm4l: 
                    self.out.fillBranch(MELA_branchname+"_ScaleUp", probVec_ScaleUp)
                    self.out.fillBranch(MELA_branchname+"_ScaleDown", probVec_ScaleDown)
                    self.out.fillBranch(MELA_branchname+"_SystUp", probVec_SystUp)
                    self.out.fillBranch(MELA_branchname+"_SystDown", probVec_SystDown)

                if MELA_computeprop and (MELA_prod or MELA_dec): 
                    self.out.fillBranch(MELA_branchname+"_prop", probPropVec)

                if MELA_addPAux:
                    self.out.fillBranch(MELA_branchname+"_aux", probVec_PAux)
                if MELA_addPmavjj:
                    self.out.fillBranch(MELA_branchname+"_mavjj", probVec_mavjj)
                if MELA_addPmavjj_true:
                    self.out.fillBranch(MELA_branchname+"_mvajj_true", probVec_mvajj_true)
                
        return True
        

                        

