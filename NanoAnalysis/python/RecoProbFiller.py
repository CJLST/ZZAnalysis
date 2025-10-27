from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from  ZZAnalysis.NanoAnalysis.initializeMELA import check_enum
from ZZAnalysis.NanoAnalysis.ZZExtraFiller import *
import Mela


class RecoProbFiller(Module):
    """Calculates proabilities with Reco-level information. 
    MELA = Pointer to MELA passed from nanoZZ4lAnalysis.py 
    """
    
    def __init__(self, MELA, NANOVERSION, settingsDict = None):
        print("***RecoProbFiller", flush=True)
        self.MELA = MELA
        self.MELAsettings = settingsDict
        self.NANOVERSION = NANOVERSION
            
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.MELAsettings != None: 
            self.sortedSettings = []
            self.denominator_name = ""
            for p, prob in enumerate(self.MELAsettings):
                if (prob["isgen"]) : continue 
                ### Sort the MELASettings dictionary so that all probabilities with divideP are last 
                if ("dividep" in prob) : 
                    self.sortedSettings.append(prob)
                    self.denominator_name = prob["dividep"]
                else: 
                    self.sortedSettings.insert(0,prob)
                self.out.branch(f"ZZCand_RecoMela_{prob['Name']}", "F", lenVar = "nZZCand", limitedPrecision=16, title="User-defined Reco-level probability")
                if prob.get("ispm4l", False): 
                    self.out.branch("ZZCand_RecoMela_"+prob["Name"]+"_ScaleUp", "F", lenVar = "nZZCand", limitedPrecision=16, title="User-defined Reco-level m4l probability with Scale uncertainties up")
                    self.out.branch("ZZCand_RecoMela_"+prob["Name"]+"_ScaleDown", "F", lenVar = "nZZCand", limitedPrecision=16, title="User-defined Reco-level m4l probability with Scale uncertainties down")
                    self.out.branch("ZZCand_RecoMela_"+prob["Name"]+"_SystUp", "F", lenVar = "nZZCand", limitedPrecision=16, title="User-defined Reco-level m4l probability with Systematic uncertainties up")
                    self.out.branch("ZZCand_RecoMela_"+prob["Name"]+"_SystDown", "F", lenVar = "nZZCand", limitedPrecision=16, title="User-defined Reco-level m4l probability with Systematic uncertainties down")
                if prob["computeprop"]: 
                    self.out.branch("ZZCand_RecoMela_"+prob["Name"]+"_prop", "F", lenVar = "nZZCand", limitedPrecision=16, title="User-defined Reco-level probability with non-default propagator scheme")
                        

        
    def analyze(self, event):
        # Only run if we have a pyfragment 
        if self.MELAsettings != None: 
            cands = Collection(event, 'ZZCand')
            electrons = Collection(event, "Electron")
            muons = Collection(event, "Muon")
            self.leps = list(electrons) + list(muons)
            fsrPhotons = Collection(event, "FsrPhoton")
            
            mothers = Mela.SimpleParticleCollection_t()
            daughters = Mela.SimpleParticleCollection_t()
            associated = Mela.SimpleParticleCollection_t()

            


            for p, prob in enumerate(self.sortedSettings):
                ### Only calculate Reco-level probabilities in this function. 
                if prob["isgen"] == False: 
                    ### Define parameters for the probability to be computed. They don't need to be changed on a per-candidate level. 
                    self.MELA.setInputEvent(daughters, None, None, 0)
                    setupInputs = {
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
                        
                    ### Parse MELA settings for the desired probability and overwrite setup for all given parameters:
                    for key, val in prob.items(): setupInputs[key] = prob[key]

                    

                    ### Define everything 
                    MELA_Name = setupInputs["Name"]
                    MELA_Process = check_enum(setupInputs["Process"], Mela.Process)
                    MELA_MatrixElement = check_enum(setupInputs["MatrixElement"], Mela.MatrixElement)
                    MELA_Production = check_enum(setupInputs["Production"], Mela.Production)
                    MELA_prod = setupInputs["Prod"]
                    MELA_dec = setupInputs["Dec"]
                    MELA_computeprop = setupInputs["computeprop"]
                    MELA_propscheme = check_enum(setupInputs["propscheme"], Mela.ResonancePropagatorScheme)
                    MELA_separatewwzz = setupInputs["separatewwzz"]
                    MELA_useconstant = setupInputs["useconstant"]
                    MELA_matchMx = setupInputs["match_mX"]
                    MELA_leptoninterference = check_enum(setupInputs["lepton_interference"], Mela.LeptonInterference)
                    MELA_ispm4l = setupInputs["ispm4l"]
                    MELA_divideP = setupInputs["dividep"]

                    ### Configure MELA for the event. 
                    self.MELA.setProcess(MELA_Process, MELA_MatrixElement, MELA_Production)
                    self.MELA.differentiate_HWW_HZZ = MELA_separatewwzz
                    self.MELA.setMelaLeptonInterference(MELA_leptoninterference)

                    

                    if MELA_matchMx: 
                        self.MELA.setMelaHiggsMassWidth(daughters.MTotal(), 0.00001, 0)
                        self.MELA.setMelaHiggsMassWidth(daughters.MTotal(), 0.00001, 1)
                    
                    for coupl, coupl_val in setupInputs["Couplings"].items(): 
                        setattr(self.MELA, coupl, coupl_val)

                    # Define an array to fill with a probability for each candidate
                    probVec = [-999.]*len(cands)
                    

                    if MELA_Name == self.denominator_name: 
                        denomVec = [-999.]*len(cands)
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
                            raise KeyError("Need to specify either production, decay, pm4l, or computeprop!")
                        
                        if MELA_computeprop and (MELA_prod or MELA_dec): 
                            probPropVec[iCand] = self.MELA.getXPropagator(MELA_propscheme)

                        elif MELA_computeprop: 
                            probVec[iCand] = self.MELA.getXPropagator(MELA_propscheme)
                        


                        # Handling divideP 
                        if MELA_Name == self.denominator_name: 
                            denomVec[iCand] = probVec[iCand]
                        
                        if MELA_divideP != None: 
                            if denomVec[iCand] != 0: 
                                probVec[iCand] /= denomVec[iCand]
                            else: 
                                print("**RecoProbFiller: Protecting against division by 0.")



                        
                        
                        
                        
                        

                            
                        ### Reset the input event on a per-candidate basis. 
                        self.MELA.resetInputEvent()

                    # Write out all branches for this probability
                    self.out.fillBranch("ZZCand_RecoMela_"+prob["Name"], probVec)

                    if MELA_ispm4l: 
                        self.out.fillBranch("ZZCand_RecoMela_"+setupInputs["Name"]+"_ScaleUp", probVec_ScaleUp)
                        self.out.fillBranch("ZZCand_RecoMela_"+setupInputs["Name"]+"_ScaleDown", probVec_ScaleDown)
                        self.out.fillBranch("ZZCand_RecoMela_"+setupInputs["Name"]+"_SystUp", probVec_SystUp)
                        self.out.fillBranch("ZZCand_RecoMela_"+setupInputs["Name"]+"_SystDown", probVec_SystDown)
                    
                    if MELA_computeprop and (MELA_prod or MELA_dec): 
                        self.out.fillBranch("ZZCand_RecoMela_"+setupInputs["Name"]+"_prop", probPropVec)



       
        
        
        return True
    
