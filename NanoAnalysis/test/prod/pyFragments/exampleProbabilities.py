#An example pyfragment for computing mela probabilities in the NanoAnalysis framework. In general, each statement of setConf("probabilities", {} ... ) should add a mela probability to the list named "probabilities". 
#Listed below are all variables that can be set when computing a probability: 
# - Name: String. The name of the probability to be calculated. When the output root file is written, the name of the branch will be "LHEMela_<Name>" for lhe-level probabilities. 
# 
# - Process: The process run by MELA. This typically denotes the spin of the resonance if a JHUGen Matrix Element is used or whether or not the process is signal or background if an MCFM Matrix Element is used.
#            The possible values for Process are found here: https://spin.pha.jhu.edu/MELA/tvar_enums.html#proc_enum 
#
# - MatrixElement: The generator from which the matrix element used to evaluate the probability is taken. The possible values are listed here: https://spin.pha.jhu.edu/MELA/tvar_enums.html#matel_enum
#
# - Production: The physics process to be evaluated in the given sample. The possible values are listed here: https://spin.pha.jhu.edu/MELA/tvar_enums.html#prod_enum
#
# - Couplings: The values of the MELA couplings to be used in the JHUGen Amplitude Basis. Possible couplings to set are listed here: https://spin.pha.jhu.edu/MELA/MELA_couplings_table.html
#             Couplings is formatted as a dictionary with keys corresponding to the desired coupling to be set. For each desired coupling, the input is a list in the format [Re, Im], where Re and Im are the real and imaginary components of the coupling. 
#
# - Prod: Boolean. Indicates whether or not the probability will be computed for the production side.
# - Dec: Boolean. Indicates whether or not the probability will be computed for the decay side.  
#   - Listed are the possible combinations for Prod and Dec and what they achieve: 
#       - Prod=True, Dec=False: calls computeProdP(), which gives the production side probability. 
#       - Prod=False, Dec=True: calls computeP(), which gives decay side probability. 
#       - Prod=True, Dec=True: calls computeProdDecP(), which gives the combined production and decay probability. 
#
# - context: String. When "LHE", will compute LHE-level probabilities, using particles from the LHEPart collection. When "Reco", will compute Reco-level probabilities using the ZZCands. When "Any", will compute the same probability at both levels.
#
# - computeprop:  Boolean. When true, will allow for a different propagator scheme to be used. This should be done in conjunction with the variable propscheme. 
# - propscheme:  The propagator scheme that defines resonancnes. By default, set to FixedWidth. The values can be found here: https://spin.pha.jhu.edu/MELA/tvar_enums.html#reso_enum
#
# - decaymode: The decay mode of the candidates being analyzed. By default, this is set to CandidateDecayMode_ZZ. All decaymodes can be found here: https://spin.pha.jhu.edu/MELA/tvar_enums.html#decmode_enum
#
# - "separatewwzz": Boolean. Treats the WW and ZZ couplings specified in MELA as separate. For instance, with this separatewwzz=False, setting the coupling ghz1 will also set ghw1 to the same value. If separatewwzz=True, then changing ghz1 would not change ghw1. 

# - "useconstant": Boolean. This turns on the calculation of a corrective constant to different probabilities through Mela::getConstant. If you would like the "pure" MELA calculation to be run, set useConstant to false. By default true.
 
# - "match_mX": Boolean. If true, will set the Higgs mass to match the invariant mass of the daughter particles in each event. 
# - "lepton_interference":"DefaultLeptonInterf",
# - "ispm4l": Boolean. When True, causes a call to computePM4l(), which is a probability useful in calculating a signal-background discrimimant. 
# - "dividep": Take the "Name" variable of one probability and divide all probabilities with divdep=True by this probability. This is typically most useful for normalizing to a native probability. Supported only for context="LHE".
# - "addPAux": also store the result of MELA.getPAux() - cf. slide 5 of https://indico.cern.ch/event/1619963/#11-hzz4l-legacy-signal-strengh
# - "addPmavjj": also store the result of computeDijetConvBW(False) - cf. slide 8 of https://indico.cern.ch/event/1619963/contributions/6826626/attachments/3190935/5678923/STXSupdate.pdf
# - "addPmavjj_true": also store the result of computeDijetConvBW(False)


#As a couple of general notes: 
# - append=True should always be set for each instance of setConf. 
# - Each variable that takes a MELA enumeration is set as a string. The string must EXACTLY match a given enumeration linked in the documentation. 


#The first probability calculates the lhe-level (isgen=True), decay-side probability (Dec=True) for a gluon fusion to Higgs sample (Production=ZZGG) decayed in JHUGen (MatrixElement=JHUGen) to ZZ {by default, decaymode=CandiateDecayMode_ZZ}.
# The final probability will be named "LHEMela_P_Native" ("Name= 'Native'")
#


from ZZAnalysis.NanoAnalysis.tools import setConf
setConf("probabilities", {'Name': "Native", 
                          "Process": "SelfDefine_spin0", 
                          "MatrixElement":"JHUGen",
                          "Production": "ZZGG", 
                          "Couplings": {'ghz1':[1,0], 
                                        'ghg2':[1,0]}, 
                          "Prod": False,
                          "Dec": True, 
                          "context": "LHE", 
                          "computeprop": False 
                          }, 
                          append=True)

setConf("probabilities", {'Name': "ggH_ghg2_1_ghz4_1", 
                        "Process": "SelfDefine_spin0", 
                        "MatrixElement":"JHUGen",
                        "Production": "ZZGG", 
                        "Couplings": {'ghz4':[1,0], 
                                        'ghg2':[1,0]}, 
                        "Prod": False,
                        "Dec": True, 
                        "context": "LHE",
                        "dividep": "Native",
                        "computeprop": False }, append=True)
