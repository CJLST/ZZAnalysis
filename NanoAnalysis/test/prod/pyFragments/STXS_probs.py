from ZZAnalysis.NanoAnalysis.tools import setConf
setConf("probabilities", {
    "Name": "m4l_BKG",
    "Process": "bkgZZ", 
    "Production": "ZZGG", 
    "MatrixElement": "JHUGen", 
    "Couplings": {}, 
    "Prod": False, 
    "Dec": False, 
    "ispm4l": True, 
    "context": "Reco", 
    "computeprop": False, 
    "useconstant": True 
    }, append=True)

# https://github.com/CJLST/ZZAnalysis/blob/e1dd72d210d05e248873389c97da5ff9c9fcb94a/AnalysisStep/test/prod/pyFragments/RecoProbabilities.py#L156
setConf("probabilities", {
    "Name": "JVBF_SIG_ghv1_1_JHUGen_JECNominal",
    "Process": "HSMHiggs", 
    "Production": "JJVBF", 
    "MatrixElement": "JHUGen", 
    "Couplings": {},
    "Prod": True, 
    "Dec": False, 
    "ispm4l": False, 
    "context": "Reco", 
    "computeprop": False,
    "useconstant": True, # This argument affects the normalization, and should generally be set to true for reconstructed events, false forgenerator level events. [from https://spin.pha.jhu.edu/Manual.pdf]
    "addPAux": True
    }, append=True)

# https://github.com/CJLST/ZZAnalysis/blob/e1dd72d210d05e248873389c97da5ff9c9fcb94a/AnalysisStep/test/prod/pyFragments/RecoProbabilities.py#L162
setConf("probabilities", {
    "Name": "JJVBF_SIG_ghv1_1_JHUGen_JECNominal",
    "Process": "HSMHiggs", 
    "Production": "JJVBF", 
    "MatrixElement": "JHUGen", 
    "Couplings": {}, 
    "Prod": True, 
    "Dec": False,
    "ispm4l": False,
    "context": "Reco", 
    "computeprop": False, # TBC
    "useconstant": True # This argument affects the normalization, and should generally be set to true for reconstructed events, false forgenerator level events. [from https://spin.pha.jhu.edu/Manual.pdf]
    }, append=True)

# https://github.com/CJLST/ZZAnalysis/blob/e1dd72d210d05e248873389c97da5ff9c9fcb94a/AnalysisStep/test/prod/pyFragments/RecoProbabilities.py#L159
setConf("probabilities", {
    "Name": "JQCD_SIG_ghg2_1_JHUGen_JECNominal",
    "Process": "HSMHiggs", 
    "Production": "JQCD", 
    "MatrixElement": "JHUGen", 
    "Couplings": {},
    "Prod": True,
    "Dec": False,  
    "ispm4l": False, 
    "context": "Reco", 
    "computeprop": False, # TBC
    "useconstant": True # This argument affects the normalization, and should generally be set to true for reconstructed events, false forgenerator level events. [from https://spin.pha.jhu.edu/Manual.pdf]
    }, append=True)

# https://github.com/CJLST/ZZAnalysis/blob/e1dd72d210d05e248873389c97da5ff9c9fcb94a/AnalysisStep/test/prod/pyFragments/RecoProbabilities.py#L171
setConf("probabilities", {
    "Name": "HadWH_SIG_ghw1_1_JHUGen_JECNominal",
    "Process": "HSMHiggs", 
    "Production": "Had_WH", 
    "MatrixElement": "JHUGen", 
    "Couplings": {},
    "Prod": True,
    "Dec": False, 
    "ispm4l": False,
    "context": "Reco", 
    "computeprop": False, # TBC
    "useconstant": True, # This argument affects the normalization, and should generally be set to true for reconstructed events, false forgenerator level events. [from https://spin.pha.jhu.edu/Manual.pdf]
    "addPmavjj": True,
    "addPmavjj_true": True
    }, append=True)

# https://github.com/CJLST/ZZAnalysis/blob/e1dd72d210d05e248873389c97da5ff9c9fcb94a/AnalysisStep/test/prod/pyFragments/RecoProbabilities.py#L168
setConf("probabilities", {
    "Name": "HadZH_SIG_ghz1_1_JHUGen_JECNominal",
    "Process": "HSMHiggs", 
    "Production": "Had_ZH", 
    "MatrixElement": "JHUGen", 
    "Couplings": {},  
    "Prod": True,
    "Dec": False, 
    "ispm4l": False,
    "context": "Reco", 
    "computeprop": False, 
    "useconstant": True, # This argument affects the normalization, and should generally be set to true for reconstructed events, false forgenerator level events. [from https://spin.pha.jhu.edu/Manual.pdf]
    "addPmavjj": True,
    "addPmavjj_true": True
    }, append=True)

# https://github.com/CJLST/ZZAnalysis/blob/e1dd72d210d05e248873389c97da5ff9c9fcb94a/AnalysisStep/test/prod/pyFragments/RecoProbabilities.py#L165
setConf("probabilities", {
    "Name": "JJQCD_SIG_ghg2_1_JHUGen_JECNominal",
    "Process": "HSMHiggs", 
    "Production": "JJQCD", 
    "MatrixElement": "JHUGen", 
    "Couplings": {}, 
    "Prod": True, 
    "Dec": False, 
    "ispm4l": False,
    "context": "Reco", 
    "computeprop": False,
    "useconstant": True # This argument affects the normalization, and should generally be set to true for reconstructed events, false forgenerator level events. [from https://spin.pha.jhu.edu/Manual.pdf]
    }, append=True)
