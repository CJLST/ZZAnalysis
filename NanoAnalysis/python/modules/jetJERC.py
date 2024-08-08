import os

def getJetCorrected(era, tag, is_mc, overwritePt=True) :
    from PhysicsTools.NATModules.modules.jetCorr import jetJERC

    if era != 2022:
        raise ValueError("getJetCorrected: Era", era, "not supported")

    if is_mc :
        if tag == "pre_EE":
            folderKey = "2022"
            L1Key = "Summer22_22Sep2023_V2_MC_L1FastJet_AK4PFPuppi"
            L2Key = "Summer22_22Sep2023_V2_MC_L2Relative_AK4PFPuppi"
            L3Key = "Summer22_22Sep2023_V2_MC_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer22_22Sep2023_V2_MC_L2L3Residual_AK4PFPuppi"
            smearKey = "JERSmear"
            JERKey = "Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi"
            JERsfKey = "Summer22_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi"
        elif tag == "post_EE":
            folderKey = "2022EE"
            L1Key = "Summer22EE_22Sep2023_V2_MC_L1FastJet_AK4PFPuppi"
            L2Key = "Summer22EE_22Sep2023_V2_MC_L2Relative_AK4PFPuppi"
            L3Key = "Summer22EE_22Sep2023_V2_MC_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer22EE_22Sep2023_V2_MC_L2L3Residual_AK4PFPuppi"
            smearKey = "JERSmear"
            JERKey = "Summer22EE_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi"
            JERsfKey = "Summer22EE_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi"
    ## Data
    ## JER are not applied to data
    else :
        if tag == "pre_EE":
            folderKey = "2022"
            L1Key = "Summer22_22Sep2023_RunCD_V2_DATA_L1FastJet_AK4PFPuppi"
            L2Key = "Summer22_22Sep2023_RunCD_V2_DATA_L2Relative_AK4PFPuppi"
            L3Key = "Summer22_22Sep2023_RunCD_V2_DATA_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer22_22Sep2023_RunCD_V2_DATA_L2L3Residual_AK4PFPuppi"
            smearKey = None
            JERKey = None
            JERsfKey = None
        elif tag == "E":
            folderKey = "2022EE"
            L1Key = "Summer22EE_22Sep2023_RunE_V2_DATA_L1FastJet_AK4PFPuppi"
            L2Key = "Summer22EE_22Sep2023_RunE_V2_DATA_L2Relative_AK4PFPuppi"
            L3Key = "Summer22EE_22Sep2023_RunE_V2_DATA_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer22EE_22Sep2023_RunE_V2_DATA_L2L3Residual_AK4PFPuppi"
            smearKey = None
            JERKey = None
            JERsfKey = None
        elif tag == "F":
            folderKey = "2022EE"
            L1Key = "Summer22EE_22Sep2023_RunF_V2_DATA_L1FastJet_AK4PFPuppi"
            L2Key = "Summer22EE_22Sep2023_RunF_V2_DATA_L2Relative_AK4PFPuppi"
            L3Key = "Summer22EE_22Sep2023_RunF_V2_DATA_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer22EE_22Sep2023_RunF_V2_DATA_L2L3Residual_AK4PFPuppi"
            smearKey = None
            JERKey = None
            JERsfKey = None
        elif tag == "G":
            folderKey = "2022EE"
            L1Key = "Summer22EE_22Sep2023_RunG_V2_DATA_L1FastJet_AK4PFPuppi"
            L2Key = "Summer22EE_22Sep2023_RunG_V2_DATA_L2Relative_AK4PFPuppi"
            L3Key = "Summer22EE_22Sep2023_RunG_V2_DATA_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer22EE_22Sep2023_RunG_V2_DATA_L2L3Residual_AK4PFPuppi"
            smearKey = None
            JERKey = None
            JERsfKey = None


    json_JERC = "%s/src/ZZAnalysis/NanoAnalysis/data/JERC/%s/jet_jerc.json" % (os.environ['CMSSW_BASE'], folderKey)
    json_JERsmear = "%s/src/ZZAnalysis/NanoAnalysis/data/JERC/jer_smear.json" % (os.environ['CMSSW_BASE'])

    print("***jetJERC: era:", era, "tag:", tag, "is MC:", is_mc, "overwritePt:", overwritePt, "json_JERC:", json_JERC, "json_JERsmear:", json_JERsmear)
    return jetJERC(json_JERC, json_JERsmear, L1Key, L2Key, L3Key, L2L3Key, smearKey, JERKey, JERsfKey, overwritePt)
