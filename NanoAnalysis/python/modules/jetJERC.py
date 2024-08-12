import os

def getJetCorrected(era, tag, is_mc, overwritePt=True) :
    from PhysicsTools.NATModules.modules.jetCorr import jetJERC

    if era != 2022:
        raise ValueError("getJetCorrected: Era", era, "not supported")

    if is_mc :
        if "pre_EE" in tag:
            folderKey = "2022_Summer22"
            L1Key = "Summer22_22Sep2023_V2_MC_L1FastJet_AK4PFPuppi"
            L2Key = "Summer22_22Sep2023_V2_MC_L2Relative_AK4PFPuppi"
            L3Key = "Summer22_22Sep2023_V2_MC_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer22_22Sep2023_V2_MC_L2L3Residual_AK4PFPuppi"
            smearKey = "JERSmear"
            JERKey = "Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi"
            JERsfKey = "Summer22_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi"
        else:
            folderKey = "2022_Summer22EE"
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
        if "pre_EE" in tag:
            folderKey = "2022_Summer22"
            L1Key = "Summer22_22Sep2023_RunCD_V2_DATA_L1FastJet_AK4PFPuppi"
            L2Key = "Summer22_22Sep2023_RunCD_V2_DATA_L2Relative_AK4PFPuppi"
            L3Key = "Summer22_22Sep2023_RunCD_V2_DATA_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer22_22Sep2023_RunCD_V2_DATA_L2L3Residual_AK4PFPuppi"
            smearKey = None
            JERKey = None
            JERsfKey = None
        elif "2022E" in tag:
            folderKey = "2022_Summer22EE"
            L1Key = "Summer22EE_22Sep2023_RunE_V2_DATA_L1FastJet_AK4PFPuppi"
            L2Key = "Summer22EE_22Sep2023_RunE_V2_DATA_L2Relative_AK4PFPuppi"
            L3Key = "Summer22EE_22Sep2023_RunE_V2_DATA_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer22EE_22Sep2023_RunE_V2_DATA_L2L3Residual_AK4PFPuppi"
            smearKey = None
            JERKey = None
            JERsfKey = None
        elif "2022F" in tag:
            folderKey = "2022_Summer22EE"
            L1Key = "Summer22EE_22Sep2023_RunF_V2_DATA_L1FastJet_AK4PFPuppi"
            L2Key = "Summer22EE_22Sep2023_RunF_V2_DATA_L2Relative_AK4PFPuppi"
            L3Key = "Summer22EE_22Sep2023_RunF_V2_DATA_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer22EE_22Sep2023_RunF_V2_DATA_L2L3Residual_AK4PFPuppi"
            smearKey = None
            JERKey = None
            JERsfKey = None
        elif "2022G" in tag:
            folderKey = "2022_Summer22EE"
            L1Key = "Summer22EE_22Sep2023_RunG_V2_DATA_L1FastJet_AK4PFPuppi"
            L2Key = "Summer22EE_22Sep2023_RunG_V2_DATA_L2Relative_AK4PFPuppi"
            L3Key = "Summer22EE_22Sep2023_RunG_V2_DATA_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer22EE_22Sep2023_RunG_V2_DATA_L2L3Residual_AK4PFPuppi"
            smearKey = None
            JERKey = None
            JERsfKey = None
        else:
            raise ValueError("getJetCorrected: tag", era, "not supported")


    json_JERC = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/%s/jet_jerc.json.gz" % (folderKey)
    json_JERsmear = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/jer_smear.json.gz"

    print("***jetJERC: era:", era, "tag:", tag, "is MC:", is_mc, "overwritePt:", overwritePt, "json_JERC:", json_JERC, "json_JERsmear:", json_JERsmear)
    return jetJERC(json_JERC, json_JERsmear, L1Key, L2Key, L3Key, L2L3Key, smearKey, JERKey, JERsfKey, overwritePt)
