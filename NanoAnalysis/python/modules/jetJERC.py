import os

def getJetCorrected(era, tag, is_mc, overwritePt=True) :
    from PhysicsTools.NATModules.modules.jetCorr import jetJERC

    if era not in [2022,2023]:
        raise ValueError("getJetCorrected: Era", era, "not supported")

    if era == 2022:
        if is_mc :
            if "pre_EE" in tag:
                folderKey = "2022_Summer22"
                L1Key = "Summer22_22Sep2023_V2_MC_L1FastJet_AK4PFPuppi"
                L2Key = "Summer22_22Sep2023_V2_MC_L2Relative_AK4PFPuppi"
                L3Key = "Summer22_22Sep2023_V2_MC_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer22_22Sep2023_V2_MC_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = "Summer22_22Sep2023_V2_MC_Total_AK4PFPuppi"
                smearKey = "JERSmear"
                JERKey = "Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi"
                JERsfKey = "Summer22_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi"
            else:
                folderKey = "2022_Summer22EE"
                L1Key = "Summer22EE_22Sep2023_V2_MC_L1FastJet_AK4PFPuppi"
                L2Key = "Summer22EE_22Sep2023_V2_MC_L2Relative_AK4PFPuppi"
                L3Key = "Summer22EE_22Sep2023_V2_MC_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer22EE_22Sep2023_V2_MC_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = "Summer22EE_22Sep2023_V2_MC_Total_AK4PFPuppi"
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
                scaleTotalKey = None
                smearKey = None
                JERKey = None
                JERsfKey = None
            elif "2022E" in tag:
                folderKey = "2022_Summer22EE"
                L1Key = "Summer22EE_22Sep2023_RunE_V2_DATA_L1FastJet_AK4PFPuppi"
                L2Key = "Summer22EE_22Sep2023_RunE_V2_DATA_L2Relative_AK4PFPuppi"
                L3Key = "Summer22EE_22Sep2023_RunE_V2_DATA_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer22EE_22Sep2023_RunE_V2_DATA_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = None
                smearKey = None
                JERKey = None
                JERsfKey = None
            elif "2022F" in tag:
                folderKey = "2022_Summer22EE"
                L1Key = "Summer22EE_22Sep2023_RunF_V2_DATA_L1FastJet_AK4PFPuppi"
                L2Key = "Summer22EE_22Sep2023_RunF_V2_DATA_L2Relative_AK4PFPuppi"
                L3Key = "Summer22EE_22Sep2023_RunF_V2_DATA_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer22EE_22Sep2023_RunF_V2_DATA_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = None
                smearKey = None
                JERKey = None
                JERsfKey = None
            elif "2022G" in tag:
                folderKey = "2022_Summer22EE"
                L1Key = "Summer22EE_22Sep2023_RunG_V2_DATA_L1FastJet_AK4PFPuppi"
                L2Key = "Summer22EE_22Sep2023_RunG_V2_DATA_L2Relative_AK4PFPuppi"
                L3Key = "Summer22EE_22Sep2023_RunG_V2_DATA_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer22EE_22Sep2023_RunG_V2_DATA_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = None
                smearKey = None
                JERKey = None
                JERsfKey = None
            else:
                raise ValueError("getJetCorrected: tag", era, "not supported")
    elif era == 2023:
        if is_mc :
            if "pre_BPix" in tag:
                folderKey = "2023_Summer23"
                L1Key = "Summer23Prompt23_V2_MC_L1FastJet_AK4PFPuppi"
                L2Key = "Summer23Prompt23_V2_MC_L2Relative_AK4PFPuppi"
                L3Key = "Summer23Prompt23_V2_MC_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer23Prompt23_V2_MC_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = "Summer23Prompt23_V2_MC_Total_AK4PFPuppi"
                smearKey = "JERSmear"
                JERKey = "Summer23Prompt23_RunCv1234_JRV1_MC_PtResolution_AK4PFPuppi"
                JERsfKey = "Summer23Prompt23_RunCv1234_JRV1_MC_ScaleFactor_AK4PFPuppi"
            else:
                folderKey = "2023_Summer23BPix"
                L1Key = "Summer23BPixPrompt23_V3_MC_L1FastJet_AK4PFPuppi"
                L2Key = "Summer23BPixPrompt23_V3_MC_L2Relative_AK4PFPuppi"
                L3Key = "Summer23BPixPrompt23_V3_MC_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer23BPixPrompt23_V3_MC_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = "Summer23BPixPrompt23_V3_MC_Total_AK4PFPuppi"
                smearKey = "JERSmear"
                JERKey = "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi"
                JERsfKey = "Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi"
        ## Data
        ## JER are not applied to data
        else :
            if "pre_BPix" in tag:
                folderKey = "2023_Summer23"
                L1Key = "Summer23Prompt23_V2_DATA_L1FastJet_AK4PFPuppi"
                L2Key = "Summer23Prompt23_V2_DATA_L2Relative_AK4PFPuppi"
                L3Key = "Summer23Prompt23_V2_DATA_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer23Prompt23_V2_DATA_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = None
                smearKey = None
                JERKey = None
                JERsfKey = None
            else:
                folderKey = "2023_Summer23BPix"
                L1Key = "Summer23BPixPrompt23_V3_DATA_L1FastJet_AK4PFPuppi"
                L2Key = "Summer23BPixPrompt23_V3_DATA_L2Relative_AK4PFPuppi"
                L3Key = "Summer23BPixPrompt23_V3_DATA_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer23BPixPrompt23_V3_DATA_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = None
                smearKey = None
                JERKey = None
                JERsfKey = None

    json_JERC = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/%s/jet_jerc.json.gz" % (folderKey)
    json_JERsmear = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/jer_smear.json.gz"

    print("***jetJERC: era:", era, "tag:", tag, "is MC:", is_mc, "overwritePt:", overwritePt, "json_JERC:", json_JERC, "json_JERsmear:", json_JERsmear)
    # Determine usePhiDependentJEC based on the tag
    usePhiDependentJEC = "post_BPix" in tag # True if "post_BPix" is in tag, False otherwise

    return jetJERC(json_JERC, json_JERsmear, L1Key, L2Key, L3Key, L2L3Key, scaleTotalKey, smearKey, JERKey, JERsfKey, overwritePt, usePhiDependentJEC)
