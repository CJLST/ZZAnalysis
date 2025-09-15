def getJetCorrected(era, tag, is_mc, overwritePt=True) :
    from PhysicsTools.NATModules.modules.jetCorr import jetJERC

    if era not in [2022, 2023, 2024]:
        raise ValueError("getJetCorrected: Era", era, "not supported")

    if era == 2022:
        if is_mc :
            if "pre_EE" in tag:
                folderKey = "2022_Summer22" # JERC file md5sum: 2fa54b84c2739e7e3c21c375d2ddcac7
                L1Key = "Summer22_22Sep2023_V2_MC_L1FastJet_AK4PFPuppi"
                L2Key = "Summer22_22Sep2023_V2_MC_L2Relative_AK4PFPuppi"
                L3Key = "Summer22_22Sep2023_V2_MC_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer22_22Sep2023_V2_MC_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = "Summer22_22Sep2023_V2_MC_Total_AK4PFPuppi"
                smearKey = "JERSmear"
                JERKey = "Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi"
                JERsfKey = "Summer22_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi"
            else:
                folderKey = "2022_Summer22EE" # JERC file md5sum: baa619665139acf05fc326aaaa0571cd
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
                folderKey = "2022_Summer22" # JERC file md5sum: 2fa54b84c2739e7e3c21c375d2ddcac7
                L1Key = "Summer22_22Sep2023_RunCD_V2_DATA_L1FastJet_AK4PFPuppi"
                L2Key = "Summer22_22Sep2023_RunCD_V2_DATA_L2Relative_AK4PFPuppi"
                L3Key = "Summer22_22Sep2023_RunCD_V2_DATA_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer22_22Sep2023_RunCD_V2_DATA_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = None
                smearKey = None
                JERKey = None
                JERsfKey = None
            elif "2022E" in tag:
                folderKey = "2022_Summer22EE" # JERC file md5sum: baa619665139acf05fc326aaaa0571cd
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
                folderKey = "2023_Summer23" # JERC file md5sum: b35c2108478e49da5d46f679c28f1111
                L1Key = "Summer23Prompt23_V2_MC_L1FastJet_AK4PFPuppi"
                L2Key = "Summer23Prompt23_V2_MC_L2Relative_AK4PFPuppi"
                L3Key = "Summer23Prompt23_V2_MC_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer23Prompt23_V2_MC_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = "Summer23Prompt23_V2_MC_Total_AK4PFPuppi"
                smearKey = "JERSmear"
                JERKey = "Summer23Prompt23_RunCv1234_JRV1_MC_PtResolution_AK4PFPuppi"
                JERsfKey = "Summer23Prompt23_RunCv1234_JRV1_MC_ScaleFactor_AK4PFPuppi"
            else:
                folderKey = "2023_Summer23BPix" # JERC file md5sum: eef8019335e15d48a7963c47f7ad306c
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
                folderKey = "2023_Summer23" # JERC file md5sum: b35c2108478e49da5d46f679c28f1111
                L1Key = "Summer23Prompt23_V2_DATA_L1FastJet_AK4PFPuppi"
                L2Key = "Summer23Prompt23_V2_DATA_L2Relative_AK4PFPuppi"
                L3Key = "Summer23Prompt23_V2_DATA_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer23Prompt23_V2_DATA_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = None
                smearKey = None
                JERKey = None
                JERsfKey = None
            else:
                folderKey = "2023_Summer23BPix" # JERC file md5sum: eef8019335e15d48a7963c47f7ad306c
                L1Key = "Summer23BPixPrompt23_V3_DATA_L1FastJet_AK4PFPuppi"
                L2Key = "Summer23BPixPrompt23_V3_DATA_L2Relative_AK4PFPuppi"
                L3Key = "Summer23BPixPrompt23_V3_DATA_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer23BPixPrompt23_V3_DATA_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = None
                smearKey = None
                JERKey = None
                JERsfKey = None

    elif era == 2024:
        if is_mc :
            raise ValueError("getJetCorrected: 2024 MC not yet supported")
        else :
            # folderKey = "2024_Winter24" # JERC file md5sum: a0c4f7f29e09162f56c07a9b5fb97d1e
            # L1Key = "Winter24Prompt24_V3_DATA_L1FastJet_AK4PFPuppi"
            # L2Key = "Winter24Prompt24_V3_DATA_L2Relative_AK4PFPuppi"
            # L3Key = "Winter24Prompt24_V3_DATA_L3Absolute_AK4PFPuppi"
            # L2L3Key = "Winter24Prompt24_V3_DATA_L2L3Residual_AK4PFPuppi"
            folderKey = "2024_Summer24" # JERC file md5sum: 754fd45c85b197ff9f7d33f68e7cd9a2
            L1Key = "Summer24Prompt24_V1_DATA_L1FastJet_AK4PFPuppi"
            L2Key = "Summer24Prompt24_V1_DATA_L2Relative_AK4PFPuppi"
            L3Key = "Summer24Prompt24_V1_DATA_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer24Prompt24_V1_DATA_L2L3Residual_AK4PFPuppi"
            scaleTotalKey = None
            smearKey = None
            JERKey = None
            JERsfKey = None

                
    json_JERC = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/%s/jet_jerc.json.gz" % (folderKey)
    json_JERsmear = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/jer_smear.json.gz" # md5sum: 390e4be4be109bb1a2d3a116f2c9386a

    # Determine usePhiDependentJEC based on the tag
    usePhiDependentJEC = era >= 2023 and not ("pre_BPix" in tag) # False up to 2023 pre_BPix, True in 2023 post_BPix and afterwards
    # Apply run-dependent JEC only for 2023 data (not MC)
    useRunDependentJEC = (era == 2023 or era == 2024 or era == 2025) and (not is_mc)

    print("***jetJERC: era:", era, "tag:", tag, "is MC:", is_mc, "overwritePt:", overwritePt, "phiDependent:", usePhiDependentJEC, "runDependent:", useRunDependentJEC, "json_JERC:", json_JERC, "json_JERsmear:", json_JERsmear)
    
    return jetJERC(json_JERC, json_JERsmear, L1Key, L2Key, L3Key, L2L3Key, scaleTotalKey, smearKey, JERKey, JERsfKey, overwritePt, usePhiDependentJEC, useRunDependentJEC)
