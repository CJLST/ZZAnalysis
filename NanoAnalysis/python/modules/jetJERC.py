def getJetCorrected(era, tag, is_mc, overwritePt=True) :
    from PhysicsTools.NATModules.modules.jetCorr import jetJERC

    if era not in [2022, 2023, 2024]:
        raise ValueError("getJetCorrected: Era", era, "not supported")

# Regrouped uncertainties (11 sources) - The {year} part indicates those uncertainties that need to be kept uncorrelated between these datasets.
# Taken from https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/merge_requests/120#9968ad259fd43c9ba2351d217c42dec468fe273f
    jes_systematics_11split = [
        "Regrouped_Absolute",
        "Regrouped_Absolute_{year}",
        "Regrouped_BBEC1",
        "Regrouped_BBEC1_{year}",
        "Regrouped_EC2",
        "Regrouped_EC2_{year}",
        "Regrouped_FlavorQCD",
        "Regrouped_HF",
        "Regrouped_HF_{year}",
        "Regrouped_RelativeBal",
        "Regrouped_RelativeSample_{year}",
    ]

    if era == 2022:
        if is_mc :
            if "pre_EE" in tag:
                folderKey = "Run3-22CDSep23-Summer22-NanoAODv12/2025-09-23"
                L1Key = "Summer22_22Sep2023_V3_MC_L1FastJet_AK4PFPuppi"
                L2Key = "Summer22_22Sep2023_V3_MC_L2Relative_AK4PFPuppi"
                L3Key = "Summer22_22Sep2023_V3_MC_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer22_22Sep2023_V3_MC_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = "Summer22_22Sep2023_V3_MC_Total_AK4PFPuppi"
                scaleKeyRegrouped11 = [
                f"Summer22_22Sep2023_V3_MC_{label.format(year=era)}_AK4PFPuppi" for label in jes_systematics_11split
                ]
                smearKey = "JERSmear"
                JERKey = "Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi"
                JERsfKey = "Summer22_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi"
            else:
                folderKey = "Run3-22EFGSep23-Summer22EE-NanoAODv12/2025-10-07"
                L1Key = "Summer22EE_22Sep2023_V3_MC_L1FastJet_AK4PFPuppi"
                L2Key = "Summer22EE_22Sep2023_V3_MC_L2Relative_AK4PFPuppi"
                L3Key = "Summer22EE_22Sep2023_V3_MC_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer22EE_22Sep2023_V3_MC_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = "Summer22EE_22Sep2023_V3_MC_Total_AK4PFPuppi"
                scaleKeyRegrouped11 = [
                f"Summer22EE_22Sep2023_V3_MC_{label.format(year='2022EE')}_AK4PFPuppi" for label in jes_systematics_11split
                ]
                smearKey = "JERSmear"
                JERKey = "Summer22EE_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi"
                JERsfKey = "Summer22EE_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi"
        ## Data
        ## JER are not applied to data
        else :
            if "pre_EE" in tag:
                folderKey = "Run3-22CDSep23-Summer22-NanoAODv12/2025-09-23"
                L1Key = "Summer22_22Sep2023_RunCD_V3_DATA_L1FastJet_AK4PFPuppi"
                L2Key = "Summer22_22Sep2023_RunCD_V3_DATA_L2Relative_AK4PFPuppi"
                L3Key = "Summer22_22Sep2023_RunCD_V3_DATA_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer22_22Sep2023_RunCD_V3_DATA_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = None
                scaleKeyRegrouped11 = None 
                smearKey = None
                JERKey = None
                JERsfKey = None
            elif "2022E" in tag:
                folderKey = "Run3-22EFGSep23-Summer22EE-NanoAODv12/2025-10-07"
                L1Key = "Summer22EE_22Sep2023_RunE_V3_DATA_L1FastJet_AK4PFPuppi"
                L2Key = "Summer22EE_22Sep2023_RunE_V3_DATA_L2Relative_AK4PFPuppi"
                L3Key = "Summer22EE_22Sep2023_RunE_V3_DATA_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer22EE_22Sep2023_RunE_V3_DATA_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = None
                scaleKeyRegrouped11 = None 
                smearKey = None
                JERKey = None
                JERsfKey = None
            elif "2022F" in tag:
                folderKey = "Run3-22EFGSep23-Summer22EE-NanoAODv12/2025-10-07"
                L1Key = "Summer22EE_22Sep2023_RunF_V3_DATA_L1FastJet_AK4PFPuppi"
                L2Key = "Summer22EE_22Sep2023_RunF_V3_DATA_L2Relative_AK4PFPuppi"
                L3Key = "Summer22EE_22Sep2023_RunF_V3_DATA_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer22EE_22Sep2023_RunF_V3_DATA_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = None
                scaleKeyRegrouped11 = None 
                smearKey = None
                JERKey = None
                JERsfKey = None
            elif "2022G" in tag:
                folderKey = "Run3-22EFGSep23-Summer22EE-NanoAODv12/2025-10-07"
                L1Key = "Summer22EE_22Sep2023_RunG_V3_DATA_L1FastJet_AK4PFPuppi"
                L2Key = "Summer22EE_22Sep2023_RunG_V3_DATA_L2Relative_AK4PFPuppi"
                L3Key = "Summer22EE_22Sep2023_RunG_V3_DATA_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer22EE_22Sep2023_RunG_V3_DATA_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = None
                scaleKeyRegrouped11 = None 
                smearKey = None
                JERKey = None
                JERsfKey = None
            else:
                raise ValueError("getJetCorrected: tag", era, "not supported")

    elif era == 2023:
        if is_mc :
            if "pre_BPix" in tag:
                folderKey = "Run3-23CSep23-Summer23-NanoAODv12/2025-10-07"
                L1Key = "Summer23Prompt23_V2_MC_L1FastJet_AK4PFPuppi"
                L2Key = "Summer23Prompt23_V2_MC_L2Relative_AK4PFPuppi"
                L3Key = "Summer23Prompt23_V2_MC_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer23Prompt23_V2_MC_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = "Summer23Prompt23_V2_MC_Total_AK4PFPuppi"
                scaleKeyRegrouped11 = [
                f"Summer23Prompt23_V2_MC_{label.format(year='2023')}_AK4PFPuppi" for label in jes_systematics_11split
                ]
                smearKey = "JERSmear"
                JERKey = "Summer23Prompt23_RunCv1234_JRV1_MC_PtResolution_AK4PFPuppi"
                JERsfKey = "Summer23Prompt23_RunCv1234_JRV1_MC_ScaleFactor_AK4PFPuppi"
            else:
                folderKey = "Run3-23DSep23-Summer23BPix-NanoAODv12/2025-10-07"
                L1Key = "Summer23BPixPrompt23_V3_MC_L1FastJet_AK4PFPuppi"
                L2Key = "Summer23BPixPrompt23_V3_MC_L2Relative_AK4PFPuppi"
                L3Key = "Summer23BPixPrompt23_V3_MC_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer23BPixPrompt23_V3_MC_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = "Summer23BPixPrompt23_V3_MC_Total_AK4PFPuppi"
                scaleKeyRegrouped11 = [
                f"Summer23BPixPrompt23_V3_MC_{label.format(year='2023BPix')}_AK4PFPuppi" for label in jes_systematics_11split
                ]
                smearKey = "JERSmear"
                JERKey = "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi"
                JERsfKey = "Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi"
        ## Data
        ## JER are not applied to data
        else :
            if "pre_BPix" in tag:
                folderKey = "Run3-23CSep23-Summer23-NanoAODv12/2025-10-07"
                L1Key = "Summer23Prompt23_V2_DATA_L1FastJet_AK4PFPuppi"
                L2Key = "Summer23Prompt23_V2_DATA_L2Relative_AK4PFPuppi"
                L3Key = "Summer23Prompt23_V2_DATA_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer23Prompt23_V2_DATA_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = None
                scaleKeyRegrouped11 = None 
                smearKey = None
                JERKey = None
                JERsfKey = None
            else:
                folderKey = "Run3-23DSep23-Summer23BPix-NanoAODv12/2025-10-07"
                L1Key = "Summer23BPixPrompt23_V3_DATA_L1FastJet_AK4PFPuppi"
                L2Key = "Summer23BPixPrompt23_V3_DATA_L2Relative_AK4PFPuppi"
                L3Key = "Summer23BPixPrompt23_V3_DATA_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer23BPixPrompt23_V3_DATA_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = None
                scaleKeyRegrouped11 = None 
                smearKey = None
                JERKey = None
                JERsfKey = None

    elif era == 2024:
        if is_mc :
            folderKey = "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2025-07-17"
            L1Key = "Summer24Prompt24_V1_MC_L1FastJet_AK4PFPuppi"
            L2Key = "Summer24Prompt24_V1_MC_L2Relative_AK4PFPuppi"
            L3Key = "Summer24Prompt24_V1_MC_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer24Prompt24_V1_MC_L2L3Residual_AK4PFPuppi"
            scaleTotalKey = "Summer24Prompt24_V1_MC_Total_AK4PFPuppi"
            scaleKeyRegrouped11 = [
                f"Summer24Prompt24_V1_MC_{label.format(year='2024')}_AK4PFPuppi" for label in jes_systematics_11split
                ]
            smearKey = "JERSmear"
            # It appears the most recent 23Bpix files are used in the following cases: 
            JERKey = "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi"
            JERsfKey = "Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi"
 
        else :
            # folderKey = "2024_Winter24" # JERC file md5sum: a0c4f7f29e09162f56c07a9b5fb97d1e
            # L1Key = "Winter24Prompt24_V3_DATA_L1FastJet_AK4PFPuppi"
            # L2Key = "Winter24Prompt24_V3_DATA_L2Relative_AK4PFPuppi"
            # L3Key = "Winter24Prompt24_V3_DATA_L3Absolute_AK4PFPuppi"
            # L2L3Key = "Winter24Prompt24_V3_DATA_L2L3Residual_AK4PFPuppi"
            folderKey = "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2025-07-17"
            L1Key = "Summer24Prompt24_V1_DATA_L1FastJet_AK4PFPuppi"
            L2Key = "Summer24Prompt24_V1_DATA_L2Relative_AK4PFPuppi"
            L3Key = "Summer24Prompt24_V1_DATA_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer24Prompt24_V1_DATA_L2L3Residual_AK4PFPuppi"
            scaleTotalKey = None
            scaleKeyRegrouped11 = None 
            smearKey = None
            JERKey = None
            JERsfKey = None

    json_JERC = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/%s/jet_jerc.json.gz" % (folderKey)
    json_JERsmear = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/jer_smear.json.gz" # md5sum: 390e4be4be109bb1a2d3a116f2c9386a

    # Determine usePhiDependentJEC based on the tag
    usePhiDependentJEC = era >= 2023 and not ("pre_BPix" in tag) # False up to 2023 pre_BPix, True in 2023 post_BPix and afterwards
    # Apply run-dependent JEC only for 2023 data (not MC)
    useRunDependentJEC = (era == 2023 or era == 2024 or era == 2025) and (not is_mc)
    # Use Splittign scheme for Jets uncertainties (11 sources)
    useJesSplittingScheme11 = False # Default

    scaleKey = scaleKeyRegrouped11 if useJesSplittingScheme11 else scaleTotalKey

    print("***jetJERC: era:", era, "tag:", tag, "is MC:", is_mc, "overwritePt:", overwritePt, "phiDependent:", usePhiDependentJEC, "runDependent:", useRunDependentJEC, "JesSplittingScheme11:", useJesSplittingScheme11,"json_JERC:", json_JERC, "json_JERsmear:", json_JERsmear)
    
    return jetJERC(json_JERC, json_JERsmear, L1Key, L2Key, L3Key, L2L3Key, scaleKey, smearKey, JERKey, JERsfKey, overwritePt, usePhiDependentJEC, useRunDependentJEC)
