def getJetCorrected(era, tag, is_mc, useJesSplittingScheme11, overwritePt=True) :
    from PhysicsTools.NATModules.modules.jetCorr import jetJERC


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


# FIXME: old nanoAODv9 corrections, option currently unsupported!!!
#    if era == 2016 and "UL" in tag:
#         folderKey = "Run2-2018-UL-NanoAODv9/2025-04-11"
#         if is_mc :
#             L1Key = "Summer19UL18_V5_MC_L1FastJet_AK4PFchs"
#             L2Key = "Summer19UL18_V5_MC_L2Relative_AK4PFchs"
#             L3Key = "Summer19UL18_V5_MC_L3Absolute_AK4PFchs"
#             L2L3Key = "Summer19UL18_V5_MC_L2L3Residual_AK4PFchs"
#             scaleTotalKey = "Summer19UL18_V5_MC_Total_AK4PFchs"
#             scaleKeyRegrouped11 = [
#                 f"Summer19UL18_V5_MC_{label.format(year='2018')}_AK4PFchs" for label in jes_systematics_11split
#                 ]
#             smearKey = "JERSmear"
#             # It appears the most recent 23Bpix files are used in the following cases: 
#             JERKey = "Summer19UL18_JRV2_MC_PtResolution_AK4PFchs"
#             JERsfKey = "Summer19UL18_JRV2_MC_ScaleFactor_AK4PFchs"
#             JERsfUncKey = None
 
#         else :
#             L1Key = "Summer19UL18_RunA_V5_DATA_L1FastJet_AK4PFchs"
#             L2Key = "Summer19UL18_RunA_V5_DATA_L2Relative_AK4PFchs"
#             L3Key = "Summer19UL18_RunA_V5_DATA_L3Absolute_AK4PFchs"
#             L2L3Key = "Summer19UL18_RunA_V5_DATA_L2L3Residual_AK4PFchs"
#             scaleTotalKey = None
#             scaleKeyRegrouped11 = None 
#             smearKey = None
#             JERKey = None
#             JERsfKey = None
#             JERsfUncKey = None


    
    ### Run 2 corrections for nanoAODv15 samples
    if era == 2016 and "UL" in tag:
        if "ULAPV" in tag:
            pass
        else:
            pass
        raise ValueError("jetJERC: 2016 to be implemented")
    
    elif era == 2017 and "UL" in tag:
        raise ValueError("jetJERC: 2017 to be implemented")
    
    elif era == 2018 and  "UL" in tag:
        folderKey = "Run2-2018-UL-NanoAODv15/2026-06-05"
        if is_mc :
            L1Key = "Summer20UL18NanoV15_V1_MC_L1FastJet_AK4PFPuppi"
            L2Key = "Summer20UL18NanoV15_V1_MC_L2Relative_AK4PFPuppi"
            L3Key = "Summer20UL18NanoV15_V1_MC_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer20UL18NanoV15_V1_MC_L2L3Residual_AK4PFPuppi"
            scaleTotalKey = "Summer20UL18NanoV15_V1_MC_Total_AK4PFPuppi"
            scaleKeyRegrouped11 = [
                f"Summer20UL18NanoV15_V1_MC_{label.format(year='2018')}_AK4PFPuppi" for label in jes_systematics_11split
                ]
            smearKey = "JERSmear"
            # It appears the most recent 23Bpix files are used in the following cases: 
            JERKey = "Summer19UL18_JRV2_MC_PtResolution_AK4PFPuppi"
            JERsfKey = "Summer19UL18_JRV2_MC_ScaleFactor_AK4PFPuppi"
            JERsfUncKey = None
 
        else :
            L1Key = "Summer20UL18NanoV15_V1_DATA_L1FastJet_AK4PFPuppi"
            L2Key = "Summer20UL18NanoV15_V1_DATA_L2Relative_AK4PFPuppi"
            L3Key = "Summer20UL18NanoV15_V1_DATA_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer20UL18NanoV15_V1_DATA_L2L3Residual_AK4PFPuppi"
            scaleTotalKey = None
            scaleKeyRegrouped11 = None 
            smearKey = None
            JERKey = None
            JERsfKey = None
            JERsfUncKey = None
        
    elif era == 2022:
        if is_mc :
            if "pre_EE" in tag:
                folderKey = "Run3-22CDSep23-Summer22-NanoAODv12/2026-06-05"
                L1Key = "Summer22_22Sep2023_V4_MC_L1FastJet_AK4PFPuppi"
                L2Key = "Summer22_22Sep2023_V4_MC_L2Relative_AK4PFPuppi"
                L3Key = "Summer22_22Sep2023_V4_MC_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer22_22Sep2023_V4_MC_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = "Summer22_22Sep2023_V4_MC_Total_AK4PFPuppi"
                scaleKeyRegrouped11 = [
                f"Summer22_22Sep2023_V4_MC_{label.format(year=era)}_AK4PFPuppi" for label in jes_systematics_11split
                ]
                smearKey = "JERSmear"
                JERKey = "Summer22_22Sep2023_JRV2_MC_PtResolution_AK4PFPuppi"
                JERsfKey = "Summer22_22Sep2023_JRV2_MC_ScaleFactor_AK4PFPuppi"
                JERsfUncKey = "Summer22_22Sep2023_JRV2_MC_SFUncertainty_AK4PFPuppi"
            else:
                folderKey = "Run3-22EFGSep23-Summer22EE-NanoAODv12/2026-06-05"
                L1Key = "Summer22EE_22Sep2023_V4_MC_L1FastJet_AK4PFPuppi"
                L2Key = "Summer22EE_22Sep2023_V4_MC_L2Relative_AK4PFPuppi"
                L3Key = "Summer22EE_22Sep2023_V4_MC_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer22EE_22Sep2023_V4_MC_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = "Summer22EE_22Sep2023_V4_MC_Total_AK4PFPuppi"
                scaleKeyRegrouped11 = [
                f"Summer22EE_22Sep2023_V4_MC_{label.format(year='2022EE')}_AK4PFPuppi" for label in jes_systematics_11split
                ]
                smearKey = "JERSmear"
                JERKey = "Summer22EE_22Sep2023_JRV2_MC_PtResolution_AK4PFPuppi"
                JERsfKey = "Summer22EE_22Sep2023_JRV2_MC_ScaleFactor_AK4PFPuppi"
                JERsfUncKey = "Summer22EE_22Sep2023_JRV2_MC_SFUncertainty_AK4PFPuppi"
        ## Data
        ## JER are not applied to data
        else :
            if "pre_EE" in tag:
                folderKey = "Run3-22CDSep23-Summer22-NanoAODv12/2026-06-05"
                L1Key = "Summer22_22Sep2023_V4_DATA_L1FastJet_AK4PFPuppi"
                L2Key = "Summer22_22Sep2023_V4_DATA_L2Relative_AK4PFPuppi"
                L3Key = "Summer22_22Sep2023_V4_DATA_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer22_22Sep2023_V4_DATA_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = None
                scaleKeyRegrouped11 = None 
                smearKey = None
                JERKey = None
                JERsfKey = None
                JERsfUncKey = None
            else:
                folderKey = "Run3-22EFGSep23-Summer22EE-NanoAODv12/2026-06-05"
                L1Key = "Summer22EE_22Sep2023_V4_DATA_L1FastJet_AK4PFPuppi"
                L2Key = "Summer22EE_22Sep2023_V4_DATA_L2Relative_AK4PFPuppi"
                L3Key = "Summer22EE_22Sep2023_V4_DATA_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer22EE_22Sep2023_V4_DATA_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = None
                scaleKeyRegrouped11 = None 
                smearKey = None
                JERKey = None
                JERsfKey = None
                JERsfUncKey = None

    elif era == 2023:
        if is_mc :
            if "pre_BPix" in tag:
                folderKey = "Run3-23CSep23-Summer23-NanoAODv12/2026-07-15"
                L1Key = "Summer23Prompt23_V4_MC_L1FastJet_AK4PFPuppi"
                L2Key = "Summer23Prompt23_V4_MC_L2Relative_AK4PFPuppi"
                L3Key = "Summer23Prompt23_V4_MC_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer23Prompt23_V4_MC_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = "Summer23Prompt23_V4_MC_Total_AK4PFPuppi"
                scaleKeyRegrouped11 = [
                f"Summer23Prompt23_V4_MC_{label.format(year='2023')}_AK4PFPuppi" for label in jes_systematics_11split
                ]
                smearKey = "JERSmear"
                JERKey = "Summer23Prompt23_RunCv1234_JRV3_MC_PtResolution_AK4PFPuppi"
                JERsfKey = "Summer23Prompt23_RunCv1234_JRV3_MC_ScaleFactor_AK4PFPuppi"
                JERsfUncKey = "Summer23Prompt23_RunCv1234_JRV3_MC_SFUncertainty_AK4PFPuppi"
            else:
                folderKey = "Run3-23DSep23-Summer23BPix-NanoAODv12/2026-07-15"
                L1Key = "Summer23BPixPrompt23_V4_MC_L1FastJet_AK4PFPuppi"
                L2Key = "Summer23BPixPrompt23_V4_MC_L2Relative_AK4PFPuppi"
                L3Key = "Summer23BPixPrompt23_V4_MC_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer23BPixPrompt23_V4_MC_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = "Summer23BPixPrompt23_V4_MC_Total_AK4PFPuppi"
                scaleKeyRegrouped11 = [
                f"Summer23BPixPrompt23_V4_MC_{label.format(year='2023BPix')}_AK4PFPuppi" for label in jes_systematics_11split
                ]
                smearKey = "JERSmear"
                JERKey = "Summer23BPixPrompt23_RunD_JRV3_MC_PtResolution_AK4PFPuppi"
                JERsfKey = "Summer23BPixPrompt23_RunD_JRV3_MC_ScaleFactor_AK4PFPuppi"
                JERsfUncKey = "Summer23BPixPrompt23_RunD_JRV3_MC_SFUncertainty_AK4PFPuppi"
        ## Data
        ## JER are not applied to data
        else :
            if "pre_BPix" in tag:
                folderKey = "Run3-23CSep23-Summer23-NanoAODv12/2026-07-15"
                L1Key = "Summer23Prompt23_V4_DATA_L1FastJet_AK4PFPuppi"
                L2Key = "Summer23Prompt23_V4_DATA_L2Relative_AK4PFPuppi"
                L3Key = "Summer23Prompt23_V4_DATA_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer23Prompt23_V4_DATA_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = None
                scaleKeyRegrouped11 = None 
                smearKey = None
                JERKey = None
                JERsfKey = None
                JERsfUncKey = None
            else:
                folderKey = "Run3-23DSep23-Summer23BPix-NanoAODv12/2026-07-15"
                L1Key = "Summer23BPixPrompt23_V4_DATA_L1FastJet_AK4PFPuppi"
                L2Key = "Summer23BPixPrompt23_V4_DATA_L2Relative_AK4PFPuppi"
                L3Key = "Summer23BPixPrompt23_V4_DATA_L3Absolute_AK4PFPuppi"
                L2L3Key = "Summer23BPixPrompt23_V4_DATA_L2L3Residual_AK4PFPuppi"
                scaleTotalKey = None
                scaleKeyRegrouped11 = None 
                smearKey = None
                JERKey = None
                JERsfKey = None
                JERsfUncKey = None

    elif era == 2024:
        folderKey = "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2026-07-16"
        if is_mc :
            L1Key = "Summer24Prompt24_V5_MC_L1FastJet_AK4PFPuppi"
            L2Key = "Summer24Prompt24_V5_MC_L2Relative_AK4PFPuppi"
            L3Key = "Summer24Prompt24_V5_MC_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer24Prompt24_V5_MC_L2L3Residual_AK4PFPuppi"
            scaleTotalKey = "Summer24Prompt24_V5_MC_Total_AK4PFPuppi"
            scaleKeyRegrouped11 = [
                f"Summer24Prompt24_V5_MC_{label.format(year='2024')}_AK4PFPuppi" for label in jes_systematics_11split
                ]
            smearKey = "JERSmear"
            JERKey = "Summer24Prompt24_JRV2_MC_PtResolution_AK4PFPuppi"
            JERsfKey = "Summer24Prompt24_JRV2_MC_ScaleFactor_AK4PFPuppi"
            JERsfUncKey = "Summer24Prompt24_JRV2_MC_SFUncertainty_AK4PFPuppi"
        else :
            L1Key = "Summer24Prompt24_V5_DATA_L1FastJet_AK4PFPuppi"
            L2Key = "Summer24Prompt24_V5_DATA_L2Relative_AK4PFPuppi"
            L3Key = "Summer24Prompt24_V5_DATA_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer24Prompt24_V5_DATA_L2L3Residual_AK4PFPuppi"
            scaleTotalKey = None
            scaleKeyRegrouped11 = None
            smearKey = None
            JERKey = None
            JERsfKey = None
            JERsfUncKey = None

    elif era == 2025:
        folderKey = "Run3-25Prompt-Summer24-NanoAODv15/2026-07-16"
        if is_mc :
            L1Key = "Summer24Prompt25_V3_MC_L1FastJet_AK4PFPuppi"
            L2Key = "Summer24Prompt25_V3_MC_L2Relative_AK4PFPuppi"
            L3Key = "Summer24Prompt25_V3_MC_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer24Prompt25_V3_MC_L2L3Residual_AK4PFPuppi"
            scaleTotalKey = "Summer24Prompt25_V3_MC_Total_AK4PFPuppi"
            scaleKeyRegrouped11 = [
                f"Summer24Prompt25_V3_MC_{label.format(year='2025')}_AK4PFPuppi" for label in jes_systematics_11split
                ]
            smearKey = "JERSmear"
            JERKey = "Summer24Prompt25_JRV2_MC_PtResolution_AK4PFPuppi"
            JERsfKey = "Summer24Prompt25_JRV2_MC_ScaleFactor_AK4PFPuppi"
            JERsfUncKey = "Summer24Prompt25_JRV2_MC_SFUncertainty_AK4PFPuppi"
        else :
            L1Key = "Summer24Prompt25_V3_DATA_L1FastJet_AK4PFPuppi"
            L2Key = "Summer24Prompt25_V3_DATA_L2Relative_AK4PFPuppi"
            L3Key = "Summer24Prompt25_V3_DATA_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer24Prompt25_V3_DATA_L2L3Residual_AK4PFPuppi"
            scaleTotalKey = None
            scaleKeyRegrouped11 = None
            smearKey = None
            JERKey = None
            JERsfKey = None
            JERsfUncKey = None

    elif era == 2026:
        folderKey = "Run3-26Prompt-Summer24-NanoAODv15/2026-07-15"
        if is_mc :
            L1Key = "Summer24Prompt26_V1_MC_L1FastJet_AK4PFPuppi"
            L2Key = "Summer24Prompt26_V1_MC_L2Relative_AK4PFPuppi"
            L3Key = "Summer24Prompt26_V1_MC_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer24Prompt26_V1_MC_L2L3Residual_AK4PFPuppi"
            scaleTotalKey = "Summer24Prompt26_V1_MC_Total_AK4PFPuppi"
            scaleKeyRegrouped11 = [
                f"Summer24Prompt26_V1_MC_{label.format(year='2026')}_AK4PFPuppi" for label in jes_systematics_11split
                ]
            smearKey = "JERSmear"
            JERKey = "Summer24Prompt26_JRV1_MC_PtResolution_AK4PFPuppi"
            JERsfKey = "Summer24Prompt26_JRV1_MC_ScaleFactor_AK4PFPuppi"
            JERsfUncKey = "Summer24Prompt26_JRV1_MC_SFUncertainty_AK4PFPuppi"
        else :
            L1Key = "Summer24Prompt26_V1_DATA_L1FastJet_AK4PFPuppi"
            L2Key = "Summer24Prompt26_V1_DATA_L2Relative_AK4PFPuppi"
            L3Key = "Summer24Prompt26_V1_DATA_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer24Prompt26_V1_DATA_L2L3Residual_AK4PFPuppi"
            scaleTotalKey = None
            scaleKeyRegrouped11 = None
            smearKey = None
            JERKey = None
            JERsfKey = None
            JERsfUncKey = None
            
    else:
        raise ValueError("getJetCorrected: Era", era, tag, "not supported")


    json_JERC = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/%s/jet_jerc.json.gz" % (folderKey)
    json_JERsmear = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/JER-Smearing/2025-11-03/jer_smear.json.gz"

    # Determine usePhiDependentJEC based on the tag
    usePhiDependentJEC = era >= 2023 and not ("pre_BPix" in tag) # False up to 2023 pre_BPix, True in 2023 post_BPix and afterwards
    # Use run-dependent L2L3Relative JEC for data (currently the case for all eras)
    useRunDependentJEC = (not is_mc)
    
    # Use Splittigng scheme for Jets uncertainties (11 sources)
    scaleKey = scaleKeyRegrouped11 if useJesSplittingScheme11 else scaleTotalKey

    print("***jetJERC: era:", era, "tag:", tag, "is MC:", is_mc, "overwritePt:", overwritePt, "phiDependent:", usePhiDependentJEC, "runDependent:", useRunDependentJEC, "JesSplittingScheme11:", useJesSplittingScheme11,"json_JERC:", json_JERC, "json_JERsmear:", json_JERsmear)
    
    return jetJERC(era, json_JERC, json_JERsmear, L1Key, L2Key, L3Key, L2L3Key, scaleKey, smearKey, JERKey, JERsfKey, JERsfUncKey, overwritePt, usePhiDependentJEC, useRunDependentJEC)
