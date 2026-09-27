def getJetVetoMap(era, tag) :
    from PhysicsTools.NATModules.modules.jetVetoMap import jetVMAP

    if era == 2016 and "UL" in tag:
        if "ULAPV" in tag:
            folderKey = "Run2-2016preVFP-UL-NanoAODv15/2026-06-05"
            corrName = "Summer19UL16_V1"
        else:
            folderKey = "Run2-2016postVFP-UL-NanoAODv15/2026-06-05"
            corrName = "Summer19UL16_V1"

    elif era == 2017 and "UL" in tag:
        folderKey = "Run2-2017-UL-NanoAODv15/2026-06-05"
        corrName = "Summer19UL17_V1"

    elif era == 2018 and "UL" in tag:
        folderKey = "Run2-2018-UL-NanoAODv15/2026-06-05"
        corrName = "Summer19UL18_V1"

    elif era == 2022:
            if "pre_EE" in tag:
                folderKey = "Run3-22CDSep23-Summer22-NanoAODv12/2026-06-05"
                corrName = "Summer22_23Sep2023_RunCD_V1"
            else:
                folderKey = "Run3-22EFGSep23-Summer22EE-NanoAODv12/2026-06-05"
                corrName = "Summer22EE_23Sep2023_RunEFG_V1"
    
    elif era == 2023:
            if "pre_BPix" in tag:
                folderKey = "Run3-23CSep23-Summer23-NanoAODv12/2026-07-15"
                corrName = "Summer23Prompt23_RunC_V1"
            else:
                folderKey = "Run3-23DSep23-Summer23BPix-NanoAODv12/2026-07-15"
                corrName = "Summer23BPixPrompt23_RunD_V1"

    elif era == 2024:
        folderKey = "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2026-07-16"
        corrName = "Summer24Prompt24_RunBCDEFGHI_V1"

    elif era == 2025:
        folderKey = "Run3-25Prompt-Summer24-NanoAODv15/2026-07-16/"
        corrName = "Summer24Prompt25_RunCDEFG_V1"

    elif era == 2026:
        folderKey = "Run3-26Prompt-Summer24-NanoAODv15/2026-07-15/"
        corrName = "Summer24Prompt26_RunBCD_V1"

    else :
        raise ValueError("getJetvetoMap: Era", era, "not supported")
        
    json_JVMAP = f"/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/{folderKey}/jetvetomaps.json.gz"
    veto_map_name= "jetvetomap"
    print("***jetJVMAP: era:", era, "tag:", tag, "corrName:", corrName, "json:", json_JVMAP)

    return jetVMAP(json_JVMAP, corrName, veto_map_name)
