def getJetVetoMap(era, tag) :
    from PhysicsTools.NATModules.modules.jetVetoMap import jetVMAP

    if era not in [2022,2023,2024]:
        raise ValueError("getJetvetoMap: Era", era, "not supported")

    if era == 2022:
            if "pre_EE" in tag:
                folderKey = "2022_Summer22" # md5sum: e35eddf2b2eb072be63c78035cba01b0
                corrName = "Summer22_23Sep2023_RunCD_V1"
            else:
                folderKey = "2022_Summer22EE" # md5sum: ea16a7d736eacd58677dd761172792aa
                corrName = "Summer22EE_23Sep2023_RunEFG_V1"
    
    elif era == 2023:
            if "pre_BPix" in tag:
                folderKey = "2023_Summer23" # md5sum: 6e5ea1c9e7303dc72485303b1a6565b2
                corrName = "Summer23Prompt23_RunC_V1"
            else:
                folderKey = "2023_Summer23BPix" # md5sum: 1f296a697ca06593ad36a7d53b56fa77
                corrName = "Summer23BPixPrompt23_RunD_V1"

    elif era == 2024:
        # folderKey = "2024_Winter24" # md5sum: 802a8be50fde0fd45d86c98ea446b16b
        # corrName = "Winter24Prompt2024BCDEFGHI_V1"
        folderKey = "2024_Summer24" # md5sum: 1cf19b46f969dbd31cb07e664dcca3cf
        corrName = "Summer24Prompt2024_RunBCDEFGHI_V1"

    json_JVMAP = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/%s/jetvetomaps.json.gz" % (folderKey)
    veto_map_name= "jetvetomap"
    print("***jetJVMAP: era:", era, "tag:", tag, "corrName:", corrName, "json:", json_JVMAP)

    return jetVMAP(json_JVMAP, corrName, veto_map_name)
