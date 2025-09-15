import os

# Set up the NATModules eleScaleRes module
def getEleScaleRes(era, tag, is_mc, overwritePt=True, EtDependent=None):
    from PhysicsTools.NATModules.modules.eleScaleRes import eleScaleRes

    # Set default behavior: Standard for 2022, EtDependent for 2023
    if EtDependent is None:
        EtDependent = (era in [2022, 2023, 2024])

    # Check for supported eras
    if era not in [2022, 2023, 2024]:
        raise ValueError(f"getEleScaleRes: Era {era} not supported")

    # localpath = "%s/src/ZZAnalysis/NanoAnalysis/data/ElectronScale/" % (os.environ['CMSSW_BASE'])
    
    if era == 2022:
            if EtDependent:
                if "pre_EE" in tag :
                    scaleKey = "EGMScale_Compound_Ele_2022preEE"
                    smearKey = "EGMSmearAndSyst_ElePTsplit_2022preEE" if is_mc else None
                    fname = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2022_Summer22/electronSS_EtDependent.json.gz" # md5sum: ccbe63a9c79802df5ff1c217d799a435
                else:
                    scaleKey = "EGMScale_Compound_Ele_2022postEE"
                    smearKey = "EGMSmearAndSyst_ElePTsplit_2022postEE" if is_mc else None
                    fname = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2022_Summer22EE/electronSS_EtDependent.json.gz" # md5sum: 24f3e4b18fcf2589050edc6b35284437  

            else: # Older "standard" version, used for early results, now superseeded by ET-dependent version
                if "pre_EE" in tag :
                    scaleKey = "2022Re-recoBCD_ScaleJSON"
                    smearKey = "2022Re-recoBCD_SmearingJSON" if is_mc else None
                    fname = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2022_Summer22/electronSS.json.gz" # md5sum: db4b5696ecf33088f1d849d95643a40f  

                else:
                    scaleKey = "2022Re-recoE+PromptFG_ScaleJSON"
                    smearKey = "2022Re-recoE+PromptFG_SmearingJSON" if is_mc else None
                    fname = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2022_Summer22EE/electronSS.json.gz" # md5sum: 471982e6e8a4bfeecc1f974058646d48  

    elif era == 2023:
        if "pre_BPix" in tag:
            scaleKey = "EGMScale_Compound_Ele_2023preBPIX"
            smearKey = "EGMSmearAndSyst_ElePTsplit_2023preBPIX" if is_mc else None
            fname = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2023_Summer23/electronSS_EtDependent.json.gz" # md5sum: e16d16f9ab0abe13cca28432417f3a48

        else:
            scaleKey = "EGMScale_Compound_Ele_2023postBPIX"
            smearKey = "EGMSmearAndSyst_ElePTsplit_2023postBPIX" if is_mc else None
            fname = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2023_Summer23BPix/electronSS_EtDependent.json.gz" #md5sum: 19be4932b71f55eb9353eb01d9fbeaca

    elif era == 2024:
        print(f"WARNING {era} electron SS - for now using 2023BPix")
        scaleKey = "EGMScale_Compound_Ele_2023postBPIX"
        smearKey = "EGMSmearAndSyst_ElePTsplit_2023postBPIX" if is_mc else None
        fname = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2023_Summer23BPix/electronSS_EtDependent.json.gz" #md5sum: 19be4932b71f55eb9353eb01d9fbeaca


    print("***eleScaleRes: era:", era, "tag:", tag, "is MC:", is_mc, "overwritePt:", overwritePt, "EtDependent:", EtDependent, "json:", fname)
    return eleScaleRes(fname, scaleKey, smearKey, overwritePt, EtDependent)
