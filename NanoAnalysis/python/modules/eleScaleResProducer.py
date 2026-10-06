import os

# Set up the NATModules eleScaleRes module
def getEleScaleRes(era, tag, is_mc, overwritePt=True):
    from PhysicsTools.NATModules.modules.eleScaleRes import eleScaleRes

    scaleKey = "Scale"
    smearKey = "SmearAndSyst" if is_mc else None

    # Note: Scale and Smearing were already applied on Run2UL NanoAODv9 samples, this is no longer the case in NanoAODv15: this module should not be called for Run2 UL v9 samples!
    if era == 2016 and "UL" in tag:
        if "ULAPV" in tag:
            fname = "/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run2-2016preVFP-UL-NanoAODv15/2025-12-05/electronSS_EtDependent.json.gz"
        else:
            fname = "/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run2-2016postVFP-UL-NanoAODv15/2025-12-05/electronSS_EtDependent.json.gz"

    elif era == 2017 and "UL" in tag:
        fname = "/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run2-2017-UL-NanoAODv15/2025-12-05/electronSS_EtDependent.json.gz"

    elif era == 2018 and "UL" in tag:
        fname = "/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run2-2018-UL-NanoAODv15/2025-12-05/electronSS_EtDependent.json.gz"
        
    elif era == 2022:
        if "pre_EE" in tag :
            fname = "/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-22CDSep23-Summer22-NanoAODv12/2025-12-15/electronSS_EtDependent.json.gz"
        else:
            fname = "/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-22EFGSep23-Summer22EE-NanoAODv12/2025-12-15/electronSS_EtDependent.json.gz"

    elif era == 2023:
        if "pre_BPix" in tag:
            fname = "/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-23CSep23-Summer23-NanoAODv12/2025-12-15/electronSS_EtDependent.json.gz"
        else:
            fname = "/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-23DSep23-Summer23BPix-NanoAODv12/2025-12-15/electronSS_EtDependent.json.gz"

    elif era == 2024:
        fname = "/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2025-12-15/electronSS_EtDependent.json.gz"

    elif era == 2025:
        fname = "/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-25Prompt-Summer24-NanoAODv15/2026-06-26/electronSS_EtDependent.json.gz"

    elif era == 2026:
        fname = "/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-25Prompt-Summer24-NanoAODv15/2026-06-08/electronSS_EtDependent.json.gz"

    else :
        raise ValueError(f"getEleScaleRes: Era {era}, tag {tag} not supported")

    print("***eleScaleRes: era:", era, "tag:", tag, "is MC:", is_mc, "overwritePt:", overwritePt, "json:", fname)
    return eleScaleRes(fname, scaleKey, smearKey, overwritePt)
