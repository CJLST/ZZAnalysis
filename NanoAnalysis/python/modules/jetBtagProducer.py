"""
Instantiate the correctionlib-based b-tagging module.
"""

import os

def getJetBtagProducer(era, tag, is_mc, is_signal, WP="M"):
    from PhysicsTools.NATModules.modules.jetBtag import jetBtag

    data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "data", "btagEff"))

    if era == 2022:
        tagger = "particleNet"
        tagger_name = "btagPNetB"
        if "pre_EE" in tag:
            json_SF = "/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-22CDSep23-Summer22-NanoAODv12/2025-08-20/btagging.json.gz"
            json_eff = os.path.join(data_dir, f"btag_{'signal' if is_signal else 'background'}_2022.json.gz")
        else:
            json_SF = "/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-22EFGSep23-Summer22EE-NanoAODv12/2025-08-20/btagging.json.gz"
            json_eff = os.path.join(data_dir, f"btag_{'signal' if is_signal else 'background'}_2022EE.json.gz")

    elif era == 2023:
        tagger = "particleNet"
        tagger_name = "btagPNetB"
        if "pre_BPix" in tag:
            json_SF = "/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-23CSep23-Summer23-NanoAODv12/2025-08-20/btagging.json.gz"
            json_eff = os.path.join(data_dir, f"btag_{'signal' if is_signal else 'background'}_2023preBPix.json.gz")
        else:
            json_SF = "/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-23DSep23-Summer23BPix-NanoAODv12/2025-08-20/btagging.json.gz"
            json_eff = os.path.join(data_dir, f"btag_{'signal' if is_signal else 'background'}_2023postBPix.json.gz")

    elif era == 2024:
        tagger = "UParTAK4"
        tagger_name = "btagUParTAK4B"
        json_SF = "/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2026-03-10/btagging.json.gz"
        json_eff = os.path.join(data_dir, f"btag_{'signal' if is_signal else 'background'}_2024.json.gz")

    elif era == 2025:
        tagger = "UParTAK4"
        tagger_name = "btagUParTAK4B"
        json_SF = "/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-25Prompt-Summer24-NanoAODv15/2026-09-21/btagging.json.gz"
        print("JetBTag: WARNING: efficiency file not yet available, using the one for 2024")
        json_eff = os.path.join(data_dir, f"btag_{'signal' if is_signal else 'background'}_2024.json.gz")
        
    elif era >= 2016 and era <= 2018:
        # FIXME: official Run 2 SF JSONs and local efficiencies should be wired when available.
        print("JetBTag: WARNING: official btagging json not available, using the one for 2022_Summer22")
        tagger = "particleNet"
        tagger_name = "btagPNetB"
        json_SF = "/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-22CDSep23-Summer22-NanoAODv12/2025-08-20/btagging.json.gz"
        json_eff = os.path.join(data_dir, f"btag_{'signal' if is_signal else 'background'}_2022.json.gz")

    else:
        raise ValueError("getJetBtagProducer: Era", era, tag, "not supported")

    if is_mc and not os.path.exists(json_eff):
        raise ValueError("getJetBtagProducer: efficiency JSON not found: %s" % json_eff)

    print(
        "***jetBtag: era:",
        era,
        "tag:",
        tag,
        "is MC:",
        is_mc,
        "tagger:",
        tagger,
        "tagger branch:",
        tagger_name,
        "WP:",
        WP,
        "json_SF:",
        json_SF,
        "json_eff:",
        json_eff,
    )

    return jetBtag(is_mc, tagger, tagger_name, WP, json_SF, json_eff, ["correlated","uncorrelated"])
