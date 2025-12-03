"""
Instantiate jetId correctionlib module (for nanoAODv13 onwards, cf. https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/jetidExample.py)
"""

def getJetIdProducer(era, tag) :
    from PhysicsTools.NATModules.modules.jetId import jetId
    
    if era == 2022:
        if "pre_EE" in tag:
            json = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2022_Summer22/jetid.json.gz" #md5sum: 2070556451837fe611d6e0b0218a5d1f
        else:
            json = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2022_Summer22EE/jetid.json.gz" # Note: file is a link to the 2022_Summer22 file
    
    elif era == 2023:
        if "pre_BPix" in tag:
            json = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2023_Summer23/jetid.json.gz" # Note: file is a link to the 2022_Summer22 file
        else:
            json = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2023_Summer23BPix/jetid.json.gz" # Note: file is a link to the 2022_Summer22 file

    elif era == 2024:
       json = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2025-07-17/jetid.json.gz" #from new versioned repository; identical to 2022_Summer22

    elif era >= 2016 and era <=2018:
        # FIXME: Assume the same as 2022_Summer22 since json file is not present in official repositories
        print("WARNING: official jetid json not available, using the one for 2022_Summer22")
        json = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2022_Summer22/jetid.json.gz"

    else: 
        raise ValueError("getJetIdProducer: get: Era", era, tag, "not supported")  

    print("***jetId: era:", era, "tag:", tag, "json:", json)
       
    return jetId(json)
