"""
Instantiate jetId correctionlib module (for nanoAODv13 onwards, cf. https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/jetidExample.py)
"""

def getJetIdProducer(era, tag) :
    from PhysicsTools.NATModules.modules.jetId import jetId
    if era not in [2022,2023,2024]:
        raise ValueError("getJetIdProducer: get: Era", era, "not supported")

    if era == 2022:
        if "pre_EE" in tag:
            folderKey = "2022_Summer22" #md5sum: 2070556451837fe611d6e0b0218a5d1f
        else:
            folderKey = "2022_Summer22EE" # Note: file is a link to the 2022_Summer22 file
    
    elif era == 2023:
        if "pre_BPix" in tag:
            folderKey = "2023_Summer23" # Note: file is a link to the 2022_Summer22 file
        else:
            folderKey = "2023_Summer23BPix" # Note: file is a link to the 2022_Summer22 file

    elif era == 2024:
       folderKey = "2024_Summer24" # Note: file is a link to the 2022_Summer22 file

    json = f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/{folderKey}/jetid.json.gz"
    print("***jetId: era:", era, "tag:", tag, "json:", json)
       
    return jetId(json)
