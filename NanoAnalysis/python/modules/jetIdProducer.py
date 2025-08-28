"""
Instantiate jetId correctionlib module (for nanoAODv13 onwards, cf. https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/jetidExample.py)
"""

def getJetIdProducer(era, tag) :
    from PhysicsTools.NATModules.modules.jetId import jetId
    if era not in [2022,2023,2024]:
        raise ValueError("getJetIdProducer: get: Era", era, "not supported")

    if era == 2022:
        if "pre_EE" in tag:
            folderKey = "2022_Summer22"
        else:
            folderKey = "2022_Summer22EE"
    
    elif era == 2023:
        if "pre_BPix" in tag:
            folderKey = "2023_Summer23"
        else:
            folderKey = "2023_Summer23BPix"

    elif era == 2024:
       folderKey = "2024_Summer24"

    json = f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/{folderKey}/jetid.json.gz"
    print("***jetId: era:", era, "tag:", tag, "json:", json)
       
    return jetId(json)
