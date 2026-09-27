"""
Instantiate jetId update/computation module according to JME recipes.
"""

def getJetIdProducer(era, tag, nanoVersion) :
    from PhysicsTools.NATModules.modules.jetId import jetId
    
    if nanoVersion <=12 :
        json = None # correctionlib not used in this case
    elif nanoVersion > 12 :
        # Note: the 2022 and 2023 jsons are intended for v13 onwards; v12 does not contain the required variables and uses a special 
        # implmentation, cf: https://cms-talk.web.cern.ch/t/should-jetid-be-recomputed-with-jetid-json-gz-also-for-nanoaodv12/145862/2
        if era == 2022: 
            if "pre_EE" in tag:
                json = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22CDSep23-Summer22-NanoAODv12/2026-06-05/jetid.json.gz"
            else:
                json = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/2026-06-05/jetid.json.gz"

        elif era == 2023:
            if "pre_BPix" in tag:
                json = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23CSep23-Summer23-NanoAODv12/2026-07-15/jetid.json.gz"
            else:
                json = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23DSep23-Summer23BPix-NanoAODv12/2026-07-15/jetid.json.gz"

        elif era == 2024:
           json = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2026-07-16/jetid.json.gz"

        elif era == 2025:
           json = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-25Prompt-Summer24-NanoAODv15/2026-07-16/jetid.json.gz"

        elif era == 2026:
           json = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-26Prompt-Summer24-NanoAODv15/2026-07-15/jetid.json.gz"
           
        elif era >= 2016 and era <=2018:
            json = None
        else: 
            raise ValueError("getJetIdProducer: get: Era:", era, "tag:", tag, "nanoVersion:", nanoVersion, "not supported")  

    print("***jetId: era:", era, "tag:", tag, "nanoVersion:", nanoVersion, "json:", json)
       
    return jetId(json, nanoVersion=nanoVersion, era=era)
