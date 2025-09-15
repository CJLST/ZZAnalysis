import os

# Set up the NATModules muonScaleRes module
def getMuonScaleRes(era, tag, is_mc, overwritePt=True) :
    from PhysicsTools.NATModules.modules.muonScaleRes import muonScaleRes 

    if era not in [2022, 2023, 2024]:
        raise ValueError(f"getMuonScaleRes: Era {era} is not supported")

    if era == 2022:
        if "pre_EE" in tag :
            fname = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2022_Summer22/muon_scalesmearing.json.gz" # md5sum: 415165703d2ca3724f1cd0f97bdf31fe

        else :
            fname = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2022_Summer22EE/muon_scalesmearing.json.gz" # md5sum: 8cb780c7e1b4507263a9edecf9b38fc2

    elif era == 2023:
        if "pre_BPix" in tag:
            fname = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2023_Summer23/muon_scalesmearing.json.gz" # md5sum: e7612461ea9416447ae9fb4a038d82cd

        else:
            fname = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2023_Summer23BPix/muon_scalesmearing.json.gz" # md5sum: 50684fe80408116aa9f59c308433b8d6

    elif era == 2024:
        print(f"WARNING {era} muonScaleRes - for now using 2023BPix")
        fname = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2023_Summer23BPix/muon_scalesmearing.json.gz"

    print("***muonScaleRes: era:", era, "tag:", tag, "is MC:", is_mc, "overwritePt:", overwritePt, "json:", fname)
    return muonScaleRes(fname, is_mc, overwritePt, minPt=3.)

