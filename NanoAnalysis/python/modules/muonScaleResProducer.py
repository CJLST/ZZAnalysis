import os

# Set up the NATModules muonScaleRes module
def getMuonScaleRes(era, tag, is_mc, overwritePt=True) :
#FIXME: module to be moved in NATModules
#    from PhysicsTools.NATModules.modules.muonScaleRes import muonScaleRes 
    from ZZAnalysis.NanoAnalysis.modules.muonScaleRes import muonScaleRes    

    if era != 2022: #FIXME add 2023
        raise ValueError("getMuonScaleRes: Era", era, "not supported")

    if "pre_EE" in tag :
        fname = "2022_schemaV2.json.gz"
    else :
        fname = "2022EE_schemaV2.json.gz"
 
    json = "%s/src/ZZAnalysis/NanoAnalysis/data/MuonScale/%s" % (os.environ['CMSSW_BASE'], fname)

    print("***muonScaleRes: era:", era, "tag:", tag, "is MC:", is_mc, "overwritePt:", overwritePt, "json:", json)
    return muonScaleRes(json, is_mc, overwritePt)

