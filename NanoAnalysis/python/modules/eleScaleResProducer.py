import os

# Set up the NATModules eleScaleRes module
def getEleScaleRes(era, tag, is_mc, overwritePt=True) :
    from PhysicsTools.NATModules.modules.eleScaleRes import eleScaleRes

    if era != 2022:
        raise ValueError("getEleScaleRes: Era", era, "not supported")

    scaleKey = "Scale" 
    if is_mc :
        smearKey = "Smearing"
    else :
        smearKey=None

    if "pre_EE" in tag :
        fname = "electronSS_preEE.json.gz"
    else:
        fname = "electronSS_postEE.json.gz"
 
    json = "%s/src/ZZAnalysis/NanoAnalysis/data/%s" % (os.environ['CMSSW_BASE'], fname)

    print("***eleScaleRes: era:", era, "tag:", tag, "is MC:", is_mc, "overwritePt:", overwritePt, "json:", json)
    return eleScaleRes(json, scaleKey, smearKey, overwritePt)
