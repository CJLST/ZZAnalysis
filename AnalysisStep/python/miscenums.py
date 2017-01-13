"""
Any enums defined in interface/miscenums.h will also be accessible here.
"""

def _():
    """put all the dirty stuff inside a function so it doesn't get imported with import *"""
    import os, re
    import ROOT

    def include(filename):
      ROOT.gROOT.ProcessLine("#include <{}>".format(filename))

    include("ZZAnalysis/AnalysisStep/interface/miscenums.h")

    with open(os.path.join(os.environ["CMSSW_BASE"], "src/ZZAnalysis/AnalysisStep/interface/miscenums.h")) as f:
        contents = f.read()

    lines = contents.split("\n")
    #remove comments and #ifdef
    lines = [line.split("//")[0] for line in lines]
    lines = [line.split("#")[0] for line in lines]
    contents = " ".join(lines)

    enums = re.findall(r'enum *(\w*) *'
                       r'[{]('
                       r'[^}]*'
                       r')[}]'
                       , contents)


    for enumname, enumcontents in enums:
        globals()[enumname] = getattr(ROOT, enumname)
        enumcontents = enumcontents.split(",")
        enumitems = [content.split("=")[0].strip() for content in enumcontents]
        for itemname in enumitems:
            globals()[itemname] = getattr(ROOT, itemname)

_()
del _
