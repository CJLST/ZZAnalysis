"""
Any enums defined in interface/miscenums.h will also be accessible here.
"""

def _():
    """put all the dirty stuff inside a function so it doesn't get imported with import *"""
    import os, re

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


    result = []
    for enumname, enumcontents in enums:
        enumcontents = enumcontents.split(",")
        for content in enumcontents:
            result.append(content.strip())
    return result

for enum in _():
    exec(enum)
del enum, _
