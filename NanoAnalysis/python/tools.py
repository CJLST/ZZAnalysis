### Various helper funcions
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

## handle configuration variables
_myConf = {}

def getConf(name,default=None):
    global _myConf
    return _myConf[name] if name in _myConf else default

def setConf(name,value=True):
    global _myConf
    _myConf[name] = value
"root://cms-xrd-global.cern.ch/"


# Insert a module in a processing sequence before the specified module
def insertBefore(sequence, moduleName, module) :
    for im, m in enumerate(sequence) :
        if type(m).__name__ == moduleName :
            sequence.insert(im, module)
            return


# Insert a module in a processing sequence after the specified module
def insertAfter(sequence, moduleName, module) :
    for im, m in enumerate(sequence) :
        if type(m).__name__ == moduleName :
            sequence.insert(im+1, module)
            return

# Get the four leptons of a ZZ or ZLL candidate
def getLeptons(aCand, event) :
    idxs = [aCand.Z1l1Idx, aCand.Z1l2Idx, aCand.Z2l1Idx, aCand.Z2l2Idx]
    electrons = Collection(event, "Electron")
    muons = Collection(event, "Muon")
    leps = list(electrons) + list(muons)
    return [leps[i] for i in idxs]
