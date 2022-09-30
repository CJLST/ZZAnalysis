### Various helper funcions

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

