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
    if type(moduleName)!=str :
        raise ValueError("insertBefore: moduleName should be a string")
    for im, m in enumerate(sequence) :
        if type(m).__name__ == moduleName :
            sequence.insert(im, module)
            return


# Insert a module in a processing sequence after the specified module
def insertAfter(sequence, moduleName, module) :
    if type(moduleName)!=str :
        raise ValueError("insertAfter: moduleName should be a string")
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

def Mother(part, gen):
    '''
        Find the ID and Idx of the mother of a given GenPart (`part`)
        amongst all the particles in GenPart (`gen`) collection.
        The function returns Idx and ID of the mother.
    '''
    idxMother= part.genPartIdxMother
    while idxMother>=0 and gen[idxMother].pdgId == part.pdgId:
        idxMother = gen[idxMother].genPartIdxMother
    idMother=0
    if idxMother >=0 : idMother = gen[idxMother].pdgId
    return idxMother, idMother

def getParentID(part, gen) :
    '''
        Return the ID of the leptons's parent:
        25 for H->Z->l; 23 for Z->l; +-15 for
        tau->l if genlep is e,mu.
    '''
    pIdx, pID = Mother(part, gen)
    if pIdx < 0 : return 0
    ppIdx = gen[pIdx].genPartIdxMother
    if pID == 23 and ppIdx>=0 and gen[ppIdx].pdgId == 25 :
        pID = 25
    return pID

def lhe_logger(genpart, LHEPart=None):
    print ("---Gen:")
    for i, gp in enumerate(genpart) :
        motherId=-1
        gmotherId=-1
        if gp.genPartIdxMother >= 0 : 
            motherId = genpart[gp.genPartIdxMother].pdgId
            if genpart[gp.genPartIdxMother].genPartIdxMother >= 0 :
                gmotherId = genpart[genpart[gp.genPartIdxMother].genPartIdxMother].pdgId
        print (i, gp.pdgId, gp.genPartIdxMother, gp.pt, gp.eta, gp.phi, gp.p4().M(), gp.status)

    if LHEPart != None :
        print("---------LHEPart---------")
        for i, Lp in enumerate(LHEPart):
            print(i, Lp.pdgId, Lp.pt, Lp.eta, Lp.status, Lp.incomingpz)

def get_genEventSumw(input_file, maxEntriesPerSample=None):
    '''
       Util function to get the sum of weights per event.
       Returns the sum of weights, similarly to what we
       stored in Counters->GetBinContent(40) in the miniAODs.
    '''
    f = input_file

    runs  = f.Runs
    event = f.Events
    nRuns = runs.GetEntries()
    nEntries = event.GetEntries()

    iRun = 0
    genEventCount = 0
    genEventSumw = 0.

    while iRun < nRuns and runs.GetEntry(iRun) :
        genEventCount += runs.genEventCount
        genEventSumw += runs.genEventSumw
        iRun +=1
    print ("gen=", genEventCount, "sumw=", genEventSumw)

    if maxEntriesPerSample is not None:
        print(f"Scaling to {maxEntriesPerSample} entries")
        if nEntries>maxEntriesPerSample :
            genEventSumw = genEventSumw*maxEntriesPerSample/nEntries
            nEntries=maxEntriesPerSample
        print("    scaled to:", nEntries, "sumw=", genEventSumw)

    return genEventSumw

# Return efficiency and asymmetric (up, down) errors for sel over tot events
def getEff(tot, sel):
    from ROOT import TEfficiency
    eff = sel/tot
    up = TEfficiency.ClopperPearson(tot, sel, 0.683, True)
    dn = TEfficiency.ClopperPearson(tot, sel, 0.683, False)
    return eff, up-eff, eff-dn
