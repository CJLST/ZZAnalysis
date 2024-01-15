### Analyze MC Truth
# -Optionally dump MC history
# -Add ZZ gen valriables:
# -GenZZFinalState: product of IDs of the four gen ZZ leptons
# -GenZZ_*Idx: ZZ lepton indices in the GenPart collection (FIXME: unsorted)
# -FsrPhoton_genFsrIdx: index of the closest gen FSR from Z->l (e, mu)
# 
###

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.HeppyCore.utils.deltar import deltaR
from ROOT import TLorentzVector

from ZZAnalysis.NanoAnalysis.tools import Mother, getParentID

class mcTruthAnalyzer(Module):
    def __init__(self, dump=False):
        print("***mcTruthAnalyzer", flush=True)
        self.writeHistFile = False
        self.printGenHist  = dump # print MC history


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("GenZZ_FinalState", "I")
        self.out.branch("GenZZ_mass", "F")
        self.out.branch("GenZZ_Z1l1Idx", "I") # Indices in the GenPart
        self.out.branch("GenZZ_Z1l2Idx", "I")
        self.out.branch("GenZZ_Z2l1Idx", "I")
        self.out.branch("GenZZ_Z2l2Idx", "I")
        self.out.branch("FsrPhoton_genFsrIdx", "I", lenVar="nFsrPhoton")

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        genpart=Collection(event,"GenPart")

        ### Print Gen history and LHE particles.
        ### See also: https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/common/hepmcDump.py
        if self.printGenHist :
            lhe_logger(genpart)

        ## Search for gen FSR from Z->ll (e, mu)
        genFSRIdxs = []
        for i, gp in enumerate(genpart) :
            midx = gp.genPartIdxMother
            if gp.pdgId==22 and gp.pt > 2. and midx >= 0 :
                mid = genpart[midx].pdgId
                if abs(mid) == 11 or abs(mid) == 13 :
                    mmidx, mmid = Mother(genpart[midx], genpart) # possibly skip intermediate rows
                    if mmid == 23 :
                        genFSRIdxs.append(i)

        # Store FsrPhoton_genFsrIdx
        fsrPhotons = Collection(event, "FsrPhoton")
        fsrPhotons_genFsrIdx = [-1]*len(fsrPhotons)
        for ifsr, fsr in enumerate(fsrPhotons) :
            dRmin = 0.3
            for igen in genFSRIdxs :
                genfsr = genpart[igen]
                dR = deltaR (fsr.eta, fsr.phi, genfsr.eta, genfsr.phi)
                if dR < dRmin :
                    dRmin = dR
                    fsrPhotons_genFsrIdx[ifsr] = igen

        self.out.fillBranch("FsrPhoton_genFsrIdx", fsrPhotons_genFsrIdx)


#        GenHLeps = filter(lambda f : 
#                           (abs(f.pdgId)==11 or abs(f.pdgId)==13 or abs(f.pdgId)==15) and
#                           f.genPartIdxMother >=0 and genpart[f.genPartIdxMother].pdgId == 23 and
#                           genpart[f.genPartIdxMother].genPartIdxMother >= 0 and
#                           genpart[ genpart[f.genPartIdxMother].genPartIdxMother].pdgId == 25, genpart) #leptons generated from H->ZZ
#        GenEl_acc = filter(lambda f: abs(f.pdgId)== 11 and abs(f.eta)<2.5 and f.pt > conf["elePt"], GenHLeps)
#        GenMu_acc = filter(lambda f: abs(f.pdgId)== 13 and abs(f.eta)<2.4 and f.pt > conf["muPt"], GenHLeps)


        # Select leptons from ZZ and associated production. 
        theGenZZLeps = []
        theAssociatedLeps = []
        theGenH = None
        for ip, p in enumerate(genpart) :
            if (p.pdgId)==25 : theGenH = p
            if (abs(p.pdgId)==11 or abs(p.pdgId)==13 or abs(p.pdgId)==15) and p.genPartIdxMother >=0 :
                mid = genpart[p.genPartIdxMother].pdgId
                pid = getParentID(p, genpart)
                if mid == 25 or (mid == 23 and pid == 25) : # Lepton from H->(Z->)ll; note that this is the first daughter in the H or Z line; ie pre-FSR
                    theGenZZLeps.append(ip)
                elif ((mid == 23 and pid == 23) or mid == 24): # Associated production, or ZZ (sorted out below).
                    # FIXME this may miss leptons from ZZTo4lamcatnlo, cf. MCHistoryTools.cc
                    theAssociatedLeps.append(ip)

        # handle ZZ samples. Note that tribosons samples are not be handled here.
        if len(theAssociatedLeps)==4 and len(theGenZZLeps)==0 :
            theGenZZLeps, theAssociatedLeps = theAssociatedLeps, theGenZZLeps
                    

        genZZFinalState=0
        genZZMass=0.
        if (len(theGenZZLeps) == 4):
            genZZFinalState=1
            gen_p4 = TLorentzVector()
            for lep in theGenZZLeps :
                gen_p4 += genpart[lep].p4()
                genZZFinalState *= genpart[lep].pdgId
                #        print("Gen: {:} leps M={:.4g} In Acc: {:} e, {:} mu,  {:} FSR".format(len(GenHLeps),gen_p4.M(), len(GenEl_acc), len(GenMu_acc), len (GenFSR)), "\n")
            genZZMass = gen_p4.M()
        else:
            theGenZZLeps = [-1]*4

#        print(theGenH.p4().M(), genZZMass, "nFSR: ", len(genFSRIdxs))

        self.out.fillBranch("GenZZ_FinalState", genZZFinalState)
        self.out.fillBranch("GenZZ_mass", genZZMass)
        self.out.fillBranch("GenZZ_Z1l1Idx", theGenZZLeps[0]) #FIXME: to be sorted with standard criteria
        self.out.fillBranch("GenZZ_Z1l2Idx", theGenZZLeps[1])
        self.out.fillBranch("GenZZ_Z2l1Idx", theGenZZLeps[2])
        self.out.fillBranch("GenZZ_Z2l2Idx", theGenZZLeps[3])

        return True
