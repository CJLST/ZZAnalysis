##
# Report genXS and genBR in the output root files. 
##
from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class genXSFiller(Module):
    def __init__(self, GENXS, GENBR):
        print("***genXSFiller: GENXS:", GENXS, "GENBR:", GENBR,flush=True)
        self.writeHistFile = False
        self.GENXS = GENXS
        self.GENBR = GENBR




    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("genxsec", "F", title="The value of the cross section as reported by the generator used to make this sample.")
        self.out.branch("genbr", "F", title="The value of the HZZ branching ratio as reported by the generator used to make this sample.")


    def analyze(self, event):
        self.out.fillBranch("genxsec", float(self.GENXS))
        self.out.fillBranch("genbr", float(self.GENBR))

        return True


