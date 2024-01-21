### 
# Clone the branches listed in "varlist" to a new tree named "treeName" in the output file,
# and stop processing event that do not satisfy the condition "continueFor".
# This is useul in particular to store gen variables for all events in a separate tree.
###

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class cloneBranches(Module):

    def __init__(self, treeName, varlist, continueFor = lambda evt : (True)) :
        self.treeName = treeName
        self.varlist = varlist
        self.continueFor = continueFor

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.inputTree = inputTree
        # Clone the requested branches into the new tree.
        eventsTree = wrappedOutputTree.tree()
        eventsTree.SetBranchStatus("*", 0)
        for var in self.varlist :
            eventsTree.SetBranchStatus(var, 1)
        self.newTree = eventsTree.CloneTree(0)
        self.newTree.SetName(self.treeName)
        eventsTree.SetBranchStatus("*", 1)
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree) :
        self.newTree.Write()
        
    def analyze(self, event) :
        self.inputTree.readAllBranches()        
        self.newTree.Fill()
        
        return self.continueFor(event) 
