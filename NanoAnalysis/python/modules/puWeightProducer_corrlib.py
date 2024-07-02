###
# Compute PU weight using correctionlib, and store in new branch.
# Load as
#  puWeight = puWeightProducer("POG/LUM/2022_Summer22EE/puWeights.json.gz", "Collisions2022_359022_362760_eraEFG_GoldenJson") 
###
from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from correctionlib import CorrectionSet

class puWeightProducer(Module):
    def __init__(self, json, key, name="puWeight", doSysVar=True) :
        """Add weights.
        Parameters:
            json: full path of json file
            key: key for the table in the file
            name: name of variable to be added 
            doSysVar: add up/dn variations
        """
        self.name = name
        self.nameUp = name+"Up"
        self.nameDn = name+"Dn"
        self.doSysVar = doSysVar

        self.evaluator = CorrectionSet.from_file(json)[key]
    
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch(self.name, "F")
        if self.doSysVar:
            self.out.branch(self.nameUp, "F")
            self.out.branch(self.nameDn, "F")
    
    def analyze(self, event):        
        self.out.fillBranch(self.name, self.evaluator.evaluate(event.Pileup_nTrueInt, "nominal"))
        if self.doSysVar:
            self.out.fillBranch(self.nameUp, self.evaluator.evaluate(event.Pileup_nTrueInt, "up"))
            self.out.fillBranch(self.nameDn, self.evaluator.evaluate(event.Pileup_nTrueInt, "down"))

        return True
    
