from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import os
import random
import numpy as np
import correctionlib

class eleScaleResProducer(Module):
    def __init__(self, dataYear, tag, is_mc, overwritePt=False):

        # TODO: this is for promptreco FG. Load correctly (using DATA_TAG should be easy).
        self.fpath = ("%s/src/ZZAnalysis/NanoAnalysis/data/SS.json" % os.environ['CMSSW_BASE'])
        self.overwritePt = overwritePt
        self.is_mc = is_mc
        print("***INIT eleScaleResProducer: dataYear:", dataYear, "tag:", tag, "overwritePt:", overwritePt)

    @property    
    def createEvaluator(self):
        evaluator = correctionlib.CorrectionSet.from_file(self.fpath)
        print(list(evaluator.keys()))

        return evaluator

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.overwritePt :
            self.out.branch("Electron_pt", "F", lenVar="nElectron")
            self.out.branch("Electron_uncorrected_pt", "F", lenVar="nElectron")
        else:
            self.out.branch("Electron_corrected_pt", "F", lenVar="nElectron")
        # self.out.branch("Electron_correctedUp_pt", "F", lenVar="nElectron")
        # self.out.branch("Electron_correctedDown_pt", "F", lenVar="nElectron")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        electrons = Collection(event, "Electron")
        evaluator = self.createEvaluator

        evaluator_smearing = evaluator["Prompt2022FG_SmearingJSON"] # TODO: Correct loading
        evaluator_scale = evaluator["Prompt2022FG_ScaleJSON"] # TODO: Correct loading

        rng = np.random.default_rng(seed=125)

        pt_corr = []
        for ele in electrons:
            if self.is_mc:
                rho = evaluator_smearing.evaluate("rho", ele.eta, ele.r9)
                smearing = rng.normal(loc=1., scale=rho)
                pt_corr.append(smearing * ele.pt)
                print("Nominal/Smeared pT for this lepton (MC):", ele.pt/(smearing * ele.pt))

                # TODO: Uncertainties

            else:
                # TODO: Check if cast to float is safe for event.run
                scale = evaluator_scale.evaluate("total_correction", ele.seedGain, float(event.run), ele.eta, ele.r9, ele.pt)
                pt_corr.append(scale * ele.pt)
                print("Nominal/Smeared pT for this lepton (DATA):", ele.pt/(scale * ele.pt))

                # TODO: Uncertainties

        if self.overwritePt :
            pt_uncorr = list(ele.pt for ele in electrons)
            self.out.fillBranch("Electron_uncorrected_pt", pt_uncorr)
            self.out.fillBranch("Electron_pt", pt_corr)
        else :
            self.out.fillBranch("Electron_corrected_pt", pt_corr)

        return True

eleScaleRes = lambda era, tag, is_mc, overwritePt=False : eleScaleResProducer(era, tag, is_mc, overwritePt)
