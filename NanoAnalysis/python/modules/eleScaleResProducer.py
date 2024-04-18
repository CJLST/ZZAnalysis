from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import os
import random
import numpy as np
import correctionlib

class eleScaleResProducer(Module):
    def __init__(self, dataYear, tag, is_mc, overwritePt=False):
        self.fname = self.get_fname(tag)
        self.fpath = ("%s/src/ZZAnalysis/NanoAnalysis/data/%s" % (os.environ['CMSSW_BASE'], self.fname))
        self.overwritePt = overwritePt
        self.is_mc = is_mc

        self.rng = np.random.default_rng()
        self.evaluator = self._get_evaluator
        self.evaluator_scale = self._get_scale_evaluator
        self.evaluator_smear = self._get_smear_evaluator

        print("***INIT eleScaleResProducer: dataYear:", dataYear, "tag:", tag, "is MC:", is_mc, "overwritePt:", overwritePt)
        print("***INIT eleScaleResProducer for: ", self.fpath)

    @property
    def _get_evaluator(self):
        evaluator = correctionlib.CorrectionSet.from_file(self.fpath)
        return evaluator

    @property    
    def _get_smear_evaluator(self):
        if "preEE" in self.fname:
            smear_evaluator = self.evaluator["2022Re-recoBCD_SmearingJSON"]
        elif "postEE" in self.fname:
            smear_evaluator = self.evaluator["2022Re-recoE+PromptFG_SmearingJSON"]

        return smear_evaluator

    @property
    def _get_scale_evaluator(self):
        if "preEE" in self.fname:
            scale_evaluator = self.evaluator["2022Re-recoBCD_ScaleJSON"]
        elif "postEE" in self.fname:
            scale_evaluator = self.evaluator["2022Re-recoE+PromptFG_ScaleJSON"]

        return scale_evaluator

    def get_fname(self, tag):
        if tag=="pre_EE":
            fname = "electronSS_preEE.json"
        elif tag=="post_EE":
            fname = "electronSS_postEE.json"
        else:
            fname = ""
            raise ValueError('[eleScaleResProducer]: specify 2022 period!')

        return fname

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
        if self.is_mc:
            self.out.branch("Electron_scaleUp_pt", "F", lenVar="nElectron")
            self.out.branch("Electron_scaleDn_pt", "F", lenVar="nElectron")
            self.out.branch("Electron_smearUp_pt", "F", lenVar="nElectron")
            self.out.branch("Electron_smearDn_pt", "F", lenVar="nElectron")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        electrons = Collection(event, "Electron")

        pt_corr = []
        pt_smear_up = []
        pt_smear_dn = []
        pt_scale_up = []
        pt_scale_dn = []

        for ele in electrons:
            if self.is_mc:
                rho = self.evaluator_smear.evaluate("rho", ele.eta, ele.r9)
                smearing = self.rng.normal(loc=1., scale=rho)
                pt_corr.append(smearing * ele.pt)

                unc_rho = self.evaluator_smear.evaluate("err_rho", ele.eta, ele.r9)
                smearing_up = self.rng.normal(loc=1., scale=rho + unc_rho)
                smearing_dn = self.rng.normal(loc=1., scale=rho - unc_rho)
                pt_smear_up.append(smearing_up * ele.pt)
                pt_smear_dn.append(smearing_dn * ele.pt)

                scale_MC_unc = self.evaluator_scale.evaluate("total_uncertainty", ele.seedGain, float(event.run), ele.eta, ele.r9, ele.pt)
                pt_scale_up.append((1+scale_MC_unc) * ele.pt)
                pt_scale_dn.append((1-scale_MC_unc) * ele.pt)
            else:
                scale = self.evaluator_scale.evaluate("total_correction", ele.seedGain, float(event.run), ele.eta, ele.r9, ele.pt)
                pt_corr.append(scale * ele.pt)

        if self.overwritePt :
            pt_uncorr = list(ele.pt for ele in electrons)
            self.out.fillBranch("Electron_uncorrected_pt", pt_uncorr)
            self.out.fillBranch("Electron_pt", pt_corr)
        else :
            self.out.fillBranch("Electron_corrected_pt", pt_corr)

        if self.is_mc:
            self.out.fillBranch("Electron_smearUp_pt", pt_smear_up)
            self.out.fillBranch("Electron_smearDn_pt", pt_smear_dn)
            self.out.fillBranch("Electron_scaleUp_pt", pt_scale_up)
            self.out.fillBranch("Electron_scaleDn_pt", pt_scale_dn)

        return True

eleScaleRes = lambda era, tag, is_mc, overwritePt=False : eleScaleResProducer(era, tag, is_mc, overwritePt)
