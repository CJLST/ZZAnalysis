from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import sys, os
import random
import numpy as np
import correctionlib

sys.path.append(f'{os.environ["CMSSW_BASE"]}/src/ZZAnalysis/NanoAnalysis/python/modules')
from MuonScaRe import pt_resol, pt_scale

class muonScaleResProducer(Module):
    def __init__(self, dataYear, tag, is_mc, overwritePt=False):
        print("***INIT muonScaleResProducer: dataYear:", dataYear, "tag:", tag, "is MC:", is_mc, "overwritePt:", overwritePt)
        self.tag = tag
        self.fname = self.get_fname(tag)
        self.fpath = ("%s/src/ZZAnalysis/NanoAnalysis/data/%s" % (os.environ['CMSSW_BASE'], self.fname))
        self.overwritePt = overwritePt
        self.is_mc = is_mc

        self.evaluator = self._get_evaluator

    @property
    def _get_evaluator(self):
        evaluator = correctionlib.CorrectionSet.from_file(self.fpath)
        return evaluator

    def get_fname(self, tag):
        # TODO: Generalize more also for 2023 corrections
        if "pre_EE" in tag :
            fname = "2022_schemaV2.json"
        else :
            fname = "2022EE_schemaV2.json"

        return fname

    def getPtCorr(self, event, muons, var = "nom"):
        is_data = 1
        if self.is_mc:
            is_data = 0

        scale_corr = pt_scale(is_data, np.array(list(event.Muon_pt)),
                              np.array(list(event.Muon_eta)),
                              np.array(list(event.Muon_phi)),
                              np.array(list(event.Muon_charge)),
                              var, self.evaluator, nested=False)
        pt_corr = scale_corr

        if self.is_mc:
            nTracks = np.array([float(mu.nTrackerLayers) for mu in muons])
            smear_corr = pt_resol(scale_corr, np.array(list(event.Muon_eta)),
                                  nTracks,
                                  var, self.evaluator, nested=False)
            pt_corr = smear_corr

        return pt_corr

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.overwritePt :
            self.out.branch("Muon_pt", "F", lenVar="nMuon")
            self.out.branch("Muon_uncorrected_pt", "F", lenVar="nMuon")
        else:
            self.out.branch("Muon_corrected_pt", "F", lenVar="nMuon")
        if self.is_mc:
            self.out.branch("Muon_syst_pt", "F", lenVar="nMuon")
            self.out.branch("Muon_stat_pt", "F", lenVar="nMuon")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        muons = Collection(event, "Muon")

        pt_corr = self.getPtCorr(event, muons, "nom")

        if self.is_mc:
            # TODO: Check. Are we assuming up=dn?
            pt_syst = self.getPtCorr(event, muons, "syst")
            pt_stat = self.getPtCorr(event, muons, "stat")

        if self.overwritePt :
            pt_uncorr = list(mu.pt for mu in muons)
            self.out.fillBranch("Muon_uncorrected_pt", pt_uncorr)
            self.out.fillBranch("Muon_pt", pt_corr)
        else :
            self.out.fillBranch("Muon_corrected_pt", pt_corr)

        if self.is_mc:
            self.out.fillBranch("Muon_syst_pt", pt_syst)
            self.out.fillBranch("Muon_stat_pt", pt_stat)

        return True

muonScaleRes = lambda era, tag, is_mc, overwritePt=False : muonScaleResProducer(era, tag, is_mc, overwritePt)