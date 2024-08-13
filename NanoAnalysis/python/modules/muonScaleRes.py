from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import os
import numpy as np
from math import pi
from ROOT import MuonScaRe

class muonScaleRes(Module):
    def __init__(self, json, is_mc, overwritePt=False):
        """Add branches for muon scale and resolution corrections.
        Parameters:
            json: full path of json file
            is_mc: True for MC (smear pt+add uncertainties), False for data (scale pt)
            overwritePt: replace value in the pt branch, and store the old one as "uncorrected_pt"
        """

        self.overwritePt = overwritePt
        self.is_mc = is_mc
        
        self.corrModule = MuonScaRe(json)


    def getPtCorr(self, event, muon, var = "nom") :
        isData = int(not self.is_mc)
        scale_corr = self.corrModule.pt_scale(isData, muon.pt, muon.eta, muon.phi, muon.charge, var)
        pt_corr = scale_corr

        if self.is_mc:
            smear_corr = self.corrModule.pt_resol(scale_corr, muon.eta, muon.nTrackerLayers, var)
            pt_corr = smear_corr

        return pt_corr


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


    def analyze(self, event):
        if event.nMuon == 0 :
            return True

        muons = Collection(event, "Muon")

        pt_corr = [0.]*len(muons)
        if self.is_mc:
            pt_syst = [0.]*len(muons)
            pt_stat = [0.]*len(muons)
                
        for imu, muon in enumerate(muons):
            # Set up a deterministic random seed.
            # The seed is unique by event and muon.
            # A fixed entropy value is also included to decorrelate different modules doing similar things.
            seedSeq = np.random.SeedSequence([event.luminosityBlock, event.event, int(abs((muon.phi/pi*100.)%1)*1e10), 351740215])
            self.corrModule.setSeed(int(seedSeq.generate_state(1,np.uint64)[0]))

            pt_corr[imu] = self.getPtCorr(event, muon, "nom")

            if self.is_mc:
                # TODO: Check. Are we assuming up=dn?
                pt_syst[imu] = self.getPtCorr(event, muon, "syst")
                pt_stat[imu] = self.getPtCorr(event, muon, "stat")

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

