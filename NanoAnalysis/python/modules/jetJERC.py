from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import os
import random
import numpy as np
import correctionlib
import envyaml

corrCfg = envyaml.EnvYAML("%s/src/ZZAnalysis/NanoAnalysis/data/JERC/jetCorrections.yaml" % (os.environ['CMSSW_BASE']))

class jetJERC(Module):
    def __init__(self, dataYear, tag, is_mc, overwritePt=False):
        print("***INIT jetJERC: dataYear:", dataYear, "tag:", tag, "is MC:", is_mc, "overwritePt:", overwritePt)
        self.tag = tag
        self.fname_JERC = "jet_jerc.json"
        self.fname_JERsmear = "jer_smear.json"
        self.fpath_JERC = ("%s/src/ZZAnalysis/NanoAnalysis/data/JERC/%s/%s" % (os.environ['CMSSW_BASE'], corrCfg["tagsMC"][self.tag]["folder"], self.fname_JERC))
        self.fpath_JERsmear = ("%s/src/ZZAnalysis/NanoAnalysis/data/JERC/%s" % (os.environ['CMSSW_BASE'], self.fname_JERsmear))
        self.overwritePt = overwritePt
        self.is_mc = is_mc

        self.evaluator_JERC, self.evaluator_jer = self._get_evaluator
        self.evaluator_L1, self.evaluator_L2, self.evaluator_L3, self.evaluator_L2L3 = self._get_JERC_evaluator
        if self.is_mc: ## JER is applied only to MC
            self.evaluator_JERsmear, self.evaluator_JER, self.evaluator_JERsf = self._get_JER_evaluator

    @property
    def _get_evaluator(self):
        evaluator_JERC = correctionlib.CorrectionSet.from_file(self.fpath_JERC)
        evaluator_JER = correctionlib.CorrectionSet.from_file(self.fpath_JERsmear)
        return evaluator_JERC, evaluator_JER

    @property
    def _get_JER_evaluator(self):
        ## JERsmear
        JERsmear_evaluator = self.evaluator_jer["JERSmear"]
        ## JER
        JER_evaluator = self.evaluator_JERC[corrCfg["tagsMC"][self.tag]["prefix_JER"] + "_PtResolution_AK4PFPuppi"]
        ## JERsf
        JERsf_evaluator = self.evaluator_JERC[corrCfg["tagsMC"][self.tag]["prefix_JER"] + "_ScaleFactor_AK4PFPuppi"]
        return JERsmear_evaluator, JER_evaluator, JERsf_evaluator

    @property
    def _get_JERC_evaluator(self):
        if self.is_mc:
            yamlTag = "tagsMC"
        else:
            yamlTag = "tagsDATA"
        ## JEC
        ## List of sequential corrections are reported in https://cms-talk.web.cern.ch/t/jes-for-2022-re-reco-cde-and-prompt-fg/32873/4
        L1_evaluator = self.evaluator_JERC[corrCfg[yamlTag][self.tag]["prefix"] + "_L1FastJet_AK4PFPuppi"]
        L2_evaluator = self.evaluator_JERC[corrCfg[yamlTag][self.tag]["prefix"] + "_L2Relative_AK4PFPuppi"]
        L3_evaluator = self.evaluator_JERC[corrCfg[yamlTag][self.tag]["prefix"] + "_L3Absolute_AK4PFPuppi"]
        L2L3_evaluator = self.evaluator_JERC[corrCfg[yamlTag][self.tag]["prefix"] + "_L2L3Residual_AK4PFPuppi"]
        return L1_evaluator, L2_evaluator, L3_evaluator, L2L3_evaluator

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.overwritePt :
            self.out.branch("Jet_pt", "F", lenVar="nJet")
            self.out.branch("Jet_mass", "F", lenVar="nJet")
            self.out.branch("Jet_uncorrected_pt", "F", lenVar="nJet")
            self.out.branch("Jet_uncorrected_mass", "F", lenVar="nJet")
        else:
            self.out.branch("Jet_corrected_pt", "F", lenVar="nJet")
            self.out.branch("Jet_corrected_mass", "F", lenVar="nJet")
        if self.is_mc:
            self.out.branch("Jet_scaleUp_pt", "F", lenVar="nJet")
            self.out.branch("Jet_scaleDn_pt", "F", lenVar="nJet")
            self.out.branch("Jet_smearUp_pt", "F", lenVar="nJet")
            self.out.branch("Jet_smearDn_pt", "F", lenVar="nJet")
            self.out.branch("Jet_smearUp_mass", "F", lenVar="nJet")
            self.out.branch("Jet_smearDn_mass", "F", lenVar="nJet")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def fixPhi(self, phi):
        if phi > np.pi:
            phi -= 2*np.pi
        elif phi < -np.pi:
            phi += 2*np.pi
        return phi

    def analyze(self, event):
        jets = Collection(event, "Jet")
        if self.is_mc: ## genJet info is necessary for JER
            gen_jets = Collection(event, "GenJet")
            gen_jets_pt = np.array([gen_jet.pt for gen_jet in gen_jets])
            gen_jets_eta = np.array([gen_jet.eta for gen_jet in gen_jets])
            gen_jets_phi = np.array([gen_jet.phi for gen_jet in gen_jets])

        pt_corr = []
        pt_uncorr = []
        mass_corr = []
        mass_uncorr = []
        pt_smear_up = []
        pt_smear_dn = []
        mass_smear_up = []
        mass_smear_dn = []
        pt_scale_up = []
        pt_scale_dn = []

        for jet in jets:
            if self.is_mc:
                #### JEC ####
                ## Jet in NanoAOD are already corrected
                ## The correction should be removed and the latest one available should be applied
                pt_raw = jet.pt * (1 - jet.rawFactor)
                mass_raw = jet.mass * (1 - jet.rawFactor)
                ## The three steps of JEC corrections are provided separately
                pt_L1 = pt_raw * self.evaluator_L1.evaluate(jet.area, jet.eta, pt_raw, event.Rho_fixedGridRhoFastjetAll)
                pt_L2 = pt_L1 * self.evaluator_L2.evaluate(jet.eta, pt_L1)
                pt_L3 = pt_L2 * self.evaluator_L3.evaluate(jet.eta, pt_L2)
                pt_JEC = pt_L3 * self.evaluator_L2L3.evaluate(jet.eta, pt_L3)
                JEC = pt_JEC / pt_raw
                mass_JEC = mass_raw * JEC

                #### JER ####
                ## Hybrid method is implemented [https://cms-jerc.web.cern.ch/JER/#smearing-procedures]
                JER = self.evaluator_JER.evaluate(jet.eta, pt_JEC, event.Rho_fixedGridRhoFastjetAll)
                ## GenMatching with genJet
                delta_eta = jet.eta - gen_jets_eta
                fixPhi = np.vectorize(self.fixPhi)
                delta_phi = fixPhi(jet.phi - gen_jets_phi)
                pt_gen = np.where((np.abs(pt_JEC - gen_jets_pt) < 3 * pt_JEC * JER) & (np.sqrt(delta_eta**2 + delta_phi**2)<0.2), gen_jets_pt, -1.0)
                pt_gen = pt_gen[pt_gen > 0][0] if np.any(pt_gen > 0) else -1. ## If no gen-matching, simply -1
                JERsf = self.evaluator_JERsf.evaluate(jet.eta, jet.pt, "nom")
                JERsf_up = self.evaluator_JERsf.evaluate(jet.eta, jet.pt, "up")
                JERsf_dn = self.evaluator_JERsf.evaluate(jet.eta, jet.pt, "down")
                JERsmear = self.evaluator_JERsmear.evaluate(pt_JEC, jet.eta, pt_gen, event.Rho_fixedGridRhoFastjetAll, event.event, JER, JERsf)
                JERsmear_up = self.evaluator_JERsmear.evaluate(pt_JEC, jet.eta, pt_gen, event.Rho_fixedGridRhoFastjetAll, event.event, JER, JERsf_up)
                JERsmear_dn = self.evaluator_JERsmear.evaluate(pt_JEC, jet.eta, pt_gen, event.Rho_fixedGridRhoFastjetAll, event.event, JER, JERsf_dn)
                pt_JEC_JER = pt_JEC * JERsmear
                pt_JEC_JER_up = pt_JEC * JERsmear_up
                pt_JEC_JER_dn = pt_JEC * JERsmear_dn
                mass_JEC_JER = mass_JEC * JERsmear
                mass_JEC_JER_up = mass_JEC * JERsmear_up
                mass_JEC_JER_dn = mass_JEC * JERsmear_dn


                pt_corr.append(pt_JEC_JER)
                pt_uncorr.append(pt_raw)
                mass_corr.append(mass_JEC_JER)
                mass_uncorr.append(mass_raw)
                pt_smear_up.append(pt_JEC_JER_up)
                pt_smear_dn.append(pt_JEC_JER_dn)
                mass_smear_up.append(mass_JEC_JER_up)
                mass_smear_dn.append(mass_JEC_JER_dn)
                pt_scale_up.append(8)
                pt_scale_dn.append(8)

            else:
                #### JEC ####
                ## Jet in NanoAOD are already corrected
                ## The correction should be removed and the latest one available should be applied
                pt_raw = jet.pt * (1 - jet.rawFactor)
                mass_raw = jet.mass * (1 - jet.rawFactor)
                ## The three steps of JEC corrections are provided separately
                pt_L1 = pt_raw * self.evaluator_L1.evaluate(jet.area, jet.eta, pt_raw, event.Rho_fixedGridRhoFastjetAll)
                pt_L2 = pt_L1 * self.evaluator_L2.evaluate(jet.eta, pt_L1)
                pt_L3 = pt_L2 * self.evaluator_L3.evaluate(jet.eta, pt_L2)
                pt_JEC = pt_L3 * self.evaluator_L2L3.evaluate(jet.eta, pt_L3)
                JEC = pt_JEC / pt_raw
                mass_JEC = mass_raw * JEC

                pt_corr.append(pt_JEC)
                pt_uncorr.append(pt_raw)
                mass_corr.append(mass_JEC)
                mass_uncorr.append(mass_raw)

        if self.overwritePt :
            self.out.fillBranch("Jet_uncorrected_pt", pt_uncorr)
            self.out.fillBranch("Jet_pt", pt_corr)
            self.out.fillBranch("Jet_uncorrected_mass", mass_uncorr)
            self.out.fillBranch("Jet_mass", mass_corr)
        else :
            self.out.fillBranch("Jet_corrected_pt", pt_corr)
            self.out.fillBranch("Jet_corrected_mass", mass_corr)

        if self.is_mc:
            self.out.fillBranch("Jet_smearUp_pt", pt_smear_up)
            self.out.fillBranch("Jet_smearDn_pt", pt_smear_dn)
            self.out.fillBranch("Jet_smearUp_mass", mass_smear_up)
            self.out.fillBranch("Jet_smearDn_mass", mass_smear_dn)
            self.out.fillBranch("Jet_scaleUp_pt", pt_scale_up)
            self.out.fillBranch("Jet_scaleDn_pt", pt_scale_dn)

        return True

jetCorrected = lambda era, tag, is_mc, overwritePt=False : jetJERC(era, tag, is_mc, overwritePt)
