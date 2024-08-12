## This is derived from  PhysicsTools/NanoAODTools/python/postprocessing/modules/common/muonScaleResProducer.py
## with modifications to:
## - add sync mode
## - replace lepton pT instead of adding a new variable (overwritePt=true)
## - take RoccoR from the CJLST library instead than recompiling it

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import os
import random
from ctypes import CDLL, c_double, c_int, c_float, c_char_p
#import ROOT
#ROOT.PyConfig.IgnoreCommandLineOptions = True

def mk_safe(fct, *args):
    try:
        return fct(*args)
    except Exception as e:
        if any('Error in function boost::math::erf_inv' in arg
               for arg in e.args):
            print('WARNING: catching exception and returning -1. Exception arguments: %s'% e.args)
            return -1.
        else:
            raise e


class muonScaleResProducer(Module):
    def __init__(self, dataYear, tag, overwritePt=False, syncMode=False):
#        p_postproc = '%s/src/PhysicsTools/NanoAODTools/python/postprocessing' % os.environ['CMSSW_BASE']
#        p_roccor = p_postproc + '/data/' + rc_dir

        if dataYear<=2018 and "UL" in tag :
            filename = 'RoccoR' + str(dataYear) + "UL.txt"
        else : 
            filename = 'RoccoR' + str(dataYear) + ".txt"
        fpath = ('%s/src/ZZAnalysis/AnalysisStep/data/RochesterCorrections/' % os.environ['CMSSW_BASE'])+filename

## Orginal implementation: load using ROOT       
#         if "/RoccoR_cc.so" not in ROOT.gSystem.GetLibraries():
#             p_helper = '%s/RoccoR.cc' % p_roccor
#             print('Loading C++ helper from ' + p_helper)
#             ROOT.gROOT.ProcessLine('.L ' + p_helper)
#         self._roccor = ROOT.RoccoR(p_roccor + '/' + rc_corrections)

        self.lib = CDLL('pluginZZAnalysisAnalysisStepPlugins.so')
        self._roccor = self.lib.get_RoccoR(c_char_p(fpath.encode("ascii")))
        self.lib.RoccoR_kSpreadMC.restype = c_double
        self.lib.RoccoR_kSmearMC.restype = c_double
        self.lib.RoccoR_kSpreadMCerror.restype = c_double
        self.lib.RoccoR_kSmearMCerror.restype = c_double
        self.lib.RoccoR_kScaleDT.restype = c_double
        self.lib.RoccoR_kScaleDTerror.restype = c_double
        self.overwritePt = overwritePt
        self.syncMode = syncMode

        print("***INIT muonScaleResProducer: dataYear:", dataYear, "tag:", tag, "overwritePt:", overwritePt, "syncMode:", syncMode)

    def beginJob(self):
        pass

    def endJob(self):
        self.lib.del_RoccoR(self._roccor)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.overwritePt :
            self.out.branch("Muon_pt", "F", lenVar="nMuon")
            self.out.branch("Muon_uncorrected_pt", "F", lenVar="nMuon")
        else:
            self.out.branch("Muon_corrected_pt", "F", lenVar="nMuon")
        self.out.branch("Muon_correctedUp_pt", "F", lenVar="nMuon")
        self.out.branch("Muon_correctedDown_pt", "F", lenVar="nMuon")
        self.is_mc = bool(inputTree.GetBranch("GenJet_pt"))

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        muons = Collection(event, "Muon")
        if self.is_mc:
            genparticles = Collection(event, "GenPart")
        roccor = self._roccor
        if self.is_mc:
            pt_corr = []
            pt_err = []
            for mu in muons:
                mu_pt = c_double(mu.pt)
                mu_eta = c_double(mu.eta)
                mu_phi = c_double(mu.phi)
                genIdx = mu.genPartIdx
                if genIdx >= 0 and genIdx < len(genparticles):
                    genMu = genparticles[genIdx]
                    pt_corr.append(mu.pt *
                                   mk_safe(self.lib.RoccoR_kSpreadMC, roccor, mu.charge, mu_pt,
                                           mu_eta, mu_phi, c_double(genMu.pt), 0, 0))
                    pt_err.append(mu.pt *
                                  mk_safe(self.lib.RoccoR_kSpreadMCerror, roccor, mu.charge,
                                          mu_pt, mu_eta, mu_phi, c_double(genMu.pt)))
                else:
                    u1 = c_double(0.5)
                    if not self.syncMode : u1 = c_double(random.uniform(0.0, 1.0))
                    pt_corr.append(
                        mu.pt * mk_safe(self.lib.RoccoR_kSmearMC, roccor, mu.charge, mu_pt,
                                        mu_eta, mu_phi, mu.nTrackerLayers, u1, 0, 0))
                    pt_err.append(
                        mu.pt * mk_safe(self.lib.RoccoR_kSmearMCerror, roccor, mu.charge, mu_pt,
                                        mu_eta, mu_phi, mu.nTrackerLayers, u1))

        else:
            pt_corr = list(
                mu.pt *
                mk_safe(self.lib.RoccoR_kScaleDT, roccor, mu.charge, c_double(mu.pt), c_double(mu.eta), c_double(mu.phi), 0, 0)
                for mu in muons)
            pt_err = list(
                mu.pt *
                mk_safe(self.lib.RoccoR_kScaleDTerror, roccor, mu.charge, c_double(mu.pt), c_double(mu.eta), c_double(mu.phi))
                for mu in muons)

        if self.overwritePt :
            pt_uncorr = list(mu.pt for mu in muons)
            self.out.fillBranch("Muon_uncorrected_pt", pt_uncorr)
            self.out.fillBranch("Muon_pt", pt_corr)
        else :
            self.out.fillBranch("Muon_corrected_pt", pt_corr)
        pt_corr_up = list(
            max(pt_corr[imu] + pt_err[imu], 0.0)
            for imu, mu in enumerate(muons))
        pt_corr_down = list(
            max(pt_corr[imu] - pt_err[imu], 0.0)
            for imu, mu in enumerate(muons))
        self.out.fillBranch("Muon_correctedUp_pt", pt_corr_up)
        self.out.fillBranch("Muon_correctedDown_pt", pt_corr_down)
        return True


muonScaleRes = lambda era, tag, overwritePt=False, syncMode=False : muonScaleResProducer(era, tag, overwritePt, syncMode)
