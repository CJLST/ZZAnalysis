# Legacy implementation based on root histograms, kept for older MC
# Imported from the standalone version of PhysicsTools.NanoAODTools.postprocessing.modules.common.
#
# NOTE: The helpers at the end set up the new correctionlib-based module from
# PhysicsTools/NATModules.

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
import os
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True

class puWeightProducer(Module):
    def __init__(self,
                 myfile,
                 targetfile,
                 myhist="pileup",
                 targethist="pileup",
                 name="puWeight",
                 norm=True,
                 verbose=False,
                 nvtx_var="Pileup_nTrueInt",
                 doSysVar=True
     ):
        self.targeth = self.loadHisto(targetfile, targethist)
        if doSysVar:
            self.targeth_plus = self.loadHisto(targetfile,
                                               targethist + "_plus")
            self.targeth_minus = self.loadHisto(targetfile,
                                                targethist + "_minus")
        self.fixLargeWeights = True  # temporary fix
        if myfile != "auto":
            self.autoPU = False
            self.myh = self.loadHisto(myfile, myhist)
        else:
            self.fixLargeWeights = False  # AR: it seems to crash with it, to be deugged
            self.autoPU = True
            ROOT.gROOT.cd()
            self.myh = self.targeth.Clone("autoPU")
            self.myh.Reset()
        self.name = name
        self.norm = norm
        self.verbose = verbose
        self.nvtxVar = nvtx_var
        self.doSysVar = doSysVar

        # Try to load module via python dictionaries
        try:
            ROOT.gSystem.Load("libZZAnalysisAnalysisStep")
            dummy = ROOT.WeightCalculatorFromHistogram
        # Load it via ROOT ACLIC. NB: this creates the object file in the
        # CMSSW directory, causing problems if many jobs are working from the
        # same CMSSW directory
        except Exception as e:
            print("Could not load module via python, trying via ROOT" + str(e))
            if "/WeightCalculatorFromHistogram_cc.so" not in ROOT.gSystem.GetLibraries(
            ):
                print("Load C++ Worker")
                ROOT.gROOT.ProcessLine(
                    ".L %s/src/PhysicsTools/NanoAODTools/src/WeightCalculatorFromHistogram.cc++"
                    % os.environ['CMSSW_BASE'])
            dummy = ROOT.WeightCalculatorFromHistogram

    def loadHisto(self, filename, hname):
        tf = ROOT.TFile.Open(filename)
        hist = tf.Get(hname)
        hist.SetDirectory(0)
        tf.Close()
        return hist

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        if self.autoPU:
            self.myh.Reset()
            print("Computing PU profile for this file")
            ROOT.gROOT.cd()
            inputFile.Get("Events").Project("autoPU",
                                            self.nvtxVar)  # doitfrom inputFile
            if outputFile:
                outputFile.cd()
                self.myh.Write()
        self._worker = ROOT.WeightCalculatorFromHistogram(
            self.myh, self.targeth, self.norm, self.fixLargeWeights,
            self.verbose)
        self.out = wrappedOutputTree
        self.out.branch(self.name, "F")
        if self.doSysVar:
            self._worker_plus = ROOT.WeightCalculatorFromHistogram(
                self.myh, self.targeth_plus, self.norm, self.fixLargeWeights,
                self.verbose)
            self._worker_minus = ROOT.WeightCalculatorFromHistogram(
                self.myh, self.targeth_minus, self.norm, self.fixLargeWeights,
                self.verbose)
            self.out.branch(self.name + "Up", "F")
            self.out.branch(self.name + "Down", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        if hasattr(event, self.nvtxVar):
            nvtx = int(getattr(event, self.nvtxVar))
            weight = self._worker.getWeight(
                nvtx) if nvtx < self.myh.GetNbinsX() else 1
            if self.doSysVar:
                weight_plus = self._worker_plus.getWeight(
                    nvtx) if nvtx < self.myh.GetNbinsX() else 1
                weight_minus = self._worker_minus.getWeight(
                    nvtx) if nvtx < self.myh.GetNbinsX() else 1
        else:
            weight = 1
        self.out.fillBranch(self.name, weight)
        if self.doSysVar:
            self.out.fillBranch(self.name + "Up", weight_plus)
            self.out.fillBranch(self.name + "Down", weight_minus)
        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

pufile_mc2016 = "%s/src/ZZAnalysis/NanoAnalysis/data/pileup_profile_Summer16.root" % os.environ[
    'CMSSW_BASE']
pufile_data2016 = "%s/src/ZZAnalysis/NanoAnalysis/data/PileupData_GoldenJSON_Full2016.root" % os.environ[
    'CMSSW_BASE']
puWeight_2016 = lambda: puWeightProducer(pufile_mc2016,
                                         pufile_data2016,
                                         "pu_mc",
                                         "pileup",
                                         verbose=False,
                                         doSysVar=True)
puAutoWeight_2016 = lambda: puWeightProducer(
    "auto", pufile_data2016, "pu_mc", "pileup", verbose=False)

pufile_data2017 = "%s/src/ZZAnalysis/NanoAnalysis/data/PileupHistogram-goldenJSON-13tev-2017-99bins_withVar.root" % os.environ[
    'CMSSW_BASE']
pufile_mc2017 = "%s/src/ZZAnalysis/NanoAnalysis/data/mcPileup2017.root" % os.environ[
    'CMSSW_BASE']
puWeight_2017 = lambda: puWeightProducer(pufile_mc2017,
                                         pufile_data2017,
                                         "pu_mc",
                                         "pileup",
                                         verbose=False,
                                         doSysVar=True)
puAutoWeight_2017 = lambda: puWeightProducer(
    "auto", pufile_data2017, "pu_mc", "pileup", verbose=False)

pufile_data2018 = "%s/src/ZZAnalysis/NanoAnalysis/data/PileupHistogram-goldenJSON-13tev-2018-100bins_withVar.root" % os.environ[
    'CMSSW_BASE']
pufile_mc2018 = "%s/src/ZZAnalysis/NanoAnalysis/data/mcPileup2018.root" % os.environ[
    'CMSSW_BASE']
puWeight_2018 = lambda: puWeightProducer(pufile_mc2018,
                                         pufile_data2018,
                                         "pu_mc",
                                         "pileup",
                                         verbose=False,
                                         doSysVar=True)
puAutoWeight_2018 = lambda: puWeightProducer(
    "auto", pufile_data2018, "pu_mc", "pileup", verbose=False)

# Ultra legacy

# 2016
pufile_dataUL2016 = "%s/src/ZZAnalysis/NanoAnalysis/data/PileupHistogram-UL2016-100bins_withVar.root" % os.environ[
    'CMSSW_BASE']
pufile_mcUL2016 = "%s/src/ZZAnalysis/NanoAnalysis/data/mcPileupUL2016.root" % os.environ[
    'CMSSW_BASE']
puWeight_UL2016 = lambda: puWeightProducer(pufile_mcUL2016,
                                           pufile_dataUL2016,
                                           "pu_mc",
                                           "pileup",
                                           verbose=False,
                                           doSysVar=True)
puAutoWeight_UL2016 = lambda: puWeightProducer(
    "auto", pufile_dataUL2016, "pu_mc", "pileup", verbose=False)

# 2017
pufile_dataUL2017 = "%s/src/ZZAnalysis/NanoAnalysis/data/PileupHistogram-UL2017-100bins_withVar.root" % os.environ[
    'CMSSW_BASE']
pufile_mcUL2017 = "%s/src/ZZAnalysis/NanoAnalysis/data/mcPileupUL2017.root" % os.environ[
    'CMSSW_BASE']
puWeight_UL2017 = lambda: puWeightProducer(pufile_mcUL2017,
                                           pufile_dataUL2017,
                                           "pu_mc",
                                           "pileup",
                                           verbose=False,
                                           doSysVar=True)
puAutoWeight_UL2017 = lambda: puWeightProducer(
    "auto", pufile_dataUL2017, "pu_mc", "pileup", verbose=False)

# 2018
pufile_dataUL2018 = "%s/src/ZZAnalysis/NanoAnalysis/data/PileupHistogram-UL2018-100bins_withVar.root" % os.environ[
    'CMSSW_BASE']
pufile_mcUL2018 = "%s/src/ZZAnalysis/NanoAnalysis/data/mcPileupUL2018.root" % os.environ[
    'CMSSW_BASE']
pufile_2022 = "%s/src/ZZAnalysis/AnalysisStep/data/PileUpWeights/pu_weights_2022.root" % os.environ['CMSSW_BASE']
puWeight_UL2018 = lambda: puWeightProducer(pufile_mcUL2018,
                                           pufile_dataUL2018,
                                           "pu_mc",
                                           "pileup",
                                           verbose=False,
                                           doSysVar=True)
puAutoWeight_UL2018 = lambda: puWeightProducer(
    "auto", pufile_dataUL2018, "pu_mc", "pileup", verbose=False)

# 2022 weights non official, merged pre and post-EE, obsolete (replaced by the split correctionlib versions, see below)
puWeight_2022 = lambda : puWeightProducer(pufile_2022,
                                          pufile_2022,
                                          "MC_out_of_the_box",
                                          "Data",
                                          verbose=False,
                                          doSysVar=True)

def puWeight(era, data_tag):
    print("***puWeight: era", era, "dataTag:", data_tag)
    if era == 2016:
        return puWeight_UL2016()
    elif era == 2017 : 
        return puWeight_UL2017()
    elif era == 2018 :
        return puWeight_UL2018()
    elif era == 2022 :
#        return puWeight_2022() # Merged pre and postEE - obsolete
        
        from PhysicsTools.NATModules.modules.puWeightProducer import puWeightProducer as puWeightProducer_corrlib
        if "pre_EE" in data_tag :
            json = "%s/src/ZZAnalysis/NanoAnalysis/data/puWeights_2022_Summer22.json.gz" % os.environ['CMSSW_BASE']
            key = "Collisions2022_355100_357900_eraBCD_GoldenJson"
        else :
            json = "%s/src/ZZAnalysis/NanoAnalysis/data/puWeights_2022_Summer22EE.json.gz" % os.environ['CMSSW_BASE']
            key = "Collisions2022_359022_362760_eraEFG_GoldenJson"
        return puWeightProducer_corrlib(json, key)
