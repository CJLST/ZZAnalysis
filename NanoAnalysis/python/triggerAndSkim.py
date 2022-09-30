####
# -Apply PV filter
# -Apply trigger requirements, and PD precedence rules for data
# FIXME: add variables to record which triggers were passed; passthrough option to check for missing triggers in data
####
from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import sys

class triggerAndSkim(Module):
    def __init__(self, isMC=True, PD="", era=2018):
        self.writeHistFile = False
        self.isMC = isMC
        self.PD = PD
        self.era = era
        print("triggerAndSkim configuration: IsMC=", self.isMC, "PD=", self.PD, "era=", self.era)
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        PD = self.PD

        ### Good PV filter
        if event.PV_npvsGood == 0 : return False


        ### Trigger requirements
        passTrigger = False
        if self.era == 2018:
            passSingleEle = event.HLT_Ele32_WPTight_Gsf
            passSingleMu = event.HLT_IsoMu24
            passDiEle = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL or event.HLT_DoubleEle25_CaloIdL_MW
            passDiMu = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8
            passMuEle = event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL or event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ or event.HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ
            passTriEle = False
            passTriMu = event.HLT_TripleMu_10_5_5_DZ or event.HLT_TripleMu_12_10_5
        else:         
            sys.exit("ERROR: era not supported: ", self.era)

        
        if self.isMC or PD == "any" :
            passTrigger = passDiEle or passDiMu or passMuEle or passTriEle or passTriMu or passSingleEle or passSingleMu

        else: # Data: ensure each event is taken only from a single PD
            if PD == "" : sys.exit("ERROR: PD must be set in data") # we may want to merge triggers for test runs 
            if ((PD=="DoubleEle" or PD=="DoubleEG"  or PD=="EGamma" ) and (passDiEle or passTriEle)) or \
               ((PD=="DoubleMu"  or PD=="DoubleMuon") and (passDiMu or passTriMu) and not passDiEle and not passTriEle) or \
               ((PD=="MuEG"      or PD=="MuonEG"    ) and passMuEle and not passDiMu and not passTriMu and not passDiEle and not passTriEle) or \
               ((PD=="SingleElectron" or PD=="EGamma") and passSingleEle and not passMuEle and not passDiMu and not passTriMu and not passDiEle and not passTriEle) or \
               ( PD=="SingleMuon" and passSingleMu and not passSingleEle and not passMuEle and not passDiMu and not passTriMu and not passDiEle and not passTriEle) :
                   passTrigger = True

        return passTrigger



