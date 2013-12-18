#$Id: ZZ4lEventContent.py,v 1.2 2012/10/17 07:59:11 argiro Exp $

import FWCore.ParameterSet.Config as cms

ZZ4lEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring(
    'drop *',
    'keep *_addPileupInfo_*_*',
    'keep *_TriggerResults_*_*',
    'keep *_EECand_*_*',
    'keep *_EEEECand_*_*',
    'keep *_EEMMCand_*_*',
    'keep *_MMCand_*_*',
    'keep *_MMMMCand_*_*',
    'keep *_ZZCand_*_*',
    'keep *_softMuons_*_*',
    'keep *_appendPhotons_*_*',
    'keep *_goodPrimaryVertices_*_*',
    'keep *_cmgPhotonSel_*_*',
    'keep *_cmgPFJetSel_*_*',
    'keep *_cmgPFMetSel_*_*',
    'keep *_preSkimCounter_*_*',
    'keep *_prePathCounter_*_*',
    )
    )
