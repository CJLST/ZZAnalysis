import FWCore.ParameterSet.Config as cms
### ----------------------------------------------------------------------
### HZZ skim (2011 version)
### ----------------------------------------------------------------------

# FIXME: to be, updated cf. ~psilva/public/HZZSkim/

# cuts
HZZSKIM_MUON_CUT=("pt > 7 && abs(eta)<2.5 && (isGlobalMuon || isTrackerMuon)")
HZZSKIM_ELECTRON_CUT=("pt > 10 && abs(eta)<2.5")
HZZSKIM_DIMUON_CUT=("mass > 40 && max(daughter(0).pt, daughter(1).pt) > 20 && min(daughter(0).pt, daughter(1).pt) > 7")
HZZSKIM_DIELECTRON_CUT=("mass > 40 && max(daughter(0).pt, daughter(1).pt) > 20 && min(daughter(0).pt, daughter(1).pt) > 10")
HZZSKIM_EMU_CUT=("mass > 40 && ((daughter(0).pt>7 && daughter(1).pt()>20) || (daughter(0).pt>20 && daughter(1).pt()>10))")

# single lepton selectors
goodHzzMuons = cms.EDFilter("PATMuonRefSelector",
                            src = cms.InputTag("patMuonsWithTrigger"),
                            cut = cms.string(HZZSKIM_MUON_CUT)
                            )
goodHzzElectrons = cms.EDFilter("PATElectronRefSelector",
                                src = cms.InputTag("patElectronsWithTrigger"),
                                cut = cms.string(HZZSKIM_ELECTRON_CUT)
                                )

# dilepton selectors
diHzzMuons = cms.EDProducer("CandViewShallowCloneCombiner",
                            decay       = cms.string("goodHzzMuons goodHzzMuons"),
                            checkCharge = cms.bool(False),
                            cut         = cms.string(HZZSKIM_DIMUON_CUT)
                            )
diHzzElectrons = cms.EDProducer("CandViewShallowCloneCombiner",
                                decay       = cms.string("goodHzzElectrons goodHzzElectrons"),
                                checkCharge = cms.bool(False),
                                cut         = cms.string(HZZSKIM_DIELECTRON_CUT)
                                )
crossHzzLeptons  = cms.EDProducer("CandViewShallowCloneCombiner",
                                  decay       = cms.string("goodHzzMuons goodHzzElectrons"),
                                  checkCharge = cms.bool(False),
                                  cut         = cms.string(HZZSKIM_EMU_CUT)
                                  )

dilepHzz = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag(cms.InputTag("diHzzMuons"), cms.InputTag("diHzzElectrons"),cms.InputTag("crossHzzLeptons"),)
)

# The actual filter
HzzSkim= cms.EDFilter("CandViewCountFilter",
                                src = cms.InputTag("dilepHzz"),
                                minNumber = cms.uint32(1)
                                )

skim2011 = cms.Sequence(
    goodHzzMuons +
    goodHzzElectrons +
    diHzzMuons +
    diHzzElectrons +
    crossHzzLeptons +
    dilepHzz +
    HzzSkim 
  )
