import FWCore.ParameterSet.Config as cms

### ----------------------------------------------------------------------
# 2012 skim (cf https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsZZ4l2012Summer), as implemented by Giovanni.
# Original selection is with MuonSelector/"muons" and ElectronSelector/gsfElectrons.
# Note that our collections are pre-filtered with a cut at 2 GeV.
### ----------------------------------------------------------------------

muons4skim = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuonsWithTrigger::PAT"),  
    cut = cms.string("(isTrackerMuon||isGlobalMuon) && pt>3 && abs(eta) < 2.4"),
)
electrons4skim = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("patElectronsWithTrigger::PAT"), 
    cut = cms.string("pt>5 && abs(eta) < 2.5"),
)


leptons4skim = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag( cms.InputTag("muons4skim"),
                         cms.InputTag("electrons4skim"), )
)
dileptons4skim = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('leptons4skim leptons4skim'),
    cut = cms.string('deltaR(daughter(0).eta,daughter(0).phi,daughter(1).eta,daughter(1).phi)> 0.01'),
    checkCharge = cms.bool(False)
)

# any (charge/flavor) pair of leptons with PT>20,10 GeV
skimL20L10 = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("dileptons4skim"),
    cut = cms.string('min(daughter(0).pt,daughter(1).pt) > 10 && max(daughter(0).pt,daughter(1).pt) > 20'),
    filter = cms.bool(True),
)

#at least one pair of SF leptons with m(ll)>40 GeV 
skim40NoOF  = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("dileptons4skim"),
    cut = cms.string('mass > 40 && abs(daughter(0).pdgId) == abs(daughter(1).pdgId)'), ## and SF only
    filter = cms.bool(True),
)


primaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string('!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2'),
    filter = cms.bool(True),
)


HZZSkim2012 = cms.Sequence(primaryVertexFilter + muons4skim + electrons4skim + leptons4skim + dileptons4skim +
                           skimL20L10 + skim40NoOF)


