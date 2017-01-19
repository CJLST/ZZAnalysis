#
# Define full paths for ZZ candidates with one trackless electron (SR and CRs)
#

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDPhotonIdProducer(process, dataFormat)
my_id_modules = [
                'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_25ns_nonTrig_V2_cff',
                'RecoEgamma.PhotonIdentification.Identification.mvaTLEID_Fall15_V1_cff',
                ]
#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)


process.appendPhotons.tleSrc = cms.InputTag("softPhotons")
process.trackless_electrons = cms.Sequence(process.bareSoftPhotons + process.egmPhotonIDSequence + process.softPhotons)
process.electrons += process.trackless_electrons

# l+l- (SFOS, both e and mu)

process.bareZCandtle = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('appendPhotons:electrons appendPhotons:electronstle'),
    cut = cms.string(''),#mass>40 && mass <180'), # protect against ghosts
    checkCharge = cms.bool(False)
)


process.ZCandtle = cms.EDProducer("ZCandidateFiller",
    src = cms.InputTag("bareZCandtle"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    bestZAmong = cms.string(BESTZ_AMONG),
    FSRMode = cms.string(FSRMODE), # "skip", "Legacy", "RunII"
    flags = cms.PSet(
        GoodLeptons = cms.string(ZLEPTONSEL),
        Z1Presel = cms.string(Z1PRESEL),
        passCombRelIsoPFFSRCorr = cms.string("abs(1)"),
        CombRelIsoPFFSRCorr = cms.string("abs(0)"),
    )
)


process.bareZZCandtle= cms.EDProducer("PATCandViewShallowCloneCombiner",
    decay = cms.string('ZCand ZCandtle'),
    cut = cms.string(LLLLPRESEL),
    checkCharge = cms.bool(False)
)


process.ZZCandtle = cms.EDProducer("ZZCandidateFiller",
    src = cms.InputTag("bareZZCandtle"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP),
    superMelaMass = cms.double(SUPERMELA_MASS),
    isMC = cms.bool(IsMC),
    bestCandAmong = cms.PSet(isBestCand = cms.string(BESTCAND_AMONG)),
    bestCandComparator = cms.string(BESTCANDCOMPARATOR),
    ZRolesByMass = cms.bool(True),
    doKinFit = cms.bool(KINREFIT),
    flags = cms.PSet(
        GoodLeptons =  cms.string(FOURGOODLEPTONS),
        Z2Mass  = cms.string(Z2MASS),
        MAllComb = cms.string(MLLALLCOMB),
        SR = cms.string(SR),
        FullSel70 = cms.string(SR), #Obsolete, use "SR"
        FullSel = cms.string(FULLSEL),
    ),
    recoProbabilities = cms.vstring(),
    muon_iso_cut = cms.double(MUISOCUT),
    electron_iso_cut = cms.double(ELEISOCUT),
)

### Trees for control regions only
# what to do in CR where charge is needed?

# Z (OSSF,both e/mu) + LL (any F/C, with no ID/iso); this is the starting point for control regions

process.bareZLLCandtle= cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('ZCand LLCandtle'),
    cut = cms.string(NOGHOST4l),
    checkCharge = cms.bool(False)
)


Z2LL_tle = "abs(daughter(1).daughter(0).pdgId * daughter(1).daughter(1).pdgId) == 242"   #Z2 = e * e w/o track
Z2LL_SS_tle = Z2LL_tle #"daughter(1).daughter(0).pdgId()==daughter(1).daughter(1).pdgId()"       #Z2 = same-sign, same-flavour
Z2LL_OS_tle = "abs(1)"

if SELSETUP == "allCutsAtOncePlusSmart" :
    CR_BESTZLLss_tle = CR_BESTCANDBASE_AA + "&&" + Z2LL_SS_tle + "&&" +CR_Z2MASS + "&&" + MLLALLCOMB + "&&" + PT20_10 + "&&" + "mass>70" + "&&" + "daughter(1).mass>12" + "&&" + SMARTMALLCOMB


# 2P2L region
CR_BESTZLLos_tle = (CR_BESTCANDBASE_AA    + "&&" +  
                CR_BASESEL            + "&&" +
                Z2LL_OS_tle               + "&&" +  
                SMARTMALLCOMB         )

CR_BESTZLL_tle = CR_BESTCANDBASE_AA + "&&" +CR_Z2MASS + "&&" + MLLALLCOMB + "&&" + PT20_10 + "&&" + "mass>70" + "&&" + "daughter(1).mass>12" + "&&" + SMARTMALLCOMB


# CR 3P1F
CR_BESTZLLos_3P1F_tle = (CR_BESTZLLos_tle + "&&" + PASSD0_OR_PASSD1)                 
CR_ZLLosSEL_3P1F_tle  = (CR_BESTZLLos_tle + "&&" + PASSD0_XOR_PASSD1) # Is the CR_BESTZLLos request redundant? 


# CR 2P2F
CR_BESTZLLos_2P2F_tle   = (CR_BESTZLLos_tle)
CR_ZLLosSEL_2P2F_tle    = (CR_BESTZLLos_tle + "&&" + BOTHFAIL)  # Is the CR_BESTZLLos request redundant? 

# ll, any combination of flavour/charge, for control regions only
process.bareLLCandtle = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('appendPhotons:electrons appendPhotons:electronstle'),
    cut = cms.string('deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)>0.02'), # protect against ghosts
    checkCharge = cms.bool(False)
)
process.LLCandtle = cms.EDProducer("ZCandidateFiller",
    src = cms.InputTag("bareLLCandtle"),
    sampleType = cms.int32(SAMPLE_TYPE),                     
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    bestZAmong = cms.string(BESTZ_AMONG),
    FSRMode = cms.string(FSRMODE), # "skip", "Legacy", "RunII"
    flags = cms.PSet(
        GoodLeptons = cms.string(ZLEPTONSEL),
        Z1Presel = cms.string(Z1PRESEL),
    )
)

process.ZLLCandtle = cms.EDProducer("ZZCandidateFiller",
    src = cms.InputTag("bareZLLCandtle"),
    sampleType = cms.int32(SAMPLE_TYPE),                    
    setup = cms.int32(LEPTON_SETUP),
    superMelaMass = cms.double(SUPERMELA_MASS),
    isMC = cms.bool(IsMC),
    bestCandComparator = cms.string(BESTCANDCOMPARATOR),
    bestCandAmong = cms.PSet(
      isBestCand    = cms.string("0"), #do not set SR best cand flag
      isBestCRZLL = cms.string(CR_BESTZLL_tle),
      isBestCRZLLss = cms.string(CR_BESTZLLss_tle),
      isBestCRZLLos_2P2F = cms.string(CR_BESTZLLos_2P2F_tle),
      isBestCRZLLos_3P1F = cms.string(CR_BESTZLLos_3P1F_tle)

    ),
    ZRolesByMass = cms.bool(False),  # daughter('Z1') = daughter(0)
    doKinFit = cms.bool(KINREFIT),
    flags = cms.PSet(
      SR = cms.string(SR),
      CRZLL = cms.string(CR_BESTZLLos_tle),
      CRZLLss = cms.string(CR_BASESEL),             #combine with proper isBestCRZLLss for AA ss/os CRss    
      CRZLLos_2P2F = cms.string(CR_ZLLosSEL_2P2F_tle),        
      CRZLLos_3P1F = cms.string(CR_ZLLosSEL_3P1F_tle),        
    ),
    recoProbabilities = cms.vstring(),
    muon_iso_cut = cms.double(MUISOCUT),
    electron_iso_cut = cms.double(ELEISOCUT),
)





process.ZlCandtle = cms.EDProducer("PATCandViewShallowCloneCombiner",
    decay = cms.string('ZCand appendPhotons:electronstle'),
    cut = cms.string("deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).eta, daughter(1).phi)>0.02 &&" + # Ghost suppression
                     "deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).eta, daughter(1).phi)>0.02 &&" +
                     ("( %s > 4) && " % ( Z_PLUS_LEP_MIJ.format(0))) + # mLL>4 for any pair cause no charge for tle OS pair (Giovanni's impl)
                     ("( %s > 4) && " % ( Z_PLUS_LEP_MIJ.format(1))) + # for trackless electrons we check all the combinations (to be validated)
                     "daughter(0).masterClone.userFloat('isBestZ') &&" +
                     "daughter(0).masterClone.userFloat('Z1Presel')"
                     ),
    checkCharge = cms.bool(False)
)


# Prepare lepton collections
process.Candidates_trackless = cms.Path(
#       process.muons             +
#       process.electrons         + process.cleanSoftElectrons +
       process.trackless_electrons+
       process.appendPhotons      + 
#       process.fsrPhotons        + process.boostedFsrPhotons +
#       process.appendPhotons     +
#       process.softLeptons       +
#       process.cleanJets         +
# Build 4-lepton candidates
#       process.bareZCand         + process.ZCand     +
       process.ZZCandSR           + ~process.ZZCandFilter + 
       process.bareZCandtle       + process.ZCandtle +  
       process.bareZZCandtle      + process.ZZCandtle
    )

process.CRtle = cms.Sequence(
    #    process.trackless_electrons +
       process.bareLLCandtle       + process.LLCandtle    +
       process.bareZLLCandtle      + process.ZLLCandtle   +
       process.ZlCandtle 
   )



process.CRlZtle = cms.Sequence(
#       process.bareZCand         + process.ZCand     +  
       process.ZlCandtle            
   )

