#
# Define full paths for ZZ candidates with one loose electron (SR and CRs)
#

SIP_LOOSE = "userFloat('SIP')<100000"
GOODLEPTON_LOOSE = "userFloat('ID') && " + SIP_LOOSE


process.softLooseElectrons = cms.EDProducer("EleFiller",
   src    = cms.InputTag("bareSoftElectrons"),
   sampleType = cms.int32(SAMPLE_TYPE),
   setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
   cut = cms.string("userFloat('dxy')<0.5 && userFloat('dz')<1"),# removed cut because variable is in Spring15 ID  && userFloat('missingHit')<=1"),
   flags = cms.PSet(
        ID = cms.string("userFloat('isBDT')"),
        isSIP = cms.string(SIP_LOOSE),
        isGood = cms.string(GOODLEPTON_LOOSE),
        isGoodRegular = cms.string(GOODLEPTON), # the "regular" (tight) selection
        isIsoFSRUncorr  = cms.string("userFloat('combRelIsoPF')<"+ELEISOCUT),
        isLoose = cms.string("abs(1)"), #FIXME: I'd set this to  (isGood&&!isGoodTight), that would be clearer I think.
#       Note: passCombRelIsoPFFSRCorr is currently set in LeptonPhotonMatcher for new FSR strategy; in ZZCandidateFiller for the old one
        ),
   mvaValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16V1Values"), # (when running VID)
   )

### ----------------------------------------------------------------------
### Lepton Cleaning (clean electrons collection from muons)
### ----------------------------------------------------------------------

process.cleanSoftLooseElectrons = cms.EDProducer("PATElectronCleaner",
    # pat electron input source
    src = cms.InputTag("softLooseElectrons"),
    # preselection (any string-based cut for pat::Electron)
    preselection = cms.string(''),
    # overlap checking configurables
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
           src       = cms.InputTag("softMuons"), # Start from loose lepton def
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string("userFloat('isGood')"),
           deltaR              = cms.double(0.05),
           checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), # overlaps don't cause the electron to be discared
        )
    ),
    # finalCut (any string-based cut for pat::Electron)
    finalCut = cms.string(''),
)

process.loose_electrons = cms.Sequence(process.bareSoftElectrons + process.softLooseElectrons + process.cleanSoftLooseElectrons )
process.electrons += process.loose_electrons

# l+l- (SFOS, both e and mu)
# Note special cuts to allow SS pairs and deactivation of checkCharge (does not combine anything if set to True)
# SKIPPERMUTATIONS is used since the loose ele collection also include electrons passing the standard (tight) selection, when both electrons pass the standard selection we end up with 2 pairings, and can skip either one. this is not the case when d1 passes the loos selection and not the tight one.
SKIPPERMUTATIONS = "((daughter(0).pt>daughter(1).pt)||(!daughter(1).masterClone.userFloat('isGoodRegular')))" 

#FIXME: in SS "standard" pairs, the second electron (loose) does not have FSR! To handle FSR correctly one would have to build
# SS "standard" pairs from 'appendPhotons:electrons appendPhotons:electrons', and merge with the SS/OS standard+loose pairs.

process.bareZCandlooseEle = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('appendPhotons:electrons appendPhotons:looseElectrons'),
    cut = cms.string(KEEPLOOSECOMB_CUT+"&&"+SKIPPERMUTATIONS),
    checkCharge = cms.bool(False)
)

TRUE_RSE_Z = "((daughter(0).pdgId*daughter(1).pdgId == +121) || (daughter(0).masterClone.userFloat('SIP') >= 4.0 || daughter(1).masterClone.userFloat('SIP') >= 4.0))"
TIGHT_RSE_Z = "mass > 60. && mass < 120."
TRUE_TIGHT_RSE_Z = TRUE_RSE_Z + '&&' + TIGHT_RSE_Z 

process.ZCandlooseEle = cms.EDProducer("ZCandidateFiller",
    src = cms.InputTag("bareZCandlooseEle"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    bestZAmong = cms.string(BESTZ_AMONG),
    FSRMode = cms.string(FSRMODE), # "skip", "Legacy", "RunII"
    flags = cms.PSet(
        GoodLeptons = cms.string(ZLEPTONSEL),
        Z1Presel = cms.string(Z1PRESEL),
        isLoose = cms.string("abs(1)"),
        isTrueTightRSEZ = cms.string(TRUE_TIGHT_RSE_Z),
    )
)

# used to be CandViewShallowCloneCombiner
process.bareZZCandlooseEle= cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('ZCand ZCandlooseEle'),
    cut = cms.string(LLLLPRESEL),
    checkCharge = cms.bool(False)
)


process.ZZCandlooseEle = cms.EDProducer("ZZCandidateFiller",
    src = cms.InputTag("bareZZCandlooseEle"),
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

###
#
# CR Need to be implemnted properly! 
#
###
### Trees for control regions only
# what to do in CR where charge is needed?

# Z (OSSF,both e/mu) + LL (any F/C, with no ID/iso); this is the starting point for control regions

# Need to add the CR where the (tight) Z1 is made of an RSE+electron and the LL-pair is from the regular e/mu colections
Z1_IS_TRUE_TIGHT_RSE = "daughter(0).masterClone.userFloat('isTrueTightRSEZ')"

process.bareZLLCandlooseEle= cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('ZCand LLCandlooseEle'),
    cut = cms.string(NOGHOST4l),
    checkCharge = cms.bool(False)
)

process.bareZLLCandlooseEleZ1RSE = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('ZCandlooseEle LLCand'),
    cut = cms.string(NOGHOST4l + '&&' + Z1_IS_TRUE_TIGHT_RSE),
    checkCharge = cms.bool(False)
)

process.ZLLCandZ1RSE = cms.EDProducer("ZZCandidateFiller",
    src = cms.InputTag("bareZLLCandlooseEleZ1RSE"),
    sampleType = cms.int32(SAMPLE_TYPE),                    
    setup = cms.int32(LEPTON_SETUP),
    superMelaMass = cms.double(SUPERMELA_MASS),
    isMC = cms.bool(IsMC),
    bestCandComparator = cms.string(BESTCANDCOMPARATOR),
    bestCandAmong = cms.PSet(
      isBestCand    = cms.string("0"), #do not set SR best cand flag
      isBestCRZLLss = cms.string(CR_BESTZLLss),
      isBestCRZLLos_2P2F = cms.string(CR_BESTZLLos_2P2F),
      isBestCRZLLos_3P1F = cms.string(CR_BESTZLLos_3P1F)

    ),
    ZRolesByMass = cms.bool(False),  # daughter('Z1') = daughter(0)
    doKinFit = cms.bool(KINREFIT),
    flags = cms.PSet(
      SR = cms.string(SR),
      CRZLLss = cms.string(CR_BASESEL),             #combine with proper isBestCRZLLss for AA ss/os CRss    
      CRZLLos_2P2F = cms.string(CR_ZLLosSEL_2P2F), 
      CRZLLos_3P1F = cms.string(CR_ZLLosSEL_3P1F),        
      number_trackless_electrons = cms.string("abs(1)"),
    ),
    recoProbabilities = cms.vstring(),
    muon_iso_cut = cms.double(MUISOCUT),
    electron_iso_cut = cms.double(ELEISOCUT),
)




Z2LL_looseEle = "abs(daughter(1).daughter(0).pdgId * daughter(1).daughter(1).pdgId) == 242"   #Z2 = e * e w/o track
Z2LL_SS_looseEle = Z2LL_looseEle #"daughter(1).daughter(0).pdgId()==daughter(1).daughter(1).pdgId()"       #Z2 = same-sign, same-flavour
Z2LL_OS_looseEle = "abs(1)"

# Need to drop the SIP cut for teh Z2 candidate
Z2SIP_looseEle = "userFloat('d1.d0.isSIP')< 4 && userFloat('d1.d1.isSIP')"  
CR_BESTCANDBASE_AA_looseEle = ("userFloat('d0.Z1Presel') && userFloat('d0.worstEleIso') <" + ELEISOCUT +
                               "&& userFloat('d0.worstMuIso') <" + MUISOCUT + "&&" +
                                Z2SIP_looseEle) # base for AA CR: # Z1 with tight leptons passing SIP and ISO, mass cuts; SIP on Z2

if SELSETUP == "allCutsAtOncePlusSmart" :
    CR_BESTZLLss_looseEle = CR_BESTCANDBASE_AA_looseEle + "&&" + Z2LL_SS_looseEle + "&&" +CR_Z2MASS + "&&" + MLLALLCOMB + "&&" + PT20_10 + "&&" + "mass>70" + "&&" + "daughter(1).mass>12" + "&&" + SMARTMALLCOMB


# 2P2L region
CR_BESTZLLos_looseEle = (CR_BESTCANDBASE_AA_looseEle    + "&&" +  
                         CR_BASESEL            + "&&" +
                         Z2LL_OS_looseEle      + "&&" + # drop OS requirement 
                         SMARTMALLCOMB         )

CR_BESTZLL_looseEle = CR_BESTCANDBASE_AA_looseEle + "&&" +CR_Z2MASS + "&&" + MLLALLCOMB + "&&" + PT20_10 + "&&" + "mass>70" + "&&" + "daughter(1).mass>12" + "&&" + SMARTMALLCOMB


# CR 3P1F
CR_BESTZLLos_3P1F_looseEle = (CR_BESTZLLos_looseEle + "&&" + PASSD0_OR_PASSD1)                 
CR_ZLLosSEL_3P1F_looseEle  = (CR_BESTZLLos_looseEle + "&&" + PASSD0_XOR_PASSD1) # Is the CR_BESTZLLos request redundant? 


# CR 2P2F
CR_BESTZLLos_2P2F_looseEle   = (CR_BESTZLLos_looseEle)
CR_ZLLosSEL_2P2F_looseEle    = (CR_BESTZLLos_looseEle + "&&" + BOTHFAIL)  # Is the CR_BESTZLLos request redundant? 
CR_ZLLosSEL_3P1F             = (CR_BESTZLLos_looseEle + "&&" + PASSD0_XOR_PASSD1)
# ll, any combination of flavour/charge, for control regions only
process.bareLLCandlooseEle = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('appendPhotons:electrons appendPhotons:looseElectrons'),
    cut = cms.string('deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)>0.02&&'+SKIPPERMUTATIONS), # protect against ghosts && skip permuations of the same 2 electrons (See above)
    checkCharge = cms.bool(False)
)
process.LLCandlooseEle = cms.EDProducer("ZCandidateFiller",
    src = cms.InputTag("bareLLCandlooseEle"),
    sampleType = cms.int32(SAMPLE_TYPE),                     
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    bestZAmong = cms.string(BESTZ_AMONG),
    FSRMode = cms.string(FSRMODE), # "skip", "Legacy", "RunII"
    flags = cms.PSet(
        GoodLeptons = cms.string(ZLEPTONSEL),
        Z1Presel = cms.string(Z1PRESEL),
    )
)

process.ZLLCandlooseEle = cms.EDProducer("ZZCandidateFiller",
    src = cms.InputTag("bareZLLCandlooseEle"),
    sampleType = cms.int32(SAMPLE_TYPE),                    
    setup = cms.int32(LEPTON_SETUP),
    superMelaMass = cms.double(SUPERMELA_MASS),
    isMC = cms.bool(IsMC),
    bestCandComparator = cms.string(BESTCANDCOMPARATOR),
    bestCandAmong = cms.PSet(
      isBestCand    = cms.string("0"), #do not set SR best cand flag
      isBestCRZLLss = cms.string(CR_BESTZLLss),
      isBestCRZLLos_2P2F = cms.string(CR_BESTZLLos_2P2F_looseEle),
      isBestCRZLLos_3P1F = cms.string(CR_BESTZLLos_3P1F_looseEle)

    ),
    ZRolesByMass = cms.bool(False),  # daughter('Z1') = daughter(0)
    doKinFit = cms.bool(KINREFIT),
    flags = cms.PSet(
      SR = cms.string(SR),
      CRZLLss = cms.string(CR_BASESEL),             #combine with proper isBestCRZLLss for AA ss/os CRss    
      CRZLLos_2P2F = cms.string(CR_ZLLosSEL_2P2F_looseEle), 
      CRZLLos_3P1F = cms.string(CR_ZLLosSEL_3P1F_looseEle),        
      number_trackless_electrons = cms.string("abs(1)"),
    ),
    recoProbabilities = cms.vstring(),
    muon_iso_cut = cms.double(MUISOCUT),
    electron_iso_cut = cms.double(ELEISOCUT),
)





process.ZlCandlooseEle = cms.EDProducer("PATCandViewShallowCloneCombiner",
    decay = cms.string('ZCand appendPhotons:looseElectrons'),
    cut = cms.string("deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).eta, daughter(1).phi)>0.02 &&" + # Ghost suppression
                     "deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).eta, daughter(1).phi)>0.02 &&" +
                     ("( %s > 4) && " % ( Z_PLUS_LEP_MIJ.format(0))) + # mLL>4 for any pair cause no charge for looseEle OS pair (Giovanni's impl)
                     ("( %s > 4) && " % ( Z_PLUS_LEP_MIJ.format(1))) + # for trackless electrons we check all the combinations (to be validated)
                     "daughter(0).masterClone.userFloat('isBestZ') &&" +
                     "daughter(0).masterClone.userFloat('Z1Presel')"
                     ),
    checkCharge = cms.bool(False)
)



# Prepare lepton collections
process.Candidates_loose = cms.Path(
#       process.muons             +
#       process.electrons         + process.cleanSoftElectrons +
       process.appendPhotons      + 
       process.ZZCandSR           + ~process.ZZCandFilter +
#       process.fsrPhotons        + process.boostedFsrPhotons +
#       process.appendPhotons     +
#       process.softLeptons       +
#       process.cleanJets         +
# Build 4-lepton candidates
#       process.bareZCand         + process.ZCand     +
       process.bareZCandlooseEle  + process.ZCandlooseEle +  
       process.bareZZCandlooseEle + process.ZZCandlooseEle
    )

process.CRlooseEle = cms.Sequence(
    #    process.trackless_electrons +
       process.bareZCand + 
       process.bareZCandlooseEle + process.ZCandlooseEle +
       process.bareLLCandlooseEle       + process.LLCandlooseEle    +
       process.bareZLLCandlooseEle       + process.ZLLCandlooseEle   +
       process.bareZLLCandlooseEleZ1RSE + process.ZLLCandZ1RSE + 

       process.ZlCandlooseEle 
   )



process.CRZllooseEle = cms.Sequence(
#       process.bareZCand         + process.ZCand     +  
       process.ZlCandlooseEle            
   )



