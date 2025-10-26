####
# Steering file for ZZ analysis starting from nanoAODs.
# Example for customization and running: test/runLocal.py
####

from __future__ import print_function
import os
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

from ZZAnalysis.NanoAnalysis.tools import setConf, getConf, insertBefore, insertAfter
from ZZAnalysis.NanoAnalysis.getEleBDTCut import *
from ZZAnalysis.NanoAnalysis.triggerAndSkim import * # Trigger requirements are defined here
from ZZAnalysis.NanoAnalysis.lepFiller import *
from ZZAnalysis.NanoAnalysis.jetFiller import *
from ZZAnalysis.NanoAnalysis.ZZFiller import *
from ZZAnalysis.NanoAnalysis.ZZExtraFiller import *


### Get processing customizations, if defined in the including .py; use defaults otherwise
DEBUG = getConf("DEBUG", False)
SAMPLENAME = getConf("SAMPLENAME", "test")
LEPTON_SETUP = getConf("LEPTON_SETUP", 2018)
DATA_TAG = getConf("DATA_TAG", "" ) # used to distinguish different subperiods/reprocessings.
                                    # Specific values currently recognized (other values->use defaults for era)
                                    # "UL" (used by muonScaleResProducer_Rochester, getEleBDTCut)
                                    # "ULAPV", (used by LeptonSFHelper)
                                    # "pre_EE" (used by LeptonSFHelper, eleScaleResProducer, muonScaleResProducer, puWeightProducer, jetJERC)
                                    # "2022E", "2022F", "2022G" (used by jetJERC)
NANOVERSION = getConf("NANOVERSION", 12)
#NANOVERSION = getConf("NANOVERSION", 15)
if not (LEPTON_SETUP == 2016 or LEPTON_SETUP == 2017 or LEPTON_SETUP == 2018 or LEPTON_SETUP == 2022 or LEPTON_SETUP == 2023 or LEPTON_SETUP == 2024 or LEPTON_SETUP == 2025) :
    print("Invalid LEPTON_SETUP", LEPTON_SETUP)
    exit(1)
IsMC = getConf("IsMC", True)
PD = getConf("PD", "")
XSEC = getConf("XSEC", 1.)
SYNCMODE = getConf("SYNCMODE", False) # fake smearing in Run2 correction modules, for synchronization purposes. No longer needed for Run3 modules.
runMELA = getConf("runMELA", True)
bestCandByMELA = getConf("bestCandByMELA", False) # requires also runMELA=True
TRIGPASSTHROUGH = getConf("TRIGPASSTHROUGH", False) # Do not filter events that do not pass triggers (HLT_passZZ4l records if they did)
PROCESS_CR = getConf("PROCESS_CR", False) # fill control regions
PROCESS_ZL = getConf("PROCESS_ZL", False) # fill ZL control region
APPLYMUCORR = getConf("APPLYMUCORR", True) # apply muon momentum scale/resolution corrections
APPLYELECORR = getConf("APPLYELECORR", True) # apply electron momentum scale/resolution corrections
APPLYJETCORR = getConf("APPLYJETCORR", True) # apply jet corrections
MUON_ID_BYMVA = getConf("MUON_ID_BYMVA", False) # if false - standard selection for muons ; if true - new WP (Muon_mvalowPt > -0.6, sip < 8, no iso)
# ggH NNLOPS weight
APPLY_QCD_GGF_UNCERT = getConf("APPLY_QCD_GGF_UNCERT", False)
# K factors for ggZZ (and old NLO ggH samples) 0:None; 1: NNLO/LO; 2: NNLO/NLO; 3: NLO/LO
APPLY_K_NNLOQCD_ZZGG = getConf("APPLY_K_NNLOQCD_ZZGG", 0)
# K factors for qqZZ
APPLY_K_NNLOQCD_ZZQQB = getConf("APPLY_K_NNLOQCD_ZZQQB", False)
APPLY_K_NNLOEW_ZZQQB  = getConf("APPLY_K_NNLOEW_ZZQQB", False)
# Add separate tree with gen info for all events
ADD_ALLEVENTS = getConf("ADD_ALLEVENTS", False)
FILTER_EVENTS = getConf("FILTER_EVENTS", 'Cands') # Filter to be applied on events. Currently supported:
                                                  # 'Cands' = any event with a SR or CR candidate (default)
                                                  # 'Z' = any event with a good Z candidate (passing the analysis Z selection criteria)
                                                  # '3L_20_10' = any event with  with 3 good leptons, pt1>20, pt2>10 (useful for trigger studies)
                                                  # 'NoFilter' = no additional filtering (besides trigger, PV filter)

CANDSTOSTORE = getConf("CANDSTOSTORE", 'BestCandOnly') # which candidates should be stored in the ZZCand collection:
                                                  # 'BestCandOnly' = only the best SR candidate in the event is saved (default)
                                                  # 'AllCands' = keep all SR candidates passing the full selection and analysis cuts
                                                  #   (including permutations of leptons).
                                                  # 'AllWithRelaxedMuId' = keep any SR candidate that can be made, even if leptons
                                                  #   don't pass ID cuts (useful for ID cut optimization studies).
                                                  # Note that this option does not affect the ZLLCand collection: for each
                                                  # CR that is activated, only the best candidate is stored.

# MELA Probabilities Dictionary:
MELAprobabilities = getConf("probabilities", None)

# Process customizations - list of functions that operate on process
customizations = getConf("customizations", [])

# Keep GenXS and GenBr for properly scaling samples with AC. 
genXS = getConf("GENXSEC", 1.)
genBR = getConf("GENBR", 1.)

from ZZAnalysis.NanoAnalysis.initializeMELA import * 
mela, melaSettings = initializeMELA(runMELA, LEPTON_SETUP, probabilities=MELAprobabilities)
                                                  
### Definition of analysis cuts
cuts = dict(
    ### lepton ID cuts
    muPt = 5.,
    elePt = 7.,
    #relIso = 0.35,
    relIso_ele = 1e9,
    relIso_mu = (1e9 if MUON_ID_BYMVA else 0.35),
    sip3d_ele = 4.,
    sip3d_mu = (8. if MUON_ID_BYMVA else 4.),
    dxy =  0.5,
    dz = 1.,
    fsr_dRET2 = 0.012,
    fsr_Iso = 1.8,
    muMva = -0.6,

    ## Relaxed ID without SIP (starting point for SIP-less CR)
    # Notes: Muon.nStations is numberOfMatchedStation, not numberOfMatches; also, muonBestTrackType!=2 is not available in nanoAODs
    muRelaxedIdNoSIP = (lambda l : (l.pt > cuts["muPt"]
                                    and abs(l.eta) < 2.4
                                    and abs(l.dxy) < cuts["dxy"]
                                    and abs(l.dz) < cuts["dz"]
                                    and (l.isGlobal or (l.isTracker and l.nStations>0)))),
    eleRelaxedIdNoSIP = (lambda l : (l.pt > cuts["elePt"]
                                     and abs(l.eta) < 2.5
                                     and abs(l.dxy) < cuts["dxy"]
                                     and abs(l.dz) < cuts["dz"])),

    passEleBDT = getEleBDTCut(LEPTON_SETUP, DATA_TAG, NANOVERSION, APPLYELECORR),

    passMuID = ((lambda l: ((l.isPFcand or (l.highPtId>0 and l.pt>200.))) and l.mvaLowPt > cuts["muMva"]) if MUON_ID_BYMVA else (lambda l: (l.isPFcand or (l.highPtId>0 and l.pt>200.)))),

    # Relaxed IDs used for CRs for fake rate method
    muRelaxedId  = (lambda l : cuts["muRelaxedIdNoSIP"](l) and abs(l.sip3d) < cuts["sip3d_mu"]),
    eleRelaxedId = (lambda l : cuts["eleRelaxedIdNoSIP"](l) and abs(l.sip3d) < cuts["sip3d_ele"]),

    # Full ID except for SIP (without isolation: FSR-corrected iso has to be applied on top, for muons)
    muFullIdNoSIP  = (lambda l, era : cuts["muRelaxedIdNoSIP"](l) and cuts["passMuID"](l)),
    eleFullIdNoSIP = (lambda l, era : cuts["eleRelaxedIdNoSIP"](l) and cuts["passEleBDT"](l)),

    # Full ID (without isolation: FSR-corrected iso has to be applied on top, for muons)
    muFullId  = (lambda l, era : cuts["muRelaxedId"](l) and cuts["passMuID"](l)),
    eleFullId = (lambda l, era : cuts["eleRelaxedId"](l) and cuts["passEleBDT"](l)),
    )

### Preselection to speed up processing.
if FILTER_EVENTS == 'NoFilter' :
    preselection = None
    postPresel =  lambda evt : (True)
else :
    if ADD_ALLEVENTS : # No preselection in the postprocessor; filter events in cloneBranches, which needs to see all events
        preselection = None
        if FILTER_EVENTS == 'Z' :
            postPresel = lambda evt : (evt.nMuon>=2 or evt.nElectron>=2)
        elif PROCESS_ZL or FILTER_EVENTS == '3L_20_10' :
            postPresel = lambda evt : (evt.nMuon+evt.nElectron>=3)
        else :
            postPresel = lambda evt : (evt.nMuon+evt.nElectron>=4)
    else : # Set a preselection for the postprocessor
        if FILTER_EVENTS == 'Z' :
            preselection = "(nMuon>=2 || nElectron>=2) && (Sum$(Muon_pt > {muPt}-2.)>=2 || Sum$(Electron_pt>{elePt}-2.)>=2)".format(**cuts)
        elif PROCESS_ZL or FILTER_EVENTS == '3L_20_10' :
            preselection = "nMuon+nElectron >= 3 && Sum$(Muon_pt > {muPt}-2.)+Sum$(Electron_pt>{elePt}-2.)>= 3".format(**cuts)
        else :
            preselection = "nMuon+nElectron >= 4 && Sum$(Muon_pt > {muPt}-2.)+Sum$(Electron_pt>{elePt}-2.)>= 4".format(**cuts)

### Input file specification
store = getConf("store","") # "/eos/cms/" for files available on eos; "root://cms-xrd-global.cern.ch/" for remote files
fileNames = getConf("fileNames", ["/store/mc/RunIIAutumn18NanoAODv7/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/BA6D7F40-ED5E-7D4E-AB14-CE8A9C5DE7EC.root",]) # to be set in calling .py
for i, file in enumerate(fileNames):
    fileNames[i] = store+file

localPath = os.environ['CMSSW_BASE']+"/src/ZZAnalysis/NanoAnalysis/"

### JSON
jsonFile = None
if not IsMC :
    if LEPTON_SETUP == 2018 :
        jsonFile = localPath+"test/prod/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt"
    elif LEPTON_SETUP == 2022 :
        jsonFile = localPath+"test/prod/Cert_Collisions2022_355100_362760_Golden.json"
    elif LEPTON_SETUP == 2023 :
        jsonFile = localPath+"test/prod/Cert_Collisions2023_366442_370790_Golden.json"
    elif LEPTON_SETUP == 2024 :
        jsonFile = localPath+"test/prod/Cert_Collisions2024_378981_386951_Golden.json"
    elif LEPTON_SETUP == 2025 :
        jsonFile = localPath+"test/prod/Cert_Collisions2025_391658_393446_Golden.json"
    else:        
        exit(1) #2016-17 to be implemented

### Modules to be run

# Standard sequence used for both data and MC
from ZZAnalysis.NanoAnalysis.RecoProbFiller import * 
reco_sequence = [lepFiller(cuts, LEPTON_SETUP, MUON_ID_BYMVA), # FSR and FSR-corrected iso; flags for passing IDs
                 ZZFiller(bestCandByMELA, mela,
                          isMC=IsMC,
                          year=LEPTON_SETUP,
                          data_tag=DATA_TAG,
                          processCR=PROCESS_CR,
                          addZL=PROCESS_ZL,
                          filter=FILTER_EVENTS,
                          candsToStore=CANDSTOSTORE,
                          debug=DEBUG), # Build ZZ candidates; choose best candidate; filter events with candidates
                 jetFiller(), # Jets cleaning with leptons
                 ZZExtraFiller(mela, IsMC, LEPTON_SETUP, DATA_TAG, PROCESS_CR), # Additional variables to selected candidates
                 ]

if MELAprobabilities != None:
    reco_sequence.append(RecoProbFiller(mela, NANOVERSION, melaSettings))  #Reco level probabilities. 

# Add muon scale corrections
if APPLYMUCORR :
    if LEPTON_SETUP < 2022 : # use Run2 Rochester corrections
        from ZZAnalysis.NanoAnalysis.modules.muonScaleResProducer_Rochester import muonScaleRes
        insertBefore(reco_sequence, 'lepFiller', muonScaleRes(LEPTON_SETUP, DATA_TAG, overwritePt=True, syncMode=SYNCMODE))
    else : # Run3 correction module
        from ZZAnalysis.NanoAnalysis.modules.muonScaleResProducer import getMuonScaleRes
        insertBefore(reco_sequence, 'lepFiller', getMuonScaleRes(LEPTON_SETUP, DATA_TAG, IsMC, overwritePt=True))
        
# Add ele scale corrections for Run 3
if APPLYELECORR and LEPTON_SETUP >=2022 :
    from ZZAnalysis.NanoAnalysis.modules.eleScaleResProducer import getEleScaleRes
    insertBefore(reco_sequence, 'lepFiller', getEleScaleRes(LEPTON_SETUP, DATA_TAG, IsMC, overwritePt=True))

# Update of JetId: manual recipe for NanoAOD v12, using JSON for v13 onwards, cf https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13p6TeV#nanoAOD_Flags
if NANOVERSION == 12 :
    from ZZAnalysis.NanoAnalysis.jetIdUpdate import *
    insertBefore(reco_sequence, 'jetFiller', jetIdUpdate())
elif NANOVERSION >=13 :
    from ZZAnalysis.NanoAnalysis.modules.jetIdProducer import getJetIdProducer
    insertBefore(reco_sequence, 'jetFiller', getJetIdProducer(LEPTON_SETUP, DATA_TAG))   

# Add jet corrections for Run 3
if APPLYJETCORR and LEPTON_SETUP >=2022 :
    from ZZAnalysis.NanoAnalysis.modules.jetJERC import getJetCorrected
    insertBefore(reco_sequence, 'jetFiller', getJetCorrected(LEPTON_SETUP, DATA_TAG, IsMC, overwritePt=True))
    from ZZAnalysis.NanoAnalysis.modules.jetVMAP import getJetVetoMap
    insertBefore(reco_sequence, 'jetFiller', getJetVetoMap(LEPTON_SETUP, DATA_TAG))

# Special modules to be applied before the reco_sequence, that may filter events
pre_sequence = [triggerAndSkim(isMC=IsMC, PD=PD, era=LEPTON_SETUP, passThru=TRIGPASSTHROUGH),  # Filter for good PV and trigger requirements; apply PD precedence rules for data
                ]
# Special modules to be applied after the reco_sequence (ie only for selected events)
post_sequence = []

if IsMC:
    from ZZAnalysis.NanoAnalysis.modules.puWeightProducer import *
    from ZZAnalysis.NanoAnalysis.lepDataMCWeight import *
    insertBefore(reco_sequence, 'ZZExtraFiller', lepDataMCWeight(LEPTON_SETUP, DATA_TAG, muonIdByMVA = MUON_ID_BYMVA))

    # Weights computation, to be placed in pre or post sequences based on the configuration
    from ZZAnalysis.NanoAnalysis.weightFiller import weightFiller
    weights = weightFiller(XSEC, APPLY_K_NNLOQCD_ZZGG, APPLY_K_NNLOQCD_ZZQQB, APPLY_K_NNLOEW_ZZQQB, APPLY_QCD_GGF_UNCERT, LEPTON_SETUP)
    
    #Protect against writing a bunch of 1's. 
    if (genXS != 1) and (genBR != 1): 
        from ZZAnalysis.NanoAnalysis.genXSFiller import * 
        post_sequence.append(genXSFiller(genXS,genBR))
    

    from ZZAnalysis.NanoAnalysis.mcTruthAnalyzer import *
    post_sequence.append(mcTruthAnalyzer()) # Gen final state etc.
    # from ZZAnalysis.NanoAnalysis.genExtraFiller import *
    # post_sequence.append(genExtraFiller(mela)) Gen-level angles (not to be confused with LHE-level angles, filled by LHEAngProbFiller

    if ADD_ALLEVENTS: # Add modules that produce the variables to be stored for all events at the beginning
        from ZZAnalysis.NanoAnalysis.genFiller import *
        from ZZAnalysis.NanoAnalysis.cloneBranches import *
        from ZZAnalysis.NanoAnalysis.LHEAngProbFiller import * 
        pre_sequence = [puWeight(LEPTON_SETUP, DATA_TAG),
                        weights, 
                        genFiller(mela, dump=False),
                        LHEAngProbFiller(mela, NANOVERSION, melaSettings),
                        cloneBranches(treeName='AllEvents',
                                      varlist=['run', 'luminosityBlock', 'event',
                                               'GenDressedLepton_*',
                                               'FidDressedLeps_*',
                                               'FidZ*',
                                               'LHE*Weight',
                                               'passedFiducial',
                                               'Generator_weight',
                                               'puWeight*',
                                               'ggH_NNLOPS_Weight',
                                               'overallEventWeight',
                                               'Pileup_nTrueInt',
                                               'LHEMela*',
                                               'GenJet*',
                                               ],
                                      #Stop further processing for events that don't have 4 reco leps
                                      continueFor = postPresel
                                      ),
                        ] + pre_sequence
        if NANOVERSION >= 15:
            from ZZAnalysis.NanoAnalysis.LHEFiller import * 
            insertBefore(pre_sequence, 'LHEAngProbFiller', LHEFiller())
        
        

    else : # Add them at the end, so that they are run only for selected events
        post_sequence.extend([puWeight(LEPTON_SETUP,DATA_TAG),
                              weights,
                              #genFiller(mela, dump=False), # Not required when ADD_ALLEVENTS = False?
                              ])
else : # Data
    post_sequence = []


ZZSequence = pre_sequence + reco_sequence + post_sequence

if CANDSTOSTORE == 'AllWithRelaxedMuId' : # Add extra variables for ID studies
    from ZZAnalysis.NanoAnalysis.ZZIDStudies import *
    insertAfter(ZZSequence, 'ZZFiller', ZZIDStudies())

### Branches to be read and written to output
branchsel_in = ['drop FatJet_*',
                'drop IsoTrack*',
                'drop L1_*',
                'drop Photon*',
                'drop SV_*',
                'drop SoftActivityJet_*',
                'drop SubJet*',
                'drop Tau*',]

branchsel_out = ['drop *',
                 'keep run',
                 'keep event',
                 'keep luminosityBlock',
                 'keep Flag*',
                 'keep Electron*',
                 'keep Muon*',
                 'keep Lepton*',
                 'keep Jet*',
                 'keep nCleanedJet*',
                 'keep FsrPhoton*',
                 # individual HLT bits are different in different data periods/eras and this causes some problems with merging data files
                 # 'keep HLT_Ele*',
                 # 'keep HLT_DoubleEle*',
                 # 'keep HLT_Mu*',
                 # 'keep HLT_DiMu*',
                 # 'keep HLT_TripleMu*',
                 # 'keep HLT_IsoMu*',
                 'keep HLT_passZZ*',
                 'keep best*', # best candidate indices
                 'keep Z*', # Z, ZZ, ZLL candidates
                 'keep MET_pt',
                 #'keep PV*',
                 #'keep Flag*',
                 ]

if IsMC:
    branchsel_in.extend(['drop GenIsolatedPhoton_*',
                         ])
    branchsel_out.extend(['keep GenPart*',
                          'keep GenZZ*',
                          'keep *eight', # Generator_weight + custom weights
                          'keep puWeight*',
                          'keep HTXS_Higgs*',
                          'keep HTXS_njets30',
                          'keep Pileup*',
                          'keep GenJet_*',
                          'keep genxsec', 
                          'keep genbr',
                          #'keep Generator*',
                          #'keep PV*',
                        ])

    if ADD_ALLEVENTS : # Gen-level variables that are relevant only for signals
        branchsel_out.extend(['keep GenDressedLepton_*',
                              'keep FidDressedLeps_*',
                              'keep FidZ*',
                              'keep passedFiducial',
                              'keep LHEPart*',
                              'keep LHEMela*'
                              ])

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
p = PostProcessor(".", fileNames,
                  prefetch=True, longTermCache=False,
                  cut=preselection, # pre-selection cuts (to speed up processing)
                  branchsel=branchsel_in, # select branches to be read
                  outputbranchsel=branchsel_out, # select branches to be written out
                  jsonInput=jsonFile, # path of json file for data
                  modules=ZZSequence,
                  noOut=False, # True = do not write out skimmed nanoAOD file
                  haddFileName="ZZ4lAnalysis.root", # name of output nanoAOD file
#                  histFileName="histos.root", histDirName="plots", # file containing histograms
                  maxEntries=0, # Number of events to be read
                  firstEntry=0, # First event to be read
                  provenance = False
                  )

for cf in customizations :
    print(f"Applying process customization: {cf.__name__}")
    cf(p)

# Print sequence to be run:
print("Sequence to be run:")
for mod in p.modules:
    print(" ", mod.__class__.__name__)
print ("", flush=True)

### Run command should be issued by the calling scripy
# p.run()


