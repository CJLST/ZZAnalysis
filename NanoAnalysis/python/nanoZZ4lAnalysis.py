####
# Steering file for ZZ analysis starting from nanoAODs.
# Example for customization and running: test/runLocal.py
####

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import *

from ZZAnalysis.NanoAnalysis.tools import setConf, getConf
from ZZAnalysis.NanoAnalysis.triggerAndSkim import * # Trigger requirements are defined here
from ZZAnalysis.NanoAnalysis.lepFiller import *
from ZZAnalysis.NanoAnalysis.jetFiller import *
from ZZAnalysis.NanoAnalysis.ZZFiller import *
from ZZAnalysis.NanoAnalysis.ZZExtraFiller import *


### Get processing customizations, if defined in the including .py; use defaults otherwise 
DEBUG = getConf("DEBUG", False)
SAMPLENAME = getConf("SAMPLENAME", "test")
LEPTON_SETUP = getConf("LEPTON_SETUP", 2018)
DATA_TAG = getConf("DATA_TAG", "" ) # flavours; at the moment only used to mark UL Run 2 samples
NANOVERSION = getConf("NANOVERSION", 11)
if not (LEPTON_SETUP == 2016 or LEPTON_SETUP == 2017 or LEPTON_SETUP == 2018 or LEPTON_SETUP == 2022 or LEPTON_SETUP == 2023) :
    print("Invalid LEPTON_SETUP", LEPTON_SETUP)
    exit(1)
IsMC = getConf("IsMC", True)
PD = getConf("PD", "")
XSEC = getConf("XSEC", 1.)
SYNCMODE = getConf("SYNCMODE", False)
runMELA = getConf("runMELA", True)
bestCandByMELA = getConf("bestCandByMELA", True) # requires also runMELA=True
TRIGPASSTHROUGH = getConf("TRIGPASSTHROUGH", True) # Do not filter events that do not pass triggers (HLT_passZZ4l records if they did)
PROCESS_CR = getConf("PROCESS_CR", False) # fill control regions
APPLYMUCORR = getConf("APPLYMUCORR", True) # apply muon momentum scale/resolution corrections
# ggH NNLOPS weight
APPLY_QCD_GGF_UNCERT = getConf("APPLY_QCD_GGF_UNCERT", False) 
# K factors for ggZZ (and old NLO ggH samples) 0:None; 1: NNLO/LO; 2: NNLO/NLO; 3: NLO/LO
APPLY_K_NNLOQCD_ZZGG = getConf("APPLY_K_NNLOQCD_ZZGG", 0) 
# K factors for qqZZ
APPLY_K_NNLOQCD_ZZQQB = getConf("APPLY_K_NNLOQCD_ZZQQB", False) 
APPLY_K_NNLOEW_ZZQQB  = getConf("APPLY_K_NNLOEW_ZZQQB", False) 


if "UL" in DATA_TAG : preUL = False # used to set the correct electron selection
else: preUL = True

### Definition of analysis cuts
cuts = dict(
    ### lepton ID cuts
    muPt = 5.,
    elePt = 7.,
    relIso = 0.35,
    sip3d = 4.,
    dxy =  0.5,
    dz = 1.,
    fsr_dRET2 = 0.012,
    fsr_Iso = 1.8,
    
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


    passEleBDT = (lambda l, era : eleBDTCut(l, era, preUL, NANOVERSION)), #actual definition in lepFiller.py

    # Relaxed IDs used for CRs for fake rate method
    muRelaxedId  = (lambda l : cuts["muRelaxedIdNoSIP"](l) and abs(l.sip3d) < cuts["sip3d"]),
    eleRelaxedId = (lambda l : cuts["eleRelaxedIdNoSIP"](l) and abs(l.sip3d) < cuts["sip3d"]),

    # Full ID except for SIP (without isolation: FSR-corrected iso has to be applied on top, for muons)
    muFullIdNoSIP  = (lambda l, era : cuts["muRelaxedIdNoSIP"](l) and (l.isPFcand or (l.highPtId>0 and l.pt>200.))),
    eleFullIdNoSIP = (lambda l, era : cuts["eleRelaxedIdNoSIP"](l) and cuts["passEleBDT"](l, era)),

    # Full ID (without isolation: FSR-corrected iso has to be applied on top, for muons)
    muFullId  = (lambda l, era : cuts["muRelaxedId"](l) and (l.isPFcand or (l.highPtId>0 and l.pt>200.))),
    eleFullId = (lambda l, era : cuts["eleRelaxedId"](l) and cuts["passEleBDT"](l, era)),
    )




### Preselection to speed up processing. Note: to be relaxed for CRs
preselection = getConf("preselection", "nMuon + nElectron >= 4 &&" +
                       "Sum$(Muon_pt > {muPt}-2.) +" + # Allow for variations due to scale calib
                       "Sum$(Electron_pt > {elePt})" +
                       ">= 4").format(**cuts)

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
        jsonFile = localPath+"test/prod/Cert_Collisions2023_366442_370355_Golden.json"
    else:        
        exit(1) #2016-17 to be implemented

### Sequence to be run
ZZSequence = [triggerAndSkim(isMC=IsMC, PD=PD, era=LEPTON_SETUP, passThru=TRIGPASSTHROUGH)] # Filter for good PV and trigger requirements; apply PD precedence rules for data

if APPLYMUCORR and LEPTON_SETUP < 2022 : # corrections not yet implemented for Run 3
    from ZZAnalysis.NanoAnalysis.muonScaleResProducer import muonScaleRes
    ZZSequence.append(muonScaleRes(LEPTON_SETUP, DATA_TAG, overwritePt=True, syncMode=SYNCMODE)) # Sets corrected muon pT and scale uncertainty


ZZSequence.extend([lepFiller(cuts, LEPTON_SETUP), # FSR and FSR-corrected iso; flags for passing IDs
                   ZZFiller(runMELA, bestCandByMELA, IsMC, LEPTON_SETUP, PROCESS_CR, debug=DEBUG), # Build ZZ candidates; choose best candidate; filter events with candidates
                   jetFiller(), # Jets cleaning with leptons
                   ZZExtraFiller('SR'),
#                  MELAFiller(), # Compute the full set of discriminants for the best candidate
                   ])

if IsMC :
    from ZZAnalysis.NanoAnalysis.mcTruthAnalyzer import *
    ZZSequence.insert(0, mcTruthAnalyzer(dump=False)) # Gen final state

    if LEPTON_SETUP < 2022 :
        from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeight_2016, puWeight_2017, puWeight_2018
        puWeight = {2016:puWeight_2016, 2017:puWeight_2017, 2018:puWeight_2018} # FIXME official weights are slightly different than the one we use (checked for 2018)
        ZZSequence.append(puWeight[LEPTON_SETUP]()) # PU reweighting
    elif LEPTON_SETUP == 2022 : #FIXME: must update for 2023 MC! 
        from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
        pufile_2022 = "%s/src/ZZAnalysis/AnalysisStep/data/PileUpWeights/pu_weights_2022.root" % os.environ['CMSSW_BASE']
        puWeight_2022 = lambda : puWeightProducer(pufile_2022,
                                                  pufile_2022,
                                                  "MC_out_of_the_box",
                                                  "Data",
                                                  verbose=False,
                                                  doSysVar=True)
        puWeight = {2022:puWeight_2022} # 2022 weights non official, still preliminary
        ZZSequence.append(puWeight[LEPTON_SETUP]()) # PU reweighting

    from ZZAnalysis.NanoAnalysis.weightFiller import weightFiller
    ZZSequence.append(weightFiller(XSEC, APPLY_K_NNLOQCD_ZZGG, APPLY_K_NNLOQCD_ZZQQB, APPLY_K_NNLOEW_ZZQQB, APPLY_QCD_GGF_UNCERT)) # total weight

### Branches to be read and written to output
branchsel_in = ""
if IsMC:
    branchsel_in  = localPath+"python/branchsel_in_MC.txt"
    branchsel_out = localPath+"python/branchsel_out_MC.txt"
else:
    branchsel_in  = localPath+"python/branchsel_in_Data.txt"
    branchsel_out = localPath+"python/branchsel_out_Data.txt"



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
                  ) 

### Run command should be issued by the calling scripy
# p.run()
