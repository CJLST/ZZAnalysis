# ----------------------------------------------------------------------
#
# Configuration for data, DoubleMu stream
#
#----------------------------------------------------------------------

LEPTON_SETUP = 2011

ELECORRTYPE   = "Paper" # "None", "Moriond", or "Paper"
ELEREGRESSION = "Paper" # "None", "Moriond", "PaperNoComb", or "Paper" 
APPLYMUCORR = True


# FORMAT OF SAMPLE STRING:
#  ('SampleName', 'Account', 'Dataset', 'Filter', FilesPerJob, 'Options')
# Account: 'cmgtools' or'cmgtools_group' =  whether pattuples are in /eos/cms/store/cmst3/user/cmgtools or group/cmgtools
# Filter: input file pattern (regexp)
# FilesPerJob: number of input files per job.
# Options:
#   DoubleMu, DoubleEle, MuEG = what PD on data
#   MCAllEvents = write all events of the proper gen final state into each tree, even if they don't pass selection (useful for signals)
#   MH126 = compute discriminants for MH=126, instead of default mass (currently 125.6)

samples = [
           
########
# DATA #
########
           
    ### DoubleMu
    ('DoubleMuA', 'cmgtools','/DoubleMu/Run2011A-16Jan2012-v1/AOD/V5/PAT_CMG_V5_15_0','cmgTuple.*root', 8, "DoubleMu"),
    ('DoubleMuB', 'cmgtools','/DoubleMu/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_15_0','cmgTuple.*root', 10, "DoubleMu"),

    ### DoubleEle
    ('DoubleEleA', 'cmgtools','/DoubleElectron/Run2011A-16Jan2012-v1/AOD/V5/PAT_CMG_V5_15_0','cmgTuple.*root', 7, "DoubleEle"),
    ('DoubleEleB', 'cmgtools','/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_15_0','cmgTuple.*root', 9, "DoubleEle"),

    ### MuEG
    ('MuEGA', 'cmgtools','/MuEG/Run2011A-13Dec2012-v1/RECO/PAT_CMG_V5_15_0', 'cmgTuple.*root', 9, "MuEG"),
    ('MuEGB', 'cmgtools','/MuEG/Run2011B-13Dec2012-v1/RECO/PAT_CMG_V5_15_0', 'cmgTuple.*root', 9, "MuEG"),


##############
# Simulation #
##############

    ### ZZ
    ('ZZTo2e2mu',   'cmgtools','/ZZTo2e2mu_mll4_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/',   'cmgTuple.*root', 5,  ""),
    ('ZZTo4mu',     'cmgtools','/ZZTo4mu_mll4_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/',     'cmgTuple.*root', 5,  ""),
    ('ZZTo4e',      'cmgtools','/ZZTo4e_mll4_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/',      'cmgTuple.*root', 5,  ""),
    ('ZZTo2mu2tau', 'cmgtools','/ZZTo2mu2tau_mll4_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 20, ""),
    ('ZZTo2e2tau',  'cmgtools','/ZZTo2e2tau_mll4_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/',  'cmgTuple.*root', 20, ""),
    ('ZZTo4tau',    'cmgtools','/ZZTo4tau_mll4_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/',    'cmgTuple.*root', 40, ""),
    ('ggZZ2l2l',    'cmgtools','/GluGluToZZTo2L2L_7TeV-gg2zz-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/',  'cmgTuple.*root', 3,  ""),
    ('ggZZ4l',      'cmgtools','/GluGluToZZTo4L_7TeV-gg2zz-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/',    'cmgTuple.*root', 1,  ""),

    # Additional ZZ samples for templates, in mass window
    ### ZZ95
#    ('ZZ95-160To2e2mu', 'cmgtools','/ZZTo2e2mu_7TeV_mll8_mZZ95-160-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 6, ""),
    ('ZZ95-160To2mu2tau', 'cmgtools','/ZZTo2mu2tau_7TeV_mll8_mZZ95-160-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 6, ""),
    ('ZZ95-160To4e', 'cmgtools','/ZZTo4e_7TeV_mll8_mZZ95-160-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 6, ""),
#    ('ZZ95-160To4mu', 'cmgtools','/ZZTo4mu_7TeV_mll8_mZZ95-160-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 6, ""),
#    ('ZZ95-160To4tau', 'cmgtools','/ZZTo4tau7TeV_mll8_mZZ95-160-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 6, ""),

    # Samples with off-shell production, continuum, and interference
    ### Off-shell
    ('ggTo2l2l_H125.6', 'cmgtools', '/GluGluTo2L2Lprime_H_M-125p6_7TeV-gg2vv315-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('ggTo2l2l_Continuum', 'cmgtools', '/GluGluTo2L2Lprime_Contin_7TeV-gg2vv315-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('ggTo2l2l_ContinuumInterfH125.6', 'cmgtools', '/GluGluTo2L2Lprime_HContinInterf_M-125p6_7TeV-gg2vv315-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('ggTo4l_H125.6', 'cmgtools', '/GluGluTo4L_H_M-125p6_7TeV-gg2vv315-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('ggTo4l_Continuum', 'cmgtools', '/GluGluTo4L_Contin_7TeV-gg2vv315-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),

    ### jhuGenV2 samples
#     ('jhuGenV2PseH126','cmgtools','/Higgs0MToZZTo4L_M-126_7TeV_ext-JHUgenV2-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0','cmgTuple.*root', 2, "MH126"),
#     ('jhuGenV2ScaHH126','cmgtools','/Higgs0PHToZZTo4L_M-126_7TeV_ext-JHUgenV2-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0','cmgTuple.*root', 2, "MH126"),
#     ('jhuGenV2Vec1MH126','cmgtools','/Vector1MToZZTo4L_M-126_7TeV-JHUgenV2-PYTHIA6_Tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0','cmgTuple.*root', 2, "MH126"),
#     ('jhuGenV2Vec1PH126','cmgtools','/Vector1PToZZTo4L_M-126_7TeV-JHUgenV2-PYTHIA6_Tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0','cmgTuple.*root', 2, "MH126"),
#     ('jhuGenV2H126','cmgtools','/SMHiggsToZZTo4L_M-126_7TeV_ext-JHUgenV2-PYTHIA6_Tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0','cmgTuple.*root', 2, "MH126"),
#     ('jhuGenV2GravH126','cmgtools','/Graviton2PMToZZTo4L_M-126_7TeV_ext-JHUgenV2-PYTHIA6_Tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0','cmgTuple.*root', 2, "MH126"),
#     ('jhuGenV2qqGravH126','cmgtools','/Graviton2PMqqbarToZZTo4L_M-126_7TeV-JHUgenV2-PYTHIA6_Tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0','cmgTuple.*root', 2, "MH126"),
#     ('jhuGenV2Grav2MHH126','cmgtools','/Graviton2MHToZZTo4L_M-126_7TeV-JHUGenV2-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MH126"),
#     ('jhuGenV2Grav2PBH126','cmgtools','/Graviton2PBToZZTo4L_M-126_7TeV-JHUGenV2-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MH126"),
#     ('jhuGenV2Grav2PHH126','cmgtools','/Graviton2PHToZZTo4L_M-126_7TeV-JHUGenV2-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MH126"),

    ### powheg15-jhuGenV3-spin/CP-MH126 samples
    ('powheg15jhuGenV3-0MH126', 'cmgtools','/Higgs0MToZZTo4L_M-126_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v2/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3-0PHH126', 'cmgtools','/Higgs0PHToZZTo4L_M-126_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents,MH126"),
    
    ('jhuGenV3Vec1MH126', 'cmgtools','/Vector1MToZZTo4L_M-126_7TeV-JHUgenV3-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents,MH126"),
    ('jhuGenV3Vec1PH126', 'cmgtools','/Vector1PToZZTo4L_M-126_7TeV-JHUgenV3-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents,MH126"),
    ('jhuGenV3Grav2PMH126', 'cmgtools','/Graviton2PMToZZTo4L_M-126_7TeV-JHUgenV3-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents,MH126"),
    ('jhuGenV3qqGravH126', 'cmgtools','/Graviton2PMqqbarToZZTo4L_M-126_7TeV-JHUgenV3-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents,MH126"),
    ('jhuGenV3Grav2PHH126', 'cmgtools','/Graviton2PHToZZTo4L_M-126_7TeV-JHUgenV3-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents,MH126"),
    ('jhuGenV3Grav2MHH126', 'cmgtools','/Graviton2MHToZZTo4L_M-126_7TeV-JHUgenV3-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents,MH126"),
    ('jhuGenV3Grav2PBH126', 'cmgtools','/Graviton2PBToZZTo4L_M-126_7TeV-JHUgenV3-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents,MH126"),
    
    ('powheg15jhuGenV3-0Mf01ph0H126', 'cmgtools','/Higgs0Mf01ph0ToZZTo4L_M-126_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3-0Mf01ph90H126', 'cmgtools','/Higgs0Mf01ph90ToZZTo4L_M-126_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3-0Mf01ph180H126', 'cmgtools','/Higgs0Mf01ph180ToZZTo4L_M-126_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3-0Mf01ph270H126', 'cmgtools','/Higgs0Mf01ph270ToZZTo4L_M-126_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3-0Mf05ph0H126', 'cmgtools','/Higgs0Mf05ph0ToZZTo4L_M-126_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3-0Mf05ph180H126', 'cmgtools','/Higgs0Mf05ph180ToZZTo4L_M-126_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3-0Mf05ph270H126', 'cmgtools','/Higgs0Mf05ph270ToZZTo4L_M-126_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3-0Mf05ph90H126', 'cmgtools','/Higgs0Mf05ph90ToZZTo4L_M-126_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents,MH126"),

    ### powheg15-jhuGenV3-spin/CP-MH126.6 samples
    ('powheg15jhuGenV3-0PMH125.6', 'cmgtools', '/Higgs0PMToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV3/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('powheg15jhuGenV3-0MH125.6',  'cmgtools', '/Higgs0MToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV3/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/',  'cmgTuple.*root', 3, "MCAllEvents"),
    ('powheg15jhuGenV3-0PHH125.6', 'cmgtools', '/Higgs0PHToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV3/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),

    ('powheg15jhuGenV4-0L1H125.6', 'cmgtools', '/Higgs0L1ToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV4/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('powheg15jhuGenV4-0L1f01ph0H125.6',   'cmgtools', '/Higgs0L1f01ph0ToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV4/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('powheg15jhuGenV4-0L1f05ph0H125.6', 'cmgtools', '/Higgs0L1f05ph0ToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV4/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('powheg15jhuGenV4-0L1f05ph180H125.6', 'cmgtools', '/Higgs0L1f05ph180ToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV4/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),

    ('powheg15jhuGenV3-0Mf01ph0H125.6',  'cmgtools', '/Higgs0Mf01ph0ToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV3/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV3-0Mf01ph90H125.6', 'cmgtools', '/Higgs0Mf01ph90ToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV3/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV3-0Mf05ph0H125.6',  'cmgtools', '/Higgs0Mf05ph0ToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV3/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('powheg15jhuGenV3-0Mf05ph90H125.6', 'cmgtools', '/Higgs0Mf05ph90ToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV3/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('powheg15jhuGenV4-0Mf05ph180H125.6', 'cmgtools', '/Higgs0Mf05ph180ToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV4/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),

    ('powheg15jhuGenV3-0PHf01ph0H125.6', 'cmgtools', '/Higgs0PHf01ph0ToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV3/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('powheg15jhuGenV3-0PHf01ph90H125.6', 'cmgtools', '/Higgs0PHf01ph90ToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV3/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('powheg15jhuGenV3-0PHf05ph0H125.6', 'cmgtools', '/Higgs0PHf05ph0ToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV3/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV3-0PHf05ph90H125.6', 'cmgtools', '/Higgs0PHf05ph90ToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV3/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('powheg15jhuGenV4-0PHf05ph180H125.6', 'cmgtools', '/Higgs0PHf05ph180ToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV4/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('powheg15jhuGenV3-0PHf01ph0Mf01ph0H125.6', 'cmgtools', '/Higgs0PHf01ph0Mf01ph0ToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV3/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('powheg15jhuGenV3-0PHf01ph0Mf01ph90H125.6', 'cmgtools', '/Higgs0PHf01ph0Mf01ph90ToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV3/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('powheg15jhuGenV3-0PHf033ph0Mf033ph0H125.6', 'cmgtools', '/Higgs0PHf033ph0Mf033ph0ToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV3/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('powheg15jhuGenV3-0PHf033ph0Mf033ph90H125.6', 'cmgtools', '/Higgs0PHf033ph0Mf033ph90ToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV3/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('powheg15jhuGenV3-0PHf05ph0Mf05ph0H125.6', 'cmgtools', '/Higgs0PHf05ph0Mf05ph0ToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV3/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('powheg15jhuGenV3-0PHf05ph0Mf05ph90H125.6', 'cmgtools', '/Higgs0PHf05ph0Mf05ph90ToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV3/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV4-0PHf05ph180Mf05ph0H125.6', 'cmgtools', '/Higgs0PHf05ph180Mf05ph0ToZZTo4L_M-125p6_7TeV-powheg15-JHUgenV4/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
                          
    ('jhuGenV3-Vec1Mf05ph01Pf05ph0H125.6', 'cmgtools', '/Vector1Mf05ph01Pf05ph0ToZZTo4L_M-125p6_7TeV-JHUGenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('jhuGenV3-Vec1Mf05ph01Pf05ph90H125.6', 'cmgtools', '/Vector1Mf05ph01Pf05ph90ToZZTo4L_M-125p6_7TeV-JHUGenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),

    ('jhuGenV3Vec1PH125.6', 'cmgtools', '/Vector1PToZZTo4L_M-125p6_7TeV-JHUGenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('jhuGenV3Vec1MH125.6', 'cmgtools', '/Vector1MToZZTo4L_M-125p6_7TeV-JHUGenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    
    ('jhuGenV3Grav2PH2H125.6', 'cmgtools', '/Graviton2PH2ToZZTo4L_M-125p6_7TeV-JHUGenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('jhuGenV3Grav2PH3H125.6', 'cmgtools', '/Graviton2PH3ToZZTo4L_M-125p6_7TeV-JHUGenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('jhuGenV3Grav2PH6H125.6', 'cmgtools', '/Graviton2PH6ToZZTo4L_M-125p6_7TeV-JHUGenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('jhuGenV3Grav2PH7H125.6', 'cmgtools', '/Graviton2PH7ToZZTo4L_M-125p6_7TeV-JHUGenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('jhuGenV3Grav2MH9H125.6', 'cmgtools', '/Graviton2MH9ToZZTo4L_M-125p6_7TeV-JHUGenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('jhuGenV3Grav2MH10H125.6', 'cmgtools', '/Graviton2MH10ToZZTo4L_M-125p6_7TeV-JHUGenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),

   
    ### MCFM67 samples
    ('ggTo4e_SMHContinInterf-MCFM67_H125.6', 'cmgtools', '/GluGluTo4e_SMHContinInterf_M-125p6_7TeV-MCFM67-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('ggTo4mu_SMH-MCFM67_H125.6', 'cmgtools', '/GluGluTo4mu_SMH_M-125p6_7TeV-MCFM67-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),


    ### VBF samples
    ('VBFH1000', 'cmgtools','/VBF_ToHToZZTo4L_M-1000_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 1, "MCAllEvents"),
    ('VBFH115', 'cmgtools','/VBF_ToHToZZTo4L_M-115_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('VBFH120', 'cmgtools','/VBF_ToHToZZTo4L_M-120_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('VBFH125', 'cmgtools','/VBF_ToHToZZTo4L_M-125_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('VBFH130', 'cmgtools','/VBF_ToHToZZTo4L_M-130_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('VBFH140', 'cmgtools','/VBF_ToHToZZTo4L_M-140_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('VBFH150', 'cmgtools','/VBF_ToHToZZTo4L_M-150_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('VBFH160', 'cmgtools','/VBF_ToHToZZTo4L_M-160_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('VBFH170', 'cmgtools','/VBF_ToHToZZTo4L_M-170_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('VBFH180', 'cmgtools','/VBF_ToHToZZTo4L_M-180_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('VBFH190', 'cmgtools','/VBF_ToHToZZTo4L_M-190_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('VBFH200', 'cmgtools','/VBF_ToHToZZTo4L_M-200_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('VBFH210', 'cmgtools','/VBF_ToHToZZTo4L_M-210_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('VBFH220', 'cmgtools','/VBF_ToHToZZTo4L_M-220_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 1, "MCAllEvents"),
    ('VBFH230', 'cmgtools','/VBF_ToHToZZTo4L_M-230_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 1, "MCAllEvents"),
    ('VBFH250', 'cmgtools','/VBF_ToHToZZTo4L_M-250_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 1, "MCAllEvents"),
    ('VBFH275', 'cmgtools','/VBF_ToHToZZTo4L_M-275_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 1, "MCAllEvents"),
    ('VBFH300', 'cmgtools','/VBF_ToHToZZTo4L_M-300_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 1, "MCAllEvents"),
    ('VBFH325', 'cmgtools','/VBF_ToHToZZTo4L_M-325_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 1, "MCAllEvents"),
    ('VBFH350', 'cmgtools','/VBF_ToHToZZTo4L_M-350_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 1, "MCAllEvents"),
    ('VBFH375', 'cmgtools','/VBF_ToHToZZTo4L_M-375_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 1, "MCAllEvents"),
    ('VBFH400', 'cmgtools','/VBF_ToHToZZTo4L_M-400_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 1, "MCAllEvents"),
    ('VBFH425', 'cmgtools','/VBF_ToHToZZTo4L_M-425_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 1, "MCAllEvents"),
    ('VBFH450', 'cmgtools','/VBF_ToHToZZTo4L_M-450_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 1, "MCAllEvents"),
    ('VBFH475', 'cmgtools','/VBF_ToHToZZTo4L_M-475_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 1, "MCAllEvents"),
    ('VBFH500', 'cmgtools','/VBF_ToHToZZTo4L_M-500_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 1, "MCAllEvents"),
    ('VBFH575', 'cmgtools','/VBF_ToHToZZTo4L_M-575_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 1, "MCAllEvents"),
    ('VBFH600', 'cmgtools','/VBF_ToHToZZTo4L_M-600_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 1, "MCAllEvents"),
    ('VBFH650', 'cmgtools','/VBF_ToHToZZTo4L_M-650_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 1, "MCAllEvents"),
    ('VBFH700', 'cmgtools','/VBF_ToHToZZTo4L_M-700_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 1, "MCAllEvents"),
    ('VBFH800', 'cmgtools','/VBF_ToHToZZTo4L_M-800_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 1, "MCAllEvents"),
    ('VBFH900', 'cmgtools','/VBF_ToHToZZTo4L_M-900_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 1, "MCAllEvents"),
    ('VBFH950', 'cmgtools','/VBF_ToHToZZTo4L_M-950_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 1, "MCAllEvents"),

    ### VBF-powheg15
    ('powheg15VBFH200', 'cmgtools','/VBFToHToZZTo4L_M-200_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH225', 'cmgtools','/VBFToHToZZTo4L_M-225_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
#    ('powheg15VBFH250', 'cmgtools','/VBFToHToZZTo4L_M-250_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH275', 'cmgtools','/VBFToHToZZTo4L_M-275_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH300', 'cmgtools','/VBFToHToZZTo4L_M-300_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH350', 'cmgtools','/VBFToHToZZTo4L_M-350_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH400', 'cmgtools','/VBFToHToZZTo4L_M-400_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH450', 'cmgtools','/VBFToHToZZTo4L_M-450_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH500', 'cmgtools','/VBFToHToZZTo4L_M-500_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH550', 'cmgtools','/VBFToHToZZTo4L_M-550_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH600', 'cmgtools','/VBFToHToZZTo4L_M-600_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH650', 'cmgtools','/VBFToHToZZTo4L_M-650_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH700', 'cmgtools','/VBFToHToZZTo4L_M-700_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
#    ('powheg15VBFH750', 'cmgtools','/VBFToHToZZTo4L_M-750_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH800', 'cmgtools','/VBFToHToZZTo4L_M-800_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
#    ('powheg15VBFH850', 'cmgtools','/VBFToHToZZTo4L_M-850_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH900', 'cmgtools','/VBFToHToZZTo4L_M-900_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
#    ('powheg15VBFH950', 'cmgtools','/VBFToHToZZTo4L_M-950_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH1000', 'cmgtools','/VBFToHToZZTo4L_M-1000_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),

    # New split ZH/WH/ttH samples
    ### ttH
    ('ttH110', 'cmgtools_group','/TTbarH_HToZZTo4L_M-110_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ttH115', 'cmgtools_group','/TTbarH_HToZZTo4L_M-115_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ttH120', 'cmgtools_group','/TTbarH_HToZZTo4L_M-120_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ttH125', 'cmgtools_group','/TTbarH_HToZZTo4L_M-125_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ttH126', 'cmgtools_group','/TTbarH_HToZZTo4L_M-126_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MH126"),
    ('ttH130', 'cmgtools_group','/TTbarH_HToZZTo4L_M-130_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ttH140', 'cmgtools_group','/TTbarH_HToZZTo4L_M-140_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ttH150', 'cmgtools_group','/TTbarH_HToZZTo4L_M-150_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ttH160', 'cmgtools_group','/TTbarH_HToZZTo4L_M-160_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ttH180', 'cmgtools_group','/TTbarH_HToZZTo4L_M-180_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ttH200', 'cmgtools_group','/TTbarH_HToZZTo4L_M-200_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),

    ### WH
    ('WH110', 'cmgtools_group','/WH_HToZZTo4L_M-110_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('WH115', 'cmgtools_group','/WH_HToZZTo4L_M-115_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('WH120', 'cmgtools_group','/WH_HToZZTo4L_M-120_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('WH125', 'cmgtools_group','/WH_HToZZTo4L_M-125_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('WH126', 'cmgtools_group','/WH_HToZZTo4L_M-126_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MH126"),
    ('WH130', 'cmgtools_group','/WH_HToZZTo4L_M-130_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('WH140', 'cmgtools_group','/WH_HToZZTo4L_M-140_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('WH150', 'cmgtools_group','/WH_HToZZTo4L_M-150_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('WH160', 'cmgtools_group','/WH_HToZZTo4L_M-160_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('WH180', 'cmgtools_group','/WH_HToZZTo4L_M-180_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('WH200', 'cmgtools_group','/WH_HToZZTo4L_M-200_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),

    ### ZH
    ('ZH110', 'cmgtools_group','/ZH_HToZZTo4L_M-110_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ZH115', 'cmgtools_group','/ZH_HToZZTo4L_M-115_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ZH120', 'cmgtools_group','/ZH_HToZZTo4L_M-120_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ZH125', 'cmgtools_group','/ZH_HToZZTo4L_M-125_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ZH126', 'cmgtools_group','/ZH_HToZZTo4L_M-126_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MH126"),
    ('ZH130', 'cmgtools_group','/ZH_HToZZTo4L_M-130_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ZH140', 'cmgtools_group','/ZH_HToZZTo4L_M-140_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ZH150', 'cmgtools_group','/ZH_HToZZTo4L_M-150_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ZH160', 'cmgtools_group','/ZH_HToZZTo4L_M-160_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ZH180', 'cmgtools_group','/ZH_HToZZTo4L_M-180_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ZH200', 'cmgtools','/ZH_HToZZTo4L_M-200_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, ""),

    ### powheg-GluGluH NEW
    ('powheg15H125', 'cmgtools','/GluGluToHToZZTo4L_M-125_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('powheg15H126', 'cmgtools','/GluGluToHToZZTo4L_M-126_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MH126"),
    ('powheg15H225', 'cmgtools','/GluGluToHToZZTo4L_M-225_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('powheg15H250', 'cmgtools','/GluGluToHToZZTo4L_M-250_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('powheg15H275', 'cmgtools','/GluGluToHToZZTo4L_M-275_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('powheg15H300', 'cmgtools','/GluGluToHToZZTo4L_M-300_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('powheg15H350', 'cmgtools','/GluGluToHToZZTo4L_M-350_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('powheg15H400', 'cmgtools','/GluGluToHToZZTo4L_M-400_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('powheg15H450', 'cmgtools','/GluGluToHToZZTo4L_M-450_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('powheg15H500', 'cmgtools','/GluGluToHToZZTo4L_M-500_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('powheg15H550', 'cmgtools','/GluGluToHToZZTo4L_M-550_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('powheg15H600', 'cmgtools','/GluGluToHToZZTo4L_M-600_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('powheg15H650', 'cmgtools','/GluGluToHToZZTo4L_M-650_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('powheg15H700', 'cmgtools','/GluGluToHToZZTo4L_M-700_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('powheg15H800', 'cmgtools','/GluGluToHToZZTo4L_M-800_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('powheg15H900', 'cmgtools','/GluGluToHToZZTo4L_M-900_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('powheg15H1000', 'cmgtools','/GluGluToHToZZTo4L_M-1000_7TeV-powheg15-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),

    ### powheg15-jhuGenV3 samples
    ('powheg15jhuGenV3H115', 'cmgtools','/SMHiggsToZZTo4L_M-115_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, ""),
    ('powheg15jhuGenV3H120', 'cmgtools','/SMHiggsToZZTo4L_M-120_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, ""),
    ('powheg15jhuGenV3H122', 'cmgtools','/SMHiggsToZZTo4L_M-122_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, ""),
    ('powheg15jhuGenV3H124', 'cmgtools','/SMHiggsToZZTo4L_M-124_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, ""),
    ('powheg15jhuGenV3H125', 'cmgtools','/SMHiggsToZZTo4L_M-125_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, ""),
    ('powheg15jhuGenV3H126', 'cmgtools','/SMHiggsToZZTo4L_M-126_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v2/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MH126"),
    ('powheg15jhuGenV3H128', 'cmgtools','/SMHiggsToZZTo4L_M-128_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, ""),
    ('powheg15jhuGenV3H130', 'cmgtools','/SMHiggsToZZTo4L_M-130_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, ""),
    ('powheg15jhuGenV3H135', 'cmgtools','/SMHiggsToZZTo4L_M-135_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, ""),
    ('powheg15jhuGenV3H140', 'cmgtools','/SMHiggsToZZTo4L_M-140_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, ""),
    ('powheg15jhuGenV3H145', 'cmgtools','/SMHiggsToZZTo4L_M-145_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, ""),
    ('powheg15jhuGenV3H150', 'cmgtools','/SMHiggsToZZTo4L_M-150_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, ""),
    ('powheg15jhuGenV3H160', 'cmgtools','/SMHiggsToZZTo4L_M-160_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, ""),
    ('powheg15jhuGenV3H170', 'cmgtools','/SMHiggsToZZTo4L_M-170_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, ""),
    ('powheg15jhuGenV3H175', 'cmgtools','/SMHiggsToZZTo4L_M-175_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, ""),
    ('powheg15jhuGenV3H180', 'cmgtools','/SMHiggsToZZTo4L_M-180_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, ""),
    ('powheg15jhuGenV3H185', 'cmgtools','/SMHiggsToZZTo4L_M-185_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, ""),
    ('powheg15jhuGenV3H190', 'cmgtools','/SMHiggsToZZTo4L_M-190_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, ""),
    ('powheg15jhuGenV3H200', 'cmgtools','/SMHiggsToZZTo4L_M-200_7TeV-powheg15-JHUgenV3-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, ""),


    # REDUCIBLE BG

    ### ZToLJ
    ('DYJetsToLLTuneZ2M50-NoB', 'cmgtools','/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 7, "DYJets_NoB"),
    ('DYJetsToLLTuneZ2M10-NoB', 'cmgtools','/DYJetsToLL_M-10To50_TuneZ2_7TeV-madgraph/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 7, "DYJets_NoB"), 
    
    ### ZToHF
    ('DYJetsToLLTuneZ2M50-B', 'cmgtools','/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 8, "DYJets_B"),
    ('DYJetsToLLTuneZ2M10-B', 'cmgtools','/DYJetsToLL_M-10To50_TuneZ2_7TeV-madgraph/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 8, "DYJets_B"),

    ### ttbar
    ('TTTo2L2Nu2B', 'cmgtools','/TTTo2L2Nu2B_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 6, ""),
    
    ### Other-bkgs - for dedicated studies
#    ('ZGMM','cmgtools_group','/ZGToMuMuG_TuneZ2_7TeV-madgraph/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0','cmgTuple.*root', 8, ""),
#    ('WZ','cmgtools_group','/WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 8, ""),
#    ('WW','cmgtools_group','/WWJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0','cmgTuple.*root', 8, ""),

    # Still running or TBD
#     #('Wjets', 'cmgtools','/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM/V5_B/PAT_CMG_V5_6_0_B','cmgTuple.*root', 8, ""),
#     #('ZGEE', 'cmgtools','/ZGToEEG_TuneZ2_7TeV-madgraph/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0','cmgTuple.*root', 8, ""),
#     #('ZGTT', 'cmgtools','/ZGToTauTauG_TuneZ2_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM/PAT_CMG_V5_15_0','cmgTuple.*root', 8, ""),
           ]

# Load deafult job config
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"
execfile(PyFilePath + "analyzer.py")        

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
