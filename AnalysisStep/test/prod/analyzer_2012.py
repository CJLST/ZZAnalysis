### ----------------------------------------------------------------------
###
### Configuration for data, DoubleMu stream
###
###----------------------------------------------------------------------

LEPTON_SETUP = 2012

ELECORRTYPE   = "Paper" # "None", "Moriond", or "Paper"
ELEREGRESSION = "Paper" # "None", "Moriond", "PaperNoComb", or "Paper" 
APPLYMUCORR = True

# Definition of samples
# PAT_CMG_V5_7_0 (regression variables implented with unofficial recipe, ptEle > 5, ptMu > 3 )
# PAT_CMG_V5_10_0 (regression variables implemented with standard recipe, ptEle > 0, ptMu > 0)


# NOTE: Sample names (first field) should not contain underscores!

samples = [

    ### DoubleMu
    ('DoubleMuA', 'cmgtools', '/DoubleMu/Run2012A-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root',        18, "DoubleMu"),
    ('DoubleMuB', 'cmgtools', '/DoubleMuParked/Run2012B-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root',  18, "DoubleMu"),
    ('DoubleMuC', 'cmgtools', '/DoubleMuParked/Run2012C-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root',  18, "DoubleMu"),
    ('DoubleMuD', 'cmgtools', '/DoubleMuParked/Run2012D-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root',  15, "DoubleMu"),


    ### DoubleEle
    ('DoubleEleA', 'cmgtools', '/DoubleElectron/Run2012A-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root', 18, "DoubleEle"),
    ('DoubleEleB', 'cmgtools', '/DoubleElectron/Run2012B-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root', 18, "DoubleEle"),
    ('DoubleEleC', 'cmgtools', '/DoubleElectron/Run2012C-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root', 18, "DoubleEle"),
    ('DoubleEleD', 'cmgtools', '/DoubleElectron/Run2012D-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root', 15, "DoubleEle"),

    
     ### MuEG
    ('MuEGA', 'cmgtools', '/MuEG/Run2012A-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root',                18, "MuEG"),
    ('MuEGB', 'cmgtools', '/MuEG/Run2012B-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root',                18, "MuEG"),
    ('MuEGC', 'cmgtools', '/MuEG/Run2012C-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root',                18, "MuEG"),
    ('MuEGD', 'cmgtools', '/MuEG/Run2012D-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root',                15, "MuEG"),
    
    

    ################ Generic MC

    ### ZZ
    ('ZZ4mu','cmgtools','/ZZTo4mu_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 4, ""),
    ('ZZ4e','cmgtools','/ZZTo4e_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 6, ""),
    ('ZZ2mu2tau','cmgtools','/ZZTo2mu2tau_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 20, ""),
    ('ZZ2e2tau','cmgtools','/ZZTo2e2tau_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 20, ""),
    ('ZZ2e2mu','cmgtools','/ZZTo2e2mu_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 12, ""),
    ('ZZ4tau','cmgtools','/ZZTo4tau_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 40, ""),
    ('ZZ4mu_ext','cmgtools','/ZZTo4mu_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 4, ""),
    ('ZZ4e_ext','cmgtools','/ZZTo4e_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 6, ""),
    ('ZZ2mu2tau_ext','cmgtools','/ZZTo2mu2tau_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 20, ""),
    ('ZZ2e2tau_ext','cmgtools','/ZZTo2e2tau_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 20, ""),
    ('ZZ2e2mu_ext','cmgtools','/ZZTo2e2mu_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 12, ""),
    ('ZZ4tau_ext','cmgtools','/ZZTo4tau_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 40, ""),
    ('ggZZ4l','cmgtools','/GluGluToZZTo4L_8TeV-gg2zz-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 3, ""),
    ('ggZZ2l2l','cmgtools','/GluGluToZZTo2L2L_TuneZ2star_8TeV-gg2zz-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 8, ""),

    ### Additional ZZ samples for templates, in mass window
    ('ZZ95-160To2e2mu','cmgtools','/ZZTo2e2mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 6, ""),
    ('ZZ95-160To2mu2tau','cmgtools','/ZZTo2mu2tau_8TeV_mll8_mZZ95-160-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 20, ""),
    ('ZZ95-160To4e','cmgtools','/ZZTo4e_8TeV_mll8_mZZ95-160-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 6, ""),
    ('ZZ95-160To4mu','cmgtools','/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 4, ""),
    ('ZZ95-160To4tau','cmgtools','/ZZTo4tau_8TeV_mll8_mZZ95-160-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 20, ""),

    ### Samples with off-shell production, continuum, and interference
    ('ggTo2l2l_H125.6',               'cmgtools','/GluGluTo2L2Lprime_H_M-125p6_8TeV-gg2vv315-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('ggTo2l2l_Continuum',            'cmgtools','/GluGluTo2L2Lprime_Contin_8TeV-gg2vv315-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('ggTo2l2l_ContinuumInterfH126.5','cmgtools','/GluGluTo2L2Lprime_HContinInterf_M-125p6_8TeV-gg2vv315-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 3, "MCAllEvents"),

    ('ggTo4l_H125.6',                 'cmgtools','/GluGluTo4L_H_M-125p6_8TeV-gg2vv315-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 3, ""),
    ('ggTo4l_Continuum',              'cmgtools','/GluGluTo4L_Contin_8TeV-gg2vv315-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 3, ""),

    ### jhuGenv2 samples
#     ('jhuGenV2H126','cmgtools',       '/SMHiggsToZZTo4L_M-126_8TeV_ext-JHUgenV2-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root',     2, "MH126"),
#     ('jhuGenV2PseH126','cmgtools',    '/Higgs0MToZZTo4L_M-126_8TeV_ext-JHUgenV2-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root',     2, "MH126"),
#     ('jhuGenV2ScaHH126','cmgtools',   '/Higgs0PHToZZTo4L_M-126_8TeV-JHUgenV2-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root',        2, "MH126"),
#     ('jhuGenV2Vec1MH126','cmgtools',  '/Vector1MToZZTo4L_M-126_8TeV-JHUgenV2-PYTHIA6_Tauola/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 4, "MH126"),
#     ('jhuGenV2Vec1PH126','cmgtools',  '/Vector1PToZZTo4L_M-126_8TeV-JHUgenV2-PYTHIA6_Tauola/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 4, "MH126"),
#     ('jhuGenV2GravH126','cmgtools',   '/Graviton2PMToZZTo4L_M-126_8TeV_ext-JHUgenV2-PYTHIA6_Tauola/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root',  2, "MH126"),
#     ('jhuGenV2qqGravH126','cmgtools', '/Graviton2PMqqbarToZZTo4L_M-126_8TeV-JHUgenV2-PYTHIA6_Tauola/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 3, "MH126"),
#     ('jhuGenV2Grav2MHH126','cmgtools','/Graviton2MHToZZTo4L_M-126_8TeV-JHUGenV2-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root',     2, "MH126"),
#     ('jhuGenV2Grav2PBH126','cmgtools','/Graviton2PBToZZTo4L_M-126_8TeV-JHUGenV2-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root',     2, "MH126"),
#     ('jhuGenV2Grav2PHH126','cmgtools','/Graviton2PHToZZTo4L_M-126_8TeV-JHUGenV2-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root',     2, "MH126"),

    ### powheg15-jhuGenv3 spin/CP samples
    ('powheg15jhuGenV3PseHH126','cmgtools','/Higgs0MToZZTo4L_M-126_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v2/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root',                3, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3ScaHH126','cmgtools','/Higgs0PHToZZTo4L_M-126_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root',               3, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3Vec1MH126','cmgtools','/Vector1MToZZTo4L_M-126_8TeV-JHUgenV3-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root',                2, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3Vec1PH126','cmgtools','/Vector1PToZZTo4L_M-126_8TeV-JHUgenV3-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root',                2, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3GravH126','cmgtools','/Graviton2PMToZZTo4L_M-126_8TeV-JHUgenV3-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root',              2, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3qqGravH126','cmgtools','/Graviton2PMqqbarToZZTo4L_M-126_8TeV-JHUgenV3-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root',       2, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3Grav2PBH126','cmgtools','/Graviton2PBToZZTo4L_M-126_8TeV-JHUGenV3-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root',           2, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3Grav2PHH126','cmgtools','/Graviton2PHToZZTo4L_M-126_8TeV-JHUGenV3-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root',           2, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3Grav2MHH126','cmgtools','/Graviton2MHToZZTo4L_M-126_8TeV-JHUGenV3-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root',           2, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3-0Mf05ph0H126','cmgtools','/Higgs0Mf05ph0ToZZTo4L_M-126_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root',     3, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3-0Mf01ph0H126','cmgtools','/Higgs0Mf01ph0ToZZTo4L_M-126_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root',     3, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3-0Mf01ph90H126','cmgtools','/Higgs0Mf01ph90ToZZTo4L_M-126_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root',   3, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3-0Mf01ph180H126','cmgtools','/Higgs0Mf01ph180ToZZTo4L_M-126_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 3, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3-0Mf01ph270H126','cmgtools','/Higgs0Mf01ph270ToZZTo4L_M-126_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 3, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3-0Mf05ph90H126','cmgtools','/Higgs0Mf05ph90ToZZTo4L_M-126_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root',   3, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3-0Mf05ph180H126','cmgtools','/Higgs0Mf05ph180ToZZTo4L_M-126_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 3, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3-0Mf05ph270H126','cmgtools','/Higgs0Mf05ph270ToZZTo4L_M-126_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 3, "MCAllEvents,MH126"),

    ### New spin/cp samples at 125.6
    ('powheg15jhuGenV3-0PHf05ph90H125.6','cmgtools','/Higgs0PHf05ph90ToZZTo4L_M-125p6_8TeV-powheg15-JHUgenV3/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('powheg15jhuGenV3-0Mf01ph90H125.6', 'cmgtools','/Higgs0Mf01ph90ToZZTo4L_M-125p6_8TeV-powheg15-JHUgenV3/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 3, "MCAllEvents"),
    ('powheg15jhuGenV3-0Mf05ph90H125.6', 'cmgtools','/Higgs0Mf05ph90ToZZTo4L_M-125p6_8TeV-powheg15-JHUgenV3/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 3, "MCAllEvents"),


    ### VBF samples
    ('VBFH116','cmgtools','/VBF_HToZZTo4L_M-116_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH117','cmgtools','/VBF_HToZZTo4L_M-117_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH118','cmgtools','/VBF_HToZZTo4L_M-118_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH119','cmgtools','/VBF_HToZZTo4L_M-119_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH120','cmgtools','/VBF_HToZZTo4L_M-120_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH121','cmgtools','/VBF_HToZZTo4L_M-121_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH122','cmgtools','/VBF_HToZZTo4L_M-122_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH123','cmgtools','/VBF_HToZZTo4L_M-123_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH124','cmgtools','/VBF_HToZZTo4L_M-124_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH125','cmgtools','/VBF_HToZZTo4L_M-125_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH126','cmgtools','/VBF_HToZZTo4L_M-126_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents,MH126"),
    ('VBFH127','cmgtools','/VBF_HToZZTo4L_M-127_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH128','cmgtools','/VBF_HToZZTo4L_M-128_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH129','cmgtools','/VBF_HToZZTo4L_M-129_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH130','cmgtools','/VBF_HToZZTo4L_M-130_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH135','cmgtools','/VBF_HToZZTo4L_M-135_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH140','cmgtools','/VBF_HToZZTo4L_M-140_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH145','cmgtools','/VBF_HToZZTo4L_M-145_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH150','cmgtools','/VBF_HToZZTo4L_M-150_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH160','cmgtools','/VBF_HToZZTo4L_M-160_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH170','cmgtools','/VBF_HToZZTo4L_M-170_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH180','cmgtools','/VBF_HToZZTo4L_M-180_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH190','cmgtools','/VBF_HToZZTo4L_M-190_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH200','cmgtools','/VBF_HToZZTo4L_M-200_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH220','cmgtools','/VBF_HToZZTo4L_M-220_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH250','cmgtools','/VBF_HToZZTo4L_M-250_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH275','cmgtools','/VBF_HToZZTo4L_M-275_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH300','cmgtools','/VBF_HToZZTo4L_M-300_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH325','cmgtools','/VBF_HToZZTo4L_M-325_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH350','cmgtools','/VBF_HToZZTo4L_M-350_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH375','cmgtools','/VBF_HToZZTo4L_M-375_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH400','cmgtools','/VBF_HToZZTo4L_M-400_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH425','cmgtools','/VBF_HToZZTo4L_M-425_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH450','cmgtools','/VBF_HToZZTo4L_M-450_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH475','cmgtools','/VBF_HToZZTo4L_M-475_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH500','cmgtools','/VBF_HToZZTo4L_M-500_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH525','cmgtools','/VBF_HToZZTo4L_M-525_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH550','cmgtools','/VBF_HToZZTo4L_M-550_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH575','cmgtools','/VBF_HToZZTo4L_M-575_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH600','cmgtools','/VBF_HToZZTo4L_M-600_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH650','cmgtools','/VBF_HToZZTo4L_M-650_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH700','cmgtools','/VBF_HToZZTo4L_M-700_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH750','cmgtools','/VBF_HToZZTo4L_M-750_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH800','cmgtools','/VBF_HToZZTo4L_M-800_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH850','cmgtools','/VBF_HToZZTo4L_M-850_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH900','cmgtools','/VBF_HToZZTo4L_M-900_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH950','cmgtools','/VBF_HToZZTo4L_M-950_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH1000','cmgtools','/VBF_HToZZTo4L_M-1000_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),

    ### VBF powheg15
    ('powheg15VBFH1000','cmgtools','/VBFToHToZZTo4L_M-1000_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH200','cmgtools','/VBFToHToZZTo4L_M-200_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH225','cmgtools','/VBFToHToZZTo4L_M-225_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH250','cmgtools','/VBFToHToZZTo4L_M-250_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH275','cmgtools','/VBFToHToZZTo4L_M-275_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH300','cmgtools','/VBFToHToZZTo4L_M-300_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH350','cmgtools','/VBFToHToZZTo4L_M-350_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH400','cmgtools','/VBFToHToZZTo4L_M-400_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH450','cmgtools','/VBFToHToZZTo4L_M-450_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH500','cmgtools','/VBFToHToZZTo4L_M-500_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH550','cmgtools','/VBFToHToZZTo4L_M-550_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH600','cmgtools','/VBFToHToZZTo4L_M-600_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH650','cmgtools','/VBFToHToZZTo4L_M-650_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH700','cmgtools','/VBFToHToZZTo4L_M-700_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH750','cmgtools','/VBFToHToZZTo4L_M-750_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH800','cmgtools','/VBFToHToZZTo4L_M-800_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH850','cmgtools','/VBFToHToZZTo4L_M-850_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH900','cmgtools','/VBFToHToZZTo4L_M-900_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15VBFH950','cmgtools','/VBFToHToZZTo4L_M-950_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),


    ### NEW split ZH/WH/ttH samples
    ('ttH110','cmgtools_group','/TTbarH_HToZZTo4L_M-110_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, ""),
    ('ttH115','cmgtools_group','/TTbarH_HToZZTo4L_M-115_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, ""),
    ('ttH120','cmgtools_group','/TTbarH_HToZZTo4L_M-120_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, ""),
    ('ttH125','cmgtools_group','/TTbarH_HToZZTo4L_M-125_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, ""),
    ('ttH126','cmgtools_group','/TTbarH_HToZZTo4L_M-126_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, "MH126"),
    ('ttH130','cmgtools_group','/TTbarH_HToZZTo4L_M-130_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, ""),
    ('ttH140','cmgtools_group','/TTbarH_HToZZTo4L_M-140_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, ""),
    ('ttH150','cmgtools_group','/TTbarH_HToZZTo4L_M-150_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, ""),
    ('ttH160','cmgtools_group','/TTbarH_HToZZTo4L_M-160_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, ""),
    ('ttH180','cmgtools_group','/TTbarH_HToZZTo4L_M-180_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, ""),
    ('ttH200','cmgtools_group','/TTbarH_HToZZTo4L_M-200_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, ""),

    ('WH110','cmgtools_group','/WH_HToZZTo4L_M-110_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('WH115','cmgtools_group','/WH_HToZZTo4L_M-115_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('WH120','cmgtools_group','/WH_HToZZTo4L_M-120_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('WH125','cmgtools_group','/WH_HToZZTo4L_M-125_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('WH126','cmgtools_group','/WH_HToZZTo4L_M-126_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MH126"),
    ('WH130','cmgtools_group','/WH_HToZZTo4L_M-130_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('WH140','cmgtools_group','/WH_HToZZTo4L_M-140_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('WH150','cmgtools_group','/WH_HToZZTo4L_M-150_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('WH160','cmgtools_group','/WH_HToZZTo4L_M-160_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('WH180','cmgtools_group','/WH_HToZZTo4L_M-180_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('WH200','cmgtools_group','/WH_HToZZTo4L_M-200_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),

    ('ZH110','cmgtools_group','/ZH_HToZZTo4L_M-110_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ZH115','cmgtools_group','/ZH_HToZZTo4L_M-115_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ZH120','cmgtools_group','/ZH_HToZZTo4L_M-120_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ZH125','cmgtools_group','/ZH_HToZZTo4L_M-125_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ZH126','cmgtools_group','/ZH_HToZZTo4L_M-126_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MH126"),
    ('ZH130','cmgtools_group','/ZH_HToZZTo4L_M-130_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ZH140','cmgtools_group','/ZH_HToZZTo4L_M-140_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ZH150','cmgtools_group','/ZH_HToZZTo4L_M-150_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ZH160','cmgtools_group','/ZH_HToZZTo4L_M-160_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ZH180','cmgtools_group','/ZH_HToZZTo4L_M-180_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ZH200','cmgtools_group','/ZH_HToZZTo4L_M-200_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),

    ### New powheg GluGluH
    ('powheg15H125','cmgtools','/GluGluToHToZZTo4L_M-125_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15H126','cmgtools','/GluGluToHToZZTo4L_M-126_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents,MH126"),
    ('powheg15H225','cmgtools','/GluGluToHToZZTo4L_M-225_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15H250','cmgtools','/GluGluToHToZZTo4L_M-250_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15H275','cmgtools','/GluGluToHToZZTo4L_M-275_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15H300','cmgtools','/GluGluToHToZZTo4L_M-300_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15H350','cmgtools','/GluGluToHToZZTo4L_M-350_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15H400','cmgtools','/GluGluToHToZZTo4L_M-400_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15H450','cmgtools','/GluGluToHToZZTo4L_M-450_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15H500','cmgtools','/GluGluToHToZZTo4L_M-500_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15H550','cmgtools','/GluGluToHToZZTo4L_M-550_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15H600','cmgtools','/GluGluToHToZZTo4L_M-600_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15H650','cmgtools','/GluGluToHToZZTo4L_M-650_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15H700','cmgtools','/GluGluToHToZZTo4L_M-700_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15H800','cmgtools','/GluGluToHToZZTo4L_M-800_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15H900','cmgtools','/GluGluToHToZZTo4L_M-900_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15H1000','cmgtools','/GluGluToHToZZTo4L_M-1000_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),


    ### powheg15-jhuGenv3 samples
    ('powheg15jhuGenV3H115','cmgtools','/SMHiggsToZZTo4L_M-115_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV3H120','cmgtools','/SMHiggsToZZTo4L_M-120_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV3H122','cmgtools','/SMHiggsToZZTo4L_M-122_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV3H124','cmgtools','/SMHiggsToZZTo4L_M-124_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV3H125','cmgtools','/SMHiggsToZZTo4L_M-125_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV3H126','cmgtools','/SMHiggsToZZTo4L_M-126_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v2/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents,MH126"),
    ('powheg15jhuGenV3H128','cmgtools','/SMHiggsToZZTo4L_M-128_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV3H130','cmgtools','/SMHiggsToZZTo4L_M-130_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV3H135','cmgtools','/SMHiggsToZZTo4L_M-135_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV3H140','cmgtools','/SMHiggsToZZTo4L_M-140_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV3H145','cmgtools','/SMHiggsToZZTo4L_M-145_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV3H150','cmgtools','/SMHiggsToZZTo4L_M-150_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV3H160','cmgtools','/SMHiggsToZZTo4L_M-160_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV3H170','cmgtools','/SMHiggsToZZTo4L_M-170_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV3H175','cmgtools','/SMHiggsToZZTo4L_M-175_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV3H180','cmgtools','/SMHiggsToZZTo4L_M-180_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV3H185','cmgtools','/SMHiggsToZZTo4L_M-185_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV3H190','cmgtools','/SMHiggsToZZTo4L_M-190_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV3H200','cmgtools','/SMHiggsToZZTo4L_M-200_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),

    # MINLO samples
    ('minloH90', 'cmgtools', '/GluGluToHToZZTo4L_M-90_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/',   'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH95', 'cmgtools', '/GluGluToHToZZTo4L_M-95_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/',   'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH100', 'cmgtools', '/GluGluToHToZZTo4L_M-100_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH105', 'cmgtools', '/GluGluToHToZZTo4L_M-105_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH110', 'cmgtools', '/GluGluToHToZZTo4L_M-110_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH115', 'cmgtools', '/GluGluToHToZZTo4L_M-115_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH120', 'cmgtools', '/GluGluToHToZZTo4L_M-120_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH124', 'cmgtools', '/GluGluToHToZZTo4L_M-124_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH125', 'cmgtools', '/GluGluToHToZZTo4L_M-125_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH126', 'cmgtools', '/GluGluToHToZZTo4L_M-126_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents,MH126"),
    ('minloH130', 'cmgtools', '/GluGluToHToZZTo4L_M-130_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH135', 'cmgtools', '/GluGluToHToZZTo4L_M-135_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH140', 'cmgtools', '/GluGluToHToZZTo4L_M-140_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH145', 'cmgtools', '/GluGluToHToZZTo4L_M-145_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH150', 'cmgtools', '/GluGluToHToZZTo4L_M-150_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH155', 'cmgtools', '/GluGluToHToZZTo4L_M-155_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH160', 'cmgtools', '/GluGluToHToZZTo4L_M-160_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH170', 'cmgtools', '/GluGluToHToZZTo4L_M-170_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH180', 'cmgtools', '/GluGluToHToZZTo4L_M-180_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH190', 'cmgtools', '/GluGluToHToZZTo4L_M-190_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH200', 'cmgtools', '/GluGluToHToZZTo4L_M-200_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH250', 'cmgtools', '/GluGluToHToZZTo4L_M-250_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH300', 'cmgtools', '/GluGluToHToZZTo4L_M-300_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH350', 'cmgtools', '/GluGluToHToZZTo4L_M-350_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH400', 'cmgtools', '/GluGluToHToZZTo4L_M-400_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH450', 'cmgtools', '/GluGluToHToZZTo4L_M-450_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH500', 'cmgtools', '/GluGluToHToZZTo4L_M-500_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH550', 'cmgtools', '/GluGluToHToZZTo4L_M-550_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH600', 'cmgtools', '/GluGluToHToZZTo4L_M-600_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH650', 'cmgtools', '/GluGluToHToZZTo4L_M-650_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH700', 'cmgtools', '/GluGluToHToZZTo4L_M-700_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH750', 'cmgtools', '/GluGluToHToZZTo4L_M-750_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH800', 'cmgtools', '/GluGluToHToZZTo4L_M-800_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH850', 'cmgtools', '/GluGluToHToZZTo4L_M-850_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH900', 'cmgtools', '/GluGluToHToZZTo4L_M-900_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH950', 'cmgtools', '/GluGluToHToZZTo4L_M-950_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),
    ('minloH1000', 'cmgtools', '/GluGluToHToZZTo4L_M-1000_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents"),


                          
    ### REDUCIBLE BG

    ### Z->LJ
    ('DYJetsToLLTuneZ2M50-NoB','cmgtools','/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0','cmgTuple.*root', 20, "DYJets_NoB"),
    ('DYJetsToLLTuneZ2M10-NoB','cmgtools','/DYJetsToLL_M-10To50filter_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 20, "DYJets_NoB"),
    ### Z->HF
    ('DYJetsToLLTuneZ2M50-B','cmgtools','/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0','cmgTuple.*root', 15, "DYJets_B"),
    ('DYJetsToLLTuneZ2M10-B','cmgtools','/DYJetsToLL_M-10To50filter_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 20, "DYJets_B"),

    ### ttbar
    ('TTTo2L2Nu2B','cmgtools','/TTTo2L2Nu2B_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 6, ""), 
        

    ### Other bkgs, for dedicated studies
    ('ZZJetsTo4L','cmgtools','/ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 6, ""), #Is this still used by anybody?
    
#      ('WWWJets','cmgtools_group','/WWWJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#      ('WWZJets','cmgtools_group','/WWZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#      ('WZZJets','cmgtools_group','/WZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""), # error in MCHistoryTool
#      ('ZZZJets','cmgtools_group','/ZZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#      ('WWGJets','cmgtools_group','/WWGJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#      ('TTWJets','cmgtools_group','/TTWJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#      ('TTZJets','cmgtools_group','/TTZJets_8TeV-madgraph_v2/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""), # error in MCHistoryTool
#      ('TTWWJets','cmgtools_group','/TTWWJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#      ('TTGJets','cmgtools_group','/TTGJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#      ('WWW_aMCatNLO','cmgtools_group','/WWW_8TeV-aMCatNLO-herwig/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#      ('TTbarW_aMCatNLO','cmgtools_group','/TTbarW_8TeV-aMCatNLO-herwig/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#      ('WZZ_aMCatNLO','cmgtools_group','/WZZ_8TeV-aMCatNLO-herwig/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),  # error in MCHistoryTool

#      ('WZ','cmgtools_group','/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0','cmgTuple.*root', 8, ""),
#      ('WWJets','cmgtools_group','/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0','cmgTuple.*root', 8, ""),
#      ('WGToLNuG','cmgtools_group','/WGToLNuG_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0','cmgTuple.*root', 8, ""),
                 
    ]

# Load deafult job config
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"
execfile(PyFilePath + "analyzer.py")        


# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


### ----------------------------------------------------------------------
###
### Change Parameters
###
###----------------------------------------------------------------------
