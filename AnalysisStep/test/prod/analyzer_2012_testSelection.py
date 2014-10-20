#----------------------------------------------------------------------
#
# Configuration for data, DoubleMu stream
#
#----------------------------------------------------------------------

SELSETUP = "conf2" # "Legacy", "conf1", "conf2", "conf3", "conf4", "conf5", "conf6"

LEPTON_SETUP = 2012

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
    ('DoubleMuA', 'cmgtools', '/DoubleMu/Run2012A-22Jan2013-v1/AOD/PAT_CMG_V5_15_0',       'cmgTuple.*root', 18, "DoubleMu"),
    ('DoubleMuB', 'cmgtools', '/DoubleMuParked/Run2012B-22Jan2013-v1/AOD/PAT_CMG_V5_15_0', 'cmgTuple.*root', 18, "DoubleMu"),
    ('DoubleMuC', 'cmgtools', '/DoubleMuParked/Run2012C-22Jan2013-v1/AOD/PAT_CMG_V5_15_0', 'cmgTuple.*root', 18, "DoubleMu"),
    ('DoubleMuD', 'cmgtools', '/DoubleMuParked/Run2012D-22Jan2013-v1/AOD/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "DoubleMu"),

    ### DoubleEle
    ('DoubleEleA', 'cmgtools', '/DoubleElectron/Run2012A-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root', 18, "DoubleEle"),
    ('DoubleEleB', 'cmgtools', '/DoubleElectron/Run2012B-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root', 18, "DoubleEle"),
    ('DoubleEleC', 'cmgtools', '/DoubleElectron/Run2012C-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root', 18, "DoubleEle"),
    ('DoubleEleD', 'cmgtools', '/DoubleElectron/Run2012D-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root', 15, "DoubleEle"),

     ### MuEG
    ('MuEGA', 'cmgtools', '/MuEG/Run2012A-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root', 18, "MuEG"),
    ('MuEGB', 'cmgtools', '/MuEG/Run2012B-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root', 18, "MuEG"),
    ('MuEGC', 'cmgtools', '/MuEG/Run2012C-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root', 18, "MuEG"),
    ('MuEGD', 'cmgtools', '/MuEG/Run2012D-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root', 15, "MuEG"),


##############
# Simulation #
##############

    ### ZZ
    ('ZZ4mu',         'cmgtools', '/ZZTo4mu_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',                    'cmgTuple.*root', 4,  ""),
    ('ZZ4e',          'cmgtools', '/ZZTo4e_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',                     'cmgTuple.*root', 6,  ""),
    ('ZZ2mu2tau',     'cmgtools', '/ZZTo2mu2tau_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',                'cmgTuple.*root', 20, ""),
    ('ZZ2e2tau',      'cmgtools', '/ZZTo2e2tau_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',                 'cmgTuple.*root', 20, ""),
    ('ZZ2e2mu',       'cmgtools', '/ZZTo2e2mu_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',                  'cmgTuple.*root', 12, ""),
    ('ZZ4tau',        'cmgtools', '/ZZTo4tau_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',                   'cmgTuple.*root', 40, ""),
#     ('ZZ4mu_ext',     'cmgtools', '/ZZTo4mu_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',                'cmgTuple.*root', 4,  ""),
#     ('ZZ4e_ext',      'cmgtools', '/ZZTo4e_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',                 'cmgTuple.*root', 6,  ""),
#     ('ZZ2mu2tau_ext', 'cmgtools', '/ZZTo2mu2tau_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',            'cmgTuple.*root', 20, ""),
##     ('ZZ2e2tau_ext',  'cmgtools', '/ZZTo2e2tau_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',             'cmgTuple.*root', 20, ""), #BUGGY sample, it was actually generated as 4tau.
#     ('ZZ2e2mu_ext',   'cmgtools', '/ZZTo2e2mu_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',              'cmgTuple.*root', 12, ""),
#     ('ZZ4tau_ext',    'cmgtools', '/ZZTo4tau_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',               'cmgTuple.*root', 40, ""),
#     ('ggZZ4l',        'cmgtools', '/GluGluToZZTo4L_8TeV-gg2zz-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',              'cmgTuple.*root', 3,  ""),
#     ('ggZZ2l2l',      'cmgtools', '/GluGluToZZTo2L2L_TuneZ2star_8TeV-gg2zz-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 8,  ""),

    ### VBF samples
#     ('VBFH116', 'cmgtools','/VBF_HToZZTo4L_M-116_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH117', 'cmgtools','/VBF_HToZZTo4L_M-117_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH118', 'cmgtools','/VBF_HToZZTo4L_M-118_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH119', 'cmgtools','/VBF_HToZZTo4L_M-119_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH120', 'cmgtools','/VBF_HToZZTo4L_M-120_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH121', 'cmgtools','/VBF_HToZZTo4L_M-121_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH122', 'cmgtools','/VBF_HToZZTo4L_M-122_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH123', 'cmgtools','/VBF_HToZZTo4L_M-123_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH124', 'cmgtools','/VBF_HToZZTo4L_M-124_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH125', 'cmgtools','/VBF_HToZZTo4L_M-125_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
    ('VBFH126', 'cmgtools','/VBF_HToZZTo4L_M-126_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents,MH126"),
#     ('VBFH127', 'cmgtools','/VBF_HToZZTo4L_M-127_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH128', 'cmgtools','/VBF_HToZZTo4L_M-128_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH129', 'cmgtools','/VBF_HToZZTo4L_M-129_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH130', 'cmgtools','/VBF_HToZZTo4L_M-130_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH135', 'cmgtools','/VBF_HToZZTo4L_M-135_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH140', 'cmgtools','/VBF_HToZZTo4L_M-140_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH145', 'cmgtools','/VBF_HToZZTo4L_M-145_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH150', 'cmgtools','/VBF_HToZZTo4L_M-150_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH160', 'cmgtools','/VBF_HToZZTo4L_M-160_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH170', 'cmgtools','/VBF_HToZZTo4L_M-170_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH180', 'cmgtools','/VBF_HToZZTo4L_M-180_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH190', 'cmgtools','/VBF_HToZZTo4L_M-190_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH200', 'cmgtools','/VBF_HToZZTo4L_M-200_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH220', 'cmgtools','/VBF_HToZZTo4L_M-220_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH250', 'cmgtools','/VBF_HToZZTo4L_M-250_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH275', 'cmgtools','/VBF_HToZZTo4L_M-275_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH300', 'cmgtools','/VBF_HToZZTo4L_M-300_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH325', 'cmgtools','/VBF_HToZZTo4L_M-325_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH350', 'cmgtools','/VBF_HToZZTo4L_M-350_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH375', 'cmgtools','/VBF_HToZZTo4L_M-375_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH400', 'cmgtools','/VBF_HToZZTo4L_M-400_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH425', 'cmgtools','/VBF_HToZZTo4L_M-425_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH450', 'cmgtools','/VBF_HToZZTo4L_M-450_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH475', 'cmgtools','/VBF_HToZZTo4L_M-475_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH500', 'cmgtools','/VBF_HToZZTo4L_M-500_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH525', 'cmgtools','/VBF_HToZZTo4L_M-525_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH550', 'cmgtools','/VBF_HToZZTo4L_M-550_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH575', 'cmgtools','/VBF_HToZZTo4L_M-575_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH600', 'cmgtools','/VBF_HToZZTo4L_M-600_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH650', 'cmgtools','/VBF_HToZZTo4L_M-650_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH700', 'cmgtools','/VBF_HToZZTo4L_M-700_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH750', 'cmgtools','/VBF_HToZZTo4L_M-750_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH800', 'cmgtools','/VBF_HToZZTo4L_M-800_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH850', 'cmgtools','/VBF_HToZZTo4L_M-850_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH900', 'cmgtools','/VBF_HToZZTo4L_M-900_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH950', 'cmgtools','/VBF_HToZZTo4L_M-950_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),
#     ('VBFH1000', 'cmgtools','/VBF_HToZZTo4L_M-1000_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "MCAllEvents"),


    # New split ZH/WH/ttH samples
    ### ttH
#     ('ttH110', 'cmgtools_group','/TTbarH_HToZZTo4L_M-110_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, ""),
#     ('ttH115', 'cmgtools_group','/TTbarH_HToZZTo4L_M-115_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, ""),
#     ('ttH120', 'cmgtools_group','/TTbarH_HToZZTo4L_M-120_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, ""),
#     ('ttH125', 'cmgtools_group','/TTbarH_HToZZTo4L_M-125_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, ""),
    ('ttH126', 'cmgtools_group','/TTbarH_HToZZTo4L_M-126_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, "MH126"),
#     ('ttH130', 'cmgtools_group','/TTbarH_HToZZTo4L_M-130_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, ""),
#     ('ttH140', 'cmgtools_group','/TTbarH_HToZZTo4L_M-140_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, ""),
#     ('ttH150', 'cmgtools_group','/TTbarH_HToZZTo4L_M-150_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, ""),
#     ('ttH160', 'cmgtools_group','/TTbarH_HToZZTo4L_M-160_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, ""),
#     ('ttH180', 'cmgtools_group','/TTbarH_HToZZTo4L_M-180_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, ""),
#     ('ttH200', 'cmgtools_group','/TTbarH_HToZZTo4L_M-200_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, ""),

    ### WH
#     ('WH110', 'cmgtools_group','/WH_HToZZTo4L_M-110_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#     ('WH115', 'cmgtools_group','/WH_HToZZTo4L_M-115_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#     ('WH120', 'cmgtools_group','/WH_HToZZTo4L_M-120_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#     ('WH125', 'cmgtools_group','/WH_HToZZTo4L_M-125_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('WH126', 'cmgtools_group','/WH_HToZZTo4L_M-126_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MH126"),
#     ('WH130', 'cmgtools_group','/WH_HToZZTo4L_M-130_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#     ('WH140', 'cmgtools_group','/WH_HToZZTo4L_M-140_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#     ('WH150', 'cmgtools_group','/WH_HToZZTo4L_M-150_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#     ('WH160', 'cmgtools_group','/WH_HToZZTo4L_M-160_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#     ('WH180', 'cmgtools_group','/WH_HToZZTo4L_M-180_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#     ('WH200', 'cmgtools_group','/WH_HToZZTo4L_M-200_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),

    ### ZH
#     ('ZH110', 'cmgtools_group','/ZH_HToZZTo4L_M-110_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#     ('ZH115', 'cmgtools_group','/ZH_HToZZTo4L_M-115_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#     ('ZH120', 'cmgtools_group','/ZH_HToZZTo4L_M-120_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#     ('ZH125', 'cmgtools_group','/ZH_HToZZTo4L_M-125_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ZH126', 'cmgtools_group','/ZH_HToZZTo4L_M-126_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MH126"),
#     ('ZH130', 'cmgtools_group','/ZH_HToZZTo4L_M-130_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#     ('ZH140', 'cmgtools_group','/ZH_HToZZTo4L_M-140_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#     ('ZH150', 'cmgtools_group','/ZH_HToZZTo4L_M-150_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#     ('ZH160', 'cmgtools_group','/ZH_HToZZTo4L_M-160_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#     ('ZH180', 'cmgtools_group','/ZH_HToZZTo4L_M-180_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
#     ('ZH200', 'cmgtools_group','/ZH_HToZZTo4L_M-200_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),

    ### powheg-GluGluH NEW
#     ('powheg15H125', 'cmgtools','/GluGluToHToZZTo4L_M-125_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15H126', 'cmgtools','/GluGluToHToZZTo4L_M-126_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents,MH126"),
#     ('powheg15H225', 'cmgtools','/GluGluToHToZZTo4L_M-225_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15H250', 'cmgtools','/GluGluToHToZZTo4L_M-250_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15H275', 'cmgtools','/GluGluToHToZZTo4L_M-275_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15H300', 'cmgtools','/GluGluToHToZZTo4L_M-300_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15H350', 'cmgtools','/GluGluToHToZZTo4L_M-350_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15H400', 'cmgtools','/GluGluToHToZZTo4L_M-400_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15H450', 'cmgtools','/GluGluToHToZZTo4L_M-450_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15H500', 'cmgtools','/GluGluToHToZZTo4L_M-500_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15H550', 'cmgtools','/GluGluToHToZZTo4L_M-550_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15H600', 'cmgtools','/GluGluToHToZZTo4L_M-600_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15H650', 'cmgtools','/GluGluToHToZZTo4L_M-650_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15H700', 'cmgtools','/GluGluToHToZZTo4L_M-700_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15H800', 'cmgtools','/GluGluToHToZZTo4L_M-800_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15H900', 'cmgtools','/GluGluToHToZZTo4L_M-900_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15H1000', 'cmgtools','/GluGluToHToZZTo4L_M-1000_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "MCAllEvents"),

    ### powheg15-jhuGenV3 samples
#     ('powheg15jhuGenV3H115', 'cmgtools','/SMHiggsToZZTo4L_M-115_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15jhuGenV3H120', 'cmgtools','/SMHiggsToZZTo4L_M-120_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15jhuGenV3H122', 'cmgtools','/SMHiggsToZZTo4L_M-122_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15jhuGenV3H124', 'cmgtools','/SMHiggsToZZTo4L_M-124_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15jhuGenV3H125', 'cmgtools','/SMHiggsToZZTo4L_M-125_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
    ('powheg15jhuGenV3H126', 'cmgtools','/SMHiggsToZZTo4L_M-126_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v2/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents,MH126"),
#     ('powheg15jhuGenV3H128', 'cmgtools','/SMHiggsToZZTo4L_M-128_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15jhuGenV3H130', 'cmgtools','/SMHiggsToZZTo4L_M-130_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15jhuGenV3H135', 'cmgtools','/SMHiggsToZZTo4L_M-135_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15jhuGenV3H140', 'cmgtools','/SMHiggsToZZTo4L_M-140_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15jhuGenV3H145', 'cmgtools','/SMHiggsToZZTo4L_M-145_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15jhuGenV3H150', 'cmgtools','/SMHiggsToZZTo4L_M-150_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15jhuGenV3H160', 'cmgtools','/SMHiggsToZZTo4L_M-160_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15jhuGenV3H170', 'cmgtools','/SMHiggsToZZTo4L_M-170_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15jhuGenV3H175', 'cmgtools','/SMHiggsToZZTo4L_M-175_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15jhuGenV3H180', 'cmgtools','/SMHiggsToZZTo4L_M-180_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15jhuGenV3H185', 'cmgtools','/SMHiggsToZZTo4L_M-185_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15jhuGenV3H190', 'cmgtools','/SMHiggsToZZTo4L_M-190_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),
#     ('powheg15jhuGenV3H200', 'cmgtools','/SMHiggsToZZTo4L_M-200_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, "MCAllEvents"),

           ]

# Load deafult job config
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"
execfile(PyFilePath + "analyzer.py")        

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
