import FWCore.ParameterSet.Config as cms

###############################################################################
#
# This MVA is a retraining of the Egamma Fall17 V1 MVA trained with TMVA on a
# more recent 94X # DY+Jets sample, which was requested at the preapproval of
# HIG-18-001:
#
# /DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
#     /RunIIFall17MiniAOD-RECOSIMstep_94X_mc2017_realistic_v10-v1/MINIAODSIM
#
# ...unlike the original Fall17 V1 ID, which is still trained with 92X DY+Jets:
#
# /DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
#     /RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10_ext1-v1/MINIAODSIM
#
# Documentation of the original Fall17 MVA:
# * https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariateElectronIdentificationRun2
# * https://rembserj.web.cern.ch/rembserj/notes/Electron_MVA_ID_2017_documentation
#
# Scale factor approval talk in Egamma:
# * https://rembserj.web.cern.ch/rembserj/slides/180216_HZZ_SF_egamma.pdf
#
# Only the scale factors for the mvaEleID-Fall17-iso-V1_94X-wpHZZ working point
# in the iso-inclusive ID are officially approved by Egamma.
#
# Also, don't confuse this ID with the Egamma Fall17V2 ID, which is trained
# with xgboost instead of TMVA, resulting in about 10 % less fakes for constant
# signal efficiency and same data/MC scale factors.
#
###############################################################################

# This MVA implementation class name
mvaFall17ClassName = "ElectronMVAEstimatorRun2Fall17NoIso"
# The tag is an extra string attached to the names of the products
# such as ValueMaps that needs to distinguish cases when the same MVA estimator
# class is used with different tuning/weights
mvaTag = "V1_94X"

# The parameters according to which the training bins are split:
ptSplit = 10.      # we have above and below 10 GeV categories
ebSplit = 0.800    # barrel is split into two regions
ebeeSplit = 1.479  # division between barrel and endcap

# There are 6 categories in this MVA. They have to be configured in this strict order
# (cuts and weight files order):
#   0   EB1 (eta<0.8)  pt 5-10 GeV     |   pt < ptSplit && |eta| < ebSplit
#   1   EB2 (eta>=0.8) pt 5-10 GeV     |   pt < ptSplit && |eta| >= ebSplit && |eta| < ebeeSplit
#   2   EE             pt 5-10 GeV     |   pt < ptSplit && |eta| >= ebeeSplit
#   3   EB1 (eta<0.8)  pt 10-inf GeV   |   pt >= ptSplit && |eta| < ebSplit
#   4   EB2 (eta>=0.8) pt 10-inf GeV   |   pt >= ptSplit && |eta| >= ebSplit && |eta| < ebeeSplit
#   5   EE             pt 10-inf GeV   |   pt >= ptSplit && |eta| >= ebeeSplit


mvaFall17WeightFiles_V1_94X = cms.vstring(
    "RecoEgamma/ElectronIdentification/data/Fall17_V1_94X/EB1_5_noIso.weights.xml.gz",
    "RecoEgamma/ElectronIdentification/data/Fall17_V1_94X/EB2_5_noIso.weights.xml.gz",
    "RecoEgamma/ElectronIdentification/data/Fall17_V1_94X/EE_5_noIso.weights.xml.gz",
    "RecoEgamma/ElectronIdentification/data/Fall17_V1_94X/EB1_10_noIso.weights.xml.gz",
    "RecoEgamma/ElectronIdentification/data/Fall17_V1_94X/EB2_10_noIso.weights.xml.gz",
    "RecoEgamma/ElectronIdentification/data/Fall17_V1_94X/EE_10_noIso.weights.xml.gz"
    )

# Load some common definitions for MVA machinery
from RecoEgamma.ElectronIdentification.Identification.mvaElectronID_tools \
    import (EleMVA_WP,
            configureVIDMVAEleID_V1)

# The locatoins of value maps with the actual MVA values and categories
# for all particles.
# The names for the maps are "<module name>:<MVA class name>Values"
# and "<module name>:<MVA class name>Categories"
mvaProducerModuleLabel = "electronMVAValueMapProducer"
mvaValueMapName        = mvaProducerModuleLabel + ":" + mvaFall17ClassName + mvaTag + "Values"
mvaCategoriesMapName   = mvaProducerModuleLabel + ":" + mvaFall17ClassName + mvaTag + "Categories"

## The working point for this MVA that is expected to have about 90% signal
# WP tuned to give about 90 and 80% signal efficiecny for electrons from Drell-Yan with pT > 25 GeV
# The working point for the low pt categories is just taken over from the high pt
idName90 = "mvaEleID-Fall17-noIso-V1_94X-wp90"
MVA_WP90 = EleMVA_WP(
    idName = idName90,
    mvaValueMapName = mvaValueMapName,           # map with MVA values for all particles
    mvaCategoriesMapName = mvaCategoriesMapName, # map with category index for all particles
    cutCategory0_C0 = 0.9194119065699909,  # EB1 low pt
    cutCategory0_C1 = 2.960512196567283,
    cutCategory0_C2 = 1.1224376343699878,
    cutCategory1_C0 = 0.8615716064570684,  # EB2 low pt
    cutCategory1_C1 = 2.5091268456912874,
    cutCategory1_C2 = 0.7471085602428349,
    cutCategory2_C0 = 0.7634510333747148,  # EE low pt
    cutCategory2_C1 = 1.3018890137217956,
    cutCategory2_C2 = 17.99527477409486,
    cutCategory3_C0 = 0.9580855123030142,  # EB1
    cutCategory3_C1 = 8.904973701688647,
    cutCategory3_C2 = 3.5017111808963404,
    cutCategory4_C0 = 0.9238323042526464,  # EB2
    cutCategory4_C1 = 9.167734241330061,
    cutCategory4_C2 = 3.64643837209186,
    cutCategory5_C0 = 0.8812593700808523,  # EE
    cutCategory5_C1 = 11.51356941995292,
    cutCategory5_C2 = 4.054192899471584
    )

idName80 = "mvaEleID-Fall17-noIso-V1_94X-wp80"
MVA_WP80 = EleMVA_WP(
    idName = idName80,
    mvaValueMapName = mvaValueMapName,           # map with MVA values for all particles
    mvaCategoriesMapName = mvaCategoriesMapName, # map with category index for all particles
    cutCategory0_C0 = 0.9563759935436272,  # EB1 low pt
    cutCategory0_C1 = 2.8841209996462385,
    cutCategory0_C2 = 0.5213141042497962,
    cutCategory1_C0 = 0.929606513571081,   # EB2 low pt
    cutCategory1_C1 = 2.0615806479950085,
    cutCategory1_C2 = 0.5720789453940258,
    cutCategory2_C0 = 0.9115136155483037,  # EE low pt
    cutCategory2_C1 = 1.2555335214477281,
    cutCategory2_C2 = 8.05816302290571,
    cutCategory3_C0 = 0.9806265216254968,  # EB1
    cutCategory3_C1 = 9.043301811060447,
    cutCategory3_C2 = 1.2220348631544082,
    cutCategory4_C0 = 0.9674974218183561,  # EB2
    cutCategory4_C1 = 8.536661593554728,
    cutCategory4_C2 = 1.7375413024707924,
    cutCategory5_C0 = 0.9472792428840157,  # EE
    cutCategory5_C1 = 8.551318188121513,
    cutCategory5_C2 = 3.3411396528324135
)

### WP tuned for HZZ analysis with very high efficiency (about 98%)
# The working points were found by requiring the same signal efficiencies in
# each category as for the Spring 16 HZZ ID
# (see RecoEgamma/ElectronIdentification/python/Identification/mvaElectronID_Spring16_HZZ_V1_94X_cff.py)
idNamewpLoose = "mvaEleID-Fall17-noIso-V1_94X-wpLoose"
MVA_WPLoose = EleMVA_WP(
    idName = idNamewpLoose,
    mvaValueMapName = mvaValueMapName,           # map with MVA values for all particles
    mvaCategoriesMapName = mvaCategoriesMapName, # map with category index for all particles
    cutCategory0 =  -0.13285867293779202, # EB1 low pt
    cutCategory1 =  -0.31765300958836074, # EB2 low pt
    cutCategory2 =  -0.0799205914718861 , # EE low pt
    cutCategory3 =  -0.856871961305474  , # EB1
    cutCategory4 =  -0.8107642141584835 , # EB2
    cutCategory5 =  -0.7179265933023059   # EE
    )

#
# Configure variable names and the values they are clipped to.
# They have to appear in the same order as in the weights xml file
#

#                Name  |  Lower clip value  | upper clip value
variablesInfo = [
                 ("ele_oldsigmaietaieta"              ,  None, None),
                 ("ele_oldsigmaiphiiphi"              ,  None, None),
                 ("ele_oldcircularity"                ,   -1.,   2.),
                 ("ele_oldr9"                         ,  None,   5.),
                 ("ele_scletawidth"                   ,  None, None),
                 ("ele_sclphiwidth"                   ,  None, None),
                 ("ele_oldhe"                         ,  None, None),
                 ("ele_kfhits"                        ,  None, None),
                 ("ele_kfchi2"                        ,  None,  10.),
                 ("ele_gsfchi2"                       ,  None, 200.),
                 ("ele_fbrem"                         ,   -1., None),
                 ("ele_gsfhits"                       ,  None, None),
                 ("ele_expected_inner_hits"           ,  None, None),
                 ("ele_conversionVertexFitProbability",  None, None),
                 ("ele_ep"                            ,  None,  20.),
                 ("ele_eelepout"                      ,  None,  20.),
                 ("ele_IoEmIop"                       ,  None, None),
                 ("ele_deltaetain"                    , -0.06, 0.06),
                 ("ele_deltaphiin"                    ,  -0.6,  0.6),
                 ("ele_deltaetaseed"                  ,  -0.2,  0.2),
                 ("rho"                               ,  None, None),
                 ("ele_psEoverEraw"                   ,  None, None), # EE only
                ]

varNames, clipLower, clipUpper = [list(l) for l in zip(*variablesInfo)]
for i, x in enumerate(clipLower):
    if x == None:
        clipLower[i] = -999999.0
for i, x in enumerate(clipUpper):
    if x == None:
        clipUpper[i] =  999999.0

#
# Finally, set up VID configuration for all cuts
#

# Create the PSet that will be fed to the MVA value map producer
mvaEleID_Fall17_noIso_V1_94X_producer_config = cms.PSet(
    mvaName            = cms.string(mvaFall17ClassName),
    mvaTag             = cms.string(mvaTag),
    # This MVA uses conversion info, so configure several data items on that
    beamSpot           = cms.InputTag('offlineBeamSpot'),
    conversionsAOD     = cms.InputTag('allConversions'),
    conversionsMiniAOD = cms.InputTag('reducedEgamma:reducedConversions'),
    # Category split parameters
    ptSplit            = cms.double(ptSplit),
    ebSplit            = cms.double(ebSplit),
    ebeeSplit          = cms.double(ebeeSplit),
    # Variable clipping parameters
    varNames           = cms.vstring(*varNames),
    clipLower          = cms.vdouble(*clipLower),
    clipUpper          = cms.vdouble(*clipUpper),
    #
    weightFileNames    = mvaFall17WeightFiles_V1_94X
    )
# Create the VPset's for VID cuts
mvaEleID_Fall17_V1_94X_wpLoose = configureVIDMVAEleID_V1( MVA_WPLoose )
mvaEleID_Fall17_V1_94X_wp90 = configureVIDMVAEleID_V1( MVA_WP90, cutName="GsfEleMVAExpoScalingCut")
mvaEleID_Fall17_V1_94X_wp80 = configureVIDMVAEleID_V1( MVA_WP80, cutName="GsfEleMVAExpoScalingCut")

mvaEleID_Fall17_V1_94X_wpLoose.isPOGApproved = cms.untracked.bool(False)
mvaEleID_Fall17_V1_94X_wp90.isPOGApproved = cms.untracked.bool(False)
mvaEleID_Fall17_V1_94X_wp80.isPOGApproved = cms.untracked.bool(False)
