# Fragment to configure nanoAOD sync jobs
from ZZAnalysis.NanoAnalysis.tools import setConf
setConf("JOBTYPE", "nanoAOD")

from PhysicsTools.NanoAODTools.postprocessing.framework.branchselection import BranchSelection
def customizeBranchselForSync_(p) :
    p.outputbranchsel=BranchSelection(['drop *',
                                       'keep run',
                                       'keep event',
                                       'keep luminosityBlock',
    #                  'keep Flag*',
                                       'keep Electron_pt',
                                       'keep Electron_uncorrected_pt',
                                       'keep Electron_eta',
                                       'keep Electron_phi',
                                       'keep Electron_charge',
                                       'keep Electron_pdgId',
                                       'keep Electron_fsrPhotonIdx',
                                       'keep Electron_sip3d',
                                       'keep Electron_mvaHZZIso',
                                       'keep Muon_pt',
                                       'keep Muon_uncorrected_pt',
                                       'keep Muon_eta',
                                       'keep Muon_phi',
                                       'keep Muon_charge',
                                       'keep Muon_pdgId',
                                       'keep Muon_fsrPhotonIdx',
                                       'keep Muon_sip3d',
    #                                   'keep Lepton*',
    #                                   'drop Lepton_ZZ*',
                                       'keep FsrPhoton*',
                                       'drop FsrPhoton_mass',
                                       'drop FsrPhoton_genFsrIdx',
                                       'keep ZZCand*',
                                       'drop ZZCand_rapidity',
                                       'drop ZZCand_Phi*',
                                       'drop ZZCand_costheta*',
                                       'drop ZZCand_nExtra*',
                                       'keep bestCandIdx',
                                       'keep Generator_weight',                                   
                                       'keep puWeight',
                                       'keep ggH_NNLOPS_Weight', # save if filled 
                                       'keep overallEventWeight',
                                       ])

setConf("customizations", customizeBranchselForSync_, append=True) 
