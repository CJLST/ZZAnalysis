### Spin-0 decay probabilities from JHUGen ###
DecayProbabilities_SpinZero_JHUGen = [
   "Name:GG_SIG_ghg2_1_ghz1_1_JHUGen Alias:<Name> Process:HSMHiggs Production:ZZGG MatrixElement:JHUGen",
]
## Production probabilities with >=1 jet(s) ##
AJetsProdProbabilities_SpinZero_JHUGen_JECNominal = [
   # JVBF
   "Name:JVBF_SIG_ghv1_1_JHUGen_JECNominal Process:HSMHiggs Production:JJVBF MatrixElement:JHUGen Cluster:J1JECNominal Options:AddPAux=1 DefaultME:-1",
   # JQCD
   "Name:JQCD_SIG_ghg2_1_JHUGen_JECNominal Process:HSMHiggs Production:JQCD MatrixElement:JHUGen Cluster:J1JECNominal DefaultME:-1",
   # JJVBF
   "Name:JJVBF_SIG_ghv1_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",
   # JJQCD
   "Name:JJQCD_SIG_ghg2_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:JJQCD MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",
   # Best-DJJ JJVBF and JJQCD
   "Name:JJVBF_SIG_ghv1_1_JHUGen_JECNominal_BestDJJ Copy:JJVBF_SIG_ghv1_1_JHUGen_JECNominal Options:MaxNumerator=JJVBF_SIG_ghv1_1_JHUGen_JECNominal;MaxDenominator=JJQCD_SIG_ghg2_1_JHUGen_JECNominal",
   "Name:JJQCD_SIG_ghg2_1_JHUGen_JECNominal_BestDJJ Copy:JJQCD_SIG_ghg2_1_JHUGen_JECNominal Options:MaxNumerator=JJVBF_SIG_ghv1_1_JHUGen_JECNominal;MaxDenominator=JJQCD_SIG_ghg2_1_JHUGen_JECNominal",
   # Hadronic ZH
   "Name:HadZH_SIG_ghz1_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",
   # Hadronic WH
   "Name:HadWH_SIG_ghw1_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:Had_WH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",
   # ttH: Undecayed MEs belong to J2-class clusters
   "Name:ttHUndecayed_SIG_kappa_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:ttH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",
   # bbH
   "Name:bbH_SIG_kappa_1_JHUGen_JECNominal Process:HSMHiggs Production:bbH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",
]
AJetsProdProbabilities_SpinZero_JHUGen_JECUp = [theME.replace("JECNominal", "JECUp") for theME in AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
AJetsProdProbabilities_SpinZero_JHUGen_JECDn = [theME.replace("JECNominal", "JECDn") for theME in AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
## Production probabilities with >=1 lepton(s) ##
ALepsProdProbabilities_SpinZero_JHUGen = [
   # Leptonic ZH
   "Name:LepZH_SIG_ghz1_1_JHUGen Alias:<Name> Process:HSMHiggs Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH DefaultME:-1",
   # Leptonic WH (CAUTION: All requiring the SM ME to be maximized)
   "Name:LepWH_SIG_ghw1_1_JHUGen Alias:<Name> Process:HSMHiggs Production:Lep_WH MatrixElement:JHUGen Cluster:LepWH Options:MaxNumerator=LepWH_SIG_ghw1_1_JHUGen DefaultME:-1",
]

### Spin-1 decay probabilities from JHUGen ###
Probabilities_SpinOne_JHUGen = []
### Spin-2 decay probabilities from JHUGen ###
Probabilities_SpinTwo_JHUGen = []

### Decay probabilities from MCFM ###
DecayProbabilities_MCFM = [
   "Name:GG_SIG_kappaTopBot_1_ghz1_1_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Options:AddPConst=1",
   "Name:GG_BSI_kappaTopBot_1_ghz1_1_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM", # Has pConst=pConst_sig+pConst_bkg
   "Name:GG_BSI_kappaTopBot_1_ghz1_i_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=0,1", # Has the same pConst
   "Name:GG_BKG_MCFM Process:bkgZZ Production:ZZGG MatrixElement:MCFM Options:AddPConst=1",
   "Name:QQB_BKG_MCFM Alias:<Name> Process:bkgZZ Production:ZZQQB MatrixElement:MCFM Options:AddPConst=1", # Aliased to construct Dbkgkin
   "Name:ZJJ_BKG_MCFM Process:bkgZJets Production:JJQCD MatrixElement:MCFM",
]
AJetsProdDecProbabilities_MCFM_JECNominal = []
AJetsProdDecProbabilities_MCFM_JECUp = [theME.replace("JECNominal", "JECUp") for theME in AJetsProdDecProbabilities_MCFM_JECNominal]
AJetsProdDecProbabilities_MCFM_JECDn = [theME.replace("JECNominal", "JECDn") for theME in AJetsProdDecProbabilities_MCFM_JECNominal]

### m4l probabilities from SuperMELA ###
PM4L_SUPERMELA = [
   "Name:m4l_SIG Process:HSMHiggs Production:ZZGG MatrixElement:JHUGen SuperMelaSyst:SMSyst_None isPM4L:1",
   "Name:m4l_BKG Process:bkgZZ Production:ZZGG MatrixElement:JHUGen SuperMelaSyst:SMSyst_None isPM4L:1",
]

### mJJ probabilities in associated production
PMAVJJ_SUPERDIJETMELA_JECNominal = [
   "Name:HadZH_mavjj_JECNominal Process:HSMHiggs Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 isPMaVJJ:1",
   "Name:HadWH_mavjj_JECNominal Process:HSMHiggs Production:Had_WH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 isPMaVJJ:1",
]
PMAVJJ_SUPERDIJETMELA_JECUp = [theME.replace("JECNominal", "JECUp") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JECDn = [theME.replace("JECNominal", "JECDn") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]


# Construct the final list
theRecoProbabilities = []
theRecoProbabilities.extend(DecayProbabilities_SpinZero_JHUGen)
theRecoProbabilities.extend(AJetsProdProbabilities_SpinZero_JHUGen_JECNominal)
theRecoProbabilities.extend(AJetsProdProbabilities_SpinZero_JHUGen_JECUp)
theRecoProbabilities.extend(AJetsProdProbabilities_SpinZero_JHUGen_JECDn)
theRecoProbabilities.extend(ALepsProdProbabilities_SpinZero_JHUGen)
theRecoProbabilities.extend(Probabilities_SpinOne_JHUGen)
theRecoProbabilities.extend(Probabilities_SpinTwo_JHUGen)
theRecoProbabilities.extend(DecayProbabilities_MCFM)
theRecoProbabilities.extend(AJetsProdDecProbabilities_MCFM_JECNominal)
theRecoProbabilities.extend(AJetsProdDecProbabilities_MCFM_JECUp)
theRecoProbabilities.extend(AJetsProdDecProbabilities_MCFM_JECDn)
theRecoProbabilities.extend(PM4L_SUPERMELA)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JECNominal)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JECUp)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JECDn)

# Append final list
for name in (
             "ZZCand",          "ZZTree",
             "ZLLCand",         "CRZLLTree",
             "ZZCandlooseEle",  "ZZTreelooseEle",
             "ZLLCandlooseEle", "CRZLLTreelooseEle",
             "ZLLCandZ1RSE",    "CRZLLTreeZ1RSE",
             "ZZCandtle",       "ZZTreetle",
             "ZLLCandtle",      "CRZLLTreetle",
            ):
    if hasattr(process, name):
        getattr(process, name).recoProbabilities.extend(theRecoProbabilities)
