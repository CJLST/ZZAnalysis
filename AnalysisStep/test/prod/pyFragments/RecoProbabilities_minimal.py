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

# Probabilities needed for the categorization (+ mJJ probabilities in associated production):
# all JEC variations for each one of the 11 sources is saved only for these MELA probabilities
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal = [
   # JVBF
   "Name:JVBF_SIG_ghv1_1_JHUGen_JECNominal Process:HSMHiggs Production:JJVBF MatrixElement:JHUGen Cluster:J1JECNominal Options:AddPAux=1 DefaultME:-1",

   # JQCD
   "Name:JQCD_SIG_ghg2_1_JHUGen_JECNominal Process:HSMHiggs Production:JQCD MatrixElement:JHUGen Cluster:J1JECNominal DefaultME:-1 Options:AddPConst=1",

   # JJVBF
   "Name:JJVBF_SIG_ghv1_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",

   # JJQCD
   "Name:JJQCD_SIG_ghg2_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:JJQCD MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",

   # Hadronic ZH
   "Name:HadZH_SIG_ghz1_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",

   # Hadronic WH
   "Name:HadWH_SIG_ghw1_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:Had_WH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",
]

AJetsProdProbabilities_SpinZero_JHUGen_JESUp = [theME.replace("JECNominal", "JESUp") for theME in AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESUp_Abs = [theME.replace("JECNominal", "JESUp_Abs") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESUp_Abs_year = [theME.replace("JECNominal", "JESUp_Abs_year") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESUp_BBEC1 = [theME.replace("JECNominal", "JESUp_BBEC1") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESUp_BBEC1_year = [theME.replace("JECNominal", "JESUp_BBEC1_year") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESUp_EC2 = [theME.replace("JECNominal", "JESUp_EC2") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESUp_EC2_year = [theME.replace("JECNominal", "JESUp_EC2_year") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESUp_FlavQCD = [theME.replace("JECNominal", "JESUp_FlavQCD") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESUp_HF = [theME.replace("JECNominal", "JESUp_HF") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESUp_HF_year = [theME.replace("JECNominal", "JESUp_HF_year") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESUp_RelBal = [theME.replace("JECNominal", "JESUp_RelBal") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESUp_RelSample_year = [theME.replace("JECNominal", "JESUp_RelSample_year") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
AJetsProdProbabilities_SpinZero_JHUGen_JESDn = [theME.replace("JECNominal", "JESDn") for theME in AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_Abs = [theME.replace("JECNominal", "JESDn_Abs") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_Abs_year = [theME.replace("JECNominal", "JESDn_Abs_year") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_BBEC1 = [theME.replace("JECNominal", "JESDn_BBEC1") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_BBEC1_year = [theME.replace("JECNominal", "JESDn_BBEC1_year") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_EC2 = [theME.replace("JECNominal", "JESDn_EC2") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_EC2_year = [theME.replace("JECNominal", "JESDn_EC2_year") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_FlavQCD = [theME.replace("JECNominal", "JESDn_FlavQCD") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_HF = [theME.replace("JECNominal", "JESDn_HF") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_HF_year = [theME.replace("JECNominal", "JESDn_HF_year") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_RelBal = [theME.replace("JECNominal", "JESDn_RelBal") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_RelSample_year = [theME.replace("JECNominal", "JESDn_RelSample_year") for theME in CAT_AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
AJetsProdProbabilities_SpinZero_JHUGen_JERUp = [theME.replace("JECNominal", "JERUp") for theME in AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
AJetsProdProbabilities_SpinZero_JHUGen_JERDn = [theME.replace("JECNominal", "JERDn") for theME in AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
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
AJetsProdDecProbabilities_MCFM_JESUp = [theME.replace("JECNominal", "JESUp") for theME in AJetsProdDecProbabilities_MCFM_JECNominal]
AJetsProdDecProbabilities_MCFM_JESDn = [theME.replace("JECNominal", "JESDn") for theME in AJetsProdDecProbabilities_MCFM_JECNominal]
AJetsProdDecProbabilities_MCFM_JERUp = [theME.replace("JECNominal", "JERUp") for theME in AJetsProdDecProbabilities_MCFM_JECNominal]
AJetsProdDecProbabilities_MCFM_JERDn = [theME.replace("JECNominal", "JERDn") for theME in AJetsProdDecProbabilities_MCFM_JECNominal]

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
PMAVJJ_SUPERDIJETMELA_JESUp = [theME.replace("JECNominal", "JESUp") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESUp_Abs = [theME.replace("JECNominal", "JESUp_Abs") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESUp_Abs_year = [theME.replace("JECNominal", "JESUp_Abs_year") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESUp_BBEC1 = [theME.replace("JECNominal", "JESUp_BBEC1") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESUp_BBEC1_year = [theME.replace("JECNominal", "JESUp_BBEC1_year") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESUp_EC2 = [theME.replace("JECNominal", "JESUp_EC2") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESUp_EC2_year = [theME.replace("JECNominal", "JESUp_EC2_year") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESUp_FlavQCD = [theME.replace("JECNominal", "JESUp_FlavQCD") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESUp_HF = [theME.replace("JECNominal", "JESUp_HF") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESUp_HF_year = [theME.replace("JECNominal", "JESUp_HF_year") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESUp_RelBal = [theME.replace("JECNominal", "JESUp_RelBal") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESUp_RelSample_year = [theME.replace("JECNominal", "JESUp_RelSample_year") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESDn = [theME.replace("JECNominal", "JESDn") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESDn_Abs = [theME.replace("JECNominal", "JESDn_Abs") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESDn_Abs_year = [theME.replace("JECNominal", "JESDn_Abs_year") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESDn_BBEC1 = [theME.replace("JECNominal", "JESDn_BBEC1") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESDn_BBEC1_year = [theME.replace("JECNominal", "JESDn_BBEC1_year") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESDn_EC2 = [theME.replace("JECNominal", "JESDn_EC2") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESDn_EC2_year = [theME.replace("JECNominal", "JESDn_EC2_year") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESDn_FlavQCD = [theME.replace("JECNominal", "JESDn_FlavQCD") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESDn_HF = [theME.replace("JECNominal", "JESDn_HF") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESDn_HF_year = [theME.replace("JECNominal", "JESDn_HF_year") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESDn_RelBal = [theME.replace("JECNominal", "JESDn_RelBal") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JESDn_RelSample_year = [theME.replace("JECNominal", "JESDn_RelSample_year") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JERUp = [theME.replace("JECNominal", "JERUp") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]
PMAVJJ_SUPERDIJETMELA_JERDn = [theME.replace("JECNominal", "JERDn") for theME in PMAVJJ_SUPERDIJETMELA_JECNominal]


# Construct the final list
theRecoProbabilities = []
theRecoProbabilities.extend(DecayProbabilities_SpinZero_JHUGen)
theRecoProbabilities.extend(AJetsProdProbabilities_SpinZero_JHUGen_JECNominal)
theRecoProbabilities.extend(AJetsProdProbabilities_SpinZero_JHUGen_JESUp)
theRecoProbabilities.extend(AJetsProdProbabilities_SpinZero_JHUGen_JESUp_Abs)
theRecoProbabilities.extend(CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESUp_Abs_year)
theRecoProbabilities.extend(CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESUp_BBEC1)
theRecoProbabilities.extend(CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESUp_BBEC1_year)
theRecoProbabilities.extend(CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESUp_EC2)
theRecoProbabilities.extend(CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESUp_EC2_year)
theRecoProbabilities.extend(CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESUp_FlavQCD)
theRecoProbabilities.extend(CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESUp_HF)
theRecoProbabilities.extend(CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESUp_HF_year)
theRecoProbabilities.extend(CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESUp_RelBal)
theRecoProbabilities.extend(CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESUp_RelSample_year)
theRecoProbabilities.extend(AJetsProdProbabilities_SpinZero_JHUGen_JESDn)
theRecoProbabilities.extend(CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_Abs)
theRecoProbabilities.extend(CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_Abs_year)
theRecoProbabilities.extend(CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_BBEC1)
theRecoProbabilities.extend(CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_BBEC1_year)
theRecoProbabilities.extend(CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_EC2)
theRecoProbabilities.extend(CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_EC2_year)
theRecoProbabilities.extend(CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_FlavQCD)
theRecoProbabilities.extend(CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_HF)
theRecoProbabilities.extend(CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_HF_year)
theRecoProbabilities.extend(CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_RelBal)
theRecoProbabilities.extend(CAT_AJetsProdProbabilities_SpinZero_JHUGen_JESDn_RelSample_year)
theRecoProbabilities.extend(AJetsProdProbabilities_SpinZero_JHUGen_JERUp)
theRecoProbabilities.extend(AJetsProdProbabilities_SpinZero_JHUGen_JERDn)
theRecoProbabilities.extend(ALepsProdProbabilities_SpinZero_JHUGen)
theRecoProbabilities.extend(Probabilities_SpinOne_JHUGen)
theRecoProbabilities.extend(Probabilities_SpinTwo_JHUGen)
theRecoProbabilities.extend(DecayProbabilities_MCFM)
theRecoProbabilities.extend(AJetsProdDecProbabilities_MCFM_JECNominal)
theRecoProbabilities.extend(AJetsProdDecProbabilities_MCFM_JESUp)
theRecoProbabilities.extend(AJetsProdDecProbabilities_MCFM_JESDn)
theRecoProbabilities.extend(AJetsProdDecProbabilities_MCFM_JERUp)
theRecoProbabilities.extend(AJetsProdDecProbabilities_MCFM_JERDn)
theRecoProbabilities.extend(PM4L_SUPERMELA)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JECNominal)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESUp)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESUp_Abs)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESUp_Abs_year)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESUp_BBEC1)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESUp_BBEC1_year)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESUp_EC2)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESUp_EC2_year)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESUp_FlavQCD)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESUp_HF)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESUp_HF_year)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESUp_RelBal)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESUp_RelSample_year)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESDn)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESDn_Abs)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESDn_Abs_year)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESDn_BBEC1)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESDn_BBEC1_year)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESDn_EC2)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESDn_EC2_year)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESDn_FlavQCD)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESDn_HF)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESDn_HF_year)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESDn_RelBal)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JESDn_RelSample_year)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JERUp)
theRecoProbabilities.extend(PMAVJJ_SUPERDIJETMELA_JERDn)

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
