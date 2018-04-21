### Spin-0 decay probabilities from JHUGen ###
DecayProbabilities_SpinZero_JHUGen = [
   "Name:GG_SIG_ghg2_1_ghz1_1_JHUGen Alias:<Name> Process:HSMHiggs Production:ZZGG MatrixElement:JHUGen Options:AddPConst=1",
   "Name:GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen Alias:<Name> Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1_prime2=10000,0",
   "Name:GG_SIG_ghg2_1_ghz2_1_JHUGen Alias:<Name> Process:H0hplus Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_ghg2_1_ghz4_1_JHUGen Alias:<Name> Process:H0minus Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen Alias:<Name> Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghzgs1_prime2=10000,0",
   "Name:GG_SIG_ghg2_1_ghza2_1_JHUGen Alias:<Name> Process:H0_Zgs Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_ghg2_1_ghza4_1_JHUGen Alias:<Name> Process:H0_Zgs_PS Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_ghg2_1_gha2_1_JHUGen Alias:<Name> Process:H0_gsgs Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_ghg2_1_gha4_1_JHUGen Alias:<Name> Process:H0_gsgs_PS Production:ZZGG MatrixElement:JHUGen",

   "Name:GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghz1_prime2=10000,0 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen",
   #"Name:GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4i_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghz1_prime2=0,10000 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen",
   "Name:GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghz2=1,0 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghz2_1_JHUGen",
   "Name:GG_SIG_ghg2_1_ghz1_1_ghz2_i_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghz2=0,1 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghz2_1_JHUGen",
   "Name:GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghz4=1,0 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghz4_1_JHUGen",
   "Name:GG_SIG_ghg2_1_ghz1_1_ghz4_i_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghz4=0,1 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghz4_1_JHUGen",

   "Name:GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghzgs1_prime2=10000,0 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen",
   "Name:GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4i_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghzgs1_prime2=0,10000 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen",
   "Name:GG_SIG_ghg2_1_ghz1_1_ghza2_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghzgs2=1,0 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghza2_1_JHUGen",
   #"Name:GG_SIG_ghg2_1_ghz1_1_ghza2_i_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghzgs2=0,1 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghza2_1_JHUGen",
   "Name:GG_SIG_ghg2_1_ghz1_1_ghza4_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghzgs4=1,0 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghza4_1_JHUGen",
   #"Name:GG_SIG_ghg2_1_ghz1_1_ghza4_i_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghzgs4=0,1 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghza4_1_JHUGen",

   "Name:GG_SIG_ghg2_1_ghz1_1_gha2_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghgsgs2=1,0 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_gha2_1_JHUGen",
   #"Name:GG_SIG_ghg2_1_ghz1_1_gha2_i_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghgsgs2=0,1 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_gha2_1_JHUGen",
   "Name:GG_SIG_ghg2_1_ghz1_1_gha4_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghgsgs4=1,0 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_gha4_1_JHUGen",
   #"Name:GG_SIG_ghg2_1_ghz1_1_gha4_i_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghgsgs4=0,1 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_gha4_1_JHUGen",

   "Name:GG_SIG_ghg2_1_ghz1prime2_1E4_ghza1prime2_1E4_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1_prime2=10000,0;ghzgs1_prime2=10000,0 Options:SubtractP=GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen,GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen",
   #"Name:GG_SIG_ghg2_1_ghz1prime2_1E4_ghza1prime2_1E4i_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1_prime2=10000,0;ghzgs1_prime2=0,10000 Options:SubtractP=GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen,GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen",
]
## Production probabilities with >=1 jet(s) ##
AJetsProdProbabilities_SpinZero_JHUGen_JECNominal = [
   # JVBF
   "Name:JVBF_SIG_ghv1_1_JHUGen_JECNominal Process:HSMHiggs Production:JJVBF MatrixElement:JHUGen Cluster:J1JECNominal Options:AddPAux=1 DefaultME:-1",

   # JQCD
   "Name:JQCD_SIG_ghg2_1_JHUGen_JECNominal Process:HSMHiggs Production:JQCD MatrixElement:JHUGen Cluster:J1JECNominal DefaultME:-1 Options:AddPConst=1",
   #"Name:JQCD_SIG_ghg4_1_JHUGen_JECNominal Process:H0minus Production:JQCD MatrixElement:JHUGen Cluster:J1JECNominal DefaultME:-1",

   # JJVBF
   "Name:JJVBF_SIG_ghv1_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",
   "Name:JJVBF_SIG_ghv1prime2_1E4_JHUGen_JECNominal Alias:<Name> Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1_prime2=10000,0 DefaultME:-1",
   "Name:JJVBF_SIG_ghv2_1_JHUGen_JECNominal Alias:<Name> Process:H0hplus Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",
   "Name:JJVBF_SIG_ghv4_1_JHUGen_JECNominal Alias:<Name> Process:H0minus Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",
   "Name:JJVBF_SIG_ghza1prime2_1E4_JHUGen_JECNominal Alias:<Name> Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghzgs1_prime2=10000,0 DefaultME:-1",
   "Name:JJVBF_SIG_ghza2_1_JHUGen_JECNominal Alias:<Name> Process:H0_Zgs Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",
   "Name:JJVBF_SIG_ghza4_1_JHUGen_JECNominal Alias:<Name> Process:H0_Zgs_PS Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",
   "Name:JJVBF_SIG_gha2_1_JHUGen_JECNominal Alias:<Name> Process:H0_gsgs Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",
   "Name:JJVBF_SIG_gha4_1_JHUGen_JECNominal Alias:<Name> Process:H0_gsgs_PS Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",

   "Name:JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghz1_prime2=10000,0 Options:SubtractP=JJVBF_SIG_ghv1_1_JHUGen_JECNominal,JJVBF_SIG_ghv1prime2_1E4_JHUGen_JECNominal DefaultME:-1",
   #"Name:JJVBF_SIG_ghv1_1_ghv1prime2_1E4i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghz1_prime2=0,10000 Options:SubtractP=JJVBF_SIG_ghv1_1_JHUGen_JECNominal,JJVBF_SIG_ghv1prime2_1E4_JHUGen_JECNominal DefaultME:-1",
   "Name:JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghz2=1,0 Options:SubtractP=JJVBF_SIG_ghv1_1_JHUGen_JECNominal,JJVBF_SIG_ghv2_1_JHUGen_JECNominal DefaultME:-1",
   #"Name:JJVBF_SIG_ghv1_1_ghv2_i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghz2=0,1 Options:SubtractP=JJVBF_SIG_ghv1_1_JHUGen_JECNominal,JJVBF_SIG_ghv2_1_JHUGen_JECNominal DefaultME:-1",
   "Name:JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghz4=1,0 Options:SubtractP=JJVBF_SIG_ghv1_1_JHUGen_JECNominal,JJVBF_SIG_ghv4_1_JHUGen_JECNominal DefaultME:-1",
   #"Name:JJVBF_SIG_ghv1_1_ghv4_i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghz4=0,1 Options:SubtractP=JJVBF_SIG_ghv1_1_JHUGen_JECNominal,JJVBF_SIG_ghv4_1_JHUGen_JECNominal DefaultME:-1",

   "Name:JJVBF_SIG_ghv1_1_ghza1prime2_1E4_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghzgs1_prime2=10000,0 Options:SubtractP=JJVBF_SIG_ghv1_1_JHUGen_JECNominal,JJVBF_SIG_ghza1prime2_1E4_JHUGen_JECNominal DefaultME:-1",
   #"Name:JJVBF_SIG_ghv1_1_ghza1prime2_1E4i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghzgs1_prime2=0,10000 Options:SubtractP=JJVBF_SIG_ghv1_1_JHUGen_JECNominal,JJVBF_SIG_ghza1prime2_1E4_JHUGen_JECNominal DefaultME:-1",
   "Name:JJVBF_SIG_ghv1_1_ghza2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghzgs2=1,0 Options:SubtractP=JJVBF_SIG_ghv1_1_JHUGen_JECNominal,JJVBF_SIG_ghza2_1_JHUGen_JECNominal DefaultME:-1",
   #"Name:JJVBF_SIG_ghv1_1_ghza2_i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghzgs2=0,1 Options:SubtractP=JJVBF_SIG_ghv1_1_JHUGen_JECNominal,JJVBF_SIG_ghza2_1_JHUGen_JECNominal DefaultME:-1",
   "Name:JJVBF_SIG_ghv1_1_ghza4_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghzgs4=1,0 Options:SubtractP=JJVBF_SIG_ghv1_1_JHUGen_JECNominal,JJVBF_SIG_ghza4_1_JHUGen_JECNominal DefaultME:-1",
   #"Name:JJVBF_SIG_ghv1_1_ghza4_i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghzgs4=0,1 Options:SubtractP=JJVBF_SIG_ghv1_1_JHUGen_JECNominal,JJVBF_SIG_ghza4_1_JHUGen_JECNominal DefaultME:-1",

   "Name:JJVBF_SIG_ghv1_1_gha2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghgsgs2=1,0 Options:SubtractP=JJVBF_SIG_ghv1_1_JHUGen_JECNominal,JJVBF_SIG_gha2_1_JHUGen_JECNominal DefaultME:-1",
   #"Name:JJVBF_SIG_ghv1_1_gha2_i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghgsgs2=0,1 Options:SubtractP=JJVBF_SIG_ghv1_1_JHUGen_JECNominal,JJVBF_SIG_gha2_1_JHUGen_JECNominal DefaultME:-1",
   "Name:JJVBF_SIG_ghv1_1_gha4_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghgsgs4=1,0 Options:SubtractP=JJVBF_SIG_ghv1_1_JHUGen_JECNominal,JJVBF_SIG_gha4_1_JHUGen_JECNominal DefaultME:-1",
   #"Name:JJVBF_SIG_ghv1_1_gha4_i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghgsgs4=0,1 Options:SubtractP=JJVBF_SIG_ghv1_1_JHUGen_JECNominal,JJVBF_SIG_gha4_1_JHUGen_JECNominal DefaultME:-1",

   # JJQCD
   "Name:JJQCD_SIG_ghg2_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:JJQCD MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",
   "Name:JJQCD_SIG_ghg4_1_JHUGen_JECNominal Alias:<Name> Process:H0minus Production:JJQCD MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",
   "Name:JJQCD_SIG_ghg2_1_ghg4_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJQCD MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghg2=1,0;ghg4=1,0 Options:SubtractP=JJQCD_SIG_ghg2_1_JHUGen_JECNominal,JJQCD_SIG_ghg4_1_JHUGen_JECNominal DefaultME:-1",
   #"Name:JJQCD_SIG_ghg2_1_ghg4_i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJQCD MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghg2=1,0;ghg4=0,1 Options:SubtractP=JJQCD_SIG_ghg2_1_JHUGen_JECNominal,JJQCD_SIG_ghg4_1_JHUGen_JECNominal DefaultME:-1",

   # Best-DJJ JJVBF and JJQCD
   "Name:JJVBF_SIG_ghv1_1_JHUGen_JECNominal_BestDJJ Copy:JJVBF_SIG_ghv1_1_JHUGen_JECNominal Options:MaxNumerator=JJVBF_SIG_ghv1_1_JHUGen_JECNominal;MaxDenominator=JJQCD_SIG_ghg2_1_JHUGen_JECNominal",
   "Name:JJQCD_SIG_ghg2_1_JHUGen_JECNominal_BestDJJ Copy:JJQCD_SIG_ghg2_1_JHUGen_JECNominal Options:MaxNumerator=JJVBF_SIG_ghv1_1_JHUGen_JECNominal;MaxDenominator=JJQCD_SIG_ghg2_1_JHUGen_JECNominal",

   # JJQCD with QQ initial state
   "Name:JJQCD_InitialQQ_SIG_ghg2_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:JJQCD MatrixElement:JHUGen Cluster:J2JECNominal ForceIncomingFlavors:-21,-21 DefaultME:-1 Options:AddPConst=1",
   "Name:JJQCD_InitialQQ_SIG_ghg4_1_JHUGen_JECNominal Alias:<Name> Process:H0minus Production:JJQCD MatrixElement:JHUGen Cluster:J2JECNominal ForceIncomingFlavors:-21,-21 DefaultME:-1",
   "Name:JJQCD_InitialQQ_SIG_ghg2_1_ghg4_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJQCD MatrixElement:JHUGen Cluster:J2JECNominal ForceIncomingFlavors:-21,-21 Couplings:ghg2=1,0;ghg4=1,0 Options:SubtractP=JJQCD_SIG_ghg2_1_JHUGen_JECNominal,JJQCD_SIG_ghg4_1_JHUGen_JECNominal DefaultME:-1",
   #"Name:JJQCD_InitialQQ_SIG_ghg2_1_ghg4_i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJQCD MatrixElement:JHUGen Cluster:J2JECNominal ForceIncomingFlavors:-21,-21 Couplings:ghg2=1,0;ghg4=0,1 Options:SubtractP=JJQCD_SIG_ghg2_1_JHUGen_JECNominal,JJQCD_SIG_ghg4_1_JHUGen_JECNominal DefaultME:-1",

   # Hadronic ZH
   "Name:HadZH_SIG_ghz1_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",
   "Name:HadZH_SIG_ghz1prime2_1E4_JHUGen_JECNominal Alias:<Name> Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1_prime2=10000,0 DefaultME:-1",
   "Name:HadZH_SIG_ghz2_1_JHUGen_JECNominal Alias:<Name> Process:H0hplus Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",
   "Name:HadZH_SIG_ghz4_1_JHUGen_JECNominal Alias:<Name> Process:H0minus Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",
   "Name:HadZH_SIG_ghza1prime2_1E4_JHUGen_JECNominal Alias:<Name> Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghzgs1_prime2=10000,0 DefaultME:-1",
   "Name:HadZH_SIG_ghza2_1_JHUGen_JECNominal Alias:<Name> Process:H0_Zgs Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",
   "Name:HadZH_SIG_ghza4_1_JHUGen_JECNominal Alias:<Name> Process:H0_Zgs_PS Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",
   "Name:HadZH_SIG_gha2_1_JHUGen_JECNominal Alias:<Name> Process:H0_gsgs Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",
   "Name:HadZH_SIG_gha4_1_JHUGen_JECNominal Alias:<Name> Process:H0_gsgs_PS Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",

   "Name:HadZH_SIG_ghz1_1_ghz1prime2_1E4_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghz1_prime2=10000,0 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghz1prime2_1E4_JHUGen_JECNominal DefaultME:-1",
   #"Name:HadZH_SIG_ghz1_1_ghz1prime2_1E4i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghz1_prime2=0,10000 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghz1prime2_1E4_JHUGen_JECNominal DefaultME:-1",
   "Name:HadZH_SIG_ghz1_1_ghz2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghz2=1,0 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghz2_1_JHUGen_JECNominal DefaultME:-1",
   #"Name:HadZH_SIG_ghz1_1_ghz2_i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghz2=0,1 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghz2_1_JHUGen_JECNominal DefaultME:-1",
   "Name:HadZH_SIG_ghz1_1_ghz4_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghz4=1,0 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghz4_1_JHUGen_JECNominal DefaultME:-1",
   #"Name:HadZH_SIG_ghz1_1_ghz4_i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghz4=0,1 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghz4_1_JHUGen_JECNominal DefaultME:-1",

   "Name:HadZH_SIG_ghz1_1_ghza1prime2_1E4_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghzgs1_prime2=10000,0 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghza1prime2_1E4_JHUGen_JECNominal DefaultME:-1",
   #"Name:HadZH_SIG_ghz1_1_ghza1prime2_1E4i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghzgs1_prime2=0,10000 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghza1prime2_1E4_JHUGen_JECNominal DefaultME:-1",
   "Name:HadZH_SIG_ghz1_1_ghza2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghzgs2=1,0 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghza2_1_JHUGen_JECNominal DefaultME:-1",
   #"Name:HadZH_SIG_ghz1_1_ghza2_i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghzgs2=0,1 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghza2_1_JHUGen_JECNominal DefaultME:-1",
   "Name:HadZH_SIG_ghz1_1_ghza4_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghzgs4=1,0 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghza4_1_JHUGen_JECNominal DefaultME:-1",
   #"Name:HadZH_SIG_ghz1_1_ghza4_i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghzgs4=0,1 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghza4_1_JHUGen_JECNominal DefaultME:-1",

   "Name:HadZH_SIG_ghz1_1_gha2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghgsgs2=1,0 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_gha2_1_JHUGen_JECNominal DefaultME:-1",
   #"Name:HadZH_SIG_ghz1_1_gha2_i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghgsgs2=0,1 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_gha2_1_JHUGen_JECNominal DefaultME:-1",
   "Name:HadZH_SIG_ghz1_1_gha4_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghgsgs4=1,0 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_gha4_1_JHUGen_JECNominal DefaultME:-1",
   #"Name:HadZH_SIG_ghz1_1_gha4_i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghgsgs4=0,1 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_gha4_1_JHUGen_JECNominal DefaultME:-1",

   # Hadronic WH
   "Name:HadWH_SIG_ghw1_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:Had_WH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",
   "Name:HadWH_SIG_ghw1prime2_1E4_JHUGen_JECNominal Alias:<Name> Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1_prime2=10000,0 DefaultME:-1",
   "Name:HadWH_SIG_ghw2_1_JHUGen_JECNominal Alias:<Name> Process:H0hplus Production:Had_WH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",
   "Name:HadWH_SIG_ghw4_1_JHUGen_JECNominal Alias:<Name> Process:H0minus Production:Had_WH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",

   "Name:HadWH_SIG_ghw1_1_ghw1prime2_1E4_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghz1_prime2=10000,0 Options:SubtractP=HadWH_SIG_ghw1_1_JHUGen_JECNominal,HadWH_SIG_ghw1prime2_1E4_JHUGen_JECNominal DefaultME:-1",
   #"Name:HadWH_SIG_ghw1_1_ghw1prime2_1E4i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghz1_prime2=0,10000 Options:SubtractP=HadWH_SIG_ghw1_1_JHUGen_JECNominal,HadWH_SIG_ghw1prime2_1E4_JHUGen_JECNominal DefaultME:-1",
   "Name:HadWH_SIG_ghw1_1_ghw2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghz2=1,0 Options:SubtractP=HadWH_SIG_ghw1_1_JHUGen_JECNominal,HadWH_SIG_ghw2_1_JHUGen_JECNominal DefaultME:-1",
   #"Name:HadWH_SIG_ghw1_1_ghw2_i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghz2=0,1 Options:SubtractP=HadWH_SIG_ghw1_1_JHUGen_JECNominal,HadWH_SIG_ghw2_1_JHUGen_JECNominal DefaultME:-1",
   "Name:HadWH_SIG_ghw1_1_ghw4_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghz4=1,0 Options:SubtractP=HadWH_SIG_ghw1_1_JHUGen_JECNominal,HadWH_SIG_ghw4_1_JHUGen_JECNominal DefaultME:-1",
   #"Name:HadWH_SIG_ghw1_1_ghw4_i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:ghz1=1,0;ghz4=0,1 Options:SubtractP=HadWH_SIG_ghw1_1_JHUGen_JECNominal,HadWH_SIG_ghw4_1_JHUGen_JECNominal DefaultME:-1",

   # ttH: Undecayed MEs belong to J2-class clusters
   "Name:ttHUndecayed_SIG_kappa_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:ttH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",
   "Name:ttHUndecayed_SIG_kappatilde_1_JHUGen_JECNominal Alias:<Name> Process:H0minus Production:ttH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",
   "Name:ttHUndecayed_SIG_kappa_1_kappatilde_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:ttH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:kappa=1,0;kappa_tilde=1,0 Options:SubtractP=ttH_SIG_kappa_1_JHUGen_JECNominal,ttH_SIG_kappatilde_1_JHUGen_JECNominal DefaultME:-1",
   #"Name:ttHUndecayed_SIG_kappa_1_kappatilde_i_JHUGen_JECNominal Process:SelfDefine_spin0 Production:ttH MatrixElement:JHUGen Cluster:J2JECNominal Couplings:kappa=1,0;kappa_tilde=0,1 Options:SubtractP=ttH_SIG_kappa_1_JHUGen_JECNominal,ttH_SIG_kappatilde_1_JHUGen_JECNominal DefaultME:-1",

   # bbH
   # Not adding kappa_tilde since it was shown bbH has close to no sensitivity any time soon
   "Name:bbH_SIG_kappa_1_JHUGen_JECNominal Process:HSMHiggs Production:bbH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1",
]
AJetsProdProbabilities_SpinZero_JHUGen_JECUp = [theME.replace("JECNominal", "JECUp") for theME in AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
AJetsProdProbabilities_SpinZero_JHUGen_JECDn = [theME.replace("JECNominal", "JECDn") for theME in AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
## Production probabilities with >=1 lepton(s) ##
ALepsProdProbabilities_SpinZero_JHUGen = [
   # Leptonic ZH
   "Name:LepZH_SIG_ghz1_1_JHUGen Alias:<Name> Process:HSMHiggs Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH DefaultME:-1",
   "Name:LepZH_SIG_ghz1prime2_1E4_JHUGen Alias:<Name> Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH Couplings:ghz1_prime2=10000,0 DefaultME:-1",
   "Name:LepZH_SIG_ghz2_1_JHUGen Alias:<Name> Process:H0hplus Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH DefaultME:-1",
   "Name:LepZH_SIG_ghz4_1_JHUGen Alias:<Name> Process:H0minus Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH DefaultME:-1",
   "Name:LepZH_SIG_ghza1prime2_1E4_JHUGen Alias:<Name> Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH Couplings:ghzgs1_prime2=10000,0 DefaultME:-1",
   "Name:LepZH_SIG_ghza2_1_JHUGen Alias:<Name> Process:H0_Zgs Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH DefaultME:-1",
   "Name:LepZH_SIG_ghza4_1_JHUGen Alias:<Name> Process:H0_Zgs_PS Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH DefaultME:-1",
   "Name:LepZH_SIG_gha2_1_JHUGen Alias:<Name> Process:H0_gsgs Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH DefaultME:-1",
   "Name:LepZH_SIG_gha4_1_JHUGen Alias:<Name> Process:H0_gsgs_PS Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH DefaultME:-1",

   "Name:LepZH_SIG_ghz1_1_ghz1prime2_1E4_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH Couplings:ghz1=1,0;ghz1_prime2=10000,0 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghz1prime2_1E4_JHUGen DefaultME:-1",
   #"Name:LepZH_SIG_ghz1_1_ghz1prime2_1E4i_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH Couplings:ghz1=1,0;ghz1_prime2=0,10000 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghz1prime2_1E4_JHUGen DefaultME:-1",
   "Name:LepZH_SIG_ghz1_1_ghz2_1_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH Couplings:ghz1=1,0;ghz2=1,0 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghz2_1_JHUGen DefaultME:-1",
   #"Name:LepZH_SIG_ghz1_1_ghz2_i_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH Couplings:ghz1=1,0;ghz2=0,1 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghz2_1_JHUGen DefaultME:-1",
   "Name:LepZH_SIG_ghz1_1_ghz4_1_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH Couplings:ghz1=1,0;ghz4=1,0 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghz4_1_JHUGen DefaultME:-1",
   #"Name:LepZH_SIG_ghz1_1_ghz4_i_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH Couplings:ghz1=1,0;ghz4=0,1 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghz4_1_JHUGen DefaultME:-1",

   "Name:LepZH_SIG_ghz1_1_ghza1prime2_1E4_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH Couplings:ghz1=1,0;ghzgs1_prime2=10000,0 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghza1prime2_1E4_JHUGen DefaultME:-1",
   #"Name:LepZH_SIG_ghz1_1_ghza1prime2_1E4i_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH Couplings:ghz1=1,0;ghzgs1_prime2=0,10000 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghza1prime2_1E4_JHUGen DefaultME:-1",
   "Name:LepZH_SIG_ghz1_1_ghza2_1_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH Couplings:ghz1=1,0;ghzgs2=1,0 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghza2_1_JHUGen DefaultME:-1",
   #"Name:LepZH_SIG_ghz1_1_ghza2_i_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH Couplings:ghz1=1,0;ghzgs2=0,1 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghza2_1_JHUGen DefaultME:-1",
   "Name:LepZH_SIG_ghz1_1_ghza4_1_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH Couplings:ghz1=1,0;ghzgs4=1,0 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghza4_1_JHUGen DefaultME:-1",
   #"Name:LepZH_SIG_ghz1_1_ghza4_i_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH Couplings:ghz1=1,0;ghzgs4=0,1 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghza4_1_JHUGen DefaultME:-1",

   "Name:LepZH_SIG_ghz1_1_gha2_1_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH Couplings:ghz1=1,0;ghgsgs2=1,0 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_gha2_1_JHUGen DefaultME:-1",
   #"Name:LepZH_SIG_ghz1_1_gha2_i_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH Couplings:ghz1=1,0;ghgsgs2=0,1 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_gha2_1_JHUGen DefaultME:-1",
   "Name:LepZH_SIG_ghz1_1_gha4_1_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH Couplings:ghz1=1,0;ghgsgs4=1,0 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_gha4_1_JHUGen DefaultME:-1",
   #"Name:LepZH_SIG_ghz1_1_gha4_i_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster:LepZH Couplings:ghz1=1,0;ghgsgs4=0,1 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_gha4_1_JHUGen DefaultME:-1",

   # Leptonic WH (CAUTION: All requiring the SM ME to be maximized)
   "Name:LepWH_SIG_ghw1_1_JHUGen Alias:<Name> Process:HSMHiggs Production:Lep_WH MatrixElement:JHUGen Cluster:LepWH Options:MaxNumerator=LepWH_SIG_ghw1_1_JHUGen DefaultME:-1",
   "Name:LepWH_SIG_ghw1prime2_1E4_JHUGen Alias:<Name> Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster:LepWH Options:MaxNumerator=LepWH_SIG_ghw1_1 Couplings:ghz1_prime2=10000,0_JHUGen DefaultME:-1",
   "Name:LepWH_SIG_ghw2_1_JHUGen Alias:<Name> Process:H0hplus Production:Lep_WH MatrixElement:JHUGen Cluster:LepWH Options:MaxNumerator=LepWH_SIG_ghw1_1_JHUGen DefaultME:-1",
   "Name:LepWH_SIG_ghw4_1_JHUGen Alias:<Name> Process:H0minus Production:Lep_WH MatrixElement:JHUGen Cluster:LepWH Options:MaxNumerator=LepWH_SIG_ghw1_1_JHUGen DefaultME:-1",

   "Name:LepWH_SIG_ghw1_1_ghw1prime2_1E4_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster:LepWH Couplings:ghz1=1,0;ghz1_prime2=10000,0 Options:SubtractP=LepWH_SIG_ghw1_1_JHUGen,LepWH_SIG_ghw1prime2_1E4_JHUGen;MaxNumerator=LepWH_SIG_ghw1_1_JHUGen DefaultME:-1",
   #"Name:LepWH_SIG_ghw1_1_ghw1prime2_1E4i_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster:LepWH Couplings:ghz1=1,0;ghz1_prime2=0,10000 Options:SubtractP=LepWH_SIG_ghw1_1_JHUGen,LepWH_SIG_ghw1prime2_1E4_JHUGen;MaxNumerator=LepWH_SIG_ghw1_1_JHUGen DefaultME:-1",
   "Name:LepWH_SIG_ghw1_1_ghw2_1_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster:LepWH Couplings:ghz1=1,0;ghz2=1,0 Options:SubtractP=LepWH_SIG_ghw1_1_JHUGen,LepWH_SIG_ghw2_1_JHUGen;MaxNumerator=LepWH_SIG_ghw1_1_JHUGen DefaultME:-1",
   #"Name:LepWH_SIG_ghw1_1_ghw2_i_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster:LepWH Couplings:ghz1=1,0;ghz2=0,1 Options:SubtractP=LepWH_SIG_ghw1_1_JHUGen,LepWH_SIG_ghw2_1_JHUGen;MaxNumerator=LepWH_SIG_ghw1_1_JHUGen DefaultME:-1",
   "Name:LepWH_SIG_ghw1_1_ghw4_1_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster:LepWH Couplings:ghz1=1,0;ghz4=1,0 Options:SubtractP=LepWH_SIG_ghw1_1_JHUGen,LepWH_SIG_ghw4_1_JHUGen;MaxNumerator=LepWH_SIG_ghw1_1_JHUGen DefaultME:-1",
   #"Name:LepWH_SIG_ghw1_1_ghw4_i_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster:LepWH Couplings:ghz1=1,0;ghz4=0,1 Options:SubtractP=LepWH_SIG_ghw1_1_JHUGen,LepWH_SIG_ghw4_1_JHUGen;MaxNumerator=LepWH_SIG_ghw1_1_JHUGen DefaultME:-1",
]

### Spin-1 decay probabilities from JHUGen ###
Probabilities_SpinOne_JHUGen = [
   "Name:QQB_SIG_ZPqqLR_1_gZPz1_1_JHUGen Process:H1minus Production:ZZQQB MatrixElement:JHUGen",
   "Name:QQB_SIG_ZPqqLR_1_gZPz2_1_JHUGen Process:H1plus Production:ZZQQB MatrixElement:JHUGen",
   "Name:INDEPENDENT_SIG_gZPz1_1_JHUGen Process:H1minus Production:ZZINDEPENDENT MatrixElement:JHUGen",
   "Name:INDEPENDENT_SIG_gZPz2_1_JHUGen Process:H1plus Production:ZZINDEPENDENT MatrixElement:JHUGen",
]
### Spin-2 decay probabilities from JHUGen ###
Probabilities_SpinTwo_JHUGen = [
   "Name:GG_SIG_gXg1_1_gXz1_1_JHUGen Process:H2_g1 Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_gXg2_1_gXz2_1_JHUGen Process:H2_g2 Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_gXg3_1_gXz3_1_JHUGen Process:H2_g3 Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_gXg4_1_gXz4_1_JHUGen Process:H2_g4 Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_gXg1_1_gXz5_1_JHUGen Process:H2_g5 Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_gXg1_1_gXz1_1_gXz5_1_JHUGen Process:H2_g1g5 Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_gXg1_1_gXz6_1_JHUGen Process:H2_g6 Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_gXg1_1_gXz7_1_JHUGen Process:H2_g7 Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_gXg5_1_gXz8_1_JHUGen Process:H2_g8 Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_gXg5_1_gXz9_1_JHUGen Process:H2_g9 Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_gXg5_1_gXz10_1_JHUGen Process:H2_g10 Production:ZZGG MatrixElement:JHUGen",

   "Name:QQB_SIG_XqqLR_1_gXz1_1_JHUGen Process:H2_g1 Production:ZZQQB MatrixElement:JHUGen",
   "Name:QQB_SIG_XqqLR_1_gXz2_1_JHUGen Process:H2_g2 Production:ZZQQB MatrixElement:JHUGen",
   "Name:QQB_SIG_XqqLR_1_gXz3_1_JHUGen Process:H2_g3 Production:ZZQQB MatrixElement:JHUGen",
   "Name:QQB_SIG_XqqLR_1_gXz4_1_JHUGen Process:H2_g4 Production:ZZQQB MatrixElement:JHUGen",
   "Name:QQB_SIG_XqqLR_1_gXz5_1_JHUGen Process:H2_g5 Production:ZZQQB MatrixElement:JHUGen",
   "Name:QQB_SIG_XqqLR_1_gXz1_1_gXz5_1_JHUGen Process:H2_g1g5 Production:ZZQQB MatrixElement:JHUGen",
   "Name:QQB_SIG_XqqLR_1_gXz6_1_JHUGen Process:H2_g6 Production:ZZQQB MatrixElement:JHUGen",
   "Name:QQB_SIG_XqqLR_1_gXz7_1_JHUGen Process:H2_g7 Production:ZZQQB MatrixElement:JHUGen",
   "Name:QQB_SIG_XqqLR_1_gXz8_1_JHUGen Process:H2_g8 Production:ZZQQB MatrixElement:JHUGen",
   "Name:QQB_SIG_XqqLR_1_gXz9_1_JHUGen Process:H2_g9 Production:ZZQQB MatrixElement:JHUGen",
   "Name:QQB_SIG_XqqLR_1_gXz10_1_JHUGen Process:H2_g10 Production:ZZQQB MatrixElement:JHUGen",

   "Name:INDEPENDENT_SIG_gXz1_1_JHUGen Process:H2_g1 Production:ZZINDEPENDENT MatrixElement:JHUGen",
   "Name:INDEPENDENT_SIG_gXz2_1_JHUGen Process:H2_g2 Production:ZZINDEPENDENT MatrixElement:JHUGen",
   "Name:INDEPENDENT_SIG_gXz3_1_JHUGen Process:H2_g3 Production:ZZINDEPENDENT MatrixElement:JHUGen",
   "Name:INDEPENDENT_SIG_gXz4_1_JHUGen Process:H2_g4 Production:ZZINDEPENDENT MatrixElement:JHUGen",
   "Name:INDEPENDENT_SIG_gXz5_1_JHUGen Process:H2_g5 Production:ZZINDEPENDENT MatrixElement:JHUGen",
   "Name:INDEPENDENT_SIG_gXz1_1_gXz5_1_JHUGen Process:H2_g1g5 Production:ZZINDEPENDENT MatrixElement:JHUGen",
   "Name:INDEPENDENT_SIG_gXz6_1_JHUGen Process:H2_g6 Production:ZZINDEPENDENT MatrixElement:JHUGen",
   "Name:INDEPENDENT_SIG_gXz7_1_JHUGen Process:H2_g7 Production:ZZINDEPENDENT MatrixElement:JHUGen",
   "Name:INDEPENDENT_SIG_gXz8_1_JHUGen Process:H2_g8 Production:ZZINDEPENDENT MatrixElement:JHUGen",
   "Name:INDEPENDENT_SIG_gXz9_1_JHUGen Process:H2_g9 Production:ZZINDEPENDENT MatrixElement:JHUGen",
   "Name:INDEPENDENT_SIG_gXz10_1_JHUGen Process:H2_g10 Production:ZZINDEPENDENT MatrixElement:JHUGen",
]

### Decay probabilities from MCFM ###
DecayProbabilities_MCFM = [
   "Name:GG_SIG_kappaTopBot_1_ghz1_1_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Options:AddPConst=1",
   "Name:GG_BSI_kappaTopBot_1_ghz1_1_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM", # Has pConst=pConst_sig+pConst_bkg
   "Name:GG_BSI_kappaTopBot_1_ghz1_i_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=0,1", # Has the same pConst
   "Name:GG_BKG_MCFM Process:bkgZZ Production:ZZGG MatrixElement:MCFM Options:AddPConst=1",
   "Name:QQB_BKG_MCFM Alias:<Name> Process:bkgZZ Production:ZZQQB MatrixElement:MCFM Options:AddPConst=1", # Aliased to construct Dbkgkin
   #"Name:INDEPENDENT_BKG_MCFM Process:bkgZZ Production:ZZINDEPENDENT MatrixElement:MCFM",
   "Name:ZJJ_BKG_MCFM Process:bkgZJets Production:JJQCD MatrixElement:MCFM Options:AddPConst=1",
]
AJetsProdDecProbabilities_MCFM_JECNominal = [
   "Name:JJEW_SIG_ghv1_1_MCFM_JECNominal Process:HSMHiggs Production:JJEW MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1",
   "Name:JJEW_BSI_ghv1_1_MCFM_JECNominal Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1",
   "Name:JJEW_BSI_ghv1_i_MCFM_JECNominal Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=0,1 Cluster:J2JECNominal DefaultME:-1",
   "Name:JJEW_BKG_MCFM_JECNominal Process:bkgZZ Production:JJEW MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1",

   "Name:JJVBF_S_SIG_ghv1_1_MCFM_JECNominal Process:HSMHiggs Production:JJVBF_S MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",
   "Name:JJVBF_S_BSI_ghv1_1_MCFM_JECNominal Process:bkgZZ_SMHiggs Production:JJVBF_S MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1",
   "Name:JJVBF_S_BSI_ghv1_i_MCFM_JECNominal Process:bkgZZ_SMHiggs Production:JJVBF_S MatrixElement:MCFM Couplings:ghz1=0,1 Cluster:J2JECNominal DefaultME:-1",
   "Name:JJVBF_SIG_ghv1_1_MCFM_JECNominal Process:HSMHiggs Production:JJVBF MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1",
   "Name:JJVBF_BSI_ghv1_1_MCFM_JECNominal Process:bkgZZ_SMHiggs Production:JJVBF MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1",
   "Name:JJVBF_BSI_ghv1_i_MCFM_JECNominal Process:bkgZZ_SMHiggs Production:JJVBF MatrixElement:MCFM Couplings:ghz1=0,1 Cluster:J2JECNominal DefaultME:-1",
   "Name:JJVBF_BKG_MCFM_JECNominal Process:bkgZZ Production:JJVBF MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",

   "Name:HadZH_S_SIG_ghz1_1_MCFM_JECNominal Process:HSMHiggs Production:Had_ZH_S MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",
   "Name:HadZH_S_BSI_ghz1_1_MCFM_JECNominal Process:bkgZZ_SMHiggs Production:Had_ZH_S MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1",
   "Name:HadZH_S_BSI_ghz1_i_MCFM_JECNominal Process:bkgZZ_SMHiggs Production:Had_ZH_S MatrixElement:MCFM Couplings:ghz1=0,1 Cluster:J2JECNominal DefaultME:-1",
   "Name:HadZH_SIG_ghz1_1_MCFM_JECNominal Process:HSMHiggs Production:Had_ZH MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1",
   "Name:HadZH_BSI_ghz1_1_MCFM_JECNominal Process:bkgZZ_SMHiggs Production:Had_ZH MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1",
   "Name:HadZH_BSI_ghz1_i_MCFM_JECNominal Process:bkgZZ_SMHiggs Production:Had_ZH MatrixElement:MCFM Couplings:ghz1=0,1 Cluster:J2JECNominal DefaultME:-1",
   "Name:HadZH_BKG_MCFM_JECNominal Process:bkgZZ Production:Had_ZH MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",

   "Name:HadWH_S_SIG_ghw1_1_MCFM_JECNominal Process:HSMHiggs Production:Had_WH_S MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",
   "Name:HadWH_S_BSI_ghw1_1_MCFM_JECNominal Process:bkgZZ_SMHiggs Production:Had_WH_S MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1",
   "Name:HadWH_S_BSI_ghw1_i_MCFM_JECNominal Process:bkgZZ_SMHiggs Production:Had_WH_S MatrixElement:MCFM Couplings:ghz1=0,1 Cluster:J2JECNominal DefaultME:-1",
   # In ZZ4l, HadWH_*_S vs HadWH_* should be the same; there is no u/t channel diagram or a valid lepton exchange
   #"Name:HadWH_SIG_ghw1_1_MCFM_JECNominal Process:HSMHiggs Production:Had_WH MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1",
   #"Name:HadWH_BSI_ghw1_1_MCFM_JECNominal Process:bkgZZ_SMHiggs Production:Had_WH MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1",
   #"Name:HadWH_BSI_ghw1_i_MCFM_JECNominal Process:bkgZZ_SMHiggs Production:Had_WH MatrixElement:MCFM Couplings:ghz1=0,1 Cluster:J2JECNominal DefaultME:-1",
   "Name:HadWH_BKG_MCFM_JECNominal Process:bkgZZ Production:Had_WH MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",

   "Name:JJQCD_BKG_MCFM_JECNominal Process:bkgZZ Production:JJQCD MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",
]
AJetsProdDecProbabilities_MCFM_JECUp = [theME.replace("JECNominal", "JECUp") for theME in AJetsProdDecProbabilities_MCFM_JECNominal]
AJetsProdDecProbabilities_MCFM_JECDn = [theME.replace("JECNominal", "JECDn") for theME in AJetsProdDecProbabilities_MCFM_JECNominal]

### m4l probabilities from SuperMELA ###
PM4L_SUPERMELA = [
   "Name:m4l_SIG Process:HSMHiggs Production:ZZGG MatrixElement:JHUGen SuperMelaSyst:SMSyst_None isPM4L:1",
   "Name:m4l_BKG Process:bkgZZ Production:ZZGG MatrixElement:JHUGen SuperMelaSyst:SMSyst_None isPM4L:1",

   "Name:m4l_SIG_ScaleDown Process:HSMHiggs Production:ZZGG MatrixElement:JHUGen SuperMelaSyst:SMSyst_ScaleDown isPM4L:1",
   "Name:m4l_BKG_ScaleDown Process:bkgZZ Production:ZZGG MatrixElement:JHUGen SuperMelaSyst:SMSyst_ScaleDown isPM4L:1",
   "Name:m4l_SIG_ResDown Process:HSMHiggs Production:ZZGG MatrixElement:JHUGen SuperMelaSyst:SMSyst_ResDown isPM4L:1",
   "Name:m4l_BKG_ResDown Process:bkgZZ Production:ZZGG MatrixElement:JHUGen SuperMelaSyst:SMSyst_ResDown isPM4L:1",

   "Name:m4l_SIG_ScaleUp Process:HSMHiggs Production:ZZGG MatrixElement:JHUGen SuperMelaSyst:SMSyst_ScaleUp isPM4L:1",
   "Name:m4l_BKG_ScaleUp Process:bkgZZ Production:ZZGG MatrixElement:JHUGen SuperMelaSyst:SMSyst_ScaleUp isPM4L:1",
   "Name:m4l_SIG_ResUp Process:HSMHiggs Production:ZZGG MatrixElement:JHUGen SuperMelaSyst:SMSyst_ResUp isPM4L:1",
   "Name:m4l_BKG_ResUp Process:bkgZZ Production:ZZGG MatrixElement:JHUGen SuperMelaSyst:SMSyst_ResUp isPM4L:1",
]

### mJJ probabilities in associated production
PMAVJJ_SUPERDIJETMELA_JECNominal = [
   "Name:HadZH_mavjj_true_JECNominal Process:HSMHiggs Production:Had_ZH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 isPMaVJJTrue:1",
   "Name:HadWH_mavjj_true_JECNominal Process:HSMHiggs Production:Had_WH MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 isPMaVJJTrue:1",
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
