### Spin-0 decay probabilities from JHUGen ###
DecayProbabilities_SpinZero_JHUGen = [
   "Name:GG_SIG_ghg2_1_ghz1_1_JHUGen Alias:<Name> Process:HSMHiggs Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_ghg2_1_ghz1prime2_1_JHUGen Alias:<Name> Process:H0_g1prime2 Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_ghg2_1_ghz2_1_JHUGen Alias:<Name> Process:H0hplus Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_ghg2_1_ghz4_1_JHUGen Alias:<Name> Process:H0minus Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_ghg2_1_ghza1prime2_1_JHUGen Alias:<Name> Process:H0_Zgsg1prime2 Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_ghg2_1_ghza2_1_JHUGen Alias:<Name> Process:H0_Zgs Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_ghg2_1_ghza4_1_JHUGen Alias:<Name> Process:H0_Zgs_PS Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_ghg2_1_gha2_1_JHUGen Alias:<Name> Process:H0_gsgs Production:ZZGG MatrixElement:JHUGen",
   "Name:GG_SIG_ghg2_1_gha4_1_JHUGen Alias:<Name> Process:H0_gsgs_PS Production:ZZGG MatrixElement:JHUGen",

   "Name:GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghz1_prime2=1,0 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghz1prime2_1_JHUGen",
   #"Name:GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1_pi2_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghz1_prime2=0,1 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghz1prime2_1_JHUGen",
   "Name:GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghz2=1,0 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghz2_1_JHUGen",
   "Name:GG_SIG_ghg2_1_ghz1_1_ghz2_1_pi2_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghz2=0,1 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghz2_1_JHUGen",
   "Name:GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghz4=1,0 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghz4_1_JHUGen",
   "Name:GG_SIG_ghg2_1_ghz1_1_ghz4_1_pi2_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghz4=0,1 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghz4_1_JHUGen",

   "Name:GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghzgs1_prime2=1,0 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghza1prime2_1_JHUGen",
   "Name:GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1_pi2_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghzgs1_prime2=0,1 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghza1prime2_1_JHUGen",
   "Name:GG_SIG_ghg2_1_ghz1_1_ghza2_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghzgs2=1,0 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghza2_1_JHUGen",
   #"Name:GG_SIG_ghg2_1_ghz1_1_ghza2_1_pi2_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghzgs2=0,1 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghza2_1_JHUGen",
   "Name:GG_SIG_ghg2_1_ghz1_1_ghza4_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghzgs4=1,0 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghza4_1_JHUGen",
   #"Name:GG_SIG_ghg2_1_ghz1_1_ghza4_1_pi2_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghzgs4=0,1 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghza4_1_JHUGen",

   "Name:GG_SIG_ghg2_1_ghz1_1_gha2_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghgsgs2=1,0 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_gha2_1_JHUGen",
   #"Name:GG_SIG_ghg2_1_ghz1_1_gha2_1_pi2_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghgsgs2=0,1 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_gha2_1_JHUGen",
   "Name:GG_SIG_ghg2_1_ghz1_1_gha4_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghgsgs4=1,0 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_gha4_1_JHUGen",
   #"Name:GG_SIG_ghg2_1_ghz1_1_gha4_1_pi2_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghgsgs4=0,1 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_gha4_1_JHUGen",
]
## Production probabilities with >=1 jet(s) ##
AJetsProdProbabilities_SpinZero_JHUGen_JECNominal = [
   # JVBF
   "Name:JVBF_SIG_ghz1_1_JHUGen_JECNominal Process:HSMHiggs Production:JJVBF MatrixElement:JHUGen Cluster=J1JECNominal Options:AddPAux",

   # JQCD
   "Name:JQCD_SIG_ghg2_1_JHUGen_JECNominal Process:HSMHiggs Production:JQCD MatrixElement:JHUGen Cluster=J1JECNominal",
   #"Name:JQCD_SIG_ghg4_1_JHUGen_JECNominal Process:H0minus Production:JQCD MatrixElement:JHUGen Cluster=J1JECNominal",

   # JJVBF
   "Name:JJVBF_SIG_ghz1_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:JJVBF_SIG_ghz1prime2_1_JHUGen_JECNominal Alias:<Name> Process:H0_g1prime2 Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:JJVBF_SIG_ghz2_1_JHUGen_JECNominal Alias:<Name> Process:H0hplus Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:JJVBF_SIG_ghz4_1_JHUGen_JECNominal Alias:<Name> Process:H0minus Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:JJVBF_SIG_ghza1prime2_1_JHUGen_JECNominal Alias:<Name> Process:H0_Zgsg1prime2 Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:JJVBF_SIG_ghza2_1_JHUGen_JECNominal Alias:<Name> Process:H0_Zgs Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:JJVBF_SIG_ghza4_1_JHUGen_JECNominal Alias:<Name> Process:H0_Zgs_PS Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:JJVBF_SIG_gha2_1_JHUGen_JECNominal Alias:<Name> Process:H0_gsgs Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:JJVBF_SIG_gha4_1_JHUGen_JECNominal Alias:<Name> Process:H0_gsgs_PS Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal",

   "Name:JJVBF_SIG_ghz1_1_ghz1prime2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghz1_prime2=1,0 Options:SubtractP=JJVBF_SIG_ghz1_1_JHUGen_JECNominal,JJVBF_SIG_ghz1prime2_1_JHUGen_JECNominal",
   #"Name:JJVBF_SIG_ghz1_1_ghz1prime2_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghz1_prime2=0,1 Options:SubtractP=JJVBF_SIG_ghz1_1_JHUGen_JECNominal,JJVBF_SIG_ghz1prime2_1_JHUGen_JECNominal",
   "Name:JJVBF_SIG_ghz1_1_ghz2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghz2=1,0 Options:SubtractP=JJVBF_SIG_ghz1_1_JHUGen_JECNominal,JJVBF_SIG_ghz2_1_JHUGen_JECNominal",
   #"Name:JJVBF_SIG_ghz1_1_ghz2_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghz2=0,1 Options:SubtractP=JJVBF_SIG_ghz1_1_JHUGen_JECNominal,JJVBF_SIG_ghz2_1_JHUGen_JECNominal",
   "Name:JJVBF_SIG_ghz1_1_ghz4_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghz4=1,0 Options:SubtractP=JJVBF_SIG_ghz1_1_JHUGen_JECNominal,JJVBF_SIG_ghz4_1_JHUGen_JECNominal",
   #"Name:JJVBF_SIG_ghz1_1_ghz4_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghz4=0,1 Options:SubtractP=JJVBF_SIG_ghz1_1_JHUGen_JECNominal,JJVBF_SIG_ghz4_1_JHUGen_JECNominal",

   "Name:JJVBF_SIG_ghz1_1_ghza1prime2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghzgs1_prime2=1,0 Options:SubtractP=JJVBF_SIG_ghz1_1_JHUGen_JECNominal,JJVBF_SIG_ghza1prime2_1_JHUGen_JECNominal",
   #"Name:JJVBF_SIG_ghz1_1_ghza1prime2_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghzgs1_prime2=0,1 Options:SubtractP=JJVBF_SIG_ghz1_1_JHUGen_JECNominal,JJVBF_SIG_ghza1prime2_1_JHUGen_JECNominal",
   "Name:JJVBF_SIG_ghz1_1_ghza2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghzgs2=1,0 Options:SubtractP=JJVBF_SIG_ghz1_1_JHUGen_JECNominal,JJVBF_SIG_ghza2_1_JHUGen_JECNominal",
   #"Name:JJVBF_SIG_ghz1_1_ghza2_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghzgs2=0,1 Options:SubtractP=JJVBF_SIG_ghz1_1_JHUGen_JECNominal,JJVBF_SIG_ghza2_1_JHUGen_JECNominal",
   "Name:JJVBF_SIG_ghz1_1_ghza4_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghzgs4=1,0 Options:SubtractP=JJVBF_SIG_ghz1_1_JHUGen_JECNominal,JJVBF_SIG_ghza4_1_JHUGen_JECNominal",
   #"Name:JJVBF_SIG_ghz1_1_ghza4_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghzgs4=0,1 Options:SubtractP=JJVBF_SIG_ghz1_1_JHUGen_JECNominal,JJVBF_SIG_ghza4_1_JHUGen_JECNominal",

   "Name:JJVBF_SIG_ghz1_1_gha2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghgsgs2=1,0 Options:SubtractP=JJVBF_SIG_ghz1_1_JHUGen_JECNominal,JJVBF_SIG_gha2_1_JHUGen_JECNominal",
   #"Name:JJVBF_SIG_ghz1_1_gha2_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghgsgs2=0,1 Options:SubtractP=JJVBF_SIG_ghz1_1_JHUGen_JECNominal,JJVBF_SIG_gha2_1_JHUGen_JECNominal",
   "Name:JJVBF_SIG_ghz1_1_gha4_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghgsgs4=1,0 Options:SubtractP=JJVBF_SIG_ghz1_1_JHUGen_JECNominal,JJVBF_SIG_gha4_1_JHUGen_JECNominal",
   #"Name:JJVBF_SIG_ghz1_1_gha4_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJVBF MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghgsgs4=0,1 Options:SubtractP=JJVBF_SIG_ghz1_1_JHUGen_JECNominal,JJVBF_SIG_gha4_1_JHUGen_JECNominal",

   # JJQCD
   "Name:JJQCD_SIG_ghg2_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:JJQCD MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:JJQCD_SIG_ghg4_1_JHUGen_JECNominal Alias:<Name> Process:H0minus Production:JJQCD MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:JJQCD_SIG_ghg2_1_ghg4_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJQCD MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghg2=1,0;ghg4=1,0 Options:SubtractP=JJQCD_SIG_ghg2_1_JHUGen_JECNominal,JJQCD_SIG_ghg4_1_JHUGen_JECNominal",
   #"Name:JJQCD_SIG_ghg2_1_ghg4_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:JJQCD MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghg2=1,0;ghg4=0,1 Options:SubtractP=JJQCD_SIG_ghg2_1_JHUGen_JECNominal,JJQCD_SIG_ghg4_1_JHUGen_JECNominal",

   # Hadronic ZH
   "Name:HadZH_SIG_ghz1_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:HadZH_SIG_ghz1prime2_1_JHUGen_JECNominal Alias:<Name> Process:H0_g1prime2 Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:HadZH_SIG_ghz2_1_JHUGen_JECNominal Alias:<Name> Process:H0hplus Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:HadZH_SIG_ghz4_1_JHUGen_JECNominal Alias:<Name> Process:H0minus Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:HadZH_SIG_ghza1prime2_1_JHUGen_JECNominal Alias:<Name> Process:H0_Zgsg1prime2 Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:HadZH_SIG_ghza2_1_JHUGen_JECNominal Alias:<Name> Process:H0_Zgs Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:HadZH_SIG_ghza4_1_JHUGen_JECNominal Alias:<Name> Process:H0_Zgs_PS Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:HadZH_SIG_gha2_1_JHUGen_JECNominal Alias:<Name> Process:H0_gsgs Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:HadZH_SIG_gha4_1_JHUGen_JECNominal Alias:<Name> Process:H0_gsgs_PS Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal",

   "Name:HadZH_SIG_ghz1_1_ghz1prime2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghz1_prime2=1,0 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghz1prime2_1_JHUGen_JECNominal",
   #"Name:HadZH_SIG_ghz1_1_ghz1prime2_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghz1_prime2=0,1 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghz1prime2_1_JHUGen_JECNominal",
   "Name:HadZH_SIG_ghz1_1_ghz2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghz2=1,0 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghz2_1_JHUGen_JECNominal",
   #"Name:HadZH_SIG_ghz1_1_ghz2_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghz2=0,1 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghz2_1_JHUGen_JECNominal",
   "Name:HadZH_SIG_ghz1_1_ghz4_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghz4=1,0 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghz4_1_JHUGen_JECNominal",
   #"Name:HadZH_SIG_ghz1_1_ghz4_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghz4=0,1 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghz4_1_JHUGen_JECNominal",

   "Name:HadZH_SIG_ghz1_1_ghza1prime2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghzgs1_prime2=1,0 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghza1prime2_1_JHUGen_JECNominal",
   #"Name:HadZH_SIG_ghz1_1_ghza1prime2_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghzgs1_prime2=0,1 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghza1prime2_1_JHUGen_JECNominal",
   "Name:HadZH_SIG_ghz1_1_ghza2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghzgs2=1,0 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghza2_1_JHUGen_JECNominal",
   #"Name:HadZH_SIG_ghz1_1_ghza2_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghzgs2=0,1 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghza2_1_JHUGen_JECNominal",
   "Name:HadZH_SIG_ghz1_1_ghza4_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghzgs4=1,0 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghza4_1_JHUGen_JECNominal",
   #"Name:HadZH_SIG_ghz1_1_ghza4_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghzgs4=0,1 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_ghza4_1_JHUGen_JECNominal",

   "Name:HadZH_SIG_ghz1_1_gha2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghgsgs2=1,0 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_gha2_1_JHUGen_JECNominal",
   #"Name:HadZH_SIG_ghz1_1_gha2_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghgsgs2=0,1 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_gha2_1_JHUGen_JECNominal",
   "Name:HadZH_SIG_ghz1_1_gha4_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghgsgs4=1,0 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_gha4_1_JHUGen_JECNominal",
   #"Name:HadZH_SIG_ghz1_1_gha4_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_ZH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghgsgs4=0,1 Options:SubtractP=HadZH_SIG_ghz1_1_JHUGen_JECNominal,HadZH_SIG_gha4_1_JHUGen_JECNominal",

   # Hadronic WH
   "Name:HadWH_SIG_ghz1_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:HadWH_SIG_ghz1prime2_1_JHUGen_JECNominal Alias:<Name> Process:H0_g1prime2 Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:HadWH_SIG_ghz2_1_JHUGen_JECNominal Alias:<Name> Process:H0hplus Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:HadWH_SIG_ghz4_1_JHUGen_JECNominal Alias:<Name> Process:H0minus Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:HadWH_SIG_ghza1prime2_1_JHUGen_JECNominal Alias:<Name> Process:H0_Zgsg1prime2 Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:HadWH_SIG_ghza2_1_JHUGen_JECNominal Alias:<Name> Process:H0_Zgs Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:HadWH_SIG_ghza4_1_JHUGen_JECNominal Alias:<Name> Process:H0_Zgs_PS Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:HadWH_SIG_gha2_1_JHUGen_JECNominal Alias:<Name> Process:H0_gsgs Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:HadWH_SIG_gha4_1_JHUGen_JECNominal Alias:<Name> Process:H0_gsgs_PS Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal",

   "Name:HadWH_SIG_ghz1_1_ghz1prime2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghz1_prime2=1,0 Options:SubtractP=HadWH_SIG_ghz1_1_JHUGen_JECNominal,HadWH_SIG_ghz1prime2_1_JHUGen_JECNominal",
   #"Name:HadWH_SIG_ghz1_1_ghz1prime2_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghz1_prime2=0,1 Options:SubtractP=HadWH_SIG_ghz1_1_JHUGen_JECNominal,HadWH_SIG_ghz1prime2_1_JHUGen_JECNominal",
   "Name:HadWH_SIG_ghz1_1_ghz2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghz2=1,0 Options:SubtractP=HadWH_SIG_ghz1_1_JHUGen_JECNominal,HadWH_SIG_ghz2_1_JHUGen_JECNominal",
   #"Name:HadWH_SIG_ghz1_1_ghz2_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghz2=0,1 Options:SubtractP=HadWH_SIG_ghz1_1_JHUGen_JECNominal,HadWH_SIG_ghz2_1_JHUGen_JECNominal",
   "Name:HadWH_SIG_ghz1_1_ghz4_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghz4=1,0 Options:SubtractP=HadWH_SIG_ghz1_1_JHUGen_JECNominal,HadWH_SIG_ghz4_1_JHUGen_JECNominal",
   #"Name:HadWH_SIG_ghz1_1_ghz4_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghz4=0,1 Options:SubtractP=HadWH_SIG_ghz1_1_JHUGen_JECNominal,HadWH_SIG_ghz4_1_JHUGen_JECNominal",

   "Name:HadWH_SIG_ghz1_1_ghza1prime2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghzgs1_prime2=1,0 Options:SubtractP=HadWH_SIG_ghz1_1_JHUGen_JECNominal,HadWH_SIG_ghza1prime2_1_JHUGen_JECNominal",
   #"Name:HadWH_SIG_ghz1_1_ghza1prime2_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghzgs1_prime2=0,1 Options:SubtractP=HadWH_SIG_ghz1_1_JHUGen_JECNominal,HadWH_SIG_ghza1prime2_1_JHUGen_JECNominal",
   "Name:HadWH_SIG_ghz1_1_ghza2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghzgs2=1,0 Options:SubtractP=HadWH_SIG_ghz1_1_JHUGen_JECNominal,HadWH_SIG_ghza2_1_JHUGen_JECNominal",
   #"Name:HadWH_SIG_ghz1_1_ghza2_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghzgs2=0,1 Options:SubtractP=HadWH_SIG_ghz1_1_JHUGen_JECNominal,HadWH_SIG_ghza2_1_JHUGen_JECNominal",
   "Name:HadWH_SIG_ghz1_1_ghza4_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghzgs4=1,0 Options:SubtractP=HadWH_SIG_ghz1_1_JHUGen_JECNominal,HadWH_SIG_ghza4_1_JHUGen_JECNominal",
   #"Name:HadWH_SIG_ghz1_1_ghza4_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghzgs4=0,1 Options:SubtractP=HadWH_SIG_ghz1_1_JHUGen_JECNominal,HadWH_SIG_ghza4_1_JHUGen_JECNominal",

   "Name:HadWH_SIG_ghz1_1_gha2_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghgsgs2=1,0 Options:SubtractP=HadWH_SIG_ghz1_1_JHUGen_JECNominal,HadWH_SIG_gha2_1_JHUGen_JECNominal",
   #"Name:HadWH_SIG_ghz1_1_gha2_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghgsgs2=0,1 Options:SubtractP=HadWH_SIG_ghz1_1_JHUGen_JECNominal,HadWH_SIG_gha2_1_JHUGen_JECNominal",
   "Name:HadWH_SIG_ghz1_1_gha4_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghgsgs4=1,0 Options:SubtractP=HadWH_SIG_ghz1_1_JHUGen_JECNominal,HadWH_SIG_gha4_1_JHUGen_JECNominal",
   #"Name:HadWH_SIG_ghz1_1_gha4_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:ghz1=1,0;ghgsgs4=0,1 Options:SubtractP=HadWH_SIG_ghz1_1_JHUGen_JECNominal,HadWH_SIG_gha4_1_JHUGen_JECNominal",


   # ttH: Undecayed MEs belong to J2-class clusters
   "Name:ttHUndecayed_SIG_kappa_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:ttH MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:ttHUndecayed_SIG_kappatilde_1_JHUGen_JECNominal Alias:<Name> Process:H0minus Production:ttH MatrixElement:JHUGen Cluster=J2JECNominal",
   "Name:ttHUndecayed_SIG_kappa_1_kappatilde_1_JHUGen_JECNominal Process:SelfDefine_spin0 Production:ttH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:kappa=1,0;kappa_tilde=1,0 Options:SubtractP=ttH_SIG_kappa_1_JHUGen_JECNominal,ttH_SIG_kappatilde_1_JHUGen_JECNominal",
   #"Name:ttHUndecayed_SIG_kappa_1_kappatilde_1_pi2_JHUGen_JECNominal Process:SelfDefine_spin0 Production:ttH MatrixElement:JHUGen Cluster=J2JECNominal Couplings:kappa=1,0;kappa_tilde=0,1 Options:SubtractP=ttH_SIG_kappa_1_JHUGen_JECNominal,ttH_SIG_kappatilde_1_JHUGen_JECNominal",

   # bbH
   # Not adding kappa_tilde since it was shown bbH has close to no sensitivity sny time soon
   "Name:bbH_SIG_kappa_1_JHUGen_JECNominal Process:HSMHiggs Production:bbH MatrixElement:JHUGen Cluster=J2JECNominal",
]
AJetsProdProbabilities_SpinZero_JHUGen_JECUp = [theME.replace("JECNominal", "JECUp") for theME in AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
AJetsProdProbabilities_SpinZero_JHUGen_JECDn = [theME.replace("JECNominal", "JECDn") for theME in AJetsProdProbabilities_SpinZero_JHUGen_JECNominal]
## Production probabilities with >=1 lepton(s) ##
ALepsProdProbabilities_SpinZero_JHUGen = [
   # Leptonic ZH
   "Name:LepZH_SIG_ghz1_1_JHUGen Alias:<Name> Process:HSMHiggs Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH",
   "Name:LepZH_SIG_ghz1prime2_1_JHUGen Alias:<Name> Process:H0_g1prime2 Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH",
   "Name:LepZH_SIG_ghz2_1_JHUGen Alias:<Name> Process:H0hplus Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH",
   "Name:LepZH_SIG_ghz4_1_JHUGen Alias:<Name> Process:H0minus Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH",
   "Name:LepZH_SIG_ghza1prime2_1_JHUGen Alias:<Name> Process:H0_Zgsg1prime2 Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH",
   "Name:LepZH_SIG_ghza2_1_JHUGen Alias:<Name> Process:H0_Zgs Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH",
   "Name:LepZH_SIG_ghza4_1_JHUGen Alias:<Name> Process:H0_Zgs_PS Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH",
   "Name:LepZH_SIG_gha2_1_JHUGen Alias:<Name> Process:H0_gsgs Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH",
   "Name:LepZH_SIG_gha4_1_JHUGen Alias:<Name> Process:H0_gsgs_PS Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH",

   "Name:LepZH_SIG_ghz1_1_ghz1prime2_1_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH Couplings:ghz1=1,0;ghz1_prime2=1,0 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghz1prime2_1_JHUGen",
   #"Name:LepZH_SIG_ghz1_1_ghz1prime2_1_pi2_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH Couplings:ghz1=1,0;ghz1_prime2=0,1 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghz1prime2_1_JHUGen",
   "Name:LepZH_SIG_ghz1_1_ghz2_1_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH Couplings:ghz1=1,0;ghz2=1,0 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghz2_1_JHUGen",
   #"Name:LepZH_SIG_ghz1_1_ghz2_1_pi2_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH Couplings:ghz1=1,0;ghz2=0,1 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghz2_1_JHUGen",
   "Name:LepZH_SIG_ghz1_1_ghz4_1_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH Couplings:ghz1=1,0;ghz4=1,0 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghz4_1_JHUGen",
   #"Name:LepZH_SIG_ghz1_1_ghz4_1_pi2_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH Couplings:ghz1=1,0;ghz4=0,1 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghz4_1_JHUGen",

   "Name:LepZH_SIG_ghz1_1_ghza1prime2_1_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH Couplings:ghz1=1,0;ghzgs1_prime2=1,0 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghza1prime2_1_JHUGen",
   #"Name:LepZH_SIG_ghz1_1_ghza1prime2_1_pi2_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH Couplings:ghz1=1,0;ghzgs1_prime2=0,1 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghza1prime2_1_JHUGen",
   "Name:LepZH_SIG_ghz1_1_ghza2_1_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH Couplings:ghz1=1,0;ghzgs2=1,0 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghza2_1_JHUGen",
   #"Name:LepZH_SIG_ghz1_1_ghza2_1_pi2_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH Couplings:ghz1=1,0;ghzgs2=0,1 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghza2_1_JHUGen",
   "Name:LepZH_SIG_ghz1_1_ghza4_1_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH Couplings:ghz1=1,0;ghzgs4=1,0 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghza4_1_JHUGen",
   #"Name:LepZH_SIG_ghz1_1_ghza4_1_pi2_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH Couplings:ghz1=1,0;ghzgs4=0,1 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_ghza4_1_JHUGen",

   "Name:LepZH_SIG_ghz1_1_gha2_1_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH Couplings:ghz1=1,0;ghgsgs2=1,0 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_gha2_1_JHUGen",
   #"Name:LepZH_SIG_ghz1_1_gha2_1_pi2_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH Couplings:ghz1=1,0;ghgsgs2=0,1 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_gha2_1_JHUGen",
   "Name:LepZH_SIG_ghz1_1_gha4_1_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH Couplings:ghz1=1,0;ghgsgs4=1,0 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_gha4_1_JHUGen",
   #"Name:LepZH_SIG_ghz1_1_gha4_1_pi2_JHUGen Process:SelfDefine_spin0 Production:Lep_ZH MatrixElement:JHUGen Cluster=LepZH Couplings:ghz1=1,0;ghgsgs4=0,1 Options:SubtractP=LepZH_SIG_ghz1_1_JHUGen,LepZH_SIG_gha4_1_JHUGen",

   # Leptonic WH (CAUTION: All requiring the SM ME to be maximized)
   "Name:LepWH_SIG_ghz1_1_JHUGen Alias:<Name> Process:HSMHiggs Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Options:MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
   "Name:LepWH_SIG_ghz1prime2_1_JHUGen Alias:<Name> Process:H0_g1prime2 Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Options:MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
   "Name:LepWH_SIG_ghz2_1_JHUGen Alias:<Name> Process:H0hplus Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Options:MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
   "Name:LepWH_SIG_ghz4_1_JHUGen Alias:<Name> Process:H0minus Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Options:MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
   "Name:LepWH_SIG_ghza1prime2_1_JHUGen Alias:<Name> Process:H0_Zgsg1prime2 Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Options:MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
   "Name:LepWH_SIG_ghza2_1_JHUGen Alias:<Name> Process:H0_Zgs Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Options:MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
   "Name:LepWH_SIG_ghza4_1_JHUGen Alias:<Name> Process:H0_Zgs_PS Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Options:MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
   "Name:LepWH_SIG_gha2_1_JHUGen Alias:<Name> Process:H0_gsgs Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Options:MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
   "Name:LepWH_SIG_gha4_1_JHUGen Alias:<Name> Process:H0_gsgs_PS Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Options:MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",

   "Name:LepWH_SIG_ghz1_1_ghz1prime2_1_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Couplings:ghz1=1,0;ghz1_prime2=1,0 Options:SubtractP=LepWH_SIG_ghz1_1_JHUGen,LepWH_SIG_ghz1prime2_1_JHUGen;MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
   #"Name:LepWH_SIG_ghz1_1_ghz1prime2_1_pi2_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Couplings:ghz1=1,0;ghz1_prime2=0,1 Options:SubtractP=LepWH_SIG_ghz1_1_JHUGen,LepWH_SIG_ghz1prime2_1_JHUGen;MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
   "Name:LepWH_SIG_ghz1_1_ghz2_1_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Couplings:ghz1=1,0;ghz2=1,0 Options:SubtractP=LepWH_SIG_ghz1_1_JHUGen,LepWH_SIG_ghz2_1_JHUGen;MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
   #"Name:LepWH_SIG_ghz1_1_ghz2_1_pi2_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Couplings:ghz1=1,0;ghz2=0,1 Options:SubtractP=LepWH_SIG_ghz1_1_JHUGen,LepWH_SIG_ghz2_1_JHUGen;MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
   "Name:LepWH_SIG_ghz1_1_ghz4_1_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Couplings:ghz1=1,0;ghz4=1,0 Options:SubtractP=LepWH_SIG_ghz1_1_JHUGen,LepWH_SIG_ghz4_1_JHUGen;MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
   #"Name:LepWH_SIG_ghz1_1_ghz4_1_pi2_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Couplings:ghz1=1,0;ghz4=0,1 Options:SubtractP=LepWH_SIG_ghz1_1_JHUGen,LepWH_SIG_ghz4_1_JHUGen;MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",

   "Name:LepWH_SIG_ghz1_1_ghza1prime2_1_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Couplings:ghz1=1,0;ghzgs1_prime2=1,0 Options:SubtractP=LepWH_SIG_ghz1_1_JHUGen,LepWH_SIG_ghza1prime2_1_JHUGen;MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
   #"Name:LepWH_SIG_ghz1_1_ghza1prime2_1_pi2_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Couplings:ghz1=1,0;ghzgs1_prime2=0,1 Options:SubtractP=LepWH_SIG_ghz1_1_JHUGen,LepWH_SIG_ghza1prime2_1_JHUGen;MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
   "Name:LepWH_SIG_ghz1_1_ghza2_1_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Couplings:ghz1=1,0;ghzgs2=1,0 Options:SubtractP=LepWH_SIG_ghz1_1_JHUGen,LepWH_SIG_ghza2_1_JHUGen;MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
   #"Name:LepWH_SIG_ghz1_1_ghza2_1_pi2_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Couplings:ghz1=1,0;ghzgs2=0,1 Options:SubtractP=LepWH_SIG_ghz1_1_JHUGen,LepWH_SIG_ghza2_1_JHUGen;MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
   "Name:LepWH_SIG_ghz1_1_ghza4_1_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Couplings:ghz1=1,0;ghzgs4=1,0 Options:SubtractP=LepWH_SIG_ghz1_1_JHUGen,LepWH_SIG_ghza4_1_JHUGen;MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
   #"Name:LepWH_SIG_ghz1_1_ghza4_1_pi2_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Couplings:ghz1=1,0;ghzgs4=0,1 Options:SubtractP=LepWH_SIG_ghz1_1_JHUGen,LepWH_SIG_ghza4_1_JHUGen;MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",

   "Name:LepWH_SIG_ghz1_1_gha2_1_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Couplings:ghz1=1,0;ghgsgs2=1,0 Options:SubtractP=LepWH_SIG_ghz1_1_JHUGen,LepWH_SIG_gha2_1_JHUGen;MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
   #"Name:LepWH_SIG_ghz1_1_gha2_1_pi2_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Couplings:ghz1=1,0;ghgsgs2=0,1 Options:SubtractP=LepWH_SIG_ghz1_1_JHUGen,LepWH_SIG_gha2_1_JHUGen;MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
   "Name:LepWH_SIG_ghz1_1_gha4_1_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Couplings:ghz1=1,0;ghgsgs4=1,0 Options:SubtractP=LepWH_SIG_ghz1_1_JHUGen,LepWH_SIG_gha4_1_JHUGen;MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
   #"Name:LepWH_SIG_ghz1_1_gha4_1_pi2_JHUGen Process:SelfDefine_spin0 Production:Lep_WH MatrixElement:JHUGen Cluster=LepWH Couplings:ghz1=1,0;ghgsgs4=0,1 Options:SubtractP=LepWH_SIG_ghz1_1_JHUGen,LepWH_SIG_gha4_1_JHUGen;MaxNumerator=LepWH_SIG_ghz1_1_JHUGen",
]

# Construct the final list
theRecoProbabilities = []
theRecoProbabilities.extend(DecayProbabilities_SpinZero_JHUGen)
theRecoProbabilities.extend(AJetsProdProbabilities_SpinZero_JHUGen_JECNominal)
theRecoProbabilities.extend(AJetsProdProbabilities_SpinZero_JHUGen_JECUp)
theRecoProbabilities.extend(AJetsProdProbabilities_SpinZero_JHUGen_JECDn)
theRecoProbabilities.extend(ALepsProdProbabilities_SpinZero_JHUGen)

# Append final list
process.ZZCand.recoProbabilities.extend(theRecoProbabilities)
process.ZZTree.recoProbabilities.extend(theRecoProbabilities)
process.CRZLLTree.recoProbabilities.extend(theRecoProbabilities)
process.ZZTreelooseEle.recoProbabilities.extend(theRecoProbabilities)
process.CRZLLTreelooseEle.recoProbabilities.extend(theRecoProbabilities)
process.CRZLLTreeZ1RSE.recoProbabilities.extend(theRecoProbabilities)
process.ZZTreetle.recoProbabilities.extend(theRecoProbabilities)
process.CRZLLTreetle.recoProbabilities.extend(theRecoProbabilities)
