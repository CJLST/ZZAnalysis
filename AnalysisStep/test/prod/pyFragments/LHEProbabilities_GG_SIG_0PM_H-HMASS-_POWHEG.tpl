LHE_PropagatorRewgt = [
   "Name:SamplePropagator Alias:<Name> PropScheme:CPS hmass:<HMASS> isGen:1 NoBranch:1 isProp:1",
   "Name:CPStoBWPropRewgt PropScheme:FixedWidth hmass:<HMASS> Options:DivideP=SamplePropagator isGen:1 isProp:1",
]
LHE_DecayProbabilities_MCFM = [
   "Name:SampleHypothesisMCFM Alias:<Name> Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=1,0 Options:DivideP=SampleHypothesisMCFM hmass:<HMASS> Cluster:NoInitialQ isGen:1 NoBranch:1",

   "Name:GG_SIG_kappaTopBot_1_ghz1_1_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_kappaTopBot_1_ghz1prime2_1E4_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1_prime2=10000,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_kappaTopBot_1_ghz2_1_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz2=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_kappaTopBot_1_ghz4_1_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz4=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_kappaTopBot_1_ghza1prime2_1E4_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghzgs1_prime2=10000,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_kappaTopBot_1_ghza2_1_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghzgs2=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_kappaTopBot_1_ghza4_1_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghzgs4=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_kappaTopBot_1_gha2_1_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghgsgs2=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_kappaTopBot_1_gha4_1_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghgsgs4=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",

   "Name:GG_SIG_kappaTopBot_1_ghz1_1_ghz1prime2_1E4_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=1,0;ghz1_prime2=10000,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   #"Name:GG_SIG_kappaTopBot_1_ghz1_1_ghz1prime2_1E4i_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=1,0;ghz1_prime2=0,10000 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_kappaTopBot_1_ghz1_1_ghz2_1_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=1,0;ghz2=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   #"Name:GG_SIG_kappaTopBot_1_ghz1_1_ghz2_i_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=1,0;ghz2=0,1 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_kappaTopBot_1_ghz1_1_ghz4_1_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=1,0;ghz4=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   #"Name:GG_SIG_kappaTopBot_1_ghz1_1_ghz4_i_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=1,0;ghz4=0,1 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_kappaTopBot_1_ghz1_1_ghza1prime2_1E4_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=1,0;ghzgs1_prime2=10000,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   #"Name:GG_SIG_kappaTopBot_1_ghz1_1_ghza1prime2_1E4i_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=1,0;ghzgs1_prime2=0,10000 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_kappaTopBot_1_ghz1_1_ghza2_1_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=1,0;ghzgs2=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   #"Name:GG_SIG_kappaTopBot_1_ghz1_1_ghza2_i_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=1,0;ghzgs2=0,1 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_kappaTopBot_1_ghz1_1_ghza4_1_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=1,0;ghzgs4=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   #"Name:GG_SIG_kappaTopBot_1_ghz1_1_ghza4_i_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=1,0;ghzgs4=0,1 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_kappaTopBot_1_ghz1_1_gha2_1_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=1,0;ghgsgs2=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   #"Name:GG_SIG_kappaTopBot_1_ghz1_1_gha2_i_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=1,0;ghgsgs2=0,1 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_kappaTopBot_1_ghz1_1_gha4_1_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=1,0;ghgsgs4=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   #"Name:GG_SIG_kappaTopBot_1_ghz1_1_gha4_i_MCFM Process:HSMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=1,0;ghgsgs4=0,1 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",

   "Name:GG_BSI_kappaTopBot_1_ghz1_1_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   #"Name:GG_BSI_kappaTopBot_1_ghz1_i_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1=0,1 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_BSI_kappaTopBot_1_ghz1prime2_1E4_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1_prime2=10000,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   #"Name:GG_BSI_kappaTopBot_1_ghz1prime2_1E4i_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz1_prime2=0,10000 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_BSI_kappaTopBot_1_ghz2_1_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz2=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   #"Name:GG_BSI_kappaTopBot_1_ghz2_i_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz2=0,1 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_BSI_kappaTopBot_1_ghz4_1_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz4=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   #"Name:GG_BSI_kappaTopBot_1_ghz4_i_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghz4=0,1 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_BSI_kappaTopBot_1_ghza1prime2_1E4_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghzgs1_prime2=10000,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   #"Name:GG_BSI_kappaTopBot_1_ghza1prime2_1E4i_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghzgs1_prime2=0,10000 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_BSI_kappaTopBot_1_ghza2_1_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghzgs2=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   #"Name:GG_BSI_kappaTopBot_1_ghza2_i_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghzgs2=0,1 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_BSI_kappaTopBot_1_ghza4_1_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghzgs4=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   #"Name:GG_BSI_kappaTopBot_1_ghza4_i_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghzgs4=0,1 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_BSI_kappaTopBot_1_gha2_1_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghgsgs2=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   #"Name:GG_BSI_kappaTopBot_1_gha2_i_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghgsgs2=0,1 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   "Name:GG_BSI_kappaTopBot_1_gha4_1_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghgsgs4=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",
   #"Name:GG_BSI_kappaTopBot_1_gha4_i_MCFM Process:bkgZZ_SMHiggs Production:ZZGG MatrixElement:MCFM Couplings:kappa_top=1,0;kappa_bot=1,0;ghgsgs4=0,1 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:NoInitialQ isGen:1",

   "Name:GG_BKG_MCFM Process:bkgZZ Production:ZZGG MatrixElement:MCFM Options:DivideP=SampleHypothesisMCFM Cluster:NoInitialQ isGen:1",
   "Name:QQB_BKG_MCFM Process:bkgZZ Production:ZZQQB MatrixElement:MCFM Options:DivideP=SampleHypothesisMCFM Cluster:NoInitialG isGen:1",
]
### Spin-0 decay probabilities from JHUGen ###
LHE_DecayProbabilities_SpinZero_JHUGen = [
   # The corresponding MCFM branch is present, and one can still compute the VA reweighting as (a1/bkg)_MCFM * (VA/a1)_JHUGen for the time being. No need for DivideP when NoBranch:1
   "Name:GG_SIG_ghg2_1_ghz1_1_JHUGen Alias:SampleHypothesisJHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0 hmass:<HMASS> Cluster:NoInitialQ isGen:1 NoBranch:1",

#   "Name:GG_SIG_ghg2_1_ghz1_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1_prime2=10000,0 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_ghz2_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz2=1,0 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_ghz4_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz4=1,0 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghzgs1_prime2=10000,0 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_ghza2_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghzgs2=1,0 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_ghza4_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghzgs4=1,0 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_gha2_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghgsgs2=1,0 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_gha4_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghgsgs4=1,0 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",

#   "Name:GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghz1_prime2=10000,0 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4i_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghz1_prime2=0,10000 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghz2=1,0 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_ghz1_1_ghz2_i_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghz2=0,1 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghz4=1,0 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_ghz1_1_ghz4_i_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghz4=0,1 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghzgs1_prime2=10000,0 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4i_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghzgs1_prime2=0,10000 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_ghz1_1_ghza2_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghzgs2=1,0 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_ghz1_1_ghza2_i_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghzgs2=0,1 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_ghz1_1_ghza4_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghzgs4=1,0 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_ghz1_1_ghza4_i_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghzgs4=0,1 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_ghz1_1_gha2_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghgsgs2=1,0 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_ghz1_1_gha2_i_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghgsgs2=0,1 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_ghz1_1_gha4_1_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghgsgs4=1,0 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
#   "Name:GG_SIG_ghg2_1_ghz1_1_gha4_i_JHUGen Process:SelfDefine_spin0 Production:ZZGG MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0;ghgsgs4=0,1 Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
]
### Spin-2 decay probabilities from JHUGen ###
LHE_Probabilities_SpinTwo_JHUGen = [
   "Name:GG_SIG_gXg1_1_gXz1_1_JHUGen Process:H2_g1 Production:ZZGG MatrixElement:JHUGen Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_gXg2_1_gXz2_1_JHUGen Process:H2_g2 Production:ZZGG MatrixElement:JHUGen Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_gXg3_1_gXz3_1_JHUGen Process:H2_g3 Production:ZZGG MatrixElement:JHUGen Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_gXg4_1_gXz4_1_JHUGen Process:H2_g4 Production:ZZGG MatrixElement:JHUGen Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_gXg1_1_gXz5_1_JHUGen Process:H2_g5 Production:ZZGG MatrixElement:JHUGen Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_gXg1_1_gXz1_1_gXz5_1_JHUGen Process:H2_g1g5 Production:ZZGG MatrixElement:JHUGen Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_gXg1_1_gXz6_1_JHUGen Process:H2_g6 Production:ZZGG MatrixElement:JHUGen Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_gXg1_1_gXz7_1_JHUGen Process:H2_g7 Production:ZZGG MatrixElement:JHUGen Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_gXg5_1_gXz8_1_JHUGen Process:H2_g8 Production:ZZGG MatrixElement:JHUGen Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_gXg5_1_gXz9_1_JHUGen Process:H2_g9 Production:ZZGG MatrixElement:JHUGen Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
   "Name:GG_SIG_gXg5_1_gXz10_1_JHUGen Process:H2_g10 Production:ZZGG MatrixElement:JHUGen Options:DivideP=SampleHypothesisJHUGen hmass:<HMASS> Cluster:NoInitialQ isGen:1",
]

# Construct the final list
theLHEProbabilities = []
theLHEProbabilities.extend(LHE_PropagatorRewgt)
theLHEProbabilities.extend(LHE_DecayProbabilities_MCFM)
theLHEProbabilities.extend(LHE_DecayProbabilities_SpinZero_JHUGen)
theLHEProbabilities.extend(LHE_Probabilities_SpinTwo_JHUGen)

# Append final list
for name in (
             "ZZTree",
             "CRZLLTree",
             "ZZTreelooseEle",
             "CRZLLTreelooseEle",
             "CRZLLTreeZ1RSE",
             "ZZTreetle",
             "CRZLLTreetle",
            ):
    if hasattr(process, name):
        tree = getattr(process, name)
        #turn on failedTree keeping the most relevant information
        tree.lheProbabilities.extend(theLHEProbabilities)
        if name == "ZZTree" and tree.skipEmptyEvents:
            tree.failedTreeLevel = max(tree.failedTreeLevel.value(), LHEFailedTree)

