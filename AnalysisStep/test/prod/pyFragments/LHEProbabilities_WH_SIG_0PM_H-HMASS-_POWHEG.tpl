LHE_PropagatorRewgt = [
   "Name:SamplePropagator Alias:<Name> PropScheme:CPS hmass:<HMASS> isGen:1 NoBranch:1 isProp:1",
   "Name:CPStoBWPropRewgt PropScheme:FixedWidth hmass:<HMASS> Options:DivideP=SamplePropagator isGen:1 isProp:1",
]
LHE_Probabilities_MCFM = [
   "Name:SampleHypothesisMCFM Alias:<Name> Process:HSMHiggs Production:Had_WH_S MatrixElement:MCFM Couplings:ghz1=1,0 hmass:<HMASS> Cluster:BestNLOWHApproximation isGen:1 NoBranch:1",

   "Name:JJEW_SIG_ghv1_1_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1prime2_1E4_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1_prime2=10000,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv2_1_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz2=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv4_1_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz4=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghza1prime2_1E4_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghzgs1_prime2=10000,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghza2_1_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghzgs2=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghza4_1_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghzgs4=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_gha2_1_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghgsgs2=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_gha4_1_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghgsgs4=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",

   "Name:JJEW_SIG_ghv1_1_ghv1prime2_25E2_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz1_prime2=2500,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghv1prime2_25E2i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz1_prime2=0,2500 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv2_0p251_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz2=0.25,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghv2_0p25i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz2=0,0.25 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv4_0p25_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz4=0.25,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghv4_0p25i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz4=0,0.25 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza1prime2_25E2_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs1_prime2=2500,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghza1prime2_25E2i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs1_prime2=0,2500 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza2_0p25_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs2=0.25,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghza2_0p25i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs2=0,0.25 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza4_0p25_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs4=0.25,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghza4_0p25i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs4=0,0.25 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1_1_gha2_0p25_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs2=0.25,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_gha2_0p25i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs2=0,0.25 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1_1_gha4_0p25_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs4=0.25,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_gha4_0p25i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs4=0,0.25 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",

   "Name:JJEW_SIG_ghv1_1_ghv1prime2_50E2_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz1_prime2=5000,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghv1prime2_50E2i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz1_prime2=0,5000 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghv1prime2_50E2_50E2i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz1_prime2=5000,5000 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv2_0p5_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz2=0.5,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghv2_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz2=0,0.5 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghv2_0p5_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz2=0.5,0.5 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv4_0p5_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz4=0.5,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghv4_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz4=0,0.5 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghv4_0p5_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz4=0.5,0.5 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza1prime2_50E2_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs1_prime2=5000,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghza1prime2_50E2i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs1_prime2=0,5000 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghza1prime2_50E2_50E2i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs1_prime2=5000,5000 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza2_0p5_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs2=0.5,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghza2_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs2=0,0.5 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghza2_0p5_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs2=0.5,0.5 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza4_0p5_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs4=0.5,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghza4_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs4=0,0.5 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghza4_0p5_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs4=0.5,0.5 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1_1_gha2_0p5_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs2=0.5,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_gha2_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs2=0,0.5 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_gha2_0p5_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs2=0.5,0.5 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1_1_gha4_0p5_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs4=0.5,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_gha4_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs4=0,0.5 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_gha4_0p5_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs4=0.5,0.5 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",

   "Name:JJEW_SIG_ghv1_1_ghv1prime2_75E2_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz1_prime2=7500,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghv1prime2_75E2i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz1_prime2=0,7500 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv2_0p75_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz2=0.75,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghv2_0p75i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz2=0,0.75 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv4_0p75_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz4=0.75,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghv4_0p75i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz4=0,0.75 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza1prime2_75E2_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs1_prime2=7500,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghza1prime2_75E2i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs1_prime2=0,7500 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza2_0p75_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs2=0.75,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghza2_0p75i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs2=0,0.75 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza4_0p75_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs4=0.75,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_ghza4_0p75i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs4=0,0.75 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1_1_gha2_0p75_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs2=0.75,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_gha2_0p75i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs2=0,0.75 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_SIG_ghv1_1_gha4_0p75_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs4=0.75,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_SIG_ghv1_1_gha4_0p75i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs4=0,0.75 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",

   "Name:JJEW_BSI_ghv1_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_BSI_ghv1prime2_1E4_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1_prime2=10000,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_BSI_ghv2_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz2=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_BSI_ghv4_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz4=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_BSI_ghza1prime2_1E4_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghzgs1_prime2=10000,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_BSI_ghza2_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghzgs2=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_BSI_ghza4_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghzgs4=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_BSI_gha2_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghgsgs2=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_BSI_gha4_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghgsgs4=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",

   "Name:JJEW_BSI_ghv1_0p5_0p5i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=0.5,0.5 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_BSI_ghv1prime2_50E2_50E2i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1_prime2=5000,5000 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_BSI_ghv2_0p5_0p5i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz2=0.5,0.5 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_BSI_ghv4_0p5_0p5i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz4=0.5,0.5 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_BSI_ghza1prime2_50E2_50E2i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghzgs1_prime2=5000,5000 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_BSI_ghza2_0p5_0p5i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghzgs2=0.5,0.5 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_BSI_ghza4_0p5_0p5i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghzgs4=0.5,0.5 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_BSI_gha2_0p5_0p5i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghgsgs2=0.5,0.5 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_BSI_gha4_0p5_0p5i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghgsgs4=0.5,0.5 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",

   "Name:JJEW_BSI_ghv1_1_ghv1prime2_1E4_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz1_prime2=10000,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_BSI_ghv1_1_ghv1prime2_1E4i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz1_prime2=0,10000 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_BSI_ghv1_1_ghv2_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz2=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_BSI_ghv1_1_ghv2_i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz2=0,1 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_BSI_ghv1_1_ghv4_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz4=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_BSI_ghv1_1_ghv4_i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz4=0,1 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_BSI_ghv1_1_ghza1prime2_1E4_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs1_prime2=10000,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_BSI_ghv1_1_ghza1prime2_1E4i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs1_prime2=0,10000 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_BSI_ghv1_1_ghza2_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs2=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_BSI_ghv1_1_ghza2_i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs2=0,1 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_BSI_ghv1_1_ghza4_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs4=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_BSI_ghv1_1_ghza4_i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs4=0,1 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_BSI_ghv1_1_gha2_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs2=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_BSI_ghv1_1_gha2_i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs2=0,1 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJEW_BSI_ghv1_1_gha4_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs4=1,0 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",
   #"Name:JJEW_BSI_ghv1_1_gha4_i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs4=0,1 Options:DivideP=SampleHypothesisMCFM hmass:125 hwidth:0.00407 Cluster:BestNLOWHApproximation isGen:1",

   "Name:JJEW_BKG_MCFM Process:bkgZZ Production:JJEW MatrixElement:MCFM Options:DivideP=SampleHypothesisMCFM Cluster:BestNLOWHApproximation isGen:1",
   "Name:JJQCD_BKG_MCFM Process:bkgZZ Production:JJQCD MatrixElement:MCFM Options:DivideP=SampleHypothesisMCFM Cluster:CommonLast isGen:1",
]
### Spin-0 decay probabilities from JHUGen ###
LHE_DecayProbabilities_SpinZero_JHUGen = [
   "Name:Dec_SIG_ghz1_1_JHUGen Alias:SampleDecayHypothesisJHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_ghz1prime2_1E4_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1_prime2=10000,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_ghz2_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz2=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_ghz4_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz4=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_ghza1prime2_1E4_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghzgs1_prime2=10000,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_ghza2_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghzgs2=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_ghza4_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghzgs4=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_gha2_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghgsgs2=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_gha4_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghgsgs4=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",

   "Name:Dec_SIG_ghz1_1_ghz1prime2_1E4_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghz1_prime2=10000,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_ghz1_1_ghz1prime2_1E4i_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghz1_prime2=0,10000 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_ghz1_1_ghz2_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghz2=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_ghz1_1_ghz2_i_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghz2=0,1 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_ghz1_1_ghz4_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghz4=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_ghz1_1_ghz4_i_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghz4=0,1 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_ghz1_1_ghza1prime2_1E4_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghzgs1_prime2=10000,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_ghz1_1_ghza1prime2_1E4i_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghzgs1_prime2=0,10000 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_ghz1_1_ghza2_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghzgs2=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_ghz1_1_ghza2_i_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghzgs2=0,1 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_ghz1_1_ghza4_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghzgs4=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_ghz1_1_ghza4_i_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghzgs4=0,1 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_ghz1_1_gha2_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghgsgs2=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_ghz1_1_gha2_i_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghgsgs2=0,1 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_ghz1_1_gha4_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghgsgs4=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
   "Name:Dec_SIG_ghz1_1_gha4_i_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghgsgs4=0,1 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
]
LHE_ProdProbabilities_SpinZero_JHUGen = [
   "Name:WH_SIG_ghw1_1_JHUGen Alias:SampleProductionHypothesisJHUGen Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Couplings:ghz1=1,0 Options:DivideP=SampleProductionHypothesisJHUGen hmass:<HMASS> Cluster:BestNLOWHApproximation isGen:1",
   "Name:WH_SIG_ghw1prime2_1E4_JHUGen Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Couplings:ghz1_prime2=10000,0 Options:DivideP=SampleProductionHypothesisJHUGen hmass:<HMASS> Cluster:BestNLOWHApproximation isGen:1",
   "Name:WH_SIG_ghw2_1_JHUGen Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Couplings:ghz2=1,0 Options:DivideP=SampleProductionHypothesisJHUGen hmass:<HMASS> Cluster:BestNLOWHApproximation isGen:1",
   "Name:WH_SIG_ghw4_1_JHUGen Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Couplings:ghz4=1,0 Options:DivideP=SampleProductionHypothesisJHUGen hmass:<HMASS> Cluster:BestNLOWHApproximation isGen:1",

   "Name:WH_SIG_ghw1_1_ghw1prime2_1E4_JHUGen Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Couplings:ghz1=1,0;ghz1_prime2=10000,0 Options:DivideP=SampleProductionHypothesisJHUGen hmass:<HMASS> Cluster:BestNLOWHApproximation isGen:1",
   "Name:WH_SIG_ghw1_1_ghw1prime2_1E4i_JHUGen Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Couplings:ghz1=1,0;ghz1_prime2=0,10000 Options:DivideP=SampleProductionHypothesisJHUGen hmass:<HMASS> Cluster:BestNLOWHApproximation isGen:1",
   "Name:WH_SIG_ghw1_1_ghw2_1_JHUGen Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Couplings:ghz1=1,0;ghz2=1,0 Options:DivideP=SampleProductionHypothesisJHUGen hmass:<HMASS> Cluster:BestNLOWHApproximation isGen:1",
   "Name:WH_SIG_ghw1_1_ghw2_i_JHUGen Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Couplings:ghz1=1,0;ghz2=0,1 Options:DivideP=SampleProductionHypothesisJHUGen hmass:<HMASS> Cluster:BestNLOWHApproximation isGen:1",
   "Name:WH_SIG_ghw1_1_ghw4_1_JHUGen Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Couplings:ghz1=1,0;ghz4=1,0 Options:DivideP=SampleProductionHypothesisJHUGen hmass:<HMASS> Cluster:BestNLOWHApproximation isGen:1",
   "Name:WH_SIG_ghw1_1_ghw4_i_JHUGen Process:SelfDefine_spin0 Production:Had_WH MatrixElement:JHUGen Couplings:ghz1=1,0;ghz4=0,1 Options:DivideP=SampleProductionHypothesisJHUGen hmass:<HMASS> Cluster:BestNLOWHApproximation isGen:1",
]

# Construct the final list
theLHEProbabilities = []
theLHEProbabilities.extend(LHE_PropagatorRewgt)
theLHEProbabilities.extend(LHE_Probabilities_MCFM)
#theLHEProbabilities.extend(LHE_DecayProbabilities_SpinZero_JHUGen)
#theLHEProbabilities.extend(LHE_ProdProbabilities_SpinZero_JHUGen)

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
