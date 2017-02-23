LHE_PropagatorRewgt = [
   "Name:SamplePropagator Alias:<Name> PropScheme:CPS hmass:<HMASS> isGen:1 NoBranch:1 isProp:1",
   "Name:CPStoBWPropRewgt PropScheme:FixedWidth hmass:<HMASS> Options:DivideP=SamplePropagator isGen:1 isProp:1",
]
### Spin-0 decay probabilities from JHUGen ###
LHE_DecayProbabilities_SpinZero_JHUGen = [
   "Name:SampleDecayHypothesisJHUGen Alias:<Name> Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0 hmass:<HMASS> isGen:1 NoBranch:1",

#   "Name:Dec_SIG_ghz1_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_ghz1prime2_1E4_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1_prime2=10000,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_ghz2_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz2=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_ghz4_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz4=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_ghza1prime2_1E4_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghzgs1_prime2=10000,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_ghza2_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghzgs2=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_ghza4_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghzgs4=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_gha2_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghgsgs2=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_gha4_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghgsgs4=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",

#   "Name:Dec_SIG_ghz1_1_ghz1prime2_1E4_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghz1_prime2=10000,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_ghz1_1_ghz1prime2_1E4i_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghz1_prime2=0,10000 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_ghz1_1_ghz2_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghz2=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_ghz1_1_ghz2_i_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghz2=0,1 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_ghz1_1_ghz4_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghz4=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_ghz1_1_ghz4_i_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghz4=0,1 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_ghz1_1_ghza1prime2_1E4_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghzgs1_prime2=10000,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_ghz1_1_ghza1prime2_1E4i_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghzgs1_prime2=0,10000 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_ghz1_1_ghza2_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghzgs2=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_ghz1_1_ghza2_i_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghzgs2=0,1 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_ghz1_1_ghza4_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghzgs4=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_ghz1_1_ghza4_i_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghzgs4=0,1 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_ghz1_1_gha2_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghgsgs2=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_ghz1_1_gha2_i_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghgsgs2=0,1 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_ghz1_1_gha4_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghgsgs4=1,0 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
#   "Name:Dec_SIG_ghz1_1_gha4_i_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0;ghgsgs4=0,1 Options:DivideP=SampleDecayHypothesisJHUGen hmass:<HMASS> isGen:1",
]
LHE_ProdProbabilities_SpinZero_JHUGen = [
#   "Name:SampleProductionHypothesisJHUGen Alias:<Name> Process:SelfDefine_spin0 Production:ttH MatrixElement:JHUGen Couplings:kappa=1,0 Cluster:NoAssociatedG hmass:<HMASS> isGen:1 NoBranch:1",

#   "Name:ttH_SIG_kappa_1_JHUGen Alias:SampleProductionHypothesisJHUGen Process:SelfDefine_spin0 Production:ttH MatrixElement:JHUGen Couplings:kappa=1,0 Options:DivideP=SampleProductionHypothesisJHUGen Cluster:NoAssociatedG hmass:<HMASS> isGen:1",
#   "Name:ttH_SIG_kappa_tilde_1_JHUGen Process:SelfDefine_spin0 Production:ttH MatrixElement:JHUGen Couplings:kappa_tilde=1,0 Options:DivideP=SampleProductionHypothesisJHUGen Cluster:NoAssociatedG hmass:<HMASS> isGen:1",
#   "Name:ttH_SIG_kappa_1_kappa_tilde_1_JHUGen Process:SelfDefine_spin0 Production:ttH MatrixElement:JHUGen Couplings:kappa=1,0;kappa_tilde=1,0 Options:DivideP=SampleProductionHypothesisJHUGen Cluster:NoAssociatedG hmass:<HMASS> isGen:1",
#   "Name:ttH_SIG_kappa_1_kappa_tilde_i_JHUGen Process:SelfDefine_spin0 Production:ttH MatrixElement:JHUGen Couplings:kappa=1,0;kappa_tilde=0,1 Options:DivideP=SampleProductionHypothesisJHUGen Cluster:NoAssociatedG hmass:<HMASS> isGen:1",
]

# Construct the final list
theLHEProbabilities = []
theLHEProbabilities.extend(LHE_PropagatorRewgt)
theLHEProbabilities.extend(LHE_DecayProbabilities_SpinZero_JHUGen)
theLHEProbabilities.extend(LHE_ProdProbabilities_SpinZero_JHUGen)

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
