###################################

###If you just need to produce the final ntuples!

Once you checked out the ZZAnalysis package and you cmsenv'ed, in ZZAnalysis/AnalysisStep/test/Macros

$ gmake all

In order to run the code to produce the final ntuples, two scripts are provided:

$ ./executeRun_7TeV.sh <input-dir> <output-dir>
$ ./executeRun_8TeV.sh <input-dir> <output-dir>


These scripts run both data and MC for all the 3 final states, for all MC samples. They also produce the CR region plots too.

Notice that the scripts rely on the fact that you have a directory with the input files called rootuples/<input-dir>/PRODFSR(_8TeV) and an output directory called trees.

Also, the scripts will create the final output directory if it doesn't exist, or they will delete all previous output files if the directory already exists

###################################

###Some more info on the code and the output


The output files will have the final "official" selection implemented. They'll also contain a variable (MC_weight) that provides for MC event-by-event weight that takes into account the proper cross sections, the data/MC efficiency corrections and the PU reweight.

If you want to change something, most of the action happens in HZZ4l.C. This is simple looper over all the events in a given input file. It inherits from the HZZ4lBase class, where all the variables contained in the input trees are datamembers of the class.

If you look inside the Loop method, you can see how the selection is done. The "sel" variable is used to implement the current "official" selection.


###Notes on event weights
MC_weight: such that sum(weights) = yield per fb
defined as (sigma*BR)/Nevt_Gen * PU_w*Powheg_w*DataMC_w
This is used to extract ZZ yields

MC_weight_norm: such that sum(weights) = efficiency (relative to the "proper" gen final state)
defined as: 1/Nevt_Gen_FS * PU_w*Powheg_w*DataMC_w
Nevt_Gen_FS is taken directly from the proper bin of the Counters histogram for H samples, while for ZZ samples it is taken from pre-computed tables in order to mix correctly the different samples.

This is used to extract signal efficiencies relative to the generated events in a given final state.


MC_weight_noxsec: such that  sum(weights) = efficiency (relative to all gen events)
defined as: 1/Nevt_Gen * PU_w*Powheg_w*DataMC_w
It is practically = MC_weight/xsection; 
it differs from MC_weight_norm because the denominator is all generated events and not only those of the "proper" final state.
This is used to fill templates. Can be used to get efficiency for filtered samples (eg VH)
