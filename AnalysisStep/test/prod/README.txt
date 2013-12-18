INTRODUCTION

The selection is driven by the master .py file in
ZZAnalysis/AnalysisStep/test/MasterPy/ZZ4lAnalysis.py. 
This controls bulding of objects with combinatorics, cuts, selection,
FSR recovery etc. This .py includes only the modules to prepare
objects with attached variables and cut information. Additional
analysis modules can be added separately to, for example, fill plots,
print debug info, and write trees.

The modules that we currently use for analysis are under ZZAnalysis/AnalysisStep/test/Ntuplizers.
As an example, a python file to run the synchronization excercise of:
https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsZZ4l2012SummerSync#Table_Summer12_NEW_ISO_including
is present in:
ZZAnalysis/AnalysisStep/analyzer_May21Sync_Summer12.py
This runs an analyzer that directly fills histograms
(test/Ntuplizers/ZZ4lAnalyzer.cc), a tree maker
(test/Ntuplizers/HZZ4lNtupleMaker.h) and a debug module that dumps all
user information attached to objects and candidates (test/dumpUserData.cc).

For the production of trees, we use the same modules, except that a
more complete set of trees and plot is made, including those for
CRs. The standard configuration for the production of these trees is in:
ZZAnalysis/AnalysisStep/test/analyzer.py. 

Here are some more technical instructions to produce these trees.

NTUPLES
Root files contains several TDirectories, in particular one for each final state and CR.
Each folder contains a tree and an histogram that stores normalization information.

The tree format is defined in test/Ntuplizers/HZZ4lNtupleFactory.cc. The tree contains 
per-event float or integer variables as well as vectors to store variables about candidates.
The variable iBC represent the index of the chosen candidate in the event if present, -1 otherwise.
The variable ZZsel stores info about the level of selection of the candidate, conventinally a value 
of 100 or larger indicates a passing event. So for example to draw the mass spectrum for passing events:

candTree->Draw("ZZMass[iBC]","iBC>=0&&ZZsel[iBC]>=100")



CHECKING UP AND COMPILING
Some additional packages are needed to compile.
Checkout instructions for 444 (2011 data/MC) and 525 (2012 data/MC) are at:
http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/CJLST/ZZAnalysis/README.txt?view=markup


RUNNING INTERACTIVELY
Set the proper input files in analyzer.py.
To run on data, set IsMC = False and PD="DoubleEle", "DoubleMu", or "MuEG".
LEPTON_SETUP has to be set to 2011 or 2012 (it controls which triggers
are requested, the effective areas for rho correction, etc.)


JOB SUBMISSION
To run on all data and MC, the processing is splitted in jobs submitted to batch queues.
Submission tools are under test/prod.

A list of samples to be processed is kept in these files:
 test/prod/analyzer_2011.py (2011 data and MC)
 test/prod/analyzer_2012.py (2012 data and MC)

These files are used to generate jobs with the proper configuration, using the script batch.py in the prod directory.
This script works only on lxplus (it needs bsub to be present).

Example: to create jobs for all samples defined in analyzer_2012.py (wihout submitting):
./batch.py -n -o PROD analyzer_2012.py
This generates a subfolder for each job in a directory called PROD. 
(Note that generating all jobs for all samples can take quite some
time: just be patient. for experiments it's better to comment out all samples but one in analyzer_2012.py).

To generate jobs AND SUBMIT them:

./batch.py -o PROD analyzer_2012.py

(To submit in a different queue, add this command line option: -b 'bsub -q 8nh < batchScript.sh' )

The job output (root file and log) is written in each job's folder as soon as the job finishes.

It can happen that a small fraction of jobs fails for technical reasons. It is important to identify these before 
trying to merge root files. For this purpose some simple scripts are available:

cd PROD
../checkProd.csh
This checks each output file and, if no problem is found, moves the job folder to a temporary folder called 
AAAOK, so to separate them from failed jobs. Once this is done, failed jobs can be resubmitted, if needed (be sure that 
no job is running or pending before resubmitting it!). Still from the PROD dir:
../cleanup.csh
../resubmit.csh

The root files for jobs completed succesfully can then be merged with this command:
haddChunks.py AAAOK

This creates one directory per sample containing the merged root file.

