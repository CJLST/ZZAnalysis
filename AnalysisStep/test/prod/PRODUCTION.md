Summary of steps for a full standard production
-----------------------------------------------
-----------------------------------------------

PRIMARY trees
-------------

1. create the jobs:
```
./batch.py -o PT13TeV samples_2015.csv 

```

2. Submit the jobs (from lxplus; does not work from private machines):

```
cd PT13TeV
resubmit.csh
``

Once jobs are done, from the same folder:

```
checkProd.csh
```
This checks all jobs and moves all those which finished correctly to a subfolder named AAAOK.


To resubmit failed jobs (repeat until all jobs succeed):

```
cleanup.csh
resubmit.csh
# wait all jobs to finish
checkProd.csh
```

Notes:
*   Sometimes jobs get stuck in the processing nodes. This can be checked with bjobs -l <jobid> - if job is running since a lot of time but appears to be idle (check e.g. IDLE_FACTOR), it's probably best to kill it; wait that bjobs reports it as finished (can take a while), and resubmit it. 
*   Sometimes failures are permanent (job keeps failing), in this case it must be investigated looking at the log.txt.gz file. In rare cases, this happens to be due to an EOS failure. If a file is reported as unaccessible: 
    *     try opening it in plain root with xrootd; that should give a reproducible error message 
    *     check =eos fileinfo <file>= : what usually the problem is that servers holding both replicas are down, or that one of the two replicas is reported with incorrect information. Open an EOS ticket with a clear pointer to the file and to the eos fileinfo message; it is normally addressed within few hours.


Once all jobs are finished, copy all trees to the standard path on lxcms03:
```
export TREESET=<YYMMDD>
mkdir /data3/Higgs/${TREESET}
mv AAAOK /data3/Higgs/${TREESET}/PT13TeV
```

Add all chunks:

```
cd /data3/Higgs/${TREESET}
haddChunks.py PT13TeV | tee haddlog_PRODFSR.txt
```
This takes a while. 

Now move the job directories aside:
```
mkdir Chunks
mv *Chunk* Chunks/
```

<!---
Now create links to the merged files into the HZZ_root folder
```
mkdir -p /data3/2014/HZZ_root/${TREESET}/PRODFSR
cd /data3/2014/HZZ_root/${TREESET}/PRODFSR
$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/test/prod/copyMergedFiles.sh /data3/2014/HZZ_out/${TREESET}/PRODFSR
-->
