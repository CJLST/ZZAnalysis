Summary of steps for a full standard production
-----------------------------------------------
-----------------------------------------------


PRIMARY trees
-------------

```
./batch.py -o PRODFSR analyzer_2011.py
./batch.py -o PRODFSR_8TeV analyzer_2012.py
```

Once jobs are done:

```
cd PRODFSR
cd PRODFSR_8TeV
../checkProd.csh
```

To resubmit failed jobs (repeat until all jobs succeed)

```
../cleanup.csh
../resubmit.csh
# wait all jobs to finish
../checkProd.csh
```

Notes:
*   Sometimes jobs get stuck in the processing nodes. This can be checked with bjobs -l <jobid> - if job is running since a lot of time but appears to be idle (check e.g. IDLE_FACTOR), it's probably best to kill it; wait that bjobs reports it as finished (can take a while), and resubmit it. 
*   Sometimes failures are permanent (job keeps failing), in this case it must be investigated looking at the log.txt.gz file. In rare cases, this happens to be due to an EOS failure. If a file is reported as unaccessible: 
    *     try opening it in plain root with xrootd; that should give a reproducible error message 
    *     check eos fileinfo <file> : what usually the problem is that servers holding both replicas are down, or that one of the two replicas is reported with incorrect information. Open an EOS ticket with a clear pointer to the file and to the eos fileinfo message; it is normally addressed within few hours.


Once all jobs are finished, copy all trees to the standard path on lxcms00:
```
export TREESET=<YYMMDD>
mkdir /data3/2014/HZZ_out/${TREESET}
mv AAAOK /data3/2014/HZZ_out/${TREESET}/PRODFSR
mv AAAOK /data3/2014/HZZ_out/${TREESET}/PRODFSR_8TeV
```

hadd all chunks:

```
cd /data3/2014/HZZ_out/${TREESET}
haddChunks.py PRODFSR | tee haddlog_PRODFSR.txt
haddChunks.py PRODFSR_8TeV | tee haddlog_PRODFSR_8TeV.txt
```
This takes a while. 

Now create links to the merged files into the HZZ_root folder:
```
mkdir -p /data3/2014/HZZ_root/${TREESET}/PRODFSR_8TeV
cd /data3/2014/HZZ_root/${TREESET}/PRODFSR_8TeV
$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/test/prod/copyMergedFiles.sh /data3/2014/HZZ_out/${TREESET}/PRODFSR_8TeV
```

Once ready, we still have to merge few datasets: data for 7TeV and data, ZZ+ZZext for 8TeV
```
$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/test/prod/hadd_7TeV.csh
$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/test/prod/hadd_8TeV.csh
```


SECONDARY trees
---------------
```
export TREESET=<YYMMDD>
cd $CMSSW_BASE/src/ZZAnalysis/AnalysisStep/test/Macros
ln -s /data3/2014/HZZ_root/ rootuples
ln -s /data3/2014/HZZ_stat/ trees
mkdir -p /data3/2014/HZZ_stat/${TREESET}/PRODFSR_8TeV
./executeRun_7TeV.sh ${TREESET} | tee trees/${TREESET}/PRODFSR/secondary_log.txt
./executeRun_8TeV.sh ${TREESET} | tee trees/${TREESET}/PRODFSR_8TeV/secondary_log.txt
```

Please store the logs as suggested, as they contain useulf information to understand what happens if some file has a problem.
