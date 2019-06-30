#!/bin/bash

set -euo pipefail

startnumber=0

for filename in AAAFAIL*; do
  if [[ $filename =~ ^AAAFAIL[0-9]+$ ]]; then
    newnumber=${filename/AAAFAIL/}
    if [ $newnumber -gt $startnumber ]; then
      startnumber=$newnumber
    fi
  fi
done

for i in {1..192}; do
  foldernumber=$(expr $startnumber + $i)
  echo -n "it's currently "; date
  echo "running for time ${i}/192 using AAAFAIL$foldernumber"

  checkProd.csh lf AAAOK AAAFAIL$foldernumber |& tee checkprodlog$foldernumber
  #the jobs that have finished are now moved to AAAOK
  #that means those symlinks are now broken
  #so we can delete them
  for folder in AAAFAIL*; do
    find $folder -xtype l -delete
    if not ls $folder/*Chunk* >& /dev/null; then
      rm -r $folder
    fi
  done

  date >> checkprodlog$foldernumber
  echo "Chunks remaining:"
  ls -d *Chunk* || break
  if [ -d AAAFAIL$foldernumber ]; then
    (
      cd AAAFAIL$foldernumber
      echo "resubmitting $(ls -d *Chunk* 2>/dev/null | wc -l) failed chunks"
      if ls -d *Chunk* >& /dev/null; then
        cleanup.csh
        resubmit_Condor.csh
      fi
    )
  fi
  echo -n "it's currently "; date
  echo "will run again in 15 minutes, which will be time $((${i}+1))/192"
  sleep 15m
done
