#!/bin/bash

logdir="/afs/cern.ch/work/l/lfeng/public_write/usarica/Production/7TeV/logs/"

while read i; do 
    if [[ "$i" != *"#"* ]] && [[ "$i" != "" ]];
    then
	bsub -q 1nd -C 0 -o $logdir"lsflog_$i.txt" -e $logdir"lsferr_$i.err" submit7TeV_Bkg.lsf.sh $i
    fi
done < Samples_gg2VV.txt