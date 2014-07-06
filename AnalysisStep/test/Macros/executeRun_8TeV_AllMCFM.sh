#!/bin/bash

logdir="/afs/cern.ch/work/l/lfeng/public_write/usarica/Production/8TeV/logs/"

while read i; do 
    if [[ "$i" != *"#"* ]] && [[ "$i" != "" ]];
    then
	bsub -q 8nh -o $logdir"lsflog_$i.txt" -e $logdir"lsferr_$i.err" submit8TeV_Bkg.lsf.sh $i
#	echo $logdir"lsflog_$i.txt" -e $logdir"lsferr_$i.err" $i
    fi
done < Samples_MCFM.txt