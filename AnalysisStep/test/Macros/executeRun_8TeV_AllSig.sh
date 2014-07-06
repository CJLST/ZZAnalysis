#!/bin/bash

#logdir="./logs/8TeV/"
logdir="/afs/cern.ch/work/l/lfeng/public_write/usarica/Production/data/logs/"

#echo $prodname
#echo $dirname

while read i; do 
    if [[ "$i" != *"#"* ]] && [[ "$i" != "" ]];
    then
	bsub -q 1nd -o $logdir"lsflog_$i.txt" -e $logdir"lsferr_$i.err" submit8TeV_Sig.lsf.sh $i
    fi
done < Samples_8TeV_Sig.txt