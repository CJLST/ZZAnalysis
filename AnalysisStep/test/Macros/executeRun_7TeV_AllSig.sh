#!/bin/bash

#prodname="./data/files"
#dirname="/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_7TeV"
logdir="/afs/cern.ch/work/g/gritsan/public/backup/data/logs/7TeV/"

#echo $prodname
#echo $dirname

while read i; do 
    if [[ "$i" != *"#"* ]] && [[ "$i" != "" ]];
    then
	bsub -q 1nd -o $logdir"lsflog_$i.txt" -e $logdir"lsferr_$i.err" submit7TeV_Sig.lsf.sh $i
    fi
done < Samples_7TeV_Sig.txt