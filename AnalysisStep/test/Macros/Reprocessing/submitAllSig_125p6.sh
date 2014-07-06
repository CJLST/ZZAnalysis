#!/bin/bash

logdir="./Logs/"

ERGTEV=7
FIRST=0
LAST=24
while [  $ERGTEV -lt 9 ];
do
	COUNTER=$FIRST
	while [  $COUNTER -lt $LAST ];
	do
		bsub -q 1nd -C 0 -o $logdir"lsflog_"$COUNTER"_"$ERGTEV".txt" -e $logdir"lsferr_"$COUNTER"_"$ERGTEV".err" submitBatchSig_125p6.lsf.sh $ERGTEV $COUNTER
		let COUNTER=$COUNTER+1
	done
	let ERGTEV=$ERGTEV+1
done
