#!/bin/bash

cd ${LS_SUBCWD}

echo "LSF job running in: " `pwd`

eval `scram runtime -sh`

echo $CMSSW_VERSION

ERGTEV=$1
COUNTER=$2
let NEXT=$COUNTER+1
COMMAND="reweightIdeal_125p6.c+($ERGTEV,$COUNTER,$NEXT)"
root -b -l -q "loadLib.C" $COMMAND


#source submitReweightIdeal_125p6.sh $1 $2
