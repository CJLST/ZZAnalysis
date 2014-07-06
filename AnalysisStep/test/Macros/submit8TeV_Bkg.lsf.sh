#!/bin/bash

cd ${LS_SUBCWD}

echo "LSF job running in: " `pwd`

eval `scram runtime -sh`

echo $CMSSW_VERSION

source executeRun_8TeV_Bkg.sh $1
