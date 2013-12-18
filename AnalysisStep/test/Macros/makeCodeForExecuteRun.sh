#!/bin/bash

# Loops over the merged files in ${inputDir} and makes the code to be pasted into execute_run_7(8)TeV.(c)sh scripts

if [ $# -lt 2 ]
then
  echo "Need to pass as input the run period and input directory!"
  exit 0
fi

energy=$1  # 0 = 7TeV, 1 = 8TeV
inputDir=$2

#4mu
for i in `ls ${inputDir}`; do echo "./run_HZZ4l 0 ${energy} \$prodname/${i} \$dirname/4mu/HZZ4lTree_`echo $i | awk -F_ '{print $2}'`"; done
echo
#4e
for i in `ls ${inputDir}`; do echo "./run_HZZ4l 1 ${energy} \$prodname/${i} \$dirname/4e/HZZ4lTree_`echo $i | awk -F_ '{print $2}'`"; done
echo
#2mu2e
for i in `ls ${inputDir}`; do echo "./run_HZZ4l 2 ${energy} \$prodname/${i} \$dirname/2mu2e/HZZ4lTree_`echo $i | awk -F_ '{print $2}'`"; done
echo
#CR
for i in `ls ${inputDir}`; do echo "./run_HZZ4l_CR ${energy} \$prodname/${i} \$dirname/CR/HZZ4lTree_`echo $i | awk -F_ '{print $2}'`"; done