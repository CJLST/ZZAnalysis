#!/bin/bash

# This script makes links of all merged .root files the input folder into the current folder and renames them into ZZ4lAnalysis_<dirname>.root

# Full path to the folder where AAAOK folder is
if [ "$1" == "" ]; then
 echo "Specify input folder, eg /data3/2014/HZZ_out/140125/PRODFSR_8TeV/"
 exit 1;
fi

INPUT_FOLDER=$1

echo Picking files from $1

for i in `ls ${INPUT_FOLDER}/ | grep -F -v '_Chunk'`; do ln -sf ${INPUT_FOLDER}/$i/ZZ4lAnalysis.root ZZ4lAnalysis_${i}.root; done
