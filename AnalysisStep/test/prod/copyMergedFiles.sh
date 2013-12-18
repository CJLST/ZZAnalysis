#!/bin/bash

# This script make a copy of merged .root files from AAAOK folder into destination folder and renames them into ZZ4lAnalysis_<dirname>.root
# Must be run from destination folder

# Full path to the folder where AAAOK folder is
if [ "$1" == "" ]; then
 echo "Specify input folder, eg /data3/2013/HZZ_out/131102/PRODFSR_8TeV/"
 exit 1;
fi

INPUT_FOLDER=$1

echo $1

for i in `ls ${INPUT_FOLDER}/ | grep -F -v '_Chunk'`; do ln -sf ${INPUT_FOLDER}/$i/ZZ4lAnalysis.root ZZ4lAnalysis_${i}.root; done
