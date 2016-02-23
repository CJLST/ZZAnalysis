#!/bin/bash

if [ $# -ne 1 ]
then
    echo "usage: moveTrees.csh newFolderName"
    exit 1
fi

## Move all trees to the standard path on lxcms03
export TREEPATH=/data3/Higgs/$1
mkdir ${TREEPATH}
echo "Moving trees ..."
mv AAAOK ${TREEPATH}/PT13TeV

## Add all chunks
cmsenv
cd ${TREEPATH}
echo "Adding chunks ..."
haddChunks.py PT13TeV > haddlog.txt 2>&1

## Put some order
mkdir Chunks
mv PT13TeV/*Chunk* Chunks/
mv PT13TeV/* .
rm -r PT13TeV/
chmod -R g+w .

## Prepare a merged file of all data
mkdir AllData
echo "Merging data trees ..."
hadd AllData/ZZ4lAnalysis.root *2015*/ZZ4lAnalysis.root > haddlog_mergingData.txt 2>&1

cd - > /dev/null
echo "Done. All trees are now in ${TREEPATH}/"
