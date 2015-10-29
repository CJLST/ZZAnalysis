#!/bin/bash

if [ $# -ne 1 ]
then
    echo "usage: moveTrees.csh newFolderName"
    exit 1
fi

## Move all trees to the standard path on lxcms03
export TREEPATH=/data3/Higgs/$1
mkdir ${TREEPATH}
mv AAAOK ${TREEPATH}/PT13TeV

## Add all chunks
cmsenv
cd ${TREEPATH}
haddChunks.py PT13TeV | tee haddlog_PRODFSR.txt

## Put some order
mkdir Chunks
mv PT13TeV/*Chunk* Chunks/
mv PT13TeV/* .
rm -r PT13TeV/
cd -
