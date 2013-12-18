#!/bin/tcsh

cmsenv



#2012 ABC
mkdir -p DoubleMu
cd DoubleMu
cmsRun ../events_8TeV_2012ABC_DoubleMu.py >&! doubleMu.txt
gzip doubleMu.txt &
cd -

mkdir -p DoubleEle
cd DoubleEle
cmsRun ../events_8TeV_2012ABC_DoubleEle.py >&! doubleEle.txt
gzip doubleEle.txt &
cd -

mkdir -p MuEG
cd MuEG
cmsRun ../events_8TeV_2012ABC_MuEG.py >&! MuEG.txt
gzip MuEG.txt
cd -
