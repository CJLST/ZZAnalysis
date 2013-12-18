#!/bin/tcsh

cmsenv

mkdir -p DoubleMu
cd DoubleMu
cmsRun ../events_7TeV_DoubleMu.py >&! doubleMu.txt
gzip doubleMu.txt &
cd -

mkdir -p DoubleEle
cd DoubleEle
cmsRun ../events_7TeV_DoubleEle.py >&! doubleEle.txt
gzip doubleEle.txt &
cd -

