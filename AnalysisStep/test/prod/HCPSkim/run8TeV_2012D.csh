#!/bin/tcsh

cmsenv

# mkdir -p DoubleMu
# cd DoubleMu
# cmsRun ../events_8TeV_DoubleMu.py >&! doubleMu.txt
# gzip doubleMu.txt &
# cd -

# mkdir -p DoubleEle
# cd DoubleEle
# cmsRun ../events_8TeV_DoubleEle.py >&! doubleEle.txt
# gzip doubleEle.txt &
# cd -

# mkdir -p MuEG
# cd MuEG
# cmsRun ../events_8TeV_MuEG.py >&! MuEG.txt
# gzip MuEG.txt
# cd -


#2012 D
mkdir -p DoubleMu_2012D
cd DoubleMu_2012D
cmsRun ../events_8TeV_2012D_DoubleMu.py >&! doubleMu.txt
gzip doubleMu.txt &
cd -

mkdir -p DoubleEle_2012D
cd DoubleEle_2012D
cmsRun ../events_8TeV_2012D_DoubleEle.py >&! doubleEle.txt
gzip doubleEle.txt &
cd -

mkdir -p MuEG_2012D
cd MuEG_2012D
cmsRun ../events_8TeV_2012D_MuEG.py >&! MuEG.txt
gzip MuEG.txt
cd -
