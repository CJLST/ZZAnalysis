#!/bin/bash

#prodname="./data/files"
#dirname="/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_7TeV"
#dirname="./data"
#prodname="/afs/cern.ch/work/g/gritsan/public/backup/data/files"
#dirname="/afs/cern.ch/work/g/gritsan/public/backup/data"
prodname="/afs/cern.ch/work/l/lfeng/public_write/usarica/Production/7TeV/files"
dirname="/afs/cern.ch/work/l/lfeng/public_write/usarica/Production/7TeV"

	    ./run_HZZ4l 0  0 DoubleMu  $prodname/ZZ4lAnalysis_DoubleMu.root  $dirname/data/HZZ4lTree_DoubleMu.root
	    ./run_HZZ4l 1  0 DoubleEle $prodname/ZZ4lAnalysis_DoubleEle.root $dirname/data/HZZ4lTree_DoubleEle.root
	    ./run_HZZ4l 2  0 DoubleOr  $prodname/ZZ4lAnalysis_DoubleOr.root  $dirname/data/HZZ4lTree_DoubleOr.root
	    ./run_HZZ4l_CR 0 DoubleOr  $prodname/ZZ4lAnalysis_DoubleOr.root  $dirname/CR/HZZ4lTree_DoubleOr_CRZLLTree.root
