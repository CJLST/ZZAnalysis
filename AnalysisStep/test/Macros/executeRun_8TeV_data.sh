#!/bin/bash

#prodname="./data/files"
#dirname="/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_8TeV"
prodname="/afs/cern.ch/work/l/lfeng/public_write/usarica/Production/8TeV/files"
dirname="/afs/cern.ch/work/l/lfeng/public_write/usarica/Production/8TeV"

	    ./run_HZZ4l 0  1 DoubleMu  $prodname/ZZ4lAnalysis_DoubleMu_1963.root  $dirname/data/HZZ4lTree_DoubleMu.root
	    ./run_HZZ4l 1  1 DoubleEle $prodname/ZZ4lAnalysis_DoubleEle_1963.root $dirname/data/HZZ4lTree_DoubleEle.root
	    ./run_HZZ4l 2  1 DoubleOr  $prodname/ZZ4lAnalysis_DoubleOr_1963.root  $dirname/data/HZZ4lTree_DoubleOr.root
	    ./run_HZZ4l_CR 1 DoubleOr  $prodname/ZZ4lAnalysis_DoubleOr_1963.root  $dirname/CR/HZZ4lTree_DoubleOr_CRZLLTree.root
