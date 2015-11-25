#!/bin/bash

if [ $# -ne 1 ]
then
    echo "usage: haddDataTrees.sh TreeSet"
else

    mkdir DataTrees_$1
    
    for PD in DoubleEG DoubleMu MuonEG SingleEle
    do
	for period in 2015B 2015C 2015C_50ns 2015D 2015Dv4
	do
	    echo "Copying data trees from /data3/Higgs/$1/ ..."
	    eos cp root://lxcms03://data3/Higgs/$1/${PD}${period}/ZZ4lAnalysis.root DataTrees_$1/ZZ4lAnalysis_${PD}${period}.root
	done
    done

    echo "Merging all trees ..."
    hadd DataTrees_$1/ZZ4lAnalysis_allData.root DataTrees_$1/ZZ4lAnalysis_*.root
    echo "Done. Merged data trees are now in ./DataTrees_$1/ZZ4lAnalysis_allData.root"
    
fi