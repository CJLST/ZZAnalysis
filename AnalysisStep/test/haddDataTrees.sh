#!/bin/bash

if [ $# -ne 1 ]
then
    echo "usage: haddDataTrees.sh TreeSet"
else

    mkdir DataTrees_$1
    
    for PD in DoubleEG DoubleMu MuonEG SingleEle
    do
	for period in 2015B 2015C_50ns 2015D 2015Dv4
	do
	    eos cp root://lxcms03://data3/Higgs/$1/${PD}${period}/ZZ4lAnalysis.root DataTrees_$1/ZZ4lAnalysis_${PD}${period}.root
	done
    done
    
    hadd DataTrees_$1/ZZ4lAnalysis_allData.root DataTrees_$1/ZZ4lAnalysis_*.root

fi