#!/bin/bash

if [ "$1" == "" ]; then
 echo "Specify set name, e.g. 140125"
 exit 1;
fi

if [ "$2" == "" ]; then
    DESTSET=$1
else
    DESTSET=$2
fi

prodname="rootuples/${1}/PRODFSR_8TeV/"
dirname="trees/${DESTSET}/PRODFSR_8TeV/"



#echo $prodname
#echo $dirname

if [ -e $dirname ]; then
    rm $dirname/data/*
    rm $dirname/4mu/*
    rm $dirname/4e/*
    rm $dirname/2mu2e/*
    rm $dirname/CR/*
fi

mkdir -p $dirname/data
mkdir -p $dirname/4mu
mkdir -p $dirname/4e
mkdir -p $dirname/2mu2e
mkdir -p $dirname/CR

while read i; do 
    if [[ "$i" != *"#"* ]] && [[ "$i" != "" ]];
    then
	if [ "$i" == "data" ]; then
	    ./run_HZZ4l 0  1 DoubleMu  $prodname/ZZ4lAnalysis_DoubleMu_1963.root $dirname/data/HZZ4lTree_DoubleMu.root
	    ./run_HZZ4l 1  1 DoubleEle $prodname/ZZ4lAnalysis_DoubleEle_1963.root $dirname/data/HZZ4lTree_DoubleEle.root
	    ./run_HZZ4l 2  1 DoubleOr  $prodname/ZZ4lAnalysis_DoubleOr_1963.root $dirname/data/HZZ4lTree_DoubleOr.root
	    ./run_HZZ4l_CR 1 DoubleOr  $prodname/ZZ4lAnalysis_DoubleOr_1963.root $dirname/CR/HZZ4lTree_DoubleOr_CRZLLTree.root
            #./run_HZZ4l_CR 1 Doublemu  $prodname/ZZ4lAnalysis_DoubleMu_1963.root $dirname/CR/HZZ4lTree_DoubleMu_CRZLLTree.root
            #./run_HZZ4l_CR 1 DoubleEle  $prodname/ZZ4lAnalysis_DoubleEle_1963.root $dirname/CR/HZZ4lTree_DoubleEle_CRZLLTree.root
	else
	    FullName=`echo $i | sed s/"ZZ4lAnalysis_"//`
	    ./run_HZZ4l 0  1 $i $prodname/${i}.root $dirname/4mu/HZZ4lTree_$FullName.root
	    ./run_HZZ4l 1  1 $i $prodname/${i}.root $dirname/4e/HZZ4lTree_$FullName.root
	    ./run_HZZ4l 2  1 $i $prodname/${i}.root $dirname/2mu2e/HZZ4lTree_$FullName.root
#	    ./run_HZZ4l_CR 1 $i $prodname/${i}.root $dirname/CR/HZZ4lTree_`echo $i | awk -F_ '{print $2}'`_CRZLLTree.root
	fi
    fi
done < Samples_8TeV.txt
