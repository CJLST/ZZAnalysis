#!/bin/bash

lstfile=$1

indicators=()
couplings=()
while IFS='' read -r line || [[ -n "$line" ]]; do
	indicator=${line%%|*}
	coupling=${line#*|}
	indicators+=($indicator)
	couplings+=($coupling)
done < $lstfile
let nACs=${#indicators[@]}

csvprodstr=("ggH" "VBFH" "WH" "ZH" "HJJ" "ttH") # Comes from the csv/lst file
tplprodstr=("GG" "VBF" "WH" "ZH" "HJJ" "ttH") # Actual tpl production indicator
let nProds=${#csvprodstr[@]}
for (( i=0; i<${nACs}; i++ ));
do
	indicator=${indicators[$i]}
	coupling=${couplings[$i]}
	mass=${indicator##*_M}
	hypo=${indicator%%_*}
	prod=""
	csvprod=""
	for (( p=0; p<${nProds}; p++ ));
	do
		if [[ "$indicator" == "${csvprodstr[$p]}"* ]];then
			prod=${tplprodstr[$p]}
			csvprod=${csvprodstr[$p]}
		fi
	done
	if [[ "$prod" != "" ]];then
		hypo=${hypo##*$csvprod}
		tplfile="LHEProbabilities_"$prod"_SIG_-HYPO-_H-HMASS-_ACJHUGen.tpl"
		if [ -f pyFragments/$tplfile ];then
			outfile=$tplfile
			outfile=${outfile/"-HMASS-"/$mass}
			outfile=${outfile/"-HYPO-"/$hypo}
			outfile=${outfile/"tpl"/"py"}
			echo $outfile" <COUPLINGS> -> "$coupling
			tplToFragment.py --outdir=pyFragments --outname=$outfile --template="pyFragments/"$tplfile --indicators="<COUPLINGS>:"$coupling
		else
			echo $tplfile" does not exist. tplToFragment.py will not proceed."
		fi
	fi
done

