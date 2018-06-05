#!/bin/bash

##
# This script reads list of couplings (see Couplings_ACJHUGen.lst for formatting) to create JHUGen anomalous couplings pyFragments.
# It replaces the sample hypothesis <COUPLINGS> in the pyFragment template with whichever actual couplings are put into the lst file.
# Specify the lst file of couplings as argument 1 to this script.
##

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
			if [[ "$prod" == "GG" ]];then
				deccoupling=$coupling
				deccoupling=${deccoupling#*"ghg2=1,0;"}
				deccoupling=${deccoupling#*"ghg4=1,0;"}
				echo $outfile" <DEC_COUPLINGS> -> "$deccoupling
				tplToFragment.py --outdir=pyFragments --outname=$outfile --template="pyFragments/"$tplfile --indicators="<COUPLINGS>:"$coupling --indicators="<DEC_COUPLINGS>:"$deccoupling
			else
				tplToFragment.py --outdir=pyFragments --outname=$outfile --template="pyFragments/"$tplfile --indicators="<COUPLINGS>:"$coupling
			fi
		else
			echo $tplfile" does not exist. tplToFragment.py will not proceed."
		fi
	fi
done

