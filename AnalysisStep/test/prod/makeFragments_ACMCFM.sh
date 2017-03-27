#!/bin/bash

##
# This script reads list of couplings (see Couplings_ACMCFM.lst for formatting) to create MCFM anomalous couplings pyFragments.
# It replaces the sample hypothesis <COUPLINGS> in the pyFragment template with whichever actual couplings are put into the lst file.
# Specify the lst file of couplings as argument 1 to this script.
##

lstfile=$1

hypoindicators=()
bsiindicators=()
couplings=()
while IFS='' read -r line || [[ -n "$line" ]]; do
	coupling=${line##*|}
	indicator=${line%|$coupling}
	hypoindicator=${indicator%|*}
	bsiindicator=${indicator#*|}
	hypoindicators+=($hypoindicator)
	bsiindicators+=($bsiindicator)
	couplings+=($coupling)
done < $lstfile
let nACs=${#couplings[@]}

for (( i=0; i<${nACs}; i++ ));
do
	hypo=${hypoindicators[$i]}
	bsi=${bsiindicators[$i]}
	coupling=${couplings[$i]}

	tplfile="LHEProbabilities_GG_"$bsi"_-HYPO-_H125_MCFM.tpl"

	if [ -f pyFragments/$tplfile ];then
		outfile=$tplfile
		outfile=${outfile/"-HYPO-"/$hypo}
		outfile=${outfile/"tpl"/"py"}
		echo $outfile" <COUPLINGS> -> "$coupling
		tplToFragment.py --outdir=pyFragments --outname=$outfile --template="pyFragments/"$tplfile --indicators="<COUPLINGS>:"$coupling
	else
		echo $tplfile" does not exist. tplToFragment.py will not proceed."
	fi
done

