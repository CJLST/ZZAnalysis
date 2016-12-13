#!/bin/bash

csvfile=$1
indicator="WminusH"

masslist=()
while IFS=',' read -ra line || [[ -n "$line" ]]; do
	if [[ "$line" == "$indicator"* ]];then
		mass="${line/$indicator/}"
		mass=${mass%%_*}
		masslist+=($mass)
	fi
done < $csvfile

indicator="WplusH"
while IFS=',' read -ra line || [[ -n "$line" ]]; do
	if [[ "$line" == "$indicator"* ]];then
		mass="${line/$indicator/}"
		mass=${mass%%_*}
		let found=0
		for x in ${masslist[@]}
		do
			if [[ "$x" == "$mass" ]];then
				let found=1
			fi
		done
		if [ $found -eq 0 ];then
			masslist+=($mass)
		fi
	fi
done < $csvfile

for mass in ${masslist[@]}
do
	tplToFragment.py --outdir=pyFragments --outname="LHEProbabilities_WH_SIG_0PM_H"$mass"_POWHEG.py" --template=pyFragments/LHEProbabilities_WH_SIG_0PM_H-HMASS-_POWHEG.tpl --indicators="<HMASS>:"$mass
done

