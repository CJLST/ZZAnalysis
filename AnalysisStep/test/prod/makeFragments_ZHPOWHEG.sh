#!/bin/bash

csvfile=$1
indicator="ZH"

masslist=()
while IFS=',' read -ra line || [[ -n "$line" ]]; do
	if [[ "$line" == "$indicator"* ]];then
		mass="${line/$indicator/}"
		mass=${mass%%_*}
		#echo $mass
		masslist+=($mass)
	fi
done < $csvfile

for mass in ${masslist[@]};
do
	tplToFragment.py --outdir=pyFragments --outname="LHEProbabilities_ZH_SIG_0PM_H"$mass"_POWHEG.py" --template=pyFragments/LHEProbabilities_ZH_SIG_0PM_H-HMASS-_POWHEG.tpl --indicators="<HMASS>:"$mass
done

