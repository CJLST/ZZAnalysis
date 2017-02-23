#!/bin/bash

csvfile=$1
indicator="ttH"

masslist=()
while IFS=',' read -ra line || [[ -n "$line" ]]; do
	if [[ "$line" == "$indicator"* ]] || [[ "$line" == "#$indicator"* ]];then
		mass="${line/$indicator/}"
		mass=${mass%%_*}
		mass=${mass#\#}
		#echo $mass
		masslist+=($mass)
	fi
done < $csvfile

for mass in ${masslist[@]};
do
	tplToFragment.py --outdir=pyFragments --outname="LHEProbabilities_ttH_SIG_0PM_H"$mass"_POWHEG.py" --template=pyFragments/LHEProbabilities_ttH_SIG_0PM_H-HMASS-_POWHEG.tpl --indicators="<HMASS>:"$mass
done

