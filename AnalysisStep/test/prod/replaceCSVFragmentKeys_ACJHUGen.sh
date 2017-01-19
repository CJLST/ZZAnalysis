#!/bin/bash

csvfile=$1
newcsvfile=$csvfile"tmp"
rm -f $newcsvfile

csvprodstr=("ggH" "VBFH" "WH" "ZH" "HJJ" "ttH") # Comes from the csv/lst file
tplprodstr=("GG" "VBF" "WH" "ZH" "HJJ" "ttH") # Actual tpl production indicator
let nProds=${#csvprodstr[@]}
while IFS='' read -r line || [[ -n "$line" ]]; do
	let hasindicator=0
	for (( p=0; p<${nProds}; p++ ));
	do
		prod=""
		csvprod=""
		if [[ "$line" == "${csvprodstr[$p]}"* ]];then
			prod=${tplprodstr[$p]}
			csvprod=${csvprodstr[$p]}
		fi
		if [[ "$prod" != "" ]];then
			indicator=${line%%,*}
			indicator=${indicator#\#}
			mass=${indicator##*_M}
			hypo=${indicator%%_*}
			hypo=${hypo##*$csvprod}
			defstr="LHEProbabilities_-PROD-_SIG_-HYPO-_H-HMASS-_ACJHUGen.py"
			repstr=$defstr
			repstr=${repstr/"-PROD-"/$prod}
			repstr=${repstr/"-HMASS-"/$mass}
			repstr=${repstr/"-HYPO-"/$hypo}
			newline=${line/$defstr/$repstr}
			echo $newline >> $newcsvfile
			let hasindicator=1
		fi
	done
	if [ $hasindicator -eq 0 ];then
		echo $line >> $newcsvfile
	fi
done < $csvfile
mv $newcsvfile $csvfile
