#!/bin/bash

csvfile=$1
newcsvfile=$csvfile"tmp"
rm -f $newcsvfile
indicatorlist=("ggH" "VBFH" "ZH" "WminusH" "WplusH")
while IFS='' read -r line || [[ -n "$line" ]]; do
	let hasindicator=0
	for indicator in "${indicatorlist[@]}"
	do
	if [[ "$line" == "$indicator"* ]] || [[ "$line" == "#$indicator"* ]];then
		mass="${line/$indicator/}"
		mass=${mass%%,*}
		mass=${mass%%_*}
		mass=${mass#\#}
		newline=${line/"-HMASS-"/$mass}
		let hasindicator=1
		echo $newline >> $newcsvfile
	fi
	done
	if [ $hasindicator -eq 0 ];then
		echo $line >> $newcsvfile
	fi
done < $csvfile
mv $newcsvfile $csvfile
