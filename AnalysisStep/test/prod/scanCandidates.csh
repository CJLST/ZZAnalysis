#!/bin/tcsh -f

set SET=$1

if ( $1 == "" ) then
 echo "Usage: scanCandidates.csh <name> [full] [file]"
 echo "name = name for the set of files to be created"
 echo "full = grep the original files (otherwise, use the cache file)"
 echo "file = single file to be scanned; if omitted, scan al Double* and MuEG subirectories"
 exit
endif

if ( $2 == full ) then
  if ( $3 != "" ) then
    grep -a "^%% " $3 | sort -s -g -k 2 > _candidates_${SET}.txt
  else 
    zgrep "^%% " Double*/*.gz MuEG*/*.gz | sort -s -g -k 2 > _candidates_${SET}.txt
  endif
endif



cat _candidates_${SET}.txt    | grep -v -e "Candidate" | sed s/".*% "// > candidates_${SET}.txt
#cat _candidates_${SET}.txt    | grep -v -e "Candidate" | grep "b%"| sed s/".*% "// > candidates_${SET}_blinded.txt
#cat _candidates_${SET}.txt    | grep -v -e "Candidate" -e "b%"| sed s/".*% "// > candidates_${SET}_unblinded.txt

cat _candidates_${SET}.txt    | grep -v -e "Candidate" | grep "%100" | sed s/".*% "// > candidates100_${SET}.txt
#cat _candidates_${SET}.txt    | grep -v -e "Candidate" -e "b%"        | grep "%100" | sed s/".*% "// > candidates100_${SET}_unblinded.txt

grep run candidates_${SET}.txt | sort -s -g -k 6 | sort -s -g -k 4 | sort -s -g -k 2 > candidates_${SET}_short.txt
#grep run candidates_${SET}_blinded.txt | sort -s -g -k 6 | sort -s -g -k 4 | sort -s -g -k 2 > candidates_${SET}_blinded_short.txt
#grep run candidates_${SET}_unblinded.txt | sort -s -g -k 6 | sort -s -g -k 4 | sort -s -g -k 2 > candidates_${SET}_unblinded_short.txt
grep run candidates100_${SET}.txt | sort -s -g -k 6 | sort -s -g -k 4 | sort -s -g -k 2 > candidates100_${SET}_short.txt
#grep run candidates100_${SET}_unblinded.txt | sort -s -g -k 6 | sort -s -g -k 4 | sort -s -g -k 2 > candidates100_${SET}_unblinded_short.txt

# echo
# echo "UnBlinded:"
# set file=candidates_${SET}_unblinded_short.txt
# echo "in m4l>70:   " `grep -c run $file`/ `grep -c 4e $file`/`grep -c 4mu $file`/`grep -c 2e2mu $file`
# set file=candidates100_${SET}_unblinded_short.txt
# echo "in m4l>100:  " `grep -c run $file`/ `grep -c 4e $file`/`grep -c 4mu $file`/`grep -c 2e2mu $file`

# echo
# set file=candidates_${SET}_blinded_short.txt
# echo "Blinded:     " `grep -c run $file`/ `grep -c 4e $file`/`grep -c 4mu $file`/`grep -c 2e2mu $file`

echo
echo "All:"
set file=candidates_${SET}_short.txt
echo "in m4l>70:   " `grep -c run $file`/ `grep -c 4e $file`/`grep -c 4mu $file`/`grep -c 2e2mu $file`
set file=candidates100_${SET}_short.txt
echo "in m4l>100:  " `grep -c run $file`/ `grep -c 4e $file`/`grep -c 4mu $file`/`grep -c 2e2mu $file`
echo
cat candidates_${SET}_short.txt | awk '{print $2 ":" $6 ":" $4}' > re${SET}.txt
#cat candidates_${SET}_short.txt | awk '{print $2 ":" $6 ":" $4 ":" $8}' > eventlist${SET}.txt
cat candidates_${SET}_short.txt | awk '{print $2 ":" $6 ":" $4 ":" $8 ":" $10 ":" $12 ":" $14 ":" $16 ":" $18 ":" $20 ":" $22 ":" $24 ":" $26 ":" $28 ":" $30 ":" $32'} > eventlist_${SET}.txt
cat candidates_${SET}_short.txt | awk '{print $2 ":" $6 ":" $4 ":" $8 ":" $10 ":" $12 ":" $14 ":" $16 ":" $18 ":" $20 ":" $22 ":" $24 ":" $26 ":" $28 ":" $30 ":" $32 ":" $34 ":" $36 ":" $38 ":" $40 ":" $42 ":" $44 '} > eventlist_${SET}_extended.txt
#cat candidates_${SET}_unblinded_short.txt | awk '{print $2 ":" $6 ":" $4 ":" $8 ":" $10 ":" $12 ":" $14 ":" $16 ":" $18 ":" $20 ":" $22 ":" $24 ":" $26 ":" $28 ":" $30 ":" $32'} > eventlist${SET}_unblinded.txt
#cat candidates_${SET}_unblinded_short.txt | awk '{print $2 ":" $6 ":" $4 ":" $8 ":" $10 ":" $12 ":" $14 ":" $16 ":" $18 ":" $20 ":" $22 ":" $24 ":" $26 ":" $28 ":" $30 ":" $32 ":" $34 ":" $36 ":" $38 ":" $40 ":" $42 ":" $44 '} > eventlist_${SET}_unblinded_extended.txt
#cat candidates_${SET}_short.txt | awk '{print $2 ":" $6 ":" $4 ":" $8 ":" $10 ":" $12 ":" $16 ":" $22 ":" $24 ":" $26}' > eventlist${SET}.txt
#cat candidates_${SET}_short.txt | awk '{print $2 " " $4 " " $7 " " $8 " " $10 " " $12 " " $16}' | sed s/"+g"// |sed -e s/"m4mu="/0/ -e s/"m4e="/1/ -e s/"m2e2mu="/2/> dump${SET}.txt
#cat candidates_${SET}_unblinded_short.txt | awk '{print $2 ":" $6 ":" $4 ":" $8 ":" $10 ":" $12 ":" $14 ":" $16 ":" $22 ":" $24 ":" $26}' > eventlist${SET}_unblinded.txt


exit

zcat DoubleMu*/*.gz  | grep -e "^%%.*run= " -e "Closed file"  | grep -A1 "run= " >! filesDoubleMu.txt
zcat DoubleEle*/*.gz | grep -e "^%%.*run= " -e "Closed file"  | grep -A1 "run= " >! filesDoubleEle.txt
zcat MuEG*/*.gz      | grep -e "^%%.*run= " -e "Closed file"  | grep -A1 "run= " >! filesMuEG.txt


grep Closed filesDoubleMu.txt | awk '{print $NF}' > filelistDoubleMu.txt
grep Closed filesDoubleEle.txt | awk '{print $NF}' > filelistDoubleEle.txt
grep Closed filesMuEG.txt | awk '{print $NF}' > filelistMuEG.txt
