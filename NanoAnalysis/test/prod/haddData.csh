#!/bin/tcsh
#
# After production is complete, the Chunks of each sample  should be merged with:
#
# haddChunks.py AAAOK 
#
# once this is done, the root files for the different PDs of each data taking years can be merged into a signle file per year with this script.
# It also archives the Chunks in a subfolder.
# 
# usage: 
# ../haddData.csh
#

if ( ! -d AAAOK) then
    echo "AAAOK folder not found!"
    exit 1
endif


cd AAAOK
# Archive Chunk folders, if any
if (`find . -maxdepth 1 -name "*_Chunk*" -print -quit` != "") then
    echo "Archiving Chunks..."
    mkdir -p Chunks
    mv *_Chunk* Chunks
endif

set dataYears=("2016" "2017" "2018" "2022" "2023")

foreach year ( $dataYears )
  if (`find . -maxdepth 1 -name "*${year}*" -print -quit` != "") then
    echo "hadding ${year}..."
    mkdir -p Data${year}
    haddnano.py Data${year}/ZZ4lAnalysis.root *${year}*/ZZ4lAnalysis.root
    if ( $status != 0 ) then
	echo "hadd failed!"
	exit 1 
    endif
  endif
end

  

  
