#!/bin/tcsh

if ( $1 == "" ) then
    echo "\nusage: moveTrees.csh newFolderName\n"
    exit 1
endif

if ( `uname -n` != "lxcms03.cern.ch" ) then
    echo "\nThis script is designed to run on lxcms03. Aborting.\n"
    exit 1
endif

if ( ! -d AAAOK) then
    echo "\nSource directory AAAOK not found. Aborting.\n"
    exit 1    
endif

eval `scram runtime -csh`

set SETNAME=`basename $PWD`

set STORAGEPATH=/data3/Higgs
set TREEPATH=${STORAGEPATH}/$1


if (-d $TREEPATH) then
    if (-w $TREEPATH) then
	echo "INFO: Target path $TREEPATH already existing; adding files therein."
    else 
	echo "Target path $TREEPATH not writable. Aborting.\n"
	exit 1
    endif
else
    mkdir $TREEPATH
    if ($? != 0) then
	echo "\nError creating destination folder $TREEPATH. Aborting.\n"
	exit 1    
    endif
    #Make the directory writable by others
    chmod g+w $TREEPATH
endif


# Check if there is a data target directort already, if we have data files
set MERGEDATA16=0
set MERGEDATA17=0
set MERGEDATA18=0
if (`find AAAOK -maxdepth 1 -name "*2016*" -print -quit` != "") then
    set MERGEDATA16=1
endif
if (`find AAAOK -maxdepth 1 -name "*2017*" -print -quit` != "") then
    set MERGEDATA17=1
endif
if (`find AAAOK -maxdepth 1 -name "*2018*" -print -quit` != "") then
    set MERGEDATA18=1
endif


if ($MERGEDATA16 || $MERGEDATA17 || $MERGEDATA18) then
    if (-d ${TREEPATH}/AllData) then
	echo "\n${TREEPATH}/AllData existing, cannot overwrite. Aborting.\n"
	exit 1
    endif
endif



## Add all chunks (in the local directory
echo "Adding chunks ..."
haddChunks.py AAAOK >&! AAAOK/haddlog_${SETNAME}.txt

## Move all trees to the standard path on lxcms03
echo "Moving trees ..." #FIXME
mv AAAOK ${TREEPATH}/${SETNAME}

## Put some order
cd ${TREEPATH}
set CHUNKPATH=Chunks_${SETNAME}
mkdir ${CHUNKPATH}
mv ${SETNAME}/*Chunk* ${CHUNKPATH}
mv ${SETNAME}/* .
rmdir ${SETNAME} # If any of the above fails (eg because a sample was already present), ${SETNAME} will not be empty and will not be deleted


## Prepare a merged file of all data


if ($MERGEDATA16 || $MERGEDATA17 || $MERGEDATA18) then
    mkdir AllData
    echo "Merging data trees ..."
    if ($MERGEDATA16) hadd AllData/ZZ4lAnalysis2016.root *2016*/ZZ4lAnalysis.root >&! haddlog_${SETNAME}_mergingData.txt
    if ($MERGEDATA17) hadd AllData/ZZ4lAnalysis.root *2017*/ZZ4lAnalysis.root >&! haddlog_${SETNAME}_mergingData.txt
    if ($MERGEDATA18) hadd AllData/ZZ4lAnalysis.root *2018*/ZZ4lAnalysis.root >&! haddlog_${SETNAME}_mergingData.txt
endif

cd - > /dev/null
echo "Done. All trees are now in ${TREEPATH}/\n"
