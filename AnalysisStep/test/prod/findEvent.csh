#!/bin/tcsh -f

set EVENT=$1 
set SET=160203

if ( $2 == "full" ) then 
    zgrep -a -E ${EVENT}\|Closed /data3/Higgs/$SET/Chunks/Double*/*.gz /data3/Higgs/$SET/Chunks/MuEG*/Double*/*.gz | grep -A1 $EVENT | tee info_${EVENT}.txt
endif

set FILE = `grep Closed info_${EVENT}.txt | awk '{print $NF}'`
set EVENTID = `grep $EVENT info_${EVENT}.txt | awk '{print $1}'`
set PATH = `echo $EVENTID | sed s/':.*$'//`
set PATH = `dirname $PATH`
set EVENTID = `echo $EVENTID | sed s/'^.*gz:'//`

echo
echo "event ID:  $EVENTID"
echo "job path:  $PATH"
echo "rooth file: $FILE" 

mkdir -p ${EVENTID}
cp ${PATH}/run_cfg.py ${EVENTID}/run_cfg.py


echo process.source.fileNames = cms.untracked.vstring\(\"$FILE\"\) >> ${EVENTID}/run_cfg.py


echo process.source.eventsToProcess = cms.untracked.VEventRange\(\"$EVENTID\"\) >> ${EVENTID}/run_cfg.py



#echo process.dumpUserData.candidateSrcs = cms.PSet\(ZZ=cms.InputTag\(\"ZZCand\"\),ZLL=\cms.InputTag\(\"\ZLLCand\"\)\) >> ${EVENTID}/run_cfg.py
echo process.dump = cms.Path\(process.dumpUserData\) >> ${EVENTID}/run_cfg.py


#echo $PATH >> ${EVENTID}/run_cfg.py



