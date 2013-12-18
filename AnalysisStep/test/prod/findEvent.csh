#!/bin/tcsh -f

#Requires run:ls:event to avoid same events in different paths
set EVENT=$1 

set SET=130715

if ( $2 == "full" ) then 
  zgrep -E ${EVENT}\|Closed /data3/2013/HZZ_out/${SET}/PRODFSR_8TeV/Double*/*.gz | grep -A1 $EVENT | tee info_${EVENT}.txt
  ##should look also in:   /data3/2013/HZZ_out/${SET}/PRODFSR_8TeV/MuEG*/*.gz
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



echo process.dumpUserData.candidateSrcs = cms.PSet\(Zmm=cms.InputTag\(\"MMCand\"\),Zee=cms.InputTag\(\"\EECand\"\),ZL=\cms.InputTag\(\"\ZlCand\"\),ZLL=\cms.InputTag\(\"\ZLLCand\"\)\) >> ${EVENTID}/run_cfg.py
echo process.dump = cms.Path\(process.dumpUserData\) >> ${EVENTID}/run_cfg.py


#echo $PATH >> ${EVENTID}/run_cfg.py



