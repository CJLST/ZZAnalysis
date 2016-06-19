#!/bin/bash

# This script edits run files produced by batch.py, applying two things
# 1) apply cern re-director for xrootd instead of default (higher chances of successful running)
# 2) make Z+jets filters active a second of the sample

# Full path to the folder where AAAOK folder is
if [ "$1" == "" ]; then
 echo "Specify input folder, eg /data3/2014/HZZ_out/140125/PRODFSR_8TeV/"
 exit 1;
fi

INPUT_FOLDER=$1

echo Picking files from $1

for i in `ls ${INPUT_FOLDER}/ | grep '_Chunk'`; do sed s#\\/store#root:\\/\\/cms-xrd-global.cern.ch\\/\\/store#g ${INPUT_FOLDER}/$i/run_cfg.py > ${INPUT_FOLDER}/$i/newrun_cfg.py ; mv ${INPUT_FOLDER}/$i/newrun_cfg.py ${INPUT_FOLDER}/$i/run_cfg.py ; done
for i in `ls ${INPUT_FOLDER}/ | grep 'DYJets'`; do sed s#"option = cms.untracked.int32(0)"#"option = cms.untracked.int32(1)"#g ${INPUT_FOLDER}/$i/run_cfg.py > ${INPUT_FOLDER}/$i/newrun_cfg.py ; mv ${INPUT_FOLDER}/$i/newrun_cfg.py ${INPUT_FOLDER}/$i/run_cfg.py ; done
for i in `ls ${INPUT_FOLDER}/ | grep 'DY1Jets'`; do sed s#"option = cms.untracked.int32(0)"#"option = cms.untracked.int32(1)"#g ${INPUT_FOLDER}/$i/run_cfg.py > ${INPUT_FOLDER}/$i/newrun_cfg.py ; mv ${INPUT_FOLDER}/$i/newrun_cfg.py ${INPUT_FOLDER}/$i/run_cfg.py ; done
for i in `ls ${INPUT_FOLDER}/ | grep 'DY2Jets'`; do sed s#"option = cms.untracked.int32(0)"#"option = cms.untracked.int32(1)"#g ${INPUT_FOLDER}/$i/run_cfg.py > ${INPUT_FOLDER}/$i/newrun_cfg.py ; mv ${INPUT_FOLDER}/$i/newrun_cfg.py ${INPUT_FOLDER}/$i/run_cfg.py ; done
for i in `ls ${INPUT_FOLDER}/ | grep 'DY3Jets'`; do sed s#"option = cms.untracked.int32(0)"#"option = cms.untracked.int32(1)"#g ${INPUT_FOLDER}/$i/run_cfg.py > ${INPUT_FOLDER}/$i/newrun_cfg.py ; mv ${INPUT_FOLDER}/$i/newrun_cfg.py ${INPUT_FOLDER}/$i/run_cfg.py ; done
for i in `ls ${INPUT_FOLDER}/ | grep 'DY4Jets'`; do sed s#"option = cms.untracked.int32(0)"#"option = cms.untracked.int32(1)"#g ${INPUT_FOLDER}/$i/run_cfg.py > ${INPUT_FOLDER}/$i/newrun_cfg.py ; mv ${INPUT_FOLDER}/$i/newrun_cfg.py ${INPUT_FOLDER}/$i/run_cfg.py ; done
for i in `ls ${INPUT_FOLDER}/ | grep 'DYBJets'`; do sed s#"option = cms.untracked.int32(0)"#"option = cms.untracked.int32(2)"#g ${INPUT_FOLDER}/$i/run_cfg.py > ${INPUT_FOLDER}/$i/newrun_cfg.py ; mv ${INPUT_FOLDER}/$i/newrun_cfg.py ${INPUT_FOLDER}/$i/run_cfg.py ; done
for i in `ls ${INPUT_FOLDER}/ | grep 'DYBFiltJets'`; do sed s#"option = cms.untracked.int32(0)"#"option = cms.untracked.int32(3)"#g ${INPUT_FOLDER}/$i/run_cfg.py > ${INPUT_FOLDER}/$i/newrun_cfg.py ; mv ${INPUT_FOLDER}/$i/newrun_cfg.py ${INPUT_FOLDER}/$i/run_cfg.py ; done

