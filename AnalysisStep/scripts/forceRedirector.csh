#!/bin/tcsh


foreach f (./*Chunk*/run_cfg.py)
echo "Forcing global redirector in " $f

cat >> $f << EOF 

for i,f in enumerate(process.source.fileNames) :
  process.source.fileNames[i] = "root://cms-xrd-global.cern.ch/"+f

EOF

end
