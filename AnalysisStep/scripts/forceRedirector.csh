#!/bin/tcsh
# Tweak miniAOD Condor jobs to use the global redirector for xrootd acces.
# Other options, depending of the region where the files are located:
# root://cmsxrootd.fnal.gov/ 
# root://xrootd-cms.infn.it/
# Note: this is relevant only for jobs running on miniAODs, not for nanoAODs.
# Usage: call forceRedirector.csh from the folder containing Chunks.

foreach f (./*Chunk*/run_cfg.py)
echo "Forcing global redirector in " $f

cat >> $f << EOF 

for i,f in enumerate(process.source.fileNames) :
  process.source.fileNames[i] = "root://cms-xrd-global.cern.ch/"+f

EOF

end
