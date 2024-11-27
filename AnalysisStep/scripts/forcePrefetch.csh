#!/bin/tcsh
# Tweak miniAOD Condor jobs to pre-fetch the input files.
# The prefetch is done with xrdcp using root://cms-xrd-global.cern.ch/
# as a redirector.
# Other options, depending of the region where the files are located:
# root://cmsxrootd.fnal.gov/
# root://xrootd-cms.infn.it/
# Note: this is relevant only for jobs running on miniAODs, not for nanoAODs.
# IMPORTANT: this should be used only for jobs created with the EOS transfer option (-t)
# Usage: call forcePrefetch.csh from the folder containing Chunks.


foreach f (./*Chunk*/run_cfg.py)
echo "Forcing prefetch in " $f


cat >> $f << EOF 

# Remove the global redirector from the file name, in case it was specified
for i,f in enumerate(process.source.fileNames) :
  process.source.fileNames[i] = f.replace("root://cms-xrd-global.cern.ch/",'')

import os
for i,f in enumerate(process.source.fileNames) :
  tempdir = "."
  print("prefetching", f)
  os.system("xrdcp -d 1 -f root://cms-xrd-global.cern.ch/"+f+" "+tempdir)
  process.source.fileNames[i] = "file:"+tempdir+"/"+os.path.basename(f)
print("New fileNames:")
print(process.source.fileNames)

EOF

end


