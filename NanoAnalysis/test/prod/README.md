NanoAOD production
==================

Job submission for nanoAODs folows the same recipe as for the miniAOD framework (see [documentation](https://github.com/CJLST/ZZAnalysis/wiki/SubmittingJobs).

**IMPORTANT**: contrary to what happens with the miniAOD tree maker, they python executable in the batch scripts depends directly on the .py code in the working area. The py code should therefore not be modified while production is being run

The script `haddData.csh` can be run from the production folder after chunks have been hadded with `haddChunks.py`, in order to merge the different data PDs into a single root file.
