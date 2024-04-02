Z+X with NanoAOD
================

All the scripts contained in this folder are meant to be used for the computation of Z+X using the fake-rate method and were written for the output of the MiniAOD framework.
To use the same scripts, the output of the NanoAOD framework should be processed with `NanoConverter.py`. This output of this script is a root file in the same format as if the sample were processed with the MiniAOD framework. At present, it only contains information for the computation of Z+X. It can be expanded to be a universal converter.

The syntax is the following

```
python3 NanoConverter.py --input="[path to the input file]" --output="[path to the output file]"
```
If you are processing MC, the flag `--mc` should be added. When processing samples that do not have the Z+L CR, like qqZZ, the flag `--skipZL` should be activated.
