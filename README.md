ZZAnalysis
==========

This is the CJLST framework for analysis of full Run 2 data for Legacy paper.

To install a complete CMSSW 10X area (including this package)
------------------------------
Used for analysis of 2016, 2017, and 2018 data

Please use **CMSSW_10_6_30**. 

Download and execute the setup script in a el7 Singularity container:
```
cmssw-el7
cmsrel CMSSW_10_6_30
cd CMSSW_10_6_30/src
cmsenv
wget -O ${TMPDIR}/checkout_10X.csh https://raw.githubusercontent.com/CJLST/ZZAnalysis/Run2UL_22/checkout.csh
chmod u+x ${TMPDIR}/checkout.csh
${TMPDIR}/checkout.csh
scramv1 b -j 4
```

Please use `ZZAnalysis/start_el7.sh` to open a slc7 container in new shells, in order to have access to the Condor queues.


To update this package from the release
------------------------------------------
In the package directory, simply issue
```
git pull
```

To commit and push new changes
------------------------------
To commit directly (you need write access to the repository):
```
git pull
[edit files]
```
Once you are ready to commit
```
git pull
git add [files to be added]
git commit -m ["commit message"] [files to be added]
git push origin Run2Legacy
```

Otherwise you can make a fork of the repository, develop therein, and make a pull request in the same way as for CMSSW.

Code documentation
------------------
Please see the [gitHub wiki](https://github.com/CJLST/ZZAnalysis) for more documentation.
