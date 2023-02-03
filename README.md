ZZAnalysis
==========

This is the CJLST framework for analysis of full Run 2 data for Legacy paper.

To install a complete CMSSW 10X area (including this package)
------------------------------
Used for analysis of 2016, 2017, and 2018 data

Please use **CMSSW_10_6_26**. 

For the NanoAnalysis package, before downloading the setup script, clone the nanoAOD repo as follows:
```
cd $CMSSW_BASE/src
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
cd PhysicsTools/NanoAODTools
cmsenv
scram b
```

Download and execute the setup script:
```
wget -O ${TMPDIR}/checkout_10X.csh https://raw.githubusercontent.com/acappati/ZZAnalysis/Run2UL_22_nano/checkout_10X.csh
cd $CMSSW_BASE/src
cmsenv
chmod u+x ${TMPDIR}/checkout_10X.csh
${TMPDIR}/checkout_10X.csh
scramv1 b -j 4
```


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
