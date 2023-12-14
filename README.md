ZZAnalysis
==========

This is the CJLST framework for analysis of Run 3 data.

To install a complete CMSSW 13X area (including this package)
------------------------------
Used for analysis of 2022 data and beyond

Please use **CMSSW_13_0_16**. 

Download and execute the setup script:
```
cmsrel CMSSW_13_0_16
cd CMSSW_13_0_16/src
cmsenv
wget -O ${TMPDIR}/checkout.csh https://raw.githubusercontent.com/CJLST/ZZAnalysis/Run3/checkout_13X.csh
chmod u+x ${TMPDIR}/checkout.csh
${TMPDIR}/checkout.csh
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
