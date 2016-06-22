ZZAnalysis
==========

To install a complete CMSSW area (including this package)
------------------------------
Please use a CMSSW_8_0_0 version >= 8_0_6.

Download and execute the setup script:
```

wget -O /tmp/checkout_80X.csh https://raw.githubusercontent.com/CJLST/ZZAnalysis/miniAOD_80X_with_loose_ele/checkout_80X.csh

cd $CMSSW_BASE/src
cmsenv
source /tmp/checkout_80X.csh
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
git push origin miniAOD_76X
```

Otherwise you can make a fork of the repository, develop therein, and make a pull request in the same way as for CMSSW.

Code documentation
------------------
Please see the [gitHub wiki](https://github.com/CJLST/ZZAnalysis) for more documentation.
