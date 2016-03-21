ZZAnalysis
==========

To install a complete CMSSW area (including this package)
------------------------------
Please use a CMSSW_7_6_X version >= 7_6_3_patch2.

Download and execute the setup script:
```
wget -O /tmp/checkout_70X.csh https://raw.githubusercontent.com/CJLST/ZZAnalysis/2l2q/checkout_70X.csh
cd $CMSSW_BASE/src
cmsenv
source /tmp/checkout_70X.csh
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
git push origin 2l2q
```

Otherwise you can make a fork of the repository, develop therein, and make a pull request in the same way as for CMSSW.

Code documentation
------------------
Please see the [gitHub wiki](https://github.com/CJLST/ZZAnalysis) for more documentation.
