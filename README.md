ZZAnalysis
==========

This package can be checked out with:

```
git clone https://github.com/CJLST/ZZAnalysis.git ZZAnalysis
```

To install a complete CMSSW area
------------------------------
Download and execute the setup script for the given release:

*   For V5_15_0/CMSSW_5_3_9 POST-LEGACY:

    ```
    ./checkout_539.csh
    ```


To update this package from the release
------------------------------------------
In the package directory, simply issue
```
git pull
```


To commit and push new changes:
------------------------------
To commit directly (you need write access to the repository):
```
git pull
[edit files]
git add [files to be added]
git commit -m ["commit message"] [files to be added]
git push origin master
```

Otherwise you can make a fork of the repository, develop therein, and make a pull request in the same way as for CMSSW.
