# Example of a fragment to be added in the "pyFragments" field in the csv. 

setConf("testvar", 1) # set or replace "testvar" 
setConf("testlist", 2, append=True) # set "testlist" as a list, append if already existing
setConf("testlist", 3, append=True)
