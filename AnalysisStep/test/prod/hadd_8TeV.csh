#!/bin/tcsh -fe

### Merge ZZ samples
hadd ZZ4lAnalysis_ZZTo4mu.root     ZZ4lAnalysis_ZZ4mu.root     ZZ4lAnalysis_ZZ4mu_ext.root
hadd ZZ4lAnalysis_ZZTo4e.root      ZZ4lAnalysis_ZZ4e.root      ZZ4lAnalysis_ZZ4e_ext.root
# Skip the buggy ZZ2e2tau_ext sample
#hadd ZZ4lAnalysis_ZZTo2e2tau.root  ZZ4lAnalysis_ZZ2e2tau.root  ZZ4lAnalysis_ZZ2e2tau_ext.root
ln -s ZZ4lAnalysis_ZZ2e2tau.root ZZ4lAnalysis_ZZTo2e2tau.root 
hadd ZZ4lAnalysis_ZZTo2mu2tau.root ZZ4lAnalysis_ZZ2mu2tau.root ZZ4lAnalysis_ZZ2mu2tau_ext.root  
hadd ZZ4lAnalysis_ZZTo2e2mu.root   ZZ4lAnalysis_ZZ2e2mu.root   ZZ4lAnalysis_ZZ2e2mu_ext.root  
hadd ZZ4lAnalysis_ZZTo4tau.root    ZZ4lAnalysis_ZZ4tau.root    ZZ4lAnalysis_ZZ4tau_ext.root 


### merge data

hadd ZZ4lAnalysis_DoubleMu_1963.root ZZ4lAnalysis_DoubleMuA.root ZZ4lAnalysis_DoubleMuB.root ZZ4lAnalysis_DoubleMuC.root ZZ4lAnalysis_DoubleMuD.root 
hadd ZZ4lAnalysis_DoubleEle_1963.root ZZ4lAnalysis_DoubleEleA.root ZZ4lAnalysis_DoubleEleB.root ZZ4lAnalysis_DoubleEleC.root ZZ4lAnalysis_DoubleEleD.root
hadd ZZ4lAnalysis_MuEG_1963.root ZZ4lAnalysis_MuEGA.root ZZ4lAnalysis_MuEGB.root ZZ4lAnalysis_MuEGC.root ZZ4lAnalysis_MuEGD.root
hadd ZZ4lAnalysis_DoubleOr_1963.root ZZ4lAnalysis_DoubleMu_1963.root ZZ4lAnalysis_DoubleEle_1963.root ZZ4lAnalysis_MuEG_1963.root
