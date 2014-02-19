#!/bin/tcsh -fe

### merge data
hadd ZZ4lAnalysis_DoubleMu.root ZZ4lAnalysis_DoubleMuA.root ZZ4lAnalysis_DoubleMuB.root
hadd ZZ4lAnalysis_DoubleEle.root ZZ4lAnalysis_DoubleEleA.root ZZ4lAnalysis_DoubleEleB.root
hadd ZZ4lAnalysis_MuEG.root ZZ4lAnalysis_MuEGA.root ZZ4lAnalysis_MuEGB.root
hadd ZZ4lAnalysis_DoubleOr.root ZZ4lAnalysis_DoubleEle.root ZZ4lAnalysis_DoubleMu.root ZZ4lAnalysis_MuEG.root
