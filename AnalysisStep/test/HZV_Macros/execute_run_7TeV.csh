#!/bin/csh

set prodname = rootuples/$1/PRODFSR/
set dirname = trees/$2/PRODFSR/

if( -e $dirname ) then
   rm $dirname/data/*
   rm $dirname/4mu/*
   rm $dirname/4e/*
   rm $dirname/2mu2e/*
   rm $dirname/CR/*
else
   mkdir trees
   mkdir trees/$2
   mkdir $dirname
   mkdir $dirname/data
   mkdir $dirname/4mu
   mkdir $dirname/4e
   mkdir $dirname/2mu2e
   mkdir $dirname/CR
endif


##data

./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_DoubleMu_NewJSON.root $dirname/data/HZZ4lTree_DoubleMu.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_DoubleEle_NewJSON.root $dirname/data/HZZ4lTree_DoubleEle.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_DoubleOr_NewJSON.root $dirname/data/HZZ4lTree_DoubleOr.root


##4mu

./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH125.root $dirname/4mu/HZZ4lTree_VBFH125.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH1000.root $dirname/4mu/HZZ4lTree_VBFH1000.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH115.root $dirname/4mu/HZZ4lTree_VBFH115.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH120.root $dirname/4mu/HZZ4lTree_VBFH120.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH130.root $dirname/4mu/HZZ4lTree_VBFH130.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH140.root $dirname/4mu/HZZ4lTree_VBFH140.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH150.root $dirname/4mu/HZZ4lTree_VBFH150.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH160.root $dirname/4mu/HZZ4lTree_VBFH160.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH170.root $dirname/4mu/HZZ4lTree_VBFH170.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH180.root $dirname/4mu/HZZ4lTree_VBFH180.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH190.root $dirname/4mu/HZZ4lTree_VBFH190.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH200.root $dirname/4mu/HZZ4lTree_VBFH200.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH210.root $dirname/4mu/HZZ4lTree_VBFH210.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH220.root $dirname/4mu/HZZ4lTree_VBFH220.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH230.root $dirname/4mu/HZZ4lTree_VBFH230.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH250.root $dirname/4mu/HZZ4lTree_VBFH250.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH275.root $dirname/4mu/HZZ4lTree_VBFH275.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH300.root $dirname/4mu/HZZ4lTree_VBFH300.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH325.root $dirname/4mu/HZZ4lTree_VBFH325.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH350.root $dirname/4mu/HZZ4lTree_VBFH350.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH375.root $dirname/4mu/HZZ4lTree_VBFH375.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH400.root $dirname/4mu/HZZ4lTree_VBFH400.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH425.root $dirname/4mu/HZZ4lTree_VBFH425.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH450.root $dirname/4mu/HZZ4lTree_VBFH450.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH475.root $dirname/4mu/HZZ4lTree_VBFH475.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH500.root $dirname/4mu/HZZ4lTree_VBFH500.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH575.root $dirname/4mu/HZZ4lTree_VBFH575.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH600.root $dirname/4mu/HZZ4lTree_VBFH600.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH650.root $dirname/4mu/HZZ4lTree_VBFH650.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH700.root $dirname/4mu/HZZ4lTree_VBFH700.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH800.root $dirname/4mu/HZZ4lTree_VBFH800.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH900.root $dirname/4mu/HZZ4lTree_VBFH900.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VBFH950.root $dirname/4mu/HZZ4lTree_VBFH950.root

./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH115.root $dirname/4mu/HZZ4lTree_VH115.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH120.root $dirname/4mu/HZZ4lTree_VH120.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH130.root $dirname/4mu/HZZ4lTree_VH130.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH135.root $dirname/4mu/HZZ4lTree_VH135.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH140.root $dirname/4mu/HZZ4lTree_VH140.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH150.root $dirname/4mu/HZZ4lTree_VH150.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH160.root $dirname/4mu/HZZ4lTree_VH160.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH170.root $dirname/4mu/HZZ4lTree_VH170.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH180.root $dirname/4mu/HZZ4lTree_VH180.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH190.root $dirname/4mu/HZZ4lTree_VH190.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH200.root $dirname/4mu/HZZ4lTree_VH200.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH210.root $dirname/4mu/HZZ4lTree_VH210.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH220.root $dirname/4mu/HZZ4lTree_VH220.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH250.root $dirname/4mu/HZZ4lTree_VH250.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH275.root $dirname/4mu/HZZ4lTree_VH275.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH300.root $dirname/4mu/HZZ4lTree_VH300.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH325.root $dirname/4mu/HZZ4lTree_VH325.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH350.root $dirname/4mu/HZZ4lTree_VH350.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH375.root $dirname/4mu/HZZ4lTree_VH375.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH425.root $dirname/4mu/HZZ4lTree_VH425.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH450.root $dirname/4mu/HZZ4lTree_VH450.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH475.root $dirname/4mu/HZZ4lTree_VH475.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH525.root $dirname/4mu/HZZ4lTree_VH525.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH550.root $dirname/4mu/HZZ4lTree_VH550.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH575.root $dirname/4mu/HZZ4lTree_VH575.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_VH600.root $dirname/4mu/HZZ4lTree_VH600.root

./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H120.root $dirname/4mu/HZZ4lTree_H120.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H124.root $dirname/4mu/HZZ4lTree_H124.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H125.root $dirname/4mu/HZZ4lTree_H125.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H126.root $dirname/4mu/HZZ4lTree_H126.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H130.root $dirname/4mu/HZZ4lTree_H130.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H140.root $dirname/4mu/HZZ4lTree_H140.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H150.root $dirname/4mu/HZZ4lTree_H150.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H160.root $dirname/4mu/HZZ4lTree_H160.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H170.root $dirname/4mu/HZZ4lTree_H170.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H180.root $dirname/4mu/HZZ4lTree_H180.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H190.root $dirname/4mu/HZZ4lTree_H190.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H200.root $dirname/4mu/HZZ4lTree_H200.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H210.root $dirname/4mu/HZZ4lTree_H210.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H220.root $dirname/4mu/HZZ4lTree_H220.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H250.root $dirname/4mu/HZZ4lTree_H250.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H275.root $dirname/4mu/HZZ4lTree_H275.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H300.root $dirname/4mu/HZZ4lTree_H300.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H325.root $dirname/4mu/HZZ4lTree_H325.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H350.root $dirname/4mu/HZZ4lTree_H350.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H400.root $dirname/4mu/HZZ4lTree_H400.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H425.root $dirname/4mu/HZZ4lTree_H425.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H450.root $dirname/4mu/HZZ4lTree_H450.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H475.root $dirname/4mu/HZZ4lTree_H475.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H525.root $dirname/4mu/HZZ4lTree_H525.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H550.root $dirname/4mu/HZZ4lTree_H550.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H575.root $dirname/4mu/HZZ4lTree_H575.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H600.root $dirname/4mu/HZZ4lTree_H600.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H650.root $dirname/4mu/HZZ4lTree_H650.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H700.root $dirname/4mu/HZZ4lTree_H700.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H750.root $dirname/4mu/HZZ4lTree_H750.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H800.root $dirname/4mu/HZZ4lTree_H800.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H900.root $dirname/4mu/HZZ4lTree_H900.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H950.root $dirname/4mu/HZZ4lTree_H950.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_H1000.root $dirname/4mu/HZZ4lTree_H1000.root

./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_ggZZ4l.root $dirname/4mu/HZZ4lTree_ggZZ4l.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_ggZZ2l2l.root $dirname/4mu/HZZ4lTree_ggZZ2l2l.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_ZZTo4mu.root $dirname/4mu/HZZ4lTree_ZZTo4mu.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_ZZTo4e.root $dirname/4mu/HZZ4lTree_ZZTo4e.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_ZZTo2e2mu.root $dirname/4mu/HZZ4lTree_ZZTo2e2mu.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_ZZTo2mu2tau.root $dirname/4mu/HZZ4lTree_ZZTo2mu2tau.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_ZZTo2e2tau.root $dirname/4mu/HZZ4lTree_ZZTo2e2tau.root
./run_HZZ4l 0 0 $prodname/ZZ4lAnalysis_ZZTo4tau.root $dirname/4mu/HZZ4lTree_ZZTo4tau.root


##4e

./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH125.root $dirname/4e/HZZ4lTree_VBFH125.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH1000.root $dirname/4e/HZZ4lTree_VBFH1000.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH115.root $dirname/4e/HZZ4lTree_VBFH115.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH120.root $dirname/4e/HZZ4lTree_VBFH120.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH130.root $dirname/4e/HZZ4lTree_VBFH130.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH140.root $dirname/4e/HZZ4lTree_VBFH140.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH150.root $dirname/4e/HZZ4lTree_VBFH150.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH160.root $dirname/4e/HZZ4lTree_VBFH160.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH170.root $dirname/4e/HZZ4lTree_VBFH170.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH180.root $dirname/4e/HZZ4lTree_VBFH180.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH190.root $dirname/4e/HZZ4lTree_VBFH190.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH200.root $dirname/4e/HZZ4lTree_VBFH200.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH210.root $dirname/4e/HZZ4lTree_VBFH210.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH220.root $dirname/4e/HZZ4lTree_VBFH220.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH230.root $dirname/4e/HZZ4lTree_VBFH230.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH250.root $dirname/4e/HZZ4lTree_VBFH250.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH275.root $dirname/4e/HZZ4lTree_VBFH275.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH300.root $dirname/4e/HZZ4lTree_VBFH300.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH325.root $dirname/4e/HZZ4lTree_VBFH325.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH350.root $dirname/4e/HZZ4lTree_VBFH350.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH375.root $dirname/4e/HZZ4lTree_VBFH375.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH400.root $dirname/4e/HZZ4lTree_VBFH400.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH425.root $dirname/4e/HZZ4lTree_VBFH425.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH450.root $dirname/4e/HZZ4lTree_VBFH450.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH475.root $dirname/4e/HZZ4lTree_VBFH475.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH500.root $dirname/4e/HZZ4lTree_VBFH500.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH575.root $dirname/4e/HZZ4lTree_VBFH575.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH600.root $dirname/4e/HZZ4lTree_VBFH600.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH650.root $dirname/4e/HZZ4lTree_VBFH650.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH700.root $dirname/4e/HZZ4lTree_VBFH700.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH800.root $dirname/4e/HZZ4lTree_VBFH800.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH900.root $dirname/4e/HZZ4lTree_VBFH900.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VBFH950.root $dirname/4e/HZZ4lTree_VBFH950.root

./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH115.root $dirname/4e/HZZ4lTree_VH115.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH120.root $dirname/4e/HZZ4lTree_VH120.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH130.root $dirname/4e/HZZ4lTree_VH130.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH135.root $dirname/4e/HZZ4lTree_VH135.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH140.root $dirname/4e/HZZ4lTree_VH140.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH150.root $dirname/4e/HZZ4lTree_VH150.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH160.root $dirname/4e/HZZ4lTree_VH160.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH170.root $dirname/4e/HZZ4lTree_VH170.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH180.root $dirname/4e/HZZ4lTree_VH180.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH190.root $dirname/4e/HZZ4lTree_VH190.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH200.root $dirname/4e/HZZ4lTree_VH200.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH210.root $dirname/4e/HZZ4lTree_VH210.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH220.root $dirname/4e/HZZ4lTree_VH220.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH250.root $dirname/4e/HZZ4lTree_VH250.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH275.root $dirname/4e/HZZ4lTree_VH275.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH300.root $dirname/4e/HZZ4lTree_VH300.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH325.root $dirname/4e/HZZ4lTree_VH325.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH350.root $dirname/4e/HZZ4lTree_VH350.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH375.root $dirname/4e/HZZ4lTree_VH375.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH425.root $dirname/4e/HZZ4lTree_VH425.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH450.root $dirname/4e/HZZ4lTree_VH450.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH475.root $dirname/4e/HZZ4lTree_VH475.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH525.root $dirname/4e/HZZ4lTree_VH525.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH550.root $dirname/4e/HZZ4lTree_VH550.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH575.root $dirname/4e/HZZ4lTree_VH575.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_VH600.root $dirname/4e/HZZ4lTree_VH600.root

./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H120.root $dirname/4e/HZZ4lTree_H120.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H124.root $dirname/4e/HZZ4lTree_H124.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H125.root $dirname/4e/HZZ4lTree_H125.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H126.root $dirname/4e/HZZ4lTree_H126.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H130.root $dirname/4e/HZZ4lTree_H130.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H140.root $dirname/4e/HZZ4lTree_H140.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H150.root $dirname/4e/HZZ4lTree_H150.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H160.root $dirname/4e/HZZ4lTree_H160.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H170.root $dirname/4e/HZZ4lTree_H170.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H180.root $dirname/4e/HZZ4lTree_H180.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H190.root $dirname/4e/HZZ4lTree_H190.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H200.root $dirname/4e/HZZ4lTree_H200.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H210.root $dirname/4e/HZZ4lTree_H210.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H220.root $dirname/4e/HZZ4lTree_H220.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H250.root $dirname/4e/HZZ4lTree_H250.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H275.root $dirname/4e/HZZ4lTree_H275.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H300.root $dirname/4e/HZZ4lTree_H300.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H325.root $dirname/4e/HZZ4lTree_H325.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H350.root $dirname/4e/HZZ4lTree_H350.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H400.root $dirname/4e/HZZ4lTree_H400.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H425.root $dirname/4e/HZZ4lTree_H425.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H450.root $dirname/4e/HZZ4lTree_H450.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H475.root $dirname/4e/HZZ4lTree_H475.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H525.root $dirname/4e/HZZ4lTree_H525.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H550.root $dirname/4e/HZZ4lTree_H550.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H575.root $dirname/4e/HZZ4lTree_H575.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H600.root $dirname/4e/HZZ4lTree_H600.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H650.root $dirname/4e/HZZ4lTree_H650.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H700.root $dirname/4e/HZZ4lTree_H700.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H750.root $dirname/4e/HZZ4lTree_H750.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H800.root $dirname/4e/HZZ4lTree_H800.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H900.root $dirname/4e/HZZ4lTree_H900.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H950.root $dirname/4e/HZZ4lTree_H950.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_H1000.root $dirname/4e/HZZ4lTree_H1000.root


./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_ggZZ4l.root $dirname/4e/HZZ4lTree_ggZZ4l.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_ggZZ2l2l.root $dirname/4e/HZZ4lTree_ggZZ2l2l.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_ZZTo4mu.root $dirname/4e/HZZ4lTree_ZZTo4mu.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_ZZTo4e.root $dirname/4e/HZZ4lTree_ZZTo4e.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_ZZTo2e2mu.root $dirname/4e/HZZ4lTree_ZZTo2e2mu.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_ZZTo2mu2tau.root $dirname/4e/HZZ4lTree_ZZTo2mu2tau.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_ZZTo2e2tau.root $dirname/4e/HZZ4lTree_ZZTo2e2tau.root
./run_HZZ4l 1 0 $prodname/ZZ4lAnalysis_ZZTo4tau.root $dirname/4e/HZZ4lTree_ZZTo4tau.root


##2mu2e

./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH125.root $dirname/2mu2e/HZZ4lTree_VBFH125.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH1000.root $dirname/2mu2e/HZZ4lTree_VBFH1000.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH115.root $dirname/2mu2e/HZZ4lTree_VBFH115.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH120.root $dirname/2mu2e/HZZ4lTree_VBFH120.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH130.root $dirname/2mu2e/HZZ4lTree_VBFH130.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH140.root $dirname/2mu2e/HZZ4lTree_VBFH140.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH150.root $dirname/2mu2e/HZZ4lTree_VBFH150.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH160.root $dirname/2mu2e/HZZ4lTree_VBFH160.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH170.root $dirname/2mu2e/HZZ4lTree_VBFH170.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH180.root $dirname/2mu2e/HZZ4lTree_VBFH180.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH190.root $dirname/2mu2e/HZZ4lTree_VBFH190.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH200.root $dirname/2mu2e/HZZ4lTree_VBFH200.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH210.root $dirname/2mu2e/HZZ4lTree_VBFH210.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH220.root $dirname/2mu2e/HZZ4lTree_VBFH220.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH230.root $dirname/2mu2e/HZZ4lTree_VBFH230.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH250.root $dirname/2mu2e/HZZ4lTree_VBFH250.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH275.root $dirname/2mu2e/HZZ4lTree_VBFH275.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH300.root $dirname/2mu2e/HZZ4lTree_VBFH300.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH325.root $dirname/2mu2e/HZZ4lTree_VBFH325.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH350.root $dirname/2mu2e/HZZ4lTree_VBFH350.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH375.root $dirname/2mu2e/HZZ4lTree_VBFH375.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH400.root $dirname/2mu2e/HZZ4lTree_VBFH400.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH425.root $dirname/2mu2e/HZZ4lTree_VBFH425.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH450.root $dirname/2mu2e/HZZ4lTree_VBFH450.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH475.root $dirname/2mu2e/HZZ4lTree_VBFH475.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH500.root $dirname/2mu2e/HZZ4lTree_VBFH500.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH575.root $dirname/2mu2e/HZZ4lTree_VBFH575.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH600.root $dirname/2mu2e/HZZ4lTree_VBFH600.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH650.root $dirname/2mu2e/HZZ4lTree_VBFH650.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH700.root $dirname/2mu2e/HZZ4lTree_VBFH700.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH800.root $dirname/2mu2e/HZZ4lTree_VBFH800.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH900.root $dirname/2mu2e/HZZ4lTree_VBFH900.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VBFH950.root $dirname/2mu2e/HZZ4lTree_VBFH950.root

./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH115.root $dirname/2mu2e/HZZ4lTree_VH115.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH120.root $dirname/2mu2e/HZZ4lTree_VH120.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH130.root $dirname/2mu2e/HZZ4lTree_VH130.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH135.root $dirname/2mu2e/HZZ4lTree_VH135.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH140.root $dirname/2mu2e/HZZ4lTree_VH140.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH150.root $dirname/2mu2e/HZZ4lTree_VH150.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH160.root $dirname/2mu2e/HZZ4lTree_VH160.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH170.root $dirname/2mu2e/HZZ4lTree_VH170.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH180.root $dirname/2mu2e/HZZ4lTree_VH180.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH190.root $dirname/2mu2e/HZZ4lTree_VH190.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH200.root $dirname/2mu2e/HZZ4lTree_VH200.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH210.root $dirname/2mu2e/HZZ4lTree_VH210.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH220.root $dirname/2mu2e/HZZ4lTree_VH220.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH250.root $dirname/2mu2e/HZZ4lTree_VH250.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH275.root $dirname/2mu2e/HZZ4lTree_VH275.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH300.root $dirname/2mu2e/HZZ4lTree_VH300.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH325.root $dirname/2mu2e/HZZ4lTree_VH325.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH350.root $dirname/2mu2e/HZZ4lTree_VH350.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH375.root $dirname/2mu2e/HZZ4lTree_VH375.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH425.root $dirname/2mu2e/HZZ4lTree_VH425.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH450.root $dirname/2mu2e/HZZ4lTree_VH450.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH475.root $dirname/2mu2e/HZZ4lTree_VH475.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH525.root $dirname/2mu2e/HZZ4lTree_VH525.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH550.root $dirname/2mu2e/HZZ4lTree_VH550.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH575.root $dirname/2mu2e/HZZ4lTree_VH575.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_VH600.root $dirname/2mu2e/HZZ4lTree_VH600.root
							   
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H120.root $dirname/2mu2e/HZZ4lTree_H120.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H124.root $dirname/2mu2e/HZZ4lTree_H124.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H125.root $dirname/2mu2e/HZZ4lTree_H125.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H126.root $dirname/2mu2e/HZZ4lTree_H126.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H130.root $dirname/2mu2e/HZZ4lTree_H130.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H140.root $dirname/2mu2e/HZZ4lTree_H140.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H150.root $dirname/2mu2e/HZZ4lTree_H150.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H160.root $dirname/2mu2e/HZZ4lTree_H160.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H170.root $dirname/2mu2e/HZZ4lTree_H170.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H180.root $dirname/2mu2e/HZZ4lTree_H180.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H190.root $dirname/2mu2e/HZZ4lTree_H190.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H200.root $dirname/2mu2e/HZZ4lTree_H200.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H210.root $dirname/2mu2e/HZZ4lTree_H210.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H220.root $dirname/2mu2e/HZZ4lTree_H220.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H250.root $dirname/2mu2e/HZZ4lTree_H250.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H275.root $dirname/2mu2e/HZZ4lTree_H275.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H300.root $dirname/2mu2e/HZZ4lTree_H300.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H325.root $dirname/2mu2e/HZZ4lTree_H325.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H350.root $dirname/2mu2e/HZZ4lTree_H350.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H400.root $dirname/2mu2e/HZZ4lTree_H400.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H425.root $dirname/2mu2e/HZZ4lTree_H425.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H450.root $dirname/2mu2e/HZZ4lTree_H450.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H475.root $dirname/2mu2e/HZZ4lTree_H475.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H525.root $dirname/2mu2e/HZZ4lTree_H525.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H550.root $dirname/2mu2e/HZZ4lTree_H550.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H575.root $dirname/2mu2e/HZZ4lTree_H575.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H600.root $dirname/2mu2e/HZZ4lTree_H600.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H650.root $dirname/2mu2e/HZZ4lTree_H650.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H700.root $dirname/2mu2e/HZZ4lTree_H700.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H750.root $dirname/2mu2e/HZZ4lTree_H750.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H800.root $dirname/2mu2e/HZZ4lTree_H800.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H900.root $dirname/2mu2e/HZZ4lTree_H900.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H950.root $dirname/2mu2e/HZZ4lTree_H950.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_H1000.root $dirname/2mu2e/HZZ4lTree_H1000.root

./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_ggZZ4l.root $dirname/2mu2e/HZZ4lTree_ggZZ4l.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_ggZZ2l2l.root $dirname/2mu2e/HZZ4lTree_ggZZ2l2l.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_ZZTo4mu.root $dirname/2mu2e/HZZ4lTree_ZZTo4mu.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_ZZTo4e.root $dirname/2mu2e/HZZ4lTree_ZZTo4e.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_ZZTo2e2mu.root $dirname/2mu2e/HZZ4lTree_ZZTo2e2mu.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_ZZTo2mu2tau.root $dirname/2mu2e/HZZ4lTree_ZZTo2mu2tau.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_ZZTo2e2tau.root $dirname/2mu2e/HZZ4lTree_ZZTo2e2tau.root
./run_HZZ4l 2 0 $prodname/ZZ4lAnalysis_ZZTo4tau.root $dirname/2mu2e/HZZ4lTree_ZZTo4tau.root


##CR

./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH125.root $dirname/CR/HZZ4lTree_VBFH125.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH1000.root $dirname/CR/HZZ4lTree_VBFH1000.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH120.root $dirname/CR/HZZ4lTree_VBFH120.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH130.root $dirname/CR/HZZ4lTree_VBFH130.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH140.root $dirname/CR/HZZ4lTree_VBFH140.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH150.root $dirname/CR/HZZ4lTree_VBFH150.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH160.root $dirname/CR/HZZ4lTree_VBFH160.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH170.root $dirname/CR/HZZ4lTree_VBFH170.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH180.root $dirname/CR/HZZ4lTree_VBFH180.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH190.root $dirname/CR/HZZ4lTree_VBFH190.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH200.root $dirname/CR/HZZ4lTree_VBFH200.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH210.root $dirname/CR/HZZ4lTree_VBFH210.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH220.root $dirname/CR/HZZ4lTree_VBFH220.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH230.root $dirname/CR/HZZ4lTree_VBFH230.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH250.root $dirname/CR/HZZ4lTree_VBFH250.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH275.root $dirname/CR/HZZ4lTree_VBFH275.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH300.root $dirname/CR/HZZ4lTree_VBFH300.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH325.root $dirname/CR/HZZ4lTree_VBFH325.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH350.root $dirname/CR/HZZ4lTree_VBFH350.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH375.root $dirname/CR/HZZ4lTree_VBFH375.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH400.root $dirname/CR/HZZ4lTree_VBFH400.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH425.root $dirname/CR/HZZ4lTree_VBFH425.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH450.root $dirname/CR/HZZ4lTree_VBFH450.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH475.root $dirname/CR/HZZ4lTree_VBFH475.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH500.root $dirname/CR/HZZ4lTree_VBFH500.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH575.root $dirname/CR/HZZ4lTree_VBFH575.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH600.root $dirname/CR/HZZ4lTree_VBFH600.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH650.root $dirname/CR/HZZ4lTree_VBFH650.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH700.root $dirname/CR/HZZ4lTree_VBFH700.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH800.root $dirname/CR/HZZ4lTree_VBFH800.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH900.root $dirname/CR/HZZ4lTree_VBFH900.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VBFH950.root $dirname/CR/HZZ4lTree_VBFH950.root


./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH115.root $dirname/CR/HZZ4lTree_VH115.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH120.root $dirname/CR/HZZ4lTree_VH120.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH130.root $dirname/CR/HZZ4lTree_VH130.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH135.root $dirname/CR/HZZ4lTree_VH135.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH140.root $dirname/CR/HZZ4lTree_VH140.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH150.root $dirname/CR/HZZ4lTree_VH150.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH160.root $dirname/CR/HZZ4lTree_VH160.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH170.root $dirname/CR/HZZ4lTree_VH170.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH180.root $dirname/CR/HZZ4lTree_VH180.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH190.root $dirname/CR/HZZ4lTree_VH190.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH200.root $dirname/CR/HZZ4lTree_VH200.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH210.root $dirname/CR/HZZ4lTree_VH210.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH220.root $dirname/CR/HZZ4lTree_VH220.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH250.root $dirname/CR/HZZ4lTree_VH250.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH275.root $dirname/CR/HZZ4lTree_VH275.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH300.root $dirname/CR/HZZ4lTree_VH300.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH325.root $dirname/CR/HZZ4lTree_VH325.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH350.root $dirname/CR/HZZ4lTree_VH350.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH375.root $dirname/CR/HZZ4lTree_VH375.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH425.root $dirname/CR/HZZ4lTree_VH425.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH450.root $dirname/CR/HZZ4lTree_VH450.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH475.root $dirname/CR/HZZ4lTree_VH475.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH525.root $dirname/CR/HZZ4lTree_VH525.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH550.root $dirname/CR/HZZ4lTree_VH550.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH575.root $dirname/CR/HZZ4lTree_VH575.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_VH600.root $dirname/CR/HZZ4lTree_VH600.root

./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H120.root $dirname/CR/HZZ4lTree_H120.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H124.root $dirname/CR/HZZ4lTree_H124.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H125.root $dirname/CR/HZZ4lTree_H125.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H126.root $dirname/CR/HZZ4lTree_H126.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H130.root $dirname/CR/HZZ4lTree_H130.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H140.root $dirname/CR/HZZ4lTree_H140.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H150.root $dirname/CR/HZZ4lTree_H150.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H160.root $dirname/CR/HZZ4lTree_H160.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H170.root $dirname/CR/HZZ4lTree_H170.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H180.root $dirname/CR/HZZ4lTree_H180.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H190.root $dirname/CR/HZZ4lTree_H190.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H200.root $dirname/CR/HZZ4lTree_H200.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H210.root $dirname/CR/HZZ4lTree_H210.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H220.root $dirname/CR/HZZ4lTree_H220.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H250.root $dirname/CR/HZZ4lTree_H250.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H275.root $dirname/CR/HZZ4lTree_H275.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H300.root $dirname/CR/HZZ4lTree_H300.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H325.root $dirname/CR/HZZ4lTree_H325.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H350.root $dirname/CR/HZZ4lTree_H350.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H400.root $dirname/CR/HZZ4lTree_H400.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H425.root $dirname/CR/HZZ4lTree_H425.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H450.root $dirname/CR/HZZ4lTree_H450.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H475.root $dirname/CR/HZZ4lTree_H475.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H525.root $dirname/CR/HZZ4lTree_H525.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H550.root $dirname/CR/HZZ4lTree_H550.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H575.root $dirname/CR/HZZ4lTree_H575.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H600.root $dirname/CR/HZZ4lTree_H600.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H650.root $dirname/CR/HZZ4lTree_H650.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H700.root $dirname/CR/HZZ4lTree_H700.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H750.root $dirname/CR/HZZ4lTree_H750.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H800.root $dirname/CR/HZZ4lTree_H800.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H900.root $dirname/CR/HZZ4lTree_H900.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H950.root $dirname/CR/HZZ4lTree_H950.root
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_H1000.root $dirname/CR/HZZ4lTree_H1000.root


./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_ggZZ4l.root $dirname/CR/HZZ4lTree_ggZZ4l
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_ggZZ2l2l.root $dirname/CR/HZZ4lTree_ggZZ2l2l
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_ZZTo4mu.root $dirname/CR/HZZ4lTree_ZZTo4mu
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_ZZTo4e.root $dirname/CR/HZZ4lTree_ZZTo4e
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_ZZTo2e2mu.root $dirname/CR/HZZ4lTree_ZZTo2e2mu
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_ZZTo2mu2tau.root $dirname/CR/HZZ4lTree_ZZTo2mu2tau
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_ZZTo2e2tau.root $dirname/CR/HZZ4lTree_ZZTo2e2tau
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_ZZTo4tau.root $dirname/CR/HZZ4lTree_ZZTo4tau
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_DYJetsToLLTuneZ2M10B.root $dirname/CR/HZZ4lTree_DYJetsToLLTuneZ2M10B
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_DYJetsToLLTuneZ2M10NoB.root $dirname/CR/HZZ4lTree_DYJetsToLLTuneZ2M10NoB
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_DYJetsToLLTuneZ2M50B.root $dirname/CR/HZZ4lTree_DYJetsToLLTuneZ2M50B
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_DYJetsToLLTuneZ2M50NoB.root $dirname/CR/HZZ4lTree_DYJetsToLLTuneZ2M50NoB

./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_DoubleMu_NewJSON.root $dirname/CR/HZZ4lTree_DoubleMu
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_DoubleEle_NewJSON.root $dirname/CR/HZZ4lTree_DoubleEle
./run_HZZ4l_CR 0 $prodname/ZZ4lAnalysis_DoubleOr_NewJSON.root $dirname/CR/HZZ4lTree_DoubleOr
