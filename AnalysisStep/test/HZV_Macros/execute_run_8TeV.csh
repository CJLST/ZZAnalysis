#!/bin/csh

set prodname = rootuples/$1/PRODFSR_8TeV/
set dirname = trees/$2/PRODFSR_8TeV/

if( -e $dirname ) then
   rm $dirname/data/*
   rm $dirname/4mu/*
   rm $dirname/4e/*
   rm $dirname/2mu2e/*
   rm $dirname/CR/*
else
   mkdir trees/$2
   mkdir $dirname
   mkdir $dirname/data
   mkdir $dirname/4mu
   mkdir $dirname/4e
   mkdir $dirname/2mu2e
   mkdir $dirname/CR
endif


##data

./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_DoubleMu_1963.root $dirname/data/HZZ4lTree_DoubleMu.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_DoubleEle_1963.root $dirname/data/HZZ4lTree_DoubleEle.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_DoubleOr_1963.root $dirname/data/HZZ4lTree_DoubleOr.root

#./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_DoubleMu_1210.root $dirname/data/HZZ4lTree_DoubleMu_1210.root
#./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_DoubleEle_1210.root $dirname/data/HZZ4lTree_DoubleEle_1210.root
#./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_DoubleOr_1210.root $dirname/data/HZZ4lTree_DoubleOr_1210.root

#ichep data
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_DoubleMu_5300.root $dirname/data/HZZ4lTree_DoubleMu_5300.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_DoubleEle_5300.root $dirname/data/HZZ4lTree_DoubleEle_5300.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_DoubleOr_5300.root $dirname/data/HZZ4lTree_DoubleOr_5300.root


##4mu
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H115.root $dirname/4mu/HZZ4lTree_H115.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H116.root $dirname/4mu/HZZ4lTree_H116.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H117.root $dirname/4mu/HZZ4lTree_H117.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H118.root $dirname/4mu/HZZ4lTree_H118.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H119.root $dirname/4mu/HZZ4lTree_H119.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H120.root $dirname/4mu/HZZ4lTree_H120.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H121.root $dirname/4mu/HZZ4lTree_H121.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H122.root $dirname/4mu/HZZ4lTree_H122.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H123.root $dirname/4mu/HZZ4lTree_H123.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H124.root $dirname/4mu/HZZ4lTree_H124.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H125.root $dirname/4mu/HZZ4lTree_H125.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H126.root $dirname/4mu/HZZ4lTree_H126.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H127.root $dirname/4mu/HZZ4lTree_H127.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H128.root $dirname/4mu/HZZ4lTree_H128.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H129.root $dirname/4mu/HZZ4lTree_H129.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H130.root $dirname/4mu/HZZ4lTree_H130.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H135.root $dirname/4mu/HZZ4lTree_H135.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H140.root $dirname/4mu/HZZ4lTree_H140.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H145.root $dirname/4mu/HZZ4lTree_H145.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H150.root $dirname/4mu/HZZ4lTree_H150.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H160.root $dirname/4mu/HZZ4lTree_H160.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H170.root $dirname/4mu/HZZ4lTree_H170.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H180.root $dirname/4mu/HZZ4lTree_H180.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H190.root $dirname/4mu/HZZ4lTree_H190.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H200.root $dirname/4mu/HZZ4lTree_H200.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H220.root $dirname/4mu/HZZ4lTree_H220.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H250.root $dirname/4mu/HZZ4lTree_H250.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H275.root $dirname/4mu/HZZ4lTree_H275.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H300.root $dirname/4mu/HZZ4lTree_H300.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H325.root $dirname/4mu/HZZ4lTree_H325.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H350.root $dirname/4mu/HZZ4lTree_H350.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H375.root $dirname/4mu/HZZ4lTree_H375.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H400.root $dirname/4mu/HZZ4lTree_H400.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H425.root $dirname/4mu/HZZ4lTree_H425.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H450.root $dirname/4mu/HZZ4lTree_H450.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H475.root $dirname/4mu/HZZ4lTree_H475.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H500.root $dirname/4mu/HZZ4lTree_H500.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H525.root $dirname/4mu/HZZ4lTree_H525.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H550.root $dirname/4mu/HZZ4lTree_H550.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H575.root $dirname/4mu/HZZ4lTree_H575.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H600.root $dirname/4mu/HZZ4lTree_H600.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H650.root $dirname/4mu/HZZ4lTree_H650.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H700.root $dirname/4mu/HZZ4lTree_H700.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H750.root $dirname/4mu/HZZ4lTree_H750.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H800.root $dirname/4mu/HZZ4lTree_H800.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H850.root $dirname/4mu/HZZ4lTree_H850.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H900.root $dirname/4mu/HZZ4lTree_H900.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H950.root $dirname/4mu/HZZ4lTree_H950.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_H1000.root $dirname/4mu/HZZ4lTree_H1000.root

./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_ggZZ4l.root $dirname/4mu/HZZ4lTree_ggZZ4l.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_ggZZ2l2l.root $dirname/4mu/HZZ4lTree_ggZZ2l2l.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_ZZTo4mu.root $dirname/4mu/HZZ4lTree_ZZTo4mu.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_ZZTo4e.root $dirname/4mu/HZZ4lTree_ZZTo4e.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_ZZTo2e2mu.root $dirname/4mu/HZZ4lTree_ZZTo2e2mu.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_ZZTo2mu2tau.root $dirname/4mu/HZZ4lTree_ZZTo2mu2tau.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_ZZTo2e2tau.root $dirname/4mu/HZZ4lTree_ZZTo2e2tau.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_ZZTo4tau.root $dirname/4mu/HZZ4lTree_ZZTo4tau.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_ZZJetsTo4L.root $dirname/4mu/HZZ4lTree_ZZJetsTo4L.root

./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_DYJetsToLLTuneZ2M50B.root $dirname/4mu/HZZ4lTree_DYJetsToLLTuneZ2M50B.root
./run_HZZ4l 0 1 $prodname/ZZ4lAnalysis_DYJetsToLLTuneZ2M50NoB.root $dirname/4mu/HZZ4lTree_DYJetsToLLTuneZ2M50NoB.root


##4e

./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H115.root $dirname/4e/HZZ4lTree_H115.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H116.root $dirname/4e/HZZ4lTree_H116.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H117.root $dirname/4e/HZZ4lTree_H117.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H118.root $dirname/4e/HZZ4lTree_H118.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H119.root $dirname/4e/HZZ4lTree_H119.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H120.root $dirname/4e/HZZ4lTree_H120.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H121.root $dirname/4e/HZZ4lTree_H121.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H122.root $dirname/4e/HZZ4lTree_H122.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H123.root $dirname/4e/HZZ4lTree_H123.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H124.root $dirname/4e/HZZ4lTree_H124.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H125.root $dirname/4e/HZZ4lTree_H125.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H126.root $dirname/4e/HZZ4lTree_H126.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H127.root $dirname/4e/HZZ4lTree_H127.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H128.root $dirname/4e/HZZ4lTree_H128.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H129.root $dirname/4e/HZZ4lTree_H129.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H130.root $dirname/4e/HZZ4lTree_H130.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H135.root $dirname/4e/HZZ4lTree_H135.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H140.root $dirname/4e/HZZ4lTree_H140.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H145.root $dirname/4e/HZZ4lTree_H145.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H150.root $dirname/4e/HZZ4lTree_H150.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H160.root $dirname/4e/HZZ4lTree_H160.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H170.root $dirname/4e/HZZ4lTree_H170.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H180.root $dirname/4e/HZZ4lTree_H180.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H190.root $dirname/4e/HZZ4lTree_H190.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H200.root $dirname/4e/HZZ4lTree_H200.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H220.root $dirname/4e/HZZ4lTree_H220.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H250.root $dirname/4e/HZZ4lTree_H250.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H275.root $dirname/4e/HZZ4lTree_H275.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H300.root $dirname/4e/HZZ4lTree_H300.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H325.root $dirname/4e/HZZ4lTree_H325.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H350.root $dirname/4e/HZZ4lTree_H350.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H375.root $dirname/4e/HZZ4lTree_H375.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H400.root $dirname/4e/HZZ4lTree_H400.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H425.root $dirname/4e/HZZ4lTree_H425.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H450.root $dirname/4e/HZZ4lTree_H450.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H475.root $dirname/4e/HZZ4lTree_H475.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H500.root $dirname/4e/HZZ4lTree_H500.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H525.root $dirname/4e/HZZ4lTree_H525.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H550.root $dirname/4e/HZZ4lTree_H550.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H575.root $dirname/4e/HZZ4lTree_H575.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H600.root $dirname/4e/HZZ4lTree_H600.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H650.root $dirname/4e/HZZ4lTree_H650.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H700.root $dirname/4e/HZZ4lTree_H700.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H750.root $dirname/4e/HZZ4lTree_H750.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H800.root $dirname/4e/HZZ4lTree_H800.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H850.root $dirname/4e/HZZ4lTree_H850.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H900.root $dirname/4e/HZZ4lTree_H900.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H950.root $dirname/4e/HZZ4lTree_H950.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_H1000.root $dirname/4e/HZZ4lTree_H1000.root


./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_ggZZ4l.root $dirname/4e/HZZ4lTree_ggZZ4l.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_ggZZ2l2l.root $dirname/4e/HZZ4lTree_ggZZ2l2l.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_ZZTo4mu.root $dirname/4e/HZZ4lTree_ZZTo4mu.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_ZZTo4e.root $dirname/4e/HZZ4lTree_ZZTo4e.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_ZZTo2e2mu.root $dirname/4e/HZZ4lTree_ZZTo2e2mu.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_ZZTo2mu2tau.root $dirname/4e/HZZ4lTree_ZZTo2mu2tau.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_ZZTo2e2tau.root $dirname/4e/HZZ4lTree_ZZTo2e2tau.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_ZZTo4tau.root $dirname/4e/HZZ4lTree_ZZTo4tau.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_ZZJetsTo4L.root $dirname/4e/HZZ4lTree_ZZJetsTo4L.root

./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_DYJetsToLLTuneZ2M50B.root $dirname/4e/HZZ4lTree_DYJetsToLLTuneZ2M50B.root
./run_HZZ4l 1 1 $prodname/ZZ4lAnalysis_DYJetsToLLTuneZ2M50NoB.root $dirname/4e/HZZ4lTree_DYJetsToLLTuneZ2M50NoB.root


##2mu2e
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H115.root $dirname/2mu2e/HZZ4lTree_H115.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H116.root $dirname/2mu2e/HZZ4lTree_H116.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H117.root $dirname/2mu2e/HZZ4lTree_H117.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H118.root $dirname/2mu2e/HZZ4lTree_H118.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H119.root $dirname/2mu2e/HZZ4lTree_H119.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H120.root $dirname/2mu2e/HZZ4lTree_H120.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H121.root $dirname/2mu2e/HZZ4lTree_H121.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H122.root $dirname/2mu2e/HZZ4lTree_H122.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H123.root $dirname/2mu2e/HZZ4lTree_H123.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H124.root $dirname/2mu2e/HZZ4lTree_H124.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H125.root $dirname/2mu2e/HZZ4lTree_H125.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H126.root $dirname/2mu2e/HZZ4lTree_H126.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H127.root $dirname/2mu2e/HZZ4lTree_H127.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H128.root $dirname/2mu2e/HZZ4lTree_H128.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H129.root $dirname/2mu2e/HZZ4lTree_H129.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H130.root $dirname/2mu2e/HZZ4lTree_H130.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H135.root $dirname/2mu2e/HZZ4lTree_H135.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H140.root $dirname/2mu2e/HZZ4lTree_H140.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H145.root $dirname/2mu2e/HZZ4lTree_H145.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H150.root $dirname/2mu2e/HZZ4lTree_H150.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H160.root $dirname/2mu2e/HZZ4lTree_H160.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H170.root $dirname/2mu2e/HZZ4lTree_H170.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H180.root $dirname/2mu2e/HZZ4lTree_H180.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H190.root $dirname/2mu2e/HZZ4lTree_H190.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H200.root $dirname/2mu2e/HZZ4lTree_H200.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H220.root $dirname/2mu2e/HZZ4lTree_H220.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H250.root $dirname/2mu2e/HZZ4lTree_H250.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H275.root $dirname/2mu2e/HZZ4lTree_H275.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H300.root $dirname/2mu2e/HZZ4lTree_H300.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H325.root $dirname/2mu2e/HZZ4lTree_H325.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H350.root $dirname/2mu2e/HZZ4lTree_H350.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H375.root $dirname/2mu2e/HZZ4lTree_H375.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H400.root $dirname/2mu2e/HZZ4lTree_H400.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H425.root $dirname/2mu2e/HZZ4lTree_H425.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H450.root $dirname/2mu2e/HZZ4lTree_H450.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H475.root $dirname/2mu2e/HZZ4lTree_H475.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H500.root $dirname/2mu2e/HZZ4lTree_H500.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H525.root $dirname/2mu2e/HZZ4lTree_H525.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H550.root $dirname/2mu2e/HZZ4lTree_H550.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H575.root $dirname/2mu2e/HZZ4lTree_H575.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H600.root $dirname/2mu2e/HZZ4lTree_H600.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H650.root $dirname/2mu2e/HZZ4lTree_H650.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H700.root $dirname/2mu2e/HZZ4lTree_H700.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H750.root $dirname/2mu2e/HZZ4lTree_H750.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H800.root $dirname/2mu2e/HZZ4lTree_H800.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H850.root $dirname/2mu2e/HZZ4lTree_H850.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H900.root $dirname/2mu2e/HZZ4lTree_H900.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H950.root $dirname/2mu2e/HZZ4lTree_H950.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_H1000.root $dirname/2mu2e/HZZ4lTree_H1000.root


./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_ggZZ4l.root $dirname/2mu2e/HZZ4lTree_ggZZ4l.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_ggZZ2l2l.root $dirname/2mu2e/HZZ4lTree_ggZZ2l2l.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_ZZTo4mu.root $dirname/2mu2e/HZZ4lTree_ZZTo4mu.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_ZZTo4e.root $dirname/2mu2e/HZZ4lTree_ZZTo4e.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_ZZTo2e2mu.root $dirname/2mu2e/HZZ4lTree_ZZTo2e2mu.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_ZZTo2mu2tau.root $dirname/2mu2e/HZZ4lTree_ZZTo2mu2tau.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_ZZTo2e2tau.root $dirname/2mu2e/HZZ4lTree_ZZTo2e2tau.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_ZZTo4tau.root $dirname/2mu2e/HZZ4lTree_ZZTo4tau.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_ZZJetsTo4L.root $dirname/2mu2e/HZZ4lTree_ZZJetsTo4L.root

./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_DYJetsToLLTuneZ2M50B.root $dirname/2mu2e/HZZ4lTree_DYJetsToLLTuneZ2M50B.root
./run_HZZ4l 2 1 $prodname/ZZ4lAnalysis_DYJetsToLLTuneZ2M50NoB.root $dirname/2mu2e/HZZ4lTree_DYJetsToLLTuneZ2M50NoB.root

##CR
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H115.root $dirname/CR/HZZ4lTree_H115.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H116.root $dirname/CR/HZZ4lTree_H116.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H117.root $dirname/CR/HZZ4lTree_H117.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H118.root $dirname/CR/HZZ4lTree_H118.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H119.root $dirname/CR/HZZ4lTree_H119.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H120.root $dirname/CR/HZZ4lTree_H120.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H121.root $dirname/CR/HZZ4lTree_H121.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H122.root $dirname/CR/HZZ4lTree_H122.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H123.root $dirname/CR/HZZ4lTree_H123.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H124.root $dirname/CR/HZZ4lTree_H124.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H125.root $dirname/CR/HZZ4lTree_H125.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H126.root $dirname/CR/HZZ4lTree_H126.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H127.root $dirname/CR/HZZ4lTree_H127.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H128.root $dirname/CR/HZZ4lTree_H128.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H129.root $dirname/CR/HZZ4lTree_H129.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H130.root $dirname/CR/HZZ4lTree_H130.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H135.root $dirname/CR/HZZ4lTree_H135.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H140.root $dirname/CR/HZZ4lTree_H140.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H145.root $dirname/CR/HZZ4lTree_H145.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H150.root $dirname/CR/HZZ4lTree_H150.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H160.root $dirname/CR/HZZ4lTree_H160.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H170.root $dirname/CR/HZZ4lTree_H170.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H180.root $dirname/CR/HZZ4lTree_H180.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H190.root $dirname/CR/HZZ4lTree_H190.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H200.root $dirname/CR/HZZ4lTree_H200.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H220.root $dirname/CR/HZZ4lTree_H220.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H250.root $dirname/CR/HZZ4lTree_H250.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H275.root $dirname/CR/HZZ4lTree_H275.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H300.root $dirname/CR/HZZ4lTree_H300.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H325.root $dirname/CR/HZZ4lTree_H325.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H350.root $dirname/CR/HZZ4lTree_H350.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H375.root $dirname/CR/HZZ4lTree_H375.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H400.root $dirname/CR/HZZ4lTree_H400.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H425.root $dirname/CR/HZZ4lTree_H425.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H450.root $dirname/CR/HZZ4lTree_H450.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H475.root $dirname/CR/HZZ4lTree_H475.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H500.root $dirname/CR/HZZ4lTree_H500.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H525.root $dirname/CR/HZZ4lTree_H525.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H550.root $dirname/CR/HZZ4lTree_H550.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H575.root $dirname/CR/HZZ4lTree_H575.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H600.root $dirname/CR/HZZ4lTree_H600.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H650.root $dirname/CR/HZZ4lTree_H650.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H700.root $dirname/CR/HZZ4lTree_H700.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H750.root $dirname/CR/HZZ4lTree_H750.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H800.root $dirname/CR/HZZ4lTree_H800.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H850.root $dirname/CR/HZZ4lTree_H850.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H900.root $dirname/CR/HZZ4lTree_H900.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H950.root $dirname/CR/HZZ4lTree_H950.root
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_H1000.root $dirname/CR/HZZ4lTree_H1000.root


./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_ggZZ4l.root $dirname/CR/HZZ4lTree_ggZZ4l
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_ggZZ2l2l.root $dirname/CR/HZZ4lTree_ggZZ2l2l
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_ZZTo4mu.root $dirname/CR/HZZ4lTree_ZZTo4mu
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_ZZTo4e.root $dirname/CR/HZZ4lTree_ZZTo4e
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_ZZTo2e2mu.root $dirname/CR/HZZ4lTree_ZZTo2e2mu
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_ZZTo2mu2tau.root $dirname/CR/HZZ4lTree_ZZTo2mu2tau
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_ZZTo2e2tau.root $dirname/CR/HZZ4lTree_ZZTo2e2tau
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_ZZTo4tau.root $dirname/CR/HZZ4lTree_ZZTo4tau
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_ZZJetsTo4L.root $dirname/CR/HZZ4lTree_ZZJetsTo4L.root

./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_DYJetsToLLTuneZ2M10B.root $dirname/CR/HZZ4lTree_DYJetsToLLTuneZ2M10B
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_DYJetsToLLTuneZ2M10NoB.root $dirname/CR/HZZ4lTree_DYJetsToLLTuneZ2M10NoB
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_DYJetsToLLTuneZ2M50B.root $dirname/CR/HZZ4lTree_DYJetsToLLTuneZ2M50B
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_DYJetsToLLTuneZ2M50NoB.root $dirname/CR/HZZ4lTree_DYJetsToLLTuneZ2M50NoB

./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_DoubleMu_1963.root $dirname/CR/HZZ4lTree_DoubleMu
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_DoubleEle_1963.root $dirname/CR/HZZ4lTree_DoubleEle
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_DoubleOr_1963.root $dirname/CR/HZZ4lTree_DoubleOr

./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_DoubleMu_5300.root $dirname/CR/HZZ4lTree_DoubleMu_5300
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_DoubleEle_5300.root $dirname/CR/HZZ4lTree_DoubleEle_5300
./run_HZZ4l_CR 1 $prodname/ZZ4lAnalysis_DoubleOr_5300.root $dirname/CR/HZZ4lTree_DoubleOr_5300
