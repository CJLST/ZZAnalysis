#!/bin/bash

csvfiles=( samples_2015.csv samples_2016_MC.csv samples_2016ICHEP_MC.csv )

for csv in "${csvfiles[@]}"
do
  bash makeFragments_GGPOWHEG.sh $csv
  bash makeFragments_VBFPOWHEG.sh $csv
  bash makeFragments_WHPOWHEG.sh $csv
  bash makeFragments_ZHPOWHEG.sh $csv
  bash makeFragments_ttHPOWHEG.sh $csv
done

bash makeFragments_ACJHUGen.sh Couplings_ACJHUGen.lst