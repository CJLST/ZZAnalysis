#!/bin/bash

##
# This is a script to loop over a list of csv files to make POWHEG pyFragments out of the tpl files present,
# and then makes the fragments for the anomalous couplings JHUGen samples from the lst file.
# One needs to update the csvfiles list to maintain the most comprehensive set of fragments.
# Temporary note specific to mH=125 GeV: The anomalous coupling probabilities from MCFM are currently commented out in the tpl, so one needs to comment them back in manually.
##

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