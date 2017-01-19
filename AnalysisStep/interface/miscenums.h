//This header is meant to contain enums that can also be set
//in python via the variables column of the csv file.
//Any enum defined here will also be a module variable
//of python/miscenums.py

#ifndef MISCENUMS_H
#define MISCENUMS_H

enum FailedTreeLevel {
    //only used if skipEmptyEvents == false
    noFailedTree = 0,      //events with no candidate are skipped entirely
    minimalFailedTree = 1, //events with no candidate are written to a separate tree with minimal information: Gen Higgs momentum, lepton flavors, and LHE level ME's
    LHEFailedTree = 2,     //same as (1) + full LHE level information
    fullFailedTree = 3,    //same as (2) + full pythia level information + reco information that is available without a candidate, e.g. jets
};

#endif
