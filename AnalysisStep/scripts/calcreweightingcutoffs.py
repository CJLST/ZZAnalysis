#!/usr/bin/env python
import argparse
import array
import numpy
import os
import ROOT

parser = argparse.ArgumentParser(description="Find the quantiles of each of the reweightingweights.  Example if you run from AAAOK: calcreweightingcutoffs.py ggH0PM_M125 99.99.  Should not be run on a sample that was already produced with a cutoff.")
parser.add_argument("directory")
parser.add_argument("percentile", nargs="+", type=float)
group = parser.add_mutually_exclusive_group()
group.add_argument("--root", action="store_true", help="fill a histogram with TTree::Draw and get the quantiles using TH1::GetQuantiles.  About 3x faster but less precise than using numpy which is the default.")
args = parser.parse_args()

filename = os.path.join(args.directory, "ZZ4lAnalysis.root")
f = ROOT.TFile(filename)
if not f:
    raise IOError("No file {}".format(filename))
t = f.Get("ZZTree/candTree")
if not t:
    raise IOError("No tree ZZTree/candTree in {}".format(filename))

t.GetEntry(0)
try:
    nweights = len(t.reweightingweights)
except AttributeError:
    raise IOError("Tree does not have reweightingweights.")

if args.root:
    cutofflists = [[] for p in args.percentile]
    for i in range(nweights):
        weightname = "reweightingweights[{}]".format(i)
        hname = "h{}".format(i)
        t.Draw("{}>>{}".format(weightname, hname))
        h = getattr(ROOT, hname)
        probsum = array.array('d', [p/100 for p in args.percentile])
        q = array.array('d', [0 for p in args.percentile])
        h.GetQuantiles(len(args.percentile), q, probsum)

        if h.GetRMS() == 0:
            q = [-1 for p in args.percentile]

        for cutoff, cutofflist in zip(q, cutofflists):
            cutofflist.append(cutoff)
else:
    weightarrays = [[] for w in t.reweightingweights]
    for entry in t:
        for w, weightarray in zip(t.reweightingweights, weightarrays):
            weightarray.append(w)
    weightarrays = numpy.array([weightarray for weightarray in weightarrays])
    percentile = numpy.percentile(weightarrays, [0]+args.percentile+[100], axis=1)
    for i, (minimum, maximum) in enumerate(zip(percentile[0], percentile[-1])):
        if minimum == maximum:
            for cutofflist in percentile:
                cutofflist[i] = -1
    cutofflists = percentile[1:-1]

for p, cutofflist in zip(args.percentile, cutofflists):
    print "{:.3f}%: REWEIGHTING_CUTOFFS={}".format(p, "|".join(str(c) for c in cutofflist))
