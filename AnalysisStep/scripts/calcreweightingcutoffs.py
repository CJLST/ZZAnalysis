#!/usr/bin/env python
import argparse
import array
from collections import OrderedDict
import numpy
import os
import ROOT

def reweightingcutoffs(directory, percentiles, root=False):

  filename = os.path.join(directory, "ZZ4lAnalysis.root")
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

  if root:
    cutofflists = [[] for p in percentiles]
    for i in range(nweights):
      weightname = "reweightingweights[{}]".format(i)
      hname = "h{}{}".format(directory, i)
      cut = "{} != 0".format(weightname)
      t.Draw("{}>>{}".format(weightname, hname), cut)
      h = getattr(ROOT, hname)
      probsum = array.array('d', [p/100 for p in percentiles])
      q = array.array('d', [0 for p in percentiles])
      h.GetQuantiles(len(percentiles), q, probsum)

      if h.GetRMS() == 0:
        q = [-1 for p in percentiles]

      for cutoff, cutofflist in zip(q, cutofflists):
        cutofflist.append(cutoff)
  else:
    weightarrays = [[] for w in t.reweightingweights]
    for entry in t:
      for w, weightarray in zip(t.reweightingweights, weightarrays):
        if w != 0:
          weightarray.append(w)
    weightarrays = numpy.array([weightarray for weightarray in weightarrays])
    percentile = numpy.percentile(weightarrays, [0]+percentiles+[100], axis=1)
    for i, (minimum, maximum) in enumerate(zip(percentile[0], percentile[-1])):
      if minimum == maximum:
        for cutofflist in percentile:
          cutofflist[i] = -1
    cutofflists = percentile[1:-1]

  return OrderedDict((p, "REWEIGHTING_CUTOFFS={}".format("|".join(str(c) for c in cutofflist)))
                           for p, cutofflist in zip(percentiles, cutofflists))

def putreweightingcutoffsincsvfile(csvfile, basedirectory, percentile, root=False):
  outfile = csvfile.replace(".csv", "_withcutoffs.csv")
  with open(csvfile) as f, open(outfile, "w") as newf:
    header = None
    for line in f:
      data = line.split(",")
      if header is None and line.strip():
        header = line
        variablesindex = data.index("::variables")
      elif line.strip():
        variablessection = data[variablesindex]
        samplename = data[0]
        directory = os.path.join(basedirectory, samplename)
        if "REWEIGHTING_CUTOFFS" not in variablessection and os.path.exists(directory):
          variables = variablessection.split(";")
          variables.append(reweightingcutoffs(directory, [percentile], root)[percentile])
          data[variablesindex] = ";".join(v for v in variables if v)
          line = ",".join(data)
      print line,
      newf.write(line)

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Find the quantiles of each of the reweightingweights.  Example if you run from AAAOK: calcreweightingcutoffs.py ggH0PM_M125 99.99.  Should not be run on a sample that was already produced with a cutoff.")
  parser.add_argument("directory", help="if no csvfile is provided, a directory containing ZZ4lAnalysis.root.  If a csvfile is provided, the directory containing those directories (e.g. AAAOK).")
  parser.add_argument("percentile", nargs="+", type=float)
  parser.add_argument("--root", action="store_true", help="fill a histogram with TTree::Draw and get the quantiles using TH1::GetQuantiles.  About 3x faster but less precise than using numpy which is the default.")
  parser.add_argument("--csvfile")
  args = parser.parse_args()

  if args.csvfile:
    if len(args.percentile) > 1:
      raise ValueError("Can only use one percentile to add cutoffs in a csv file!")
    putreweightingcutoffsincsvfile(args.csvfile, args.directory, args.percentile[0], args.root)
  else:
    for k, v in reweightingcutoffs(args.directory, args.percentile, args.root).iteritems():
      print "{:.3f}%: {}".format(k, v)
