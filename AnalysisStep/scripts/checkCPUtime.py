#!/usr/bin/env python

import argparse, glob, gzip, os, re, sys
import ROOT

def checkCPUtime(folder):
  if not os.path.isdir(folder):
    raise IOError("folder {} does not exist")

  try:
    with open(os.path.join(folder, "exitStatus.txt")) as f:
      contents = f.read()
      if contents.strip() and int(contents):
        return "exited with exit status "+contents.strip()
  except IOError:
    return "no exitStatus.txt"

  jobouts = sorted(glob.glob(os.path.join(folder, "job_*.txt")))
  if not jobouts: return
  jobout = jobouts[-1]

  with open(jobout) as f:
    contents = f.read()
    if "CPU time limit [NOT] exceeded" in contents: return "CPU time limit is reported as exceeded, but we already saw this file and changed it to [NOT]"
    if "CPU time limit exceeded" not in contents: return "CPU time limit not exceeded"  #nothing more to do

  try:
    with gzip.GzipFile(os.path.join(folder, "log.txt.gz")) as f:
      contents = f.read()
      if "Begin Fatal Exception" in contents:
        return "log has fatal exception"
      match = re.search("\nNevt_Gen: *([0-9]*)\n", contents)
      if not match:
        return "can't find nevents in log"
      nevents = int(match.group(1))
  except IOError:
    return "no log.txt.gz"

  f = ROOT.TFile(os.path.join(folder, "ZZ4lAnalysis.root"))
  if not f: return "no ZZ4lAnalysis.root"
  if f.TestBit(ROOT.TFile.kRecovered): return "ZZ4lAnalysis.root is recovered (probably not all there)"
  if not f.ZZTree: return "no ZZTree"
  t = f.ZZTree.Get("candTree")
  if not t: return "no candTree"
  t2 = f.ZZTree.Get("candTree_failed")
  if t2:
    if t.GetEntries() + t2.GetEntries() != nevents: return "not all events recorded in candTree or candTree_failed"

  with open(jobout) as f:
    contents = f.read()
  contents = contents.replace("CPU time limit exceeded", "CPU time limit [NOT] exceeded")
  with open(jobout, "w") as f:
    f.write(contents)

  return "CPU time limit is reported as exceeded, but looks like everything is there :).  Marking it as [NOT] exceeded."

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Check folders that exceeded the CPU time limit.')
  parser.add_argument('folders', metavar='folder', nargs='*', help='folders to check')
  parser.add_argument('--quiet', '-q', action='store_true')
  args = parser.parse_args()

  if not args.folders: args.folders = glob.glob("*Chunk*/")

  for folder in args.folders:
    result = checkCPUtime(folder)
    if not args.quiet:
      print "check CPU time for {}: {}".format(folder, result)
