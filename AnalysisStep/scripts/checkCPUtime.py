#!/usr/bin/env python

"""
Check folders that exceeded the CPU time limit.
If it looks like everything is there, edit the
output file to say that they didn't exceed the
limit.  This can happen if LSF tries to kill the
job but doesn't manage to do it until the job
finishes.

Because there might be corner cases, and it's
extra important to have all the data, this script
doesn't run on data unless the --force-data flag
is used.
"""

import argparse, glob, gzip, os, re, sys
import ROOT

def checkCPUtime(folder, forcedata=False):
  if not os.path.isdir(folder):
    raise IOError("folder {} does not exist")

  try:
    with open(os.path.join(folder, "exitStatus.txt")) as f:
      contents = f.read()
      if contents.strip() and int(contents):
        return "exited with exit status "+contents.strip()
  except IOError:
    return "no exitStatus.txt"

  try:
    with open(os.path.join(folder, "run_cfg.py")) as f:
      matches = set()
      for line in f:
        match = re.search(r"""\bPD = cms.string[(]["'](.*)["'][)]""")
        if match: matches.add(match.group(1))
      if matches != {""} and not forcedata:
        return "this folder is data, so we want to be extra careful with failed jobs.  If you're sure, you can run with --force-data."
  except IOError:
    return "no run_cfg.py"

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
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('folders', metavar='folder', nargs='*', help='folders to check')
  parser.add_argument('--quiet', '-q', action='store_true', help="don't print anything")
  parser.add_argument('--force-data', action='store_true', help="run even on data folders (dangerous!)")
  args = parser.parse_args()

  if not args.folders: args.folders = glob.glob("*Chunk*/")

  for folder in args.folders:
    result = checkCPUtime(folder, forcedata=args.force_data)
    if not args.quiet:
      print "check CPU time for {}: {}".format(folder, result)
