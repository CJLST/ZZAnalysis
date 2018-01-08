#!/usr/bin/env python

"""
This script prints the branching ratio for each production mode and mass point.
Below 300 GeV, where interference effects are important, this is taken from YR3.
At 300 GeV and above, it's taken from the weights that JHUGen writes to the LHE file.
"""

from __future__ import print_function
from utilities import cache, cd, TFile
import argparse, contextlib, csv, functools, itertools, math, numpy, os, re, subprocess, tempfile, urllib
try:
  import yellowhiggs
except ImportError:
  subprocess.check_call("pip install --user -e git+git://github.com/hroskes/yellowhiggs.git@master#egg=yellowhiggs".split())
  import yellowhiggs

basemass = 300


@contextlib.contextmanager
def cd(newdir):
  """http://stackoverflow.com/a/24176022/5228524"""
  prevdir = os.getcwd()
  os.chdir(os.path.expanduser(newdir))
  try:
    yield
  finally:
    os.chdir(prevdir)

class TFile(object):
  def __init__(self, *args, **kwargs):
    self.args, self.kwargs = args, kwargs
  def __enter__(self):
    import ROOT
    self.__tfile = ROOT.TFile.Open(*self.args, **self.kwargs)
    if not self.__tfile or self.__tfile.IsZombie(): return None
    return self.__tfile
  def __exit__(self, *err):
    if self.__tfile:
      self.__tfile.Close()

def cache(function):
  cache = {}
  @functools.wraps(function)
  def newfunction(*args, **kwargs):
    try:
      return cache[args, tuple(sorted(kwargs.iteritems()))]
    except TypeError:
      print(args, tuple(sorted(kwargs.iteritems())))
      raise
    except KeyError:
      cache[args, tuple(sorted(kwargs.iteritems()))] = function(*args, **kwargs)
      return newfunction(*args, **kwargs)
  return newfunction

def BR_YR3(mass):
  try:
    return yellowhiggs.br(mass, "llll(e,mu,tau)")[0]
  except ValueError:
    return .5*(yellowhiggs.br(mass+1, "llll(e,mu,tau)")[0]+yellowhiggs.br(mass-1, "llll(e,mu,tau)")[0])

def GammaHZZ_YR3(mass):
  return yellowhiggs.width(mass)[0] * yellowhiggs.br(mass, "llll(e,mu,tau)")[0]

def sgn(number):
  return math.copysign(1, number)

@cache
def setupJHUGen():
  if not os.path.exists("JHUGen"):
    os.mkdir("JHUGen")
  with open("JHUGen/.gitignore", "w") as f:
    f.write("*")
  if not os.path.exists("JHUGen/JHUGenerator"):
    with cd("JHUGen"):
      subprocess.check_call(["wget", "http://spin.pha.jhu.edu/Generator/JHUGenerator.v7.0.11.tar.gz"])
      subprocess.check_call(["tar" ,"xvzf", "JHUGenerator.v7.0.11.tar.gz"])
  with cd("JHUGen/JHUGenerator"):
    if not os.path.exists("PMZZdistribution.out"):
      subprocess.check_call(["wget", "https://github.com/cms-sw/genproductions/raw/be7e09d66bb79807c5fe12cbb7b2504424f76079/bin/JHUGen/Pdecay/PMZZdistribution.out"])

    changed = False

    with open("makefile") as f:
      makefile = f.read()
    newmakefile = makefile.replace("linkMELA = Yes", "linkMELA = No")
    if newmakefile != makefile:
      changed = True
      with open("makefile", "w") as f:
        f.write(newmakefile)

    with open("mod_PMZZ.F90") as f:
      modPMZZ = f.read()
    newmodPMZZ = modPMZZ.replace(" call HTO_gridHt(EHat/GeV,BigGamma)", " BigGamma = 1\n             !call HTO_gridHt(EHat/GeV,BigGamma)")
    if newmodPMZZ != modPMZZ:
      changed = True
      with open("mod_PMZZ.F90", "w") as f:
        f.write(newmodPMZZ)

    if changed or not os.path.exists("JHUGen"):
      os.system("make")

@cache
def setupBigGamma():
  setupJHUGen()
  with cd("BigGamma"):
    subprocess.check_call(["gfortran", "-o", "BigGamma", "-J", "modules", "BigGamma.F90", "../JHUGen/JHUGenerator/CPS/CALLING_cpHTO.f"])

@cache
def GammaHZZ_JHU(mass):
  setupJHUGen()
  with cd("JHUGen/JHUGenerator"):
    output = subprocess.check_output("./JHUGen ReadPMZZ PrintPMZZ={mass},{mass} ReweightDecay WidthScheme=3 PrintPMZZIntervals=0".format(mass=mass).split())
    for line in output.split("\n"):
      match = re.match(" *([0-9.]+) *([0-9.+-E]+)", line)
      if match and "%" not in line:  #why is % not in line needed?
        if float("{:.3f}".format(float(match.group(1)))) == float("{:.3f}".format(mass)):
          return float(match.group(2))
    raise RuntimeError("Couldn't find {}\n\n\n".format(mass) + output)

@cache
def GammaH_YR2(mass):
  setupBigGamma()
  with cd("BigGamma"):
    return float(subprocess.check_output(["./BigGamma", str(mass)]))

def averageBR(productionmode, mass):
  production = "171217"
  if productionmode == "VBF":
    productionmode = "VBFH"
    production = "171005"
  folder = "{}{:d}".format(productionmode, mass)

  if not os.path.exists("/data3/Higgs"): raise RuntimeError("Have to run this on lxcms03")

  with TFile("/data3/Higgs/"+production+"/"+folder+"/ZZ4lAnalysis.root") as f:
    if not f: return float("nan"), float("nan")
    t = f.ZZTree.Get("candTree")
    failedt = f.ZZTree.Get("candTree_failed")
    if not failedt: failedt = []
    if not t and not failedt: return float("nan"), float("nan")
    for _ in t, failedt:
      if not _: continue
      _.SetBranchStatus("*", 0)
      _.SetBranchStatus("GenHMass", 1)
      _.SetBranchStatus("genHEPMCweight", 1)
      _.SetBranchStatus("p_Gen_CPStoBWPropRewgt", 1)
    t.GetEntry(0)
    multiplyweight = GammaHZZ_YR3(basemass) * GammaHZZ_JHU(t.GenHMass) / (GammaHZZ_JHU(basemass) * abs(t.genHEPMCweight) * GammaH_YR2(t.GenHMass))
    t.GetEntry(1)
#    print("test should be equal:", multiplyweight, GammaHZZ_YR3(basemass) * GammaHZZ_JHU(t.GenHMass) / (GammaHZZ_JHU(basemass) * t.genHEPMCweight * GammaH_YR2(t.GenHMass)))

    bothtrees = itertools.chain(t, failedt)
#    bothtrees = itertools.islice(bothtrees, 5)

    BR, weights, weights_rwttoBW = \
      zip(*([multiplyweight * abs(entry.genHEPMCweight), sgn(entry.genHEPMCweight), sgn(entry.genHEPMCweight)*entry.p_Gen_CPStoBWPropRewgt if productionmode not in ("ZH", "WplusH", "WminusH") else float("nan")] for entry in bothtrees))

    return numpy.average(BR, weights=weights), numpy.average(BR, weights=weights_rwttoBW)

def getmasses(productionmode):
  if productionmode in ("ggH", "VBF", "ZH", "WplusH", "WminusH"):
    return 115, 120, 124, 125, 126, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 190, 200, 210, 230, 250, 270, 300, 350, 400, 450, 500, 550, 600, 700, 750, 800, 900, 1000, 1500, 2000, 2500, 3000
  if productionmode == "ttH":
    return 115, 120, 124, 125, 126, 130, 135, 140, 145

def update_spreadsheet(filename, p, m, BR):
  if p == "VBF": p = "VBFH"
  if (p, m) == ("ggH", 125): requirednmatches = 11
  elif (p, m) == ("ggH", 300): requirednmatches = 6
  elif m == 125: requirednmatches = 5
  else: requirednmatches = 1
  nmatches = 0
  m = str(int(m))
  with tempfile.NamedTemporaryFile(bufsize=0) as newf:
    with open(filename) as f, open(filename) as f2:
      reader = csv.DictReader(f)
      next(f2)
      writer = csv.DictWriter(newf, fieldnames=reader.fieldnames)
      writer.writeheader()
      for row, line in itertools.izip(reader, f2):
        while not line.strip():
          newf.write("\n")
          line = next(f2)
        if p+m in row["identifier"]:
          if re.match("#*"+p+m+"[0-9]+", row["identifier"]):
            pass
          elif re.match("#*"+p+m+"[_scaletuneupdownminloHJJNNLOPS]*", row["identifier"]):
            nmatches += 1
            regex = "(GENBR=)[0-9.eE+-]+(;)"
            match = re.search(regex, row["::variables"])
            assert match, row["::variables"]
            row["::variables"] = re.sub(regex, r"\g<1>{:g}\g<2>".format(BR), row["::variables"])
          else:
            assert False, row["identifier"]
        writer.writerow(row)

      if nmatches != requirednmatches: assert False, (p, m, nmatches, requirednmatches)

    with open(filename, "w") as finalf, open(newf.name) as newf2:
      finalf.write(newf2.read().replace('\r\n', '\n'))


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-p", action="append", choices="ggH VBF ZH WplusH WminusH ttH".split())
  parser.add_argument("--update-spreadsheet")
  args = parser.parse_args()

  line = "{:3} {:4d} {:8.4%}"
  for p in "ggH", "VBF", "ZH", "WplusH", "WminusH", "ttH":
    if args.p and p not in args.p: continue
    for m in getmasses(p):
      if m < 300:
        BR = BR_YR3(min(m, 1000))
      else:
        BR = averageBR(p, m)[0]

      print(line.format(p, m, BR))
      if numpy.isnan(BR): BR = 999
      if args.update_spreadsheet:
        update_spreadsheet(args.update_spreadsheet, p, m, BR)