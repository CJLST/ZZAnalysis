#!/usr/bin/env python

"""
This script prints the branching ratio for each production mode and mass point.
Below 300 GeV, where interference effects are important, this is taken from YR3.
At 300 GeV and above, it's taken from the weights that JHUGen writes to the LHE file.
"""

from __future__ import print_function
import argparse, contextlib, csv, functools, itertools, math, numpy, os, re, subprocess, tempfile, urllib
try:
  import yellowhiggs
except ImportError:
  subprocess.check_call("pip install --user -e git+git://github.com/hroskes/yellowhiggs.git@master#egg=yellowhiggs".split())
  import yellowhiggs
if os.path.exists("src") and not os.path.exists("src/.gitignore"):
  with open("src/.gitignore", "w") as f:
    f.write("yellowhiggs\n.gitignore\npip-delete-this-directory.txt")

basemass = 300
oldfilterefficiencyZH = 0.15038
filterefficiencyZH = 0.1483
filterefficiencyttH = 0.1544
CJLSTproduction = "180121"

YR4data_BR_4l = {
  120: 0.0001659,
  124: 0.0002502,
  125: 0.0002745,
  126: 0.0003001,
  130: 0.0004124,
}

YR4data_BR_ZZ = {
  120: 0.01572,
  124: 0.02383,
  125: 0.02619,
  126: 0.02866,
  130: 0.03955,
}

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

def BR_YR3(mass, productionmode):
  if productionmode in ("ZH", "ttH") and mass < 300:
    result = BR_fromspreadsheet(productionmode, mass)
    if 120 <= m <= 130:
      result *= YR4data_BR_ZZ[m] / yellowhiggs.br(m, "ZZ")[0]
    if productionmode == "ZH":
      result *= filterefficiencyZH / oldfilterefficiencyZH
    return result

  mass = int(str(mass))
  try:
    result = yellowhiggs.br(mass, "llll(e,mu,tau)")[0]
  except ValueError:
    return .5*(BR_YR3(mass+1, productionmode)+BR_YR3(mass-1, productionmode))

  if productionmode in ("ggH", "VBF", "WplusH", "WminusH"):
    if 120 <= m <= 130:
      result *= YR4data_BR_4l[m] / yellowhiggs.br(m, "llll(e,mu,tau)")[0]
    return result

  if 120 <= m <= 130: assert False

  #need to get ZZ2l2x BR
  #from the comment here: https://github.com/CJLST/ZZAnalysis/blob/454fffb5842f470e71e89348a6de7a8d30f4a813/AnalysisStep/test/prod/samples_2016_MC.csv#L5-L6
  #it doesn't say how to get BR(H->ZZ*->2L2nu) or BR(H->WW*->2L2nu)
  #llnunu  = ZZllnunu + WWlnulnu
  #        = ZZllnunu + 9*WWemununu
  #        = ZZllnunu + 9*emununu

  result += yellowhiggs.br(mass, "llqq(e,mu,tau)")[0] + yellowhiggs.br(mass, "llnunu(e,mu,tau)")[0] - 9*yellowhiggs.br(mass, "enuemunumu")[0]
  if productionmode == "ZH":
    result *= filterefficiencyZH
  elif productionmode == "ttH":
    result *= filterefficiencyttH
  else:
    assert False, productionmode

  

  return result

def GammaHZZ_YR3(mass, productionmode):
  return yellowhiggs.width(mass)[0] * BR_YR3(mass, productionmode)

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

def averageBR(productionmode, mass, spreadsheet=None):
  production = CJLSTproduction
  if productionmode == "VBF":
    productionmode = "VBFH"
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
    multiplyweight = GammaHZZ_YR3(basemass, p) * GammaHZZ_JHU(t.GenHMass) / (GammaHZZ_JHU(basemass) * abs(t.genHEPMCweight) * GammaH_YR2(t.GenHMass))
    t.GetEntry(1)
#    print("test should be equal:", multiplyweight, GammaHZZ_YR3(basemass, p) * GammaHZZ_JHU(t.GenHMass) / (GammaHZZ_JHU(basemass) * t.genHEPMCweight * GammaH_YR2(t.GenHMass)))

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
        if p+str(m) in row["identifier"]:
          if re.match("#*"+p+str(m)+"[0-9]+", row["identifier"]):
            pass
          elif re.match("#*"+p+str(m)+"[_scaletuneupdownminloHJJNNLOPS]*", row["identifier"]):
            nmatches += 1
            regex = "(GENBR=)[0-9.eE+-]+(;)"
            match = re.search(regex, row["::variables"])
            assert match, row["::variables"]
            row["::variables"] = re.sub(regex, r"\g<1>{:g}\g<2>".format(BR), row["::variables"])
            if p in ("ZH", "ttH") and int(m) >= 300:
              row["BR=1"] = BR_YR3(min(m, 1000), p)
            elif p in ("ZH", "ttH"):
              row["BR=1"] = BR
              if "_scale" in row["identifier"] or "_tune" in row["identifier"]:
                if p == "ZH": row["BR=1"] /= filterefficiencyZH
                if p == "ttH": row["BR=1"] /= filterefficiencyttH
            elif 120 <= m <= 130:
              row["BR=1"] = YR4data_BR_4l[m]
          else:
            assert False, row["identifier"]
        writer.writerow(row)

      if nmatches != requirednmatches: assert False, (p, m, nmatches, requirednmatches)

    with open(filename, "w") as finalf, open(newf.name) as newf2:
      finalf.write(newf2.read().replace('\r\n', '\n'))

def BR_fromspreadsheet(p, m):
  with contextlib.closing(urllib.urlopen("https://raw.githubusercontent.com/CJLST/ZZAnalysis/87bce99aa845936454ce302a70095797b81c194c/AnalysisStep/test/prod/samples_2016_MC.csv")) as f:
    reader = csv.DictReader(f)
    for row in reader:
      if re.match("#*"+p+str(m), row["identifier"]):
        return float(row["BR=1"])

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-p", action="append", choices="ggH VBF ZH WplusH WminusH ttH".split())
  parser.add_argument("--minimum-mass", type=float, default=0)
  parser.add_argument("--maximum-mass", type=float, default=float("inf"))
  parser.add_argument("--update-spreadsheet")
  args = parser.parse_args()

  line = "{:3} {:4d} {:8.4%}"
  for p in "ggH", "VBF", "ZH", "WplusH", "WminusH", "ttH":
    if args.p and p not in args.p: continue
    for m in getmasses(p):
      if not args.minimum_mass <= m <= args.maximum_mass: continue
      if m < 300:
        BR = BR_YR3(m, p)
      else:
        BR = averageBR(p, m)[0]

      print(line.format(p, m, BR))
      if numpy.isnan(BR): BR = 999
      if args.update_spreadsheet:
        update_spreadsheet(args.update_spreadsheet, p, m, BR)
