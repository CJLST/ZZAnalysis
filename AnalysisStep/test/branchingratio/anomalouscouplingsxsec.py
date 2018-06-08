#!/usr/bin/env python

import argparse, csv, itertools, os, re, subprocess, tempfile

from utilities import cache, cd

if not os.path.exists("anomalouscouplingsconstants"):
  subprocess.check_call(["git", "clone", "git@github.com:hroskes/anomalouscouplingsconstants"])

commit = "36d952846228fc7bb555ed935c7b2954eb634fe8"

with cd("anomalouscouplingsconstants"):
  if subprocess.check_output(["git", "rev-parse", "HEAD"]).strip() != commit:
    subprocess.check_call(["git", "fetch", "--all"])
    subprocess.check_call(["git", "checkout", commit])

import anomalouscouplingsconstants

from uncertainties import nominal_value, ufloat

def genxsecfromrow(row):
  return float(row["::variables"].split("GENXSEC=")[1].split(";")[0])
def genBRfromrow(row):
  return float(row["::variables"].split("GENBR=")[1].split(";")[0])

"""
These numbers are from 2017.
However, they only depend on the Z branching fractions
and lepton interference, neither of which depend on the PDF
or anything like that.
"""
filterefficiencyZHpowheg = ufloat(0.14836, 0.00018)
filterefficienciesZH = {
  "0L1":            ufloat(0.1490, 0.0011),
  "0L1Zg":          ufloat(0.1877, 0.0012),
  "0L1Zgf05ph0":    ufloat(0.1502, 0.0011),
  "0L1f05ph0":      ufloat(0.1520, 0.0012),
  "0PM":            ufloat(0.1520, 0.0011),
  "0PH":            ufloat(0.1472, 0.0011),
  "0PHf05ph0":      ufloat(0.1529, 0.0011),
  "0M":             ufloat(0.1483, 0.0011),
  "0Mf05ph0":       ufloat(0.1502, 0.0011),
}
filterefficiencyttHpowheg = ufloat(0.1541, 0.0004)
filterefficienciesttH = {
  "0M":             ufloat(0.1320, 0.0011),
  "0Mf05ph0":       ufloat(0.1338, 0.0011),
  "0PM":            ufloat(0.1339, 0.0011),
}

@cache
def SMgenxsecBR(year, process):
  if process == "WH": return SMgenxsec(year, "WplusH") + SMgenxsec(year, "WminusH")
  if process == "HJJ": process = "ggH"
  if process == "VBF": process += "H"
  process += "125"

  SMfilename =  os.path.join(os.path.dirname(__file__), "../prod/samples_{}_MC.csv".format(year))
  with open(SMfilename) as f:
    for row in csv.DictReader(f):
      if row["identifier"] == process: return genxsecfromrow(row), genBRfromrow(row)
  raise ValueError("Didn't find " + process)

def SMgenxsec(year, process):
  if process == "WH": return SMgenxsec(year, "WplusH") + SMgenxsec(year, "WminusH")
  return SMgenxsecBR(year, process)[0]
def SMgenBR(year, process):
  if process == "WH":
    assert SMgenBR(year, "WplusH") == SMgenBR(year, "WminusH")
    return SMgenBR(year, "WplusH")
  return SMgenxsecBR(year, process)[1]

def JHUxsec(year, process, coupling):
  if process == "VBFH": process = "VBF"

  gdecayname = {
    "0M": "g4",
    "0Mf05ph0": "g4",
    "0PH": "g2",
    "0PHf05ph0": "g2",
    "0L1": "g1prime2",
    "0L1f05ph0": "g1prime2",
    "0L1Zg": "ghzgs1prime2",
    "0L1Zgf05ph0": "ghzgs1prime2",
    "0PM": "g4",  #doesn't actually matter but have to set it to something
  }[coupling]

  if process in ("ggH", "VBF", "ZH", "WH"):
    decaycoupling = {
      "0PM":         {"a1": 0                        },
      "0M":          {                      "a3":   2},
      "0Mf05ph0":    {"a1": 0, "a1a3":   1, "a3":   2},
      "0PH":         {                      "a2":   2},
      "0PHf05ph0":   {"a1": 0, "a1a2":   1, "a2":   2},
      "0L1":         {                      "L1":   2},
      "0L1f05ph0":   {"a1": 0, "a1L1":   1, "L1":   2},
      "0L1Zg":       {                      "L1Zg": 2},
      "0L1Zgf05ph0": {"a1": 0, "a1L1Zg": 1, "L1Zg": 2},
    }[coupling]
  else:
    decaycoupling = {"a1": 0}

  if process == "ggH":
    prodcoupling = None
  elif process in ("VBF", "ZH", "WH"):
    prodSMcoupling = "a1"
    prodcoupling = decaycoupling
    gprodname = gdecayname
  elif process == "HJJ":
    prodSMcoupling = "a2"
    prodcoupling = {
      "0PM":         {"a2": 0                    },
      "0M":          {                    "a3": 2},
      "0Mf05ph0":    {"a2": 0, "a2a3": 1, "a3": 2},
    }[coupling]
    gprodname = "ghg4"
  elif process == "ttH":
    prodSMcoupling = "kappa"
    prodcoupling = {
      "0PM":         {"kappa": 0                                       },
      "0M":          {                                  "kappatilde": 2},
      "0Mf05ph0":    {"kappa": 0, "kappakappatilde": 1, "kappatilde": 2},
    }[coupling]
    gprodname = "kappa_tilde_"

  gdecay = getattr(anomalouscouplingsconstants, gdecayname+"HZZ")
  decaypart = sum(gdecay**power * getattr(anomalouscouplingsconstants, "JHUXSHZZ4l"+xs) for xs, power in decaycoupling.iteritems()) / anomalouscouplingsconstants.JHUXSHZZ4la1
  if process == "ggH":
    if year == 2016:
      prodpart = 0.9901
    else:
      prodpart = 1
  else:
    pdf = {2016: "NNPDF30_lo_as_0130", 2017: "NNPDF31_lo_as_0130"}[year]
    pdfmodule = getattr(anomalouscouplingsconstants, pdf)
    gprod = getattr(anomalouscouplingsconstants, gprodname+process)
    prodpart = sum(gprod**power * getattr(pdfmodule, "JHUXS"+process+xs) for xs, power in prodcoupling.iteritems()) / getattr(pdfmodule, "JHUXS"+process+prodSMcoupling)

  return prodpart * decaypart

def update_spreadsheet(year):
  filename =  os.path.join(os.path.dirname(__file__), "../prod/samples_{}_MC_anomalous.csv".format(year))
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
        identifier = row["identifier"].lstrip("#")
        if identifier.strip() and not identifier.strip().startswith("ggGrav"):
          process, coupling = re.match("(.*H|HJJ)(0.*)_M125", identifier).groups()
          genxsec = SMgenxsec(year, process) * JHUxsec(year, process, coupling)
          regex = "(GENXSEC=)[0-9.eE+-]+(;)"
          match = re.search(regex, row["::variables"])
          assert match, row["::variables"]
          row["::variables"] = re.sub(regex, r"\g<1>{:g}\g<2>".format(genxsec.nominal_value), row["::variables"])

          genBR = SMgenBR(year, process)
          if process == "ZH":
            genBR *= filterefficienciesZH[coupling] / filterefficiencyZHpowheg
          if process == "ttH":
            genBR *= filterefficienciesttH[coupling] / filterefficiencyttHpowheg

          regex = "(GENBR=)[0-9.eE+-]+(;)"
          match = re.search(regex, row["::variables"])
          assert match, row["::variables"]
          row["::variables"] = re.sub(regex, r"\g<1>{:g}\g<2>".format(nominal_value(genBR)), row["::variables"])

          print "{:4} {:13} {:10.3g} {:10.3g}".format(process, coupling, genxsec, genBR)

        writer.writerow(row)

    with open(filename, "w") as finalf, open(newf.name) as newf2:
      finalf.write(newf2.read().replace('\r\n', '\n'))


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("year", type=int)
  args = parser.parse_args()
  update_spreadsheet(args.year)

