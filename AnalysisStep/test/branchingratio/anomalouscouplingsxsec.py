#!/usr/bin/env python

import argparse, csv, itertools, os, re, subprocess, tempfile

from utilities import cache, cd

if not os.path.exists("anomalouscouplingsconstants"):
  subprocess.check_call(["git", "clone", "git@github.com:hroskes/anomalouscouplingsconstants"])

commit = "f7b1311913ced69cfdfdc4ebae6179674872b63b"

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
  with open("../prod/Couplings_ACJHUGen.lst") as f:
    for line in f:
      split = line.split("|")
      if split[0] == process+coupling+"_M125":
        couplingsfromlst = split[1].split(";")
        break
    else:
      raise ValueError("Didn't find "+process+coupling+"_M125 in ../prod/Couplings_ACJHUGen.lst")

  if process == "VBFH": process = "VBF"

  couplingsbrokendown = gglist, ttlist, VVlist = [[], [], []]
  for coupling in couplingsfromlst:
    couplingname, couplingval = coupling.split("=")
    couplingval = complex(*(float(_) for _ in couplingval.split(",")))
    assert couplingval.imag == 0
    couplingval = couplingval.real
    coupling = couplingname, couplingval
    if couplingname.startswith("ghz"): VVlist.append(coupling)
    elif couplingname.startswith("ghg"): gglist.append(coupling)
    elif couplingname.startswith("kappa"): ttlist.append(coupling)
    else: assert 0, coupling

  decaypart = 0
  for (coupling1name, coupling1val), (coupling2name, coupling2val) in itertools.combinations_with_replacement(VVlist, 2):
    if coupling1name == coupling2name:
      xs = coupling1name
    else:
      xs = coupling1name + coupling2name
    xs = xs.replace("ghz1_prime2", "L1").replace("ghzgs1_prime2", "L1Zg").replace("ghz1", "a1").replace("ghz2", "a2").replace("ghz4", "a3")
    decaypart += coupling1val * coupling2val * getattr(anomalouscouplingsconstants, "JHUXSHZZ4l"+xs)
  decaypart /= anomalouscouplingsconstants.JHUXSHZZ4la1

  if process == "ggH":
    if year == 2016:
      prodpart = 0.9901
    else:
      prodpart = 1
  else:
    pdf = {2016: "NNPDF30_lo_as_0130", 2017: "NNPDF31_lo_as_0130", 2018: "NNPDF31_lo_as_0130"}[year]
    pdfmodule = getattr(anomalouscouplingsconstants, pdf)

    if process == "ttH":
      lst = ttlist
      assert not gglist
    elif process == "HJJ":
      lst = gglist
      assert not ttlist
    else:
      lst = VVlist
      assert not gglist and not ttlist

    prodpart = 0
    for (coupling1name, coupling1val), (coupling2name, coupling2val) in itertools.combinations_with_replacement(lst, 2):
      if coupling1name == coupling2name:
        xs = coupling1name
      else:
        xs = coupling1name + coupling2name
      xs = xs.replace("ghz1_prime2", "L1").replace("ghzgs1_prime2", "L1Zg").replace("ghz1", "a1").replace("ghz2", "a2").replace("ghz4", "a3").replace("ghg2", "a2").replace("ghg4", "a3").replace("kappa_tilde", "kappatilde")
      prodpart += coupling1val * coupling2val * getattr(pdfmodule, "JHUXS"+process+xs)

    if process == "ttH":
      prodpart /= getattr(pdfmodule, "JHUXS"+process+"kappa")
    elif process == "HJJ":
      prodpart /= getattr(pdfmodule, "JHUXS"+process+"a2")
    else:
      prodpart /= getattr(pdfmodule, "JHUXS"+process+"a1")

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

