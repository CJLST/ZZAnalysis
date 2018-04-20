#!/usr/bin/env python

import argparse, contextlib, csv, itertools, numpy, os, re, tempfile, urllib
from utilities import CJLSTproduction, TFile

def averageNLOweight(identifier, year):
  if year == 2016: return 1
  assert year == 2017

  production = CJLSTproduction[year]

  if not os.path.exists("/data3/Higgs"): raise RuntimeError("Have to run this on lxcms03")

  with TFile("/data3/Higgs/"+production+"/"+identifier+"/ZZ4lAnalysis.root") as f:
    if not f: return float("nan")
    try:
      t = f.ZZTree.Get("candTree")
      failedt = f.ZZTree.Get("candTree_failed")
    except AttributeError:
      t = failedt = None
    if not failedt: failedt = []
    if not t and not failedt: return float("nan")
    for _ in t, failedt:
      if not _: continue
      _.SetBranchStatus("*", 0)
      _.SetBranchStatus("genHEPMCweight", 1)
      _.SetBranchStatus("genHEPMCweight_NNLO", 1)

    bothtrees = itertools.chain(t, failedt)

    return numpy.average([entry.genHEPMCweight / entry.genHEPMCweight_NNLO for entry in bothtrees])

def rawgenxsec(identifier, year):
  with contextlib.closing(urllib.urlopen("https://raw.githubusercontent.com/CJLST/ZZAnalysis/0c7f32434c0c549d15955936570e16e29e95855d/AnalysisStep/test/prod/samples_{}_MC.csv".format(year))) as f:
    reader = csv.DictReader(f)
    for row in reader:
      if re.match("#*"+identifier+"$", row["identifier"]):
        variables = {_.split("=")[0]: _.split("=")[1] for _ in row["::variables"].split(";")}
        return float(variables["GENXSEC"])
  assert False, (identifier, year)

def update_spreadsheet(year):
  filename =  os.path.join(os.path.dirname(__file__), "../prod/samples_{}_MC.csv".format(year))
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
        if identifier.strip() and identifier.strip() != "See comment field":
          multiply = averageNLOweight(identifier, year)
          print identifier, multiply
          if multiply != 1 and multiply == multiply:
            genxsec = rawgenxsec(identifier, year) * multiply
            regex = "(GENXSEC=)[0-9.eE+-]+(;)"
            match = re.search(regex, row["::variables"])
            assert match, row["::variables"]
            row["::variables"] = re.sub(regex, r"\g<1>{:g}\g<2>".format(genxsec), row["::variables"])
        writer.writerow(row)

    with open(filename, "w") as finalf, open(newf.name) as newf2:
      finalf.write(newf2.read().replace('\r\n', '\n'))

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("year", type=int)
  args = parser.parse_args()
  update_spreadsheet(args.year)
