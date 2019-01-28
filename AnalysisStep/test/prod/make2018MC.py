#!/usr/bin/env python

import argparse, cPickle, csv, datetime, functools, itertools, logging, os, re

import Utilities.General.cmssw_das_client as das_client

logging.getLogger().setLevel(logging.INFO)

def error(self, message, *args, **kws):
  if self.isEnabledFor(logging.ERROR):
    self._log(logging.ERROR, message, args, **kws)
    raise RuntimeError(message)

logging.Logger.error = error

p = argparse.ArgumentParser()
p.add_argument("--clear-cache", action="store_true")
args = p.parse_args()

if args.clear_cache:
  try:
    os.remove("dasqueries.pkl")
  except OSError:
    pass

def cache_file(filename, *argkeys, **kwargkeys):
  def inner_cache_file(function):
    try:
      with open(filename, "rb") as f:
        cache = cPickle.load(f)
    except IOError:
      cache = {}
    @functools.wraps(function)
    def newfunction(*args, **kwargs):
      argsforcache = tuple(key(arg) for key, arg in itertools.izip_longest(argkeys, args, fillvalue=lambda x: x))
      kwargsforcache = {kw: kwargkeys.get(kw, lambda x: x)(kwargs[kw]) for kw in kwargs}
      keyforcache = argsforcache, tuple(sorted(kwargsforcache.iteritems()))
      try:
        return cache[keyforcache]
      except TypeError:
        print keyforcache
        raise
      except KeyError:
        result = function(*args, **kwargs)
        try:
          with open(filename, "rb") as f:
            cache.update(cPickle.load(f))
        except IOError:
          pass
        cache[keyforcache] = result
        with open(filename, "wb") as f:
          cPickle.dump(cache, f)
      return newfunction(*args, **kwargs)
    return newfunction
  return inner_cache_file

@cache_file("dasqueries.pkl")
def das(query):
  data = das_client.get_data(query, 0)
  if data["status"] == "error": raise RuntimeError(data["reason"])
  return data["data"]

for filename in "samples_{year}_MC.csv", "samples_{year}_MC_anomalous.csv":
  with open(filename.format(year=2018), "w") as newf, open(filename.format(year=2017)) as f, open(filename.format(year=2017)) as f2:
    reader = csv.DictReader(f)
    next(f2)
    writer = csv.DictWriter(newf, fieldnames=reader.fieldnames)
    writer.writeheader()
    seenPDs = set()

    for row, line in itertools.izip(reader, f2):
      while not line.strip():
        newf.write("\n")
        line = next(f2)

      if row["dataset"]:
        row["::variables"] = row["::variables"].replace("LEPTON_SETUP=2017", "LEPTON_SETUP=2018")
        assert "LEPTON_SETUP=2018" in row["::variables"]

        row["identifier"] = row["identifier"].lstrip("#")
        identifier = row["identifier"].replace("ext", "")
        datasetparts = list(row["dataset"].split("/"))
        datasetparts[2] = "RunIIAutumn18*"

        if "PSweights" in datasetparts[1]: continue

        logging.info("Found PD: "+datasetparts[1])
        datasetparts[1] = datasetparts[1].replace("JHUgenV700", "JHUGenV7011").replace("JHUgenV698", "JHUGenV7011")
        if datasetparts[1].startswith("ZH"): datasetparts[1] = datasetparts[1].replace("HWJ", "HZJ").replace("powheg-", "powheg2-")
        if datasetparts[1].startswith("ZZ"): datasetparts[1] = datasetparts[1].replace("13TeV", "TuneCP5_13TeV")

        if datasetparts[1] in seenPDs:
          logging.info("Already dealt with that one")
          logging.info("\n")
          continue


        seenPDs.add(datasetparts[1])
        dasoutput = das("dataset dataset="+("/".join(datasetparts)) + " status=*")
        if not dasoutput:
          isok = False

          if datetime.date.today() < datetime.date(2019, 2, 1):
            #GEN approved but not submitted yet
            match = re.match("GluGluHToZZTo4L_M([0-9]+)_13TeV_powheg2_JHUGenV7011_pythia8", datasetparts[1])
            if match and int(match.group(1)) >= 145: isok = True

            #Qianying is working on it
            if "NNLOPS" in datasetparts[1]: isok = True

            #these are not in DAS yet but should be soon
            if re.match("GluGluHToZZTo4L_M300_13TeV(?:_tune(?:up|down))?_powheg2_minloHJJ_JHUGenV7011_pythia8", datasetparts[1]): isok = True
            if datasetparts[1] in (
              "ttH_HToZZ_4LFilter_M145_13TeV_powheg2_JHUGenV7011_pythia8",
              "ttH_HToZZ_4LFilter_M125_13TeV_tunedown_powheg2_JHUGenV7011_pythia8",
              "WminusH_HToZZTo4L_M125_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV7011_pythia8",
              "ZZJJTo4L_EWK_TuneCP5_13TeV-madgraph-pythia8",
              "ZZJJTo4L_QCD_TuneCP5_13TeV-madgraph-pythia8",
              "Higgs0PMToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8",
            ): isok = True

            match = re.match("(?:VBF|WplusH|WminusH|ZH)_HToZZ(?:To4L|_4LFilter)_M([0-9]+)_13TeV_powheg2(?:-minlo-H[ZW]J)?_JHUGenV7011_pythia8", datasetparts[1])
            if match and int(match.group(1)) >= 145: isok = True

            #these don't exist at all.  do we need them?
            if datasetparts[1] in (
              "ZZTo4L_TuneCP5_13TeV-amcatnloFXFX-pythia8",
              "ZZTo4L_TuneCP5_13TeV-sherpa",
              "WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8",
              "TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8", #CUETP8M1???
              "ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8", #CUETP8M1???
            ): isok = True

            if "phantom" in datasetparts[1]: isok = True

          (logging.warning if isok else logging.error)("Didn't find "+datasetparts[1])
          row["identifier"] = "#" + identifier
          datasetparts[2] = "?????????"
          row["dataset"] = "/".join(datasetparts)
          writer.writerow(row)
          
        for result in dasoutput:
          newdataset = row["dataset"] = result["dataset"][0]["name"]

          row["identifier"] = identifier

          newdatasetparts = newdataset.split("/")
          if "FlatPU" in newdatasetparts[2]: continue
          if "ForMUOVal" in newdatasetparts[2]: continue

          match = re.match(r"RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15(?:_(ext[0-9]+))?-v[0-9]+$", newdatasetparts[2])
          if not match:
            raise ValueError("Didn't recognize second part of the dataset:\n"+newdataset+"\nIf this is a new version of the MiniAOD, update the script.")
          if match.group(1):
            row["identifier"] += match.group(1)

          status = result["dataset"][0]["status"]
          row["identifier"] = {
            "VALID": "",
            "PRODUCTION": "#",
          }[status] + row["identifier"]

          logging.info("Found dataset: "+newdataset)
          logging.info("status: "+status)

          writer.writerow(row)

        logging.info("\n")

      else:  #comment rows
        logging.debug("writing a comment row")
        writer.writerow(row)

