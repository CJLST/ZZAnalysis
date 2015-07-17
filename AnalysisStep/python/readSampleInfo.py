# Imported from:
# https://raw.githubusercontent.com/bellan/VVXAnalysis/master/Producers/python/readSampleInfo.py
#
import sys, os, commands, math


def checkBool(val):
  if val == 'False' or val == 'FALSE' or val == 'false' or val == 'NO' or val == 'no' or val == 'No':
    return False
  elif val == 'True' or val == 'TRUE' or val == 'true' or val == 'YES' or val == 'yes' or val == 'Yes':
    return True
  else:
    return val
  

def isFloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False


def readSamplesInfo(infoFilePath = 'samples_8TeV.csv', indexBy = 'identifier'):
  """
  Loads the sample information database from the given comma-separated-values
  (csv) file.
  """

  import string

  infoFile                  = open(infoFilePath, "r")
  database                  = {}
  defaults                  = {}
  header                    = None
  for line in infoFile:
    line                    = line.strip()
    if not line:              continue
    if line.startswith("#"):  continue

    data                    = line.split(",")

    # Assign info to the entry indexed by the dataset identifier
    if header:            
      if len(data) != len(header):
        raise ValueError, "Inconsistent number of columns in data '" + line + "', expected header = " + str(header)

      info                  = {}
      for (key, value) in zip(header, data):
        if key.startswith("::"):
          dictionary = {}
          if len(value.strip()) == 0:
            value           = []
          else:
            value           = map(string.strip, value.split(";"))
            for v in value:
              bareelement = map(string.strip, v.split("="))
              if len(bareelement) == 2:
                if isFloat(bareelement[1]):
                  dictionary[bareelement[0]] = float(bareelement[1])
                else:
                  dictionary[bareelement[0]] = checkBool(bareelement[1])
              else:
                dictionary[bareelement[0]] = bareelement[0]

          info[key]           = dictionary
        else:
          info[key]           = value.strip()

      index                 = info[indexBy]
      if index in database:
        raise ValueError, "Duplicate entries encountered for %d" % index
      del info[indexBy]
      database[index]       = info

    # Read header information (first line in file)
    else:
      header                = []
      for datum in data:
        if not datum:       break
        if "=" in datum:
          (datum, default)  = map(str.strip, datum.split("="))
          defaults[datum]   = checkBool(default)
        header.append(datum)


  if len(database) == 0:
    raise ValueError, "Invalid information file '" + infoFilePath + "', no entries found."
  return (database, defaults)



def readSampleInfo(sample, infoFilePath = 'samples_8TeV.csv', indexBy = 'identifier'):
  db,defaults = readSamplesInfo(infoFilePath, indexBy)

  if sample in db:
    return db[sample]
  else:
    print "Unknown sample", sample
    sys.exit(2)


def crossSection(sample, infoFilePath = 'samples_8TeV.csv', indexBy = 'identifier'):
  return float(readSampleInfo(sample, infoFilePath, indexBy)['crossSection'])

#merge together db and defaults
def readSampleDB(infoFilePath = 'samples_8TeV.csv', indexBy = 'identifier'):
  db,defaults = readSamplesInfo(infoFilePath, indexBy)
  for sample in db:
    for key,val in db[sample].iteritems():
      if key in defaults:
        if val == "":
          db[sample][key] = defaults[key]
          #print "setting default for ", key, "=",db[sample][key]

  return db


def typeOfSamples(infoFilePath = 'samples_8TeV.csv', indexBy = 'identifier'):
  db,defaults = readSamplesInfo(infoFilePath, indexBy)
  typeofsamples = []
  for sample in db:
    typeofsample = db[sample]['process']
    if not typeofsample in typeofsamples and not typeofsample == "":
      typeofsamples.append(typeofsample)
      
  return typeofsamples

def getSamplesBy(category, categoryvalue, infoFilePath = 'samples_8TeV.csv', indexBy = 'identifier'):
  DB = readSampleDB(infoFilePath, indexBy)
  samples = []
  for sample in DB:

    if(category == indexBy):
      if(sample == categoryvalue):
        samples.append(sample)
    else:
      if not category in DB[sample]: 
        print "ERROR unknown category!"
        return samples
      if DB[sample][category] == categoryvalue:
        samples.append(sample)
  
  return samples

