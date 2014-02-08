file = open('../prod/analyzer_2011.py', 'r')

b = "'(,"

print "#################################################################"
print "# DATA: the string \"data\" will trigger processing of data trees #"
print "#################################################################"
print "data"
print
print "######"
print "# MC #"
print "######"

for line in file:
   if "### " in line and "Double" not in line and "MuEG" not in line:
      name = line.split()[1]
      print
      print "#", name

   if not line.startswith('#'):
      if "('" in line and "Double" not in line and "MuEG" not in line:
         sample = line.split()[0]
         for i in range(0, len(b)):
            sample = sample.replace(b[i], "")
         print "ZZ4lAnalysis_"+sample