#!/usr/bin/python

filename = 'Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_Silver.txt'

json = eval( open(filename).read() )

print "# from", filename
print 'process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange( *('
for run,lumisections in sorted(json.iteritems()):
 for ls in lumisections:
   print "    '%s:%s-%s:%s'," % (run, ls[0], run, ls[1])
print '))'


