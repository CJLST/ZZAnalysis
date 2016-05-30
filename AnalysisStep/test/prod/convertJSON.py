#!/usr/bin/python

filename = 'json_2016_Silver.txt'

json = eval( open(filename).read() )

print "# from", filename
print 'process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange( *('
for run,lumisections in sorted(json.iteritems()):
 for ls in lumisections:
   print "    '%s:%s-%s:%s'," % (run, ls[0], run, ls[1])
print '))'


