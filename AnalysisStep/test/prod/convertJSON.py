#!/usr/bin/python

#json = eval( open('/afs/cern.ch/user/a/anlevin/public/moriond_2012_remove_pixel_ecal_recovered.txt').read() ) # 2012 bulk
#json = eval( open('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_201191-201191_8TeV_11Dec2012ReReco-recover_Collisions12_JSON.txt').read() ) #for Run2012C-EcalRecover_11Dec2012
#json = eval( open('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt').read() ) #for Run2012D-16Jan2013-v1

json = eval( open('Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON_v2.txt').read() ) #for Run2012D-16Jan2013-v1


print 'process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange( *('
for run,lumisections in sorted(json.iteritems()):
 for ls in lumisections:
   print "    '%s:%s-%s:%s'," % (run, ls[0], run, ls[1])
print '))'


