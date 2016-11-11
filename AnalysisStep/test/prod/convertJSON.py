#!/usr/bin/python

#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-273730_13TeV_PromptReco_Collisions16_JSON.txt'
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-275783_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt' #Jul 1, 2016; 5.76/fb; https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2675.html
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276097_13TeV_PromptReco_Collisions16_JSON_NoL1T_v2.txt' #Jul 9, 2016; 7.65/fb; https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2681.html
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276384_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt' # Jul 15, 2016; 9.235/fb, https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2684/1.html
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt' #Jul 22, 2016; 12.9/fb; https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2689.html
filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-283685_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt' #Oct 28, 2016; 33.59/fb; https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2741.html

json = eval( open(filename).read() )

print "# from", filename
print 'process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange( *('
for run,lumisections in sorted(json.iteritems()):
 for ls in lumisections:
   print "    '%s:%s-%s:%s'," % (run, ls[0], run, ls[1])
print '))'


