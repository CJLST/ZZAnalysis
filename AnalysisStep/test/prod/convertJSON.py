#!/usr/bin/python

# Running: set the correct flename (below) and run:
# ./convertJSON.py > pyFragments/json_2018.py 


#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-273730_13TeV_PromptReco_Collisions16_JSON.txt'
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-275783_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt' #Jul 1, 2016; 5.76/fb; https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2675.html
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276097_13TeV_PromptReco_Collisions16_JSON_NoL1T_v2.txt' #Jul 9, 2016; 7.65/fb; https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2681.html
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276384_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt' # Jul 15, 2016; 9.235/fb, https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2684/1.html
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt' #Jul 22, 2016; 12.9/fb; https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2689.html
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-283685_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt' #Oct 28, 2016; 33.59/fb; https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2741.html
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt' #05 Dec, 2016; 36.8/fb;

#filename ='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/DCSOnly/json_DCSONLY.txt' # DCS-only, updated regularly
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-297723_13TeV_PromptReco_Collisions17_JSON.txt' # Jul 11, 2017; 3.785/fb; https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2863.html
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-299042_13TeV_PromptReco_Collisions17_JSON.txt' # Jul 21, 2017; 3.912/fb; https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2873.html
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-301141_13TeV_PromptReco_Collisions17_JSON.txt' # Aug 23, 2017; 10.090/fb;  
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-301567_13TeV_PromptReco_Collisions17_JSON.txt' # Sep 1, 2017; 12.136/fb;
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-301997_13TeV_PromptReco_Collisions17_JSON.txt' # Sep 9, 2017; 13.88/fb;
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-305185_13TeV_PromptReco_Collisions17_JSON.txt' # Oct 26, 2017; 29.02/fb; https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2949.html
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-305364_13TeV_PromptReco_Collisions17_JSON.txt' # Nov 2, 2017; 32.501/fb;
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-305636_13TeV_PromptReco_Collisions17_JSON.txt' # Nov 10, 2017; 35.88/fb; https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2963.html
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt' # Dec 15, 2017; 41.86/fb; https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2989/1/1/1.html
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt' # Feb 02, 2018; 41.37 fb; https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/3009.html

#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-317591_13TeV_PromptReco_Collisions18_JSON.txt' #June 22, 2018
filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-324209_13TeV_PromptReco_Collisions18_JSON.txt' #Oct 19, 2018; 50.98/fb; https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/3226.html

json = eval( open(filename).read() )

print "# from", filename
print 'process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange( *('
for run,lumisections in sorted(json.iteritems()):
# if int(run) in goodPhysics :
    for ls in lumisections:
        print "    '%s:%s-%s:%s'," % (run, ls[0], run, ls[1])
print '))'


