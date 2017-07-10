#!/usr/bin/python

# Running: set the correct flename (below) and run:
# ./convertJSON.py > pyFragments/json_2017.py 


#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-273730_13TeV_PromptReco_Collisions16_JSON.txt'
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-275783_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt' #Jul 1, 2016; 5.76/fb; https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2675.html
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276097_13TeV_PromptReco_Collisions16_JSON_NoL1T_v2.txt' #Jul 9, 2016; 7.65/fb; https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2681.html
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276384_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt' # Jul 15, 2016; 9.235/fb, https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2684/1.html
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt' #Jul 22, 2016; 12.9/fb; https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2689.html
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-283685_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt' #Oct 28, 2016; 33.59/fb; https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2741.html
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt' #05 Dec, 2016; 36.8/fb;

filename ='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/DCSOnly/json_DCSONLY.txt' # DCS-only, updated regularly


goodPhysics = [297047,297050,297056,297057,297099,297100,297101,297113,297114,297168,297169,297170,297171,297175,297176,297177,297178,297179,297180,297181,297211,297215,297218,297224,297225,297227,297281,297286,297292,297293,297296,297308,297046,297219,297359,297411,297424,297425,297426,297429,297430,297431,297432,297433,297434,297435,297467,297468,297469,297474,297483,297484,297485,297486,297488,297503,297504,297505,297557,297558,297559,297560,297562,297563,297598,297603,297606,297620,297659,297660,297664,297665,297666,297670,297671,297672,297674,297675,297678] # FIXME: temporary list from https://indico.cern.ch/event/651440/contributions/2651026/attachments/1489063/2313823/17-07-06_PPD_News.pdf;


json = eval( open(filename).read() )

print "# from", filename
print 'process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange( *('
for run,lumisections in sorted(json.iteritems()):

    if int(run) in goodPhysics :
        for ls in lumisections:
            print "    '%s:%s-%s:%s'," % (run, ls[0], run, ls[1])
print '))'


