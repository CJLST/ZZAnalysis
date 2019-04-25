#!/usr/bin/python

# Running: uncomment the correct filename (below) and run:
# ./convertJSON.py > pyFragments/json_201X.py

#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt' #for 2016 taken from https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2016Analysis#DATA
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt' #for 2017 taken from https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2017Analysis#13_TeV_pp_runs_ReReco
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt' #for 2018 taken from https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2018Analysis#Early2018Re_reco_17Sep2018_datas

json = eval( open(filename).read() )

print "# from", filename
print 'process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange( *('
for run,lumisections in sorted(json.iteritems()):
# if int(run) in goodPhysics :
    for ls in lumisections:
        print "    '%s:%s-%s:%s'," % (run, ls[0], run, ls[1])
print '))'


