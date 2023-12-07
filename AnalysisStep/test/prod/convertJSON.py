#!/usr/bin/env python3

# Running: uncomment the correct filename (below) and run:
# ./convertJSON.py > pyFragments/json_201X.py

#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt' #for 2016 taken from https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2016Analysis#Re_reco_datasets_07Aug17
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt' #for 2017 taken from https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2017Analysis#13_TeV_pp_runs_ReReco
#filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt' #for 2018 taken from https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2018Analysis#Early2018Re_reco_17Sep2018_datas
filename = '../../../NanoAnalysis/test/prod/Cert_Collisions2022_355100_362760_Golden.json'

json = eval( open(filename).read() )

print("# from", filename)
print('process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange( *(')
for run,lumisections in sorted(json.items()):
# if int(run) in goodPhysics :
    for ls in lumisections:
        print("    '%s:%s-%s:%s'," % (run, ls[0], run, ls[1]))
print('))')


