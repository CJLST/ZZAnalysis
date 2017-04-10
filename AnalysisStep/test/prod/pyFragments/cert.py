import FWCore.PythonUtilities.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(filename = 'pyFragments/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt').getVLuminosityBlockRange()

