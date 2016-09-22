import FWCore.PythonUtilities.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(filename = 'pyFragments/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.json').getVLuminosityBlockRange()
