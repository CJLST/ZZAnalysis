#!/usr/bin/env python
# Spring2016MC_M17PUscenario taken from https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_20_patchX/SimGeneral/MixingModule/python/mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi.py#L25
# To obtain DataPileup root histogram run pileupCalc.py:

# For 2016 data (Moriond 2017 setup:)
# pileupCalc.py -i JSONfile.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 75 --numPileupBins 75 DataPileupHistogram2016.root

#Up and down variation histograms are created varying the minBiasXsec by 4.6% (72383, 66017)

# For 2017 data:
#pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-305636_13TeV_PromptReco_Collisions17_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 75 --numPileupBins 75 DataPileupHistogram2017_69200_75bins.root
#pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-305636_13TeV_PromptReco_Collisions17_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 66017 --maxPileupBin 75 --numPileupBins 75 DataPileupHistogram2017_66017_75bins.root
#pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-305636_13TeV_PromptReco_Collisions17_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 72383 --maxPileupBin 75 --numPileupBins 75 DataPileupHistogram2017_72383_75bins.root

# Afetrwards run this script to produce root file which contains PU weights


import ROOT as rt
puMC = {
    'Spring2016MC_PUscenarioV1' : [ 0.000829312873542, 0.00124276120498, 0.00339329181587, 0.00408224735376, 0.00383036590008,
                                    0.00659159288946,  0.00816022734493, 0.00943640833116, 0.0137777376066,  0.017059392038,
                                    0.0213193035468,   0.0247343174676,  0.0280848773878,  0.0323308476564,  0.0370394341409,
                                    0.0456917721191,   0.0558762890594,  0.0576956187107,  0.0625325287017,  0.0591603758776,
                                    0.0656650815128,   0.0678329011676,  0.0625142146389,  0.0548068448797,  0.0503893295063,
                                    0.040209818868,    0.0374446988111,  0.0299661572042,  0.0272024759921,  0.0219328403791,
                                    0.0179586571619,   0.0142926728247,  0.00839941654725, 0.00522366397213, 0.00224457976761,
                                    0.000779274977993, 0.000197066585944,7.16031761328e-05,0.0             , 0.0,
                                    0.0,        0.0,        0.0,        0.0,        0.0,
                                    0.0,        0.0,        0.0,        0.0,        0.0],

    'Spring2016MC_M17PUscenario' : [1.78653e-05 ,2.56602e-05 ,5.27857e-05 ,8.88954e-05 ,0.000109362 ,0.000140973 ,0.000240998 ,
                                    0.00071209 , 0.00130121 , 0.00245255 , 0.00502589 , 0.00919534 , 0.0146697 ,  0.0204126 ,
                                    0.0267586 ,  0.0337697 ,  0.0401478 ,  0.0450159 ,  0.0490577 ,  0.0524855 ,  0.0548159 ,
                                    0.0559937 ,  0.0554468 ,  0.0537687 ,  0.0512055 ,  0.0476713 ,  0.0435312 ,  0.0393107 ,
                                    0.0349812 ,  0.0307413 ,  0.0272425 ,  0.0237115 ,  0.0208329 ,  0.0182459 ,  0.0160712 ,
                                    0.0142498 ,  0.012804 ,   0.011571 ,   0.010547 ,   0.00959489 , 0.00891718 , 0.00829292 ,
                                    0.0076195 ,  0.0069806 ,  0.0062025 ,  0.00546581 , 0.00484127 , 0.00407168 , 0.00337681 ,
                                    0.00269893 , 0.00212473 , 0.00160208 , 0.00117884 , 0.000859662 ,0.000569085 ,0.000365431 ,
                                    0.000243565 ,0.00015688 , 9.88128e-05 ,6.53783e-05 ,3.73924e-05 ,2.61382e-05 ,2.0307e-05 ,
                                    1.73032e-05 ,1.435e-05 ,  1.36486e-05 ,1.35555e-05 ,1.37491e-05 ,1.34255e-05 ,1.33987e-05 ,
                                    1.34061e-05 ,1.34211e-05 ,1.34177e-05 ,1.32959e-05 ,1.33287e-05],
}

### MC pu scenario to be used
puMCscenario = puMC['Spring2016MC_M17PUscenario']
len_mc = len(puMCscenario)

#data_file_name = 'DataPileupHistogram69200_75bins.root'
#data_file_name_varUp = 'DataPileupHistogram72383_75bins.root'
#data_file_name_varDn = 'DataPileupHistogram66017_75bins.root'

data_file_name       = 'DataPileupHistogram2017_69200_75bins.root'
data_file_name_varUp = 'DataPileupHistogram2017_72383_75bins.root'
data_file_name_varDn = 'DataPileupHistogram2017_66017_75bins.root'

rt.TH1.SetDefaultSumw2(True)

h_d = rt.TH1F('Data', '', len_mc , 0, len_mc) 
h_d_varUp = rt.TH1F('Data_varUp', '', len_mc , 0, len_mc) 
h_d_varDn = rt.TH1F('Data_varDn', '', len_mc , 0, len_mc) 


fpu = rt.TFile.Open(data_file_name,'read')
h_din = fpu.Get('pileup')
for i in range(1, len_mc + 1) :
    h_d.SetBinContent(i, h_din.GetBinContent(i))
h_d.Scale(1./h_d.Integral())
fpu.Close()

fpu = rt.TFile.Open(data_file_name_varUp,'read')
h_din = fpu.Get('pileup')
for i in range(1, len_mc + 1) :
    h_d_varUp.SetBinContent(i, h_din.GetBinContent(i))
h_d_varUp.Scale(1./h_d_varUp.Integral())
fpu.Close()

fpu = rt.TFile.Open(data_file_name_varDn,'read')
h_din = fpu.Get('pileup')
for i in range(1, len_mc + 1) :
    h_d_varDn.SetBinContent(i, h_din.GetBinContent(i))
h_d_varDn.Scale(1./h_d_varDn.Integral())
fpu.Close()


h_mc = rt.TH1F('MC out-of-the-box', ';true number of interactions;normalized to unity', len_mc , 0, len_mc)
for ipu in range(len(puMCscenario)) :
    puMCscenario[ipu]
    h_mc.SetBinContent(ipu + 1, puMCscenario[ipu])

h_mc.Scale(1./h_mc.Integral())

h_w = h_d.Clone('weights')
h_w.Divide(h_mc)


h_w_varUp = h_d_varUp.Clone('weights_varUp')
h_w_varUp.Divide(h_mc)

h_w_varDn = h_d_varDn.Clone('weights_varDn')
h_w_varDn.Divide(h_mc)

h_mc_rw = h_mc.Clone('MC reweighted')
h_mc_rw_varUp = h_mc.Clone('MC_up')
h_mc_rw_varUp.SetTitle("MC reweighted +1#sigma")
h_mc_rw_varDn = h_mc.Clone('MC reweighted')
h_mc_rw_varDn.SetTitle("MC reweighted -1#sigma")


for i in range(1, len_mc + 1) :
    h_mc_rw.SetBinContent(i, h_mc.GetBinContent(i)*h_w.GetBinContent(i))
    h_mc_rw_varUp.SetBinContent(i, h_mc.GetBinContent(i)*h_w_varUp.GetBinContent(i))
    h_mc_rw_varDn.SetBinContent(i, h_mc.GetBinContent(i)*h_w_varDn.GetBinContent(i))

can = rt.TCanvas('can', 'can', 400, 400)
h_mc.SetLineColor(rt.kBlue)
h_mc.Draw()
#h_mc.GetYaxis().SetRangeUser
h_mc.SetMaximum(0.12)
h_mc.Draw()
h_d.SetLineColor(rt.kRed-2)
h_d.SetFillColor(rt.kRed-2)
h_d.SetFillStyle(3004)
h_d.Draw('HISTSAME')
h_mc_rw.SetLineColor(rt.kGreen)
h_mc_rw.Draw('SAME')

h_mc_rw_varUp.SetLineColor(rt.kGreen + 2)
h_mc_rw_varUp.SetLineStyle(2)
h_mc_rw_varUp.Draw('SAME')

h_mc_rw_varDn.SetLineColor(rt.kGreen + 2)
h_mc_rw_varDn.SetLineStyle(3)
h_mc_rw_varDn.Draw('SAME')

leg = can.BuildLegend()
leg.Draw('SAME')

f_out = rt.TFile.Open('pu_weights.root', 'recreate')
h_w.Write()
h_w_varDn.Write()
h_w_varUp.Write()

h_mc.Write()
h_d.Write()
h_mc_rw.Write()
can.Write()
f_out.Close()
