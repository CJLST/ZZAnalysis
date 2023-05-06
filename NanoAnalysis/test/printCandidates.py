#!/usr/bin/env python
##
# Print a list of selected candidates from nanoAOD trees.
##
from __future__ import print_function
import math
import argparse
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

parser = argparse.ArgumentParser(description='Print selected 4l events in the specified nanoAOD files')

parser.add_argument('file', help='filename to read', nargs='?', default="ggH125_fixedFSR_Skim.root")
parser.add_argument('-r', dest ='region', help='which region to print out: "SR"=signal (default), or CRs: "2P2F", "3P1F", "SS", "SIP"', default="SR")
parser.add_argument('-f', dest ='outputFormat', help='format: "full", "eventID"', default="full")

args = parser.parse_args()

nanoFile = args.file
region = args.region
outputFormat = args.outputFormat
checkTrigger = True

if region != 'SR' and region != 'SS' and region != '2P2F' and region != '3P1F' and region != 'SIP' :
    print ('Region', region, 'not supported, must be"SR"=signal (default), or CRs: "2P2F", "3P1F", "SS", "SIP"')
    exit(1)

t = ROOT.TChain("Events")
t.Add(nanoFile)

iEntry = 0
n4e=0
n4mu=0
n2e2mu=0
n2mu2e=0

isMC = False
if (t.GetBranch("overallEventWeight")) :
    isMC = True

while t.GetEntry(iEntry):
    iEntry+=1

    eventId='{}:{}:{}'.format(t.run,t.luminosityBlock,t.event)

    flavZ1 = 0
    flavZ2 = 0

    if region == 'SR' :
        iZZ = t.bestCandIdx
        if iZZ < 0 : continue # no candidate passes selection in this event, or
            # the event does not pass the required triggers (for samples processed
            # with TRIGPASSTHROUGH=True)
        if checkTrigger and not t.HLT_passZZ4l : continue 
        if outputFormat == "eventID" :
            print(eventId)
        else : 
            flavZ1 = t.ZZCand_Z1flav[iZZ]
            flavZ2 = t.ZZCand_Z2flav[iZZ]
            print('{}:{:.2f}:{:.2f}:{:.2f}:{:.3f}:{:.3f}:{:.3f}'.format(eventId,
                                                                        t.ZZCand_mass[iZZ],
                                                                        t.ZZCand_Z1mass[iZZ],
                                                                        t.ZZCand_Z2mass[iZZ],
                                                                        t.ZZCand_KD[iZZ],
                                                                        t.overallEventWeight if isMC else 1.,
                                                                        t.ZZCand_dataMCWeight[iZZ] if isMC else 1.,
                                                                        ))
    else:
        iZLL = -1
        if   region == 'SS':   iZLL = t.ZLLbestSSIdx
        elif region == '2P2F': iZLL = t.ZLLbest2P2FIdx
        elif region == '3P1F': iZLL = t.ZLLbest3P1FIdx
        elif region == 'SIP':  iZLL = t.ZLLbestSIPCRIdx
        if iZLL < 0 : continue
        if outputFormat == "eventID" :
            print(eventId)
        else : 
            flavZ1 = t.ZLLCand_Z1flav[iZLL]
            flavZ2 = t.ZLLCand_Z2flav[iZLL]
            print('{}:{:.2f}:{:.2f}:{:.2f}:{:.3f}:{:.3f}'.format(eventId,
                                                                 t.ZLLCand_mass[iZLL],
                                                                 t.ZLLCand_Z1mass[iZLL],
                                                                 t.ZLLCand_Z2mass[iZLL],
                                                                 t.ZLLCand_KD[iZLL],
                                                                 t.overallEventWeight if isMC else 1.,
                                                                 # t.ZLLCand_dataMCWeight[iZLL],
                                                                 ))

    flav = abs(flavZ1*flavZ2)
    if flav == 11**4 :
        n4e = n4e+1
    elif flav == 13**4 :
        n4mu = n4mu+1
    elif abs(flavZ1) == (11*11) :        
        n2e2mu = n2e2mu+1
    elif abs(flavZ1) == (13*13) :        
        n2mu2e = n2mu2e+1




if region == 'SR' :    
    print ('\n## Selected events all/4e/4mu/2e2mu : ', n4e+n4mu+n2e2mu+n2mu2e, "/", n4e, "/", n4mu, "/", n2e2mu+n2mu2e, sep="")
else:
    print ('\n## ', region, " CR : ", n4e+n4mu+n2e2mu+n2mu2e, sep='')


    
